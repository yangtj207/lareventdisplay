////////////////////////////////////////////////////////////////////////
/// \file   DrawLArVoxel3D_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lareventdisplay/EventDisplay/Style.h"
#include "larsim/Simulation/LArVoxelData.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimListUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"

#include "TDatabasePDG.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

namespace evdb_tool {

  class DrawLArVoxel3D : public ISim3DDrawer {
  public:
    explicit DrawLArVoxel3D(const fhicl::ParameterSet&);

    void Draw(const art::Event&, evdb::View3D*) const override;

  private:
    int GetMCTruth(const art::Event&, std::vector<const simb::MCTruth*>&) const;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawLArVoxel3D::DrawLArVoxel3D(const fhicl::ParameterSet& pset) {}

  void
  DrawLArVoxel3D::Draw(const art::Event& evt, evdb::View3D* view) const
  {
    art::ServiceHandle<evd::SimulationDrawingOptions const> drawOpt;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowSimChannelInfo) return;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowMCTruthTrajectories) return;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    art::ServiceHandle<geo::Geometry const> geom;

    // Recover a handle to the collection of MCParticles
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;

    evt.getByLabel(drawOpt->fG4ModuleLabel, mcParticleHandle);

    if (!mcParticleHandle.isValid()) return;

    // Define a couple of colors for neutrals and if we gray it out...
    int neutralColor(12);
    int grayedColor(15);
    int neutrinoColor(38);

    // Use the LArVoxelList to get the true energy deposition locations as opposed to using MCTrajectories
    const sim::LArVoxelList voxels =
      sim::SimListUtils::GetLArVoxelList(evt, drawOpt->fSimChannelLabel.encode());

    mf::LogDebug("SimulationDrawer")
      << "Starting loop over " << mcParticleHandle->size() << " McParticles, voxel list size is "
      << voxels.size() << std::endl;

    // Using the voxel information can be slow (see previous implementation of this code).
    // In order to speed things up we have modified the strategy:
    // 1) Make one pass through the list of voxels
    // 2) For each voxel, keep track of the MCParticle contributing energy to it and it's position
    //    which is done by keeping a map between the MCParticle and a vector of positions
    // 3) Then loop through the map to draw the particle trajectories.
    // One caveat is the need for MCParticles... and the voxels contain the track ids. So we'll need one
    // more loop to make a map of track id's and MCParticles.

    // First up is to build the map between track id's and associated MCParticles so we can recover when looping over voxels
    std::map<int, const simb::MCParticle*> trackToMcParticleMap;

    // Should we display the trajectories too?
    double minPartEnergy(0.01);

    for (size_t p = 0; p < mcParticleHandle->size(); ++p) {
      art::Ptr<simb::MCParticle> mcParticle(mcParticleHandle, p);

      trackToMcParticleMap[mcParticle->TrackId()] = mcParticle.get();

      // Quick loop through to draw trajectories...
      if (drawOpt->fShowMCTruthTrajectories) {
        // Is there an associated McTrajectory?
        const simb::MCTrajectory& mcTraj = mcParticle->Trajectory();

        int pdgCode(mcParticle->PdgCode());
        int colorIdx(evd::Style::ColorFromPDG(mcParticle->PdgCode()));
        TParticlePDG* partPDG(TDatabasePDG::Instance()->GetParticle(pdgCode));
        double partCharge = partPDG ? partPDG->Charge() : 0.;
        double partEnergy = mcParticle->E();

        if (!drawOpt->fShowMCTruthColors) colorIdx = grayedColor;

        if (!mcTraj.empty() && partEnergy > minPartEnergy && mcParticle->TrackId() < 100000000) {
          double g4Ticks(clockData.TPCG4Time2Tick(mcParticle->T()) - trigger_offset(clockData));
          double xOffset(0.);
          double xPosMinTick = 0.;
          double xPosMaxTick = std::numeric_limits<double>::max();

          // collect the points from this particle
          int numTrajPoints = mcTraj.size();

          std::unique_ptr<double[]> hitPositions(new double[3 * numTrajPoints]);
          int hitCount(0);

          for (int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++) {
            double xPos = mcTraj.X(hitIdx);
            double yPos = mcTraj.Y(hitIdx);
            double zPos = mcTraj.Z(hitIdx);

            // If we have cosmic rays then we need to get the offset which allows translating from
            // when they were generated vs when they were tracked.
            // Note that this also explicitly checks that they are in a TPC volume
            geo::Point_t hitLocation(xPos, yPos, zPos);

            try {
              geo::TPCID tpcID = geom->PositionToTPCID(hitLocation);
              geo::PlaneID planeID(tpcID, 0);

              xPosMinTick = detProp.ConvertTicksToX(0, planeID);
              xPosMaxTick = detProp.ConvertTicksToX(detProp.NumberTimeSamples(), planeID);
              xOffset = detProp.ConvertTicksToX(g4Ticks, planeID) - xPosMinTick;

              if (xPosMaxTick < xPosMinTick) std::swap(xPosMinTick, xPosMaxTick);
            }
            catch (...) {
              continue;
            }

            // Now move the hit position to correspond to the timing
            xPos += xOffset;

            // Check fiducial limits
            if (xPos > xPosMinTick && xPos < xPosMaxTick) {
              hitPositions[3 * hitCount] = xPos;
              hitPositions[3 * hitCount + 1] = yPos;
              hitPositions[3 * hitCount + 2] = zPos;
              hitCount++;
            }
          }

          TPolyLine3D& pl(view->AddPolyLine3D(1, colorIdx, 1, 1));

          // Draw neutrals as a gray dotted line to help fade into background a bit...
          if (partCharge == 0.) {
            pl.SetLineColor(neutralColor);
            pl.SetLineStyle(3);
            pl.SetLineWidth(1);
          }
          pl.SetPolyLine(hitCount, hitPositions.get(), "");
        }
      }
    }

    // Now we set up and build the map between MCParticles and a vector of positions obtained from the voxels
    std::map<const simb::MCParticle*, std::vector<std::vector<double>>> partToPosMap;

    sim::LArVoxelList::const_iterator vxitr;
    for (vxitr = voxels.begin(); vxitr != voxels.end(); vxitr++) {
      const sim::LArVoxelData& vxd = (*vxitr).second;

      for (size_t partIdx = 0; partIdx < vxd.NumberParticles(); partIdx++) {
        if (vxd.Energy(partIdx) > drawOpt->fMinEnergyDeposition) {
          int trackId = vxd.TrackID(partIdx);

          // It can be in some instances that mcPart here could be zero.
          const simb::MCParticle* mcPart = trackToMcParticleMap[trackId];

          partToPosMap[mcPart].push_back(std::vector<double>(3));

          partToPosMap[mcPart].back()[0] = vxd.VoxelID().X();
          partToPosMap[mcPart].back()[1] = vxd.VoxelID().Y();
          partToPosMap[mcPart].back()[2] = vxd.VoxelID().Z();
        }
      } // end if this track id is in the current voxel
    }   // end loop over voxels

    // Finally ready for the main event! Simply loop through the map between MCParticle and positions to
    // draw the trajectories
    std::map<const simb::MCParticle*, std::vector<std::vector<double>>>::iterator partToPosMapItr;

    for (partToPosMapItr = partToPosMap.begin(); partToPosMapItr != partToPosMap.end();
         partToPosMapItr++) {
      // Recover the McParticle, we'll need to access several data members so may as well dereference it
      const simb::MCParticle* mcPart = partToPosMapItr->first;

      // Apparently, it can happen that we get a null pointer here or maybe no points to plot
      if (!mcPart || partToPosMapItr->second.empty()) continue;

      double g4Ticks(clockData.TPCG4Time2Tick(mcPart->T()) - trigger_offset(clockData));
      double xOffset = 0.;
      double xPosMinTick = 0.;
      double xPosMaxTick = std::numeric_limits<double>::max();

      int colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
      int markerIdx(kFullDotSmall);
      int markerSize(2);

      if (!drawOpt->fShowMCTruthFullSize) {
        colorIdx = grayedColor;
        markerIdx = kDot;
        markerSize = 1;
      }

      std::unique_ptr<double[]> hitPositions(new double[3 * partToPosMapItr->second.size()]);
      int hitCount(0);

      // Now loop over points and add to trajectory
      for (size_t posIdx = 0; posIdx < partToPosMapItr->second.size(); posIdx++) {
        const std::vector<double>& posVec = partToPosMapItr->second[posIdx];

        // Check xOffset state and set if necessary
        geo::Point_t hitLocation(posVec[0], posVec[1], posVec[2]);

        try {
          geo::TPCID tpcID = geom->PositionToTPCID(hitLocation);
          geo::PlaneID planeID(tpcID, 0);

          xPosMinTick = detProp.ConvertTicksToX(0, planeID);
          xPosMaxTick = detProp.ConvertTicksToX(detProp.NumberTimeSamples(), planeID);
          xOffset = detProp.ConvertTicksToX(g4Ticks, planeID) - xPosMinTick;

          if (xPosMaxTick < xPosMinTick) std::swap(xPosMinTick, xPosMaxTick);
        }
        catch (...) {
          continue;
        }

        double xCoord = posVec[0] + xOffset;

        // If a voxel records an energy deposit then must have been in the TPC
        // But because things get shifted still need to cut off if outside drift
        if (xCoord > xPosMinTick && xCoord < xPosMaxTick) {
          hitPositions[3 * hitCount] = xCoord;
          hitPositions[3 * hitCount + 1] = posVec[1];
          hitPositions[3 * hitCount + 2] = posVec[2];
          hitCount++;
        }
      }

      TPolyMarker3D& pm = view->AddPolyMarker3D(1, colorIdx, markerIdx, markerSize);
      pm.SetPolyMarker(hitCount, hitPositions.get(), markerIdx);
    }

    // Finally, let's see if we can draw the incoming particle from the MCTruth information
    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);

    // Loop through the MCTruth vector
    for (unsigned int idx = 0; idx < mctruth.size(); idx++) {
      // Go through each MCTruth object in the list
      for (int particleIdx = 0; particleIdx < mctruth[idx]->NParticles(); particleIdx++) {
        const simb::MCParticle& mcPart = mctruth[idx]->GetParticle(particleIdx);

        // A negative mother id indicates the "primary" particle
        if (mcPart.Mother() == -1 && mcPart.StatusCode() == 0) {
          mf::LogDebug("SimulationDrawer") << mcPart << std::endl;

          // Get position vector
          TVector3 particlePosition(mcPart.Vx(), mcPart.Vy(), mcPart.Vz());

          // Get direction vector (in opposite direction)
          TVector3 oppPartDir(-mcPart.Px(), -mcPart.Py(), -mcPart.Pz());

          if (oppPartDir.Mag2() > 0.) oppPartDir.SetMag(1.);

          double arcLenToDraw = -particlePosition.Z() / oppPartDir.CosTheta();

          // No point in drawing if arc length is zero (e.g. Ar nucleus)
          if (arcLenToDraw > 0.) {
            // Draw the line, use an off color to be unique
            TPolyLine3D& pl(view->AddPolyLine3D(2, neutrinoColor, 1, 2));

            pl.SetPoint(0, particlePosition.X(), particlePosition.Y(), particlePosition.Z());

            particlePosition += std::min(arcLenToDraw + 10., 1000.) * oppPartDir;

            pl.SetPoint(1, particlePosition.X(), particlePosition.Y(), particlePosition.Z());
          }
        }
        // The particles we want to draw will be early in the list so break out if we didn't find them
        else
          break;
      } // loop on particles in list
    }

    return;
  }

  //......................................................................
  int
  DrawLArVoxel3D::GetMCTruth(const art::Event& evt, std::vector<const simb::MCTruth*>& mcvec) const
  {
    mcvec.clear();

    if (evt.isRealData()) return 0;

    std::vector<const simb::MCTruth*> temp;

    std::vector<art::Handle<std::vector<simb::MCTruth>>> mctcol;

    // use get by Type because there should only be one collection of these in the event
    try {
      //evt.getManyByType(mctcol);
      mctcol = evt.getMany<std::vector<simb::MCTruth>>();

      for (size_t mctc = 0; mctc < mctcol.size(); ++mctc) {
        art::Handle<std::vector<simb::MCTruth>> mclistHandle = mctcol[mctc];

        for (size_t i = 0; i < mclistHandle->size(); ++i) {
          temp.push_back(&(mclistHandle->at(i)));
        }
      }
      temp.swap(mcvec);
    }
    catch (cet::exception& e) {
      mf::LogWarning("DrawLArVoxel3D") << "GetMCTruth:"
                                       << " failed with message:\n"
                                       << e;
    }

    return mcvec.size();
  }

  DEFINE_ART_CLASS_TOOL(DrawLArVoxel3D)
}
