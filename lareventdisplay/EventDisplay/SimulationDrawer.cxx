///
/// \file    SimulationDrawer.cxx
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
///

#include <algorithm>
#include <iomanip>

#include "TDatabasePDG.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "lareventdisplay/EventDisplay/Style.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/LArVoxelData.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimListUtils.h"
#include "nuevdb/EventDisplayBase/View2D.h"
#include "nuevdb/EventDisplayBase/View3D.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void
  writeErrMsg(const char* fcn, cet::exception const& e)
  {
    mf::LogWarning("SimulationDrawer") << "SimulationDrawer::" << fcn << " failed with message:\n"
                                       << e;
  }
}

namespace evd {

  SimulationDrawer::SimulationDrawer()
  {
    // For now only draw cryostat=0.
    art::ServiceHandle<geo::Geometry const> geom;
    minx = 1e9;
    maxx = -1e9;
    miny = 1e9;
    maxy = -1e9;
    minz = 1e9;
    maxz = -1e9;

    // This is looking to find the range of the complete active volume... this may not be the
    // best way to do this...
    for (size_t cryoIdx = 0; cryoIdx < geom->Ncryostats(); cryoIdx++) {
      const geo::CryostatGeo& cryoGeo = geom->Cryostat(cryoIdx);

      for (size_t tpcIdx = 0; tpcIdx < cryoGeo.NTPC(); tpcIdx++) {
        const geo::TPCGeo& tpc = cryoGeo.TPC(tpcIdx);

        mf::LogDebug("SimulationDrawer")
          << "Cryo/TPC idx: " << cryoIdx << "/" << tpcIdx << ", TPC center: " << tpc.GetCenter()[0]
          << ", " << tpc.GetCenter()[1] << ", " << tpc.GetCenter()[2] << std::endl;
        mf::LogDebug("SimulationDrawer")
          << "         TPC Active center: " << tpc.GetActiveVolumeCenter()[0] << ", "
          << tpc.GetActiveVolumeCenter()[1] << ", " << tpc.GetActiveVolumeCenter()[2]
          << ", H/W/L: " << tpc.ActiveHalfHeight() << "/" << tpc.ActiveHalfWidth() << "/"
          << tpc.ActiveLength() << std::endl;

        if (minx > tpc.GetCenter()[0] - tpc.HalfWidth())
          minx = tpc.GetCenter()[0] - tpc.HalfWidth();
        if (maxx < tpc.GetCenter()[0] + tpc.HalfWidth())
          maxx = tpc.GetCenter()[0] + tpc.HalfWidth();
        if (miny > tpc.GetCenter()[1] - tpc.HalfHeight())
          miny = tpc.GetCenter()[1] - tpc.HalfHeight();
        if (maxy < tpc.GetCenter()[1] + tpc.HalfHeight())
          maxy = tpc.GetCenter()[1] + tpc.HalfHeight();
        if (minz > tpc.GetCenter()[2] - tpc.Length() / 2.)
          minz = tpc.GetCenter()[2] - tpc.Length() / 2.;
        if (maxz < tpc.GetCenter()[2] + tpc.Length() / 2.)
          maxz = tpc.GetCenter()[2] + tpc.Length() / 2.;

        mf::LogDebug("SimulationDrawer")
          << "        minx/maxx: " << minx << "/" << maxx << ", miny/maxy: " << miny << "/" << maxy
          << ", minz/miny: " << minz << "/" << maxz << std::endl;
      }
    }
  }

  //......................................................................

  SimulationDrawer::~SimulationDrawer() {}

  //......................................................................

  void
  SimulationDrawer::MCTruthShortText(const art::Event& evt, evdb::View2D* view)
  {

    if (evt.isRealData()) return;

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;

    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);

    for (unsigned int i = 0; i < mctruth.size(); ++i) {
      std::string mctext;
      bool firstin = true;
      bool firstout = true;
      std::string origin;
      std::string incoming;
      std::string outgoing;
      // Label cosmic rays -- others are pretty obvious
      if (mctruth[i]->Origin() == simb::kCosmicRay) origin = "c-ray: ";
      int jmax = TMath::Min(20, mctruth[i]->NParticles());
      for (int j = 0; j < jmax; ++j) {
        const simb::MCParticle& p = mctruth[i]->GetParticle(j);
        char buff[1024];
        if (p.P() > 0.05) {
          sprintf(buff,
                  "#color[%d]{%s #scale[0.75]{[%.1f GeV/c]}}",
                  Style::ColorFromPDG(p.PdgCode()),
                  Style::LatexName(p.PdgCode()),
                  p.P());
        }
        else {
          sprintf(buff,
                  "#color[%d]{%s}",
                  Style::ColorFromPDG(p.PdgCode()),
                  Style::LatexName(p.PdgCode()));
        }
        if (p.StatusCode() == 0) {
          if (firstin == false) incoming += " + ";
          incoming += buff;
          firstin = false;
        }
        if (p.StatusCode() == 1) {
          if (firstout == false) outgoing += " + ";
          outgoing += buff;
          firstout = false;
        }
      } // loop on j particles
      if (origin == "" && incoming == "") { mctext = outgoing; }
      else {
        mctext = origin + incoming + " #rightarrow " + outgoing;
      }
      TLatex& latex = view->AddLatex(0.03, 0.2, mctext.c_str());
      latex.SetTextSize(0.6);

    } // loop on i mctruth objects
  }

  //......................................................................

  void
  SimulationDrawer::MCTruthLongText(const art::Event& evt, evdb::View2D* /*view*/)
  {
    if (evt.isRealData()) return;

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;

    std::vector<const simb::MCTruth*> mctruth;
    this->GetMCTruth(evt, mctruth);
    std::cout << "\nMCTruth Ptcl trackID            PDG      P      T   Moth  Process\n";
    for (unsigned int i = 0; i < mctruth.size(); ++i) {
      for (int j = 0; j < mctruth[i]->NParticles(); ++j) {
        const simb::MCParticle& p = mctruth[i]->GetParticle(j);
        if (p.StatusCode() == 0 || p.StatusCode() == 1) {
          int KE = 1000 * (p.E() - p.Mass());
          std::cout << std::right << std::setw(7) << i << std::setw(5) << j << std::setw(8)
                    << p.TrackId() << " " << std::setw(14) << Style::LatexName(p.PdgCode())
                    << std::setw(7) << int(1000 * p.P()) << std::setw(7) << KE << std::setw(7)
                    << p.Mother() << " " << p.Process() << "\n";
        }
        /*
                std::cout << Style::LatexName(p.PdgCode())
                  << "\t\t" << p.E() << " GeV"
                  << "\t"   << "(" << p.P() << " GeV/c)"
                  << std::endl;
*/
      } // loop on j particles in list
    }
    std::cout << "Note: Momentum, P, and kinetic energy, T, in MeV/c\n";
  } // MCTruthLongText

  //......................................................................
  //this is the method you would use to color code hits by the MC truth pdg value
  void
  SimulationDrawer::MCTruthVectors2D(const art::Event& evt, evdb::View2D* view, unsigned int plane)
  {
    if (evt.isRealData()) return;

    const spacecharge::SpaceCharge* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;
    // If the option is turned off, there's nothing to do
    if (drawopt->fShowMCTruthVectors > 3) {
      std::cout << "Unsupported ShowMCTruthVectors option (> 2)\n";
      return;
    }

    art::ServiceHandle<geo::Geometry const> geo;
    art::ServiceHandle<evd::RawDrawingOptions const> rawopt;
    // get the x position of the plane in question
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    // shift the truth by a fixed amount so it doesn't overlay the reco
    double xShift = -2;

    static bool first = true;
    if (first) {
      std::cout
        << "******** Show MCTruth (Genie) particles when ShowMCTruthVectors = 1 or 3 ******** \n";
      std::cout << "  MCTruth vectors corrected for space charge? " << sce->EnableCorrSCE()
                << " and shifted by " << xShift << " cm in X\n";
      std::cout << "  Neutrons and photons drawn with dotted lines. \n";
      std::cout << "  Red = e+/-, nue, nuebar. Blue = mu+/-, numu, numubar. Green = tau+/-, nutau, "
                   "nutaubar.\n";
      std::cout << "  Yellow = photons. Magenta = pions, protons and nuetrons.\n";
      std::cout << "******** Show MCParticle (Geant) decay photons (e.g. from pizeros) when "
                   "ShowMCTruthVectors = 2 or 3 ******** \n";
      std::cout << "  Photons > 50 MeV are drawn as dotted lines corrected for space charge and "
                   "are not shifted.\n";
      std::cout << "  Decay photon end points are drawn at 2 interaction lengths (44 cm) from the "
                   "start position.\n";
      std::cout << "  Color: Green = (50 < E_photon < 100 MeV), Blue = (100 MeV < E_photon < 200 "
                   "MeV), Red = (E_photon > 300 MeV).\n";
    }

    bool showTruth = (drawopt->fShowMCTruthVectors == 1 || drawopt->fShowMCTruthVectors == 3);
    bool showPhotons = (drawopt->fShowMCTruthVectors > 1);

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    if (showTruth) {
      // Unpack and draw the MC vectors
      std::vector<const simb::MCTruth*> mctruth;
      this->GetMCTruth(evt, mctruth);

      for (size_t i = 0; i < mctruth.size(); ++i) {
        if (mctruth[i]->Origin() == simb::kCosmicRay) continue;
        for (int j = 0; j < mctruth[i]->NParticles(); ++j) {
          const simb::MCParticle& p = mctruth[i]->GetParticle(j);

          // Skip all but incoming and out-going particles
          if (!(p.StatusCode() == 0 || p.StatusCode() == 1)) continue;

          double r = p.P() * 10.0; // Scale length so 10 cm = 1 GeV/c

          if (p.StatusCode() == 0) r = -r; // Flip for incoming particles

          auto gptStart = geo::Point_t(p.Vx(), p.Vy(), p.Vz());
          geo::Point_t sceOffset{0, 0, 0};
          if (sce->EnableCorrSCE()) sceOffset = sce->GetPosOffsets(gptStart);
          xyz1[0] = p.Vx() - sceOffset.X();
          xyz1[1] = p.Vy() + sceOffset.Y();
          xyz1[2] = p.Vz() + sceOffset.Z();
          xyz2[0] = xyz1[0] + r * p.Px() / p.P();
          xyz2[1] = xyz1[1] + r * p.Py() / p.P();
          xyz2[2] = xyz1[2] + r * p.Pz() / p.P();

          double w1 =
            geo->WireCoordinate(xyz1[1], xyz1[2], (int)plane, rawopt->fTPC, rawopt->fCryostat);
          double w2 =
            geo->WireCoordinate(xyz2[1], xyz2[2], (int)plane, rawopt->fTPC, rawopt->fCryostat);

          double time =
            detProp.ConvertXToTicks(xyz1[0] + xShift, (int)plane, rawopt->fTPC, rawopt->fCryostat);
          double time2 =
            detProp.ConvertXToTicks(xyz2[0] + xShift, (int)plane, rawopt->fTPC, rawopt->fCryostat);

          if (rawopt->fAxisOrientation < 1) {
            TLine& l = view->AddLine(w1, time, w2, time2);
            evd::Style::FromPDG(l, p.PdgCode());
          }
          else {
            TLine& l = view->AddLine(time, w1, time2, w2);
            evd::Style::FromPDG(l, p.PdgCode());
          }
          //

        } // loop on j particles in list
      }   // loop on truths
    }     // showTruth

    if (showPhotons) {
      // draw pizero photons with T > 30 MeV
      art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
      sim::ParticleList const& plist = pi_serv->ParticleList();
      if (plist.empty()) return;
      // photon interaction length approximate
      double r = 44;
      for (sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
        simb::MCParticle* p = (*ipart).second;
        int trackID = p->TrackId();
        art::Ptr<simb::MCTruth> theTruth = pi_serv->TrackIdToMCTruth_P(trackID);
        if (theTruth->Origin() == simb::kCosmicRay) continue;
        if (p->PdgCode() != 22) continue;
        if (p->Process() != "Decay") continue;
        int TMeV = 1000 * (p->E() - p->Mass());
        if (TMeV < 30) continue;
        auto gptStart = geo::Point_t(p->Vx(), p->Vy(), p->Vz());
        geo::Point_t sceOffset{0, 0, 0};
        if (sce->EnableCorrSCE()) sceOffset = sce->GetPosOffsets(gptStart);
        xyz1[0] = p->Vx() - sceOffset.X();
        xyz1[1] = p->Vy() + sceOffset.Y();
        xyz1[2] = p->Vz() + sceOffset.Z();
        xyz2[0] = xyz1[0] + r * p->Px() / p->P();
        xyz2[1] = xyz1[1] + r * p->Py() / p->P();
        xyz2[2] = xyz1[2] + r * p->Pz() / p->P();
        double w1 =
          geo->WireCoordinate(xyz1[1], xyz1[2], (int)plane, rawopt->fTPC, rawopt->fCryostat);
        double t1 = detProp.ConvertXToTicks(xyz1[0], (int)plane, rawopt->fTPC, rawopt->fCryostat);
        double w2 =
          geo->WireCoordinate(xyz2[1], xyz2[2], (int)plane, rawopt->fTPC, rawopt->fCryostat);
        double t2 = detProp.ConvertXToTicks(xyz2[0], (int)plane, rawopt->fTPC, rawopt->fCryostat);
        TLine& l = view->AddLine(w1, t1, w2, t2);
        l.SetLineWidth(2);
        l.SetLineStyle(kDotted);
        if (TMeV < 100) { l.SetLineColor(kGreen); }
        else if (TMeV < 200) {
          l.SetLineColor(kBlue);
        }
        else {
          l.SetLineColor(kRed);
        }
      } // ipart
    }   // showPhotons

    first = false;

  } // MCTruthVectors2D

  //......................................................................
  //this method draws the true particle trajectories in 3D
  void
  SimulationDrawer::MCTruth3D(const art::Event& evt, evdb::View3D* view)
  {
    if (evt.isRealData()) return;

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;
    // If the option is turned off, there's nothing to do
    if (!drawopt->fShowMCTruthTrajectories) return;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    art::ServiceHandle<geo::Geometry const> geom;

    // get the particles from the Geant4 step
    std::vector<const simb::MCParticle*> plist;
    this->GetParticle(evt, plist);

    // Define a couple of colors for neutrals and if we gray it out...
    int neutralColor(12);
    int grayedColor(15);
    int neutrinoColor(38);

    // Use the LArVoxelList to get the true energy deposition locations as opposed to using MCTrajectories
    const sim::LArVoxelList voxels =
      sim::SimListUtils::GetLArVoxelList(evt, drawopt->fG4ModuleLabel.label());

    mf::LogDebug("SimulationDrawer")
      << "Starting loop over " << plist.size() << " McParticles, voxel list size is "
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

    for (size_t p = 0; p < plist.size(); ++p) {
      trackToMcParticleMap[plist[p]->TrackId()] = plist[p];

      // Quick loop through to draw trajectories...
      if (drawopt->fShowMCTruthTrajectories) {
        // Is there an associated McTrajectory?
        const simb::MCParticle* mcPart = plist[p];
        const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

        int pdgCode(mcPart->PdgCode());
        int colorIdx(evd::Style::ColorFromPDG(mcPart->PdgCode()));
        TParticlePDG* partPDG(TDatabasePDG::Instance()->GetParticle(pdgCode));
        double partCharge = partPDG ? partPDG->Charge() : 0.;
        double partEnergy = mcPart->E();

        if (!drawopt->fShowMCTruthColors) colorIdx = grayedColor;

        if (!mcTraj.empty() && partEnergy > minPartEnergy && mcPart->TrackId() < 100000000) {
          double g4Ticks(clockData.TPCG4Time2Tick(mcPart->T()) - trigger_offset(clockData));
          double xOffset = 0.;
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
              // Check for space charge offsets
              //                        if (spaceCharge->EnableSimEfieldSCE())
              //                        {
              //                            std::vector<double> offsetVec = spaceCharge->GetPosOffsets(xPos,yPos,zPos);
              //                            xPos += offsetVec[0] - 0.7;
              //                            yPos -= offsetVec[1];
              //                            zPos -= offsetVec[2];
              //                        }

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
        if (vxd.Energy(partIdx) > drawopt->fMinEnergyDeposition) {
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

      if (!drawopt->fShowMCTruthFullSize) {
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
  //this method draws the true particle trajectories in 3D Ortho view.
  void
  SimulationDrawer::MCTruthOrtho(const art::Event& evt,
                                 evd::OrthoProj_t proj,
                                 double msize,
                                 evdb::View2D* view)
  {
    if (evt.isRealData()) return;

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;

    // If the option is turned off, there's nothing to do
    if (!drawopt->fShowMCTruthTrajectories) return;

    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    // get the particles from the Geant4 step
    std::vector<const simb::MCParticle*> plist;
    this->GetParticle(evt, plist);

    // Useful variables

    double xMinimum(-1. * (maxx - minx));
    double xMaximum(2. * (maxx - minx));

    // Use the LArVoxelList to get the true energy deposition locations as
    // opposed to using MCTrajectories
    // GetLarVoxelList should be adjusted to take an InputTag instead of strange.
    const sim::LArVoxelList voxels =
      sim::SimListUtils::GetLArVoxelList(evt, drawopt->fSimChannelLabel.encode());

    mf::LogDebug("SimulationDrawer")
      << "Starting loop over " << plist.size() << " McParticles, voxel list size is "
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
    bool displayMcTrajectories(true);
    double minPartEnergy(0.025);

    double tpcminx = 1.0;
    double tpcmaxx = -1.0;
    double xOffset = 0.0;
    double g4Ticks = 0.0;
    double coeff = 0.0;
    double readoutwindowsize = 0.0;
    double vtx[3] = {0.0, 0.0, 0.0};
    for (size_t p = 0; p < plist.size(); ++p) {
      trackToMcParticleMap[plist[p]->TrackId()] = plist[p];

      // Quick loop through to drawn trajectories...
      if (displayMcTrajectories) {
        // Is there an associated McTrajectory?
        const simb::MCParticle* mcPart = plist[p];
        const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

        int pdgCode(mcPart->PdgCode());
        TParticlePDG* partPDG(TDatabasePDG::Instance()->GetParticle(pdgCode));
        double partCharge = partPDG ? partPDG->Charge() : 0.;
        double partEnergy = mcPart->E();

        if (!mcTraj.empty() && partEnergy > minPartEnergy && mcPart->TrackId() < 100000000) {
          // collect the points from this particle
          int numTrajPoints = mcTraj.size();

          std::unique_ptr<double[]> hitPosX(new double[numTrajPoints]);
          std::unique_ptr<double[]> hitPosY(new double[numTrajPoints]);
          std::unique_ptr<double[]> hitPosZ(new double[numTrajPoints]);
          int hitCount(0);

          double xPos = mcTraj.X(0);
          double yPos = mcTraj.Y(0);
          double zPos = mcTraj.Z(0);

          tpcminx = 1.0;
          tpcmaxx = -1.0;
          xOffset = 0.0;
          g4Ticks = 0.0;
          vtx[0] = 0.0;
          vtx[1] = 0.0;
          vtx[2] = 0.0;
          coeff = 0.0;
          readoutwindowsize = 0.0;
          for (int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++) {
            xPos = mcTraj.X(hitIdx);
            yPos = mcTraj.Y(hitIdx);
            zPos = mcTraj.Z(hitIdx);

            // If the original simulated hit did not occur in the TPC volume then don't draw it
            if (xPos < minx || xPos > maxx || yPos < miny || yPos > maxy || zPos < minz ||
                zPos > maxz)
              continue;

            if ((xPos < tpcminx) || (xPos > tpcmaxx)) {
              vtx[0] = xPos;
              vtx[1] = yPos;
              vtx[2] = zPos;
              geo::TPCID tpcid = geom->FindTPCAtPosition(vtx);
              unsigned int cryo = geom->FindCryostatAtPosition(vtx);

              if (tpcid.isValid) {
                unsigned int tpc = tpcid.TPC;
                const geo::TPCGeo& tpcgeo = geom->GetElement(tpcid);
                tpcminx = tpcgeo.MinX();
                tpcmaxx = tpcgeo.MaxX();

                coeff = detProp.GetXTicksCoefficient(tpc, cryo);
                readoutwindowsize =
                  detProp.ConvertTicksToX(detProp.ReadOutWindowSize(), 0, tpc, cryo);

                // The following is meant to get the correct offset for drawing
                // the particle trajectory In particular, the cosmic rays will
                // not be correctly placed without this
                g4Ticks = clockData.TPCG4Time2Tick(mcPart->T()) +
                          detProp.GetXTicksOffset(0, tpc, cryo) - trigger_offset(clockData);

                xOffset = detProp.ConvertTicksToX(g4Ticks, 0, tpc, cryo);
              }
              else {
                xOffset = 0;
                tpcminx = 1.0;
                tpcmaxx = -1.0;
                coeff = 0.0;
                readoutwindowsize = 0.0;
              }
            }

            // Now move the hit position to correspond to the timing
            xPos += xOffset;

            bool inreadoutwindow = false;
            if (coeff < 0) {
              if ((xPos > readoutwindowsize) && (xPos < tpcmaxx)) inreadoutwindow = true;
            }
            else if (coeff > 0) {
              if ((xPos > tpcminx) && (xPos < readoutwindowsize)) inreadoutwindow = true;
            }

            if (!inreadoutwindow) continue;

            // Check fiducial limits
            if (xPos > xMinimum && xPos < xMaximum) {
              hitPosX[hitCount] = xPos;
              hitPosY[hitCount] = yPos;
              hitPosZ[hitCount] = zPos;
              hitCount++;
            }
          }

          TPolyLine& pl = view->AddPolyLine(
            1, evd::Style::ColorFromPDG(mcPart->PdgCode()), 1, 1); //kFullCircle, msize);

          // Draw neutrals as a gray dotted line to help fade into background a bit...
          if (partCharge == 0.) {
            pl.SetLineColor(13);
            pl.SetLineStyle(3);
            pl.SetLineWidth(1);
          }

          if (proj == evd::kXY)
            pl.SetPolyLine(hitCount, hitPosX.get(), hitPosY.get(), "");
          else if (proj == evd::kXZ)
            pl.SetPolyLine(hitCount, hitPosZ.get(), hitPosX.get(), "");
          else if (proj == evd::kYZ)
            pl.SetPolyLine(hitCount, hitPosZ.get(), hitPosY.get(), "");
        }
      }
    }

    // Now we set up and build the map between MCParticles and a vector of positions obtained from the voxels
    std::map<const simb::MCParticle*, std::vector<std::vector<double>>> partToPosMap;

    sim::LArVoxelList::const_iterator vxitr;
    for (vxitr = voxels.begin(); vxitr != voxels.end(); vxitr++) {
      const sim::LArVoxelData& vxd = (*vxitr).second;

      for (size_t partIdx = 0; partIdx < vxd.NumberParticles(); partIdx++) {
        if (vxd.Energy(partIdx) > drawopt->fMinEnergyDeposition) {
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

      tpcminx = 1.0;
      tpcmaxx = -1.0;
      xOffset = 0.0;
      g4Ticks = 0.0;
      std::vector<std::array<double, 3>> posVecCorr;
      posVecCorr.reserve(partToPosMapItr->second.size());
      coeff = 0.0;
      readoutwindowsize = 0.0;

      // Now loop over points and add to trajectory
      for (size_t posIdx = 0; posIdx < partToPosMapItr->second.size(); posIdx++) {
        const std::vector<double>& posVec = partToPosMapItr->second[posIdx];

        if ((posVec[0] < tpcminx) || (posVec[0] > tpcmaxx)) {
          vtx[0] = posVec[0];
          vtx[1] = posVec[1];
          vtx[2] = posVec[2];
          geo::TPCID tpcid = geom->FindTPCAtPosition(vtx);
          unsigned int cryo = geom->FindCryostatAtPosition(vtx);

          if (tpcid.isValid) {
            unsigned int tpc = tpcid.TPC;

            const geo::TPCGeo& tpcgeo = geom->GetElement(tpcid);
            tpcminx = tpcgeo.MinX();
            tpcmaxx = tpcgeo.MaxX();

            coeff = detProp.GetXTicksCoefficient(tpc, cryo);
            readoutwindowsize = detProp.ConvertTicksToX(detProp.ReadOutWindowSize(), 0, tpc, cryo);
            // The following is meant to get the correct offset for drawing the
            // particle trajectory In particular, the cosmic rays will not be
            // correctly placed without this
            g4Ticks = clockData.TPCG4Time2Tick(mcPart->T()) +
                      detProp.GetXTicksOffset(0, tpc, cryo) - trigger_offset(clockData);

            xOffset = detProp.ConvertTicksToX(g4Ticks, 0, tpc, cryo);
          }
          else {
            xOffset = 0;
            tpcminx = 1.0;
            tpcmaxx = -1.0;
            coeff = 0.0;
            readoutwindowsize = 0.0;
          }
        }

        double xCoord = posVec[0] + xOffset;

        bool inreadoutwindow = false;
        if (coeff < 0) {
          if ((xCoord > readoutwindowsize) && (xCoord < tpcmaxx)) inreadoutwindow = true;
        }
        else if (coeff > 0) {
          if ((xCoord > tpcminx) && (xCoord < readoutwindowsize)) inreadoutwindow = true;
        }

        if (inreadoutwindow && (xCoord > xMinimum && xCoord < xMaximum)) {
          posVecCorr.push_back({{xCoord, posVec[1], posVec[2]}});
        }
      }

      TPolyMarker& pm = view->AddPolyMarker(posVecCorr.size(),
                                            evd::Style::ColorFromPDG(mcPart->PdgCode()),
                                            kFullDotMedium,
                                            2); //kFullCircle, msize);

      for (size_t p = 0; p < posVecCorr.size(); ++p) {
        if (proj == evd::kXY)
          pm.SetPoint(p, posVecCorr[p][0], posVecCorr[p][1]);
        else if (proj == evd::kXZ)
          pm.SetPoint(p, posVecCorr[p][2], posVecCorr[p][0]);
        else if (proj == evd::kYZ)
          pm.SetPoint(p, posVecCorr[p][2], posVecCorr[p][1]);
      }
    }

    return;
  }

  //......................................................................
  int
  SimulationDrawer::GetParticle(const art::Event& evt, std::vector<const simb::MCParticle*>& plist)
  {
    plist.clear();

    if (evt.isRealData()) return 0;

    art::ServiceHandle<evd::SimulationDrawingOptions const> drawopt;

    std::vector<const simb::MCParticle*> temp;

    art::View<simb::MCParticle> plcol;
    // use get by Type because there should only be one collection of these in the event
    try {
      evt.getView(drawopt->fG4ModuleLabel, plcol);
      for (unsigned int i = 0; i < plcol.vals().size(); ++i) {
        temp.push_back(plcol.vals().at(i));
      }
      temp.swap(plist);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetRawDigits", e);
    }

    return plist.size();
  }

  //......................................................................

  int
  SimulationDrawer::GetMCTruth(const art::Event& evt, std::vector<const simb::MCTruth*>& mcvec)
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
      writeErrMsg("GetMCTruth", e);
    }

    return mcvec.size();
  }

  //......................................................................

  void
  SimulationDrawer::HiLite(int trkId, bool dohilite)
  {
    fHighlite[trkId] = dohilite;
  }

} // namespace
////////////////////////////////////////////////////////////////////////
