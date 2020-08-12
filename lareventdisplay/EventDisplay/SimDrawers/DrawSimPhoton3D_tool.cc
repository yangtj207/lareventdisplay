////////////////////////////////////////////////////////////////////////
/// \file   DrawSimPhoton3D_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TPolyLine3D.h"

// Eigen
#include <Eigen/Core>

namespace evdb_tool {

  class DrawSimPhoton3D : public ISim3DDrawer {
  public:
    explicit DrawSimPhoton3D(const fhicl::ParameterSet& pset);

    ~DrawSimPhoton3D();

    void Draw(const art::Event&, evdb::View3D*) const override;

  private:
    void DrawRectangularBox(evdb::View3D*,
                            const Eigen::Vector3f&,
                            const Eigen::Vector3f&,
                            int,
                            int,
                            int) const;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawSimPhoton3D::DrawSimPhoton3D(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);

    return;
  }

  DrawSimPhoton3D::~DrawSimPhoton3D() {}

  void
  DrawSimPhoton3D::Draw(const art::Event& evt, evdb::View3D* view) const
  {
    art::ServiceHandle<evd::SimulationDrawingOptions> drawOpt;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowSimPhotonInfo) return;

    // Recover a handle to the collection of MCParticles
    // We need these so we can determine the offset (if any)
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;

    evt.getByLabel(drawOpt->fG4ModuleLabel, mcParticleHandle);

    if (!mcParticleHandle.isValid()) return;

    // Create a mapping between track ID's and MCParticles
    using TrackToMcParticleMap = std::map<int, const simb::MCParticle*>;

    TrackToMcParticleMap trackToMcParticleMap;

    for (const auto& mcParticle : *mcParticleHandle)
      trackToMcParticleMap[mcParticle.TrackId()] = &mcParticle;

    // Now recover the simphotons
    art::Handle<std::vector<sim::SimPhotons>> simPhotonsHandle;

    evt.getByLabel(drawOpt->fSimPhotonLabel, simPhotonsHandle);

    if (simPhotonsHandle.isValid() && simPhotonsHandle->size() > 0) {
      mf::LogDebug("SimPhoton3DDrawer")
        << "Starting loop over " << simPhotonsHandle->size() << " SimPhotons, " << std::endl;

      // Get the detector properties, clocks...
      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<evd::ColorDrawingOptions> cst;

      // First step is to create a map between MCParticle and SimEnergyDeposit objects...
      using MCPartToOnePhotonMap =
        std::map<const simb::MCParticle*, std::vector<const sim::OnePhoton*>>;
      using ChanToMCPartToOnePhotonMap = std::map<int, MCPartToOnePhotonMap>;

      ChanToMCPartToOnePhotonMap chanToMCPartToOnePhotonMap;

      // Go through the SimEnergyDeposits and populate the map
      for (const auto& simPhoton : *simPhotonsHandle) {
        MCPartToOnePhotonMap& mcPartToOnePhotonMap =
          chanToMCPartToOnePhotonMap[simPhoton.OpChannel()];

        for (const auto& onePhoton : simPhoton) {
          TrackToMcParticleMap::const_iterator trackMCItr =
            trackToMcParticleMap.find(onePhoton.MotherTrackID);

          if (trackMCItr == trackToMcParticleMap.end()) continue;

          mcPartToOnePhotonMap[trackMCItr->second].push_back(&onePhoton);
        }
      }

      // Mapping of energy deposited per channel...
      std::map<int, float> channelToEnergyMap;

      // Keep track of mininum and maximum
      float maxEnergy = std::numeric_limits<float>::lowest();
      float minEnergy = std::numeric_limits<float>::max();

      // Go through everything and get cenergy desposited per channel...
      for (const auto& chanToMCPartToOnePhoton : chanToMCPartToOnePhotonMap) {
        float totalE(0.);

        // Go through all contributors to this channel
        for (const auto& mcPartToOnePhoton : chanToMCPartToOnePhoton.second) {
          // Current scheme will ignore displacement in time... need to come back to this at some point...
          //                // The first task we need to take on is to find the offset for the energy deposit
          //                // This is for the case of "out of time" particles... (e.g. cosmic rays)
          //                double g4Ticks(detClocks->TPCG4Time2Tick(mcPartToOnePhoton.first->T())-theDetector->TriggerOffset());
          //                double xOffset(0.);
          //                double xPosMinTick(0.);
          //                double xPosMaxTick(std::numeric_limits<double>::max());

          for (const auto& onePhoton : mcPartToOnePhoton.second) {
            //                    Eigen::Vector3f point(onePhoton->InitialPosition.X(),onePhoton->InitialPosition.Y(),onePhoton->InitialPosition.Z());
            //
            //                    std::cout << "    - Initial: " << onePhoton->InitialPosition.X() << "/" << onePhoton->InitialPosition.Y() << "/" << onePhoton->InitialPosition.Z() << ", final : " //<<   onePhoton->FinalLocalPosition.X() << "/" << onePhoton->FinalLocalPosition.Y() << "/" << onePhoton->FinalLocalPosition.Z() << std::endl;
            //
            //                    // If we have cosmic rays then we need to get the offset which allows translating from
            //                    // when they were generated vs when they were tracked.
            //                    // Note that this also explicitly checks that they are in a TPC volume
            //                    try
            //                    {
            //                        geo::TPCID   tpcID = geom->PositionToTPCID(geo::Point_t(point[0],point[1],point[2]));
            //
            //                        if (tpcID.Cryostat == geo::CryostatID::InvalidID || tpcID.TPC == geo::TPCID::InvalidID) continue;
            //
            //                        geo::PlaneID planeID(tpcID,0);
            //
            //                        xPosMinTick = theDetector->ConvertTicksToX(0,planeID);
            //                        xPosMaxTick = theDetector->ConvertTicksToX(theDetector->NumberTimeSamples(),planeID);
            //                        xOffset     = theDetector->ConvertTicksToX(g4Ticks, planeID) - xPosMinTick;
            //
            //                        if (xPosMaxTick < xPosMinTick) std::swap(xPosMinTick,xPosMaxTick);
            //                    }
            //                    catch(...) {continue;}

            // Recover the deposited energy
            totalE += onePhoton->Energy;
          }
        }

        channelToEnergyMap[chanToMCPartToOnePhoton.first] = totalE;

        maxEnergy = std::max(maxEnergy, totalE);
        minEnergy = std::min(minEnergy, totalE);
      }

      // Get the scale factor from energy deposit range
      float yzWidthScale(1. / (maxEnergy - minEnergy));
      float energyDepositScale(
        (cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) * yzWidthScale);

      // Go through the channels and draw the objects
      for (const auto& channelToEnergy : channelToEnergyMap) {
        // Recover the color index based on energy
        float widthFactor =
          0.95 * std::max(float(0.), std::min(float(1.), yzWidthScale* channelToEnergy.second));
        float energyFactor =
          cst->fRecoQLow[geo::kCollection] + energyDepositScale * channelToEnergy.second;

        // Recover the position for this channel
        const geo::OpDetGeo& opHitGeo = geom->OpDetGeoFromOpChannel(channelToEnergy.first);
        const geo::Point_t& opHitPos = opHitGeo.GetCenter();
        float xWidth = 0.01;
        float zWidth = widthFactor * opHitGeo.HalfW();
        float yWidth = widthFactor * opHitGeo.HalfH();

        // Get widths of box to draw
        Eigen::Vector3f coordsLo(
          opHitPos.X() - xWidth, opHitPos.Y() - yWidth, opHitPos.Z() - zWidth);
        Eigen::Vector3f coordsHi(
          opHitPos.X() + xWidth, opHitPos.Y() + yWidth, opHitPos.Z() + zWidth);

        int energyColorIdx = cst->CalQ(geo::kCollection).GetColor(energyFactor);

        DrawRectangularBox(view, coordsLo, coordsHi, energyColorIdx, 1, 1);
      }
    }

    return;
  }

  void
  DrawSimPhoton3D::DrawRectangularBox(evdb::View3D* view,
                                      const Eigen::Vector3f& coordsLo,
                                      const Eigen::Vector3f& coordsHi,
                                      int color,
                                      int width,
                                      int style) const
  {
    TPolyLine3D& top = view->AddPolyLine3D(5, color, width, style);
    top.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(1, coordsHi[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(2, coordsHi[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& side = view->AddPolyLine3D(5, color, width, style);
    side.SetPoint(0, coordsHi[0], coordsHi[1], coordsLo[2]);
    side.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    side.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    side.SetPoint(3, coordsHi[0], coordsHi[1], coordsHi[2]);
    side.SetPoint(4, coordsHi[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& side2 = view->AddPolyLine3D(5, color, width, style);
    side2.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    side2.SetPoint(1, coordsLo[0], coordsLo[1], coordsLo[2]);
    side2.SetPoint(2, coordsLo[0], coordsLo[1], coordsHi[2]);
    side2.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    side2.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& bottom = view->AddPolyLine3D(5, color, width, style);
    bottom.SetPoint(0, coordsLo[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(3, coordsLo[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(4, coordsLo[0], coordsLo[1], coordsLo[2]);

    return;
  }

  DEFINE_ART_CLASS_TOOL(DrawSimPhoton3D)
}
