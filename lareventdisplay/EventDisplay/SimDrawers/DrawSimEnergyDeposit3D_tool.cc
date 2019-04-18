////////////////////////////////////////////////////////////////////////
/// \file   DrawSimEnergyDeposit3D_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lareventdisplay/EventDisplay/Style.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

#include "TF1.h"
#include "TPolyLine.h"

#include <fstream>
#include <cmath>

namespace evdb_tool
{

class DrawSimEnergyDeposit3D : public ISim3DDrawer
{
public:
    explicit DrawSimEnergyDeposit3D(const fhicl::ParameterSet& pset);

    ~DrawSimEnergyDeposit3D();

    void Draw(const art::Event&, evdb::View3D*) const override;

private:

//    double fMinX;
//    double fMaxX;
//    double fMinY;
//    double fMaxY;
//    double fMinZ;
//    double fMaxZ;

};

//----------------------------------------------------------------------
// Constructor.
DrawSimEnergyDeposit3D::DrawSimEnergyDeposit3D(const fhicl::ParameterSet& pset)
{
//    fNumPoints     = pset.get< int>("NumPoints",     1000);
//    fFloatBaseline = pset.get<bool>("FloatBaseline", false);

    return;
}

DrawSimEnergyDeposit3D::~DrawSimEnergyDeposit3D()
{

}


void DrawSimEnergyDeposit3D::Draw(const art::Event& evt, evdb::View3D* view) const
{
    art::ServiceHandle<evd::SimulationDrawingOptions const> drawOpt;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowSimEnergyInfo) return;

    // Recover a handle to the collection of MCParticles
    // We need these so we can determine the offset (if any)
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;

    evt.getByLabel(drawOpt->fG4ModuleLabel, mcParticleHandle);

    if (!mcParticleHandle.isValid()) return;

    // Create a mapping between track ID's and MCParticles
    using TrackToMcParticleMap = std::map<int, const simb::MCParticle*>;

    TrackToMcParticleMap trackToMcParticleMap;

    for(const auto& mcParticle : *mcParticleHandle) trackToMcParticleMap[mcParticle.TrackId()] = &mcParticle;

    // Now recover the simchannels
    art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyDepositHandle;

    evt.getByLabel(drawOpt->fSimEnergyLabel, simEnergyDepositHandle);

    if (simEnergyDepositHandle.isValid() && simEnergyDepositHandle->size() > 0)
    {
        mf::LogDebug("SimEnergyDeposit3DDrawer") << "Starting loop over " << simEnergyDepositHandle->size() << " SimEnergyDeposits, " << std::endl;

        // Get the detector properties, clocks...
        detinfo::DetectorProperties const*      theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
        detinfo::DetectorClocks     const*      detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
        art::ServiceHandle<geo::Geometry const> geom;

        // First step is to create a map between MCParticle and SimEnergyDeposit objects...
        using MCPartToSimEnergyMap = std::map<const simb::MCParticle*, std::vector<const sim::SimEnergyDeposit*>>;

        MCPartToSimEnergyMap mcPartToSimEnergyMap;

        // Go through the SimEnergyDeposits and populate the map
        for(const auto& simEnergyDeposit : *simEnergyDepositHandle)
        {
            TrackToMcParticleMap::const_iterator trackMCItr = trackToMcParticleMap.find(simEnergyDeposit.TrackID());

            if (trackMCItr == trackToMcParticleMap.end()) continue;

            mcPartToSimEnergyMap[trackMCItr->second].push_back(&simEnergyDeposit);
        }

        // Would like to draw the deposits as markers with colors given by particle id
        // So we make two passes, first to fill a map with color the key and positions for the markers
        std::map<int,std::vector<sim::SimEnergyDeposit::Point_t>> colorToPositionMap;

        // Now we loop through and build the mapping of color to positions
        for(const auto& mcPartToSimEnergy : mcPartToSimEnergyMap)
        {
            // The first task we need to take on is to find the offset for the energy deposit
            // This is for the case of "out of time" particles... (e.g. cosmic rays)
            double g4Ticks(detClocks->TPCG4Time2Tick(mcPartToSimEnergy.first->T())-theDetector->TriggerOffset());
            double xOffset(0.);
            double xPosMinTick(0.);
            double xPosMaxTick(std::numeric_limits<double>::max());

            for(const auto& simEnergyDeposit : mcPartToSimEnergy.second)
            {
                sim::SimEnergyDeposit::Point_t point = simEnergyDeposit->MidPoint();

                // If we have cosmic rays then we need to get the offset which allows translating from
                // when they were generated vs when they were tracked.
                // Note that this also explicitly checks that they are in a TPC volume
                try
                {
                    geo::TPCID   tpcID = geom->PositionToTPCID(point);
                    geo::PlaneID planeID(tpcID,0);

                    xPosMinTick = theDetector->ConvertTicksToX(0,planeID);
                    xPosMaxTick = theDetector->ConvertTicksToX(theDetector->NumberTimeSamples(),planeID);
                    xOffset     = theDetector->ConvertTicksToX(g4Ticks, planeID) - xPosMinTick;

                    if (xPosMaxTick < xPosMinTick) std::swap(xPosMinTick,xPosMaxTick);
                }
                catch(...) {continue;}

                colorToPositionMap[evd::Style::ColorFromPDG(simEnergyDeposit->PdgCode())].emplace_back(sim::SimEnergyDeposit::Point_t(point.X() + xOffset, point.Y(), point.Z()));
            }
        }

        // Now we can do some drawing
        for(const auto& pair : colorToPositionMap)
        {
            int colorIdx(pair.first);
            int markerIdx(kFullDotMedium);
            int markerSize(2);

            TPolyMarker3D& pm = view->AddPolyMarker3D(1, colorIdx, markerIdx, markerSize);

            // Import positions into an array
            std::vector<double> posArrayVec;
            int    hitCount(0);

            posArrayVec.resize(3 * pair.second.size());

            for(const auto& point : pair.second)
            {
                posArrayVec[3*hitCount    ] = point.X();
                posArrayVec[3*hitCount + 1] = point.Y();
                posArrayVec[3*hitCount + 2] = point.Z();
                hitCount++;
            }

            pm.SetPolyMarker(hitCount, posArrayVec.data(), markerIdx);
        }
    }

    return;
}

DEFINE_ART_CLASS_TOOL(DrawSimEnergyDeposit3D)
}
