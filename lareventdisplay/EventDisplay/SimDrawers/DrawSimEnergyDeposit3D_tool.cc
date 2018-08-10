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
    art::ServiceHandle<evd::SimulationDrawingOptions> drawOpt;
    
    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowMCTruthTrajectories) return;
    
    //  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
//    detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
//    detinfo::DetectorClocks     const* detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    art::ServiceHandle<geo::Geometry>  geom;
    
    // Recover a handle to the collection of MCParticles
    art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyDepositHandle;
    
    evt.getByLabel(drawOpt->fSimEnergyLabel, simEnergyDepositHandle);
    
    if (simEnergyDepositHandle.isValid() && simEnergyDepositHandle->size() > 0)
    {
        mf::LogDebug("SimEnergyDeposit3DDrawer") << "Starting loop over " << simEnergyDepositHandle->size() << " SimEnergyDeposits, " << std::endl;
        
        // Would like to draw the deposits as markers with colors given by particle id
        // So we make two passes, first to fill a map with color the key and positions for the markers
        std::map<int,std::vector<sim::SimEnergyDeposit::Point_t>> colorToPositionMap;
        
        // First loop through energy deposits
        for(const auto& simEnergyDeposit : *simEnergyDepositHandle)
        {
            colorToPositionMap[evd::Style::ColorFromPDG(simEnergyDeposit.PdgCode())].emplace_back(simEnergyDeposit.MidPoint());
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
