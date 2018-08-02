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
    art::Handle< std::vector<sim::SimEnergyDeposit>> simEnergyDepositHandle;
    
    evt.getByLabel(drawOpt->fSimEnergyLabel, simEnergyDepositHandle);
    
    if (simEnergyDepositHandle.isValid() && simEnergyDepositHandle->size() > 0)
    {
        // Define a couple of colors for neutrals and if we gray it out...
//    int neutralColor(12);
//    int grayedColor(15);
//    int neutrinoColor(38);
    
        mf::LogDebug("SimEnergyDeposit3DDrawer") << "Starting loop over " << simEnergyDepositHandle->size() << " SimEnergyDeposits, " << std::endl;
    
        // The simplest approach is to loop over the SimEnergyDeposits and draw their centers
        std::unique_ptr<double[]> hitPositions(new double[3*simEnergyDepositHandle->size()]);
        int                       hitCount(0);
        
        for(const auto& simEnergyDeposit : *simEnergyDepositHandle)
        {
            geo::Point_t startPoint = simEnergyDeposit.Start();
            geo::Point_t stopPoint  = simEnergyDeposit.End();
        
            hitPositions[3*hitCount    ] = 0.5 * (startPoint.X() + stopPoint.X());
            hitPositions[3*hitCount + 1] = 0.5 * (startPoint.Y() + stopPoint.Y());
            hitPositions[3*hitCount + 2] = 0.5 * (startPoint.Z() + stopPoint.Z());
            hitCount++;
    
            mf::LogDebug("SimEnergyDeposit3DDrawer") << "--> Hit: " << hitCount << " x,y,z: " << 0.5 * (startPoint.X() + stopPoint.X()) << "," << 0.5 * (startPoint.Y() + stopPoint.Y()) << "," << 0.5 * (startPoint.Z() + stopPoint.Z()) << ", # e: " <<     simEnergyDeposit.NumElectrons() << ", # gamma: " << simEnergyDeposit.NumPhotons() << ", edep: " << simEnergyDeposit.Energy() << std::endl;
        }
        
        int colorIdx(3);
        int markerIdx(kFullDotSmall);
        int markerSize(2);
    
        TPolyMarker3D& pm = view->AddPolyMarker3D(1, colorIdx, markerIdx, markerSize);
        pm.SetPolyMarker(hitCount, hitPositions.get(), markerIdx);
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(DrawSimEnergyDeposit3D)
}
