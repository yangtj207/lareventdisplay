////////////////////////////////////////////////////////////////////////
/// \file   DrawSimEnergyDeposit3D_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "larcore/Geometry/Geometry.h"

#include "TF1.h"
#include "TPolyLine.h"

#include <fstream>

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
//    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    
    // If the option is turned off, there's nothing to do
//    if (!drawopt->fShowMCTruthTrajectories) return;

//    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod)
//    {
//        // Step one is to recover the hits for this label that match the input channel
//        art::InputTag const which = recoOpt->fHitLabels[imod];
//
//        art::Handle< std::vector<recob::Hit> > hitVecHandle;
//        event->getByLabel(which, hitVecHandle);
//    }//end loop over HitFinding modules

    return;
}
    
DEFINE_ART_CLASS_TOOL(DrawSimEnergyDeposit3D)
}
