////////////////////////////////////////////////////////////////////////
/// \file   PFParticle3DDrawer_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lareventdisplay/EventDisplay/Style.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

#include <fstream>

namespace evdb_tool
{

class PFParticle3DDrawer : public I3DDrawer
{
public:
    explicit PFParticle3DDrawer(const fhicl::ParameterSet&);

    ~PFParticle3DDrawer();

    void Draw(const art::Event&, evdb::View3D*) const override;

private:
};

//----------------------------------------------------------------------
// Constructor.
PFParticle3DDrawer::PFParticle3DDrawer(const fhicl::ParameterSet& pset)
{
//    fNumPoints     = pset.get< int>("NumPoints",     1000);
//    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.

    return;
}

PFParticle3DDrawer::~PFParticle3DDrawer()
{
}

void PFParticle3DDrawer::Draw(const art::Event& evt, evdb::View3D* view) const
{
/*
    art::ServiceHandle<evd::SimulationDrawingOptions const> drawOpt;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowMCTruthTrajectories) return;

    //  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    detinfo::DetectorProperties const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks     const* detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    art::ServiceHandle<geo::Geometry const>  geom;

    // Recover a handle to the collection of MCParticles
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;

    evt.getByLabel(drawOpt->fG4ModuleLabel, mcParticleHandle);

    // Define a couple of colors for neutrals and if we gray it out...
    int neutralColor(12);
    int grayedColor(15);
    int neutrinoColor(38);
*/
    return;
}

DEFINE_ART_CLASS_TOOL(PFParticle3DDrawer)
}
