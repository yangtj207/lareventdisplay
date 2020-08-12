////////////////////////////////////////////////////////////////////////
/// \file   PCA3DDrawer_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/fwd.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/fwd.h"
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"
namespace evdb {
  class View3D;
}

namespace evdb_tool {

  class PCA3DDrawer : public I3DDrawer {
  public:
    explicit PCA3DDrawer(const fhicl::ParameterSet&);

    ~PCA3DDrawer();

    void Draw(const art::Event&, evdb::View3D*) const override;

  private:
  };

  //----------------------------------------------------------------------
  // Constructor.
  PCA3DDrawer::PCA3DDrawer(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.

    return;
  }

  PCA3DDrawer::~PCA3DDrawer() {}

  void
  PCA3DDrawer::Draw(const art::Event& evt, evdb::View3D* view) const
  {
    /*
    art::ServiceHandle<evd::SimulationDrawingOptions const> drawOpt;

    // If the option is turned off, there's nothing to do
    if (!drawOpt->fShowMCTruthTrajectories) return;

    //  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
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

  DEFINE_ART_CLASS_TOOL(PCA3DDrawer)
}
