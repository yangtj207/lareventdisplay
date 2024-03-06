////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawerStandard_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lareventdisplay/EventDisplay/3DDrawers/ISpacePoints3D.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include "nuevdb/EventDisplayBase/View3D.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TPolyMarker3D.h"

namespace evdb_tool {

  class SpacePoint3DDrawerStandard : public ISpacePoints3D {
  public:
    explicit SpacePoint3DDrawerStandard(const fhicl::ParameterSet&);

    ~SpacePoint3DDrawerStandard();

    void Draw(const std::vector<art::Ptr<recob::SpacePoint>>&, // Space points
              evdb::View3D*,                                   // 3D display
              int,                                             // Color
              int,                                             // Marker
              float,                                           // Size) const override;
              const art::FindManyP<recob::Hit>*                // pointer to associated hits
    ) const;

  private:
  };

  //----------------------------------------------------------------------
  // Constructor.
  SpacePoint3DDrawerStandard::SpacePoint3DDrawerStandard(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.

    return;
  }

  SpacePoint3DDrawerStandard::~SpacePoint3DDrawerStandard()
  {
    return;
  }

  void SpacePoint3DDrawerStandard::Draw(const std::vector<art::Ptr<recob::SpacePoint>>& spts,
                                        evdb::View3D* view,
                                        int color,
                                        int marker,
                                        float size,
                                        const art::FindManyP<recob::Hit>*) const
  {
    // Get services.

    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    // Organize space points into separate collections according to the color
    // we want them to be.
    // If option If option fColorSpacePointsByChisq is false, this means
    // having a single collection with color inherited from the prong
    // (specified by the argument color).

    std::map<int, std::vector<const recob::SpacePoint*>> spmap; // Indexed by color.
    int spcolor = color;

    for (auto& pspt : spts) {
      //std::cout<<pspt<<std::endl;
      //if(pspt == 0) throw cet::exception("RecoBaseDrawer:DrawSpacePoint3D") << "space point is null\n";

      // For rainbow effect, choose root colors in range [51,100].
      // We are using 100=best (red), 51=worst (blue).

      //        if (pspt->Chisq() > -100.) continue;

      spcolor = 12;

      if (pspt->Chisq() < -100.) spcolor = 16;

      if (recoOpt->fColorSpacePointsByChisq) {
        spcolor = 100 - 2.5 * pspt->Chisq();

        if (spcolor < 51) spcolor = 51;
        if (spcolor > 100) spcolor = 100;
      }
      else
        spcolor = color;
      //if (pspt->Chisq() < -1.) spcolor += 6;

      spmap[spcolor].push_back(&*pspt);
    }

    // Loop over colors.
    // Note that larger (=better) space points are plotted on
    // top for optimal visibility.

    for (auto const& [spcolor, psps] : spmap) {

      // Make and fill a polymarker.

      TPolyMarker3D& pm = view->AddPolyMarker3D(psps.size(), spcolor, marker, size);

      for (size_t s = 0; s < psps.size(); ++s) {
        const recob::SpacePoint& spt = *psps[s];

        const double* xyz = spt.XYZ();
        pm.SetPoint(s, xyz[0], xyz[1], xyz[2]);
      }
    }

    return;
  }

  DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawerStandard)
}
