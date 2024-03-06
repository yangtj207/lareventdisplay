////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawerChiSquare_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lareventdisplay/EventDisplay/3DDrawers/ISpacePoints3D.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"

#include "nuevdb/EventDisplayBase/View3D.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TPolyMarker3D.h"

namespace evdb_tool {

  class SpacePoint3DDrawerChiSquare : public ISpacePoints3D {
  public:
    explicit SpacePoint3DDrawerChiSquare(const fhicl::ParameterSet&);

    ~SpacePoint3DDrawerChiSquare();

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
  SpacePoint3DDrawerChiSquare::SpacePoint3DDrawerChiSquare(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.

    return;
  }

  SpacePoint3DDrawerChiSquare::~SpacePoint3DDrawerChiSquare()
  {
    return;
  }

  void SpacePoint3DDrawerChiSquare::Draw(const std::vector<art::Ptr<recob::SpacePoint>>& hitsVec,
                                         evdb::View3D* view,
                                         int color,
                                         int marker,
                                         float size,
                                         const art::FindManyP<recob::Hit>* hitAssns) const
  {
    // Get services.
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;

    using HitPosition = std::array<double, 6>;
    std::map<int, std::vector<HitPosition>> colorToHitMap;

    float minHitChiSquare(0.);
    float maxHitChiSquare(2.);
    float hitChiSqScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) /
                        (maxHitChiSquare - minHitChiSquare));

    for (const auto& spacePoint : hitsVec) {
      const double* pos = spacePoint->XYZ();
      const double* err = spacePoint->ErrXYZ();

      int chargeColorIdx(0);
      float spacePointChiSq(spacePoint->Chisq());

      float hitChiSq = std::max(minHitChiSquare, std::min(maxHitChiSquare, spacePointChiSq));

      float chgFactor = cst->fRecoQHigh[geo::kCollection] - hitChiSqScale * hitChiSq;

      chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(chgFactor);

      colorToHitMap[chargeColorIdx].push_back(
        HitPosition() = {{pos[0], pos[1], pos[2], err[3], err[3], err[5]}});
    }

    for (auto& hitPair : colorToHitMap) {
      TPolyMarker3D& pm =
        view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotLarge, 0.17);
      for (const auto& hit : hitPair.second)
        pm.SetNextPoint(hit[0], hit[1], hit[2]);
    }

    return;
  }

  DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawerChiSquare)
}
