////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawerHitAsymmetry_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lareventdisplay/EventDisplay/3DDrawers/ISpacePoints3D.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"

#include "nuevdb/EventDisplayBase/View3D.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TMath.h"
#include "TPolyMarker3D.h"

namespace evdb_tool {

  class SpacePoint3DDrawerHitAsymmetry : public ISpacePoints3D {
  public:
    explicit SpacePoint3DDrawerHitAsymmetry(const fhicl::ParameterSet&);

    ~SpacePoint3DDrawerHitAsymmetry();

    void Draw(const std::vector<art::Ptr<recob::SpacePoint>>&, // Space points
              evdb::View3D*,                                   // 3D display
              int,                                             // Color
              int,                                             // Marker
              float,                                           // Size) const override;
              const art::FindManyP<recob::Hit>*                // pointer to associated hits
    ) const;

  private:
    float fMinAsymmetry;
    float fMaxAsymmetry;
  };

  //----------------------------------------------------------------------
  // Constructor.
  SpacePoint3DDrawerHitAsymmetry::SpacePoint3DDrawerHitAsymmetry(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.
    fMinAsymmetry = pset.get<float>("MinAsymmetry", -1.);
    fMaxAsymmetry = pset.get<float>("MaxAsymmetry", 1.);

    return;
  }

  SpacePoint3DDrawerHitAsymmetry::~SpacePoint3DDrawerHitAsymmetry()
  {
    return;
  }

  void SpacePoint3DDrawerHitAsymmetry::Draw(const std::vector<art::Ptr<recob::SpacePoint>>& hitsVec,
                                            evdb::View3D* view,
                                            int color,
                                            int marker,
                                            float size,
                                            const art::FindManyP<recob::Hit>* hitAssnVec) const
  {
    // Let's not crash
    if (hitsVec.empty() || !hitAssnVec) return;

    // Get services.
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;

    using HitPosition = std::array<double, 6>;
    std::map<int, std::vector<HitPosition>> colorToHitMap;

    // Get the scale factor
    float asymmetryScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) /
                         (fMaxAsymmetry - fMinAsymmetry));

    for (const auto& spacePoint : hitsVec) {
      float hitAsymmetry = spacePoint->ErrXYZ()[3] - fMinAsymmetry;

      if (std::abs(hitAsymmetry) <= fMaxAsymmetry - fMinAsymmetry) {
        float chgFactor = cst->fRecoQLow[geo::kCollection] + asymmetryScale * hitAsymmetry;
        int chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(chgFactor);
        const double* pos = spacePoint->XYZ();
        const double* err = spacePoint->ErrXYZ();

        colorToHitMap[chargeColorIdx].push_back(
          HitPosition() = {{pos[0], pos[1], pos[2], err[2], err[2], err[5]}});
      }
    }

    for (auto& hitPair : colorToHitMap) {
      TPolyMarker3D& pm =
        view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotLarge, 0.25);
      for (const auto& hit : hitPair.second)
        pm.SetNextPoint(hit[0], hit[1], hit[2]);
    }

    return;
  }

  DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawerHitAsymmetry)
}
