////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawerHitCharge_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
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

  class SpacePoint3DDrawerHitCharge : public ISpacePoints3D {
  public:
    explicit SpacePoint3DDrawerHitCharge(const fhicl::ParameterSet&);

    ~SpacePoint3DDrawerHitCharge();

    void Draw(const std::vector<art::Ptr<recob::SpacePoint>>&, // Space points
              evdb::View3D*,                                   // 3D display
              int,                                             // Color
              int,                                             // Marker
              float,                                           // Size) const override;
              const art::FindManyP<recob::Hit>*                // pointer to associated hits
    ) const;

  private:
    double getSpacePointCharge(const art::Ptr<recob::SpacePoint>&,
                               const art::FindManyP<recob::Hit>*) const;
    double chargeIntegral(double, double, double, double, int, int) const;

    bool fUseAbsoluteScale;
    float fMinHitCharge;
    float fMaxHitCharge;
  };

  //----------------------------------------------------------------------
  // Constructor.
  SpacePoint3DDrawerHitCharge::SpacePoint3DDrawerHitCharge(const fhicl::ParameterSet& pset)
  {
    //    fNumPoints     = pset.get< int>("NumPoints",     1000);
    //    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.
    fUseAbsoluteScale = pset.get<bool>("UseAbsoluteScale", false);
    fMinHitCharge = pset.get<float>("MinHitCharge", 0.);
    fMaxHitCharge = pset.get<float>("MaxHitCharge", 2500.);

    return;
  }

  SpacePoint3DDrawerHitCharge::~SpacePoint3DDrawerHitCharge()
  {
    return;
  }

  void SpacePoint3DDrawerHitCharge::Draw(const std::vector<art::Ptr<recob::SpacePoint>>& hitsVec,
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

    float minHitCharge(std::numeric_limits<float>::max());
    float maxHitCharge(std::numeric_limits<float>::lowest());

    if (fUseAbsoluteScale) {
      minHitCharge = fMinHitCharge;
      maxHitCharge = fMaxHitCharge;
    }
    else
    // Find the range in the input space point list
    {
      for (const auto& spacePoint : hitsVec) {
        //float hitCharge = getSpacePointCharge(spacePoint, hitAssnVec);
        float hitCharge = spacePoint->ErrXYZ()[1];

        minHitCharge = std::min(minHitCharge, hitCharge);
        maxHitCharge = std::max(maxHitCharge, hitCharge);
      }
    }

    // Make sure we really have something here
    if (maxHitCharge > minHitCharge) {
      float hitChiSqScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) /
                          (maxHitCharge - minHitCharge));

      for (const auto& spacePoint : hitsVec) {
        float hitCharge = getSpacePointCharge(spacePoint, hitAssnVec);

        if (hitCharge > 0.) {
          float chgFactor = cst->fRecoQLow[geo::kCollection] + hitChiSqScale * hitCharge;
          int chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(chgFactor);
          const double* pos = spacePoint->XYZ();
          const double* err = spacePoint->ErrXYZ();

          colorToHitMap[chargeColorIdx].push_back(
            HitPosition() = {{pos[0], pos[1], pos[2], err[3], err[3], err[5]}});
        }
      }

      for (auto& hitPair : colorToHitMap) {
        TPolyMarker3D& pm =
          view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotLarge, 0.25);
        for (const auto& hit : hitPair.second)
          pm.SetNextPoint(hit[0], hit[1], hit[2]);
      }
    }

    return;
  }

  double SpacePoint3DDrawerHitCharge::getSpacePointCharge(
    const art::Ptr<recob::SpacePoint>& spacePoint,
    const art::FindManyP<recob::Hit>* hitAssnVec) const
  {
    double totalCharge(0.);

    // Need to recover the integrated charge from the collection plane, so need to loop through associated hits
    const std::vector<art::Ptr<recob::Hit>>& hit2DVec(hitAssnVec->at(spacePoint.key()));

    float hitCharge(0.);
    int lowIndex(std::numeric_limits<int>::min());
    int hiIndex(std::numeric_limits<int>::max());

    for (const auto& hit2D : hit2DVec) {
      int hitStart = hit2D->PeakTime() - 2. * hit2D->RMS() - 0.5;
      int hitStop = hit2D->PeakTime() + 2. * hit2D->RMS() + 0.5;

      lowIndex = std::max(hitStart, lowIndex);
      hiIndex = std::min(hitStop + 1, hiIndex);

      hitCharge += hit2D->Integral();
    }

    if (!hit2DVec.empty()) hitCharge /= float(hit2DVec.size());

    if (hitCharge > 0.) {
      if (hiIndex > lowIndex) {
        for (const auto& hit2D : hit2DVec)
          totalCharge += chargeIntegral(
            hit2D->PeakTime(), hit2D->PeakAmplitude(), hit2D->RMS(), 1., lowIndex, hiIndex);

        totalCharge /= float(hit2DVec.size());
      }
    }

    return totalCharge;
  }

  double SpacePoint3DDrawerHitCharge::chargeIntegral(double peakMean,
                                                     double peakAmp,
                                                     double peakWidth,
                                                     double areaNorm,
                                                     int low,
                                                     int hi) const
  {
    double integral(0);

    for (int sigPos = low; sigPos < hi; sigPos++)
      integral += peakAmp * TMath::Gaus(double(sigPos) + 0.5, peakMean, peakWidth);

    return integral;
  }

  DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawerHitCharge)
}
