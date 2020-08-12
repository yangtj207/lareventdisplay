////////////////////////////////////////////////////////////////////////
/// \file   DrawGausHits_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWFHitDrawer.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"

#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TF1.h"
#include "TPolyLine.h"

#include <fstream>

namespace evdb_tool {

  class DrawGausHits : public IWFHitDrawer {
  public:
    explicit DrawGausHits(const fhicl::ParameterSet& pset);

    ~DrawGausHits();

    void configure(const fhicl::ParameterSet& pset) override;
    void Draw(evdb::View2D&, raw::ChannelID_t&) const override;

  private:
    using HitParams_t = struct HitParams_t {
      float hitCenter;
      float hitSigma;
      float hitHeight;
      float hitStart;
      float hitEnd;
    };

    using ROIHitParamsVec = std::vector<HitParams_t>;
    using HitParamsVec = std::vector<ROIHitParamsVec>;

    int fNumPoints;
    bool fFloatBaseline;
    std::vector<int> fColorVec;
    mutable std::vector<std::unique_ptr<TF1>> fHitFunctionVec;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawGausHits::DrawGausHits(const fhicl::ParameterSet& pset) { configure(pset); }

  DrawGausHits::~DrawGausHits() {}

  void
  DrawGausHits::configure(const fhicl::ParameterSet& pset)
  {
    fNumPoints = pset.get<int>("NumPoints", 1000);
    fFloatBaseline = pset.get<bool>("FloatBaseline", false);

    fColorVec.clear();
    fColorVec.push_back(kRed);
    fColorVec.push_back(kGreen);
    fColorVec.push_back(kMagenta);
    fColorVec.push_back(kCyan);

    fHitFunctionVec.clear();

    return;
  }

  void
  DrawGausHits::Draw(evdb::View2D& view2D, raw::ChannelID_t& channel) const
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if (!event) return;

    fHitFunctionVec.clear();

    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      // Step one is to recover the hits for this label that match the input channel
      art::InputTag const which = recoOpt->fHitLabels[imod];

      art::Handle<std::vector<recob::Hit>> hitVecHandle;
      event->getByLabel(which, hitVecHandle);

      // Get a container for the subset of hits we are drawing
      art::PtrVector<recob::Hit> hitPtrVec;

      try {
        for (size_t hitIdx = 0; hitIdx < hitVecHandle->size(); hitIdx++) {
          art::Ptr<recob::Hit> hit(hitVecHandle, hitIdx);

          if (hit->Channel() == channel) hitPtrVec.push_back(hit);
        }
      }
      catch (cet::exception& e) {
        mf::LogWarning("DrawGausHits") << "DrawGausHits"
                                       << " failed with message:\n"
                                       << e;
      }

      if (hitPtrVec.empty()) continue;

      // Apparently you cannot trust some hit producers to put the hits in the correct order!
      std::sort(hitPtrVec.begin(), hitPtrVec.end(), [](const auto& left, const auto& right) {
        return left->PeakTime() < right->PeakTime();
      });

      // Get associations to wires
      art::FindManyP<recob::Wire> wireAssnsVec(hitPtrVec, *event, which);
      std::vector<float> wireDataVec;

      // Recover the full (zero-padded outside ROI's) deconvolved waveform for this wire
      if (wireAssnsVec.isValid() && wireAssnsVec.size() > 0) {
        auto hwafp = wireAssnsVec.at(0).front();
        if (!hwafp.isNull() && hwafp.isAvailable()) { wireDataVec = hwafp->Signal(); }
      }

      // Now go through and process the hits back into the hit parameters
      using ROIHitParamsVec = std::vector<HitParams_t>;
      using HitParamsVec = std::vector<ROIHitParamsVec>;

      // Get an initial container for common hits on ROI
      HitParamsVec hitParamsVec;
      ROIHitParamsVec roiHitParamsVec;
      raw::TDCtick_t lastEndTick(10000);

      for (const auto& hit : hitPtrVec) {
        // check roi end condition
        if (hit->PeakTime() - 3. * hit->RMS() > lastEndTick) {
          if (!roiHitParamsVec.empty()) hitParamsVec.push_back(roiHitParamsVec);
          roiHitParamsVec.clear();
        }

        HitParams_t hitParams;

        hitParams.hitCenter = hit->PeakTime();
        hitParams.hitSigma = hit->RMS();
        hitParams.hitHeight = hit->PeakAmplitude();
        hitParams.hitStart = hit->StartTick(); //PeakTime() - 3. * hit->RMS();
        hitParams.hitEnd = hit->EndTick();     //PeakTime() + 3. * hit->RMS();

        lastEndTick = hitParams.hitEnd;

        roiHitParamsVec.emplace_back(hitParams);

      } //end loop over reco hits

      // Just in case (probably never called...)
      if (!roiHitParamsVec.empty()) hitParamsVec.push_back(roiHitParamsVec);

      size_t roiCount(0);

      for (const auto& roiHitParamsVec : hitParamsVec) {
        // Create a histogram here...
        double roiStart = roiHitParamsVec.front().hitStart;
        double roiStop = roiHitParamsVec.back().hitEnd;

        std::string funcString = "gaus(0)";
        std::string funcName = Form("hitshape_%05zu_c%02zu", size_t(channel), roiCount++);

        for (size_t idx = 1; idx < roiHitParamsVec.size(); idx++)
          funcString += "+gaus(" + std::to_string(3 * idx) + ")";

        // Include a baseline
        float baseline(0.);

        if (fFloatBaseline && !wireDataVec.empty()) baseline = wireDataVec.at(roiStart);

        funcString += "+" + std::to_string(baseline);

        fHitFunctionVec.emplace_back(
          std::make_unique<TF1>(funcName.c_str(), funcString.c_str(), roiStart, roiStop));
        TF1& hitFunc = *fHitFunctionVec.back().get();

        hitFunc.SetLineColor(fColorVec[imod % fColorVec.size()]);

        size_t idx(0);
        for (const auto& hitParams : roiHitParamsVec) {
          hitFunc.SetParameter(idx + 0, hitParams.hitHeight);
          hitFunc.SetParameter(idx + 1, hitParams.hitCenter);
          hitFunc.SetParameter(idx + 2, hitParams.hitSigma);

          TPolyLine& hitHeight = view2D.AddPolyLine(2, kBlack, 1, 1);

          hitHeight.SetPoint(0, hitParams.hitCenter, baseline);
          hitHeight.SetPoint(1, hitParams.hitCenter, hitParams.hitHeight + baseline);

          hitHeight.Draw("same");

          TPolyLine& hitSigma = view2D.AddPolyLine(2, kGray, 1, 1);

          hitSigma.SetPoint(
            0, hitParams.hitCenter - hitParams.hitSigma, 0.6 * hitParams.hitHeight + baseline);
          hitSigma.SetPoint(
            1, hitParams.hitCenter + hitParams.hitSigma, 0.6 * hitParams.hitHeight + baseline);

          hitSigma.Draw("same");

          idx += 3;
        }

        hitFunc.Draw("same");
      }
    } //end loop over HitFinding modules

    return;
  }

  DEFINE_ART_CLASS_TOOL(DrawGausHits)
}
