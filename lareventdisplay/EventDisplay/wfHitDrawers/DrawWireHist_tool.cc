////////////////////////////////////////////////////////////////////////
/// \file   DrawWireHist_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWaveformDrawer.h"

#include "nuevdb/EventDisplayBase/EventHolder.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "TH1F.h"

namespace evdb_tool {

  class DrawWireHist : public IWaveformDrawer {
  public:
    explicit DrawWireHist(const fhicl::ParameterSet& pset);

    ~DrawWireHist();

    void configure(const fhicl::ParameterSet& pset) override;
    void Fill(evdb::View2D&, raw::ChannelID_t&, float, float) override;
    void Draw(const std::string&, float, float) override;

    float getMaximum() const override { return fMaximum; };
    float getMinimum() const override { return fMinimum; };

  private:
    void BookHistogram(raw::ChannelID_t&, float, float);

    float fMaximum;
    float fMinimum;

    std::vector<int> fColorMap;
    std::unordered_map<std::string, std::unique_ptr<TH1F>> fRecoHistMap;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawWireHist::DrawWireHist(const fhicl::ParameterSet& pset) { configure(pset); }

  DrawWireHist::~DrawWireHist() {}

  void DrawWireHist::configure(const fhicl::ParameterSet& pset)
  {
    fColorMap.push_back(kBlue);
    fColorMap.push_back(kMagenta);
    fColorMap.push_back(kBlack);
    fColorMap.push_back(kRed);

    fRecoHistMap.clear();

    return;
  }

  void DrawWireHist::Fill(evdb::View2D& view2D,
                          raw::ChannelID_t& channel,
                          float lowBin,
                          float numTicks)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    // Check if we're supposed to draw raw hits at all
    if (rawOpt->fDrawRawDataOrCalibWires == 0) return;

    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if (!event) return;

    // Handle histograms
    BookHistogram(channel, lowBin, numTicks);

    fMinimum = std::numeric_limits<float>::max();
    fMaximum = std::numeric_limits<float>::lowest();

    int nWireLabels = 0;
    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      // Step one is to recover the hits for this label that match the input channel
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::Handle<std::vector<recob::Wire>> wireVecHandle;
      if (!event->getByLabel(which, wireVecHandle)) continue;
      ++nWireLabels;

      for (size_t wireIdx = 0; wireIdx < wireVecHandle->size(); wireIdx++) {
        art::Ptr<recob::Wire> wire(wireVecHandle, wireIdx);

        if (wire->Channel() != channel) continue;

        const std::vector<float>& signalVec = wire->Signal();

        TH1F* histPtr = fRecoHistMap.at(which.encode()).get();

        for (size_t idx = 0; idx < signalVec.size(); idx++) {
          histPtr->Fill(float(idx) + 0.5, signalVec[idx]);

          fMinimum = std::min(fMinimum, signalVec[idx]);
          fMaximum = std::max(fMaximum, signalVec[idx]);
        }

        histPtr->SetLineColor(fColorMap.at((nWireLabels - 1) % recoOpt->fWireLabels.size()));

        // There is only one channel displayed so if here we are done
        break;
      }
    } //end loop over HitFinding modules

    return;
  }

  void DrawWireHist::Draw(const std::string& options, float maxLowVal, float maxHiVal)
  {
    for (const auto& histMap : fRecoHistMap) {
      TH1F* histPtr = histMap.second.get();

      // Set the limits
      histPtr->SetMaximum(maxHiVal);
      histPtr->SetMinimum(maxLowVal);

      histPtr->Draw(options.c_str());
    }

    return;
  }

  //......................................................................
  void DrawWireHist::BookHistogram(raw::ChannelID_t& channel, float startTick, float numTicks)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;
    art::ServiceHandle<evd::RawDrawingOptions const> drawopt;
    art::ServiceHandle<geo::Geometry const> geo;

    // Get rid of the previous histograms
    fRecoHistMap.clear();

    // Now add a histogram for each of the wire labels
    for (auto& tag : recoOpt->fWireLabels) {
      // figure out the signal type for this plane, assume that
      // plane n in each TPC/cryostat has the same type
      geo::SigType_t sigType = geo->SignalType(channel);
      std::string tagString(tag.encode());
      int numBins = numTicks;

      fRecoHistMap[tagString] = std::make_unique<TH1F>(
        "fCALTQHisto", ";t [ticks];q [ADC]", numBins, startTick, startTick + numTicks);

      TH1F* histPtr = fRecoHistMap.at(tagString).get();

      histPtr->SetMaximum(cst->fRecoQHigh[(size_t)sigType]);
      histPtr->SetMinimum(cst->fRecoQLow[(size_t)sigType]);

      histPtr->SetLineColor(kBlue);
      histPtr->SetLineWidth(1);

      histPtr->GetXaxis()->SetLabelSize(0.10);   // was 0.15
      histPtr->GetXaxis()->SetLabelOffset(0.01); // was 0.00
      histPtr->GetXaxis()->SetTitleSize(0.10);   // was 0.15
      histPtr->GetXaxis()->SetTitleOffset(0.60); // was 0.80

      histPtr->GetYaxis()->SetLabelSize(0.10);    // was 0.15
      histPtr->GetYaxis()->SetLabelOffset(0.002); // was 0.00
      histPtr->GetYaxis()->SetTitleSize(0.10);    // was 0.15
      histPtr->GetYaxis()->SetTitleOffset(0.16);  // was 0.80
    }
  }

  DEFINE_ART_CLASS_TOOL(DrawWireHist)
}
