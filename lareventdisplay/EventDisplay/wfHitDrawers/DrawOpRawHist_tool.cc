////////////////////////////////////////////////////////////////////////
/// \file   DrawOpRawHist_tool.cc
/// \author tjyang@fnal.gov
/// Based on DrawRawHist_tool.cc
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWaveformDrawer.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"

#include "nuevdb/EventDisplayBase/EventHolder.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"

namespace evdb_tool {

  class DrawOpRawHist : public IWaveformDrawer {
  public:
    explicit DrawOpRawHist(const fhicl::ParameterSet& pset);

    ~DrawOpRawHist();

    void configure(const fhicl::ParameterSet& pset) override;
    void Fill(evdb::View2D&, raw::ChannelID_t&, float, float) override;
    void Draw(const std::string&, float, float) override;

    float getMaximum() const override { return fMaximum; };
    float getMinimum() const override { return fMinimum; };

  private:
    void BookHistogram(raw::ChannelID_t&, float, float);

    float fMaximum;
    float fMinimum;

    std::unique_ptr<TH1F> fRawDigitHist;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawOpRawHist::DrawOpRawHist(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }

  DrawOpRawHist::~DrawOpRawHist() {}

  void DrawOpRawHist::configure(const fhicl::ParameterSet& pset)
  {
    return;
  }

  void DrawOpRawHist::Fill(evdb::View2D& view2D,
                           raw::ChannelID_t& channel,
                           float lowBin,
                           float numTicks)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;

    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if (!event) return;

    // Handle histograms
    BookHistogram(channel, lowBin, numTicks);

    fMinimum = std::numeric_limits<float>::max();
    fMaximum = std::numeric_limits<float>::lowest();

    // Loop over the possible producers of RawDigits
    for (const auto& rawDataLabel : rawOpt->fOpRawDataLabels) {
      art::Handle<std::vector<raw::OpDetWaveform>> rawDigitVecHandle;
      event->getByLabel(rawDataLabel, rawDigitVecHandle);

      if (!rawDigitVecHandle.isValid()) continue;

      for (size_t rawDigitIdx = 0; rawDigitIdx < rawDigitVecHandle->size(); rawDigitIdx++) {
        art::Ptr<raw::OpDetWaveform> rawDigit(rawDigitVecHandle, rawDigitIdx);

        if (rawDigit->ChannelNumber() != channel) continue;

        auto const & waveform = rawDigit->Waveform();

        // recover the pedestal
        float pedestal = 0;
        
        TH1D *basehelp= new TH1D("basehelp","basehelp",2000, 1300,1800);
        for(size_t j=0; j<waveform.size(); j++){
          if(j<1000){
            basehelp->Fill(waveform[j]);
          }
        }
        int basebinmax = basehelp->GetMaximumBin();
        pedestal = basehelp->GetXaxis()->GetBinCenter(basebinmax);
        basehelp->Delete();   

        TH1F* histPtr = fRawDigitHist.get();

        for (size_t idx = 0; idx < waveform.size(); ++idx) {
          float signalVal = float(waveform[idx]) - pedestal;

          histPtr->Fill(float(idx) + 0.5, signalVal);
        }

        short minimumVal = *std::min_element(waveform.begin(), waveform.end());
        short maximumVal = *std::max_element(waveform.begin(), waveform.end());

        fMinimum = float(minimumVal) - pedestal;
        fMaximum = float(maximumVal) - pedestal;

        histPtr->SetLineColor(kBlack);

        // There is only one channel displayed so if here we are done
        break;
      }
    }

    return;
  }

  void DrawOpRawHist::Draw(const std::string& options, float maxLowVal, float maxHiVal)
  {
    TH1F* histPtr = fRawDigitHist.get();

    // Do we have valid limits to set?
    histPtr->SetMaximum(maxHiVal);
    histPtr->SetMinimum(maxLowVal);

    histPtr->Draw(options.c_str());

    return;
  }

  //......................................................................
  void DrawOpRawHist::BookHistogram(raw::ChannelID_t& channel, float startTick, float numTicks)
  {
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;
    art::ServiceHandle<geo::Geometry const> geo;

    // Get rid of the previous histograms
    if (fRawDigitHist.get()) fRawDigitHist.reset();

    // figure out the signal type for this plane, assume that
    // plane n in each TPC/cryostat has the same type
    geo::SigType_t sigType = geo->SignalType(channel);
    int numBins = numTicks;

    fRawDigitHist = std::make_unique<TH1F>(
                                           "fRAWQHisto", ";t [ticks];q [ADC]", numBins, startTick, startTick + numTicks);

    TH1F* histPtr = fRawDigitHist.get();

    histPtr->SetMaximum(cst->fRawQHigh[(size_t)sigType]);
    histPtr->SetMinimum(cst->fRawQLow[(size_t)sigType]);

    histPtr->SetLineColor(kBlack);
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

  DEFINE_ART_CLASS_TOOL(DrawOpRawHist)
}
