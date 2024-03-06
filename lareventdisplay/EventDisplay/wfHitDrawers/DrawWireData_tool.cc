////////////////////////////////////////////////////////////////////////
/// \file   DrawWireData_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Wire.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWaveformDrawer.h"

#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "nuevdb/EventDisplayBase/View2D.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "TPolyLine.h"

namespace evdb_tool {

  class DrawWireData : public IWaveformDrawer {
  public:
    explicit DrawWireData(const fhicl::ParameterSet& pset);

    ~DrawWireData();

    void configure(const fhicl::ParameterSet& pset) override;
    void Fill(evdb::View2D&, raw::ChannelID_t&, float, float) override;
    void Draw(const std::string&, float, float) override;

    float getMaximum() const override { return fMaximum; };
    float getMinimum() const override { return fMinimum; };

  private:
    float fMaximum;
    float fMinimum;

    std::vector<int> fColorMap;
  };

  //----------------------------------------------------------------------
  // Constructor.
  DrawWireData::DrawWireData(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }

  DrawWireData::~DrawWireData() {}

  void DrawWireData::configure(const fhicl::ParameterSet& pset)
  {
    fColorMap.push_back(kBlue);
    fColorMap.push_back(kMagenta);
    fColorMap.push_back(kBlack);
    fColorMap.push_back(kRed);

    return;
  }

  void DrawWireData::Fill(evdb::View2D& view2D,
                          raw::ChannelID_t& channel,
                          float lowBin,
                          float hiBin)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if (!event) return;

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      // Step one is to recover the hits for this label that match the input channel
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::Handle<std::vector<recob::Wire>> wireVecHandle;
      event->getByLabel(which, wireVecHandle);

      for (size_t wireIdx = 0; wireIdx < wireVecHandle->size(); wireIdx++) {
        art::Ptr<recob::Wire> wire(wireVecHandle, wireIdx);

        if (wire->Channel() != channel) continue;

        // Recover a full wire version of the deconvolved wire data
        // (the ROIs don't tend to display well)
        std::vector<float> signal = wire->Signal();

        TPolyLine& wireWaveform =
          view2D.AddPolyLine(signal.size(), fColorMap[imod % fColorMap.size()], 2, 1);

        for (size_t idx = 0; idx < signal.size(); idx++) {
          float bin = float(idx) + 0.5;

          if (bin >= lowBin && bin <= hiBin) wireWaveform.SetPoint(idx, bin, signal[idx]);
        }

        wireWaveform.Draw("same");
      }
    } //end loop over HitFinding modules

    return;
  }

  void DrawWireData::Draw(const std::string&, float, float)
  {
    return;
  }

  DEFINE_ART_CLASS_TOOL(DrawWireData)
}
