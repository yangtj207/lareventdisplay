///
/// \file    TQPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
///
#include "lareventdisplay/EventDisplay/TQPad.h"
#include "TH1F.h"
#include "TPad.h"

#include "art/Utilities/make_tool.h"
#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWFHitDrawer.h"
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWaveformDrawer.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "nuevdb/EventDisplayBase/View2D.h"

// C/C++ standard libraries
#include <algorithm> // std::min(), std::max()

namespace evd {

  static const int kRAW = 0;
  static const int kCALIB = 1;
  //static const int kRAWCALIB = 2;
  static const int kQ = 0;
  static const int kTQ = 1;

  //......................................................................

  TQPad::TQPad(const char* nm,
               const char* ti,
               double x1,
               double y1,
               double x2,
               double y2,
               const char* opt,
               unsigned int plane,
               unsigned int wire)
    : DrawingPad(nm, ti, x1, y1, x2, y2), fWire(wire), fPlane(plane), fFrameHist(0)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    unsigned int planes = geo->Nplanes();

    this->Pad()->cd();

    this->Pad()->SetLeftMargin(0.050);
    this->Pad()->SetRightMargin(0.050);

    this->Pad()->SetTopMargin(0.005);
    this->Pad()->SetBottomMargin(0.110);

    // there has to be a better way of doing this that does
    // not have a case for each number of planes in a detector
    if (planes == 2 && fPlane > 0) {
      this->Pad()->SetTopMargin(0.110);
      this->Pad()->SetBottomMargin(0.010);
    }
    else if (planes > 2) {
      if (fPlane == 1) {
        this->Pad()->SetTopMargin(0.005);
        this->Pad()->SetBottomMargin(0.010);
      }
      else if (fPlane == 2) {
        this->Pad()->SetTopMargin(0.110);
        this->Pad()->SetBottomMargin(0.010);
      }
    }

    std::string opts(opt);
    if (opts == "TQ") {
      fTQ = kTQ;
      // BB adjust the vertical spacing
      this->Pad()->SetTopMargin(0);
      this->Pad()->SetBottomMargin(0.2);
    }
    if (opts == "Q") { fTQ = kQ; }

    this->BookHistogram();
    fView = new evdb::View2D();

    art::ServiceHandle<evd::RawDrawingOptions const> rawOptions;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOptions;

    fHitDrawerTool = art::make_tool<evdb_tool::IWFHitDrawer>(recoOptions->fHitDrawerParams);
    fRawDigitDrawerTool =
      art::make_tool<evdb_tool::IWaveformDrawer>(rawOptions->fRawDigitDrawerParams);
    fWireDrawerTool = art::make_tool<evdb_tool::IWaveformDrawer>(recoOptions->fWireDrawerParams);
  }

  //......................................................................

  TQPad::~TQPad()
  {
    if (fView) {
      delete fView;
      fView = 0;
    }
    if (fFrameHist) {
      delete fFrameHist;
      fFrameHist = 0;
    }
  }

  //......................................................................
  void
  TQPad::Draw()
  {
    art::ServiceHandle<evd::RawDrawingOptions const> drawopt;

    //grab the singleton with the event
    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    if (!evt) return;

    art::ServiceHandle<geo::Geometry const> geoSvc;

    fPad->Clear();
    fPad->cd();

    // Note this handles drawing waveforms for both SP and DP where the difference is handled by the tools
    if (fTQ == kTQ) {
      // Recover a channel number from current information
      raw::ChannelID_t channel =
        geoSvc->PlaneWireToChannel(fPlane, fWire, drawopt->fTPC, drawopt->fCryostat);

      // Call the tools to fill the histograms for RawDigits and Wire data
      fRawDigitDrawerTool->Fill(
        *fView, channel, this->RawDataDraw()->StartTick(), this->RawDataDraw()->TotalClockTicks());
      fWireDrawerTool->Fill(
        *fView, channel, this->RawDataDraw()->StartTick(), this->RawDataDraw()->TotalClockTicks());

      // Vertical limits set for the enclosing histogram, then draw it with axes only
      float maxLowVal = std::min(fRawDigitDrawerTool->getMinimum(), fWireDrawerTool->getMinimum());
      float maxHiVal = std::max(fRawDigitDrawerTool->getMaximum(), fWireDrawerTool->getMaximum());

      if (drawopt->fDrawRawDataOrCalibWires == kCALIB) {
        maxLowVal = fWireDrawerTool->getMinimum();
        maxHiVal = fWireDrawerTool->getMaximum();
      }

      if (maxLowVal < std::numeric_limits<float>::max())
        maxLowVal -= 5.;
      else
        maxLowVal = -10.;
      if (maxHiVal > std::numeric_limits<float>::lowest())
        maxHiVal += 5.;
      else
        maxHiVal = 10.;

      fFrameHist->SetMaximum(maxHiVal);
      fFrameHist->SetMinimum(maxLowVal);
      fFrameHist->Draw("AXIS");

      // draw with histogram style, only (square) lines, no errors
      static const std::string defaultDrawOptions = "HIST same";

      // Draw the desired histograms
      // If its not just the raw hists then we output the wire histograms
      if (drawopt->fDrawRawDataOrCalibWires != kRAW) {
        fWireDrawerTool->Draw(defaultDrawOptions.c_str(), maxLowVal, maxHiVal);

        fHitDrawerTool->Draw(*fView, channel);
      }

      // Likewise, if it is not just the calib hists then we output the raw histogram
      if (drawopt->fDrawRawDataOrCalibWires != kCALIB)
        fRawDigitDrawerTool->Draw(defaultDrawOptions.c_str(), maxLowVal, maxHiVal);

      // This is a remnant from a time long past...
      fFrameHist->SetTitleOffset(0.2, "Y");
    } // end if fTQ == kTQ

    // I am not sure what the block below is trying to do... I don't see where the hists are actually filled.
    // ** remove this for now until someone can explain what it is **
    //    else if(fTQ == kQ && fTQ == -1)
    //    {
    //        // figure out the signal type for this plane, assume that
    //        // plane n in each TPC/cryostat has the same type
    //        geo::PlaneID planeid(drawopt->CurrentTPC(), fPlane);
    //        geo::SigType_t sigType = geoSvc->SignalType(planeid);
    //
    //        art::ServiceHandle<evd::ColorDrawingOptions const> cst;
    //
    //        TH1F *hist;
    //
    //        int ndiv = 0;
    //        if(drawopt->fDrawRawDataOrCalibWires != kCALIB){
    //            hist = fRawHisto;
    //            hist->SetMinimum(cst->fRawQLow [(size_t)sigType]);
    //            hist->SetMaximum(cst->fRawQHigh[(size_t)sigType]);
    //            ndiv = cst->fRawDiv[(size_t)sigType];
    //        }
    //        if(drawopt->fDrawRawDataOrCalibWires == kCALIB){
    //            hist = fRecoHisto;
    //            hist->SetMinimum(cst->fRecoQLow [(size_t)sigType]);
    //            hist->SetMaximum(cst->fRecoQHigh[(size_t)sigType]);
    //            ndiv = cst->fRecoDiv[(size_t)sigType];
    //        }
    //
    //        hist->SetLabelSize(0, "X");
    //        hist->SetLabelSize(0, "Y");
    //        hist->SetTickLength(0, "X");
    //        hist->SetTickLength(0, "Y");
    //        hist->Draw("pY+");
    //
    //        //
    //        // Use this to fill the histogram with colors from the color scale
    //        //
    //        double x1, x2, y1, y2;
    //        x1 = 0.;
    //        x2 = 1.;
    //
    //        for(int i = 0; i < ndiv; ++i){
    //            y1 = hist->GetMinimum() + i*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);
    //            y2 = hist->GetMinimum() + (i + 1)*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);
    //
    //            int c = 1;
    //            if (drawopt->fDrawRawDataOrCalibWires==kRAW) {
    //                c = cst->RawQ(sigType).GetColor(0.5*(y1+y2));
    //            }
    //            if (drawopt->fDrawRawDataOrCalibWires!=kRAW) {
    //                c= cst->CalQ(sigType).GetColor(0.5*(y1+y2));
    //            }
    //
    //            TBox& b = fView->AddBox(x1,y1,x2,y2);
    //            b.SetFillStyle(1001);
    //            b.SetFillColor(c);
    //            b.Draw();
    //        } // end loop over Q histogram bins
    //
    //        hist->Draw("same");
    //    } // end if fTQ == kQ

    return;
  }

  //......................................................................
  void
  TQPad::BookHistogram()
  {
    if (fFrameHist) {
      delete fFrameHist;
      fFrameHist = 0;
    }

    art::ServiceHandle<evd::ColorDrawingOptions const> cst;
    art::ServiceHandle<evd::RawDrawingOptions const> drawopt;

    // figure out the signal type for this plane, assume that
    // plane n in each TPC/cryostat has the same type
    geo::PlaneID planeid(drawopt->CurrentTPC(), fPlane);
    art::ServiceHandle<geo::Geometry const> geo;
    geo::SigType_t sigType = geo->SignalType(planeid);

    /// \todo decide if ndivraw and ndivreco are useful
    double qxloraw = cst->fRawQLow[(size_t)sigType];
    double qxhiraw = cst->fRawQHigh[(size_t)sigType];
    double tqxlo = 1. * this->RawDataDraw()->StartTick();
    double tqxhi = 1. * this->RawDataDraw()->TotalClockTicks();

    switch (fTQ) {
    case kQ:
      fFrameHist = new TH1F("fFrameHist", ";t [ticks];[ADC]", 2, 0., 1.);
      fFrameHist->SetMaximum(qxhiraw);
      fFrameHist->SetMinimum(qxloraw);
      break; // kQ
    case kTQ:
      fFrameHist = new TH1F("fFrameHist", ";t [ticks];q [ADC]", (int)tqxhi, tqxlo, tqxhi + tqxlo);
      break;
    default: throw cet::exception("TQPad") << __func__ << ": unexpected quantity #" << fTQ << "\n";
    } //end if fTQ == kTQ

    // Set the label, title size and offsets
    // Note this is the base histogram so this control these for both the raw and wire histograms
    fFrameHist->GetXaxis()->SetLabelSize(0.10);
    fFrameHist->GetXaxis()->SetLabelOffset(0.00);
    fFrameHist->GetXaxis()->SetTitleSize(0.10);
    fFrameHist->GetXaxis()->SetTitleOffset(0.80);

    fFrameHist->GetYaxis()->SetLabelSize(0.10);
    fFrameHist->GetYaxis()->SetLabelOffset(0.01);
    fFrameHist->GetYaxis()->SetTitleSize(0.10);
    fFrameHist->GetYaxis()->SetTitleOffset(0.80);
  }

}
//////////////////////////////////////////////////////////////////////////
