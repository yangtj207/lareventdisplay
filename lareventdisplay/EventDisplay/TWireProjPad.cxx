////////////////////////////////////////////////////////////////////////
///
/// \file    TWireProjPad.cxx
/// \brief   Drawing pad for X-Z or Y-Z projections of events
/// \author  messier@indiana.edu
///
////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "TCanvas.h"
#include "TClass.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TList.h"
#include "TPad.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/PxUtils.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/Style.h"
#include "lareventdisplay/EventDisplay/TWireProjPad.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "nuevdb/EventDisplayBase/View2D.h"

#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {

  template <typename Stream>
  void
  DumpPad(Stream&& log, TVirtualPad* pPad)
  {
    if (!pPad) {
      log << "pad not available";
      return;
    }

    log << pPad->IsA()->GetName() << "[" << ((void*)pPad) << "](\"" << pPad->GetName() << "\")";
    TFrame const* pFrame = pPad->GetFrame();
    if (pFrame) {
      double const low_wire = pFrame->GetX1(), high_wire = pFrame->GetX2();
      double const low_tdc = pFrame->GetY1(), high_tdc = pFrame->GetY2();
      double const wire_pixels = pPad->XtoAbsPixel(high_wire) - pPad->XtoAbsPixel(low_wire);
      double const tdc_pixels = -(pPad->YtoAbsPixel(high_tdc) - pPad->YtoAbsPixel(low_tdc));
      log << " has frame spanning wires " << low_wire << "-" << high_wire << " and TDC " << low_tdc
          << "-" << high_tdc << " in a window " << wire_pixels << "x" << tdc_pixels << " pixel big";
    }
    else {
      log << " has no frame";
    }

  } // DumpPad()

  [[maybe_unused]] void
  DumpPadsInCanvas(TVirtualPad* pPad, std::string caller, std::string msg = "")
  {
    mf::LogDebug log(caller);
    if (!msg.empty()) log << msg << ": ";
    if (!pPad) {
      log << "pad not available";
      return;
    }

    DumpPad(log, pPad);

    TCanvas const* pCanvas = pPad->GetCanvas();
    log << "\nCanvas is: (TCanvas*) (" << ((void*)pPad->GetCanvas()) << ") with "
        << pCanvas->GetListOfPrimitives()->GetSize() << " primitives and the following pads:";
    TIterator* pIter = pCanvas->GetListOfPrimitives()->MakeIterator();
    TObject const* pObject;
    while ((pObject = pIter->Next())) {
      if (!pObject->InheritsFrom(TVirtualPad::Class())) continue;
      log << "\n  " << ((pObject == pPad) ? '*' : '-') << "  ";
      DumpPad(log, (TVirtualPad*)pObject);
    }
    log << "\n";
    delete pIter;
  } // DumpPadsInCanvas()

} // local namespace

namespace evd {

  ///
  /// Create a pad showing a single X-Z or Y-Z projection of the detector
  /// \param nm : Name of the pad
  /// \param ti : Title of the pad
  /// \param x1 : Location of left  edge of pad (0-1)
  /// \param x2 : Location of right edge of pad (0-1)
  /// \param y1 : Location of bottom edge of pad (0-1)
  /// \param y2 : Location of top    edge of pad (0-1)
  /// \param plane : plane number of view
  ///
  TWireProjPad::TWireProjPad(const char* nm,
                             const char* ti,
                             double x1,
                             double x2,
                             double y1,
                             double y2,
                             unsigned int plane)
    : DrawingPad(nm, ti, x1, x2, y1, y2), fPlane(plane)
  {
    fCurrentZoom.resize(4);

    art::ServiceHandle<geo::Geometry const> geo;

    this->Pad()->cd();

    this->Pad()->SetLeftMargin(0.070);
    this->Pad()->SetRightMargin(0.010);

    // how many planes in the detector and
    // which plane is this one?

    unsigned int planes = geo->Nplanes();
    this->Pad()->SetTopMargin(0.005);
    this->Pad()->SetBottomMargin(0.110);

    // there has to be a better way of doing this that does
    // not have a case for each number of planes in a detector
    if (planes == 2 && fPlane > 0) {
      this->Pad()->SetTopMargin(0.110);
      this->Pad()->SetBottomMargin(0.005);
    }
    else if (planes > 2) {
      if (fPlane == 1) {
        this->Pad()->SetTopMargin(0.055);
        this->Pad()->SetBottomMargin(0.055);
      }
      else if (fPlane == 2) {
        this->Pad()->SetTopMargin(0.110);
        this->Pad()->SetBottomMargin(0.005);
      }
    }

    TString planeNo = "fTWirePlane";
    planeNo += fPlane;

    // picking the information from the current TPC
    art::ServiceHandle<evd::RawDrawingOptions const> rawopt;
    auto const signalType = geo->SignalType({rawopt->CurrentTPC(), fPlane});
    TString xtitle = ";Induction Wire;t (tdc)";
    if (signalType == geo::kCollection) xtitle = ";Collection Wire;t (tdc)";

    unsigned int const nWires = geo->Nwires(fPlane);
    unsigned int const nTicks = RawDataDraw()->TotalClockTicks();

    fXLo = -0.005 * (nWires - 1);
    fXHi = 1.005 * (nWires - 1);
    fYLo = 0.990 * (unsigned int)(this->RawDataDraw()->StartTick());
    fYHi = 1.005 * std::min((unsigned int)(this->RawDataDraw()->StartTick() + nTicks), nTicks);

    fOri = rawopt->fAxisOrientation;
    if (fOri > 0) {
      fYLo = -0.005 * (nWires - 1);
      fYHi = 1.005 * (nWires - 1);
      fYLo = 0.990 * (unsigned int)(this->RawDataDraw()->StartTick());
      fYHi = 1.005 * std::min((unsigned int)(this->RawDataDraw()->StartTick() + nTicks), nTicks);
      fXLo = -0.005 * nTicks;
      fXHi = 1.010 * nTicks;
      xtitle = ";t (tdc);InductionWire";
      if (signalType == geo::kCollection) xtitle = ";t (tdc);Collection Wire";
    }

    // make the range of the histogram be the biggest extent
    // in both directions and then use SetRangeUser() to shrink it down
    // that will allow us to change the axes on the fly
    double min = std::min(fXLo, fYLo);
    double max = std::max(fXHi, fYHi);

    fHisto = new TH1F(*(fPad->DrawFrame(min, min, max, max)));

    fHisto->SetTitleOffset(0.5, "Y");
    fHisto->SetTitleOffset(0.75, "X");
    SetZoomRange(fXLo, fXHi, fYLo, fYHi);
    fHisto->GetYaxis()->SetLabelSize(0.05);
    fHisto->GetYaxis()->CenterTitle();
    fHisto->GetXaxis()->SetLabelSize(0.05);
    fHisto->GetXaxis()->CenterTitle();
    fHisto->Draw("AB");

    fView = new evdb::View2D();
  }

  //......................................................................
  TWireProjPad::~TWireProjPad()
  {
    if (fHisto) {
      delete fHisto;
      fHisto = 0;
    }
    if (fView) {
      delete fView;
      fView = 0;
    }
  }

  //......................................................................
  void
  TWireProjPad::Draw(const char* opt)
  {
    // DumpPadsInCanvas(fPad, "TWireProjPad", "Draw()");
    MF_LOG_DEBUG("TWireProjPad") << "Started to draw plane " << fPlane;

    ///\todo: Why is kSelectedColor hard coded?
    int kSelectedColor = 4;
    fView->Clear();

    // grab the singleton holding the art::Event
    art::Event const* evtPtr = evdb::EventHolder::Instance()->GetEvent();
    if (evtPtr) {
      auto const& evt = *evtPtr;
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
      art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

      this->SimulationDraw()->MCTruthVectors2D(evt, fView, fPlane);

      // the 2D pads have too much detail to be rendered on screen;
      // to act smarter, RawDataDrawer needs to know the range being plotted
      this->RawDataDraw()->ExtractRange(fPad, &GetCurrentZoom());
      this->RawDataDraw()->RawDigit2D(
        evt, detProp, fView, fPlane, GetDrawOptions().bZoom2DdrawToRoI);

      this->RecoBaseDraw()->Wire2D(evt, fView, fPlane);
      this->RecoBaseDraw()->Hit2D(evt, detProp, fView, fPlane);

      if (recoOpt->fUseHitSelector)
        this->RecoBaseDraw()->Hit2D(
          this->HitSelectorGet()->GetSelectedHits(fPlane), kSelectedColor, fView, true);

      this->RecoBaseDraw()->Slice2D(evt, detProp, fView, fPlane);
      this->RecoBaseDraw()->Cluster2D(evt, clockData, detProp, fView, fPlane);
      this->RecoBaseDraw()->EndPoint2D(evt, fView, fPlane);
      this->RecoBaseDraw()->Prong2D(evt, clockData, detProp, fView, fPlane);
      this->RecoBaseDraw()->Vertex2D(evt, detProp, fView, fPlane);
      this->RecoBaseDraw()->Seed2D(evt, detProp, fView, fPlane);
      this->RecoBaseDraw()->OpFlash2D(evt, clockData, detProp, fView, fPlane);
      this->RecoBaseDraw()->Event2D(evt, fView, fPlane);
      this->RecoBaseDraw()->DrawTrackVertexAssns2D(evt, clockData, detProp, fView, fPlane);

      UpdatePad();
    } // if (evt)

    ClearandUpdatePad();

    // check if we need to swap the axis ranges
    art::ServiceHandle<evd::RawDrawingOptions const> rawopt;
    if (fOri != rawopt->fAxisOrientation) {
      fOri = rawopt->fAxisOrientation;
      double max = fXHi;
      double min = fXLo;
      fXHi = fYHi;
      fXLo = fYLo;
      fYHi = max;
      fYLo = min;

      SetZoomRange(fXLo, fXHi, fYLo, fYHi);

      TString xtitle = fHisto->GetXaxis()->GetTitle();
      fHisto->GetXaxis()->SetTitle(fHisto->GetYaxis()->GetTitle());
      fHisto->GetYaxis()->SetTitle(xtitle);
    }

    if (fPlane > 0)
      fHisto->Draw("X+");
    else
      fHisto->Draw("");

    // Check if we should zoom the displays;
    // if there is no event, we have no clue about the region of interest
    // and therefore we don't touch anything
    if (opt == 0 && evtPtr) { this->ShowFull(); }

    MF_LOG_DEBUG("TWireProjPad") << "Started rendering plane " << fPlane;

    fView->Draw();

    MF_LOG_DEBUG("TWireProjPad") << "Drawing of plane " << fPlane << " completed";
  }

  //......................................................................
  void
  TWireProjPad::ClearHitList()
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    if (recoOpt->fUseHitSelector) {
      this->HitSelectorGet()->ClearHitList(fPlane);
      this->Draw();
    }
  }

  //......................................................................
  // the override parameter is needed to unzoom to full range when the fAutoZoomInterest is on.

  void
  TWireProjPad::ShowFull(int override)
  {
    // x values are wire numbers, y values are ticks of the clock
    int xmin = fXLo;
    int xmax = fXHi;
    int ymax = fYHi;
    int ymin = fYLo;

    art::ServiceHandle<evd::EvdLayoutOptions const> evdlayoutopt;
    art::ServiceHandle<evd::RawDrawingOptions const> rawopt;

    if (GetDrawOptions().bZoom2DdrawToRoI && !override) {
      int test = 0;
      if (rawopt->fDrawRawDataOrCalibWires == 0)
        test = RawDataDraw()->GetRegionOfInterest((int)fPlane, xmin, xmax, ymin, ymax);
      else
        test = RecoBaseDraw()->GetRegionOfInterest((int)fPlane, xmin, xmax, ymin, ymax);

      if (test != 0) return;
    }

    SetZoomRange(xmin, xmax, ymin, ymax);
  }

  //......................................................................
  void
  TWireProjPad::GetWireRange(int* i1, int* i2) const
  {
    if (fOri < 1) {
      *i1 = fHisto->GetXaxis()->GetFirst();
      *i2 = fHisto->GetXaxis()->GetLast();
    }
    else {
      *i1 = fHisto->GetYaxis()->GetFirst();
      *i2 = fHisto->GetYaxis()->GetLast();
    }
  }

  //......................................................................
  // Set the X axis range only
  //
  void
  TWireProjPad::SetWireRange(int i1, int i2)
  {
    if (fOri < 1) { fHisto->GetXaxis()->SetRange(i1, i2); }
    else {
      fHisto->GetYaxis()->SetRange(i1, i2);
    }
    fCurrentZoom[0] = i1;
    fCurrentZoom[1] = i2;
  }

  //......................................................................
  // Set the visible range of the wire / time view
  //
  void
  TWireProjPad::SetZoomRange(int i1, int i2, int y1, int y2)
  {
    MF_LOG_DEBUG("TWireProjPad") << "SetZoomRange(" << i1 << ", " << i2 << ", " << y1 << ", " << y2
                                 << ") on plane #" << fPlane;

    fHisto->GetXaxis()->SetRangeUser(i1, i2);
    fHisto->GetYaxis()->SetRangeUser(y1, y2);
    fCurrentZoom[0] = i1;
    fCurrentZoom[1] = i2;
    fCurrentZoom[2] = y1;
    fCurrentZoom[3] = y2;
  }

  //......................................................................
  // Set the visible range of the wire / time view from the view
  //
  void
  TWireProjPad::SetZoomFromView()
  {
    TAxis const& xaxis = *(fHisto->GetXaxis());
    fCurrentZoom[0] = xaxis.GetBinLowEdge(xaxis.GetFirst());
    fCurrentZoom[1] = xaxis.GetBinUpEdge(xaxis.GetLast());
    fCurrentZoom[2] = fHisto->GetMinimum();
    fCurrentZoom[3] = fHisto->GetMaximum();
    MF_LOG_DEBUG("TWireProjPad") << "Zoom set to wires (" << fCurrentZoom[0] << "; "
                                 << fCurrentZoom[1] << " ), tick (" << fCurrentZoom[2] << "; "
                                 << fCurrentZoom[3] << ") for plane #" << fPlane;
  } // TWireProjPad::SetZoomFromView()
  //......................................................................
  void
  TWireProjPad::SaveHitList(double i1,
                            double i2,
                            double y1,
                            double y2,
                            double distance,
                            const char* zoom_opt,
                            bool good_plane)
  {
    const art::Event* evtPtr = evdb::EventHolder::Instance()->GetEvent();
    if (evtPtr) {
      auto const& evt = *evtPtr;
      art::ServiceHandle<evd::RecoDrawingOptions const> recoopt;
      if (recoopt->fUseHitSelector) {
        this->HitSelectorGet()->SaveHits(evt, fPlane, i1, i2, y1, y2, distance, good_plane);
        this->Draw(zoom_opt);
      }
    }
  }

  /////////////////////////////////////////////////
  // Pass the seed list onwards to InfoTransfer
  //
  double
  TWireProjPad::SaveSeedList(std::vector<util::PxLine> seedlines, double distance)
  {
    double KineticEnergy = util::kBogusD;
    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    if (evt) {
      art::ServiceHandle<evd::RecoDrawingOptions const> recoopt;
      if (recoopt->fUseHitSelector)
        KineticEnergy = this->HitSelectorGet()->SaveSeedLines(*evt, seedlines, distance);
    }
    return KineticEnergy;
  }

  //......................................................................
  void
  TWireProjPad::SelectOneHit(double x, double y, const char* zoom_opt)
  {

    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    if (evt) {
      art::ServiceHandle<evd::RecoDrawingOptions const> recoopt;
      if (recoopt->fUseHitSelector) {
        this->HitSelectorGet()->ChangeHit(*evt, fPlane, x, y);
        this->Draw(zoom_opt);
      }
    }
  }

  //......................................................................
  void
  TWireProjPad::ClearandUpdatePad()
  {
    fPad->Clear();
    this->UpdatePad();

    return;
  }

  //......................................................................
  void
  TWireProjPad::UpdatePad()
  {
    fPad->cd();
    fPad->Modified();
    fPad->Update();
    fPad->GetFrame()->SetBit(TPad::kCannotMove, true);
    fPad->SetBit(TPad::kCannotMove, true);

    return;
  }

  //......................................................................
  void
  TWireProjPad::DrawLinesinView(std::vector<util::PxLine> lines,
                                bool deleting,
                                const char* zoom_opt)
  {
    fPad->cd();
    if (deleting) {
      fPad->Clear();
      this->Draw(zoom_opt);
    }
    else {
      fView->Clear();
      fView->Draw();
    }

    mf::LogVerbatim("TWireProjPad") << "Drawing " << lines.size() << " lines";

    for (size_t is = 0; is < lines.size(); ++is) {
      if (fPlane != lines[is].plane) continue;

      TLine& l = fView->AddLine(lines[is].w0, lines[is].t0, lines[is].w1, lines[is].t1);

      fView->Draw();
      evd::Style::FromPDG(l, 11);
    }

    fView->Draw();
    UpdatePad();
    fView->Draw();

    return;
  }

} // namespace
////////////////////////////////////////////////////////////////////////
