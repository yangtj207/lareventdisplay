///
/// \file    Ortho3DPad.cxx
/// \brief   Drawing pad showing an orthographic rendering of 3D objects
///          in the detector
/// \author  greenlee@fnal.gov
///

#include "TPad.h"
#include "TH1F.h"
#include "TBox.h"
#include "TPolyMarker.h"
#include "TGNumberEntry.h"

#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "EventDisplay/Ortho3DPad.h"
#include "EventDisplay/SimulationDrawer.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/EventHolder.h"
#include "Geometry/Geometry.h"

///
/// Create a pad to show an orthographic rendering of 3D objcts.
/// @param name : Name of the pad
/// @param title : Title of the pad
/// @param proj : Choose orthographic projection
/// @param x1 : Location of left  edge of pad (0-1)
/// @param x2 : Location of right edge of pad (0-1)
/// @param y1 : Location of bottom edge of pad (0-1)
/// @param y2 : Location of top    edge of pad (0-1)
///
evd::Ortho3DPad::Ortho3DPad(const char* name, const char* title,
			    evd::OrthoProj_t proj,
			    double x1, double y1,
			    double x2, double y2) :
  DrawingPad(name, title, x1, y1, x2, y2),
  fHisto(0),
  fProj(proj),
  fXLo(0.),
  fXHi(0.),
  fYLo(0.),
  fYHi(0.),
  fMSize(0.25),
  fMSizeEntry(0),
  fPress(false),
  fBoxDrawn(false),
  fPressPx(0),
  fPressPy(0),
  fCurrentPx(0),
  fCurrentPy(0),
  fPressX(0.),
  fPressY(0.),
  fReleaseX(0.),
  fReleaseY(0.)
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;

  // Set up pad.

  Pad()->SetBit(kCannotPick);
  Pad()->SetBit(TPad::kCannotMove);
  Pad()->Draw();
  Pad()->cd();
  Pad()->SetLeftMargin (0.080);
  Pad()->SetRightMargin (0.010);
  Pad()->SetTopMargin (0.010);
  Pad()->SetBottomMargin (0.10);

  // Define histogram boundaries (cm).
  // For now only draw tpc=0.
  switch (proj) {
    case evd::kXY:
      fXLo = 0.;
      fXHi = 2.*geo->DetHalfWidth();
      fYLo = -geo->DetHalfHeight();
      fYHi = geo->DetHalfHeight();
      break;
    case evd::kXZ:
      fXLo = 0.;
      fXHi = geo->DetLength();
      fYLo = 0.;
      fYHi = 2.*geo->DetHalfWidth();
      break;
    case evd::kYZ:
      fXLo = 0.;
      fXHi = geo->DetLength();
      fYLo = -geo->DetHalfHeight();
      fYHi = geo->DetHalfHeight();
      break;
    default:
      throw cet::exception("Ortho3DPad")
        << __func__ << ": unwknow projection " << ((int) proj);
  } // switch

  // Make enclosing histogram.

  fHisto = new TH1F(*(Pad()->DrawFrame(fXLo, fYLo, fXHi, fYHi)));
  fHisto->SetBit(kCannotPick);
  fHisto->SetBit(TPad::kCannotMove);
  fHisto->SetTitleOffset(1.,"Y");
  fHisto->SetTitleOffset(1.,"X");
  fHisto->GetXaxis()->SetLabelSize(0.04);
  fHisto->GetXaxis()->SetTitleSize(0.04);
  switch (proj) {
    case evd::kXY:
      fHisto->GetXaxis()->SetTitle("x (cm)");
      fHisto->GetYaxis()->SetTitle("y (cm)");
      break;
    case evd::kXZ:
      fHisto->GetXaxis()->SetTitle("z (cm)");
      fHisto->GetYaxis()->SetTitle("x (cm)");
      break;
    case evd::kYZ:
      fHisto->GetXaxis()->SetTitle("z (cm)");
      fHisto->GetYaxis()->SetTitle("y (cm)");
      break;
    default:
      throw cet::exception("Ortho3DPad")
        << __func__ << ": unexpected flow (projection: " << ((int) proj) << ")";
  } // switch
  fHisto->GetXaxis()->CenterTitle();
  fHisto->GetYaxis()->SetLabelSize(0.04);
  fHisto->GetYaxis()->SetTitleSize(0.04);
  fHisto->GetYaxis()->CenterTitle();
  fHisto->Draw("AB");

  fView = new evdb::View2D();

  // Install mouse event handler.

  std::ostringstream ostr;
  ostr << "evd::Ortho3DPad::MouseEvent((evd::Ortho3DPad*)" << this << ")";
  fPad->AddExec("getmousezoom", ostr.str().c_str());
}

//......................................................................
// Destructor.

evd::Ortho3DPad::~Ortho3DPad() 
{
  if (fHisto) { delete fHisto; fHisto = nullptr; }
  if (fView) { delete fView; fView = nullptr; }
}

//......................................................................
// Draw selected objects.

void evd::Ortho3DPad::Draw(const char* /*opt*/)
{
  fPad->Clear();
  fView->Clear();

  // Remove zoom.

  UnZoom(false);

  // grab the event from the singleton

  const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();

  // Insert graphic objects into fView collection.

  if(evt)
    {
      SimulationDraw()->MCTruthOrtho(*evt, fProj, fMSize, fView);
      RecoBaseDraw()->SpacePointOrtho(*evt, fProj, fMSize, fView);
      RecoBaseDraw()->ProngOrtho(*evt, fProj, fMSize, fView);
      RecoBaseDraw()->SeedOrtho(*evt, fProj, fView);
    }
  // Draw objects on pad.

  fPad->cd();
  fHisto->Draw("X-");
  fView->Draw();
  fPad->Modified();
  fPad->Update();
}

//......................................................................
// Set zoom region.  Limits are specified in user coordinates.

void evd::Ortho3DPad::SetZoom(double xlo, double ylo,
			      double xhi, double yhi,
			      bool update)
{
  fHisto->GetXaxis()->SetRangeUser(xlo, xhi);
  fHisto->GetYaxis()->SetRangeUser(ylo, yhi);  
  fPad->Modified();
  if(update)
    fPad->Update();
}

//......................................................................
// Set zoom region to full range.

void evd::Ortho3DPad::UnZoom(bool update)
{
  fHisto->GetXaxis()->SetRangeUser(fXLo, fXHi);
  fHisto->GetYaxis()->SetRangeUser(fYLo, fYHi);
  fPad->Modified();

  // Also set marker size to default one pixel.

  SetMarkerSize(1., false);

  if(update)
    fPad->Update();
}

//......................................................................
// Set the marker size (measured in pixels).

void evd::Ortho3DPad::SetMarkerSize(double size, bool update)
{
  // Update marker size.

  if(fMSize != size/4.) {

    // Update marker size attribute.

    fMSize = size/4.;    // Scale to pixels.

    // Update widget.

    if(fMSizeEntry)
      fMSizeEntry->SetNumber(size);

    // Loop over graphic objects that are currently drawn on 
    // pad, and update any that are polymarkers.

    TIter next(fPad->GetListOfPrimitives());
    while(TObject* obj = next()) {
      if(obj->InheritsFrom(TPolyMarker::Class())) {
	TPolyMarker* pm = (TPolyMarker*)obj;
	pm->SetMarkerSize(fMSize);
      }
    }

    fPad->Modified();
    if(update)
      fPad->Update();
  }
}

//......................................................................
// Save a reference to the marker size number entry widget.
// This class gets the marker size from the widget if it is changed
// via the gui.  It also sets the number in the widget if it is changed
// not from the gui.

void evd::Ortho3DPad::SetMSizeEntry(TGNumberEntry* p)
{
  fMSizeEntry = p;
  if(fMSizeEntry)
    fMSizeEntry->SetNumber(4.*fMSize);
}

//......................................................................
// Slot method for marker size number entry widget.
// This method is called when the user changes the marker size via the gui.

void evd::Ortho3DPad::SetMSize()
{

  // Get marker size from number entry widget.

  if (!fMSizeEntry)
    throw cet::exception("Ortho3DPad") << __func__ << ": no MSize entry";
  double val = fMSizeEntry->GetNumber();

  // Scale the marker size such that the displayed marker size
  // is measured in pixels.

  SetMarkerSize(val, true);
}

//......................................................................
// Static mouse event handler.
// This method is called by the gui for mouse events in the graphics pad.

void evd::Ortho3DPad::MouseEvent(evd::Ortho3DPad* p)
{
  TObject *select = gPad->GetSelected();
  if(!select)
    return;
  if(!select->InheritsFrom("TBox"))
    return;
  ((TBox*)select)->SetBit(TBox::kCannotMove);
  p->MouseEvent();
}

//......................................................................
// Mouse event handler (called from static handler).

void evd::Ortho3DPad::MouseEvent()
{
  // Get event type and location.

  int event = gPad->GetEvent();
  int px = gPad->GetEventX();        // pixels.
  int py = gPad->GetEventY();        // pixels.
  double x = gPad->AbsPixeltoX(px);  // User coordinates.
  double y = gPad->AbsPixeltoY(py);  // User coordinates.

  // Handle different events.

  switch(event) {
  case kMouseEnter:

    // Main purpose of this case is to set the cursor shape.

    gPad->SetCursor(kCross);
    fCurrentPx = px;
    fCurrentPy = py;
    break;

  case kMouseMotion:

    // Not really needed...

    gPad->SetCursor(kCross);
    fCurrentPx = px;
    fCurrentPy = py;
    break;

  case kMouseLeave:

    // Undraw box.

    if(fBoxDrawn) {
      double pxlo = std::min(fPressPx, fCurrentPx);
      double pxhi = std::max(fPressPx, fCurrentPx);
      double pylo = std::min(fPressPy, fCurrentPy);
      double pyhi = std::max(fPressPy, fCurrentPy);
      gVirtualX->DrawBox(pxlo, pylo, pxhi, pyhi, TVirtualX::kHollow);
      fBoxDrawn = false;
    }

    // Set everything to default.

    fPress = false;
    fPressPx = 0;
    fPressPy = 0;
    fCurrentPx = 0;
    fCurrentPy = 0;
    fPressX = 0.;
    fPressY = 0.;
    fReleaseX = 0.;
    fReleaseY = 0.;
    break;

  case kButton1Motion:

    // Undraw old selection box.

    if(fBoxDrawn) {
      double pxlo = std::min(fPressPx, fCurrentPx);
      double pxhi = std::max(fPressPx, fCurrentPx);
      double pylo = std::min(fPressPy, fCurrentPy);
      double pyhi = std::max(fPressPy, fCurrentPy);
      gVirtualX->DrawBox(pxlo, pylo, pxhi, pyhi, TVirtualX::kHollow);
      fBoxDrawn = false;
    }

    // Update cursor location.

    gPad->SetCursor(kCross);
    fCurrentPx = px;
    fCurrentPy = py;

    // Draw new selection box.

    {
      double pxlo = std::min(fPressPx, fCurrentPx);
      double pxhi = std::max(fPressPx, fCurrentPx);
      double pylo = std::min(fPressPy, fCurrentPy);
      double pyhi = std::max(fPressPy, fCurrentPy);
      gVirtualX->DrawBox(pxlo, pylo, pxhi, pyhi, TVirtualX::kHollow);
      fBoxDrawn = true;
    }
    break;

  case kButton1Down:

    // Record the location of the button press event, which will be
    // one corner of zoom region.

    gPad->SetCursor(kCross);
    fPress = true;
    fPressPx = px;
    fPressPy = py;
    fCurrentPx = px;
    fCurrentPy = py;
    fPressX = x;
    fPressY = y;
    fReleaseX = 0.;
    fReleaseY = 0.;
    break;

  case kButton1Up:

    // Get the location of button release event, then zoom.

    gPad->SetCursor(kCross);
    fPress = false;
    fCurrentPx = px;
    fCurrentPy = py;
    fReleaseX = x;
    fReleaseY = y;
    {
      double xlo = std::min(fPressX, fReleaseX);
      double xhi = std::max(fPressX, fReleaseX);
      double ylo = std::min(fPressY, fReleaseY);
      double yhi = std::max(fPressY, fReleaseY);
      SetZoom(xlo, ylo, xhi, yhi, true);
    }
    break;
  }
}

////////////////////////////////////////////////////////////////////////
