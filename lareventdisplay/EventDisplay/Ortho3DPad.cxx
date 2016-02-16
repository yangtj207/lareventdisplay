///
/// \file    Ortho3DPad.cxx
/// \brief   Drawing pad showing an orthographic rendering of 3D objects
///          in the detector
/// \author  greenlee@fnal.gov
///

#include "TPad.h"
#include "TFrame.h"
#include "TPadPainter.h"
#include "TH1F.h"
#include "TBox.h"
#include "TPolyMarker.h"
#include "TGNumberEntry.h"
#include "TLatex.h"

#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lareventdisplay/EventDisplay/Ortho3DPad.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/EventHolder.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"

/// Define static data members.

evd::Ortho3DPad* evd::Ortho3DPad::fMousePad = 0;

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
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;

    switch (proj) {
    case evd::kXY:
      TPCBox.push_back(TBox(world[0]-geo->DetHalfWidth(i),
			    world[1]-geo->DetHalfHeight(i),
			    world[0]+geo->DetHalfWidth(i),
			    world[1]+geo->DetHalfHeight(i)));
      break;
    case evd::kXZ:
      TPCBox.push_back(TBox(world[2]-geo->DetLength(i)/2.,
			    world[0]-geo->DetHalfWidth(i),
			    world[2]+geo->DetLength(i)/2.,
			    world[0]+geo->DetHalfWidth(i)));
      break;
    case evd::kYZ:
      TPCBox.push_back(TBox(world[2]-geo->DetLength(i)/2.,
			    world[1]-geo->DetHalfHeight(i),
			    world[2]+geo->DetLength(i)/2.,
			    world[1]+geo->DetHalfHeight(i)));
      break;
    default:
      throw cet::exception("Ortho3DPad")
        << __func__ << ": unwknow projection " << ((int) proj) << "\n";
  } // switch
    TPCBox.back().SetFillStyle(0);
    TPCBox.back().SetLineStyle(2);
    TPCBox.back().SetLineWidth(2);
    TPCBox.back().SetLineColor(16);
  }
  
  switch (proj) {
    case evd::kXY:
      fXLo = minx;
      fXHi = maxx;
      fYLo = miny;
      fYHi = maxy;
      break;
    case evd::kXZ:
      fXLo = minz;
      fXHi = maxz;
      fYLo = minx;
      fYHi = maxx;
      break;
    case evd::kYZ:
      fXLo = minz;
      fXHi = maxz;
      fYLo = miny;
      fYHi = maxy;
      break;
    default:
      throw cet::exception("Ortho3DPad")
        << __func__ << ": unwknow projection " << ((int) proj) << "\n";
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
        << __func__ << ": unexpected flow (projection: " << ((int) proj) << ")\n";
  } // switch

  fHisto->GetXaxis()->CenterTitle();
  fHisto->GetYaxis()->SetLabelSize(0.04);
  fHisto->GetYaxis()->SetTitleSize(0.04);
  fHisto->GetYaxis()->CenterTitle();
  fHisto->SetFillColor(18);
  fHisto->Draw("AB");

  fView = new evdb::View2D();
    
  // Set pad fill color
  Pad()->SetFillColor(18);
  Pad()->SetFrameFillColor(18);
  Pad()->GetPainter()->SetFillColor(18);
  Pad()->Modified();
  Pad()->Update();

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
      RecoBaseDraw()->PFParticleOrtho(*evt, fProj, fMSize, fView);
      RecoBaseDraw()->ProngOrtho(*evt, fProj, fMSize, fView);
      RecoBaseDraw()->SeedOrtho(*evt, fProj, fView);
      RecoBaseDraw()->OpFlashOrtho(*evt, fProj, fView);
      RecoBaseDraw()->VertexOrtho(*evt, fProj, fView);
    }
  // Draw objects on pad.

  fPad->cd();
  fPad->GetPainter()->SetFillColor(18);
  fHisto->Draw("X-");
  fView->Draw();
  TLatex latex;
  latex.SetTextColor(16);
  latex.SetTextSize(0.05);
  for (size_t i = 0; i<TPCBox.size(); ++i){
    TPCBox[i].Draw();
    double x1 = TPCBox[i].GetX2() - 0.02*(fXHi-fXLo);
    double y1 = TPCBox[i].GetY2() - 0.05*(fYHi-fYLo);
    for (size_t j = 0; j<i; ++j){
      if (std::abs(x1-(TPCBox[j].GetX2() - 0.02*(fXHi-fXLo)))<1e-6&&
	  std::abs(y1-(TPCBox[j].GetY2() - 0.05*(fYHi-fYLo)))<1e-6){
	y1 -= 0.05*(fYHi-fYLo);
      }
    }
    latex.DrawLatex(x1,y1,Form("%lu",i));
  }
  fPad->Modified();
  fPad->Update();
  fBoxDrawn = false;
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
  if(update) {
    fPad->Update();
    fBoxDrawn = false;
  }
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

  if(update) {
    fPad->Update();
    fBoxDrawn = false;
  }
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
    if(update) {
      fPad->Update();
      fBoxDrawn = false;
    }
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
    throw cet::exception("Ortho3DPad") << __func__ << ": no MSize entry\n";
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

  // This line is a workaround for a root bug that sends mouse events
  // to the wrong pad.

  if(fMousePad != 0)
    p = fMousePad;

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
    fMousePad = this;
    break;

  case kMouseMotion:

    // Not really needed...

    gPad->SetCursor(kCross);
    fCurrentPx = px;
    fCurrentPy = py;
    break;

  case kMouseLeave:

    // Undraw box.

    //if(fBoxDrawn) {
    //  double pxlo = std::min(fPressPx, fCurrentPx);
    //  double pxhi = std::max(fPressPx, fCurrentPx);
    //  double pylo = std::min(fPressPy, fCurrentPy);
    //  double pyhi = std::max(fPressPy, fCurrentPy);
    //  gVirtualX->DrawBox(pxlo, pylo, pxhi, pyhi, TVirtualX::kHollow);
    //  fBoxDrawn = false;
    //}

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
    fMousePad = 0;
    break;

  case kButton1Motion:

    // Undraw old selection box.

    //if(fBoxDrawn) {
    //  double pxlo = std::min(fPressPx, fCurrentPx);
    //  double pxhi = std::max(fPressPx, fCurrentPx);
    //  double pylo = std::min(fPressPy, fCurrentPy);
    //  double pyhi = std::max(fPressPy, fCurrentPy);
    //  gVirtualX->DrawBox(pxlo, pylo, pxhi, pyhi, TVirtualX::kHollow);
    //  fBoxDrawn = false;
    //}

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
    gVirtualX->SetLineColor(-1);
    gVirtualX->SetLineStyle(0);
    gVirtualX->SetLineWidth(1);

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
    fMousePad = this;
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
    fMousePad = 0;
    break;
  }
}

////////////////////////////////////////////////////////////////////////
