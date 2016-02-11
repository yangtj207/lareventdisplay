//
/// \file    Ortho3DView.cxx
/// \brief   Orthographic view display window
/// \author  greenlee@fnal.gov
///
#include <string>
#include <cmath>
#include "TCanvas.h"
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TVirtualX.h"
#include "TRootEmbeddedCanvas.h"
#include "lareventdisplay/EventDisplay/Ortho3DView.h"
#include "lareventdisplay/EventDisplay/Ortho3DPad.h"

#include "cetlib/exception.h"
#include "art/Framework/Principal/Event.h"

//......................................................................
// Consgtructor.

evd::Ortho3DView::Ortho3DView(TGMainFrame* mf) :
  evdb::Canvas(mf)
{
  // Remove everything that's in the main frame by default.

  mf->RemoveFrame(fEmbCanvas);
  mf->RemoveFrame(fFrame);

  // Make meta frame to hold graphics pad and our widgets.

  fMetaFrame = new TGHorizontalFrame(mf);

  // Refill main frame, replacing root canvas with meta frame.

  mf->AddFrame(fMetaFrame, fLayout);
  mf->AddFrame(fFrame);

  // Add two frames inside the meta frame.
  // Add a vertical frame on the left to hold our widgets.
  // Add root canvas on the right.

  fWidgetFrame = new TGVerticalFrame(fMetaFrame);
  fEmbCanvas->ReparentWindow(fMetaFrame);
  fMetaFrame->AddFrame(fWidgetFrame, new TGLayoutHints(kLHintsTop | 
						       kLHintsLeft | 
						       kLHintsExpandY));
  fMetaFrame->AddFrame(fEmbCanvas, fLayout);

  // Make vertically stacked subpads and widget subframes.
    
  int npad = 2;
  for(int ipad = 0; ipad < npad; ++ipad) {
    evdb::Canvas::fCanvas->cd();
    OrthoProj_t proj = evd::kNoProj;
    std::string projname;
    switch (ipad) {
      case 0:
        proj = kXZ;
        projname = "XZ";
        break;
      case 1:
        proj = kYZ;
        projname = "YZ";
        break;
      default:
        throw cet::exception("Ortho3DView")
          << __func__ << ": unknown projection pad " << ipad << "\n";
    } // switch
    
    std::string padname = std::string("Ortho3DPad") + projname;
    std::string padtitle = projname + std::string(" View");
    double ylo = double(npad - ipad - 1) / double(npad);
    double yhi = double(npad - ipad) / double(npad);
    Ortho3DPad* pad = new Ortho3DPad(padname.c_str(), padtitle.c_str(),
				     proj, 0.0, ylo, 1.0, yhi);
    fOrtho3DPads.push_back(pad);

    // Add subframe for this pad's widgets.

    TGCompositeFrame* wframe = new TGVerticalFrame(fWidgetFrame);
    fWidgetSubFrames.push_back(wframe);
    fWidgetFrame->AddFrame(wframe, new TGLayoutHints(kLHintsTop | 
						     kLHintsLeft | 
						     kLHintsExpandY));

    // Add widgets.

    // Label.

    TGLabel* label = new TGLabel(wframe, padtitle.c_str());
    wframe->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft,
					      5, 5, 5, 1));

    // Unzoom button.

    TGTextButton* unzoom = new TGTextButton(wframe, "&Unzoom");
    wframe->AddFrame(unzoom, new TGLayoutHints(kLHintsTop | kLHintsLeft,
					       5, 5, 5, 1));
    unzoom->Connect("Clicked()", "evd::Ortho3DPad", pad,
		    "UnZoom(=true)");

    // Marker size entry.

    TGCompositeFrame* msize_frame = new TGHorizontalFrame(wframe);
    wframe->AddFrame(msize_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft,
						    5, 5, 5, 1));
    int val = pad->GetMarkerSize();
    TGNumberEntry* msize_entry = new TGNumberEntry(msize_frame, val, 3, -1,
						   TGNumberFormat::kNESInteger,
						   TGNumberFormat::kNEANonNegative,
						   TGNumberFormat::kNELLimitMin, 1.);
    msize_frame->AddFrame(msize_entry);
    pad->SetMSizeEntry(msize_entry);

    TGLabel* msize_label = new TGLabel(msize_frame, "Marker Size");
    msize_frame->AddFrame(msize_label,
			  new TGLayoutHints(kLHintsTop | kLHintsLeft,
					    5, 0, 0, 1));
    msize_entry->Connect("ValueSet(Long_t)", "evd::Ortho3DPad", pad,
			 "SetMSize()");
  }

  // Draw everything and update canvas.

  Draw();
  evdb::Canvas::fCanvas->Update();
}
  
//......................................................................
// Destructor.
evd::Ortho3DView::~Ortho3DView() 
{
}

//......................................................................
// Draw object in graphics pads.
void evd::Ortho3DView::Draw(const char* /*opt*/) 
{
  for(std::vector<Ortho3DPad*>::const_iterator i = fOrtho3DPads.begin();
      i != fOrtho3DPads.end(); ++i) {
    Ortho3DPad* pad = *i;
    pad->Draw();
  }
}

////////////////////////////////////////////////////////////////////////
