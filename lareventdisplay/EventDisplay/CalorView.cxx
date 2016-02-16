//
/// \file    CalorView.cxx
/// \brief   Calorimetric view display window
/// \author  msoderbe@syr.edu
///
#include <iostream>
#include <sstream>
#include <cmath>

#include "TCanvas.h"
#include "TVirtualX.h"
#include "TRootEmbeddedCanvas.h"
#include "lareventdisplay/EventDisplay/CalorView.h"
#include "lareventdisplay/EventDisplay/CalorPad.h"
#include "lareventdisplay/EventDisplay/AnalysisDrawingOptions.h"

#include "art/Framework/Principal/Event.h"

//......................................................................
// Constructor.

evd::CalorView::CalorView(TGMainFrame* mf) : evdb::Canvas(mf)
{

  art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;
  
  evdb::Canvas::fCanvas->cd();
  if (anaOpt->fDrawShowerCalor){
    fDeDxPad = new CalorPad("fDeDxPad","DeDx Pad",0.0,0.5,1.0,1.0,2);
  }
  else{
    fDeDxPad = new CalorPad("fDeDxPad","DeDx Pad",0.0,0.5,1.0,1.0,1);
  }
  evdb::Canvas::fCanvas->cd();
  fKEPad = new CalorPad("fKEPad","Kinetic Energy Pad",0.0,0.0,1.0,0.5,0);
  
  this->Connect("CloseWindow()","evd::CalorView",this,"CloseWindow()");

  evdb::Canvas::fCanvas->Update(); 
}
  
//......................................................................
// Destructor.
evd::CalorView::~CalorView() 
{
  //if(fDeDxPad){ delete fDeDxPad; fDeDxPad = 0;}
  //if(fKEPad){ delete fKEPad; fKEPad = 0;}
}

//......................................................................
void evd::CalorView::CloseWindow()
{
  delete this;
}

//......................................................................
// Draw object in graphics pads.
void evd::CalorView::Draw(const char* /*opt*/) 
{

  //evdb::Canvas::fCanvas->ls();
  fDeDxPad->Pad()->cd();
  fDeDxPad->Draw();
  
  fKEPad->Pad()->cd();
  fKEPad->Draw();
  
  evdb::Canvas::fCanvas->Update();
}

////////////////////////////////////////////////////////////////////////
