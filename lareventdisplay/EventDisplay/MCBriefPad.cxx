///
/// \file    MCBriefPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
///
#include "TPad.h"

#include "art/Framework/Principal/fwd.h"

#include "lareventdisplay/EventDisplay/MCBriefPad.h"
#include "nuevdb/EventDisplayBase/View2D.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"

namespace evd{

  //......................................................................

  MCBriefPad::MCBriefPad(const char* nm, const char* ti,
			 double x1, double y1,
			 double x2, double y2,
			 const char* /*opt*/) :
    DrawingPad(nm, ti, x1, y1, x2, y2)
  {
    this->Pad()->cd();

    fView = new evdb::View2D();
  }

  //......................................................................

  MCBriefPad::~MCBriefPad()
  {
    if (fView) { delete fView; fView = 0; }
  }

  //......................................................................

  void MCBriefPad::Draw()
  {
    fView->Clear();
    this->Pad()->Clear();

    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      this->SimulationDraw()->MCTruthShortText(*evt, fView);
      this->SimulationDraw()->MCTruthLongText (*evt, fView);
    }
    fPad->cd();
    fView->Draw();
  }
}//namespace
//////////////////////////////////////////////////////////////////////////
