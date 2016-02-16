///
/// \file    Display3DPad.cxx
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
/// \version $Id: Display3DPad.cxx,v 1.5 2011/04/12 22:07:16 bckhouse Exp $
///
#include <iostream>
#include "TPad.h"
#include "TView3D.h"
#include "TGLViewer.h"
#include "TPolyLine3D.h"

#include "lareventdisplay/EventDisplay/Display3DPad.h"
#include "EventDisplayBase/View3D.h"
#include "EventDisplayBase/EventHolder.h"
#include "larcore/Geometry/Geometry.h"
#include "lareventdisplay/EventDisplay/GeometryDrawer.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lardata/RecoBase/Seed.h"
#include "lardata/RecoObjects/BezierTrack.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd{

  ///
  /// Create a pad to show a 3D rendering of the detector and events
  /// @param nm : Name of the pad
  /// @param ti : Title of the pad
  /// @param x1 : Location of left  edge of pad (0-1)
  /// @param x2 : Location of right edge of pad (0-1)
  /// @param y1 : Location of bottom edge of pad (0-1)
  /// @param y2 : Location of top    edge of pad (0-1)
  /// @param opt: Options. Currently just a place holder
  ///
  Display3DPad::Display3DPad(const char* nm, const char* ti,
			     double x1, double y1,
			     double x2, double y2,
			     const char* /*opt*/) :
    DrawingPad(nm, ti, x1, y1, x2, y2)
  {
    this->Pad()->SetFillColor(kBlack);
    this->Pad()->Draw();
    this->Pad()->cd();
    fView = new evdb::View3D();
  }

  //......................................................................

  Display3DPad::~Display3DPad() 
  {
    if (fView) { delete fView; fView = 0; }
  }

  //......................................................................

  void Display3DPad::Draw() 
  {
    fView->Clear();

    art::ServiceHandle<geo::Geometry> geo;

    // grab the event from the singleton
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();

    if(evt){
      this->GeometryDraw()->  DetOutline3D(fView);
      this->SimulationDraw()->MCTruth3D   (*evt, fView);
      this->RecoBaseDraw()->  PFParticle3D(*evt, fView);
      this->RecoBaseDraw()->  SpacePoint3D(*evt, fView);
      this->RecoBaseDraw()->  Prong3D     (*evt, fView);
      this->RecoBaseDraw()->  Seed3D      (*evt, fView);
      this->RecoBaseDraw()->  BezierTrack3D (*evt, fView);
      this->RecoBaseDraw()->  Vertex3D    (*evt, fView);
      this->RecoBaseDraw()->  Event3D     (*evt, fView);
    
      art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutoptions;
      if(evdlayoutoptions->fMakeSeeds)
	UpdateSeedCurve();
    }

    this->Pad()->Clear();
    this->Pad()->cd();
    if (fPad->GetView()==0) {
      int irep;
      double rmin[]={-2.1*geo->DetHalfWidth(),-2.1*geo->DetHalfHeight(),-0.5*geo->DetLength()};
      double rmax[]={ 2.1*geo->DetHalfWidth(), 2.1*geo->DetHalfHeight(), 2.1*geo->DetLength()};
      TView3D* v = new TView3D(1,rmin,rmax);
      v->SetPerspective();
      v->SetView(0.0,260.0,270.0,irep);
      fPad->SetView(v); // ROOT takes ownership of object *v
    }
    fView->Draw();
    fPad->Update();
  }


  ////////////////////////////////////////////////////////////////////////


  void Display3DPad::UpdateSeedCurve()
  {
    std::vector<recob::Seed> SeedVec = this->HitSelectorGet()->SeedVector();
    std::cout<<"Seeds available to Display3DPad : " << SeedVec.size()<<std::endl;
    trkf::BezierTrack BTrack(SeedVec);
  
    int N=100;
    TPolyLine3D& pl = fView->AddPolyLine3D(N,kOrange+5,2,0);
    fView->Draw();
  

    double  xyzpt[3] ;
  
    for(int i=0; i!=N; i++)
      {
	BTrack.GetTrackPoint(float(i)/N,xyzpt );
	double x = xyzpt[0];
	double y = xyzpt[1];
	double z = xyzpt[2];

	pl.SetPoint(i,x,y,z);
      
      }
  
    fView->Draw();
    //  UpdatePad();
  }

}//namespace
