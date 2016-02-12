//
/// \file    TWQMultiTPCProjectionView.cxx
/// \brief   The "main" event display view that most people will want to use
/// \author  brebel@fnal.gov
/// \version $Id: XZYZProjectionsView.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "TCanvas.h"
#include "TFrame.h"
#include "TGFrame.h"  // For TGMainFrame, TGHorizontalFrame
#include "TGLayout.h" // For TGLayoutHints
#include "TGDimension.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TMath.h"
#include "TString.h"
#include "TRootEmbeddedCanvas.h"
#include "TLine.h"
#include "Buttons.h"
#include "TGTextView.h"
#include "TROOT.h"
#include "sstream"

#include "EventDisplayBase/View2D.h"
#include "lareventdisplay/EventDisplay/TWQMultiTPCProjection.h"
#include "lareventdisplay/EventDisplay/HeaderPad.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/TWireProjPad.h"
#include "lareventdisplay/EventDisplay/TQPad.h"
#include "lareventdisplay/EventDisplay/MCBriefPad.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lareventdisplay/EventDisplay/Style.h"

#include "lardata/RecoBase/Seed.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//#include "EventDisplay/InfoTransfer.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace evd{

  static unsigned int kPlane;
  static unsigned int kWire;
  static double       kDistance;
  static int          curr_zooming_plane;
  static const char*  zoom_opt=0;
  
  static int          shift_lock;
  

  //......................................................................
  TWQMultiTPCProjectionView::TWQMultiTPCProjectionView(TGMainFrame* mf) 
    : evdb::Canvas(mf)
  {  
    
    art::ServiceHandle<geo::Geometry> geo;

    // first make pads for things that don't depend on the number of 
    // planes in the detector
    // bottom left corner is (0.,0.), top right is  (1., 1.)

    evdb::Canvas::fCanvas->cd();  
    fHeaderPad = new HeaderPad("fHeaderPadMultiTPC","Header",0.0,0.0,0.15,0.13,"");  
    fHeaderPad->Draw();  

    evdb::Canvas::fCanvas->cd();  
    fMC = new MCBriefPad("fMCPadMultiTPC","MC Info.",0.15,0.13,1.0,0.17,"");  
    fMC->Draw();  

    evdb::Canvas::fCanvas->cd();  
    fWireQ = new TQPad("fWireQPadMultiTPC", "ADCvsTime",0.15,0.0,1.0,0.13,"TQ", 0, 0);  
    fWireQ->Pad()->SetBit(TPad::kCannotMove,true);
    fWireQ->Draw();  


    // add new "meta frame" to hold the GUI Canvas and a side frame (vframe)
    fMetaFrame  = new TGCompositeFrame(mf, 60, 60, kHorizontalFrame);
    fMetaFrame->SetBit(TPad::kCannotMove,true);

    //new frame organizing the buttons on the left of the canvas.
    fVFrame  = new TGCompositeFrame(fMetaFrame, 60, 60, kVerticalFrame); 
    // Define a layout for placing the canvas within the frame.
    fLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
				kLHintsExpandY, 5, 5, 5, 5);

    mf->RemoveFrame((TGFrame *)fEmbCanvas);
    mf->RemoveFrame(fFrame);

    fEmbCanvas->ReparentWindow( fMetaFrame, fXsize, fYsize);
	

    fMetaFrame->AddFrame(fVFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY));
    fMetaFrame->AddFrame(fEmbCanvas, fLayout);


    mf->AddFrame(fMetaFrame,fLayout);
    mf->AddFrame(fFrame);

    // plane number entry  
    fPlaneEntry = new TGNumberEntry(fFrame,
				    0,2,-1,
				    TGNumberFormat::kNESInteger, 
				    TGNumberFormat::kNEAAnyNumber, 
				    TGNumberFormat::kNELLimitMinMax, 
				    0, geo->Nplanes()-1 );

    kPlane = 0;
    kWire = TMath::Nint(0.5*geo->Nwires(0));
    kDistance=1.5;
    fWireQ->SetPlaneWire(kPlane, kWire);

    // Initial value
    fPlaneEntry->SetNumber( kPlane );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fPlaneEntry->Connect("ValueSet(Long_t)", "evd::TWQMultiTPCProjectionView", this, "SetPlane()");
    fPlaneEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQMultiTPCProjectionView", this, "SetPlane()");
    // Text label for this numeric field.
    fPlaneLabel= new TGLabel(fFrame,"Plane");

    // wire number entry 
    fWireEntry = new TGNumberEntry(fFrame,0,6,-1,
				   TGNumberFormat::kNESInteger, 
				   TGNumberFormat::kNEAAnyNumber, 
				   TGNumberFormat::kNELLimitMinMax, 
				   0, geo->Nwires(0)-1 );
    // Initial value
    fWireEntry->SetNumber( kWire );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fWireEntry->Connect("ValueSet(Long_t)", "evd::TWQMultiTPCProjectionView", this, "SetWire()");
    fWireEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQMultiTPCProjectionView", this, "SetWire()");

    // Text label for this numeric field.
    fWireLabel= new TGLabel(fFrame,"Wire");

    // adc threshold number entry 
    fThresEntry = new TGNumberEntry(fFrame,0,6,-1,
				    TGNumberFormat::kNESInteger, 
				    TGNumberFormat::kNEAAnyNumber, 
				    TGNumberFormat::kNELLimitMinMax, 
				    0 , geo->Nwires(0)-1 );
    // Initial value
    art::ServiceHandle<evd::ColorDrawingOptions>      cst;
    art::ServiceHandle<evd::SimulationDrawingOptions> sdo;
    art::ServiceHandle<evd::RawDrawingOptions>        rawopt;
    art::ServiceHandle<evd::EvdLayoutOptions>         evdlayoutopt;


    fThresEntry->SetNumber( rawopt->fMinSignal );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fThresEntry->Connect("ValueSet(Long_t)", "evd::TWQMultiTPCProjectionView", this, "SetThreshold()");
    fThresEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQMultiTPCProjectionView", this, "SetThreshold()");

    // Text label for this numeric field.
    fThresLabel= new TGLabel(fFrame,"ADC Threshold");

    // check button to toggle color vs grey
    fGreyScale = new TGCheckButton(fFrame,"Grayscale",1);
    fGreyScale->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetGreyscale()");
    if(cst->fColorOrGray == 1) fGreyScale->SetState(kButtonDown);

    // check button to toggle MC information
    if(evdlayoutopt->fEnableMCTruthCheckBox){
      fMCOn = new TGCheckButton(fFrame,"MC Truth",5);
      fMCOn->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetMCInfo()");
      if(sdo->fShowMCTruthText == 1) fMCOn->SetState(kButtonDown);
    }

    // radio buttons to toggle drawing raw vs calibrated information
    fRawCalibDraw = new TGRadioButton(fFrame,"Both",          2);
    fCalibDraw    = new TGRadioButton(fFrame,"Reconstructed", 3);
    fRawDraw      = new TGRadioButton(fFrame,"Raw",           4);
    fRawDraw     ->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetRawCalib()");
    fCalibDraw   ->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetRawCalib()");
    fRawCalibDraw->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetRawCalib()");
    if(rawopt->fDrawRawDataOrCalibWires == 0)      fRawDraw->SetState(kButtonDown);
    else if(rawopt->fDrawRawDataOrCalibWires == 1) fCalibDraw->SetState(kButtonDown);
    else if(rawopt->fDrawRawDataOrCalibWires == 2) fRawCalibDraw->SetState(kButtonDown);

    // Put all these widgets into the frame.  The last
    // four numbers in each TGLayoutHint are padleft, padright,
    // padtop, padbottom.
   if(evdlayoutopt->fEnableMCTruthCheckBox){
     fFrame->AddFrame(fMCOn,          new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   }
   fFrame->AddFrame(fGreyScale,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fRawCalibDraw,  new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fCalibDraw,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fRawDraw,       new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fPlaneEntry,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
   fFrame->AddFrame(fPlaneLabel,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fWireEntry,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
   fFrame->AddFrame(fWireLabel,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
   fFrame->AddFrame(fThresEntry,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
   fFrame->AddFrame(fThresLabel,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );

    if(evdlayoutopt->fShowSideBar)
      SetUpSideBar();    
    else
      evdlayoutopt->fShowEndPointSection=0;  // zero it to avoid a misconfiguration in the fcl file. 

    //zero the ppoints queue.
    ppoints.clear();
    pline.clear();

    // geometry to figure out the number of TPCs
    unsigned int ntpc = geo->NTPC();

    // geometry to figure out the number of planes
    unsigned int nplanes = geo->Nplanes();

    // now determine the positions of all the time vs wire number 
    // and charge histograms for the planes
    for(unsigned int t = 0; t < ntpc; ++t){
      for(unsigned int i = 0; i < nplanes; ++i){
	double twx1 = 0. + t*0.97/(1.*ntpc);
	double twx2 = t*0.97/(1.*ntpc);
	double twx3 = 1.0;
	double twy1 = 0.17 +   (i)*(1.0-0.171)/(1.*nplanes);
	double twy2 = 0.17 + (i+1)*(1.0-0.171)/(1.*nplanes);
      
	TString padname = "fWireProjTPC";
	padname += t;
	padname += "Plane";
	padname += i;
      
	TString padtitle = "TPC";
	padtitle += t;
	padtitle += "Plane";
	padtitle += i;
      
	evdb::Canvas::fCanvas->cd();

	mf::LogVerbatim("MultiTPC") << "make new plane ";
	fPlanes.push_back(new TWireProjPad(padname, padtitle, twx1, twy1, twx2, twy2, i+t*nplanes));
	fPlanes.back()->Draw();
//	fPlanes.back()->Pad()->AddExec("mousedispatch",
//				       Form("evd::TWQMultiTPCProjectionView::MouseDispatch(%d, (void*)%d)", 
//					    i+t*nplanes, this)
//				       ); 
	fPlanes.back()->Pad()->AddExec("mousedispatch",
				       Form("evd::TWQMultiTPCProjectionView::MouseDispatch(%d, (void*)%lu)", 
					    i+t*nplanes, (unsigned long) this)
				       ); 
      
	mf::LogVerbatim("MultiTPC") << "size of planes vec is now " << fPlanes.size();

	if(t+1 == ntpc){
	  padname = "fQPadTPC";
	  padname += t;
	  padname += "Plane";
	  padname += i;
	  
	  padtitle = "QTPC";
	  padtitle += t;
	  padname += "Plane";
	  padname += i;
	  
	  evdb::Canvas::fCanvas->cd();
	  fPlaneQ.push_back(new TQPad(padname, padtitle, twx2, twy1, twx3, twy2, "Q", i, 0));
	  fPlaneQ[i]->Draw();
	}
      }// end loop to draw pads
    }

    evdb::Canvas::fCanvas->Update();

  }

  //......................................................................
  TWQMultiTPCProjectionView::~TWQMultiTPCProjectionView() 
  {  
    if (fHeaderPad) { delete fHeaderPad;  fHeaderPad  = 0; }  
    if (fMC)        { delete fMC;         fMC         = 0; }  
    if (fWireQ)     { delete fWireQ;      fWireQ      = 0; }  
    if (fPlaneEntry){ delete fPlaneEntry; fPlaneEntry = 0; }
    if (fWireEntry) { delete fWireEntry;  fWireEntry  = 0; }
    if (fPlaneLabel){ delete fPlaneLabel; fPlaneLabel = 0; }
    if (fWireLabel) { delete fWireLabel;  fWireLabel  = 0; }
    for(unsigned int i = 0; i < fPlanes.size(); ++i){
      if(fPlanes[i]){ delete fPlanes[i];  fPlanes[i]  = 0; }
      if(fPlaneQ[i]){ delete fPlaneQ[i];  fPlaneQ[i]  = 0; }
    }
    fPlanes.clear();
    fPlaneQ.clear();
  }

  //......................................................................
  void TWQMultiTPCProjectionView::DrawPads(const char* /*opt*/)
  {
    for(unsigned int i=0; i<fPlanes.size();++i){
      fPlanes[i]->Draw();
      fPlanes[i]->Pad()->Update();
      fPlanes[i]->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);
    }
    for(unsigned int j=0;j<fPlaneQ.size();++j){
      fPlaneQ[j]->Draw();
      fPlaneQ[j]->Pad()->Update();
      fPlaneQ[j]->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);
    }
  }
  //......................................................................
  void TWQMultiTPCProjectionView::Draw(const char* opt) 
  {  
    art::ServiceHandle<geo::Geometry> geo;

    fPrevZoomOpt.clear();
  
    evdb::Canvas::fCanvas->cd();    
    zoom_opt=0;
    fHeaderPad->Draw();    
    fMC       ->Draw(); 
    fWireQ->Draw();
  
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;

    if(evdlayoutopt->fPrintTotalCharge) PrintCharge();
  
    //clear queue of selected points
    ppoints.clear();
    pline.clear();
    // Reset current zooming plane - since it's not currently zooming.
    curr_zooming_plane=-1;
  
    //  double Charge=0, ConvCharge=0;
    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw(opt);
      fPlanes[i]->Pad()->Update();
      fPlanes[i]->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);
      fPlaneQ[i]->Draw();
      std::vector<double> ZoomParams = fPlanes[i]->GetCurrentZoom();
      fZoomOpt.wmin[i] = ZoomParams[0];
      fZoomOpt.wmax[i] = ZoomParams[1];
      fZoomOpt.tmin[i] = ZoomParams[2];
      fZoomOpt.tmax[i] = ZoomParams[3];
    }
  
    // Reset any text boxes which are enabled
    if(fXYZPosition)
      fXYZPosition->SetForegroundColor(kBlack);

    if(fAngleInfo)
      fAngleInfo->SetForegroundColor(kBlack);

    evdb::Canvas::fCanvas->Update();
  }

  // comment out this method as for now we don't want to change every
  // plane to have the same range in wire number because wire numbers
  // don't necessarily overlap from plane to plane, ie the same range
  // isn't appropriate for every plane
  //......................................................................
  //   void TWQMultiTPCProjectionView::RangeChanged() 
  //   {  
  //     static int ilolast = -1;  
  //     static int ihilast = -1;    
  // 
  //     int ilo; 
  //     int ihi;
  //     std::vector<int> lo;
  //     std::vector<int> hi;
  //     std::vector<bool> axischanged;
  //     for(unsigned int i = 0; i < fPlanes.size(); ++i){
  //       fPlanes[i]->GetWireRange(&ilo, &ihi);  
  //       lo.push_back(ilo);
  //       hi.push_back(ihi);
  //       axischanged.push_back((ilo != ilolast) || (ihi != ihilast));
  //     }
  //       
  //     TVirtualPad* ori = gPad;  
  // 
  //     // loop over the bools to see which axes need to change
  //     for(unsigned int i = 0; i < axischanged.size(); ++i){
  //       if (axischanged[i]) {    
  // 	fPlanes[i]->SetWireRange(ilo, ihi);    
  // 	fPlanes[i]->Pad()->cd();    
  // 	fPlanes[i]->Pad()->Modified();    
  // 	fPlanes[i]->Pad()->Update();    
  // 
  // 	ilolast = ilo;    
  // 	ihilast = ihi;  
  //       }  
  //     }
  // 
  //     evdb::Canvas::fCanvas->cd();  
  //     evdb::Canvas::fCanvas->Modified();  
  //     evdb::Canvas::fCanvas->Update();  
  //     ori->cd();
  //   }
  //......................................................................
  //......................................................................
  void 	TWQMultiTPCProjectionView::PrintCharge()
  {

    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;

    for(size_t iplane = 0; iplane < fPlanes.size(); ++iplane){ 
      if(geo->Plane(iplane).SignalType()==geo::kCollection){
	double ch=0,convch=0;
	if(rawopt->fDrawRawDataOrCalibWires == 0){
	  fPlanes[iplane]->RawDataDraw()->GetChargeSum(iplane,ch,convch);
	  mf::LogVerbatim("TWQMultiTPCProjectionView") << "Warning! Calculating for RawData! ";
	}
	else{  
	  fPlanes[iplane]->RecoBaseDraw()->GetChargeSum(iplane,ch,convch);  
	}    

	mf::LogVerbatim("TWQMultiTPCProjectionView") << "\ncharge collected at collection plane: " 
					     << iplane << " " << ch << " " << convch;
      }
    }


  }

  //-------------------------------------------------------------------
  //......................................................................
  void TWQMultiTPCProjectionView::MouseDispatch(int plane, void * wqpv)
  {
    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();
    evd::TWQMultiTPCProjectionView *wqpp = (evd::TWQMultiTPCProjectionView*)wqpv;
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;

    switch (event){
	
    case kButton1Shift:           	
      shift_lock=1;
      // 	TWQMultiTPCProjectionView::SelectHit() is undefined
      //if(evdlayoutopt->fMakeClusters==1){wqpp->SelectHit(plane);}
      //else {wqpp->SelectPoint(plane);}
      wqpp->SelectPoint(plane);
      break;
    case kButton1Up:    
      if(shift_lock==1) break;
      if(evdlayoutopt-> fChangeWire==1) wqpp->ChangeWire(plane);		
    case kButton1Down: shift_lock=0;
    case kButton1Motion:	
      wqpp->SetMouseZoomRegion(plane);
      break;
      //  default:		
    }
  }


  //......................................................................
  void TWQMultiTPCProjectionView::ChangeWire(int plane)
  {
    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();
    int px = gPad->GetEventX();
    if(event!=11) return;
    TObject *select = gPad->GetSelected();
    if(!select) return;
    if(!select->InheritsFrom("TBox")) return;

    //now find wire that was clicked on
    float xx = gPad->AbsPixeltoX(px);
    float x = gPad->PadtoX(xx);
	

    kPlane = plane;
    kWire  = (unsigned int)TMath::Nint(x);

    this->SetPlaneWire();

    return;

  }

  //......................................................................
  void TWQMultiTPCProjectionView::SelectPoint(int plane)
  {
    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();

    if(event!=7) return;

    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    if(evdlayoutopt->fShowEndPointSection!=1)
      return;
    //struct planepoint;
    int px = gPad->GetEventX();
    double w0 = gPad->AbsPixeltoX(px);
    double x = gPad->PadtoX(w0);

    int py = gPad->GetEventY();
    double t0 = gPad->AbsPixeltoY(py);
    double y = gPad->PadtoY(t0);

    util::PxPoint ppx(plane,x,y);
    curr_zooming_plane=-1;

    // check if not clicking on a plane that is already in the ppoints list:
    int repeat_plane=-1;
    for(size_t ii = 0; ii < this->ppoints.size(); ++ii)
      if(ppx.plane==this->ppoints[ii].plane){
	this->ppoints[ii]=ppx;
	//clear View and draw new Marker
	this->fPlanes[this->ppoints[ii].plane]->View()->Clear();
	if(evdlayoutopt->fShowEndPointMarkers)
	  this->fPlanes[this->ppoints[ii].plane]->View()->AddMarker(ppx.w, ppx.t, kRed, 29, 2.0);
	else
	  this->fPlanes[plane]->View()->AddMarker(0.0,0.0,2,1,0.1);
	this->fPlanes[this->ppoints[ii].plane]->View()->Draw();
	repeat_plane=this->ppoints[ii].plane;  
	break;
      }
	
    //if plane does not repeat and size of list is larger than 2 pop_front 
    // and delete its marker. Otherwise just push_back.
    if(repeat_plane==-1){
      if( this->ppoints.size()>=2){
	this->fPlanes[this->ppoints[0].plane]->Pad()->cd();
	this->fPlanes[this->ppoints[0].plane]->View()->Clear();
	this->fPlanes[this->ppoints[0].plane]->View()->Draw();
	this->ppoints.pop_front();
      }
      this->ppoints.push_back(ppx);
      this->fPlanes[plane]->Pad()->cd();
      this->fPlanes[plane]->View()->Clear();
      if(evdlayoutopt->fShowEndPointMarkers)
	this->fPlanes[plane]->View()->AddMarker(ppx.w, ppx.t, kRed, 29, 2.0);
      else
	this->fPlanes[plane]->View()->AddMarker(0.0,0.0,2,1,0.1);
      this->fPlanes[plane]->View()->Draw();
    }

	
    return;

  }

  //......................................................................
  void TWQMultiTPCProjectionView::ClearEndPoints()
  {
    for (size_t x = 0; x < fPlanes.size(); ++x){
      fPlanes[x]->Pad()->cd();
      fPlanes[x]->View()->Clear();
      fPlanes[x]->View()->AddMarker(0.0,0.0,2,1,0.1);
      fPlanes[x]->Pad()->Update(); 
      fPlanes[x]->View()->Draw();
    }
    ppoints.clear();
    gPad->Modified();
    gPad->Update();
    gPad->cd();

  }


  //......................................................................
  double TWQMultiTPCProjectionView::FindLineLength()
  {
    // if list is larger than or equal to two, can project to XYZ and extrapolate to third plane (if exists)

    ////for now leaving commented. At some point be useful to display this information somewhere
    // for(unsigned int ix=0;ix<ppoints.size();ix++)
    //   std::cout << "ppoints, planes,x,y :" << ix << " " << ppoints[ix].plane << " " << ppoints[ix].x << " " << ppoints[ix].y << std::endl;

    if(pline.size() >= 2){

      double xyz_vertex_fit[3];
      double second_time;
      double pos[3];
      const double origin[3] = {0.,0.,0.};
      double xx0 = 0., yy0 = 0., zz0 = 0.;
      double xx1 = 0., yy1 = 0., zz1 = 0.;
      double length;

      double y,z;

      art::ServiceHandle<geo::Geometry> geom;
      const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
      art::ServiceHandle<evd::RawDrawingOptions> rawOpt;
      double ftimetick = detp->SamplingRate()/1000.;
      double larv = detp->DriftVelocity(detp->Efield(), detp->Temperature());
		
      //find channels corresponding to found wires.
      int chan1 = geom->PlaneWireToChannel(pline[0].plane,pline[0].w0, rawOpt->fTPC, rawOpt->fCryostat);
      int chan2 = geom->PlaneWireToChannel(pline[1].plane,pline[1].w0, rawOpt->fTPC, rawOpt->fCryostat);

      bool wires_cross=false;
      bool time_good=false;
	
      if(fabs(pline[0].t0-pline[1].t0) < 200){
	wires_cross= geom->ChannelsIntersect(chan1,chan2,y,z);
	time_good=true;
      }
      else{
	TGText *tt=new TGText("too big");
	tt->InsLine(1,"time distance");  
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	// return; //not returning, because may need to delete marker from wplane
      }

      if(wires_cross){
	TGText *tt=new TGText("wires cross");
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	geom->Plane(pline[0].plane).LocalToWorld(origin, pos);
	xyz_vertex_fit[0]=(pline[0].t0-detp->TriggerOffset())*larv*ftimetick+pos[0];
	geom->Plane(pline[1].plane).LocalToWorld(origin, pos);
	second_time=(pline[1].t0-detp->TriggerOffset())*larv*ftimetick+pos[0];
	
	xx0=(xyz_vertex_fit[0]+second_time)/2;
	yy0=y;
	zz0=z;

	//////////// the xyz vertex is found. Can proceed to calulate distance from edge
      }
      else{
	if(time_good){    //otherwise the wires_cross are false by default
	  TGText *tt=new TGText("cross");
	  tt->InsLine(1,"wires do not");
	  fXYZPosition->SetText(tt);
	  fXYZPosition->Update();
	}
	// return; //not returning, because may need to delete marker from wplanereturn;
      }
      //find channels corresponding to found wires AT END OF LINE.
      chan1 = geom->PlaneWireToChannel(pline[0].plane,pline[0].w1, rawOpt->fTPC, rawOpt->fCryostat);
      chan2 = geom->PlaneWireToChannel(pline[1].plane,pline[1].w1, rawOpt->fTPC, rawOpt->fCryostat);

      wires_cross=false;
      time_good=false;
	
      if(fabs(pline[0].t1-pline[1].t1) < 200){
	wires_cross= geom->ChannelsIntersect(chan1,chan2,y,z);
	time_good=true;
      }
      else{
	TGText *tt=new TGText("too big");
	tt->InsLine(1,"time distance");  
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	// return; //not returning, because may need to delete marker from wplane
      }

      if(wires_cross){
	TGText *tt=new TGText("wires do cross");
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	geom->Plane(pline[0].plane).LocalToWorld(origin, pos);
	xyz_vertex_fit[0]=(pline[0].t1-detp->TriggerOffset())*larv*ftimetick+pos[0];
	geom->Plane(pline[1].plane).LocalToWorld(origin, pos);
	second_time=(pline[1].t1-detp->TriggerOffset())*larv*ftimetick+pos[0];

	xx1=(xyz_vertex_fit[0]+second_time)/2;
	yy1=y;
	zz1=z;

      }
      else{
	if(time_good){    //otherwise the wires_cross are false by default
	  TGText *tt=new TGText("cross");
	  tt->InsLine(1,"wires do not");
	  fXYZPosition->SetText(tt);
	  fXYZPosition->Update();
	}
	// return; //not returning, because may need to delete marker from wplanereturn;
      }
      //update pad?
      gPad->Modified();
      gPad->Update();
      gPad->cd();

      length=pow(xx0-xx1,2)+pow(yy0-yy1,2)+pow(zz0-zz1,2);
      length=pow(length,0.5);
      return length;
    } // end if( this->ppoints.size()>=2)

    else{
      TGText *tt=new TGText("selected points");
      tt->InsLine(1,"not enough");
      fXYZPosition->SetText(tt);
      fXYZPosition->Update();
    }



    return -99;
  }

  //......................................................................
  void TWQMultiTPCProjectionView::FindEndPoint()
  {
    // if list is larger than or equal to two, can project to XYZ and extrapolate to third plane (if exists)

    if(ppoints.size()>=2 ){

      double xyz_vertex_fit[3] = {0.};
      double second_time = 0.;
      double pos[3] = {0.};
      const double origin[3] = {0.,0.,0.};
      double y = 0.;
      double z = 0.;

      art::ServiceHandle<geo::Geometry> geom;
      const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
      art::ServiceHandle<evd::RawDrawingOptions> rawOpt;
      double ftimetick = detp->SamplingRate()/1000.;
      double larv = detp->DriftVelocity(detp->Efield(), detp->Temperature());
		
      //find channels corresponding to found wires.
      int chan1 = geom->PlaneWireToChannel(ppoints[0].plane,ppoints[0].w, rawOpt->fTPC, rawOpt->fCryostat);
      int chan2 = geom->PlaneWireToChannel(ppoints[1].plane,ppoints[1].w, rawOpt->fTPC, rawOpt->fCryostat);

      bool wires_cross=false;
      bool time_good=false;
	
      if(fabs(ppoints[0].t-ppoints[1].t) < 200){
	wires_cross= geom->ChannelsIntersect(chan1,chan2,y,z);
	time_good=true;
      }
      else{
	TGText *tt=new TGText("too big");
	tt->InsLine(1,"time distance");  
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	// return; //not returning, because may need to delete marker from wplane
      }

      if(wires_cross){
	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	geom->Plane(ppoints[0].plane).LocalToWorld(origin, pos);
	xyz_vertex_fit[0]=(ppoints[0].t-detp->TriggerOffset())*larv*ftimetick+pos[0];
	geom->Plane(ppoints[1].plane).LocalToWorld(origin, pos);
	second_time=(ppoints[1].t-detp->TriggerOffset())*larv*ftimetick+pos[0];
		
	TGText *tt=new TGText(Form("z:%4.1f",z));
	tt->InsLine(1,Form("x:%4.1f,",(xyz_vertex_fit[0]+second_time)/2)); 
	tt->InsLine(1,Form("y:%4.1f,",y));  
	fXYZPosition->SetText(tt);
	fXYZPosition->Update();
	//////////// the xyz vertex is found. Can proceed to calulate distance from edge
	
      }
      else{
	if(time_good){    //otherwise the wires_cross are false by default
	  TGText *tt=new TGText("cross");
	  tt->InsLine(1,"wires do not");
	  fXYZPosition->SetText(tt);
	  fXYZPosition->Update();
	}
	// return; //not returning, because may need to delete marker from wplanereturn;
      }
      // extrapolate third point only if there are enough planes
      if(fPlanes.size() > 2){
	
	unsigned int wplane = 0;
	unsigned int wirevertex = 0;
	art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
	
	for(size_t xx = 0; xx < fPlanes.size(); ++xx){
	  wplane = 0;
	  for(int yy = 0; yy < 2; ++yy)
	    if(ppoints[yy].plane == xx)
	      ++wplane;
	
	  if(!wplane){ 
	    wplane = xx;
	    break;
	  }
	}
	
	
	geom->Plane(wplane).LocalToWorld(origin, pos);
	pos[1]=xyz_vertex_fit[1];
	pos[2]=xyz_vertex_fit[2];
	
	wirevertex = geom->NearestWire(pos, wplane, rawOpt->fTPC, rawOpt->fCryostat);
	
	double drifttick=((xyz_vertex_fit[0])/detp->DriftVelocity(detp->Efield(),detp->Temperature()))*(1./ftimetick);
	double timestart=drifttick-(pos[0]/detp->DriftVelocity(detp->Efield(),detp->Temperature()))*(1./ftimetick)+detp->TriggerOffset();
	
	fPlanes[wplane]->Pad()->cd();
	fPlanes[wplane]->View()->Clear();
	if(wires_cross && evdlayoutopt->fShowEndPointMarkers)  //only Draw if it makes sense
	  fPlanes[wplane]->View()->AddMarker(wirevertex, timestart, kMagenta, 29, 2.0);
	else  //draw dummy marker to delete old one
	  fPlanes[wplane]->View()->AddMarker(0.0,0.0,2,1,0.1);
	fPlanes[wplane]->Pad()->Update(); 
	fPlanes[wplane]->View()->Draw();	
      }// end if(fPlanes.size()>2)
      //update pad?
      gPad->Modified();
      gPad->Update();
      gPad->cd();
    } // end if( this->ppoints.size()>=2)
    else{
      TGText *tt=new TGText("selected points");
      tt->InsLine(1,"not enough");
      fXYZPosition->SetText(tt);
      fXYZPosition->Update();
    }

    return;
  }

  //.......................................................................
  void TWQMultiTPCProjectionView::SetMouseZoomRegion(int plane)
  {
    //*-*-*-*-*-*-*-*-*-*-*Create a new arrow in this pad*-*-*-*-*-*-*-*-*-*-*-*-*
    //*-*                  ==============================
    //
    TObject *select = gPad->GetSelected();
    if(!select) return;
    if(!select->InheritsFrom("TBox")) return;

    static Float_t w0=-1, t0=-1, w1=-1, t1=-1;

    static Int_t pxold, pyold;
    static Int_t pw0, pt0;
    static Int_t linedrawn;
    //static int curr_plane;
    //TLine *line;

    static int wstart,wend;
    static float tstart,tend;

    int event = gPad->GetEvent();
    int px = gPad->GetEventX();
    int py = gPad->GetEventY();

    switch (event){

    case kButton1Down:{
      gVirtualX->SetLineColor(-1);
      w0 = gPad->AbsPixeltoX(px);
      t0 = gPad->AbsPixeltoY(py);
      pw0   = px; pt0   = py;
      pxold = px; pyold = py;
      linedrawn = 0;
      float x = gPad->PadtoX(w0);
      tstart = gPad->PadtoY(t0);

      wstart  = (unsigned int)TMath::Nint(x);
      curr_zooming_plane=plane;
      break;
    }
    case kButton1Motion:{ 
      int lx,hx,ly,hy;
      if (pw0 < pxold){
	lx=pw0; 
	hx=pxold; 
      }
      else{
	lx=pxold;
	hx=pw0;
      }

      if (pt0 < pyold){
	ly=pt0; 
	hy=pyold;
      }
      else{
	ly=pyold;
	hy=pt0;
      }

      if (linedrawn) gVirtualX->DrawBox(lx, ly, hx, hy,TVirtualX::kHollow);
      pxold = px;
      pyold = py;
      linedrawn = 1;

      if (pw0 < pxold){
	lx=pw0; 
	hx=pxold; 
      }
      else{
	lx=pxold;
	hx=pw0;
      }

      if (pt0 < pyold){
	ly=pt0; 
	hy=pyold;
      }
      else{
	ly=pyold;
	hy=pt0;
      }

      gVirtualX->DrawBox(lx, ly, hx, hy,TVirtualX::kHollow);
      break;
    }
    case kButton1Up:{
      if (px == pw0 && py == pt0) break;
      w1 = gPad->AbsPixeltoX(px);
      t1 = gPad->AbsPixeltoY(py);
      gPad->Modified(kTRUE);
		  
      //   line = new TLine(w0,t0,w1,t1);
      //   line->Draw();
		  
      float x = gPad->PadtoX(w1);
      tend = gPad->PadtoY(t1);
      wend  = (unsigned int)TMath::Nint(x);
		  
      gROOT->SetEditorMode();
		  
      //make sure the box is significantly big to avoid accidental zooms on nothing.
      double xx1,yy1,xx2,yy2;
		  
      gPad->GetRangeAxis(xx1, yy1, xx2, yy2);
		  
      if(wstart != 0 && tstart != 0 && 
	 ( fabs(wend-wstart ) > 0.01*(xx2-xx1) ) && 
	 ( fabs(tend-tstart ) > 0.01*(yy2-yy1)   &&  
	   curr_zooming_plane==plane ) ){
		    		    
	this->SetZoom(plane,wstart,wend,tstart,tend);
	wstart=-1;
	tstart=-1;
      }
      break;
    }
    }// end switch
  }


  //......................................................................
  // if flag is true then zoom. If flag is false then unzoom.
  void 	TWQMultiTPCProjectionView::ZoomInterest(bool flag)
  {
    mf::LogVerbatim("TWQMultiTPCProjectionView") <<"ZoomInterest called";
  
    if(flag==true) zoom_opt="1";
    else zoom_opt="0";
  
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt; 
 
    ZoomOptionsMultiTPC zo;
    //  mf::LogVerbatim("TWQMultiTPCProjectionView") <<"Zoom interest pushing back zoom options"<<std::endl;
    fPrevZoomOpt.push_back(fZoomOpt);

 
    for(size_t iplane = 0; iplane < fPlanes.size(); ++iplane){
      int minw,maxw,mint,maxt;
      if(flag){
	int test=0;
	if(rawopt->fDrawRawDataOrCalibWires == 0)
	  fPlanes[iplane]->RawDataDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
	else
	  fPlanes[iplane]->RecoBaseDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
	  
	if(test==-1)
	  continue;
      }
      else{
	minw = -0.005*(geo->Nwires(iplane)-1);
	maxw =  1.005*(geo->Nwires(iplane)-1);
	mint = -0.005*fPlanes[iplane]->RawDataDraw()->TotalClockTicks();
	maxt =  1.01*fPlanes[iplane]->RawDataDraw()->TotalClockTicks();
      }
      
      SetZoom(iplane,minw,maxw,mint,maxt,false);
      zo.wmin[iplane]=minw;
      zo.tmin[iplane]=mint;
      zo.wmax[iplane]=maxw;
      zo.tmax[iplane]=maxt;
      zo.OnlyPlaneChanged=-1;
    }
    fZoomOpt=zo;

  }

  //......................................................................
  void TWQMultiTPCProjectionView::SetUpSideBar()
  {  
    SetUpZoomButtons();
    SetUpPositionFind();
  }

  //......................................................................
  void TWQMultiTPCProjectionView::SetZoomInterest()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fAutoZoomInterest = fToggleAutoZoom->GetState();
  }

  //......................................................................   
  void TWQMultiTPCProjectionView::ToggleEndPointMarkers()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fShowEndPointMarkers= fToggleShowMarkers->GetState();
  }

  //......................................................................
  void TWQMultiTPCProjectionView::SetUpZoomButtons()
  {
    // enter zoom buttons
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  

    fZoomInterest=new TGTextButton(fVFrame,"&Zoom Interest",150);
    fZoomInterest->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "ZoomInterest()");



    fUnZoomInterest=new TGTextButton(fVFrame,"&UnZoom Interest",150);
    fUnZoomInterest->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "ZoomInterest(=false)");


    fZoomBack=new TGTextButton(fVFrame,"&Zoom Back",150);
    fZoomBack->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "ZoomBack()");


    fToggleAutoZoom=new TGCheckButton(fVFrame,"AutoZoom",0);;       ///< Toggle the autozoom setting 
    fToggleAutoZoom->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "SetZoomInterest()");
    if(evdlayoutopt->fAutoZoomInterest == 1) fToggleAutoZoom->SetState(kButtonDown);

    fVFrame->AddFrame(fZoomInterest,     new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fUnZoomInterest,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );

    fVFrame->AddFrame(fZoomBack,         new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );

    fVFrame->AddFrame(fToggleAutoZoom,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
  }


  //----------------------------------------------------------------------------
  void	TWQMultiTPCProjectionView::RadioButtonsDispatch(int parameter)
  {
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt; 
    if(parameter==1 || parameter == 2){
      fToggleZoom->SetState(kButtonUp); 	 
    }
 
  }

  //......................................................................
  void TWQMultiTPCProjectionView::SetUpPositionFind()
  {
    // enter zoom buttons
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  
    if(!evdlayoutopt->fShowEndPointSection)            
      return;

    // int    	 fShowEndPointMarkers;             ///< Draw EndPoint Markers if clicked. 

    fFindEndpoint=new TGTextButton(fVFrame,"&Find XYZ",150);
    fFindEndpoint->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "FindEndPoint()");

    fXYZPosition=new TGTextView(fVFrame,100,55,999,TGView::kNoHSB | TGView::kNoVSB); ///< Display the xyz position
    fXYZPosition->SetEditable("false");
    TGText *tt=new TGText("x,y,z");
    fXYZPosition->SetText(tt);


    fClearPPoints=new TGTextButton(fVFrame,"&Clear Points",150);
    fClearPPoints->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "ClearEndPoints()");  // ?

    fToggleShowMarkers=new TGCheckButton(fVFrame,"ShowMarkers",0);       ///< Toggle the ShowEndPointMarkers Setting 
    fToggleShowMarkers->Connect("Clicked()", "evd::TWQMultiTPCProjectionView", this, "ToggleEndPointMarkers()");
    if(evdlayoutopt->fShowEndPointMarkers == 1) fToggleShowMarkers->SetState(kButtonDown);

    fVFrame->AddFrame(fFindEndpoint,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fXYZPosition,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fClearPPoints,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fToggleShowMarkers,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
  }


  /////////////////////////////////////////
  //  Go back one step in zoom

  void    TWQMultiTPCProjectionView::ZoomBack()
  {
    if(fPrevZoomOpt.size()>0){
      ZoomOptionsMultiTPC ThePrevZoomOpt = fPrevZoomOpt.at(fPrevZoomOpt.size()-1);
      int plane = fZoomOpt.OnlyPlaneChanged;
      if(plane != -1){
	SetZoom(plane, 
		ThePrevZoomOpt.wmin[plane],
		ThePrevZoomOpt.wmax[plane],
		ThePrevZoomOpt.tmin[plane],
		ThePrevZoomOpt.tmax[plane],
		false);
      }
      else{
	for( size_t iplane = 0; iplane != fPlanes.size(); ++iplane){
	  SetZoom(iplane, 
		  ThePrevZoomOpt.wmin[iplane],
		  ThePrevZoomOpt.wmax[iplane],
		  ThePrevZoomOpt.tmin[iplane],
		  ThePrevZoomOpt.tmax[iplane],
		  false);
	      
	}
      }

      fPrevZoomOpt.pop_back();
    }
    else
      mf::LogVerbatim("TWQMultiTPCProjectionView") <<"unable to unzoom further - no zoom settings left on stack"<<std::endl;
  }

  //------------------------------------
  void    TWQMultiTPCProjectionView::SetZoom(int plane,
					     int wirelow,
					     int wirehi,
					     int timelow,
					     int timehi, 
					     bool StoreZoom)
  {

    if(StoreZoom){
      fPrevZoomOpt.push_back(fZoomOpt);
      fZoomOpt.OnlyPlaneChanged = plane;
    }

    fZoomOpt.wmin[plane]    = wirelow;
    fZoomOpt.wmax[plane]    = wirehi;
    fZoomOpt.tmin[plane]    = timelow;
    fZoomOpt.tmax[plane]    = timehi;


    TVirtualPad *ori = gPad;
    zoom_opt="1";

    // error checking - useful for the mouse zoom.
    if(wirehi<wirelow){
      int temp=wirelow;
      wirelow=wirehi;
      wirehi=temp;
    }
  
    if(timehi<timelow){
      int temp=timelow;
      timelow=timehi;
      timehi=temp;
    }
  
    //if drawing, then currently not zooming
    curr_zooming_plane=-1;
  
    fPlanes[plane]->SetZoomRange(wirelow, wirehi,timelow,timehi);
    fPlanes[plane]->Draw("1");
    fPlanes[plane]->UpdatePad();
  
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    //  UpdateSeedCurve();
  
    ori->cd();
  
    return;  
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetPlaneWire()
  {
    TVirtualPad *ori = gPad;

    fWireQ->SetPlaneWire(kPlane, kWire);

    fWireQ->Draw();
    fWireQ->Pad()->cd();
    fWireQ->Pad()->Modified();
    fWireQ->Pad()->Update();
    fWireQ->Pad()->SetBit(TPad::kCannotMove,true);
    fWireQ->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);

    fPlaneEntry->SetNumber(kPlane);
    fWireEntry->SetNumber(kWire);

    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetPlane()
  {
    kPlane = (unsigned int)fPlaneEntry->GetNumberEntry()->GetNumber();

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetWire()
  {
    kWire = (unsigned int)fWireEntry->GetNumberEntry()->GetNumber();

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetDistance()
  {
    kDistance = (double)fDistance->GetNumberEntry()->GetNumber();
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetThreshold()
  {
    double threshold = fThresEntry->GetNumberEntry()->GetNumber();

    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    rawopt->fMinSignal = threshold;

    TVirtualPad *ori = gPad;
    this->DrawPads(zoom_opt);
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetGreyscale()
  {
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    TGButton *b = (TGButton *)gTQSender;
    if(b->GetState() == kButtonDown){
      cst->fColorOrGray = 1;
    }
    else{
      cst->fColorOrGray = 0;
    }

    TVirtualPad *ori = gPad;
    this->DrawPads(zoom_opt);
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetRawCalib()
  {
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;

    TGButton *b = (TGButton *)gTQSender;
    int id = b->WidgetId();

    // id values are set in lines 125 - 127
    if(id == 4){
      rawopt->fDrawRawDataOrCalibWires = 0;
      fRawDraw->SetState(kButtonDown);
      fCalibDraw->SetState(kButtonUp);
      fRawCalibDraw->SetState(kButtonUp);
    }
    else if(id == 3){
      rawopt->fDrawRawDataOrCalibWires = 1;
      fRawDraw->SetState(kButtonUp);
      fCalibDraw->SetState(kButtonDown);
      fRawCalibDraw->SetState(kButtonUp);
    }
    else if(id == 2){
      rawopt->fDrawRawDataOrCalibWires = 2;
      fRawDraw->SetState(kButtonUp);
      fCalibDraw->SetState(kButtonUp);
      fRawCalibDraw->SetState(kButtonDown);
    }

    TVirtualPad *ori = gPad;

    fWireQ->Draw();
    fWireQ->Pad()->cd();
    fWireQ->Pad()->Modified();
    fWireQ->Pad()->Update();

    this->DrawPads(zoom_opt);
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQMultiTPCProjectionView::SetMCInfo()
  {
    art::ServiceHandle<evd::SimulationDrawingOptions> sdo;

    TGButton *b = (TGButton *)gTQSender;
    if(b->GetState() == kButtonDown){
      sdo->fShowMCTruthText    = 1;
      sdo->fShowMCTruthVectors = 1;
    }
    else{
      sdo->fShowMCTruthText    = 0;
      sdo->fShowMCTruthVectors = 0;
    }

    TVirtualPad *ori = gPad;

    fMC->Draw();
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();
  }


}// namespace
