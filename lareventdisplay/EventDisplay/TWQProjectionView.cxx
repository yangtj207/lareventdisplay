//
/// \file    TWQProjectionView.cxx
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
#include "EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/ChangeTrackers.h" // util::DataProductChangeTracker_t
#include "lareventdisplay/EventDisplay/TWQProjectionView.h"
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
#include "lareventdisplay/EventDisplay/InfoTransfer.h"

#include "lardata/RecoBase/Seed.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/Utilities/GeometryUtilities.h"

namespace evd{

  static unsigned int kPlane;
  static unsigned int kWire;
  static double       kDistance;
  static int          curr_zooming_plane;
  static const char*  zoom_opt=0;
  
  static int          shift_lock;
  

  //......................................................................
  TWQProjectionView::TWQProjectionView(TGMainFrame* mf) 
    : evdb::Canvas(mf)
    , fRedraw(nullptr)
    , fCryoInput(nullptr), fTPCInput(nullptr), fTotalTPCLabel(nullptr)
    , isZoomAutomatic
        (art::ServiceHandle<evd::EvdLayoutOptions>()->fAutoZoomInterest)
    , fLastEvent(new util::DataProductChangeTracker_t)
  {  
    
    art::ServiceHandle<geo::Geometry> geo;

    // first make pads for things that don't depend on the number of 
    // planes in the detector
    // bottom left corner is (0.,0.), top right is  (1., 1.)
    fAngleInfo=NULL;
    fXYZPosition=NULL;
      
    fLastThreshold = -1.;
    
    evdb::Canvas::fCanvas->cd();  
    fHeaderPad = new HeaderPad("fHeaderPad","Header",0.0,0.0,0.15,0.13,"");  
    fHeaderPad->Draw();  

    evdb::Canvas::fCanvas->cd();  
    fMC = new MCBriefPad("fMCPad","MC Info.",0.15,0.13,1.0,0.17,"");  
    fMC->Draw();  

    evdb::Canvas::fCanvas->cd();  
//    fWireQ = new TQPad("fWireQPad", "ADCvsTime",0.15,0.0,1.0,0.13,"TQ", 0, 0);  
    fWireQ = new TQPad("fWireQPad", "ADCvsTime",0.15,0.0,1.0,0.14,"TQ", 0, 0);  
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
				    0 , geo->Nplanes()-1 );

    kPlane = 0;
    kWire = TMath::Nint(0.5*geo->Nwires(0));
    kDistance=1.5;
    fWireQ->SetPlaneWire(kPlane, kWire);

    // Initial value
    fPlaneEntry->SetNumber( kPlane );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fPlaneEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetPlane()");
    fPlaneEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetPlane()");
    // Text label for this numeric field.
    fPlaneLabel= new TGLabel(fFrame,"Plane");

    // wire number entry 
    unsigned int maxwire=0;
    for (unsigned int ip=0;ip<geo->Nplanes();ip++)
      maxwire = (geo->Nwires(ip)-1 > maxwire) ? geo->Nwires(ip)-1 : maxwire;
    
    
    fWireEntry = new TGNumberEntry(fFrame,0,6,-1,
				   TGNumberFormat::kNESInteger, 
				   TGNumberFormat::kNEAAnyNumber, 
				   TGNumberFormat::kNELLimitMinMax, 
				   0 , maxwire );
    // Initial value
    fWireEntry->SetNumber( kWire );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fWireEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetWire()");
    fWireEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetWire()");

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
    fThresEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetThreshold()");
    fThresEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetThreshold()");

    // Text label for this numeric field.
    fThresLabel= new TGLabel(fFrame,"ADC Threshold");

    // check button to toggle color vs grey
    fGreyScale = new TGCheckButton(fFrame,"Grayscale",1);
    fGreyScale->Connect("Clicked()", "evd::TWQProjectionView", this, "SetGreyscale()");
    if(cst->fColorOrGray == 1) fGreyScale->SetState(kButtonDown);

    // check button to toggle MC information
    if(evdlayoutopt->fEnableMCTruthCheckBox){
      fMCOn = new TGCheckButton(fFrame,"MC Truth",5);
      fMCOn->Connect("Clicked()", "evd::TWQProjectionView", this, "SetMCInfo()");
      if(sdo->fShowMCTruthText == 1) fMCOn->SetState(kButtonDown);
    }

    // radio buttons to toggle drawing raw vs calibrated information
    fRawCalibDraw = new TGRadioButton(fFrame,"Both",          2);
    fCalibDraw    = new TGRadioButton(fFrame,"Reconstructed", 3);
    fRawDraw      = new TGRadioButton(fFrame,"Raw",           4);
    fRawDraw     ->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
    fCalibDraw   ->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
    fRawCalibDraw->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
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

    // geometry to figure out the number of planes
    unsigned int nplanes = geo->Nplanes();

    if(evdlayoutopt-> fShowSideBar)
      SetUpSideBar();    
    else
      evdlayoutopt->fShowEndPointSection=0;  // zero it to avoid a misconfiguration in the fcl file. 

    //zero the ppoints queue.
    ppoints.clear();
    pline.clear();
    seedlines.clear();
       

    // now determine the positions of all the time vs wire number 
    // and charge histograms for the planes
    for(unsigned int i = 0; i < nplanes; ++i){
      double twx1 = 0.;
      double twx2 = 0.97;
      double twx3 = 1.0;
      double twy1 = 0.17 +   (i)*(1.0-0.171)/(1.*nplanes);
      double twy2 = 0.17 + (i+1)*(1.0-0.171)/(1.*nplanes);



      TString padname = "fWireProjP";
      padname += i;

      TString padtitle = "Plane";
      padtitle += i;

      evdb::Canvas::fCanvas->cd();
      fPlanes.push_back(new TWireProjPad(padname, padtitle, twx1, twy1, twx2, twy2, i));
      fPlanes[i]->Draw();
//      fPlanes[i]->Pad()->AddExec("mousedispatch",Form("evd::TWQProjectionView::MouseDispatch(%d, (void*)%d)", i, this)); 
      fPlanes[i]->Pad()->AddExec("mousedispatch",Form("evd::TWQProjectionView::MouseDispatch(%d, (void*)%lu)", i, (unsigned long) this)); 

      padname = "fQPadPlane";
      padname += i;

      padtitle = "QPlane";
      padtitle += i;

      evdb::Canvas::fCanvas->cd();
      fPlaneQ.push_back(new TQPad(padname, padtitle, twx2, twy1, twx3, twy2, "Q", i, 0));
      fPlaneQ[i]->Draw();

    } 

    // propagate the zoom setting
    SetAutomaticZoomMode(isZoomAutomatic);

    evdb::Canvas::fCanvas->Update();

  }

  //......................................................................
  TWQProjectionView::~TWQProjectionView() 
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
    
    delete fLastEvent;
  }

  //......................................................................
  void TWQProjectionView::ResetRegionsOfInterest() {
    for(TWireProjPad* planePad: fPlanes)
      planePad->RawDataDraw()->ResetRegionOfInterest();
  } // TWQProjectionView::ResetRegionsOfInterest()
  
  //......................................................................
  void TWQProjectionView::DrawPads(const char* /*opt*/)
  {
    
    OnNewEvent(); // if the current event is a new one, we need some resetting
    
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
  void TWQProjectionView::SetAutomaticZoomMode(bool bSet /* = true */) {
    isZoomAutomatic = bSet;
    for (TWireProjPad* pPlane: fPlanes)
      pPlane->SetZoomToRoI(isZoomAutomatic);
  } // TWQProjectionView::SetAutomaticZoomMode()

  
  //......................................................................
  void TWQProjectionView::Draw(const char* opt) 
  {  
    mf::LogDebug("TWQProjectionView") << "Starting to draw";
    
    OnNewEvent(); // if the current event is a new one, we need some resetting
    
    art::ServiceHandle<geo::Geometry> geo;

    ClearAllSeeds();
    fPrevZoomOpt.clear();
  
  
    evdb::Canvas::fCanvas->cd();    
    zoom_opt=0;
    fHeaderPad->Draw();    
    fMC       ->Draw(); 
    fWireQ    ->Draw();
  
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
    if(evdlayoutopt->fPrintTotalCharge)
      PrintCharge();
  
    //clear queue of selected points
    ppoints.clear();
    pline.clear();
    seedlines.clear();
    // Reset current zooming plane - since it's not currently zooming.
    curr_zooming_plane=-1;
  
    unsigned int const nPlanes = fPlanes.size();
    LOG_DEBUG("TWQProjectionView") << "Start drawing " << nPlanes << " planes";
    //  double Charge=0, ConvCharge=0;
    for(unsigned int i=0;i<nPlanes;++i){
      TWireProjPad* planePad = fPlanes[i];
      planePad->Draw(opt);
      planePad->Pad()->Update();
      planePad->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);
      fPlaneQ[i]->Draw();
      std::vector<double> ZoomParams = planePad->GetCurrentZoom();
      fZoomOpt.wmin[i] = ZoomParams[0];
      fZoomOpt.wmax[i] = ZoomParams[1];
      fZoomOpt.tmin[i] = ZoomParams[2];
      fZoomOpt.tmax[i] = ZoomParams[3];
      // Charge deposit feature - not working yet
      // 
      //  if(geo->Plane(i).SignalType()==geo::kCollection)
      //  {
      //	planePad->RecoBaseDraw()->GetChargeSum(i,Charge,ConvCharge);  
      //   }
   
    }
    mf::LogDebug("TWQProjectionView")
      << "Done drawing " << nPlanes << " planes";
  

    // Charge deposit feature - not working yet
    //  std::stringstream ss;
    // if(ConvCharge!=0)
    // {
    //  ss << ConvCharge << "MeV"<<std::endl;
    //  }
    // else
    //  {
    //    ss<<" no reco info";
    //  }
    //	 	
    //  TGText * tt = new TGText(ss.str().c_str());
    //  tt->InsLine(1, "Approx EDep:");
    //  fAngleInfo->SetText(tt);
    //
    //  ss.flush();
    //

    // Reset any text boxes which are enabled

    if(fXYZPosition)
      fXYZPosition->SetForegroundColor(kBlack);

    if(fAngleInfo)
      fAngleInfo->SetForegroundColor(kBlack);

    evdb::Canvas::fCanvas->Update();
    mf::LogDebug("TWQProjectionView") << "Done drawing";
  }

  // comment out this method as for now we don't want to change every
  // plane to have the same range in wire number because wire numbers
  // don't necessarily overlap from plane to plane, ie the same range
  // isn't appropriate for every plane
  //......................................................................
  //   void TWQProjectionView::RangeChanged() 
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
  void 	TWQProjectionView::SetTestFlag(int number)
  {
    art::ServiceHandle<evd::InfoTransfer>   infot; 
    infot->SetTestFlag(number);  
  }
  
  
  //......................................................................
  void 	TWQProjectionView::PrintCharge()
  {

    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;

    for(size_t iplane = 0; iplane < fPlanes.size(); ++iplane){ 
      if(geo->Plane(iplane).SignalType()==geo::kCollection){
	double ch=0,convch=0;
	if(rawopt->fDrawRawDataOrCalibWires == 0){
	  fPlanes[iplane]->RawDataDraw()->GetChargeSum(iplane,ch,convch);
	  mf::LogVerbatim("TWQProjectionView") << "Warning! Calculating for RawData! ";
	}
	else{  
	  fPlanes[iplane]->RecoBaseDraw()->GetChargeSum(iplane,ch,convch);  
	}    

	mf::LogVerbatim("TWQProjectionView") << "\ncharge collected at collection plane: " 
					     << iplane << " " << ch << " " << convch;
      }
    }


  }

  //-------------------------------------------------------------------
  //......................................................................
  void TWQProjectionView::MouseDispatch(int plane, void * wqpv)
  {
    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();
    evd::TWQProjectionView *wqpp = (evd::TWQProjectionView*)wqpv;
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;

    switch (event){
	
    case kButton1Shift:           	
      shift_lock=1;	
      if(evdlayoutopt->fMakeClusters==1){wqpp->SelectHit(plane);}
      else {wqpp->SelectPoint(plane);}
      break;
    case kButton1Up:    
      if(shift_lock==1) break;
      if(evdlayoutopt-> fChangeWire==1) wqpp->ChangeWire(plane);		
    case kButton1Down: shift_lock=0;
    case kButton1Motion:	
      if(evdlayoutopt->fMakeClusters == 1 && !evdlayoutopt->fMakeSeeds == 1){ wqpp->SetClusters(plane);}
      else if(evdlayoutopt->fMakeClusters == 0 && evdlayoutopt->fMakeSeeds == 1){ wqpp->SetSeeds(plane);}
      else { wqpp->SetMouseZoomRegion(plane);}
      break;
      //  default:		
    }
  }


  //......................................................................
  void TWQProjectionView::ChangeWire(int plane)
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
  void TWQProjectionView::SelectPoint(int plane)
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
  void TWQProjectionView::ClearEndPoints()
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
  double TWQProjectionView::FindLineLength()
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
		
      //find wireIDs corresponding to found wires.
      geo::WireID wire1(rawOpt->fCryostat,rawOpt->fTPC,pline[0].plane,pline[0].w0);
      geo::WireID wire2(rawOpt->fCryostat,rawOpt->fTPC,pline[1].plane,pline[1].w0);
      
      bool wires_cross=false;
      bool time_good=false;
	
      if(std::abs(pline[0].t0-pline[1].t0) < 200){
	geo::WireIDIntersection widIntersect;
	wires_cross = geom->WireIDsIntersect(wire1,wire2,widIntersect);
	y = widIntersect.y;
	z = widIntersect.z;
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
      //find wireIDs corresponding to found wires AT END OF LINE.
      wire1.Wire = pline[0].w1;
      wire2.Wire = pline[1].w1;

      wires_cross=false;
      time_good=false;
	
      if(std::abs(pline[0].t1-pline[1].t1) < 200){
	geo::WireIDIntersection widIntersect;
	wires_cross = geom->WireIDsIntersect(wire1,wire2,widIntersect);
	y = widIntersect.y;
	z = widIntersect.z;
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
  void TWQProjectionView::FindEndPoint()
  {
    // if list is larger than or equal to two, can project to XYZ and extrapolate to third plane (if exists)

    ////for now leaving commented. At some point be useful to display this information somewhere
    // for(unsigned int ix=0;ix<ppoints.size();ix++)
    //   std::cout << "ppoints, planes,x,y :" << ix << " " << ppoints[ix].plane << " " << ppoints[ix].x << " " << ppoints[ix].y << std::endl;

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
      //double ftimetick = detp->SamplingRate()/1000.;
      //double larv = detp->DriftVelocity(detp->Efield(), detp->Temperature());
		
      //find channels corresponding to found wires.
      geo::WireID wire1(rawOpt->fCryostat,rawOpt->fTPC,ppoints[0].plane,ppoints[0].w);
      geo::WireID wire2(rawOpt->fCryostat,rawOpt->fTPC,ppoints[1].plane,ppoints[1].w);
 
      bool wires_cross=false;
      bool time_good=false;
	
      if(std::abs(ppoints[0].t-ppoints[1].t) < 200){
	geo::WireIDIntersection widIntersect;
	geom->WireIDsIntersect(wire1,wire2,widIntersect);
	y = widIntersect.y;
	z = widIntersect.z;
	wires_cross=true;
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
	
	xyz_vertex_fit[0]=detp->ConvertTicksToX(ppoints[0].t,
						ppoints[0].plane,
						rawOpt->fTPC,
						rawOpt->fCryostat);
	second_time=detp->ConvertTicksToX(ppoints[1].t,
					  ppoints[1].plane,
					  rawOpt->fTPC,
					  rawOpt->fCryostat);
	
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
	
	double timestart = detp->ConvertXToTicks(xyz_vertex_fit[0],
						 wplane,
						 rawOpt->fTPC,
						 rawOpt->fCryostat);

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

  //......................................................................
  void TWQProjectionView::RefitSeeds(){
  
  
  }

  //......................................................................
  // SaveSelection
  void TWQProjectionView::SaveSelection()
  {
    util::GeometryUtilities gser;
	
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutoptions;
    if(evdlayoutoptions->fMakeSeeds){
      double KineticEnergy = fPlanes[0]->SaveSeedList( seedlines, kDistance);
      std::stringstream ss;
      ss<< " " << std::setprecision(1) << std::fixed << KineticEnergy << " MeV";
	    
      TGText * tt = new TGText(fAngleInfo->GetText());
      tt->InsLine(1,ss.str().c_str());
      fAngleInfo->SetText(tt);
      fAngleInfo->Update();
    }
    else if(evdlayoutoptions->fMakeClusters){
      //only calculating in 3 planes now, needs to be generalized eventually
      //	evd::TWQProjectionView tt;
      double omx[3];
      //double omx0,omx1;
      double xphi;
      double xtheta;
      unsigned int ii;
      if(pline.size()<2) {
	TGText *tt=new TGText("not enough lines selected");
	fAngleInfo->SetText(tt);
	fAngleInfo->Update();
	return;
      }
      double deltawire;
      double deltatime;
      for(ii = 0; ii < pline.size(); ++ii){
	//if(pline[ii].plane==0) break;
	deltawire=pline[ii].w1-pline[ii].w0;
	deltatime=pline[ii].t1-pline[ii].t0;
	omx[ii]=gser.Get2Dangle(deltawire,deltatime);
      }

      for(size_t ii = 0; ii < pline.size(); ++ii){
	fPlanes[pline[ii].plane]->SaveHitList(pline[ii].w0,
					      pline[ii].t0,
					      pline[ii].w1,
					      pline[ii].t1,
					      kDistance, 
					      zoom_opt);
      }
      if(fPlanes.size()>pline.size() && pline.size() >=2 ){   // need to project to third plane

		
	util::PxPoint p00(pline[0].plane,pline[0].w0,pline[0].t0);
	util::PxPoint p01(pline[1].plane,pline[1].w0,pline[1].t0);
	util::PxPoint p0N(0,0,0);
	int error1=gser.GetProjectedPoint(&p00,&p01,p0N);
		
				
	util::PxPoint p10(pline[0].plane,pline[0].w1,pline[0].t1);
	util::PxPoint p11(pline[1].plane,pline[1].w1,pline[1].t1);
	util::PxPoint p1N(0,0,0);
	int error2=gser.GetProjectedPoint(&p10,&p11,p1N);
	if(error1!=-1 && error2!=-1)
	  fPlanes[p0N.plane]->SaveHitList(p0N.w,p0N.t,p1N.w,p1N.t,kDistance, zoom_opt,false);
		
		
				
      }
	
      for(size_t jj = 0;jj < fPlanes.size(); ++jj){
	fPlanes[jj]->UpdatePad();
      }

	
      gser.Get3DaxisN(pline[0].plane,pline[1].plane,omx[0],omx[1],xphi,xtheta);

//       double xtheta2;
//       if(xphi < 2 && xphi > -2)
// 	xtheta2=gser.Get3DSpecialCaseTheta(pline[0].plane,
// 					   pline[1].plane,
// 					   pline[0].w1-pline[0].w0, 
// 					   pline[1].w1-pline[1].w0);
	
      double length=FindLineLength();
      TGText *tt=new TGText(Form("Length:%4.1f",length));
      tt->InsLine(1,Form("Omega P%d:%4.1f,",pline[0].plane,omx[0])); 
      tt->InsLine(2,Form("Omega P%d:%4.1f,",pline[1].plane,omx[1]));  
      tt->InsLine(3,Form("Phi: %4.1f,",xphi));

      tt->InsLine(4, Form("Theta: %4.1f",xtheta));
      fAngleInfo->SetText(tt);
      fAngleInfo->Update();
    } // end else if
  }

  //.......................................................................
  void TWQProjectionView::ClearSelection()
  {
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
  
    if(!evdlayoutopt->fMakeSeeds && !evdlayoutopt->fMakeClusters)
      ppoints.clear();
    else if(evdlayoutopt->fMakeSeeds)
      ClearAllSeeds();
    else if(evdlayoutopt->fMakeClusters){
      if(this->pline.size()==0) return;
      for(size_t i = 0;i < fPlanes.size(); ++i){
	fPlanes[i]->ClearHitList();
	fPlanes[i]->UpdatePad();
      }
      pline.clear();
	
    }	
  }

  //.......................................................................
  void TWQProjectionView::ClearAllSeeds()
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if(!recoopt->fUseHitSelector) return;

    if(this->seedlines.size()==0) return;
    //this->DrawPads();
    seedlines.clear();
    fPlanes[0]->HitSelectorGet()->SeedVector().clear();
    // only for last plane
    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->DrawLinesinView(seedlines,true);
      fPlanes[i]->UpdatePad();
    }
	
    return;
  }



  //.......................................................................
  // If we are mid seed creation, stop and clear.
  // Otherwise, remove the last seed.
  //
  void TWQProjectionView::ClearLastSeed()
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if(!recoopt->fUseHitSelector) return;

    if(this->seedlines.size()==0) return;

    // If the seed lines vector is an event seed number, remove
    //  one seed from the hit selector seed collection  
    if(seedlines.size()%3==0)
      fPlanes[0]->HitSelectorGet()->SeedVector().pop_back();

    // Figure out how many display lines need to be rejected
    int SeedsToRemove = seedlines.size() % fPlanes.size();
    if(SeedsToRemove==0) SeedsToRemove=3;

    mf::LogVerbatim("TWQProjectionView") <<"Removing  " 
					 << SeedsToRemove 
					 << " of " << seedlines.size() <<" seed lines";
  
    // Throw them out
    for(int i=0; i!=SeedsToRemove; ++i)
      seedlines.pop_back();
  
    // Then update the drawing pads
    for(size_t i=0; i!=fPlanes.size(); ++i){
      fPlanes[i]->UpdatePad();
      fPlanes[i]->DrawLinesinView(seedlines,true);
      fPlanes[i]->UpdatePad();
    }
  
  }

  //.......................................................................
  void TWQProjectionView::SetMouseZoomRegion(int plane)
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
	 ( std::abs(wend-wstart ) > 0.01*(xx2-xx1) ) && 
	 ( std::abs(tend-tstart ) > 0.01*(yy2-yy1)   &&  
	   curr_zooming_plane==plane ) ){
	
	SetAutomaticZoomMode(false);
	this->SetZoom(plane,wstart,wend,tstart,tend);
	wstart=-1;
	tstart=-1;
      }
      break;
    }
    }// end switch
  }

  //......................................................................
  int TWQProjectionView::DrawLine(int plane,util::PxLine &pline)
  {
    art::ServiceHandle<evd::EvdLayoutOptions>     evdlayoutopt;
    const detinfo::DetectorProperties* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
    static Float_t w0=-1, t0=-1, w1=-1, t1=-1;

    static Int_t pxold, pyold;
    static Int_t pw0, pt0;
  
    static Int_t linedrawn;
  
    // These static variables keep track of seed lines which 
    //  were drawn
    static Int_t LastSeedt0;
    static Int_t LastSeedt1;
    static Int_t LastSeedPlane;


    //Temporary evil - eventually from pset
    double MaxSeedTimeDiff=100;

    int event = gPad->GetEvent();
    int px = gPad->GetEventX();
    int py = gPad->GetEventY();
  
    int SeedCounter = seedlines.size();


    int linefinished=0;
  
    switch (event){

    case kButton1Down:{
      //not doing anything right now
      w0 = gPad->AbsPixeltoX(px);
      t0 = gPad->AbsPixeltoY(py);
      if(evdlayoutopt->fMakeSeeds){
	if((SeedCounter%3)==1){
	  t0=det->ConvertXToTicks(det->ConvertTicksToX(LastSeedt0,LastSeedPlane,0,0),plane,0,0);				
	  py=gPad->YtoAbsPixel(t0);
	}
	else{
	  LastSeedt0=t0;
	  LastSeedPlane=plane;
	}
      }
      pw0   = px; pt0   = py;
      pxold = px; pyold = py;
      linedrawn = 0;
      curr_zooming_plane=plane;
      break;
    }

    case kButton1Motion:{ 
      int lx,hx,ly,hy;
		  
      // If we are in seed mode, and one seed line
      //  was already placed, constrain head of next line
      //  to be at same t0

      if((evdlayoutopt->fMakeSeeds)&&((SeedCounter%3)==1)){
	t0 = det->ConvertXToTicks(det->ConvertTicksToX(LastSeedt0, LastSeedPlane,0,0),plane,0,0);
	t1 = det->ConvertXToTicks(det->ConvertTicksToX(LastSeedt1, LastSeedPlane,0,0),plane,0,0);
	pt0=gPad->YtoAbsPixel(t0);
      }

      lx=pxold;
      hx=pw0;
		  
      ly=pyold;
      hy=pt0;
		  
		  

      if (linedrawn) gVirtualX->DrawLine(lx, ly, hx, hy);
		  
      pxold = px;
      pyold = py;
      linedrawn = 1;
		  
      lx=pxold;
      hx=pw0;
		  
      ly=pyold;
      hy=pt0;
		  
      if (linedrawn) gVirtualX->DrawLine(lx, ly, hx, hy);
      break;
    }

    case kButton1Up:{
      if (px == pw0 && py == pt0) break;
      w1 = gPad->AbsPixeltoX(px);
      t1 = gPad->AbsPixeltoY(py);
			
      gPad->Modified(kTRUE);

      // if we are making seeds, process this mouse operation
      // in context of how many unconnected seeds we have
      if(evdlayoutopt->fMakeSeeds){
	LastSeedt1=t1;
			    
	// we keep track of how many seed lines we placed
	SeedCounter++;

	// If this is a second seed line, it should match an
	// existing line.  \todo: Add wire check
	if((SeedCounter%3)==2){
	  if(LastSeedPlane==plane){
	    ClearLastSeed();
	    mf::LogVerbatim("TWQProjectionView") <<"Cannot draw seed with 2 lines in same plane!";
	    linefinished=0;
	    return 0;
	  }
	  if(((LastSeedt0-t0)>MaxSeedTimeDiff)||(LastSeedt1-t1)>MaxSeedTimeDiff){
	    mf::LogVerbatim("TWQProjectionView")<<"This seed line does not fit the last proj. not adding it";
	    linefinished=0;
	    SeedCounter--;
	  }
	  else{
	    // This proj is compatible to make a 3D seed.  Indicate this with the linefinished status
	    mf::LogVerbatim("TWQProjectionView")<<"This is a good seed. adding line and returing 2";
	    pline=util::PxLine(plane,w0,t0,w1,t1);
	    linefinished=2;
	  }
				  

	}
	// If a first seed line, can be drawn with no constraints.
	// Put down guide lines in other 2 views for t
			    
	else if((SeedCounter%3)==1){
	  pline=util::PxLine(plane,w0,t0,w1,t1);
	  linefinished=1;
				
	}
	// Of it is a thid seed line, something is terribly
	//  wrong.  Shout a warning.
	else if((SeedCounter%3)==0){
	  mf::LogVerbatim("TWQProjectionView") <<"You have been eaten by a Grue. And seed mode is doing something wrong.";
	}
      }
      else{
			  
	pline=util::PxLine(plane,w0,t0,w1,t1);
	linefinished=1;
      }
      //curr_zooming_plane=-1;
      //gROOT->SetEditorMode();
    }
    }   //end switch
  
    return linefinished;
  }

  //.......................................................................
  void TWQProjectionView::SetClusters(int plane)
  {

    TObject *select = gPad->GetSelected();
    if(!select) return;
    if(!select->InheritsFrom("TBox")) return;


    util::PxLine ppx;
	 
    if(!DrawLine(plane,ppx))
      return;
	 
    curr_zooming_plane=-1;
    gROOT->SetEditorMode();

    // check if not clicking on a plane that is already in the ppoints list:
    int repeat_plane=-1;

    for(size_t ii = 0; ii < this->pline.size(); ++ii){
      if(ppx.plane==this->pline[ii].plane){
	this->pline[ii]=ppx;
		
	//clear View and draw new Marker
	this->fPlanes[plane]->Pad()->cd();
	this->fPlanes[this->pline[ii].plane]->View()->Clear();
	this->fPlanes[this->pline[ii].plane]->View()->Draw();
	TLine& l = this->fPlanes[this->pline[ii].plane]->View()->AddLine(pline[ii].w0,
									 pline[ii].t0,
									 pline[ii].w1,
									 pline[ii].t1);
	evd::Style::FromPDG(l,11);
	this->fPlanes[this->pline[ii].plane]->View()->Draw();
	repeat_plane=this->pline[ii].plane;
	break;
      }
    }
	
    //if plane does not repeat and size of list is larger than 2 pop_front and delete its marker. Otherwise just push_back.
    if(repeat_plane==-1){
      if( this->pline.size()>=2){
	this->fPlanes[this->pline[0].plane]->Pad()->cd();
	this->fPlanes[this->pline[0].plane]->View()->Clear();
	this->fPlanes[this->pline[0].plane]->View()->Draw();
	this->fPlanes[this->pline[0].plane]->Pad()->Clear();
	this->fPlanes[this->pline[0].plane]->Pad()->Draw();
	this->fPlanes[this->pline[0].plane]->Draw();
	this->fPlanes[this->pline[0].plane]->Pad()->Modified();
	this->fPlanes[this->pline[0].plane]->Pad()->Update();
	this->fPlanes[this->pline[0].plane]->Pad()->GetFrame()->SetBit(TPad::kCannotMove,true);
	this->fPlanes[this->pline[0].plane]->Pad()->SetBit(TPad::kCannotMove,true);
	this->pline.pop_front();
      }

      this->pline.push_back(ppx);
      this->fPlanes[plane]->Pad()->cd();
      this->fPlanes[plane]->View()->Clear();
      this->fPlanes[plane]->View()->Draw();

      int size=pline.size()-1;
      TLine& l = this->fPlanes[this->pline[size].plane]->View()->AddLine(pline[size].w0,
									 pline[size].t0,
									 pline[size].w1,
									 pline[size].t1);
      this->fPlanes[this->pline[size].plane]->View()->Draw();
      evd::Style::FromPDG(l,11);


    }
	
  }

  //.......................................................................
  void TWQProjectionView::SetSeeds(int plane)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if(!recoopt->fUseHitSelector) return;

    TObject *select = gPad->GetSelected();
    if(!select) return;
    if(!select->InheritsFrom("TBox")) return;
  
  
    util::PxLine ppx;
  
    int DrawStatus = DrawLine(plane,ppx);
  
    //  std::cout<<"Draw status return value for this line " <<DrawStatus<<std::endl;

    if(DrawStatus==0)
      return;
  
    curr_zooming_plane=-1;
    gROOT->SetEditorMode();
  
  
    this->seedlines.push_back(ppx);
  
    if(DrawStatus==1){
      // We drew the first of seed projection - go ahead and pass it on, and draq guide lines
      for(size_t p = 0; p != fPlanes.size(); ++p)
	this->fPlanes[p]->DrawLinesinView(seedlines);
    }
    else if(DrawStatus==2){
      mf::LogVerbatim("TWQProjectionView") << "Adding third seed projection";
      // We drew the second seed projection - make the 3D seed
      
      util::PxLine l0 = seedlines.at(seedlines.size()-1);
      util::PxLine l1 = seedlines.at(seedlines.size()-2);
      
      const detinfo::DetectorProperties* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
      
      double x0 = det->ConvertTicksToX(l0.t0, l0.plane, 0,0);
      double x1 = det->ConvertTicksToX(l0.t1, l0.plane, 0,0);

      

      util::PxPoint thirdp0(0,0,0), thirdp1(0,0,0);
      
      util::GeometryUtilities geomutil;
      
      double yz0[2];
      double yz1[2];
      util::PxPoint l0_pt0 = l0.pt0();
      util::PxPoint l0_pt1 = l0.pt1();
      util::PxPoint l1_pt0 = l1.pt0();
      util::PxPoint l1_pt1 = l1.pt1();
      geomutil.GetProjectedPoint(&l0_pt0, &l1_pt0, thirdp0);
      geomutil.GetProjectedPoint(&l0_pt1, &l1_pt1, thirdp1);
      geomutil.GetYZ(&l0_pt0, &l1_pt0, yz0);
      geomutil.GetYZ(&l0_pt1, &l1_pt1, yz1);
	     
      thirdp0.t = det->ConvertXToTicks(x0, thirdp0.plane,0,0);
      thirdp1.t = det->ConvertXToTicks(x1, thirdp1.plane,0,0);

      util::PxLine thirdline(thirdp0.plane, thirdp0.w, thirdp0.t, thirdp1.w, thirdp1.t);
    
      this->seedlines.push_back(thirdline);
    
      mf::LogVerbatim("TWQProjectionView") <<"Third projection goes " 
					   << thirdline.w0<<" " 
					   <<thirdline.t0<< " to " 
					   << thirdline.w1 <<  " " 
					   << thirdline.t1;


      // Update all the views, since we are also adding
      //  guide lines
      for(size_t p=0; p!=fPlanes.size(); ++p){
	this->fPlanes[p]->DrawLinesinView(seedlines);
      }
      

      // Now make an actual seed object to store
      //  in the hit selector
      double seedpt[3], seeddir[3], seederr[3];
      
      for(int i=0; i!=3; ++i){
	seederr[i] = 0;
      }
      
    
      seedpt[0]  = (x0+x1)/2.;
      seeddir[0] = (x1-x0)/2.;

      seedpt[1]  = (yz0[0]+yz1[0])/2.;
      seeddir[1] = (yz1[0]-yz0[0])/2.;
      seedpt[2]  = (yz0[1]+yz1[1])/2.;
      seeddir[2] = (yz1[1]-yz0[1])/2.;
	  
      recob::Seed NewSeed(seedpt,seeddir,seederr,seederr); 
      fPlanes[0]->HitSelectorGet()->SeedVector().push_back(NewSeed);
      mf::LogVerbatim("TWQProjectionView") <<"Adding seed, vector now of size " 
					   << fPlanes[0]->HitSelectorGet()->SeedVector().size();
    

    }  
  
    // Update the bezier curve on the screen, and the track params
    //  text box
    mf::LogVerbatim("TWQProjectionView") << "size of seed vector: " 
					 << fPlanes[0]->HitSelectorGet()->SeedVector().size();
    
    int NSeeds = fPlanes[0]->HitSelectorGet()->SeedVector().size();

    double Length=0;
    std::stringstream ss;
    if(NSeeds>1)
      Length = UpdateSeedCurve();
    else if(NSeeds==1)
      Length = fPlanes[0]->HitSelectorGet()->SeedVector().at(0).GetLength()*2.; 
    
    if(Length>0){
      ss<<" " << Length<< " cm";   
      TGText * tt = new TGText(ss.str().c_str());
      tt->InsLine(1,"3D Track : ");
      fAngleInfo->SetText(tt);
      fAngleInfo->Update();
      ss.clear();
      ss.str("");
      
      double FirstPt[3], LastPt[3], FirstDir[3], LastDir[3], Err[3];
      fPlanes[0]->HitSelectorGet()->SeedVector().at(0).GetPoint(FirstPt,Err);
      fPlanes[0]->HitSelectorGet()->SeedVector().at(0).GetDirection(FirstDir,Err);
      fPlanes[0]->HitSelectorGet()->SeedVector().at(NSeeds-1).GetPoint(LastPt,Err);
      fPlanes[0]->HitSelectorGet()->SeedVector().at(NSeeds-1).GetDirection(LastDir,Err);
      for(int i=0; i!=3; ++i){
	FirstPt[i] -= FirstDir[i];
	LastPt[i]  += LastDir[i];
      }
      ss<<"z| " <<std::setprecision(0)<<std::fixed<< FirstPt[2] <<"," <<"\t" <<LastPt[2];
      tt = new TGText(ss.str().c_str());
      ss.clear();
      ss.str("");
      
      ss<<"x| " <<std::setprecision(0)<<std::fixed<< FirstPt[0] <<","<<"\t" <<LastPt[0];
      
      
      tt->InsLine(1,ss.str().c_str());    
      ss.clear();
      ss.str("");    
      
      ss<<"y| " <<std::setprecision(0)<<std::fixed<< FirstPt[1] <<","<<"\t" <<LastPt[1];
      
      tt->InsLine(1,ss.str().c_str());    
      ss.clear();
      ss.str("");    
      
      fXYZPosition->SetText(tt);
      
      art::ServiceHandle<geo::Geometry> geo;
      
      TVector3 TPCHighBoundary = TVector3(2 * geo->DetHalfWidth(),
					  geo->DetHalfHeight(),
					  geo->DetLength());
      TVector3 TPCLowBoundary = TVector3(0,
					 0-geo->DetHalfHeight(),
					 0);
      mf::LogVerbatim("TWQProjectionView") << "Boundaries " 
					   << TPCLowBoundary[0]  << " " 
					   << TPCLowBoundary[1]  << " " 
					   << TPCLowBoundary[2]  << "\n"
					   << TPCHighBoundary[0] << " " 
					   << TPCHighBoundary[1] << " " 
					   << TPCHighBoundary[2];
      bool OverPointFound=false;
      for(int i=0; i!=3; ++i){
	if((FirstPt[i]>(TPCHighBoundary[i])) ||
	   (FirstPt[i]<(TPCLowBoundary[i]))  ||
	   (LastPt[i]>(TPCHighBoundary[i]))  ||
	   (LastPt[i]<(TPCLowBoundary[i]))){
	  mf::LogVerbatim("TWQProjectionView") <<"One point over TPC boundary";
	  fXYZPosition->SetForegroundColor(0xff0000);
	  OverPointFound=true;
	}
	if(!OverPointFound){
	  if((FirstPt[i]>(TPCHighBoundary[i]-5))||(FirstPt[i]<(TPCLowBoundary[i]+5)) ||
	     (LastPt[i]>(TPCHighBoundary[i]-5))||(LastPt[i]<(TPCLowBoundary[i]+5))){
	    mf::LogVerbatim("TWQProjectionView") <<"One point 5cm from TPC boundary";
	    fXYZPosition->SetForegroundColor(0x0000ff);
	    OverPointFound=true;
	  }
	      
	}
      }
    
      fXYZPosition->Update();
  
    }

    return;
  }

  //......................................................................
  void TWQProjectionView::SelectHit(int plane)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if(!recoopt->fUseHitSelector) return;

    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();

    ///\todo What the heck is this all about???
    if(event!=7) return;

    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    if(evdlayoutopt->fMakeClusters != 1)
      return;
    //struct planepoint;
    int px = gPad->GetEventX();
    double w0 = gPad->AbsPixeltoX(px);
    double x = gPad->PadtoX(w0);

    int py = gPad->GetEventY();
    double t0 = gPad->AbsPixeltoY(py);
    double y = gPad->PadtoY(t0);

    fPlanes[plane]->SelectOneHit(x,y,zoom_opt);
    fPlanes[plane]->UpdatePad();
    return;

  }

  //......................................................................
  // if flag is true then zoom. If flag is false then unzoom.
  void 	TWQProjectionView::ZoomInterest(bool flag)
  {
    mf::LogVerbatim("TWQProjectionView") <<"ZoomInterest called";
  
    if(flag==true) zoom_opt="1";
    else zoom_opt="0";
  
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt; 
 
    ZoomOptions zo;
    //  mf::LogVerbatim("TWQProjectionView") <<"Zoom interest pushing back zoom options"<<std::endl;
    fPrevZoomOpt.push_back(fZoomOpt);

 
    for(size_t iplane = 0; iplane < fPlanes.size(); ++iplane){
      int minw,maxw,mint,maxt;
      if(flag){
	int test=0;
	if(rawopt->fDrawRawDataOrCalibWires == 0)
	  test = fPlanes[iplane]->RawDataDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
	else
	  fPlanes[iplane]->RecoBaseDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
	  
	if(test != 0)
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
  void TWQProjectionView::SetUpSideBar()
  {  
    SetUpZoomButtons();
    SetUpPositionFind();
    SetUpClusterButtons();
    SetUpDrawingButtons();
    SetUpTPCselection();
  }

  //......................................................................
  void TWQProjectionView::SetZoomInterest()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fAutoZoomInterest = fToggleAutoZoom->GetState();
    SetAutomaticZoomMode(evdlayoutopt->fAutoZoomInterest == 1);
  }

  //......................................................................
  void TWQProjectionView::SetZoomFromView() {
    for (TWireProjPad* pPlane: fPlanes)
      pPlane->SetZoomFromView();
  }
  

  //......................................................................
  void TWQProjectionView::SetClusterInterest()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fMakeClusters = fToggleClusters->GetState();
  }


  //......................................................................
  void TWQProjectionView::SetSeedInterest()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fMakeSeeds = fToggleSeeds->GetState();
    if(evdlayoutopt->fMakeSeeds==true && evdlayoutopt->fMakeClusters==false ){
      evdlayoutopt->fMakeClusters = fToggleSeeds->GetState();
      fToggleClusters->SetState(kButtonDown);
    }
  }

  //......................................................................   
  void TWQProjectionView::ToggleEndPointMarkers()
  {  
    art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
    evdlayoutopt->fShowEndPointMarkers= fToggleShowMarkers->GetState();
  }
  
  //......................................................................
  void TWQProjectionView::ForceRedraw() {
    LOG_DEBUG("TWQProjectionView") << "Explicit request for redrawing";
    
    // for now, bother only about redrawing the plane pads
    SetZoomFromView();
    DrawPads();
    
  } // TWQProjectionView::ForceRedraw()
  
  //......................................................................
  void TWQProjectionView::SetUpZoomButtons()
  {
    // enter zoom buttons
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  

    fZoomInterest=new TGTextButton(fVFrame,"&Zoom Interest",150);
    fZoomInterest->Connect("Clicked()", "evd::TWQProjectionView", this, "ZoomInterest()");



    fUnZoomInterest=new TGTextButton(fVFrame,"&UnZoom Interest",150);
    fUnZoomInterest->Connect("Clicked()", "evd::TWQProjectionView", this, "ZoomInterest(=false)");


    fZoomBack=new TGTextButton(fVFrame,"&Zoom Back",150);
    fZoomBack->Connect("Clicked()", "evd::TWQProjectionView", this, "ZoomBack()");


    fToggleAutoZoom=new TGCheckButton(fVFrame,"AutoZoom",0);;       ///< Toggle the autozoom setting 
    fToggleAutoZoom->Connect("Clicked()", "evd::TWQProjectionView", this, "SetZoomInterest()");
    if(evdlayoutopt->fAutoZoomInterest == 1) fToggleAutoZoom->SetState(kButtonDown);

    fVFrame->AddFrame(fZoomInterest,     new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fUnZoomInterest,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );

    fVFrame->AddFrame(fZoomBack,         new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );

    fVFrame->AddFrame(fToggleAutoZoom,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
  }

  //......................................................................
  void TWQProjectionView::SetUpClusterButtons()
  {
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  
    if(!evdlayoutopt->fShowClusterSection)            
      return;
    // enter zoom buttons
		
    fToggleZoom = new TGRadioButton(fVFrame,"Use Zoom",          2);
    fToggleClusters    = new TGRadioButton(fVFrame,"Select Clusters", 3);
    fToggleSeeds      = new TGRadioButton(fVFrame,"Select Seeds", 4);	
	
    fToggleZoom->Connect("Clicked()", "evd::TWQProjectionView", this, "RadioButtonsDispatch(=0)");
    fToggleClusters->Connect("Clicked()", "evd::TWQProjectionView", this, "RadioButtonsDispatch(=1)");
    fToggleSeeds->Connect("Clicked()", "evd::TWQProjectionView", this, "RadioButtonsDispatch(=2)");
	
	
	
    //fToggleClusters=new TGCheckButton(fVFrame,"&Make Clusters",0);      ///< Toggle the make cluster setting 
	
    fCalcAngle=new TGTextButton(fVFrame, "&Save Selection", 150);
    fCalcAngle->Connect("Clicked()", "evd::TWQProjectionView",this,"SaveSelection()");
	
    fClear=new TGTextButton(fVFrame, "&Clear Selection",0);
    fClear->Connect("Clicked()","evd::TWQProjectionView", this, "ClearSelection()");
	
    fClearLastSeed=new TGTextButton(fVFrame, "&Clear Last Seed",0);
    fClearLastSeed->Connect("Clicked()","evd::TWQProjectionView", this, "ClearLastSeed()");

    //	fRefitSeeds=new TGTextButton(fVFrame, "&Refit Seeds",0);
    //	fRefitSeeds->Connect("Clicked()","evd::TWQProjectionView", this, "RefitSeeds()");

    //fClearSeeds=new TGTextButton(fVFrame, "&Clear All Seeds",0);
    //fClearSeeds->Connect("Clicked()","evd::TWQProjectionView", this, "ClearAllSeeds()");
	
    if(evdlayoutopt->fMakeClusters == 1) fToggleClusters->SetState(kButtonDown);  
    else if(evdlayoutopt->fMakeSeeds == 1) fToggleSeeds->SetState(kButtonDown); 
    else fToggleZoom->SetState(kButtonDown); 
	
	
	
    fAngleInfo=new TGTextView(fVFrame,115,75,999,TGView::kNoHSB | TGView::kNoVSB); ///< Display the calculated angles
    fAngleInfo->SetEditable("false");
    TGText *tt=new TGText("...");
    fAngleInfo->SetText(tt);
	

    fDistance = new TGNumberEntry(fVFrame,0,6,-1,
				  TGNumberFormat::kNESReal, 
				  TGNumberFormat::kNEAPositive, 
				  TGNumberFormat::kNELLimitMinMax, 
				  0 , 100);
    // Initial value
    fDistance->SetNumber( 1.5 );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fDistance->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetDistance()");
    fDistance->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetDistance()");

    // Text label for this numeric field.
    fDistanceLabel= new TGLabel(fVFrame,"Distance");

	
    fVFrame->AddFrame(fToggleZoom, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
    fVFrame->AddFrame(fToggleClusters, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
    fVFrame->AddFrame(fToggleSeeds, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
	
    fVFrame->AddFrame(fCalcAngle,new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
    fVFrame->AddFrame(fClear,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,0,5,1));
    fVFrame->AddFrame(fClearLastSeed, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
    //	fVFrame->AddFrame(fRefitSeeds, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
	
    fVFrame->AddFrame(fDistance,     new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fDistanceLabel,     new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
	
    fVFrame->AddFrame(fAngleInfo,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,0,5,1));
	
	fVFrame->AddFrame(fDistance,     new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
	
	
    //fVFrame->AddFrame(fClearSeeds, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5,1));
	
  }
  
  //......................................................................
  void TWQProjectionView::SetUpDrawingButtons()
  {
    // enter zoom buttons
  //  art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  

    fRedraw = new TGTextButton(fVFrame, "&Redraw", 120);
    fRedraw->Connect("Clicked()", "evd::TWQProjectionView", this, "ForceRedraw()");

    fVFrame->AddFrame(fRedraw,           new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    
  } // SetUpDrawingButtons()


  //......................................................................
  std::string TWQProjectionView::TotalElementsString(unsigned int NElements)
    { return "(" + std::to_string(NElements) + " total)"; }
  
  void TWQProjectionView::SetUpTPCselection()
  {
    geo::GeometryCore const& geom = *(art::ServiceHandle<geo::Geometry>());
    art::ServiceHandle<evd::RawDrawingOptions> rawOpt;
    
    TGHorizontalFrame* pRow = nullptr;
    //
    // Cryostat selection line
    //
    // this is the subframe with horizontal alignment where we place our widgets:
    pRow = new TGHorizontalFrame(fVFrame, 216, 32, kHorizontalFrame);
    
    geo::CryostatID::CryostatID_t const CurrentCryo = rawOpt->fCryostat;
    unsigned int const NCryo = geom.Ncryostats();
    if (NCryo > 1) { // allow a selector
      unsigned int const NCryoDigits = std::to_string(NCryo - 1).length(); // a silly way fast to code...

      geo::CryostatID::CryostatID_t const CurrentCryo = rawOpt->fCryostat;
      
      // label
      TGLabel* pLabel = new TGLabel(pRow, "Cryo #");
      pLabel->SetTextJustify(kTextRight | kTextCenterY);
      
      // numerical input
      fCryoInput = new TGNumberEntry(
        pRow, (Double_t) CurrentCryo, // parent widget; starting value;
        NCryoDigits, -1, // number of digits for the input field; ID;
        TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, // type of number; number attributes;
        TGNumberFormat::kNELLimitMinMax, -1, NCryo // limits
        );

      TGLabel* pTotalCryoLabel = new TGLabel(pRow, TotalElementsString(NCryo).c_str());
      pTotalCryoLabel->SetTextJustify(kTextLeft | kTextCenterY);
      // the numbers are padding on the four sides
      pRow->AddFrame(pLabel,          new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
      pRow->AddFrame(fCryoInput,      new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
      pRow->AddFrame(pTotalCryoLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
     
      fCryoInput->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SelectTPC()");
    }
    else { // just place a static label
      TGLabel* pLabel = new TGLabel(pRow, "Cryo #0 (1 total)");
      // the numbers are padding on the four sides
      pRow->AddFrame(pLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
    }
    
    fVFrame->AddFrame(pRow, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    
    //
    // TPC selection line
    //
    // this is the subframe with horizontal alignment where we place our widgets:
    pRow = new TGHorizontalFrame(fVFrame, 216, 32, kHorizontalFrame);
    
    unsigned int MaxTPC = geom.MaxTPCs();
    if (MaxTPC > 1) { // label, numeric input, then total
      unsigned int const NTPCDigits = std::to_string(MaxTPC - 1).length(); // a silly way fast to code...

      geo::TPCID::TPCID_t const CurrentTPC = rawOpt->fTPC;
      unsigned int const NTPCs = geom.NTPC(geo::CryostatID(CurrentCryo));
      
      // label
      TGLabel* pLabel = new TGLabel(pRow, "TPC  #");
      pLabel->SetTextJustify(kTextRight | kTextCenterY);
      
      // numerical input
      fTPCInput = new TGNumberEntry(
        pRow, (Double_t) CurrentTPC, // parent widget; starting value;
        NTPCDigits, -1, // number of digits for the input field; ID;
        TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, // type of number; number attributes;
        TGNumberFormat::kNELLimitMinMax, -1, MaxTPC // limits
        );

      fTotalTPCLabel = new TGLabel(pRow, TotalElementsString(NTPCs).c_str());
      fTotalTPCLabel->SetTextJustify(kTextRight | kTextCenterY);
      // the numbers are padding on the four sides
      pRow->AddFrame(pLabel,         new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
      pRow->AddFrame(fTPCInput,      new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
      pRow->AddFrame(fTotalTPCLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
     
      fTPCInput->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SelectTPC()");
    }
    else { // just place another static label
      TGLabel* pLabel = new TGLabel(pRow, "TPC  #0 (1 total)");
      // the numbers are padding on the four sides
      pRow->AddFrame(pLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 5, 5));
    }
    
    fVFrame->AddFrame(pRow, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    
  } // TWQProjectionView::SetUpTPCselection()
  
  //----------------------------------------------------------------------------
  void	TWQProjectionView::SelectTPC()
  {
    /*
     * This function takes care of the input in the cryostat and TPC fields.
     * It is called whenever any of the two values (cryostat and TPC) change;
     * it can perform a number of tasks:
     * - first checks if the new input is valid
     * - if it's not, tries to wrap it to a valid TPC ID
     * - if it's not possible, goes back to the current TPC
     * - if the resulting TPC ID is different from the original one,
     *   updates the content of the current TPC and cryostat in the service,
     *   and asks for redrawing
     * - if the cryostat is changed, updates the TPC count
     * - if it changed the values of the TPC ID respect to what was in input
     *   at the beginning, then updates the input fields
     *
     */
    evd::RawDrawingOptions& rawOpt = *(art::ServiceHandle<evd::RawDrawingOptions>());
    geo::GeometryCore const& geom = *(art::ServiceHandle<geo::Geometry>());
    
    geo::TPCID CurrentTPC(rawOpt.fCryostat, rawOpt.fTPC);
    geo::TPCID RequestedTPC(
      fCryoInput? (unsigned int) fCryoInput->GetIntNumber(): 0U,
      fTPCInput? (unsigned int) fTPCInput->GetIntNumber(): 0U
      );
    geo::TPCID NewTPC(RequestedTPC);
    
    // if the input ends up being invalid, try to fix it somehow;
    // we give a special meaning to negative values;
    // since they are not supported by the container we store them in
    // (that is, TPCID) we have to handle negative values as special case:
    if (fCryoInput && (fCryoInput->GetIntNumber() < 0)) {
      // wrap back to the last cryostat, last TPC
      NewTPC.Cryostat = (geo::CryostatID::CryostatID_t) (geom.Ncryostats() - 1);
      NewTPC.TPC = (geo::TPCID::TPCID_t) (geom.NTPC(NewTPC) - 1);
    }
    else if (fTPCInput && (fTPCInput->GetIntNumber() < 0)) {
      // wrap back to the previous cryostat, last TPC
      if (NewTPC.Cryostat-- == 0)
        NewTPC.Cryostat = (geo::CryostatID::CryostatID_t) (geom.Ncryostats() - 1);
      NewTPC.TPC = (geo::TPCID::TPCID_t) (geom.NTPC(NewTPC) - 1);
    }
    else if (!geom.HasTPC(NewTPC)) {
      // both elements are positive: it must be overflow
      unsigned int const NCryos = geom.Ncryostats();
      // what's wrong?
      if (NewTPC.Cryostat >= NCryos) { // cryostat wrap
        NewTPC.Cryostat = (geo::CryostatID::CryostatID_t) 0;
        NewTPC.TPC = (geo::TPCID::TPCID_t) 0;
      }
      else { // TPC wrap
        if (++NewTPC.Cryostat >= NCryos) 
          NewTPC.Cryostat = (geo::CryostatID::CryostatID_t) 0;
        NewTPC.TPC = (geo::TPCID::TPCID_t) 0;
      }
      
      LOG_DEBUG("TWQProjectionView") << __func__ << ": invalid TPC "
        << RequestedTPC << ", corrected as " << NewTPC << " instead";
    }
    
    if (!geom.HasTPC(NewTPC)) { // weird...
      LOG_ERROR("TWQProjectionView") << __func__ << ": internal error: "
        << RequestedTPC << " turned into an invalid TPC " << NewTPC;
    }
    else if (NewTPC != CurrentTPC) { // do we need to change after all? 
      LOG_DEBUG("TWQProjectionView") << __func__ << ": switching from "
        << CurrentTPC << " to " << NewTPC;
    
      // new cryostat?
      if (rawOpt.fCryostat != NewTPC.Cryostat) { // update the TPC count
        unsigned int const NTPCs = geom.NTPC(NewTPC);
        fTotalTPCLabel->SetText(TotalElementsString(NTPCs).c_str());
      //  fTotalTPCLabel->Modified();
      }
      // update the current TPC in the service
      // (that is the thing everything else refers to)
      rawOpt.fCryostat = NewTPC.Cryostat;
      rawOpt.fTPC = NewTPC.TPC;
      
      // redraw the content
      ResetRegionsOfInterest();
      DrawPads();
    //  evdb::Canvas::fCanvas->cd();
    //  evdb::Canvas::fCanvas->Modified();
    //  evdb::Canvas::fCanvas->Update();
    }
    
    // if we have changed the requested TPC, we need to update the input fields
    if (NewTPC != RequestedTPC) {
      if (fCryoInput) fCryoInput->SetIntNumber(NewTPC.Cryostat);
      if (fTPCInput) fTPCInput->SetIntNumber(NewTPC.TPC);
    }
    
  } // TWQProjectionView::SelectTPC()
  
  
  //----------------------------------------------------------------------------
  void	TWQProjectionView::RadioButtonsDispatch(int parameter)
  {
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt; 
    if(parameter==0){evdlayoutopt->fMakeClusters=0; 
      evdlayoutopt->fMakeSeeds=0; 
      fToggleClusters->SetState(kButtonUp); 	 
      fToggleSeeds->SetState(kButtonUp);  
    }
    else if(parameter==1){
      evdlayoutopt->fMakeClusters=1; 
      evdlayoutopt->fMakeSeeds=0; 
      fToggleZoom->SetState(kButtonUp); 	 
      fToggleSeeds->SetState(kButtonUp);  
    }
    else if(parameter==2){
      evdlayoutopt->fMakeClusters=0; 
      evdlayoutopt->fMakeSeeds=1; 
      fToggleZoom->SetState(kButtonUp); 	 
      fToggleClusters->SetState(kButtonUp);  
    }
 
  }

  //......................................................................
  void TWQProjectionView::SetUpPositionFind()
  {
    // enter zoom buttons
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  
    if(!evdlayoutopt->fShowEndPointSection)            
      return;

    // int    	 fShowEndPointMarkers;             ///< Draw EndPoint Markers if clicked. 

    fFindEndpoint=new TGTextButton(fVFrame,"&Find XYZ",150);
    fFindEndpoint->Connect("Clicked()", "evd::TWQProjectionView", this, "FindEndPoint()");

    fXYZPosition=new TGTextView(fVFrame,100,55,999,TGView::kNoHSB | TGView::kNoVSB); ///< Display the xyz position
    fXYZPosition->SetEditable("false");
    TGText *tt=new TGText("x,y,z");
    fXYZPosition->SetText(tt);


    fClearPPoints=new TGTextButton(fVFrame,"&Clear Points",150);
    fClearPPoints->Connect("Clicked()", "evd::TWQProjectionView", this, "ClearEndPoints()");  // ?

    fToggleShowMarkers=new TGCheckButton(fVFrame,"ShowMarkers",0);       ///< Toggle the ShowEndPointMarkers Setting 
    fToggleShowMarkers->Connect("Clicked()", "evd::TWQProjectionView", this, "ToggleEndPointMarkers()");
    if(evdlayoutopt->fShowEndPointMarkers == 1) fToggleShowMarkers->SetState(kButtonDown);

    fVFrame->AddFrame(fFindEndpoint,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fXYZPosition,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fClearPPoints,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fToggleShowMarkers,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
  }


  /////////////////////////////////////////
  //  Go back one step in zoom

  void    TWQProjectionView::ZoomBack()
  {
    if(fPrevZoomOpt.size()>0){
      ZoomOptions ThePrevZoomOpt = fPrevZoomOpt.at(fPrevZoomOpt.size()-1);
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
      mf::LogVerbatim("TWQProjectionView") <<"unable to unzoom further - no zoom settings left on stack"<<std::endl;
  }

  //------------------------------------
  void    TWQProjectionView::SetZoom(int plane,
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
  void TWQProjectionView::SetPlaneWire()
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
  void TWQProjectionView::SetPlane()
  {
    kPlane = (unsigned int)fPlaneEntry->GetNumberEntry()->GetNumber();

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetWire()
  {
    art::ServiceHandle<geo::Geometry> geo;
    kWire = (geo->Nwires(kPlane)-1 > (unsigned int)fWireEntry->GetNumberEntry()->GetNumber() ) ?  (unsigned int)fWireEntry->GetNumberEntry()->GetNumber() : geo->Nwires(kPlane)-1;

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetDistance()
  {
    kDistance = (double)fDistance->GetNumberEntry()->GetNumber();
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetThreshold()
  {
    double threshold = fThresEntry->GetNumberEntry()->GetNumber();
      
    if (threshold != fLastThreshold)
    {
        art::ServiceHandle<evd::RawDrawingOptions> rawopt;
        rawopt->fMinSignal = threshold;

        TVirtualPad *ori = gPad;
        this->DrawPads(zoom_opt);
        evdb::Canvas::fCanvas->cd();
        evdb::Canvas::fCanvas->Modified();
        evdb::Canvas::fCanvas->Update();

        ori->cd();
    }
      
    fLastThreshold = threshold;

    return;
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetGreyscale()
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
  void TWQProjectionView::SetRawCalib()
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
  void TWQProjectionView::SetMCInfo()
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

  //-----------------------------------------------------------------
  double TWQProjectionView::UpdateSeedCurve()
  {
    double BezLength = util::kBogusD;

    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if(!recoopt->fUseHitSelector) return BezLength;

    std::vector<recob::Seed> SeedVec=fPlanes[0]->HitSelectorGet()->SeedVector();
    for(size_t p = 0; p != fPlanes.size(); ++p){
      // This method draws the existing seed lines on every
      //   plane, as well as guide lines mid seeding.
      
      fPlanes[p]->DrawLinesinView(seedlines);

      //  The update method draws the curve onto each 2D pad
      //   and returns the track length.
      // note - the length is the same in every view

      BezLength = fPlanes[p]->UpdateSeedCurve(SeedVec, p);
    }
    return BezLength;
  }


  //-----------------------------------------------------------------
  bool TWQProjectionView::OnNewEvent() {
    
    // first check if it's new...
    art::Event const* pEvent = evdb::EventHolder::Instance()->GetEvent();
    if (!pEvent) {
      if (!fLastEvent->isValid()) return false; // no event before, nor now
      fLastEvent->clear();
      return true;
    }
    
    // do we have a new event?
    if (!fLastEvent->update
      ({*pEvent, art::ServiceHandle<evd::RawDrawingOptions>()->fRawDataLabel})
      )
      return false;
    
    LOG_DEBUG("TWQProjectionView") << "New event or product: " << *fLastEvent;
    
    art::ServiceHandle<evd::EvdLayoutOptions> drawopt;
    SetAutomaticZoomMode(drawopt->fAutoZoomInterest == 1);
    
    return true; // yes, a new event is here
  } // TWQProjectionView::OnNewEvent()

  //-----------------------------------------------------------------
  
}// namespace
