////////////////////////////////////////////////////////////////////////
///
/// \file    TWireProjPad.cxx
/// \brief   Drawing pad for X-Z or Y-Z projections of events
/// \author  messier@indiana.edu
/// \version $Id: TZProjPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include "lareventdisplay/EventDisplay/TWireProjPad.h"
#include "TPad.h"
#include "TH1F.h"
#include "TString.h"
#include "TMarker.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TClass.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TList.h"

#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lareventdisplay/EventDisplay/GeometryDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/Style.h"
#include "TFrame.h"


#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace {
  
  template <typename Stream>
  void DumpPad(Stream&& log, TVirtualPad* pPad) {
    if (!pPad) {
      log << "pad not available";
      return;
    }
    
    log << pPad->IsA()->GetName() << "[" << ((void*) pPad) << "](\"" << pPad->GetName() << "\")";
    TFrame const* pFrame = pPad->GetFrame();
    if (pFrame) {
      double const low_wire = pFrame->GetX1(), high_wire = pFrame->GetX2();
      double const low_tdc = pFrame->GetY1(), high_tdc = pFrame->GetY2();
      double const wire_pixels = pPad->XtoAbsPixel(high_wire) - pPad->XtoAbsPixel(low_wire);
      double const tdc_pixels = -(pPad->YtoAbsPixel(high_tdc) - pPad->YtoAbsPixel(low_tdc));
      log << " has frame spanning wires "
        << low_wire << "-" << high_wire << " and TDC " << low_tdc << "-" << high_tdc
        << " in a window " << wire_pixels << "x" << tdc_pixels << " pixel big";
    }
    else {
      log << " has no frame";
    }

  } // DumpPad()
  
  
  [[gnu::unused]] void DumpPadsInCanvas
    (TVirtualPad* pPad, std::string caller, std::string msg = "")
  {
    mf::LogDebug log(caller);
    if (!msg.empty()) log << msg << ": ";
    if (!pPad) {
      log << "pad not available";
      return;
    }
    
    DumpPad(log, pPad);
    
    TCanvas const* pCanvas = pPad->GetCanvas();
    log << "\nCanvas is: (TCanvas*) (" << ((void*)pPad->GetCanvas())
       << ") with " << pCanvas->GetListOfPrimitives()->GetSize() << " primitives and the following pads:";
    TIterator* pIter = pCanvas->GetListOfPrimitives()->MakeIterator();
    TObject const* pObject;
    while ((pObject = pIter->Next())) {
      if (!pObject->InheritsFrom(TVirtualPad::Class())) continue;
      log << "\n  " << ((pObject == pPad)? '*': '-') << "  ";
      DumpPad(log, (TVirtualPad*) pObject);
    }
    log << "\n";
    delete pIter;
  } // DumpPadsInCanvas()

} // local namespace

namespace evd{

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
			     double x1, double x2,
			     double y1, double y2,
			     unsigned int plane)
    : DrawingPad(nm, ti, x1, x2, y1, y2)
    , fPlane(plane)
  {
    
    fCurrentZoom.resize(4);
    
    art::ServiceHandle<geo::Geometry> geo;
    
    this->Pad()->SetBit(kCannotPick);
    // this->Pad()->SetBit(TPad::kCannotMove);
    this->Pad()->cd();
    
    this->Pad()->SetLeftMargin  (0.070);
    this->Pad()->SetRightMargin (0.010);
    
    // how many planes in the detector and 
    // which plane is this one?
    
    unsigned int planes = geo->Nplanes();
    this->Pad()->SetTopMargin   (0.005);
    this->Pad()->SetBottomMargin(0.110);
    
    // there has to be a better way of doing this that does
    // not have a case for each number of planes in a detector
    if(planes == 2 && fPlane > 0){
      this->Pad()->SetTopMargin   (0.110);
      this->Pad()->SetBottomMargin(0.005);
    }
    else if(planes > 2){
      if(fPlane == 1){
	this->Pad()->SetTopMargin   (0.055);
	this->Pad()->SetBottomMargin(0.055);
      }
      else if(fPlane == 2){
	this->Pad()->SetTopMargin   (0.110);
	this->Pad()->SetBottomMargin(0.005);
      }
    }
    
    TString planeNo = "fTWirePlane";
    planeNo += fPlane;
    
    
    TString xtitle = ";Induction Wire;t (tdc)";
    if(geo->Plane(fPlane).SignalType() == geo::kCollection) xtitle = ";Collection Wire;t (tdc)";
    
    unsigned int const nWires = geo->Nwires(fPlane);
    unsigned int const nTicks = RawDataDraw()->TotalClockTicks();
    
    fXLo = -0.005 * (nWires - 1);
    fXHi =  1.005 * (nWires - 1);
    fYLo =  0.990*(unsigned int)(this->RawDataDraw()->StartTick());
    fYHi =  1.005*std::min((unsigned int)(this->RawDataDraw()->StartTick()+nTicks), nTicks);
    
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    fOri = rawopt->fAxisOrientation;
    if(fOri > 0){
      fYLo = -0.005 * (nWires - 1);
      fYHi =  1.005 * (nWires - 1);
      fYLo =  0.990*(unsigned int)(this->RawDataDraw()->StartTick());
      fYHi =  1.005*std::min((unsigned int)(this->RawDataDraw()->StartTick()+nTicks), nTicks);
      fXLo = -0.005 * nTicks;
      fXHi =  1.010 * nTicks;
      xtitle = ";t (tdc);InductionWire";
      if(geo->Plane(fPlane).SignalType() == geo::kCollection) xtitle = ";t (tdc);Collection Wire";
    }      
    
    // make the range of the histogram be the biggest extent 
    // in both directions and then use SetRangeUser() to shrink it down
    // that will allow us to change the axes on the fly
    double min = std::min(fXLo, fYLo);
    double max = std::max(fXHi, fYHi);
    
    fHisto = new TH1F(*(fPad->DrawFrame(min, min, max, max)));
    
    fHisto->SetTitleOffset(0.5,"Y");
    fHisto->SetTitleOffset(0.75,"X");
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
    if (fHisto) { delete fHisto; fHisto = 0; }
    if (fView)  { delete fView;  fView  = 0; }
  }
  
  //......................................................................
  void TWireProjPad::Draw(const char* opt) 
  {
    // DumpPadsInCanvas(fPad, "TWireProjPad", "Draw()");
    LOG_DEBUG("TWireProjPad") << "Started to draw plane " << fPlane;
    
    ///\todo: Why is kSelectedColor hard coded?
    int kSelectedColor = 4;
    fView->Clear();
    
    // grab the singleton holding the art::Event
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

      this->SimulationDraw()->MCTruthVectors2D(*evt, fView, fPlane);
      
      // the 2D pads have too much detail to be rendered on screen;
      // to act smarter, RawDataDrawer needs to know the range being plotted
      this->RawDataDraw()->   ExtractRange    (fPad, &GetCurrentZoom());
      this->RawDataDraw()->   RawDigit2D
        (*evt, fView, fPlane, GetDrawOptions().bZoom2DdrawToRoI);
      
      this->RecoBaseDraw()->  Wire2D          (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Hit2D           (*evt, fView, fPlane);
      
      if(recoOpt->fUseHitSelector)
        this->RecoBaseDraw()->Hit2D(this->HitSelectorGet()->GetSelectedHits(fPlane),
                                    kSelectedColor,
                                    fView);
   
      this->RecoBaseDraw()->  Cluster2D             (*evt, fView, fPlane);
      this->RecoBaseDraw()->  EndPoint2D            (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Prong2D               (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Vertex2D              (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Seed2D                (*evt, fView, fPlane);
      this->RecoBaseDraw()->  BezierTrack2D         (*evt, fView, fPlane);
      this->RecoBaseDraw()->  OpFlash2D             (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Event2D               (*evt, fView, fPlane);
      this->RecoBaseDraw()->  DrawTrackVertexAssns2D(*evt, fView, fPlane);
      
    //  DumpPadsInCanvas(fPad, "TWireProjPad", "Before UpdatePad()");
      UpdatePad();
    } // if (evt)

    // DumpPadsInCanvas(fPad, "TWireProjPad", "Before ClearandUpdatePad()");
    ClearandUpdatePad();
    
    // DumpPadsInCanvas(fPad, "TWireProjPad", "After ClearandUpdatePad()");
       
    // check if we need to swap the axis ranges
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if(fOri != rawopt->fAxisOrientation){
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

    if (fPlane > 0) fHisto->Draw("X+");
    else            fHisto->Draw("");

      
    // Check if we should zoom the displays;
    // if there is no event, we have no clue about the region of interest
    // and therefore we don't touch anything
    if (opt==0 && evt) {
      //       if (drawopt->fAutoZoom) this->AutoZoom();
      //       else                    this->ShowFull();
      this->ShowFull();
    }
    
    LOG_DEBUG("TWireProjPad") << "Started rendering plane " << fPlane;
    
    fView->Draw();

    LOG_DEBUG("TWireProjPad") << "Drawing of plane " << fPlane << " completed";
  }

  //......................................................................
  void TWireProjPad::ClearHitList()
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    if(recoOpt->fUseHitSelector){
      this->HitSelectorGet()->ClearHitList(fPlane);
      this->Draw();
    }

    return;
  }

  //......................................................................

  ///
  /// Automatically zoom the view to a size just larger than the
  /// events. Also ensures that the aspect ratio is the same for the XZ
  /// and YZ projections.
  ///
  //   void TWireProjPad::AutoZoom()
  //   {
  //     double xmin, ymin, zmin;
  //     double xmax, ymax, zmax;
  //     this->RawDataDraw()->GetLimits(&xmin, &xmax, 
  // 				   &ymin, &ymax, 
  // 				   &zmin, &zmax);
  //     double dx = xmax-xmin;
  //     double dy = ymax-ymin;
  //     double dz = zmax-zmin;
  
  //     if (dx<dy) dx = dy;
  //     else       dy = dx;
  //     xmin = 0.5*(xmin+xmax)-0.6*dx;
  //     xmax = 0.5*(xmin+xmax)+0.6*dx;
  //     ymin = 0.5*(ymin+ymax)-0.6*dy;
  //     ymax = 0.5*(ymin+ymax)+0.6*dy;
  //     zmin -= 0.1*dz;
  //     zmax += 0.1*dz;
  
  //     fHisto->GetXaxis()->SetRangeUser(zmin,zmax);
  //     if (fXorY==kX) fHisto->GetYaxis()->SetRangeUser(xmin,xmax);
  //     else           fHisto->GetYaxis()->SetRangeUser(ymin,ymax);
  //   }

  //......................................................................
  // the override parameter is needed to unzoom to full range when the fAutoZoomInterest is on. 

  void TWireProjPad::ShowFull(int override)
  {
    art::ServiceHandle<geo::Geometry> g;
    
    // x values are wire numbers, y values are ticks of the clock
    int xmin = fXLo;
    int xmax = fXHi;
    int ymax = fYHi;
    int ymin = fYLo;
    
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    
    if(GetDrawOptions().bZoom2DdrawToRoI && !override){
      int test=0;
      if(rawopt->fDrawRawDataOrCalibWires == 0)
	test=RawDataDraw()->GetRegionOfInterest((int)fPlane,xmin,xmax,ymin,ymax);
      else
	test=RecoBaseDraw()->GetRegionOfInterest((int)fPlane,xmin,xmax,ymin,ymax);
    
      if(test != 0) return;     
    }
	
    SetZoomRange(xmin, xmax, ymin, ymax);

    return;
  }

  //......................................................................
  void TWireProjPad::GetWireRange(int* i1, int* i2) const 
  {
    if(fOri < 1){
      *i1 = fHisto->GetXaxis()->GetFirst();
      *i2 = fHisto->GetXaxis()->GetLast();
    }
    else{
      *i1 = fHisto->GetYaxis()->GetFirst();
      *i2 = fHisto->GetYaxis()->GetLast();
    }
    
    return;
  }

  //......................................................................
  // Set the X axis range only
  //
  void TWireProjPad::SetWireRange(int i1, int i2)
  {
    if(fOri < 1){
      fHisto->GetXaxis()->SetRange(i1,i2);
    }
    else{
      fHisto->GetYaxis()->SetRange(i1,i2);
    }
    fCurrentZoom[0] = i1;
    fCurrentZoom[1] = i2;

    return;
  }

  //......................................................................
  // Set the visible range of the wire / time view
  //
  void TWireProjPad::SetZoomRange(int i1, int i2,int y1, int y2)
  {
    LOG_DEBUG("TWireProjPad")
      << "SetZoomRange(" << i1 << ", " << i2 << ", " << y1 << ", " << y2
      << ") on plane #" << fPlane;
    
    fHisto->GetXaxis()->SetRangeUser(i1,i2);
    fHisto->GetYaxis()->SetRangeUser(y1,y2);
    fCurrentZoom[0]=i1;
    fCurrentZoom[1]=i2;
    fCurrentZoom[2]=y1;
    fCurrentZoom[3]=y2;
  }
  //......................................................................
  // Set the visible range of the wire / time view from the view
  //
  void TWireProjPad::SetZoomFromView() {
    TAxis const& xaxis = *(fHisto->GetXaxis());
    fCurrentZoom[0] = xaxis.GetBinLowEdge(xaxis.GetFirst());
    fCurrentZoom[1] = xaxis.GetBinUpEdge(xaxis.GetLast());
    fCurrentZoom[2] = fHisto->GetMinimum();
    fCurrentZoom[3] = fHisto->GetMaximum();
    LOG_DEBUG("TWireProjPad") << "Zoom set to wires ("
      << fCurrentZoom[0] << "; " << fCurrentZoom[1] << " ), tick ("
      << fCurrentZoom[2] << "; " << fCurrentZoom[3] << ") for plane #"
      << fPlane;
  } // TWireProjPad::SetZoomFromView()
  //......................................................................
  void TWireProjPad::SaveHitList(double i1, 
				 double i2,
				 double y1, 
				 double y2, 
				 double distance, 
				 const char* zoom_opt,
				 bool good_plane)
  {  
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
      if(recoopt->fUseHitSelector){
	this->HitSelectorGet()->SaveHits(*evt, fView, fPlane, i1, i2, y1, y2, distance, good_plane);
	this->Draw(zoom_opt);
      }
    }

    return;
  }

  /////////////////////////////////////////////////
  // Pass the seed list onwards to InfoTransfer
  //
  double TWireProjPad::SaveSeedList(std::vector< util::PxLine > seedlines, 
				    double distance)
  {
    double KineticEnergy = util::kBogusD;
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
      if(recoopt->fUseHitSelector)
	KineticEnergy = this->HitSelectorGet()->SaveSeedLines(*evt, fView,seedlines, distance);
    }
    return KineticEnergy;
  }

  //......................................................................
  void TWireProjPad::SelectOneHit(double x, 
				  double y, 
				  const char* zoom_opt)
  {
  
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
      if(recoopt->fUseHitSelector){
	this->HitSelectorGet()->ChangeHit(*evt, fView, fPlane,x,y);
	this->Draw(zoom_opt);
      }
    }
    
    return;
  }
  
  //......................................................................
  void TWireProjPad::ClearandUpdatePad()
  {
    fPad->Clear(); 
    this->UpdatePad(); 

    return;
  }

  //......................................................................
  void TWireProjPad::UpdatePad(){
    fPad->cd();
    fPad->Modified();
    fPad->Update();
    fPad->GetFrame()->SetBit(TPad::kCannotMove,true);
    fPad->SetBit(TPad::kCannotMove,true);

    return;
  }
  
  //......................................................................
  void TWireProjPad::DrawLinesinView(std::vector< util::PxLine > lines,
				     bool deleting, 
				     const char * zoom_opt)
  {

    art::ServiceHandle<evd::EvdLayoutOptions>    evdlayoutopt;
    detinfo::DetectorProperties const* det = lar::providerFrom<detinfo::DetectorPropertiesService>();

    fPad->cd();
    if(deleting) {
      fPad->Clear();
      this->Draw(zoom_opt);
    }
    else {
      fView->Clear();
      fView->Draw();
    }
  
    mf::LogVerbatim("TWireProjPad") << "Drawing " << lines.size() <<" lines";
  
    for(size_t is = 0; is < lines.size(); ++is){
      if(fPlane!=lines[is].plane)
	continue;
      
      TLine& l = fView->AddLine(lines[is].w0,lines[is].t0,lines[is].w1,lines[is].t1);
      
      fView->Draw();
      evd::Style::FromPDG(l,11);
     
      // In Seed mode, colour of "sealed" seeds to green
      if(evdlayoutopt->fMakeSeeds){
	if( ( (lines.size()%3)==0 ) ||
	    ( is < ( lines.size()-(lines.size()%3) )  )) {
	  l.SetLineColor(kGreen);	      
	}	  
	else
	  l.SetLineColor(kRed);
      }
    } 
  
    // Seed mode guide lines
  
    if(evdlayoutopt->fMakeSeeds){
      TLine &lg1 = fView->AddLine(0,0,0,0);
      fView->Draw();
      TLine &lg0 = fView->AddLine(0,0,0,0);
      fView->Draw();
      lg1.SetLineStyle(kDashed);
      lg0.SetLineStyle(kDashed);      
      lg0.SetLineWidth(1);
      lg1.SetLineWidth(1);
      lg0.SetLineColor(kGray);
      lg1.SetLineColor(kGray);
      lg0.SetBit(kCannotPick);
      lg1.SetBit(kCannotPick);
	
      if((lines.size()%3)==1){
	mf::LogVerbatim("TWireProjPad") << "adding guide lines";
	util::PxLine TopLine = lines.at(lines.size()-1);
	lg0.SetX1(1);
	lg0.SetX2(5000);

	double TopT0 = det->ConvertXToTicks(det->ConvertTicksToX(TopLine.t0,TopLine.plane,0,0),fPlane,0,0); 
	double TopT1 = det->ConvertXToTicks(det->ConvertTicksToX(TopLine.t1,TopLine.plane,0,0),fPlane,0,0); 

	lg0.SetY1(TopT0);
	lg0.SetY2(TopT0);

	lg1.SetX1(1);
	lg1.SetX2(5000);
	  
	lg1.SetY1(TopT1);
	lg1.SetY2(TopT1);
	  
	
      }
    }
    fView->Draw();
    UpdatePad();
    fView->Draw();
  
    return;
  }

  ////////////////////////////////////////////
  //
  // This method creates and draws a curve on 
  // each view, given the seeds selected by
  // the user in the HitSelector.
  //
  // The return value is the length of the 3D
  // track
  //
  double TWireProjPad::UpdateSeedCurve(std::vector<recob::Seed> SeedVec, int plane)
  {
    mf::LogVerbatim("TWireProjPad") <<"running updateseedcurve for plane " << plane;
    fView->Draw();
    UpdatePad();
    int N=100;

    double ReturnVal=0;
  
    // For some reason, this line is needed to prevent lines being drawn twice
    //  TPolyLine& pldummy1 = fView->AddPolyLine(2,kBlue,1,0);
    // pldummy1.SetPoint(0,0,0);
    // pldummy1.SetPoint(1,0,0);
 
    int c=0; int t=0;
    int LastGoodValue=0;
    double ticks[3];
    double projpt[3]; 

    if(SeedVec.size() > 1){
      TPolyLine& pl = fView->AddPolyLine(N,kOrange+9,2,0);
      fView->Draw();
        
      trkf::BezierTrack BTrack(SeedVec);
      
      for(int i = 0; i != N; ++i){
	try{
	  BTrack.GetProjectedPointUVWT(float(i)/N,projpt,ticks,c,t );
	  LastGoodValue=i;
	  LOG_DEBUG("TWireProjPad") << i << " ";
	}
	catch(cet::exception excp){
	  BTrack.GetProjectedPointUVWT(float(LastGoodValue)/N, projpt, ticks, c, t);
	}
	  
	double x = projpt[plane];
	double y = ticks[plane];
	pl.SetPoint(i,x,y);
	  
	if(LastGoodValue!=i){
	  TMarker& mrk = fView->AddMarker(x, y, 3, 34, 1.5);
	  mrk.SetMarkerColor(3);
	}
      }
      ReturnVal =  BTrack.GetLength();
    }
    else
      ReturnVal=0;

    fView->Draw();
    UpdatePad();
  
    return ReturnVal;
  }
  
}// namespace
////////////////////////////////////////////////////////////////////////





