///
/// \file    TQPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: TQPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "lareventdisplay/EventDisplay/TQPad.h"

#include "TBox.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPad.h"
#include "TPolyLine.h"
#include "TText.h"

#include "cetlib/exception.h"

#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"

// C/C++ standard libraries
#include <algorithm> // std::min(), std::max()


namespace evd{

   static const int kRAW      = 0;
   static const int kCALIB    = 1;
   static const int kRAWCALIB = 2;
   static const int kQ        = 0;
   static const int kTQ       = 1;

   //......................................................................

   TQPad::TQPad(const char* nm, const char* ti,
                double x1, double y1,
                double x2, double y2,
                const char *opt,
                unsigned int plane,
                unsigned int wire) :
      DrawingPad(nm, ti, x1, y1, x2, y2),
      fWire(wire),
      fPlane(plane),
      fRawHisto(0),
      fRecoHisto(0)
   {

      art::ServiceHandle<geo::Geometry> geo;
      unsigned int planes = geo->Nplanes();

      this->Pad()->cd();

      this->Pad()->SetLeftMargin  (0.050);
      this->Pad()->SetRightMargin (0.050);

      this->Pad()->SetTopMargin   (0.005);
      this->Pad()->SetBottomMargin(0.110);

      // there has to be a better way of doing this that does
      // not have a case for each number of planes in a detector
      if(planes == 2 && fPlane > 0){
         this->Pad()->SetTopMargin   (0.110);
         this->Pad()->SetBottomMargin(0.010);
      }
      else if(planes > 2){
         if(fPlane == 1){
            this->Pad()->SetTopMargin   (0.005);
            this->Pad()->SetBottomMargin(0.010);
         }
         else if(fPlane == 2){
            this->Pad()->SetTopMargin   (0.110);
            this->Pad()->SetBottomMargin(0.010);
         }
      }


      std::string opts(opt);
      if(opts == "TQ") {
        fTQ = kTQ;
        // BB adjust the vertical spacing
        this->Pad()->SetTopMargin   (0);
        this->Pad()->SetBottomMargin(0.2);
      }
      if(opts == "Q" ){
         fTQ = kQ;
      }

      this->BookHistogram();
      fView = new evdb::View2D();
   }

   //......................................................................

   TQPad::~TQPad() 
   {
      if (fView)      { delete fView;      fView      = 0; }
      if (fRawHisto)  { delete fRawHisto;  fRawHisto  = 0; }
      if (fRecoHisto) { delete fRecoHisto; fRecoHisto = 0; }
   }

   //......................................................................
   void TQPad::Draw() 
   {
      art::ServiceHandle<evd::RawDrawingOptions> drawopt;

      if (!fRawHisto || !fRecoHisto) return;

      //grab the singleton with the event
      const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
      if(!evt) return;

      fPad->Clear();
      fPad->cd();

      std::vector<double> hstart;
      std::vector<double> hend;
      std::vector<double> hamplitudes;
      std::vector<double> hpeaktimes;

      if(fTQ == kTQ){
         fRawHisto->Reset("ICEM");
         fRecoHisto->Reset("ICEM");

         this->RawDataDraw()->FillTQHisto(*evt,
                                          fPlane,
                                          fWire, 
                                          fRawHisto);

         this->RecoBaseDraw()->FillTQHisto(*evt,
                                           fPlane,
                                           fWire, 
                                           fRecoHisto,
                                           hstart,
                                           hend,
                                           hamplitudes,
                                           hpeaktimes);
         
         // draw with histogram style, only (square) lines, no errors
         static const std::string defaultDrawOptions = "HIST";
         
         switch (drawopt->fDrawRawDataOrCalibWires) {
           case kRAW:
             fRawHisto->Draw(defaultDrawOptions.c_str());
             break;
           case kCALIB:
             fRecoHisto->Draw(defaultDrawOptions.c_str());
             break;
           case kRAWCALIB:
             fRawHisto->SetMaximum(1.1*std::max(fRawHisto->GetMaximum(), fRecoHisto->GetMaximum()));
             fRawHisto->SetMinimum(1.1*std::min(fRawHisto->GetMinimum(), fRecoHisto->GetMinimum()));
             fRawHisto->Draw(defaultDrawOptions.c_str());
             fRecoHisto->Draw((defaultDrawOptions + " same").c_str());
             break;
         } // switch

         // this loop draws Gaussian shapes for identified hits in the reco histo
         for (size_t i = 0; i < hstart.size() && drawopt->fDrawRawDataOrCalibWires != kRAW; ++i) {
            //hend and hstart are 1-sigma away from peak
            double width = (hend[i]-hstart[i])/2;

            //create a function corresponding to the Gaussian shape
            TF1 *f1 = new TF1("hitshape","gaus(0)",hstart[i]-1.5*width,hend[i]+1.5*width);//5-sigma wide window
            f1->SetParameters(hamplitudes[i],hpeaktimes[i],width);

            //create TPolyLine that actually gets drawn
            TPolyLine& p1 = fView->AddPolyLine(1001, 
                                               kOrange+7,
                                               3,
                                               1);

            //set coordinates of TPolyLine based on Gaussian function
            for(int j = 0; j<1001; ++j){ 
               double x = hstart[i]-1.5*width+j*5*width/1000;
               double y = f1->Eval(x); 
               p1.SetPoint(j, x, y);
            }
            p1.Draw("same");
            if(f1) delete f1;
         }
         
/* This code needs additional work to draw the text on the pad
        // BB: draw plane and wire number on the histogram
        std::string pw = "P:W = " + std::to_string(this->fPlane) +
          ":" + std::to_string(this->fWire);
        char const* txt = pw.c_str();
        std::cout<<" my text "<<txt<<"\n";
        double xp = 0.1, yp = 0.9;
        this->cd();
        TText& plnwir = fView->AddText(xp, yp, txt);xxx
        plnwir.SetTextColor(kBlack);
        plnwir.Draw("same");
*/
         if     (drawopt->fDrawRawDataOrCalibWires == kCALIB) fRecoHisto->Draw((defaultDrawOptions + " same").c_str());
         else if(drawopt->fDrawRawDataOrCalibWires == kRAWCALIB){
            fRawHisto->Draw((defaultDrawOptions + " same").c_str());
            fRecoHisto->Draw((defaultDrawOptions + " same").c_str());
         }

         fRawHisto->SetTitleOffset(0.2, "Y");
         fRecoHisto->SetLabelSize(0.2, "Y");

      } // end if fTQ == kTQ
    
      if(fTQ == kQ){

	// figure out the signal type for this plane, assume that
	// plane n in each TPC/cryostat has the same type
	/// \todo: May need to figure out cryostat and tpc number to be displayed
	art::ServiceHandle<geo::Geometry> geo;
	geo::SigType_t sigType = geo->Cryostat(0).TPC(0).Plane(fPlane).SignalType();

	art::ServiceHandle<evd::ColorDrawingOptions> cst;
	
	TH1F *hist;
	
	int ndiv = 0;
	if(drawopt->fDrawRawDataOrCalibWires != kCALIB){
	  hist = fRawHisto;
	  hist->SetMinimum(cst->fRawQLow [(size_t)sigType]);
	  hist->SetMaximum(cst->fRawQHigh[(size_t)sigType]);
	  ndiv = cst->fRawDiv[(size_t)sigType];
	}
	if(drawopt->fDrawRawDataOrCalibWires == kCALIB){
	  hist = fRecoHisto;
	  hist->SetMinimum(cst->fRecoQLow [(size_t)sigType]);
	  hist->SetMaximum(cst->fRecoQHigh[(size_t)sigType]);
	  ndiv = cst->fRecoDiv[(size_t)sigType];
	}
	
	hist->SetLabelSize(0, "X");
	hist->SetLabelSize(0, "Y");
	hist->SetTickLength(0, "X");
	hist->SetTickLength(0, "Y");
	hist->Draw("pY+");
	
	//
	// Use this to fill the histogram with colors from the color scale
	//
	double x1, x2, y1, y2;
	x1 = 0.;
	x2 = 1.;
	
	for(int i = 0; i < ndiv; ++i){
	  y1 = hist->GetMinimum() + i*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);
	  y2 = hist->GetMinimum() + (i + 1)*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);
	  
	  int c = 1;
	  if (drawopt->fDrawRawDataOrCalibWires==kRAW) {
	    c = cst->RawQ(sigType).GetColor(0.5*(y1+y2));
	  }
	  if (drawopt->fDrawRawDataOrCalibWires!=kRAW) {
	    c= cst->CalQ(sigType).GetColor(0.5*(y1+y2));
	  }
	  
	  TBox& b = fView->AddBox(x1,y1,x2,y2);
	  b.SetFillStyle(1001);
	  b.SetFillColor(c);      
	  b.Draw();
	} // end loop over Q histogram bins
        

	hist->Draw("same");
	
      } // end if fTQ == kQ

      return;
   }


   //......................................................................
   void TQPad::BookHistogram() 
   {

     if (fRawHisto) {
       delete fRawHisto; 
       fRawHisto = 0;
     }
     if (fRecoHisto) {
       delete fRecoHisto; 
       fRecoHisto = 0;
     }
     
     art::ServiceHandle<evd::ColorDrawingOptions> cst;
     art::ServiceHandle<evd::RawDrawingOptions> drawopt;

     // figure out the signal type for this plane, assume that
     // plane n in each TPC/cryostat has the same type
     /// \todo: May need to figure out cryostat and tpc number to be displayed
     art::ServiceHandle<geo::Geometry> geo;
     geo::SigType_t sigType = geo->Cryostat(0).TPC(0).Plane(fPlane).SignalType();
     
     /// \todo decide if ndivraw and ndivreco are useful
     //     int    ndivraw   = cst->fRawDiv;
     //     int    ndivreco  = cst->fRecoDiv;
     double qxloraw   = cst->fRawQLow[(size_t)sigType];
     double qxhiraw   = cst->fRawQHigh[(size_t)sigType];
     double qxloreco  = cst->fRecoQLow[(size_t)sigType];
     double qxhireco  = cst->fRecoQHigh[(size_t)sigType];
     double tqxlo     = 1.*this->RawDataDraw()->StartTick();
     double tqxhi     = 1.*this->RawDataDraw()->TotalClockTicks();
     
     switch (fTQ) {
       case kQ:
         fRawHisto = new TH1F("fRAWQHisto", ";;", 2,0.,1.);
         fRawHisto->SetMaximum(qxhiraw);
         fRawHisto->SetMinimum(qxloraw);
         
         fRecoHisto = new TH1F("fCALQHisto", ";;", 1,0.,1.);
         fRecoHisto->SetMaximum(qxhireco);
         fRecoHisto->SetMinimum(qxloreco);
         break; // kQ
       case kTQ:
         fRawHisto = new TH1F("fRAWTQHisto", ";t [ticks];q [ADC]", (int)tqxhi,tqxlo,tqxhi+tqxlo);
         fRecoHisto = new TH1F("fCALTQHisto", ";t [ticks];q [ADC]", (int)tqxhi,tqxlo,tqxhi+tqxlo);
         fRecoHisto->SetLineColor(kBlue);
         break;
       default:
         throw cet::exception("TQPad") << __func__ << ": unexpected quantity #" << fTQ << "\n";
     }//end if fTQ == kTQ
     
     fRawHisto->SetLabelSize  (0.15,"X");
     fRawHisto->SetLabelOffset(0.00,"X");
     fRawHisto->SetTitleSize  (0.15,"X");
     fRawHisto->SetTitleOffset(0.80,"X");
     
     fRawHisto->SetLabelSize  (0.10,"Y");
     fRawHisto->SetLabelOffset(0.01,"Y");
     fRawHisto->SetTitleSize  (0.15,"Y");
     fRawHisto->SetTitleOffset(0.80,"Y");
     
     fRecoHisto->SetLabelSize  (0.15,"X");
     fRecoHisto->SetLabelOffset(0.00,"X");
     fRecoHisto->SetTitleSize  (0.15,"X");
     fRecoHisto->SetTitleOffset(0.80,"X");
     
     fRecoHisto->SetLabelSize  (0.15,"Y");
     fRecoHisto->SetLabelOffset(0.00,"Y");
     fRecoHisto->SetTitleSize  (0.15,"Y");
     fRecoHisto->SetTitleOffset(0.80,"Y");
   }

}
//////////////////////////////////////////////////////////////////////////
