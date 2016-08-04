/// \file    AnalysisBaseDrawer.cxx
/// \brief   Class to aid in the rendering of AnalysisBase objects
/// \author  msoderbe@syr.edu
#include "TMarker.h"
#include "TBox.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TVector3.h"
#include "TLatex.h"

#include "nutools/EventDisplayBase/View2D.h"
#include "lareventdisplay/EventDisplay/eventdisplay.h"
#include "lareventdisplay/EventDisplay/Style.h"
#include "lareventdisplay/EventDisplay/AnalysisBaseDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/AnalysisDrawingOptions.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/RecoObjects/BezierTrack.h"
#include "lardata/AnalysisAlg/CalorimetryAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/FindMany.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <math.h>

/* unused function
namespace {
   // Utility function to make uniform error messages.
   void writeErrMsg(const char* fcn,
                    cet::exception const& e)
   {
      mf::LogWarning("AnalysisBaseDrawer") << "AnalysisBaseDrawer::" << fcn
                                           << " failed with message:\n"
                                           << e;
   }
}
*/

namespace evd{

   //......................................................................
   AnalysisBaseDrawer::AnalysisBaseDrawer() 
   {
    
   }

   //......................................................................
   AnalysisBaseDrawer::~AnalysisBaseDrawer() 
   {
 
   }

   //......................................................................
   void AnalysisBaseDrawer::DrawDeDx(const art::Event& evt,
                                  evdb::View2D* view)
   {
      art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
      art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;
      art::ServiceHandle<geo::Geometry>               geom;

      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
       
         //Get Track collection
         std::string which = recoOpt->fTrackLabels[imod];
         art::Handle<std::vector<recob::Track> > trackListHandle;
         evt.getByLabel(which,trackListHandle);
         std::vector<art::Ptr<recob::Track> > tracklist;
         art::fill_ptr_vector(tracklist, trackListHandle);
       
         //Loop over Calorimetry collections
         for(size_t cmod = 0; cmod < anaOpt->fCalorimetryLabels.size(); ++cmod) {
            std::string const callabel = anaOpt->fCalorimetryLabels[cmod];
            //Association between Tracks and Calorimetry
            art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, callabel);
	    if (!fmcal.isValid()) continue;
            //Loop over PID collections
            for(size_t pmod = 0; pmod < anaOpt->fParticleIDLabels.size(); ++pmod) {
               std::string const pidlabel = anaOpt->fParticleIDLabels[pmod];
               //Association between Tracks and PID
               art::FindMany<anab::ParticleID> fmpid(trackListHandle, evt, pidlabel);
	       if (!fmpid.isValid()) continue;
         
               //Loop over Tracks
               for(size_t trkIter = 0; trkIter<tracklist.size(); ++trkIter){
		 int color = tracklist[trkIter].key()%evd::kNCOLS;
		 std::vector<const anab::Calorimetry*> calos = fmcal.at(trkIter);
		 std::vector<const anab::ParticleID*> pids = fmpid.at(trkIter);
		 if (!calos.size()) continue;
		 if (calos.size()!=pids.size()) continue;
		 size_t bestplane = 0;
                 size_t calopl = 0;
                 size_t pidpl = 0;
		 size_t nmaxhits = 0;
		 for (size_t icalo = 0; icalo < calos.size(); ++icalo){
		   if (calos[icalo]->dEdx().size() > nmaxhits){
		     nmaxhits = calos[icalo]->dEdx().size();
		     bestplane = calos[icalo]->PlaneID().Plane;
		   }
		 }
		 if (anaOpt->fCaloPlane>=0 and anaOpt->fCaloPlane<int(geom->Nplanes())){
                   for (size_t icalo = 0; icalo < calos.size(); ++icalo){
                     if (int(calos[icalo]->PlaneID().Plane)==anaOpt->fCaloPlane&&
                         calos[icalo]->dEdx().size())
                       bestplane = calos[icalo]->PlaneID().Plane;
                   }
                 }

                 for (size_t icalo = 0; icalo < calos.size(); ++icalo){
                   if (calos[icalo]->PlaneID().Plane==bestplane){
                     calopl = icalo;
                   }
                 }
                 for (size_t ipid = 0; ipid < pids.size(); ++ipid){
                   if (pids[ipid]->PlaneID().Plane==bestplane){
                     pidpl = ipid;
                   }
                 }

		 TPolyMarker& pm = view->AddPolyMarker(calos[calopl]->dEdx().size(),evd::kColor[color],8,0.8);
		 for(size_t h = 0; h<calos[calopl]->dEdx().size();++h){
		   double xvalue = calos[calopl]->ResidualRange().at(h);
		   double yvalue = calos[calopl]->dEdx().at(h);
		   pm.SetPoint(h,xvalue,yvalue);
                   
		   double error = yvalue*(0.04231 + 0.0001783*(yvalue*yvalue));
		   TLine& l = view->AddLine(xvalue,yvalue-error,xvalue,yvalue+error);
		   l.SetLineColor(evd::kColor[color]);
		 }
		 
		 char trackinfo[80];
		 char pida[80];
		 char proton[80];
		 //char kaon[80];
		 char pion[80];
		 //char muon[80];
		 sprintf(trackinfo,"Track #%d: K.E. = %.1f MeV , Range = %.1f cm",
			 int(tracklist[trkIter].key()),
			 calos[calopl]->KineticEnergy(),
			 calos[calopl]->Range());
		 sprintf(proton,"Proton Chi2 = %.1f, Kaon Chi2 = %.1f", 
			 pids[pidpl]->Chi2Proton(),
			 pids[pidpl]->Chi2Kaon());
//		 sprintf(kaon,"Kaon Chi2 = %.1f", 
//			 pids[pidpl]->Chi2Kaon());
		 sprintf(pion,"Pion Chi2 = %.1f, Muon Chi2 = %.1f", 
			 pids[pidpl]->Chi2Pion(),
			 pids[pidpl]->Chi2Muon());
//		 sprintf(muon,"Muon Chi2 = %.1f", 
//			 pids[pidpl]->Chi2Muon());
		 sprintf(pida,"Plane %d, PIDA = %.1f, NHits = %d",
			 calos[calopl]->PlaneID().Plane,
			 pids[pidpl]->PIDA(),
                         int(calos[calopl]->dEdx().size()));

		 double offset = ((double)trkIter)*10.0;
		 TLatex& track_tex  = view->AddLatex(13.0, (46.0)     - offset,trackinfo);
		 TLatex& pida_tex   = view->AddLatex(13.0, (46.0-2.5) - offset,pida);
		 TLatex& proton_tex = view->AddLatex(13.0, (46.0-5.0) - offset,proton);
		 //TLatex& kaon_tex   = view->AddLatex(13.0, (46.0-4.0) - offset,kaon);
		 TLatex& pion_tex   = view->AddLatex(13.0, (46.0-7.5) - offset,pion);
		 //TLatex& muon_tex   = view->AddLatex(13.0, (46.0-8.0) - offset,muon);
		 track_tex.SetTextColor(evd::kColor[color]);
		 proton_tex.SetTextColor(evd::kColor[color]);
		 //kaon_tex.SetTextColor(evd::kColor[color]);
		 pion_tex.SetTextColor(evd::kColor[color]);
		 //muon_tex.SetTextColor(evd::kColor[color]);
		 pida_tex.SetTextColor(evd::kColor[color]);
		 track_tex.SetTextSize(0.05);
		 proton_tex.SetTextSize(0.05);
		 //kaon_tex.SetTextSize(0.05);
		 pion_tex.SetTextSize(0.05);
		 //muon_tex.SetTextSize(0.05);
		 pida_tex.SetTextSize(0.05);
	       }
            }
         }
      }
   }
  
  //......................................................................
   void AnalysisBaseDrawer::DrawKineticEnergy(const art::Event& evt,
                                              evdb::View2D* view)
   {
      art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
      art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;
      art::ServiceHandle<geo::Geometry>               geom;
      //add some legend-like labels with appropriate grayscale
      char proton[80];
      char kaon[80];
      char pion[80];
      char muon[80];
      sprintf(proton,"proton");
      sprintf(kaon,"kaon");
      sprintf(pion,"pion");
      sprintf(muon,"muon");
      TLatex& proton_tex = view->AddLatex(2.0, 180.0,proton);
      TLatex& kaon_tex   = view->AddLatex(2.0, 165.0,kaon);
      TLatex& pion_tex   = view->AddLatex(2.0, 150.0,pion);
      TLatex& muon_tex   = view->AddLatex(2.0, 135.0,muon);
      proton_tex.SetTextColor(kBlack);
      kaon_tex.SetTextColor(kGray+2);
      pion_tex.SetTextColor(kGray+1);
      muon_tex.SetTextColor(kGray);
      proton_tex.SetTextSize(0.075);
      kaon_tex.SetTextSize(0.075);
      pion_tex.SetTextSize(0.075);
      muon_tex.SetTextSize(0.075);
      
      //now get the actual data
      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
         //Get Track collection
         std::string which = recoOpt->fTrackLabels[imod];
         art::Handle<std::vector<recob::Track> > trackListHandle;
         evt.getByLabel(which,trackListHandle);
         std::vector<art::Ptr<recob::Track> > tracklist;
         art::fill_ptr_vector(tracklist, trackListHandle);
       
         //Loop over Calorimetry collections
         for(size_t cmod = 0; cmod < anaOpt->fCalorimetryLabels.size(); ++cmod) {
            std::string const callabel = anaOpt->fCalorimetryLabels[cmod];
            //Association between Tracks and Calorimetry
            art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, callabel);
	    if (!fmcal.isValid()) continue;

            //Loop over PID collections
            for(size_t pmod = 0; pmod < anaOpt->fParticleIDLabels.size(); ++pmod) {
               std::string const pidlabel = anaOpt->fParticleIDLabels[pmod];
               //Association between Tracks and PID
               art::FindMany<anab::ParticleID> fmpid(trackListHandle, evt, pidlabel);
	       if (!fmpid.isValid()) continue;

               //Loop over Tracks
               for(size_t trkIter = 0; trkIter<tracklist.size(); ++trkIter){
		 int color = tracklist[trkIter]->ID()%evd::kNCOLS;

		 std::vector<const anab::Calorimetry*> calos = fmcal.at(trkIter);
		 if (!calos.size()) continue;
		 size_t bestplane = 0;
		 size_t nmaxhits = 0;
		 for (size_t icalo = 0; icalo < calos.size(); ++icalo){
		   if (calos[icalo]->dEdx().size() > nmaxhits){
		     nmaxhits = calos[icalo]->dEdx().size();
		     bestplane = icalo;
		   }
		 }
		 if (anaOpt->fCaloPlane>=0 and anaOpt->fCaloPlane<int(geom->Nplanes())){
		   for (size_t i = 0; i<geom->Nplanes(); ++i){
		     if (int(calos[i]->PlaneID().Plane)==anaOpt->fCaloPlane)
		       bestplane = i;
		   }
		 }

		 double xvalue = calos[bestplane]->Range();
		 double yvalue = calos[bestplane]->KineticEnergy();
		 view->AddMarker(xvalue,yvalue,evd::kColor[color],8,0.8);
		 if(yvalue>0.0){
		   double error = yvalue*(0.6064/std::sqrt(yvalue));
		   TLine& l = view->AddLine(xvalue,yvalue-error,xvalue,yvalue+error);
		   l.SetLineColor(evd::kColor[color]);
		 }

               }
            }
         }
      }
  
   }

 //......................................................................
   void AnalysisBaseDrawer::CalorShower(const art::Event& evt,
					     evdb::View2D* view)
   {
     art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
     art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;

      for(size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
       
         //Get Track collection
         std::string which = recoOpt->fShowerLabels[imod];
         art::Handle<std::vector<anab::Calorimetry> > caloListHandle;
         evt.getByLabel(which,caloListHandle);
         std::vector<art::Ptr<anab::Calorimetry> > calolist;
         art::fill_ptr_vector(calolist, caloListHandle);
       
         
            //Loop over PID collections
            for(size_t pmod = 0; pmod < anaOpt->fParticleIDLabels.size(); ++pmod) {
               std::string const pidlabel = anaOpt->fParticleIDLabels[pmod];
               //Association between Tracks and PID
               
               //Loop over Tracks
               for(size_t shwIter = 0; shwIter<calolist.size(); ++shwIter){
		 int color = kRed;

		    TPolyMarker& pm = view->AddPolyMarker((*calolist.at(shwIter)).dEdx().size(),color,8,0.8);
                     for(size_t h = 0; h<(*calolist.at(shwIter)).dEdx().size();++h){
                        pm.SetPoint(h,(*calolist.at(shwIter)).ResidualRange().at(h),(*calolist.at(shwIter)).dEdx().at(h));
                     }
                  
                  
                
               }
            }
      }

      char mip[80];
      char mip2[80];

     sprintf(mip,"1 MIP");
     sprintf(mip2,"2 MIP");   
     double offset = 0;

     double MIP = 2.12;  // This is one mip in LAr, taken from uboone docdb #414 
     TLine & Line1Mip = view->AddLine(0, MIP, 100, MIP);
     TLine & Line2Mip = view->AddLine(0, 2*MIP, 100, 2*MIP);

     TLatex& mip_tex   = view->AddLatex(40.0, (23.0-20.0) - offset,mip);
     TLatex& mip2_tex  = view->AddLatex(40.0, (23.0-18.0) - offset,mip2);     

     mip_tex.SetTextColor(kGray+3);
     mip2_tex.SetTextColor(kGray+2);
     mip_tex.SetTextSize(0.02);
     mip2_tex.SetTextSize(0.02);
   
     Line1Mip.SetLineStyle(kDashed);
     Line1Mip.SetLineColor(kGray+3);
     Line2Mip.SetLineStyle(kDashed);
     Line2Mip.SetLineColor(kGray+2);
   }

   //......................................................................
   void AnalysisBaseDrawer::CalorInteractive(const art::Event& /*evt*/,
					     evdb::View2D* view,
					     trkf::BezierTrack BTrack,
					     trkf::HitPtrVec HitHider
					     )
   {
     
     // slightly dirty workaround to get the hits passed on
     std::vector<art::Ptr<recob::Hit> > hits = HitHider.Hits;
     
     art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
     art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;
     art::ServiceHandle<geo::Geometry>               geom;
     
     

     calo::CalorimetryAlg calalg(recoOpt->fCaloPSet);

     

     anab::Calorimetry cal = BTrack.GetCalorimetryObject(hits, geo::kCollection, calalg);
     
     TPolyMarker& pm = view->AddPolyMarker(cal.dEdx().size(),kOrange,8,0.8);
     for(size_t h = 0; h!=cal.dEdx().size();++h){
       pm.SetPoint(h,cal.ResidualRange().at(h), cal.dEdx().at(h));
       


     }

     char proton[80];
     char kaon[80];
     char pion[80];
     char muon[80];
     char mip[80];
     char mip2[80];

     sprintf(proton,"Proton");
     sprintf(kaon,"Kaon");   
     sprintf(pion,"Pion");
     sprintf(muon,"Muon");   

     sprintf(mip,"1 MIP");
     sprintf(mip2,"2 MIP");   
     double offset = 0;
     TLatex& proton_tex = view->AddLatex(40.0, (23.0-1.0) - offset,proton);
     TLatex& kaon_tex   = view->AddLatex(40.0, (23.0-2.0) - offset,kaon);
     TLatex& pion_tex   = view->AddLatex(40.0, (23.0-3.0) - offset,pion);
     TLatex& muon_tex   = view->AddLatex(40.0, (23.0-4.0) - offset,muon);
     TLatex& mip_tex   = view->AddLatex(40.0, (23.0-20.0) - offset,mip);
     TLatex& mip2_tex  = view->AddLatex(40.0, (23.0-18.0) - offset,mip2);
     proton_tex.SetTextColor(Style::ColorFromPDG(2212));
     kaon_tex.SetTextColor(Style::ColorFromPDG(321));
     pion_tex.SetTextColor(Style::ColorFromPDG(211));
     muon_tex.SetTextColor(Style::ColorFromPDG(13));


     proton_tex.SetTextSize(0.03);
     kaon_tex.SetTextSize(0.03);
     pion_tex.SetTextSize(0.03);
     muon_tex.SetTextSize(0.03);

     mip_tex.SetTextColor(kGray+3);
     mip2_tex.SetTextColor(kGray+2);
     mip_tex.SetTextSize(0.02);
     mip2_tex.SetTextSize(0.02);

     double MIP = 1.5 * 1.4;
     TLine & Line1Mip = view->AddLine(0, MIP, 100, MIP);
     TLine & Line2Mip = view->AddLine(0, 2*MIP, 100, 2*MIP);
   
     Line1Mip.SetLineStyle(kDashed);
     Line1Mip.SetLineColor(kGray+3);
     Line2Mip.SetLineStyle(kDashed);
     Line2Mip.SetLineColor(kGray+2);
     
     
   }
  
 

}// namespace
////////////////////////////////////////////////////////////////////////
