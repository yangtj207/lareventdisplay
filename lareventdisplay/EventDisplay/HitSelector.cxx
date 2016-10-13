/// \file    HitSelector.cxx
/// \brief   /// Class to perform operations needed to select hits and pass them to InfoTransfer.
/// \author  andrzej.szelc@yale.edu
/// \author ellen.klein@yale.edu
////////////////////////////////////////////////


#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nutools/EventDisplayBase/evdb.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/PxHitConverter.h"
#include "lareventdisplay/EventDisplay/InfoTransfer.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lardataobj/RecoBase/Seed.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "TMath.h"
namespace{
  void WriteMsg(const char* fcn)
  {
    mf::LogVerbatim("HitSelector") << "HitSelector::" << fcn << " \n";
  }
}

namespace evd{

  //----------------------------------------------------------------------------
  HitSelector::HitSelector() 
  {
    starthitout.resize(3);
    endhitout.resize(3);
  }
  
  //----------------------------------------------------------------------------
  HitSelector::~HitSelector() 
  {
 
  }


  //......................................................................
  ///
  /// Save SeedLineList to infoTransfer service
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// 
  /// Return value is the kinetic enegry of the track
  //
  double HitSelector::SaveSeedLines(const art::Event& evt,
				    evdb::View2D*     /*view*/,
				    std::vector< util::PxLine > seedlines,
				    double distance
				    ) 
  {
    art::ServiceHandle<evd::InfoTransfer>   infot;
    art::ServiceHandle<geo::Geometry>       geom;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    std::map<int, std::vector < art::Ptr < recob::Hit> > > hits_to_save;

    double KineticEnergy=0;
    art::Handle< std::vector<recob::Hit> > HitListHandle;


    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      std::string const which = recoOpt->fHitLabels[imod];

      evt.getByLabel(which,HitListHandle);


      art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutoptions;
      if(evdlayoutoptions->fMakeSeeds)
	{
	  // Code to collect hits around curve

	  // First, put the hits into a PtrVector
	  art::PtrVector<recob::Hit> HitVec;
	  for(unsigned int i=0; i < HitListHandle->size(); ++i)
	    {
	      art::Ptr<recob::Hit> hit(HitListHandle,i);
	      HitVec.push_back(hit);
	    }

	  // Get the bezier curve which was created
	  if(SeedVector().size()>1)
	    {
	      mf::LogVerbatim("HitSelector") <<"Hit selector passing bezier track " 
					     << SeedVector().size() << " seeds";
	    }
	  else if(SeedVector().size()==1)
	    {
	      mf::LogVerbatim("HitSelector") <<"One seed defined - splitting in half to make track";
	      double Pt[3], Dir[3], Err[3];
	      SeedVector().at(0).GetPoint(     Pt,  Err);
	      SeedVector().at(0).GetDirection( Dir, Err);
	     
	      double NewPt1[3], NewPt2[3], NewDir[3];
	      for(int i=0; i!=3; ++i)
		{
		  NewPt1[i] = Pt[i] - Dir[i]/2;
		  NewPt2[i] = Pt[i] + Dir[i]/2;
		  NewDir[i] = Dir[i] / 2;
		}
	      SeedVector().clear();
	      
	      recob::Seed TheSeed1(NewPt1, NewDir, Err, Err);
	      SeedVector().push_back(TheSeed1);
	      recob::Seed TheSeed2(NewPt2, NewDir, Err, Err);
	      SeedVector().push_back(TheSeed2);
	    }
	  else
	    {
	      mf::LogWarning("HitSelector") << "Cannot select hits, no seeds selected!";
	      return 0;
	    }
	   
	  trkf::BezierTrack BTrack(SeedVector());

	  // These hold output,
	  //  s = fractional distance along curve of close approach
	  //  d = min distance from line
	  std::vector<double> s;
	  std::vector<double> d;

	  // Run through the hits, getting their distances
	  BTrack.GetClosestApproaches(HitVec, s, d);
	 
	  // Put them into the vector for output
	  for(size_t ih=0; ih!=HitVec.size(); ++ih){
	    if(d.at(ih)<distance){
	      hits_to_save[HitVec.at(ih)->WireID().Plane].push_back(HitVec.at(ih));
	    } 
	  }
	  art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
	  calo::CalorimetryAlg calalg(recoOpt->fCaloPSet);
	  for(std::map<int, std::vector<art::Ptr<recob::Hit> > >::iterator it=hits_to_save.begin(); it!=hits_to_save.end(); ++it){
	    KineticEnergy += BTrack.GetCalorimetryObject(it->second, geo::kCollection,calalg).KineticEnergy();
	    
	  }
	  mf::LogVerbatim("HitSelector") <<"Track kinetic energy : " << KineticEnergy;
	}
    }
   
   
    infot->SetSeedList(seedlines);
    for(std::map<int, std::vector<art::Ptr<recob::Hit> > >::iterator it=hits_to_save.begin(); it!=hits_to_save.end(); ++it)
      {
	infot->SetHitList(it->first, it->second);
      }
    infot->SetTestFlag(1);
    infot->SetEvtNumber(evt.id().event());
   
    return KineticEnergy;
  }

  //......................................................................
  ///
  /// Save HitList to infoTransfer service
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void HitSelector::SaveHits(const art::Event& evt,
			     evdb::View2D*     /*view*/,
			     unsigned int      plane,
			     double x,  double y,
			     double x1, double y1,
			     double distance,
			     bool good_plane
			     ) 
  {
    // unsigned int w  = 0;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>  geo;
    util::GeometryUtilities  gser;
    std::vector < art::Ptr < recob::Hit> > hits_to_save; //std::vector < azart::Ptr >


    //std::vector < recob::Hit> hits_to_point; //needed to convert art::Ptr<recob:: Hit> to Hit* to pass to Hit2D
   // std::vector< const recob::Hit*> hits_to_draw; //draw selected hits in a different color

    art::Handle< std::vector<recob::Hit> > HitListHandle;
    //  recob::Hit * startHit;
    //  recob::Hit * endHit;
    //preprocess event - load up all the hits with std::vector, as in BackTracker
    
	  
		  
    starthitout[plane].clear();
    endhitout[plane].clear();
 
    starthitout[plane].resize(2);
    endhitout[plane].resize(2);
	  
	
   
    
    
    double lslope=0;
    
    if((x-x1)!=0)
      lslope=(y-y1)/(x-x1);
    
     for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      std::string const which = recoOpt->fHitLabels[imod];
    
     std::vector<art::Ptr<recob::Hit>> hitlist; 
     hitlist.clear();
     evt.getByLabel(which,HitListHandle);
     
    for(unsigned int ii = 0; ii < HitListHandle->size(); ++ii){
       art::Ptr<recob::Hit> hit(HitListHandle, ii);
       if(hit->WireID().Plane==plane)
	  hitlist.push_back(hit);
      }
    
    
    // Select Local Hit List
     util::PxHitConverter PxC;
     std::vector <util::PxHit> pxhitlist;
     PxC.GeneratePxHit(hitlist,pxhitlist);
     std::vector< unsigned int > pxhitlist_local_index;
     std::vector< util::PxHit > pxhitlist_local;
     pxhitlist_local.clear();
     
     util::PxPoint startHit;
     startHit.plane=pxhitlist.at(0).plane;
     startHit.w=(x+x1)/2;
     startHit.t=(y+y1)/2;
     
     double orttemp=TMath::Sqrt((y1-y)*(y1-y)*gser.TimeToCm()*gser.TimeToCm() + (x1-x)*(x1-x)*gser.WireToCm()*gser.WireToCm())/2;
     
     gser.SelectLocalHitlistIndex(pxhitlist, pxhitlist_local_index, startHit, 
			     orttemp, 
			     distance, lslope);
    

     for(unsigned int idx=0;idx<pxhitlist_local_index.size();idx++)
       {
        hits_to_save.push_back(hitlist.at(pxhitlist_local_index.at(idx)));
	pxhitlist_local.push_back(pxhitlist.at(pxhitlist_local_index.at(idx)));
       }
    
    
   /*  gser.SelectLocalHitlist( hitlist, 
			 hits_to_save,
			   (x+x1)/2,
			   (y+y1)/2, 
			   TMath::Sqrt((y1-y)*(y1-y)*gser.TimeToCm()*gser.TimeToCm() + (x1-x)*(x1-x)*gser.WireToCm()*gser.WireToCm())/2,   
			   distance, 
			   lslope);*/		
    //const_cast<recob::Hit *> (nearHit.get())
     recob::Hit * hit=  const_cast<recob::Hit *> ((hits_to_save[gser.FindClosestHitIndex(pxhitlist_local,x,y)]).get());
     //gser.FindClosestHit(hits_to_save,
	//		      x, y);
      
       starthitout[plane][1] = hit->PeakTime() ;  
       starthitout[plane][0] = hit->WireID().Wire;
 
//        // this obviously not correct, the fact that x,y are used for both start and end point, A.S. -> to debug
       recob::Hit * endhit= const_cast<recob::Hit *> ((hits_to_save[gser.FindClosestHitIndex(pxhitlist_local,x,y)]).get());
       //FindClosestHit(hits_to_save,
	//		      x, y);
      
        endhitout[plane][1] = endhit->PeakTime() ;  
        endhitout[plane][0] = endhit->WireID().Wire;
    
     }
     
    art::ServiceHandle<evd::InfoTransfer>   infot;
    infot->SetHitList(plane,hits_to_save);
    infot->SetTestFlag(1);
    if(good_plane)
      {
	infot->SetStartHitCoords(plane,starthitout[plane]);
	infot->SetEndHitCoords(plane,endhitout[plane]);
      } 
    infot->SetEvtNumber(evt.id().event());
   
  }

  //......................................................................
  ///
  /// Save HitList to infoTransfer service
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void HitSelector::ChangeHit(const art::Event& evt,
			      evdb::View2D*     /*view*/,
			      unsigned int      plane,
			      double x,  double y
			      ) 
  {
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>  geo;
    util::GeometryUtilities  gser;
    
    
    //get hits from info transfer, see if our selected hit is in it
    std::vector < art::Ptr<recob::Hit> > hits_saved;
      std::vector < art::Ptr<recob::Hit> > hitlist;
    std::vector < art::Ptr < recob::Hit> > hits_to_save; 
     hitlist.clear();


    art::Handle< std::vector<recob::Hit> > HitListHandle;
   
    //preprocess event - load up all the hits with std::vector, as in BackTracker
    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      std::string const which = recoOpt->fHitLabels[imod];
      evt.getByLabel(which,HitListHandle);
    
     for(unsigned int ii = 0; ii < HitListHandle->size(); ++ii){
       art::Ptr<recob::Hit> hit(HitListHandle, ii);
        if(hit->WireID().Plane==plane)
	  hitlist.push_back(hit);
      }
      
     util::PxHitConverter PxC;
     std::vector <util::PxHit> pxhitlist;
     PxC.GeneratePxHit(hitlist,pxhitlist);
      
      art::Ptr < recob::Hit > selected_hit= hits_to_save[gser.FindClosestHitIndex(pxhitlist,x,y)];
      //FindClosestHitPtr(hitlist, x, y); 
      
      // find selected hit in evD
     //art::Ptr<recob::Hit> selected_hit;
//       for(unsigned int ii = 0; ii < HitListHandle->size(); ++ii){
// 	art::Ptr<recob::Hit> hit(HitListHandle, ii);
// 	double time=hit->PeakTime();
// 	if(hit->WireID().Plane != plane)
// 	  continue;
// 	if(hit->WireID().Wire < (x+1) &&
// 	   hit->WireID().Wire > (x-1) ){
// 	  if(y>(time-5) && y<(time+5)){
// 	    selected_hit=hit;
// 	  }
// 	}
		  
    //}
      if(selected_hit.isNull()){
	WriteMsg("no luck finding hit in evd, please try again");
	break;
      }
   
    
      art::ServiceHandle<evd::InfoTransfer>   infot;
      hits_saved=infot->GetSelectedHitList(plane);
      int found_it=0;
      for(unsigned int jj=0; jj<hits_saved.size(); ++jj){
	if(selected_hit->PeakTime() == hits_saved[jj]->PeakTime()){
	  if(selected_hit->Channel() == hits_saved[jj]->Channel()){
	    found_it=1;
	    hits_saved.erase(hits_saved.begin()+jj);
	  }
	}
      }
    
      //if didn't find it, add it
      if(found_it!=1){
	hits_saved.push_back(selected_hit);
      }
    
      //update the info transfer list
      infot->SetHitList(plane,hits_saved);
      infot->SetTestFlag(1);
      infot->SetEvtNumber(evt.id().event());
    }
  }
  


  //......................................................................
  ///
  /// Save HitList to infoTransfer service
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  
  trkf::HitPtrVec  HitSelector::GetSelectedHitPtrs(unsigned int plane)
  {
    art::ServiceHandle<evd::InfoTransfer>   infot;
    trkf::HitPtrVec ToReturn;
    ToReturn.Hits = infot->GetSelectedHitList(plane);
    return ToReturn;
  }

  std::vector< const recob::Hit*>  HitSelector::GetSelectedHits(unsigned int plane)
  {
    std::vector < art::Ptr<recob::Hit> > hits_saved;
    //std::vector < art::Ptr < recob::Hit> > hits_to_save; 

    //std::vector < recob::Hit> hits_to_point; //needed to convert art::Ptr<recob:: Hit> to Hit* to pass to Hit2D
    std::vector< const recob::Hit*> hits_to_draw; //draw selected hits in a different color
    
    art::ServiceHandle<evd::InfoTransfer>   infot;
    hits_saved=infot->GetSelectedHitList(plane);
	
    // if(hits_saved.size()!=0){
    for(unsigned int i=0;i<hits_saved.size();i++)
      hits_to_draw.push_back(hits_saved[i].get());
    //}
    
    return hits_to_draw;     
    
  }
  
  
  //......................................................................
  void HitSelector::ClearHitList(unsigned int plane)
  {
    art::ServiceHandle<evd::InfoTransfer>   infot;
    infot->ClearSelectedHitList(plane);
  }
  
  //----------------------------------------------------------------------------
  std::vector<recob::Seed>& HitSelector::SeedVector()
  {
    return fSeedVector;
  }
  
  
  
} //end namespace
