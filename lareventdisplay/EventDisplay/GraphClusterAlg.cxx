////////////////////////////////////////////////////////////////////////
//
// GraphClusterAlg class
//
// andrzej.szelc@yale.edu
// ellen.klein@yale.edu
//
//  Methods to use by a dummy producer
////////////////////////////////////////////////////////////////////////
#include "lareventdisplay/EventDisplay/GraphClusterAlg.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"


 //-------------------------------------------------
  evd::GraphClusterAlg::GraphClusterAlg(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset);
    art::ServiceHandle<geo::Geometry>  geo;
	
	
    
	
	
    fNPlanes = geo->Nplanes();
//     starthit.resize(fNPlanes);
//     endhit.resize(fNPlanes);
    
   /* 
    swire.resize(fNPlanes);
    ewire.resize(fNPlanes);
    stime.resize(fNPlanes);
    etime.resize(fNPlanes);*/
  }
  
  
 
  //-------------------------------------------------
  void evd::GraphClusterAlg::reconfigure(fhicl::ParameterSet const& /*pset*/)
  {


    return;
  }


  
  void evd::GraphClusterAlg::GetHitListAndEndPoints(unsigned int plane, art::PtrVector <recob::Hit>  &ptrhitlist,util::PxLine &startendpoints)
  {
   GetHitList(plane,ptrhitlist);
   GetStartEndHits(plane,startendpoints);
    
       
  }
  
 
  
  void evd::GraphClusterAlg::GetStartEndHits(unsigned int plane,util::PxLine &startendpoints)
  {
    std::vector < double > starthit;
    std::vector < double > endhit;
    art::ServiceHandle<evd::InfoTransfer> intr;
    starthit=intr->GetStartHitCoords(plane);
    endhit=intr->GetEndHitCoords(plane);
  
    startendpoints.w0=starthit[0];
    startendpoints.t0=starthit[1];
    startendpoints.w1=endhit[0];
    startendpoints.t1=endhit[1];
    startendpoints.plane=plane;
    
  }
  

//   //----------------------------------------------------------------------------
//   void evd::GraphClusterAlg::GetStartEndHits(unsigned int plane)
//   {
//     std::vector < double > starthit;
//     std::vector < double > endhit;
//     art::ServiceHandle<evd::InfoTransfer> intr;
//     starthit=intr->GetStartHitCoords(plane);
//     endhit=intr->GetEndHitCoords(plane);
//   
//    
//     swire[plane]=starthit[0];
//     stime[plane]=starthit[1];
//     ewire[plane]=endhit[0];
//     etime[plane]=endhit[1];
//   
//   }



//   //----------------------------------------------------------------------------
//   void evd::GraphClusterAlg::GetStartEndHits(unsigned int plane, 
// 				     recob::Hit *starthit,
// 				     recob::Hit *endhit) 
//   {
//    
//     art::ServiceHandle<evd::InfoTransfer> intr;
//     art::ServiceHandle<geo::Geometry>  geo;
//    
//     starthit=intr->GetStartHit(plane);
//     endhit=intr->GetEndHit(plane);
// 		
//     // error checking for bogus transfers		
//     if(starthit!=NULL){
//       stime[plane] = starthit->PeakTime() ;  
//       try{
// 	if(starthit->Wire() != NULL){
// 	  swire[plane] = starthit->WireID().Wire;
// 	}
// 	else{
// 	  swire[plane]=0;
// 	}
//       }
//       catch(cet::exception e) {
// 	mf::LogWarning("GraphClusterAlg") << "caught exception \n"
// 				       << e;
// 	swire[plane]=0;			   
//       }
//     }
//     else{
//       stime[plane]=0.;
//     }
// 			
//     if(endhit!=NULL){
//       etime[plane] = endhit->PeakTime(); 
//       if(endhit->Wire()!=NULL){
// 	ewire[plane] = endhit->WireID().Wire;
//       }
//       else{
// 	ewire[plane]=0;
//       }
//     }
//     else{
//       etime[plane]=0.;
//     }
//    
//   }
 
  //---------------------------------------------------------------------------- 
  //  void evd::GraphClusterAlg::GetHitList(unsigned int plane,std::vector< art::Ptr <recob::Hit> > ptrhitlist)
  //  {
  //   art::ServiceHandle<evd::InfoTransfer> intr;
  //    
  //   
  //   ptrhitlist=intr->GetHitList(plane);
  //   //std::vector <recob::Hit *> hitlist_out;
  //   
  //   if(ptrhitlist.size()==0) {
  // 	WriteMsg("hit list of zero size, please select some hits");
  // 		return;
  // 	}
  //     
  //    for(art::PtrVector<recob::Hit>::const_iterator hitIter = ptrhitlist.begin(); hitIter != ptrhitlist.end();  hitIter++){
  // // 	art::Ptr<recob::Hit> theHit = (*hitIter);
  // // 	unsigned int plane,cstat,tpc,wire;
  // //	hitlist_out.push_back((*hitIter)->Get());
  //  	} 
  //    
  //    
  //    
  //    
  //    return;// hitlist_out;
  //  }
  
  //----------------------------------------------------------------------------
  void evd::GraphClusterAlg::GetHitList(unsigned int plane, art::PtrVector <recob::Hit>  &ptrhitlist)
  {
    art::ServiceHandle<evd::InfoTransfer> intr;
   
    std::vector< art::Ptr <recob::Hit> > ptlist=intr->GetHitList(plane);
  
 
    //std::vector <recob::Hit *> hitlist_out;
  
    if(ptlist.size()==0) {
      mf::LogVerbatim("GraphClusterAlg") << ("hit list of zero size, please select some hits");
      return;
    }
    
    for(art::PtrVector<recob::Hit>::const_iterator hitIter = ptlist.begin(); hitIter != ptlist.end();  hitIter++){
      // 	art::Ptr<recob::Hit> theHit = (*hitIter);
      // 	unsigned int plane,cstat,tpc,wire;
      ptrhitlist.push_back((*hitIter));
    } 
      
    return;// hitlist_out;
  }
 
 
  //----------------------------------------------------------------------------
  std::vector < util::PxLine > evd::GraphClusterAlg::GetSeedLines()
  {
   
    art::ServiceHandle<evd::InfoTransfer> intr;
    //////////////////////////////////////////////////
    //this is where you could create Bezier Tracks if you wanted to do it inside a producer	 
    //////////////////////////////////////////////////
    std::vector < util::PxLine > plines = intr->GetSeedList();
	
    std::cout << " Received Seed List of Size: " << plines.size() << std::endl; 
    
    return plines;
  }


 int evd::GraphClusterAlg::CheckValidity(art::Event& evt){
    art::ServiceHandle<evd::InfoTransfer> intr;
    TestFlag=intr->GetTestFlag();
  
    fEvent=intr->GetEvtNumber();
    fRun=intr->GetRunNumber();
    fSubRun=intr->GetSubRunNumber();
      
   
    
    if(TestFlag==-1)
      return -1;
    
    if(fEvent!=(int)evt.id().event() || fRun!=(int)evt.id().run() || fSubRun!=(int)evt.id().subRun() ) {
      mf::LogVerbatim("GraphClusterAlg") << (" old event ");
      return -1;
    }
   
    return TestFlag;
  }
