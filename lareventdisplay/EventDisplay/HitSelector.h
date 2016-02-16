/// \file    HitSelector.h
/// \brief   Class to perform operations needed to select hits and pass them to a cluster.
/// \author  andrzej.szelc@yale.edu

#ifndef EVD_HITSELECTOR_H
#define EVD_HITSELECTOR_H
#include <vector>
#ifdef __CINT__
namespace art { 
  class Event;
  class PtrVector;
  class Ptr;
  class ServiceHandle;
  class View;
}
namespace trkf
{
 class HitPtrVec;
}
#else
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "lardata/RecoObjects/BezierTrack.h"
#endif
#include "lareventdisplay/EventDisplay/OrthoProj.h"

class TH1F;
namespace evdb { 
  class View2D;   
  class View3D;   
}

namespace recob { 
   class Hit;     
   class Cluster;
   class Seed;
}

 namespace util {
   class PxPoint;
   class PxLine;
 }

/// Class to perform operations needed to select hits and pass them to InfoTransfer.
namespace evd {
  
  class HitSelector {
  public:
    HitSelector();
    ~HitSelector();

    void SaveHits(const art::Event& evt,
		     evdb::View2D*     view,
		     unsigned int      plane,
		      double x,  double y,
		      double x1, double y1,
		      double distance,
		      bool good_plane=true
 		      ); 
    
    double SaveSeedLines(const art::Event& evt,
		       evdb::View2D*     view,
		       std::vector < util::PxLine > seedline,
		       double distance
		       ); 
  
    
     void ChangeHit(const art::Event& evt,
		     evdb::View2D*     view,
		     unsigned int      plane,
		       double x,  double y
 		      ); 
    
    std::vector< const recob::Hit*> GetSelectedHits(unsigned int plane);
    trkf::HitPtrVec GetSelectedHitPtrs(unsigned int plane);
     
    void ClearHitList(unsigned int plane);

    std::vector<recob::Seed>& SeedVector();
    
   private:
      int test;
      std::vector<recob::Seed> fSeedVector;   
    
      std::vector < std::vector <double > > starthitout;
      std::vector < std::vector <double > > endhitout;
      
  };
}

#endif
////////////////////////////////////////////////////////////////////////

