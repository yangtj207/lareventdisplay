/// \file    HitSelector.h
/// \brief   Class to perform operations needed to select hits and pass them to a cluster.
/// \author  andrzej.szelc@yale.edu

#ifndef EVD_HITSELECTOR_H
#define EVD_HITSELECTOR_H

#include <vector>

#include "art/Framework/Principal/fwd.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/Seed.h"

namespace evdb {
  class View2D;
}

namespace recob {
   class Hit;
}

 namespace util {
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
    std::vector<art::Ptr<recob::Hit>> GetSelectedHitPtrs(unsigned int plane);

    void ClearHitList(unsigned int plane);

    std::vector<recob::Seed>& SeedVector();

   private:
      //int test;
      std::vector<recob::Seed> fSeedVector;

      std::vector < std::vector <double > > starthitout;
      std::vector < std::vector <double > > endhitout;

  };
}

#endif
////////////////////////////////////////////////////////////////////////
