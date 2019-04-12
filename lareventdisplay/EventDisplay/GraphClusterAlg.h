////////////////////////////////////////////////////////////////////////
// GraphClusterAlg.h
//
// GraphClusterAlg class
//
// Andrzej Szelc (andrzej.szelc@yale.edu)
//
////////////////////////////////////////////////////////////////////////
#ifndef GRAPHCLUSTERALG_H
#define GRAPHCLUSTERALG_H

#include <vector>

#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace art {
  class Event;
}

#ifdef __ROOTCLING__
namespace fhicl {
  class ParameterSet;
}

namespace recob {
 class Hit;
}


#else
#include "lareventdisplay/EventDisplay/InfoTransfer.h"
#include "larcore/Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"
#endif


namespace util {
 class PxLine;
 class PxPoint;
}

namespace geo {
  class Geometry;
}


namespace recob {
  class Hit;
  class Cluster;
}

namespace evd {

  class InfoTransfer;


  class GraphClusterAlg {

  public:

  GraphClusterAlg(fhicl::ParameterSet const& pset);

  void reconfigure(fhicl::ParameterSet const& pset);

//   void GetStartEndHits(unsigned int plane, recob::Hit * starthit,recob::Hit * endhit);
//   void GetStartEndHits(unsigned int plane);
  void GetStartEndHits(unsigned int plane,util::PxLine &startendpoints);



    //void GetHitList(unsigned int plane,std::vector< art::Ptr <recob::Hit> > ptrhitlist);
  void GetHitList(unsigned int plane, art::PtrVector <recob::Hit>  &ptrhitlist);

  void GetHitListAndEndPoints(unsigned int plane, art::PtrVector <recob::Hit>  &ptrhitlist,util::PxLine &startendpoints);

  int CheckValidity(art::Event& evt);

  private:
    std::vector < util::PxLine > GetSeedLines();



    unsigned int fNPlanes;

    int TestFlag;
    int fRun;
    int fSubRun;
    int fEvent;



    /*
    std::vector< recob::Hit * > starthit;
    std::vector< recob::Hit * > endhit;
    */
//     std::vector < std::vector< recob::Hit * > > hitlist;

//     std::vector < util::PxLine > plines;
//
//     std::vector <unsigned int> swire;
//     std::vector <unsigned int> ewire;
//     std::vector <double> stime;
//     std::vector <double> etime;
//



  }; //class GraphClusterAlg

} //namespace evd





#endif
