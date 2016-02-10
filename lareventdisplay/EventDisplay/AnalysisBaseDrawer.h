/// \file    AnalysisBaseDrawer.h
/// \brief   Class to aid in the rendering of AnalysisBase objects
#ifndef EVD_ANALYSISBASEDRAWER_H
#define EVD_ANALYSISBASEDRAWER_H
#include <vector>
#ifdef __CINT__
namespace art { 
  class Event;
  class PtrVector;
  class Ptr;
  class ServiceHandle;
  class View;
}
namespace recob{
  class Hit;
}
#else
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#endif

namespace evdb{
   class View2D;
}

namespace recob{
  class Hit;
}

namespace trkf
{
  class BezierTrack;
  class HitPtrVec;
}

namespace evd {

  /// Aid in the rendering of AnalysisBase objects
  class AnalysisBaseDrawer {
  public:
    AnalysisBaseDrawer();
    ~AnalysisBaseDrawer();

  public:

    void DrawDeDx(const art::Event& evt,
               evdb::View2D* view);

    void DrawKineticEnergy(const art::Event& evt,
                           evdb::View2D* view);

    void CalorShower(const art::Event& evt,
					     evdb::View2D* view   ); 
    
    void CalorInteractive(const art::Event& evt,
			  evdb::View2D* view,
			  trkf::BezierTrack BTrack,
			  trkf::HitPtrVec Hits );

    
  private:
    
  };
}

#endif
////////////////////////////////////////////////////////////////////////
