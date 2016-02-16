///
/// \file    SimulationDrawer.h
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: SimulationDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
///
#ifndef EVD_SIMULATIONDRAWER_H
#define EVD_SIMULATIONDRAWER_H
#include <string>
#include <vector>
#include <map>
#include "lareventdisplay/EventDisplay/OrthoProj.h"

namespace art  { class Event;       }
namespace evdb { class View2D;      }
namespace evdb { class View3D;      }
namespace geo  { class Geometry;    }
namespace simb { 
  class MCTruth;     
  class MCParticle;
}

namespace evd {
  class SimulationDrawer {
  public:
    SimulationDrawer();
    ~SimulationDrawer();

  public:
    // Drawing functions
    void MCTruthShortText(const art::Event& evt,
			  evdb::View2D*     view);
    void MCTruthLongText(const art::Event& evt,
			 evdb::View2D*     view);
    void MCTruthVectors2D(const art::Event& evt,
			  evdb::View2D*     view,
			  unsigned int      plane);
    void MCTruth3D(const art::Event& evt,
		   evdb::View3D*     view);
    void MCTruthOrtho(const art::Event& evt,
		      evd::OrthoProj_t  proj,
		      double            msize,
		      evdb::View2D*     view);

    void HiLite(int trkId, bool hlt=true);

    double minx;
    double maxx;
    double miny;
    double maxy;
    double minz;
    double maxz;

  private:
    int GetMCTruth(const art::Event&                  evt,
		   std::vector<const simb::MCTruth*>& mctruth);

    int GetParticle(const art::Event&                     evt,
		    std::vector<const simb::MCParticle*>& plist);

  private:
    std::map<int,bool>       fHighlite;

  };
}

#endif
////////////////////////////////////////////////////////////////////////
