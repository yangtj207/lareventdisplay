///
/// \file    MCBriefPad.h
/// \brief   Drawing pad for short summary of an MC event
/// \author  messier@indiana.edu
///
#ifndef EVD_MCBRIEF_H
#define EVD_MCBRIEF_H

#include "lareventdisplay/EventDisplay/DrawingPad.h"
namespace evdb {
  class View2D;
}

namespace evd {
  class MCBriefPad : public DrawingPad {
  public:
    MCBriefPad(const char* nm,
               const char* ti,
               double x1,
               double y1,
               double x2,
               double y2,
               const char* opt);
    ~MCBriefPad();
    void Draw();

  private:
    evdb::View2D* fView; ///< Superimpose scale on 1D histo
  };
}

#endif
////////////////////////////////////////////////////////////////////////
