///
/// \file    Display3DPad.h
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
///
#ifndef EVD_DISPLAY3DPAD_H
#define EVD_DISPLAY3DPAD_H

#include "lareventdisplay/EventDisplay/DrawingPad.h"

class TH3F;
namespace evdb      { class View3D;       }
namespace evdb_tool { class ISim3DDrawer; class I3DDrawer;}


namespace evd {
    
class RawDataDrawer;
class RecoBaseDrawer;

/// A drawing pad showing a 3D rendering of the detector
class Display3DPad : public DrawingPad {
public:
    Display3DPad(const char* nm, const char* ti,
                 double x1, double y1,
                 double x2, double y2,
                 const char* opt);
    ~Display3DPad();
  
  
  
    void Draw();
  
    void UpdateSeedCurve();
private:
    evdb::View3D* fView;  ///< Collection of graphics objects to render
    
    std::vector<std::unique_ptr<evdb_tool::ISim3DDrawer>> fSim3DDrawerVec;
    std::vector<std::unique_ptr<evdb_tool::I3DDrawer>>    fReco3DDrawerVec;
};
}

#endif
////////////////////////////////////////////////////////////////////////
