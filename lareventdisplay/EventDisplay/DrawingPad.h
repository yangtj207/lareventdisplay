///
/// \file    DrawingPad.h
/// \brief   Base class for all event display drawing pads
/// \author  messier@indiana.edu
///
#ifndef EVD_DRAWINGPAD_H
#define EVD_DRAWINGPAD_H

#include <memory>

class TPad;

namespace evd_tool {
  class IExperimentDrawer;
}

namespace evd {
  class HeaderDrawer;
  class SimulationDrawer;
  class RawDataDrawer;
  class RecoBaseDrawer;
  class AnalysisBaseDrawer;
  class HitSelector;

  /// Base class for event display drawing pads
  class DrawingPad {
  public:
    DrawingPad(const char* nm, const char* ti, double x1, double y1, double y2, double x2);
    ~DrawingPad();
    TPad* Pad() { return fPad; }

    // Access to the drawing utilities
    HeaderDrawer* HeaderDraw();
    evd_tool::IExperimentDrawer* GeometryDraw();
    SimulationDrawer* SimulationDraw();
    RawDataDrawer* RawDataDraw();
    RecoBaseDrawer* RecoBaseDraw();
    AnalysisBaseDrawer* AnalysisBaseDraw();
    HitSelector* HitSelectorGet();

  protected:
    using IExperimentDrawerPtr = std::unique_ptr<evd_tool::IExperimentDrawer>;

    TPad* fPad;                            ///< The ROOT graphics pad
    HeaderDrawer* fHeaderDraw;             ///< Drawer for event header info
    IExperimentDrawerPtr fGeometryDraw;    ///< Drawer for detector geometry
    SimulationDrawer* fSimulationDraw;     ///< Drawer for simulation objects
    RawDataDrawer* fRawDataDraw;           ///< Drawer for raw data
    RecoBaseDrawer* fRecoBaseDraw;         ///< Drawer for recobase objects
    AnalysisBaseDrawer* fAnalysisBaseDraw; ///< Drawer for analysisbase objects
  };
}
#endif
////////////////////////////////////////////////////////////////////////
