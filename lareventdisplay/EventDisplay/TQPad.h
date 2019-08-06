///
/// \file    TQPad.h
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
///
#ifndef EVD_TQPAD_H
#define EVD_TQPAD_H

#include "lareventdisplay/EventDisplay/DrawingPad.h"

#include <memory>

class TH1F;

namespace evdb { class View2D; }
namespace evdb_tool
{
    class IWaveformDrawer;
    class IWFHitDrawer;
}

namespace evd
{
class TQPad : public DrawingPad
{
public:
    TQPad(const char*  nm,
          const char*  ti,
          double       x1,
          double       y1,
          double       x2,
          double       y2,
          const char*  opt,
          unsigned int plane,
          unsigned int wire);

    ~TQPad();

    void Draw();

    void SetPlaneWire(unsigned int plane=0, unsigned int wire=0) { fPlane = plane; fWire = wire; }

private:
    void BookHistogram();

    using IWFHitDrawerPtr    = std::unique_ptr<evdb_tool::IWFHitDrawer>;
    using IWaveformDrawerPtr = std::unique_ptr<evdb_tool::IWaveformDrawer>;

    unsigned int                      fWire;
    unsigned int                      fPlane;               ///< Which plane in the detector
    int                               fTQ;                  ///< 0 = plot shows charge only, 1 = plot shows charge vs time for a wire
    TH1F*                             fFrameHist;           ///< A dummy histogram to define the axes
    evdb::View2D*                     fView;                ///< Superimpose scale on 1D histo
    IWFHitDrawerPtr                   fHitDrawerTool;       ///< An instance of the tool to draw hits
    IWaveformDrawerPtr                fRawDigitDrawerTool;  ///< An instance of the tool to draw hits
    IWaveformDrawerPtr                fWireDrawerTool;      ///< An instance of the tool to draw hits
  };
}

#endif
////////////////////////////////////////////////////////////////////////
