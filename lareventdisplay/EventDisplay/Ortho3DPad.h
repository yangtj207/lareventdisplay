///
/// \file    Ortho3DPad.h
/// \brief   Drawing pad showing an orthographic projection of 3D objects 
///          in the detector
/// \author  greenlee@fnal.gov
///

#ifndef EVD_ORTHO3DPAD_H
#define EVD_ORTHO3DPAD_H

#include "lareventdisplay/EventDisplay/DrawingPad.h"
#include "lareventdisplay/EventDisplay/OrthoProj.h"
#include "RQ_OBJECT.h"
#include <vector>

class TH1F;
class TGNumberEntry;
class TBox;
namespace evdb { class View2D; }

namespace evd {
  
  /// A drawing pad showing an orthographic rendering of 3D objects.

  class Ortho3DPad : public DrawingPad {
  public:

    // So this class can receive gui signals.
    
    RQ_OBJECT("evd::Ortho3DPad")

  public:

    // Constructor, destructor.

    Ortho3DPad(const char* nm, const char* ti,
 	     evd::OrthoProj_t proj,
	     double x1, double y1,
	     double x2, double y2);
    ~Ortho3DPad();

    // Accessors.

    double GetMarkerSize() const {return 4.*fMSize;}  // Size in pixels.

    // Methods.

    void Draw(const char* opt=0);
    void SetZoom(double xlo, double ylo, double xhi, double yhi, bool update);
    void UnZoom(bool update);
    void SetMarkerSize(double size, bool update);  // Size in pixels.

    // Widget-related methods.

    void SetMSizeEntry(TGNumberEntry* p);   // Add number entry widget.
    void SetMSize();                        // Slot for marker size signals.

    // Handler for mouse events.

    static void MouseEvent(evd::Ortho3DPad* p);
    void MouseEvent();

  private:

    // Static attributes.

    static Ortho3DPad* fMousePad;  ///< Selected pad for mouse action.

    // Attributes.

    TH1F* fHisto;             ///< Enclosing histogram.
    evd::OrthoProj_t fProj;   ///< Projection.
    double fXLo;              ///< Low x value.
    double fXHi;              ///< High x value.
    double fYLo;              ///< Low y value.
    double fYHi;              ///< High y value.
    double fMSize;            ///< Marker size.
    std::vector<TBox> TPCBox; ///< TPC box
    evdb::View2D* fView;      ///< Collection of graphics objects to render

    // Widgets.

    TGNumberEntry* fMSizeEntry;   ///< For changing marker size.

    // Mouse/zoom status attributes.

    bool fPress;            ///< Is button 1 pressed?
    bool fBoxDrawn;         ///< Is selection box drawn?
    int fPressPx;           ///< Pixel location where button 1 was pressed.
    int fPressPy;           ///< Poxel location where button 1 was pressed.
    int fCurrentPx;         ///< Current pixel location of mouse.
    int fCurrentPy;         ///< Current pixel location of mouse.
    double fPressX;         ///< User location where button 1 was pressed.
    double fPressY;         ///< User location where button 1 was pressed.
    double fReleaseX;       ///< User location where button 1 was released.
    double fReleaseY;       ///< User location where button 1 was released.
  };
}

#endif
////////////////////////////////////////////////////////////////////////
