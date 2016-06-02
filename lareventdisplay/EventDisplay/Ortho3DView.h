///
/// \file    Ortho3DView.h
/// \brief   A view showing an orthographic projection of 3D objects.
/// \author  greenlee@fnal.gov
///
/// Only 3D objects (prongs and space points) can be displayed in this view.
///

#ifndef EVD_ORTHO3DVIEW_H
#define EVD_ORTHO3DVIEW_H
#include <vector>

#include "nutools/EventDisplayBase/Canvas.h"

namespace evd {

  class Ortho3DPad;

  /// View of event shoing orthographic view of 3D objects.

  class Ortho3DView : public evdb::Canvas {

  public:

    // Constructor, destructor.

    Ortho3DView(TGMainFrame* mf);
    virtual ~Ortho3DView();

    // Required methods.
    
    const char* Description() const { return "Orthographic 3D Detector Display"; }
    const char* PrintTag()    const { return "larortho3d";               }
    void Draw(const char* opt="");

  private:

    // Attributes.

    // Graphics pads.

    std::vector<Ortho3DPad*> fOrtho3DPads; ///< Graphics pads.

    // Frames.

    TGCompositeFrame* fMetaFrame;    ///< Frame holding root canvas and widget frame.
    TGCompositeFrame* fWidgetFrame;  ///< Frame holding widgets.
    std::vector<TGCompositeFrame*> fWidgetSubFrames; // Frame holding widgets for one pad.
  };
}

#endif
////////////////////////////////////////////////////////////////////////
