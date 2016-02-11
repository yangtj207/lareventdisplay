/////////////////////////////////////////////////////////////////////////////
///
/// \file    TWireProjPad.h
/// \brief   Drawing pad showing a single X-Z or Y-Z projection of an event
/// \author  messier@indiana.edu
/// \version $Id: TWireProjPad.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
/////////////////////////////////////////////////////////////////////////////

#ifndef EVD_TWIREPROJPAD_H
#define EVD_TWIREPROJPAD_H
#include "lareventdisplay/EventDisplay/DrawingPad.h"
#include <vector>


class TH1F;

namespace evdb { class View2D;   }

namespace util {
  class PxLine;
  class PxPoint;
}

namespace recob
{
  class Seed;
}

namespace evd {

  /// A drawing pad for time vs wire
  class TWireProjPad : public DrawingPad {
  public:
    struct DrawOptions_t {
      bool bZoom2DdrawToRoI = false; ///< whether to force zoom to RoI or not
    }; // DrawOptions_t
    
    TWireProjPad(const char* nm, const char* ti,
		 double x1, double y1,
		 double x2, double y2,
		 unsigned int plane);
    ~TWireProjPad();
    void Draw(const char* opt=0);
    void GetWireRange(int *i1, int *i2) const;
    void SetWireRange(int i1, int i2);

    void SetZoomRange(int i1, int i2,int y1, int y2);
    
    /// Return the current draw options
    DrawOptions_t const& GetDrawOptions() const { return fDrawOpts; }
    /// Receive the full set of draw options
    void SetDrawOptions(DrawOptions_t const& opt) { fDrawOpts = opt; }
    
    /// Sets the draw option about zooming to the region of interest
    void SetZoomToRoI(bool bZoomToRoI)
      { fDrawOpts.bZoom2DdrawToRoI = bZoomToRoI; }
    
    /// Sets the zoom parameters from the current histogram view
    void SetZoomFromView();
    
    void SaveHitList(double i1, double i2,double y1, double y2, double distance, const char* zoom_opt,bool good_plane=true);
			
    double SaveSeedList(std::vector < util::PxLine > seedlines, double distance);
			
    void ClearHitList();
    void SelectOneHit(double x, double y, const char* zoom_opt);
    
    unsigned int GetPlane() const { return fPlane; }

    void ClearandUpdatePad();
    void UpdatePad();
    void DrawLinesinView(std::vector< util::PxLine > lines,bool deleting=false,const char * zoom_opt=0);
			
    void ShowFull(int override=0);

    double UpdateSeedCurve(std::vector<recob::Seed> SeedVec, int plane);
    evdb::View2D*  View() const { return fView; }

    std::vector<double> const& GetCurrentZoom() const {return fCurrentZoom;}
    
  private:
    /*     void AutoZoom(); */
   

  private:

    std::vector<double> fCurrentZoom;
    DrawOptions_t fDrawOpts; ///< set of current draw options


    unsigned int  fPlane; ///< Which plane in the detector
    TH1F*         fHisto; ///< Histogram to draw object on
    evdb::View2D* fView;  ///< Collection of graphics objects to render

    double        fXLo;   ///< Low  value of x axis
    double        fXHi;   ///< High value of x axis
    double        fYLo;   ///< Low  value of y axis
    double        fYHi;   ///< High value of y axis
    int           fOri;   ///< Orientation of the axes - see RawDrawingOptions for values
  };
}

#endif
////////////////////////////////////////////////////////////////////////
