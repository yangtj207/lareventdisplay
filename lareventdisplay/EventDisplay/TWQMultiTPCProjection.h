////////////////////////////////////////////////////////////////////////
///
/// \file    TWQProjectionView.h
/// \brief   A view showing the time vs wire, charge and charge vs time information for an event
/// \author  brebel@fnal.gov
/// \version $Id: TQPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
// clustering code written by ellen.klein@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef EVD_TWQMULTITPCPROJECTION_H
#define EVD_TWQMULTITPCPROJECTION_H

#include "TGButton.h" // For TGCheckButton

#include "nutools/EventDisplayBase/Canvas.h"

#include "RQ_OBJECT.h"

#include <vector>
#include <deque>

// Forward declarations
class TGNumberEntry;
class TGTextView;
class TGLabel;
class TRootEmbeddedCanvas;

namespace util {
class PxLine;
class PxPoint;

}

namespace evd {

  class MCBriefPad;
  class TQPad;
  class TWireProjPad;
  class HeaderPad;


  // Helper class, to store zoom settings in each view

  class ZoomOptionsMultiTPC
  {
    
  public:
    ZoomOptionsMultiTPC() {OnlyPlaneChanged=-1;}
    ~ZoomOptionsMultiTPC() {}
    std::map<int, double>  wmin;
    std::map<int, double>  wmax;
    std::map<int, double>  tmin;
    std::map<int, double>  tmax;
    int OnlyPlaneChanged;
  };

  
   
  class TWQMultiTPCProjectionView : public evdb::Canvas {

  public:
    
    RQ_OBJECT("evd::TWQMultiTPCProjectionView")

  public:
    TWQMultiTPCProjectionView(TGMainFrame* mf);
    ~TWQMultiTPCProjectionView();

    const char* Description() const { return "Time/Wire/Charge Projections"; }
    const char* PrintTag()    const { return "twq-multiTPC-proj";            }

    void Draw(const char* opt="");

 //   void        RangeChanged();
    static void MouseDispatch(int plane, void * wqpv);
    
    
    void 	ChangeWire(int plane); 
    
    void 	SetMouseZoomRegion(int plane);
    
    void 	SelectPoint(int plane);
    // 	TWQMultiTPCProjectionView::SelectHit() is undefined
    //void 	SelectHit(int plane);
    void        SetPlaneWire();
    void        SetPlane();
    void        SetWire();
    void 	SetDistance();
    void        SetThreshold();
    void        SetGreyscale();
    void        SetMCInfo();
    void        SetRawCalib();
    void 	SetUpSideBar();
    void        SetUpZoomButtons();
    void 	SetUpPositionFind();
    void        SetZoom(int plane,int wirelow,int wirehi,int timelo,int timehi, bool StoreZoom=true); 
    void 	ZoomInterest(bool flag=true);

    void        ZoomBack();  // Revert to the previous zoom setting
        
    void 	SetZoomInterest();
    void 	PrintCharge();
    void 	DrawPads(const char* opt="");
    
    void 	FindEndPoint();
    double 	FindLineLength();
    void 	ClearEndPoints();
    void 	ToggleEndPointMarkers();

    void	RadioButtonsDispatch(int parameter);
  
  private:

    HeaderPad*  fHeaderPad;              ///< Show header information    
    TQPad*      fWireQ;                  ///< Histogram of charge vs time on selected wire
    MCBriefPad* fMC;                     ///< Short summary of MC event    
    std::vector<TQPad* >        fPlaneQ; ///< charge on each plane
    std::vector<TWireProjPad*>  fPlanes; ///< time vs wire projection for each plane
    
    
    TGCompositeFrame*    fVFrame;      ///< needed for the side frame
    TGCompositeFrame*    fMetaFrame;   ///< needed for the side frame

    TGLabel* fWireLabel;
    TGLabel* fDistanceLabel;
    TGLabel* fPlaneLabel;
    TGLabel* fThresLabel;
    TGLabel* fGreyLabel;

    TGNumberEntry* fWireEntry;     ///< Wire number displayed.	
    TGNumberEntry* fPlaneEntry;    ///< Plane number displayed.	
    TGNumberEntry* fThresEntry;    ///< ADC threshold to display.
    TGNumberEntry* fDistance;      ///<Distance from line to find hits in cluster 	
    TGCheckButton* fGreyScale;     ///< Display gray or color scale
    TGCheckButton* fMCOn;          ///< Display MC truth information
    TGRadioButton* fRawDraw;       ///< Draw Raw information only
    TGRadioButton* fCalibDraw;     ///< Draw calibrated information only
    TGRadioButton* fRawCalibDraw;  ///< Draw raw and calibrated information
   
    TGTextButton*  fZoomInterest;   ///< Zoom on iteresting region
    TGTextButton*  fUnZoomInterest; ///< Unzoom on iteresting region
    TGTextButton*  fZoomBack;       ///< Unzoom on iteresting region
    TGCheckButton* fToggleAutoZoom; ///< Toggle the autozoom setting 
    TGRadioButton* fToggleZoom;     ///< Use zoom setting
    TGTextView*    fAngleInfo;
    
    TGTextButton*  fFindEndpoint;      ///< Calculate XYZ position of two points in wire planes
    TGTextButton*  fClearPPoints;      ///< Clear current list of End Points
    TGCheckButton* fToggleShowMarkers; ///< Toggle the ShowEndPointMarkersSetting
    TGTextView*    fXYZPosition;       ///< Display the xyz position
    
    TGTextButton* fCalcAngle; ///< Calculate the 2D & 3D angles between lines
    TGTextButton* fClear;     ///< Clears the selected points in an event
    
    int DrawLine(int plane,util::PxLine &pline);
    
    std::deque<util::PxPoint>  ppoints; ///< list of points in each WireProjPad used for x,y,z finding
    std::deque<util::PxLine>   pline;   ///< list of lines in each WireProjPad used for calculating 2d and 3d angles

    ZoomOptionsMultiTPC fZoomOpt;
    std::vector<ZoomOptionsMultiTPC> fPrevZoomOpt;
    
  };

}// namespace

#endif //EVD_TWQPROJECTION_H
