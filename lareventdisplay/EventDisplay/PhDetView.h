////////////////////////////////////////////////////////////////////////
///
/// \file    PhDetView.h
/// \brief   A view showing the time vs wire, charge and charge vs time information for an event
/// \author  brebel@fnal.gov
// clustering code written by ellen.klein@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef EVD_PHDET_H
#define EVD_PHDET_H

#include "RQ_OBJECT.h"
#include "lardata/Utilities/PxUtils.h"
#include "nuevdb/EventDisplayBase/Canvas.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include <deque>
#include <map>
#include <vector>

// Forward declarations
class TGCheckButton;
class TGCompositeFrame;
class TGLabel;
class TGMainFrame;
class TGNumberEntry;
class TGRadioButton;
class TGTextButton;
class TGTextView;

namespace util {
  class DataProductChangeTracker_t;
}

namespace evd {

  class HeaderPad;
  class MCBriefPad;
  class TQPad;
  class TWireProjPad;

  // Helper class, to store zoom settings in each view

  struct ZoomOptionsPD {
    std::map<int, double> wmin;
    std::map<int, double> wmax;
    std::map<int, double> tmin;
    std::map<int, double> tmax;
    int OnlyPlaneChanged{-1};
  };

  class PhDetView : public evdb::Canvas {
  public:
    RQ_OBJECT("evd::PhDetView")

  public:
    PhDetView(TGMainFrame* mf);
    ~PhDetView();

    const char* Description() const { return "Photon Detector"; }
    const char* PrintTag() const { return "phdet"; }

    void Draw(const char* opt = "");

    //   void        RangeChanged();
    static void MouseDispatch(int plane, void* wqpv);

    void ChangeWire(int plane);
    void SelectTPC(); ///< select TPC from GUI

    void SetMouseZoomRegion(int plane);

    void SelectPoint(int plane);
    void SelectHit(int plane);
    void SetPlaneWire();
    void SetPlane();
    void SetWire();
    void SetDistance();
    void SetThreshold();
    void SetGreyscale();
    void SetMCInfo();
    void SetRawCalib();
    void SetUpSideBar();
    void ForceRedraw(); ///< Forces a redraw of the window
    void SetUpZoomButtons();
    void SetUpClusterButtons();
    void SetUpDrawingButtons();
    void SetUpTPCselection();
    void SetUpPositionFind();
    void SetZoom(int plane, int wirelow, int wirehi, int timelo, int timehi, bool StoreZoom = true);
    void ZoomInterest(bool flag = true);
    /// Clear all the regions of interest
    void ResetRegionsOfInterest();

    void ZoomBack(); // Revert to the previous zoom setting
    void SetClusterInterest();

    void SetClusters(int plane);

    void SetZoomInterest();
    void SetZoomFromView();
    void PrintCharge();
    void DrawPads(const char* opt = "");

    void FindEndPoint();
    double FindLineLength(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp);
    void ClearEndPoints();
    void ToggleEndPointMarkers();

    void RadioButtonsDispatch(int parameter);
    void SaveSelection();
    void ClearSelection();

    /// Returns if a new event is detected; if so, it also resets accordingly
    bool OnNewEvent();

    void SetTestFlag(int number = 1);

  protected:
    TGCompositeFrame* fVFrame; ///< needed for the side frame

  private:
    HeaderPad* fHeaderPad;              ///< Show header information
    TQPad* fWireQ;                      ///< Histogram of charge vs time on selected wire
    MCBriefPad* fMC;                    ///< Short summary of MC event
    std::vector<TQPad*> fPlaneQ;        ///< charge on each plane
    std::vector<TWireProjPad*> fPlanes; ///< time vs wire projection for each plane

    double fLastThreshold; ///< Kludge to prevent double drawing when changing threshold

    TGCompositeFrame* fMetaFrame; ///< needed for the side frame

    TGLabel* fWireLabel;
    TGLabel* fDistanceLabel;
    TGLabel* fPlaneLabel;
    TGLabel* fThresLabel;

    TGNumberEntry* fWireEntry;    ///< Wire number displayed.
    TGNumberEntry* fPlaneEntry;   ///< Plane number displayed.
    TGNumberEntry* fThresEntry;   ///< ADC threshold to display.
    TGNumberEntry* fDistance;     ///< Distance from line to find hits in cluster
    TGCheckButton* fGreyScale;    ///< Display gray or color scale
    TGCheckButton* fMCOn;         ///< Display MC truth information
    TGRadioButton* fRawDraw;      ///< Draw Raw information only
    TGRadioButton* fCalibDraw;    ///< Draw calibrated information only
    TGRadioButton* fRawCalibDraw; ///< Draw raw and calibrated information

    TGTextButton* fZoomInterest;    ///< Zoom on iteresting region
    TGTextButton* fUnZoomInterest;  ///< Unzoom on iteresting region
    TGTextButton* fZoomBack;        ///< Unzoom on iteresting region
    TGCheckButton* fToggleAutoZoom; ///< Toggle the autozoom setting
    TGRadioButton* fToggleZoom;     ///< Use zoom setting
    TGRadioButton* fToggleClusters; ///< Use make cluster setting
    TGTextView* fAngleInfo;

    TGTextButton* fFindEndpoint;       ///< Calculate XYZ position of two points in wire planes
    TGTextButton* fClearPPoints;       ///< Clear current list of End Points
    TGCheckButton* fToggleShowMarkers; ///< Toggle the ShowEndPointMarkersSetting
    TGTextView* fXYZPosition;          ///< Display the xyz position

    TGTextButton* fCalcAngle; ///<Calculate the 2D & 3D angles between lines
    TGTextButton* fClear;     ///<Clears the selected points in an event

    TGTextButton* fRedraw; ///< Button to force redraw

    TGNumberEntry* fCryoInput; ///< current cryostat
    TGNumberEntry* fTPCInput;  ///< current TPC
    TGLabel* fTotalTPCLabel;   ///< total TPCs in the current cryostat

    int DrawLine(int plane, util::PxLine& pline);

    std::deque<util::PxPoint>
      ppoints; ///< list of points in each WireProjPad used for x,y,z finding
    std::deque<util::PxLine>
      pline; ///< list of lines in each WireProjPad used for calculating 2d and 3d angles, also making seeds (eventually)

    std::vector<util::PxLine>
      seedlines; ///< list of lines in each WireProjPad used for calculating 2d and 3d angles, also making seeds (eventually)

    ZoomOptionsPD fZoomOpt;
    std::vector<ZoomOptionsPD> fPrevZoomOpt;

    bool isZoomAutomatic; ///< true if user did not ask for custom zoom

    util::DataProductChangeTracker_t* fLastEvent; ///< keeps track of latest event

    /// Records whether we are automatically zooming to the region of interest
    void SetAutomaticZoomMode(bool bSet = true);

    /// Returns a string visualizing the total number of elements
    static std::string TotalElementsString(unsigned int NElements);
  };

} // namespace

#endif //EVD_PHDET_H
