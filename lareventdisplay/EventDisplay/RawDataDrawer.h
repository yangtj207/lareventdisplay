/// \file    RawDataDrawer.h
/// \brief   Class to aid in the rendering of RawData objects
/// \author  messier@indiana.edu
/// \version $Id: RawDataDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
#ifndef EVD_RAWDATADRAWER_H
#define EVD_RAWDATADRAWER_H

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::PlaneID

#include <vector>
#ifndef __CINT__

#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h" // lariov::ChannelStatusProvider::Status_t

#endif

class TH1F;
class TVirtualPad;
namespace art    { class Event;     }
namespace evdb   { class View2D;    }
namespace evdb   { class View3D;    }
namespace raw    { class RawDigit;  }
namespace lar {
  namespace util {
    template <typename T> class MinMaxCollector;
  }
}
namespace util {
  class PlaneDataChangeTracker_t;
} // namespace util

namespace evd {
  
  namespace details {
    class RawDigitCacheDataClass;
    class CellGridClass;
    typedef ::util::PlaneDataChangeTracker_t CacheID_t;
  } // namespace details

  /// Aid in the rendering of RawData objects
  class RawDataDrawer {
  public:
    
    RawDataDrawer();
    ~RawDataDrawer();
    
    /**
     * @brief Draws raw digit content in 2D wire plane representation
     * @param evt source for raw digits
     * @param view target rendered object
     * @param plane number of the plane to be drawn
     * @param bZoomToRoI whether to render only te region of interest
     * 
     * This function performs pre-rendering of the raw digit content into a
     * 2D view of the wire plane as TDC vs. wire number.
     * The material for rendering is created and sent to view object for actual
     * rendering.
     * The pre-rendering result currently depends on information from the
     * current rendering canvas, and in particular on its viewport in the
     * (wire, TDC) space.
     * 
     * If no zoom to the region of interest is required, the region itself is
     * computed (only if not known yet) while the rendering is performed.
     * If the zoom is required instead, rendering is performed in two steps;
     * in the first, run only of no region of interest is known yet, the region
     * is extracted. In the second, that information is used for rendering.
     */
    void RawDigit2D(
      art::Event const& evt, evdb::View2D* view, unsigned int plane,
      bool bZoomToRoI = false
      );

/*     void RawDigit3D(const art::Event& evt, */
/* 		    evdb::View3D*     view); */

    void FillQHisto(const art::Event& evt, 
		    unsigned int      plane,
		    TH1F*             histo);

    void FillTQHisto(const art::Event& evt, 
		     unsigned int      plane,
		     unsigned int      wire,
		     TH1F*             histo);

    double StartTick()       const { return fStartTick; }
    double TotalClockTicks() const { return fTicks; }
    
    /// Fills the viewport information from the specified pad
    void ExtractRange
      (TVirtualPad* pPad, std::vector<double> const* zoom = nullptr);
   
    /// Fills the viewport borders from the specified extremes
    void SetDrawingLimits
      (float low_wire, float high_wire, float low_tdc, float high_tdc);
    
    int GetRegionOfInterest(int plane,
			    int& minw,
			    int& maxw,
			    int& mint,
			    int& maxt);
    
    /// Forgets about the current region of interest
    void ResetRegionOfInterest();
    
    /// Returns whether there is currently a valid region of interest
    /// for the specified plane
    bool hasRegionOfInterest(geo::PlaneID::PlaneID_t plane) const;
    
    void GetChargeSum(int plane,
		      double& charge,
		      double& convcharge);
   
  private:
    typedef struct {
      int    adc   = 0;  ///< total ADC count in this box
      bool   good  = false; ///< whether the channel is not bad
    } BoxInfo_t;
    
    typedef struct {
      unsigned int width = 0; // width of pad in pixels
      unsigned int height = 0; // heigt of pad in pixels
      
      /// Returns whether the stored value is valid
      bool isFilled() const { return (width != 0) && (height != 0); }
      
      /// Returns whether the stored value is valid
      operator bool() const { return isFilled(); }
      
    } PadResolution_t; ///< Stores the information about the drawing area
    
    /// Helper class to be used with ChannelLooper()
    class OperationBaseClass;
    class ManyOperations;
    class BoxDrawer;
    class RoIextractorClass;
    
    // Since this is a private facility, we indulge in non-recommended practises
    // like friendship; these classes have the ability to write their findings
    // directly into RawDataDrawer; the alternative would be to have an
    // interface writing that information exposed to the public, that has the
    // draw back that we completely lose control on who can use it.
    // Other solutions might be also possible...
    // but this is already more complicated than needed.
    friend class BoxDrawer;
    friend class RoIextractorClass;
    
    /// Cache of raw digits
    // Never use raw pointers. Unless you are dealing with CINT, that is.
    evd::details::RawDigitCacheDataClass* digit_cache;
    
#ifndef __CINT__
    /// Prepares for a new event (if somebody tells it to)
    void Reset(art::Event const& event);
    
    /// Reads raw::RawDigits; also triggers Reset()
    void GetRawDigits(art::Event const& evt);
    
    /// Returns whether a channel with the specified status should be processed
    bool ProcessChannelWithStatus
      (lariov::ChannelStatusProvider::Status_t channel_status) const;
#endif // __CINT__

    double fStartTick;                       ///< low tick
    double fTicks;                           ///< number of ticks of the clock
   
    std::vector<int> fWireMin;     ///< lowest wire in interesting region for each plane
    std::vector<int> fWireMax;     ///< highest wire in interesting region for each plane
    std::vector<int> fTimeMin;     ///< lowest time in interesting region for each plane
    std::vector<int> fTimeMax;     ///< highest time in interesting region for each plane
    
    std::vector<double> fRawCharge;     ///< Sum of Raw Charge
    std::vector<double> fConvertedCharge;     ///< Sum of Charge Converted using Birks' formula
    
    PadResolution_t PadResolution; ///< stored pad resolution
    
    details::CacheID_t* fCacheID; ///< information about the last processed plane
    
    // TODO with ROOT 6, turn this into a std::unique_ptr()
    details::CellGridClass* fDrawingRange; ///< information about the viewport
    
    /// Performs the 2D wire plane drawing
    void DrawRawDigit2D
      (art::Event const& evt, evdb::View2D* view, unsigned int plane);
    
    /**
     * @brief Makes sure raw::RawDigit's are available for the current settings
     * @param evt event to read the digits from
     * @param ts a cache ID assessing the new state the cache should move to
     * 
     * The function will ask the data cache for an update
     * (RawDigitCacheDataClass::Update()).
     * The cache will evaluate whether it is already in a state compatible with
     * ts or if cache needs to be invalidated, in which case it will fill with
     * new data.
     * This method also triggers a Reset() if the target state differs from the
     * old one.
     */
    void GetRawDigits
      (art::Event const& evt, details::CacheID_t const& new_timestamp);
    
    // Helper functions for drawing
    bool RunOperation(art::Event const& evt, OperationBaseClass* operation);
    void QueueDrawingBoxes(
      evdb::View2D* view,
      geo::PlaneID const& pid,
      std::vector<BoxInfo_t> const& BoxInfo
      );
    void RunDrawOperation
      (art::Event const& evt, evdb::View2D* view, unsigned int plane);
    void RunRoIextractor(art::Event const& evt, unsigned int plane);
    void SetDrawingLimitsFromRoI(geo::PlaneID::PlaneID_t plane);
    void SetDrawingLimitsFromRoI(geo::PlaneID const pid)
      { SetDrawingLimitsFromRoI(pid.Plane); }
    
    /// Empty collection, used in return value of invalid digits
    static std::vector<raw::RawDigit> const EmptyRawDigits;
    
  }; // class RawDataDrawer
  
  
}

#endif
////////////////////////////////////////////////////////////////////////
