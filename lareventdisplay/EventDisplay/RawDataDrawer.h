/// \file    RawDataDrawer.h
/// \brief   Class to aid in the rendering of RawData objects
/// \author  messier@indiana.edu
/// \version $Id: RawDataDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
#ifndef EVD_RAWDATADRAWER_H
#define EVD_RAWDATADRAWER_H
#include <vector>
#ifndef __CINT__

#include "CalibrationDBI/Interface/ChannelStatusProvider.h" // lariov::ChannelStatusProvider::Status_t

#endif

class TH1F;
class TVirtualPad;
namespace art    { class Event;     }
namespace evdb   { class View2D;    }
namespace evdb   { class View3D;    }
namespace raw    { class RawDigit;  }

namespace evd {
  
  namespace details {
    class RawDigitCacheClass;
    class CellGridClass;
  } // namespace details

  /// Aid in the rendering of RawData objects
  class RawDataDrawer {
  public:
    RawDataDrawer();
    ~RawDataDrawer();

    void RawDigit2D(const art::Event& evt,
		    evdb::View2D*     view,
		    unsigned int      plane);

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
   
    int GetRegionOfInterest(int plane,
			    int& minw,
			    int& maxw,
			    int& mint,
			    int& maxt);
    
    void GetChargeSum(int plane,
		      double& charge,
		      double& convcharge);
   
  private:
    /// Cache of raw digits
    // Never use raw pointers. Unless you are dealing with CINT, that is.
    evd::details::RawDigitCacheClass* digit_cache;
    
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
    
    // TODO with ROOT 6, turn this into a std::unique_ptr()
    details::CellGridClass* fDrawingRange; ///< information about the viewport
    
    /// Empty collection, used in return value of invalid digits
    static std::vector<raw::RawDigit> const EmptyRawDigits;
    
  }; // class RawDataDrawer
  
  
}

#endif
////////////////////////////////////////////////////////////////////////
