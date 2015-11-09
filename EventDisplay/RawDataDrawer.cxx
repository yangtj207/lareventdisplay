/**
 * @file    RawDataDrawer.cxx
 * @brief   Class to aid in the rendering of RawData objects
 * @author  messier@indiana.edu
 * 
 * This class prepares the rendering of the raw digits content for the 2D view
 * of the display. In particular, it fills the 2D view of detected charge
 * vs. time and wire number for a single TPC.
 * The code is not ready to support an arbitrary TPC from the detector;
 * it can be fixed to support that, but a good deal of the calling code also
 * has to change.
 * 
 * Notes from Gianluca Petrillo (petrillo@fnal.gov) on August 19, 2015
 * --------------------------------------------------------------------
 * 
 * As of August 2015, this class performs some preprocessing for the data to
 * be displayed.
 * The main argument here is that the display assigns to each plane a space
 * of 600x250 pixels, while the amound of data is typically 500x5000 elements.
 * These numbers vary wildly: while the first number is close to the current
 * default size for a new window, that window can be resized and windows many
 * times larger can be obtained; for the data content, a MicroBooNE plane
 * contains 2400x9600 elements (that's roughly 25 millions, and there is three
 * of them), each of the two LArIAT planes contains 240x3200 (less than 1
 * million), and a single TPC of DUNE's 35t prototype about 350x4800 (but it
 * is larger when extended readout window modes are used).
 * The code produces TBox'es to be rendered by the nutools event display
 * infrastructure. The code before August 2015 would produce one box for each
 * of the aggregated charge; charge could be aggregated in time by FHiCL
 * configuration to merge TDC ticks. This typically bloats rendering time,
 * since rendering of all boxes is attempted.
 * The new code performs dynamic aggregation after discovering the actual size
 * of the graphical viewport, and it submits at most one TBox per pixel.
 * Additional improvement is caching of the uncompressed raw data, so that
 * following zooming is faster, and especially a way to bypass the decompression
 * when the original data is not compressed in the first place, that saves
 * a bit of time and quite some memory.
 * 
 * @todo There are probably a number of glitches and shortcuts in the current
 * preprocessing implementation. If they become a problem, they can probably be
 * fixed on demand. Examples include:
 * - no alignment of the boxes with the wire and tick numbers
 * - possible border effects
 * - the fact that only the maximum charge is displayed on each box (an option
 *   is almost implemented to replace that with a sum, and an average is also
 *   possible), and no dynamic charge range detection is in place
 * - the first drawing is performed with a grid that is not the final one
 *   (because for some reason the frame starts larger than it should and is
 *   resized later)
 * - the drawing honours the zoom region the user selected, even if the viewport
 *   ends up being larger (that is, if upstream decides that the zoom region is
 *   too small and a larger area will be drawn, the additional area will be
 *   blank)
 */
#include <cmath> // std::abs(), ...
#include <utility> // std::pair<>, std::move()
#include <memory> // std::unique_ptr()
#include <tuple>
#include <limits> // std::numeric_limits<>
#include <type_traits> // std::add_const_t<>, ...
#include <algorithm> // std::fill(), std::find_if(), ...
#include <cstddef> // std::ptrdiff_t

#include "TH1F.h"
// #include "TPolyLine3D.h"
#include "TBox.h"
#include "TFrame.h"
#include "TVirtualPad.h"

#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/EventHolder.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/StatCollector.h" // lar::util::MinMaxCollector<>
#include "Geometry/Geometry.h"
#include "Utilities/LArPropertiesService.h"
#include "Utilities/IDetectorPropertiesService.h"


#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"


// internal use classes declaration;
// it can't live in the header because it uses C++11/14
namespace details {
  /// @todo Document this code and make it into a library
  template <typename T>
  class PointerToData_t {
      public:
    using value_type = T;
    using const_value_type = std::add_const_t<value_type>;
    using reference = std::add_lvalue_reference_t<value_type>;
    using const_reference
      = std::add_lvalue_reference_t<const_value_type>;
    using pointer = std::add_pointer_t<value_type>;
    using const_pointer = std::add_pointer_t<const_value_type>;
    
    /// Destructor: gets rid of the owned data
    ~PointerToData_t() { Clear(); }
    
    /// @name Dereferencing
    /// @{
    const_reference operator* () const { return *pData; }
    reference operator* () { return *pData; }
    
    const_pointer operator-> () const { return pData; }
    pointer operator-> () { return pData; }
    /// @}
    
    /// Returns whether we point to something
    operator bool() const { return hasData(); }
    
    /// Returns whether we point to nothing
    bool operator! () const { return !hasData(); }
    
    /// Returns whether we have data
    bool hasData() const { return bool(pData); }
    
    /// Returns whether we have data and we own it
    bool owned() const { return bOwned && hasData(); }
    
    /// Sets the data and the ownership
    void SetData(pointer data, bool owned)
      {
        Clear();
        bOwned = owned;
        pData = data;
      }
    /// Acquire ownership of the specified data
    void AcquireData(pointer data) { SetData(data, true); }
    /// Point to the specified data, not acquiring ownership
    void PointToData(pointer data) { SetData(data, false); }
    /// Point to the specified data, not acquiring ownership
    void PointToData(reference data) { SetData(&data, false); }
    /// Move data from the specified object, and own it
    void StealData(std::remove_const_t<T>&& data)
      { AcquireData(new T(std::move(data))); }
    /// Create a owned copy of the specified object
    void NewData(T const& data) { AcquireData(new T(data)); }
    /// Stop pointing to the data; if owned, delete it
    void Clear()
      {
        if (bOwned) delete pData;
        pData = nullptr;
        bOwned = false;
      } // Clear()
    
      protected:
    bool bOwned = false;   ///< whether we own our data
    pointer pData = nullptr; ///< pointer to data
  }; // class PointerToData_t<>
}

namespace evd {
  namespace details {
    
    /// Information about a RawDigit; may contain uncompressed duplicate of data
    class RawDigitInfo_t {
        public:
      /// Returns an art pointer to the actual digit
      art::Ptr<raw::RawDigit> DigitPtr() const { return digit; }
      
      /// Returns an art pointer to the actual digit
      raw::RawDigit const& Digit() const { return *digit; }
      
      /// Returns the channel of this digit (for convenience)
      raw::ChannelID_t Channel() const
        { return digit? digit->Channel(): raw::InvalidChannelID; }
      
      /// minimum charge
      short MinCharge() const { return SampleInfo().min_charge; }
      
      /// maximum charge
      short MaxCharge() const { return SampleInfo().max_charge; }
    
      /// average charge
    //  short AverageCharge() const { return SampleInfo().average_charge; }
      
      /// Returns the uncompressed data
      raw::RawDigit::ADCvector_t const& Data() const;
      
      /// Parses the specified digit
      void Fill(art::Ptr<raw::RawDigit> const& src);
      
      /// Deletes the data
      void Clear();
      
      /// Dumps the content of the digit info
      template <typename Stream>
      void Dump(Stream&& out) const;
      
        private:
      typedef struct {
        short min_charge = std::numeric_limits<short>::max(); ///< minimum charge
        short max_charge = std::numeric_limits<short>::max(); ///< maximum charge
      //  float average_charge = 0.; ///< average charge
      } SampleInfo_t; // SampleInfo_t
      
      art::Ptr<raw::RawDigit> digit; ///< a pointer to the actual digit
      
      /// Uncompressed data
      mutable ::details::PointerToData_t<raw::RawDigit::ADCvector_t const> data;
      
      /// Information collected from the uncompressed data
      mutable std::unique_ptr<SampleInfo_t> sample_info;
      
      /// Fills the uncompressed data cache
      void UncompressData() const;
      
      /// Fills the sample info cache
      void CollectSampleInfo() const;
      
      /// Returns the uncompressed data
      SampleInfo_t const& SampleInfo() const;
      
    }; // class RawDigitInfo_t
    
    
    /// Cached set of RawDigitInfo_t
    class RawDigitCacheClass {
        public:
     
      /// Returns the list of digit info
      std::vector<RawDigitInfo_t> const& Digits() const { return digits; }
      
      /// Returns a pointer to the digit info of given channel, nullptr if none
      RawDigitInfo_t const* FindChannel(raw::ChannelID_t channel) const;

      
      /// Returns the largest number of samples in the unpacked raw digits
      size_t MaxSamples() const { return max_samples; }
      
      /// Returns whether the cache is empty()
      bool empty() const { return digits.empty(); }
      
      /// Updates the cache if needed; returns whether it was needed
      bool Update(art::Event const& evt, art::InputTag const& label);
      /// Empties the cache
      void Clear();
      /// Declare the cache invalid
      void Invalidate();
      
      /// Dump the content of the cache
      template <typename Stream>
      void Dump(Stream&& out) const;
      
        private:
      struct CacheTimestamp_t {
        art::EventID event_id;
        art::InputTag input_label;
        
        void clear()
          { event_id = art::EventID(); input_label = art::InputTag(); }
        
        bool operator== (CacheTimestamp_t const& ts) const
          {
            return (ts.event_id == event_id) && (ts.input_label == input_label);
          } // operator==
        bool operator!= (CacheTimestamp_t const& ts) const
          {
            return (ts.event_id != event_id) || (ts.input_label != input_label);
          } // operator!=
        
        operator std::string() const
          {
            return "R: " + std::to_string(event_id.run())
              + " S: " + std::to_string(event_id.subRun())
              + " E: " + std::to_string(event_id.event())
              + " I: " + input_label.encode();
          }
      }; // CacheTimestamp_t
      
      std::vector<RawDigitInfo_t> digits; ///< vector of raw digit information
      
      size_t max_samples = 0; ///< the largest number of ticks in any digit
      
      CacheTimestamp_t timestamp; ///< time stamp for the cache
      
      /// Fills the cache
      void Refill(art::Event const& evt, CacheTimestamp_t const& ts);
      /// Returns whether the cache is updated to this event
      bool isUpToDate(CacheTimestamp_t const& ts) const;
      
    }; // struct RawDigitCacheClass
    

    std::vector<evd::details::RawDigitInfo_t>::const_iterator begin
      (evd::details::RawDigitCacheClass const& cache)
      { return cache.Digits().cbegin(); }
    std::vector<evd::details::RawDigitInfo_t>::const_iterator end
      (evd::details::RawDigitCacheClass const& cache)
      { return cache.Digits().cend(); }


    
    /// Manages a cell-like division of a coordinate
    class GridAxisClass {
        public:
      /// Default constructor: an invalid range
      GridAxisClass() { Init(0, 0., 0.); }
      
      /// Constructor: sets the limits and the number of cells
      GridAxisClass(size_t nDiv, float new_min, float new_max)
        { Init(nDiv, new_min, new_max); }
      
      //@{
      /// Returns the index of the specified cell
      std::ptrdiff_t GetCell(float coord) const;
      std::ptrdiff_t operator() (float coord) const { return GetCell(coord); }
      //@}
      
      /// Returns whether the cell is present or not
      bool hasCell(std::ptrdiff_t iCell) const
        { return (iCell >= 0) && ((size_t) iCell < NCells()); }
      
      /// Returns whether the coordinate is included in the range or not
      bool hasCoord(float coord) const
        { return (coord >= Min()) && (coord < Max()); }
      
      
      //@{
      /// Returns the extremes of the axis
      float Min() const { return min; }
      float Max() const { return max; }
      //@}
      
      /// Returns the length of the axis
      float Length() const { return max - min; }
      
      /// Returns the length of the axis
      size_t NCells() const { return n_cells; }
      
      /// Returns whether minimum and maximum match
      bool isEmpty() const { return max == min; }
      
      /// Returns the cell size
      float CellSize() const { return cell_size; }
      
      /// Returns the lower edge of the cell
      float LowerEdge(std::ptrdiff_t iCell) const
        { return Min() + CellSize() * iCell; }
      
      /// Returns the upper edge of the cell
      float UpperEdge(std::ptrdiff_t iCell) const
        { return LowerEdge(iCell + 1); }
      
      /// Initialize the axis
      bool Init(size_t nDiv, float new_min, float new_max);
      
      /// Expands the cell (at fixed range) to meet minimum cell size
      /// @return Whether the cell size was changed
      bool SetMinCellSize(float min_size);
      
      /// Expands the cell (at fixed range) to meet maximum cell size
      /// @return Whether the cell size was changed
      bool SetMaxCellSize(float max_size);
      
      /// Expands the cell (at fixed range) to meet maximum cell size
      /// @return Whether the cell size was changed
      bool SetCellSizeBoundary(float min_size, float max_size)
        { return SetMinCellSize(min_size) || SetMaxCellSize(max_size); }
      
      
        private:
      size_t n_cells; ///< number of cells in the axis
      float min, max; ///< extremes of the axis
      
      float cell_size; ///< size of each cell
      
    }; // GridAxisClass
    
    
    /// Manages a grid-like division of 2D space
    class CellGridClass {
        public:
      
      /// Default constructor: invalid ranges
      CellGridClass(): wire_axis(), tdc_axis() {}
      
      /// Constructor: sets the extremes and assumes one cell for each element
      CellGridClass(unsigned int nWires, unsigned int nTDC);
      
      /// Constructor: sets the wire and TDC ranges in detail
      CellGridClass(
        float min_wire, float max_wire, unsigned int nWires,
        float min_tdc, float max_tdc, unsigned int nTDC
        );
      
      /// Returns the total number of cells in the grid
      size_t NCells() const { return wire_axis.NCells() * tdc_axis.NCells(); }
      
      /// Return the information about the wires
      GridAxisClass const& WireAxis() const { return wire_axis; }
      
      /// Return the information about the TDCs
      GridAxisClass const& TDCAxis() const { return tdc_axis; }
      
      
      /// Returns the index of specified cell, or -1 if out of range
      std::ptrdiff_t GetCell(float wire, float tick) const;
      
      /// Returns the coordinates { w1, t1, w2, t2 } of specified cell
      std::tuple<float, float, float, float> GetCellBox
        (std::ptrdiff_t iCell) const;
      
      //@{
      /// Returns whether the range includes the specified wire
      bool hasWire(float wire) const { return wire_axis.hasCoord(wire); }
      bool hasWire(int wire) const { return hasWire((float) wire); }
      //@}
      
      //@{
      /// Returns whether the range includes the specified wire
      bool hasTick(float tick) const { return tdc_axis.hasCoord(tick); }
      bool hasTick(int tick) const { return hasTick((float) tick); }
      //@}
      
      
      /// Increments the specified cell of cont with the value v
      /// @return whether there was such a cell
      template <typename CONT>
      bool Add(CONT& cont, float wire, float tick, typename CONT::value_type v)
        {
          std::ptrdiff_t cell = GetCell(wire, tick);
          if (cell < 0) return false;
          cont[(size_t) cell] += v;
          return true;
        } // Add()
      
      
      /// @name Setters
      /// @{
      /// Sets a simple wire range: all the wires, one cell per wire
      void SetWireRange(unsigned int nWires)
        { SetWireRange(0., (float) nWires, nWires); }
      
      /// Sets the complete wire range
      void SetWireRange(float min_wire, float max_wire, unsigned int nWires)
        { wire_axis.Init(nWires, min_wire, max_wire); }
      
      /// Sets the complete wire range, with minimum cell size
      void SetWireRange
        (float min_wire, float max_wire, unsigned int nWires, float min_size)
        {
          wire_axis.Init(nWires, min_wire, max_wire);
          wire_axis.SetMinCellSize(min_size);
        }
      
      /// Sets a simple TDC range: all the ticks, one cell per tick
      void SetTDCRange(unsigned int nTDC)
        { SetTDCRange(0., (float) nTDC, nTDC); }
      
      /// Sets the complete TDC range
      void SetTDCRange(float min_tdc, float max_tdc, unsigned int nTDC)
        { tdc_axis.Init(nTDC, min_tdc, max_tdc); }
      
      /// Sets the complete TDC range, with minimum cell size
      void SetTDCRange
        (float min_tdc, float max_tdc, unsigned int nTDC, float min_size)
        {
          tdc_axis.Init(nTDC, min_tdc, max_tdc);
          tdc_axis.SetMinCellSize(min_size);
        }
      
      /// @}
      
      /// Sets the minimum size for wire cells
      bool SetMinWireCellSize(float min_size)
        { return wire_axis.SetMinCellSize(min_size); }
      
      /// Sets the minimum size for TDC cells
      bool SetMinTDCCellSize(float min_size)
        { return tdc_axis.SetMinCellSize(min_size); }
      
        protected:
      GridAxisClass wire_axis;
      GridAxisClass tdc_axis;
    }; // CellGridClass
    
    
    //--------------------------------------------------------------------------
    /// Applies Birks correction
    class ADCCorrectorClass {
        public:
      /// Default constructor: awaits for update()
      ADCCorrectorClass() {}
      
      /// Constructor: update()s with the specified information
      ADCCorrectorClass(geo::PlaneID const& pid) { update(pid); }
      
      /// Applies Birks correction to the specified pedestal-subtracted charge
      double Correct(float adc) const
        {
          if (adc < 0.) return 0.;
          register double const dQdX = adc / wirePitch / electronsToADC;
          dataprov::IDetectorProperties const* detp = lar::providerFrom<util::IDetectorPropertiesService>();
          return detp->BirksCorrection(dQdX);
        } // Correct()
      double operator() (float adc) const { return Correct(adc); }
      
      void update(geo::PlaneID const& pid)
        {
          art::ServiceHandle<geo::Geometry> geo;
          wirePitch = geo->WirePitch(pid);

          dataprov::IDetectorProperties const* detp = lar::providerFrom<util::IDetectorPropertiesService>();
          electronsToADC = detp->ElectronsToADC();
        } // update()
      
        protected:
      float wirePitch; ///< wire pitch
      float electronsToADC; ///< conversion constant
      
    }; // ADCCorrectorClass
    //--------------------------------------------------------------------------
  } // namespace details
} // namespace evd

namespace evd {
  
  // empty vector
  std::vector<raw::RawDigit> const RawDataDrawer::EmptyRawDigits;
  
  //......................................................................
  RawDataDrawer::RawDataDrawer()
    : digit_cache(new details::RawDigitCacheClass)
    , fStartTick(0),fTicks(2048)
    , fDrawingRange(new details::CellGridClass)
  { 
    art::ServiceHandle<geo::Geometry> geo;

    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    geo::TPCID tpcid(drawopt->fCryostat, drawopt->fTPC);
    
    fStartTick = drawopt->fStartTick;
    fTicks = drawopt->fTicks;

    // set the list of bad channels in this detector
    unsigned int nplanes=geo->Nplanes(tpcid);
    fWireMin.resize(nplanes,-1);   
    fWireMax.resize(nplanes,-1);    
    fTimeMin.resize(nplanes,-1);    
    fTimeMax.resize(nplanes,-1);    
    fRawCharge.resize(nplanes,0);   
    fConvertedCharge.resize(nplanes,0);
  }

  //......................................................................
  RawDataDrawer::~RawDataDrawer()
  {
    delete digit_cache;
    delete fDrawingRange;
  }


  //......................................................................
  void RawDataDrawer::ExtractRange
    (TVirtualPad* pPad, std::vector<double> const* zoom /* = nullptr */)
  {
    mf::LogDebug log("RawDataDrawer");
    log << "ExtractRange() on pad '" << pPad->GetName() << "'";
    
    TFrame const* pFrame = pPad->GetFrame();
    if (pFrame) {
      // these coordinates are used to find the actual extent of pad in pixels
      double low_wire = pFrame->GetX1(), high_wire = pFrame->GetX2();
      double low_tdc = pFrame->GetY1(), high_tdc = pFrame->GetY2();
      double const wire_pixels = pPad->XtoAbsPixel(high_wire) - pPad->XtoAbsPixel(low_wire);
      double const tdc_pixels = -(pPad->YtoAbsPixel(high_tdc) - pPad->YtoAbsPixel(low_tdc));
      
      log << "\n frame window is " << wire_pixels << "x" << tdc_pixels << " pixel big and";
      // those coordinates also are a (unreliable) estimation of the zoom;
      // if we have a better one, let's use it
      // (this does not change the size of the window in terms of pixels)
      if (zoom) {
        log << ", from external source,";
        low_wire = (*zoom)[0];
        high_wire = (*zoom)[1];
        low_tdc = (*zoom)[2];
        high_tdc = (*zoom)[3];
      }
      
      log << " spans wires "
        << low_wire << "-" << high_wire << " and TDC " << low_tdc << "-" << high_tdc;
      
      fDrawingRange->SetWireRange(low_wire, high_wire, wire_pixels, 1.0);
      fDrawingRange->SetTDCRange(low_tdc, high_tdc, tdc_pixels, 1.0);
    }
    else {
      // keep the old frame (if any)
      log << "\n  no frame!";
    }
    
  } // RawDataDrawer::ExtractRange()
  
  
  //......................................................................
  void RawDataDrawer::RawDigit2D(const art::Event& evt,
				                 evdb::View2D*     view,
				                 unsigned int      plane)
  {
    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires == 1) return;

    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    //get pedestal conditions
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();  
    
    geo::PlaneID pid(drawopt->fCryostat, drawopt->fTPC, plane);

    fRawCharge[plane]       = 0.;
    fConvertedCharge[plane] = 0.;
    
    // set the minimum cell in ticks to at least match fTicksPerPoint
    fDrawingRange->SetMinTDCCellSize((float) drawopt->fTicksPerPoint);

    GetRawDigits(evt);

    if(digit_cache->empty()) return;
    
    // TODO turn this into a user option
    constexpr bool bDrawSum = false;
    
    typedef struct {
      int    adc   = 0;  ///< total ADC count in this box
      bool   good  = false; ///< whether the channel is not bad
    } BoxInfo_t;
    
    
    // set up the size of the grid to be visualized;
    // the information on the size has to be already there:
    // caller should have user ExtractRange() first.
    std::vector<BoxInfo_t> BoxInfo(fDrawingRange->NCells());
    
    lariov::IChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();
    
    // collect the interesting range
    lar::util::MinMaxCollector<float> WireRange, TDCrange;
    
    // loop over all the channels/raw digits
    for (evd::details::RawDigitInfo_t const& digit_info: *digit_cache) {
      raw::RawDigit const& hit = digit_info.Digit();
      raw::ChannelID_t const channel = hit.Channel();
      
      // skip the bad channels
      if (!channelStatus.IsPresent(channel)) continue;
      // The following test is meant to be temporary until the "correct" solution is implemented
      if (!ProcessChannelWithStatus(channelStatus.Status(channel))) continue;
      
      // we have a list of all channels, but we are drawing only on one plane;
      // most of the channels will not contribute to this plane,
      // and before we start querying databases, unpacking data etc.
      // we want to know it's for something
      
      std::vector<geo::WireID> WireIDs = geo->ChannelToWire(channel);
      
      bool bDrawChannel = false;
      for (geo::WireID const& wireID: WireIDs){
        if (wireID.planeID() != pid) continue; // not us!
        bDrawChannel = true;
        break;
      } // for wires
      if (!bDrawChannel) continue;
      
      // collect bad channels
      bool const bGood = !channelStatus.IsBad(channel);
      
      // nothing else to be done if the channel is not good:
      // cells are marked bad by default and if any good channel falls in any of
      // them, they become good
      if (!bGood) continue;
      
      // at this point we know we have to process this channel
      raw::RawDigit::ADCvector_t const& uncompressed = digit_info.Data();
      
      details::ADCCorrectorClass ADCCorrector(pid);
      
      // recover the pedestal
      float const pedestal = pedestalRetrievalAlg.PedMean(channel);
      
      // loop over all the wires that are covered by this channel;
      // without knowing better, we have to draw into all of them
      for (geo::WireID const& wireID: WireIDs){
        // check that the plane and tpc are the correct ones to draw
        if (wireID.planeID() != pid) continue; // not us!
        
        float const wire = float(wireID.Wire);
        
        // if this wire number is out of range, we ignore it
        if (!fDrawingRange->hasWire(wire)) continue;
        
        // get an iterator over the adc values
        // accumulate all the data of this wire in our "cells"
        // TODO we could extract the number of ticks falling in the current cell
        // rather than checking all the ticks in sequence
        // TODO support back ticksPerPoint (automatic accumulation of ticks)
        size_t const max_tick = std::min({
          size_t(fDrawingRange->TDCAxis().Max()),
          uncompressed.size(),
          size_t(fTicks)
          });
        
        for (size_t iTick = fStartTick; iTick < max_tick; ++iTick) {
          
          // check if we are out of range
          const float tick = float(iTick);
          
          std::ptrdiff_t cell = fDrawingRange->GetCell(wire, tick);
          if (cell < 0) continue;
          
          float const adc = uncompressed[iTick] - pedestal;
          
          BoxInfo_t& info = BoxInfo[cell];
          info.good = true; // if in range, we mark this cell as good
          
          fRawCharge[plane] += adc;
          fConvertedCharge[plane] += ADCCorrector(adc);
          
          if (bDrawSum) { // draw total charge in cell
            info.adc += adc;
          }
          else { // draw maximum digit in the cell
            if (std::abs(info.adc) > std::abs(adc)) continue;
            info.adc = adc;
          }
          
          // this range is only within the currently selected region
          WireRange.add(wire);
          TDCrange.add(tick);
          
        } // if good
      } // for wires
    } // for channels
    
    //
    // All the information is now collected in BoxInfo.
    // Make boxes out of it.
    //
    
    // drawing options:
    float const MinSignal = drawopt->fMinSignal;
    bool const bScaleDigitsByCharge = drawopt->fScaleDigitsByCharge;

    geo::SigType_t const sigType = geo->SignalType(pid);
    evdb::ColorScale const& ColorSet = cst->RawQ(sigType);
    size_t const nBoxes = BoxInfo.size();
    unsigned int nDrawnBoxes = 0;
    for (size_t iBox = 0; iBox < nBoxes; ++iBox) {
      BoxInfo_t const& info = BoxInfo[iBox];
      
      // too little signal, don't bother drawing
      if(info.good && (std::abs(info.adc) < MinSignal)) continue;
      
      // skip the bad cells
      if (!info.good) continue;
      
      // box color, proportional to the ADC count
      int const color = ColorSet.GetColor(info.adc);
      
      // scale factor, proportional to ADC count (optional)
      constexpr float q0 = 1000.;
      float const sf = bScaleDigitsByCharge
        ? std::min(std::sqrt((float) info.adc / q0), 1.0F)
        : 1.;
      
      // coordinates of the cell box
      float min_wire, max_wire, min_tick, max_tick;
      std::tie(min_wire, min_tick, max_wire, max_tick)
        = fDrawingRange->GetCellBox(iBox);
      
      if (sf != 1.) { // need to shrink the box
        float const nsf = 1. - sf; // negation of scale factor
        float const half_box_wires = (max_wire - min_wire) / 2.,
          half_box_ticks = (max_tick - min_tick) / 2.;
        
        // shrink the box:
        min_wire += nsf * half_box_wires;
        max_wire -= nsf * half_box_wires;
        min_tick += nsf * half_box_ticks;
        max_tick -= nsf * half_box_ticks;
      } // if scaling
      
      // allocate the box on the view;
      // the order of the coordinates depends on the orientation
      TBox* pBox;
      if (drawopt->fAxisOrientation < 1)
        pBox = &(view->AddBox(min_wire, min_tick, max_wire, max_tick));
      else
        pBox = &(view->AddBox(min_tick, min_wire, max_tick, max_wire));
      
      pBox->SetFillStyle(1001);
      pBox->SetFillColor(color);
      pBox->SetBit(kCannotPick);
      
      ++nDrawnBoxes;
    } // for (iBox)
    
    // now, when do we want to fill the interesting region?
    // this region as computed in the current implementation is within
    // the selected region for visualization (fDrawingRange);
    // so we update it only if it has not been initialized yet
    // (that is, it's empty)
    if ((fWireMin[plane] == fWireMax[plane]) && WireRange.has_data()) {
      mf::LogInfo("RawDataDrawer") << "Region of interest for "
        << std::string(pid) << " detected to be within wires "
        << WireRange.min() << " and " << WireRange.max();
      fWireMax[plane] = WireRange.max() + 1;
      fWireMin[plane] = WireRange.min();
    }
    if ((fTimeMin[plane] == fTimeMax[plane]) && TDCrange.has_data()) {
      mf::LogInfo("RawDataDrawer") << "Region of interest for "
        << std::string(pid) << " detected to be within ticks "
        << TDCrange.min() << " and " << TDCrange.max();
      fTimeMax[plane] = TDCrange.max() + 1;
      fTimeMin[plane] = TDCrange.min();
    }
    
    
#if 0
    const dataprov::IDetectorProperties* detp = art::ServiceHandle<util::IDetectorPropertiesService>()->provider();
    
    for (evd::details::RawDigitInfo_t& digit_info: digit_cache->digits) {
      raw::RawDigit const& hit = digit_info.Digit();
      raw::ChannelID_t const channel = hit.Channel();
      
      // The following test is meant to be temporary until the "correct" solution is implemented
      auto const channel_status = channelFilter.GetChannelStatus(channel);
      if (channel_status == filter::ChannelFilter::NOTPHYSICAL) continue;
      if (channel_status >  drawopt->fMaxChannelStatus)         continue;
        
      geo::SigType_t sigType = geo->SignalType(channel);
      geo::View_t    v       = geo->View(channel);
      double const wirePitch = geo->WirePitch(v); // FIXME this assumes all TPC are the same
      std::vector<geo::WireID> wireids = geo->ChannelToWire(channel);
      bool skipchan = true;
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid){
	  skipchan = false;
	}
      }
      if (skipchan) continue;
      
      raw::RawDigit::ADCvector_t const& uncompressed = digit_info.Data();
      
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid){
        
          // recover the pedestal
          float pedestal = pedestalRetrievalAlg.PedMean(channel);
	  
	  double wire = 1.*wid.Wire;
	  double tick = 0;
	  // get an iterator over the adc values
	  std::vector<short>::const_iterator itr = uncompressed.cbegin(),
	    dend = uncompressed.end();
	  while( itr != dend ){
	    int ticksUsed = 0;
	    double tdcsum = 0.;
	    double adcsum = 0.;
	    while(ticksUsed < ticksPerPoint && itr != dend){
	      tdcsum  += tick;
//          adcsum  += (1.*(*itr)) - fPedestalVec[channel]; // pedestals[plane]; //hit.GetPedestal();
          adcsum  += (1.*(*itr)) - pedestal;
	      ++ticksUsed;
	      tick += 1.;
	      itr++; // this advance of the iterator is sufficient for the external loop too
	    }
	    double adc = adcsum/ticksPerPoint;
	    double tdc = tdcsum/ticksPerPoint;
	    
	    if(std::abs(adc) < drawopt->fMinSignal) continue;
	
	    fRawCharge[plane] += std::abs(adc);

	    double const dQdX = std::abs(adc)/wirePitch/detp->ElectronsToADC();
	    fConvertedCharge[plane] += detp->BirksCorrection(dQdX);
	
	    int    co = 0;
	    double sf = 1.;
	    constexpr double q0 = 1000.0;
	    
	    co = cst->RawQ(sigType).GetColor(adc);
	    if (drawopt->fScaleDigitsByCharge) {
	      sf = std::sqrt(adc/q0);
	      if (sf>1.0) sf = 1.0;
	    }
	    if(wire < minw)
	      minw = wire;
	    if(wire > maxw)
	      maxw = wire;
	    if(tdc < mint)
	      mint = tdc;
	    if(tdc > maxt)
	      maxt = tdc;
	
	    // don't draw boxes for tdc values that don't exist
	    if(tdc > fStartTick+fTicks || tdc < fStartTick) continue;
	    
	    /* FIXME need to find which is the right tdc index
	    BoxInfo_t& boxInfo = BoxInfo[wire * max_ticks + tdc];
	    boxInfo.color = co;
	    boxInfo.scale = sf;
	    */
	    
	    TBox* b1;
	    if(drawopt->fAxisOrientation < 1){
	      b1 = &(view->AddBox(wire-sf*0.5,
				      tdc-sf*0.5*ticksPerPoint,
				      wire+sf*0.5,
				      tdc+sf*0.5*ticksPerPoint));
	      ++nBoxes;
	    }
	    else{
	      b1 = &(view->AddBox(tdc-sf*0.5*ticksPerPoint,
				      wire-sf*0.5,
				      tdc+sf*0.5*ticksPerPoint,
				      wire+sf*0.5));
	      ++nBoxes;
	    }
	    b1->SetFillStyle(1001);
	    b1->SetFillColor(co);    
	    b1->SetBit(kCannotPick);
	    
	    // 	  TBox& b2 = view->AddBox(wire-0.1,tdc-0.1,wire+0.1,tdc+0.1);
	    // 	  b2.SetFillStyle(0);
	    // 	  b2.SetLineColor(15);
	    // 	  b2.SetBit(kCannotPick);
	 
	  }// end loop over samples 
	}// end if in the right plane
      }// end loopo over wireids
    }//end loop over raw hits

    fWireMin[plane] = minw;   
    fWireMax[plane] = maxw;    
    fTimeMin[plane] = mint;    
    fTimeMax[plane] = maxt; 
    
    // now loop over all the bad channels and set them to 0 adc
    for(size_t bc = 0; bc < fBadChannels.size(); ++bc){
      
      geo::SigType_t sigType = geo->SignalType(fBadChannels[bc]);
	
      std::vector<geo::WireID> wireids = geo->ChannelToWire(fBadChannels[bc]);
      
      // check this is the correct plane and tpc
      for( auto const& wid : wireids){
	if(wid.planeID() == pid){
	
	  if(drawopt->fMinSignal > 0) continue;
	
	  int      co = cst->RawQ(sigType).GetColor(0);
	  double wire = 1.*w;
	  
	  for(int iTick = fStartTick; iTick < fStartTick+fTicks;
            iTick += ticksPerPoint)
          {
	    double const tdc = iTick + 0.5*ticksPerPoint;
	  /* FIXME
	    BoxInfo_t& boxInfo = BoxInfo[wire * max_ticks + iTick];
	    boxInfo.color = co;
	    boxInfo.scale = 1.;
	  */
	    if(drawopt->fAxisOrientation < 1){
	      TBox& b1 = view->AddBox(wire-0.5,tdc-0.5*ticksPerPoint,wire+0.5,tdc+0.5*ticksPerPoint);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }
	    else{
	      TBox &b1 = view->AddBox(tdc-0.5*ticksPerPoint,wire-0.5,tdc+0.5*ticksPerPoint,wire+0.5);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }	  
	  }
	}// end if in the right plane
      }// end loop over wireids
    }// end loop over bad channels    
    
#endif // 0
    
    mf::LogDebug("RawDataDrawer") << __func__ << "(" << std::string(pid)
      << ") added " << nDrawnBoxes
      << " " << fDrawingRange->WireAxis().CellSize()
      << "x" << fDrawingRange->TDCAxis().CellSize()
      << " boxes (out of a grid of " << fDrawingRange->WireAxis().NCells()
      << "x" << fDrawingRange->TDCAxis().NCells() << " cells) covering ( "
      << fDrawingRange->WireAxis().Min() << " -- " << fDrawingRange->WireAxis().Max()
      << " ) wire and ( "
      << fDrawingRange->TDCAxis().Min() << " -- " << fDrawingRange->TDCAxis().Max()
      << " ) tick ranges";
    
  } // RawDataDrawer::RawDigit2D()

  //........................................................................
  int RawDataDrawer::GetRegionOfInterest(int plane,int& minw,int& maxw,int& mint,int& maxt)
  {
    art::ServiceHandle<geo::Geometry> geo;
 
    if((unsigned int)plane>fWireMin.size())
      {mf::LogWarning  ("RawDataDrawer") << " Requested plane " << plane <<" is larger than those available " << std::endl;
	return -1;
      }
  
    minw=fWireMin[plane];
    maxw=fWireMax[plane];
    mint=fTimeMin[plane];
    maxt=fTimeMax[plane];
  
    //make values a bit larger, but make sure they don't go out of bounds 
    minw= (minw-30<0) ? 0 : minw-30;
    mint= (mint-10<0) ? 0 : mint-10;

    maxw= (maxw+10>(int)geo->Nwires(plane)) ? geo->Nwires(plane) : maxw+10;
    maxt= (maxt+10>TotalClockTicks()) ? TotalClockTicks() : maxt+10;
    
    return 0;
  }

  //......................................................................
  void RawDataDrawer::GetChargeSum(int plane,double& charge,double& convcharge)
  {
    charge=fRawCharge[plane]; 
    convcharge=fConvertedCharge[plane];    
    
  }

  //......................................................................
  void RawDataDrawer::FillQHisto(const art::Event& evt,
				                 unsigned int      plane,
				                 TH1F*             histo)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires==1) return;

    art::ServiceHandle<geo::Geometry> geo;
  
    GetRawDigits(evt);
    
    geo::PlaneID pid(drawopt->fCryostat, drawopt->fTPC, plane);
      
    lariov::IChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();
    
    //get pedestal conditions
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();

    for (evd::details::RawDigitInfo_t const& digit_info: *digit_cache) {
      raw::RawDigit const& hit = digit_info.Digit();
      raw::ChannelID_t const channel = hit.Channel();
        
      if (!channelStatus.IsPresent(channel)) continue;
      
      // The following test is meant to be temporary until the "correct" solution is implemented
      if (!ProcessChannelWithStatus(channelStatus.Status(channel))) continue;
      
      // to be explicit: we don't cound bad channels in
      if (channelStatus.IsBad(channel)) continue;
      
      std::vector<geo::WireID> wireids = geo->ChannelToWire(channel);
      for(auto const& wid : wireids){
        // check that the plane and tpc are the correct ones to draw
        if (wid.planeID() != pid) continue;

        raw::RawDigit::ADCvector_t const& uncompressed = digit_info.Data();
        
        float const pedestal = pedestalRetrievalAlg.PedMean(channel);

        for(short d: uncompressed)
          histo->Fill(float(d) - pedestal); //pedestals[plane]); //hit.GetPedestal());
        
        // this channel is on the correct plane, don't double count the raw signal
        // if there are more than one wids for the channel
        break;
      }// end loop over wids
    }//end loop over raw hits

  }

  //......................................................................
  void RawDataDrawer::FillTQHisto(const art::Event& evt,
				  unsigned int      plane,
				  unsigned int      wire,
				  TH1F*             histo)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires==1) return;
    
    // make sure we have the raw digits cached
    GetRawDigits(evt);

    if (digit_cache->empty()) return;

    geo::WireID const wireid(drawopt->fCryostat, drawopt->fTPC, plane, wire);
    
    // find the channel
    art::ServiceHandle<geo::Geometry> geom;
    raw::ChannelID_t const channel = geom->PlaneWireToChannel(wireid);
    if (!raw::isValidChannelID(channel)) { // no channel, empty histogram
      mf::LogError("RawDataDrawer") << __func__ << ": no channel associated to "
        << std::string(wireid);
      return;
    } // if no channel
    
    // check the channel status; bad channels are still ok.
    lariov::IChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();
    
    if (!channelStatus.IsPresent(channel)) return;
    
    // The following test is meant to be temporary until the "correct" solution is implemented
    if (!ProcessChannelWithStatus(channelStatus.Status(channel))) return;
    
    
    // we accept to see the content of a bad channel, so this is commented out:
    // if (channelStatus.IsBad()) return;
    
    //get pedestal conditions
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();
    
    // find the raw digit
    // (iDigit is an iterator to a evd::details::RawDigitInfo_t)
    evd::details::RawDigitInfo_t const* pDigit
      = digit_cache->FindChannel(channel);
    if (!pDigit) { // this is weird...
      mf::LogWarning("RawDataDrawer") << __func__
        << ": can't find raw digit for channel #" << channel
        << " (" << std::string(wireid) << ")";
      return;
    }
    
    raw::RawDigit::ADCvector_t const& uncompressed = pDigit->Data();
    
    float const pedestal = pedestalRetrievalAlg.PedMean(channel);
    
    for(size_t j = 0; j < uncompressed.size(); ++j)
      histo->Fill(float(j), float(uncompressed[j]) - pedestal); //pedestals[plane]); //hit.GetPedestal());
    
  } // RawDataDrawer::FillTQHisto()

  //......................................................................

  //   void RawDataDrawer::RawDigit3D(const art::Event& evt,
  // 				 evdb::View3D*     view)
  //   {
  //     // Check if we're supposed to draw raw hits at all
  //     art::ServiceHandle<evd::RawDrawingOptions> drawopt;
  //     if (drawopt->fDrawRawOrCalibHits!=0) return;


  //     art::ServiceHandle<geo::Geometry> geom;

  //     HitTower tower;
  //     tower.fQscale = 0.01;

  //     for (unsigned int imod=0; imod<drawopt->fRawDigitModules.size(); ++imod) {
  //       const char* which = drawopt->fRawDigitModules[imod].c_str();

  //       std::vector<raw::RawDigit> rawhits;
  //       GetRawDigits(evt, which, rawhits);

  //       for (unsigned int i=0; i<rawhits.size(); ++i) {
  // 	double t = 0;
  // 	double q = 0;
  // 	t = rawhits[i]->fTDC[0];
  // 	for (unsigned int j=0; j<rawhits[i]->NADC(); ++j) {
  // 	  q += rawhits[i]->ADC(j);
  // 	}
  // 	// Hack for now...
  // 	if (q<=0.0) q = 1+i%10;
      
  // 	// Get the cell geometry for the hit
  // 	int         iplane = cmap->GetPlane(rawhits[i].get());
  // 	int         icell  = cmap->GetCell(rawhits[i].get());
  // 	double      xyz[3];
  // 	double      dpos[3];
  // 	geo::View_t v;
  // 	geom->CellInfo(iplane, icell, &v, xyz, dpos);
      
  // 	switch (drawopt->fRawDigit3DStyle) {
  // 	case 1:
  // 	  //
  // 	  // Render digits as towers
  // 	  //
  // 	  if (v==geo::kX) {
  // 	    tower.AddHit(v, iplane, icell, xyz[0], xyz[2], q,  t);
  // 	  }
  // 	  else if (v==geo::kY) {
  // 	    tower.AddHit(v, iplane, icell, xyz[1], xyz[2], q, t);
  // 	  }
  // 	  else abort();
  // 	  break;
  // 	default:
  // 	  //
  // 	  // Render Digits as boxes 
  // 	  //
  // 	  TPolyLine3D& p = view->AddPolyLine3D(5,kGreen+1,1,2);
  // 	  double sf = std::sqrt(0.01*q);
  // 	  if (v==geo::kX) {
  // 	    double x1 = xyz[0] - sf*dpos[0];
  // 	    double x2 = xyz[0] + sf*dpos[0];
  // 	    double z1 = xyz[2] - sf*dpos[2];
  // 	    double z2 = xyz[2] + sf*dpos[2];
  // 	    p.SetPoint(0, x1, geom->DetHalfHeight(), z1);
  // 	    p.SetPoint(1, x2, geom->DetHalfHeight(), z1);
  // 	    p.SetPoint(2, x2, geom->DetHalfHeight(), z2);
  // 	    p.SetPoint(3, x1, geom->DetHalfHeight(), z2);
  // 	    p.SetPoint(4, x1, geom->DetHalfHeight(), z1);
  // 	  }
  // 	  else if (v==geo::kY) {
  // 	    double y1 = xyz[1] - sf*dpos[1];
  // 	    double y2 = xyz[1] + sf*dpos[1];
  // 	    double z1 = xyz[2] - sf*dpos[2];
  // 	    double z2 = xyz[2] + sf*dpos[2];
  // 	    p.SetPoint(0, geom->DetHalfWidth(), y1, z1);
  // 	    p.SetPoint(1, geom->DetHalfWidth(), y2, z1);
  // 	    p.SetPoint(2, geom->DetHalfWidth(), y2, z2);
  // 	    p.SetPoint(3, geom->DetHalfWidth(), y1, z2);
  // 	    p.SetPoint(4, geom->DetHalfWidth(), y1, z1);
  // 	  }
  // 	  else abort();
  // 	  break;
  // 	} // switch fRawDigit3DStyle    
  //       }//end loop over raw digits
  //     }// end loop over RawDigit modules
  
  //     // Render the towers for that style choice
  //     if (drawopt->fRawDigit3DStyle==1) tower.Draw(view);
  //   }

  //......................................................................    

  void RawDataDrawer::Reset(art::Event const& event) {
    // Prepares for a new event.
    
    // Reset the cache of the raw digits;
    // This looks very much circular, since in the current implementation
    // it is precisely GetRawDigits() that gives the trigger for Reset().
    // In practise the following operation is no-op, since the raw digit cache
    // is by now already aligned with the event, so no update is performed
    // and Reset() is not triggered again.
    // This call is here for the future, if the trigger gets moved from
    // GetRawDigits() to somewhere else, raw digits will be properly updated
    // anyway. Not that necessary anyway.
  //  GetRawDigits(event);
    
    // reset the region of interest
    std::fill(fWireMin.begin(), fWireMin.end(), -1);
    std::fill(fWireMax.begin(), fWireMax.end(), -1);
    std::fill(fTimeMin.begin(), fTimeMin.end(), -1);
    std::fill(fTimeMax.begin(), fTimeMax.end(), -1);
    
  } // RawDataDrawer::GetRawDigits()
  
  
  //......................................................................    

  void RawDataDrawer::GetRawDigits(art::Event const& evt) {
  //  digit_cache->Dump(mf::LogDebug("RawDataDrawer"));
    // update cache
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if (digit_cache->Update(evt, rawopt->fRawDataLabel)) Reset(evt);
  } // RawDataDrawer::GetRawDigits()
  
  
  //......................................................................    
  bool RawDataDrawer::ProcessChannelWithStatus
    (lariov::IChannelStatusProvider::Status_t channel_status) const
  {
    // if we don't have a valid status, we can't reject the channel
    if (!lariov::IChannelStatusProvider::IsValidStatus(channel_status))
      return true;
    
    // is the status "too bad"?
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (channel_status > drawopt->fMaxChannelStatus) return false;
    if (channel_status < drawopt->fMinChannelStatus) return false;
    
    // no reason to reject it...
    return true;
  } // RawDataDrawer::ProcessChannel()
  
  //----------------------------------------------------------------------------
  namespace details {
    
    //--------------------------------------------------------------------------
    //--- RawDigitInfo_t
    //---
    raw::RawDigit::ADCvector_t const& RawDigitInfo_t::Data() const {
      if (!data.hasData()) UncompressData();
      return *data;
    } // RawDigitInfo_t::Data()
    
    
    void RawDigitInfo_t::Fill(art::Ptr<raw::RawDigit> const& src) {
      data.Clear();
      digit = src;
    } // RawDigitInfo_t::Fill()
    
    
    void RawDigitInfo_t::Clear() {
      data.Clear();
      sample_info.reset();
    }
    
    
    void RawDigitInfo_t::UncompressData() const {
      data.Clear();
      
      if (!digit) return; // no original data, can't do anything
      
      if (digit->Compression() == kNone) {
        // no compression, we can refer to the original data directly
        data.PointToData(digit->ADCs());
      }
      else {
        // data is compressed, need to do the real work
        raw::RawDigit::ADCvector_t samples;
        Uncompress(digit->ADCs(), samples, digit->Compression());
        data.StealData(std::move(samples));
      }
    } // RawDigitInfo_t::UncompressData()
    
    
    void RawDigitInfo_t::CollectSampleInfo() const {
      raw::RawDigit::ADCvector_t const& samples = Data();
      
    //  lar::util::StatCollector<double> stat;
    //  stat.add(samples.begin(), samples.end());
      
      lar::util::MinMaxCollector<raw::RawDigit::ADCvector_t::value_type> stat
        (samples.begin(), samples.end());
      
      sample_info.reset(new SampleInfo_t);
      sample_info->min_charge = stat.min();
      sample_info->max_charge = stat.max();
    //  sample_info->average_charge = stat.Average();
      
    } // RawDigitInfo_t::CollectSampleInfo()
    
    
    RawDigitInfo_t::SampleInfo_t const& RawDigitInfo_t::SampleInfo() const {
      if (!sample_info) CollectSampleInfo();
      return *sample_info;
    } // SampleInfo()
    
    template <typename Stream>
    void RawDigitInfo_t::Dump(Stream&& out) const {
      out << "  digit at " << ((void*) digit.get()) << " on channel #" << digit->Channel()
        << " with " << digit->NADC();
      if (digit->Compression() == kNone) out << " uncompressed data";
      else out << " data items compressed with <" << digit->Compression() << ">";
      if (data.hasData())
        out << " with data (" << data->size() << " samples)";
      else out << " without data";
    } // RawDigitInfo_t::Dump()
    
    //--------------------------------------------------------------------------
    //--- RawDigitCacheClass
    //---
   
    RawDigitInfo_t const* RawDigitCacheClass::FindChannel
      (raw::ChannelID_t channel) const
    {
      auto iDigit = std::find_if(
        digits.cbegin(), digits.cend(), 
        [channel](evd::details::RawDigitInfo_t const& digit)
          { return digit.Channel() == channel; }
        );
      return (iDigit == digits.cend())? nullptr: &*iDigit;
    } // RawDigitCacheClass::FindChannel()
    
    bool RawDigitCacheClass::Update
      (art::Event const& evt, art::InputTag const& label)
    {
      CacheTimestamp_t ts = { evt.id(), label };
      if (isUpToDate(ts)) return false;
      Refill(evt, ts);
      return true;
    } // RawDigitCacheClass::Update()
    
    void RawDigitCacheClass::Refill
      (art::Event const& evt, CacheTimestamp_t const& ts)
    {
      Clear();
      
      art::Handle< std::vector<raw::RawDigit>> rdcol;
      if (!evt.getByLabel(ts.input_label, rdcol)) {
        mf::LogWarning("RawDataDrawer")
          << "no RawDigit collection '" << ts.input_label << "' found";
        return;
      }
      
      digits.resize(rdcol->size());
      for(size_t iDigit = 0; iDigit < rdcol->size(); ++iDigit) {
        art::Ptr<raw::RawDigit> pDigit(rdcol, iDigit);
        digits[iDigit].Fill(pDigit);
        size_t samples = pDigit->Samples();
        if (samples > max_samples) max_samples = samples;
      } // for
      
      timestamp = ts;
      
    } // RawDigitCacheClass::Refill()
    
    
    void RawDigitCacheClass::Invalidate() {
      timestamp.clear();
    } // RawDigitCacheClass::Invalidate()
    
    
    void RawDigitCacheClass::Clear() {
      Invalidate();
      digits.clear();
      max_samples = 0;
    } // RawDigitCacheClass::Clear()
    
    bool RawDigitCacheClass::isUpToDate(CacheTimestamp_t const& ts) const {
      /*
        I am disabling the cache.
        The reason is that I can't reliably determine whether it became invalid.
        My method was to check whether the event ID or the input tag have changed.
        Fine, except that there are situations where the event display triggers
        a complete reread of the event, and therefore the event data (including
        RawDigit) is deleted. Event the art pointers become "corrupted", because
        they seem to cache the pointed address.
        And how can I know that happened? the event ID is not changed, just all
        the event content is gone.
        Need to find out if EventDisplayBase, in its infinite wisdom, has a mean
        to tell me that this slaughter happened.
      */
      return false;
    // return timestamp == ts;
    } // RawDigitCacheClass::isUpToDate()
    
    template <typename Stream>
    void RawDigitCacheClass::Dump(Stream&& out) const {
      out << "Cache at " << ((void*) this)
        << " with time stamp " << std::string(timestamp)
        << " and " << digits.size()
        << " entries (maximum sample: " << max_samples << ");"
        << " data at " << ((void*) digits.data());
      for (RawDigitInfo_t const& digitInfo: digits) {
        out << "\n  ";
        digitInfo.Dump(out);
      } // for
      out << "\n";
    } // RawDigitCacheClass::Dump()
    
    
    //--------------------------------------------------------------------------
    //--- GridAxisClass
    //---
    std::ptrdiff_t GridAxisClass::GetCell(float coord) const {
      return std::ptrdiff_t((coord - min) / cell_size); // truncate
    } // GridAxisClass::GetCell()
    
    
    //--------------------------------------------------------------------------
    bool GridAxisClass::Init(size_t nDiv, float new_min, float new_max) {
      
      min = new_min;
      max = new_max;
      n_cells = std::max(nDiv, size_t(1));
      
      cell_size = Length() / float(n_cells);
      
      return std::isnormal(cell_size);
    } // GridAxisClass::Init()
    
    
    //--------------------------------------------------------------------------
    bool GridAxisClass::SetMinCellSize(float min_size) {
      if (cell_size >= min_size) return false;
      
      // n_cells gets truncated
      n_cells = (size_t) std::max(std::floor(Length() / min_size), 1.0F);
      
      // reevaluate cell size, that might be different than min_size
      // because of n_cells truncation or minimum value
      cell_size = Length() / float(n_cells);
      return true;
    } // GridAxisClass::SetMinCellSize()
    
    
    //--------------------------------------------------------------------------
    bool GridAxisClass::SetMaxCellSize(float max_size) {
      if (cell_size <= max_size) return false;
      
      // n_cells gets rounded up
      n_cells = (size_t) std::max(std::ceil(Length() / max_size), 1.0F);
      
      // reevaluate cell size, that might be different than max_size
      // because of n_cells rounding or minimum value
      cell_size = Length() / float(n_cells);
      return true;
    } // GridAxisClass::SetMaxCellSize()
    
    
    //--------------------------------------------------------------------------
    //--- CellGridClass
    //---
    CellGridClass::CellGridClass(unsigned int nWires, unsigned int nTDC)
      : wire_axis((size_t) nWires, 0., float(nWires))
      , tdc_axis((size_t) nTDC, 0., float(nTDC))
    {
    } // CellGridClass::CellGridClass(int, int)
    
    
    //--------------------------------------------------------------------------
    CellGridClass::CellGridClass(
      float min_wire, float max_wire, unsigned int nWires,
      float min_tdc, float max_tdc, unsigned int nTDC
      )
      : wire_axis((size_t) nWires, min_wire, max_wire)
      , tdc_axis((size_t) nTDC, min_tdc, max_tdc)
    {
    } // CellGridClass::CellGridClass({ float, float, int } x 2)
    
    
    //--------------------------------------------------------------------------
    std::ptrdiff_t CellGridClass::GetCell(float wire, float tick) const {
      std::ptrdiff_t iWireCell = wire_axis.GetCell(wire);
      if (!wire_axis.hasCell(iWireCell)) return std::ptrdiff_t(-1);
      std::ptrdiff_t iTDCCell = tdc_axis.GetCell(tick);
      if (!tdc_axis.hasCell(iTDCCell)) return std::ptrdiff_t(-1);
      return iWireCell * TDCAxis().NCells() + iTDCCell;
    } // CellGridClass::GetCell()
    
    
    //--------------------------------------------------------------------------
    std::tuple<float, float, float, float> CellGridClass::GetCellBox
      (std::ptrdiff_t iCell) const
    {
      // { w1, t1, w2, t2 }
      register size_t const nWireCells = TDCAxis().NCells();
      std::ptrdiff_t iWireCell = (std::ptrdiff_t) (iCell / nWireCells),
        iTDCCell = (std::ptrdiff_t) (iCell % nWireCells);
      
      
      return std::tuple<float, float, float, float>(
        WireAxis().LowerEdge(iWireCell), TDCAxis().LowerEdge(iTDCCell),
        WireAxis().UpperEdge(iWireCell), TDCAxis().UpperEdge(iTDCCell)
        );
    } // CellGridClass::GetCellBox()
    
    
    //--------------------------------------------------------------------------
    
  } // details

} // namespace evd

////////////////////////////////////////////////////////////////////////
