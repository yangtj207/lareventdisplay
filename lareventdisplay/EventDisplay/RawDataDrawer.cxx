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
 * - the fact that only the maximum charge is displayed on each box , and no
 *   dynamic charge range detection is in place
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
#include <typeinfo> // to use typeid()
#include <algorithm> // std::fill(), std::find_if(), ...
#include <cstddef> // std::ptrdiff_t

#include "TH1F.h"
// #include "TPolyLine3D.h"
#include "TBox.h"
#include "TFrame.h"
#include "TVirtualPad.h"

#include "lareventdisplay/EventDisplay/ChangeTrackers.h" // util::PlaneDataChangeTracker_t
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/Utilities/StatCollector.h" // lar::util::MinMaxCollector<>
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/demangle.h"


namespace {
  template <typename Stream, typename T>
  void PrintRange
    (Stream&& out, std::string header, lar::util::MinMaxCollector<T> const& range)
  {
    out << header << ": " << range.min() << " -- " << range.max()
      << " (" << (range.has_data()? "valid": "invalid") << ")";
  } // PrintRange()
} // local namespace


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
    class RawDigitCacheDataClass {
        public:
     
      /// Returns the list of digit info
      std::vector<RawDigitInfo_t> const& Digits() const { return digits; }
      
      /// Returns a pointer to the digit info of given channel, nullptr if none
      RawDigitInfo_t const* FindChannel(raw::ChannelID_t channel) const;
      
      /// Returns the largest number of samples in the unpacked raw digits
      size_t MaxSamples() const { return max_samples; }
      
      /// Returns whether the cache is empty() (STL-like interface)
      bool empty() const { return digits.empty(); }
      
      /// Empties the cache
      void Clear();
      
      /// Fills the cache from the specified raw digits product handle
      void Refill(art::Handle<std::vector<raw::RawDigit>>& rdcol);
      
      /// Clears the cache and marks it as invalid (use Update() to fill it)
      void Invalidate();
      
      /// Updates the cache for new_timestamp using the specified event
      /// @return true if it needed to update (that might have failed)
      bool Update(art::Event const& evt, CacheID_t const& new_timestamp);
      
      /// Dump the content of the cache
      template <typename Stream>
      void Dump(Stream&& out) const;
      
        private:
      
      struct BoolWithUpToDateMetadata {
        
        bool bUpToDate = false;
        std::vector<raw::RawDigit> const* digits = nullptr;
        
        BoolWithUpToDateMetadata() = default;
        BoolWithUpToDateMetadata
          (bool uptodate, std::vector<raw::RawDigit> const* newdigits):
          bUpToDate(uptodate), digits(newdigits)
          {}
        
        operator bool() const { return bUpToDate; }
      }; // struct BoolWithUpToDateMetadata
    
      std::vector<RawDigitInfo_t> digits; ///< vector of raw digit information
      
      CacheID_t timestamp; ///< object expressing validity range of cached data
      
      size_t max_samples = 0; ///< the largest number of ticks in any digit
      
      /// Checks whether an update is needed; can load digits in the process
      BoolWithUpToDateMetadata CheckUpToDate
        (CacheID_t const& ts, art::Event const* evt = nullptr) const;
      
      static std::vector<raw::RawDigit> const* ReadProduct
        (art::Event const& evt, art::InputTag label);
      
    }; // struct RawDigitCacheDataClass
    

    std::vector<evd::details::RawDigitInfo_t>::const_iterator begin
      (RawDigitCacheDataClass const& cache)
      { return cache.Digits().cbegin(); }
    std::vector<evd::details::RawDigitInfo_t>::const_iterator end
      (RawDigitCacheDataClass const& cache)
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
      
      /// Initialize the axis, returns whether cell size is finite
      bool Init(size_t nDiv, float new_min, float new_max);
      
      /// Initialize the axis limits, returns whether cell size is finite
      bool SetLimits(float new_min, float new_max);
      
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
      
      template <typename Stream>
      void Dump(Stream&& out) const;
      
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
      
      /// Sets the wire range, leaving the number of wire cells unchanged
      void SetWireRange(float min_wire, float max_wire)
        { wire_axis.SetLimits(min_wire, max_wire); }
      
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
      
      /// Sets the TDC range, leaving the number of ticks unchanged
      void SetTDCRange(float min_tdc, float max_tdc)
        { tdc_axis.SetLimits(min_tdc, max_tdc); }
      
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
      
      /// Prints the current axes on the specified stream
      template <typename Stream>
      void Dump(Stream&& out) const;
      
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
          detinfo::DetectorProperties const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
          return detp->BirksCorrection(dQdX);
        } // Correct()
      double operator() (float adc) const { return Correct(adc); }
      
      void update(geo::PlaneID const& pid)
        {
          art::ServiceHandle<geo::Geometry> geo;
          wirePitch = geo->WirePitch(pid);

          detinfo::DetectorProperties const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
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
    : digit_cache(new details::RawDigitCacheDataClass)
    , fStartTick(0),fTicks(2048)
    , fCacheID(new details::CacheID_t)
    , fDrawingRange(new details::CellGridClass)
  { 
    art::ServiceHandle<geo::Geometry> geo;

    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    geo::TPCID tpcid(rawopt->fCryostat, rawopt->fTPC);
    
    fStartTick = rawopt->fStartTick;
    fTicks = rawopt->fTicks;

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
    delete fCacheID;
  }
  
  
  //......................................................................
  void RawDataDrawer::SetDrawingLimits
    (float low_wire, float high_wire, float low_tdc, float high_tdc)
  {
    LOG_DEBUG("RawDataDrawer") << __func__
      << "() setting drawing range as wires ( "
      << low_wire << " - " << high_wire
      << " ), ticks ( " << low_tdc << " - " << high_tdc << " )";
    
    // we need to set the minimum cell size to 1, otherwise some cell will not
    // cover any wire/tick and they will be always empty
    if (PadResolution) {
      // TODO implement support for swapping axes here
      unsigned int wire_pixels = PadResolution.width;
      unsigned int tdc_pixels = PadResolution.height;
      fDrawingRange->SetWireRange(low_wire, high_wire, wire_pixels, 1.F);
      fDrawingRange->SetTDCRange(low_tdc, high_tdc, tdc_pixels, 1.F);
    }
    else {
      LOG_DEBUG("RawDataDrawer")
        << "Pad size not available -- using existing cell size";
      fDrawingRange->SetWireRange(low_wire, high_wire);
      fDrawingRange->SetMinWireCellSize(1.F);
      fDrawingRange->SetTDCRange(low_tdc, high_tdc);
      fDrawingRange->SetMinTDCCellSize(1.F);
    }
    
  } // RawDataDrawer::SetDrawingLimits()
  
  
  void RawDataDrawer::SetDrawingLimitsFromRoI
    (geo::PlaneID::PlaneID_t plane)
  {
    SetDrawingLimits
      (fWireMin[plane], fWireMax[plane], fTimeMin[plane], fTimeMax[plane]);
  } // RawDataDrawer::SetDrawingLimitsFromRoI()
  
  
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
      
      PadResolution.width = (unsigned int) wire_pixels;
      PadResolution.height = (unsigned int) tdc_pixels;
      
      log << "\n frame window is "
        << PadResolution.width << "x" << PadResolution.height
        << " pixel big and";
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
      
      // TODO support swapping axes here:
      // if (rawopt.fAxisOrientation < 1) { normal ; } else { swapped; }
      
      fDrawingRange->SetWireRange(low_wire, high_wire, (unsigned int) wire_pixels, 1.0);
      fDrawingRange->SetTDCRange(low_tdc, high_tdc, (unsigned int) tdc_pixels, 1.0);
    }
    else {
      // keep the old frame (if any)
      log << "\n  no frame!";
    }
    
  } // RawDataDrawer::ExtractRange()
  
  
  //......................................................................
  class RawDataDrawer::OperationBaseClass {
      public:
    OperationBaseClass
      (geo::PlaneID const& pid, RawDataDrawer* data_drawer = nullptr)
      : pRawDataDrawer(data_drawer), planeID(pid) {}
    
    virtual ~OperationBaseClass() = default;
    
    virtual bool Initialize() { return true; }
    
    virtual bool ProcessWire(geo::WireID const&) { return true; }
    virtual bool ProcessTick(size_t) { return true; }
    
    virtual bool Operate(geo::WireID const& wireID, size_t tick, float adc) = 0;
    
    virtual bool Finish() { return true; }
    
    virtual std::string Name() const
      { return cet::demangle_symbol(typeid(*this).name()); }
    
    bool operator() (geo::WireID const& wireID, size_t tick, float adc)
      { return Operate(wireID, tick, adc); }
    
    geo::PlaneID const& PlaneID() const { return planeID; }
    RawDataDrawer* RawDataDrawerPtr() const { return pRawDataDrawer; }
    
      protected:
    RawDataDrawer* pRawDataDrawer = nullptr;
      
      private:
    geo::PlaneID planeID;
  }; // class RawDataDrawer::OperationBaseClass
  
  //......................................................................
  class RawDataDrawer::ManyOperations: public RawDataDrawer::OperationBaseClass
  {
    std::vector<std::unique_ptr<RawDataDrawer::OperationBaseClass>> operations;
      public:
    
    ManyOperations
      (geo::PlaneID const& pid, RawDataDrawer* data_drawer = nullptr):
      OperationBaseClass(pid, data_drawer)
      {}
    
    virtual bool Initialize() override
      {
        bool bAllOk = true;
        for (std::unique_ptr<OperationBaseClass> const& op: operations)
          if (!op->Initialize()) bAllOk = false;
        return bAllOk;
      }
    
    virtual bool ProcessWire(geo::WireID const& wireID) override
      {
        for (std::unique_ptr<OperationBaseClass> const& op: operations)
          if (op->ProcessWire(wireID)) return true;
        return false;
      }
    virtual bool ProcessTick(size_t tick) override
      {
        for (std::unique_ptr<OperationBaseClass> const& op: operations)
          if (op->ProcessTick(tick)) return true;
        return false;
      }
    
    virtual bool Operate
      (geo::WireID const& wireID, size_t tick, float adc) override
      {
        for (std::unique_ptr<OperationBaseClass> const& op: operations)
          if (!op->Operate(wireID, tick, adc)) return false;
        return true;
      }
    
    virtual bool Finish() override
      {
        bool bAllOk = true;
        for (std::unique_ptr<OperationBaseClass> const& op: operations)
          if (!op->Finish()) bAllOk = false;
        return bAllOk;
      }
    
    virtual std::string Name() const
      {
        std::string msg = cet::demangle_symbol(typeid(*this).name());
        msg +=
          (" [running " + std::to_string(operations.size()) + " operations:");
        for (auto const& op: operations) {// it's unique_ptr<OperationBaseClass>
          if (op) msg += " " + op->Name();
          else msg += " <invalid>";
        }
        return msg + " ]";
      }
    
    OperationBaseClass* Operator(size_t iOp)
      { return operations.at(iOp).get(); }
    OperationBaseClass const* Operator(size_t iOp) const
      { return operations.at(iOp).get(); }
    
    void AddOperation(std::unique_ptr<OperationBaseClass> new_op)
      {
        if (!new_op) return;
        if (PlaneID() != new_op->PlaneID()) {
          throw art::Exception(art::errors::LogicError)
            << "RawDataDrawer::ManyOperations(): trying to run operations on "
            << std::string(PlaneID()) << " and "
            << std::string(new_op->PlaneID()) << " at the same time";
        }
        if (RawDataDrawerPtr()
          && (RawDataDrawerPtr() != new_op->RawDataDrawerPtr())
          )
        {
          throw art::Exception(art::errors::LogicError)
            << "RawDataDrawer::ManyOperations(): "
            "trying to run operations on different RawDataDrawer"
            ; // possible, but very unlikely
        }
        operations.emplace_back(std::move(new_op));
      }
    
  }; // class RawDataDrawer::ManyOperations
  
  //......................................................................
  bool RawDataDrawer::RunOperation
    (art::Event const& evt, OperationBaseClass* operation)
  {
    geo::PlaneID const& pid = operation->PlaneID();
    
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    details::CacheID_t NewCacheID(evt, rawopt->fRawDataLabel, pid);
    GetRawDigits(evt, NewCacheID);
    
    if(digit_cache->empty()) return true;
    
    LOG_DEBUG("RawDataDrawer") << "RawDataDrawer::RunOperation() running "
      << operation->Name();
    
    // if we have an initialization failure, return false immediately;
    // but it's way better if the failure throws an exception
    if (!operation->Initialize()) return false;
    
    lariov::ChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = *(lar::providerFrom<lariov::DetPedestalService>());
    
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
    
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
      
      std::vector<geo::WireID> WireIDs = geom.ChannelToWire(channel);
      
      bool bDrawChannel = false;
      for (geo::WireID const& wireID: WireIDs){
        if (wireID.planeID() != pid) continue; // not us!
        bDrawChannel = true;
        break;
      } // for wires
      if (!bDrawChannel) continue;
      
      // collect bad channels
      bool const bGood = rawopt->fSeeBadChannels || !channelStatus.IsBad(channel);
      
      // nothing else to be done if the channel is not good:
      // cells are marked bad by default and if any good channel falls in any of
      // them, they become good
      if (!bGood) continue;
      
      // at this point we know we have to process this channel
      raw::RawDigit::ADCvector_t const& uncompressed = digit_info.Data();
      
      // recover the pedestal
      float const pedestal = pedestalRetrievalAlg.PedMean(channel);
      
      // loop over all the wires that are covered by this channel;
      // without knowing better, we have to draw into all of them
      for (geo::WireID const& wireID: WireIDs){
        // check that the plane and tpc are the correct ones to draw
        if (wireID.planeID() != pid) continue; // not us!
        
        // do we have anything to do with this wire?
        if (!operation->ProcessWire(wireID)) continue;
        
        // get an iterator over the adc values
        // accumulate all the data of this wire in our "cells"
        size_t const max_tick = std::min({
          uncompressed.size(),
          size_t(fStartTick + fTicks)
          });
        
        for (size_t iTick = fStartTick; iTick < max_tick; ++iTick) {
          
          // do we have anything to do with this wire?
          if (!operation->ProcessTick(iTick)) continue;
          
          float const adc = uncompressed[iTick] - pedestal;
          
          if (!operation->Operate(wireID, iTick, adc)) return false;
          
        } // if good
      } // for wires
    } // for channels
    
    return operation->Finish();
  } // ChannelLooper()
  
  
  //......................................................................
  class RawDataDrawer::BoxDrawer: public RawDataDrawer::OperationBaseClass {
      public:
    
    BoxDrawer(
      geo::PlaneID const& pid,
      RawDataDrawer* dataDrawer,
      evdb::View2D* new_view
      )
      : OperationBaseClass(pid, dataDrawer)
      , view(new_view)
      , rawCharge(0.), convertedCharge(0.)
      , drawingRange(*(dataDrawer->fDrawingRange))
      , ADCCorrector(PlaneID())
      {}
    
    virtual bool Initialize() override
      {
        art::ServiceHandle<evd::RawDrawingOptions> rawopt;
        
        // set up the size of the grid to be visualized;
        // the information on the size has to be already there:
        // caller should have user ExtractRange(), or similar, first.
        // set the minimum cell in ticks to at least match fTicksPerPoint
        drawingRange.SetMinTDCCellSize((float) rawopt->fTicksPerPoint);
        // also set the minimum wire cell size to 1,
        // otherwise there will be cells represented by no wire.
        drawingRange.SetMinWireCellSize(1.F);
        boxInfo.clear();
        boxInfo.resize(drawingRange.NCells());
        return true;
      }
    
    virtual bool ProcessWire(geo::WireID const& wire) override
      { return drawingRange.hasWire((int) wire.Wire); }
    
    virtual bool ProcessTick(size_t tick) override
      { return drawingRange.hasTick((float) tick); }
    
    virtual bool Operate
      (geo::WireID const& wireID, size_t tick, float adc) override
      {
        geo::WireID::WireID_t const wire = wireID.Wire;
        std::ptrdiff_t cell = drawingRange.GetCell(wire, tick);
        if (cell < 0) return true;
        
        BoxInfo_t& info = boxInfo[cell];
        info.good = true; // if in range, we mark this cell as good
        
        rawCharge += adc;
        convertedCharge += ADCCorrector(adc);
        
        // draw maximum digit in the cell
        if (std::abs(info.adc) <= std::abs(adc)) info.adc = adc;
        
        return true;
      }
    
    virtual bool Finish() override
      {
        // write the information back
        geo::PlaneID::PlaneID_t const plane = PlaneID().Plane;
        RawDataDrawerPtr()->fRawCharge[plane] = rawCharge;
        RawDataDrawerPtr()->fConvertedCharge[plane] = convertedCharge;
        
        // the cell size might have changed because of minimum size settings
        // from configuration (see Initialize())
        *(RawDataDrawerPtr()->fDrawingRange) = drawingRange;
        
        // complete the drawing
        RawDataDrawerPtr()->QueueDrawingBoxes(view, PlaneID(), boxInfo);
        
        return true;
      }
    
      private:
    evdb::View2D* view;
    
    double rawCharge = 0., convertedCharge = 0.;
    details::CellGridClass drawingRange;
    std::vector<BoxInfo_t> boxInfo;
    details::ADCCorrectorClass ADCCorrector;
  }; // class RawDataDrawer::BoxDrawer
  
  
  void RawDataDrawer::QueueDrawingBoxes(
    evdb::View2D* view,
    geo::PlaneID const& pid,
    std::vector<BoxInfo_t> const& BoxInfo
    )
  {
    //
    // All the information is now collected in BoxInfo.
    // Make boxes out of it.
    //
    evd::RawDrawingOptions const& rawopt
      = *art::ServiceHandle<evd::RawDrawingOptions>();
    
    LOG_DEBUG("RawDataDrawer")
      << "Filling " << BoxInfo.size() << " boxes to be rendered";
    
    // drawing options:
    float const MinSignal = rawopt.fMinSignal;
    bool const bScaleDigitsByCharge = rawopt.fScaleDigitsByCharge;

    art::ServiceHandle<evd::ColorDrawingOptions> cst;
    
    geo::GeometryCore const& geom = *art::ServiceHandle<geo::Geometry>();
    geo::SigType_t const sigType = geom.SignalType(pid);
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
    /*
      LOG_TRACE("RawDataDrawer")
        << "Wires ( " << min_wire << " - " << max_wire << " ) ticks ("
        << min_tick << " - " << max_tick << " ) for cell " << iBox;
    */
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
      if (rawopt.fAxisOrientation < 1)
        pBox = &(view->AddBox(min_wire, min_tick, max_wire, max_tick));
      else
        pBox = &(view->AddBox(min_tick, min_wire, max_tick, max_wire));
      
      pBox->SetFillStyle(1001);
      pBox->SetFillColor(color);
      pBox->SetBit(kCannotPick);
      
      ++nDrawnBoxes;
    } // for (iBox)
    
    LOG_DEBUG("RawDataDrawer")
      << "Sent " << nDrawnBoxes << "/" << BoxInfo.size()
      << " boxes to be rendered";
  } // RawDataDrawer::QueueDrawingBoxes()
  
  
  void RawDataDrawer::RunDrawOperation
    (art::Event const& evt, evdb::View2D* view, unsigned int plane)
  {
    
    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if (rawopt->fDrawRawDataOrCalibWires == 1) return;
    
    geo::PlaneID const pid(rawopt->CurrentTPC(), plane);
    BoxDrawer drawer(pid, this, view);
    if (!RunOperation(evt, &drawer)) {
      throw art::Exception(art::errors::Unknown)
        << "RawDataDrawer::RunDrawOperation(): "
        "somewhere something went somehow wrong";
    }
    
  } // RawDataDrawer::RunDrawOperation()
  
  
  //......................................................................
  class RawDataDrawer::RoIextractorClass:
    public RawDataDrawer::OperationBaseClass
  {
      public:
    
    float const RoIthreshold;
    
    RoIextractorClass(geo::PlaneID const& pid, RawDataDrawer* data_drawer)
      : OperationBaseClass(pid, data_drawer)
      , RoIthreshold
        (art::ServiceHandle<evd::RawDrawingOptions>()->RoIthreshold(PlaneID()))
      {}
    
    virtual bool Operate
      (geo::WireID const& wireID, size_t tick, float adc) override
      {
        if (std::abs(adc) < RoIthreshold) return true;
        WireRange.add(wireID.Wire);
        TDCrange.add(tick);
        return true;
      } // Operate()
    
    virtual bool Finish() override
      {
        geo::PlaneID::PlaneID_t const plane = PlaneID().Plane;
        int& WireMin = pRawDataDrawer->fWireMin[plane];
        int& WireMax = pRawDataDrawer->fWireMax[plane];
        int& TimeMin = pRawDataDrawer->fTimeMin[plane];
        int& TimeMax = pRawDataDrawer->fTimeMax[plane];
        
        if ((WireMin == WireMax) && WireRange.has_data()) {
          geo::GeometryCore const& geom = *art::ServiceHandle<geo::Geometry>();
          mf::LogInfo("RawDataDrawer") << "Region of interest for "
            << std::string(PlaneID()) << " detected to be within wires "
            << WireRange.min() << " to " << WireRange.max()
            << " (plane has " << geom.Nwires(PlaneID()) << " wires)";
          WireMax = WireRange.max() + 1;
          WireMin = WireRange.min();
        }
        if ((TimeMin == TimeMax) && TDCrange.has_data()) {
          mf::LogInfo("RawDataDrawer") << "Region of interest for "
            << std::string(PlaneID()) << " detected to be within ticks "
            << TDCrange.min() << " to " << TDCrange.max();
          TimeMax = TDCrange.max() + 1;
          TimeMin = TDCrange.min();
        }
        return true;
      } // Finish()
    
      private:
    lar::util::MinMaxCollector<float> WireRange, TDCrange;
  }; // class RawDataDrawer::RoIextractorClass
  
  
  void RawDataDrawer::RunRoIextractor
    (art::Event const& evt, unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    geo::PlaneID const pid(rawopt->CurrentTPC(), plane);
    
    // if we have no region of interest, prepare to extract it
    bool const bExtractRoI = !hasRegionOfInterest(plane);
    LOG_TRACE("RawDataDrawer") << "Region of interest for " << pid
      << (bExtractRoI? " extracted": " not extracted") << " on this draw";
    
    if (!bExtractRoI) return;
    
    RoIextractorClass Extractor(pid, this);
    if (!RunOperation(evt, &Extractor)) {
      throw std::runtime_error
        ("RawDataDrawer::RunRoIextractor(): somewhere something went somehow wrong");
    }
    
  } // RawDataDrawer::RunRoIextractor()
  
  
  //......................................................................
  
  void RawDataDrawer::RawDigit2D(
    art::Event const& evt, evdb::View2D* view, unsigned int plane,
    bool bZoomToRoI /* = false */
    )
  {
    
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    geo::PlaneID const pid(rawopt->CurrentTPC(), plane);
    
    bool const bDraw = (rawopt->fDrawRawDataOrCalibWires != 1);
    // if we don't need to draw, don't bother doing anything;
    // if the region of interest is required, RunRoIextractor() should be called
    // (ok, now it's private, but it could be exposed)
    if (!bDraw) return;
    
    // make sure we reset what needs to be reset
    // before the operations are initialized;
    // we call for reading raw digits; they will be cached, so it's not a waste
    details::CacheID_t NewCacheID(evt, rawopt->fRawDataLabel, pid);
    GetRawDigits(evt, NewCacheID);
    
    bool const hasRoI = hasRegionOfInterest(plane);
    
    // - if we don't have a RoI yet, we want to get it while we draw
    //   * if we are zooming into it now, we have to extract it first, then draw
    //   * if we are not zooming, we can do both at the same time
    // - if we have a RoI, we don't want to extract it again
    if (!bZoomToRoI) { // we are not required to zoom to the RoI
      
      std::unique_ptr<OperationBaseClass> operation;
      
      // we will do the drawing in one pass
      LOG_DEBUG("RawDataDrawer")
        << __func__ << "() setting up one-pass drawing";
      operation.reset(new BoxDrawer(pid, this, view));
      
      if (!hasRoI) { // we don't have any RoI; since it's cheap, let's get it
        LOG_DEBUG("RawDataDrawer") << __func__ << "() adding RoI extraction";
        
        // swap cards: operation becomes a multiple operation:
        // - prepare the two operations (one is there already, somehow)
        std::unique_ptr<OperationBaseClass> drawer(std::move(operation));
        std::unique_ptr<OperationBaseClass> extractor
          (new RoIextractorClass(pid, this));
        // - create a new composite operation and give it the sub-ops
        operation.reset(new ManyOperations(pid, this));
        ManyOperations* pManyOps
          = static_cast<ManyOperations*>(operation.get());
        pManyOps->AddOperation(std::move(drawer));
        pManyOps->AddOperation(std::move(extractor));
      }
      
      if (!RunOperation(evt, operation.get())) {
        throw art::Exception(art::errors::Unknown)
          << "RawDataDrawer::RunDrawOperation(): "
          "somewhere something went somehow wrong";
      }
    }
    else { // we are zooming to RoI
      // first, we want the RoI extracted; the extractor will update this object
      if (!hasRoI) {
        LOG_DEBUG("RawDataDrawer") << __func__
          << "() setting up RoI extraction for " << pid;
        RoIextractorClass extractor(pid, this);
        if (!RunOperation(evt, &extractor)) {
          throw art::Exception(art::errors::Unknown)
            << "RawDataDrawer::RunDrawOperation():"
            " something went somehow wrong while extracting RoI";
        }
      }
      else {
        LOG_DEBUG("RawDataDrawer") << __func__
          << "() using existing RoI for " << pid
          << ": wires ( " << fWireMin[plane] << " - " << fWireMax[plane]
          << " ), ticks ( " << fTimeMin[plane] << " - " << fTimeMax[plane]
          << " )";
      }
      
      // adopt the drawing limits information from the wire/time limits
      SetDrawingLimitsFromRoI(pid);
      
      // then we draw
      LOG_DEBUG("RawDataDrawer") << __func__ << "() setting up drawing";
      BoxDrawer drawer(pid, this, view);
      if (!RunOperation(evt, &drawer)) {
        throw art::Exception(art::errors::Unknown)
          << "RawDataDrawer::RunDrawOperation():"
          " something went somehow wrong while drawing";
      }
    }
  } // RawDataDrawer::RawDigit2D()
  
  
  //........................................................................
  int RawDataDrawer::GetRegionOfInterest(int plane,int& minw,int& maxw,int& mint,int& maxt)
  {
    art::ServiceHandle<geo::Geometry> geo;
 
    if((unsigned int)plane>=fWireMin.size())
      {mf::LogWarning  ("RawDataDrawer") << " Requested plane " << plane <<" is larger than those available " << std::endl;
	return -1;
      }
    
    minw=fWireMin[plane];
    maxw=fWireMax[plane];
    mint=fTimeMin[plane];
    maxt=fTimeMax[plane];
    
    if ((minw == maxw) || (mint == maxt)) return 1;
    
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
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if (rawopt->fDrawRawDataOrCalibWires==1) return;

    art::ServiceHandle<geo::Geometry> geo;
    
    geo::PlaneID const pid(rawopt->CurrentTPC(), plane);
    details::CacheID_t NewCacheID(evt, rawopt->fRawDataLabel, pid);
    GetRawDigits(evt, NewCacheID);
      
    lariov::ChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();

    for (evd::details::RawDigitInfo_t const& digit_info: *digit_cache) {
      raw::RawDigit const& hit = digit_info.Digit();
      raw::ChannelID_t const channel = hit.Channel();
        
      if (!channelStatus.IsPresent(channel)) continue;
      
      // The following test is meant to be temporary until the "correct" solution is implemented
      if (!ProcessChannelWithStatus(channelStatus.Status(channel))) continue;
      
      // to be explicit: we don't cound bad channels in
      if (!rawopt->fSeeBadChannels && channelStatus.IsBad(channel)) continue;
      
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
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if (rawopt->fDrawRawDataOrCalibWires==1) return;
    
    // make sure we have the raw digits cached
    geo::PlaneID const pid(rawopt->CurrentTPC(), plane);
    details::CacheID_t NewCacheID(evt, rawopt->fRawDataLabel, pid);
    GetRawDigits(evt, NewCacheID);

    if (digit_cache->empty()) return;

    geo::WireID const wireid(pid, wire);
    
    // find the channel
    art::ServiceHandle<geo::Geometry> geom;
    raw::ChannelID_t const channel = geom->PlaneWireToChannel(wireid);
    if (!raw::isValidChannelID(channel)) { // no channel, empty histogram
      mf::LogError("RawDataDrawer") << __func__ << ": no channel associated to "
        << std::string(wireid);
      return;
    } // if no channel
    
    // check the channel status; bad channels are still ok.
    lariov::ChannelStatusProvider const& channelStatus
      = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    if (!channelStatus.IsPresent(channel)) return;
    
    // The following test is meant to be temporary until the "correct" solution is implemented
    if (!ProcessChannelWithStatus(channelStatus.Status(channel))) return;
    
    
    // we accept to see the content of a bad channel, so this is commented out:
    if (!rawopt->fSeeBadChannels && channelStatus.IsBad(channel)) return;
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
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
  //     art::ServiceHandle<evd::RawDrawingOptions> rawopt;
  //     if (rawopt->fDrawRawOrCalibHits!=0) return;


  //     art::ServiceHandle<geo::Geometry> geom;

  //     HitTower tower;
  //     tower.fQscale = 0.01;

  //     for (unsigned int imod=0; imod<rawopt->fRawDigitModules.size(); ++imod) {
  //       const char* which = rawopt->fRawDigitModules[imod].c_str();

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
      
  // 	switch (rawopt->fRawDigit3DStyle) {
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
  //     if (rawopt->fRawDigit3DStyle==1) tower.Draw(view);
  //   }

  //......................................................................    
  bool RawDataDrawer::hasRegionOfInterest(geo::PlaneID::PlaneID_t plane) const {
    
    return (fWireMax[plane] != -1) && (fTimeMax[plane] != -1);
    
  } // RawDataDrawer::hasRegionOfInterest()
  
  
  //......................................................................    
  void RawDataDrawer::ResetRegionOfInterest() {
    
    LOG_DEBUG("RawDataDrawer") << "RawDataDrawer[" << ((void*) this)
      << "]: resetting the region of interest";
    
    std::fill(fWireMin.begin(), fWireMin.end(), -1);
    std::fill(fWireMax.begin(), fWireMax.end(), -1);
    std::fill(fTimeMin.begin(), fTimeMin.end(), -1);
    std::fill(fTimeMax.begin(), fTimeMax.end(), -1);
    
  } // RawDataDrawer::ResetRegionOfInterest()
  
  
  //......................................................................    

  void RawDataDrawer::GetRawDigits
    (art::Event const& evt, details::CacheID_t const& new_timestamp)
  {
    LOG_DEBUG("RawDataDrawer") << "GetRawDigits() for " << new_timestamp
      << " (last for: " << *fCacheID << ")";
    
    // update cache
    digit_cache->Update(evt, new_timestamp);
    
    // if time stamp is changing, we want to reconsider which region is
    // interesting
    if (!fCacheID->sameTPC(new_timestamp)) ResetRegionOfInterest();
    
    // all the caches have been properly updated or invalidated;
    // we are now on a new cache state
    *fCacheID = new_timestamp;
    
  } // RawDataDrawer::GetRawDigits()
  
  
  //......................................................................    
  bool RawDataDrawer::ProcessChannelWithStatus
    (lariov::ChannelStatusProvider::Status_t channel_status) const
  {
    // if we don't have a valid status, we can't reject the channel
    if (!lariov::ChannelStatusProvider::IsValidStatus(channel_status))
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

      art::ServiceHandle<evd::RawDrawingOptions> drawopt;      

      if (digit->Compression() == kNone) {
        // no compression, we can refer to the original data directly
        data.PointToData(digit->ADCs());
      }
      else {
        // data is compressed, need to do the real work
	if (drawopt->fUncompressWithPed){//Use pedestal in uncompression
	  int pedestal = (int)digit->GetPedestal();
	  raw::RawDigit::ADCvector_t samples;
	  Uncompress(digit->ADCs(), samples, pedestal, digit->Compression());
	  data.StealData(std::move(samples));
	}
	else{
	  raw::RawDigit::ADCvector_t samples;
	  Uncompress(digit->ADCs(), samples, digit->Compression());
	  data.StealData(std::move(samples));
	}
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
    //--- RawDigitCacheDataClass
    //---
   
    RawDigitInfo_t const* RawDigitCacheDataClass::FindChannel
      (raw::ChannelID_t channel) const
    {
      auto iDigit = std::find_if(
        digits.cbegin(), digits.cend(), 
        [channel](evd::details::RawDigitInfo_t const& digit)
          { return digit.Channel() == channel; }
        );
      return (iDigit == digits.cend())? nullptr: &*iDigit;
    } // RawDigitCacheDataClass::FindChannel()
    
    std::vector<raw::RawDigit> const* RawDigitCacheDataClass::ReadProduct
      (art::Event const& evt, art::InputTag label)
    {
      art::Handle< std::vector<raw::RawDigit>> rdcol;
      if (!evt.getByLabel(label, rdcol)) return nullptr;
      return &*rdcol;
    } // RawDigitCacheDataClass::ReadProduct()
    
    
    void RawDigitCacheDataClass::Refill
      (art::Handle<std::vector<raw::RawDigit>>& rdcol)
    {
      digits.resize(rdcol->size());
      for(size_t iDigit = 0; iDigit < rdcol->size(); ++iDigit) {
        art::Ptr<raw::RawDigit> pDigit(rdcol, iDigit);
        digits[iDigit].Fill(pDigit);
        size_t samples = pDigit->Samples();
        if (samples > max_samples) max_samples = samples;
      } // for
    } // RawDigitCacheDataClass::Refill()
    
    
    void RawDigitCacheDataClass::Invalidate() {
      timestamp.clear();
    } // RawDigitCacheDataClass::Invalidate()
    
    
    void RawDigitCacheDataClass::Clear() {
      Invalidate();
      digits.clear();
      max_samples = 0;
    } // RawDigitCacheDataClass::Clear()
    
    
    RawDigitCacheDataClass::BoolWithUpToDateMetadata
    RawDigitCacheDataClass::CheckUpToDate
      (CacheID_t const& ts, art::Event const* evt /* = nullptr */) const
    {
      BoolWithUpToDateMetadata res{ false, nullptr };
      
      // normally only if either the event or the product label have changed,
      // cache becomes invalid:
      if (!ts.sameProduct(timestamp)) return res; // outdated cache
      
      // But: our cache stores pointers to the original data, and on a new TPC
      // the event display may reload the event anew, removing the "old" data
      // from memory.
      // Since TPC can change with or without the data being invalidated,
      // a more accurate verification is needed.
      
      // if the cache is empty, well, it does not make much difference;
      // we invalidate it to be sure
      if (empty())
        return res; // outdated cache
      
      if (!evt)
        return res; // outdated, since we can't know better without the event
      
      // here we force reading of the product
      res.digits = ReadProduct(*evt, ts.inputLabel());
      if (!res.digits)
        return res; // outdated cache; this is actually an error
      
      if (res.digits->empty())
        return res; // outdated; no digits (strange!), invalidate just in case
      
      // use the first digit as test
      raw::ChannelID_t channel = res.digits->front().Channel();
      RawDigitInfo_t const* pInfo = FindChannel(channel);
      if (!pInfo)
        return res; // outdated: we don't even have this channel in cache!
      
      if (&(pInfo->Digit()) != &(res.digits->front()))
        return res; // outdated: different memory address for data
      
      res.bUpToDate = true;
      return res; // cache still valid
    } // RawDigitCacheDataClass::CheckUpToDate()
    
    
    bool RawDigitCacheDataClass::Update
      (art::Event const& evt, CacheID_t const& new_timestamp)
    {
      BoolWithUpToDateMetadata update_info = CheckUpToDate(new_timestamp, &evt);
      
      if (update_info) return false; // already up to date: move on!
      
      LOG_DEBUG("RawDataDrawer")
        << "Refilling raw digit cache RawDigitCacheDataClass["
        << ((void*) this ) << "] for " << new_timestamp;
      
      Clear();
      
      art::Handle< std::vector<raw::RawDigit>> rdcol;
      if (!evt.getByLabel(new_timestamp.inputLabel(), rdcol)) {
        mf::LogWarning("RawDataDrawer") << "no RawDigit collection '"
          << new_timestamp.inputLabel() << "' found";
        return true;
      }
      
      Refill(rdcol);
      
      timestamp = new_timestamp;
      return true;
    } // RawDigitCacheDataClass::Update()
    
    
    template <typename Stream>
    void RawDigitCacheDataClass::Dump(Stream&& out) const {
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
    } // RawDigitCacheDataClass::Dump()
    
    
    //--------------------------------------------------------------------------
    //--- GridAxisClass
    //---
    std::ptrdiff_t GridAxisClass::GetCell(float coord) const {
      return std::ptrdiff_t((coord - min) / cell_size); // truncate
    } // GridAxisClass::GetCell()
    
    
    //--------------------------------------------------------------------------
    bool GridAxisClass::Init(size_t nDiv, float new_min, float new_max) {
      
      n_cells = std::max(nDiv, size_t(1));
      return SetLimits(new_min, new_max);
      
    } // GridAxisClass::Init()
    
    
    //--------------------------------------------------------------------------
    bool GridAxisClass::SetLimits(float new_min, float new_max) {
      min = new_min;
      max = new_max;
      cell_size = Length() / float(n_cells);
      
      return std::isnormal(cell_size);
    } // GridAxisClass::SetLimits()
    
    
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
    template <typename Stream>
    void GridAxisClass::Dump(Stream&& out) const {
      out << NCells() << " cells from " << Min() << " to " << Max()
        << " (length: " << Length() << ")";
    } // GridAxisClass::Dump()
    
    
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
      register size_t const nTDCCells = TDCAxis().NCells();
      std::ptrdiff_t iWireCell = (std::ptrdiff_t) (iCell / nTDCCells),
        iTDCCell = (std::ptrdiff_t) (iCell % nTDCCells);
      
      
      return std::tuple<float, float, float, float>(
        WireAxis().LowerEdge(iWireCell), TDCAxis().LowerEdge(iTDCCell),
        WireAxis().UpperEdge(iWireCell), TDCAxis().UpperEdge(iTDCCell)
        );
    } // CellGridClass::GetCellBox()
    
    
    //--------------------------------------------------------------------------
    template <typename Stream>
    void CellGridClass::Dump(Stream&& out) const {
      out << "Wire axis: ";
      WireAxis().Dump(out);
      out << "; time axis: ";
      TDCAxis().Dump(out);
    } // CellGridClass::Dump()
    
    
    //--------------------------------------------------------------------------
    
  } // details

} // namespace evd

////////////////////////////////////////////////////////////////////////
