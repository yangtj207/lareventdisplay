////////////////////////////////////////////////////////////////////////
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef RAWDRAWINGOPTIONS_H
#define RAWDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::TPCID

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nuevdb/EventDisplayBase/Reconfigurable.h"

namespace evd {
  /**
   * @brief Display parameters for the raw data
   *
   * Configuration parameters
   * -------------------------
   *
   * This is an incomplete list of the supported parameters:
   * - *RoIthresholds* (list of real numbers, default: empty): threshold in ADC
   *   counts for a tick on a wire to be "interesting", thus extending the
   *   region of interest to include it. The thresholds are specified one per
   *   plane; if no threshold is specified for a plane, the threshold for the
   *   last plane is used instead (therefore specifying just one threshold will
   *   apply the same threshold to all planes). If no threshold is specified
   *   at all, the value of 'MinSignal' parameter is used as threshold for all
   *   planes
   *
   */
  class RawDrawingOptions : public evdb::Reconfigurable
  {
  public:
      explicit RawDrawingOptions(fhicl::ParameterSet const& pset);

      void reconfigure(fhicl::ParameterSet const& pset) override;

      int                        fDrawRawDataOrCalibWires;                 ///< 0 for raw
      int    	                 fTicksPerPoint;                           ///< number of ticks to include in one point
      int    	                 fScaleDigitsByCharge;                     ///< scale the size of the digit by the charge
      double 	                 fMinSignal;                               ///< minimum ADC count to display a time bin
      double                     fStartTick;                               ///< Starting tick for the display
      double 	                 fTicks;                                   ///< number of TDC ticks to display, ie # fTicks past fStartTick
      int    	                 fAxisOrientation;                         ///< 0 = TDC values on y-axis, wire number on x-axis, 1 = swapped
      unsigned int               fTPC;                                     ///< TPC number to draw, typically set by TWQProjectionView
      unsigned int               fCryostat;                                ///< Cryostat number to draw, typically set by TWQProjectionView
      unsigned int               fMinChannelStatus;                        ///< Display channels with this status and above
      unsigned int               fMaxChannelStatus;                        ///< Display channels with this status and below
      std::vector<art::InputTag> fRawDataLabels;                           ///< module label that made the raw digits, default is daq

      bool                       fUncompressWithPed;                       ///< Option to uncompress with pedestal. Turned off by default
      bool                       fSeeBadChannels;                          ///< Allow "bad" channels to be viewed
       
      std::vector<float>         fRoIthresholds;                           ///< region of interest thresholds, per plane
      
      int                        fPedestalOption;                          ///< 0: use DetPedestalService;   1:  Use pedestal in raw::RawDigt;   2:  no ped subtraction
      
      fhicl::ParameterSet        fRawDigitDrawerParams;                    ///< FHICL parameters for the RawDigit waveform display

      /// Returns the current TPC as a TPCID
      geo::TPCID   CurrentTPC() const { return geo::TPCID(fCryostat, fTPC); }

      /// Returns the region of interest threshold for the specified wire plane
      double RoIthreshold(geo::PlaneID const& planeID) const
      { return RoIthreshold(planeID.Plane); }

      /// Returns the region of interest threshold for the specified wire plane
      double RoIthreshold(geo::PlaneID::PlaneID_t plane) const
      {
        return (plane < fRoIthresholds.size())?
          fRoIthresholds[plane]: fRoIthresholds.back();
      } // RoIthreshold(plane number)

  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::RawDrawingOptions, LEGACY)
#endif
