////////////////////////////////////////////////////////////////////////
/// \file RawDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

/// LArSoft includes
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"

namespace evd {

  //......................................................................
  RawDrawingOptions::RawDrawingOptions(fhicl::ParameterSet const& pset)
  : evdb::Reconfigurable{pset}
  {
    this->reconfigure(pset);
  }

  //......................................................................
  void RawDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
      fDrawRawDataOrCalibWires    = pset.get< int                        >("DrawRawDataOrCalibWires"    );
      fScaleDigitsByCharge        = pset.get< int                        >("ScaleDigitsByCharge"        );
      fTicksPerPoint              = pset.get< int                        >("TicksPerPoint"              );
      fMinSignal                  = pset.get< double                     >("MinimumSignal"              );
      fStartTick                  = pset.get< double                     >("StartTick",            0    );
      fTicks                      = pset.get< double                     >("TotalTicks",           2048 );
      fAxisOrientation         	  = pset.get< int                        >("AxisOrientation",      0    );
      fRawDataLabels              = pset.get< std::vector<art::InputTag> >("RawDataLabels",        std::vector<art::InputTag>() = {"daq"} );
      fTPC                        = pset.get< unsigned int               >("TPC",                  0    );
      fCryostat                   = pset.get< unsigned int               >("Cryostat",             0    );
      fMinChannelStatus           = pset.get< unsigned int               >("MinChannelStatus",     0    );
      fMaxChannelStatus           = pset.get< unsigned int               >("MaxChannelStatus",     lariov::ChannelStatusProvider::InvalidStatus - 1);
      fUncompressWithPed          = pset.get< bool                       >("UncompressWithPed",    false);
      fSeeBadChannels             = pset.get< bool                       >("SeeBadChannels",       false);
      fRoIthresholds              = pset.get< std::vector<float>         >("RoIthresholds",        std::vector<float>());
      fPedestalOption             = pset.get< int                        >("PedestalOption",       0    );

      if (fRoIthresholds.empty()) fRoIthresholds.push_back((float) fMinSignal);

      fRawDigitDrawerParams       = pset.get< fhicl::ParameterSet >("RawDigitDrawer"             );
  }
}

DEFINE_ART_SERVICE(evd::RawDrawingOptions)
