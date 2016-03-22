////////////////////////////////////////////////////////////////////////
/// \file RawDrawingOption_plugin.cc
///
/// \version $Id: RecoDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// LArSoft includes
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "larevt/Filters/ChannelFilter.h"

#include <iostream>

namespace evd {

  //......................................................................
  RawDrawingOptions::RawDrawingOptions(fhicl::ParameterSet const& pset, 
                                       art::ActivityRegistry& /* reg */) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  RawDrawingOptions::~RawDrawingOptions() 
  {
  }

  //......................................................................
  void RawDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawRawDataOrCalibWires    = pset.get< int         >("DrawRawDataOrCalibWires"    );
    fScaleDigitsByCharge     	= pset.get< int         >("ScaleDigitsByCharge"        );
    fTicksPerPoint              = pset.get< int         >("TicksPerPoint"              );
    fMinSignal                  = pset.get< double      >("MinimumSignal"              );
    fStartTick                  = pset.get< double      >("StartTick",            0    );
    fTicks                      = pset.get< double      >("TotalTicks",           2048 );
    fAxisOrientation         	= pset.get< int         >("AxisOrientation",      0    );
    fRawDataLabel               = pset.get< art::InputTag >("RawDataLabel",        "daq");
    fTPC                        = pset.get< unsigned int>("TPC",                  0    );
    fCryostat                   = pset.get< unsigned int>("Cryostat",             0    );
    fMinChannelStatus           = pset.get< unsigned int>("MinChannelStatus",     0    );
    fMaxChannelStatus           = pset.get< unsigned int>("MaxChannelStatus",     lariov::ChannelStatusProvider::InvalidStatus - 1);
    fUncompressWithPed          = pset.get< bool >("UncompressWithPed",false);
    fSeeBadChannels             = pset.get< bool >("SeeBadChannels", false);
    fRoIthresholds              = pset.get< std::vector<float> >("RoIthresholds", std::vector<float>());
    
    if (fRoIthresholds.empty()) fRoIthresholds.push_back((float) fMinSignal);
  }
}

namespace evd {

  DEFINE_ART_SERVICE(RawDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
