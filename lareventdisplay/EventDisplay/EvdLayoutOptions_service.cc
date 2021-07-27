////////////////////////////////////////////////////////////////////////
/// \file EvdLayoutOptions_service.cc
///
/// \author  andrzejs@fnal.gov

/// LArSoft includes
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

namespace evd {

  //......................................................................
  EvdLayoutOptions::EvdLayoutOptions(fhicl::ParameterSet const& pset) :
                                     evdb::Reconfigurable{pset}, fParameterSet(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  void EvdLayoutOptions::reconfigure(fhicl::ParameterSet const& pset){

      fShowSideBar    	       = pset.get<     int     >("ShowSideBar");
      fAutoZoomInterest        = pset.get<     int     >("AutoZoomInterest");
      fPrintTotalCharge  	   = pset.get<     int     >("PrintTotalCharge");
      fShowEndPointSection     = pset.get<     int     >("ShowEndPointSection");
      fShowEndPointMarkers     = pset.get<     int     >("ShowEndPointMarkers");
      fShowClusterSection      = pset.get<     int     >("ShowClusterSection");
      fMakeClusters 	       = pset.get<     int     >("MakeClusters");
      fChangeWire		       = pset.get<     int     >("ChangeWire");
      fEnableMCTruthCheckBox   = pset.get<     int     >("EnableMCTruthCheckBox");

      fThreeWindow             = pset.get<     bool    >("DrawThreeWindow", true);
      fDrawGrid                = pset.get<     bool    >("DrawGrid",        true);
      fDrawAxes                = pset.get<     bool    >("DrawAxes",        true);
      fDrawBadChannels         = pset.get<     bool    >("DrawBadChannels", true);

      fDisplayName             = pset.get< std::string >("DisplayName",     "LArSoft");
  }
}

DEFINE_ART_SERVICE(evd::EvdLayoutOptions)
