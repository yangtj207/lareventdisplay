////////////////////////////////////////////////////////////////////////
/// \file EvdLayoutOptions_service.cc
///
/// \author  andrzejs@fnal.gov

// Framework includes

/// LArSoft includes
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"

#include <iostream>

namespace evd {

  //......................................................................
  EvdLayoutOptions::EvdLayoutOptions(fhicl::ParameterSet const& pset, 
				       art::ActivityRegistry& /* reg */) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  EvdLayoutOptions::~EvdLayoutOptions() 
  {
  }

  //......................................................................
  void EvdLayoutOptions::reconfigure(fhicl::ParameterSet const& pset){

    fShowSideBar    	       = pset.get< int  >("ShowSideBar");
    fAutoZoomInterest          = pset.get< int  >("AutoZoomInterest");
    fPrintTotalCharge  	       = pset.get< int  >("PrintTotalCharge");
    fShowEndPointSection       = pset.get< int  >("ShowEndPointSection");
    fShowEndPointMarkers       = pset.get< int  >("ShowEndPointMarkers");
    fShowClusterSection        = pset.get< int  >("ShowClusterSection");
    fMakeClusters 	           = pset.get< int  >("MakeClusters");
    fMakeSeeds	 	           = pset.get< int  >("MakeSeeds");
    fChangeWire		           = pset.get< int  >("ChangeWire");
    fEnableMCTruthCheckBox     = pset.get< int  >("EnableMCTruthCheckBox");
      
    fThreeWindow               = pset.get< bool >("DrawThreeWindow", true);
    fDrawGrid                  = pset.get< bool >("DrawGrid",        true);
    fDrawAxes                  = pset.get< bool >("DrawAxes",        true);
    fDrawBadChannels           = pset.get< bool >("DrawBadChannels", true);
  }
}

namespace evd {

  DEFINE_ART_SERVICE(EvdLayoutOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
