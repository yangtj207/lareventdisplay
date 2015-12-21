////////////////////////////////////////////////////////////////////////
/// \file EvdLayoutOptions_service.cc
///
/// \version $Id: EvdLayoutOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
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

    fShowSideBar    	       = pset.get< int         >("ShowSideBar");
    fAutoZoomInterest          = pset.get< int         >("AutoZoomInterest");  
    fPrintTotalCharge  	       = pset.get< int         >("PrintTotalCharge"); 
    fShowEndPointSection       = pset.get< int         >("ShowEndPointSection"); 
    fShowEndPointMarkers       = pset.get< int         >("ShowEndPointMarkers");   
    fShowClusterSection        = pset.get< int	       >("ShowClusterSection");
    fMakeClusters 	       = pset.get< int	       >("MakeClusters");
    fMakeSeeds	 	       = pset.get< int	       >("MakeSeeds");
    fChangeWire		       = pset.get< int	       >("ChangeWire");
    fEnableMCTruthCheckBox     = pset.get< int         >("EnableMCTruthCheckBox");
  }   
}

namespace evd {

  DEFINE_ART_SERVICE(EvdLayoutOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
