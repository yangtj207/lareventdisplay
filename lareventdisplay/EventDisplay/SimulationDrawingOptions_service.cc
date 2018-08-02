////////////////////////////////////////////////////////////////////////
/// \file SimulationDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes

/// LArSoft includes
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"

#include <iostream>

namespace evd {

  //......................................................................
  SimulationDrawingOptions::SimulationDrawingOptions(fhicl::ParameterSet const& pset, 
						     art::ActivityRegistry& /* reg */)
  : evdb::Reconfigurable{pset}
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  SimulationDrawingOptions::~SimulationDrawingOptions() 
  {
  }

  //......................................................................
  void SimulationDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fShowMCTruthText         = pset.get< bool          >      ("ShowMCTruthText",         true);
    try {
      fShowMCTruthVectors = pset.get< unsigned short >("ShowMCTruthVectors", 0);
    } // try
    catch (...) {
      std::cout<<"ShowMCTruthVectors changed to unsigned short. Please update your fcl configuration\n";
      fShowMCTruthVectors      = pset.get< bool >("ShowMCTruthVectors",0);
    } // catch
    fShowMCTruthTrajectories = pset.get< bool          >      ("ShowMCTruthTrajectories",  true);
    fShowMCTruthColors       = pset.get< bool          >      ("ShowMCTruthColors",        true);
    fShowMCTruthFullSize     = pset.get< bool          >      ("ShowMCTruthFullSize",      true);
    fShowScintillationLight  = pset.get< bool          >      ("ShowScintillationLight",  false);
    fMinEnergyDeposition     = pset.get< double        >      ("MinimumEnergyDeposition"       );
    fG4ModuleLabel           = pset.get< art::InputTag >      ("G4ModuleLabel"                 );
    fSimEnergyLabel          = pset.get< art::InputTag >      ("SimEnergyLabel"                );
      
    f3DDrawerParams          = pset.get< fhicl::ParameterSet >("Draw3DTools"                   );
  }
  
}

namespace evd {

  DEFINE_ART_SERVICE(SimulationDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
