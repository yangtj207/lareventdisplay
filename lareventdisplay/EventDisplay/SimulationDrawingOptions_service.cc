////////////////////////////////////////////////////////////////////////
/// \file SimulationDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

/// LArSoft includes
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

namespace evd {

  //......................................................................
  SimulationDrawingOptions::SimulationDrawingOptions(fhicl::ParameterSet const& pset)
  : evdb::Reconfigurable{pset}
  {
    this->reconfigure(pset);
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
    fShowSimChannelInfo      = pset.get< bool          >      ("ShowSimChannelInfo",       true);
    fShowSimEnergyInfo       = pset.get< bool          >      ("ShowSimEnergyInfo",        true);
    fShowSimPhotonInfo       = pset.get< bool          >      ("ShowSimPhotonInfo",        true);
    fShowMCTruthColors       = pset.get< bool          >      ("ShowMCTruthColors",        true);
    fShowMCTruthFullSize     = pset.get< bool          >      ("ShowMCTruthFullSize",      true);
    fShowScintillationLight  = pset.get< bool          >      ("ShowScintillationLight",  false);
    fMinEnergyDeposition     = pset.get< double        >      ("MinimumEnergyDeposition"       );
    fG4ModuleLabel           = pset.get< art::InputTag >      ("G4ModuleLabel"                 );
    fSimChannelLabel         = pset.get< art::InputTag >      ("SimChannelLabel"               );
    fSimEnergyLabel          = pset.get< art::InputTag >      ("SimEnergyLabel"                );
    fSimPhotonLabel          = pset.get< art::InputTag >      ("SimPhotonLabel"                );

    f3DDrawerParams          = pset.get< fhicl::ParameterSet >("Draw3DTools"                   );
  }

}

DEFINE_ART_SERVICE(evd::SimulationDrawingOptions)
