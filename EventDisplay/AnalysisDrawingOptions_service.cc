////////////////////////////////////////////////////////////////////////
/// \file AnalysisDrawingOption_plugin.cc
///
/// \version $Id: AnalysisDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// LArSoft includes
#include "EventDisplay/AnalysisDrawingOptions.h"

#include <iostream>

namespace evd {

  //......................................................................
  AnalysisDrawingOptions::AnalysisDrawingOptions(fhicl::ParameterSet const& pset, 
					 art::ActivityRegistry& reg) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  AnalysisDrawingOptions::~AnalysisDrawingOptions() 
  {
  }

  //......................................................................
  void AnalysisDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawCalorimetry           = pset.get< int >("DrawCalorimetry"        );
    fDrawParticleID    	       = pset.get< int >("DrawParticleID"   	    );

    fCalorimetryLabels         = pset.get< std::vector<std::string> >("CalorimetryModuleLabels" );
    fParticleIDLabels          = pset.get< std::vector<std::string> >("ParticleIDModuleLabels"  );

    fCalorTemplateFileNames    = pset.get< std::vector<std::string> >("CalorTemplateFileNames"  );
  }
  
}

namespace evd {

  DEFINE_ART_SERVICE(AnalysisDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
