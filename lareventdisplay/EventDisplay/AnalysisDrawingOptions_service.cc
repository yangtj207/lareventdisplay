////////////////////////////////////////////////////////////////////////
/// \file AnalysisDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

/// LArSoft includes
#include "lareventdisplay/EventDisplay/AnalysisDrawingOptions.h"


namespace evd {

  //......................................................................
  AnalysisDrawingOptions::AnalysisDrawingOptions(fhicl::ParameterSet const& pset)
  : evdb::Reconfigurable{pset}
  {
    this->reconfigure(pset);
  }

  //......................................................................
  void AnalysisDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawCalorimetry           = pset.get< int >("DrawCalorimetry"        );
    fDrawParticleID    	       = pset.get< int >("DrawParticleID"   	  );
    fDrawShowerCalor           = pset.get< int >("DrawShowerCalor"        );
    fCaloPlane                 = pset.get< int >("CaloPlane"              );
    fTrackID                   = pset.get< int >("TrackID"                );
    fCalorimetryLabels         = pset.get< std::vector<std::string> >("CalorimetryModuleLabels" );
    fParticleIDLabels          = pset.get< std::vector<std::string> >("ParticleIDModuleLabels"  );

    fCalorTemplateFileName     = pset.get< std::string >("CalorTemplateFileName"  );
  }

}

namespace evd {

  DEFINE_ART_SERVICE(AnalysisDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
