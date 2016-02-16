#ifndef ANALYSISDRAWINGOPTIONS_H
#define ANALYSISDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace evd {
  class AnalysisDrawingOptions 
  {
  public:
    AnalysisDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~AnalysisDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset);

    int fDrawCalorimetry;
    int fDrawParticleID;
    int fDrawShowerCalor;
    int fCaloPlane;

    std::vector<std::string> fCalorimetryLabels;         ///< module labels that produced calorimetry 
    std::vector<std::string> fParticleIDLabels;     	   ///< module labels that produced particleid

    std::string fCalorTemplateFileName;    ///< files that have calorimetry template curves

  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::AnalysisDrawingOptions, LEGACY)
#endif

