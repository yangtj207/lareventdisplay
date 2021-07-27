#ifndef ANALYSISDRAWINGOPTIONS_H
#define ANALYSISDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nuevdb/EventDisplayBase/Reconfigurable.h"

namespace evd {
  class AnalysisDrawingOptions : public evdb::Reconfigurable
  {
  public:
    explicit AnalysisDrawingOptions(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset) ;

    int fDrawCalorimetry;
    int fDrawParticleID;
    int fDrawShowerCalor;
    int fCaloPlane;
    int fTrackID;

    std::vector<std::string> fCalorimetryLabels;         ///< module labels that produced calorimetry
    std::vector<std::string> fParticleIDLabels;     	   ///< module labels that produced particleid

    std::string fCalorTemplateFileName;    ///< files that have calorimetry template curves

  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::AnalysisDrawingOptions, LEGACY)
#endif
