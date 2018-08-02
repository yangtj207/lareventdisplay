////////////////////////////////////////////////////////////////////////
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SIMULATIONDRAWINGOPTIONS_H
#define SIMULATIONDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "nutools/EventDisplayBase/Reconfigurable.h"

namespace evd {
    
class SimulationDrawingOptions : public evdb::Reconfigurable
{
public:
    explicit SimulationDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~SimulationDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset) ;
  
    bool                fShowMCTruthText;
    unsigned short      fShowMCTruthVectors;
    bool                fShowMCTruthTrajectories;
    bool                fShowMCTruthColors;
    bool                fShowMCTruthFullSize;
    bool                fShowScintillationLight = false; ///< Whether to draw low energy light (default: no).
    double              fMinEnergyDeposition;
    art::InputTag       fG4ModuleLabel;                  ///< module label producing sim::SimChannel objects
    art::InputTag       fSimEnergyLabel;
    
    fhicl::ParameterSet f3DDrawerParams;                  ///< FHICL paramegers for the 3D drawers
};
    
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::SimulationDrawingOptions, LEGACY)
#endif

