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
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nuevdb/EventDisplayBase/Reconfigurable.h"

namespace evd {
    
class SimulationDrawingOptions : public evdb::Reconfigurable
{
public:
    explicit SimulationDrawingOptions(fhicl::ParameterSet const& pset);
    
    void reconfigure(fhicl::ParameterSet const& pset) ;
  
    bool                fShowMCTruthText;
    unsigned short      fShowMCTruthVectors;
    bool                fShowMCTruthTrajectories;
    bool                fShowSimChannelInfo;
    bool                fShowSimEnergyInfo;
    bool                fShowSimPhotonInfo;              ///< Display SimPhoton info in 3D display
    bool                fShowMCTruthColors;
    bool                fShowMCTruthFullSize;
    bool                fShowScintillationLight = false; ///< Whether to draw low energy light (default: no).
    double              fMinEnergyDeposition;
    art::InputTag       fG4ModuleLabel;                  ///< module label producing sim::SimChannel objects
    art::InputTag       fSimChannelLabel;                ///< SimChannels may be independent of MC stuff
    art::InputTag       fSimEnergyLabel;                 ///< Also for SimEnergyDeposits
    art::InputTag       fSimPhotonLabel;                 ///< and for SimPhotons
    
    fhicl::ParameterSet f3DDrawerParams;                  ///< FHICL paramegers for the 3D drawers
};
    
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::SimulationDrawingOptions, LEGACY)
#endif
