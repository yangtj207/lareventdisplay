////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawer_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lareventdisplay/EventDisplay/Style.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

#include <fstream>

namespace evdb_tool
{

class SpacePoint3DDrawer : public I3DDrawer
{
public:
    explicit SpacePoint3DDrawer(const fhicl::ParameterSet&);
    
    ~SpacePoint3DDrawer();
    
    void Draw(const art::Event&, evdb::View3D*) const override;
    
private:
};
    
//----------------------------------------------------------------------
// Constructor.
SpacePoint3DDrawer::SpacePoint3DDrawer(const fhicl::ParameterSet& pset)
{
//    fNumPoints     = pset.get< int>("NumPoints",     1000);
//    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.

    return;
}

SpacePoint3DDrawer::~SpacePoint3DDrawer()
{
    return;
}

void SpacePoint3DDrawer::Draw(const art::Event& evt, evdb::View3D* view) const
{
/*
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    
    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;
    
    std::vector<art::InputTag> labels;
    if(recoOpt->fDrawSpacePoints != 0)
    {
        for(size_t imod = 0; imod < recoOpt->fSpacePointLabels.size(); ++imod)
            labels.push_back(recoOpt->fSpacePointLabels[imod]);
    }
    
    for(size_t imod = 0; imod < labels.size(); ++imod)
    {
        art::InputTag const which = labels[imod];
        
        std::vector<art::Ptr<recob::SpacePoint>> spts;
        this->GetSpacePoints(evt, which, spts);
        int color = 10*imod + 11;
        
        color = 0;
        
        //        std::vector<const recob::SpacePoint*> sptsVec;
        //
        //        sptsVec.resize(spts.size());
        //        for(const auto& spt : spts){
        //          std::cout<<spt<<" "<<*spt<<" "<<&*spt<<std::endl;
        //          sptsVec.push_back(&*spt);
        //          std::cout<<sptsVec.back()<<std::endl;
        //        }
        this->DrawSpacePoint3D(spts, view, color, kFullDotMedium, 1);
    }
*/
   return;
}

DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawer)
}
