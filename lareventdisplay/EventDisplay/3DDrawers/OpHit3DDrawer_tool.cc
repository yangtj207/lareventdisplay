////////////////////////////////////////////////////////////////////////
/// \file   OpHit3DDrawer_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"

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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

// Eigen
#include <Eigen/Dense>

#include <fstream>

namespace evdb_tool
{

class OpHit3DDrawer : public I3DDrawer
{
public:
    explicit OpHit3DDrawer(const fhicl::ParameterSet&);
    
    ~OpHit3DDrawer();
    
    void Draw(const art::Event&, evdb::View3D*) const override;
    
private:
    void DrawRectangularBox(evdb::View3D*, const Eigen::Vector3f&, const Eigen::Vector3f&, int, int, int) const;
};
    
//----------------------------------------------------------------------
// Constructor.
OpHit3DDrawer::OpHit3DDrawer(const fhicl::ParameterSet& pset)
{
    return;
}

OpHit3DDrawer::~OpHit3DDrawer()
{
}

void OpHit3DDrawer::Draw(const art::Event& event, evdb::View3D* view) const
{
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    
    if (recoOpt->fDrawOpHits == 0) return;
    
    // Service recovery
    art::ServiceHandle<geo::Geometry>  geo;
    art::ServiceHandle<evd::ColorDrawingOptions> cst;
    
    art::Handle<std::vector<recob::OpHit>> opHitHandle;
    
    // This seems like a time waster but we want to get the full color scale for all OpHits... so loops away...
    std::vector<float> opHitPEVec;
    
    // This is almost identically the same for loop we will re-excute below... sigh...
    // But the idea is to get a min/max range for the drawing colors...
    for(size_t idx = 0; idx < recoOpt->fOpHitLabels.size(); idx++)
    {
        art::InputTag opHitProducer = recoOpt->fOpHitLabels[idx];
        
        event.getByLabel(opHitProducer, opHitHandle);
        
        if (!opHitHandle.isValid()   ) continue;
        if ( opHitHandle->size() == 0) continue;
        
        // Start the loop over flashes
        for(const auto& opHit : *opHitHandle)
        {
            // Make some selections...
            if (opHit.PE()       < recoOpt->fFlashMinPE) continue;
            if (opHit.PeakTime() < recoOpt->fFlashTMin)  continue;
            if (opHit.PeakTime() > recoOpt->fFlashTMax)  continue;
            
            opHitPEVec.push_back(opHit.PE());
        }
    }
    
    // Do we have any flashes and hits?
    if (!opHitPEVec.empty())
    {
        // Sorting is good for mind and body...
        std::sort(opHitPEVec.begin(),opHitPEVec.end());
        
        float minTotalPE = opHitPEVec.front();
        float maxTotalPE = opHitPEVec[0.9 * opHitPEVec.size()];
        
        // Now we can set the scaling factor for PE
        float opHitPEScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) / (maxTotalPE - minTotalPE));
    
        // We are meant to draw the flashes/hits, so loop over the list of input flashes
        for(size_t idx = 0; idx < recoOpt->fOpHitLabels.size(); idx++)
        {
            art::InputTag opHitProducer = recoOpt->fOpHitLabels[idx];
            
            event.getByLabel(opHitProducer, opHitHandle);
            
            if (!opHitHandle.isValid()   ) continue;
            if ( opHitHandle->size() == 0) continue;
        
            // Start the loop over flashes
            for(const auto& opHit : *opHitHandle)
            {
                // Make some selections...
                if (opHit.PE()       < recoOpt->fFlashMinPE) continue;
                if (opHit.PeakTime() < recoOpt->fFlashTMin)  continue;
                if (opHit.PeakTime() > recoOpt->fFlashTMax)  continue;

                unsigned int         opChannel = opHit.OpChannel();
                const geo::OpDetGeo& opHitGeo  = geo->OpDetGeoFromOpChannel(opChannel);
                const geo::Point_t&  opHitPos  = opHitGeo.GetCenter();
                float                xWidth    = opHit.Width();
                float                zWidth    = opHitGeo.HalfW();
                float                yWidth    = opHitGeo.HalfH();
        
                Eigen::Vector3f opHitLo(opHitPos.X() - xWidth, opHitPos.Y() - yWidth, opHitPos.Z() - zWidth);
                Eigen::Vector3f opHitHi(opHitPos.X() + xWidth, opHitPos.Y() + yWidth, opHitPos.Z() + zWidth);
        
                float peFactor = cst->fRecoQLow[geo::kCollection] + opHitPEScale * std::min(maxTotalPE,float(opHit.PE()));
        
                int chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(peFactor);
        
                DrawRectangularBox(view, opHitLo, opHitHi, chargeColorIdx, 2, 1);
            }
        }
    }

    return;
}

void OpHit3DDrawer::DrawRectangularBox(evdb::View3D* view, const Eigen::Vector3f& coordsLo, const Eigen::Vector3f& coordsHi, int color, int width, int style) const
{
    TPolyLine3D& top = view->AddPolyLine3D(5, color, width, style);
    top.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(1, coordsHi[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(2, coordsHi[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);
    
    TPolyLine3D& side = view->AddPolyLine3D(5, color, width, style);
    side.SetPoint(0, coordsHi[0], coordsHi[1], coordsLo[2]);
    side.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    side.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    side.SetPoint(3, coordsHi[0], coordsHi[1], coordsHi[2]);
    side.SetPoint(4, coordsHi[0], coordsHi[1], coordsLo[2]);
    
    TPolyLine3D& side2 = view->AddPolyLine3D(5, color, width, style);
    side2.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    side2.SetPoint(1, coordsLo[0], coordsLo[1], coordsLo[2]);
    side2.SetPoint(2, coordsLo[0], coordsLo[1], coordsHi[2]);
    side2.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    side2.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);
    
    TPolyLine3D& bottom = view->AddPolyLine3D(5, color, width, style);
    bottom.SetPoint(0, coordsLo[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(3, coordsLo[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(4, coordsLo[0], coordsLo[1], coordsLo[2]);
    
    return;
}


DEFINE_ART_CLASS_TOOL(OpHit3DDrawer)
}
