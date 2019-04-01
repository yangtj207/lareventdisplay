////////////////////////////////////////////////////////////////////////
/// \file   OpFlash3DDrawer_tool.cc
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

class OpFlash3DDrawer : public I3DDrawer
{
public:
    explicit OpFlash3DDrawer(const fhicl::ParameterSet&);
    
    ~OpFlash3DDrawer();
    
    void Draw(const art::Event&, evdb::View3D*) const override;
    
private:
    void DrawRectangularBox(evdb::View3D*, const Eigen::Vector3f&, const Eigen::Vector3f&, int, int, int) const;
};
    
//----------------------------------------------------------------------
// Constructor.
OpFlash3DDrawer::OpFlash3DDrawer(const fhicl::ParameterSet& pset)
{
    return;
}

OpFlash3DDrawer::~OpFlash3DDrawer()
{
}

void OpFlash3DDrawer::Draw(const art::Event& event, evdb::View3D* view) const
{
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    
    if (recoOpt->fDrawOpFlashes == 0) return;
    
    // Service recovery
    art::ServiceHandle<geo::Geometry>  geo;
    detinfo::DetectorProperties const* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<evd::ColorDrawingOptions> cst;
    
    std::vector<geo::PlaneID> planeIDVec;
    
    planeIDVec.push_back(geo::PlaneID(0,0,0));
    planeIDVec.push_back(geo::PlaneID(0,1,0));
    planeIDVec.push_back(geo::PlaneID(1,0,0));
    planeIDVec.push_back(geo::PlaneID(1,1,0));
    
    art::Handle<std::vector<recob::OpFlash>> opFlashHandle;
    
    // This seems like a time waster but we want to get the full color scale for all OpHits... so loops away...
    std::vector<float> opHitPEVec;
    
    // This is almost identically the same for loop we will re-excute below... sigh...
    for(size_t idx = 0; idx < recoOpt->fOpFlashLabels.size(); idx++)
    {
        art::InputTag opFlashProducer = recoOpt->fOpFlashLabels[idx];
        
        event.getByLabel(opFlashProducer, opFlashHandle);
        
        if (!opFlashHandle.isValid()   ) continue;
        if ( opFlashHandle->size() == 0) continue;
        
        // To get associations we'll need an art ptr vector...
        art::PtrVector<recob::OpFlash> opFlashVec;
        
        for(size_t idx = 0; idx < opFlashHandle->size(); idx++) opFlashVec.push_back(art::Ptr<recob::OpFlash>(opFlashHandle,idx));
        
        // Recover the associations to op hits
        art::FindManyP<recob::OpHit> opHitAssnVec(opFlashVec, event, opFlashProducer);
        
        if (opHitAssnVec.size() == 0) continue;
        
        // Start the loop over flashes
        for(const auto& opFlashPtr : opFlashVec)
        {
            std::cout << "--> opFlash PE: " << opFlashPtr->TotalPE() << ", Time: " << opFlashPtr->Time() << ", width: " << opFlashPtr->TimeWidth() << ", y/w: " << opFlashPtr->YCenter() << "/" << opFlashPtr->YWidth() << ", Z/w: " << opFlashPtr->ZCenter() << "/" << opFlashPtr->ZWidth() << std::endl;
            // Make some selections...
            if (opFlashPtr->TotalPE() < recoOpt->fFlashMinPE) continue;
            if (opFlashPtr->Time()    < recoOpt->fFlashTMin) continue;
            if (opFlashPtr->Time()    > recoOpt->fFlashTMax) continue;
            
            // Start by going through the associated OpHits
            const std::vector<art::Ptr<recob::OpHit>> opHitVec = opHitAssnVec.at(opFlashPtr.key());

            for(const auto& opHit : opHitVec) opHitPEVec.push_back(opHit->PE());
        }
    }
    
    // Sorting is good for mind and body...
    std::sort(opHitPEVec.begin(),opHitPEVec.end());
    
    float minTotalPE = opHitPEVec.front();
    float maxTotalPE = opHitPEVec[0.9 * opHitPEVec.size()];
    
    // Now we can set the scaling factor for PE
    float opHitPEScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) / (maxTotalPE - minTotalPE));

    // We are meant to draw the flashes/hits, so loop over the list of input flashes
    for(size_t idx = 0; idx < recoOpt->fOpFlashLabels.size(); idx++)
    {
        art::InputTag opFlashProducer = recoOpt->fOpFlashLabels[idx];
        
        event.getByLabel(opFlashProducer, opFlashHandle);
        
        if (!opFlashHandle.isValid()   ) continue;
        if ( opFlashHandle->size() == 0) continue;
        
        // To get associations we'll need an art ptr vector...
        art::PtrVector<recob::OpFlash> opFlashVec;
        
        for(size_t idx = 0; idx < opFlashHandle->size(); idx++) opFlashVec.push_back(art::Ptr<recob::OpFlash>(opFlashHandle,idx));
        
        // Recover the associations to op hits
        art::FindManyP<recob::OpHit> opHitAssnVec(opFlashVec, event, opFlashProducer);
        
        if (opHitAssnVec.size() == 0) continue;
        
        // Start the loop over flashes
        for(const auto& opFlashPtr : opFlashVec)
        {
            std::cout << "--> opFlash PE: " << opFlashPtr->TotalPE() << ", Time: " << opFlashPtr->Time() << ", width: " << opFlashPtr->TimeWidth() << ", y/w: " << opFlashPtr->YCenter() << "/" << opFlashPtr->YWidth() << ", Z/w: " << opFlashPtr->ZCenter() << "/" << opFlashPtr->ZWidth() << std::endl;
            // Make some selections...
            if (opFlashPtr->TotalPE() < recoOpt->fFlashMinPE) continue;
            if (opFlashPtr->Time()    < recoOpt->fFlashTMin) continue;
            if (opFlashPtr->Time()    > recoOpt->fFlashTMax) continue;
            
            // Start by going through the associated OpHits
            const std::vector<art::Ptr<recob::OpHit>> opHitVec = opHitAssnVec.at(opFlashPtr.key());
            
            // We use the flash time to give us an x position (for now... will need a better way eventually)
            float flashTick  = opFlashPtr->Time()/det->SamplingRate()*1e3 + det->GetXTicksOffset(planeIDVec[idx]);
            float flashWidth = opFlashPtr->TimeWidth()/det->SamplingRate()*1e3 + det->GetXTicksOffset(planeIDVec[idx]);
            
            // Now convert from time to distance...
            float flashXpos = det->ConvertTicksToX(flashTick,  planeIDVec[idx]);
            float flashXWid = det->ConvertTicksToX(flashWidth, planeIDVec[idx]);
            
            // Loop through the OpHits here
            for(const auto& opHit : opHitVec)
            {
                unsigned int         opChannel = opHit->OpChannel();
                const geo::OpDetGeo& opHitGeo  = geo->OpDetGeoFromOpChannel(opChannel);
                const geo::Point_t&  opHitPos  = opHitGeo.GetCenter();
                float                zWidth    = opHitGeo.HalfW();
                float                yWidth    = opHitGeo.HalfH();
                
                Eigen::Vector3f opHitLo(opHitPos.X() - flashXWid, opHitPos.Y() - yWidth, opHitPos.Z() - zWidth);
                Eigen::Vector3f opHitHi(opHitPos.X() + flashXWid, opHitPos.Y() + yWidth, opHitPos.Z() + zWidth);
                
                // Temporary kludge...
                flashXpos = opHitPos.X();
                
                float peFactor = cst->fRecoQLow[geo::kCollection] + opHitPEScale * std::min(maxTotalPE,float(opHit->PE()));
                
                int chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(peFactor);
                
                DrawRectangularBox(view, opHitLo, opHitHi, chargeColorIdx, 2, 1);
            }

            std::cout << "     == flashtick: " << flashTick << ", flashwidth: " << flashWidth << ", flashXpos: " << flashXpos << ", wid: " << flashXWid << std::endl;
            
            //            std::vector<Eigen::Vector3f>
            Eigen::Vector3f coordsLo(flashXpos - flashXWid,opFlashPtr->YCenter() - opFlashPtr->YWidth(),opFlashPtr->ZCenter() - opFlashPtr->ZWidth());
            Eigen::Vector3f coordsHi(flashXpos + flashXWid,opFlashPtr->YCenter() + opFlashPtr->YWidth(),opFlashPtr->ZCenter() + opFlashPtr->ZWidth());
            
            DrawRectangularBox(view, coordsLo, coordsHi, kRed, 2, 1);
            
        }
    }
    
    /*
     
     art::ServiceHandle<geo::Geometry>            geo;
     detinfo::DetectorProperties const* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
     geo::PlaneID pid(rawOpt->fCryostat, rawOpt->fTPC, plane);
     
     
     for(size_t imod = 0; imod < recoOpt->fOpFlashLabels.size(); ++imod) {
     const art::InputTag which = recoOpt->fOpFlashLabels[imod];
     
     art::PtrVector<recob::OpFlash> opflashes;
     this->GetOpFlashes(evt, which, opflashes);
     
     if(opflashes.size() < 1) continue;
     
     int NFlashes = opflashes.size();
     //double TopCoord = 1000;
     
     MF_LOG_VERBATIM("RecoBaseDrawer") <<"Total "<<NFlashes<<" flashes.";
     
     // project each seed into this view
     for (size_t iof = 0; iof < opflashes.size(); ++iof) {
     if (opflashes[iof]->TotalPE() < recoOpt->fFlashMinPE) continue;
     if (opflashes[iof]->Time() < recoOpt->fFlashTMin) continue;
     if (opflashes[iof]->Time() > recoOpt->fFlashTMax) continue;
     int Color = evd::kColor[(iof)%evd::kNCOLS];
     MF_LOG_VERBATIM("RecoBaseDrawer") << "Flash t: "
     << opflashes[iof]->Time()
     << "\t y,z : "
     << opflashes[iof]->YCenter()
     << ", "
     << opflashes[iof]->ZCenter()
     << " \t PE :"
     << opflashes[iof]->TotalPE();
     
     float flashtick = opflashes[iof]->Time()/det->SamplingRate()*1e3 + det->GetXTicksOffset(pid);
     float wire0 = FLT_MAX;
     float wire1 = FLT_MIN;
     
     //Find the 4 corners and convert them to wire numbers
     std::vector<TVector3> points;
     points.push_back(TVector3(0, opflashes[iof]->YCenter()-opflashes[iof]->YWidth(), opflashes[iof]->ZCenter()-opflashes[iof]->ZWidth()));
     points.push_back(TVector3(0, opflashes[iof]->YCenter()-opflashes[iof]->YWidth(), opflashes[iof]->ZCenter()+opflashes[iof]->ZWidth()));
     points.push_back(TVector3(0, opflashes[iof]->YCenter()+opflashes[iof]->YWidth(), opflashes[iof]->ZCenter()-opflashes[iof]->ZWidth()));
     points.push_back(TVector3(0, opflashes[iof]->YCenter()+opflashes[iof]->YWidth(), opflashes[iof]->ZCenter()+opflashes[iof]->ZWidth()));
     
     for (size_t i = 0; i<points.size(); ++i){
     geo::WireID wireID;
     try{
     wireID = geo->NearestWireID(points[i], pid);
     }
     catch(geo::InvalidWireError const& e) {
     wireID = e.suggestedWireID(); // pick the closest valid wire
     }
     if (wireID.Wire < wire0) wire0 = wireID.Wire;
     if (wireID.Wire > wire1) wire1 = wireID.Wire;
     }
     if(rawOpt->fAxisOrientation > 0){
     TLine& line = view->AddLine(flashtick, wire0, flashtick, wire1);
     line.SetLineWidth(2);
     line.SetLineStyle(2);
     line.SetLineColor(Color);
     }
     else{
     TLine& line = view->AddLine(wire0, flashtick, wire1, flashtick);
     line.SetLineWidth(2);
     line.SetLineStyle(2);
     line.SetLineColor(Color);
     }
     } // loop on opflashes
     } // loop on imod folders
     */

    return;
}

void OpFlash3DDrawer::DrawRectangularBox(evdb::View3D* view, const Eigen::Vector3f& coordsLo, const Eigen::Vector3f& coordsHi, int color, int width, int style) const
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


DEFINE_ART_CLASS_TOOL(OpFlash3DDrawer)
}
