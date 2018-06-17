////////////////////////////////////////////////////////////////////////
/// \file   DrawGausHits_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWFHitDrawer.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TF1.h"
#include "TPolyLine.h"

#include <fstream>

namespace evdb_tool
{

class DrawGausHits : public IWFHitDrawer
{
public:
    explicit DrawGausHits(const fhicl::ParameterSet& pset);
    
    ~DrawGausHits();
    
    void configure(const fhicl::ParameterSet& pset)                                       override;
    void Draw(evdb::View2D&, std::vector<std::unique_ptr<TF1>>&, raw::ChannelID_t&) const override;
    
private:
    
    using HitParams_t = struct HitParams_t
    {
        float hitCenter;
        float hitSigma;
        float hitHeight;
        float hitStart;
        float hitEnd;
    };
    
    using ROIHitParamsVec = std::vector<HitParams_t>;
    using HitParamsVec    = std::vector<ROIHitParamsVec>;
    
    int     fNumPoints;
};
    
//----------------------------------------------------------------------
// Constructor.
DrawGausHits::DrawGausHits(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
DrawGausHits::~DrawGausHits()
{
}
    
void DrawGausHits::configure(const fhicl::ParameterSet& pset)
{
    fNumPoints = pset.get<int>("NumPoints", 1000);
    return;
}

    
void DrawGausHits::Draw(evdb::View2D&                      view2D,
                        std::vector<std::unique_ptr<TF1>>& hitFuncVec,
                        raw::ChannelID_t&                  channel) const
{
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    
    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if(!event) return;
    
    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod)
    {
        // Step one is to recover the hits for this label that match the input channel
        art::InputTag const which = recoOpt->fHitLabels[imod];
        
        std::vector<const recob::Hit*> hits;
        std::vector<const recob::Hit*> temp;
        
        event->getView(which, temp);
        
        for(const auto& hit : temp)
        {
            if (hit->Channel() == channel) hits.push_back(hit);
        }
        
        // Now go through and process the hits back into the hit parameters
        
        using HitParams_t = struct HitParams_t
        {
            float hitCenter;
            float hitSigma;
            float hitHeight;
            float hitStart;
            float hitEnd;
        };
        
        using ROIHitParamsVec = std::vector<HitParams_t>;
        using HitParamsVec    = std::vector<ROIHitParamsVec>;
        
        // Get an initial container for common hits on ROI
        HitParamsVec hitParamsVec;
        ROIHitParamsVec roiHitParamsVec;
        raw::TDCtick_t  lastEndTick(6400);
        
        for (size_t i = 0; i < hits.size(); ++i)
        {
            // check roi end condition
            if (hits[i]->EndTick() > lastEndTick)
            {
                if (!roiHitParamsVec.empty()) hitParamsVec.push_back(roiHitParamsVec);
                roiHitParamsVec.clear();
            }
            
            HitParams_t hitParams;
            
            hitParams.hitCenter = hits[i]->PeakTime();
            hitParams.hitSigma  = hits[i]->RMS();
            hitParams.hitHeight = hits[i]->PeakAmplitude();
            hitParams.hitStart  = hits[i]->StartTick();
            hitParams.hitEnd    = hits[i]->EndTick();
            
            roiHitParamsVec.emplace_back(hitParams);
            
            lastEndTick = hits[i]->EndTick();
        }//end loop over reco hits
        
        // Just in case (probably never called...)
        if (!roiHitParamsVec.empty()) hitParamsVec.push_back(roiHitParamsVec);
        
        size_t roiCount(0);
        
        for(const auto& roiHitParamsVec : hitParamsVec)
        {
            // Create a histogram here...
            double roiStart = roiHitParamsVec.front().hitStart; //roiHitParamsVec.front().hitStart - 3. * roiHitParamsVec.front().hitSigma;
            double roiStop  = roiHitParamsVec.back().hitEnd;    //roiHitParamsVec.back().hitEnd    + 3. * roiHitParamsVec.back().hitSigma;
            //                        double width    = roiStop - roiStart;
            
            std::string funcString = "gaus(0)";
            std::string funcName   = Form("hitshape_%05zu_c%02zu",size_t(channel),roiCount++);
            
            for(size_t idx = 1; idx < roiHitParamsVec.size(); idx++) funcString += "+gaus(" + std::to_string(3*idx) + ")";
            
            hitFuncVec.emplace_back(std::make_unique<TF1>(TF1(funcName.c_str(),funcString.c_str(),roiStart,roiStop)));
            
            TF1* f1 = hitFuncVec.back().get();
            
            size_t idx(0);
            for(const auto& hitParams : roiHitParamsVec)
            {
                f1->SetParameter(idx + 0, hitParams.hitHeight);
                f1->SetParameter(idx + 1, hitParams.hitCenter);
                f1->SetParameter(idx + 2, hitParams.hitSigma);
                
                TPolyLine& hitHeight = view2D.AddPolyLine(2, kBlack, 1, 1);
                
                hitHeight.SetPoint(0, hitParams.hitCenter, 0.);
                hitHeight.SetPoint(1, hitParams.hitCenter, hitParams.hitHeight);
                
                hitHeight.Draw("same");
                
                TPolyLine& hitSigma = view2D.AddPolyLine(2, kGray, 1, 1);
                
                hitSigma.SetPoint(0, hitParams.hitCenter - hitParams.hitSigma, 0.6 * hitParams.hitHeight);
                hitSigma.SetPoint(1, hitParams.hitCenter + hitParams.hitSigma, 0.6 * hitParams.hitHeight);
                
                hitSigma.Draw("same");
                
                idx += 3;
            }
        }
    }//end loop over HitFinding modules

    return;
}
    
DEFINE_ART_CLASS_TOOL(DrawGausHits)
}
