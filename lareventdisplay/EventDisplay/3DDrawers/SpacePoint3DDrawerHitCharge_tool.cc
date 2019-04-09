////////////////////////////////////////////////////////////////////////
/// \file   SpacePoint3DDrawerHitCharge_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/3DDrawers/ISpacePoints3D.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lareventdisplay/EventDisplay/Style.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "TMath.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TDatabasePDG.h"

#include <fstream>

namespace evdb_tool
{

class SpacePoint3DDrawerHitCharge : public ISpacePoints3D
{
public:
    explicit SpacePoint3DDrawerHitCharge(const fhicl::ParameterSet&);
    
    ~SpacePoint3DDrawerHitCharge();
    
    void Draw(const std::vector<art::Ptr<recob::SpacePoint>>&,  // Space points
              evdb::View3D*,                                    // 3D display
              int,                                              // Color
              int,                                              // Marker
              float,                                            // Size) const override;
              const art::FindManyP<recob::Hit>*                 // pointer to associated hits
             ) const;
    
private:
    double chargeIntegral(double,double,double,double,int,int) const;
    
    float          fMinHitCharge;
    float          fMaxHitCharge;
};
    
//----------------------------------------------------------------------
// Constructor.
SpacePoint3DDrawerHitCharge::SpacePoint3DDrawerHitCharge(const fhicl::ParameterSet& pset)
{
//    fNumPoints     = pset.get< int>("NumPoints",     1000);
//    fFloatBaseline = pset.get<bool>("FloatBaseline", false);
    // For now only draw cryostat=0.
    fMinHitCharge = pset.get<float>("MinHitCharge",    0.);
    fMaxHitCharge = pset.get<float>("MaxHitCharge", 2500.);

    return;
}

SpacePoint3DDrawerHitCharge::~SpacePoint3DDrawerHitCharge()
{
    return;
}

void SpacePoint3DDrawerHitCharge::Draw(const std::vector<art::Ptr<recob::SpacePoint>>& hitsVec,
                                       evdb::View3D*                                   view,
                                       int                                             color,
                                       int                                             marker,
                                       float                                           size,
                                       const art::FindManyP<recob::Hit>*               hitAssnVec) const
{
    // Let's not crash
    if (hitsVec.empty() || !hitAssnVec) return;
    
    // Get services.
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;
    
    using HitPosition = std::array<double,6>;
    std::map<int,std::vector<HitPosition>> colorToHitMap;
    
    float hitChiSqScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) / (fMaxHitCharge - fMinHitCharge));
    
    for(const auto& spacePoint : hitsVec)
    {
        const double* pos = spacePoint->XYZ();
        const double* err = spacePoint->ErrXYZ();
        
        // Need to recover the integrated charge from the collection plane, so need to loop through associated hits
        const std::vector<art::Ptr<recob::Hit>>& hit2DVec(hitAssnVec->at(spacePoint.key()));
        
        float hitCharge(0.);
        int   lowIndex(std::numeric_limits<int>::min());
        int   hiIndex(std::numeric_limits<int>::max());
        
        for(const auto& hit2D : hit2DVec)
        {
            int hitStart = hit2D->PeakTime() - 2. * hit2D->RMS() - 0.5;
            int hitStop  = hit2D->PeakTime() + 2. * hit2D->RMS() + 0.5;
            
            lowIndex = std::max(hitStart,    lowIndex);
            hiIndex  = std::min(hitStop + 1, hiIndex);
            
            hitCharge += hit2D->Integral();
        }
        
        if (!hit2DVec.empty()) hitCharge /= float(hit2DVec.size());
        
        if (hitCharge > 0.)
        {
            int   chargeColorIdx(0);
            float integral(0.);
            
            if (hiIndex > lowIndex)
            {
                for(const auto& hit2D : hit2DVec)
                    integral += chargeIntegral(hit2D->PeakTime(),hit2D->PeakAmplitude(),hit2D->RMS(),1.,lowIndex,hiIndex);
                
                integral /= float(hit2DVec.size());
            }
            
//            hitCharge = std::min(hitCharge, fMaxHitCharge);
            integral = std::min(integral, fMaxHitCharge);
        
            float chgFactor = cst->fRecoQLow[geo::kCollection] + hitChiSqScale * integral;
        
            chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(chgFactor);
        
            colorToHitMap[chargeColorIdx].push_back(HitPosition()={{pos[0],pos[1],pos[2],err[3],err[3],err[5]}});
        }
    }
    
    for(auto& hitPair : colorToHitMap)
    {
        TPolyMarker3D& pm = view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotLarge, 0.25);
        for (const auto& hit : hitPair.second) pm.SetNextPoint(hit[0],hit[1],hit[2]);
    }

    return;
}

double SpacePoint3DDrawerHitCharge::chargeIntegral(double peakMean,
                                                   double peakAmp,
                                                   double peakWidth,
                                                   double areaNorm,
                                                   int    low,
                                                   int    hi) const
{
    double integral(0);
    
    for(int sigPos = low; sigPos < hi; sigPos++) integral += peakAmp * TMath::Gaus(double(sigPos)+0.5,peakMean,peakWidth);
    
    return integral;
}
    
DEFINE_ART_CLASS_TOOL(SpacePoint3DDrawerHitCharge)
}
