////////////////////////////////////////////////////////////////////////
// $Id: RecoDrawingOption.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef RECODRAWINGOPTIONS_H
#define RECODRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace evd {
  
class RecoDrawingOptions
{
public:
    RecoDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~RecoDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset);

    int  fDrawHits;
    int  fDrawClusters;
    int  fDrawPFParticles;
    int  fDraw2DSlopeEndPoints;
    int  fDrawSpacePoints;
    int  fDrawProngs;
    int  fDrawTracks;
    int  fDrawOpFlashes;
    int  fDrawTrackTrajectoryPoints;
    int  fDrawTrackSpacePoints;
    int  fDrawShowers;
    int  fDrawVertices;
    int  fDrawEvents;
    int  fDraw2DEndPoints;
    int  fDrawSeeds;
    int  fDrawBezierTracks;
    int  fDrawCosmicTags;
    int  fSelectedHitColor;
    bool fUseHitSelector;
    bool fSkeletonOnly;
    bool fBestPCAAxisOnly;
    bool fDrawTrackVertexAssns;
    
    std::vector<std::string> fWireLabels;           ///< module labels that produced wires
    std::vector<std::string> fHitLabels;     		///< module labels that produced hits  		       
    std::vector<std::string> fEndPoint2DLabels; 	///< module labels that produced end point 2d objects  
    std::vector<std::string> fClusterLabels;    	///< module labels that produced clusters
    std::vector<std::string> fPFParticleLabels;     ///< module labels that produced PFParticles
    std::vector<std::string> fSpacePointLabels; 	///< module labels that produced space points
    std::vector<std::string> fProngLabels;   		///< module labels that produced prongs   	       
    std::vector<std::string> fTrackLabels;   		///< module labels that produced tracks   	       
    std::vector<std::string> fShowerLabels;  		///< module labels that produced showers  	       
    std::vector<std::string> fVertexLabels;  		///< module labels that produced vertices 	       
    std::vector<std::string> fEventLabels;   		///< module labels that produced events		       
    std::vector<std::string> fOpFlashLabels;        ///< module labels that produced events
    std::vector<std::string> fSeedLabels;       	///< module labels that produced events                
    std::vector<std::string> fBezierTrackLabels;    ///< module labels that produced events
    std::vector<std::string> fCosmicTagLabels;	    ///< module labels that produced cosmic tags
    std::vector<std::string> fTrkVtxTrackLabels;    ///< module labels that produced tracks (Track/Vertex module)
    std::vector<std::string> fTrkVtxCosmicLabels;   ///< module labels that tagged track as CR (Track/Vertex module)
    std::vector<std::string> fTrkVtxFilterLabels;   ///< module labels that filtered event (Track/Vertex module)


    ///\todo Why are calorimetry related drawing options in RecoDrawingOptions instead of a separate service?
    fhicl::ParameterSet      fCaloPSet;                 /// < parameterset for calorimetry algorithm 
    fhicl::ParameterSet      fSeedPSet;                 /// < parameterset for seed algorithm        


    int                      fColorProngsByLabel;       ///< Generate prong colors by label or id?
    int                      fColorSpacePointsByChisq;  ///< Generate space point colors by chisquare?

    double fFlashMinPE;                                 ///< Minimal PE for a flash to be displayed. 
    double fFlashTMin;                                  ///< Minimal time for a flash to be displayed.
    double fFlashTMax;                                  ///< Maximum time for a flash to be displayed. 
  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::RecoDrawingOptions, LEGACY)
#endif

