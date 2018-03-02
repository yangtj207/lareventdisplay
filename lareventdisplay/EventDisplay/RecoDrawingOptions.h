////////////////////////////////////////////////////////////////////////
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
#include "nutools/EventDisplayBase/Reconfigurable.h"
#include "canvas/Utilities/InputTag.h"

namespace evd {
  
  class RecoDrawingOptions : public evdb::Reconfigurable
{
public:
    explicit RecoDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~RecoDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset) ;

    int  fDrawHits;
    int  fDrawClusters;
    int  fDrawPFParticles;
    int  fDrawEdges;
    int  fDraw2DSlopeEndPoints;
    int  fDrawSpacePoints;
    int  fDrawProngs;
    int  fDrawTracks;
    int  fDrawOpFlashes;
    int  fDrawTrackTrajectoryPoints;
    int  fDrawTrackSegments;
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
    bool fDraw3DSpacePointHeatMap;
    
    std::vector<art::InputTag> fWireLabels;         ///< module labels that produced wires
    std::vector<art::InputTag> fHitLabels;     		///< module labels that produced hits
    std::vector<art::InputTag> fEndPoint2DLabels; 	///< module labels that produced end point 2d objects
    std::vector<art::InputTag> fClusterLabels;    	///< module labels that produced clusters
    std::vector<art::InputTag> fPFParticleLabels;   ///< module labels that produced PFParticles
    std::vector<art::InputTag> fEdgeLabels;         ///< module labels that produced Edge objects
    std::vector<art::InputTag> fSpacePointLabels; 	///< module labels that produced space points
    std::vector<art::InputTag> fProngLabels;   		///< module labels that produced prongs
    std::vector<art::InputTag> fTrackLabels;   		///< module labels that produced tracks
    std::vector<art::InputTag> fShowerLabels;  		///< module labels that produced showers
    std::vector<art::InputTag> fVertexLabels;  		///< module labels that produced vertices
    std::vector<art::InputTag> fEventLabels;   		///< module labels that produced events
    std::vector<art::InputTag> fOpFlashLabels;      ///< module labels that produced events
    std::vector<art::InputTag> fSeedLabels;       	///< module labels that produced events
    std::vector<art::InputTag> fBezierTrackLabels;  ///< module labels that produced events
    std::vector<art::InputTag> fCosmicTagLabels;	///< module labels that produced cosmic tags
    std::vector<art::InputTag> fTrkVtxTrackLabels;  ///< module labels that produced tracks (Track/Vertex module)
    std::vector<art::InputTag> fTrkVtxCosmicLabels; ///< module labels that tagged track as CR (Track/Vertex module)
    std::vector<art::InputTag> fTrkVtxFilterLabels; ///< module labels that filtered event (Track/Vertex module)

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

