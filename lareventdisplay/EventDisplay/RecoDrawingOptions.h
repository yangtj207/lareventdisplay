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
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nuevdb/EventDisplayBase/Reconfigurable.h"
#include "canvas/Utilities/InputTag.h"

namespace evd {
  
  class RecoDrawingOptions : public evdb::Reconfigurable
{
public:
    explicit RecoDrawingOptions(fhicl::ParameterSet const& pset);
    
    void reconfigure(fhicl::ParameterSet const& pset) ;

    int  fDrawHits;
    int  fDrawClusters;
    int  fDrawSlices;
    int  fDrawSliceSpacePoints;
    int  fDrawPFParticles;
    int  fDrawEdges;
    int  fDraw2DSlopeEndPoints;
    int  fDrawSpacePoints;
    int  fDrawProngs;
    int  fDrawTracks;
    int  fDrawOpHits;
    int  fDrawOpFlashes;
    int  fDrawTrackTrajectoryPoints;
    int  fDrawTrackSegments;
    int  fDrawTrackSpacePoints;
    int  fDrawShowers;
    int  fDrawVertices;
    int  fDrawEvents;
    int  fDraw2DEndPoints;
    int  fDrawSeeds;
    int  fDrawCosmicTags;
    int  fSelectedHitColor;
    bool fUseHitSelector;
    bool fSkeletonOnly;
    bool fBestPCAAxisOnly;
    bool fDrawTrackVertexAssns;
    bool fDraw3DSpacePoints;
    bool fDraw3DSpacePointHeatMap;
    bool fDraw3DEdges;
    bool fDraw3DPCAAxes;
    bool fDrawAllWireIDs;
    
    std::vector<art::InputTag> fWireLabels;                 ///< module labels that produced wires
    std::vector<art::InputTag> fHitLabels;                  ///< module labels that produced hits
    std::vector<art::InputTag> fSliceLabels;                ///< module labels that produced slices
    std::vector<art::InputTag> fEndPoint2DLabels;           ///< module labels that produced end point 2d objects
    std::vector<art::InputTag> fClusterLabels;              ///< module labels that produced clusters
    std::vector<art::InputTag> fPFParticleLabels;           ///< module labels that produced PFParticles
    std::vector<art::InputTag> fEdgeLabels;                 ///< module labels that produced Edge objects
    std::vector<art::InputTag> fExtremePointLabels;         ///< module labels that produced Extreme Points
    std::vector<art::InputTag> fSpacePointLabels;           ///< module labels that produced space points
    std::vector<art::InputTag> fProngLabels;                ///< module labels that produced prongs
    std::vector<art::InputTag> fTrackLabels;                ///< module labels that produced tracks
    std::vector<art::InputTag> fShowerLabels;               ///< module labels that produced showers
    std::vector<art::InputTag> fVertexLabels;               ///< module labels that produced vertices
    std::vector<art::InputTag> fEventLabels;                ///< module labels that produced events
    std::vector<art::InputTag> fOpHitLabels;                ///< module labels that produced events
    std::vector<art::InputTag> fOpFlashLabels;              ///< module labels that produced events
    std::vector<art::InputTag> fSeedLabels;                 ///< module labels that produced events
    std::vector<art::InputTag> fCosmicTagLabels;            ///< module labels that produced cosmic tags
    std::vector<art::InputTag> fTrkVtxTrackLabels;          ///< module labels that produced tracks (Track/Vertex module)
    std::vector<art::InputTag> fTrkVtxCosmicLabels;         ///< module labels that tagged track as CR (Track/Vertex module)
    std::vector<art::InputTag> fTrkVtxFilterLabels;         ///< module labels that filtered event (Track/Vertex module)

    ///\todo Why are calorimetry related drawing options in RecoDrawingOptions instead of a separate service?
    fhicl::ParameterSet        fCaloPSet;                   /// < parameterset for calorimetry algorithm
    fhicl::ParameterSet        fSeedPSet;                   /// < parameterset for seed algorithm
   
    int                        fColorProngsByLabel;         ///< Generate prong colors by label or id?
    int                        fColorSpacePointsByChisq;    ///< Generate space point colors by chisquare?
   
    double                     fFlashMinPE;                 ///< Minimal PE for a flash to be displayed.
    double                     fFlashTMin;                  ///< Minimal time for a flash to be displayed.
    double                     fFlashTMax;                  ///< Maximum time for a flash to be displayed.
    
    fhicl::ParameterSet        fHitDrawerParams;            ///< FHICL parameters for the hit drawing
    fhicl::ParameterSet        fWireDrawerParams;           ///< FHICL parameters for the wire drawing
    
    fhicl::ParameterSet        fSpacePointDrawerParams;     ///< FHICL parameters for SpacePoint drawing
    fhicl::ParameterSet        fAllSpacePointDrawerParams;  ///< FHICL parameters for SpacePoint drawing
    
    fhicl::ParameterSet        f3DDrawerParams;             ///< FHICL paramegers for the 3D drawers
  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::RecoDrawingOptions, LEGACY)
#endif
