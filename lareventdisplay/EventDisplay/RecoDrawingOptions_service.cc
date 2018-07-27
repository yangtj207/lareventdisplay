////////////////////////////////////////////////////////////////////////
/// \file RecoDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes

/// LArSoft includes
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include <iostream>

namespace evd {

//......................................................................
RecoDrawingOptions::RecoDrawingOptions(fhicl::ParameterSet const& pset, 
                                       art::ActivityRegistry& /* reg */)
  : evdb::Reconfigurable{pset}
{
    this->reconfigure(pset);
}

//......................................................................
RecoDrawingOptions::~RecoDrawingOptions() 
{
}

//......................................................................
void RecoDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
{
    fDrawHits                  = pset.get< int                        >("DrawHits"                 );
    fDrawClusters    	       = pset.get< int                        >("DrawClusters"   	       );
    fDrawSlices = pset.get<int>("DrawSlices", 0);
    fDrawPFParticles           = pset.get< int                        >("DrawPFParticles"          );
    fDrawSpacePoints   	       = pset.get< int                        >("DrawSpacePoints"          );
    fDrawProngs      	       = pset.get< int                        >("DrawProngs"     	       );
    fDrawTracks      	       = pset.get< int                        >("DrawTracks"     	       );
    fDrawTrackTrajectoryPoints = pset.get< int                        >("DrawTrackTrajectoryPoints");
    fDrawTrackSegments         = pset.get< int                        >("DrawTrackSegments"        );
    fDrawTrackSpacePoints      = pset.get< int                        >("DrawTrackSpacePoints"     );
    fDrawShowers     	       = pset.get< int                        >("DrawShowers"    	       );
    fDrawVertices    	       = pset.get< int                        >("DrawVertices"   	       );
    fDrawOpFlashes             = pset.get< int                        >("DrawOpFlashes"            );
    fDrawSeeds     	           = pset.get< int                        >("DrawSeeds"    	     	   );
    fDrawBezierTracks          = pset.get< int                        >("DrawBezierTracks"         );
    fDrawEvents      	       = pset.get< int                        >("DrawEvents"     	       );
    fDraw2DEndPoints           = pset.get< int                        >("Draw2DEndPoints"	       );
    fDraw2DSlopeEndPoints      = pset.get< int                        >("Draw2DSlopeEndPoints"     );
    fSelectedHitColor	       = pset.get< int                        >("SelectedHitColor"         );
    fUseHitSelector            = pset.get< bool                       >("UseHitSelector"           );
    fSkeletonOnly              = pset.get< bool                       >("DrawSkeleton3DHitsOnly"   );
    fBestPCAAxisOnly           = pset.get< bool                       >("DrawBestPCAAxisOnly"      );
    fDrawTrackVertexAssns      = pset.get< bool                       >("DrawTrackVertexAssns"     );
    fDraw3DSpacePointHeatMap   = pset.get< bool                       >("Draw3DSpacePointHeatMap"  );
    fHitLabels                 = pset.get< std::vector<art::InputTag> >("HitModuleLabels"          );
    if(pset.has_key("SliceModuleLabels")) fSliceLabels = pset.get< std::vector<art::InputTag> >("SliceModuleLabels");
    fSpacePointLabels 	       = pset.get< std::vector<art::InputTag> >("SpacePointModuleLabels"   );
    fProngLabels      	       = pset.get< std::vector<art::InputTag> >("ProngModuleLabels"        );
    fEndPoint2DLabels 	       = pset.get< std::vector<art::InputTag> >("EndPoint2DModuleLabels"   );
    fClusterLabels    	       = pset.get< std::vector<art::InputTag> >("ClusterModuleLabels"      );
    fPFParticleLabels          = pset.get< std::vector<art::InputTag> >("PFParticleModuleLabels"   );
    fTrackLabels      	       = pset.get< std::vector<art::InputTag> >("TrackModuleLabels"        );
    fShowerLabels     	       = pset.get< std::vector<art::InputTag> >("ShowerModuleLabels"       );
    fVertexLabels     	       = pset.get< std::vector<art::InputTag> >("VertexModuleLabels"       );
    fOpFlashLabels             = pset.get< std::vector<art::InputTag> >("OpFlashModuleLabels"      );
    fSeedLabels       	       = pset.get< std::vector<art::InputTag> >("SeedModuleLabels"         );
    fBezierTrackLabels         = pset.get< std::vector<art::InputTag> >("BezierTrackModuleLabels"  );
    fTrkVtxTrackLabels         = pset.get< std::vector<art::InputTag> >("TrkVtxTrackLabels"        );
    fTrkVtxCosmicLabels        = pset.get< std::vector<art::InputTag> >("TrkVtxCosmicLabels"       );
    fTrkVtxFilterLabels        = pset.get< std::vector<art::InputTag> >("TrkVtxFilterLabels"       );

    fEventLabels      	       = pset.get< std::vector<art::InputTag> >("EventModuleLabels"        );
    fWireLabels       	       = pset.get< std::vector<art::InputTag> >("WireModuleLabels"         );
    fColorProngsByLabel        = pset.get< int                        >("ColorProngsByLabel"       );
    fColorSpacePointsByChisq   = pset.get< int                        >("ColorSpacePointsByChisq"  );
    fCaloPSet                  = pset.get< fhicl::ParameterSet        >("CalorimetryAlgorithm"     );
    //   fSeedPSet = pset.get< fhicl::ParameterSet >("SeedAlgorithm");
    
    fCosmicTagLabels           = pset.get< std::vector<art::InputTag> >("CosmicTagLabels", std::vector<art::InputTag>() );
    fDrawCosmicTags            = pset.get< int                        >("DrawCosmicTags"           );
    fFlashMinPE                = pset.get< double                     >("FlashMinPE", 0.0          );
    fFlashTMin                 = pset.get< double                     >("FlashTMin", -1e9          );
    fFlashTMax                 = pset.get< double                     >("FlashTMax", 1e9           );
  }
  
}

namespace evd {

  DEFINE_ART_SERVICE(RecoDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
