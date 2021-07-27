////////////////////////////////////////////////////////////////////////
/// \file RecoDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

/// LArSoft includes
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"


namespace evd {

//......................................................................
RecoDrawingOptions::RecoDrawingOptions(fhicl::ParameterSet const& pset)
  : evdb::Reconfigurable{pset}
{
    this->reconfigure(pset);
}

//......................................................................
void RecoDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
{
    fDrawHits                  = pset.get< int                        >("DrawHits"                 );
    fDrawClusters    	       = pset.get< int                        >("DrawClusters"   	       );
    fDrawSlices                = pset.get< int                        >("DrawSlices",             0);
    fDrawSliceSpacePoints      = pset.get< int                        >("DrawSliceSpacePoints",   0);
    fDrawPFParticles           = pset.get< int                        >("DrawPFParticles"          );
    fDrawEdges                 = pset.get< int                        >("DrawEdges"                );
    fDrawSpacePoints   	       = pset.get< int                        >("DrawSpacePoints"          );
    fDrawProngs      	       = pset.get< int                        >("DrawProngs"     	       );
    fDrawTracks      	       = pset.get< int                        >("DrawTracks"     	       );
    fDrawTrackTrajectoryPoints = pset.get< int                        >("DrawTrackTrajectoryPoints");
    fDrawTrackSegments         = pset.get< int                        >("DrawTrackSegments"        );
    fDrawTrackSpacePoints      = pset.get< int                        >("DrawTrackSpacePoints"     );
    fDrawShowers     	       = pset.get< int                        >("DrawShowers"    	       );
    fDrawVertices    	       = pset.get< int                        >("DrawVertices"   	       );
    fDrawOpHits                = pset.get< int                        >("DrawOpHits"               );
    fDrawOpFlashes             = pset.get< int                        >("DrawOpFlashes"            );
    fDrawSeeds     	           = pset.get< int                        >("DrawSeeds"    	     	   );
    fDrawEvents      	       = pset.get< int                        >("DrawEvents"     	       );
    fDraw2DEndPoints           = pset.get< int                        >("Draw2DEndPoints"	       );
    fDraw2DSlopeEndPoints      = pset.get< int                        >("Draw2DSlopeEndPoints"     );
    fSelectedHitColor	       = pset.get< int                        >("SelectedHitColor"         );
    fUseHitSelector            = pset.get< bool                       >("UseHitSelector"           );
    fSkeletonOnly              = pset.get< bool                       >("DrawSkeleton3DHitsOnly"   );
    fBestPCAAxisOnly           = pset.get< bool                       >("DrawBestPCAAxisOnly"      );
    fDrawTrackVertexAssns      = pset.get< bool                       >("DrawTrackVertexAssns"     );
    fDraw3DSpacePoints         = pset.get< bool                       >("Draw3DSpacePoints"        );
    fDraw3DSpacePointHeatMap   = pset.get< bool                       >("Draw3DSpacePointHeatMap"  );
    fDraw3DEdges               = pset.get< bool                       >("Draw3DEdges"              );
    fDraw3DPCAAxes             = pset.get< bool                       >("Draw3DPCAAxes"            );
    fDrawAllWireIDs            = pset.get< bool                       >("DrawAllWireIDs"           );
    fHitLabels                 = pset.get< std::vector<art::InputTag> >("HitModuleLabels"          );
    if(pset.has_key("SliceModuleLabels")) fSliceLabels = pset.get< std::vector<art::InputTag> >("SliceModuleLabels");
    fSpacePointLabels 	       = pset.get< std::vector<art::InputTag> >("SpacePointModuleLabels"   );
    fProngLabels      	       = pset.get< std::vector<art::InputTag> >("ProngModuleLabels"        );
    fEndPoint2DLabels 	       = pset.get< std::vector<art::InputTag> >("EndPoint2DModuleLabels"   );
    fClusterLabels    	       = pset.get< std::vector<art::InputTag> >("ClusterModuleLabels"      );
    fPFParticleLabels          = pset.get< std::vector<art::InputTag> >("PFParticleModuleLabels"   );
    fEdgeLabels                = pset.get< std::vector<art::InputTag> >("EdgeModuleLabels"         );
    fExtremePointLabels        = pset.get< std::vector<art::InputTag> >("ExtremeModuleLabels"      );
    fTrackLabels      	       = pset.get< std::vector<art::InputTag> >("TrackModuleLabels"        );
    fShowerLabels     	       = pset.get< std::vector<art::InputTag> >("ShowerModuleLabels"       );
    fVertexLabels     	       = pset.get< std::vector<art::InputTag> >("VertexModuleLabels"       );
    fOpHitLabels               = pset.get< std::vector<art::InputTag> >("OpHitModuleLabels"        );
    fOpFlashLabels             = pset.get< std::vector<art::InputTag> >("OpFlashModuleLabels"      );
    fSeedLabels       	       = pset.get< std::vector<art::InputTag> >("SeedModuleLabels"         );
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

    fHitDrawerParams           = pset.get< fhicl::ParameterSet        >("HitDrawer"                );
    fWireDrawerParams          = pset.get< fhicl::ParameterSet        >("WireDrawer"               );

    fSpacePointDrawerParams    = pset.get< fhicl::ParameterSet        >("SpacePointDrawer"         );
    fAllSpacePointDrawerParams = pset.get< fhicl::ParameterSet        >("AllSpacePointDrawer"      );
    
    f3DDrawerParams            = pset.get< fhicl::ParameterSet        >("Reco3DDrawers"            );
  }

  DEFINE_ART_SERVICE(RecoDrawingOptions)
} // namespace evd
////////////////////////////////////////////////////////////////////////
