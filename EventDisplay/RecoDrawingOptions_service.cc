////////////////////////////////////////////////////////////////////////
/// \file RecoDrawingOption_plugin.cc
///
/// \version $Id: RecoDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// LArSoft includes
#include "EventDisplay/RecoDrawingOptions.h"

#include <iostream>

namespace evd {

//......................................................................
RecoDrawingOptions::RecoDrawingOptions(fhicl::ParameterSet const& pset, 
                                       art::ActivityRegistry& /* reg */)
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
    fDrawHits                  = pset.get< int                      >("DrawHits"                 );
    fDrawClusters    	       = pset.get< int                      >("DrawClusters"   	         );
    fDrawPFParticles           = pset.get< int                      >("DrawPFParticles"          );
    fDrawSpacePoints   	       = pset.get< int                      >("DrawSpacePoints"          );
    fDrawProngs      	       = pset.get< int                      >("DrawProngs"     	         );
    fDrawTracks      	       = pset.get< int                      >("DrawTracks"     	         );
    fDrawTrackTrajectoryPoints = pset.get< int                      >("DrawTrackTrajectoryPoints");
    fDrawTrackSpacePoints      = pset.get< int                      >("DrawTrackSpacePoints"     );
    fDrawShowers     	       = pset.get< int                      >("DrawShowers"    	         );
    fDrawVertices    	       = pset.get< int                      >("DrawVertices"   	     	 );
    fDrawOpFlashes             = pset.get< int                      >("DrawOpFlashes"            );
    fDrawSeeds     	           = pset.get< int                      >("DrawSeeds"    	     	 );
    fDrawBezierTracks          = pset.get< int                      >("DrawBezierTracks"         );
    fDrawEvents      	       = pset.get< int                      >("DrawEvents"     	     	 );
    fDraw2DEndPoints           = pset.get< int                      >("Draw2DEndPoints"	     	 );
    fDraw2DSlopeEndPoints      = pset.get< int                      >("Draw2DSlopeEndPoints"     );
    fSelectedHitColor	       = pset.get< int                      >("SelectedHitColor"         );
    fUseHitSelector            = pset.get< bool                     >("UseHitSelector"           );
    fSkeletonOnly              = pset.get< bool                     >("DrawSkeleton3DHitsOnly"   );
    fBestPCAAxisOnly           = pset.get< bool                     >("DrawBestPCAAxisOnly"      );
    fDrawTrackVertexAssns      = pset.get< bool                     >("DrawTrackVertexAssns"     );
    fHitLabels                 = pset.get< std::vector<std::string> >("HitModuleLabels"          );
    fSpacePointLabels 	       = pset.get< std::vector<std::string> >("SpacePointModuleLabels"	 );
    fProngLabels      	       = pset.get< std::vector<std::string> >("ProngModuleLabels"     	 );
    fEndPoint2DLabels 	       = pset.get< std::vector<std::string> >("EndPoint2DModuleLabels"	 );
    fClusterLabels    	       = pset.get< std::vector<std::string> >("ClusterModuleLabels"   	 );
    fPFParticleLabels          = pset.get< std::vector<std::string> >("PFParticleModuleLabels"   );
    fTrackLabels      	       = pset.get< std::vector<std::string> >("TrackModuleLabels"     	 );
    fShowerLabels     	       = pset.get< std::vector<std::string> >("ShowerModuleLabels"    	 );
    fVertexLabels     	       = pset.get< std::vector<std::string> >("VertexModuleLabels"    	 );
    fOpFlashLabels             = pset.get< std::vector<std::string> >("OpFlashModuleLabels"      );
    fSeedLabels       	       = pset.get< std::vector<std::string> >("SeedModuleLabels"      	 );
    fBezierTrackLabels         = pset.get< std::vector<std::string> >("BezierTrackModuleLabels"  );
    fTrkVtxTrackLabels         = pset.get< std::vector<std::string> >("TrkVtxTrackLabels"        );
    fTrkVtxCosmicLabels        = pset.get< std::vector<std::string> >("TrkVtxCosmicLabels"       );
    fTrkVtxFilterLabels        = pset.get< std::vector<std::string> >("TrkVtxFilterLabels"       );

    fEventLabels      	       = pset.get< std::vector<std::string> >("EventModuleLabels"     	 );
    fWireLabels       	       = pset.get< std::vector<std::string> >("WireModuleLabels"      	 );
    fColorProngsByLabel        = pset.get< int                      >("ColorProngsByLabel"       );
    fColorSpacePointsByChisq   = pset.get< int                      >("ColorSpacePointsByChisq"  );
    fCaloPSet                  = pset.get< fhicl::ParameterSet      >("CalorimetryAlgorithm"     );
    //   fSeedPSet = pset.get< fhicl::ParameterSet >("SeedAlgorithm");
    
    fCosmicTagLabels           = pset.get< std::vector<std::string> >("CosmicTagLabels", std::vector<std::string>() );
    fDrawCosmicTags            = pset.get< int                      >("DrawCosmicTags"           );
  }
  
}

namespace evd {

  DEFINE_ART_SERVICE(RecoDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
