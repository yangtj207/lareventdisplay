/// \file    RecoBaseDrawer.h
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  messier@indiana.edu
/// \version $Id: RecoBaseDrawer.h,v 1.3 2010/11/11 22:47:20 p-novaart Exp $
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOBASEDRAWER_H

#include <vector>

#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Principal/View.h"

#ifdef __ROOTCLING__

namespace art { 
    class Event;
    class ServiceHandle;
}

#else
#include "art/Framework/Services/Registry/ServiceHandle.h"
#endif

#include "lareventdisplay/EventDisplay/OrthoProj.h"
#include "lardataobj/RecoBase/SpacePoint.h"

class TVector3;
class TH1F;
namespace evdb { 
    class View2D;
    class View3D;
}

namespace geo   { class Geometry; }

namespace recob { 
    class Hit;
    class Cluster;
    class PCAxis;
    class PFParticle;
    class EndPoint2D;
    class Prong;
    class Track;
    class Shower;
    class Wire;
    class Vertex;
    class Seed;
    class Event;
    class OpFlash;
}

namespace anab {
    class CosmicTag;
}

namespace util {
    class LArPropertiesService;
    class DetectorProperties;
}

namespace evd {

class ColorDrawingOptions;
class RawDrawingOptions;
class RecoDrawingOptions;

/// Aid in the rendering of RecoBase objects
class RecoBaseDrawer
{
public:
    RecoBaseDrawer();
    ~RecoBaseDrawer();

public:

    void Wire2D(const art::Event& evt,
                evdb::View2D*     view,
		        unsigned int      plane);
    void Hit2D(const art::Event& evt,
	           evdb::View2D*     view,
	           unsigned int      plane);
    void Hit2D(std::vector<const recob::Hit*> hits,
	           int                            color,
	           evdb::View2D*                  view,
               bool                           drawConnectingLines = false,
               int                            lineWidth = 1);
    void EndPoint2D(const art::Event& evt,
		            evdb::View2D*     view,
                    unsigned int      plane);
    void OpFlash2D(const art::Event& evt,
		           evdb::View2D*     view,
		           unsigned int      plane);

    void Seed2D(const art::Event& evt,
		        evdb::View2D*     view,
		        unsigned int      plane);

    void BezierTrack2D(const art::Event& evt,
		               evdb::View2D*     view,
		               unsigned int      plane);


    void Draw2DSlopeEndPoints(double        xStart,
                              double        yStart,
                              double        xEnd,
                              double        yEnd,
                              double        slope,
                              int           color,
                              evdb::View2D* view);
    void Draw2DSlopeEndPoints(double        x,
                              double        y,
                              double        slope,
                              int           color,
                              evdb::View2D* view);
    void Draw2DSlopeEndPoints(double        x,
                              double        y,
                              double        cosx,
			      double        cosy,
                              int           color,
                              evdb::View2D* view);
    void Cluster2D(const art::Event& evt,
		           evdb::View2D*     view,
		           unsigned int      plane);
    void Prong2D(const art::Event&        evt,
		         evdb::View2D*            view,
		         unsigned int             plane);
    void DrawTrackVertexAssns2D(const art::Event& evt,
                                evdb::View2D*     view,
                                unsigned int      plane);
    void DrawProng2D(std::vector<const recob::Hit*>& hits,
                     evdb::View2D*                   view,
                     unsigned int                    plane,
                     TVector3                 const& startPos,
                     TVector3                 const& startDir,
                     int                             id,
		     float cscore = -5);
    void DrawTrack2D(std::vector<const recob::Hit*>& hits,
                     evdb::View2D*                   view,
                     unsigned int                    plane,
                     const recob::Track*             track,
                     int                             color,
                     int                             lineWidth);
    void Vertex2D(const art::Event& evt,
		          evdb::View2D*     view,
		          unsigned int      plane);

    void Event2D(const art::Event& evt,
		         evdb::View2D*     view,
		         unsigned int      plane);

    void SpacePoint3D(const art::Event& evt,
		              evdb::View3D*     view);
    void PFParticle3D(const art::Event& evt,
		              evdb::View3D*     view);
    void DrawPFParticle3D(const art::Ptr<recob::PFParticle>&       pfPart,
                          const art::PtrVector<recob::PFParticle>& pfParticleVec,
                          const art::FindMany<recob::SpacePoint>&  spacePointAssnsVec,
                          const art::FindMany<recob::Track>&       trackAssnVec,
                          const art::FindMany<recob::PCAxis>&      pcAxisAssnVec,
                          const art::FindMany<anab::CosmicTag>&    cosmicTagAssnVec,
                          int                                      depth,
                          evdb::View3D*                            view);
    void Prong3D(const art::Event& evt,
		         evdb::View3D*     view);
    void DrawSpacePoint3D(const std::vector<const recob::SpacePoint*>& spts,
                          evdb::View3D*                                view,
                          int                                          color,
                          int                                          marker = 3,
                          int                                          size = 1);
    void DrawTrack3D(const recob::Track& track,
		             evdb::View3D*       view,
                     int                 color,
                     int                 marker = 1,
                     int                 size = 2);
    void DrawShower3D(const recob::Shower& shower,
		              int                  color,
                      evdb::View3D*        view);
    void Seed3D(const art::Event& evt,
		        evdb::View3D*     view);
    void BezierTrack3D(const art::Event& evt,
		               evdb::View3D*     view);
    void Vertex3D(const art::Event& evt,
		          evdb::View3D*     view);
    void Event3D(const art::Event& evt,
                evdb::View3D*     view);
    void OpFlashOrtho(const art::Event& evt,
		      evd::OrthoProj_t  proj,
		      evdb::View2D*     view);
    void VertexOrtho(const art::PtrVector<recob::Vertex>& vertex,
				  evd::OrthoProj_t proj,
				  evdb::View2D* view,
				  int marker);
    void VertexOrtho(const art::Event& evt,
		     evd::OrthoProj_t  proj,
		     evdb::View2D*     view);
    void SpacePointOrtho(const art::Event& evt,
			             evd::OrthoProj_t  proj,
			             double            msize,
			             evdb::View2D*     view);
    void PFParticleOrtho(const art::Event& evt,
                         evd::OrthoProj_t  proj,
                         double            msize,
		                 evdb::View2D*     view);
    void DrawPFParticleOrtho(const art::Ptr<recob::PFParticle>&       pfPart,
                             const art::PtrVector<recob::PFParticle>& pfParticleVec,
                             const art::FindMany<recob::SpacePoint>&  spacePointAssnsVec,
                             const art::FindMany<recob::PCAxis>&      pcAxisAssnVec,
                             int                                      depth,
                             evd::OrthoProj_t                         proj,
                             evdb::View2D*                            view);
    void ProngOrtho(const art::Event& evt,
		            evd::OrthoProj_t  proj,
		            double            msize,
		            evdb::View2D*     view);
    void DrawSpacePointOrtho(const std::vector<const recob::SpacePoint*>& spts,
			                 int                 color,
                             evd::OrthoProj_t    proj,
			                 double              msize,
                             evdb::View2D*       view,
			     int mode = 0); ///< 0: track, 1: shower
    void DrawProngOrtho(const recob::Prong& prong,
			            int                 color,
			            evd::OrthoProj_t    proj,
			            double              msize,
			            evdb::View2D*       view);
    void DrawTrackOrtho(const recob::Track& track,
			            int                 color,
			            evd::OrthoProj_t    proj,
			            double              msize,
			            evdb::View2D*       view);
    void DrawShowerOrtho(const recob::Shower& shower,
			             int                 color,
			             evd::OrthoProj_t    proj,
			             double              msize,
                         evdb::View2D*       view);
    void SeedOrtho(const art::Event& evt,
		           evd::OrthoProj_t  proj,
		           evdb::View2D*     view);

    void FillTQHisto(const art::Event&    evt, 
		             unsigned int         plane,
                     unsigned int         wire,
		             TH1F*                histo,
		             std::vector<double>& hstart,
                     std::vector<double>& hend,
                     std::vector<double>& hitamplitudes,
                     std::vector<double>& hpeaktimes);

    void FillQHisto(const art::Event& evt,
		            unsigned int      plane,
		            TH1F*             histo);

    int GetRegionOfInterest(int  plane,
			                int& minw,
                            int& maxw,
			                int& mint,
			                int& maxt);
    
    void GetChargeSum(int     plane,
		              double& charge,
		              double& convcharge);
    
  private:
    void GetClusterOutlines(std::vector<const recob::Hit*>& hits,
			                std::vector<double>&            tpts,
			                std::vector<double>&      	    wpts,
			                unsigned int              	    plane);
    int GetWires(const art::Event&            evt,
                 const std::string&           which,
		         art::PtrVector<recob::Wire>& wires);
    int GetHits(const art::Event&               evt,
                const std::string&              which,
		        std::vector<const recob::Hit*>& hits,
		        unsigned int                    plane);
    int GetClusters(const art::Event&               evt,
		            const std::string&              which,
		            art::PtrVector<recob::Cluster>& clust);
    int GetPFParticles(const art::Event&                  evt,
		               const std::string&                 which,
		               art::PtrVector<recob::PFParticle>& pfpart);
    int GetEndPoint2D(const art::Event&                  evt, 
		              const std::string&                 which,
		              art::PtrVector<recob::EndPoint2D>& ep2d);
    int GetSpacePoints(const art::Event&               evt,
		               const std::string&              which,
		               std::vector<const recob::SpacePoint*>& spts);

    int GetTracks(const art::Event&        evt,
		          const std::string&       which,
		          art::View<recob::Track>& track);

    int GetShowers(const art::Event&        evt,
		           const std::string&        which,
		           art::View<recob::Shower>& shower);

    int GetVertices(const art::Event&              evt,
		            const art::InputTag&           which,
		            art::PtrVector<recob::Vertex>& vertex);

    int GetSeeds(const art::Event&            evt,
                 const std::string&           which,
		         art::PtrVector<recob::Seed>& seed);

    int GetBezierTracks(const art::Event&             evt,
			            const std::string&            which,
			            art::PtrVector<recob::Track>& btbs);


    int GetOpFlashes(const art::Event&               evt,
                     const std::string&              which,
		             art::PtrVector<recob::OpFlash>& opflash);

    int GetEvents(const art::Event&             evt,
		          const std::string&            which,
		          art::PtrVector<recob::Event>& event);

        
  private:

    std::vector<int>          fWireMin;         ///< lowest wire in interesting region for each plane
    std::vector<int>          fWireMax;         ///< highest wire in interesting region for each plane
    std::vector<int>          fTimeMin;         ///< lowest time in interesting region for each plane
    std::vector<int>          fTimeMax;         ///< highest time in interesting region for each plane
    
    std::vector<double>       fRawCharge;       ///< Sum of Raw Charge
    std::vector<double>       fConvertedCharge; ///< Sum of Charge Converted using Birks' formula
    
  };
}

#endif
////////////////////////////////////////////////////////////////////////
