/// \file    RecoBaseDrawer.h
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  messier@indiana.edu
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOBASEDRAWER_H

#include <vector>

#include "art/Framework/Principal/DataViewImpl.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Principal/fwd.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lareventdisplay/EventDisplay/OrthoProj.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

class TVector3;
class TH1F;
namespace evdb {
  class View2D;
  class View3D;
}

namespace recob {
  class Edge;
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

namespace evdb_tool {
  class ISpacePoints3D;
}

namespace evd {

  /// Aid in the rendering of RecoBase objects
  class RecoBaseDrawer {
  public:
    RecoBaseDrawer();
    ~RecoBaseDrawer();

  public:
    void Wire2D(const art::Event& evt, evdb::View2D* view, unsigned int plane);
    int Hit2D(const art::Event& evt,
              detinfo::DetectorPropertiesData const& detProp,
              evdb::View2D* view,
              unsigned int plane);
    int Hit2D(std::vector<const recob::Hit*> hits,
              int color,
              evdb::View2D* view,
              bool allWireIds,
              bool drawConnectingLines = false,
              int lineWidth = 1);
    int Hit2D(std::vector<const recob::Hit*> hits, evdb::View2D* view, float cosmicscore);

    void EndPoint2D(const art::Event& evt, evdb::View2D* view, unsigned int plane);
    void OpFlash2D(const art::Event& evt,
                   detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   evdb::View2D* view,
                   unsigned int plane);

    void Seed2D(const art::Event& evt,
                detinfo::DetectorPropertiesData const& detProp,
                evdb::View2D* view,
                unsigned int plane);

    void Draw2DSlopeEndPoints(double xStart,
                              double yStart,
                              double xEnd,
                              double yEnd,
                              double slope,
                              int color,
                              evdb::View2D* view);
    void Draw2DSlopeEndPoints(double x, double y, double slope, int color, evdb::View2D* view);
    void Draw2DSlopeEndPoints(double x,
                              double y,
                              double cosx,
                              double cosy,
                              int color,
                              evdb::View2D* view);
    void Slice2D(const art::Event& evt,
                 detinfo::DetectorPropertiesData const& detProp,
                 evdb::View2D* view,
                 unsigned int plane);
    void Cluster2D(const art::Event& evt,
                   detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   evdb::View2D* view,
                   unsigned int plane);
    void Prong2D(const art::Event& evt,
                 detinfo::DetectorClocksData const& clockData,
                 detinfo::DetectorPropertiesData const& detProp,
                 evdb::View2D* view,
                 unsigned int plane);
    void DrawTrackVertexAssns2D(const art::Event& evt,
                                detinfo::DetectorClocksData const& clockData,
                                detinfo::DetectorPropertiesData const& detProp,
                                evdb::View2D* view,
                                unsigned int plane);
    void DrawProng2D(detinfo::DetectorPropertiesData const& detProp,
                     std::vector<const recob::Hit*>& hits,
                     evdb::View2D* view,
                     unsigned int plane,
                     TVector3 const& startPos,
                     TVector3 const& startDir,
                     int id,
                     float cscore = -5);
    void DrawTrack2D(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     std::vector<const recob::Hit*>& hits,
                     evdb::View2D* view,
                     unsigned int plane,
                     const recob::Track* track,
                     int color,
                     int lineWidth);
    void Vertex2D(const art::Event& evt,
                  detinfo::DetectorPropertiesData const& detProp,
                  evdb::View2D* view,
                  unsigned int plane);

    void Event2D(const art::Event& evt, evdb::View2D* view, unsigned int plane);

    void SpacePoint3D(const art::Event& evt, evdb::View3D* view);
    void PFParticle3D(const art::Event& evt, evdb::View3D* view);
    void DrawPFParticle3D(const art::Ptr<recob::PFParticle>& pfPart,
                          const art::PtrVector<recob::PFParticle>& pfParticleVec,
                          const std::vector<art::Ptr<recob::SpacePoint>>& spacePointVec,
                          const art::FindManyP<recob::Edge>& edgeAssnsVec,
                          const art::FindManyP<recob::SpacePoint>& spacePointAssnsVec,
                          const art::FindManyP<recob::SpacePoint>& edgeSPAssnVec,
                          const art::FindManyP<recob::Hit>& spHitAssnVec,
                          const art::FindMany<recob::Track>& trackAssnVec,
                          const art::FindMany<recob::PCAxis>& pcAxisAssnVec,
                          const art::FindMany<anab::CosmicTag>& cosmicTagAssnVec,
                          int depth,
                          evdb::View3D* view);
    void Edge3D(const art::Event& evt, evdb::View3D* view);
    void Prong3D(const art::Event& evt, evdb::View3D* view);
    void DrawTrack3D(const recob::Track& track,
                     evdb::View3D* view,
                     int color,
                     int marker = 1,
                     float size = 2.);
    void DrawShower3D(const recob::Shower& shower, int color, evdb::View3D* view);
    void Seed3D(const art::Event& evt, evdb::View3D* view);
    void Vertex3D(const art::Event& evt, evdb::View3D* view);
    void Event3D(const art::Event& evt, evdb::View3D* view);
    void Slice3D(const art::Event& evt, evdb::View3D* view);
    void OpFlashOrtho(const art::Event& evt,
                      detinfo::DetectorClocksData const& clockData,
                      detinfo::DetectorPropertiesData const& detProp,
                      evd::OrthoProj_t proj,
                      evdb::View2D* view);
    void VertexOrtho(const art::PtrVector<recob::Vertex>& vertex,
                     evd::OrthoProj_t proj,
                     evdb::View2D* view,
                     int marker);
    void VertexOrtho(const art::Event& evt, evd::OrthoProj_t proj, evdb::View2D* view);
    void SpacePointOrtho(const art::Event& evt,
                         evd::OrthoProj_t proj,
                         double msize,
                         evdb::View2D* view);
    void PFParticleOrtho(const art::Event& evt,
                         evd::OrthoProj_t proj,
                         double msize,
                         evdb::View2D* view);
    void DrawPFParticleOrtho(const art::Ptr<recob::PFParticle>& pfPart,
                             const art::PtrVector<recob::PFParticle>& pfParticleVec,
                             const art::FindMany<recob::SpacePoint>& spacePointAssnsVec,
                             const art::FindMany<recob::PCAxis>& pcAxisAssnVec,
                             int depth,
                             evd::OrthoProj_t proj,
                             evdb::View2D* view);
    void ProngOrtho(const art::Event& evt, evd::OrthoProj_t proj, double msize, evdb::View2D* view);
    void DrawSpacePointOrtho(std::vector<art::Ptr<recob::SpacePoint>>& spts,
                             int color,
                             evd::OrthoProj_t proj,
                             double msize,
                             evdb::View2D* view,
                             int mode = 0); ///< 0: track, 1: shower
    void DrawProngOrtho(const recob::Prong& prong,
                        int color,
                        evd::OrthoProj_t proj,
                        double msize,
                        evdb::View2D* view);
    void DrawTrackOrtho(const recob::Track& track,
                        int color,
                        evd::OrthoProj_t proj,
                        double msize,
                        evdb::View2D* view);
    void DrawShowerOrtho(const recob::Shower& shower,
                         int color,
                         evd::OrthoProj_t proj,
                         double msize,
                         evdb::View2D* view);
    void SeedOrtho(const art::Event& evt, evd::OrthoProj_t proj, evdb::View2D* view);

    void FillTQHisto(const art::Event& evt, unsigned int plane, unsigned int wire, TH1F* histo);

    void FillQHisto(const art::Event& evt, unsigned int plane, TH1F* histo);

    void FillTQHistoDP(const art::Event& evt,
                       unsigned int plane,
                       unsigned int wire,
                       TH1F* histo,
                       std::vector<double>& htau1,
                       std::vector<double>& htau2,
                       std::vector<double>& hitamplitudes,
                       std::vector<double>& hpeaktimes,
                       std::vector<int>& hstartT,
                       std::vector<int>& hendT,
                       std::vector<int>& hNMultiHit,
                       std::vector<int>& hLocalHitIndex);

    int GetRegionOfInterest(int plane, int& minw, int& maxw, int& mint, int& maxt);

    void GetChargeSum(int plane, double& charge, double& convcharge);

    //double EvalExpoFit(double x,
    //	       double tau1,
    //	       double tau2,
    //	       double amplitude,
    //	       double peaktime);

    //double EvalMultiExpoFit(double x,
    //		    int HitNumber,
    //		    int NHits,
    //		    std::vector<double> tau1,
    //		    std::vector<double> tau2,
    //		    std::vector<double> amplitude,
    //		    std::vector<double> peaktime);

  private:
    void GetClusterOutlines(std::vector<const recob::Hit*>& hits,
                            std::vector<double>& tpts,
                            std::vector<double>& wpts,
                            unsigned int plane);
    int GetWires(const art::Event& evt,
                 const art::InputTag& which,
                 art::PtrVector<recob::Wire>& wires);
    int GetHits(const art::Event& evt,
                const art::InputTag& which,
                std::vector<const recob::Hit*>& hits,
                unsigned int plane);
    int GetSlices(const art::Event& evt,
                  const art::InputTag& which,
                  art::PtrVector<recob::Slice>& slices);
    int GetClusters(const art::Event& evt,
                    const art::InputTag& which,
                    art::PtrVector<recob::Cluster>& clust);
    int GetPFParticles(const art::Event& evt,
                       const art::InputTag& which,
                       art::PtrVector<recob::PFParticle>& pfpart);
    int GetEndPoint2D(const art::Event& evt,
                      const art::InputTag& which,
                      art::PtrVector<recob::EndPoint2D>& ep2d);
    int GetSpacePoints(const art::Event& evt,
                       const art::InputTag& which,
                       std::vector<art::Ptr<recob::SpacePoint>>& spts);
    int GetEdges(const art::Event& evt,
                 const art::InputTag& which,
                 std::vector<art::Ptr<recob::Edge>>& edges);

    int GetTracks(const art::Event& evt,
                  const art::InputTag& which,
                  art::View<recob::Track>& track);

    int GetShowers(const art::Event& evt,
                   const art::InputTag& which,
                   art::View<recob::Shower>& shower);

    int GetVertices(const art::Event& evt,
                    const art::InputTag& which,
                    art::PtrVector<recob::Vertex>& vertex);

    int GetSeeds(const art::Event& evt,
                 const art::InputTag& which,
                 art::PtrVector<recob::Seed>& seed);

    int GetOpFlashes(const art::Event& evt,
                     const art::InputTag& which,
                     art::PtrVector<recob::OpFlash>& opflash);

    int GetEvents(const art::Event& evt,
                  const art::InputTag& which,
                  art::PtrVector<recob::Event>& event);

    float SpacePointChiSq(const std::vector<art::Ptr<recob::Hit>>&) const;

    std::vector<std::array<double, 3>> Circle3D(const TVector3& pos,
                                                const TVector3& axisDir,
                                                const double& radius);

    int CountHits(const art::Event& evt,
                  const art::InputTag& which,
                  unsigned int cryostat,
                  unsigned int tpc,
                  unsigned int plane);

  private:
    using ISpacePointDrawerPtr = std::unique_ptr<evdb_tool::ISpacePoints3D>;

    ISpacePointDrawerPtr fAllSpacePointDrawer;
    ISpacePointDrawerPtr fSpacePointDrawer;

    std::vector<int> fWireMin; ///< lowest wire in interesting region for each plane
    std::vector<int> fWireMax; ///< highest wire in interesting region for each plane
    std::vector<int> fTimeMin; ///< lowest time in interesting region for each plane
    std::vector<int> fTimeMax; ///< highest time in interesting region for each plane

    std::vector<double> fRawCharge;       ///< Sum of Raw Charge
    std::vector<double> fConvertedCharge; ///< Sum of Charge Converted using Birks' formula
  };
}

#endif
////////////////////////////////////////////////////////////////////////
