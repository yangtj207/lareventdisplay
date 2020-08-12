////////////////////////////////////////////////////////////////////////
/// \file   OpFlash3DDrawer_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TPolyLine3D.h"

// Eigen
#include <Eigen/Core>

namespace evdb_tool {

  class OpFlash3DDrawer : public I3DDrawer {
  public:
    explicit OpFlash3DDrawer(const fhicl::ParameterSet&);

    void Draw(const art::Event&, evdb::View3D*) const override;

  private:
    void DrawRectangularBox(evdb::View3D*,
                            const Eigen::Vector3f&,
                            const Eigen::Vector3f&,
                            int,
                            int,
                            int) const;
  };

  //----------------------------------------------------------------------
  // Constructor.
  OpFlash3DDrawer::OpFlash3DDrawer(const fhicl::ParameterSet& pset) {}

  void
  OpFlash3DDrawer::Draw(const art::Event& event, evdb::View3D* view) const
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

    if (recoOpt->fDrawOpFlashes == 0) return;

    // Service recovery
    art::ServiceHandle<geo::Geometry> geo;
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clock_data);
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    std::vector<geo::PlaneID> planeIDVec;

    planeIDVec.push_back(geo::PlaneID(0, 0, 0));
    planeIDVec.push_back(geo::PlaneID(0, 1, 0));
    planeIDVec.push_back(geo::PlaneID(1, 0, 0));
    planeIDVec.push_back(geo::PlaneID(1, 1, 0));

    art::Handle<std::vector<recob::OpFlash>> opFlashHandle;

    // This seems like a time waster but we want to get the full color scale for all OpHits... so loops away...
    std::vector<float> opHitPEVec;

    // This is almost identically the same for loop we will re-excute below... sigh...
    for (size_t idx = 0; idx < recoOpt->fOpFlashLabels.size(); idx++) {
      art::InputTag opFlashProducer = recoOpt->fOpFlashLabels[idx];

      event.getByLabel(opFlashProducer, opFlashHandle);

      if (!opFlashHandle.isValid()) continue;
      if (opFlashHandle->size() == 0) continue;

      // To get associations we'll need an art ptr vector...
      art::PtrVector<recob::OpFlash> opFlashVec;

      for (size_t idx = 0; idx < opFlashHandle->size(); idx++)
        opFlashVec.push_back(art::Ptr<recob::OpFlash>(opFlashHandle, idx));

      // Recover the associations to op hits
      art::FindManyP<recob::OpHit> opHitAssnVec(opFlashVec, event, opFlashProducer);

      if (opHitAssnVec.size() == 0) continue;

      // Start the loop over flashes
      for (const auto& opFlashPtr : opFlashVec) {
        std::cout << "--> opFlash PE: " << opFlashPtr->TotalPE() << ", Time: " << opFlashPtr->Time()
                  << ", width: " << opFlashPtr->TimeWidth() << ", y/w: " << opFlashPtr->YCenter()
                  << "/" << opFlashPtr->YWidth() << ", Z/w: " << opFlashPtr->ZCenter() << "/"
                  << opFlashPtr->ZWidth() << std::endl;
        // Make some selections...
        if (opFlashPtr->TotalPE() < recoOpt->fFlashMinPE) continue;
        if (opFlashPtr->Time() < recoOpt->fFlashTMin) continue;
        if (opFlashPtr->Time() > recoOpt->fFlashTMax) continue;

        // Start by going through the associated OpHits
        const std::vector<art::Ptr<recob::OpHit>>& opHitVec = opHitAssnVec.at(opFlashPtr.key());

        for (const auto& opHit : opHitVec)
          opHitPEVec.push_back(opHit->PE());
      }
    }

    // Do we have any flashes and hits?
    if (!opHitPEVec.empty()) {
      // Sorting is good for mind and body...
      std::sort(opHitPEVec.begin(), opHitPEVec.end());

      float minTotalPE = opHitPEVec.front();
      float maxTotalPE = opHitPEVec[0.9 * opHitPEVec.size()];

      // Now we can set the scaling factor for PE
      float opHitPEScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) /
                         (maxTotalPE - minTotalPE));

      // We are meant to draw the flashes/hits, so loop over the list of input flashes
      for (size_t idx = 0; idx < recoOpt->fOpFlashLabels.size(); idx++) {
        art::InputTag opFlashProducer = recoOpt->fOpFlashLabels[idx];

        event.getByLabel(opFlashProducer, opFlashHandle);

        if (!opFlashHandle.isValid()) continue;
        if (opFlashHandle->size() == 0) continue;

        // To get associations we'll need an art ptr vector...
        art::PtrVector<recob::OpFlash> opFlashVec;

        for (size_t idx = 0; idx < opFlashHandle->size(); idx++)
          opFlashVec.push_back(art::Ptr<recob::OpFlash>(opFlashHandle, idx));

        // Recover the associations to op hits
        art::FindManyP<recob::OpHit> opHitAssnVec(opFlashVec, event, opFlashProducer);

        if (opHitAssnVec.size() == 0) continue;

        // Start the loop over flashes
        for (const auto& opFlashPtr : opFlashVec) {
          // Make some selections...
          if (opFlashPtr->TotalPE() < recoOpt->fFlashMinPE) continue;
          if (opFlashPtr->Time() < recoOpt->fFlashTMin) continue;
          if (opFlashPtr->Time() > recoOpt->fFlashTMax) continue;

          // Start by going through the associated OpHits
          const std::vector<art::Ptr<recob::OpHit>> opHitVec = opHitAssnVec.at(opFlashPtr.key());

          // We use the flash time to give us an x position (for now... will
          // need a better way eventually)
          float flashTick = opFlashPtr->Time() / sampling_rate(clock_data) * 1e3 +
                            det_prop.GetXTicksOffset(planeIDVec[idx]);
          float flashWidth = opFlashPtr->TimeWidth() / sampling_rate(clock_data) * 1e3 +
                             det_prop.GetXTicksOffset(planeIDVec[idx]);

          // Now convert from time to distance...
          float flashXpos = det_prop.ConvertTicksToX(flashTick, planeIDVec[idx]);
          float flashXWid = det_prop.ConvertTicksToX(flashWidth, planeIDVec[idx]);

          // Loop through the OpHits here
          for (const auto& opHit : opHitVec) {
            unsigned int opChannel = opHit->OpChannel();
            const geo::OpDetGeo& opHitGeo = geo->OpDetGeoFromOpChannel(opChannel);
            const geo::Point_t& opHitPos = opHitGeo.GetCenter();
            float zWidth = opHitGeo.HalfW();
            float yWidth = opHitGeo.HalfH();

            Eigen::Vector3f opHitLo(
              opHitPos.X() - flashXWid, opHitPos.Y() - yWidth, opHitPos.Z() - zWidth);
            Eigen::Vector3f opHitHi(
              opHitPos.X() + flashXWid, opHitPos.Y() + yWidth, opHitPos.Z() + zWidth);

            // Temporary kludge...
            flashXpos = opHitPos.X();

            float peFactor = cst->fRecoQLow[geo::kCollection] +
                             opHitPEScale * std::min(maxTotalPE, float(opHit->PE()));

            int chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(peFactor);

            DrawRectangularBox(view, opHitLo, opHitHi, chargeColorIdx, 2, 1);
          }

          std::cout << "     == flashtick: " << flashTick << ", flashwidth: " << flashWidth
                    << ", flashXpos: " << flashXpos << ", wid: " << flashXWid
                    << ", opHitPEScale: " << opHitPEScale << std::endl;

          //            std::vector<Eigen::Vector3f>
          Eigen::Vector3f coordsLo(flashXpos - flashXWid,
                                   opFlashPtr->YCenter() - opFlashPtr->YWidth(),
                                   opFlashPtr->ZCenter() - opFlashPtr->ZWidth());
          Eigen::Vector3f coordsHi(flashXpos + flashXWid,
                                   opFlashPtr->YCenter() + opFlashPtr->YWidth(),
                                   opFlashPtr->ZCenter() + opFlashPtr->ZWidth());

          DrawRectangularBox(view, coordsLo, coordsHi, kRed, 2, 1);
        }
      }
    }

    return;
  }

  void
  OpFlash3DDrawer::DrawRectangularBox(evdb::View3D* view,
                                      const Eigen::Vector3f& coordsLo,
                                      const Eigen::Vector3f& coordsHi,
                                      int color,
                                      int width,
                                      int style) const
  {
    TPolyLine3D& top = view->AddPolyLine3D(5, color, width, style);
    top.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(1, coordsHi[0], coordsHi[1], coordsLo[2]);
    top.SetPoint(2, coordsHi[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    top.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& side = view->AddPolyLine3D(5, color, width, style);
    side.SetPoint(0, coordsHi[0], coordsHi[1], coordsLo[2]);
    side.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    side.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    side.SetPoint(3, coordsHi[0], coordsHi[1], coordsHi[2]);
    side.SetPoint(4, coordsHi[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& side2 = view->AddPolyLine3D(5, color, width, style);
    side2.SetPoint(0, coordsLo[0], coordsHi[1], coordsLo[2]);
    side2.SetPoint(1, coordsLo[0], coordsLo[1], coordsLo[2]);
    side2.SetPoint(2, coordsLo[0], coordsLo[1], coordsHi[2]);
    side2.SetPoint(3, coordsLo[0], coordsHi[1], coordsHi[2]);
    side2.SetPoint(4, coordsLo[0], coordsHi[1], coordsLo[2]);

    TPolyLine3D& bottom = view->AddPolyLine3D(5, color, width, style);
    bottom.SetPoint(0, coordsLo[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(1, coordsHi[0], coordsLo[1], coordsLo[2]);
    bottom.SetPoint(2, coordsHi[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(3, coordsLo[0], coordsLo[1], coordsHi[2]);
    bottom.SetPoint(4, coordsLo[0], coordsLo[1], coordsLo[2]);

    return;
  }

  DEFINE_ART_CLASS_TOOL(OpFlash3DDrawer)
}
