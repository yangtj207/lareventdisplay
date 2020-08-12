/// \file    RecoBaseDrawer.cxx
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  brebel@fnal.gov

#include <cmath>
#include <limits>
#include <map>
#include <stdint.h>

#include "TBox.h"
#include "TH1.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TRotation.h"
#include "TText.h"
#include "TVector3.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Edge.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lareventdisplay/EventDisplay/3DDrawers/ISpacePoints3D.h"
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/eventdisplay.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "nuevdb/EventDisplayBase/View2D.h"
#include "nuevdb/EventDisplayBase/View3D.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void
  writeErrMsg(const char* fcn, cet::exception const& e)
  {
    mf::LogWarning("RecoBaseDrawer") << "RecoBaseDrawer::" << fcn << " failed with message:\n" << e;
  }
} // namespace

namespace evd {

  //......................................................................
  RecoBaseDrawer::RecoBaseDrawer()
  {
    art::ServiceHandle<geo::Geometry const> geo;
    art::ServiceHandle<evd::RawDrawingOptions const> rawOptions;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOptions;

    fWireMin.resize(0);
    fWireMax.resize(0);
    fTimeMin.resize(0);
    fTimeMax.resize(0);
    fRawCharge.resize(0);
    fConvertedCharge.resize(0);

    // set the list of channels in this detector
    for (size_t t = 0; t < geo->NTPC(); ++t) {
      unsigned int nplanes = geo->Nplanes(t);
      fWireMin.resize(nplanes, -1);
      fWireMax.resize(nplanes, -1);
      fTimeMin.resize(nplanes, -1);
      fTimeMax.resize(nplanes, -1);
      fRawCharge.resize(nplanes, 0);
      fConvertedCharge.resize(nplanes, 0);
      for (size_t p = 0; p < geo->Nplanes(t); ++p) {
        fWireMin[p] = 0;
        fWireMax[p] = geo->TPC(t).Plane(p).Nwires();
        fTimeMin[p] = 0;
        fTimeMax[p] = rawOptions->fTicks;
      } // end loop over planes
    }   // end loop over TPCs

    fAllSpacePointDrawer =
      art::make_tool<evdb_tool::ISpacePoints3D>(recoOptions->fAllSpacePointDrawerParams);
    fSpacePointDrawer =
      art::make_tool<evdb_tool::ISpacePoints3D>(recoOptions->fSpacePointDrawerParams);
  }

  //......................................................................
  RecoBaseDrawer::~RecoBaseDrawer() {}

  //......................................................................
  void
  RecoBaseDrawer::Wire2D(const art::Event& evt, evdb::View2D* view, unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    lariov::ChannelStatusProvider const& channelStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    int ticksPerPoint = rawOpt->fTicksPerPoint;

    // to make det independent later:
    double mint = 5000;
    double maxt = 0;
    double minw = 5000;
    double maxw = 0;

    geo::PlaneID pid(rawOpt->fCryostat, rawOpt->fTPC, plane);

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      if (wires.size() < 1) continue;

      for (size_t i = 0; i < wires.size(); ++i) {

        uint32_t channel = wires[i]->Channel();

        if (!rawOpt->fSeeBadChannels && channelStatus.IsBad(channel)) continue;

        std::vector<geo::WireID> wireids = geo->ChannelToWire(channel);

        geo::SigType_t sigType = geo->SignalType(channel);

        for (auto const& wid : wireids) {
          if (wid.planeID() != pid) continue;

          double wire = 1. * wid.Wire;
          double tick = 0;
          // get the unpacked ROIs
          std::vector<float> wirSig = wires[i]->Signal();
          if (wirSig.size() == 0) continue;
          // get an iterator over the adc values
          std::vector<float>::const_iterator itr = wirSig.begin();
          while (itr != wirSig.end()) {
            int ticksUsed = 0;
            double tdcsum = 0.;
            double adcsum = 0.;
            while (ticksUsed < ticksPerPoint && itr != wirSig.end()) {
              tdcsum += tick;
              adcsum += (1. * (*itr));
              ++ticksUsed;
              tick += 1.;
              itr++; // this advance of the iterator is sufficient for the external loop too
            }
            double adc = adcsum / ticksPerPoint;
            double tdc = tdcsum / ticksPerPoint;

            if (TMath::Abs(adc) < rawOpt->fMinSignal) continue;
            if (tdc > rawOpt->fTicks) continue;

            int co = 0;
            double sf = 1.;
            double q0 = 1000.0;

            co = cst->CalQ(sigType).GetColor(adc);
            if (rawOpt->fScaleDigitsByCharge) {
              sf = sqrt(adc / q0);
              if (sf > 1.0) sf = 1.0;
            }

            if (wire < minw) minw = wire;
            if (wire > maxw) maxw = wire;
            if (tdc < mint) mint = tdc;
            if (tdc > maxt) maxt = tdc;

            if (rawOpt->fAxisOrientation < 1) {
              TBox& b1 = view->AddBox(wire - sf * 0.5,
                                      tdc - sf * 0.5 * ticksPerPoint,
                                      wire + sf * 0.5,
                                      tdc + sf * 0.5 * ticksPerPoint);
              b1.SetFillStyle(1001);
              b1.SetFillColor(co);
              b1.SetBit(kCannotPick);
            }
            else {
              TBox& b1 = view->AddBox(tdc - sf * 0.5 * ticksPerPoint,
                                      wire - sf * 0.5,
                                      tdc + sf * 0.5 * ticksPerPoint,
                                      wire + sf * 0.5);
              b1.SetFillStyle(1001);
              b1.SetFillColor(co);
              b1.SetBit(kCannotPick);
            }
          } // end loop over samples
        }   //end loop over wire segments
      }     //end loop over wires
    }       // end loop over wire module labels

    fWireMin[plane] = minw;
    fWireMax[plane] = maxw;
    fTimeMin[plane] = mint;
    fTimeMax[plane] = maxt;

    // Add a loop to draw dead wires in 2D display
    double startTick(50.);
    double endTick((rawOpt->fTicks - 50.) * ticksPerPoint);

    for (size_t wireNo = 0; wireNo < geo->Nwires(pid); wireNo++) {
      raw::ChannelID_t channel =
        geo->PlaneWireToChannel(geo::WireID(rawOpt->fCryostat, rawOpt->fTPC, plane, wireNo));

      if (!rawOpt->fSeeBadChannels && channelStatus.IsBad(channel)) {
        double wire = 1. * wireNo;
        TLine& line = view->AddLine(wire, startTick, wire, endTick);
        line.SetLineColor(kGray);
        line.SetLineWidth(1.0);
        line.SetBit(kCannotPick);
      }
    }
  }

  //......................................................................
  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  int
  RecoBaseDrawer::Hit2D(const art::Event& evt,
                        detinfo::DetectorPropertiesData const& detProp,
                        evdb::View2D* view,
                        unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    int nHitsDrawn(0);

    if (recoOpt->fDrawHits == 0) return nHitsDrawn;
    if (rawOpt->fDrawRawDataOrCalibWires < 1) return nHitsDrawn;

    fRawCharge[plane] = 0;
    fConvertedCharge[plane] = 0;

    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fHitLabels[imod];

      std::vector<const recob::Hit*> hits;
      this->GetHits(evt, which, hits, plane);

      // Display all hits on the two 2D views provided
      for (auto itr : hits) {

        if (itr->WireID().TPC != rawOpt->fTPC || itr->WireID().Cryostat != rawOpt->fCryostat)
          continue;

        // Try to get the "best" charge measurement, ie. the one last in
        // the calibration chain
        fRawCharge[itr->WireID().Plane] += itr->PeakAmplitude();
        double dQdX = itr->PeakAmplitude() / geo->WirePitch() / detProp.ElectronsToADC();
        fConvertedCharge[itr->WireID().Plane] += detProp.BirksCorrection(dQdX);
      } // loop on hits

      nHitsDrawn = this->Hit2D(hits, kBlack, view, recoOpt->fDrawAllWireIDs);

    } // loop on imod folders

    return nHitsDrawn;
  }

  //......................................................................
  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param hits   : vector of hits for the veiw
  /// @param color  : color of associated cluster/prong
  /// @param view   : Pointer to view to draw on
  ///
  /// assumes the hits are all from the correct plane for the given view
  int
  RecoBaseDrawer::Hit2D(std::vector<const recob::Hit*> hits,
                        int color,
                        evdb::View2D* view,
                        bool allWireIDs,
                        bool drawConnectingLines,
                        int lineWidth)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    unsigned int w = 0;
    unsigned int wold = 0;
    float timeold = 0.;

    if (color == -1) color = recoOpt->fSelectedHitColor;

    int nHitsDrawn(0);

    for (const auto& hit : hits) {
      // Note that the WireID in the hit object is useless for those detectors where a channel can correspond to
      // more than one plane/wire. So our plan is to recover the list of wire IDs from the channel number and
      // loop over those (if there are any)
      // However, we need to preserve the option for drawing hits only associated to the wireID it contains
      std::vector<geo::WireID> wireIDs;

      if (allWireIDs)
        wireIDs = geo->ChannelToWire(hit->Channel());
      else
        wireIDs.push_back(hit->WireID());

      // Loop to find match
      for (const auto& wireID : wireIDs) {
        if (wireID.TPC != rawOpt->fTPC || wireID.Cryostat != rawOpt->fCryostat) continue;

        if (std::isnan(hit->PeakTime()) || std::isnan(hit->Integral())) {
          std::cout << "====>> Found hit with a NAN, channel: " << hit->Channel()
                    << ", start/end: " << hit->StartTick() << "/" << hit->EndTick()
                    << ", chisquare: " << hit->GoodnessOfFit() << std::endl;
        }

        if (hit->PeakTime() > rawOpt->fTicks) continue;

        w = wireID.Wire;

        // Try to get the "best" charge measurement, ie. the one last in
        // the calibration chain
        float time = hit->PeakTime();
        float rms = 0.5 * hit->RMS();

        if (rawOpt->fAxisOrientation < 1) {
          TBox& b1 = view->AddBox(w - 0.5, time - rms, w + 0.5, time + rms);
          if (drawConnectingLines && nHitsDrawn > 0) {
            TLine& l = view->AddLine(w, time, wold, timeold);
            l.SetLineColor(color);
            l.SetBit(kCannotPick);
          }
          b1.SetFillStyle(0);
          b1.SetBit(kCannotPick);
          b1.SetLineColor(color);
          b1.SetLineWidth(lineWidth);
        }
        else {
          TBox& b1 = view->AddBox(time - rms, w - 0.5, time + rms, w + 0.5);
          if (drawConnectingLines && nHitsDrawn > 0) {
            TLine& l = view->AddLine(time, w, timeold, wold);
            l.SetLineColor(color);
            l.SetBit(kCannotPick);
          }
          b1.SetFillStyle(0);
          b1.SetBit(kCannotPick);
          b1.SetLineColor(color);
          b1.SetLineWidth(lineWidth);
        }
        wold = w;
        timeold = time;
        nHitsDrawn++;
      }
    } // loop on hits

    return nHitsDrawn;
  }

  //........................................................................
  int
  RecoBaseDrawer::Hit2D(std::vector<const recob::Hit*> hits, evdb::View2D* view, float cosmicscore)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    unsigned int w(0);
    unsigned int wold(0);
    float timeold(0.);
    int nHitsDrawn(0);

    for (const auto& hit : hits) {
      // check that we are in the correct TPC
      // the view should tell use we are in the correct plane
      if (hit->WireID().TPC != rawOpt->fTPC || hit->WireID().Cryostat != rawOpt->fCryostat)
        continue;

      w = hit->WireID().Wire;

      // Try to get the "best" charge measurement, ie. the one last in
      // the calibration chain
      float time = hit->PeakTime();

      if (rawOpt->fAxisOrientation < 1) {
        if (nHitsDrawn > 0) {
          TLine& l = view->AddLine(w, time + 100, wold, timeold + 100);
          l.SetLineWidth(3);
          l.SetLineColor(1);
          if (cosmicscore > 0.5) l.SetLineColor(kMagenta);
          l.SetBit(kCannotPick);
        }
      }
      else {
        if (nHitsDrawn > 0) {
          TLine& l = view->AddLine(time + 20, w, timeold + 20, wold);
          l.SetLineColor(1);
          if (cosmicscore > 0.5) l.SetLineStyle(2);
          l.SetBit(kCannotPick);
        }
      }

      wold = w;
      timeold = time;
      nHitsDrawn++;
    } // loop on hits

    return nHitsDrawn;
  }

  //........................................................................
  int
  RecoBaseDrawer::GetRegionOfInterest(int plane, int& minw, int& maxw, int& mint, int& maxt)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if ((unsigned int)plane > fWireMin.size()) {
      mf::LogWarning("RecoBaseDrawer")
        << " Requested plane " << plane << " is larger than those available ";
      return -1;
    }

    minw = fWireMin[plane];
    maxw = fWireMax[plane];
    mint = fTimeMin[plane];
    maxt = fTimeMax[plane];

    //make values a bit larger, but make sure they don't go out of bounds
    minw = (minw - 30 < 0) ? 0 : minw - 30;
    mint = (mint - 10 < 0) ? 0 : mint - 10;

    int fTicks = rawOpt->fTicks;

    maxw = (maxw + 10 > (int)geo->Nwires(plane)) ? geo->Nwires(plane) : maxw + 10;
    maxt = (maxt + 10 > fTicks) ? fTicks : maxt + 10;

    return 0;
  }

  //......................................................................
  void
  RecoBaseDrawer::GetChargeSum(int plane, double& charge, double& convcharge)
  {
    charge = fRawCharge[plane];
    convcharge = fConvertedCharge[plane];

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::EndPoint2D(const art::Event& evt, evdb::View2D* view, unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDraw2DEndPoints == 0) return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    for (size_t imod = 0; imod < recoOpt->fEndPoint2DLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fEndPoint2DLabels[imod];

      art::PtrVector<recob::EndPoint2D> ep2d;
      this->GetEndPoint2D(evt, which, ep2d);

      for (size_t iep = 0; iep < ep2d.size(); ++iep) {
        // only worry about end points with the correct view
        if (ep2d[iep]->View() != gview) continue;

        ///\todo - have to verify that we are in the right TPC, but to do that we
        // need to be sure that all EndPoint2D objects have filled the required information

        // draw cluster with unique marker
        // Place this cluster's unique marker at the hit's location
        int color = evd::kColor[ep2d[iep]->ID() % evd::kNCOLS];

        double x = ep2d[iep]->WireID().Wire;
        double y = ep2d[iep]->DriftTime();

        if (rawOpt->fAxisOrientation > 0) {
          x = ep2d[iep]->DriftTime();
          y = ep2d[iep]->WireID().Wire;
        }

        TMarker& strt = view->AddMarker(x, y, color, 30, 2.0);
        strt.SetMarkerColor(color);
        // BB: draw the ID
        if (recoOpt->fDraw2DEndPoints > 1) {
          std::string s = "2V" + std::to_string(ep2d[iep]->ID());
          char const* txt = s.c_str();
          TText& vtxID = view->AddText(x, y + 20, txt);
          vtxID.SetTextColor(color);
          vtxID.SetTextSize(0.05);
        }

      } // loop on iep end points
    }   // loop on imod folders

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::OpFlash2D(const art::Event& evt,
                            detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            evdb::View2D* view,
                            unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawOpFlashes == 0) return;

    art::ServiceHandle<geo::Geometry const> geo;
    geo::PlaneID pid(rawOpt->fCryostat, rawOpt->fTPC, plane);

    for (size_t imod = 0; imod < recoOpt->fOpFlashLabels.size(); ++imod) {
      const art::InputTag which = recoOpt->fOpFlashLabels[imod];

      art::PtrVector<recob::OpFlash> opflashes;
      this->GetOpFlashes(evt, which, opflashes);

      if (opflashes.size() < 1) continue;

      int NFlashes = opflashes.size();
      //double TopCoord = 1000;

      MF_LOG_VERBATIM("RecoBaseDrawer") << "Total " << NFlashes << " flashes.";

      // project each seed into this view
      for (size_t iof = 0; iof < opflashes.size(); ++iof) {
        if (opflashes[iof]->TotalPE() < recoOpt->fFlashMinPE) continue;
        if (opflashes[iof]->Time() < recoOpt->fFlashTMin) continue;
        if (opflashes[iof]->Time() > recoOpt->fFlashTMax) continue;
        int Color = evd::kColor[(iof) % evd::kNCOLS];
        MF_LOG_VERBATIM("RecoBaseDrawer")
          << "Flash t: " << opflashes[iof]->Time() << "\t y,z : " << opflashes[iof]->YCenter()
          << ", " << opflashes[iof]->ZCenter() << " \t PE :" << opflashes[iof]->TotalPE();

        float flashtick =
          opflashes[iof]->Time() / sampling_rate(clockData) * 1e3 + detProp.GetXTicksOffset(pid);
        float wire0 = FLT_MAX;
        float wire1 = FLT_MIN;

        //Find the 4 corners and convert them to wire numbers
        std::vector<TVector3> points;
        points.push_back(TVector3(0,
                                  opflashes[iof]->YCenter() - opflashes[iof]->YWidth(),
                                  opflashes[iof]->ZCenter() - opflashes[iof]->ZWidth()));
        points.push_back(TVector3(0,
                                  opflashes[iof]->YCenter() - opflashes[iof]->YWidth(),
                                  opflashes[iof]->ZCenter() + opflashes[iof]->ZWidth()));
        points.push_back(TVector3(0,
                                  opflashes[iof]->YCenter() + opflashes[iof]->YWidth(),
                                  opflashes[iof]->ZCenter() - opflashes[iof]->ZWidth()));
        points.push_back(TVector3(0,
                                  opflashes[iof]->YCenter() + opflashes[iof]->YWidth(),
                                  opflashes[iof]->ZCenter() + opflashes[iof]->ZWidth()));

        for (size_t i = 0; i < points.size(); ++i) {
          geo::WireID wireID;
          try {
            wireID = geo->NearestWireID(points[i], pid);
          }
          catch (geo::InvalidWireError const& e) {
            wireID = e.suggestedWireID(); // pick the closest valid wire
          }
          if (wireID.Wire < wire0) wire0 = wireID.Wire;
          if (wireID.Wire > wire1) wire1 = wireID.Wire;
        }
        if (rawOpt->fAxisOrientation > 0) {
          TLine& line = view->AddLine(flashtick, wire0, flashtick, wire1);
          line.SetLineWidth(2);
          line.SetLineStyle(2);
          line.SetLineColor(Color);
        }
        else {
          TLine& line = view->AddLine(wire0, flashtick, wire1, flashtick);
          line.SetLineWidth(2);
          line.SetLineStyle(2);
          line.SetLineColor(Color);
        }
      } // loop on opflashes
    }   // loop on imod folders

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Seed2D(const art::Event& evt,
                         detinfo::DetectorPropertiesData const& detProp,
                         evdb::View2D* view,
                         unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawSeeds == 0) return;

    for (size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fSeedLabels[imod];

      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);

      if (seeds.size() < 1) continue;

      // project each seed into this view
      for (size_t isd = 0; isd < seeds.size(); ++isd) {
        double SeedPoint[3];
        double SeedDir[3];
        double SeedPointErr[3];
        double SeedDirErr[3];
        double SeedEnd1[3];
        double SeedEnd2[3];

        seeds[isd]->GetPoint(SeedPoint, SeedPointErr);
        seeds[isd]->GetDirection(SeedDir, SeedDirErr);

        SeedEnd1[0] = SeedPoint[0] + SeedDir[0];
        SeedEnd1[1] = SeedPoint[1] + SeedDir[1];
        SeedEnd1[2] = SeedPoint[2] + SeedDir[2];

        SeedEnd2[0] = SeedPoint[0] - SeedDir[0];
        SeedEnd2[1] = SeedPoint[1] - SeedDir[1];
        SeedEnd2[2] = SeedPoint[2] - SeedDir[2];

        // Draw seed on evd
        // int color  = kColor[seeds[isd]->ID()%kNCOLS];
        int color = evd::kColor[0];
        unsigned int wirepoint = 0;
        unsigned int wireend1 = 0;
        unsigned int wireend2 = 0;
        try {
          wirepoint = geo->NearestWire(SeedPoint, plane, rawOpt->fTPC, rawOpt->fCryostat);
        }
        catch (cet::exception& e) {
          wirepoint = atoi(e.explain_self().substr(e.explain_self().find("#") + 1, 5).c_str());
        }
        try {
          wireend1 = geo->NearestWire(SeedEnd1, plane, rawOpt->fTPC, rawOpt->fCryostat);
        }
        catch (cet::exception& e) {
          wireend1 = atoi(e.explain_self().substr(e.explain_self().find("#") + 1, 5).c_str());
        }
        try {
          wireend2 = geo->NearestWire(SeedEnd2, plane, rawOpt->fTPC, rawOpt->fCryostat);
        }
        catch (cet::exception& e) {
          wireend2 = atoi(e.explain_self().substr(e.explain_self().find("#") + 1, 5).c_str());
        }

        double x = wirepoint;
        double y = detProp.ConvertXToTicks(SeedPoint[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
        double x1 = wireend1;
        double y1 = detProp.ConvertXToTicks(SeedEnd1[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
        double x2 = wireend2;
        double y2 = detProp.ConvertXToTicks(SeedEnd2[0], plane, rawOpt->fTPC, rawOpt->fCryostat);

        if (rawOpt->fAxisOrientation > 0) {
          x = detProp.ConvertXToTicks(SeedPoint[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
          y = wirepoint;
          x1 = detProp.ConvertXToTicks(SeedEnd1[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
          y1 = wireend1;
          x2 = detProp.ConvertXToTicks(SeedEnd2[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
          y2 = wireend2;
        }

        TMarker& strt = view->AddMarker(x, y, color, 4, 1.5);
        TLine& line = view->AddLine(x1, y1, x2, y2);
        strt.SetMarkerColor(color);
        line.SetLineColor(color);
        line.SetLineWidth(2.0);
      } // loop on seeds
    }   // loop on imod folders

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Slice2D(const art::Event& evt,
                          detinfo::DetectorPropertiesData const& detProp,
                          evdb::View2D* view,
                          unsigned int plane)
  {
    // Color code hits associated with Slices
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawSlices == 0) return;

    art::ServiceHandle<geo::Geometry const> geo;

    static bool first = true;
    if (first) {
      std::cout
        << "******** DrawSlices: 0 = none, 1 = color coded, 2 = color coded + ID at slice center\n";
      std::cout << "  3 = open circle at slice center with size proportional to the AspectRatio. "
                   "Closed circles";
      std::cout << "      at the slice ends with connecting dotted lines\n";
      first = false;
    }
    unsigned int c = rawOpt->fCryostat;
    unsigned int t = rawOpt->fTPC;
    geo::PlaneID planeID(c, t, plane);

    for (size_t imod = 0; imod < recoOpt->fSliceLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fSliceLabels[imod];
      art::PtrVector<recob::Slice> slices;
      this->GetSlices(evt, which, slices);
      if (slices.size() < 1) continue;
      art::FindMany<recob::Hit> fmh(slices, evt, which);
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        int slcID(std::abs(slices[isl]->ID()));
        int color(evd::kColor[slcID % evd::kNCOLS]);
        if (recoOpt->fDrawSlices < 3) {
          // draw color-coded hits
          std::vector<const recob::Hit*> hits = fmh.at(isl);
          std::vector<const recob::Hit*> hits_on_plane;
          for (auto hit : hits) {
            if (hit->WireID().Plane == plane) { hits_on_plane.push_back(hit); }
          }
          if (this->Hit2D(hits_on_plane, color, view, false, false) < 1) continue;
          if (recoOpt->fDrawSlices == 2) {
            geo::Point_t slicePos(
              slices[isl]->Center().X(), slices[isl]->Center().Y(), slices[isl]->Center().Z());
            double tick = detProp.ConvertXToTicks(slices[isl]->Center().X(), planeID);
            double wire = geo->WireCoordinate(slicePos, planeID);
            std::string s = std::to_string(slcID);
            char const* txt = s.c_str();
            TText& slcID = view->AddText(wire, tick, txt);
            slcID.SetTextSize(0.05);
            slcID.SetTextColor(color);
          } // draw ID
        }
        else {
          // draw the center, end points and direction vector
          geo::Point_t slicePos(
            slices[isl]->Center().X(), slices[isl]->Center().Y(), slices[isl]->Center().Z());
          double tick = detProp.ConvertXToTicks(slices[isl]->Center().X(), planeID);
          double wire = geo->WireCoordinate(slicePos, planeID);
          float markerSize = 1;
          if (slices[isl]->AspectRatio() > 0) {
            markerSize = 1 / slices[isl]->AspectRatio();
            if (markerSize > 3) markerSize = 3;
          }
          TMarker& ctr = view->AddMarker(wire, tick, color, 24, markerSize);
          ctr.SetMarkerColor(color);
          // npts, color, width, style
          TPolyLine& pline = view->AddPolyLine(2, color, 2, 3);
          geo::Point_t slicePos0(
            slices[isl]->End0Pos().X(), slices[isl]->End0Pos().Y(), slices[isl]->End0Pos().Z());
          tick = detProp.ConvertXToTicks(slices[isl]->End0Pos().X(), planeID);
          wire = geo->WireCoordinate(slicePos0, planeID);
          TMarker& end0 = view->AddMarker(wire, tick, color, 20, 1.0);
          end0.SetMarkerColor(color);
          pline.SetPoint(0, wire, tick);
          geo::Point_t slicePos1(
            slices[isl]->End1Pos().X(), slices[isl]->End1Pos().Y(), slices[isl]->End1Pos().Z());
          tick = detProp.ConvertXToTicks(slices[isl]->End1Pos().X(), plane, t, c);
          wire = geo->WireCoordinate(slicePos1, planeID);
          TMarker& end1 = view->AddMarker(wire, tick, color, 20, 1.0);
          end1.SetMarkerColor(color);
          pline.SetPoint(1, wire, tick);
        }
      } // isl

    } // imod

  } // Slice2D
  //......................................................................
  void
  RecoBaseDrawer::Cluster2D(const art::Event& evt,
                            detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            evdb::View2D* view,
                            unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawClusters == 0) return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    // if user sets "DrawClusters" to 2, draw the clusters differently:
    //    bool drawAsMarkers = (recoOpt->fDrawClusters == 1 ||
    //                          recoOpt->fDrawClusters == 3);
    bool drawAsMarkers = recoOpt->fDrawClusters != 2;

    // draw connecting lines between cluster hits?
    bool drawConnectingLines = (recoOpt->fDrawClusters >= 3);

    static bool first = true;
    if (first) {
      std::cout << "******** DrawClusters: 0 = none, 1 = cluster hits, 2 = unique marker, 3 = "
                   "cluster hits with connecting lines.\n";
      std::cout << " 4 = with T<cluster or trajectory ID> P<PFParticle ID> color-matched. "
                   "Unmatched cluster IDs shown in black.\n";
      std::cout << " Color scheme: By cluster ID in each plane or by PFParticle ID (Self) if a "
                   "PFParticle - Cluster association exists.\n";
      first = false;
    }

    for (size_t imod = 0; imod < recoOpt->fClusterLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fClusterLabels[imod];

      art::PtrVector<recob::Cluster> clust;
      this->GetClusters(evt, which, clust);

      if (clust.size() < 1) continue;

      // We want to draw the hits that are associated to "free" space points (non clustered)
      // This is done here, before drawing the hits on clusters so they will be "under" the cluster
      // hits (since spacepoints could be made from a used 2D hit but then not used themselves)
      // Get the space points created by the PFParticle producer
      std::vector<art::Ptr<recob::SpacePoint>> spacePointVec;
      this->GetSpacePoints(evt, which, spacePointVec);

      // No space points no continue
      if (spacePointVec.size() > 0) {
        // Add the relations to recover associations cluster hits
        art::FindManyP<recob::Hit> spHitAssnVec(spacePointVec, evt, which);

        if (spHitAssnVec.isValid()) {
          // Create a local hit vector...
          std::vector<const recob::Hit*> freeHitVec;

          // loop through space points looking for those that are free
          for (const auto& spacePointPtr : spacePointVec) {
            if (spacePointPtr->Chisq() < -99.) {
              // Recover associated hits
              const std::vector<art::Ptr<recob::Hit>>& hitVec =
                spHitAssnVec.at(spacePointPtr.key());

              for (const auto& hitPtr : hitVec) {
                if (hitPtr.get()->WireID().Plane != plane) continue;

                freeHitVec.push_back(hitPtr.get());
              }
            }
          }

          // Draw the free hits in gray
          this->Hit2D(freeHitVec, kGray, view, false, false, false);
        }
      }

      // Ok, now proceed with our normal processing of hits on clusters
      art::FindMany<recob::Hit> fmh(clust, evt, which);
      art::FindManyP<recob::PFParticle> fmc(clust, evt, which);

      for (size_t ic = 0; ic < clust.size(); ++ic) {
        // only worry about clusters with the correct view
        //            if(clust[ic]->View() != gview) continue;
        if (clust[ic]->Plane().Plane != plane) continue;

        // see if we can set the color index in a sensible fashion
        int clusterIdx(std::abs(clust[ic]->ID()));
        int colorIdx(clusterIdx % evd::kNCOLS);
        bool pfpAssociation = false;
        int pfpIndex = INT_MAX;
        float cosmicscore = FLT_MIN;

        if (fmc.isValid()) {
          std::vector<art::Ptr<recob::PFParticle>> pfplist = fmc.at(ic);
          // Use the first one
          if (!pfplist.empty()) {
            clusterIdx = pfplist[0]->Self();
            colorIdx = clusterIdx % evd::kNCOLS;
            pfpAssociation = true;
            pfpIndex = pfplist[0]->Self();
            //Get cosmic score
            if (recoOpt->fDrawCosmicTags) {
              art::FindManyP<anab::CosmicTag> fmct(pfplist, evt, which);
              if (fmct.isValid()) {
                std::vector<art::Ptr<anab::CosmicTag>> ctlist = fmct.at(0);
                if (!ctlist.empty()) {
                  //std::cout<<"cosmic tag "<<ctlist[0]->CosmicScore()<<std::endl;
                  cosmicscore = ctlist[0]->CosmicScore();
                }
              }
            }
          } // pfplist is not empty
        }

        std::vector<const recob::Hit*> hits = fmh.at(ic);

        if (drawAsMarkers) {
          // draw cluster with unique marker
          // Place this cluster's unique marker at the hit's location
          int color = evd::kColor[colorIdx];

          // If there are no hits in this cryostat/TPC then we skip the rest
          // That no hits were drawn is the sign for this
          if (this->Hit2D(hits, color, view, false, drawConnectingLines) < 1) continue;

          if (recoOpt->fDrawCosmicTags && cosmicscore != FLT_MIN)
            this->Hit2D(hits, view, cosmicscore);

          if (recoOpt->fDrawClusters > 3) {
            // BB: draw the cluster ID
            //std::string s = std::to_string(clusterIdx);
            // TY: change to draw cluster id instead of index
            //                std::string s = std::to_string(clusterIdx) + "," + std::to_string(clust[ic]->ID());
            // BB: Put a T in front to denote a trajectory ID
            std::string s = "T" + std::to_string(clust[ic]->ID());
            // append the PFP index + 1 (sort of the ID)
            if (pfpIndex != INT_MAX) s = s + " P" + std::to_string(pfpIndex + 1);
            char const* txt = s.c_str();
            double wire = 0.5 * (clust[ic]->StartWire() + clust[ic]->EndWire());
            double tick = 20 + 0.5 * (clust[ic]->StartTick() + clust[ic]->EndTick());
            TText& clID = view->AddText(wire, tick, txt);
            clID.SetTextSize(0.05);
            if (pfpAssociation) { clID.SetTextColor(color); }
            else {
              clID.SetTextColor(kBlack);
            }
          } // recoOpt->fDrawClusters > 3
        }
        else {

          // default "outline" method:
          std::vector<double> tpts, wpts;

          this->GetClusterOutlines(hits, tpts, wpts, plane);

          int lcolor = 9; // line color
          int fcolor = 9; // fill color
          int width = 2;  // line width
          int style = 1;  // 1=solid line style
          if (view != 0) {
            TPolyLine& p1 = view->AddPolyLine(wpts.size(), lcolor, width, style);
            TPolyLine& p2 = view->AddPolyLine(wpts.size(), lcolor, width, style);
            p1.SetOption("f");
            p1.SetFillStyle(3003);
            p1.SetFillColor(fcolor);
            for (size_t i = 0; i < wpts.size(); ++i) {
              if (rawOpt->fAxisOrientation < 1) {
                p1.SetPoint(i, wpts[i], tpts[i]);
                p2.SetPoint(i, wpts[i], tpts[i]);
              }
              else {
                p1.SetPoint(i, tpts[i], wpts[i]);
                p2.SetPoint(i, tpts[i], wpts[i]);
              }
            } // loop on i points in ZX view
          }   // if we have a cluster in the ZX view
        }     // end if outline mode

        // draw the direction cosine of the cluster as well as it's starting point
        // (average of the start and end angle -- by default they are the same value)
        // thetawire is the angle measured CW from +z axis to wire
        //double thetawire = geo->TPC(t).Plane(plane).Wire(0).ThetaZ();
        double wirePitch = geo->WirePitch(gview);
        double driftvelocity = detProp.DriftVelocity();    // cm/us
        double timetick = sampling_rate(clockData) * 1e-3; // time sample in us
        // rotate coord system CCW around x-axis by pi-thetawire
        //   new yprime direction is perpendicular to the wire direction
        //   in the same plane as the wires and in the direction of
        //   increasing wire number
        //use yprime-component of dir cos in rotated coord sys to get
        //   dTdW (number of time ticks per unit of wire pitch)
        //double rotang = 3.1416-thetawire;
        this->Draw2DSlopeEndPoints(
          clust[ic]->StartWire(),
          clust[ic]->StartTick(),
          clust[ic]->EndWire(),
          clust[ic]->EndTick(),
          std::tan((clust[ic]->StartAngle() + clust[ic]->EndAngle()) / 2.) * wirePitch /
            driftvelocity / timetick,
          evd::kColor[colorIdx],
          view);

      } // loop on ic clusters
    }   // loop on imod folders

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Draw2DSlopeEndPoints(double xStart,
                                       double yStart,
                                       double xEnd,
                                       double yEnd,
                                       double slope,
                                       int color,
                                       evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDraw2DSlopeEndPoints < 1) return;

    double x1 = xStart;
    double y1 = yStart;
    double x2 = xEnd;
    double slope1 = slope;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (rawOpt->fAxisOrientation > 0) {
      x1 = yStart;
      y1 = xStart;
      x2 = yEnd;
      if (std::abs(slope) > 0.)
        slope1 = 1. / slope;
      else
        slope1 = 1.e6;
    }

    double deltaX = 0.5 * (x2 - x1);
    double xm = x1 + deltaX;
    double ym = y1 + deltaX * slope;

    TMarker& strt = view->AddMarker(xm, ym, color, kFullCircle, 1.0);
    strt.SetMarkerColor(color); // stupid line to shut up compiler warning

    //    double stublen = 50.0 ;
    double stublen = 2. * deltaX;
    TLine& l = view->AddLine(x1, y1, x1 + stublen, y1 + slope1 * stublen);
    l.SetLineColor(color);
    l.SetLineWidth(1); //2);

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Draw2DSlopeEndPoints(double x,
                                       double y,
                                       double slope,
                                       int color,
                                       evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDraw2DSlopeEndPoints < 1) return;

    double x1 = x;
    double y1 = y;
    double slope1 = slope;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (rawOpt->fAxisOrientation > 0) {
      x1 = y;
      y1 = x;
      if (std::abs(slope) > 0.)
        slope1 = 1. / slope;
      else
        slope1 = 1.e6;
    }

    TMarker& strt = view->AddMarker(x1, y1, color, kFullStar, 2.0);
    strt.SetMarkerColor(color); // stupid line to shut up compiler warning

    //    double stublen = 50.0 ;
    double stublen = 300.0;
    TLine& l = view->AddLine(x1, y1, x1 + stublen, y1 + slope1 * stublen);
    l.SetLineColor(color);
    l.SetLineWidth(2);
    l.SetLineStyle(2);

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Draw2DSlopeEndPoints(double x,
                                       double y,
                                       double cosx,
                                       double cosy,
                                       int color,
                                       evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDraw2DSlopeEndPoints < 1) return;

    double x1 = x;
    double y1 = y;
    double cosx1 = cosx;
    double cosy1 = cosy;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (rawOpt->fAxisOrientation > 0) {
      x1 = y;
      y1 = x;
      cosx1 = cosy;
      cosy1 = cosx;
    }

    TMarker& strt = view->AddMarker(x1, y1, color, kFullStar, 2.0);
    strt.SetMarkerColor(color); // stupid line to shut up compiler warning

    //    double stublen = 50.0 ;
    double stublen = 300.0;
    TLine& l = view->AddLine(x1, y1, x1 + stublen * cosx1, y1 + stublen * cosy1);
    l.SetLineColor(color);
    l.SetLineWidth(2);
    l.SetLineStyle(2);

    return;
  }

  //......................................................................
  ///
  /// Make a set of points which outline a cluster
  ///
  /// @param c      : Reco base cluster to outline
  /// @param wpts   : wire values of the outlines
  /// @param tpts   : tdc values of the outlines
  /// @param plane  : plane number
  ///
  void
  RecoBaseDrawer::GetClusterOutlines(std::vector<const recob::Hit*>& hits,
                                     std::vector<double>& wpts,
                                     std::vector<double>& tpts,
                                     unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;

    // Map wire numbers to highest and lowest in the plane
    std::map<unsigned int, double> wlo, whi;
    // On first pass, initialize
    for (size_t j = 0; j < hits.size(); ++j) {
      // check that we are on the correct plane and TPC
      if (hits[j]->WireID().Plane != plane || hits[j]->WireID().TPC != rawOpt->fTPC ||
          hits[j]->WireID().Cryostat != rawOpt->fCryostat)
        continue;

      wlo[hits[j]->WireID().Wire] = hits[j]->PeakTime();
      whi[hits[j]->WireID().Wire] = hits[j]->PeakTime();
    }

    double t = 0.;

    // Finalize on second pass
    for (size_t j = 0; j < hits.size(); ++j) {
      t = hits[j]->PeakTime();

      if (t < wlo[hits[j]->WireID().Wire]) wlo[hits[j]->WireID().Wire] = t;
      if (t > whi[hits[j]->WireID().Wire]) whi[hits[j]->WireID().Wire] = t;
    }

    // Loop over wires and low times to make lines along bottom
    // edge. Work from upstream edge to downstream edge
    std::map<unsigned int, double>::iterator itr(wlo.begin());
    std::map<unsigned int, double>::iterator itrEnd(wlo.end());
    for (; itr != itrEnd; ++itr) {
      unsigned int w = itr->first;
      t = itr->second;

      wpts.push_back(1. * w - 0.1);
      tpts.push_back(t - 0.1);
      wpts.push_back(1. * w + 0.1);
      tpts.push_back(t - 0.1);
    }

    // Loop over planes and high cells to make lines along top
    // edge. Work from downstream edge toward upstream edge
    std::map<unsigned int, double>::reverse_iterator ritr(whi.rbegin());
    std::map<unsigned int, double>::reverse_iterator ritrEnd(whi.rend());
    for (; ritr != ritrEnd; ++ritr) {
      unsigned int w = ritr->first;
      t = ritr->second;

      wpts.push_back(1. * w + 0.1);
      tpts.push_back(t + 0.1);
      wpts.push_back(1. * w - 0.1);
      tpts.push_back(t + 0.1);
    }

    // Add link to starting point to close the box
    wpts.push_back(wpts[0]);
    tpts.push_back(tpts[0]);

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawProng2D(detinfo::DetectorPropertiesData const& detProp,
                              std::vector<const recob::Hit*>& hits,
                              evdb::View2D* view,
                              unsigned int plane,
                              TVector3 const& startPos,
                              TVector3 const& startDir,
                              int id,
                              float cscore)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    unsigned int c = rawOpt->fCryostat;
    unsigned int t = rawOpt->fTPC;
    geo::PlaneID planeID(c, t, plane);
    geo::Point_t localPos(startPos.X(), startPos.Y(), startPos.Z());

    int color(evd::kColor2[id % evd::kNCOLS]);
    int lineWidth(1);

    if (cscore > 0.1 && recoOpt->fDrawCosmicTags) {
      color = kRed;
      if (cscore < 0.6) color = kMagenta;
      lineWidth = 3;
    }
    else if (cscore < -10000) { //shower hits
      lineWidth = 3;
    }

    // first draw the hits
    if (cscore < -1000) { //shower
      this->Hit2D(hits, color, view, false, false, lineWidth);
      if (recoOpt->fDrawShowers >= 1) {
        //draw the shower ID at the beginning of shower
        std::string s = std::to_string(id);
        char const* txt = s.c_str();
        double tick = 30 + detProp.ConvertXToTicks(startPos.X(), planeID);
        double wire = geo->WireCoordinate(localPos, planeID);
        TText& shwID = view->AddText(wire, tick, txt);
        shwID.SetTextColor(evd::kColor2[id % evd::kNCOLS]);
        shwID.SetTextSize(0.1);
      }
    }
    else
      this->Hit2D(hits, color, view, false, false, lineWidth);

    double tick0 = detProp.ConvertXToTicks(startPos.X(), planeID);
    double wire0 = geo->WireCoordinate(localPos, planeID);

    localPos = geo::Point_t(startPos + startDir); // Huh? what is this?

    double tick1 = detProp.ConvertXToTicks((startPos + startDir).X(), planeID);
    double wire1 = geo->WireCoordinate(localPos, planeID);
    double cost = 0;
    double cosw = 0;
    double ds = sqrt(pow(tick0 - tick1, 2) + pow(wire0 - wire1, 2));

    if (ds) {
      cost = (tick1 - tick0) / ds;
      cosw = (wire1 - wire0) / ds;
    }

    this->Draw2DSlopeEndPoints(wire0, tick0, cosw, cost, evd::kColor[id % evd::kNCOLS], view);

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawTrack2D(detinfo::DetectorClocksData const& clockData,
                              detinfo::DetectorPropertiesData const& detProp,
                              std::vector<const recob::Hit*>& hits,
                              evdb::View2D* view,
                              unsigned int plane,
                              const recob::Track* track,
                              int color,
                              int lineWidth)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;
    unsigned int c = rawOpt->fCryostat;
    unsigned int t = rawOpt->fTPC;

    // first draw the hits
    this->Hit2D(hits, color, view, false, true, lineWidth);

    const auto& startPos = track->Vertex();
    const auto& startDir = track->VertexDirection();

    // prepare to draw prongs
    double local[3] = {0.};
    double world[3] = {0.};
    geo->Cryostat(c).TPC(t).Plane(plane).LocalToWorld(local, world);
    world[1] = startPos.Y();
    world[2] = startPos.Z();

    // convert the starting position and direction from 3D to 2D coordinates
    double tick = detProp.ConvertXToTicks(startPos.X(), plane, t, c);
    double wire = 0.;
    try {
      wire = 1. * geo->NearestWire(world, plane, t, c);
    }
    catch (cet::exception& e) {
      wire = 1. * atoi(e.explain_self().substr(e.explain_self().find("#") + 1, 5).c_str());
    }

    // thetawire is the angle measured CW from +z axis to wire
    double thetawire = geo->TPC(t).Plane(plane).Wire(0).ThetaZ();
    double wirePitch = geo->WirePitch(hits[0]->View());
    double driftvelocity = detProp.DriftVelocity();    // cm/us
    double timetick = sampling_rate(clockData) * 1e-3; // time sample in us
    // rotate coord system CCW around x-axis by pi-thetawire
    //   new yprime direction is perpendicular to the wire direction
    //   in the same plane as the wires and in the direction of
    //   increasing wire number
    //use yprime-component of dir cos in rotated coord sys to get
    //   dTdW (number of time ticks per unit of wire pitch)
    double rotang = 3.1416 - thetawire;
    double yprime = std::cos(rotang) * startDir.Y() + std::sin(rotang) * startDir.Z();
    double dTdW = startDir.X() * wirePitch / driftvelocity / timetick / yprime;

    this->Draw2DSlopeEndPoints(wire, tick, dTdW, color, view);

    // Draw a line to the hit positions, starting from the vertex
    size_t nTrackHits = track->NumberTrajectoryPoints();
    //TPolyLine& pl         = view->AddPolyLine(track->CountValidPoints(),1,1,0); //kColor[id%evd::kNCOLS],1,0);
    TPolyLine& pl = view->AddPolyLine(0, 1, 1, 0); //kColor[id%evd::kNCOLS],1,0);

    size_t vidx = 0;
    for (size_t idx = 0; idx < nTrackHits; idx++) {
      if (track->HasValidPoint(idx) == 0) continue;
      const auto& hitPos = track->LocationAtPoint(idx);

      // Use "world" from above
      world[1] = hitPos.Y();
      world[2] = hitPos.Z();

      // convert the starting position and direction from 3D to 2D coordinates
      double tickHit = detProp.ConvertXToTicks(hitPos.X(), plane, t, c);
      double wireHit = 0.;
      try {
        wireHit = 1. * geo->NearestWire(world, plane, t, c);
      }
      catch (cet::exception& e) {
        wireHit = 1. * atoi(e.explain_self().substr(e.explain_self().find("#") + 1, 5).c_str());
      }
      const size_t tpc = geo->FindTPCAtPosition(hitPos).TPC;
      const size_t cryo = geo->FindCryostatAtPosition(hitPos);
      if (tpc == t && cryo == c) { pl.SetPoint(vidx++, wireHit, tickHit); }
    }
    //pl.SetPolyLine(vidx);

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Prong2D(const art::Event& evt,
                          detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          evdb::View2D* view,
                          unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    unsigned int cstat = rawOpt->fCryostat;
    unsigned int tpc = rawOpt->fTPC;
    geo::PlaneID planeID(cstat, tpc, plane);
    int tid = 0;

    if (recoOpt->fDrawTracks != 0) {
      for (size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
        art::InputTag const which = recoOpt->fTrackLabels[imod];

        art::View<recob::Track> track;
        this->GetTracks(evt, which, track);

        if (track.vals().size() < 1) continue;

        art::FindMany<recob::Hit> fmh(track, evt, which);

        art::InputTag const whichTag(
          recoOpt->fCosmicTagLabels.size() > imod ? recoOpt->fCosmicTagLabels[imod] : "");
        art::FindManyP<anab::CosmicTag> cosmicTrackTags(track, evt, whichTag);

        auto tracksProxy = proxy::getCollection<proxy::Tracks>(evt, which);

        // loop over the prongs and get the clusters and hits associated with
        // them.  only keep those that are in this view
        for (size_t t = 0; t < track.vals().size(); ++t) {
          // Check for possible issue
          if (track.vals().at(t)->NumberTrajectoryPoints() == 0) {
            std::cout << "***** Track with no trajectory points ********" << std::endl;
            continue;
          }

          if (recoOpt->fDrawTracks > 1) {
            // BB: draw the track ID at the end of the track
            geo::Point_t trackPos(track.vals().at(t)->End().X(),
                                  track.vals().at(t)->End().Y(),
                                  track.vals().at(t)->End().Z());
            double tick = 30 + detProp.ConvertXToTicks(trackPos.X(), plane, tpc, cstat);
            double wire = geo->WireCoordinate(trackPos, geo::PlaneID(cstat, tpc, plane));
            tid =
              track.vals().at(t)->ID() &
              65535; //this is a hack for PMA track id which uses the 16th bit to identify shower-like track.;
            std::string s = std::to_string(tid);
            char const* txt = s.c_str();
            TText& trkID = view->AddText(wire, tick, txt);
            trkID.SetTextColor(evd::kColor[tid % evd::kNCOLS]);
            trkID.SetTextSize(0.1);
          }

          float Score = -999;
          if (cosmicTrackTags.isValid()) {
            if (cosmicTrackTags.at(t).size() > 0) {
              art::Ptr<anab::CosmicTag> currentTag = cosmicTrackTags.at(t).at(0);
              Score = currentTag->CosmicScore();
            }
          }

          std::vector<const recob::Hit*> hits;
          if (track.vals().at(t)->NumberTrajectoryPoints() == fmh.at(t).size()) {
            auto tp = tracksProxy[t];
            for (auto point : tp.points()) {
              if (!point.isPointValid()) continue;
              hits.push_back(point.hit());
            }
          }
          else {
            hits = fmh.at(t);
          }
          // only get the hits for the current view
          std::vector<const recob::Hit*>::iterator itr = hits.begin();
          while (itr < hits.end()) {
            if ((*itr)->View() != gview)
              hits.erase(itr);
            else
              itr++;
          }

          const recob::Track* aTrack(track.vals().at(t));
          int color(evd::kColor[(aTrack->ID() & 65535) % evd::kNCOLS]);
          int lineWidth(1);

          if (Score > 0.1 && recoOpt->fDrawCosmicTags) {
            color = kRed;
            if (Score < 0.6) color = kMagenta;
            lineWidth = 3;
          }
          else if (Score < -10000) { //shower hits
            lineWidth = 3;
          }

          this->DrawTrack2D(clockData, detProp, hits, view, plane, aTrack, color, lineWidth);
        } // end loop over prongs
      }   // end loop over labels
    }     // end draw tracks

    if (recoOpt->fDrawShowers != 0) {
      static bool first = true;

      if (first) {
        std::cout << "DrawShower options: \n";
        std::cout << " 1 = Hits in shower color-coded by the shower ID\n";
        std::cout << " 2 = Same as 1 + shower axis and circle representing the shower cone\n";
        std::cout << "     Black cone = shower start dE/dx < 1 MeV/cm (< 1/2 MIP)\n";
        std::cout << "     Blue cone = shower start dE/dx < 3 MeV/cm (~1 MIP)\n";
        std::cout << "     Green cone = shower start 3 MeV/cm < dE/dx < 5 MeV/cm (~2 MIP)\n";
        std::cout << "     Red cone = shower start 5 MeV/cm < dE/dx (>2 MIP)\n";
        first = false;
      }
      for (size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
        art::InputTag const which = recoOpt->fShowerLabels[imod];

        art::View<recob::Shower> shower;
        this->GetShowers(evt, which, shower);
        if (shower.vals().size() < 1) continue;

        art::FindMany<recob::Hit> fmh(shower, evt, which);

        // loop over the prongs and get the clusters and hits associated with
        // them.  only keep those that are in this view
        for (size_t s = 0; s < shower.vals().size(); ++s) {

          std::vector<const recob::Hit*> hits = fmh.at(s);
          // only get the hits for the current view
          std::vector<const recob::Hit*>::iterator itr = hits.begin();
          while (itr < hits.end()) {
            if ((*itr)->View() != gview)
              hits.erase(itr);
            else
              itr++;
          }
          if (recoOpt->fDrawShowers > 1) {
            // BB draw a line between the start and end points and a "circle" that represents
            // the shower cone angle at the end point
            if (!shower.vals().at(s)->has_length()) continue;
            if (!shower.vals().at(s)->has_open_angle()) continue;

            TVector3 startPos = shower.vals().at(s)->ShowerStart();
            TVector3 dir = shower.vals().at(s)->Direction();
            double length = shower.vals().at(s)->Length();
            double openAngle = shower.vals().at(s)->OpenAngle();

            // Find the center of the cone base
            TVector3 endPos = startPos + length * dir;

            geo::Point_t localStart(startPos);
            geo::Point_t localEnd(endPos);

            double swire = geo->WireCoordinate(localStart, planeID);
            double stick = detProp.ConvertXToTicks(startPos.X(), planeID);
            double ewire = geo->WireCoordinate(localEnd, planeID);
            double etick = detProp.ConvertXToTicks(endPos.X(), planeID);
            TLine& coneLine = view->AddLine(swire, stick, ewire, etick);
            // color coding by dE/dx
            std::vector<double> dedxVec = shower.vals().at(s)->dEdx();
            //                      float dEdx = shower.vals().at(s)->dEdx()[plane];
            // use black for too-low dE/dx
            int color = kBlack;
            if (plane < dedxVec.size()) {
              if (dedxVec[plane] > 1 && dedxVec[plane] < 3) {
                // use blue for ~1 MIP
                color = kBlue;
              }
              else if (dedxVec[plane] < 5) {
                // use green for ~2 MIP
                color = kGreen;
              }
              else {
                // use red for >~ 2 MIP
                color = kRed;
              }
            }
            coneLine.SetLineColor(color);

            // Now find the 3D circle that represents the base of the cone
            double radius = length * openAngle;
            auto coneRim = Circle3D(endPos, dir, radius);
            TPolyLine& pline = view->AddPolyLine(coneRim.size(), color, 2, 0);
            // project these points into the plane
            for (unsigned short ipt = 0; ipt < coneRim.size(); ++ipt) {
              geo::Point_t localPos(coneRim[ipt][0], coneRim[ipt][1], coneRim[ipt][2]);

              double wire = geo->WireCoordinate(localPos, planeID);
              double tick = detProp.ConvertXToTicks(coneRim[ipt][0], planeID);
              pline.SetPoint(ipt, wire, tick);
            } // ipt
          }
          this->DrawProng2D(detProp,
                            hits,
                            view,
                            plane,
                            shower.vals().at(s)->ShowerStart(),
                            shower.vals().at(s)->Direction(),
                            s,
                            -10001); //use -10001 to increase shower hit size

        } // end loop over prongs
      }   // end loop over labels
    }     // end draw showers

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawTrackVertexAssns2D(const art::Event& evt,
                                         detinfo::DetectorClocksData const& clockData,
                                         detinfo::DetectorPropertiesData const& detProp,
                                         evdb::View2D* view,
                                         unsigned int plane)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    if (!recoOpt->fDrawTrackVertexAssns) return;

    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    unsigned int cstat = rawOpt->fCryostat;
    unsigned int tpc = rawOpt->fTPC;
    geo::PlaneID planeID(cstat, tpc, plane);
    int tid = 0;

    for (size_t imod = 0; imod < recoOpt->fTrkVtxTrackLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fTrkVtxTrackLabels[imod];

      art::View<recob::Track> trackCol;
      this->GetTracks(evt, which, trackCol);

      if (trackCol.vals().size() < 1) continue;

      // Recover associations output from the filter
      std::unique_ptr<art::Assns<recob::Vertex, recob::Track>> vertexTrackAssociations(
        new art::Assns<recob::Vertex, recob::Track>);

      // Recover a handle to the collection of associations between vertices and tracks
      // This is a bit non-standard way to do this but trying to avoid complications
      art::Handle<art::Assns<recob::Vertex, recob::Track>> vertexTrackAssnsHandle;

      evt.getByLabel(recoOpt->fTrkVtxFilterLabels[imod], vertexTrackAssnsHandle);

      if (vertexTrackAssnsHandle->size() < 1) continue;

      // Get the rest of the associations in the standard way
      art::FindMany<recob::Hit> fmh(trackCol, evt, which);

      art::FindManyP<anab::CosmicTag> cosmicTrackTags(
        trackCol, evt, recoOpt->fTrkVtxCosmicLabels[imod]);

      auto tracksProxy = proxy::getCollection<proxy::Tracks>(evt, which);

      // Need to keep track of vertices unfortunately
      int lastVtxIdx(-1);
      int color(kRed);

      std::cout << "==> Neutrino Candidate drawing for tagger "
                << recoOpt->fTrkVtxFilterLabels[imod] << std::endl;

      // Now we can iterate over the vertex/track associations and do some drawing
      for (const auto& vertexTrackAssn : *vertexTrackAssnsHandle) {
        // Start by drawing the vertex
        art::Ptr<recob::Vertex> vertex = vertexTrackAssn.first;

        if (vertex->ID() != lastVtxIdx) {
          // BB: draw polymarker at the vertex position in this plane
          double xyz[3];

          vertex->XYZ(xyz);

          geo::Point_t localXYZ(xyz[0], xyz[1], xyz[2]);

          double wire = geo->WireCoordinate(localXYZ, planeID);
          double time = detProp.ConvertXToTicks(xyz[0], planeID);

          TMarker& strt = view->AddMarker(wire, time, color, 24, 3.0);
          strt.SetMarkerColor(color);

          std::cout << "    --> Drawing vertex id: " << vertex->ID() << std::endl;
        }

        lastVtxIdx = vertex->ID();

        const art::Ptr<recob::Track>& track = vertexTrackAssn.second;

        // BB: draw the track ID at the end of the track
        double x = track->End().X();
        geo::Point_t trackEnd(track->End());
        double tick = 30 + detProp.ConvertXToTicks(x, planeID);
        double wire = geo->WireCoordinate(trackEnd, planeID);

        tid = track->ID() & 65535;

        std::cout << "        --> Drawing Track id: " << tid << std::endl;

        std::string s = std::to_string(tid);
        char const* txt = s.c_str();

        TText& trkID = view->AddText(wire, tick, txt);
        trkID.SetTextColor(color);
        trkID.SetTextSize(0.1);

        float cosmicScore = -999;
        if (cosmicTrackTags.isValid()) {
          if (cosmicTrackTags.at(track.key()).size() > 0) {
            art::Ptr<anab::CosmicTag> currentTag = cosmicTrackTags.at(track.key()).at(0);
            cosmicScore = currentTag->CosmicScore();
          }
        }

        std::vector<const recob::Hit*> hits;
        if (track->NumberTrajectoryPoints() == fmh.at(track.key()).size()) {
          auto tp = tracksProxy[track.key()];
          for (auto point : tp.points()) {
            if (!point.isPointValid()) continue;
            hits.push_back(point.hit());
          }
        }
        else {
          hits = fmh.at(track.key());
        }
        // only get the hits for the current view
        std::vector<const recob::Hit*>::iterator itr = hits.begin();
        while (itr < hits.end()) {
          if ((*itr)->View() != gview)
            hits.erase(itr);
          else
            itr++;
        }

        int lineWidth(1);

        if (cosmicScore > 0.1) {
          color = kRed;
          if (cosmicScore < 0.6) color = kMagenta;
          lineWidth = 3;
        }
        else if (cosmicScore < -10000) { //shower hits
          lineWidth = 3;
        }

        this->DrawTrack2D(clockData, detProp, hits, view, plane, track.get(), color, lineWidth);

      } // end loop over vertex/track associations

    } // end loop over labels
  }

  //......................................................................
  void
  RecoBaseDrawer::Vertex2D(const art::Event& evt,
                           detinfo::DetectorPropertiesData const& detProp,
                           evdb::View2D* view,
                           unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawVertices == 0) return;

    art::ServiceHandle<geo::Geometry const> geo;
    static bool first = true;

    if (first) {
      std::cout << "******** DrawVertices: Open circles color coded across all planes. Set "
                   "DrawVertices > 1 to display the vertex ID\n";
      first = false;
    }

    for (size_t imod = 0; imod < recoOpt->fVertexLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fVertexLabels[imod];

      art::PtrVector<recob::Vertex> vertex;
      this->GetVertices(evt, which, vertex);

      if (vertex.size() < 1) continue;

      double local[3] = {0., 0., 0.};
      double world[3] = {0., 0., 0.};
      const geo::TPCGeo& tpc = geo->TPC(rawOpt->fTPC);
      tpc.LocalToWorld(local, world);
      double minxyz[3], maxxyz[3];
      minxyz[0] = world[0] - geo->DetHalfWidth(rawOpt->fTPC, rawOpt->fCryostat);
      maxxyz[0] = world[0] + geo->DetHalfWidth(rawOpt->fTPC, rawOpt->fCryostat);
      minxyz[1] = world[1] - geo->DetHalfWidth(rawOpt->fTPC, rawOpt->fCryostat);
      maxxyz[1] = world[1] + geo->DetHalfWidth(rawOpt->fTPC, rawOpt->fCryostat);
      minxyz[2] = world[2] - geo->DetLength(rawOpt->fTPC, rawOpt->fCryostat) / 2;
      maxxyz[2] = world[2] + geo->DetLength(rawOpt->fTPC, rawOpt->fCryostat) / 2;

      for (size_t v = 0; v < vertex.size(); ++v) {
        // ensure the vertex is inside the current tpc
        double xyz[3];
        vertex[v]->XYZ(xyz);
        if (xyz[0] < minxyz[0] || xyz[0] > maxxyz[0]) continue;
        if (xyz[1] < minxyz[1] || xyz[1] > maxxyz[1]) continue;
        if (xyz[2] < minxyz[2] || xyz[2] > maxxyz[2]) continue;

        geo::Point_t localPos(xyz[0], xyz[1], xyz[2]);

        // BB: draw polymarker at the vertex position in this plane
        double wire =
          geo->WireCoordinate(localPos, geo::PlaneID(rawOpt->fCryostat, rawOpt->fTPC, plane));
        double time = detProp.ConvertXToTicks(xyz[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
        int color = evd::kColor[vertex[v]->ID() % evd::kNCOLS];
        TMarker& strt = view->AddMarker(wire, time, color, 24, 1.0);
        strt.SetMarkerColor(color);

        // BB: draw the vertex ID
        if (recoOpt->fDrawVertices > 1) {
          std::string s = "3V" + std::to_string(vertex[v]->ID());
          char const* txt = s.c_str();
          TText& vtxID = view->AddText(wire, time + 30, txt);
          vtxID.SetTextColor(color);
          vtxID.SetTextSize(0.05);
        }
      } // end loop over vertices to draw from this label
    }   // end loop over vertex module lables

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Event2D(const art::Event& evt, evdb::View2D* view, unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    if (recoOpt->fDrawEvents != 0) {
      geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

      for (unsigned int imod = 0; imod < recoOpt->fEventLabels.size(); ++imod) {
        art::InputTag const which = recoOpt->fEventLabels[imod];

        art::PtrVector<recob::Event> event;
        this->GetEvents(evt, which, event);

        if (event.size() < 1) continue;

        art::FindMany<recob::Hit> fmh(event, evt, which);

        for (size_t e = 0; e < event.size(); ++e) {
          std::vector<const recob::Hit*> hits;

          hits = fmh.at(e);

          // only get the hits for the current view
          std::vector<const recob::Hit*>::iterator itr = hits.begin();
          while (itr < hits.end()) {
            if ((*itr)->View() != gview)
              hits.erase(itr);
            else
              itr++;
          }

          this->Hit2D(hits, evd::kColor[event[e]->ID() % evd::kNCOLS], view, false, true);
        } // end loop over events
      }   // end loop over event module lables
    }     // end if we are drawing events

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Seed3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    std::vector<art::InputTag> labels;
    if (recoOpt->fDrawSeeds != 0)
      for (size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod)
        labels.push_back(recoOpt->fSeedLabels[imod]);

    for (size_t imod = 0; imod < labels.size(); ++imod) {
      art::InputTag const which = labels[imod];

      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);

      int color = 0;

      if (seeds.size() < 1) continue;

      TPolyMarker3D& pmrk = view->AddPolyMarker3D(seeds.size(), color, 4, 1);

      for (size_t iseed = 0; iseed != seeds.size(); ++iseed) {
        double pt[3], pterr[3], dir[3], direrr[3];
        seeds.at(iseed)->GetPoint(pt, pterr);
        seeds.at(iseed)->GetDirection(dir, direrr);

        double end1[3], end2[3];
        for (int i = 0; i != 3; ++i) {
          end1[i] = pt[i] + dir[i];
          end2[i] = pt[i] - dir[i];
        }

        TPolyLine3D& pline = view->AddPolyLine3D(2, color, 2, 0);

        pmrk.SetPoint(iseed, pt[0], pt[1], pt[2]);
        pline.SetPoint(0, end1[0], end1[1], end1[2]);
        pline.SetPoint(1, end2[0], end2[1], end2[2]);
      } // end loop over seeds
    }   // end loop over module labels

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::SeedOrtho(const art::Event& evt, evd::OrthoProj_t proj, evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    std::vector<art::InputTag> labels;
    if (recoOpt->fDrawSeeds != 0)
      for (size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod)
        labels.push_back(recoOpt->fSeedLabels[imod]);

    for (size_t imod = 0; imod < labels.size(); ++imod) {
      art::InputTag const which = labels[imod];

      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);

      int color = 0;

      for (size_t iseed = 0; iseed != seeds.size(); ++iseed) {
        double pt[3], pterr[3], dir[3], direrr[3];
        seeds.at(iseed)->GetPoint(pt, pterr);
        seeds.at(iseed)->GetDirection(dir, direrr);

        double end1[3], end2[3];
        for (int i = 0; i != 3; ++i) {
          end1[i] = pt[i] + dir[i];
          end2[i] = pt[i] - dir[i];
        }

        if (proj == evd::kXY) {
          TMarker& strt = view->AddMarker(pt[1], pt[0], color, 4, 1.5);
          TLine& line = view->AddLine(end1[1], end1[0], end2[1], end2[0]);
          strt.SetMarkerColor(evd::kColor[color]);
          line.SetLineColor(evd::kColor[color]);
          line.SetLineWidth(2.0);
        }
        else if (proj == evd::kXZ) {
          TMarker& strt = view->AddMarker(pt[2], pt[0], color, 4, 1.5);
          TLine& line = view->AddLine(end1[2], end1[0], end2[2], end2[0]);
          strt.SetMarkerColor(evd::kColor[color]);
          line.SetLineColor(evd::kColor[color]);
          line.SetLineWidth(2.0);
        }
        else {
          if (proj != evd::kYZ)
            throw cet::exception("RecoBaseDrawer:SeedOrtho")
              << "projection is not YZ as expected\n";

          TMarker& strt = view->AddMarker(pt[2], pt[1], color, 4, 1.5);
          TLine& line = view->AddLine(end1[2], end1[1], end2[2], end2[1]);
          strt.SetMarkerColor(evd::kColor[color]);
          line.SetLineColor(evd::kColor[color]);
          line.SetLineWidth(2.0);
        }
      }
    }
  }

  //......................................................................
  void
  RecoBaseDrawer::SpacePoint3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    std::vector<art::InputTag> labels;
    if (recoOpt->fDrawSpacePoints != 0) {
      for (size_t imod = 0; imod < recoOpt->fSpacePointLabels.size(); ++imod)
        labels.push_back(recoOpt->fSpacePointLabels[imod]);
    }

    for (size_t imod = 0; imod < labels.size(); ++imod) {
      art::InputTag const which = labels[imod];

      std::vector<art::Ptr<recob::SpacePoint>> spts;
      this->GetSpacePoints(evt, which, spts);
      int color = 10 * imod + 11;

      color = 0;

      //        std::vector<const recob::SpacePoint*> sptsVec;
      //
      //        sptsVec.resize(spts.size());
      //        for(const auto& spt : spts){
      //          std::cout<<spt<<" "<<*spt<<" "<<&*spt<<std::endl;
      //          sptsVec.push_back(&*spt);
      //          std::cout<<sptsVec.back()<<std::endl;
      //        }
      fAllSpacePointDrawer->Draw(spts, view, color, kFullDotMedium, 1);
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::PFParticle3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawPFParticles < 1) return;

    // The plan is to loop over the list of possible particles
    for (size_t imod = 0; imod < recoOpt->fPFParticleLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fPFParticleLabels[imod];
      art::InputTag const assns = recoOpt->fSpacePointLabels[imod];

      // Start off by recovering our 3D Clusters for this label
      art::PtrVector<recob::PFParticle> pfParticleVec;
      this->GetPFParticles(evt, which, pfParticleVec);

      mf::LogDebug("RecoBaseDrawer")
        << "RecoBaseDrawer: number PFParticles to draw: " << pfParticleVec.size() << std::endl;

      // Make sure we have some clusters
      if (pfParticleVec.size() < 1) continue;

      // Get the space points created by the PFParticle producer
      std::vector<art::Ptr<recob::SpacePoint>> spacePointVec;
      this->GetSpacePoints(evt, assns, spacePointVec);

      // Recover the edges
      std::vector<art::Ptr<recob::Edge>> edgeVec;
      if (recoOpt->fDrawEdges) this->GetEdges(evt, assns, edgeVec);

      // No space points no continue
      if (spacePointVec.empty()) continue;

      // Add the relations to recover associations cluster hits
      art::FindManyP<recob::SpacePoint> edgeSpacePointAssnsVec(edgeVec, evt, assns);
      art::FindManyP<recob::SpacePoint> spacePointAssnVec(pfParticleVec, evt, assns);
      art::FindManyP<recob::Hit> spHitAssnVec(spacePointVec, evt, assns);
      art::FindManyP<recob::Edge> edgeAssnsVec(pfParticleVec, evt, assns);

      // If no valid space point associations then nothing to do
      if (!spacePointAssnVec.isValid()) continue;

      // Need the PCA info as well
      art::FindMany<recob::PCAxis> pcAxisAssnVec(pfParticleVec, evt, which);

      // Want CR tagging info
      // Note the cosmic tags come from a different producer - we assume that the producers are
      // matched in the fcl label vectors!
      art::InputTag cosmicTagLabel =
        imod < recoOpt->fCosmicTagLabels.size() ? recoOpt->fCosmicTagLabels[imod] : "";
      art::FindMany<anab::CosmicTag> pfCosmicAssns(pfParticleVec, evt, cosmicTagLabel);

      // We also want to drive display of tracks but have the same issue with production... so follow the
      // same prescription.
      art::InputTag trackTagLabel =
        imod < recoOpt->fTrackLabels.size() ? recoOpt->fTrackLabels[imod] : "";
      art::FindMany<recob::Track> pfTrackAssns(pfParticleVec, evt, trackTagLabel);

      // Commence looping over possible clusters
      for (size_t idx = 0; idx < pfParticleVec.size(); idx++) {
        // Recover cluster
        const art::Ptr<recob::PFParticle> pfParticle(pfParticleVec.at(idx));

        // Drawing will be done recursively down the chain of hieirarchy... so we want to begin
        // with only "primary" particles, if we find one that isn't then we skip
        if (!pfParticle->IsPrimary()) continue;

        // Call the recursive drawing routine
        DrawPFParticle3D(pfParticle,
                         pfParticleVec,
                         spacePointVec,
                         edgeAssnsVec,
                         spacePointAssnVec,
                         edgeSpacePointAssnsVec,
                         spHitAssnVec,
                         pfTrackAssns,
                         pcAxisAssnVec,
                         pfCosmicAssns,
                         0,
                         view);
      }
    }

    return;
  }

  float
  RecoBaseDrawer::SpacePointChiSq(const std::vector<art::Ptr<recob::Hit>>& hitVec) const
  {
    float hitChiSq(0.);

    bool usePlane[] = {false, false, false};
    float peakTimeVec[] = {0., 0., 0.};
    float peakSigmaVec[] = {0., 0., 0.};
    float aveSum(0.);
    float weightSum(0.);

    // Temp ad hoc correction to investigate...
    std::map<size_t, double> planeOffsetMap;

    planeOffsetMap[0] = 0.;
    planeOffsetMap[1] = 4.;
    planeOffsetMap[2] = 8.;

    for (const auto& hit : hitVec) {
      if (!hit) continue;

      float peakTime = hit->PeakTime() - planeOffsetMap[hit->WireID().Plane];
      float peakRMS = hit->RMS();

      aveSum += peakTime / (peakRMS * peakRMS);
      weightSum += 1. / (peakRMS * peakRMS);

      peakTimeVec[hit->WireID().Plane] = peakTime;
      peakSigmaVec[hit->WireID().Plane] = peakRMS;
      usePlane[hit->WireID().Plane] = true;
    }

    aveSum /= weightSum;

    for (int idx = 0; idx < 3; idx++) {
      if (usePlane[idx]) {
        float deltaTime = peakTimeVec[idx] - aveSum;
        float sigmaPeakTimeSq = peakSigmaVec[idx] * peakSigmaVec[idx];

        hitChiSq += deltaTime * deltaTime / sigmaPeakTimeSq;
      }
    }

    return hitChiSq;
  }

  void
  RecoBaseDrawer::DrawPFParticle3D(const art::Ptr<recob::PFParticle>& pfPart,
                                   const art::PtrVector<recob::PFParticle>& pfParticleVec,
                                   const std::vector<art::Ptr<recob::SpacePoint>>& spacePointVec,
                                   const art::FindManyP<recob::Edge>& edgeAssnsVec,
                                   const art::FindManyP<recob::SpacePoint>& spacePointAssnVec,
                                   const art::FindManyP<recob::SpacePoint>& edgeSPAssnVec,
                                   const art::FindManyP<recob::Hit>& spHitAssnVec,
                                   const art::FindMany<recob::Track>& trackAssnVec,
                                   const art::FindMany<recob::PCAxis>& pcAxisAssnVec,
                                   const art::FindMany<anab::CosmicTag>& cosmicTagAssnVec,
                                   int depth,
                                   evdb::View3D* view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<evd::ColorDrawingOptions const> cst;

    // First let's draw the hits associated to this cluster
    const std::vector<art::Ptr<recob::SpacePoint>>& hitsVec(spacePointAssnVec.at(pfPart.key()));

    // Use the particle ID to determine the color to draw the points
    // Ok, this is what we would like to do eventually but currently all particles are the same...
    bool isCosmic(false);
    int colorIdx(evd::kColor[pfPart->Self() % evd::kNCOLS]);

    // Recover cosmic tag info if any
    if (cosmicTagAssnVec.isValid() && recoOpt->fDrawPFParticles > 3) {
      std::vector<const anab::CosmicTag*> pfCosmicTagVec = cosmicTagAssnVec.at(pfPart.key());

      if (!pfCosmicTagVec.empty()) {
        const anab::CosmicTag* cosmicTag = pfCosmicTagVec.front();

        if (cosmicTag->CosmicScore() > 0.6) isCosmic = true;
      }
    }

    // Reset color index if a cosmic
    if (isCosmic) colorIdx = 12;

    if (!hitsVec.empty() && recoOpt->fDraw3DSpacePoints)
      fSpacePointDrawer->Draw(hitsVec, view, 1, kFullDotLarge, 0.25, &spHitAssnVec);
    /*
    {
        using HitPosition = std::array<double,6>;
        std::map<int,std::vector<HitPosition>> colorToHitMap;

        float minHitChiSquare(0.);
        float maxHitChiSquare(2.);
        float hitChiSqScale((cst->fRecoQHigh[geo::kCollection] - cst->fRecoQLow[geo::kCollection]) / (maxHitChiSquare - minHitChiSquare));

        for(const auto& spacePoint : hitsVec)
        {
            const double* pos = spacePoint->XYZ();
            const double* err = spacePoint->ErrXYZ();

            bool  storeHit(false);
            int   chargeColorIdx(0);
            float spacePointChiSq(spacePoint->Chisq());

            if (recoOpt->fDraw3DSpacePointHeatMap)
            {
                storeHit = true;

                float hitChiSq = std::max(minHitChiSquare, std::min(maxHitChiSquare, spacePointChiSq));

                float chgFactor = cst->fRecoQHigh[geo::kCollection] - hitChiSqScale * hitChiSq;
                //float chgFactor = delTScaleFctr * hitChiSq + cst->fRecoQLow[geo::kCollection];

                chargeColorIdx = cst->CalQ(geo::kCollection).GetColor(chgFactor);
            }
            else
            {
                if (spacePointChiSq > 0. && !recoOpt->fSkeletonOnly)         // All cluster hits which are unmarked
                {
                    storeHit = true;
                }
                else if (spacePointChiSq > -2.)                              // Skeleton hits
                {
                    chargeColorIdx = 5;
                    storeHit = true;
                }
                else if (spacePointChiSq > -3.)                              // Pure edge hits
                {
                    if (chargeColorIdx < 0) chargeColorIdx = !isCosmic ? 3 : colorIdx + 3;
                    storeHit = true;
                }
                else if (spacePointChiSq > -4.)                              // Skeleton hits which are also edge hits
                {
                    if (chargeColorIdx < 0) chargeColorIdx = !isCosmic ? 0 : colorIdx + 3;
                    storeHit = true;
                }
                else if (spacePoint->Chisq() > -5.)                             // Hits which form seeds for tracks
                {
                    if (chargeColorIdx < 0) chargeColorIdx = !isCosmic ? 5 : colorIdx + 3;
                    storeHit = true;
                }
                else if (!recoOpt->fSkeletonOnly)                                // hits made from pairs
                {
                    chargeColorIdx = 15;
                    storeHit       = true;
                }

                if (chargeColorIdx < 0) chargeColorIdx = 0;
            }

            if (storeHit) colorToHitMap[chargeColorIdx].push_back(HitPosition()={{pos[0],pos[1],pos[2],err[3],err[3],err[5]}});
        }

        size_t nHitsDrawn(0);

        for(auto& hitPair : colorToHitMap)
        {
            //TPolyMarker3D& pm = view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotMedium, 3);
            TPolyMarker3D& pm = view->AddPolyMarker3D(hitPair.second.size(), hitPair.first, kFullDotLarge, 0.25); //kFullDotLarge, 0.3);
            for (const auto& hit : hitPair.second) pm.SetNextPoint(hit[0],hit[1],hit[2]);

            nHitsDrawn += hitPair.second.size();
        }
    }
*/

    // Now try to draw any associated edges
    if (edgeAssnsVec.isValid() && recoOpt->fDraw3DEdges) {
      const std::vector<art::Ptr<recob::Edge>>& edgeVec(edgeAssnsVec.at(pfPart.key()));

      if (!edgeVec.empty()) {
        TPolyMarker3D& pm = view->AddPolyMarker3D(
          2 * edgeVec.size(), colorIdx, kFullDotMedium, 1.25); //kFullDotLarge, 0.5);

        for (const auto& edge : edgeVec) {
          try {
            const std::vector<art::Ptr<recob::SpacePoint>>& spacePointVec(
              edgeSPAssnVec.at(edge.key()));

            if (spacePointVec.size() != 2) {
              std::cout << "Space Point vector associated to edge is not of size 2: "
                        << spacePointVec.size() << std::endl;
              continue;
            }

            const recob::SpacePoint* firstSP = spacePointVec[0].get();
            const recob::SpacePoint* secondSP = spacePointVec[1].get();

            TVector3 startPoint(firstSP->XYZ()[0], firstSP->XYZ()[1], firstSP->XYZ()[2]);
            TVector3 endPoint(secondSP->XYZ()[0], secondSP->XYZ()[1], secondSP->XYZ()[2]);
            TVector3 lineVec(endPoint - startPoint);

            pm.SetNextPoint(startPoint[0], startPoint[1], startPoint[2]);
            pm.SetNextPoint(endPoint[0], endPoint[1], endPoint[2]);

            double length = lineVec.Mag();

            if (length == 0.) {
              std::cout << "Edge length is zero, index 1: " << edge->FirstPointID()
                        << ", index 2: " << edge->SecondPointID() << std::endl;
              continue;
            }

            double minLen = std::max(2.01, length);

            if (minLen > length) {
              lineVec.SetMag(1.);

              startPoint += -0.5 * (minLen - length) * lineVec;
              endPoint += 0.5 * (minLen - length) * lineVec;
            }

            // Get a polyline object to draw from the first to the second space point
            TPolyLine3D& pl = view->AddPolyLine3D(2, colorIdx, 4, 1);

            pl.SetPoint(0, startPoint[0], startPoint[1], startPoint[2]);
            pl.SetPoint(1, endPoint[0], endPoint[1], endPoint[2]);
          }
          catch (...) {
            continue;
          }
        }
      }
    }

    // Draw associated tracks
    if (trackAssnVec.isValid()) {
      std::vector<const recob::Track*> trackVec(trackAssnVec.at(pfPart.key()));

      if (!trackVec.empty()) {
        for (const auto& track : trackVec)
          DrawTrack3D(*track, view, colorIdx, kFullDotLarge, 0.5);
      }
    }

    // Look up the PCA info
    if (pcAxisAssnVec.isValid() && recoOpt->fDraw3DPCAAxes) {
      std::vector<const recob::PCAxis*> pcaVec(pcAxisAssnVec.at(pfPart.key()));

      if (!pcaVec.empty()) {
        // For each axis we are going to draw a solid line between two points
        int numPoints(2);
        //int lineWidth[2] = {       3,  1};
        int lineWidth[2] = {2, 1};
        int lineStyle[2] = {1, 13};
        int lineColor[2] = {colorIdx, 18};
        //int markStyle[2] = {       4,  4};
        int markStyle[2] = {kFullDotLarge, kFullDotLarge};
        double markSize[2] = {0.5, 0.2};
        int pcaIdx(0);

        if (!isCosmic) lineColor[1] = colorIdx;

        // The order of axes in the returned association vector is arbitrary... the "first" axis is
        // better and we can divine that by looking at the axis id's (the best will have been made first)
        if (pcaVec.size() > 1 && pcaVec.front()->getID() > pcaVec.back()->getID())
          std::reverse(pcaVec.begin(), pcaVec.end());

        for (const auto& pca : pcaVec) {
          // We need the mean position
          const double* avePosition = pca->getAvePosition();

          // Let's draw a marker at the interesting points
          int pmrkIdx(0);
          TPolyMarker3D& pmrk =
            view->AddPolyMarker3D(7, lineColor[pcaIdx], markStyle[pcaIdx], markSize[pcaIdx]);

          pmrk.SetPoint(pmrkIdx++, avePosition[0], avePosition[1], avePosition[2]);

          // Loop over pca dimensions
          for (int dimIdx = 0; dimIdx < 3; dimIdx++) {
            // Oh please oh please give me an instance of a poly line...
            TPolyLine3D& pl = view->AddPolyLine3D(
              numPoints, lineColor[pcaIdx], lineWidth[pcaIdx], lineStyle[pcaIdx]);

            // We will use the eigen value to give the length of the line we're going to plot
            double eigenValue = pca->getEigenValues()[dimIdx];

            // Make sure a valid eigenvalue
            if (eigenValue > 0) {
              // Really want the root of the eigen value
              eigenValue = 3. * sqrt(eigenValue);

              // Recover the eigenvector
              const std::vector<double>& eigenVector = pca->getEigenVectors()[dimIdx];

              // Set the first point
              double xl = avePosition[0] - 0.5 * eigenValue * eigenVector[0];
              double yl = avePosition[1] - 0.5 * eigenValue * eigenVector[1];
              double zl = avePosition[2] - 0.5 * eigenValue * eigenVector[2];

              pl.SetPoint(0, xl, yl, zl);
              pmrk.SetPoint(pmrkIdx++, xl, yl, zl);

              // Set the second point
              double xu = avePosition[0] + 0.5 * eigenValue * eigenVector[0];
              double yu = avePosition[1] + 0.5 * eigenValue * eigenVector[1];
              double zu = avePosition[2] + 0.5 * eigenValue * eigenVector[2];

              pl.SetPoint(1, xu, yu, zu);
              pmrk.SetPoint(pmrkIdx++, xu, yu, zu);
            }
          }

          // By convention we will have drawn the "best" axis first
          if (recoOpt->fBestPCAAxisOnly) break;

          pcaIdx++;
        }
      }
    }

    // Now let's loop over daughters and call drawing routine for them
    if (pfPart->NumDaughters() > 0) {
      depth++;

      //        std::string indent(depth, ' ');

      //        std::cout << indent << "-drawPFParticle3D--> pfPart idx: " << pfPart->Self() << ", daughter list size: " << pfPart->Daughters().size() << std::endl;

      for (const auto& daughterIdx : pfPart->Daughters()) {
        DrawPFParticle3D(pfParticleVec.at(daughterIdx),
                         pfParticleVec,
                         spacePointVec,
                         edgeAssnsVec,
                         spacePointAssnVec,
                         edgeSPAssnVec,
                         spHitAssnVec,
                         trackAssnVec,
                         pcAxisAssnVec,
                         cosmicTagAssnVec,
                         depth,
                         view);
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Edge3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDrawEdges < 1) return;

    // The plan is to loop over the list of possible particles
    for (size_t imod = 0; imod < recoOpt->fEdgeLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fEdgeLabels[imod];

      // Start off by recovering our 3D Clusters for this label
      std::vector<art::Ptr<recob::Edge>> edgeVec;
      this->GetEdges(evt, which, edgeVec);

      mf::LogDebug("RecoBaseDrawer")
        << "RecoBaseDrawer: number Edges to draw: " << edgeVec.size() << std::endl;

      if (!edgeVec.empty()) {
        // Get the space points created by the PFParticle producer
        std::vector<art::Ptr<recob::SpacePoint>> spacePointVec;
        this->GetSpacePoints(evt, which, spacePointVec);

        // First draw the space points (all of them), then circle back on the edges...
        int colorIdx(41); //2);

        TPolyMarker3D& pm = view->AddPolyMarker3D(
          spacePointVec.size(), colorIdx, kFullDotMedium, 0.5); //kFullDotLarge, 0.5);

        for (const auto& spacePoint : spacePointVec) {
          TVector3 spPosition(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]);

          pm.SetNextPoint(spPosition[0], spPosition[1], spPosition[2]);
        }

        // Now draw the edges
        for (const auto& edge : edgeVec) {
          art::Ptr<recob::SpacePoint> firstSP = spacePointVec.at(edge->FirstPointID());
          art::Ptr<recob::SpacePoint> secondSP = spacePointVec.at(edge->SecondPointID());

          if (firstSP->ID() != edge->FirstPointID() || secondSP->ID() != edge->SecondPointID()) {
            mf::LogDebug("RecoBaseDrawer")
              << "Edge: Space point index mismatch, first: " << firstSP->ID() << ", "
              << edge->FirstPointID() << ", second: " << secondSP->ID() << ", "
              << edge->SecondPointID() << std::endl;
            continue;
          }

          TVector3 startPoint(firstSP->XYZ()[0], firstSP->XYZ()[1], firstSP->XYZ()[2]);
          TVector3 endPoint(secondSP->XYZ()[0], secondSP->XYZ()[1], secondSP->XYZ()[2]);
          TVector3 lineVec(endPoint - startPoint);

          double length = lineVec.Mag();

          if (length == 0.) {
            //                    std::cout << "Edge length is zero, index 1: " << edge->FirstPointID() << ", index 2: " << edge->SecondPointID() << std::endl;
            continue;
          }

          // Get a polyline object to draw from the first to the second space point
          //                TPolyLine3D& pl = view->AddPolyLine3D(2, colorIdx, 1, 1); //4, 1);
          //
          //                pl.SetPoint(0, startPoint[0], startPoint[1], startPoint[2]);
          //                pl.SetPoint(1, endPoint[0],   endPoint[1],   endPoint[2]);
          TPolyMarker3D& fakeLine = view->AddPolyMarker3D(10, 5, kFullDotMedium, 1.0);

          lineVec.SetMag(1.);

          for (int idx = 1; idx <= 10; idx++) {
            TVector3 plotPoint = startPoint + 0.1 * double(idx) * length * lineVec;

            fakeLine.SetNextPoint(plotPoint[0], plotPoint[1], plotPoint[2]);
          }
        }
      }
    }

    // Draw any associated Extreme Points
    for (size_t imod = 0; imod < recoOpt->fExtremePointLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fExtremePointLabels[imod];

      // Start off by recovering our 3D Clusters for this label
      std::vector<art::Ptr<recob::SpacePoint>> spacePointVec;
      this->GetSpacePoints(evt, which, spacePointVec);

      mf::LogDebug("RecoBaseDrawer")
        << "RecoBaseDrawer: number Extreme points to draw: " << spacePointVec.size() << std::endl;

      if (!spacePointVec.empty()) {
        // First draw the space points (all of them), then circle back on the edges...
        int colorIdx(kYellow);

        TPolyMarker3D& pm = view->AddPolyMarker3D(
          spacePointVec.size(), colorIdx, kFullDotLarge, 1.0); //kFullDotLarge, 0.5);

        for (const auto& spacePoint : spacePointVec) {
          TVector3 spPosition(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]);

          pm.SetNextPoint(spPosition[0], spPosition[1], spPosition[2]);
        }
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Prong3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    // Tracks.

    if (recoOpt->fDrawTracks > 2) {
      for (size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
        art::InputTag which = recoOpt->fTrackLabels[imod];
        art::View<recob::Track> trackView;
        this->GetTracks(evt, which, trackView);
        if (!trackView.isValid())
          continue; //Prevent potential segmentation fault if no tracks found. aoliv23@lsu.edu

        art::PtrVector<recob::Track> trackVec;

        trackView.fill(trackVec);

        art::InputTag const cosmicTagLabel(
          recoOpt->fCosmicTagLabels.size() > imod ? recoOpt->fCosmicTagLabels[imod] : "");
        art::FindMany<anab::CosmicTag> cosmicTagAssnVec(trackVec, evt, cosmicTagLabel);

        for (const auto& track : trackVec) {
          int color = evd::kColor[track.key() % evd::kNCOLS];
          int marker = kFullDotLarge;
          float size = 2.0;

          // Check if a CosmicTag object is available

          // Recover cosmic tag info if any
          if (cosmicTagAssnVec.isValid()) {
            std::vector<const anab::CosmicTag*> tkCosmicTagVec = cosmicTagAssnVec.at(track.key());

            if (!tkCosmicTagVec.empty()) {
              const anab::CosmicTag* cosmicTag = tkCosmicTagVec.front();

              // If tagged as Cosmic then neutralize the color
              if (cosmicTag->CosmicScore() > 0.6) {
                color = 14;
                size = 0.5;
              }
            }
          }

          // Draw track using only embedded information.

          DrawTrack3D(*track, view, color, marker, size);
        }
      }
    }

    // Showers.

    if (recoOpt->fDrawShowers != 0) {
      for (size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
        art::InputTag which = recoOpt->fShowerLabels[imod];
        art::View<recob::Shower> shower;
        this->GetShowers(evt, which, shower);

        for (size_t s = 0; s < shower.vals().size(); ++s) {
          const recob::Shower* pshower = shower.vals().at(s);
          int color = pshower->ID();
          DrawShower3D(*pshower, color, view);
        }
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawTrack3D(const recob::Track& track,
                              evdb::View3D* view,
                              int color,
                              int marker,
                              float size)
  {
    // Get options.
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDrawTrackSpacePoints) {
      // Use brute force to find the module label and index of this
      // track, so that we can find associated space points and draw
      // them.
      const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
      std::vector<art::Handle<std::vector<recob::Track>>> handles;
      evt->getManyByType(handles);

      for (auto ih : handles) {
        const art::Handle<std::vector<recob::Track>> handle = ih;

        if (handle.isValid()) {
          const std::string& which = handle.provenance()->moduleLabel();
          art::FindManyP<recob::SpacePoint> fmsp(handle, *evt, which);

          if (fmsp.isValid() && fmsp.size() > 0) {
            int n = handle->size();
            float spSize = 0.5 * size;

            for (int i = 0; i < n; ++i) {
              art::Ptr<recob::Track> p(handle, i);
              if (&*p == &track) {
                std::vector<art::Ptr<recob::SpacePoint>> spts = fmsp.at(i);
                fSpacePointDrawer->Draw(spts, view, color, marker, spSize);
              }
            }
          }
        }
      }
    }

    if (recoOpt->fDrawTrackTrajectoryPoints) {
      // Draw trajectory points.
      int np = track.NumberTrajectoryPoints();

      int lineSize = size;

      if (lineSize < 1) lineSize = 1;

      // Make and fill a special polymarker for the head of the track
      //TPolyMarker3D& pmStart = view->AddPolyMarker3D(1, color, 4, 3);
      TPolyMarker3D& pmStart = view->AddPolyMarker3D(1, 0, marker, 2. * size);

      const auto& firstPos = track.LocationAtPoint(0);
      pmStart.SetPoint(0, firstPos.X(), firstPos.Y(), firstPos.Z());

      // Make and fill a polymarker.
      TPolyMarker3D& pm = view->AddPolyMarker3D(track.CountValidPoints(), color, 1, 3);

      for (int p = 0; p < np; ++p) {
        if (!track.HasValidPoint(p)) continue;

        const auto& pos = track.LocationAtPoint(p);
        pm.SetPoint(p, pos.X(), pos.Y(), pos.Z());
      }

      // As we are a track, should we not be drawing a line here?
      TPolyLine3D& pl = view->AddPolyLine3D(track.CountValidPoints(), color, lineSize, 7);

      for (int p = 0; p < np; ++p) {
        if (!track.HasValidPoint(p)) continue;

        const auto pos = track.LocationAtPoint(p);

        pl.SetPoint(p, pos.X(), pos.Y(), pos.Z());
      }

      if (recoOpt->fDrawTrackTrajectoryPoints > 4) {
        // Can we add the track direction at each point?
        // This won't work for the last point... but let's try
        auto startPos = track.LocationAtPoint(0);
        auto startDir = track.DirectionAtPoint(0);

        for (int p = 1; p < np; ++p) {
          if (!track.HasValidPoint(p)) continue;

          TPolyLine3D& pl = view->AddPolyLine3D(2, (color + 1) % evd::kNCOLS, size, 7); //1, 3);

          auto nextPos = track.LocationAtPoint(p);
          auto deltaPos = nextPos - startPos;
          double arcLen = deltaPos.Dot(
            startDir); // arc len to plane containing next point perpendicular to current point dir

          if (arcLen < 0.) arcLen = 3.;

          auto endPoint = startPos + arcLen * startDir;

          pl.SetPoint(0, startPos.X(), startPos.Y(), startPos.Z());
          pl.SetPoint(1, endPoint.X(), endPoint.Y(), endPoint.Z());

          startPos = nextPos;
          startDir = track.DirectionAtPoint(p);
        }
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawShower3D(const recob::Shower& shower, int color, evdb::View3D* view)
  {
    // Use brute force to find the module label and index of this
    // shower, so that we can find associated space points and draw
    // them.
    // B. Baller: Catch an exception if there are no space points and draw a cone instead.

    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    std::vector<art::Handle<std::vector<recob::Shower>>> handles;
    evt->getManyByType(handles);

    bool noSpts = false;

    for (auto ih : handles) {
      const art::Handle<std::vector<recob::Shower>> handle = ih;

      if (handle.isValid()) {

        const std::string& which = handle.provenance()->moduleLabel();
        art::FindManyP<recob::SpacePoint> fmsp(handle, *evt, which);

        int n = handle->size();
        for (int i = 0; i < n; ++i) {
          art::Ptr<recob::Shower> p(handle, i);
          if (&*p == &shower) {
            // BB catch if no space points
            std::vector<art::Ptr<recob::SpacePoint>> spts;
            try {
              spts = fmsp.at(i);
              fSpacePointDrawer->Draw(spts, view, color);
            }
            catch (...) {
              noSpts = true;
              continue;
            } // catch
          }   // shower
        }     // i
      }       // ih
    }

    if (noSpts && shower.has_length() && shower.has_open_angle()) {
      std::cout << "No space points associated with the shower. Drawing a cone instead\n";
      color = kRed;
      auto& dedx = shower.dEdx();
      if (!dedx.empty()) {
        double dedxAve = 0;
        for (auto& dedxInPln : dedx)
          dedxAve += dedxInPln;
        dedxAve /= (double)dedx.size();
        // Use blue for ~1 MIP
        color = kBlue;
        // use green for ~2 MIP
        if (dedxAve > 3 && dedxAve < 5) color = kGreen;
      }
      double radius = shower.Length() * shower.OpenAngle();
      TVector3 startPos = shower.ShowerStart();
      TVector3 endPos = startPos + shower.Length() * shower.Direction();
      auto coneRim = Circle3D(endPos, shower.Direction(), radius);
      TPolyLine3D& pl = view->AddPolyLine3D(coneRim.size(), color, 2, 0);
      for (unsigned short ipt = 0; ipt < coneRim.size(); ++ipt) {
        auto& pt = coneRim[ipt];
        pl.SetPoint(ipt, pt[0], pt[1], pt[2]);
      }
      // Draw a line from the start position to each point on the rim
      for (unsigned short ipt = 0; ipt < coneRim.size(); ++ipt) {
        TPolyLine3D& panel = view->AddPolyLine3D(2, color, 2, 0);
        panel.SetPoint(0, startPos.X(), startPos.Y(), startPos.Z());
        panel.SetPoint(1, coneRim[ipt][0], coneRim[ipt][1], coneRim[ipt][2]);
      } //  ipt

    } // no space points

    return;
  }

  //......................................................................
  std::vector<std::array<double, 3>>
  RecoBaseDrawer::Circle3D(const TVector3& centerPos, const TVector3& axisDir, const double& radius)
  {
    // B. Baller Create a polyline circle in 3D

    // Make the rotation matrix to transform into the circle coordinate system
    TRotation r;
    r.RotateX(axisDir.X());
    r.RotateY(axisDir.Y());
    r.RotateZ(axisDir.Z());
    constexpr unsigned short nRimPts = 16;
    std::vector<std::array<double, 3>> rimPts(nRimPts + 1);
    for (unsigned short iang = 0; iang < nRimPts; ++iang) {
      double rimAngle = iang * 2 * M_PI / (float)nRimPts;
      TVector3 rim = {0, 0, 1};
      rim.SetX(radius * cos(rimAngle));
      rim.SetY(radius * sin(rimAngle));
      rim.SetZ(0);
      rim.Transform(r);
      rim += centerPos;
      for (unsigned short ixyz = 0; ixyz < 3; ++ixyz)
        rimPts[iang][ixyz] = rim[ixyz];
    } // iang
    // close the circle
    rimPts[nRimPts] = rimPts[0];
    return rimPts;
  } // PolyLineCircle

  //......................................................................
  void
  RecoBaseDrawer::Vertex3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    if (recoOpt->fDrawVertices != 0) {

      for (size_t imod = 0; imod < recoOpt->fVertexLabels.size(); ++imod) {
        art::InputTag const which = recoOpt->fVertexLabels[imod];

        art::PtrVector<recob::Vertex> vertex;
        this->GetVertices(evt, which, vertex);

        art::FindManyP<recob::Track> fmt(vertex, evt, which);
        art::FindManyP<recob::Shower> fms(vertex, evt, which);

        for (size_t v = 0; v < vertex.size(); ++v) {

          if (fmt.isValid()) {
            std::vector<art::Ptr<recob::Track>> tracks = fmt.at(v);

            // grab the Prongs from the vertex and draw those
            for (size_t t = 0; t < tracks.size(); ++t)
              this->DrawTrack3D(*(tracks[t]), view, vertex[v]->ID());
          }

          if (fms.isValid()) {
            std::vector<art::Ptr<recob::Shower>> showers = fms.at(v);

            for (size_t s = 0; s < showers.size(); ++s)
              this->DrawShower3D(*(showers[s]), vertex[v]->ID(), view);
          }

          double xyz[3] = {0.};
          vertex[v]->XYZ(xyz);
          TPolyMarker3D& pm =
            view->AddPolyMarker3D(1, evd::kColor[vertex[v]->ID() % evd::kNCOLS], 29, 6);
          pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);

        } // end loop over vertices to draw from this label
      }   // end loop over vertex module lables
    }     // end if we are drawing vertices

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::Event3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawEvents != 0) {

      for (size_t imod = 0; imod < recoOpt->fEventLabels.size(); ++imod) {
        art::InputTag const which = recoOpt->fEventLabels[imod];

        art::PtrVector<recob::Event> event;
        this->GetEvents(evt, which, event);

        if (event.size() < 1) continue;

        art::FindManyP<recob::Vertex> fmvp(event, evt, which);
        art::FindMany<recob::Vertex> fmv(event, evt, which);

        for (size_t e = 0; e < event.size(); ++e) {

          // grab the vertices for this event
          std::vector<art::Ptr<recob::Vertex>> vertex = fmvp.at(e);

          if (vertex.size() < 1) continue;

          art::FindManyP<recob::Track> fmt(vertex, evt, recoOpt->fVertexLabels[0]);
          art::FindManyP<recob::Shower> fms(vertex, evt, recoOpt->fVertexLabels[0]);

          for (size_t v = 0; v < vertex.size(); ++v) {

            /// \todo need a better way to grab the vertex module labels,
            // right now assume there is only 1 in the list
            std::vector<art::Ptr<recob::Track>> tracks = fmt.at(v);
            std::vector<art::Ptr<recob::Shower>> showers = fms.at(v);

            // grab the Prongs from the vertex and draw those
            for (size_t t = 0; t < tracks.size(); ++t)
              this->DrawTrack3D(*(tracks[t]), view, event[e]->ID());

            for (size_t s = 0; s < showers.size(); ++s)
              this->DrawShower3D(*(showers[s]), event[e]->ID(), view);

          } // end loop over vertices from this event

          double xyz[3] = {0.};
          std::vector<const recob::Vertex*> vts = fmv.at(e);

          event[e]->PrimaryVertex(vts)->XYZ(xyz);
          TPolyMarker3D& pm =
            view->AddPolyMarker3D(1, evd::kColor[event[e]->ID() % evd::kNCOLS], 29, 6);
          pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);

        } // end loop over events
      }   // end loop over event module lables
    }     // end if we are drawing events

    return;
  }
  //......................................................................
  void
  RecoBaseDrawer::Slice3D(const art::Event& evt, evdb::View3D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawSlices < 1) return;
    if (recoOpt->fDrawSliceSpacePoints < 1) return;
    for (size_t imod = 0; imod < recoOpt->fSliceLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fSliceLabels[imod];
      art::PtrVector<recob::Slice> slices;
      this->GetSlices(evt, which, slices);
      if (slices.size() < 1) continue;
      art::FindManyP<recob::SpacePoint> fmsp(slices, evt, which);
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        int slcID = std::abs(slices[isl]->ID());
        int color = evd::kColor[slcID % evd::kNCOLS];
        std::vector<art::Ptr<recob::SpacePoint>> spts = fmsp.at(isl);
        fSpacePointDrawer->Draw(spts, view, color, kFullDotLarge, 2);
      }
    }
  }
  //......................................................................
  void
  RecoBaseDrawer::OpFlashOrtho(const art::Event& evt,
                               detinfo::DetectorClocksData const& clockData,
                               detinfo::DetectorPropertiesData const& detProp,
                               evd::OrthoProj_t proj,
                               evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawOpFlashes == 0) return;

    geo::PlaneID pid(rawOpt->fCryostat, rawOpt->fTPC, 0);

    double minx = 1e9;
    double maxx = -1e9;
    for (size_t i = 0; i < geo->NTPC(); ++i) {
      double local[3] = {0., 0., 0.};
      double world[3] = {0., 0., 0.};
      const geo::TPCGeo& tpc = geo->TPC(i);
      tpc.LocalToWorld(local, world);
      if (minx > world[0] - geo->DetHalfWidth(i)) minx = world[0] - geo->DetHalfWidth(i);
      if (maxx < world[0] + geo->DetHalfWidth(i)) maxx = world[0] + geo->DetHalfWidth(i);
    }

    for (size_t imod = 0; imod < recoOpt->fOpFlashLabels.size(); ++imod) {
      const art::InputTag which = recoOpt->fOpFlashLabels[imod];

      art::PtrVector<recob::OpFlash> opflashes;
      this->GetOpFlashes(evt, which, opflashes);

      if (opflashes.size() < 1) continue;

      int NFlashes = opflashes.size();

      // project each seed into this view
      for (int iof = 0; iof < NFlashes; ++iof) {

        if (opflashes[iof]->TotalPE() < recoOpt->fFlashMinPE) continue;
        if (opflashes[iof]->Time() < recoOpt->fFlashTMin) continue;
        if (opflashes[iof]->Time() > recoOpt->fFlashTMax) continue;

        double YCentre = opflashes[iof]->YCenter();
        double YHalfWidth = opflashes[iof]->YWidth();
        double ZCentre = opflashes[iof]->ZCenter();
        double ZHalfWidth = opflashes[iof]->ZWidth();

        int Colour = evd::kColor[(iof) % evd::kNCOLS];

        if (proj == evd::kXY) {
          TBox& b1 = view->AddBox(YCentre - YHalfWidth, minx, YCentre + YHalfWidth, maxx);
          b1.SetFillStyle(3004 + (iof % 3));
          b1.SetFillColor(Colour);
        }
        else if (proj == evd::kXZ) {
          float xflash = detProp.ConvertTicksToX(
            opflashes[iof]->Time() / sampling_rate(clockData) * 1e3 + detProp.GetXTicksOffset(pid),
            pid);
          TLine& line = view->AddLine(ZCentre - ZHalfWidth, xflash, ZCentre + ZHalfWidth, xflash);
          line.SetLineWidth(2);
          line.SetLineStyle(2);
          line.SetLineColor(Colour);
        }
        else if (proj == evd::kYZ) {
          TBox& b1 = view->AddBox(
            ZCentre - ZHalfWidth, YCentre - YHalfWidth, ZCentre + ZHalfWidth, YCentre + YHalfWidth);
          b1.SetFillStyle(3004 + (iof % 3));
          b1.SetFillColor(Colour);
          view->AddMarker(ZCentre, YCentre, Colour, 4, 1.5);
        }

      } // Flashes with this label
    }   // Vector of OpFlash labels
  }
  //......................................................................
  void
  RecoBaseDrawer::VertexOrtho(const art::PtrVector<recob::Vertex>& vertex,
                              evd::OrthoProj_t proj,
                              evdb::View2D* view,
                              int marker)
  {
    for (size_t v = 0; v < vertex.size(); ++v) {

      double xyz[3] = {0.};
      vertex[v]->XYZ(xyz);

      int color = evd::kColor[vertex[v]->ID() % evd::kNCOLS];

      if (proj == evd::kXY) {
        TMarker& strt = view->AddMarker(xyz[1], xyz[0], color, marker, 1.0);
        strt.SetMarkerColor(color);
      }
      else if (proj == evd::kXZ) {
        TMarker& strt = view->AddMarker(xyz[2], xyz[0], color, marker, 1.0);
        strt.SetMarkerColor(color);
      }
      else if (proj == evd::kYZ) {
        TMarker& strt = view->AddMarker(xyz[2], xyz[1], color, marker, 1.0);
        strt.SetMarkerColor(color);
      }
    }
  }
  void
  RecoBaseDrawer::VertexOrtho(const art::Event& evt, evd::OrthoProj_t proj, evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawVertices == 0) return;

    for (size_t imod = 0; imod < recoOpt->fVertexLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fVertexLabels[imod];

      art::PtrVector<recob::Vertex> vertex;
      this->GetVertices(evt, which, vertex);
      this->VertexOrtho(vertex, proj, view, 24);

      this->GetVertices(evt, art::InputTag(which.label(), "kink", which.process()), vertex);
      this->VertexOrtho(vertex, proj, view, 27);

      this->GetVertices(evt, art::InputTag(which.label(), "node", which.process()), vertex);
      this->VertexOrtho(vertex, proj, view, 22);
    }
    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::SpacePointOrtho(const art::Event& evt,
                                  evd::OrthoProj_t proj,
                                  double msize,
                                  evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    std::vector<art::InputTag> labels;
    if (recoOpt->fDrawSpacePoints != 0) {
      for (size_t imod = 0; imod < recoOpt->fSpacePointLabels.size(); ++imod)
        labels.push_back(recoOpt->fSpacePointLabels[imod]);
    }

    for (size_t imod = 0; imod < labels.size(); ++imod) {
      art::InputTag const which = labels[imod];

      std::vector<art::Ptr<recob::SpacePoint>> spts;
      this->GetSpacePoints(evt, which, spts);
      int color = imod;

      this->DrawSpacePointOrtho(spts, color, proj, msize, view);
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::PFParticleOrtho(const art::Event& evt,
                                  evd::OrthoProj_t proj,
                                  double msize,
                                  evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawPFParticles < 1) return;

    // The plan is to loop over the list of possible particles
    for (size_t imod = 0; imod < recoOpt->fPFParticleLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fPFParticleLabels[imod];

      // Start off by recovering our 3D Clusters for this label
      art::PtrVector<recob::PFParticle> pfParticleVec;
      this->GetPFParticles(evt, which, pfParticleVec);

      // Make sure we have some clusters
      if (pfParticleVec.size() < 1) continue;

      // Add the relations to recover associations cluster hits
      art::FindMany<recob::SpacePoint> spacePointAssnVec(pfParticleVec, evt, which);

      // If no valid space point associations then nothing to do
      if (!spacePointAssnVec.isValid()) continue;

      // Need the PCA info as well
      art::FindMany<recob::PCAxis> pcAxisAssnVec(pfParticleVec, evt, which);

      if (!pcAxisAssnVec.isValid()) continue;

      // Commence looping over possible clusters
      for (size_t idx = 0; idx < pfParticleVec.size(); idx++) {
        // Recover cluster
        const art::Ptr<recob::PFParticle> pfParticle(pfParticleVec.at(idx));

        // Drawing will be done recursively down the chain of hieirarchy... so we want to begin
        // with only "primary" particles, if we find one that isn't then we skip
        if (!pfParticle->IsPrimary()) continue;

        // Call the recursive drawing routine
        DrawPFParticleOrtho(
          pfParticle, pfParticleVec, spacePointAssnVec, pcAxisAssnVec, 0, proj, view);
      }
    }

    return;
  }

  void
  RecoBaseDrawer::DrawPFParticleOrtho(const art::Ptr<recob::PFParticle>& pfPart,
                                      const art::PtrVector<recob::PFParticle>& pfParticleVec,
                                      const art::FindMany<recob::SpacePoint>& spacePointAssnVec,
                                      const art::FindMany<recob::PCAxis>& pcAxisAssnVec,
                                      int depth,
                                      evd::OrthoProj_t proj,
                                      evdb::View2D* view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    // First let's draw the hits associated to this cluster
    const std::vector<const recob::SpacePoint*>& hitsVec(spacePointAssnVec.at(pfPart->Self()));

    // Use the particle ID to determine the color to draw the points
    // Ok, this is what we would like to do eventually but currently all particles are the same...
    //        int colorIdx = evd::Style::ColorFromPDG(pfPart->PdgCode());
    int colorIdx = evd::kColor[pfPart->Self() % evd::kNCOLS];

    if (!hitsVec.empty()) {
      std::vector<const recob::SpacePoint*> hitPosVec;
      std::vector<const recob::SpacePoint*> skeletonPosVec;
      std::vector<const recob::SpacePoint*> skelEdgePosVec;
      std::vector<const recob::SpacePoint*> edgePosVec;
      std::vector<const recob::SpacePoint*> seedPosVec;
      std::vector<const recob::SpacePoint*> pairPosVec;

      for (const auto& spacePoint : hitsVec) {
        if (spacePoint->Chisq() > 0.)
          hitPosVec.push_back(spacePoint);
        else if (spacePoint->Chisq() == -1.)
          skeletonPosVec.push_back(spacePoint);
        else if (spacePoint->Chisq() == -3.)
          skelEdgePosVec.push_back(spacePoint);
        else if (spacePoint->Chisq() == -4.)
          seedPosVec.push_back(spacePoint);
        else if (spacePoint->Chisq() > -10.)
          edgePosVec.push_back(spacePoint);
        else
          pairPosVec.push_back(spacePoint);
      }

      int hitIdx(0);

      if (!recoOpt->fSkeletonOnly) {
        TPolyMarker& pm1 = view->AddPolyMarker(hitPosVec.size(), colorIdx, kFullDotMedium, 1);
        for (const auto* spacePoint : hitPosVec) {
          const double* pos = spacePoint->XYZ();

          if (proj == evd::kXY)
            pm1.SetPoint(hitIdx++, pos[0], pos[1]);
          else if (proj == evd::kXZ)
            pm1.SetPoint(hitIdx++, pos[2], pos[0]);
          else if (proj == evd::kYZ)
            pm1.SetPoint(hitIdx++, pos[2], pos[1]);
        }

        hitIdx = 0;

        TPolyMarker& pm2 = view->AddPolyMarker(edgePosVec.size(), 28, kFullDotMedium, 1);
        for (const auto* spacePoint : edgePosVec) {
          const double* pos = spacePoint->XYZ();

          if (proj == evd::kXY)
            pm2.SetPoint(hitIdx++, pos[0], pos[1]);
          else if (proj == evd::kXZ)
            pm2.SetPoint(hitIdx++, pos[2], pos[0]);
          else if (proj == evd::kYZ)
            pm2.SetPoint(hitIdx++, pos[2], pos[1]);
        }

        hitIdx = 0;

        TPolyMarker& pm3 = view->AddPolyMarker(pairPosVec.size(), 2, kFullDotMedium, 1);
        for (const auto* spacePoint : pairPosVec) {
          const double* pos = spacePoint->XYZ();

          if (proj == evd::kXY)
            pm3.SetPoint(hitIdx++, pos[0], pos[1]);
          else if (proj == evd::kXZ)
            pm3.SetPoint(hitIdx++, pos[2], pos[0]);
          else if (proj == evd::kYZ)
            pm3.SetPoint(hitIdx++, pos[2], pos[1]);
        }
      }

      hitIdx = 0;

      TPolyMarker& pm4 = view->AddPolyMarker(skeletonPosVec.size(), 1, kFullDotMedium, 1);
      for (const auto* spacePoint : skeletonPosVec) {
        const double* pos = spacePoint->XYZ();

        if (proj == evd::kXY)
          pm4.SetPoint(hitIdx++, pos[0], pos[1]);
        else if (proj == evd::kXZ)
          pm4.SetPoint(hitIdx++, pos[2], pos[0]);
        else if (proj == evd::kYZ)
          pm4.SetPoint(hitIdx++, pos[2], pos[1]);
      }

      hitIdx = 0;

      TPolyMarker& pm5 = view->AddPolyMarker(skelEdgePosVec.size(), 3, kFullDotMedium, 1);
      for (const auto* spacePoint : skelEdgePosVec) {
        const double* pos = spacePoint->XYZ();

        if (proj == evd::kXY)
          pm5.SetPoint(hitIdx++, pos[0], pos[1]);
        else if (proj == evd::kXZ)
          pm5.SetPoint(hitIdx++, pos[2], pos[0]);
        else if (proj == evd::kYZ)
          pm5.SetPoint(hitIdx++, pos[2], pos[1]);
      }

      hitIdx = 0;

      TPolyMarker& pm6 = view->AddPolyMarker(seedPosVec.size(), 6, kFullDotMedium, 1);
      for (const auto* spacePoint : seedPosVec) {
        const double* pos = spacePoint->XYZ();

        if (proj == evd::kXY)
          pm6.SetPoint(hitIdx++, pos[0], pos[1]);
        else if (proj == evd::kXZ)
          pm6.SetPoint(hitIdx++, pos[2], pos[0]);
        else if (proj == evd::kYZ)
          pm6.SetPoint(hitIdx++, pos[2], pos[1]);
      }
    }

    // Look up the PCA info
    if (pcAxisAssnVec.isValid()) {
      std::vector<const recob::PCAxis*> pcaVec(pcAxisAssnVec.at(pfPart->Self()));

      if (!pcaVec.empty()) {
        // For each axis we are going to draw a solid line between two points
        int numPoints(2);
        int lineWidth[2] = {3, 1};
        int lineStyle[2] = {1, 13};
        int lineColor[2] = {colorIdx, 18};
        int markStyle[2] = {4, 4};
        int pcaIdx(0);

        // The order of axes in the returned association vector is arbitrary... the "first" axis is
        // better and we can divine that by looking at the axis id's (the best will have been made first)
        if (pcaVec.size() > 1 && pcaVec.front()->getID() > pcaVec.back()->getID())
          std::reverse(pcaVec.begin(), pcaVec.end());

        for (const auto& pca : pcaVec) {
          // We need the mean position
          const double* avePosition = pca->getAvePosition();

          // Let's draw a marker at the interesting points
          int pmrkIdx(0);
          TPolyMarker& pmrk = view->AddPolyMarker(7, lineColor[pcaIdx], markStyle[pcaIdx], 1);

          if (proj == evd::kXY)
            pmrk.SetPoint(pmrkIdx++, avePosition[0], avePosition[1]);
          else if (proj == evd::kXZ)
            pmrk.SetPoint(pmrkIdx++, avePosition[2], avePosition[0]);
          else if (proj == evd::kYZ)
            pmrk.SetPoint(pmrkIdx++, avePosition[2], avePosition[1]);

          // Loop over pca dimensions
          for (int dimIdx = 0; dimIdx < 3; dimIdx++) {
            // Oh please oh please give me an instance of a poly line...
            TPolyLine& pl =
              view->AddPolyLine(numPoints, lineColor[pcaIdx], lineWidth[pcaIdx], lineStyle[pcaIdx]);

            // We will use the eigen value to give the length of the line we're going to plot
            double eigenValue = pca->getEigenValues()[dimIdx];

            // Make sure a valid eigenvalue
            if (eigenValue > 0) {
              // Really want the root of the eigen value
              eigenValue = 3. * sqrt(eigenValue);

              // Recover the eigenvector
              const std::vector<double>& eigenVector = pca->getEigenVectors()[dimIdx];

              // Set the first point
              double xl = avePosition[0] - 0.5 * eigenValue * eigenVector[0];
              double yl = avePosition[1] - 0.5 * eigenValue * eigenVector[1];
              double zl = avePosition[2] - 0.5 * eigenValue * eigenVector[2];

              if (proj == evd::kXY) {
                pl.SetPoint(0, xl, yl);
                pmrk.SetPoint(pmrkIdx++, xl, yl);
              }
              else if (proj == evd::kXZ) {
                pl.SetPoint(0, zl, xl);
                pmrk.SetPoint(pmrkIdx++, zl, xl);
              }
              else if (proj == evd::kYZ) {
                pl.SetPoint(0, zl, yl);
                pmrk.SetPoint(pmrkIdx++, zl, yl);
              }

              // Set the second point
              double xu = avePosition[0] + 0.5 * eigenValue * eigenVector[0];
              double yu = avePosition[1] + 0.5 * eigenValue * eigenVector[1];
              double zu = avePosition[2] + 0.5 * eigenValue * eigenVector[2];

              if (proj == evd::kXY) {
                pl.SetPoint(1, xu, yu);
                pmrk.SetPoint(pmrkIdx++, xu, yu);
              }
              else if (proj == evd::kXZ) {
                pl.SetPoint(1, zu, xu);
                pmrk.SetPoint(pmrkIdx++, zu, xu);
              }
              else if (proj == evd::kYZ) {
                pl.SetPoint(1, zu, yu);
                pmrk.SetPoint(pmrkIdx++, zu, yu);
              }
            }
          }

          // By convention we will have drawn the "best" axis first
          if (recoOpt->fBestPCAAxisOnly) break;

          pcaIdx++;
        }
      }
    }

    // Now let's loop over daughters and call drawing routine for them
    if (pfPart->NumDaughters() > 0) {
      depth++;

      for (const auto& daughterIdx : pfPart->Daughters()) {
        DrawPFParticleOrtho(pfParticleVec.at(daughterIdx),
                            pfParticleVec,
                            spacePointAssnVec,
                            pcAxisAssnVec,
                            depth,
                            proj,
                            view);
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::ProngOrtho(const art::Event& evt,
                             evd::OrthoProj_t proj,
                             double msize,
                             evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    // Tracks.

    if (recoOpt->fDrawTracks != 0) {
      for (size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
        art::InputTag which = recoOpt->fTrackLabels[imod];
        art::View<recob::Track> track;
        this->GetTracks(evt, which, track);

        for (size_t t = 0; t < track.vals().size(); ++t) {
          const recob::Track* ptrack = track.vals().at(t);
          int color = ptrack->ID() & 65535;

          // Draw track using only embedded information.

          DrawTrackOrtho(*ptrack, color, proj, msize, view);
        }
      }
    }

    // Showers.

    if (recoOpt->fDrawShowers != 0) {
      for (size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
        art::InputTag which = recoOpt->fShowerLabels[imod];
        art::View<recob::Shower> shower;
        this->GetShowers(evt, which, shower);

        for (size_t s = 0; s < shower.vals().size(); ++s) {
          const recob::Shower* pshower = shower.vals().at(s);
          int color = pshower->ID();
          DrawShowerOrtho(*pshower, color, proj, msize, view);
        }
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawSpacePointOrtho(std::vector<art::Ptr<recob::SpacePoint>>& spts,
                                      int color,
                                      evd::OrthoProj_t proj,
                                      double msize,
                                      evdb::View2D* view,
                                      int mode)
  {
    // Get services.

    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    // Organize space points into separate collections according to the color
    // we want them to be.
    // If option If option fColorSpacePointsByChisq is false, this means
    // having a single collection with color inherited from the prong
    // (specified by the argument color).

    std::map<int, std::vector<art::Ptr<recob::SpacePoint>>> spmap; // Indexed by color.

    for (auto& pspt : spts) {

      // By default use event display palette.

      int spcolor = evd::kColor[color % evd::kNCOLS];
      if (mode == 1) { //shower hits
        spcolor = evd::kColor2[color % evd::kNCOLS];
      }
      // For rainbow effect, choose root colors in range [51,100].
      // We are using 100=best (red), 51=worst (blue).

      if (recoOpt->fColorSpacePointsByChisq) {
        spcolor = 100 - 2.5 * pspt->Chisq();
        if (spcolor < 51) spcolor = 51;
        if (spcolor > 100) spcolor = 100;
      }
      spmap[spcolor].push_back(pspt);
    }

    // Loop over colors.
    // Note that larger (=better) space points are plotted on
    // top for optimal visibility.

    for (auto icolor : spmap) {
      int spcolor = icolor.first;
      std::vector<art::Ptr<recob::SpacePoint>>& psps = icolor.second;

      // Make and fill a polymarker.

      TPolyMarker& pm = view->AddPolyMarker(psps.size(), spcolor, kFullCircle, msize);
      for (size_t s = 0; s < psps.size(); ++s) {
        const recob::SpacePoint& spt = *psps[s];
        const double* xyz = spt.XYZ();
        switch (proj) {
        case evd::kXY: pm.SetPoint(s, xyz[0], xyz[1]); break;
        case evd::kXZ: pm.SetPoint(s, xyz[2], xyz[0]); break;
        case evd::kYZ: pm.SetPoint(s, xyz[2], xyz[1]); break;
        default:
          throw cet::exception("RecoBaseDrawer")
            << __func__ << ": unknown projection #" << ((int)proj) << "\n";
        } // switch
      }
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawTrackOrtho(const recob::Track& track,
                                 int color,
                                 evd::OrthoProj_t proj,
                                 double msize,
                                 evdb::View2D* view)
  {
    // Get options.

    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;

    if (recoOpt->fDrawTrackSpacePoints) {

      // Use brute force to find the module label and index of this
      // track, so that we can find associated space points and draw
      // them.

      const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
      std::vector<art::Handle<std::vector<recob::Track>>> handles;
      evt->getManyByType(handles);
      for (auto ih : handles) {
        const art::Handle<std::vector<recob::Track>> handle = ih;
        if (handle.isValid()) {
          const std::string& which = handle.provenance()->moduleLabel();
          art::FindManyP<recob::SpacePoint> fmsp(handle, *evt, which);

          int n = handle->size();
          for (int i = 0; i < n; ++i) {
            art::Ptr<recob::Track> p(handle, i);
            if (&*p == &track) {
              std::vector<art::Ptr<recob::SpacePoint>> spts = fmsp.at(i);
              DrawSpacePointOrtho(spts, color, proj, msize, view);
            }
          }
        }
      }
    }
    if (recoOpt->fDrawTrackTrajectoryPoints) {

      // Draw trajectory points.

      int np = track.NumberTrajectoryPoints();
      int vp = track.CountValidPoints();

      // Make and fill a polymarker.

      TPolyMarker& pm =
        view->AddPolyMarker(vp, evd::kColor[color % evd::kNCOLS], kFullCircle, msize);
      TPolyLine& pl = view->AddPolyLine(vp, evd::kColor[color % evd::kNCOLS], 2, 0);
      for (int p = 0; p < np; ++p) {
        if (track.HasValidPoint(p) == 0) continue;
        const auto& pos = track.LocationAtPoint(p);
        switch (proj) {
        case evd::kXY:
          pm.SetPoint(p, pos.X(), pos.Y());
          pl.SetPoint(p, pos.X(), pos.Y());
          break;
        case evd::kXZ:
          pm.SetPoint(p, pos.Z(), pos.X());
          pl.SetPoint(p, pos.Z(), pos.X());
          break;
        case evd::kYZ:
          pm.SetPoint(p, pos.Z(), pos.Y());
          pl.SetPoint(p, pos.Z(), pos.Y());
          break;
        default:
          throw cet::exception("RecoBaseDrawer")
            << __func__ << ": unknown projection #" << ((int)proj) << "\n";
        } // switch
      }   // p
      // BB: draw the track ID at the end of the track
      if (recoOpt->fDrawTracks > 1) {
        int tid =
          track.ID() &
          65535; //this is a hack for PMA track id which uses the 16th bit to identify shower-like track.
        std::string s = std::to_string(tid);
        char const* txt = s.c_str();
        double x = track.End().X();
        double y = track.End().Y();
        double z = track.End().Z();
        if (proj == evd::kXY) {
          TText& trkID = view->AddText(x, y, txt);
          trkID.SetTextColor(evd::kColor[tid % evd::kNCOLS]);
          trkID.SetTextSize(0.03);
        }
        else if (proj == evd::kXZ) {
          TText& trkID = view->AddText(z, x, txt);
          trkID.SetTextColor(evd::kColor[tid % evd::kNCOLS]);
          trkID.SetTextSize(0.03);
        }
        else if (proj == evd::kYZ) {
          TText& trkID = view->AddText(z, y, txt);
          trkID.SetTextColor(evd::kColor[tid % evd::kNCOLS]);
          trkID.SetTextSize(0.03);
        } // proj
      }   // recoOpt->fDrawTracks > 1
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::DrawShowerOrtho(const recob::Shower& shower,
                                  int color,
                                  evd::OrthoProj_t proj,
                                  double msize,
                                  evdb::View2D* view)
  {
    // Use brute force to find the module label and index of this
    // shower, so that we can find associated space points and draw
    // them.

    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    std::vector<art::Handle<std::vector<recob::Shower>>> handles;
    evt->getManyByType(handles);
    for (auto ih : handles) {
      const art::Handle<std::vector<recob::Shower>> handle = ih;
      if (handle.isValid()) {
        const std::string& which = handle.provenance()->moduleLabel();
        art::FindManyP<recob::SpacePoint> fmsp(handle, *evt, which);
        if (!fmsp.isValid()) continue;
        int n = handle->size();
        for (int i = 0; i < n; ++i) {
          art::Ptr<recob::Shower> p(handle, i);
          if (&*p == &shower) {
            switch (proj) {
            case evd::kXY:
              view->AddMarker(p->ShowerStart().X(),
                              p->ShowerStart().Y(),
                              evd::kColor2[color % evd::kNCOLS],
                              5,
                              2.0);
              break;
            case evd::kXZ:
              view->AddMarker(p->ShowerStart().Z(),
                              p->ShowerStart().X(),
                              evd::kColor2[color % evd::kNCOLS],
                              5,
                              2.0);
              break;
            case evd::kYZ:
              view->AddMarker(p->ShowerStart().Z(),
                              p->ShowerStart().Y(),
                              evd::kColor2[color % evd::kNCOLS],
                              5,
                              2.0);
              break;
            default:
              throw cet::exception("RecoBaseDrawer")
                << __func__ << ": unknown projection #" << ((int)proj) << "\n";
            } // switch

            if (fmsp.isValid()) {
              std::vector<art::Ptr<recob::SpacePoint>> spts = fmsp.at(i);
              DrawSpacePointOrtho(spts, color, proj, msize, view, 1);
            }
          }
        }
      }
    }

    return;
  }

  //......................................................................
  int
  RecoBaseDrawer::GetWires(const art::Event& evt,
                           const art::InputTag& which,
                           art::PtrVector<recob::Wire>& wires)
  {
    wires.clear();

    art::Handle<std::vector<recob::Wire>> wcol;
    art::PtrVector<recob::Wire> temp;

    try {
      evt.getByLabel(which, wcol);

      for (unsigned int i = 0; i < wcol->size(); ++i) {
        art::Ptr<recob::Wire> w(wcol, i);
        temp.push_back(w);
      }
      temp.swap(wires);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetWires", e);
    }

    return wires.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetHits(const art::Event& evt,
                          const art::InputTag& which,
                          std::vector<const recob::Hit*>& hits,
                          unsigned int plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    hits.clear();

    std::vector<const recob::Hit*> temp;

    try {
      evt.getView(which, temp);
      for (const auto& hit : temp) {
        // Note that the WireID in the hit object is useless for those detectors where a channel can correspond to
        // more than one plane/wire. So our plan is to recover the list of wire IDs from the channel number and
        // loop over those (if there are any)
        const std::vector<geo::WireID>& wireIDs = geo->ChannelToWire(hit->Channel());

        // Loop to find match
        for (const auto& wireID : wireIDs) {
          if (wireID.Plane == plane && wireID.TPC == rawOpt->fTPC &&
              wireID.Cryostat == rawOpt->fCryostat)
            hits.push_back(hit);
        }
      }
    }
    catch (cet::exception& e) {
      writeErrMsg("GetHits", e);
    }

    return hits.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetSlices(const art::Event& evt,
                            const art::InputTag& which,
                            art::PtrVector<recob::Slice>& slices)
  {
    slices.clear();
    art::PtrVector<recob::Slice> temp;

    art::Handle<std::vector<recob::Slice>> slcCol;

    try {
      evt.getByLabel(which, slcCol);
      temp.reserve(slcCol->size());
      for (unsigned int i = 0; i < slcCol->size(); ++i) {
        art::Ptr<recob::Slice> slc(slcCol, i);
        temp.push_back(slc);
      }
      temp.swap(slices);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetSlices", e);
    }

    return slices.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetClusters(const art::Event& evt,
                              const art::InputTag& which,
                              art::PtrVector<recob::Cluster>& clust)
  {
    clust.clear();
    art::PtrVector<recob::Cluster> temp;

    art::Handle<std::vector<recob::Cluster>> clcol;

    try {
      evt.getByLabel(which, clcol);
      temp.reserve(clcol->size());
      for (unsigned int i = 0; i < clcol->size(); ++i) {
        art::Ptr<recob::Cluster> cl(clcol, i);
        temp.push_back(cl);
      }
      temp.swap(clust);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetClusters", e);
    }

    return clust.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetPFParticles(const art::Event& evt,
                                 const art::InputTag& which,
                                 art::PtrVector<recob::PFParticle>& clust)
  {
    clust.clear();
    art::PtrVector<recob::PFParticle> temp;

    art::Handle<std::vector<recob::PFParticle>> clcol;

    try {
      evt.getByLabel(which, clcol);
      for (unsigned int i = 0; i < clcol->size(); ++i) {
        art::Ptr<recob::PFParticle> cl(clcol, i);
        temp.push_back(cl);
      }
      temp.swap(clust);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetPFParticles", e);
    }

    return clust.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetEndPoint2D(const art::Event& evt,
                                const art::InputTag& which,
                                art::PtrVector<recob::EndPoint2D>& ep2d)
  {
    ep2d.clear();
    art::PtrVector<recob::EndPoint2D> temp;

    art::Handle<std::vector<recob::EndPoint2D>> epcol;

    try {
      evt.getByLabel(which, epcol);
      for (unsigned int i = 0; i < epcol->size(); ++i) {
        art::Ptr<recob::EndPoint2D> ep(epcol, i);
        temp.push_back(ep);
      }
      temp.swap(ep2d);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetEndPoint2D", e);
    }

    return ep2d.size();
  }

  //......................................................................

  int
  RecoBaseDrawer::GetOpFlashes(const art::Event& evt,
                               const art::InputTag& which,
                               art::PtrVector<recob::OpFlash>& opflashes)
  {
    opflashes.clear();
    art::PtrVector<recob::OpFlash> temp;

    art::Handle<std::vector<recob::OpFlash>> opflashcol;

    try {
      evt.getByLabel(which, opflashcol);
      for (unsigned int i = 0; i < opflashcol->size(); ++i) {
        art::Ptr<recob::OpFlash> opf(opflashcol, i);
        temp.push_back(opf);
      }
      temp.swap(opflashes);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetOpFlashes", e);
    }

    return opflashes.size();
  }

  //......................................................................

  int
  RecoBaseDrawer::GetSeeds(const art::Event& evt,
                           const art::InputTag& which,
                           art::PtrVector<recob::Seed>& seeds)
  {
    seeds.clear();
    art::PtrVector<recob::Seed> temp;

    art::Handle<std::vector<recob::Seed>> seedcol;

    try {
      evt.getByLabel(which, seedcol);
      for (unsigned int i = 0; i < seedcol->size(); ++i) {
        art::Ptr<recob::Seed> sd(seedcol, i);
        temp.push_back(sd);
      }
      temp.swap(seeds);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetSeeds", e);
    }

    return seeds.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetSpacePoints(const art::Event& evt,
                                 const art::InputTag& which,
                                 std::vector<art::Ptr<recob::SpacePoint>>& spts)
  {
    spts.clear();
    art::Handle<std::vector<recob::SpacePoint>> spcol;
    if (evt.getByLabel(which, spcol)) art::fill_ptr_vector(spts, spcol);

    return spts.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetEdges(const art::Event& evt,
                           const art::InputTag& which,
                           std::vector<art::Ptr<recob::Edge>>& edges)
  {
    edges.clear();

    art::Handle<std::vector<recob::Edge>> edgeCol;

    evt.getByLabel(which, edgeCol);

    for (unsigned int i = 0; i < edgeCol->size(); ++i)
      edges.emplace_back(edgeCol, i);

    return edges.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetTracks(const art::Event& evt,
                            const art::InputTag& which,
                            art::View<recob::Track>& track)
  {
    try {
      evt.getView(which, track);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetTracks", e);
    }

    return track.vals().size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetShowers(const art::Event& evt,
                             const art::InputTag& which,
                             art::View<recob::Shower>& shower)
  {
    try {
      evt.getView(which, shower);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetShowers", e);
    }

    return shower.vals().size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetVertices(const art::Event& evt,
                              const art::InputTag& which,
                              art::PtrVector<recob::Vertex>& vertex)
  {
    vertex.clear();
    art::PtrVector<recob::Vertex> temp;

    art::Handle<std::vector<recob::Vertex>> vcol;

    try {
      evt.getByLabel(which, vcol);
      for (size_t i = 0; i < vcol->size(); ++i) {
        art::Ptr<recob::Vertex> v(vcol, i);
        temp.push_back(v);
      }
      temp.swap(vertex);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetVertices", e);
    }

    return vertex.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::GetEvents(const art::Event& evt,
                            const art::InputTag& which,
                            art::PtrVector<recob::Event>& event)
  {
    event.clear();
    art::PtrVector<recob::Event> temp;

    art::Handle<std::vector<recob::Event>> ecol;

    try {
      evt.getByLabel(which, ecol);
      for (size_t i = 0; i < ecol->size(); ++i) {
        art::Ptr<recob::Event> e(ecol, i);
        temp.push_back(e);
      }
      temp.swap(event);
    }
    catch (cet::exception& e) {
      writeErrMsg("GetEvents", e);
    }

    return event.size();
  }

  //......................................................................
  int
  RecoBaseDrawer::CountHits(const art::Event& evt,
                            const art::InputTag& which,
                            unsigned int cryostat,
                            unsigned int tpc,
                            unsigned int plane)
  {
    std::vector<const recob::Hit*> temp;
    int NumberOfHitsBeforeThisPlane = 0;
    evt.getView(
      which,
      temp); //temp.size() = total number of hits for this event (number of all hits in all Cryostats, TPC's, planes and wires)
    for (size_t t = 0; t < temp.size(); ++t) {
      if (temp[t]->WireID().Cryostat == cryostat && temp[t]->WireID().TPC == tpc &&
          temp[t]->WireID().Plane == plane)
        break;
      NumberOfHitsBeforeThisPlane++;
    }
    return NumberOfHitsBeforeThisPlane;
  }

  //......................................................................
  void
  RecoBaseDrawer::FillTQHisto(const art::Event& evt,
                              unsigned int plane,
                              unsigned int wire,
                              TH1F* histo)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    float minSig(std::numeric_limits<float>::max());
    float maxSig(std::numeric_limits<float>::lowest());
    bool setLimits(false);

    // Check if we're supposed to draw raw hits at all
    if (rawOpt->fDrawRawDataOrCalibWires == 0) return;

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (size_t i = 0; i < wires.size(); ++i) {

        std::vector<geo::WireID> wireids = geo->ChannelToWire(wires[i]->Channel());

        bool goodWID = false;
        for (auto const& wid : wireids) {
          // check for correct plane, wire and tpc
          if (wid.Plane == plane && wid.Wire == wire && wid.TPC == rawOpt->fTPC &&
              wid.Cryostat == rawOpt->fCryostat)
            goodWID = true;
        }
        if (!goodWID) continue;

        std::vector<float> wirSig = wires[i]->Signal();
        for (unsigned int ii = 0; ii < wirSig.size(); ++ii) {
          //                histo->SetLineColor(imod+4);
          //                histo->Fill(1.*ii, wirSig[ii]);
          minSig = std::min(minSig, wirSig[ii]);
          maxSig = std::max(maxSig, wirSig[ii]);
        }

        setLimits = true;
      } //end loop over wires
    }   //end loop over wire modules

    if (setLimits) {
      histo->SetMaximum(1.2 * maxSig);
      histo->SetMinimum(1.2 * minSig);
    }

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::FillQHisto(const art::Event& evt, unsigned int plane, TH1F* histo)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    // Check if we're supposed to draw raw hits at all
    if (rawOpt->fDrawRawDataOrCalibWires == 0) return;

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (unsigned int i = 0; i < wires.size(); ++i) {

        std::vector<geo::WireID> wireids = geo->ChannelToWire(wires[i]->Channel());

        bool goodWID = false;
        for (auto const& wid : wireids) {
          // check for correct plane, wire and tpc
          if (wid.Plane == plane && wid.TPC == rawOpt->fTPC && wid.Cryostat == rawOpt->fCryostat)
            goodWID = true;
        }
        if (!goodWID) continue;
        std::vector<float> wirSig = wires[i]->Signal();
        for (unsigned int ii = 0; ii < wirSig.size(); ++ii)
          histo->Fill(wirSig[ii]);
        /*
        for(size_t s = 0; s < wires[i]->NSignal(); ++s)
          histo->Fill(wires[i]->Signal()[s]);
*/

      } //end loop over raw hits
    }   //end loop over Wire modules

    return;
  }

  //......................................................................
  void
  RecoBaseDrawer::FillTQHistoDP(const art::Event& evt,
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
                                std::vector<int>& hLocalHitIndex)
  {
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions const> recoOpt;
    art::ServiceHandle<geo::Geometry const> geo;

    // Check if we're supposed to draw raw hits at all
    if (rawOpt->fDrawRawDataOrCalibWires == 0) return;

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (size_t i = 0; i < wires.size(); ++i) {

        std::vector<geo::WireID> wireids = geo->ChannelToWire(wires[i]->Channel());

        bool goodWID = false;
        for (auto const& wid : wireids) {
          if (wid.Plane == plane && wid.Wire == wire && wid.TPC == rawOpt->fTPC &&
              wid.Cryostat == rawOpt->fCryostat)
            goodWID = true;
        }

        if (!goodWID) continue;

        std::vector<float> wirSig = wires[i]->Signal();
        for (unsigned int ii = 0; ii < wirSig.size(); ++ii)
          histo->Fill(1. * ii, wirSig[ii]);
        break;
      } //end loop over wires
    }   //end loop over wire modules

    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      art::InputTag const which = recoOpt->fHitLabels[imod];

      std::vector<const recob::Hit*> hits;
      this->GetHits(evt, which, hits, plane);

      auto hitResults = anab::FVectorReader<recob::Hit, 4>::create(evt, "dprawhit");
      const auto& fitParams = hitResults->vectors();

      int FitParamsOffset = CountHits(evt, which, rawOpt->fCryostat, rawOpt->fTPC, plane);

      for (size_t i = 0; i < hits.size(); ++i) {
        // check for correct wire. Plane, cryostat and tpc were checked in GetHits
        if (hits[i]->WireID().Wire != wire) continue;

        hpeaktimes.push_back(fitParams[FitParamsOffset + i][0]);
        htau1.push_back(fitParams[FitParamsOffset + i][1]);
        htau2.push_back(fitParams[FitParamsOffset + i][2]);
        hitamplitudes.push_back(fitParams[FitParamsOffset + i][3]);
        hstartT.push_back(hits[i]->StartTick());
        hendT.push_back(hits[i]->EndTick());
        hNMultiHit.push_back(hits[i]->Multiplicity());
        hLocalHitIndex.push_back(hits[i]->LocalIndex());
      } //end loop over reco hits
    }   //end loop over HitFinding modules

    return;
  }

  //......................................................................
  //double RecoBaseDrawer::EvalExpoFit(double x,
  //				   double tau1,
  //				   double tau2,
  //				   double amplitude,
  //				   double peaktime)
  //{
  //return (amplitude * exp(0.4*(x-peaktime)/tau1) / ( 1 + exp(0.4*(x-peaktime)/tau2) ) );
  //}

  //......................................................................
  //double RecoBaseDrawer::EvalMultiExpoFit(double x,
  //					int HitNumber,
  //					int NHits,
  //                                    std::vector<double> tau1,
  //                                    std::vector<double> tau2,
  //                                    std::vector<double> amplitude,
  //                                    std::vector<double> peaktime)
  //{
  //    double x_sum = 0.;
  //
  //    for(int i = HitNumber; i < HitNumber+NHits; i++)
  //    {
  //    x_sum += (amplitude[i] * exp(0.4*(x-peaktime[i])/tau1[i]) / ( 1 + exp(0.4*(x-peaktime[i])/tau2[i]) ) );
  //    }
  //
  //return x_sum;
  //}

} // namespace evd
////////////////////////////////////////////////////////////////////////
