////////////////////////////////////////////////////////////////////////
/// \file   ICARUSDrawer.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"

#include "larcore/Geometry/Geometry.h"
#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "nuevdb/EventDisplayBase/View3D.h"

#include "TPolyLine3D.h"

namespace evd_tool {

  class ICARUSDrawer : IExperimentDrawer {
  public:
    explicit ICARUSDrawer(const fhicl::ParameterSet& pset);

    void DetOutline3D(evdb::View3D* view) override;

    ~ICARUSDrawer() {}

  private:
    void configure(const fhicl::ParameterSet& pset);
    void DrawRectangularBox(evdb::View3D* view,
                            double* coordsLo,
                            double* coordsHi,
                            int color = kGray,
                            int width = 1,
                            int style = 1);
    void DrawGrids(evdb::View3D* view,
                   double* coordsLo,
                   double* coordsHi,
                   bool verticalGrid,
                   int color = kGray,
                   int width = 1,
                   int style = 1);
    void DrawAxes(evdb::View3D* view,
                  double* coordsLo,
                  double* coordsHi,
                  int color = kGray,
                  int width = 1,
                  int style = 1);
    void DrawBadChannels(evdb::View3D* view, double* coords, int color, int width, int style);

    // Member variables from the fhicl file
    bool fDrawGrid;        ///< true to draw backing grid
    bool fDrawAxes;        ///< true to draw coordinate axes
    bool fDrawBadChannels; ///< true to draw bad channels
  };

  //----------------------------------------------------------------------
  // Constructor.
  ICARUSDrawer::ICARUSDrawer(const fhicl::ParameterSet& pset) { configure(pset); }

  void ICARUSDrawer::configure(const fhicl::ParameterSet& pset)
  {
    // Start by recovering the parameters
    fDrawGrid = pset.get<bool>("DrawGrid", true);
    fDrawAxes = pset.get<bool>("DrawAxes", true);
    fDrawBadChannels = pset.get<bool>("DrawBadChannels", true);

    return;
  }

  //......................................................................
  void ICARUSDrawer::DetOutline3D(evdb::View3D* view)
  {
    art::ServiceHandle<geo::Geometry const> geo;

    bool axesNotDrawn(true);

    double xl, xu, yl, yu, zl, zu;

    geo->WorldBox(&xl, &xu, &yl, &yu, &zl, &zu);

    std::cout << "--- building ICARUS 3D display, low coord: " << xl << ", " << yl << ", " << zl
              << ", hi coord: " << xu << ", " << yu << ", " << zu << std::endl;

    // Loop over the number of cryostats
    for (geo::cryostat_iterator cryoItr = geo->begin_cryostat(); cryoItr != geo->end_cryostat();
         cryoItr++) {
      const geo::CryostatGeo& cryoGeo = *cryoItr;

      double cryoCoordsLo[] = {cryoGeo.MinX(), cryoGeo.MinY(), cryoGeo.MinZ()};
      double cryoCoordsHi[] = {cryoGeo.MaxX(), cryoGeo.MaxY(), cryoGeo.MaxZ()};

      std::cout << "    - cryostat: " << cryoGeo.ID() << ", low coord: " << cryoCoordsLo[0] << ", "
                << cryoCoordsLo[1] << ", " << cryoCoordsLo[2] << ", hi coord: " << cryoCoordsHi[0]
                << ", " << cryoCoordsHi[1] << ", " << cryoCoordsHi[2] << std::endl;

      DrawRectangularBox(view, cryoCoordsLo, cryoCoordsHi, kWhite, 2, 1);

      if (fDrawAxes && axesNotDrawn) {
        DrawAxes(view, cryoCoordsLo, cryoCoordsHi, kBlue, 1, 1);
        axesNotDrawn = true;
      }

      // Now draw the TPC's associated to this cryostat
      for (size_t tpcIdx = 0; tpcIdx < cryoGeo.NTPC(); tpcIdx++) {
        const geo::TPCGeo& tpcGeo = cryoGeo.TPC(tpcIdx);

        // Find the center of the current TPC
        TVector3 tpcCenter = tpcGeo.GetCenter();

        // Now draw the standard volume
        double coordsLo[] = {tpcCenter.X() - tpcGeo.HalfWidth(),
                             tpcCenter.Y() - tpcGeo.HalfHeight(),
                             tpcCenter.Z() - 0.5 * tpcGeo.Length()};
        double coordsHi[] = {tpcCenter.X() + tpcGeo.HalfWidth(),
                             tpcCenter.Y() + tpcGeo.HalfHeight(),
                             tpcCenter.Z() + 0.5 * tpcGeo.Length()};

        std::cout << "     - TPC: " << tpcGeo.ID() << ", low coord: " << coordsLo[0] << ", "
                  << coordsLo[1] << ", " << coordsLo[2] << ", hi coord: " << coordsHi[0] << ", "
                  << coordsHi[1] << ", " << coordsHi[2] << std::endl;

        DrawRectangularBox(view, coordsLo, coordsHi, kRed, 2, 1);

        // It could be that we don't want to see the grids
        if (fDrawGrid) DrawGrids(view, coordsLo, coordsHi, tpcIdx > 0, kGray + 2, 1, 1);

        if (fDrawBadChannels) DrawBadChannels(view, coordsHi, kGray, 1, 1);
      }
    }

    return;
  }

  void ICARUSDrawer::DrawRectangularBox(evdb::View3D* view,
                                        double* coordsLo,
                                        double* coordsHi,
                                        int color,
                                        int width,
                                        int style)
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

  void ICARUSDrawer::DrawGrids(evdb::View3D* view,
                               double* coordsLo,
                               double* coordsHi,
                               bool verticalGrid,
                               int color,
                               int width,
                               int style)
  {
    double z = coordsLo[2];
    // Grid running along x and y at constant z
    while (1) {
      TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
      gridt.SetPoint(0, coordsLo[0], coordsLo[1], z);
      gridt.SetPoint(1, coordsHi[0], coordsLo[1], z);

      if (verticalGrid) {
        TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
        grids.SetPoint(0, coordsHi[0], coordsLo[1], z);
        grids.SetPoint(1, coordsHi[0], coordsHi[1], z);
      }

      z += 10.0;
      if (z > coordsHi[2]) break;
    }

    // Grid running along z at constant x
    double x = coordsLo[0];
    while (1) {
      TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
      gridt.SetPoint(0, x, coordsLo[1], coordsLo[2]);
      gridt.SetPoint(1, x, coordsLo[1], coordsHi[2]);
      x += 10.0;
      if (x > coordsHi[0]) break;
    }

    // Grid running along z at constant y
    if (verticalGrid) {
      double y = coordsLo[1];
      while (1) {
        TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
        grids.SetPoint(0, coordsHi[0], y, coordsLo[2]);
        grids.SetPoint(1, coordsHi[0], y, coordsHi[2]);
        y += 10.0;
        if (y > coordsHi[1]) break;
      }
    }

    return;
  }

  void ICARUSDrawer::DrawAxes(evdb::View3D* view,
                              double* coordsLo,
                              double* coordsHi,
                              int color,
                              int width,
                              int style)
  {

    // Indicate coordinate system
    double x0 = -0.20;               // Center location of the key
    double y0 = 1.10 * coordsLo[1];  // Center location of the key
    double z0 = -0.10 * coordsHi[2]; // Center location of the key
    double sz = 0.20 * coordsHi[2];  // Scale size of the key in z direction

    TPolyLine3D& xaxis = view->AddPolyLine3D(2, color, style, width);
    TPolyLine3D& yaxis = view->AddPolyLine3D(2, color, style, width);
    TPolyLine3D& zaxis = view->AddPolyLine3D(2, color, style, width);
    xaxis.SetPoint(0, x0, y0, z0);
    xaxis.SetPoint(1, sz + x0, y0, z0);

    yaxis.SetPoint(0, x0, y0, z0);
    yaxis.SetPoint(1, x0, y0 + sz, z0);

    zaxis.SetPoint(0, x0, y0, z0);
    zaxis.SetPoint(1, x0, y0, z0 + sz);

    TPolyLine3D& xpoint = view->AddPolyLine3D(3, color, style, width);
    TPolyLine3D& ypoint = view->AddPolyLine3D(3, color, style, width);
    TPolyLine3D& zpoint = view->AddPolyLine3D(3, color, style, width);

    xpoint.SetPoint(0, 0.95 * sz + x0, y0, z0 - 0.05 * sz);
    xpoint.SetPoint(1, 1.00 * sz + x0, y0, z0);
    xpoint.SetPoint(2, 0.95 * sz + x0, y0, z0 + 0.05 * sz);

    ypoint.SetPoint(0, x0, 0.95 * sz + y0, z0 - 0.05 * sz);
    ypoint.SetPoint(1, x0, 1.00 * sz + y0, z0);
    ypoint.SetPoint(2, x0, 0.95 * sz + y0, z0 + 0.05 * sz);

    zpoint.SetPoint(0, x0 - 0.05 * sz, y0, 0.95 * sz + z0);
    zpoint.SetPoint(1, x0 + 0.00 * sz, y0, 1.00 * sz + z0);
    zpoint.SetPoint(2, x0 + 0.05 * sz, y0, 0.95 * sz + z0);

    TPolyLine3D& zleg = view->AddPolyLine3D(4, color, style, width);
    zleg.SetPoint(0, x0 - 0.05 * sz, y0 + 0.05 * sz, z0 + 1.05 * sz);
    zleg.SetPoint(1, x0 + 0.05 * sz, y0 + 0.05 * sz, z0 + 1.05 * sz);
    zleg.SetPoint(2, x0 - 0.05 * sz, y0 - 0.05 * sz, z0 + 1.05 * sz);
    zleg.SetPoint(3, x0 + 0.05 * sz, y0 - 0.05 * sz, z0 + 1.05 * sz);

    TPolyLine3D& yleg = view->AddPolyLine3D(5, color, style, width);
    yleg.SetPoint(0, x0 - 0.05 * sz, y0 + 1.15 * sz, z0);
    yleg.SetPoint(1, x0 + 0.00 * sz, y0 + 1.10 * sz, z0);
    yleg.SetPoint(2, x0 + 0.00 * sz, y0 + 1.05 * sz, z0);
    yleg.SetPoint(3, x0 + 0.00 * sz, y0 + 1.10 * sz, z0);
    yleg.SetPoint(4, x0 + 0.05 * sz, y0 + 1.15 * sz, z0);

    TPolyLine3D& xleg = view->AddPolyLine3D(7, color, style, width);
    xleg.SetPoint(0, x0 + 1.05 * sz, y0 + 0.05 * sz, z0 - 0.05 * sz);
    xleg.SetPoint(1, x0 + 1.05 * sz, y0 + 0.00 * sz, z0 - 0.00 * sz);
    xleg.SetPoint(2, x0 + 1.05 * sz, y0 + 0.05 * sz, z0 + 0.05 * sz);
    xleg.SetPoint(3, x0 + 1.05 * sz, y0 + 0.00 * sz, z0 - 0.00 * sz);
    xleg.SetPoint(4, x0 + 1.05 * sz, y0 - 0.05 * sz, z0 - 0.05 * sz);
    xleg.SetPoint(5, x0 + 1.05 * sz, y0 + 0.00 * sz, z0 - 0.00 * sz);
    xleg.SetPoint(6, x0 + 1.05 * sz, y0 - 0.05 * sz, z0 + 0.05 * sz);

    return;
  }

  void ICARUSDrawer::DrawBadChannels(evdb::View3D* view,
                                     double* coords,
                                     int color,
                                     int width,
                                     int style)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    art::ServiceHandle<evd::RawDrawingOptions const> rawOpt;

    lariov::ChannelStatusProvider const& channelStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    // We want to translate the wire position to the opposite side of the TPC...
    for (size_t viewNo = 0; viewNo < geo->Nviews(); viewNo++) {
      for (size_t wireNo = 0; wireNo < geo->Nwires(viewNo); wireNo++) {
        geo::WireID wireID = geo::WireID(rawOpt->fCryostat, rawOpt->fTPC, viewNo, wireNo);

        raw::ChannelID_t channel = geo->PlaneWireToChannel(wireID);

        if (channelStatus.IsBad(channel)) {
          const geo::WireGeo* wireGeo = geo->WirePtr(wireID);

          double wireStart[3];
          double wireEnd[3];

          wireGeo->GetStart(wireStart);
          wireGeo->GetEnd(wireEnd);

          TPolyLine3D& pl = view->AddPolyLine3D(2, color, style, width);
          pl.SetPoint(0, coords[0] - 0.5, wireStart[1], wireStart[2]);
          pl.SetPoint(1, coords[0] - 0.5, wireEnd[1], wireEnd[2]);
        }
      }
    }

    return;
  }

  DEFINE_ART_CLASS_TOOL(ICARUSDrawer)
}
