////////////////////////////////////////////////////////////////////////
/// \file   StandardDrawer.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"

#include "art/Utilities/ToolMacros.h"

#include "fhiclcpp/ParameterSet.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "nuevdb/EventDisplayBase/View3D.h"

#include "TPolyLine3D.h"

#include <algorithm> // std::min()
#include <array>
#include <cmath> // std::abs()

namespace evd_tool {

  class StandardDrawer : public IExperimentDrawer {
  public:
    explicit StandardDrawer(const fhicl::ParameterSet& pset);

    virtual void DetOutline3D(evdb::View3D* view) override;

  protected:
    /// Draw the outline of an object bounded by a box.
    void DrawBoxBoundedGeoOutline(evdb::View3D* view,
                                  geo::BoxBoundedGeo const& bb,
                                  Color_t color,
                                  Width_t width,
                                  Style_t style) const;

    /// Draw the outline of the TPC volume.
    void DrawTPCoutline(evdb::View3D* view,
                        geo::TPCGeo const& TPC,
                        Color_t color,
                        Width_t width,
                        Style_t style) const
    {
      DrawBoxBoundedGeoOutline(view, TPC, color, width, style);
    }

    /// Draw the outline of the TPC active volume.
    void DrawActiveTPCoutline(evdb::View3D* view,
                              geo::TPCGeo const& TPC,
                              Color_t color,
                              Width_t width,
                              Style_t style) const;

    void DrawRectangularBox(evdb::View3D* view,
                            double const* coordsLo,
                            double const* coordsHi,
                            int color = kGray,
                            int width = 1,
                            int style = 1) const;
    void DrawGrids(evdb::View3D* view,
                   double const* coordsLo,
                   double const* coordsHi,
                   int color = kGray,
                   int width = 1,
                   int style = 1) const;
    void DrawAxes(evdb::View3D* view,
                  double const* coordsLo,
                  double const* coordsHi,
                  int color = kGray,
                  int width = 1,
                  int style = 1) const;

  private:
    void configure(const fhicl::ParameterSet& pset);
    // Member variables from the fhicl file
    bool fDrawGrid;   ///< true to draw backing grid
    bool fDrawAxes;   ///< true to draw coordinate axes
    bool fDrawActive; ///< true to outline TPC sensitive volumes
  };

  //----------------------------------------------------------------------
  // Constructor.
  StandardDrawer::StandardDrawer(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }

  void StandardDrawer::configure(const fhicl::ParameterSet& pset)
  {
    // Start by recovering the parameters
    fDrawGrid = pset.get<bool>("DrawGrid", true);
    fDrawAxes = pset.get<bool>("DrawAxes", true);
    fDrawActive = pset.get<bool>("DrawActive", true);

    return;
  }

  //......................................................................
  void StandardDrawer::DetOutline3D(evdb::View3D* view)
  {
    auto const& geom = *(lar::providerFrom<geo::Geometry>());

    // we compute the total volume of the detector, to be used for the axes;
    // we do include the origin by choice
    geo::BoxBoundedGeo detector({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});

    // Draw a box for each cryostat, and, within it, for each TPC;
    // the outlined volumes are the ones from the geometry boxes
    for (geo::CryostatGeo const& cryo : geom.Iterate<geo::CryostatGeo>()) {

      // include this cryostat in the detector volume
      detector.ExtendToInclude(cryo);

      // draw the cryostat box
      DrawBoxBoundedGeoOutline(view, cryo.Boundaries(), kRed + 2, 1, kSolid);

      // draw all TPC boxes
      for (geo::TPCGeo const& TPC : cryo.IterateTPCs()) {

        DrawTPCoutline(view, TPC, kRed, 2, kSolid);

        // BUG the double brace syntax is required to work around clang bug 21629
        // optionally draw the grid
        if (fDrawGrid) {
          std::array<double, 3U> const tpcLow{{TPC.MinX(), TPC.MinY(), TPC.MinZ()}},
            tpcHigh{{TPC.MaxX(), TPC.MaxY(), TPC.MaxZ()}};
          DrawGrids(view, tpcLow.data(), tpcHigh.data(), kGray + 2, 1, kSolid);
        }

        // optionally draw the active volume
        if (fDrawActive) DrawActiveTPCoutline(view, TPC, kCyan + 2, 1, kDotted);

      } // for TPCs in cryostat

    } // for cryostats

    // draw axes if requested
    if (fDrawAxes) {
      // BUG the double brace syntax is required to work around clang bug 21629
      std::array<double, 3U> const detLow = {{detector.MinX(), detector.MinY(), detector.MinZ()}},
                                   detHigh = {{detector.MaxX(), detector.MaxY(), detector.MaxZ()}};
      DrawAxes(view, detLow.data(), detHigh.data(), kBlue, 1, kSolid);
    } // if draw axes
  }

  void StandardDrawer::DrawBoxBoundedGeoOutline(evdb::View3D* view,
                                                geo::BoxBoundedGeo const& bb,
                                                Color_t color,
                                                Width_t width,
                                                Style_t style) const
  {
    // BUG the double brace syntax is required to work around clang bug 21629
    std::array<double, 3U> const low{{bb.MinX(), bb.MinY(), bb.MinZ()}},
      high{{bb.MaxX(), bb.MaxY(), bb.MaxZ()}};
    ;
    DrawRectangularBox(view, low.data(), high.data(), color, width, style);
  } // StandardDrawer::DrawBoxBoundedGeoOutline()

  void StandardDrawer::DrawActiveTPCoutline(evdb::View3D* view,
                                            geo::TPCGeo const& TPC,
                                            Color_t color,
                                            Width_t width,
                                            Style_t style) const
  {
    auto const& activeCenter = TPC.GetActiveVolumeCenter();
    DrawBoxBoundedGeoOutline(view,
                             {{activeCenter.X() - TPC.ActiveHalfWidth(),
                               activeCenter.Y() - TPC.ActiveHalfHeight(),
                               activeCenter.Z() - TPC.ActiveHalfLength()},
                              {activeCenter.X() + TPC.ActiveHalfWidth(),
                               activeCenter.Y() + TPC.ActiveHalfHeight(),
                               activeCenter.Z() + TPC.ActiveHalfLength()}},
                             color,
                             width,
                             style);
  }

  void StandardDrawer::DrawRectangularBox(evdb::View3D* view,
                                          double const* coordsLo,
                                          double const* coordsHi,
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

  void StandardDrawer::DrawGrids(evdb::View3D* view,
                                 double const* coordsLo,
                                 double const* coordsHi,
                                 int color,
                                 int width,
                                 int style) const
  {
    // uniform step size, each 25 cm except that at least 5 per plane
    double const gridStep = std::min(25.0,
                                     std::min({std::abs(coordsHi[0] - coordsLo[0]),
                                               std::abs(coordsHi[1] - coordsLo[1]),
                                               std::abs(coordsHi[2] - coordsLo[2])}) /
                                       5);

    // Grid running along x and y at constant z
    for (double z = coordsLo[2]; z <= coordsHi[2]; z += gridStep) {

      // across x, on bottom plane, fixed z
      TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
      gridt.SetPoint(0, coordsLo[0], coordsLo[1], z);
      gridt.SetPoint(1, coordsHi[0], coordsLo[1], z);

      // on right plane, across y, fixed z
      TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
      grids.SetPoint(0, coordsHi[0], coordsLo[1], z);
      grids.SetPoint(1, coordsHi[0], coordsHi[1], z);
    }

    // Grid running along z at constant x
    for (double x = coordsLo[0]; x <= coordsHi[0]; x += gridStep) {
      // fixed x, on bottom plane, across z
      TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
      gridt.SetPoint(0, x, coordsLo[1], coordsLo[2]);
      gridt.SetPoint(1, x, coordsLo[1], coordsHi[2]);
    }

    // Grid running along z at constant y
    for (double y = coordsLo[1]; y <= coordsHi[1]; y += gridStep) {
      // on right plane, fixed y, across z
      TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
      grids.SetPoint(0, coordsHi[0], y, coordsLo[2]);
      grids.SetPoint(1, coordsHi[0], y, coordsHi[2]);
    }

    return;
  }

  void StandardDrawer::DrawAxes(evdb::View3D* view,
                                double const* coordsLo,
                                double const* coordsHi,
                                int color,
                                int width,
                                int style) const
  {
    /*
     * Axes are drawn encompassing the whole detector volume,
     * the axis length being a fraction of the detector dimensions
     */
    double const vertexMargin = 0.06;
    double const axisLength = 0.40; // 20% of the shortest

    double const dx = (coordsHi[0] - coordsLo[0]);
    double const dy = (coordsHi[1] - coordsLo[1]);
    double const dz = (coordsHi[2] - coordsLo[2]);

    // axes origin
    double const x0 = coordsLo[0] - dx * vertexMargin;
    double const y0 = coordsLo[1] - dy * vertexMargin;
    double const z0 = coordsLo[2] - dz * vertexMargin;
    // axis length
    double const sz = axisLength * std::min({std::abs(dx), std::abs(dy), std::abs(dz)});

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

  DEFINE_ART_CLASS_TOOL(StandardDrawer)
}
