////////////////////////////////////////////////////////////////////////
/// \file   MicroBooNEDrawer.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"

#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"

#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/RawDrawingOptions.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"

#include <fstream>

namespace evd_tool
{

class MicroBooNEDrawer : IExperimentDrawer
{
public:
    explicit MicroBooNEDrawer(const fhicl::ParameterSet& pset);
    
    void DetOutline3D(evdb::View3D* view) override;
    
    ~MicroBooNEDrawer() {}
    
private:
    void configure(const fhicl::ParameterSet& pset);
    void DrawRectangularBox(evdb::View3D* view, double* coordsLo, double* coordsHi, int color=kGray, int width = 1, int style = 1);
    void DrawGrids(evdb::View3D* view, double* coordsLo, double* coordsHi, int color=kGray, int width = 1, int style = 1);
    void DrawAxes(evdb::View3D* view, double* coordsLo, double* coordsHi, int color=kGray, int width = 1, int style = 1);
    void DrawBadChannels(evdb::View3D* view, double* coords, int color, int width, int style);
    
    // Member variables from the fhicl file
    bool fThreeWindow;                 ///< true to draw rectangular box representing 3 windows
    bool fDrawGrid;                    ///< true to draw backing grid
    bool fDrawAxes;                    ///< true to draw coordinate axes
    bool fDrawBadChannels;             ///< true to draw bad channels
};
    
//----------------------------------------------------------------------
// Constructor.
MicroBooNEDrawer::MicroBooNEDrawer(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void MicroBooNEDrawer::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fThreeWindow     = pset.get< bool >("DrawThreeWindow", true);
    fDrawGrid        = pset.get< bool >("DrawGrid",        true);
    fDrawAxes        = pset.get< bool >("DrawAxes",        true);
    fDrawBadChannels = pset.get< bool >("DrawBadChannels", true);
    
    return;
}
    
//......................................................................
void MicroBooNEDrawer::DetOutline3D(evdb::View3D* view)
{
    art::ServiceHandle<geo::Geometry>         geo;
    
    // If requested, draw the outer three window volume first
    if (fThreeWindow)
    {
        double threeWinCoordsLo[] = {-2.*geo->DetHalfWidth(), -geo->DetHalfHeight(),               0.};
        double threeWinCoordsHi[] = { 4.*geo->DetHalfWidth(),  geo->DetHalfHeight(), geo->DetLength()};
        
        DrawRectangularBox(view, threeWinCoordsLo, threeWinCoordsHi, kGray);
    }
    
    // Now draw the standard volume
    double coordsLo[] = {                    0., -geo->DetHalfHeight(),               0.};
    double coordsHi[] = {2.*geo->DetHalfWidth(),  geo->DetHalfHeight(), geo->DetLength()};
    
    DrawRectangularBox(view, coordsLo, coordsHi, kRed, 2, 1);
    
    // It could be that we don't want to see the grids
    if (fDrawGrid)        DrawGrids(view, coordsLo, coordsHi, kGray+2, 1, 1);
    
    if (fDrawAxes)        DrawAxes(view, coordsLo, coordsHi, kBlue, 1, 1);
    
    if (fDrawBadChannels) DrawBadChannels(view, coordsHi, kGray, 1, 1);
    
    return;
}

void MicroBooNEDrawer::DrawRectangularBox(evdb::View3D* view, double* coordsLo, double* coordsHi, int color, int width, int style)
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

void MicroBooNEDrawer::DrawGrids(evdb::View3D* view, double* coordsLo, double* coordsHi, int color, int width, int style)
{
    double z = coordsLo[2];
    // Grid running along x and y at constant z
    for (;;) {
        TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
        gridt.SetPoint(0, coordsLo[0], coordsLo[1], z);
        gridt.SetPoint(1, coordsHi[0], coordsLo[1], z);
        
        TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
        grids.SetPoint(0, coordsHi[0], coordsLo[1], z);
        grids.SetPoint(1, coordsHi[0], coordsHi[1], z);
        
        z += 10.0;
        if (z>coordsHi[2]) break;
    }
    
    // Grid running along z at constant x
    double x = 0.0;
    for (;;) {
        TPolyLine3D& gridt = view->AddPolyLine3D(2, color, style, width);
        gridt.SetPoint(0, x, coordsLo[1], coordsLo[2]);
        gridt.SetPoint(1, x, coordsLo[1], coordsHi[2]);
        x += 10.0;
        if (x>coordsHi[0]) break;
    }
    
    // Grid running along z at constant y
    double y = 0.0;
    for (;;) {
        TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
        grids.SetPoint(0, coordsHi[0], y, coordsLo[2]);
        grids.SetPoint(1, coordsHi[0], y, coordsHi[2]);
        y += 10.0;
        if (y>coordsHi[1]) break;
    }
    y = -10.0;
    for (;;) {
        TPolyLine3D& grids = view->AddPolyLine3D(2, color, style, width);
        grids.SetPoint(0, coordsHi[0], y, coordsLo[2]);
        grids.SetPoint(1, coordsHi[0], y, coordsHi[2]);
        y -= 10.0;
        if (y<coordsLo[1]) break;
    }
    
    return;
}

void MicroBooNEDrawer::DrawAxes(evdb::View3D* view, double* coordsLo, double* coordsHi, int color, int width, int style)
{
    
    // Indicate coordinate system
    double x0 = -0.20;             // Center location of the key
    double y0 =  1.10*coordsLo[1]; // Center location of the key
    double z0 = -0.10*coordsHi[2]; // Center location of the key
    double sz =  0.20*coordsHi[2]; // Scale size of the key in z direction
    
    TPolyLine3D& xaxis = view->AddPolyLine3D(2, color, style, width);
    TPolyLine3D& yaxis = view->AddPolyLine3D(2, color, style, width);
    TPolyLine3D& zaxis = view->AddPolyLine3D(2, color, style, width);
    xaxis.SetPoint(0, x0,    y0, z0);
    xaxis.SetPoint(1, sz+x0, y0, z0);
    
    yaxis.SetPoint(0, x0, y0,     z0);
    yaxis.SetPoint(1, x0, y0+sz,  z0);
    
    zaxis.SetPoint(0, x0, y0, z0);
    zaxis.SetPoint(1, x0, y0, z0+sz);
    
    TPolyLine3D& xpoint = view->AddPolyLine3D(3, color, style, width);
    TPolyLine3D& ypoint = view->AddPolyLine3D(3, color, style, width);
    TPolyLine3D& zpoint = view->AddPolyLine3D(3, color, style, width);
    
    xpoint.SetPoint(0, 0.95*sz+x0, y0, z0-0.05*sz);
    xpoint.SetPoint(1, 1.00*sz+x0, y0, z0);
    xpoint.SetPoint(2, 0.95*sz+x0, y0, z0+0.05*sz);
    
    ypoint.SetPoint(0, x0, 0.95*sz+y0, z0-0.05*sz);
    ypoint.SetPoint(1, x0, 1.00*sz+y0, z0);
    ypoint.SetPoint(2, x0, 0.95*sz+y0, z0+0.05*sz);
    
    zpoint.SetPoint(0, x0-0.05*sz, y0, 0.95*sz+z0);
    zpoint.SetPoint(1, x0+0.00*sz, y0, 1.00*sz+z0);
    zpoint.SetPoint(2, x0+0.05*sz, y0, 0.95*sz+z0);
    
    TPolyLine3D& zleg = view->AddPolyLine3D(4, color, style, width);
    zleg.SetPoint(0,  x0-0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(1,  x0+0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(2,  x0-0.05*sz, y0-0.05*sz, z0+1.05*sz);
    zleg.SetPoint(3,  x0+0.05*sz, y0-0.05*sz, z0+1.05*sz);
    
    TPolyLine3D& yleg = view->AddPolyLine3D(5, color, style, width);
    yleg.SetPoint(0,  x0-0.05*sz, y0+1.15*sz, z0);
    yleg.SetPoint(1,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(2,  x0+0.00*sz, y0+1.05*sz, z0);
    yleg.SetPoint(3,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(4,  x0+0.05*sz, y0+1.15*sz, z0);
    
    TPolyLine3D& xleg = view->AddPolyLine3D(7, color, style, width);
    xleg.SetPoint(0,  x0+1.05*sz, y0+0.05*sz, z0-0.05*sz);
    xleg.SetPoint(1,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(2,  x0+1.05*sz, y0+0.05*sz, z0+0.05*sz);
    xleg.SetPoint(3,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(4,  x0+1.05*sz, y0-0.05*sz, z0-0.05*sz);
    xleg.SetPoint(5,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(6,  x0+1.05*sz, y0-0.05*sz, z0+0.05*sz);
    
    return;
}

void MicroBooNEDrawer::DrawBadChannels(evdb::View3D* view, double* coords, int color, int width, int style)
{
    art::ServiceHandle<geo::Geometry>          geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawOpt;
    
    lariov::ChannelStatusProvider const& channelStatus
    = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    // We want to translate the wire position to the opposite side of the TPC...
    for(size_t viewNo = 0; viewNo < geo->Nviews(); viewNo++)
    {
        for(size_t wireNo = 0; wireNo < geo->Nwires(viewNo); wireNo++)
        {
            geo::WireID wireID = geo::WireID(rawOpt->fCryostat, rawOpt->fTPC, viewNo, wireNo);
            
            raw::ChannelID_t channel = geo->PlaneWireToChannel(wireID);
            
            if (channelStatus.IsBad(channel))
            {
                const geo::WireGeo* wireGeo = geo->WirePtr(wireID);
                
                double wireStart[3];
                double wireEnd[3];
                
                wireGeo->GetStart(wireStart);
                wireGeo->GetEnd(wireEnd);
                
                TPolyLine3D& pl = view->AddPolyLine3D(2, color, style, width);
                pl.SetPoint(0, coords[0]-0.5, wireStart[1], wireStart[2]);
                pl.SetPoint(1, coords[0]-0.5,   wireEnd[1],   wireEnd[2]);
            }
        }
    }
    
    return;
}
    
    
DEFINE_ART_CLASS_TOOL(MicroBooNEDrawer)
}
