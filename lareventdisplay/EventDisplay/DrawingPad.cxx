///
/// \file    DrawingPad.cxx
/// \brief   Base class for all event display drawing pads
/// \author  messier@indiana.edu
///
#include "lareventdisplay/EventDisplay/DrawingPad.h"
#include "TPad.h"
#include "nutools/EventDisplayBase/evdb.h"
#include "lareventdisplay/EventDisplay/HeaderDrawer.h"
#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/RawDataDrawer.h"
#include "lareventdisplay/EventDisplay/AnalysisBaseDrawer.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/make_tool.h"

namespace evd{

  // Declare singleton HitSelector

  HitSelector * gTheHitSelector;

  /// Create a drawing pad for the event display
  ///
  /// @param nm : Name of the TPad
  /// @param ti : Title of the TPad
  /// @param x1 : Relative x position (0-1) of lower left  corner
  /// @param y1 : Relative y position (0-1) of lower left  corner
  /// @param x2 : Relative x position (0-1) of upper right corner
  /// @param y2 : Relative y position (0-1) of upper right corner
  ///
  DrawingPad::DrawingPad(const char* nm,
			 const char* ti,
			 double x1, double y1,
			 double x2, double y2)
    : fPad(0)
    , fHeaderDraw(0) //Every pointer checked for a 0 value in the destructor should be set to 0 here.  aoliv23@lsu.edu
    , fSimulationDraw(0)
    , fRawDataDraw(0)
    , fRecoBaseDraw(0)
    , fAnalysisBaseDraw(0)
  {
    fPad = new TPad(nm,ti,x1,y1,x2,y2);
    fPad->Draw();
    fPad->cd();
  }

  //......................................................................

  DrawingPad::~DrawingPad()
  {
    if (fHeaderDraw)       { delete fHeaderDraw;       fHeaderDraw       = 0; }
    if (fSimulationDraw)   { delete fSimulationDraw;   fSimulationDraw   = 0; }
    if (fRawDataDraw)      { delete fRawDataDraw;      fRawDataDraw      = 0; }
    if (fRecoBaseDraw)     { delete fRecoBaseDraw;     fRecoBaseDraw     = 0; }
    if (fAnalysisBaseDraw) { delete fAnalysisBaseDraw; fAnalysisBaseDraw = 0; }
    //  if (fHitSelector)   { delete fHitSelector;     fHitSelector   = 0; }
    if (fPad)              { delete fPad;              fPad = 0;              }
  }

  // //......................................................................

  //......................................................................

  ///
  /// Provide access to the drawer for the detector geometry
  ///
  HeaderDrawer* DrawingPad::HeaderDraw()
  {
    if (fHeaderDraw==0) fHeaderDraw = new HeaderDrawer();
    return fHeaderDraw;
  }

  ///
  /// Provide access to the drawer for the detector geometry
  ///
  evd_tool::IExperimentDrawer* DrawingPad::GeometryDraw()
  {
    if (fGeometryDraw==0)
    {
        art::ServiceHandle<evd::EvdLayoutOptions const> layoutOptions;
        const fhicl::ParameterSet&                pset = layoutOptions->fParameterSet;

        fGeometryDraw = art::make_tool<evd_tool::IExperimentDrawer>(pset.get<fhicl::ParameterSet>("Experiment3DDrawer"));
    }
    return fGeometryDraw.get();
  }

  ///
  /// Provide access to the drawer for the Simulation classes
  ///
  SimulationDrawer* DrawingPad::SimulationDraw()
  {
    if (fSimulationDraw==0) fSimulationDraw = new SimulationDrawer();
    return fSimulationDraw;

  }

  ///
  /// Provide access to the drawer for the RawData classes
  ///
  RawDataDrawer* DrawingPad::RawDataDraw()
  {
    if (fRawDataDraw==0) fRawDataDraw = new RawDataDrawer();
    return fRawDataDraw;
  }

  //......................................................................

  ///
  /// Provide access to the drawer for RecoBase classes
  ///
  RecoBaseDrawer* DrawingPad::RecoBaseDraw()
  {
    if (fRecoBaseDraw==0) fRecoBaseDraw = new RecoBaseDrawer();
    return fRecoBaseDraw;

  }

  //......................................................................

  ///
  /// Provide access to the drawer for AnalysisBase classes
  ///
  AnalysisBaseDrawer* DrawingPad::AnalysisBaseDraw()
  {
    if (fAnalysisBaseDraw==0) fAnalysisBaseDraw = new AnalysisBaseDrawer();
    return fAnalysisBaseDraw;
  }

  //......................................................................
  //......................................................................

  ///
  /// Provide access to the HitSelector
  ///
  HitSelector* DrawingPad::HitSelectorGet()
  {
    if (gTheHitSelector==0) gTheHitSelector = new HitSelector();
    return gTheHitSelector;
  }

}// namespace
////////////////////////////////////////////////////////////////////////
