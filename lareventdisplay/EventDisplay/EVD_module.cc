////////////////////////////////////////////////////////////////////////
/// \file EVD_module.cc
///
/// \author  jpaley@anl.gov

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

//LArSoft includes
#include "lareventdisplay/EventDisplay/CalorView.h"
#include "lareventdisplay/EventDisplay/Display3DView.h"
#include "lareventdisplay/EventDisplay/Ortho3DView.h"
#include "lareventdisplay/EventDisplay/TWQMultiTPCProjection.h"
#include "lareventdisplay/EventDisplay/TWQProjectionView.h"
#include "lareventdisplay/EventDisplay/PhDetView.h"
#include "nuevdb/EventDisplayBase/DisplayWindow.h"

// Framework includes
namespace art {
  class Event;
}
namespace fhicl {
  class ParameterSet;
}

#if defined __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-private-field"
#endif

/// The Event Display
namespace evd {

  class Canvas;

  /// a class for transporting photons in a roughly realistic way
  class EVD : public art::EDAnalyzer {
  public:
    explicit EVD(fhicl::ParameterSet const& pset);
    virtual ~EVD();

    void analyze(art::Event const& evt);
    void beginJob();

  private:
    bool fWindowsDrawn; ///< flag for whether windows are already drawn
  };
}

// Builder for the TWQProjectionView canvas
static evdb::Canvas* mk_twqprojectionview_canvas(TGMainFrame* mf)
{
  return new evd::TWQProjectionView(mf);
}

// Builder for the TWQMultiTPCProjectionView canvas
static evdb::Canvas* mk_twqmtpcprojectionview_canvas(TGMainFrame* mf)
{
  return new evd::TWQMultiTPCProjectionView(mf);
}

// Builder for the Display3D view
static evdb::Canvas* mk_display3d_canvas(TGMainFrame* mf)
{
  return new evd::Display3DView(mf);
}

// Builder for the Ortho view
static evdb::Canvas* mk_ortho3d_canvas(TGMainFrame* mf)
{
  return new evd::Ortho3DView(mf);
}

// Builder for the Calor view
static evdb::Canvas* mk_calor_canvas(TGMainFrame* mf)
{
  return new evd::CalorView(mf);
}

// Builder for the Photon Detector view
static evdb::Canvas* mk_phdet_canvas(TGMainFrame* mf)
{
  return new evd::PhDetView(mf);
}

// // Builder for the MCTruth view
// static evdb::ObjListCanvas* mk_mctrue_canvas(TGMainFrame* mf)
// {
//   return new evd::MCTrueView(mf);
// }

namespace evd {

  //----------------------------------------------------
  EVD::EVD(fhicl::ParameterSet const& pset) : EDAnalyzer(pset), fWindowsDrawn(false) {}

  //----------------------------------------------------
  EVD::~EVD() {}

  //----------------------------------------------------
  void EVD::beginJob()
  {
    // Register the list of windows used by the event display
    evdb::DisplayWindow::Register("Time vs Wire, Charge View",
                                  "Time vs Wire, Charge View",
                                  700,
                                  700,
                                  mk_twqprojectionview_canvas);

    evdb::DisplayWindow::Register("Time vs Wire, Charge View, Multi-TPC",
                                  "Time vs Wire, Charge View, Multi-TPC",
                                  700,
                                  700,
                                  mk_twqmtpcprojectionview_canvas);

    evdb::DisplayWindow::Register("Display3D", "Display3D", 700, 700, mk_display3d_canvas);

    evdb::DisplayWindow::Register("Ortho3D", "Ortho3D", 700, 700, mk_ortho3d_canvas);

    evdb::DisplayWindow::Register("Calorimetry", "Calorimetry", 700, 700, mk_calor_canvas);

    evdb::DisplayWindow::Register("Photon Detector", "Photon Detector", 700, 700, mk_phdet_canvas);

    //     evdb::ListWindow::Register("MC Particle List",
    // 			       "MC Particle List",
    // 			       400,
    // 			       800,
    // 			       mk_mctrue_canvas);

    // Open up the main display window and run
    evdb::DisplayWindow::OpenWindow(0);
  }

  //----------------------------------------------------
  void EVD::analyze(const art::Event& /*evt*/) {}

} //namespace

namespace evd {

  DEFINE_ART_MODULE(EVD)

} // namespace evd
#if defined __clang__
#pragma clang diagnostic pop
#endif
