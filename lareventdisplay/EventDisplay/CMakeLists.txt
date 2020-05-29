include(CetRootCint)

set(PACKAGE EventDisplay)
FILE(GLOB src_files *.cxx)

cet_rootcint(${PACKAGE})

add_subdirectory(ExptDrawers)
add_subdirectory(wfHitDrawers)
add_subdirectory(SimDrawers)
add_subdirectory(3DDrawers)

art_make_library(
  SOURCE
    ${src_files}
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
  LIBRARIES
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_PERSISTENCY_PROVENANCE}
    ${MF_MESSAGELOGGER}
    ROOT::Core
    ROOT::EG
    ROOT::Gpad
    ROOT::Graf3d
    ROOT::Graf
    ROOT::Gui
    ROOT::Hist
    ROOT::Physics
    ROOT::RIO
    canvas
    cetlib_except
    lardataalg_DetectorInfo
    larcorealg_Geometry
    lardata_ArtDataHelper
    lardata_Utilities
    lardataobj_RawData
    lareventdisplay_EventDisplay_AnalysisDrawingOptions_service
    lareventdisplay_EventDisplay_ColorDrawingOptions_service
    lareventdisplay_EventDisplay_EvdLayoutOptions_service
    lareventdisplay_EventDisplay_InfoTransfer_service
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    lareventdisplay_EventDisplay_RecoDrawingOptions_service
    lareventdisplay_EventDisplay_SimulationDrawingOptions_service
    larsim_MCCheater_ParticleInventoryService_service
    larsim_Simulation
    nusimdata_SimulationBase
    nuevdb_EventDisplayBase
)

simple_plugin(GraphCluster "module" lareventdisplay_EventDisplay)
simple_plugin(EVD "module" lareventdisplay_EventDisplay)

simple_plugin(AnalysisDrawingOptions "service" nuevdb_EventDisplayBase)
simple_plugin(EvdLayoutOptions "service" nuevdb_EventDisplayBase)
simple_plugin(RawDrawingOptions "service" nuevdb_EventDisplayBase)
simple_plugin(RecoDrawingOptions "service" nuevdb_EventDisplayBase)
simple_plugin(SimulationDrawingOptions "service" nuevdb_EventDisplayBase)

simple_plugin(ColorDrawingOptions "service"
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ROOT::Core
              ROOT::Graf3d
              ROOT::Hist
              ROOT::MathCore
              larcorealg_Geometry
              lardataobj_RecoBase
              nuevdb_EventDisplayBase)

simple_plugin(InfoTransfer "service"
              larcorealg_Geometry
              nuevdb_EventDisplayBase
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ROOT::Core)

install_headers()
install_fhicl()
install_source()
