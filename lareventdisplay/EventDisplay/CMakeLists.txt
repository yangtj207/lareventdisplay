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
    ${ROOT_CORE}
    ${ROOT_EG}
    ${ROOT_GPAD}
    ${ROOT_GRAF3D}
    ${ROOT_GRAF}
    ${ROOT_GUI}
    ${ROOT_HIST}
    ${ROOT_PHYSICS}
    ${ROOT_RIO}
    canvas
    cetlib_except
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
    nutools_EventDisplayBase
)

simple_plugin(GraphCluster "module" lareventdisplay_EventDisplay)
simple_plugin(EVD "module" lareventdisplay_EventDisplay)

simple_plugin(AnalysisDrawingOptions "service" nutools_EventDisplayBase)
simple_plugin(EvdLayoutOptions "service" nutools_EventDisplayBase)
simple_plugin(RawDrawingOptions "service" nutools_EventDisplayBase)
simple_plugin(RecoDrawingOptions "service" nutools_EventDisplayBase)
simple_plugin(SimulationDrawingOptions "service" nutools_EventDisplayBase)

simple_plugin(ColorDrawingOptions "service"
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ${ROOT_CORE}
              ${ROOT_GRAF3D}
              ${ROOT_HIST}
              ${ROOT_MATHCORE}
              larcorealg_Geometry
              lardataobj_RecoBase
              nutools_EventDisplayBase)

simple_plugin(InfoTransfer "service"
              nutools_EventDisplayBase
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_PERSISTENCY_PROVENANCE}
              ${MF_MESSAGELOGGER}
              ${ROOT_CORE})

install_headers()
install_fhicl()
install_source()
