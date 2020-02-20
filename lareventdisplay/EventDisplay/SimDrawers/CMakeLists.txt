cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ROOT::Core
    ROOT::EG
    ROOT::Physics
    canvas
    larcorealg_Geometry
    lareventdisplay_EventDisplay
    lareventdisplay_EventDisplay_ColorDrawingOptions_service
    lareventdisplay_EventDisplay_SimulationDrawingOptions_service
    larsim_Simulation
    nusimdata_SimulationBase
    nuevdb_EventDisplayBase
)

install_headers()
install_source()
install_fhicl()
