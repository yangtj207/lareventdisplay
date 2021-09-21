cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    messagefacility::MF_MessageLogger
    ROOT::Core
    ROOT::EG
    ROOT::Physics
    canvas::canvas
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
