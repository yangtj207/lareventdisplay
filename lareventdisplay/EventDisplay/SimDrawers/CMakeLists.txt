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
    nusimdata::SimulationBase
    nuevdb::EventDisplayBase
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.13.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
install_fhicl()
