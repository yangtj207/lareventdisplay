cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    art::Framework_Services_Registry
    ROOT::Core
    ROOT::Physics
    canvas::canvas
    larcorealg_Geometry
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    nuevdb::EventDisplayBase
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.13.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
install_fhicl()
