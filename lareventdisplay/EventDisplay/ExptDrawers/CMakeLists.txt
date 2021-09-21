cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    art::Framework_Services_Registry
    ROOT::Core
    ROOT::Physics
    canvas::canvas
    larcorealg_Geometry
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    nuevdb_EventDisplayBase
  )

install_headers()
install_source()
install_fhicl()
