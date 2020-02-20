cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ROOT::Core
    ROOT::Physics
    canvas
    larcorealg_Geometry
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    nuevdb_EventDisplayBase
  )

install_headers()
install_source()
install_fhicl()
