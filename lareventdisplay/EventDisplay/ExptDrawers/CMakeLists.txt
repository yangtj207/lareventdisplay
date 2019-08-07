cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ROOT_CORE}
    ${ROOT_PHYSICS}
    canvas
    larcorealg_Geometry
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    nutools_EventDisplayBase
  )

install_headers()
install_source()
install_fhicl()
