cet_enable_asserts()

art_make(TOOL_LIBRARIES
           art::Framework_Principal
           canvas::canvas
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           lareventdisplay_EventDisplay_ColorDrawingOptions_service
           lareventdisplay_EventDisplay_RecoDrawingOptions_service)

install_headers()
install_source()
install_fhicl()
