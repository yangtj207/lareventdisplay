cet_enable_asserts()

art_make(TOOL_LIBRARIES
           art::Framework_Principal
           canvas::canvas
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           lareventdisplay_EventDisplay_ColorDrawingOptions_service
           lareventdisplay_EventDisplay_RecoDrawingOptions_service) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.13.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
install_fhicl()
