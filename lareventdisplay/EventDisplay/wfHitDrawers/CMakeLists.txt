cet_enable_asserts()

art_make(
  TOOL_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    ROOT::Core
    canvas::canvas
    lardata_ArtDataHelper
    lardataobj_RawData
    lareventdisplay_EventDisplay_ColorDrawingOptions_service
    lareventdisplay_EventDisplay_RawDrawingOptions_service
    lareventdisplay_EventDisplay_RecoDrawingOptions_service
    nuevdb_EventDisplayBase
)

install_headers()
install_source()
install_fhicl()
