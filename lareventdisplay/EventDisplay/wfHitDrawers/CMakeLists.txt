cet_enable_asserts()

cet_write_plugin_builder(lar::WFHitDrawer art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  lareventdisplay::WFHitDrawer
)

cet_write_plugin_builder(lar::WaveformDrawer art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  lareventdisplay::WaveformDrawer
)

include(lar::WFHitDrawer)
include(lar::WaveformDrawer)

cet_build_plugin(DrawGausHits lar::WFHitDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::Graf3d
  ROOT::Hist
)

cet_build_plugin(DrawRawHist lar::WaveformDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_RawDrawingOptions_service
  larevt::DetPedestalProvider
  larevt::DetPedestalService
  lardataobj::RawData
  larcorealg::Geometry
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
)

cet_build_plugin(DrawSkewHits lar::WFHitDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_RawDrawingOptions_service
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Graf3d
)

cet_build_plugin(DrawWireData lar::WaveformDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Graf3d
)

cet_build_plugin(DrawWireHist lar::WaveformDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_RawDrawingOptions_service
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Hist
)

install_headers()
install_source()
install_fhicl()
