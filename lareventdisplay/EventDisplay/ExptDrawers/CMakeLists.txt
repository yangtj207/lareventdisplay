cet_enable_asserts()

cet_build_plugin(ICARUSDrawer art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_RawDrawingOptions_service
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  nuevdb::EventDisplayBase
  art::Framework_Services_Registry
  ROOT::Graf3d
)

cet_build_plugin(MicroBooNEDrawer art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_RawDrawingOptions_service
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  nuevdb::EventDisplayBase
  art::Framework_Services_Registry
  ROOT::Graf3d
)

cet_build_plugin(ProtoDUNEDrawer art::tool
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  nuevdb::EventDisplayBase
  fhiclcpp::fhiclcpp
  ROOT::Graf3d
)

cet_build_plugin(StandardDrawer art::tool
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  larcorealg::Geometry
  nuevdb::EventDisplayBase
  fhiclcpp::fhiclcpp
  ROOT::Graf3d
)

install_headers()
install_source()
install_fhicl()
