cet_enable_asserts()

cet_write_plugin_builder(lar::SpacePoint3DDrawer art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  lareventdisplay::SpacePoint3DDrawer
)

include(lar::SpacePoint3DDrawer)

cet_build_plugin(Edge3DDrawer art::tool
  IMPL_SOURCE
)

cet_build_plugin(OpFlash3DDrawer art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Graf3d
  Eigen3::Eigen
)

cet_build_plugin(OpHit3DDrawer art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nuevdb::EventDisplayBase
  art::Framework_Principal
  art::Framework_Services_Registry
  ROOT::Graf3d
  Eigen3::Eigen
)

cet_build_plugin(PCA3DDrawer art::tool
  IMPL_SOURCE
)

cet_build_plugin(PFParticle3DDrawer art::tool
  LIBRARIES PRIVATE
  lardataalg::DetectorInfo
)

cet_build_plugin(SpacePoint3DDrawerAsymmetry lar::SpacePoint3DDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Graf3d
  ROOT::MathCore
)

cet_build_plugin(SpacePoint3DDrawerChiSquare lar::SpacePoint3DDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Graf3d
)

cet_build_plugin(SpacePoint3DDrawerHitCharge lar::SpacePoint3DDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Graf3d
  ROOT::MathCore
)

cet_build_plugin(SpacePoint3DDrawerStandard lar::SpacePoint3DDrawer
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_RecoDrawingOptions_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Graf3d
)

install_headers()
install_source()
install_fhicl()
