cet_enable_asserts()

cet_build_plugin(DrawLArVoxel3D art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay
  lareventdisplay::EventDisplay_SimulationDrawingOptions_service
  larsim::Simulation
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  nuevdb::EventDisplayBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::EG
  ROOT::Graf3d
)

cet_build_plugin(DrawSimEnergyDeposit3D art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay
  lareventdisplay::EventDisplay_SimulationDrawingOptions_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::Graf3d
)

cet_build_plugin(DrawSimPhoton3D art::tool
  LIBRARIES PRIVATE
  lareventdisplay::EventDisplay_ColorDrawingOptions_service
  lareventdisplay::EventDisplay_SimulationDrawingOptions_service
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::Simulation
  nuevdb::EventDisplayBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::Graf3d
  Eigen3::Eigen
)

install_headers()
install_source()
install_fhicl()
