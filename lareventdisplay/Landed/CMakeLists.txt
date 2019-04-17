#cet_find_library( SQLITE3  NAMES sqlite3  PATHS ENV NO_DEFAULT_PATH )
#cet_find_library( ZLIB  NAMES z  PATHS ENV NO_DEFAULT_PATH )

art_make( 
          MODULE_LIBRARIES 
            larcore_Geometry_Geometry_service
	    lareventdisplay_Landed_LandedSocket_service
            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            larcorealg_Geometry
            larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            nutools_NuReweight_art
            nusimdata_SimulationBase
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
            art_Persistency_Common
            art_Persistency_Provenance
            art_Utilities
	    canvas
            ${MF_MESSAGELOGGER}
            ${FHICLCPP}
            cetlib cetlib_except
            ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
	    ${SQLITE3}
	    ${G4VRML}
	    ${ZLIB_LIBRARIES}
	    nutools_G4Base
	    ${G4GLOBAL}   
	    ${G4INTERFACES}   
	    ${G4INTERCOMS}   
	    ${G4RUN}   
	    ${G4VIS_MANAGEMENT}   
	    ${G4MODELING}
	    ${G4TREE}
	    ${G4FR}
	    ${G4VISHEPREP}
	    ${G4RAYTRACER}
	    ${G4GMOCREN}
	    ${G4PHYSICSLISTS}
	  SERVICE_LIBRARIES
    	    ${ART_FRAMEWORK_CORE}
    	    art_Utilities canvas
    	    ${ART_ROOT_IO_ROOTINPUT_SOURCE}
    	    ${ART_ROOT_IO}
            art_Framework_Principal
    	    art_Persistency_Provenance
	    canvas
    	    ${MF_MESSAGELOGGER}
            cetlib cetlib_except
            pthread
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
install_fhicl()
install_source()

include_directories ( $ENV{G4INCLUDE} )
