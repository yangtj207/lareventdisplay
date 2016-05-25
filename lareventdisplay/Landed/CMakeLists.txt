cet_find_library( SQLITE3  NAMES sqlite3  PATHS ENV NO_DEFAULT_PATH )
cet_find_library( ZLIB  NAMES z  PATHS ENV NO_DEFAULT_PATH )
cet_find_library( G4BASE  NAMES G4Base  PATHS ENV NO_DEFAULT_PATH )

art_make( BASENAME_ONLY
          LIBRARY_NAME Landed
          MODULE_LIBRARIES 
            larcore_Geometry_Geometry_service
	    LandedSocket_service
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            larcore_Geometry
            larsim_Simulation
            lardata_Utilities
            larevt_Filters
            lardata_RawData
            NuReweightArt
            ${SIMULATIONBASE}
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_PERSISTENCY_COMMON}
            ${ART_PERSISTENCY_PROVENANCE}
            ${ART_UTILITIES}
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${CLHEP}
            ${ROOT_BASIC_LIB_LIST}
	    ${SQLITE3}
	    ${ZLIB}
	    ${G4VRML}
	    ${G4BASE}
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
    	    ${ART_UTILITIES}
    	    ${ART_FRAMEWORK_IO_ROOT_ROOTINPUT_SOURCE}
    	    ${ART_FRAMEWORK_IO_ROOT}
    	    ${ART_PERSISTENCY_PROVENANCE}
    	    ${MF_MESSAGELOGGER}
            ${ROOT_BASIC_LIB_LIST}
          )

install_headers()
install_fhicl()
# install_source()

include_directories ( $ENV{G4INCLUDE} )
