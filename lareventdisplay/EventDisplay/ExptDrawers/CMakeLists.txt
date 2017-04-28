cet_enable_asserts()

set( evd_tool_lib_list
			larcore_Geometry
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			larcore_Geometry_Geometry_service
			lardata_Utilities
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common
			art_Persistency_Provenance
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			art_Utilities
			canvas
            nutools_EventDisplayBase
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
            ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
            ${ROOT_EVE_LIB_LIST}
            ${ROOT_X3d}
    )


simple_plugin( StandardDrawer  "tool"
               ${evd_tool_lib_list}
               BASENAME_ONLY
             )

simple_plugin( MicroBooNEDrawer  "tool"
               ${evd_tool_lib_list}
	           BASENAME_ONLY
             )

simple_plugin( ICARUSDrawer  "tool"
               ${evd_tool_lib_list}
               BASENAME_ONLY
             )

install_headers()
install_source()
install_fhicl()
