
art_make( BASENAME_ONLY
          LIBRARY_NAME  EventDisplay
          LIB_LIBRARIES ColorDrawingOptions_service
	                InfoTransfer_service
	                Filters
                        LArProperties_service
                        DetectorProperties_service
 			Simulation
			RawData
			RecoBase
			RecoObjects
			Utilities
			AnalysisBase
			AnalysisAlg
	                ${SIMULATIONBASE}
			${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_EVE_LIB_LIST}
			${ROOT_X3d}
          MODULE_LIBRARIES EventDisplay
			${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
          SERVICE_LIBRARIES 
			Geometry_service
			Geometry
			RecoBase
			${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
     )

install_headers()
install_fhicl()
install_source()
