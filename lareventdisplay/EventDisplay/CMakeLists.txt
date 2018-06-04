
include(CetRootCint)

set( PACKAGE EventDisplay )
FILE( GLOB src_files *.cxx )

cet_rootcint( ${PACKAGE} )

add_subdirectory(ExptDrawers)

art_make_library( 
        	  SOURCE ${src_files} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
                  LIBRARIES lareventdisplay_EventDisplay_ColorDrawingOptions_service
                            lareventdisplay_EventDisplay_RecoDrawingOptions_service
	                    lareventdisplay_EventDisplay_InfoTransfer_service
	                    lareventdisplay_EventDisplay_AnalysisDrawingOptions_service
	                    lareventdisplay_EventDisplay_EvdLayoutOptions_service
	                    lareventdisplay_EventDisplay_RawDrawingOptions_service
	                    lareventdisplay_EventDisplay_SimulationDrawingOptions_service
	                    larevt_Filters
          larsim_MCCheater_BackTrackerService_service
          larsim_MCCheater_ParticleInventoryService_service
 			    larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			    lardataobj_RawData
			    lardataobj_RecoBase
			    lardata_RecoObjects larreco_Deprecated
			    lardata_ArtDataHelper
			    lardata_Utilities
			    lardataobj_AnalysisBase
			    larreco_Calorimetry
		 	    larevt_CalibrationDBI_IOVData
			    larevt_CalibrationDBI_Providers
	                    nusimdata_SimulationBase
			    nutools_EventDisplayBase
                            ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    art_Persistency_Common
			    art_Utilities
			    canvas
			    ${MF_MESSAGELOGGER}
			    
	                    ${FHICLCPP}
			    cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
			    ${ROOT_X3d}
                )

simple_plugin ( GraphCluster module
                        lareventdisplay_EventDisplay
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                
     )

simple_plugin ( EVD module
                        lareventdisplay_EventDisplay
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                
     )

simple_plugin ( AnalysisDrawingOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( ColorDrawingOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( EvdLayoutOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( InfoTransfer service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( RawDrawingOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( RecoDrawingOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

simple_plugin ( SimulationDrawingOptions service
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			nutools_EventDisplayBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                
     )

install_headers()
install_fhicl()
install_source()
