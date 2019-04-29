art_make(LIB_LIBRARIES larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CLHEP_LIBRARIES}
                       cetlib_except
                       ${ROOT_EG_LIBRARY}
                       ${ROOT_RIO_LIBRARY}
                       ${ROOT_Hist_LIBRARY}
                       ${ROOT_Matrix_LIBRARY}
                       ${ROOT_Tree_LIBRARY}
                       ${ROOT_LIBRARIES}
          SERVICE_LIBRARIES larsim_PhotonPropagation
                       larsim_Simulation
                       nutools_ParticleNavigation
                       lardataobj_Simulation
                       larevt_Filters
                       lardataobj_RawData
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_LIBRARIES}
                       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                       ${MF_MESSAGELOGGER}
          MODULE_LIBRARIES larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        larsim_IonizationScintillation
                        lardataobj_Simulation
                        larsim_Simulation
                        nutools_ParticleNavigation
                        lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        nutools_RandomUtils_NuRandomService_service
                        art_Persistency_Provenance
                        ${MF_MESSAGELOGGER}
                        ${CLHEP_LIBRARIES}
                        ${ROOT_Gpad_LIBRARY}
                        ${ROOT_GenVector_LIBRARY}
                        ${ROOT_LIBRARIES}
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)
add_subdirectory(LibraryMappingTools)
