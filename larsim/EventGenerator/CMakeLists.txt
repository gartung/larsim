art_make(LIB_LIBRARIES nusimdata_SimulationBase
                       cetlib
                       cetlib_except
                       ${FHICLCPP}
                       ${CLHEP_LIBRARIES}
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${ROOT_Core_LIBRARY}
                       ${ROOT_RIO_LIBRARY}
                       ${ROOT_Hist_LIBRARY}
                       ${ROOT_Physics_LIBRARY}
          MODULE_LIBRARIES larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_Simulation
                           lardataobj_Simulation
                           larcoreobj_SummaryData
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           nutools_RandomUtils_NuRandomService_service
                           larsim_EventGenerator
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ${CLHEP_LIBRARIES}
                           ${ROOT_EG_LIBRARY}
                           ${ROOT_Geom_LIBRARY}
                           ${ROOT_LIBRARIES}
        )


install_headers()
install_fhicl()
install_source()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(MARLEY)
add_subdirectory(MuonPropagation)
add_subdirectory(CORSIKA)

# See ubooneoffline repository for a functioning NuWro implementation.
# Leave it commented out here, as below.
##add_subdirectory(NuWro)
