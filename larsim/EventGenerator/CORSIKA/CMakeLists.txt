include_directories(${IFDHC_INC})
find_library(IFDH name ifdh PATHS ${CMAKE_INSTALL_PREFIX} NO DEFAULTS)
art_make(MODULE_LIBRARIES
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nutools_RandomUtils_NuRandomService_service
         larsim_Simulation lardataobj_Simulation
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${CLHEP_LIBRARIES}
         ${SQLITE3}
         ${ROOT_EG_LIBRARY}
         ${IFDH})

install_headers()
install_fhicl()
install_source()

