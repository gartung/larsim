
include_directories ( ${CRY_DIR}/src )

art_make(MODULE_LIBRARIES
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nutools_RandomUtils_NuRandomService_service
         nutools_EventGeneratorBase_CRY
         ${CRY}
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${CLHEP_LIBRARIES}
         ${MF_MESSAGELOGGER}
         ${ROOT_Hist_LIBRARY}
         ${ROOT_Physics_LIBRARY})

install_headers()
install_fhicl()
install_source()

