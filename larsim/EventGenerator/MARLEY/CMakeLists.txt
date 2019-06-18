
# for c2: the problem is in fhiclcpp/types/detail/Sequence.icc
#cet_add_compiler_flags( CXX -Wno-missing-braces )

cet_find_simple_package(MARLEY HEADERS marley/Generator.hh)
string(REPLACE libMARLEY libMARLEY_ROOT MARLEY_ROOT ${MARLEY})

art_make(
      LIB_LIBRARIES
                        larsim_Simulation
                        lardataobj_Simulation
                        lardataobj_Simulation_dict
                        nusimdata_SimulationBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${MARLEY}
                        ${MARLEY_ROOT}
                        ${FHICLCPP}
                        ${MF_MESSAGELOGGER}
                        cetlib_except
                        ${ROOT_Geom_LIBRARY}
      MODULE_LIBRARIES
                        larsim_EventGenerator_MARLEY
                        larsim_Simulation
                        lardata_Utilities
                        lardataobj_Simulation
                        lardataobj_Simulation_dict
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        ${MARLEY}
                        ${MARLEY_ROOT}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Provenance
                        ${CLHEP_LIBRARIES}
                        ${ROOT_Geom_LIBRARY}
                        ${ROOT_Hist_LIBRARY}
                        ${ROOT_Tree_LIBRARY}
                        ${ROOT_RIO_LIBRARY}
                        ${Geant4_LIBRARIES}
        )

install_headers()
install_fhicl()
install_source()
