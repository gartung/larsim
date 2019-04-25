
add_subdirectory(job)

art_make(LIB_LIBRARIES larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       larsim_Simulation
                       nutools_ParticleNavigation
                       lardataobj_Simulation
                       nusimdata_SimulationBase
                       larcoreobj_SummaryData
                       lardata_Utilities
                       lardataobj_MCBase
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${ROOT_Geom_LIBRARY}
                       ${ROOT_XMLIO_LIBRARY}
                       ${ROOT_Gdml_LIBRARY}
          MODULE_LIBRARIES larsim_MCSTReco
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation
                           nutools_ParticleNavigation
                           lardataobj_Simulation
                           nusimdata_SimulationBase
                           larcoreobj_SummaryData
                           lardata_Utilities
                           lardataobj_MCBase
                           nutools_RandomUtils_NuRandomService_service
                           art_Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ${ROOT_Geom_LIBRARY}
                           ${ROOT_XMLIO_LIBRARY}
                           ${ROOT_Gdml_LIBRARY}
                           ${ROOT_EG_LIBRARY}
                           ${ROOT_Hist_LIBRARY}
                           ${CLHEP_LIBRARIES}
                        )

install_headers()
install_fhicl()
install_source()
