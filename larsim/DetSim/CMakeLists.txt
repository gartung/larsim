simple_plugin(SimWireAna "module"
                         lardataobj_RawData
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                         art_Persistency_Provenance
                         ${MF_MESSAGELOGGER}
                         ${ROOT_Hist_LIBRARY})

simple_plugin(SimWire "module"
                      larsim_Simulation lardataobj_Simulation
                      lardataobj_RawData
                      larcorealg_Geometry
                      larcore_Geometry_Geometry_service
                      lardata_Utilities_LArFFT_service
                      lardata_Utilities
                      nutools_RandomUtils_NuRandomService_service
                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
                      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                      art_Persistency_Provenance
                      ${MF_MESSAGELOGGER}
                      ${CLHEP_LIBRARIES}
                      ${ROOT_MathCore_LIBRARY}
                      ${ROOT_Hist_LIBRARY})


simple_plugin(WienerFilterAna "module"
                              lardataobj_RawData
                              larcorealg_Geometry
                              larcore_Geometry_Geometry_service
                              ${ART_FRAMEWORK_CORE}
                              ${ART_FRAMEWORK_SERVICES_REGISTRY}
                              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                              art_Persistency_Provenance
                              ${MF_MESSAGELOGGER}
                              ${ROOT_Hist_LIBRARY})

simple_plugin(RecompressRawDigits "module"
              lardataobj_RawData
              art_Persistency_Provenance)

install_headers()
install_fhicl()
install_source()
