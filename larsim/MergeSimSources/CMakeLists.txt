include_directories ( ${GEANT4_INCLUDE_DIRS} )
include_directories ( ${GENIE_INC}/GENIE )
include_directories ( ${LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( ${XERCES_C_INC} )

art_make(LIB_LIBRARIES  larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        lardata_Utilities
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib_except
                        ROOT::GenVector
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        ROOT::Geom
                        ROOT::MathCore
                        ${ROOT_LIBRARIES}
                        ${CLHEP_LIBRARIES}
                        ${Geant4_LIBRARIES}
           MODULE_LIBRARIES larsim_MergeSimSources
                        larsim_MCCheater_ParticleInventoryService_service
                        lardataalg_MCDumpers
                        larsim_Simulation  lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        ROOT::GenVector
                        ${ROOT_LIBRARIES}
                        ${CLHEP_LIBRARIES}
                        ${Geant4_LIBRARIES}
                )

install_headers()
install_fhicl()
