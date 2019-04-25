include(${Geant4_USE_FILE})
include_directories ( ${GENIE_INC}/GENIE )
include_directories ( ${LIBXML2_INC})
include_directories ( ${GEANT4_INCLUDE_DIRS} )
include_directories ( ${XERCES_C_INC} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS} -Wno-shadow -Wno-unused-parameter")

art_make(LIB_LIBRARIES  larsim_PhotonPropagation
                        larsim_PhotonPropagation_PhotonVisibilityService_service
                        lardata_Utilities
                        larsim_Simulation
                        nutools_ParticleNavigation
                        lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        nutools_RandomUtils_NuRandomService_service
                        nutools_G4Base
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib_except
                        ${ROOT_GenVector_LIBRARY}
                        ${ROOT_EG_LIBRARY}
                        ${ROOT_TreePlayer_LIBRARY}
                        ${ROOT_FFTW_LIBRARY}
                        ${ROOT_Reflex_LIBRARY}
                        ${ROOT_EGPythia6_LIBRARY}
                        ${ROOT_Gui_LIBRARY}
                        ${ROOT_LIBRARIES}
                        ${ROOT_Geom_LIBRARY}
                        ${CLHEP_LIBRARIES}
                        ${Geant4_LIBRARIES}
           MODULE_LIBRARIES larsim_LArG4
                        larsim_MCCheater_ParticleInventoryService_service
                        lardataalg_MCDumpers
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nutools_RandomUtils_NuRandomService_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_Core_LIBRARY}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        ${ROOT_GenVector_LIBRARY}
                        nutools_G4Base
                        ${CLHEP_LIBRARIES}
                        ${Geant4_LIBRARIES}
                )

install_headers()
install_fhicl()
install_source()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
