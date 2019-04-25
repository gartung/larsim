include_directories ( ${GENIE_INC}/GENIE )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGENIE_PRE_R3")

art_make(LIB_LIBRARIES larsim_EventWeight_Base
                       lardata_Utilities
                       nutools_NuReweight_art
                       nutools_NuReweight
                       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
                       nusimdata_SimulationBase
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${FHICLCPP}
                       ${CLHEP_LIBRARIES})

install_headers()
install_fhicl()
install_source()
