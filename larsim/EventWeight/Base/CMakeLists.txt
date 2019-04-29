include_directories(${Boost_INCLUDE_DIRS})
art_make(LIB_LIBRARIES
         nutools_RandomUtils_NuRandomService_service
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         canvas
         cetlib_except
         ${FHICLCPP}
         ${CLHEP_LIBRARIES}
         ${ROOT_Core_LIBRARY}
         ${ROOT_Matrix_LIBRARY}
         Boost::system
         )

install_headers()
install_fhicl()
install_source()
