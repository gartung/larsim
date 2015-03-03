set(HEADERS
	PhotonLibrary.h
	PhotonVisibilityService.h
	)

add_library(PhotonPropagation SHARED
	${HEADERS}
	PhotonLibrary.cxx
	)

target_link_libraries(PhotonPropagation
     larsoft::RawData
     larsoft::Geometry
     art::art_Framework_Core
     art::art_Framework_IO_Sources
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     FNALCore::FNALCore
     ${ROOT_BASIC_LIB_LIST}
     ${ROOT_EG}
)


art_add_module(PhotonLibraryAnalyzer_module
	PhotonLibraryAnalyzer_module.cc
	)

art_add_service(PhotonVisibilityService_service
	PhotonVisibilityService_service.cc
	)

install(FILES ${HEADERS} DESTINATION
     ${CMAKE_INSTALL_INCLUDEDIR}/PhotonPropagation COMPONENT Development 
)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     PhotonPropagation
     PhotonLibraryAnalyzer_module
     PhotonVisibilityService_service
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )


add_subdirectory(LibraryBuildTools)
