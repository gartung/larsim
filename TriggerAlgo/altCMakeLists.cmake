set(TriggerAlgo_HEADERS
	TriggerAlgoBase.h
	TriggerAlgoMicroBoone.h
	TriggerTypes.hh
	)

art_add_service(TriggerAlgoBase_service SHARED
	TriggerAlgoBase_service.cc
	)

art_add_service(TriggerAlgoMicroBoone_service SHARED
	TriggerAlgoMicroBoone_service.cc
	)

install(TARGETS
     TriggerAlgoBase_service
     TriggerAlgoMicroBoone_service
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${TriggerAlgo_HEADERS} DESTINATION
     ${CMAKE_INSTALL_INCLUDEDIR}/TriggerAlgo COMPONENT Development )

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)
