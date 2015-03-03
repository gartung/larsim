art_add_service(BackTracker_service SHARED
	BackTracker_service.cc
	)

art_add_module(BackTrackerLoader_module
	BackTrackerLoader_module.cc
	)


art_add_module(CheckBackTracking_module 
	CheckBackTracking_module.cc
	)

art_add_module(RecoCheckAna_module
	BackTracker.h
	RecoCheckAna_module.cc
	)

install(TARGETS 
     BackTracker_service
     BackTrackerLoader_module 
     CheckBackTracking_module
     RecoCheckAna_module 
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(FILES BackTracker.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MCCheater
     COMPONENT Development)
