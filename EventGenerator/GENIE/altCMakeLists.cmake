include_directories ( "${GENIE_INC}/GENIE" )
 
art_add_module(GENIEGen_module
	GENIEGen_module.cc
	)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     GENIEGen_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )
