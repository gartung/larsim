art_add_module(FileMuons_module
	FileMuons_module.cc
	)

art_add_module(LightSource_module
	LightSource_module.cc
	)

art_add_module(NDKGen_module
	NDKGen_module.cc
	)

art_add_module(NUANCEGen_module
	NUANCEGen_module.cc
	)

art_add_module(RadioGen_module
	RadioGen_module.cc
	)

art_add_module(SingleGen_module
	SingleGen_module.cc
	)

art_add_module(TextFileGen_module
	TextFileGen_module.cc
	)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     FileMuons_module
     LightSource_module
     NDKGen_module
     NUANCEGen_module
     RadioGen_module
     SingleGen_module
     TextFileGen_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )



add_subdirectory(CRY)
add_subdirectory(GENIE)
