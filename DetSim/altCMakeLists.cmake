art_add_module(DumpRawDigits_module
     DumpRawDigits_module.cc
     )

art_add_module(SimWireAna_module
	SimWireAna_module.cc
	)

art_add_module(SimWire_module
	SimWire_module.cc
	)

art_add_module(SimWireT962_module
	SimWireT962_module.cc
	)

art_add_module(WienerFilterAna_module
	WienerFilterAna_module.cc
	)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     DumpRawDigits_module
     SimWireAna_module
     SimWire_module
     SimWireT962_module
     WienerFilterAna_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )
