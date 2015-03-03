art_add_module(FilterNoDirtNeutrinos_module
	FilterNoDirtNeutrinos_module.cc
)


art_add_module(FilterPrimaryPDG_module
	FilterPrimaryPDG_module.cc
)

install(TARGETS
     FilterNoDirtNeutrinos_module
     FilterPrimaryPDG_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )
