set(HEADERS
     MCRecoEdep.h
     MCRecoPart.h
     MCShowerRecoAlg.h
     MCShowerRecoPart.h
     MCTrackRecoAlg.h
     )

add_library(MCSTReco SHARED
     ${HEADERS}
     MCRecoEdep.cxx
     MCRecoPart.cxx
     MCShowerRecoAlg.cxx
     MCShowerRecoPart.cxx
     MCTrackRecoAlg.cxx
     )


target_link_libraries(MCSTReco
     larsoft::Simulation
     larsoft::SimulationBase
     larsoft::Utilities
     larsoft::Geometry
     FNALCore::FNALCore
     )

art_add_module(MCReco_module MCReco_module.cc)

art_add_module(ToyOneShowerGen_module ToyOneShowerGen_module.cc)

install(TARGETS
     MCSTReco
     MCReco_module
     ToyOneShowerGen_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

add_subdirectory(job)
