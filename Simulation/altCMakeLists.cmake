include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set(HEADERS
	AuxDetSimChannel.h
	BeamGateInfo.h
	BeamTypes.h
	EmEveIdCalculator.h
	EveIdCalculator.h
	LArG4Parameters.h
	LArVoxelCalculator.h
	LArVoxelData.h
	LArVoxelID.h
	LArVoxelList.h
	ParticleHistory.h
	ParticleList.h
	PhotonVoxels.h
	SimChannel.h
	SimListUtils.h
	SimPhotons.h
	sim.h
	)

add_library(Simulation SHARED
	${HEADERS}
	AuxDetSimChannel.cxx
	EmEveIdCalculator.cxx
	EveIdCalculator.cxx
	LArVoxelData.cxx
	LArVoxelID.cxx
	LArVoxelList.cxx
	ParticleHistory.cxx
	ParticleList.cxx
	PhotonVoxels.cxx
	SimChannel.cxx
	SimListUtils.cxx
	SimPhotons.cxx
	)

target_link_libraries(Simulation
     larsoft::Utilities
     ${SIMULATIONBASE}
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     ${ROOT_BASIC_LIB_LIST}
     )

art_add_dictionary(DICTIONARY_LIBRARIES art::art_Framework_Core)

art_add_service(LArVoxelCalculator_service SHARED
	LArVoxelCalculator_service.cc
	)

art_add_service(LArG4Parameters_service SHARED
	LArG4Parameters_service.cc
	)

install(FILES ${HEADERS} DESTINATION
     ${CMAKE_INSTALL_INCLUDEDIR}/Simulation COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     Simulation
     Simulation_map
     Simulation_dict
     LArVoxelCalculator_service 
     LArG4Parameters_service
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )
