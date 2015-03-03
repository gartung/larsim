include_directories ( ${GENIE_INC}/GENIE )
include_directories ( ${LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( ${GEANT4_FQ_DIR}/include )
include_directories ( ${XERCES_C_INC} )
include_directories ( ${CRYHOME}/src )

set(HEADERS
     AuxDetReadout.h
     AuxDetReadoutGeometry.h
     FastOpticalPhysics.h
     G4BadIdeaAction.h
     ISCalculation.h
     ISCalculationNEST.h
     ISCalculationSeparate.h
     IonizationAndScintillation.h
     IonizationAndScintillationAction.h
     LArStackingAction.h
     LArVoxelReadout.h
     LArVoxelReadoutGeometry.h
     MaterialPropertyLoader.h
     MuNuclearSplittingProcess.h
     MuNuclearSplittingProcessXSecBias.h
     NestAlg.h
     OpDetLookup.h
     OpDetPhotonTable.h
     OpDetReadoutGeometry.h
     OpDetSensitiveDetector.h
     OpParamAction.h
     OpParamSD.h
     ParticleListAction.h
     PhysicsList.h
     VisualizationAction.h
     ConfigurablePhysicsList.hh
     CustomPhysicsBuiltIns.hh
     CustomPhysicsFactory.hh
     CustomPhysicsTable.hh
     G4ThermalElectron.hh
     OpBoundaryProcessSimple.hh
     OpFastScintillation.hh
     OpticalPhysics.hh
     ConfigurablePhysicsList.icc
     )

add_library(LArG4 SHARED
     ${HEADERS}
     AuxDetReadout.cxx
     AuxDetReadoutGeometry.cxx
     CustomPhysicsBuiltIns.cxx
     CustomPhysicsTable.cxx
     FastOpticalPhysics.cxx
     G4BadIdeaAction.cxx
     G4ThermalElectron.cxx
     ISCalculation.cxx
     ISCalculationNEST.cxx
     ISCalculationSeparate.cxx
     IonizationAndScintillation.cxx
     IonizationAndScintillationAction.cxx
     LArStackingAction.cxx
     LArVoxelReadout.cxx
     LArVoxelReadoutGeometry.cxx
     MaterialPropertyLoader.cxx
     MuNuclearSplittingProcess.cxx
     MuNuclearSplittingProcessXSecBias.cxx
     NestAlg.cxx
     OpBoundaryProcessSimple.cxx
     OpDetLookup.cxx
     OpDetPhotonTable.cxx
     OpDetReadoutGeometry.cxx
     OpDetSensitiveDetector.cxx
     OpFastScintillation.cxx
     OpParamAction.cxx
     OpParamSD.cxx
     OpticalPhysics.cxx
     ParticleListAction.cxx
     PhysicsList.cxx
     VisualizationAction.cxx
     )

target_link_libraries(LArG4
     PhotonPropagation
     Simulation
     larsoft::Utilities
     larsoft::Geometry
     FNALCore::FNALCore
     ${SIMULATIONBASE}
     ${ROOT_EG}
     ${ROOT_TREEPLAYER} 
     ${ROOT_FFTW}
     ${ROOT_REFLEX}
     ${ROOT_BASIC_LIB_LIST}
     ${G4_LIB_LIST}
     ${ROOT_BASIC_LIB_LIST}
     )

art_add_module(LArG4Ana_module
     LArG4Ana_module.cc
     )


art_add_module(LArG4_module
     LArG4_module.cc
     )

art_add_module(LArSimChannelAna_module
     LArSimChannelAna_module.cc
     )

install(FILES ${HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/LArG4 COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


install(TARGETS
     LArG4
     LArG4_module
     LArG4Ana_module
     LArSimChannelAna_module
     EXPORT larsimLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

set(LArG4_MAC_FILES
	LArG4.mac
	atree.mac
	)	

install(FILES ${LArG4_MAC_FILES} DESTINATION G4 COMPONENT Runtime)
