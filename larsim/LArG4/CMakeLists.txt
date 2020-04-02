include(${Geant4_USE_FILE})
include_directories ( ${GENIE_INC}/GENIE )
include_directories ( ${LIBXML2_INC})
include_directories ( ${GEANT4_INCLUDE_DIRS} )
include_directories ( ${XERCES_C_INC} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS} -Wno-shadow -Wno-unused-parameter")

art_make(LIB_LIBRARIES
           larsim_PhotonPropagation_PhotonVisibilityService_service
           lardataobj_Simulation
           larcorealg_Geometry
           nug4_G4Base
           nug4_ParticleNavigation
           nurandom_RandomUtils_NuRandomService_service
           ${ART_UTILITIES}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ROOT::Geom
           ROOT::MathCore
           ${CLHEP_LIBRARIES}
           G4digits_hits
           G4event
           G4geometry
           G4global
           G4graphics_reps
           G4materials
           G4intercoms
           G4particles
           G4physicslists
           G4processes
           G4track
           G4run
         MODULE_LIBRARIES
           larsim_LArG4
           larsim_MCCheater_ParticleInventoryService_service
           lardataalg_MCDumpers
           larcorealg_Geometry
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ROOT::Hist
           ROOT::Core
           ROOT::Tree
         )

install_headers()
install_fhicl()

# install macro files in lib
file(GLOB g4_macros *.mac )
install_fw( LIST ${g4_macros} )
