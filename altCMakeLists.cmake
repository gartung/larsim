cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0025)
	cmake_policy(SET CMP0025 OLD)
endif()

if(POLICY CMP0042)
	cmake_policy(SET CMP0042 NEW)
endif()


project(larsim)


set(larsim_VERSION "04.00.01")
set(larsim_VERSION_MAJOR 04)
set(larsim_VERSION_MINOR 00)
set(larsim_VERSION_PATCH 01)

set(larsim_SOVERSION "1.0.0")


set(larsim_DEBUG_POSTFIX "d")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Modules ${CMAKE_MODULE_PATH})

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(BASE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/BuildProducts")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")

find_package(nutools 1.07.00 REQUIRED)

find_package(art 1.11.3 REQUIRED)

find_package(FNALCore 0.1.0 REQUIRED)

find_package( larcore 1.0.0 )

find_package( lardata 1.0.0 )

find_package( larevt 1.0.0 )

find_package( Geant4 9.6.3 REQUIRED )

find_package( Boost 1.53.0  REQUIRED
      serialization
      date_time
     )

find_package(GCCXML 0.9.0 REQUIRED)

find_package(CLHEP 2.2.0.3 REQUIRED)

find_package(SQLite3 3.8.5 REQUIRED)

find_package(ROOT 5.34.20 REQUIRED
     Core
     Cint
     Cintex
     Hist
     Matrix
     Reflex
     RIO
     Thread
     Tree
)

if(NOT ROOT_python_FOUND)
     message(FATAL_ERROR "art requires ROOT with Python support")
endif()

find_package(TBB 4.1.0 REQUIRED)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FNALCore_CXX_FLAGS} -O3 -g -DNDEBUG -fno-omit-frame-pointer -pedantic -Wno-unused-local-typedefs")

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})
include_directories(${art_INCLUDE_DIRS})
include_directories(${FNALCore_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${larcore_INCLUDE_DIRS})
include_directories(${lardata_INCLUDE_DIRS})
include_directories(${larevt_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${Geant4_INCLUDE_DIRS})
include_directories( ${nutools_INCLUDE_DIRS} )
include_directories( ${PostgreSQL_INCLUDE_DIRS})
include_directories( ${GEANT4_FQ_DIR}/include )

set( PQ ${PostgreSQL_LIBRARIES})
find_library( SIMULATIONBASE NAMES SimulationBase PATHS ${NUTOOLS_LIB} NO_DEFAULT_PATH )
#find_library( EVENTGENERATORBASECRY NAMES EventGeneratorBaseCRY PATHS ${NUTOOLS_LIB} NO_DEFAULT_PATH )
#find_library( EVENTGENERATORBASEGENIE NAMES EventGeneratorBaseGENIE PATHS ${NUTOOLS_LIB} NO_DEFAULT_PATH )

find_library( CRY NAMES CRY PATHS ${CRY_HOME}/lib NO_DEFAULT_PATH )
# genie
find_library( GALGORITHM NAMES GAlgorithm PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GBARYONRESONANCE NAMES GBaryonResonance PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GBASE NAMES GBase PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GBODEKYANG NAMES GBodekYang PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GCHARM NAMES GCharm PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GCOH NAMES GCoh PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GDFRC NAMES GDfrc PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GDIS NAMES GDIS PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GCROSSSECTIONS NAMES GCrossSections PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GDECAY NAMES GDecay PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GELAS NAMES GElas PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GELFF NAMES GElFF PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GHEP NAMES GHEP PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GEVGCORE NAMES GEVGCore  PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GEVGMODULES NAMES GEVGModules PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GEVGDRIVERS NAMES GEVGDrivers PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GGIBUU NAMES GGiBUU PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GHADRONTRANSP NAMES GHadronTransp PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GFRAGMENTATION NAMES GFragmentation PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GINTERACTION NAMES GInteraction PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GLLEWELLYNSMITH NAMES GLlewellynSmith  PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GMEC NAMES GMEC PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GMESSENGER NAMES GMessenger PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GNUGAMMA NAMES GNuGamma PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GNUE NAMES GNuE PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GNTUPLE NAMES GNtuple PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GNUCLEAR NAMES GNuclear PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GNUMERICAL NAMES GNumerical PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GQPM NAMES GQPM PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GPDG NAMES GPDG PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GPDF NAMES GPDF PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GQEL NAMES GQEL PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GRES NAMES GRES PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GREGISTRY NAMES GRegistry PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GREINSEGHAL NAMES GReinSeghal PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GUTILS NAMES GUtils PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GGEO NAMES GGeo PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GFLUXDRIVERS NAMES GFluxDrivers PATHS $(GENIE_LIB) NO_DEFAULT_PATH )
find_library( GMUELOSS NAMES GMuELoss PATHS $(GENIE_LIB) NO_DEFAULT_PATH )

add_subdirectory(DetSim)
add_subdirectory(EventGenerator)
add_subdirectory(LArG4)
add_subdirectory(MCCheater)
add_subdirectory(PhotonPropagation)
add_subdirectory(Simulation)
add_subdirectory(TriggerAlgo)
add_subdirectory(SimFilters)
 
configure_package_config_file(
  Modules/larsimConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/larsimConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/larsim-${larsim_VERSION}
  PATH_VARS
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_LIBDIR
  )

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/larsimConfigVersion.cmake
  VERSION ${larsim_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/larsimConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/larsimConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/larsim-${larsim_VERSION}
  COMPONENT Development
  )

install(EXPORT larsimLibraries
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/larsim-${larsim_VERSION}
  NAMESPACE larsoft::
  COMPONENT Development
  )



include(ArtCPack)
