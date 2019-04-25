# ======================================================================
#  larsim main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsim/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.11 or later
cmake_minimum_required (VERSION 3.11)

project(larsim VERSION 8.03.03 LANGUAGES CXX)

set(fcl_dir job)
set(gdml_dir gdml)
set(fw_dir fw)
find_package( cetmodules )
include(CetCMakeEnv)
cet_cmake_env()
include(ParseVersion)

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)
add_definitions(-DGENIE_PRE_R3)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( larevt REQUIRED )
find_package( marley )
find_package( larsoft-data )

find_library( CRY NAMES CRY )
find_library( SQLITE3 NAMES sqlite3 )
# genie
cet_find_library( GALGORITHM NAMES GAlgorithm )
cet_find_library( GBARYONRESONANCE NAMES GBaryonResonance )
cet_find_library( GBASE NAMES GBase )
cet_find_library( GBODEKYANG NAMES GBodekYang )
cet_find_library( GCHARM NAMES GCharm )
cet_find_library( GCOH NAMES GCoh )
cet_find_library( GDFRC NAMES GDfrc )
cet_find_library( GDIS NAMES GDIS )
cet_find_library( GCROSSSECTIONS NAMES GCrossSections )
cet_find_library( GDECAY NAMES GDecay )
cet_find_library( GELAS NAMES GElas )
cet_find_library( GELFF NAMES GElFF )
cet_find_library( GHEP NAMES GHEP )
cet_find_library( GEVGCORE NAMES GEVGCore  )
cet_find_library( GEVGMODULES NAMES GEVGModules )
cet_find_library( GEVGDRIVERS NAMES GEVGDrivers )
cet_find_library( GGIBUU NAMES GGiBUU )
cet_find_library( GHADRONTRANSP NAMES GHadronTransp )
cet_find_library( GFRAGMENTATION NAMES GFragmentation )
cet_find_library( GINTERACTION NAMES GInteraction )
cet_find_library( GLLEWELLYNSMITH NAMES GLlewellynSmith  )
cet_find_library( GMEC NAMES GMEC )
cet_find_library( GMESSENGER NAMES GMessenger )
cet_find_library( GNUGAMMA NAMES GNuGamma )
cet_find_library( GNUE NAMES GNuE )
cet_find_library( GNTUPLE NAMES GNtuple )
cet_find_library( GNUCLEAR NAMES GNuclear )
cet_find_library( GNUMERICAL NAMES GNumerical )
cet_find_library( GQPM NAMES GQPM )
cet_find_library( GPDG NAMES GPDG )
cet_find_library( GPDF NAMES GPDF )
cet_find_library( GQEL NAMES GQEL )
cet_find_library( GRES NAMES GRES )
cet_find_library( GREGISTRY NAMES GRegistry )
cet_find_library( GREINSEGHAL NAMES GReinSeghal )
cet_find_library( GUTILS NAMES GUtils )
cet_find_library( GGEO NAMES GGeo )
cet_find_library( GFLUXDRIVERS NAMES GFluxDrivers )
cet_find_library( GMUELOSS NAMES GMuELoss )
# marley libraries
cet_find_library( MARLEY NAMES MARLEY )
cet_find_library( MARLEY_ROOT NAMES MARLEY_ROOT )

# temporarily needed since this is an unexpected header path
include_directories( ${IFDHC_DIR}/inc )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(larsim)

# tests
add_subdirectory(test)
add_subdirectory(cmakeConfig)
# ups - table and config files
#add_subdirectory(ups)

# packaging utility
#include(UseCPack)
