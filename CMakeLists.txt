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


# use cmake 3.12 or later
cmake_minimum_required (VERSION 3.12)

project(larsim VERSION 8.19.03 LANGUAGES CXX)

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
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable -Wno-shadow -Wno-unused-parameter -Wno-implicit-fallthrough -Wno-ignored-qualifiers
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( nurandom REQUIRED )
find_package( larevt REQUIRED )
find_package( Geant4 REQUIRED )
include(${Geant4_USE_FILE})
find_package( ROOT REQUIRED )
find_package( CLHEP REQUIRED )
find_package( marley )
find_package( larsoft-data )

find_library( CRY NAMES CRY )
find_library( SQLITE3 NAMES sqlite3 )

# genie (copied from nugen)
cet_find_library( LOG4CPP NAMES log4cpp PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWMSG       NAMES GFwMsg       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWREG       NAMES GFwReg       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWALG       NAMES GFwAlg       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWINT       NAMES GFwInt       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWGHEP      NAMES GFwGHEP      PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWNUM       NAMES GFwNum       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWUTL       NAMES GFwUtl       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWPARDAT    NAMES GFwParDat    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWEG        NAMES GFwEG        PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GFWNTP       NAMES GFwNtp       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHXSIG      NAMES GPhXSIg      PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHPDF       NAMES GPhPDF       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHNUCLST    NAMES GPhNuclSt    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHCMN       NAMES GPhCmn       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDCY       NAMES GPhDcy       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHHADTRANSP NAMES GPhHadTransp PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHHADNZ     NAMES GPhHadnz     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDEEX      NAMES GPhDeEx      PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHAMNGXS    NAMES GPhAMNGXS    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHAMNGEG    NAMES GPhAMNGEG    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHCHMXS     NAMES GPhChmXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHCOHXS     NAMES GPhCohXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHCOHEG     NAMES GPhCohEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDISXS     NAMES GPhDISXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDISEG     NAMES GPhDISEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDFRCXS    NAMES GPhDfrcXS    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHDFRCEG    NAMES GPhDfrcEG    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHGLWRESXS  NAMES GPhGlwResXS  PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHGLWRESEG  NAMES GPhGlwResEG  PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHIBDXS     NAMES GPhIBDXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHIBDEG     NAMES GPhIBDEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHMNUCXS    NAMES GPhMNucXS    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHMNUCEG    NAMES GPhMNucEG    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHMEL       NAMES GPhMEL       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHNUELXS    NAMES GPhNuElXS    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHNUELEG    NAMES GPhNuElEG    PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHQELXS     NAMES GPhQELXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHQELEG     NAMES GPhQELEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHRESXS     NAMES GPhResXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHRESEG     NAMES GPhResEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHSTRXS     NAMES GPhStrXS     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHSTREG     NAMES GPhStrEG     PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GPHNDCY      NAMES GPhNDcy      PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GTLGEO       NAMES GTlGeo       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GTLFLX       NAMES GTlFlx       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GRWFWK       NAMES GRwFwk       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GRWIO        NAMES GRwIO        PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )
cet_find_library( GRWCLC       NAMES GRwClc       PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )

# GENIE pre-R-3 used GSL via ROOT's Math libs, now makes direct call not available via that route
cet_find_library( GSL NAMES gsl PATHS ${CMAKE_INSTALL_RPATH} NO_DEFAULT_PATH )

set (GENIE_LIB_LIST ${GFWMSG}
                    ${GFWREG}
                    ${GFWALG}
                    ${GFWINT}
                    ${GFWGHEP}
                    ${GFWNUM}
                    ${GSL}  # FWNUM relies on GSL
                    ${GFWUTL}
                    ${GFWPARDAT}
                    ${GFWEG}
                    ${GFWNTP}
                    ${GPHXSIG}
                    ${GPHPDF}
                    ${GPHNUCLST}
                    ${GPHCMN}
                    ${GPHDCY}
                    ${GPHHADTRANSP}
                    ${GPHHADNZ}
                    ${GPHDEEX}
                    ${GPHAMNGXS}
                    ${GPHAMNGEG}
                    ${GPHCHMXS}
                    ${GPHCOHXS}
                    ${GPHCOHEG}
                    ${GPHDISXS}
                    ${GPHDISEG}
                    ${GPHDFRCXS}
                    ${GPHDFRCEG}
                    ${GPHGLWRESXS}
                    ${GPHGLWRESEG}
                    ${GPHIBDXS}
                    ${GPHIBDEG}
                    ${GPHMNUCXS}
                    ${GPHMNUCEG}
                    ${GPHMEL}
                    ${GPHNUELXS}
                    ${GPHNUELEG}
                    ${GPHQELXS}
                    ${GPHQELEG}
                    ${GPHRESXS}
                    ${GPHRESEG}
                    ${GPHSTRXS}
                    ${GPHSTREG}
                    ${GPHNDCY}
                    ${GTLGEO}
                    ${GTLFLX}
                    ${GRWFWK}
                    ${GRWIO}
                    ${GRWCLC}
                   ${LOG4CPP}
                    )


# marley libraries
cet_find_library( MARLEY NAMES MARLEY)
cet_find_library( MARLEY_ROOT NAMES MARLEY_ROOT)

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
