
include_directories ( ${GENIE_INC}/GENIE )

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
#cet_find_library( GREINSEGHAL NAMES GReinSeghal )
cet_find_library( GUTILS NAMES GUtils )
cet_find_library( GGEO NAMES GGeo )
cet_find_library( GFLUXDRIVERS NAMES GFluxDrivers )
cet_find_library( GMUELOSS NAMES GMuELoss )
cet_find_library( GREWEIGHT NAMES GReWeight )
cet_find_library( GNUCLEONDECAY NAMES GNucleonDecay )
cet_find_library( GNEUTRONOSC NAMES GNeutronOsc )
# prior to GENIE R-2_9_0 Sehgal's name was spelled incorrectly
# so if minor version # is 8 or less use GReinSeghal, otherwise GReinSehgal
# parse_ups_version returns VMAJ, VMIN, and VPRJ
parse_underscored_version( ${GENIE_VERSION} )
set ( GRSSPELLING GReinSehgal )
if ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
  set ( GRSSPELLING GReinSeghal )
endif ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
message(STATUS "GENIE version ${GENIE_VERSION} major ${VMAJ} minor ${VMIN} ${GRSSPELLING}")
cet_find_library( GREINSEHGAL NAMES ${GRSSPELLING} )

art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         larsim_Simulation lardataobj_Simulation
         nusimdata_SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nutools_RandomUtils_NuRandomService_service
         nutools_EventGeneratorBase_GENIE
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${CLHEP_LIBRARIES}
#	 ${GMESSENGER}
         ${GALGORITHM}
         ${GBARYONRESONANCE}
         ${GBASE}
         ${GBODEKYANG}
         ${GCHARM}
         ${GCOH}
         ${GDFRC}
         ${GDIS}
         ${GCROSSSECTIONS}
         ${GDECAY}
         ${GELAS}
         ${GELFF}
         ${GHEP}
         ${GEVGCORE}
         ${GEVGMODULES}
         ${GEVGDRIVERS}
         ${GGIBUU}
         ${GHADRONTRANSP}
         ${GFRAGMENTATION}
         ${GINTERACTION}
         ${GLLEWELLYNSMITH}
         ${GMEC}
         ${GNUGAMMA}
         ${GNUE}
         ${GNTUPLE}
         ${GNUCLEAR}
         ${GNUMERICAL}
         ${GQPM}
         ${GPDG}
         ${GPDF}
         ${GQEL}
         ${GRES}
         ${GREGISTRY}
         ${GREINSEHGAL}
         ${GUTILS}
         ${GGEO}
         ${GFLUXDRIVERS}
         ${GMUELOSS}
         ${GREWEIGHT}
         ${GNUCLEONDECAY}
         ${GNEUTRONOSC}
         ${GMESSENGER}
         ${ROOT_EGPythia6_LIBRARY}    # FIXME!!! - resolving genie run time reference
         ${ROOT_EG_LIBRARY}
         ${ROOT_Hist_LIBRARY}
         ${ROOT_MathCore_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
