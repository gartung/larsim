include_directories ( ${GEANT4_INCLUDE_DIRS} )

art_make_exec(NAME simplifyGDML
              LIBRARIES
              ${CLHEP_LIBRARIES}
	      ${Geant4_LIBRARIES}
             )
