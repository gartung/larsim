include_directories ( ${GEANT4_INCLUDE_DIRS} )
include_directories ( ${XERCES_C_INC} )

art_make_exec(NAME simplifyGDML
  LIBRARIES
    ${CLHEP_LIBRARIES}
    ${Geant4_LIBRARIES}
  )

