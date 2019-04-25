art_make(MODULE_LIBRARIES
         lardataalg_MCDumpers
         nusimdata_SimulationBase
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${ROOT_EG_LIBRARY}
         ${ROOT_Core_LIBRARY}
         ${ROOT_Physics_LIBRARY})

install_headers()
install_fhicl()
install_source()
