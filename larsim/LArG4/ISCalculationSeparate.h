////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.h
/// \brief Interface to algorithm class for a specific calculation of 
///        ionization electrons and scintillation photons assuming there
///        is no correlation between the two
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATIONSEPARATE_H
#define LARG4_ISCALCULATIONSEPARATE_H

#include <map>

#include "Geant4/G4EmSaturation.hh"

#include "larsim/LArG4/ISCalculation.h"

// forward declaration
namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

 class ISCalculationSeparate : public ISCalculation {

 public:

   ISCalculationSeparate(CLHEP::HepRandomEngine&);
   virtual ~ISCalculationSeparate();

   void   Initialize();
   void   Reset();
   void   CalculateIonizationAndScintillation(const G4Step* step);
   double StepSizeLimit()              const { return fStepSize;            }

  CLHEP::HepRandomEngine& fEngine;     ///< random engine
  int BinomFluct ( int N0, double prob );

 private:

   double                fStepSize;            ///< maximum step to take				  
   double                fEfield;              ///< value of electric field from LArProperties service
   double 	   	 fGeVToElectrons;      ///< conversion factor from LArProperties service	  
   double 	   	 fRecombA;             ///< from LArG4Parameters service			  
   double 	   	 fRecombk;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxA;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxB;             ///< from LArG4Parameters service			  
   bool   	   	 fUseModBoxRecomb;     ///< from LArG4Parameters service			  
   bool   	   	 fScintByParticleType; ///< from LArProperties service			  
   double 	   	 fScintYieldFactor;    ///< scintillation yield factor                             
   G4EmSaturation* 	 fEMSaturation;        ///< pointer to EM saturation                            

   double                fEfieldNominal;
   int                   fCachedTrackID;  
     
 };
}
#endif // LARG4_ISCALCULATIONSEPARATE_H

