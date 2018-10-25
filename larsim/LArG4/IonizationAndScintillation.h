////////////////////////////////////////////////////////////////////////
/// \file  IonizationAndScintillation.h
/// \brief Singleton to access a unified treatment of ionization and 
///        scintillation in LAr
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_IONIZATIONANDSCINTILLATION_H
#define LARG4_IONIZATIONANDSCINTILLATION_H

#include <cstring>

#include "larsim/LArG4/ISCalculation.h"

#include "Geant4/G4Step.hh"

#include "TH1.h"
#include "TH2.h"

namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

  // The Ionization and Scintillation singleton
  class IonizationAndScintillation
  {
  public:

    static IonizationAndScintillation* CreateInstance(CLHEP::HepRandomEngine& engine);
    static IonizationAndScintillation* Instance();

    // Method to reset the internal variables held in the ISCalculation
    // This method should be called at the start of any G4Step
    void Reset(const G4Step* step); 

    double EnergyDeposit()              const { return fISCalc->EnergyDeposit();              } 
    double NumberIonizationElectrons()  const { return fISCalc->NumberIonizationElectrons();  } 
    double NumberIons()                 const { return fISCalc->NumberIons();  } 
    double NumberScintillationPhotons() const { return fISCalc->NumberScintillationPhotons(); } 
    double StepSizeLimit()              const { return fISCalc->StepSizeLimit();              }
    double dEdx()                       const { return fISCalc->dEdx();                       }

  private:

    IonizationAndScintillation(CLHEP::HepRandomEngine& engine);
    ~IonizationAndScintillation();

    larg4::ISCalculation* fISCalc;             ///< object to calculate ionization and scintillation
                                               ///< produced by an energy deposition
    std::string           fISCalculator;       ///< name of calculator to use, NEST or Separate
    G4Step const*         fStep;               ///< pointer to the current G4 step
    int                   fStepNumber;         ///< last StepNumber checked
    int                   fTrkID;              ///< last TrkID checked
    double                fNumIon;

    TH1F*                 fElectronsPerStep;   ///< histogram of electrons per step
    TH1F*                 fStepSize;           ///< histogram of the step sizes
    TH2F*                 fStepSizeCompare;   
    TH1F*                 fPhotonsPerStep;     ///< histogram of the photons per step
    TH1F*                 fEnergyPerStep;      ///< histogram of the energy deposited per step
    TH1F*                 fElectronsPerLength; ///< histogram of electrons per cm
    TH1F*                 fPhotonsPerLength;   ///< histogram of photons per cm
    TH1F*                 fElectronsPerEDep;   ///< histogram of electrons per MeV deposited
    TH1F*                 fPhotonsPerEDep;     ///< histogram of photons per MeV deposited
    TH2F*                 fElectronsVsPhotons; ///< histogram of electrons vs photons per step 
    TH1F*                 fRecombSurvivalFrac; ///< histogram of recombination factor
    TH1F*                 fRecombSurvivalFrac_el; ///< histogram of recombination factor
    TH1F*                 fEnergyPerLength;
    TH1F*                 fEnergyPerLength_el;
    TH1F*                 fEnergyPerLength_mu;
    TH2F*                 fdEdxVsRecomb;
    TH2F*                 fdEdxVsStepLength;
    TH1F*                 fRecomb_dEdxBins_2;
    TH1F*                 fRecomb_dEdxBins_5;
    TH1F*                 fRecomb_dEdxBins_15;
    TH1F*                 fRecomb_dEdxBins_30;
    TH1F*                 fdQdx_dEdxBins_2;
    TH1F*                 fdQdx_dEdxBins_5;
    TH1F*                 fdQdx_dEdxBins_15;
    TH1F*                 fdQdx_dEdxBins_30;
    CLHEP::HepRandomEngine& fEngine;           ///< random engine
  };

} // namespace larg4


#endif // LARG4_IONIZATIONANDSCINTILLATION
