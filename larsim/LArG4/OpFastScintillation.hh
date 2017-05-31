// Class adapted for LArSoft by Ben Jones, MIT 10/10/12
//
// This class is a physics process based on the standard Geant4
// scintillation process.
//
// It has been stripped down and adapted to form the backbone of
// the LArG4 fast optical simulation.  Photons, instead of being
// produced and added to the geant4 particle stack, are logged
// and used to predict the visibility of this step to each PMT in
// the detector.
//
// The photonvisibilityservice looks up the visibility of the relevant
// xyz point, and if a photon is detected at a given PMT, one OnePhoton
// object is logged in the OpDetPhotonTable
//
// At the end of the event, the OpDetPhotonTable is read out
// by LArG4, and detected photons are stored in the event.
//
// This process can be used alongside the standard Cerenkov process,
// which does step geant4 opticalphotons.  Both the fast scintillation
// table and the geant4 sensitive detectors are read out by LArG4 to
// produce a combined SimPhoton collection.
//
//
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: OpFastScintillation.hh,v 1.21 2010-10-28 23:29:21 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        OpFastScintillation.hh
// Description:	Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//                         of energy deposited by particle type
//                         Thanks to Zach Hartwig (Department of Nuclear
//                         Science and Engineeering - MIT)
//              2005-07-28 add G4ProcessType to constructor
//              2002-11-21 change to user G4Poisson for small MeanNumPotons
//              2002-11-07 allow for fast and slow scintillation
//              2002-11-05 make use of constant material properties
//              2002-05-16 changed to inherit from VRestDiscreteProcess
//              2002-05-09 changed IsApplicable method
//              1999-10-29 add method and class descriptors
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef OpFastScintillation_h
#define OpFastScintillation_h 1

/////////////
// Includes
/////////////

#include "Geant4/globals.hh"
#include "Geant4/templates.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleMomentum.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4VRestDiscreteProcess.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4PhysicsTable.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4PhysicsOrderedFreeVector.hh"
#include "Geant4/G4EmSaturation.hh"

#include "TVector3.h"

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

namespace larg4{

class OpFastScintillation : public G4VRestDiscreteProcess
{
public:
  OpFastScintillation(const G4String& processName = "Scintillation",
                      G4ProcessType type = fElectromagnetic);

  OpFastScintillation(const OpFastScintillation &right);

  ~OpFastScintillation();

public:
  // OpFastScintillation Process has both PostStepDoIt (for energy
  // deposition of particles in flight) and AtRestDoIt (for energy
  // given to the medium by particles at rest)

  // Returns true -> 'is applicable', for any particle type except
  // for an 'opticalphoton' and for short-lived particles
  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.
  G4double GetMeanFreePath(const G4Track& aTrack,
                           G4double ,
                           G4ForceCondition* );

  // Returns infinity; i. e. the process does not limit the time,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.
  G4double GetMeanLifeTime(const G4Track& aTrack,
                           G4ForceCondition* );

  // These are the methods implementing the scintillation process.
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step&  aStep);
  virtual G4VParticleChange* AtRestDoIt (const G4Track& aTrack,
                                         const G4Step& aStep);

  // If set, the primary particle tracking is interrupted and any produced
  // scintillation photons are tracked next. When all have been tracked, the
  // tracking of the primary resumes.
  void SetTrackSecondariesFirst(const G4bool state)
  {
    fTrackSecondariesFirst = state;
  }

  // If set, the OpFastScintillation process expects the user to have
  // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.
  void SetFiniteRiseTime(const G4bool state)
  {
    fFiniteRiseTime = state;
  }

  // Returns the boolean flag for tracking secondaries first.
  G4bool GetTrackSecondariesFirst() const {return fTrackSecondariesFirst;}

  // Returns the boolean flag for a finite scintillation rise time.
  G4bool GetFiniteRiseTime() const {return fFiniteRiseTime;}

  // Called to set the scintillation exciation ratio, needed when
  // the scintillation level excitation is different for different
  // types of particles. This overwrites the YieldRatio obtained
  // from the G4MaterialPropertiesTable.
  void SetScintillationExcitationRatio(const G4double excitationratio);

  G4double GetScintillationExcitationRatio() const {return fExcitationRatio;}

  G4PhysicsTable* GetFastIntegralTable() const {return fFastIntegralTable;}
  G4PhysicsTable* GetSlowIntegralTable() const {return fSlowIntegralTable;}

  // Adds Birks Saturation to the process.
  void AddSaturation(G4EmSaturation* sat) { fEMSaturation = sat; }

  // Removes the Birks Saturation from the process.
  void RemoveSaturation() { fEMSaturation = NULL; }

  // Returns the Birks Saturation.
  G4EmSaturation* GetSaturation() const { return fEMSaturation; }

  // Called by the user to set the scintillation yield as a function
  // of energy deposited by particle type
  void SetScintillationByParticleType(const G4bool );

  // Return the boolean that determines the method of scintillation
  // production
  G4bool GetScintillationByParticleType() const { return fScintillationByParticleType; }

  // Prints the fast and slow scintillation integral tables.
  void DumpPhysicsTable() const;

protected:
  // It builds either the fast or slow scintillation integral table;
  // or both.
  void BuildPhysicsTable();

  // Note the production of N photons in at point xyz.
  // pass on to generate detector response, etc.
  bool RecordPhotonsProduced(const G4Step& aStep, double N);

  double GetYieldRatio(const G4DynamicParticle* aParticle,
                       G4MaterialPropertiesTable* aMaterialPropertiesTable) const;

  void PropogatePhoton(TVector3 r0);

  G4PhysicsTable* fSlowIntegralTable;
  G4PhysicsTable* fFastIntegralTable;

  G4bool fTrackSecondariesFirst;
  G4bool fFiniteRiseTime;

  G4double fExcitationRatio;

  G4bool fScintillationByParticleType;

  G4double single_exp(G4double t, G4double tau2);
  G4double bi_exp(G4double t, G4double tau1, G4double tau2);

  // emission time distribution when there is a finite rise time
  G4double sample_time(G4double tau1, G4double tau2);

  static TVector3 random_unit();

  static void ortho_basis(TVector3 p0, TVector3& p1, TVector3& p2);

  G4EmSaturation* fEMSaturation;

};

inline
G4bool OpFastScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if (aParticleType.GetParticleName() == "opticalphoton") return false;
  if (aParticleType.IsShortLived()) return false;

  return true;
}

inline
void OpFastScintillation::DumpPhysicsTable() const
{
  if (fFastIntegralTable) {
    G4int PhysicsTableSize = fFastIntegralTable->entries();
    G4PhysicsOrderedFreeVector *v;

    for (G4int i = 0 ; i < PhysicsTableSize ; i++ ){
      v = (G4PhysicsOrderedFreeVector*)(*fFastIntegralTable)[i];
      v->DumpValues();
    }
  }

  if (fSlowIntegralTable) {
    G4int PhysicsTableSize = fSlowIntegralTable->entries();
    G4PhysicsOrderedFreeVector *v;

    for (G4int i = 0 ; i < PhysicsTableSize ; i++ ){
      v = (G4PhysicsOrderedFreeVector*)(*fSlowIntegralTable)[i];
      v->DumpValues();
    }
  }
}

inline
G4double OpFastScintillation::single_exp(G4double t, G4double tau2)
{
  return std::exp(-1.0*t/tau2)/tau2;
}

inline
G4double OpFastScintillation::bi_exp(G4double t, G4double tau1, G4double tau2)
{
  return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
}

} //namespace

#endif /* OpFastScintillation_h */
