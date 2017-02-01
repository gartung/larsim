////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.cxx
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.cxx,v 1.4 2009/08/12 20:52:14 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
///
///
///  This version has undergone changes from the original LArVoxelReadout to
///  remove post G4 energy deposition steps. Ionization/Scintillation,
///  drifitng, and other elements are now separated.
///
///  Wes, February 2017
///
////////////////////////////////////////////////////////////////////////

// C/C++ standard library
#include <cstdio> // std::sscanf()
#include <cmath> // std::ceil()
#include <string>
#include <map>
#include <utility> // std::move()
#include <cassert>

// GEANT
#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4TouchableHistory.hh"
#include "Geant4/G4TouchableHandle.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"

// framework libraries
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft code
#include "larsim/LArG4/LArVoxelReadout.h"
#include "larsim/LArG4/ParticleListAction.h"

// CLHEP
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"

namespace larg4 {
  
  
  //---------------------------------------------------------------------------------------
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout(std::string const& name)
    : G4VSensitiveDetector(name)
  {
  } // LArVoxelReadout::LArVoxelReadout()

  //---------------------------------------------------------------------------------------
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout
  (std::string const& name, unsigned int, unsigned int)
    : LArVoxelReadout(name)
  { }
    
  //---------------------------------------------------------------------------------------
  LArVoxelReadout::~LArVoxelReadout() {}

  //---------------------------------------------------------------------------------------
  // Called at the start of each event.
  void LArVoxelReadout::Initialize(G4HCofThisEvent*)
  {
    fSimEDepCol.clear();
    fSimEDepCol.reserve(5000000);
  }  

  //---------------------------------------------------------------------------------------
  // Called at the end of each event.
  void LArVoxelReadout::EndOfEvent(G4HCofThisEvent*)
  {
    fSimEDepCol.clear();
  }// LArVoxelReadout::EndOfEvent()

  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::clear()
  {
  }

  
  //---------------------------------------------------------------------------------------
  //This is just a function for initializing the internals once we set the step we are looking at.
  void LArVoxelReadout::ProcessStep( G4Step* step)
  {
    if(step->GetTotalEnergyDeposit() <= 0) return;

    fStepMidPoint = 0.5*( step->GetPreStepPoint()->GetPosition() + 
			  step->GetPostStepPoint()->GetPosition() );
    
    fSimEDepCol.emplace_back(fStepMidPoint.x(),fStepMidPoint.y(),fStepMidPoint.z(),
			     step->GetPreStepPoint()->GetGlobalTime(),
			     step->GetTotalEnergyDeposit()/CLHEP::MeV,
			     ParticleListAction::GetCurrentTrackID());    
  }

  
  //---------------------------------------------------------------------------------------
  // Called for each step.
  G4bool LArVoxelReadout::ProcessHits( G4Step* step, G4TouchableHistory* pHistory)
  {
    // All work done for the "parallel world" "box of voxels" in
    // LArVoxelReadoutGeometry makes this a fairly simple routine.
    // First, the usual check for non-zero energy:

    // Only process the hit if the step is inside the active volume and
    // it has deposited energy.  The hit being inside the active volume
    // is virtually sure to happen because the LArVoxelReadoutGeometry
    // that this class makes use of only has voxels for inside the TPC.

    // The step can be no bigger than the size of the voxel,
    // because of the geometry set up in LArVoxelGeometry and the
    // transportation set up in PhysicsList.  Find the mid-point
    // of the step.

    ProcessStep(step);

    return true;
  }

  //---------------------------------------------------------------------------------------
  // Never used but still have to be defined for G4
  void LArVoxelReadout::DrawAll()  {}
  void LArVoxelReadout::PrintAll() {}

} // namespace larg4
