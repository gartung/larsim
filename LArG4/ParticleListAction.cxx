////////////////////////////////////////////////////////////////////////
/// \file  ParticleListAction.cxx
/// \brief Use Geant4's user "hooks" to maintain a list of particles generated by Geant4.
///
/// \version $Id: ParticleListAction.cxx,v 1.8 2010/04/29 15:39:33 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "LArG4/ParticleListAction.h"
#include "G4Base/PrimaryParticleInformation.h"
#include "Simulation/sim.h"
#include "Simulation/ParticleList.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4Event.hh"
//#include <G4RunManager.hh>
#include "Geant4/G4Track.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PrimaryParticle.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4VUserPrimaryParticleInformation.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4VProcess.hh"
#include "Geant4/G4String.hh"

#include <TLorentzVector.h>
#include <TString.h>

#include <algorithm>

const G4bool debug = false;

namespace larg4 {

  // Initialize static members.
  simb::MCParticle* ParticleListAction::fparticle       = 0;
  int               ParticleListAction::fCurrentTrackID = sim::NoParticleId;
  int               ParticleListAction::fTrackIDOffset  = 0;

  //----------------------------------------------------------------------------
  // Constructor.
  ParticleListAction::ParticleListAction(double energyCut, 
                                         bool   storeTrajectories,
                                         bool   keepEMShowerDaughters)
    : fenergyCut(energyCut * GeV)
    , fparticleList(0)
    , fstoreTrajectories(storeTrajectories)
    , fKeepEMShowerDaughters(keepEMShowerDaughters)
  {
    // Create the particle list that we'll (re-)use during the course
    // of the Geant4 simulation.
    fparticleList = new sim::ParticleList;
    fParentIDMap.clear();
  }

  //----------------------------------------------------------------------------
  // Destructor.
  ParticleListAction::~ParticleListAction()
  {
    // Delete anything that we created with "new'.
    delete fparticleList;
  }

  //----------------------------------------------------------------------------
  // Begin the event
  void ParticleListAction::BeginOfEventAction(const G4Event*)
  {
    // Clear any previous particle information.
    fparticle = 0;
    fparticleList->clear();
    fParentIDMap.clear();
    fTrackIDToMCTruthIndex.clear();
    fCurrentTrackID = sim::NoParticleId;
  }

  //-------------------------------------------------------------
  // figure out the ultimate parentage of the particle with track ID
  // trackid
  // assume that the current track id has already been added to
  // the fParentIDMap
  int ParticleListAction::GetParentage(int trackid) const
  {
    int parentid = sim::NoParticleId;
    
    // search the fParentIDMap recursively until we have the parent id 
    // of the first EM particle that led to this one
    std::map<int,int>::const_iterator itr = fParentIDMap.find(trackid);
    while( itr != fParentIDMap.end() ){
      LOG_DEBUG("ParticleListAction") << "parentage for " << trackid
				      << " " << (*itr).second;
      // set the parentid to the current parent ID, when the loop ends
      // this id will be the first EM particle 
      parentid = (*itr).second;
      itr = fParentIDMap.find(parentid);
    }
    LOG_DEBUG("ParticleListAction") << "final parent ID " << parentid; 

    return parentid;
  }
  
  //----------------------------------------------------------------------------
  // Create our initial simb::MCParticle object and add it to the sim::ParticleList.
  void ParticleListAction::PreTrackingAction(const G4Track* track)
  {
    // Particle type.
    G4ParticleDefinition* particleDefinition = track->GetDefinition();
    G4int pdgCode = particleDefinition->GetPDGEncoding();
    
    // Get Geant4's ID number for this track.  This will be the same
    // ID number that we'll use in the ParticleList.
    G4int trackID = track->GetTrackID() + fTrackIDOffset;
    fCurrentTrackID = trackID;
    size_t mcTruthIndex = 0;
    // And the particle's parent:
    G4int parentID = track->GetParentID();
    
    std::string process_name = "unknown";
    
    // Is there an MCTruth object associated with this G4Track?  We
    // have to go up a "chain" of information to find out:
    const G4DynamicParticle* dynamicParticle = track->GetDynamicParticle();
    const G4PrimaryParticle* primaryParticle = dynamicParticle->GetPrimaryParticle();
    if ( primaryParticle != 0 ){
      const G4VUserPrimaryParticleInformation* gppi = primaryParticle->GetUserInformation();
      const g4b::PrimaryParticleInformation* ppi = dynamic_cast<const g4b::PrimaryParticleInformation*>(gppi);
      if ( ppi != 0 ){
        // If we've made it this far, a PrimaryParticleInformation
        // object exists and we are using a primary particle, set the
        // process name accordingly
        process_name = "primary";
        
        // primary particles should have parentID = 0, even if there
        // are multiple MCTruths for this event
        parentID = 0;
        
        mcTruthIndex = ppi->MCTruthIndex();
      } // end else no primary particle information
    } // Is there a G4PrimaryParticle?
    // If this is not a primary particle...
    else{
      // check if this particle was made in an EM shower, don't put it in the particle
      // list as we don't care about secondaries, tertiaries, etc for these showers
      // figure out what process is making this track - skip it if it is
      // one of pair production, compton scattering, photoelectric effect
      // bremstrahlung, annihilation, any ionization - who wants to save
      // a buttload of electrons that arent from a CC interaction?
      process_name = track->GetCreatorProcess()->GetProcessName();
      if( !fKeepEMShowerDaughters
         && (process_name.find("conv")               != std::string::npos
             || process_name.find("LowEnConversion") != std::string::npos
             || process_name.find("Pair")            != std::string::npos
             || process_name.find("compt")           != std::string::npos
             || process_name.find("Compt")           != std::string::npos
             || process_name.find("Brem")            != std::string::npos
             || process_name.find("phot")            != std::string::npos
             || process_name.find("Photo")           != std::string::npos
             || process_name.find("Ion")             != std::string::npos
             || process_name.find("annihil")         != std::string::npos)
         ){
        
        // figure out the ultimate parentage of this particle
        // first add this track id and its parent to the fParentIDMap
        fParentIDMap[trackID] = parentID;
        
        fCurrentTrackID = -1 * this->GetParentage(trackID);
        
        // check that fCurrentTrackID is in the particle list - it is possible
        // that this particle's parent is a particle that did not get tracked.
        // An example is a partent that was made due to muMinusCaptureAtRest
        // and the daughter was made by the phot process.  The parent likely
        // isn't saved in the particle list because it is below the energy cut
        // which will put a bogus track id value into the sim::IDE object for
        // the sim::SimChannel if we don't check it.
        if(fparticleList->find(fCurrentTrackID) == fparticleList->end() )
        fCurrentTrackID = sim::NoParticleId;
        
        // set fparticle to 0 as we are not stepping this particle and
        // adding trajectory points to it
        fparticle = 0;
        return;
        
      } // end if keeping EM shower daughters
      
      // Check the energy of the particle.  If it falls below the energy
      // cut, don't add it to our list.
      G4double energy = track->GetKineticEnergy();
      if( energy < fenergyCut ){
        fparticle = 0;
        
        // do add the particle to the parent id map though
        // and set the current track id to be it's ultimate parent
        fParentIDMap[trackID] = parentID;
        
        fCurrentTrackID = -1*this->GetParentage(trackID);
        
        return;
      }
      
      // check to see if the parent particle has been stored in the particle navigator
      // if not, then see if it is possible to walk up the fParentIDMap to find the
      // ultimate parent of this particle.  Use that ID as the parent ID for this
      // particle
      if( fparticleList->find(parentID) == fparticleList->end() ){
        // do add the particle to the parent id map
        // just in case it makes a daughter that we have to track as well
        fParentIDMap[trackID] = parentID;
        int pid = this->GetParentage(parentID);
        
        // if we still can't find the parent in the particle navigator,
        // we have to give up
        if( fparticleList->find(pid) == fparticleList->end() ){
          mf::LogWarning("ParticleListAction") << "can't find parent id: "
          << parentID
          << " in the particle list, or fParentIDMap."
          << " Make " << parentID << " the mother ID for"
          << " track ID " << fCurrentTrackID
          << " in the hope that it will aid debugging.";
        }
        else
          parentID = pid;
      }
      
    }// end if not a primary particle
    
    // Attempt to find the MCTruth index corresponding to the
    // current particle.  If the fCurrentTrackID is not in the
    // map try the parent ID, if that is not there, throw an
    // exception
    if(fTrackIDToMCTruthIndex.count(fCurrentTrackID) > 0 )
      mcTruthIndex = fTrackIDToMCTruthIndex.at(fCurrentTrackID);
    else if(fTrackIDToMCTruthIndex.count(parentID) > 0 )
      mcTruthIndex = fTrackIDToMCTruthIndex.at(parentID);
    else
      throw cet::exception("ParticleListAction") << "Cannot find MCTruth index for track id "
      << fCurrentTrackID << " or " << parentID;

    // This is probably the PDG mass, but just in case:
    double mass = dynamicParticle->GetMass()/GeV;
    
    // Create the sim::Particle object.
    fparticle = new simb::MCParticle( trackID, pdgCode, process_name, parentID, mass);
    
    // Polarization.
    const G4ThreeVector& polarization = track->GetPolarization();
    fparticle->SetPolarization( TVector3(polarization.x(),
                                         polarization.y(),
                                         polarization.z() ) );
    
    // Save the particle in the ParticleList.
    fparticleList->Add( fparticle );
    
    if(fTrackIDToMCTruthIndex.count(fCurrentTrackID) > 0)
    LOG_WARNING("ParticleListAction") << "attempting to put " << fCurrentTrackID
    << " into fTrackIDToMCTruthIndex map "
    << " particle is\n" << *fParticle;
    
    fTrackIDToMCTruthIndex[fCurrentTrackID] = mcTruthIndex;

  }
  
  //----------------------------------------------------------------------------
  void ParticleListAction::PostTrackingAction( const G4Track* aTrack)
  {
    if(fparticle && aTrack){
      fparticle->SetWeight(aTrack->GetWeight());
      G4String process = aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
      fparticle->SetEndProcess(process);
      
      
    }
    
    return;
  }
  
  //----------------------------------------------------------------------------
  // With every step, add to the particle's trajectory.
  void ParticleListAction::SteppingAction(const G4Step* step)
  {
    
    if ( fparticle == 0 ) {
      return;
    }
    
    // For the most part, we just want to add the post-step
    // information to the particle's trajectory.  There's one
    // exception: In PreTrackingAction, the correct time information
    // is not available.  So add the correct vertex information here.
    
    if ( fparticle->NumberTrajectoryPoints() == 0 ){
      
      // Get the pre/along-step information from the G4Step.
      const G4StepPoint* preStepPoint = step->GetPreStepPoint();
      
      const G4ThreeVector position = preStepPoint->GetPosition();
      G4double time = preStepPoint->GetGlobalTime();
      
      // Remember that LArSoft uses cm, ns, GeV.
      TLorentzVector fourPos(position.x() / cm,
                             position.y() / cm,
                             position.z() / cm,
                             time / ns);
      
      const G4ThreeVector momentum = preStepPoint->GetMomentum();
      const G4double energy = preStepPoint->GetTotalEnergy();
      TLorentzVector fourMom(momentum.x() / GeV,
                             momentum.y() / GeV,
                             momentum.z() / GeV,
                             energy / GeV);
      
      // Add the first point in the trajectory.
      fparticle->AddTrajectoryPoint( fourPos, fourMom );
    } // end if this is the first step
    
    // At this point, the particle is being transported through the
    // simulation. This method is being called for every voxel that
    // the track passes through, but we don't want to update the
    // trajectory information if we're just updating voxels. To check
    // for this, look at the process name for the step, and compare it
    // against the voxelization process name (set in PhysicsList.cxx).
    G4String process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4bool ignoreProcess = process.contains("LArVoxel") || process.contains("OpDetReadout");
    
    LOG_DEBUG("ParticleListAction::SteppingAction") << ": DEBUG - process='"
    << process << "'"
    << " ignoreProcess=" << ignoreProcess
    << " fstoreTrajectories="
    << fstoreTrajectories;
    
    // We store the initial creation point of the particle
    // and its final position (ie where it has no more energy, or at least < 1 eV) no matter
    // what, but whether we store the rest of the trajectory depends
    // on the process, and on a user switch.
    if ( fstoreTrajectories  &&  !ignoreProcess ){
      // Get the post-step information from the G4Step.
      const G4StepPoint* postStepPoint = step->GetPostStepPoint();
      
      const G4ThreeVector position = postStepPoint->GetPosition();
      G4double time = postStepPoint->GetGlobalTime();
      
      // Remember that LArSoft uses cm, ns, GeV.
      TLorentzVector fourPos(position.x() / cm,
                             position.y() / cm,
                             position.z() / cm,
                             time / ns );
      
      const G4ThreeVector momentum = postStepPoint->GetMomentum();
      const G4double energy = postStepPoint->GetTotalEnergy();
      TLorentzVector fourMom(momentum.x() / GeV,
                             momentum.y() / GeV,
                             momentum.z() / GeV,
                             energy / GeV );
      
      // Add another point in the trajectory.
      fparticle->AddTrajectoryPoint( fourPos, fourMom );
    }
  }
  
  //-------------------------------------------------------------
  std::map<int, size_t> ParticleListAction::TrackIDToMCTruthIndexMap() const
  {
    return fTrackIDToMCTruthIndex;
  }
  
  //----------------------------------------------------------------------------
  /// Utility class for the EndOfEventAction method: update the
  /// daughter relationships in the particle list.
  class UpdateDaughterInformation
  : public std::unary_function<sim::ParticleList::value_type, void>
  {
  public:
    UpdateDaughterInformation()
    : particleList(0)
    {}
    void SetParticleList( sim::ParticleList* p ) { particleList = p; }
    void operator()( sim::ParticleList::value_type& particleListEntry )
    {
      // We're looking at this Particle in the list.
      simb::MCParticle* particle = particleListEntry.second;
      
      // The parent ID of this particle.
      int parentID = particle->Mother();
      
      // If the parentID <= 0, this is a primary particle.
      if ( parentID <= 0 ) return;
      
      // If we get here, this particle is somebody's daughter.  Add
      // it to the list of daughter particles for that parent.
      
      // Get the parent particle from the list.
      sim::ParticleList::iterator parentEntry = particleList->find( parentID );
      
      if ( parentEntry == particleList->end() ){
        // We have an "orphan": a particle whose parent isn't
        // recorded in the particle list.  This is not signficant;
        // it's possible for a particle not to be saved in the list
        // because it failed an energy cut, but for it to have a
        // daughter that passed the cut (e.g., a nuclear decay).
        return;
      }
      
      // Add the current particle to the daughter list of the
      // parent.
      simb::MCParticle* parent = (*parentEntry).second;
      parent->AddDaughter( particle->TrackId() );
    }
  private:
    sim::ParticleList* particleList;
  };
  
  //----------------------------------------------------------------------------
  // There's one last thing to do: All the particles have their
  // parent IDs set (in PostTrackingAction), but we haven't set the
  // daughters yet.  That's done in this method.
  void ParticleListAction::EndOfEventAction(const G4Event*)
  {
    // Set up the utility class for the "for_each" algorithm.  (We only
    // need a separate set-up for the utility class because we need to
    // give it the pointer to the particle list.  We're using the STL
    // "for_each" instead of the C++ "for loop" because it's supposed
    // to be faster.
    UpdateDaughterInformation updateDaughterInformation;
    updateDaughterInformation.SetParticleList( fparticleList );
    
    // Update the daughter information for each particle in the list.
    std::for_each(fparticleList->begin(), 
                  fparticleList->end(), 
                  updateDaughterInformation);
  }
  
  //----------------------------------------------------------------------------
  // Returns the ParticleList accumulated during the current event.
  const sim::ParticleList* ParticleListAction::GetList() const
  {
    // check if the ParticleNavigator has entries, and if
    // so grab the highest track id value from it to
    // add to the fTrackIDOffset
    int highestID = 0;
    for( auto pn = fparticleList->begin(); pn != fparticleList->end(); pn++)
    if( (*pn).first > highestID ) highestID = (*pn).first;
    
    fTrackIDOffset = highestID + 1;
    
    return fparticleList;
  }
  
} // namespace LArG4
