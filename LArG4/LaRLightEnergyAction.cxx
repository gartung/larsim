////////////////////////////////////////////////////////////////////////
/// \file  LaRLightEnergyAction.cxx
/// \brief Use Geant4's user "hooks" to maintain a list of particles generated by Geant4.
///
/// Based on G4BadIdeaAction written by Brian Rebel (brebel@fnal.gov)
// Name is a reminescence of LA(small cap. R)TF case. Capitalization will change accordingly
////////////////////////////////////////////////////////////////////////

#include "LArG4/LaRLightEnergyAction.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4Track.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4TransportationManager.hh" 
#include "Geant4/G4Navigator.hh" 
#include "Geant4/G4VProcess.hh" 
#include "Geant4/G4ProcessManager.hh" 
#include "Geant4/G4OpBoundaryProcess.hh" 

#include <algorithm>
#include <vector>

//added to write energy deposit in histogram

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

//#include "art/Framework/Principal/Event.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include <TTree.h>
G4int ktrack1=0;
G4int ttracks1=0;
G4int tracknumber=0;//>current track number 
G4int tracknumber_old=0;//> used to store track number to check whether it's changed
G4bool change_track=false;//>is track changed?
std::vector<std::string> volnames ;//>vector of volume names - used for debugging
std::vector<std::vector<double>> volpos ;//>vector of volume positions - used for debugging
std::vector<double> tmppos;
G4bool change_event=false;//>is event nr. changed?
double energydep1=0.;//>deposited energy
G4int evnr=0;//>temp event number
G4int refl=0;//>number of reflections - used for debugging - may be accesible in output data(planned)
G4int refltpb=0;//> reflections from TPB - as above
G4bool test2=false;
G4bool test1=false;
G4bool tracktest=false;
G4double refl_all=0;//> all reflections count
G4String processn;//> number of physical process
G4String thePrePVname1;
G4String thePostPVname1;

// Photon variables defined at each step, for use 
// in temporary velocity bug fix. -wforeman
double globalTime, velocity_G4, velocity_step;

namespace larg4 {

  //----------------------------------------------------------------------------
  // Constructor.
  LaRLightEnergyAction::LaRLightEnergyAction(int trkOption)
  {
  }

  //----------------------------------------------------------------------------
  // Destructor.
  LaRLightEnergyAction::~LaRLightEnergyAction()
  {
  }

  //----------------------------------------------------------------------------
  // With every step, add to the particle's trajectory.
  void LaRLightEnergyAction::SteppingAction(const G4Step* step)
  {


	change_track=false;//>new track
	change_event=false;//>new event
	tracknumber =step->GetTrack()->GetTrackID();
	if(tracknumber!=0){
		if(tracknumber!=tracknumber_old){
			ttracks1++;
			tracknumber_old=tracknumber;
			change_track=true;
		}
	}
	else{
		ttracks1++;
		tracknumber_old=tracknumber;
		change_track=true;

	}
	if (change_track==true){
		//volnames.clear();//used for debugging
		//volpos.clear();//used for debugging
		energydep1=0.;
		test1=false;
		test2=false;
		
	}

	if(eventnumber_fast!=0){
		if(eventnumber_fast!=evnr){
			change_event=true;
			evnr=eventnumber_fast;
	//		mf::LogWarning("LaRLightEnergyAction")<<"--------------EVENT NUMBER"<<evnr<<std::endl;
		}

	}
	else{
		change_event=true;
		evnr=eventnumber_fast;
	//	mf::LogWarning("LaRLightEnergyAction")<<" EVENT NUMBER laRlight energyaction "<<evnr<<std::endl;	
	}
	if(change_event==true){
		energy_deposit_step=0.;


	}

        // Temporary fix for problem where  DeltaTime on the first step
        // of optical photon propagation is calculated incorrectly. -wforeman
        globalTime = step->GetTrack()->GetGlobalTime();
        velocity_G4 = step->GetTrack()->GetVelocity();
        velocity_step = step->GetStepLength() / step->GetDeltaTime();
        if ( (step->GetTrack()->GetDefinition()->GetPDGEncoding()==0) && 
             fabs(velocity_G4 - velocity_step) > 0.0001 ) {
          // Subtract the faulty step time from the global time,
          // and add the correct step time based on G4 velocity.
          step->GetPostStepPoint()->SetGlobalTime( 
           globalTime - step->GetDeltaTime() + step->GetStepLength()/velocity_G4);
        }

	if (step->GetTrack()->GetDefinition()->GetPDGEncoding()==0){
		tmppos.clear();

		//check if reflected from wall && thePrePVname1.contains("refl") double(step->GetTrack()->GetKineticEnergy()/eV) <4.2) &&
		thePrePVname1 = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
		thePostPVname1 = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
		if (step->GetTrack()) tracktest=true;
		else tracktest=false;


		if( tracktest && (double(step->GetTrack()->GetKineticEnergy()/eV) >7.2) && (step->GetTrack()->GetParentID()==0) &&!test2 && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="volBeamBox_PV") && step->GetPostStepPoint()->GetPhysicalVolume()->GetName().contains("voltpb")){
			refltpb++;
			test2=true;
		}

		if(  !test1 && step->GetPreStepPoint()->GetPhysicalVolume()->GetName().contains("voltpb") && step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="volBeamBox_PV"){
			refl++;
			//std::cout<<"reflected photons!!!!!!!! "<<refl<<" dirceted into tpb "<<refltpb<<" reflectivity  "<<double(refl)/double(refltpb)<<" energy "<<step->GetTrack()->GetKineticEnergy()/eV<<std::endl;	
			test1=true;

		}


 
		if ( step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="volArgon_cap_L_PV" && ((step->GetPostStepPoint()->GetPosition()[2]/cm)>90.0 ||(step->GetPostStepPoint()->GetPosition()[2]/cm)<0.0)){
			step->GetTrack()->SetTrackStatus(fStopAndKill);
			ktrack1++;
			}//if in argon cap behind TPC - killing track
	}//if opticalphoton - check


if(step->GetTrack()->GetVolume()->GetName()=="volBeamBox_PV" || step->GetTrack()->GetVolume()->GetName()=="volTPCActive_PV"){ 
		

			 energy_deposit_step+=step->GetTotalEnergyDeposit();
			//mf::LogWarning("LaRLightEnergyAction")<<"energy deposit step filled!!!!!!!!! "<<energy_deposit_step<<" change event "<<change_event<<std::endl;

	}
			
    
    return;
  }

} // namespace LArG4

