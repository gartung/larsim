////////////////////////////////////////////////////////////////////////
/// \file  G4BadIdeaAction.cxx
/// \brief Use Geant4's user "hooks" to maintain a list of particles generated by Geant4.
///
/// \version $Id: G4BadIdeaAction.cxx,v 1.8 2010/04/29 15:39:33 seligman Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "LArG4/G4BadIdeaAction.h"

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
G4int tracknumber=0;
G4int tracknumber_old=0;
G4bool change_track=false;
std::vector<std::string> volnames ;
std::vector<std::vector<double>> volpos ;
std::vector<double> tmppos;
G4bool change_event=false;
double energydep1=0.;
G4int evnr=0;
G4int refl=0;
G4int refltpb=0;
G4bool test2=false;
G4bool test3=false;
G4bool test1=false;
G4bool tracktest=false;
G4double refl_all=0;
G4String processn;
namespace larg4 {

  //----------------------------------------------------------------------------
  // Constructor.
  G4BadIdeaAction::G4BadIdeaAction(int trkOption) :
    fNoIncomingMuons(trkOption)
  {
    // trkOption comes from LArG4's fSmartStacking
    // Negative values effect action in this routine. Positive values
    // effect action in LArStackingAction.
    mf::LogWarning("G4BadIdeaAction") << "instantiating the G4BadIdeaAction \n"
				      << "This UserAction is only to be used with "
				      << "Geant4 v4.9.4.p02 to solve a stepping bug.\n"
				      << "If you are using a different version of G4, "
				      << "remove this UserAction from your list in LArG4.cxx";
  }

  //----------------------------------------------------------------------------
  // Destructor.
  G4BadIdeaAction::~G4BadIdeaAction()
  {
  }

  //----------------------------------------------------------------------------
  // With every step, add to the particle's trajectory.
  void G4BadIdeaAction::SteppingAction(const G4Step* step)
  {


change_track=false;
change_event=false;
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
 volnames.clear();
volpos.clear();
energydep1=0.;
test1=false;
test2=false;
test3=false;
}

if(eventnumber_fast!=0){
	if(eventnumber_fast!=evnr){
		change_event=true;
		evnr=eventnumber_fast;
std::cout<<"--------------EVENT NUMBER"<<evnr<<std::endl;
	}

}
else{
	change_event=true;
	evnr=eventnumber_fast;
std::cout<<" EVENT NUMBER badideaaction "<<evnr<<std::endl;	
}
if(change_event==true){
energy_deposit_step=0.;
//phot_refl_ev=0.;

}
  //G4Navigator* theNavigator =
    //G4TransportationManager::GetTransportationManager()->
    //GetNavigatorForTracking();

    ////////////////////////////////////////////////////////////////////
    ///\todo Do not copy the code below. Contact Brian Rebel, Eric Church, 
    ///\todo Bill Seligman and Andrzej Szelc for reasons why not to!
    ////////////////////////////////////////////////////////////////////////

    // If the step size is such that the particle appears to be "stuck"
    // in its trajectory, give it a kick.
    const double epsilon   = 1000. * std::numeric_limits<double>::epsilon();
    const double stepSize  = step->GetStepLength();
    G4Track* nonConstTrack = const_cast<G4Track*>( step->GetTrack() );
//G4bool valid;
  //theNavigator->GetLocalExitNormal(&valid);
 const G4StepPoint* preStepPoint = step->GetPreStepPoint();
      const G4StepPoint* postStepPoint = step->GetPostStepPoint();
      const G4ThreeVector position = preStepPoint->GetPosition();

//PRINT OUT PRIMARY PARTICLE INFO HERE it's proton PDG 2212

if (step->GetTrack()->GetDefinition()->GetPDGEncoding()==0){
tmppos.clear();

//check if reflected from wall && thePrePVname1.contains("refl") double(step->GetTrack()->GetKineticEnergy()/eV) <4.2) &&
	G4String thePrePVname1 = preStepPoint->GetPhysicalVolume()->GetName();
	G4String thePostPVname1 = postStepPoint->GetPhysicalVolume()->GetName();
	if (step->GetTrack()) tracktest=true;
	else tracktest=false;
//std::cout<<"coming out of argon "<<tracktest<<" energy "<<double(step->GetTrack()->GetKineticEnergy()/eV)<<" parent ID, 0=primary "<<step->GetTrack()->GetParentID()<<" track change "<<test2<<" pre volume "<<preStepPoint->GetPhysicalVolume()->GetName()<<" post volume "<<postStepPoint->GetPhysicalVolume()->GetName()<<std::endl;

if( tracktest && (double(step->GetTrack()->GetKineticEnergy()/eV) >7.2) && (step->GetTrack()->GetParentID()==0) &&!test2 && (preStepPoint->GetPhysicalVolume()->GetName()=="volBeamBox_PV") && postStepPoint->GetPhysicalVolume()->GetName().contains("voltpb")){
	refltpb++;


	test2=true;
		}

if(  !test1 && preStepPoint->GetPhysicalVolume()->GetName().contains("voltpb") && postStepPoint->GetPhysicalVolume()->GetName()=="volBeamBox_PV"){
	refl++;
	//std::cout<<"reflected photons!!!!!!!! "<<refl<<" dirceted into tpb "<<refltpb<<" reflectivity  "<<double(refl)/double(refltpb)<<" energy "<<step->GetTrack()->GetKineticEnergy()/eV<<std::endl;	
test1=true;

		}






//tmppos.push_back(double(postStepPoint->GetPosition()[0]/cm));
//tmppos.push_back(double(postStepPoint->GetPosition()[1]/cm));
//tmppos.push_back(double(postStepPoint->GetPosition()[2]/cm));

//volpos.push_back(tmppos);      

//volnames.push_back(postStepPoint->GetPhysicalVolume()->GetName());
//std::cout<<"111111111111111111111111 process "<<postStepPoint->GetProcessDefinedStep()->GetProcessName()<<" ID "<<step->GetTrack()->GetTrackID()<<" "<<postStepPoint->GetPhysicalVolume()->GetName()<<" "<<postStepPoint->GetStepStatus()<<std::endl;
//   G4OpBoundaryProcessStatus boundaryStatus=Undefined;
//G4OpBoundaryProcess* boundary;
  
    //find the boundary process only once
 /*
      G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
     G4ProcessVector* pv = pm->GetProcessList();

  //  std::cout<<nprocesses<<"is the process list length "<<std::endl;
   for(int ij=0;ij<nprocesses;ij++){

std::cout<<(*pv)[ij]->GetProcessName()<<std::endl;
if((*pv)[ij]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[ij];
//std::cout<<boundary->GetStatus()<<std::endl;
		} 
	}

*/
if ( postStepPoint->GetPhysicalVolume()->GetName()=="volArgon_cap_L_PV" && ((postStepPoint->GetPosition()[2]/cm)>90.0 ||(postStepPoint->GetPosition()[2]/cm)<0.0)){
 step->GetTrack()->SetTrackStatus(fStopAndKill);
ktrack1++;
//std::cout<<" >>>>>>>>>>>>>>>>>>>>>>>>>bad idea action track killed in volArgon_cap "<<ktrack1<<"ID "<<tracknumber<<" total track number is "<<ttracks1<<" volume "<<preStepPoint->GetPhysicalVolume()->GetName()<<" particle "<<step->GetTrack()->GetDefinition()->GetName()<<std::endl;
//for(int ii=0;ii<int(volnames.size());ii++) std::cout<<ii<<" previous volume "<<volnames[ii]<<std::endl;
//for(int ij=0;ij<int(volpos.size());ij++) std::cout<<ij<<" previous pos x "<<volpos[ij][0]<<" previous pos y "<<volpos[ij][1]<<" previous pos z "<<volpos[ij][2]<<std::endl;
		}//if in argon cap behind TPC - killing track
	}//if opticalphoton - check



    if (step->GetTrack()->GetCurrentStepNumber() > 5e4 
	&&  stepSize < epsilon 
	&& (step->GetTrack()->GetCurrentStepNumber()%1000) == 1  ){

      // Cast away the const-ness of the pointer to G4Step.
      // This is dangerous. Don't do this at home. We're
      // only doing this because we're desperate.
      // The need to do this is the result of a bug in Geant4 v4.9.4.p02
      // We should no longer call this code when we move beyond that version - but you still do it in 4.9.6.. 
      const double kick = 0.001;
      
      G4ThreeVector aValue = nonConstTrack->GetPosition();
      
      mf::LogWarning("G4BadIdeaAction") << "##### In endless loop. Kicking particle by "
					<< " (+0.001,+0.001,+0.001) --- "
					<< " PDG and encoding "
					<< step->GetTrack()->GetDynamicParticle()->GetPDGcode()  
					<< " "
					<<  step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() 
					<<  " current step number: " 
					<< step->GetTrack()->GetCurrentStepNumber() 
					<<  " stepsize: "  << stepSize 
					<<  " x,y,z  " 
					<< aValue.x() << " " << aValue.y() << " " << aValue.z();
      
      G4ThreeVector  translate(kick,kick,kick);
      aValue+=translate;
      
      nonConstTrack->SetPosition(aValue);
    }


    if (fNoIncomingMuons<0)
      {
	// This is for overlays of, say, rock muons, which we have
	G4StepPoint * thePrePoint = step->GetPreStepPoint(); 
	G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume(); 
	G4String thePrePVname = thePrePV->GetName();
	G4StepPoint * thePostPoint = step->GetPostStepPoint(); 
	G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume(); 
	G4String thePostPVname("null");
	if (thePostPV) thePostPVname = thePostPV->GetName();
	
	if (abs(step->GetTrack()->GetDynamicParticle()->GetPDGcode()) == 13
	    && thePostPVname.contains("volTPCActive")
	    && !thePrePVname.contains("volTPCActive")
	    )
	  ((G4Track*)nonConstTrack)->SetTrackStatus(fStopAndKill);
      }





if(step->GetTrack()->GetVolume()->GetName()=="volBeamBox_PV" || step->GetTrack()->GetVolume()->GetName()=="volTPCActive_PV"){ 
			//std::cout<<"in box!!!!!!!!! "<<std::endl;

			 energy_deposit_step+=step->GetTotalEnergyDeposit();
		
// std::cout<<"energy deposit step filled!!!!!!!!! "<<energy_deposit_step<<" change event "<<change_event<<std::endl;

	}
			




 
 



/*
std::cout<<" >>>>>>>>>>>>>>>>>>>>>>>>>bad idea action primary particle PDG -- track nr --Track ID "<<step->GetTrack()->GetDefinition()->GetPDGEncoding()<<" "<<ktrack1<<" "<<step->GetTrack()->GetTrackID()<<" volume "<<preStepPoint->GetPhysicalVolume()->GetName()<<" pre step point position x "<<preStepPoint->GetPosition()[0]/cm<<" position y "<<preStepPoint->GetPosition()[1]/cm<<" position z "<<preStepPoint->GetPosition()[2]/cm<<" ENERGY DEPOSIT IN THIS STEP "<<step->GetTotalEnergyDeposit()/MeV<<" TOTAL "<<energy_deposit_step<< std::endl;
*/
 //else std::cout<<" in volWorld or no step "<<std::endl;
    
    return;
  }

} // namespace LArG4

