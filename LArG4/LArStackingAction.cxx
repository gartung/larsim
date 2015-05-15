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

#include "LArG4/LArStackingAction.h"
#include "Geometry/Geometry.h"

#include "Geant4/G4SDManager.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Event.hh"
#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4TrackStatus.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ios.hh"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

LArStackingAction::LArStackingAction(G4int stack, G4int alg /*=0*/)
 : fStage(0)
 , fStack(stack)
 , fAlg(alg)
 , freqMuon(2)
 , freqIsoMuon(0)
 , freqIso(10)
 , fangRoI(30.*deg)
{ 
  //theMessenger = new LArStackingActionMessenger(this);
}

LArStackingAction::~LArStackingAction()
{ //delete theMessenger; 
}

G4ClassificationOfNewTrack 
LArStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  
  G4ClassificationOfNewTrack classification = fUrgent;  
  
  art::ServiceHandle<geo::Geometry> geom; 
  
  if (fAlg > 0 && aTrack->GetDefinition()->GetPDGEncoding()==11) {
    
    double bounds[6] = {0};
    bool insideCryo = false;
    const G4ThreeVector tr4Pos = aTrack->GetPosition();
    const TVector3 trPos(tr4Pos.x()/cm,tr4Pos.y()/cm,tr4Pos.z()/cm);
    for(unsigned int c = 0; c < geom->Ncryostats(); ++c) {
      geom->CryostatBoundaries(bounds, c);
      if ( trPos.X() - bounds[0] > 0.0 && trPos.X() - bounds[1] < 0.0 &&
           trPos.Y() - bounds[2] > 0.0 && trPos.Y() - bounds[3] < 0.0 &&
	   trPos.Z() - bounds[4] > 0.0 && trPos.Z() - bounds[5] < 0.0 ) {
        insideCryo = true;
	break;
      }
    }
    
    if (!insideCryo && aTrack->GetKineticEnergy() < 5.0) return fKill;
  
    //reject electrons that point away from cryostat
    if (fAlg > 1 && !insideCryo) {

      //unit direction vector
      TVector3 p_hat(aTrack->GetMomentumDirection.x(), aTrack->GetMomentumDirection.y(), aTrack->GetMomentumDirection.z());
      p_hat = p_hat.Unit();

      //Start position vector
      TVector3 p0(aTrack->GetPosition().x()/cm, aTrack->GetPosition().y()/cm, aTrack->GetPosition().z()/cm);

      double bounds[6] = {0};
      double min_deflection = 4.0;
      for (unsigned int c = 0; c < geom->Ncryostats(); ++c) {
	geom->CryostatBoundaries(bounds, c);

	//does particle point at the cryostat?
	if (PointsAtCryoStat(p0, p_hat, bounds)) {
          min_deflection = 0.0;
	  break;
	}

	//make edges
	std::vector< std::pair<TVector3, TVector3> > edges;
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[2], bounds[4]), TVector3(bounds[1]-bounds[0], 0.0, 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[3], bounds[4]), TVector3(bounds[1]-bounds[0], 0.0, 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[2], bounds[5]), TVector3(bounds[1]-bounds[0], 0.0, 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[3], bounds[5]), TVector3(bounds[1]-bounds[0], 0.0, 0.0) ) );

	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[2], bounds[4]), TVector3(0.0, bounds[3]-bounds[2], 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[1], bounds[2], bounds[4]), TVector3(0.0, bounds[3]-bounds[2], 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[2], bounds[5]), TVector3(0.0, bounds[3]-bounds[2], 0.0) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[1], bounds[2], bounds[5]), TVector3(0.0, bounds[3]-bounds[2], 0.0) ) );

	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[2], bounds[4]), TVector3(0.0, 0.0, bounds[5]-bounds[4]) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[1], bounds[2], bounds[4]), TVector3(0.0, 0.0, bounds[5]-bounds[4]) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[0], bounds[3], bounds[4]), TVector3(0.0, 0.0, bounds[5]-bounds[4]) ) );
	edges.push_back( std::pair<TVector3,TVector3>( TVector3(bounds[1], bounds[3], bounds[4]), TVector3(0.0, 0.0, bounds[5]-bounds[4]) ) );

	//loop over edges and find minimum deflection
	for (auto itE = edges.begin(); itE != edges.end(); ++itE) {

	  double A = (itE->first - p0).Dot(p_hat);
	  double B = (itE->second).Dot(p_hat);

	  double D, E;
	  if ( (itE->second).x() != 0.0 ) {
	    D = itE->second.x();
	    E = -p0.x();
	  }
	  else if ( (itE->second).y() != 0.0 ) {
	    D = itE->second.y();
	    E = -p0.y();
	  }
	  else if ( (itE->second).z() != 0.0 ) {
	    D = itE->second.z();
	    E = -p0.z();
	  }
	  else { std::cout << "PROBLEM" <<std::endl; }
	  double C2 = (itE->first - p0).Mag2() - E*E;

	  double tmp_min = 99.9;
	  if (B==0.0) {
	    if (A==0.0) {
	      std::cout<<"DEGENERATE CASE: A=B=0"<<std::endl;
	      tmp_min = 3.14145 / 2.0; //90 degrees
	    }
	    else if ( 0.0 <= -E/D && 1.0 >= -E/D) {
	      std::cout<<"DEGENERATE CASE: -E/D = t"<<std::endl;
	      tmp_min = acos(A/sqrt(C2));
	    }
	    else {
	      std::cout<<"ODD THIRD CASE"<<std::endl;
            }

	    if (tmp_min < min_deflection) {
	      min_deflection = tmp_min;
	      if (min_deflection > 3.5) std::cout<<"BAD DEFLECTION"<<std::endl;
	      continue;
	    }
	    else continue;
	  }

	  double tmp_minlocation = B*C2/(D*(A*D-B*E)) - E/D;
	  if (tmp_minlocation < 0.0 || tmp_minlocation > 1.0) {
	    double arg_at_0 = A/sqrt(C2 + E*E);
	    double arg_at_1 = (A+B)/sqrt(C2+ (D+E)*(D+E));

	    double val_at_0 = acos(arg_at_0);
	    double val_at_1 = acos(arg_at_1);
	    if (val_at_0 < min_deflection) min_deflection = val_at_0;
	    if (val_at_1 < min_deflection) min_deflection = val_at_1;

	    if (min_deflection > 3.5) std::cout<<"BAD DEFLECTION OB "<<std::endl;
	  }
	  else {
	    tmp_min = acos( (A+B*tmp_minlocation)/sqrt(C2 + (D*tmp_minlocation + E)*(D*tmp_minlocation+E)  ) );
	    if (tmp_min < min_deflection) min_deflection = tmp_min;
	    if (min_deflection > 3.5) std::cout<<"BAD DEFLECTION NORM"<<std::endl;
	  }
	}//end loop over edges
      }//end loop over cryostats
      
      if (min_deflection > 1.4) return fKill; //hardcoded cut for now
    }//end if fAlg > 1
    
    return classification;
  }//end if fAlg > 0
  
  //The following only runs for fAlg==0  
  
  classification = fWaiting;
  TString volName(InsideTPC(aTrack));
  Double_t buffer = 500; // Keep muNucl neutrals within 5m (for now) of larTPC.

  // These 3 for now for investigation in gdb.
  //int pdg = aTrack->GetDefinition()->GetPDGEncoding();
  int ppdg = aTrack->GetParentID();
  TString process("NA");
  if (ppdg) process = (TString)aTrack->GetCreatorProcess()->GetProcessName();


  switch(fStage){
  case 0: // Fstage 0 : Primary muons only
    if (aTrack->GetParentID()==0)
    {
      G4ParticleDefinition *particleType = aTrack->GetDefinition();
      if( ((particleType==G4MuonPlus::MuonPlusDefinition())
	  || (particleType==G4MuonMinus::MuonMinusDefinition())
	   )
	  && !volName.Contains("unknown")
	  ){ 
	classification = fUrgent; 
      }
    }
    if (volName.Contains("unknown"))   classification = fKill;
    break;
    
  case 1: // Stage 1 : K0,Lambda,n's made urgent here.
          //           Suspended tracks will be sent to the waiting stack
    if(aTrack->GetTrackStatus()==fSuspend) { break; }

    if ((aTrack->GetDefinition()->GetPDGEncoding()==2112 || aTrack->GetDefinition()->GetPDGEncoding()==130 || aTrack->GetDefinition()->GetPDGEncoding()==310 || aTrack->GetDefinition()->GetPDGEncoding()==311 || aTrack->GetDefinition()->GetPDGEncoding()==3122 ) && (aTrack->GetParentID()==1) && !volName.Contains("unknown"))
      {
	
	const G4ThreeVector tr4Pos = aTrack->GetPosition();
	// G4 returns positions in mm, have to convert to cm for LArSoft coordinate systems
	const TVector3 trPos(tr4Pos.x()/cm,tr4Pos.y()/cm,tr4Pos.z()/cm);
	//double locNeut = trPos.Mag();
	classification = fUrgent; 
	// std::cout << "LArStackingAction: DetHalfWidth, Height, FullLength: " << geom->DetHalfWidth() << ", " << geom->DetHalfHeight() << ", " << geom->DetLength() << std::endl;
	
	if (
	    trPos.X() < (geom->DetHalfWidth()*2.0 + buffer) && trPos.X() > (-buffer) &&
	    trPos.Y() < (geom->DetHalfHeight()*2.0 + buffer) && trPos.Y() > (-geom->DetHalfHeight()*2.0 - buffer) &&
	    trPos.Z() < (geom->DetLength() + buffer) && trPos.Z() > (-buffer) 
	    )

	  { classification = fUrgent; break; }
	// These tracks need to be "scored" cuz every now and then they
	// might get to the LAr.
	else
	  { classification = fKill; break; }
	
	
      }
    
	//    if(aTrack->GetDefinition()->GetPDGCharge()==0.) { break; }
    if (volName.Contains("unknown"))   classification = fKill;
    break;

  default: 
    // Track all other Primaries. Accept all secondaries in TPC.
    // Kill muon ionization electrons outside TPC
    // ignore primaries since they have no creator process
    
    if(aTrack->GetParentID() == 0 && !volName.Contains("unknown")){ 
      classification = fUrgent;
      break;
    }
       
    if(volName.Contains(geom->GetLArTPCVolumeName()) && aTrack->GetParentID()!=0)
      { 
	classification = fUrgent;
	if (fStack & 0x4 && 
	    aTrack->GetCreatorProcess()->GetProcessName().contains("muIoni")
	    ) 
	  {
	    classification = fKill;
	  }
	break;
      }
    else if (volName.Contains("unknown") ){
      classification = fKill;
      break;
    }
    // Leave this here, even though I claim we've Killed these in stage 2.
    if(aTrack->GetDefinition()->GetPDGEncoding()==11 
       && aTrack->GetCreatorProcess()->GetProcessName().contains("muIoni") )
      {
	classification = fKill;
	break;
      }
    // For now, kill every other thing, no matter where it is.
    classification = fKill;

  } // end switch

  return classification;
}


std::string LArStackingAction::InsideTPC(const G4Track * aTrack)
{

  art::ServiceHandle<geo::Geometry> geom;
  const G4ThreeVector tr4Pos = aTrack->GetPosition();

  // G4 returns positions in mm, have to convert to cm for LArSoft coordinate systems
  const TVector3 trPos(tr4Pos.x()/cm,tr4Pos.y()/cm,tr4Pos.z()/cm);

  const std::string volName(geom->VolumeName(trPos));

  return volName;
}


void LArStackingAction::NewStage()
{

  // Here when Urgent stack is empty. Waiting stack about to be made Urgent,
  // upon saying ReClassify().
  fStage++;

    // I yanked the ExN04's use here of stackManager->clear(), which clears stack
    // and prepares to end the event. Think I may wanna do something like this if
    // muon has been tracked, its doca is large, there are no hit voxels in the TPC,
    // and I'm in further need of optimization.

  if(fStage==1){
  // Stage 0->1 : check if at least "reqMuon" hits on muon chamber
  //              otherwise abort current event

    stackManager->ReClassify();
    return;
  }

  else if(fStage==2){
  // Stage 1->2 : check the isolation of muon tracks
  //              at least "reqIsoMuon" isolated muons
  //              otherwise abort current event.
  //              Isolation requires "reqIso" or less hits
  //              (including own hits) in the RoI region
  //              in the tracker layers.
    stackManager->ReClassify();
    return;
  }

  else{
  // Other stage change : just re-classify
    stackManager->ReClassify();
  }
}
    
void LArStackingAction::PrepareNewEvent()
{ 
  fStage = 0; 
  //trkHits = 0;
  //muonHits = 0;
}


