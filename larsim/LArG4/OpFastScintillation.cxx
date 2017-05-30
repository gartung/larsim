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
// Added disclaimer : This code is gross.  Thats basically because
//  it adheres to the original, gross Geant4 implementation.
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
// $Id: OpFastScintillation.cc,v 1.38 2010-12-15 07:39:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        OpFastScintillation.cc
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma)
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

//#include "Geant4/g4ios.hh"
#include "TLorentzVector.h"

#include "Geant4/globals.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4EmProcessSubType.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4Poisson.hh"

#include "larsim/LArG4/IonizationAndScintillation.h"
#include "larsim/LArG4/OpFastScintillation.hh"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/OpDetGeo.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// support libraries
#include "cetlib/exception.h"


namespace larg4
{
  //-------------------------------------------------------------
  OpFastScintillation::OpFastScintillation(const G4String& processName,
                                           G4ProcessType type)
    : G4VRestDiscreteProcess(processName, type),
      fSlowIntegralTable(0),
      fFastIntegralTable(0),
      fTrackSecondariesFirst(false),
      fFiniteRiseTime(false),
      fExcitationRatio(1),
      fEMSaturation(0)
  {
    SetProcessSubType(25);

    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
    fScintillationByParticleType = larp->ScintByParticleType();

    if (verboseLevel>0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }

    BuildPhysicsTable();
  }

  //-------------------------------------------------------------
  OpFastScintillation::OpFastScintillation(const OpFastScintillation& rhs)
    : G4VRestDiscreteProcess(rhs.GetProcessName(), rhs.GetProcessType())
  {
    fSlowIntegralTable = 0;
    fFastIntegralTable = 0;

    fTrackSecondariesFirst       = rhs.GetTrackSecondariesFirst();
    fFiniteRiseTime              = rhs.GetFiniteRiseTime();
    fExcitationRatio             = rhs.GetScintillationExcitationRatio();
    fScintillationByParticleType = rhs.GetScintillationByParticleType();
    fEMSaturation                = rhs.GetSaturation();

    BuildPhysicsTable();
  }

  //-------------------------------------------------------------
  OpFastScintillation::~OpFastScintillation()
  {
    if(fFastIntegralTable){
      fFastIntegralTable->clearAndDestroy();
      delete fFastIntegralTable;
    }
    if(fSlowIntegralTable){
      fSlowIntegralTable->clearAndDestroy();
      delete fSlowIntegralTable;
    }
  }

  //-------------------------------------------------------------
  // This routine simply calls the equivalent PostStepDoIt since all the
  // necessary information resides in aStep.GetTotalEnergyDeposit()
  G4VParticleChange*
  OpFastScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)
  {
    return OpFastScintillation::PostStepDoIt(aTrack, aStep);
  }

  // This routine is called for each tracking step of a charged particle
  // in a scintillator. A Poisson/Gauss-distributed number of photons is
  // generated according to the scintillation yield formula, distributed
  // evenly along the track segment and uniformly into 4pi.
  G4VParticleChange*
  OpFastScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  {
    aParticleChange.Initialize(aTrack);

    // Check that we are in a material with a properties table, if not
    // just return
    const G4Material* aMaterial = aTrack.GetMaterial();
    G4MaterialPropertiesTable* aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();
    if (!aMaterialPropertiesTable)
      return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();

    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4ThreeVector p0 = aStep.GetDeltaPosition().unit();

    // We don't want to produce any trackable G4 secondaries
    aParticleChange.SetNumberOfSecondaries(0);


    // Retrieve the Scintillation Integral for this material
    // new G4PhysicsOrderedFreeVector allocated to hold CII's


    // Some explanation for later improvements to scint yield code:
    //
    // What does G4 do here?
    //  It produces light in 2 steps, fast (scnt=1) then slow (scnt=2)
    //
    // The ratio of slow photons to fast photons is related	by the yieldratio
    //  parameter.  G4's poisson fluctuating scheme is a bit different to ours
    //  - we should check that they are equivalent.
    //
    // G4 poisson fluctuates the number of initial photons then divides them
    //  with a constant factor between fast + slow, whereas we poisson
    //  fluctuate separateyly the fast and slow detection numbers.
    //

    // get the number of photons produced from the IonizationAndScintillation
    // singleton
    larg4::IonizationAndScintillation::Instance()->Reset(&aStep);
    double MeanNumberOfPhotons = larg4::IonizationAndScintillation::Instance()->NumberScintillationPhotons();
    RecordPhotonsProduced(aStep, MeanNumberOfPhotons);


    if (verboseLevel>0) {
      G4cout << "\n Exiting from OpFastScintillation::DoIt -- NumberOfSecondaries = "
             << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }

    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  //-------------------------------------------------------------
  bool OpFastScintillation::RecordPhotonsProduced(const G4Step& aStep, double MeanNumberOfPhotons)
  {
    // Get the pointer to the fast scintillation table
    OpDetPhotonTable * fst = OpDetPhotonTable::Instance();
    OpDetPhotonTable* litefst = OpDetPhotonTable::Instance();

    // Get the pointer to the visibility service
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    const G4Track * aTrack = aStep.GetTrack();

    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

    const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
    const G4Material* aMaterial = aTrack->GetMaterial();

    G4int materialIndex = aMaterial->GetIndex();


    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
    //G4double      t0 = pPreStepPoint->GetGlobalTime() - fGlobalTimeOffset;
    G4double      t0 = pPreStepPoint->GetGlobalTime();


    G4MaterialPropertiesTable* aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();

    // Get the visibility vector for this point
    size_t NOpChannels = 0;
    NOpChannels = pvs->NOpChannels();


    G4MaterialPropertyVector* Fast_Intensity =
      aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
    G4MaterialPropertyVector* Slow_Intensity =
      aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

    if (!Fast_Intensity && !Slow_Intensity )
      return 1;


    G4int nscnt = 1;
    if (Fast_Intensity && Slow_Intensity) nscnt = 2;


    G4int Num = 0;
    double YieldRatio = GetYieldRatio(aParticle,
                                      aMaterialPropertiesTable);



    double const xyz[3] = { x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm };
    float const* Visibilities = pvs->GetAllVisibilities(xyz);

    for (G4int scnt = 1; scnt <= nscnt; scnt++) {

      G4double ScintillationTime = 0.*CLHEP::ns;
      G4double ScintillationRiseTime = 0.*CLHEP::ns;
      G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;

      if (scnt == 1) {
        if (nscnt == 1) {
          if(Fast_Intensity){
            ScintillationTime   = aMaterialPropertiesTable->
              GetConstProperty("FASTTIMECONSTANT");
            if (fFiniteRiseTime) {
              ScintillationRiseTime = aMaterialPropertiesTable->
                GetConstProperty("FASTSCINTILLATIONRISETIME");
            }
            ScintillationIntegral =
              (G4PhysicsOrderedFreeVector*)((*fFastIntegralTable)(materialIndex));
          }
          if(Slow_Intensity){
            ScintillationTime   = aMaterialPropertiesTable->
              GetConstProperty("SLOWTIMECONSTANT");
            if (fFiniteRiseTime) {
              ScintillationRiseTime = aMaterialPropertiesTable->
                GetConstProperty("SLOWSCINTILLATIONRISETIME");
            }
            ScintillationIntegral =
              (G4PhysicsOrderedFreeVector*)((*fSlowIntegralTable)(materialIndex));
          }
        }
        else {
          if(YieldRatio==0)
            YieldRatio = aMaterialPropertiesTable->
              GetConstProperty("YIELDRATIO");

          if ( fExcitationRatio == 1.0 ) {
            Num = G4int (std::min(YieldRatio,1.0)*MeanNumberOfPhotons);
          }
          else {
            Num = G4int (std::min(fExcitationRatio,1.0)*MeanNumberOfPhotons);
          }
          ScintillationTime   = aMaterialPropertiesTable->
            GetConstProperty("FASTTIMECONSTANT");
          if (fFiniteRiseTime) {
            ScintillationRiseTime = aMaterialPropertiesTable->
              GetConstProperty("FASTSCINTILLATIONRISETIME");
          }
          ScintillationIntegral =
            (G4PhysicsOrderedFreeVector*)((*fFastIntegralTable)(materialIndex));
        }
      }

      else {
        Num = MeanNumberOfPhotons - Num;
        ScintillationTime   =   aMaterialPropertiesTable->
          GetConstProperty("SLOWTIMECONSTANT");
        if (fFiniteRiseTime) {
          ScintillationRiseTime = aMaterialPropertiesTable->
            GetConstProperty("SLOWSCINTILLATIONRISETIME");
        }
        ScintillationIntegral =
          (G4PhysicsOrderedFreeVector*)((*fSlowIntegralTable)(materialIndex));
      }

      if (!ScintillationIntegral) continue;

      static long NN = 0;
      NN += Num;
      std::cout << " **** NUM = " << NN << " ****" << std::endl;

      if(!Visibilities){
      }else{
        std::map<int, int> DetectedNum;
        for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++){
          G4int DetThisPMT = G4int(G4Poisson(Visibilities[OpDet] * Num));
          if(DetThisPMT>0){
            DetectedNum[OpDet]=DetThisPMT;
            //   mf::LogInfo("OpFastScintillation") << "FastScint: " <<
            //   //   it->second<<" " << Num << " " << DetThisPMT;
          }
        }
        // Now we run through each PMT figuring out num of detected photons

        if(lgp->UseLitePhotons()){
          std::map<int, std::map<int, int>> StepPhotonTable;
          // And then add these to the total collection for the event
          for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
              itdetphot!=DetectedNum.end(); ++itdetphot){
            std::map<int, int>  StepPhotons;
            for (G4int i = 0; i < itdetphot->second; ++i){
              G4double deltaTime = aStep.GetStepLength() /
                ((pPreStepPoint->GetVelocity()+ pPostStepPoint->GetVelocity())/2.);


              if (ScintillationRiseTime==0.0) {
                deltaTime = deltaTime -
                  ScintillationTime * std::log( G4UniformRand() );
              }
              else {
                deltaTime = deltaTime +
                  sample_time(ScintillationRiseTime, ScintillationTime);
              }

              G4double aSecondaryTime = t0 + deltaTime;
              float Time = aSecondaryTime;
              int ticks = static_cast<int>(Time);
              StepPhotons[ticks]++;
            }
            StepPhotonTable[itdetphot->first] = StepPhotons;
            //Iterate over Step Photon Table to add photons to OpDetBacktrackerRecords.

            sim::OpDetBacktrackerRecord tmpOpDetBTRecord(itdetphot->first);
            int thisG4TrackID = (aStep.GetTrack())->GetTrackID();
            CLHEP::Hep3Vector prePoint  = (aStep.GetPreStepPoint())->GetPosition();
            CLHEP::Hep3Vector postPoint = (aStep.GetPostStepPoint())->GetPosition();
            //Note the use of xO (letter O) instead of x0. This is to differentiate the positions here with the earlier declared double* x0
            double xO = ( ( (prePoint.getX() + postPoint.getX() ) / 2.0) / CLHEP::cm );
            double yO = ( ( (prePoint.getY() + postPoint.getY() ) / 2.0) / CLHEP::cm );
            double zO = ( ( (prePoint.getZ() + postPoint.getZ() ) / 2.0) / CLHEP::cm );
            double const xyzPos[3] = {xO,yO,zO};
            double energy  = ( aStep.GetTotalEnergyDeposit() / CLHEP::GeV );

            //Loop over StepPhotons to get number of photons detected at each time for this channel and G4Step.
            for(std::map<int,int>::iterator stepPhotonsIt = StepPhotons.begin(); stepPhotonsIt != StepPhotons.end(); ++stepPhotonsIt){
              int photonTime = stepPhotonsIt->first;
              int numPhotons = stepPhotonsIt->second;
              tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, photonTime, numPhotons, xyzPos, energy);
            }
            //Add OpDetBackTrackerRecord. (opdetphotonTABLE->instance().addOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord BTRrecord)
            litefst->AddOpDetBacktrackerRecord(tmpOpDetBTRecord);
          }
          litefst->AddPhoton(&StepPhotonTable);
        }
        else{
          // And then add these to the total collection for the event
          for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
              itdetphot!=DetectedNum.end(); ++itdetphot){
            for (G4int i = 0; i < itdetphot->second; ++i){
              G4double deltaTime = aStep.GetStepLength() /
                ((pPreStepPoint->GetVelocity()+
                  pPostStepPoint->GetVelocity())/2.);

              if (ScintillationRiseTime==0.0){
                deltaTime = deltaTime -
                  ScintillationTime * std::log( G4UniformRand() );
              }
              else{
                deltaTime = deltaTime +
                  sample_time(ScintillationRiseTime, ScintillationTime);
              }

              G4double aSecondaryTime = t0 + deltaTime;

              // The sim photon in this case stores its production point and time
              TVector3 PhotonPosition(x0[0],x0[1],x0[2]);

              // We don't know anything about the momentum dir, so set it to be Z
              float Energy = 9.7*CLHEP::eV;
              float Time = aSecondaryTime;

              // Make a photon object for the collection
              sim::OnePhoton PhotToAdd;
              PhotToAdd.InitialPosition  = PhotonPosition;
              PhotToAdd.Energy           = Energy;
              PhotToAdd.Time             = Time;
              PhotToAdd.SetInSD          = false;

              fst->AddPhoton(itdetphot->first, std::move(PhotToAdd));
            }
          }
        }
      }
    }

    return 0;
  }

  // --------------------------------------------------
  double OpFastScintillation::
  GetYieldRatio(const G4DynamicParticle* aParticle,
                G4MaterialPropertiesTable* aMaterialPropertiesTable) const
  {
    if(!fScintillationByParticleType) return 0;

    // The scintillation response is a function of the energy
    // deposited by particle types.

    // Get the definition of the current particle
    G4ParticleDefinition* pDef = aParticle->GetDefinition();

    // Default for particles not enumerated below
    std::string propName = "ELECTRON";

    /* */if(pDef == G4Proton::ProtonDefinition()) propName = "PROTON";
    else if(pDef == G4MuonPlus::MuonPlusDefinition() ||
            pDef == G4MuonMinus::MuonMinusDefinition()) propName = "MUON";
    else if(pDef == G4PionPlus::PionPlusDefinition() ||
            pDef == G4PionMinus::PionMinusDefinition()) propName = "PION";
    else if(pDef == G4KaonPlus::KaonPlusDefinition() ||
            pDef == G4KaonMinus::KaonMinusDefinition()) propName = "KAON";
    else if(pDef == G4Alpha::AlphaDefinition()) propName = "ALPHA";
    // Electrons (must also account for shell-binding energy
    // attributed to gamma from standard PhotoElectricEffect)
    else if(pDef == G4Electron::ElectronDefinition() ||
	    pDef == G4Gamma::GammaDefinition()) propName = "ELECTRON";


    const double ret = aMaterialPropertiesTable->GetConstProperty((propName+"YIELDRATIO").c_str());

    // If the user has not specified yields for (p,d,t,a,carbon)
    // then these unspecified particles will default to the
    // electron's scintillation yield
    if(ret == 0) return aMaterialPropertiesTable->GetConstProperty("ELECTRONYIELDRATIO");

    return ret;
  }

  // --------------------------------------------------
  void OpFastScintillation::BuildPhysicsTable()
  {
    if (fFastIntegralTable && fSlowIntegralTable) return;

    const G4MaterialTable* theMaterialTable =
      G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    // create new physics table

    if(!fFastIntegralTable)fFastIntegralTable = new G4PhysicsTable(numOfMaterials);
    if(!fSlowIntegralTable)fSlowIntegralTable = new G4PhysicsTable(numOfMaterials);

    // loop for materials
    for (G4int i=0 ; i < numOfMaterials; i++){
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();
      G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();

      // Retrieve vector of scintillation wavelength intensity for
      // the material from the material's optical properties table.

      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
        aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {

        G4MaterialPropertyVector* fFastLightVector =
          aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");

        if (fFastLightVector) {

          // Retrieve the first intensity point in vector
          // of (photon energy, intensity) pairs

          G4double currentIN = (*fFastLightVector)[0];

          if (currentIN >= 0.0) {

            // Create first (photon energy, Scintillation
            // Integral pair

            G4double currentPM = fFastLightVector->Energy(0);

            G4double currentCII = 0.0;

            aPhysicsOrderedFreeVector->
              InsertValues(currentPM , currentCII);

            // Set previous values to current ones prior to loop

            G4double prevPM  = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN  = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material

            for (size_t i = 1;
                 i < fFastLightVector->GetVectorLength();
                 i++)
              {
                currentPM = fFastLightVector->Energy(i);
                currentIN = (*fFastLightVector)[i];

                currentCII = 0.5 * (prevIN + currentIN);

                currentCII = prevCII +
                  (currentPM - prevPM) * currentCII;

                aPhysicsOrderedFreeVector->
                  InsertValues(currentPM, currentCII);

                prevPM  = currentPM;
                prevCII = currentCII;
                prevIN  = currentIN;
              }
          }
        }

        G4MaterialPropertyVector* theSlowLightVector =
          aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

        if (theSlowLightVector) {

          // Retrieve the first intensity point in vector
          // of (photon energy, intensity) pairs

          G4double currentIN = (*theSlowLightVector)[0];

          if (currentIN >= 0.0) {

            // Create first (photon energy, Scintillation
            // Integral pair

            G4double currentPM = theSlowLightVector->Energy(0);

            G4double currentCII = 0.0;

            bPhysicsOrderedFreeVector->
              InsertValues(currentPM , currentCII);

            // Set previous values to current ones prior to loop

            G4double prevPM  = currentPM;
            G4double prevCII = currentCII;
            G4double prevIN  = currentIN;

            // loop over all (photon energy, intensity)
            // pairs stored for this material

            for (size_t i = 1;
                 i < theSlowLightVector->GetVectorLength();
                 i++)
              {
                currentPM = theSlowLightVector->Energy(i);
                currentIN = (*theSlowLightVector)[i];

                currentCII = 0.5 * (prevIN + currentIN);

                currentCII = prevCII +
                  (currentPM - prevPM) * currentCII;

                bPhysicsOrderedFreeVector->
                  InsertValues(currentPM, currentCII);

                prevPM  = currentPM;
                prevCII = currentCII;
                prevIN  = currentIN;
              }
          }
        }
      }

      // The scintillation integral(s) for a given material
      // will be inserted in the table(s) according to the
      // position of the material in the material table.
      fFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
      fSlowIntegralTable->insertAt(i,bPhysicsOrderedFreeVector);
    }
  }

  // --------------------------------------------------
  // Called by the user to set the scintillation yield as a function
  // of energy deposited by particle types
  void OpFastScintillation::SetScintillationByParticleType(const G4bool scintType)
  {
    if (fEMSaturation) {
      G4Exception("OpFastScintillation::SetScintillationByParticleType", "Scint02",
                  JustWarning, "Redefinition: Birks Saturation is replaced by ScintillationByParticleType!");
      RemoveSaturation();
    }
    fScintillationByParticleType = scintType;
  }

  // --------------------------------------------------
  G4double OpFastScintillation::GetMeanFreePath(const G4Track&,
                                                G4double ,
                                                G4ForceCondition* condition)
  {
    *condition = StronglyForced;

    return DBL_MAX;
  }

  // --------------------------------------------------
  G4double OpFastScintillation::GetMeanLifeTime(const G4Track&,
                                                G4ForceCondition* condition)
  {
    *condition = Forced;

    return DBL_MAX;
  }

  // --------------------------------------------------
  G4double OpFastScintillation::sample_time(G4double tau1, G4double tau2)
  {
    // tau1: rise time and tau2: decay time

    while(1) {
      // two random numbers
      G4double ran1 = G4UniformRand();
      G4double ran2 = G4UniformRand();
      //
      // exponential distribution as envelope function: very efficient
      //
      G4double d = (tau1+tau2)/tau2;
      // make sure the envelope function is
      // always larger than the bi-exponential
      G4double t = -1.0*tau2*std::log(1-ran1);
      G4double g = d*single_exp(t,tau2);
      if (ran2 <= bi_exp(t,tau1,tau2)/g) return t;
    }
    return -1.0;
  }

} // namespace
