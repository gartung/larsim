////////////////////////////////////////////////////////////////////////
/// \file IonizationAndScintillation.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// lar includes
#include "larsim/LArG4/IonizationAndScintillation.h"
#include "larsim/LArG4/ISCalculationNEST.h"
#include "larsim/LArG4/ISCalculationSeparate.h"
#include "larsim/Simulation/LArG4Parameters.h"

// ROOT includes

// Framework includes
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// C/C++ standard libraries
#include <cassert>

namespace larg4 {

  static IonizationAndScintillation* gInstance = 0;

  //......................................................................
  IonizationAndScintillation* IonizationAndScintillation::CreateInstance
    (CLHEP::HepRandomEngine& engine)
  {
    if(!gInstance) gInstance = new IonizationAndScintillation(engine);
    return gInstance;
  }

  //......................................................................
  IonizationAndScintillation* IonizationAndScintillation::Instance()
  {
    // the instance must have been created already by CreateInstance()
    assert(gInstance);
    return gInstance;
  }

  //......................................................................
  // Constructor.
  IonizationAndScintillation::IonizationAndScintillation
    (CLHEP::HepRandomEngine& engine)
    : fISCalc(0)
    , fStep(0)
    , fElectronsPerStep(0)
    , fStepSize(0)
    , fStepSizeCompare(0)
    , fPhotonsPerStep(0)
    , fEnergyPerStep(0)
    , fElectronsVsPhotons(0)
    , fRecombSurvivalFrac(0)
    , fEnergyPerLength(0)
    , fdEdxVsRecomb(0)
    , fEngine(engine)
  {
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    fISCalculator = lgp->IonAndScintCalculator();

    if(fISCalculator.compare("NEST") == 0)
      fISCalc = new larg4::ISCalculationNEST(fEngine);
    else if(fISCalculator.compare("Separate") == 0)
      fISCalc = new larg4::ISCalculationSeparate(fEngine);
    else
      mf::LogWarning("IonizationAndScintillation") << "No ISCalculation set, this can't be good.";

    // Reset the values for the electrons, photons, and energy to 0
    // in the calculator
    fISCalc->Reset();
    //set the current track and step number values to bogus so that it will run the first reset:
    fStepNumber=-1;
    fTrkID=-1;

    // initialize the calculator
    fISCalc->Initialize();

    // make the histograms
    art::ServiceHandle< art::TFileService> tfs;

    fElectronsPerStep   = tfs->make<TH1F>("electronsPerStep", ";Electrons;Steps", 
					  500, 0., 5000.);			
    fPhotonsPerStep   	= tfs->make<TH1F>("photonsPerStep", ";Photons;Steps", 	
					  500, 0., 5000.);			
    fEnergyPerStep    	= tfs->make<TH1F>("energyPerStep", ";Energy (MeV);Steps", 
					  100, 0., 0.5);				
    fStepSize         	= tfs->make<TH1F>("stepSize", ";Step Size (CLHEP::cm);Steps", 	
					  500, 0., 0.2);                          
    fStepSizeCompare    = tfs->make<TH2F>("stepSizeCompare", ";G4Step->GetStepLength();PostStepPoint-PreStepPoint",
                                          200,0.,0.2,
                                          200,0.,0.2);
    fElectronsPerLength = tfs->make<TH1F>("electronsPerLength", ";Electrons #times 10^{3}/CLHEP::cm;Steps",
					  1000, 0., 1000.);
    fPhotonsPerLength   = tfs->make<TH1F>("photonsPerLength", ";Photons #times 10^{3}/CLHEP::cm;Steps",
					  1000, 0., 1000.);
    fElectronsPerEDep   = tfs->make<TH1F>("electronsPerEDep", ";Electrons #times 10^{3}/MeV;Steps",
					  1000, 0., 1000.);
    fPhotonsPerEDep     = tfs->make<TH1F>("photonsPerEDep", ";Photons #times 10^{3}/MeV;Steps",
					  1000, 0., 1000.);
					  
    fElectronsVsPhotons = tfs->make<TH2F>("electronsVsPhotons", ";Photons;Electrons",
					  500, 0., 5000., 500, 0., 5000.);
    
    fRecombSurvivalFrac = tfs->make<TH1F>("recombSurvivalFrac", ";Recombination survival fraction (N_{e}/N_{i})",
					  1000, 0., 1.2);
    fRecombSurvivalFrac_el = tfs->make<TH1F>("recombSurvivalFrac_el", ";Recombination survival fraction (electrons hower) (N_{e}/N_{i})",
					  1000, 0., 1.2);
    
    fEnergyPerLength = tfs->make<TH1F>("energyPerLength", "Step size >= 0.03 cm;dE/dx [MeV/cm]",
					  300, 0., 30.);
      fEnergyPerLength->SetOption("hist");    
    fEnergyPerLength_el = tfs->make<TH1F>("energyPerLength_el", "Step size >= 0.03 cm;dE/dx [MeV/cm]",
					  300, 0., 30.);
      fEnergyPerLength_el->SetOption("hist");    
    
    fEnergyPerLength_mu = tfs->make<TH1F>("energyPerLength_mu", "Step size >= 0.03 cm;dE/dx [MeV/cm]",
					  300, 0., 30.);
      fEnergyPerLength_mu->SetOption("hist");    
    
    fdEdxVsRecomb       = tfs->make<TH2F>("dEdxVsRecomb", ";dE/dx [MeV/cm];Recombination survival fraction",
					  1000, 0., 50.,
                                          1000, 0., 1.0);

    fdEdxVsStepLength   = tfs->make<TH2F>("dEdxVsStepLength",";dE/dx [MeV/cm];Step length in Geant4 [cm]",
                                          1000,0.,50.,
                                          1000,0.,0.05);
    fdEdxVsStepLength->SetOption("colz");
    
    fRecomb_dEdxBins_2   = tfs->make<TH1F>("Recomb_dEdxBins_2","2 MeV/cm;Recombination survival fraction",
                                          1000,0.,1.0);
    fRecomb_dEdxBins_5   = tfs->make<TH1F>("Recomb_dEdxBins_5","5 MeV/cm;Recombination survival fraction",
                                          1000,0.,1.0);
    fRecomb_dEdxBins_15   = tfs->make<TH1F>("Recomb_dEdxBins_15","15 MeV/cm;Recombination survival fraction",
                                          1000,0.,1.0);
    fRecomb_dEdxBins_30   = tfs->make<TH1F>("Recomb_dEdxBins_30","30 MeV/cm;Recombination survival fraction",
                                          1000,0.,1.0);
    fdQdx_dEdxBins_2   = tfs->make<TH1F>("dQdx_dEdxBins_2","2 MeV/cm;dQ/dx [e-/cm]",
                                          1000,0.,500e3);
    fdQdx_dEdxBins_5   = tfs->make<TH1F>("dQdx_dEdxBins_5","5 MeV/cm;dQ/dx [e-/cm]",
                                          1000,0.,500e3);
    fdQdx_dEdxBins_15   = tfs->make<TH1F>("dQdx_dEdxBins_15","15 MeV/cm;dQ/dx [e-/cm]",
                                          1000,0.,500e3);
    fdQdx_dEdxBins_30   = tfs->make<TH1F>("dQdx_dEdxBins_30","30 MeV/cm;dQ/dx [e-/cm]",
                                          1000,0.,500e3);


    return;
  }

  //......................................................................
  IonizationAndScintillation::~IonizationAndScintillation() 
  {
    if(fISCalc) delete fISCalc;
  }


  //......................................................................
  void IonizationAndScintillation::Reset(const G4Step* step)
  {

    if(fStepNumber==step->GetTrack()->GetCurrentStepNumber() && fTrkID==step->GetTrack()->GetTrackID())
      return;
    
    fStepNumber=step->GetTrack()->GetCurrentStepNumber(); 
    fTrkID=step->GetTrack()->GetTrackID();

    fStep = step;

    // reset the calculator
    fISCalc->Reset();

    // check the material for this step and be sure it is LAr
    if(step->GetTrack()->GetMaterial()->GetName() != "LAr") return;

    // double check that the energy deposit is non-zero
    // then do the calculation if it is
    if( step->GetTotalEnergyDeposit() > 0 ){
 
      fISCalc->CalculateIonizationAndScintillation(fStep);
    
      LOG_DEBUG("IonizationAndScintillation") << "Step Size: "   << fStep->GetStepLength()/CLHEP::cm
					      << "\nEnergy: "    << fISCalc->EnergyDeposit()
					      << "\nElectrons: " << fISCalc->NumberIonizationElectrons()
					      << "\nPhotons: "   << fISCalc->NumberScintillationPhotons();

      G4ThreeVector totstep = fStep->GetPostStepPoint()->GetPosition();
      totstep -= fStep->GetPreStepPoint()->GetPosition();
      double dx     = totstep.mag()/CLHEP::cm;
      //double dx     = fStep->GetStepLength()/CLHEP::cm;
     
      
      // Fill the histograms
      fStepSize          ->Fill(dx);
      fStepSizeCompare   ->Fill( fStep->GetStepLength()/CLHEP::cm, totstep.mag()/CLHEP::cm );
      fEnergyPerStep     ->Fill(fISCalc->EnergyDeposit());
      fElectronsPerStep  ->Fill(fISCalc->NumberIonizationElectrons());
      fPhotonsPerStep    ->Fill(fISCalc->NumberScintillationPhotons());
      fElectronsVsPhotons->Fill(fISCalc->NumberScintillationPhotons(), 
				fISCalc->NumberIonizationElectrons());
      if( fISCalc->Recomb() > 0. ){
        if( dx >= 0.03 ) fEnergyPerLength  ->Fill( fISCalc->dEdx() );
        if( fStep->GetTrack()->GetTrackID() == 1 ) {
          fdEdxVsStepLength->Fill( fISCalc->dEdx(), dx );
          if( dx >= 0.03 ) fEnergyPerLength_mu->Fill( fISCalc->dEdx() );
        }
        if( fStep->GetTrack()->GetTrackID() > 1 ) {
          if( dx >= 0.03 ) fEnergyPerLength_el->Fill( fISCalc->dEdx() );
          fRecombSurvivalFrac_el->Fill( fISCalc->Recomb() );
        }
        fRecombSurvivalFrac->Fill(fISCalc->Recomb(), fISCalc->EnergyDeposit() );
        fdEdxVsRecomb->Fill( fISCalc->dEdx(), fISCalc->Recomb() );
        if( fISCalc->dEdx() > 1.99 && fISCalc->dEdx() < 2.01 ) {
          fRecomb_dEdxBins_2  ->Fill( fISCalc->Recomb() );
          fdQdx_dEdxBins_2    ->Fill( fISCalc->NumberIonizationElectrons()/dx );
        } else 
        if( fISCalc->dEdx() > 4.98 && fISCalc->dEdx() < 5.02 ) {
          fRecomb_dEdxBins_5  ->Fill( fISCalc->Recomb() );
          fdQdx_dEdxBins_5    ->Fill( fISCalc->NumberIonizationElectrons()/dx );
        } else 
        if( fISCalc->dEdx() > 14.92 && fISCalc->dEdx() < 15.08 ) {
          fRecomb_dEdxBins_15 ->Fill( fISCalc->Recomb() );
          fdQdx_dEdxBins_15   ->Fill( fISCalc->NumberIonizationElectrons()/dx );
        } else 
        if( fISCalc->dEdx() > 29.85 && fISCalc->dEdx() < 30.15 ) {
          fRecomb_dEdxBins_30 ->Fill( fISCalc->Recomb() );
          fdQdx_dEdxBins_30   ->Fill( fISCalc->NumberIonizationElectrons()/dx );
        }
      }
//      double const stepSize = totstep.mag()/CLHEP::cm;
      double const stepSize = dx;
      if (stepSize > 0.0) {
        fElectronsPerLength->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/stepSize);
        fPhotonsPerLength  ->Fill(fISCalc->NumberScintillationPhotons()*1.e-3/stepSize);
      }
      double const energyDep = fISCalc->EnergyDeposit();
      if (energyDep) {
        fElectronsPerEDep  ->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/energyDep);
        fPhotonsPerEDep    ->Fill(fISCalc->NumberScintillationPhotons()*1.e-3/energyDep);
      }

    } // end if the energy deposition is non-zero

    return;
  }

} // namespace
