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

    Int_t isrun;
    Int_t isevent;
    Int_t issubrun;
    int steps;
    bool isprim;
    Double_t dLdxStep[200];
    Double_t zStep[200];	//<---number of protons per length/step
    Double_t dLdEStep[200];	//<---number of protons per edep/step
	Double_t dLdIStep[200];
	Int_t isions;
	Int_t isions_allpart;
	Int_t isphotons_allpart;
    double size_max;
    double size_sum;
    int photons_sum;
	int ions_sum;
    double energy_sum;


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
    , fPhotonsPerStep(0)
    , fEnergyPerStep(0)
    , fElectronsVsPhotons(0)
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
    evnr=0;
    steps=0;
    isphotons=0;
	isions=0;
    isphotons_allpart=0;
	isions_allpart=0;
    size_max=0.5;
    size_sum=0.0;
    photons_sum=0.0;
    energy_sum=0.0;
;
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
fisTree = tfs->make<TTree>("IonScintTree","IonScintTree");
    fisTree->Branch("isrun", &isrun1, "isrun/I");
    fisTree->Branch("issubrun", &issubrun1, "issubrun/I");
    fisTree->Branch("isevent", &isevent1, "isevent/I");
    fisTree->Branch("isphotons", &isphotons, "isphotons/I");
    fisTree->Branch("isphotons_allpart", &isphotons_allpart, "isphotons_allpart/I");
    fisTree->Branch("isions_allpart", &isions_allpart, "isions_allpart/I");
    fisTree->Branch("dLdxStep", &dLdxStep, "dLdxStep[200]/D");
   fisTree->Branch("dLdEStep", &dLdEStep, "dLdEStep[200]/D");
    fisTree->Branch("isions", &isions, "isions/I");
   fisTree->Branch("dLdIStep", &dLdIStep, "dLdIStep[200]/D");





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
    if(fTrkID==1) isprim=true;
    else isprim=false;
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
      
      // Fill the histograms
      fStepSize          ->Fill(totstep.mag()/CLHEP::cm);
      fEnergyPerStep     ->Fill(fISCalc->EnergyDeposit());
      fElectronsPerStep  ->Fill(fISCalc->NumberIonizationElectrons());
      fPhotonsPerStep    ->Fill(fISCalc->NumberScintillationPhotons());
      fElectronsVsPhotons->Fill(fISCalc->NumberScintillationPhotons(), 
				fISCalc->NumberIonizationElectrons());
      fElectronsPerLength->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/(totstep.mag()/CLHEP::cm));
      fPhotonsPerLength  ->Fill(fISCalc->NumberScintillationPhotons()*1.e-3/(totstep.mag()/CLHEP::cm));
      fElectronsPerEDep  ->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/fISCalc->EnergyDeposit());
      fPhotonsPerEDep    ->Fill(fISCalc->NumberScintillationPhotons()*1.e-3/fISCalc->EnergyDeposit());
	if(isprim) isphotons+=fISCalc->NumberScintillationPhotons();
	if(isprim) isions+=fISCalc->NumberIonizationElectrons();
	isphotons_allpart+=fISCalc->NumberScintillationPhotons();
	isions_allpart+=fISCalc->NumberIonizationElectrons();
mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction "<<fTrkID<<" "<<step->GetTrack()->GetDefinition()->GetParticleName()<<" "<<fStepNumber<<std::endl;
    if(isprim ){ //only primary particles are interesting in this step of the the NEST validation
		size_sum+=double(totstep.mag()/CLHEP::cm);
		photons_sum+=int(fISCalc->NumberScintillationPhotons());
		energy_sum+=double(fISCalc->EnergyDeposit()/CLHEP::MeV);
		ions_sum+=int(fISCalc->NumberIonizationElectrons());
		if(size_sum>0.5&&steps<200){
		dLdxStep[steps]=Double_t(photons_sum/size_sum);
		dLdIStep[steps]=Double_t(ions_sum/size_sum);
		dLdEStep[steps]=Double_t(photons_sum/energy_sum);
		zStep[steps]=Double_t(fStep->GetPreStepPoint()->GetPosition().z()/CLHEP::cm);
		++steps;
		size_sum=0.0;
		photons_sum=0.0;
		
		mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction "<<fTrkID<<" "<<fStepNumber<<" "<<steps<<" "<<dLdxStep[steps]<<std::endl;
		}
	}//end of statement for storing primary photons/e dep
	change_event=false;//>new track	change_event=false;//>new event

	if(isevent!=0){
		if(isevent!=evnr){
			change_event=true;
			evnr=isevent;
			/*for(int ii=0;ii<steps;++ii){
				mf::LogWarning("IonisationScintillation")<<" photons stepby step "<<ii<<" "<<dLdxStep[ii]<<std::endl;

				}*/
			steps=0;
			mf::LogWarning("ionistation scint")<<"--------------EVENT NUMBER"<<evnr<<std::endl;
		}
		

	}

	else{
		change_event=true;
		/*for(int ii=0;ii<steps;++ii){
			mf::LogWarning("IonisationScintillation")<<" photons step by step "<<ii<<" "<<dLdxStep[ii]<<std::endl;

		}*/
		steps=0;
		evnr=isevent;
		//mf::LogWarning("IonisationScintillation")<<" EVENT NUMBER scint energyaction "<<evnr<<std::endl;	
	}
	if (change_event==true){

		isevent1=isevent;
		isrun1=isrun;
		issubrun1=issubrun;
mf::LogWarning("IonisationScintillation")<<" EVENT NUMBER scint energyaction "<<isevent1<<std::endl;
if(steps<200&&isprim) mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction written to the tree "<<fStepNumber<<" "<<dLdxStep[steps]<<" "<<dLdxStep[1]<<std::endl;
		fisTree->Fill();
		isphotons=0;
		isions=0;
		size_sum=0.0;
	}

    } // end if the energy deposition is non-zero
else{		if(steps<200&&isprim){
		dLdxStep[steps]=0.0;
		dLdEStep[steps]=0.0;
		dLdIStep[steps]=0.0;
		mf::LogWarning("IonisationScintillation")<<" trk id scint zero energyaction "<<fTrkID<<" "<<fStepNumber<<" "<<steps<<" "<<dLdxStep[steps]<<std::endl;
		}

}
	//if(isprim) ++steps;
    return;
  }

} // namespace
