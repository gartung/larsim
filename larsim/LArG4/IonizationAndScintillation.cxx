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
        int steps1;
	    int steps2;
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
Int_t isphotons_allpartEn;
Int_t isphotonsEn;
float photons_count;
float energy_dep;
Double_t Ene;
Double_t Ene_all;  
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
    steps1=0;
    steps2=0;
    isphotons=0;
    isions=0;
    isphotons_allpart=0;
    isions_allpart=0;
    size_max=0.5;
    size_sum=0.0;
    photons_sum=0.0;
    energy_sum=0.0;
    isphotonsEn=0;
    isphotons_allpartEn=0;
    Ene=0.0;
    Ene_all=0.0;

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

    fEnergyVsPhotons = tfs->make<TH2F>("energyVsPhotons", ";Photons;Energy",
                                          500, 0., 5000., 100, 0., 0.5);


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
    fisTree->Branch("isphotons_allpartEn", &isphotons_allpartEn, "isphotons_allpartEn/I");
    fisTree->Branch("isphotonsEn", &isphotonsEn, "isphotonsEn/I");    
fisTree->Branch("Ene", &Ene, "Ene/D");
fisTree->Branch("Ene_all", &Ene_all, "Ene_all/D");




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
    energy_dep=0.0;
    photons_count=0;
    // reset the calculator
    fISCalc->Reset();

    // check the material for this step and be sure it is LAr
    if(step->GetTrack()->GetMaterial()->GetName() != "LAr") return;

    // double check that the energy deposit is non-zero
    // then do the calculation if it is
    if( step->GetTotalEnergyDeposit() > 0 ){
      //getting number of photons and energy deposit in this step

      fISCalc->CalculateIonizationAndScintillation(fStep);
    photons_count=fISCalc->NumberScintillationPhotons(); 
      LOG_DEBUG("IonizationAndScintillation") << "Step Size: "   << fStep->GetStepLength()/CLHEP::cm
					      << "\nEnergy: "    << energy_dep
					      << "\nElectrons: " << fISCalc->NumberIonizationElectrons()
					      << "\nPhotons: "   << photons_count;

      G4ThreeVector totstep = fStep->GetPostStepPoint()->GetPosition();
      totstep -= fStep->GetPreStepPoint()->GetPosition();
      
      // Fill the histograms
      fStepSize          ->Fill(totstep.mag()/CLHEP::cm);
      fEnergyPerStep     ->Fill(energy_dep);
      fElectronsPerStep  ->Fill(fISCalc->NumberIonizationElectrons());
      fPhotonsPerStep    ->Fill(photons_count);
      fElectronsVsPhotons->Fill(photons_count, 
				fISCalc->NumberIonizationElectrons());
      fEnergyVsPhotons->Fill(photons_count, 
                                energy_dep);
 
     fElectronsPerLength->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/(totstep.mag()/CLHEP::cm));
      fPhotonsPerLength  ->Fill(photons_count*1.e-3/(totstep.mag()/CLHEP::cm));
      fElectronsPerEDep  ->Fill(fISCalc->NumberIonizationElectrons()*1.e-3/energy_dep);
      fPhotonsPerEDep    ->Fill(photons_count*1.e-3/energy_dep);

	isphotons_allpart+=int(photons_count);
		energy_dep=fISCalc->EnergyDeposit()/CLHEP::MeV;
	
	isions_allpart+=int(fISCalc->NumberIonizationElectrons());
        ++steps1;
	Ene_all+=energy_dep;
mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction "<<fTrkID<<" "<<step->GetTrack()->GetDefinition()->GetParticleName()<<" "<<fStepNumber<<std::endl;
                              if(Double_t(energy_dep)!=0.0)	isphotons_allpartEn+=Int_t(photons_count); 
   if(isprim ){ //only primary particles are interesting in this step of the the NEST validation
		isphotons+=photons_count;
		isions+=fISCalc->NumberIonizationElectrons();
		size_sum+=double(totstep.mag()/CLHEP::cm);
		photons_sum+=int(photons_count);
		energy_sum+=double(energy_dep);
		Ene+=energy_dep;
		ions_sum+=int(fISCalc->NumberIonizationElectrons());
		if(Double_t(energy_dep)!=0.0) isphotonsEn+=Int_t(fISCalc->NumberScintillationPhotons());
		
		if(size_sum>0.5&&steps<200){
		dLdxStep[steps]=Double_t(photons_sum/size_sum);
		dLdIStep[steps]=Double_t(ions_sum/size_sum);
		dLdEStep[steps]=Double_t(photons_sum/energy_sum);
		zStep[steps]=Double_t(fStep->GetPreStepPoint()->GetPosition().z()/CLHEP::cm);
		++steps;
                size_sum=0.0;
		
                photons_sum=0.0;
		++steps2;//steps of primary particle
		mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction photons test "<<fTrkID<<" "<<fStepNumber<<" "<<steps<<" "<<isphotons<<" added --> "<<photons_count<<" should be "<<fISCalc->NumberScintillationPhotons()<<" to integer "<<int(fISCalc->NumberScintillationPhotons())<<" int_t "<<Int_t(fISCalc->NumberScintillationPhotons())<<std::endl;
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
		//	steps=0;
			//steps1=0;
			//	steps2=0;
			mf::LogWarning("ionistation scint")<<"--------------EVENT NUMBER"<<evnr<<std::endl;
		}
		

	}

	else{
		change_event=true;
		evnr=isevent;
	}
	if (change_event==true){

		isevent1=isevent;
		isrun1=isrun;
		issubrun1=issubrun;
		mf::LogWarning("IonisationScintillation")<<" EVENT NUMBER scint energyaction "<<isevent1<<" photons count "<<isphotons<<" all particles: "<<isphotons_allpart<<" different method "<<isphotonsEn<<std::endl;
		if(steps<200&&isprim) mf::LogWarning("IonisationScintillation")<<" trk id scint energyaction written to the tree "<<fStepNumber<<" "<<dLdxStep[steps]<<" "<<dLdxStep[1]<<std::endl;
		
		

		fisTree->Fill();
		isphotons=0;
		isions=0;
		size_sum=0.0;
		isphotons_allpart=0;
		isphotonsEn=0;
		isphotons_allpartEn=0;
		Ene=0.0;
		Ene_all=0.0;
		steps=0;
		steps1=0;
		steps2=0;
	}

    } // end if the energy deposition is non-zero
else{		if(steps<200&&isprim){
		dLdxStep[steps]=0.0;
		dLdEStep[steps]=0.0;
		dLdIStep[steps]=0.0;
		mf::LogWarning("IonisationScintillation")<<" trk id scint zero energyaction "<<fTrkID<<" "<<fStepNumber<<" "<<steps<<" "<<dLdxStep[steps]<<std::endl;
		}

}
    return;
  }

} // namespace
