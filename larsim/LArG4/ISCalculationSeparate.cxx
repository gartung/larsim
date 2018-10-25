////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "CLHEP/Vector/ThreeVector.h"

#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4EmSaturation.hh"

#include "larsim/LArG4/ISCalculationSeparate.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Simulation/LArVoxelCalculator.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

namespace larg4{

  //----------------------------------------------------------------------------
  ISCalculationSeparate::ISCalculationSeparate(CLHEP::HepRandomEngine& engine)
  : fEngine(engine)
  {
  }

  //----------------------------------------------------------------------------
  ISCalculationSeparate::~ISCalculationSeparate()
  {
  }

  //----------------------------------------------------------------------------
  void ISCalculationSeparate::Initialize()
  {
    art::ServiceHandle<sim::LArG4Parameters> lgpHandle;
    const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    double density       = detprop->Density(detprop->Temperature());
    fEfield       = detprop->Efield();
    fEfieldNominal = fEfield;
    fScintByParticleType = larp->ScintByParticleType();
    fGeVToElectrons      = lgpHandle->GeVToElectrons();

    // \todo get scintillation yield from LArG4Parameters or LArProperties
    fScintYieldFactor  = 1.;

    // the recombination coefficient is in g/(MeVcm^2), but
    // we report energy depositions in MeV/cm, need to divide
    // Recombk from the LArG4Parameters service by the density
    // of the argon we got above.
    fRecombA             = lgpHandle->RecombA();
    fRecombk             = lgpHandle->Recombk()/density;
    fModBoxA             = lgpHandle->ModBoxA();
    fModBoxB             = lgpHandle->ModBoxB()/density;
    fUseModBoxRecomb     = lgpHandle->UseModBoxRecomb();  

    // Use Birks Correction in the Scintillation process    
    fEMSaturation = G4LossTableManager::Instance()->EmSaturation();

    // determine the step size using the voxel sizes
    art::ServiceHandle<sim::LArVoxelCalculator> lvc;
    double maxsize = std::max(lvc->VoxelSizeX(), std::max(lvc->VoxelSizeY(), lvc->VoxelSizeZ())) * CLHEP::cm;

    fStepSize = 0.1 * maxsize;

    fCachedTrackID = -1;

    return;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  void ISCalculationSeparate::Reset()
  {
    fEnergyDeposit   = 0.;
    fNumScintPhotons = 0.;
    fNumIonElectrons = 0.;
    fNumIons         = 0.;
    fRecomb          = 0.; 
    fdEdx            = 0.; 
    return;
  }

  //----------------------------------------------------------------------------
  // fNumIonElectrons returns a value that is not corrected for life time effects
  void ISCalculationSeparate::CalculateIonizationAndScintillation(const G4Step* step)
  {
    
    CLHEP::RandGauss GaussGen(fEngine);
    CLHEP::RandFlat  UniformGen(fEngine);
    bool simStochasticity = false;
    bool simAngDependence = false;
    bool simEfieldFluct   = false;

    fEnergyDeposit = step->GetTotalEnergyDeposit()/CLHEP::MeV;

    // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
    // R = A/(1 + (dE/dx)*k)
    // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
    // from GeV/voxel width
    // A = 0.800 +/- 0.003
    // k = (0.097+/-0.001) g/(MeVcm^2)
    // the dx depends on the trajectory of the step
    // k should be divided by the density as our dE/dx is in MeV/cm,
    // the division is handled in the constructor when we set fRecombk
    // B.Baller: Add Modified Box recombination - ArgoNeuT result submitted to JINST
   
    // Update the cached track ID. If the ID == 1, means new event has started
    int trkID = step->GetTrack()->GetTrackID();
    if( simEfieldFluct && trkID == 1 && fCachedTrackID != 1 )
      fEfield = GaussGen.fire(fEfieldNominal,fEfieldNominal*0.15); 
    
    fCachedTrackID = trkID;

    G4ThreeVector totstep = step->GetPostStepPoint()->GetPosition();
    totstep -= step->GetPreStepPoint()->GetPosition();
    double dx     = totstep.mag()/CLHEP::cm;
    //double dx     = step->GetStepLength()/CLHEP::cm;
    double recomb = 0.;
    double dEdx   = (dx == 0.0)? 0.0: fEnergyDeposit/dx;
    double EFieldStep = EFieldAtStep(fEfield,step);
    
    /* 
    const G4Material* aMaterial = step->GetPreStepPoint()->GetMaterial();
//    G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
    G4double Density = aMaterial->GetDensity()/(CLHEP::g/CLHEP::cm3);
    std::cout
    <<"Material (post-step) = "<<step->GetPostStepPoint()->GetMaterial()->GetName()
    <<", density = "<<Density
    <<", temp = "<<step->GetPostStepPoint()->GetMaterial()->GetTemperature()
    <<", press = "<<step->GetPostStepPoint()->GetMaterial()->GetPressure()<<"\n";
    */

    /*    
    // .........................................
          if( UniformGen.fire() < 0.25 ) dEdx = 2.;
    else  if( UniformGen.fire() < 0.50 ) dEdx = 5.;
    else  if( UniformGen.fire() < 0.75 ) dEdx = 15.;
    else  if( UniformGen.fire() < 1.00 ) dEdx = 30.;
    fEnergyDeposit = dx * dEdx;
    //............................................
    */

    //std::cout<<"Using EField = "<<EFieldStep<<"\n";

    if( dx ) {
      // At low dE/dx, Modified Box becomes unphysical and 
      // turns downward. So calculate R with both models, and
      // if dE/dx is low, use whichever is higher.
      double Xi = fModBoxB * dEdx / EFieldStep;
      double recomb_ModBox = log(fModBoxA + Xi) / Xi;
      double recomb_Birks = fRecombA/(1. + dEdx * fRecombk / EFieldStep);
      if( fUseModBoxRecomb && ( dEdx > 2.0 || recomb_ModBox >= recomb_Birks ) )
        recomb = recomb_ModBox; 
      else 
        recomb = recomb_Birks;
    }


    /* 
    if( fUseModBoxRecomb ) {
      // Guard against low dE/dx where modified box isn't applicable
      if(dEdx < 1.) dEdx = 1.;
      if (dx){
	double Xi = fModBoxB * dEdx / EFieldStep;
	recomb = log(fModBoxA + Xi) / Xi;
      }
      else 
	recomb = 0;
    } 
    else{
      recomb = fRecombA/(1. + dEdx * fRecombk / EFieldStep);
    }
    */
    
    
     
    // ------------------------------------------------------------------- 
    // Calculate angle of this step rel. to the E-field
    // which is assumed to be along x axis.
    if( simAngDependence && totstep.mag() > 0 ) {
      totstep.setMag(1.);
      double dy = totstep.getY();
      double dz = totstep.getZ();
      double arg = sqrt(dy*dy+dz*dz);
      double angle = 90.;
      if( arg < 1. ) angle = (180./3.14159)*std::asin( sqrt(dy*dy+dz*dz) );

      // Scale recombination up/down from the central value
      // depending on angle, using ArgoNeuT results, arxiv:1306.1712.
      // Study done on 6/4/2018 to find this parameterization using
      // values for the parameterized fits at 15 MeV/cm.
      //double a0 = 0.9534;
      //double a1 = 1.022e-4;
      //double a2 = 1.025e-5;
      //double recombScale = a0 + a1*angle + a2*pow(angle,2);
      // Revised study on 6/7/2018 using more realistic bin center on
      // lowest-angle bin, and fit linearly
      double a0 = 0.9116;
      double a1 = 0.001435;
      double recombScale = a0 + a1*angle;

      recomb *= recombScale;
      if( recomb > 1. ) recomb = 1.;
      if( recomb < 0. ) recomb = 0.;
    }
   
    fRecomb = recomb;
    fdEdx   = dEdx;

    // 1.e-3 converts fEnergyDeposit to GeV

      // First divide up energy into total number of "quanta" (ionizations+excitations),
      //    fGeVToElectrons = W_ion^-1
      //    W_ph            = W_ion / (1+excRatio)
      //    totalQuanta     = dE / W_ph
    double excRatio = 0.21;     // LAr excitation ratio, N_exc / N_i
    int totalQuanta = floor( fGeVToElectrons * 1.e-3 * fEnergyDeposit * (1.+excRatio) );

    if( simStochasticity ) {
      // --------------------------------------------
      // Here we make some modifications that follow along with what's done 
      // by the NEST alg to more naturally simulate the fluctuations in charge.
      double FanoFactor = 0.107;  // Fano factor for LAr (Doke)
      // Applying Fano factor smearing for LAr (sub-Poisson when factoring in energy 
      // conservation constraint). Note that 
      totalQuanta = int(floor(GaussGen.fire(totalQuanta, sqrt(FanoFactor*totalQuanta))+0.5));
      // Separate into excited and ionized populations
      int numExc = BinomFluct( totalQuanta, excRatio/(1.+excRatio) );
      int numIon = totalQuanta - numExc; 
      fNumIons = numIon;

      // Simulate the recombination with stochasticity by drawing from binomial
      fNumIonElectrons  = (double)BinomFluct( numIon, recomb );
    
    } else {

      fNumIonElectrons = fGeVToElectrons * 1.e-3 * fEnergyDeposit * recomb;
    
    }
      
    fNumScintPhotons = totalQuanta - fNumIonElectrons;

    //std::cout
    //<<" -------------------------------\n"
    //<<"  Efield             : "<<fEfield<<"\n"
    //<<"  dE                 : "<<fEnergyDeposit<<"\n"
    //<<"  dE/dx              : "<<dEdx<<"\n"
    //<<"  total quanta       : "<<totalQuanta<<"\n"
    //<<"  num ion            : "<<numIon<<"\n"
    //<<"  recomb:            : "<<recomb<<"\n"
    //<<"  num free electrons : "<<fNumIonElectrons<<"\n";
    // ------------------------------------------------------------------- 

    LOG_DEBUG("ISCalculationSeparate") << " Electrons produced for " << fEnergyDeposit 
				       << " MeV deposited with "     << recomb 
				       << " recombination: "         << fNumIonElectrons;

    // Now do the scintillation
    /*
    G4MaterialPropertiesTable* mpt = step->GetTrack()->GetMaterial()->GetMaterialPropertiesTable();
    if( !mpt) 
      throw cet::exception("ISCalculationSeparate") << "Cannot find materials property table"
						    << " for this step! "
						    << step->GetTrack()->GetMaterial() << "\n";

    // if not doing the scintillation by particle type, use the saturation
    double scintYield = mpt->GetConstProperty("SCINTILLATIONYIELD");

    if(fScintByParticleType){

      LOG_DEBUG("ISCalculationSeparate") << "scintillating by particle type";

      // Get the definition of the current particle
      G4ParticleDefinition *pDef = step->GetTrack()->GetDynamicParticle()->GetDefinition();
      //G4MaterialPropertyVector *Scint_Yield_Vector = NULL;

      // Obtain the G4MaterialPropertyVectory containing the
      // scintillation light yield as a function of the deposited
      // energy for the current particle type
      
      // Protons
      if(pDef == G4Proton::ProtonDefinition()){
	scintYield = mpt->GetConstProperty("PROTONSCINTILLATIONYIELD");
      }
      // Muons
      else if(pDef == G4MuonPlus::MuonPlusDefinition() ||
	      pDef == G4MuonMinus::MuonMinusDefinition()){
	scintYield = mpt->GetConstProperty("MUONSCINTILLATIONYIELD");
      }
      // Pions
      else if(pDef == G4PionPlus::PionPlusDefinition() ||
	      pDef == G4PionMinus::PionMinusDefinition()){
	scintYield = mpt->GetConstProperty("PIONSCINTILLATIONYIELD");
      }
      // Kaons
      else if(pDef == G4KaonPlus::KaonPlusDefinition() ||
	      pDef == G4KaonMinus::KaonMinusDefinition()){
	scintYield = mpt->GetConstProperty("KAONSCINTILLATIONYIELD");
      }
      // Alphas
      else if(pDef == G4Alpha::AlphaDefinition()){
	scintYield = mpt->GetConstProperty("ALPHASCINTILLATIONYIELD");
      }
      // Electrons (must also account for shell-binding energy
      // attributed to gamma from standard PhotoElectricEffect)
      else if(pDef == G4Electron::ElectronDefinition() ||
	      pDef == G4Gamma::GammaDefinition()){
	scintYield = mpt->GetConstProperty("ELECTRONSCINTILLATIONYIELD");
      }       
      // Default for particles not enumerated/listed above
      else{
	scintYield = mpt->GetConstProperty("ELECTRONSCINTILLATIONYIELD");
      }

      // If the user has not specified yields for (p,d,t,a,carbon)
      // then these unspecified particles will default to the 
      // electron's scintillation yield
      //if(!Scint_Yield_Vector){{
	//  scintYield = mpt->GetConstProperty("ELECTRONSCINTILLATIONYIELD");
	//}
      //}
	   
      // Throw an exception if no scintillation yield is found
      if (!scintYield) 
	throw cet::exception("ISCalculationSeparate") << "Request for scintillation yield for energy "
						      << "deposit and particle type without correct "
						      << "entry in MaterialPropertiesTable\n"
						      << "ScintillationByParticleType requires at "
						      << "minimum that ELECTRONSCINTILLATIONYIELD is "
						      << "set by the user\n";

      fNumScintPhotons =  scintYield * fEnergyDeposit;
    }
    else if(fEMSaturation){
      // The default linear scintillation process
      fNumScintPhotons = fScintYieldFactor * scintYield * fEMSaturation->VisibleEnergyDepositionAtAStep(step);
    }
    else{
      fNumScintPhotons = fScintYieldFactor * scintYield * fEnergyDeposit;
    }
    */


    LOG_DEBUG("ISCalculationSeparate") << "number photons: " << fNumScintPhotons 
				       << " energy: "        << fEnergyDeposit/CLHEP::MeV
				       << " saturation: " 
				       << fEMSaturation->VisibleEnergyDepositionAtAStep(step)
				       << " step length: "   << step->GetStepLength()/CLHEP::cm;


    return;
  }

  //----------------------------------------------------------------------------
  int ISCalculationSeparate::BinomFluct ( int N0, double prob ) {
    CLHEP::RandGauss GaussGen(fEngine);
    CLHEP::RandFlat  UniformGen(fEngine);

    double mean = N0*prob;
    double sigma = sqrt(N0*prob*(1-prob));
    int N1 = 0;
    if ( prob == 0.00 ) return N1;
    if ( prob == 1.00 ) return N0;
    
    if ( N0 < 10 ) {
      for(G4int i = 0; i < N0; i++) {
        if(UniformGen.fire() < prob) N1++;
      }
    }
    else {
      N1 = G4int(floor(GaussGen.fire(mean,sigma)+0.5));
    }
    if ( N1 > N0 ) N1 = N0;
    if ( N1 < 0 ) N1 = 0;
    return N1;
  }


}// namespace
