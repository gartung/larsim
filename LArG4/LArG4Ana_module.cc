////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_LARG4ANA_H
#define LARG4_LARG4ANA_H 

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

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
#include "LaRLightEnergyAction.h"

#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ios.hh"
// LArSoft Includes
#include "MCCheater/BackTracker.h"
#include "Simulation/ParticleList.h"
#include "Simulation/sim.h"
#include "Geometry/Geometry.h"

// ROOT includes
#include "TH1.h"
#include "TMath.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include <TTree.h>

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

///Geant4 interface 
namespace larg4 {  
 
  class LArG4Ana : public art::EDAnalyzer{
  public:
 
    /// Standard constructor and destructor for an FMWK module.
    explicit LArG4Ana(fhicl::ParameterSet const& pset);
    virtual ~LArG4Ana();

    void analyze (const art::Event& evt); 
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& pset);

  private:

    std::string fG4ModuleLabel;     ///< module label for the Geant
    std::string fTruthModuleLabel;  ///< module label for the Geant


    TH1D *fPDGCodes;
    TH1D *fEdep;
    TH1D *fPi0Momentum;
    TH1D *fnEnergy;
    TH1D *fnDist;
    TH1D *fnumChannels; ///< The number of channels recieving charge per event
    TProfile *fnumIDEs; ///< Number of drift electrons per channel.
    TH1D *fEventCharge; ///< Charge collected per event
    TH1D *fEventEnergy; ///< Energy collected per event
    TProfile *fChannelCharge; ///< Charge per channel.
    TProfile *fChannelEnergy; ///< Energy per channel.
    
    //    Int_t stringDim = 35;

    TTree *fTree;
//added
//TTree *fTree2;
    Int_t fTEvt;
Float_t fTEnedep;
Float_t fTDIST;
    Int_t fTSub;
  //  Int_t fTRefl;;
    Int_t fTRun;
    Int_t fTPdg;
    Int_t fTID;
    Int_t fTNds;
    Int_t fTNdsOriginal;
    Int_t fTNds4;
    Int_t *fTDID;
    Int_t *fTDPdg;
    Float_t *fTDWt;
    Char_t fTProcess[35];
    Char_t fTVolume[35];
    Char_t fTTVolume[35]; // Termination Volume
    Char_t fTMaterial[35];
    Char_t fTDProcess[200][35];
    Int_t fTParentID; 
    Int_t fTStatus;
    Float_t fTWeight;
    Float_t* fT4Origin;
    Float_t* fT4DOrigin;
    Float_t* fT4Termination; // Termination Coordinates
    Float_t* fT4Momentum;
    Float_t* fT4DMomentum;
bool fill_tree;
  };

} // namespace larg4

namespace larg4 {

  //-----------------------------------------------------------------------
  // Constructor
  LArG4Ana::LArG4Ana(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  LArG4Ana::~LArG4Ana() 
  {
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;
    
    fPDGCodes    = tfs->make<TH1D>("pdgcodes", ";PDG Code;",               5000, -2500, 2500);
fEdep    = tfs->make<TH1D>("energy deposit", "edeposit",300, 0, 1.5);
//fTRefl    = tfs->make<TH1D>("reflected photons", "reflphot",1000, 0, 10000);
    fPi0Momentum = tfs->make<TH1D>("pi0mom",   ";#pi^{0} Momentum (GeV);", 1000, 0.,    1000.);

    fTree = tfs->make<TTree>("MCTTree","MCTTree");
    fnEnergy = tfs->make<TH1D>("nEnergy", ";n,#Lambda^{0},K^{0} Momentum (GeV);", 100, 0., 10.);
    fnDist = tfs->make<TH1D>("nDistance", ";n,#Lambda^{0},K^{0} Distance (m);", 200, -30000.0, +30000.);


    // Some histograms relating to drift electrons, active detector
    // channels and charge/energy on channels
    fnumChannels = tfs->make<TH1D>("fnumChannels", 
				   "Active channels;Active channels;# events",
				   256, 0, geo->Nchannels()); 
    fnumIDEs = tfs->make<TProfile>("fnumIDEs", 
				   "Drift Electrons per channel;Channel;Drift electrons",
				   geo->Nchannels()+1, 0, geo->Nchannels(),
				   0, 1e4); 
    fEventCharge = tfs->make<TH1D>("fEventCharge", 
				   "Charge in event;Total charge per event;# events",
				   100, 0, 2.5e8); 
    fEventEnergy = tfs->make<TH1D>("fEventEnergy", 
				   "Energy in event;Total energy per event;# events",
				   100, 0, 1e4); 
    fChannelCharge = tfs->make<TProfile>("fChannelCharge", 
					 "Charge on channel;Channel;Total charge per channel",
					 geo->Nchannels()+1,0,geo->Nchannels(),
					 0,1e5); 
    fChannelEnergy = tfs->make<TProfile>("fChannelEnergy", 
					 "Energy on channel;Channel;Total energy per channel",
					 geo->Nchannels()+1,0,geo->Nchannels(),
					 0, 1e3); 


    fT4Origin = new Float_t[4];
    fT4DOrigin = new Float_t[fTNds*4];
    fT4Termination = new Float_t[4];
    fT4Momentum = new Float_t[4];
    fT4DMomentum = new Float_t[fTNds*4];
    fTDID = new Int_t[fTNds];
    fTDPdg = new Int_t[fTNds];
    fTDWt  = new Float_t[fTNds];
    fTNds4 = fTNds*4; //  TTree/Branch requirement to store this.

    fTree->Branch("MCEvt", &fTEvt, "MCEvt/I");

    fTree->Branch("MCEnedep", &fTEnedep, "MCEnedep/F");
    //fTree->Branch("MCRefl", &fTRefl, "MCRefl/I");
fTree->Branch("MCDIST", &fTDIST, "MCDIST/F");
    fTree->Branch("MCSub", &fTSub, "MCSub/I");
    fTree->Branch("MCRun", &fTRun, "MCRun/I");
    fTree->Branch("MCWt", &fTWeight, "MCWt/F");
    fTree->Branch("MCPdg", &fTPdg, "MCPdg/I");
    fTree->Branch("MCID", &fTID, "MCID/I");
    fTree->Branch("MCParentID", &fTParentID, "MCParentID/I");
    fTree->Branch("MCNumDs", &fTNds, "MCNumDs/I");
    fTree->Branch("MCNumDs4", &fTNds4, "MCNumDs4/I");
    fTree->Branch("MCDID", fTDID, "MCDID[MCNumDs]/I");
    fTree->Branch("MCDPdg", fTDPdg, "MCDPdg[MCNumDs]/I");
    fTree->Branch("MCDWt", fTDWt, "MCDWt[MCNumDs]/I");
    fTree->Branch("MCProcess", fTProcess, "MCProcess/C");
    fTree->Branch("MCVolume", fTVolume, "MCVolume/C");
    fTree->Branch("MCTVolume", fTTVolume, "MCTVolume/C");
    fTree->Branch("MCMaterial", fTMaterial, "MCMaterial/C");
    fTree->Branch("MCDProcess", fTDProcess, "MCDProcess[MCNumDs]/C");
    fTree->Branch("MCStatus", &fTStatus, "MCStatus/I");
    fTree->Branch("MCOrigin", fT4Origin, "MCOrigin[4]/F");
    fTree->Branch("MCDOrigin", fT4DOrigin, "MCDOrigin[MCNumDs4]/F");
    fTree->Branch("MCTermination", fT4Termination, "MCTermination[4]/F");
    fTree->Branch("MCMomentum", fT4Momentum, "MCMomentum[4]/F");
    fTree->Branch("MCDMomentum", fT4DMomentum, "MCDMomentum[MCNumDs4]/F");
  

//fTree2->Branch("EvtID", &fTEvt, "EvtID/I");
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::reconfigure(fhicl::ParameterSet const& p)
  {
    fG4ModuleLabel    = p.get< std::string >("GeantModuleLabel");
    fTNdsOriginal     = p.get< int         >("Ndaughters"      );
    fTruthModuleLabel = p.get< std::string >("TruthModuleLabel"); 
    fTNds = fTNdsOriginal;
    
    return;
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::analyze(const art::Event& evt) 
  {

    //get the list of particles from this event
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    art::ServiceHandle<geo::Geometry> geom;

//-----------------------------------------------------
  // Define a "handle" to point to BackTracker.
    // BackTracker is a tool connecting the tracks with the particles.

 ////////////////////////////////////////////////////////////////////
    //Determination of the electron energy at the beginning of the TPC
    //JUST THE ELECTRON!!! IF IT MAKES BREMS,THOSE ARE NOT CONSIDERED!!!
    int tsize=plist.Particle(0)->Trajectory().size();
int list_size=plist.size();
    int tt=0;
int tts=0;
int tu=0;
   // double en_ent=0.0;
    //double en_ext=0.0;
    double en_dep=0.0;
double dist2=0.;
std::vector<double> vec1;
std::vector<double> pt1;
std::vector<double> pt2;

std::cout<<"analyzing trajectory  in list of size "<<list_size<<std::endl;

bool first_entry=false;
bool entry_test=false;
vec1.clear();
pt1.clear();
pt2.clear();
fTEnedep=0.;
//fTRefl=0;
for(tu=0;tu<list_size;tu++){
tts=0;
en_dep=0.;
entry_test=false;
    for(tt=0;tt<int(plist.Particle(tu)->Trajectory().size());tt++)
    {
std::cout<<tu<<" "<<plist.Particle(tu)->TrackId()<<" BEFORE ENTERING BOX @@@@@@@@@@@@ POSITION X "<<plist.Particle(tu)->Trajectory().X(tt)<<" @@@@@@@@@@@@ POSITION Y "<<plist.Particle(tu)->Trajectory().Y(tt)<<" @@@@@@@@@@@@ POSITION Z "<<plist.Particle(tu)->Trajectory().Z(tt)<<" ENERGY "<<plist.Particle(tu)->Trajectory().E(tt)<<std::endl;
     if ((TMath::Abs(double(plist.Particle(tu)->Trajectory().X(tt))) <= 7.5) && (TMath::Abs(double(plist.Particle(tu)->Trajectory().Y(tt))) <= 7.5) && (TMath::Abs(double(plist.Particle(tu)->Trajectory().Z(tt))) <= 12.5)){

std::cout<<tu<<"&&&&&&&& IN CHAMBER @@@@@@@@@@@@ POSITION X "<<plist.Particle(tu)->Trajectory().X(tt)<<" @@@@@@@@@@@@ POSITION Y "<<plist.Particle(tu)->Trajectory().Y(tt)<<" @@@@@@@@@@@@ POSITION Z "<<plist.Particle(tu)->Trajectory().Z(tt)<<" ENERGY "<<plist.Particle(tu)->Trajectory().E(tt)<<std::endl;
   	  if(entry_test==false){
			 if(tt>0) tts=tt-1;
			 else tts=tt;
				}
	  
	  entry_test=true;
						}

if (entry_test && ((TMath::Abs(double(plist.Particle(tu)->Trajectory().X(tt))) > 7.5) || (TMath::Abs(double(plist.Particle(tu)->Trajectory().Y(tt))) > 7.5) || (TMath::Abs(double(plist.Particle(tu)->Trajectory().Z(tt))) > 12.5))){
		en_dep=double(plist.Particle(tu)->Trajectory().E(tts))-double(plist.Particle(tu)->Trajectory().E(tt));

std::cout<<tu<<" LEAVING CHAMBER @@@@@@@@@@@@ POSITION X "<<plist.Particle(tu)->Trajectory().X(tt)<<" @@@@@@@@@@@@ POSITION Y "<<plist.Particle(tu)->Trajectory().Y(tt)<<" @@@@@@@@@@@@ POSITION Z "<<plist.Particle(tu)->Trajectory().Z(tt)<<" ENERGY "<<plist.Particle(tu)->Trajectory().E(tt)<<" starting point "<<tts<<" end "<<tt<<std::endl;
 		break;
		}
 

	    }//loop over each trajectory
std::cout<<tu<<" TRAJECTORY finished size -> "<<plist.Particle(tu)->Trajectory().size()<<" energy deposited "<<float(en_dep)*1000.<<std::endl;
	//if (entry_test==true)fTEnedep+=float(en_dep)*1000.;
	//else fTEnedep+=0.;
	//en_ent=0.;
	//en_ext=0.;
	}//loop over trajectories

//dist2=sqrt(pow((plist.Particle(0)->Trajectory().X(tt)-plist.Particle(0)->Trajectory().X(tts)),2.0)+pow((plist.Particle(0)->Trajectory().Y(tt)-plist.Particle(0)->Trajectory().Y(tts)),2.0)+pow((plist.Particle(0)->Trajectory().Z(tt)-plist.Particle(0)->Trajectory().Z(tts)),2.0));



fTEnedep=energy_deposit_step;
//fTRefl=phot_refl_ev;
fTDIST=float(dist2);
//fEdep->Fill(fTEnedep);//MeV

//evt.id().event(); 


//-----------------------------------------------------

    // loop over all sim::SimChannels in the event and make sure there are no
    // sim::IDEs with trackID values that are not in the sim::ParticleList
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);

    double totalCharge=0.0;
    double totalEnergy=0.0;
    int temp_id=0;
    fnumChannels->Fill(sccol.size());
    for(size_t sc = 0; sc < sccol.size(); ++sc){
      double numIDEs=0.0;
      double scCharge=0.0;
      double scEnergy=0.0;
      const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = sccol[sc]->TDCIDEMap();
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	const std::vector<sim::IDE> idevec = (*mapitr).second;
	numIDEs += idevec.size();
	for(size_t iv = 0; iv < idevec.size(); ++iv){
	  if(plist.find( idevec[iv].trackID ) == plist.end()
	     && idevec[iv].trackID != sim::NoParticleId) 
	  mf::LogWarning("LArG4Ana") << idevec[iv].trackID << " is not in particle list"; 
	  totalCharge +=idevec[iv].numElectrons;
	  scCharge += idevec[iv].numElectrons;
	  totalEnergy +=idevec[iv].energy;
	  scEnergy += idevec[iv].energy;
	}
      }
      fnumIDEs->Fill(sc,numIDEs);
      fChannelCharge->Fill(sc,scCharge);
      fChannelEnergy->Fill(sc,scEnergy);
    }
    fEventCharge->Fill(totalCharge);
    fEventEnergy->Fill(totalEnergy);

    // get the particles from the back tracker
    const sim::ParticleList& Particles = bt->ParticleList();
    std::vector<const simb::MCParticle*> pvec;
    pvec.reserve(Particles.size());
    for (const auto& PartPair: Particles) {
      pvec.push_back(PartPair.second);
      fPDGCodes->Fill(PartPair.second->PdgCode());
    }

    // now look for pi0's that decay to 2 gammas
    int pi0loc = -1;
    int numpi0gamma = 0;
    for(unsigned int i = 0; i < pvec.size(); ++i){
      if(pvec[i]->PdgCode() == 111) pi0loc = i;
      if(pvec[i]->Mother() == pi0loc+1 && 
	 pi0loc > 0 && 
	 pvec[i]->PdgCode() == 22){
	mf::LogInfo("LArG4Ana") << pvec[i]->E() << " gamma energy ";
	++numpi0gamma;
      }
      
      // n,Lambda,K0s,K0L,K0
      if (pvec[i]->PdgCode() == 2112 || 
	  pvec[i]->PdgCode() == 3122 || 
	  pvec[i]->PdgCode() == 130  || 
	  pvec[i]->PdgCode() == 310  || 
	  pvec[i]->PdgCode() == 311 ) {
	fnEnergy->Fill(pvec[i]->E(),pvec[i]->Weight());
	fnDist->Fill(pvec[i]->Vx(),pvec[i]->Weight());
      }

      fTPdg = pvec[i]->PdgCode();
      fTID = pvec[i]->TrackId();
      // 0 out strings, else there may be cruft in here from prev evt.
      for (unsigned int s = 0; s < 35; ++s){
	  *(fTProcess+s)  = 0; 
	  *(fTProcess+s)  = 0;
	  *(fTMaterial+s) = 0; 
	  *(fTMaterial+s) = 0;
	  *(fTVolume+s)   = 0;
	  *(fTVolume+s)   = 0;
	  *(fTTVolume+s)  = 0; 
	  *(fTTVolume+s)  = 0;
      }
      
      for(unsigned int s = 0; s < pvec[i]->Process().length(); ++s) *(fTProcess+s) = pvec[i]->Process()[s];

      TVector3 dum = pvec[i]->Position().Vect();

      for (unsigned int s = 0; s < geom->MaterialName(pvec[i]->Position().Vect()).length(); ++s) 
	*(fTMaterial+s) = geom->MaterialName(pvec[i]->Position().Vect())[s];

      for (unsigned int s = 0; s < geom->VolumeName(pvec[i]->Position().Vect()).length(); ++s) 
	*(fTVolume+s) = geom->VolumeName(pvec[i]->Position().Vect())[s];

      for (unsigned int s = 0; s < geom->VolumeName(pvec[i]->EndPosition().Vect()).length(); ++s) 
	*(fTTVolume+s) = geom->VolumeName(pvec[i]->EndPosition().Vect())[s];

      fTEvt = evt.id().event(); 
    if(fTEvt!=0){
	if(fTEvt!=temp_id){
			fill_tree=true;
			temp_id=fTEvt;
			}
	else fill_tree=false;

		}
    else{
	fill_tree=true;


		}
      fTSub = evt.subRun();
      fTRun = evt.run();
      fTParentID = pvec[i]->Mother();
      fTStatus = pvec[i]->StatusCode();
      int daughter = 9999;
      fTNds = TMath::Min(pvec[i]->NumberDaughters(),fTNdsOriginal);
      for( int d = 0; d < fTNds; d++ ){
	daughter = pvec[i]->Daughter(d);
	fTDID[d] = daughter; 
	// zero it out.
	for (unsigned int s = 0; s < 35; ++s) *(fTDProcess[d]+s) = 0; 
	
	for(unsigned int jj = i; jj < pvec.size(); ++jj){ // Don't look below i.
       
	  if (fTDID[d] == pvec[jj]->TrackId()){
	    fTDPdg[d] = pvec[jj]->PdgCode(); // get the pointer,  
	    fTDWt[d]  = pvec[jj]->Weight();
	    
	    for (unsigned int s = 0; s < pvec[jj]->Process().length(); ++s) 
	      *(fTDProcess[d]+s) = pvec[jj]->Process()[s];
		
	    for (unsigned int kk = 0; kk < 4; ++kk){
	      fT4DOrigin[d*4+kk] = pvec[jj]->Position()[kk];
	      fT4DMomentum[d*4+kk] = pvec[jj]->Momentum()[kk];
	    }
	    break;
	  }	      
	}
      }//end loop over d
      
      for (unsigned int ii = 0; ii < 4; ++ii){
	fT4Termination[ii] = 1e9;
	fT4Origin[ii] = pvec[i]->Position()[ii];
	if (ii!=3) fT4Termination[ii] = pvec[i]->EndPosition()[ii];
	if (ii==4) fT4Termination[ii] = pvec[i]->Momentum()[ii]; // yes, odd
	fT4Momentum[ii] = pvec[i]->Momentum()[ii];
      }

      fTWeight = pvec[i]->Weight();
      if(fill_tree==true) fTree->Fill();
      
    } // end loop on particles in list 
    if(numpi0gamma == 2 && pi0loc > 0){
      mf::LogInfo("LArG4Ana") << pvec[pi0loc]->E();
      fPi0Momentum->Fill(pvec[pi0loc]->E());
    }

    return;
  }



} // namespace larg4

namespace larg4 {

  DEFINE_ART_MODULE(LArG4Ana)

} // namespace LArG4

#endif // LARG4_LARG4_H

