////////////////////////////////////////////////////////////////////////
// Class:       PhotonLibraryPropagationS2
// Plugin Type: producer (art v2_05_00)
// File:        PhotonLibraryPropagationS2_module.cc
//
// Generated at Tue Mar 21 07:45:42 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "Geant4/G4Poisson.hh"

#include <memory>
#include <iostream>

#include "larsim/PhotonPropagation/PhotonVisibilityServiceS2.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "larsim/LArG4/OpDetSensitiveDetector.h"
#include "larsim/LArG4/OpDetReadoutGeometry.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/PhotonVoxels.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
#include "larsim/IonizationScintillation/ISCalculationSeparate.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
// ROOT Includes
#include "TGeoManager.h"

namespace phot {
  class PhotonLibraryPropagationS2;
}


class phot::PhotonLibraryPropagationS2 : public art::EDProducer {
public:
  explicit PhotonLibraryPropagationS2(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonLibraryPropagationS2(PhotonLibraryPropagationS2 const &) = delete;
  PhotonLibraryPropagationS2(PhotonLibraryPropagationS2 &&) = delete;
  PhotonLibraryPropagationS2 & operator = (PhotonLibraryPropagationS2 const &) = delete;
  PhotonLibraryPropagationS2 & operator = (PhotonLibraryPropagationS2 &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void Print(std::map<int, std::map<int, int>>* StepPhotonTable);

private:

  bool fUseLitePhotons;

  std::string fDriftEModuleLabel;
  double fGain;//Number of photons created per drifted electron.

};


phot::PhotonLibraryPropagationS2::PhotonLibraryPropagationS2(fhicl::ParameterSet const & p)
{

  art::ServiceHandle<sim::LArG4Parameters> lgp;

  fUseLitePhotons=lgp->UseLitePhotons();
  if(fUseLitePhotons)
  {
     produces< std::vector<sim::OpDetBacktrackerRecord> >();
     produces< std::vector<sim::SimPhotonsLite> >();
  }
  else     produces< std::vector<sim::SimPhotons> >();  

  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "photon",    p, "SeedPhoton");
  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", p, "SeedScintTime");  
  this->reconfigure(p);
}


void phot::PhotonLibraryPropagationS2::produce(art::Event & e)
{


  art::ServiceHandle<PhotonVisibilityServiceS2> pvs;

  art::ServiceHandle<geo::Geometry> geom;
  const detinfo::LArProperties* larp = lar::providerFrom<detinfo::LArPropertiesService>();
  
  art::ServiceHandle<art::RandomNumberGenerator> rng;  
  CLHEP::HepRandomEngine &engine_photon = rng->getEngine("photon");
  CLHEP::RandPoissonQ randpoisphot(engine_photon);
  CLHEP::HepRandomEngine &engine_scinttime = rng->getEngine("scinttime");
  CLHEP::RandFlat randflatscinttime(engine_scinttime);


    larg4::OpDetPhotonTable::Instance()->ClearTable(geom->NOpDets());

  // Get the pointer to the fast scintillation table
  larg4::OpDetPhotonTable * fst = larg4::OpDetPhotonTable::Instance();
  larg4::OpDetPhotonTable* litefst = larg4::OpDetPhotonTable::Instance();
  
  const size_t NOpChannels = pvs->NOpChannels();
  double nphot,nphot_fast,nphot_slow;

  std::unique_ptr< std::vector<sim::SimPhotons>  >               PhotonCol                  (new std::vector<sim::SimPhotons>);
  std::unique_ptr< std::vector<sim::SimPhotonsLite>  >           LitePhotonCol              (new std::vector<sim::SimPhotonsLite>);
  std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > >  cOpDetBacktrackerRecordCol (new std::vector<sim::OpDetBacktrackerRecord>);

  art::ValidHandle< std::vector<sim::SimDriftedElectronCluster> > ElectronClusters_handle = e.getValidHandle<std::vector<sim::SimDriftedElectronCluster>>(fDriftEModuleLabel);

  //std::cout << "... GETTING ELECTRON CLUSTER HANDLE"<<std::endl;

  int counter=0;

  if(fUseLitePhotons) //std::cout << "... USING LITE PHOTONS" << std::endl;
  else //std::cout << "... NOT USING LITE PHOTONS!!!! ---ERROR---" << std::endl;
  for (sim::SimDriftedElectronCluster const& ElectronCluster: *ElectronClusters_handle)
  {	
	//std::cout << "reading a cluster " << counter <<std::endl; counter ++;

	double const xyz[3] = { ElectronCluster.FinalPositionX(), ElectronCluster.FinalPositionY(), ElectronCluster.FinalPositionZ() };
	float const* Visibilities = pvs->GetAllVisibilities(xyz);
	if(!Visibilities)
	continue;
	const std::vector<float>* PropParameters = nullptr;
	if(pvs->IncludeParPropTime()) PropParameters = pvs->GetTimingPar(xyz);
	//std::cout <<"\tPosition> " <<ElectronCluster.FinalPositionX()<<" " << ElectronCluster.FinalPositionY()<<" " <<ElectronCluster.FinalPositionZ() <<std::endl;	
	//std::cout <<"\tTiming> " <<ElectronCluster.getTime()<<std::endl;
	//std::cout <<"\tEnergy> " <<ElectronCluster.getEnergy()<<std::endl;
	//std::cout <<"\telectrons> " <<ElectronCluster.getNumberOfElectrons()<<std::endl;
	//std::cout <<"\tTrack> " <<ElectronCluster.getTrackID()<<std::endl;	
      
	nphot =ElectronCluster.getNumberOfElectrons()*fGain; // # photons generated in the Gas Ar phase

	if(!Visibilities)
	{
	}
	else
	{

		//std::cout <<"\t\tVisibilities loaded. Let's calculate the number of detected photons. OpChannels: " <<NOpChannels<<std::endl;	

		std::map<int, int> DetectedNum;
		std::map<int, TF1> PropTimeFunction;

		for(size_t OpDet=0; OpDet!=NOpChannels; OpDet++)
		{
			G4int DetThisPMT = G4int(randpoisphot.fire(Visibilities[OpDet] * nphot));

			//std::cout <<"\t\tVisibilities loaded. Let's calculate the number of detected photons. OpChannels: " <<Visibilities[OpDet] << " "<< nphot << " " << Visibilities[OpDet] * nphot << std::endl;	
			if(DetThisPMT>0) 
        		{
			      DetectedNum[OpDet]=DetThisPMT;
			     //std::cout <<"\t\t\t#photons per pmt"<<OpDet<<": "<<DetThisPMT<<std::endl;
			      //   //   it->second<<" " << Num << " " << DetThisPMT;  

			      //det_photon_ctr += DetThisPMT; // CASE-DEBUG DO NOT REMOVE THIS COMMENT
        		}

		    	if(pvs->IncludeParPropTime() && PropParameters)
			{
		    	  double range=0;
			  for(size_t i=0; i<pvs->ParPropTimeNpar();i++) range+=300*TMath::Abs(PropParameters[OpDet][i]);
			  TF1 AuxFunction(Form("timingfunc%i",(int)((size_t)OpDet)),Form("%s",pvs->ParPropTimeFormula().c_str()),PropParameters[OpDet][0],PropParameters[OpDet][0]+range); 
			  for(size_t i=1; i<pvs->ParPropTimeNpar();i++) AuxFunction.SetParameter(i-1, PropParameters[OpDet][i]);
			  PropTimeFunction[OpDet]=AuxFunction;
	  		}
                }
	    // Now we run through each PMT figuring out num of detected photons
	
		if(fUseLitePhotons)
		{
			//std::cout << "\t\tLet's create our StepPhotonTable Map." << std::endl;	

			std::map<int, std::map<int, int>> StepPhotonTable;
			// And then add these to the total collection for the event  

			//std::cout << "before"<<std::endl;
			Print(&StepPhotonTable);

			for(std::map<int,int>::const_iterator itdetphot = DetectedNum.begin();
				itdetphot!=DetectedNum.end(); ++itdetphot)
			{

				//std::cout << "... iterating! " <<  std::endl;
				std::map<int, int>  StepPhotons;
				for (G4int i = 0; i < itdetphot->second; ++i)
				{
					G4double deltaTime = ElectronCluster.getTime();

					if(pvs->IncludeParPropTime())
					{
						deltaTime += PropTimeFunction[itdetphot->first].GetRandom(); 
					}

					G4double aSecondaryTime = deltaTime;
					float Time = aSecondaryTime;
					int ticks = static_cast<int>(Time);
					StepPhotons[ticks]++;
			 	}
				 StepPhotonTable[itdetphot->first] = StepPhotons;
				 //Iterate over Step Photon Table to add photons to OpDetBacktrackerRecords.

				 sim::OpDetBacktrackerRecord tmpOpDetBTRecord(itdetphot->first);
				 //int thisG4TrackID = (aStep.GetTrack())->GetTrackID();
				 int thisG4TrackID = ElectronCluster.getTrackID();
				 double xO = ( ElectronCluster.FinalPositionX() / CLHEP::cm );
				 double yO = ( ElectronCluster.FinalPositionY() / CLHEP::cm );
				 double zO = ( ElectronCluster.FinalPositionZ() / CLHEP::cm );
				 double const xyzPos[3] = {xO,yO,zO};
				//         double energy  = ( aStep.GetTotalEnergyDeposit() / CLHEP::MeV );
				 double energy  = ElectronCluster.getEnergy()/CLHEP::MeV;
			 	//Loop over StepPhotons to get number of photons detected at each time for this channel and G4Step.
			 	for(std::map<int,int>::iterator stepPhotonsIt = StepPhotons.begin(); stepPhotonsIt != StepPhotons.end(); ++stepPhotonsIt)
			 	{
			   		int photonTime = stepPhotonsIt->first;
			   		int numPhotons = stepPhotonsIt->second;
			   		tmpOpDetBTRecord.AddScintillationPhotons(thisG4TrackID, photonTime, numPhotons, xyzPos, energy);
				}
			 	//Add OpDetBackTrackerRecord. (larg4::OpDetPhotonTABLE->instance().addOpDetBacktrackerRecord(sim::OpDetBacktrackerRecord BTRrecord)
				 litefst->AddOpDetBacktrackerRecord(tmpOpDetBTRecord);
			}
			//std::cout << "after"<<std::endl;
			Print(&StepPhotonTable);

			litefst->AddPhoton(&StepPhotonTable);


			////std::cout << "total map....."<<std::endl;
			//Print(&(litefst->GetLitePhotons()));
		}
		else //std::cout << "SIMPHOTONS NOT SUPPORTED!!! "<< std::endl;
	}

//	OpDetSensitiveDetector *theOpDetDet = dynamic_cast<OpDetSensitiveDetector*>(sdManager->FindSensitiveDetector("OpDetSensitiveDetector"));
//	if (OpDetSensitiveDetector)
//	{
  }// endif electroncluster loop

  if(!fUseLitePhotons)
  {

	//std::cout << "SIMPHOTONS NOT SUPPORTED!!! "<< std::endl;
      /*
	LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
	std::vector<sim::SimPhotons>& ThePhotons = larg4::OpDetPhotonTable::Instance()->GetPhotons();
	PhotonCol->reserve(ThePhotons.size());
	for(auto& it : ThePhotons)
	PhotonCol->push_back(std::move(it));*/
  }
  else
  {



	LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";

//	std::map<int, std::map<int, int> > ThePhotons = larg4::OpDetPhotonTable::Instance()->GetLitePhotons();
	std::map<int, std::map<int, int> > ThePhotons = litefst->GetLitePhotons();

	//std::cout << "... Storing OpDet Hit Collection in Event. Size of map>" << ThePhotons.size() <<  std::endl;
	Print(&ThePhotons);


	if(ThePhotons.size() > 0)
	{
		LitePhotonCol->reserve(ThePhotons.size());

		for(auto const& it : ThePhotons)
		{

	                
			sim::SimPhotonsLite ph;
			ph.OpChannel = it.first;
			ph.DetectedPhotons = it.second;
			LitePhotonCol->push_back(ph);
			//std::cout << "... \t" << it.first << std::endl;
		}
	}
	//*cOpDetBacktrackerRecordCol = larg4::OpDetPhotonTable::Instance()->YieldOpDetBacktrackerRecords();
	*cOpDetBacktrackerRecordCol = litefst->YieldOpDetBacktrackerRecords();
  }

  if(!fUseLitePhotons) e.put(std::move(PhotonCol));
  else
  {
	//std::cout << "... Storing litephotons in file " << LitePhotonCol->size() <<  std::endl;
	e.put(std::move(LitePhotonCol));
	e.put(std::move(cOpDetBacktrackerRecordCol));
   }
  
}

  void phot::PhotonLibraryPropagationS2::Print(std::map<int, std::map<int, int>>* StepPhotonTable)
  {
    for(auto it = StepPhotonTable->begin(); it!=StepPhotonTable->end(); it++)
    {
      for(auto in_it = it->second.begin(); in_it!=it->second.end(); in_it++)
      {
        //std::cout << in_it->second << " ";
      }
	//std::cout << std::endl;
    }
  }


void phot::PhotonLibraryPropagationS2::reconfigure(fhicl::ParameterSet const & p)
{
  fGain = p.get<double>("Gain",500);
  fDriftEModuleLabel= p.get< std::string >("DriftEModuleLabel");
}

void phot::PhotonLibraryPropagationS2::beginJob()
{

}

DEFINE_ART_MODULE(phot::PhotonLibraryPropagationS2)
