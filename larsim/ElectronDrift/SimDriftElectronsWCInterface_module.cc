/**
 * @file SimDriftElectonsWCInterface_module.cxx
 *
 * @brief Transports energy depositions in the LAr TPC to the TPC
 * channels.
 *
 * author:
 * This module was prepared by William Seligman (me), based on code
 * that had been in
 * `LArG4::LArVoxelReadout::DriftIonizationElectrons`. However, though
 * I wrote the original LArVoxelReadout code, I have no idea who added
 * DriftIonizationElectrons. I probably will not be able to answer any
 * questions about how this code works.
 *
 * This module acts on sim::SimEnergyDeposit, the single energy
 * depositions from the detector simulation (LArG4), and simulates the
 * transport of the ensuing ionization electrons to the readout
 * channels:
 *
 * 1. the number of ionisation electrons is read from the current
 *   `larg4::IonizationAndScintillation` instance
 * 2. space charge displacement is optionally applied
 * 3. lifetime correction is applied
 * 4. charge is split in small electron clusters
 * 5. each cluster is subject to longitudinal and transverse diffusion
 * 6. each cluster is assigned to one TPC channel for each wire plane
 * 7. optionally, charge is forced to stay on the planes; otherwise charge
 *    drifting outside the plane is lost
 *
 * For each energy deposition, entries on the appropriate
 * `sim::SimChannel` are added, with the information of the position
 * where the energy deposit happened (in global coordinates,
 * centimeters), the ID of the track in the detector simulation which
 * produced the deposition, and the quantized time of arrival to the
 * channel (in global TDC tick units).  At most one entry is added for
 * each electron cluster, but entries from the same energy deposit can
 * be compacted if falling on the same TDC tick.
 *
 * Options
 * --------
 *
 * A few optional behaviours are supported:
 *
 * * lead off-plane charge to the planes: regulated by
 *   `RecoverOffPlaneDeposit()`, if charge which reaches a wire plane
 *   is actually off it by less than the chosen margin, it's accounted for by
 *   that plane; by default the margin is 0 and all the charge off the plane
 *   is lost (with a warning)
 *
 */

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
//#include "larcore/Geometry/GeometryCore.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

// External libraries
#include "CLHEP/Random/RandGauss.h"
#include "TMath.h"

// C++ includes
#include <vector>
#include <map>
#include <algorithm> // std::find
#include <cmath>

//stuff for WireCell Interface
#include "WireCellIface/SimpleDepo.h"

namespace detsim {

  // Base class for creation of raw signals on wires.
  class SimDriftElectronsWCInterface : public art::EDProducer {

  public:

    explicit SimDriftElectronsWCInterface(fhicl::ParameterSet const& pset);
    virtual ~SimDriftElectronsWCInterface() {};

    // Methods that that are available for a module derived from
    // art::EDProducer.
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    // The label of the module that created the sim::SimEnergyDeposit
    // objects (as of Oct-2017, this is probably "largeant").
    art::InputTag fSimModuleLabel;

    const detinfo::DetectorClocks* fTimeService;
    std::unique_ptr<CLHEP::RandGauss> fRandGauss;

    double fElectronLifetime;
    double fElectronClusterSize;
    int    fMinNumberOfElCluster;
    //double fSampleRate; // unused
    //int    fTriggerOffset; // unused
    double fLongitudinalDiffusion;
    double fTransverseDiffusion;

    double fLifetimeCorr_const;
    double fLDiff_const;
    double fTDiff_const;
    double fRecipDriftVel[3];

    bool fStoreDriftedElectronClusters;

    //double fOffPlaneMargin;

    // In order to create the associations, for each channel we create
    // we have to keep track of its index in the output vector, and the
    // indexes of all the steps that contributed to it.
    typedef struct {
      size_t              channelIndex;
      std::vector<size_t> stepList;
    } ChannelBookKeeping_t;

    // Define type: channel -> sim::SimChannel's bookkeeping.
    typedef std::map<raw::ChannelID_t, ChannelBookKeeping_t> ChannelMap_t;


    // Array of maps of channel data indexed by [cryostat,tpc]
    std::vector< std::vector<ChannelMap_t> > fChannelMaps;
    // The above ensemble may be thought of as a 3D array of
    // ChannelBookKeepings: e.g., SimChannel[cryostat,tpc,channel ID].

    // Save the number of cryostats, and the number of TPCs within
    // each cryostat.
    size_t fNCryostats;
    std::vector<size_t> fNTPCs;

    // Per-cluster information.
    std::vector< double > fXDiff;
    std::vector< double > fYDiff;
    std::vector< double > fZDiff;
    std::vector< double > fnElDiff;
    std::vector< double > fnEnDiff;

    double fDriftClusterPos[3];

    // Utility routine.
    //geo::vect::Vector_t RecoverOffPlaneDeposit (geo::vect::Vector_t const&, geo::PlaneGeo const&) const;

    art::ServiceHandle<geo::Geometry> fGeometry;  ///< Handle to the Geometry service
    ::detinfo::ElecClock              fClock;     ///< TPC electronics clock


    geo::TPCID GetTPCID(std::shared_ptr<WireCell::SimpleDepo> depo);
    std::shared_ptr<WireCell::SimpleDepo> CreateDriftedDepo(std::shared_ptr<WireCell::SimpleDepo> const orig_depo,
							    geo::TPCGeo const& tpcGeo);
    void FillElClusterRandomArrays(std::shared_ptr<WireCell::SimpleDepo> const drifted_idepo,
				   double edep_energy);
    void FillSimChannels(std::shared_ptr<WireCell::SimpleDepo> const drifted_idepo,
			 std::unique_ptr< std::vector<sim::SimChannel> > & channels,
			 geo::TPCGeo const& tpcGeo);
    double fXYZ[3];

  }; // class SimDriftElectronsWCInterface

  //-------------------------------------------------
  SimDriftElectronsWCInterface::SimDriftElectronsWCInterface(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<sim::SimChannel> >();

    // create a default random engine; obtain the random seed from
    // NuRandomService, unless overridden in configuration with key
    // "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, pset, "Seed");
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::reconfigure(fhicl::ParameterSet const& p)
  {
    std::string label= p.get< std::string >("SimulationLabel");
    fSimModuleLabel = art::InputTag(label);

    fStoreDriftedElectronClusters = p.get< bool >("StoreDriftedElectronClusters", false);


    return;
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::beginJob()
  {
    //needed for TPC ticks
    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock = fTimeService->TPCClock();

    // Set up the gaussian generator.
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine();
    fRandGauss = std::unique_ptr<CLHEP::RandGauss>(new CLHEP::RandGauss(engine));

    // Define the physical constants we'll use.
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fElectronLifetime      = detprop->ElectronLifetime(); // Electron lifetime as returned by the DetectorProperties service assumed to be in us;
    for (int i = 0; i<3; ++i) {
      double driftVelocity = detprop->DriftVelocity(detprop->Efield(i),
						    detprop->Temperature())*1.e-3;
      //  Drift velocity as returned by the DetectorProperties service assumed to be in cm/us.
      //  Multiply by 1.e-3 to convert into LArSoft standard velocity units, cm/ns;

      fRecipDriftVel[i] = 1./driftVelocity;
    }

    // To-do: Move the parameters we fetch from "LArG4" to detector
    // properties.
    art::ServiceHandle<sim::LArG4Parameters> paramHandle;
    fElectronClusterSize   = paramHandle->ElectronClusterSize();
    fMinNumberOfElCluster  = paramHandle->MinNumberOfElCluster();
    fLongitudinalDiffusion = paramHandle->LongitudinalDiffusion(); // cm^2/ns units
    fTransverseDiffusion   = paramHandle->TransverseDiffusion(); // cm^2/ns units
    
    LOG_DEBUG("SimDriftElectronsWCInterface")  << " e lifetime (ns): "        << fElectronLifetime
					      << "\n Temperature (K): "     << detprop->Temperature()
					      << "\n Drift velocity (cm/ns): "  << 1./fRecipDriftVel[0]
					      <<" "<<1./fRecipDriftVel[1]<<" "<<1./fRecipDriftVel[2];

    // Opposite of lifetime. Convert from us to standard LArSoft time units, ns;
    fLifetimeCorr_const = -1000. * fElectronLifetime;
    fLDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    fTDiff_const        = std::sqrt(2.*fTransverseDiffusion);

    // For this detector's geometry, save the number of cryostats and
    // the number of TPCs within each cryostat.
    fNCryostats = fGeometry->Ncryostats();
    fNTPCs.resize(fNCryostats);
    for ( size_t n = 0; n < fNCryostats; ++n )
      fNTPCs[n] = fGeometry->NTPC(n);

    return;
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::endJob()
  {
  }


  //-------------------------------------------------
  geo::TPCID SimDriftElectronsWCInterface::GetTPCID(std::shared_ptr<WireCell::SimpleDepo> const depo)
  {
    return fGeometry->PositionToTPCID({depo->pos().x(),depo->pos().y(),depo->pos().z()});
  }

  //-------------------------------------------------
  std::shared_ptr<WireCell::SimpleDepo> 
  SimDriftElectronsWCInterface::CreateDriftedDepo(std::shared_ptr<WireCell::SimpleDepo> const orig_idepo,
						  geo::TPCGeo const& tpcGeo)
  {

    // X drift distance - the drift direction can be either in
    // the positive or negative direction, so use std::abs
    /// \todo think about effects of drift between planes.
    double XDrift = std::abs(orig_idepo->pos().x() - tpcGeo.PlaneLocation(0)[0]);

    // Drift time in ns
    double TDrift = XDrift * fRecipDriftVel[0];
    if (tpcGeo.Nplanes() == 2){// special case for ArgoNeuT (plane 0 is the second wire plane)
      TDrift = ((XDrift - tpcGeo.PlanePitch(0,1)) * fRecipDriftVel[0]
		+ tpcGeo.PlanePitch(0,1) * fRecipDriftVel[1]);
    }
	
    //fill out a drifted idepo object...
    auto drifted_idepo = std::make_shared< WireCell::SimpleDepo >
      (
       orig_idepo->time() + TDrift, //simtime + drift time to first plane  in ns
       WireCell::Point({tpcGeo.PlaneLocation(0)[0],orig_idepo->pos().y(),orig_idepo->pos().z()}), //center position at first plane in cm
       orig_idepo->charge()*TMath::Exp(TDrift/fLifetimeCorr_const), //number of electrons after lifetime correction
       orig_idepo,     //pointer to previous depo
       std::sqrt(TDrift)*fLDiff_const*0.5,   //half longitudinal diffusion sigma in cm
       std::sqrt(TDrift)*fTDiff_const*0.5,   //half transverse diffusion sigma in cm
       orig_idepo->id(),
       orig_idepo->pdg()
      );
    return drifted_idepo;
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::FillElClusterRandomArrays(std::shared_ptr<WireCell::SimpleDepo> const drifted_idepo,
							       double edep_energy)
  {

    //ok, now fill the random arrays for electron clusters from the drifted idepo object
    //this is the part that could for the most part stay the same
    //with the exception that 
    double electronclsize = fElectronClusterSize;
    
    // Number of electron clusters.
    int nClus = (int) std::ceil(drifted_idepo->charge() / electronclsize);
    if (nClus < fMinNumberOfElCluster)
      {
	electronclsize = drifted_idepo->charge() / fMinNumberOfElCluster;
	if (electronclsize < 1.0)
	  electronclsize = 1.0;
	nClus = (int) std::ceil(drifted_idepo->charge() / electronclsize);
      }
    
    // Empty and resize the electron-cluster vectors.
    fXDiff.clear();
    fYDiff.clear();
    fZDiff.clear();
    fnElDiff.clear();
    fnEnDiff.clear();
    fXDiff.resize(nClus);
    fYDiff.resize(nClus);
    fZDiff.resize(nClus);
    fnElDiff.resize(nClus, electronclsize);
    fnEnDiff.resize(nClus);
    
    // fix the number of electrons in the last cluster, that has a smaller size
    fnElDiff.back() = drifted_idepo->charge() - (nClus-1)*electronclsize;
    
    for(size_t xx = 0; xx < fnElDiff.size(); ++xx){
      if(drifted_idepo->charge() > 0)
	fnEnDiff[xx] = edep_energy/drifted_idepo->charge()*fnElDiff[xx]; //todo energy from IDepo
      else
	fnEnDiff[xx] = 0.;
    }
    
    // Smear drift times by x position and drift time
    if (drifted_idepo->extent_long() > 0.0)
      fRandGauss->fireArray( nClus, &fXDiff[0], 0., 2.0*drifted_idepo->extent_long());
    else
      fXDiff.assign(nClus, 0.0);
    
    if (drifted_idepo->extent_tran() > 0.0) {
      // Smear the Y,Z position by the transverse diffusion
      fRandGauss->fireArray( nClus, &fYDiff[0], drifted_idepo->pos().y(), 2.0*drifted_idepo->extent_tran());
      fRandGauss->fireArray( nClus, &fZDiff[0], drifted_idepo->pos().z(), 2.0*drifted_idepo->extent_tran());
    }
    else {
      fYDiff.assign(nClus, drifted_idepo->pos().y());
      fZDiff.assign(nClus, drifted_idepo->pos().z());
    }
    
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::FillSimChannels(std::shared_ptr<WireCell::SimpleDepo> const drifted_idepo,
						     std::unique_ptr< std::vector<sim::SimChannel> > & channels,
						     geo::TPCGeo const& tpcGeo)
  {
    auto tpcid = tpcGeo.ID();

    // make a collection of electrons for each plane
    for(size_t p = 0; p < tpcGeo.Nplanes(); ++p){
      
      double Plane0Pitch = tpcGeo.Plane0Pitch(p);
      
      // "-" sign is because Plane0Pitch output is positive. Andrzej
      fDriftClusterPos[0] = tpcGeo.PlaneLocation(0)[0] - Plane0Pitch;
      
      // Drift nClus electron clusters to the induction plane
      for(unsigned int k = 0; k < fnElDiff.size(); ++k){
	
	// Correct drift time for longitudinal diffusion and plane
	double TDiff = drifted_idepo->time() + fXDiff[k] * fRecipDriftVel[0];
	
	// Take into account different Efields between planes
	// Also take into account special case for ArgoNeuT where Nplanes = 2.
	for (size_t ip = 0; ip<p; ++ip){
	  TDiff += tpcGeo.PlanePitch(ip,ip+1) * fRecipDriftVel[tpcGeo.Nplanes()==3?ip+1:ip+2];
	}
	fDriftClusterPos[1] = fYDiff[k];
	fDriftClusterPos[2] = fZDiff[k];
	
	/// \todo think about effects of drift between planes
	
	// grab the nearest channel to the fDriftClusterPos position
	try{
	  
	  raw::ChannelID_t channel = fGeometry->NearestChannel(fDriftClusterPos, p, tpcid.TPC, tpcid.Cryostat);
	  unsigned int tdc = fClock.Ticks(fTimeService->G4ToElecTime(TDiff));
	  
	  // Find whether we already have this channel in our map.
	  ChannelMap_t& channelDataMap = fChannelMaps[tpcid.Cryostat][tpcid.TPC];
	  auto search = channelDataMap.find(channel);
	  
	  // We will find (or create) the pointer to a
	  // sim::SimChannel.
	  //sim::SimChannel* channelPtr = NULL;
	  size_t channelIndex=0;
	  
	  // Have we created the sim::SimChannel corresponding to
	  // channel ID?
	  if (search == channelDataMap.end())
	    {
	      
	      // We haven't. Initialize the bookkeeping information
	      // for this channel.
	      ChannelBookKeeping_t bookKeeping;
	      
	      // Add a new channel to the end of the list we'll
	      // write out after we've processed this event.
	      bookKeeping.channelIndex = channels->size();
	      channels->emplace_back( channel );
	      channelIndex = bookKeeping.channelIndex;
	      
	      // Save the bookkeeping information for this channel.
	      channelDataMap[channel] = bookKeeping;
	    }
	  else {
	    // We've created this SimChannel for a previous energy
	    // deposit. Get its address.
	    
	    auto& bookKeeping = search->second;
	    channelIndex = bookKeeping.channelIndex;
	    
	  }
	  sim::SimChannel* channelPtr = &(channels->at(channelIndex));
	  
	  // Add the electron clusters and energy to the
	  // sim::SimChannel
	  fXYZ[0] = drifted_idepo->prior()->pos().x();
	  fXYZ[1] = drifted_idepo->prior()->pos().y();
	  fXYZ[2] = drifted_idepo->prior()->pos().z();
	  channelPtr->AddIonizationElectrons(drifted_idepo->id(), //original track ID
					     tdc, //time tick
					     fnElDiff[k], //n electrons
					     fXYZ, //original xyz
					     fnEnDiff[k] //energy);
					     );
	}
	catch(cet::exception &e) {
	  mf::LogWarning("SimDriftElectronsWCInterface") << "unable to drift electrons from point ("
							 << drifted_idepo->prior()->pos().x() << "," 
							 << drifted_idepo->prior()->pos().y() << "," 
							 << drifted_idepo->prior()->pos().z() << ")" 
							 << " with exception " << e;
	} // end try to determine channel
      } // end loop over clusters
    } // end loop over planes
  }

  //-------------------------------------------------
  void SimDriftElectronsWCInterface::produce(art::Event& event)
  {
    // Fetch the SimEnergyDeposit objects for this event.
    typedef art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle_t;
    energyDepositHandle_t energyDepositHandle;
    // If there aren't any energy deposits for this event, don't
    // panic. It's possible someone is doing a study with events
    // outside the TPC, or where there are only non-ionizing
    // particles, or something like that.
    if (!event.getByLabel(fSimModuleLabel, energyDepositHandle))
      return;

    // Define the container for the SimChannel objects that will be
    // transferred to the art::Event after the put statement below.
    std::unique_ptr< std::vector<sim::SimChannel> > 		channels(new std::vector<sim::SimChannel>);

    // Clear the channel maps from the last event. Remember,
    // fChannelMaps is an array[cryo][tpc] of maps.
    size_t cryo = 0;
    fChannelMaps.resize(fNCryostats);
    for (auto& cryoData: fChannelMaps) { // each, a vector of maps
      cryoData.resize(fNTPCs[cryo++]);
      for (auto& channelsMap: cryoData) channelsMap.clear(); // each, a map
    }

    // We're going through the input vector by index, rather than by
    // iterator, because we need the index number to compute the
    // associations near the end of this method.
    auto const& energyDeposits = *energyDepositHandle;
    auto energyDepositsSize = energyDeposits.size();

    // For each energy deposit in this event
    for ( size_t edIndex = 0; edIndex < energyDepositsSize; ++edIndex )
      {

	//starting point: edep
	auto const& energyDeposit = energyDeposits[edIndex];

	//if no electrons to drift, then continue
	if(energyDeposit.NumElectrons()<=0) {
	  LOG_DEBUG("SimDriftElectronsWCInterface")
	    << "step " << edIndex << ". "
	    << "No electrons drifted to readout, " << energyDeposit.Energy() << " MeV lost.";
	  continue;
	}

	//let's turn the larsoft object into a WireCell::SimpleDepo for fun
	auto orig_idepo = std::make_shared< WireCell::SimpleDepo >
	  (
	   energyDeposit.Time(), //simtime (this is in ns...)
	   WireCell::Point({energyDeposit.X(),
			    energyDeposit.Y(),
			    energyDeposit.Z()}
			   ), //center position at first plane in cm
	   energyDeposit.NumElectrons(), //number of electrons after lifetime correction
	   WireCell::IDepo::pointer(nullptr), //would be pointer to previous depo
	   0.0, //no diffusion on initial
	   0.0, //no diffusion on initial
	   energyDeposit.TrackID(),
	   energyDeposit.PdgCode()
	   );


	//need to get the TPC we are in for this depo...
	geo::TPCID tpcid = GetTPCID(orig_idepo);

	//check validity
	if(!tpcid){ //if not valid
	  mf::LogWarning("SimDriftElectronsWCInterface") << "step cannot be found in a cryostat or TPC\n";
	  continue; //move to next edep
	}
	
	//get the actual geo object
	geo::TPCGeo const& tpcGeo = fGeometry->TPC(tpcid);

	//make a check on the drift distance...
	if(tpcGeo.DriftDirection()==geo::kNegX && tpcGeo.PlaneLocation(0)[0]>orig_idepo->pos().x())
	  continue;
	else if(tpcGeo.DriftDirection()==geo::kPosX && tpcGeo.PlaneLocation(0)[0]<orig_idepo->pos().x())
	  continue;

	//and, create a new 'drifted' idepo object, based on the original one and the plane location
	auto drifted_idepo = CreateDriftedDepo(orig_idepo,tpcGeo);

	//Wes 15Aug2018
	//now the hard part ... we need to figure calculate the spread on each wire for each plane
	//analytically it can be done. but drift module now does this approximately via random numbers
	//it'd be best to match this to how WireCell calculates this to avoid confusion due to random
	//but (1) if electron clusters are small enough and (2) we don't demand perfect matching to blips
	//I think this will be ok. Importantly, we ensure the same underlying timing and diffusion
	FillElClusterRandomArrays(drifted_idepo,
				  energyDeposit.Energy()); //note, need the energy and it's not in the drited_idepo
	
	//then we can fill the sim channels
	FillSimChannels(drifted_idepo,channels,tpcGeo);

      } // for each sim::SimEnergyDeposit 

   
    
    // Write the sim::SimChannel collection.
    event.put(std::move(channels));
    
    return;
  }
  
} 

DEFINE_ART_MODULE(detsim::SimDriftElectronsWCInterface)

