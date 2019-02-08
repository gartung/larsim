/**
 * @file SimDriftElectrons_module.cxx
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
 * Update:
 * Christoph Alt, September 2018 (christoph.alt@cern.ch)
 * Break hardcoded charge drift in x to support charge drift in y and z.
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
#include "art/Utilities/make_tool.h"
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

//stuff from wes
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

// Drift tool
#include "larsim/ElectronDrift/DriftTools/IElectronDriftTool.h"

namespace detsim {
    
// Base class for creation of raw signals on wires.
class SimDriftElectrons : public art::EDProducer {
    
public:
    
    explicit SimDriftElectrons(fhicl::ParameterSet const& pset);
    virtual ~SimDriftElectrons() {};
    
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
    
    std::unique_ptr<CLHEP::RandGauss> fRandGauss;
    
    bool fStoreDriftedElectronClusters;
    
    // In order to create the associations, for each channel we create
    // we have to keep track of its index in the output vector, and the
    // indexes of all the steps that contributed to it.
    // Array of maps of channel data indexed by [cryostat,tpc]
    ChannelToSimEnergyMap fChannelToSimEnergyMap;
    
    // The above ensemble may be thought of as a 3D array of
    // ChannelBookKeepings: e.g., SimChannel[cryostat,tpc,channel ID].
    
    // Save the number of cryostats, and the number of TPCs within
    // each cryostat.
    size_t fNCryostats;
    std::vector<size_t> fNTPCs;
    
    // Utility routine.
    //geo::vect::Vector_t RecoverOffPlaneDeposit (geo::vect::Vector_t const&, geo::PlaneGeo const&) const;
    
    art::ServiceHandle<geo::Geometry> fGeometry;  ///< Handle to the Geometry service
    
    // The tool to drift individual energy deposits
    std::unique_ptr<detsim::IElectronDriftTool> fDriftTool; ///< Tool for generating noise
    
}; // class SimDriftElectrons

//-------------------------------------------------
SimDriftElectrons::SimDriftElectrons(fhicl::ParameterSet const& pset)
: art::EDProducer{pset}
{
    this->reconfigure(pset);
        
    //produces< art::Assns<sim::SimChannel, sim::SimEnergyDeposit> >();
    
    // create a default random engine; obtain the random seed from
    // NuRandomService, unless overridden in configuration with key
    // "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
    ->createEngine(*this, pset, "Seed");
    
    /**
     * @brief Sets the margin for recovery of charge drifted off-plane.
     * @param margin the extent of the margin on each frame coordinate [cm]
     *
     * This method sets the margin for the recovery of off-plane ionization
     * charge. See `RecoverOffPlaneDeposit()` for a description of that feature.
     *
     */
        //fOffPlaneMargin = pset.get< double >("ChargeRecoveryMargin",0.0);
        // Protection against a silly value.
        //fOffPlaneMargin = std::max(fOffPlaneMargin,0.0);
}

//-------------------------------------------------
void SimDriftElectrons::reconfigure(fhicl::ParameterSet const& p)
{
    std::string label= p.get< std::string >("SimulationLabel");
    fSimModuleLabel = art::InputTag(label);
    
    fStoreDriftedElectronClusters = p.get< bool >("StoreDriftedElectronClusters", false);
    
    fDriftTool = art::make_tool<detsim::IElectronDriftTool>(p.get<fhicl::ParameterSet>("ElectronDriftTool"));
    
    fDriftTool->produces(*this);

    return;
}

//-------------------------------------------------
void SimDriftElectrons::beginJob()
{
    // Set up the gaussian generator.
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(art::ScheduleID::first(),
                                                    moduleDescription().moduleLabel(),
                                                    {});
    fRandGauss = std::unique_ptr<CLHEP::RandGauss>(new CLHEP::RandGauss(engine));
    
    // For this detector's geometry, save the number of cryostats and
    // the number of TPCs within each cryostat.
    fNCryostats = fGeometry->Ncryostats();
    fNTPCs.resize(fNCryostats);
    for ( size_t n = 0; n < fNCryostats; ++n )
        fNTPCs[n] = fGeometry->NTPC(n);
    
    return;
}
    
//-------------------------------------------------
void SimDriftElectrons::endJob()
{
}

//-------------------------------------------------
void SimDriftElectrons::produce(art::Event& event)
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
    
    // Set up the output data objects in the tool
    fDriftTool->setupOutputDataProducts();
    
    // Clear the channel maps from the last event.
    fChannelToSimEnergyMap.clear();
    
    // We're going through the input vector by index, rather than by
    // iterator, because we need the index number to compute the
    // associations near the end of this method.
    auto const& energyDeposits = *energyDepositHandle;
    auto energyDepositsSize = energyDeposits.size();
    
    // For each energy deposit in this event
    for ( size_t edIndex = 0; edIndex < energyDepositsSize; ++edIndex )
    {
        auto const& energyDeposit = energyDeposits[edIndex];
        
        fDriftTool->driftElectrons(edIndex, energyDeposit, *fRandGauss, fChannelToSimEnergyMap);
    } // for each sim::SimEnergyDeposit
    /*
     // Now that we've processed the information for all the
     // sim::SimEnergyDeposit objects into sim::SimChannel objects,
     // create the associations between them.
     
     // Define the container for the associations between the channels
     // and the energy deposits (steps). Note it's possible for an
     // energy deposit to be associated with more than one channel (if
     // its electrons drift to multiple wires), and a channel will
     // almost certainly have multiple energy deposits.
     std::unique_ptr< art::Assns<sim::SimEnergyDeposit, sim::SimChannel> >
     step2channel (new art::Assns<sim::SimEnergyDeposit, sim::SimChannel>);
     
     // For every element in the 3-D fChannelMaps array...
     for ( auto i = fChannelMaps.begin(); i != fChannelMaps.end(); ++i ) {
     for ( auto j = i->begin(); j != i->end(); ++j ) {
     for ( auto k = j->begin(); k != j->end(); ++k ) {
     const ChannelBookKeeping_t& bookKeeping = (*k).second;
     const size_t channelIndex = bookKeeping.channelIndex;
     
     // Create a one-to-one association between each channel and
     // each step that created it.
     for ( size_t m = 0; m < bookKeeping.stepList.size(); ++m)
     {
     const size_t edIndex = bookKeeping.stepList[m];
     // Props to me for figuring out the following two
     // statements. You have to look deeply in the
     // documentation for art::Ptr and util::Associations to
     // put this together.
     art::Ptr<sim::SimEnergyDeposit> energyDepositPtr( energyDepositHandle, edIndex );
     util::CreateAssn(*this, event, *channels, energyDepositPtr, *step2channel, channelIndex);
     }
     }
     }
     }
     */
    // Write the sim::SimChannel collection.
    fDriftTool->put(event);
    
    // ... and its associations.
    //event.put(std::move(step2channel));
    
    return;
}
    
    //----------------------------------------------------------------------------
/*
 geo::vect::Vector_t SimDriftElectrons::RecoverOffPlaneDeposit
 (geo::vect::Vector_t const& pos, geo::PlaneGeo const& plane) const
 {
 //
 // translate the landing position on the two frame coordinates
 // ("width" and "depth")
 //
 auto const landingPos = plane.PointWidthDepthProjection(pos);
 
 //
 // compute the distance of the landing position on the two frame
 // coordinates ("width" and "depth");
 // keep the point within 10 micrometers (0.001 cm) from the border
 //
 auto const offPlane = plane.DeltaFromActivePlane(landingPos, 0.001);
 
 //
 // if both the distances are below the margin, move the point to
 // the border
 //
 
 // nothing to recover: landing is inside
 if ((offPlane.X() == 0.0) && (offPlane.Y() == 0.0)) return pos;
 
 // won't recover: too far
 if ((std::abs(offPlane.X()) > fOffPlaneMargin)
 || (std::abs(offPlane.Y()) > fOffPlaneMargin))
 return pos;
 
 // we didn't fully decompose because it might be unnecessary;
 // now we need the full thing
 auto const distance = plane.DistanceFromPlane(pos);
 
 return plane.ComposePoint(distance, landingPos + offPlane);
 
 } // SimDriftElectrons::RecoverOffPlaneDeposit()
 */
} // namespace detsim

DEFINE_ART_MODULE(detsim::SimDriftElectrons)
