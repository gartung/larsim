////////////////////////////////////////////////////////////////////////
/// \file   ElectronDriftStandard.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larsim/ElectronDrift/DriftTools/IElectronDriftTool.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larsim/ElectronDrift/DriftTools/IDiffusionTool.h"

#include "CLHEP/Random/RandGauss.h"

//stuff from wes
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "TH1F.h"
#include "TProfile.h"

#include <cmath>
#include <fstream>

namespace detsim
{

class ElectronDriftStandard : IElectronDriftTool
{
public:
    explicit ElectronDriftStandard(const fhicl::ParameterSet& pset);
    
    ~ElectronDriftStandard();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    // Search for candidate hits on the input waveform
    void driftElectrons(const size_t,
                        const sim::SimEnergyDeposit&,
                        CLHEP::RandGauss&,
                        std::vector<sim::SimChannel>&,
                        std::vector<sim::SimDriftedElectronCluster>&,
                        ChannelMapByCryoTPC&,
                        ChannelToSimEnergyMap&) override;         // output candidate hits

private:
    
    bool                              fStoreDriftedElectronClusters;

    const detinfo::DetectorClocks*    fTimeService;

    double                            fRecipDriftVel[3];
    
    // The tool to handle the diffusion
    std::unique_ptr<detsim::IDiffusionTool> fDiffusionTool; ///< Tool for handling the diffusion during the drift

    const geo::GeometryCore*          fGeometry = lar::providerFrom<geo::Geometry>();
    ::detinfo::ElecClock              fClock;     ///< TPC electronics clock
};
    
//----------------------------------------------------------------------
// Constructor.
ElectronDriftStandard::ElectronDriftStandard(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ElectronDriftStandard::~ElectronDriftStandard()
{
}
    
void ElectronDriftStandard::configure(const fhicl::ParameterSet& pset)
{
    // Recover our parameters
    fStoreDriftedElectronClusters = pset.get< bool >("StoreDriftedElectronClusters", false);

    // Define the physical constants we'll use.
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    for (int i = 0; i<3; ++i)
    {
        //  Drift velocity as returned by the DetectorProperties service assumed to be in cm/us. Multiply by 1.e-3 to convert into LArSoft standard velocity units, cm/ns;
        double driftVelocity = detprop->DriftVelocity(detprop->Efield(i), detprop->Temperature())*1.e-3;
        
        fRecipDriftVel[i] = 1./driftVelocity;
    }

    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock       = fTimeService->TPCClock();
    
    fDiffusionTool = art::make_tool<detsim::IDiffusionTool>(pset.get<fhicl::ParameterSet>("DiffusionTool"));

    return;
}
    
void ElectronDriftStandard::driftElectrons(const size_t                                 edIndex,
                                           const sim::SimEnergyDeposit&                 energyDeposit,
                                           CLHEP::RandGauss&                            randGauss,
                                           std::vector<sim::SimChannel>&                simChannelVec,
                                           std::vector<sim::SimDriftedElectronCluster>& simDriftedElectronClustersVec,
                                           ChannelMapByCryoTPC&                         channelMap,
                                           ChannelToSimEnergyMap&                       channelToSimEnergyMap)
{
    // First call the tool to do the drifting
    if (fDiffusionTool->getDiffusionVec(energyDeposit, randGauss))
    {
        // "xyz" is the position of the energy deposit in world
        // coordinates. Note that the units of distance in
        // sim::SimEnergyDeposit are supposed to be cm.
        auto const mp = energyDeposit.MidPoint();
        double const xyz[3] = { mp.X(), mp.Y(), mp.Z() };
        
        // Recover drift distance and drift time
        double driftTime = fDiffusionTool->getDriftTime();
        
        // Recover the coordinate indices
        int driftCoordinate       = fDiffusionTool->getDriftCoordinate();
        int transverseCoordinate1 = fDiffusionTool->getTransverseCoordinate1();
        int transverseCoordinate2 = fDiffusionTool->getTransverseCoordinate2();
        
        // Get a pointer to the TPC geometry
        const geo::TPCGeo* tpcGeo = fDiffusionTool->getTPCGeo();
        
        double fDriftClusterPos[3];

        // make a collection of electrons for each plane
        for(size_t p = 0; p < tpcGeo->Nplanes(); ++p)
        {
            fDriftClusterPos[driftCoordinate] = tpcGeo->PlaneLocation(p)[driftCoordinate];
            
            // Now loop over each of the clusters we have created
            for (const auto &cluster : fDiffusionTool->getDiffusionValsVec())
            {
                // Correct drift time for longitudinal diffusion and plane
                double TDiff = driftTime + cluster.longitudinalDiffusion * fRecipDriftVel[0];
                
                // Take into account different Efields between planes
                // Also take into account special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is the second wire plane
                for (size_t ip = 0; ip<p; ++ip)
                {
                    TDiff += (tpcGeo->PlaneLocation(ip+1)[driftCoordinate] - tpcGeo->PlaneLocation(ip)[driftCoordinate]) * fRecipDriftVel[(tpcGeo->Nplanes() == 2 && driftCoordinate == 0)?ip+2:ip+1];
                }
                
                fDriftClusterPos[transverseCoordinate1] = cluster.transverseDiffusion_1;
                fDriftClusterPos[transverseCoordinate2] = cluster.transverseDiffusion_2;
                
                /// \todo think about effects of drift between planes
                
                // grab the nearest channel to the fDriftClusterPos position
                try{
                    /*
                     if (fOffPlaneMargin != 0) {
                     // get the effective position where to consider the charge landed;
                     //
                     // Some optimisations are possible; in particular, this method
                     // could be extended to inform us if the point was too far.
                     // Currently, if that is the case the code will proceed, find the
                     // point is off plane, emit a warning and skip the deposition.
                     //
                     auto const landingPos
                     = RecoverOffPlaneDeposit({ fDriftClusterPos[0], fDriftClusterPos[1], fDriftClusterPos[2] }, plane);
                     fDriftClusterPos[0] = landingPos.X();
                     fDriftClusterPos[1] = landingPos.Y();
                     fDriftClusterPos[2] = landingPos.Z();
                     
                     } // if charge lands off plane
                     */
                    const geo::TPCID& tpcID = tpcGeo->ID();
                    
                    raw::ChannelID_t channel = fGeometry->NearestChannel(fDriftClusterPos, p, tpcID.TPC, tpcID.Cryostat);
                    
                    //          std::cout << "fDriftClusterPos[0]: " << fDriftClusterPos[0] << "\t fDriftClusterPos[1]: " << fDriftClusterPos[1] << "\t fDriftClusterPos[2]: " << fDriftClusterPos[2] << std::endl;
                    //          std::cout << "channel: " << channel << std::endl;
                    
                    //std::cout << "\tgot channel " << channel << " for cluster " << k << std::endl;
                    
                    /// \todo check on what happens if we allow the tdc value to be
                    /// \todo beyond the end of the expected number of ticks
                    // Add potential decay/capture/etc delay effect, simTime.
                    auto const   simTime = energyDeposit.Time();
                    unsigned int tdc     = fClock.Ticks(fTimeService->G4ToElecTime(TDiff + simTime));
                    
                    // Recover the SimChannel and handle the bookeeping
                    size_t channelIndex(0);
                    
                    ChannelToSimEnergyMap::iterator chanIdxItr = channelToSimEnergyMap.find(channel);
                    
                    // Have we already created a SimChannel for this channel?
                    if (chanIdxItr != channelToSimEnergyMap.end())
                    {
                        // Recover the info...
                        ChannelIdxSimEnergyVec& channelIdxSimEnergyVec = chanIdxItr->second;
                        
                        channelIndex = channelIdxSimEnergyVec.first;
                        
                        std::vector<size_t>& simEnergyVec = channelIdxSimEnergyVec.second;
                        
                        // Has this step contributed to this channel already?
                        std::vector<size_t>::iterator chanItr = std::find(simEnergyVec.begin(),simEnergyVec.end(),edIndex);
                        
                        // If not then keep track
                        if (chanItr == simEnergyVec.end()) simEnergyVec.push_back(edIndex);
                    }
                    // Otherwise need to create a new entry
                    else
                    {
                        // Add the new SimChannel (and get the index to it)
                        channelIndex = simChannelVec.size();
                        
                        simChannelVec.emplace_back(channel);
                        
                        ChannelIdxSimEnergyVec& channelIdxToSimEnergyVec = channelToSimEnergyMap[channel];
                        
                        channelIdxToSimEnergyVec.first = channelIndex;
                        channelIdxToSimEnergyVec.second.clear();
                    }
                    
                    sim::SimChannel* channelPtr = &(simChannelVec.at(channelIndex));
                    
                    //if(!channelPtr) std::cout << "\tUmm...ptr is NULL?" << std::endl;
                    //else std::cout << "\tChannel is " << channelPtr->Channel() << std::endl;
                    // Add the electron clusters and energy to the
                    // sim::SimChannel
                    channelPtr->AddIonizationElectrons(energyDeposit.TrackID(),
                                                       tdc,
                                                       cluster.clusterNumElectrons,
                                                       xyz,
                                                       cluster.clusterEnergy);
                    
                    if(fStoreDriftedElectronClusters)
                        simDriftedElectronClustersVec.push_back(sim::SimDriftedElectronCluster(cluster.clusterNumElectrons,
                                                                                               TDiff + simTime,        // timing
                                                                                               {mp.X(),mp.Y(),mp.Z()}, // mean position of the deposited energy
                                                                                               {fDriftClusterPos[0],fDriftClusterPos[1],fDriftClusterPos[2]}, // final position of the drifted cluster
                                                                                               {fDiffusionTool->getLongitudinalDiffusionSig(),
                                                                                                fDiffusionTool->getTransverseDiffusionSig(),
                                                                                                fDiffusionTool->getTransverseDiffusionSig()}, // Longitudinal (X) and transverse (Y,Z) diffusion
                                                                                               cluster.clusterEnergy, //deposited energy that originated this cluster
                                                                                               energyDeposit.TrackID()) );
                }
                catch(cet::exception &e) {
                    mf::LogWarning("SimDriftElectrons") << "unable to drift electrons from point ("
                    << xyz[0] << "," << xyz[1] << "," << xyz[2]
                    << ") with exception " << e;
                } // end try to determine channel
            } // end loop over clusters
        } // end loop over planes
    }
    
    return;
}

DEFINE_ART_CLASS_TOOL(ElectronDriftStandard)
}
