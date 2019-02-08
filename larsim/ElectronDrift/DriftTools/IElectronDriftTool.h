///////////////////////////////////////////////////////////////////////
///
/// \file   IElectronDriftTool.h
///
/// \brief  This handles the actual drift of electrons to the anode
///         plane(s) and their output
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IElectronDriftTool_H
#define IElectronDriftTool_H

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "CLHEP/Random/RandGauss.h"

namespace art
{
    class EDProducer;
    class Event;
}

namespace sim
{
    class SimEnergyDeposit;
    class SimChannel;
    class SimDriftedElectronCluster;
}

namespace detsim
{
    // Define a mapping between a channel and the SimEnergyDeposit objects contributing to it
    using ChannelIdxSimEnergyVec = std::pair<size_t, std::vector<size_t>>;
    using ChannelToSimEnergyMap  = std::unordered_map<raw::ChannelID_t, ChannelIdxSimEnergyVec>;
    
    class IElectronDriftTool
    {
    public:
        virtual ~IElectronDriftTool() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Allow our tools to declare data products they plan to output to event store
        virtual void produces(art::EDProducer&) = 0;
        
        // Set up output data products
        virtual void setupOutputDataProducts() = 0;
        
        // Search for candidate hits on the input waveform
        virtual void driftElectrons(const size_t,
                                    const sim::SimEnergyDeposit&,
                                    CLHEP::RandGauss&,
                                    ChannelToSimEnergyMap&) = 0;         // output candidate hits
        
        // Output the data products
        virtual void put(art::Event&) = 0;
    };
}

#endif
