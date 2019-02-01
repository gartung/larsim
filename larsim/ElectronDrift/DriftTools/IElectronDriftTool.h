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

namespace sim
{
    class SimEnergyDeposit;
    class SimChannel;
    class SimDriftedElectronCluster;
}

namespace detsim
{
    // Declare a structure for keeping track of mapping channel indcices
    // to energy deposits
    typedef struct {
        size_t              channelIndex;
        std::vector<size_t> stepList;
    } ChannelBookKeeping_t;
    
    // Define type: channel -> sim::SimChannel's bookkeeping.
    using ChannelMap_t = std::map<raw::ChannelID_t, ChannelBookKeeping_t>;
    
    // Array of maps of channel data indexed by [cryostat,tpc]
    using ChannelMapByCryoTPC = std::vector< std::vector<ChannelMap_t>>;
    
    class IElectronDriftTool
    {
    public:
        virtual ~IElectronDriftTool() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Search for candidate hits on the input waveform
        virtual void driftElectrons(const size_t,
                                    const sim::SimEnergyDeposit&,
                                    CLHEP::RandGauss&,
                                    std::vector<sim::SimChannel>&,
                                    std::vector<sim::SimDriftedElectronCluster>&,
                                    ChannelMapByCryoTPC&) = 0;         // output candidate hits
    };
}

#endif
