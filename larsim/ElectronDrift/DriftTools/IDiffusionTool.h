///////////////////////////////////////////////////////////////////////
///
/// \file   IDiffusionTool.h
///
/// \brief  This is intended to handle the longitudinal and transverse diffusion
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IDiffusionTool_H
#define IDiffusionTool_H

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcore/Geometry/Geometry.h"

#include "CLHEP/Random/RandGauss.h"

namespace sim
{
    class SimEnergyDeposit;
}

namespace detsim
{
    // Define a structure to contain the results of the diffusion on transport
    class DiffusionVals
    {
    public:
        DiffusionVals() {}
        DiffusionVals(double longDiff, double transDiff1, double transDiff2, double nElec, double clusE) :
            longitudinalDiffusion(longDiff),
            transverseDiffusion_1(transDiff1),
            transverseDiffusion_2(transDiff2),
            clusterNumElectrons(nElec),
            clusterEnergy(clusE)
        {}
        double longitudinalDiffusion;
        double transverseDiffusion_1;
        double transverseDiffusion_2;
        double clusterNumElectrons;
        double clusterEnergy;
    };
    
    // Define type: channel -> sim::SimChannel's bookkeeping.
    using DiffusionValsVec = std::vector<DiffusionVals>;
    
    class IDiffusionTool
    {
    public:
        virtual ~IDiffusionTool() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Search for candidate hits on the input waveform
        virtual bool getDiffusionVec(const sim::SimEnergyDeposit&, CLHEP::RandGauss&) = 0;
        
        // Recover vector of clusters
        virtual const DiffusionValsVec& getDiffusionValsVec() const = 0;
        
        // Recover the TPC geometry information
        virtual const geo::TPCGeo* getTPCGeo() const = 0;
        
        // recover the longitudinal and transverse sigma for diffusion
        virtual double getLongitudinalDiffusionSig() const = 0;
        virtual double getTransverseDiffusionSig()   const = 0;
        
        // Use the diffusion tool to also return drift parameters
        virtual double getDriftDistance() const = 0;
        virtual double getDriftTime()     const = 0;
        
        // Recover the indices for the drift
        virtual int getDriftCoordinate()       const = 0;
        virtual int getTransverseCoordinate1() const = 0;
        virtual int getTransverseCoordinate2() const = 0;

    };
}

#endif
