////////////////////////////////////////////////////////////////////////
/// \file   DiffusionStandard.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larsim/ElectronDrift/DriftTools/IDiffusionTool.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// LArSoft includes
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimDriftedElectronCluster.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

#include "CLHEP/Random/RandGauss.h"

//stuff from wes
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

// Eigen
#include <Eigen/Dense>

#include "TFile.h"
#include "TTree.h"
#include "TCollection.h"
#include "TKey.h"
#include "TH1F.h"
#include "TProfile.h"

#include <cmath>
#include <fstream>

namespace detsim
{

class DiffusionStandard : IDiffusionTool
{
public:
    explicit DiffusionStandard(const fhicl::ParameterSet& pset);
    
    ~DiffusionStandard();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    // Search for candidate hits on the input waveform
    bool getDiffusionVec(const sim::SimEnergyDeposit&, CLHEP::RandGauss&) override;
    
    // Recover vector of clusters
    const DiffusionValsVec& getDiffusionValsVec() const override {return fDiffusionValsVec;}
    
    // Recover the TPC geometry information
    const geo::TPCGeo* getTPCGeo() const override {return fTPCGeo;}
    
    // recover the longitudinal and transverse sigma for diffusion
    double getLongitudinalDiffusionSig() const override {return fLongitudinalDiffusionSig;}
    double getTransverseDiffusionSig()   const override {return fTransverseDiffusionSig;}
    
    // Use the diffusion tool to also return drift parameters
    double getDriftDistance() const override {return fDriftDistance;}
    double getDriftTime()     const override {return fDriftTime;}
    
    // Recover the indices for the drift
    int getDriftCoordinate()       const override {return fDriftCoordinate;}
    int getTransverseCoordinate1() const override {return fTransverseCoordinate1;}
    int getTransverseCoordinate2() const override {return fTransverseCoordinate2;}

private:

    const detinfo::DetectorClocks*    fTimeService;

    double                            fElectronLifetime;
    double                            fElectronClusterSize;
    int                               fMinNumberOfElCluster;
    double                            fLongitudinalDiffusion;
    double                            fTransverseDiffusion;
    double                            fLifetimeCorr_const;
    double                            fLDiff_const;
    double                            fTDiff_const;
    double                            fRecipDriftVel[3];
    
    // output variables
    const geo::TPCGeo*                fTPCGeo;
    double                            fDriftDistance;
    double                            fDriftTime;
    double                            fLongitudinalDiffusionSig;
    double                            fTransverseDiffusionSig;
    int                               fDriftCoordinate;
    int                               fTransverseCoordinate1;
    int                               fTransverseCoordinate2;

    DiffusionValsVec                  fDiffusionValsVec;
    
    //IS calculationg
    larg4::ISCalcSeparate             fISAlg;

    const geo::GeometryCore*          fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
DiffusionStandard::DiffusionStandard(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
DiffusionStandard::~DiffusionStandard()
{
}
    
void DiffusionStandard::configure(const fhicl::ParameterSet& pset)
{
    // Recover our parameters
//    fResponseFileName             = pset.get<std::string>("ResponseFileName",             "t600_responses_2D_v0.0.root");
//    fStoreDriftedElectronClusters = pset.get< bool      >("StoreDriftedElectronClusters", false);

    // Define the physical constants we'll use.
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    fElectronLifetime = detprop->ElectronLifetime(); // Electron lifetime as returned by the DetectorProperties service assumed to be in us;
    
    for (int i = 0; i<3; ++i)
    {
        //  Drift velocity as returned by the DetectorProperties service assumed to be in cm/us. Multiply by 1.e-3 to convert into LArSoft standard velocity units, cm/ns;
        double driftVelocity = detprop->DriftVelocity(detprop->Efield(i), detprop->Temperature())*1.e-3;
        
        fRecipDriftVel[i] = 1./driftVelocity;
    }
    
    // To-do: Move the parameters we fetch from "LArG4" to detector
    // properties.
    art::ServiceHandle<sim::LArG4Parameters> paramHandle;
    fElectronClusterSize   = paramHandle->ElectronClusterSize();
    fMinNumberOfElCluster  = paramHandle->MinNumberOfElCluster();
    fLongitudinalDiffusion = paramHandle->LongitudinalDiffusion(); // cm^2/ns units
    fTransverseDiffusion   = paramHandle->TransverseDiffusion(); // cm^2/ns units
    
    MF_LOG_DEBUG("SimDriftElectrons")  << " e lifetime (ns): " << fElectronLifetime
                                       << "\n Temperature (K): "     << detprop->Temperature()
                                       << "\n Drift velocity (cm/ns): "  << 1./fRecipDriftVel[0]
                                       <<" "<<1./fRecipDriftVel[1]<<" "<<1./fRecipDriftVel[2];

    // Opposite of lifetime. Convert from us to standard LArSoft time units, ns;
    fLifetimeCorr_const = -1000. * fElectronLifetime;
    fLDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    fTDiff_const        = std::sqrt(2.*fTransverseDiffusion);

    fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
                      detprop,
                      &(*paramHandle),
                      lar::providerFrom<spacecharge::SpaceChargeService>());

    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();

    return;
}
    
bool DiffusionStandard::getDiffusionVec(const sim::SimEnergyDeposit& energyDeposit, CLHEP::RandGauss& randGauss)
{
    // Clear the last diffusion vector
    fDiffusionValsVec.clear();
    
    // "xyz" is the position of the energy deposit in world
    // coordinates. Note that the units of distance in
    // sim::SimEnergyDeposit are supposed to be cm.
    auto const mp = energyDeposit.MidPoint();
    double const xyz[3] = { mp.X(), mp.Y(), mp.Z() };
    
    // From the position in world coordinates, determine the
    // cryostat and tpc. If somehow the step is outside a tpc
    // (e.g., cosmic rays in rock) just move on to the next one.
    unsigned int cryostat = 0;
    try {
        fGeometry->PositionToCryostat(xyz, cryostat);
    }
    catch(cet::exception &e){
        mf::LogWarning("SimDriftElectrons") << "step "// << energyDeposit << "\n"
        << "cannot be found in a cryostat\n"
        << e;
        return false;
    }
    
    unsigned int tpc = 0;
    try {
        fGeometry->PositionToTPC(xyz, tpc, cryostat);
    }
    catch(cet::exception &e){
        mf::LogWarning("SimDriftElectrons") << "step "// << energyDeposit << "\n"
        << "cannot be found in a TPC\n"
        << e;
        return false;
    }
    
    fTPCGeo = fGeometry->TPCPtr(geo::TPCID(cryostat, tpc));
    
    // The drift direction can be either in the positive
    // or negative direction in any coordinate x, y or z.
    // Charge drift in ...
    // +x: tpcGeo->DetectDriftDirection()==1
    // -x: tpcGeo->DetectDriftDirection()==-1
    // +y: tpcGeo->DetectDriftDirection()==2
    // -y tpcGeo->DetectDriftDirection()==-2
    // +z: tpcGeo->DetectDriftDirection()==3
    // -z: tpcGeo->DetectDriftDirection()==-3
    //Define charge drift direction: driftcoordinate (x, y or z) and driftsign (positive or negative). Also define coordinates perpendicular to drift direction.
    fDriftCoordinate = std::abs(fTPCGeo->DetectDriftDirection())-1;  //x:0, y:1, z:2
    
    if(fDriftCoordinate == 0)
    {
        fTransverseCoordinate1 = 1;
        fTransverseCoordinate2 = 2;
    }
    else if(fDriftCoordinate == 1)
    {
        fTransverseCoordinate1 = 0;
        fTransverseCoordinate2 = 2;
    }
    else if(fDriftCoordinate == 2)
    {
        fTransverseCoordinate1 = 0;
        fTransverseCoordinate2 = 1;
    }
    
    if(fTransverseCoordinate1 == fTransverseCoordinate2) return false; //this is the case when driftcoordinate != 0, 1 or 2
    
    int driftsign = 0; //1: +x, +y or +z, -1: -x, -y or -z
    
    if(fTPCGeo->DetectDriftDirection() > 0) driftsign = 1;
    else driftsign = -1;
    
    //Check for charge deposits behind charge readout planes
    if(driftsign ==  1 && fTPCGeo->PlaneLocation(0)[fDriftCoordinate] < xyz[fDriftCoordinate] ) return false;
    if(driftsign == -1 && fTPCGeo->PlaneLocation(0)[fDriftCoordinate] > xyz[fDriftCoordinate] ) return false;
    
    /// \todo think about effects of drift between planes.
    // Center of plane is also returned in cm units
    fDriftDistance = std::abs(xyz[fDriftCoordinate] - fTPCGeo->PlaneLocation(0)[fDriftCoordinate]);
    
    // Space-charge effect (SCE): Get SCE {x,y,z} offsets for
    // particular location in TPC
    geo::Vector_t posOffsets{0.0,0.0,0.0};
    double posOffsetxyz[3] = {0.0,0.0,0.0}; //need this array for the driftcoordinate and transversecoordinates
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    
    if (SCE->EnableSimSpatialSCE() == true)
    {
        posOffsets = SCE->GetPosOffsets(mp);
        posOffsetxyz[0] = posOffsets.X();
        posOffsetxyz[1] = posOffsets.Y();
        posOffsetxyz[2] = posOffsets.Z();
    }
    
    double avegagetransversePos1 = 0.;
    double avegagetransversePos2 = 0.;
    
    fDriftDistance         += -1.*posOffsetxyz[fDriftCoordinate];
    avegagetransversePos1  = xyz[fTransverseCoordinate1] + posOffsetxyz[fTransverseCoordinate1];
    avegagetransversePos2  = xyz[fTransverseCoordinate2] + posOffsetxyz[fTransverseCoordinate2];
    
    // Space charge distortion could push the energy deposit beyond the wire
    // plane (see issue #15131). Given that we don't have any subtlety in the
    // simulation of this region, bringing the deposit exactly on the plane
    // should be enough for the time being.
//    if (DriftDistance < 0.) DriftDistance = 0.;
    
    // Drift time in ns
    fDriftTime = 0.;
    
    if (fDriftDistance > 0.) fDriftTime = fDriftDistance * fRecipDriftVel[0];
    
    if (fTPCGeo->Nplanes() == 2 && fDriftCoordinate == 0){// special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is the second wire plane
        fDriftTime = ((fDriftDistance - fTPCGeo->PlanePitch(0,1)) * fRecipDriftVel[0] + fTPCGeo->PlanePitch(0,1) * fRecipDriftVel[1]);
    }
    
    fISAlg.Reset();
    fISAlg.CalculateIonizationAndScintillation(energyDeposit);
    //std::cout << "Got " << fISAlg.NumberIonizationElectrons() << "." << std::endl;
    
    const double lifetimecorrection = TMath::Exp(fDriftTime / fLifetimeCorr_const);
    const int    nIonizedElectrons  = fISAlg.NumberIonizationElectrons();
    const double energy             = energyDeposit.Energy();
    
    // if we have no electrons (too small energy or too large recombination)
    // we are done already here
    if (nIonizedElectrons <= 0) {
        MF_LOG_DEBUG("getDiffusionVec")
        << "step "// << energyDeposit << "\n"
        << "No electrons drifted to readout, " << energy << " MeV lost.";
        return false;
    }
    
    // includes the effect of lifetime: lifetimecorrection = exp[-tdrift/tau]
    const double nElectrons = nIonizedElectrons * lifetimecorrection;
    //std::cout << "After lifetime, " << nElectrons << " electrons." << std::endl;
    
    // Longitudinal & transverse diffusion sigma (cm)
    double SqrtT          = std::sqrt(fDriftTime);
    double electronclsize = fElectronClusterSize;
    
    fLongitudinalDiffusionSig = SqrtT * fLDiff_const;
    fTransverseDiffusionSig   = SqrtT * fTDiff_const;
    
    // Number of electron clusters.
    int nClus = (int) std::ceil(nElectrons / electronclsize);
    if (nClus < fMinNumberOfElCluster)
    {
        electronclsize = nElectrons / fMinNumberOfElCluster;
        if (electronclsize < 1.0)
        {
            electronclsize = 1.0;
        }
        nClus = (int) std::ceil(nElectrons / electronclsize);
    }
    
    // Empty and resize the electron-cluster vectors.
    std::vector< double > longDiff(nClus);
    std::vector< double > transDiff1(nClus);
    std::vector< double > transDiff2(nClus);
    std::vector< double > nElDiff(nClus, electronclsize);
    std::vector< double > nEnDiff(nClus);

    // fix the number of electrons in the last cluster, that has a smaller size
    nElDiff.back() = nElectrons - (nClus-1)*electronclsize;
    
    for(size_t xx = 0; xx < nElDiff.size(); ++xx)
    {
        if(nElectrons > 0) nEnDiff[xx] = energy/nElectrons*nElDiff[xx];
        else               nEnDiff[xx] = 0.;
    }

    // Smear drift times by longitudinal diffusion
    if (fLongitudinalDiffusionSig > 0.0)
        randGauss.fireArray( nClus, &longDiff[0], 0., fLongitudinalDiffusionSig);
    else
        longDiff.assign(nClus, 0.0);
    
    if (fTransverseDiffusionSig > 0.0)
    {
        // Smear the coordinates in plane perpendicular to drift direction by the transverse diffusion
        randGauss.fireArray( nClus, &transDiff1[0], avegagetransversePos1, fTransverseDiffusionSig);
        randGauss.fireArray( nClus, &transDiff2[0], avegagetransversePos2, fTransverseDiffusionSig);
    }
    else
    {
        transDiff1.assign(nClus, avegagetransversePos1);
        transDiff2.assign(nClus, avegagetransversePos2);
    }
    
    for(size_t idx = 0; idx < nElDiff.size(); idx++)
        fDiffusionValsVec.emplace_back(longDiff[idx], transDiff1[idx], transDiff2[idx], nElDiff[idx], nEnDiff[idx]);

    return true;
}

DEFINE_ART_CLASS_TOOL(DiffusionStandard)
}
