////////////////////////////////////////////////////////////////////////
/// \file   ElectronDrift2D.cc
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
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

#include "CLHEP/Random/RandGauss.h"

//stuff from wes
#include "larsim/IonizationScintillation/ISCalcSeparate.h"

#include "TH1F.h"
#include "TProfile.h"

#include <cmath>
#include <fstream>

namespace detsim
{

class ElectronDrift2D : IElectronDriftTool
{
public:
    explicit ElectronDrift2D(const fhicl::ParameterSet& pset);
    
    ~ElectronDrift2D();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    // Search for candidate hits on the input waveform
    void driftElectrons(const size_t,
                        const sim::SimEnergyDeposit&,
                        CLHEP::RandGauss&,
                        std::vector<sim::SimChannel>&,
                        std::vector<sim::SimDriftedElectronCluster>&,
                        ChannelMapByCryoTPC&) override;         // output candidate hits

private:
    
    bool                              fStoreDriftedElectronClusters;

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
    
    //IS calculationg
    larg4::ISCalcSeparate             fISAlg;

    const geo::GeometryCore*          fGeometry = lar::providerFrom<geo::Geometry>();
    ::detinfo::ElecClock              fClock;     ///< TPC electronics clock
};
    
//----------------------------------------------------------------------
// Constructor.
ElectronDrift2D::ElectronDrift2D(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ElectronDrift2D::~ElectronDrift2D()
{
}
    
void ElectronDrift2D::configure(const fhicl::ParameterSet& pset)
{
    // Recover our parameters
    fStoreDriftedElectronClusters = pset.get< bool >("StoreDriftedElectronClusters", false);

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
    fClock       = fTimeService->TPCClock();
    
    // If asked, define the global histograms
//    if (fOutputHistograms)
//    {
//        // Access ART's TFileService, which will handle creating and writing
//        // histograms and n-tuples for us.
//        art::ServiceHandle<art::TFileService> tfs;
//
//        fHistDirectory = tfs.get();
//
//        // Make a directory for these histograms
//        art::TFileDirectory dir = fHistDirectory->mkdir(Form("HitPlane_%1zu",fPlane));
//
//        fDStopStartHist        = dir.make<TH1F>(Form("DStopStart_%1zu",   fPlane), ";Delta Stop/Start;",    100,   0., 100.);
//        fDMaxTickMinTickHist   = dir.make<TH1F>(Form("DMaxTMinT_%1zu",    fPlane), ";Delta Max/Min Tick;",  100,   0., 100.);
//        fDMaxDerivMinDerivHist = dir.make<TH1F>(Form("DMaxDMinD_%1zu",    fPlane), ";Delta Max/Min Deriv;", 200,   0., 100.);
//        fMaxErosionHist        = dir.make<TH1F>(Form("MaxErosion_%1zu",   fPlane), ";Max Erosion;",         200, -50., 150.);
//        fMaxDilationHist       = dir.make<TH1F>(Form("MaxDilation_%1zu",  fPlane), ";Max Dilation;",        200, -50., 150.);
//        fMaxDilEroRatHist      = dir.make<TH1F>(Form("MaxDilEroRat_%1zu", fPlane), ";Max Dil/Ero;",         200,  -1.,   1.);
//    }

    return;
}
    
void ElectronDrift2D::driftElectrons(const size_t                                 edIndex,
                                           const sim::SimEnergyDeposit&                 energyDeposit,
                                           CLHEP::RandGauss&                            randGauss,
                                           std::vector<sim::SimChannel>&                simChannelVec,
                                           std::vector<sim::SimDriftedElectronCluster>& simDriftedElectronClustersVec,
                                           ChannelMapByCryoTPC&                         channelMap)
{
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
        return;
    }
    
    unsigned int tpc = 0;
    try {
        fGeometry->PositionToTPC(xyz, tpc, cryostat);
    }
    catch(cet::exception &e){
        mf::LogWarning("SimDriftElectrons") << "step "// << energyDeposit << "\n"
        << "cannot be found in a TPC\n"
        << e;
        return;
    }
    
    const geo::TPCGeo& tpcGeo = fGeometry->TPC(tpc, cryostat);
    
    // The drift direction can be either in the positive
    // or negative direction in any coordinate x, y or z.
    // Charge drift in ...
    // +x: tpcGeo.DetectDriftDirection()==1
    // -x: tpcGeo.DetectDriftDirection()==-1
    // +y: tpcGeo.DetectDriftDirection()==2
    // -y tpcGeo.DetectDriftDirection()==-2
    // +z: tpcGeo.DetectDriftDirection()==3
    // -z: tpcGeo.DetectDriftDirection()==-3
    
    
    //Define charge drift direction: driftcoordinate (x, y or z) and driftsign (positive or negative). Also define coordinates perpendicular to drift direction.
    int driftcoordinate = std::abs(tpcGeo.DetectDriftDirection())-1;  //x:0, y:1, z:2
    
    int transversecoordinate1 = 0;
    int transversecoordinate2 = 0;
    if(driftcoordinate == 0)
    {
        transversecoordinate1 = 1;
        transversecoordinate2 = 2;
    }
    else if(driftcoordinate == 1)
    {
        transversecoordinate1 = 0;
        transversecoordinate2 = 2;
    }
    else if(driftcoordinate == 2)
    {
        transversecoordinate1 = 0;
        transversecoordinate2 = 1;
    }
    
    if(transversecoordinate1 == transversecoordinate2) return; //this is the case when driftcoordinate != 0, 1 or 2
    
    int driftsign = 0; //1: +x, +y or +z, -1: -x, -y or -z
    if(tpcGeo.DetectDriftDirection() > 0) driftsign = 1;
    else driftsign = -1;
    
    //Check for charge deposits behind charge readout planes
    if(driftsign == 1 && tpcGeo.PlaneLocation(0)[driftcoordinate] < xyz[driftcoordinate] )  return;
    if(driftsign == -1 && tpcGeo.PlaneLocation(0)[driftcoordinate] > xyz[driftcoordinate] ) return;
    
    /// \todo think about effects of drift between planes.
    // Center of plane is also returned in cm units
    double DriftDistance = std::abs(xyz[driftcoordinate] - tpcGeo.PlaneLocation(0)[driftcoordinate]);
    
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
    
    DriftDistance += -1.*posOffsetxyz[driftcoordinate];
    avegagetransversePos1 = xyz[transversecoordinate1] + posOffsetxyz[transversecoordinate1];
    avegagetransversePos2 = xyz[transversecoordinate2] + posOffsetxyz[transversecoordinate2];
    
    
    // Space charge distortion could push the energy deposit beyond the wire
    // plane (see issue #15131). Given that we don't have any subtlety in the
    // simulation of this region, bringing the deposit exactly on the plane
    // should be enough for the time being.
    if (DriftDistance < 0.) DriftDistance = 0.;
    
    // Drift time in ns
    double TDrift = DriftDistance * fRecipDriftVel[0];
    
    if (tpcGeo.Nplanes() == 2 && driftcoordinate == 0){// special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is the second wire plane
        TDrift = ((DriftDistance - tpcGeo.PlanePitch(0,1)) * fRecipDriftVel[0]
                  + tpcGeo.PlanePitch(0,1) * fRecipDriftVel[1]);
    }
    
    fISAlg.Reset();
    fISAlg.CalculateIonizationAndScintillation(energyDeposit);
    //std::cout << "Got " << fISAlg.NumberIonizationElectrons() << "." << std::endl;
    
    const double lifetimecorrection = TMath::Exp(TDrift / fLifetimeCorr_const);
    const int    nIonizedElectrons  = fISAlg.NumberIonizationElectrons();
    const double energy             = energyDeposit.Energy();
    
    // if we have no electrons (too small energy or too large recombination)
    // we are done already here
    if (nIonizedElectrons <= 0) {
        MF_LOG_DEBUG("SimDriftElectrons")
        << "step "// << energyDeposit << "\n"
        << "No electrons drifted to readout, " << energy << " MeV lost.";
        return;
    }
    
    // includes the effect of lifetime: lifetimecorrection = exp[-tdrift/tau]
    const double nElectrons = nIonizedElectrons * lifetimecorrection;
    //std::cout << "After lifetime, " << nElectrons << " electrons." << std::endl;
    
    // Longitudinal & transverse diffusion sigma (cm)
    double SqrtT    = std::sqrt(TDrift);
    double LDiffSig = SqrtT * fLDiff_const;
    double TDiffSig = SqrtT * fTDiff_const;
    double electronclsize = fElectronClusterSize;
    
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
    
    double fDriftClusterPos[3];

    // fix the number of electrons in the last cluster, that has a smaller size
    nElDiff.back() = nElectrons - (nClus-1)*electronclsize;
    
    for(size_t xx = 0; xx < nElDiff.size(); ++xx)
    {
        if(nElectrons > 0) nEnDiff[xx] = energy/nElectrons*nElDiff[xx];
        else               nEnDiff[xx] = 0.;
    }

    // Smear drift times by longitudinal diffusion
    if (LDiffSig > 0.0)
        randGauss.fireArray( nClus, &longDiff[0], 0., LDiffSig);
    else
        longDiff.assign(nClus, 0.0);
    
    if (TDiffSig > 0.0)
    {
        // Smear the coordinates in plane perpendicular to drift direction by the transverse diffusion
        randGauss.fireArray( nClus, &transDiff1[0], avegagetransversePos1, TDiffSig);
        randGauss.fireArray( nClus, &transDiff2[0], avegagetransversePos2, TDiffSig);
    }
    else
    {
        transDiff1.assign(nClus, avegagetransversePos1);
        transDiff2.assign(nClus, avegagetransversePos2);
    }
    
    // make a collection of electrons for each plane
    for(size_t p = 0; p < tpcGeo.Nplanes(); ++p)
    {
        fDriftClusterPos[driftcoordinate] = tpcGeo.PlaneLocation(p)[driftcoordinate];
        
        // Drift nClus electron clusters to the induction plane
        for(int k = 0; k < nClus; ++k)
        {
            // Correct drift time for longitudinal diffusion and plane
            double TDiff = TDrift + longDiff[k] * fRecipDriftVel[0];
            
            // Take into account different Efields between planes
            // Also take into account special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is the second wire plane
            for (size_t ip = 0; ip<p; ++ip)
            {
                TDiff += (tpcGeo.PlaneLocation(ip+1)[driftcoordinate] - tpcGeo.PlaneLocation(ip)[driftcoordinate]) * fRecipDriftVel[(tpcGeo.Nplanes() == 2 && driftcoordinate == 0)?ip+2:ip+1];
            }
            
            fDriftClusterPos[transversecoordinate1] = transDiff1[k];
            fDriftClusterPos[transversecoordinate2] = transDiff2[k];
            
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
                raw::ChannelID_t channel = fGeometry->NearestChannel(fDriftClusterPos, p, tpc, cryostat);
                
                //          std::cout << "fDriftClusterPos[0]: " << fDriftClusterPos[0] << "\t fDriftClusterPos[1]: " << fDriftClusterPos[1] << "\t fDriftClusterPos[2]: " << fDriftClusterPos[2] << std::endl;
                //          std::cout << "channel: " << channel << std::endl;
                
                //std::cout << "\tgot channel " << channel << " for cluster " << k << std::endl;
                
                /// \todo check on what happens if we allow the tdc value to be
                /// \todo beyond the end of the expected number of ticks
                // Add potential decay/capture/etc delay effect, simTime.
                auto const simTime = energyDeposit.Time();
                unsigned int tdc = fClock.Ticks(fTimeService->G4ToElecTime(TDiff + simTime));
                
                // Find whether we already have this channel in our map.
                ChannelMap_t& channelDataMap = channelMap[cryostat][tpc];
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
                    bookKeeping.channelIndex = simChannelVec.size();
                    simChannelVec.emplace_back( channel );
                    channelIndex = bookKeeping.channelIndex;
                    
                    // Save the pointer to the newly-created
                    // sim::SimChannel.
                    //channelPtr = &(channels->back());
                    //bookKeeping.channelPtr = channelPtr;
                    
                    // Initialize a vector with the index of the step that
                    // created this channel.
                    bookKeeping.stepList.push_back( edIndex );
                    
                    // Save the bookkeeping information for this channel.
                    channelDataMap[channel] = bookKeeping;
                }
                else {
                    // We've created this SimChannel for a previous energy
                    // deposit. Get its address.
                    
                    //std::cout << "\tHave seen this channel before." << std::endl;
                    
                    auto& bookKeeping = search->second;
                    channelIndex = bookKeeping.channelIndex;
                    //channelPtr = bookKeeping.channelPtr;
                    
                    // Has this step contributed to this channel before?
                    auto& stepList = bookKeeping.stepList;
                    auto stepSearch = std::find(stepList.begin(), stepList.end(), edIndex );
                    if ( stepSearch == stepList.end() ) {
                        // No, so add this step's index to the list.
                        stepList.push_back( edIndex );
                    }
                }
                
                sim::SimChannel* channelPtr = &(simChannelVec.at(channelIndex));
                
                //if(!channelPtr) std::cout << "\tUmm...ptr is NULL?" << std::endl;
                //else std::cout << "\tChannel is " << channelPtr->Channel() << std::endl;
                // Add the electron clusters and energy to the
                // sim::SimChannel
                channelPtr->AddIonizationElectrons(energyDeposit.TrackID(),
                                                   tdc,
                                                   nElDiff[k],
                                                   xyz,
                                                   nEnDiff[k]);
                
                if(fStoreDriftedElectronClusters)
                    simDriftedElectronClustersVec.push_back(sim::SimDriftedElectronCluster(nElDiff[k],
                                                                                           TDiff + simTime,        // timing
                                                                                           {mp.X(),mp.Y(),mp.Z()}, // mean position of the deposited energy
                                                                                           {fDriftClusterPos[0],fDriftClusterPos[1],fDriftClusterPos[2]}, // final position of the drifted cluster
                                                                                           {LDiffSig,TDiffSig,TDiffSig}, // Longitudinal (X) and transverse (Y,Z) diffusion
                                                                                           nEnDiff[k], //deposited energy that originated this cluster
                                                                                           energyDeposit.TrackID()) );
            }
            catch(cet::exception &e) {
                mf::LogWarning("SimDriftElectrons") << "unable to drift electrons from point ("
                << xyz[0] << "," << xyz[1] << "," << xyz[2]
                << ") with exception " << e;
            } // end try to determine channel
        } // end loop over clusters
    } // end loop over planes
    return;
}

DEFINE_ART_CLASS_TOOL(ElectronDrift2D)
}
