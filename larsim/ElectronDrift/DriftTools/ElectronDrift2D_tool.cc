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
                        ChannelMapByCryoTPC&,
                        ChannelToSimEnergyMap&) override;         // output candidate hits

private:
    
    std::string                       fResponseFileName;
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
    
    // Define the structure to contain the look up tables
    // We'll need to divide by the 1D distance to the nearest wire,
    // then by the plane (we'll loop over all three), then by the
    // wire to add charge to
    using ResponseVec                 = std::vector<float>;
    using WireResponseVecMap          = std::unordered_map<int,ResponseVec>;
    using DistWireResponseVecMap      = std::unordered_map<int,WireResponseVecMap>;
    using PlaneDistWireResponseVecMap = std::unordered_map<int,DistWireResponseVecMap>;
    using DistDecoderVec              = std::vector<float>;
    
    PlaneDistWireResponseVecMap       fPlaneDistWireResponseVecMap;
    const DistDecoderVec              fDistDecoderVec = {-0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135};
    float                             fDriftPosBinSize = 0.03;
    
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
    fResponseFileName             = pset.get<std::string>("ResponseFileName",             "t600_responses_2D_v0.0.root");
    fStoreDriftedElectronClusters = pset.get< bool      >("StoreDriftedElectronClusters", false);

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
    
    // The actual work for this model is done by look up tables which recover in the form of histograms
    std::unordered_map<std::string, int> planeDecoderMap = {{"I1",0},{"I2",1},{"C",2}};

    // Recover the input field response histogram
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    
    if (!searchPath.find_file(fResponseFileName, fullFileName))
        throw cet::exception("ElectronDrift2D::configure") << "Can't find input file: '" << fullFileName << "'\n";
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("ElectronDrift2D::configure") << "Unable to open input file: " << fullFileName << std::endl;
    
    TIter next(inputFile.GetListOfKeys());
    
    while(TKey* key = (TKey*)next())
    {
        //std::cout << "key: " << key->GetName() << ", class:" << key->GetClassName() << std::endl;
        std::string className(key->GetClassName());
        
        if (className.find("TH1") < className.npos)
        {
            TH1F* histPtr = dynamic_cast<TH1F*>(inputFile.Get(key->GetName()));
            
            // Expect the histogram name to have the form
            // "FieldResponse_(plane)_(wire)_(pos)
            // Where "plane" will be a string like "C", "I1", or "I2"
            // "wire" will be 0-9
            // "pos" is a position in cm
            // So... we need to decode it here
            std::string histName(histPtr->GetName());
            
            size_t planeStartPos = histName.find("_") + 1;
            size_t wireStartPos  = histName.find("_", planeStartPos) + 1;
            size_t posStartPos   = histName.find("_", wireStartPos) + 1;
            
            std::string planeString = histName.substr(planeStartPos, wireStartPos - planeStartPos - 1);
            int         wireNum     = std::stoi(histName.substr(wireStartPos, posStartPos - wireStartPos - 1));
            float       wirePos     = std::stof(histName.substr(posStartPos, className.size() - posStartPos));
            
            std::vector<float>::const_iterator distIterator = std::find_if(fDistDecoderVec.begin(),fDistDecoderVec.end(),[wirePos](const auto& val){return std::abs(val-wirePos) < std::numeric_limits<float>::epsilon();});
            
            if (distIterator == fDistDecoderVec.end())
                throw cet::exception("ElectronDrift2D::configure") << "Unable to decode drift ditsance: " << wirePos << std::endl;

            int wirePosIdx = std::distance(fDistDecoderVec.begin(),distIterator);
            
            std::cout << "--> Recovered hist: " << histName << ", title: " << histPtr->GetTitle() << ", plane: " << planeString << ", wire: " << wireNum << ", pos: " << wirePos << ", idx: " << wirePosIdx << ", " << fDriftPosBinSize << std::endl;
            
            // The histograms will contain the charge deposited on the wire as a function of time in us
            // We need to convert this to ticks...
            int   numBins      = histPtr->GetXaxis()->GetNbins();
            float histBinWidth = histPtr->GetBinWidth(1);
            float histLowEdge  = histPtr->GetXaxis()->GetXmin();
            float firstBinLow  = histPtr->GetBinCenter(1) - 0.5 * histBinWidth;
            float histMaxEdge  = histPtr->GetXaxis()->GetXmax();
            float lastBinMax   = histPtr->GetBinCenter(numBins) + 0.5 * histBinWidth;
            
            float samplingRate = detprop->SamplingRate() * 1.e-3;                // We want this in us/bin
            float rateRatio    = samplingRate / histBinWidth;                    // This gives the ratio of time bins for full readout to response bins
            int   binsPerTick  = int(std::round(rateRatio));
            int   nBinsTickVec = numBins / binsPerTick;
            
            // Recover the vector we want for this drift position, then plane, then wire
            DistWireResponseVecMap&  distWireResponseVecMap = fPlaneDistWireResponseVecMap[planeDecoderMap[planeString]];
            WireResponseVecMap&      wireResponseVecMap      = distWireResponseVecMap[wirePosIdx];
            ResponseVec&             responseVec             = wireResponseVecMap[wireNum];
            
            responseVec.resize(nBinsTickVec,0.);

            std::cout << "    hist has " << numBins << " bins, min: " << histLowEdge << "/" << firstBinLow << ", max: " << histMaxEdge << "/" << lastBinMax << ", bin Size: " << (histMaxEdge - histLowEdge - histPtr->GetBinWidth(1))/float(numBins) << "/" << histPtr->GetBinWidth(1) << ", sampling rate: " << samplingRate << ", rat: " << rateRatio << ", binsPerTick: " << binsPerTick << std::endl;
            
            for(int binIdx = 1; binIdx <= numBins; binIdx++)
            {
                size_t responseVecIdx = (binIdx - 1) / binsPerTick;
                responseVec[responseVecIdx] += histPtr->GetBinContent(binIdx);
            }
        }
    }
    
    inputFile.Close();
        

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
                                     ChannelMapByCryoTPC&                         channelMap,
                                     ChannelToSimEnergyMap&                       channelToSimEnergyMap)
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
    
    // Now have everything we need to start filling SimChannel output
    // We need to loop over planes since we'll then need to find the doca to the wire in the plane for each cluster
    for(const auto& planeItr : fPlaneDistWireResponseVecMap)
    {
        int thisPlane = planeItr.first;
        
        fDriftClusterPos[driftcoordinate] = tpcGeo.PlaneLocation(thisPlane)[driftcoordinate];
        
        // Now loop over each of the clusters we have created
        for(int clusIdx = 0; clusIdx < nClus; clusIdx++)
        {
            // Correct drift time for longitudinal diffusion
            double TDiff = TDrift + longDiff[clusIdx] * fRecipDriftVel[0];
            
            // And correct position for transverse...
            fDriftClusterPos[transversecoordinate1] = transDiff1[clusIdx];
            fDriftClusterPos[transversecoordinate2] = transDiff2[clusIdx];
            
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
                //raw::ChannelID_t channel = fGeometry->NearestChannel(fDriftClusterPos, thisPlane, tpc, cryostat);
                geo::WireID wireID = fGeometry->NearestWireID(fDriftClusterPos, thisPlane, tpc, cryostat);
                raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(wireID);

                // Find the transverse doca to this wire...
                Eigen::Vector3d wireEndPoint1;
                Eigen::Vector3d wireEndPoint2;
                
                fGeometry->WireEndPoints(wireID, &wireEndPoint1[0], &wireEndPoint2[0]);
                
                // Want the distance to the wire in the plane transverse to the drift
                Eigen::Vector2f wireEndPoint2D1(wireEndPoint1[transversecoordinate1],wireEndPoint1[transversecoordinate2]);
                Eigen::Vector2f wireEndPoint2D2(wireEndPoint2[transversecoordinate1],wireEndPoint2[transversecoordinate2]);
                Eigen::Vector2f driftPos2D(fDriftClusterPos[transversecoordinate1],fDriftClusterPos[transversecoordinate2]);
                
                Eigen::Vector2f wireDirVec  = wireEndPoint2D2 - wireEndPoint2D1;
                Eigen::Vector2f driftPosVec = driftPos2D - wireEndPoint2D1;
                
                float wireLength   = wireDirVec.norm();
                
                wireDirVec = wireDirVec.normalized();
                
                float arcLenToDoca = wireDirVec.dot(driftPosVec);
                
                if (arcLenToDoca > wireLength || arcLenToDoca < 0.)
                    std::cout << "***** ACK! ******" << std::endl;
                
                Eigen::Vector2f docaPos = wireEndPoint2D1 + arcLenToDoca * wireDirVec;
                Eigen::Vector2f docaVec = driftPos2D - docaPos;
                
                float distWireToPos = docaVec.norm();
                
                docaVec.normalize();
                
                if (wireDirVec[0] * docaVec[1] - wireDirVec[1] * docaVec[0] < 0.) distWireToPos = -distWireToPos;
                
                int distWireToPosIdx = int(std::round((distWireToPos - fDistDecoderVec[0]) / fDriftPosBinSize));
                
                distWireToPosIdx = std::max(0,std::min(distWireToPosIdx,int(fDistDecoderVec.size())));
                
                // Now we can get the look up tables we need to complete the task
                const DistWireResponseVecMap& distWireResponseVecMap = planeItr.second;
                const WireResponseVecMap&     wireResponseVecMap     = distWireResponseVecMap.at(distWireToPosIdx);
                
                // Skip looping over all wires for the moment... take the central wire to start
                const ResponseVec& responseVec = wireResponseVecMap.at(0);

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
                
                // Loop over the ticks in the response vector
                unsigned int tdcIdx(0);
                
                for(const auto& charge : responseVec)
                {
                    channelPtr->AddIonizationElectrons(energyDeposit.TrackID(),
                                                       tdc + tdcIdx++,
                                                       nElDiff[clusIdx] * charge,
                                                       xyz,
                                                       nEnDiff[clusIdx]);
                }
                
                if(fStoreDriftedElectronClusters)
                    simDriftedElectronClustersVec.push_back(sim::SimDriftedElectronCluster(nElDiff[clusIdx],
                                                                                            TDiff + simTime,        // timing
                                                                                            {mp.X(),mp.Y(),mp.Z()}, // mean position of the deposited energy
                                                                                            {fDriftClusterPos[0],fDriftClusterPos[1],fDriftClusterPos[2]}, // final position of thedrifted  cluster
                                                                                            {LDiffSig,TDiffSig,TDiffSig}, // Longitudinal (X) and transverse (Y,Z) diffusion
                                                                                            nEnDiff[clusIdx], //deposited energy that originated this cluster
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
