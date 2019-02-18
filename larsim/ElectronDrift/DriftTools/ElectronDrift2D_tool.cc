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
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larsim/ElectronDrift/DriftTools/IDiffusionTool.h"

#include "lardataobj/Simulation/Waveform.h"

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
    
    // Allow our tools to declare data products they plan to output to event store
    void produces(art::EDProducer&) override;
    
    // Set up output data products
    void setupOutputDataProducts() override;

    // Search for candidate hits on the input waveform
    void driftElectrons(const size_t,
                        const sim::SimEnergyDeposit&,
                        CLHEP::RandGauss&,
                        ChannelToSimEnergyMap&) override;         // output candidate hits
    
    // Output the data products
    void put(art::Event&) override;

private:
    
    int getTransverseBin(const geo::WireID&, const Eigen::Vector3d&, int, int) const;
    sim::SimChannel* getSimChannel(ChannelToSimEnergyMap&, std::vector<sim::SimChannel>&, const size_t, const raw::ChannelID_t) const;
    
    std::string                                fResponseFileName;
    bool                                       fStoreDriftedElectronClusters;

    const detinfo::DetectorClocks*             fTimeService;

    double                                     fRecipDriftVel[3];
    
    // Define the structure to contain the look up tables
    // We'll need to divide by the 1D distance to the nearest wire,
    // then by the plane (we'll loop over all three), then by the
    // wire to add charge to
    using ResponseVec                 = std::vector<float>;
    using WireResponseVecMap          = std::unordered_map<int,ResponseVec>;
    using DistWireResponseVecMap      = std::unordered_map<int,WireResponseVecMap>;
    using PlaneDistWireResponseVecMap = std::unordered_map<int,DistWireResponseVecMap>;
    using DistDecoderVec              = std::vector<float>;
    using WireOffsetMap               = std::unordered_map<int, int>;
    using DistWireOffsetMap           = std::unordered_map<int, WireOffsetMap>;
    using PlaneDistWireOffsetMap      = std::unordered_map<int, DistWireOffsetMap>;
    
    PlaneDistWireResponseVecMap                 fPlaneDistWireResponseVecMap;
    PlaneDistWireOffsetMap                      fPlaneDistWireOffsetMap;
    const DistDecoderVec                        fDistDecoderVec   = {-0.135, -0.105, -0.075, -0.045, -0.015, 0.015, 0.045, 0.075, 0.105, 0.135};
    double                                      fDriftPosBinSize  = 0.03;
    double                                      fDriftPlaneOffset = 10.;
    double                                      fTicksPerNanosecond;
    
    size_t                                      fMaxWaveformTicks = 2*4096;
    
    // Define structures for holding output waveforms
    using WaveformVec             = std::vector<float>;
//    using WaveformVec             = lar::sparse_vector<float>;
    using ChannelToWaveformVecMap = std::unordered_map<raw::ChannelID_t, WaveformVec>;
    
    ChannelToWaveformVecMap                     fChannelToWaveformVecMap;
    
    // Output objects
    std::unique_ptr<std::vector<sim::SimChannel>>                fSimChannelVec;
    std::unique_ptr<std::vector<sim::Waveform>>                  fSimWaveformVec;
    std::unique_ptr<std::vector<sim::SimDriftedElectronCluster>> fSimDriftedElectronClustersVec;

    // The tool to handle the diffusion
    std::unique_ptr<detsim::IDiffusionTool>     fDiffusionTool; ///< Tool for handling the diffusion during the drift

    const geo::GeometryCore*                    fGeometry = lar::providerFrom<geo::Geometry>();
    ::detinfo::ElecClock                        fClock;     ///< TPC electronics clock
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
    
    for (int i = 0; i<3; ++i)
    {
        //  Drift velocity as returned by the DetectorProperties service assumed to be in cm/us. Multiply by 1.e-3 to convert into LArSoft standard velocity units, cm/ns;
        double driftVelocity = detprop->DriftVelocity(detprop->Efield(i), detprop->Temperature())*1.e-3;
        
        fRecipDriftVel[i] = 1./driftVelocity;
    }
    
    double samplingRate = detprop->SamplingRate() * 1.e-3;                // We want this in us/bin
    
    fTicksPerNanosecond = 1. / (1.e3 * samplingRate);   // want this in ticks/ns

    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock       = fTimeService->TPCClock();
    
    fDiffusionTool = art::make_tool<detsim::IDiffusionTool>(pset.get<fhicl::ParameterSet>("DiffusionTool"));

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
    
    float histIntegral(0.);
    
    while(TKey* key = (TKey*)next())
    {
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
            
            // Assumption: for a given set of responses for a plane and position, the hists are read in wire number order, the first will be zero and be
            // the central wire. So we want to recover the integral for that hist
            if (planeString == "C" && wireNum == 0) histIntegral = histPtr->Integral();
            
            // Now normalize the histogram by the cached integral for wire 0
            histPtr->Scale(1./histIntegral);
            
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
            
            ResponseVec::iterator firstNonZeroItr = std::find_if(responseVec.begin(),responseVec.end(),[](const auto& val){return std::abs(val) > 1.e-6;}); //std::numeric_limits<float>::epsilon();});
            
            fPlaneDistWireOffsetMap[planeDecoderMap[planeString]][wirePosIdx][wireNum] = std::distance(responseVec.begin(),firstNonZeroItr);
            
            std::cout << "    ---> Starting index for plane: " << planeString << ", dist: " << wirePos << "/" << wirePosIdx << ", wireNum: " << wireNum << ", is: " << fPlaneDistWireOffsetMap[planeDecoderMap[planeString]][wirePosIdx][wireNum] << std::endl;
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
    
void ElectronDrift2D::produces(art::EDProducer& producer)
{
    producer.produces<std::vector<sim::SimChannel>>();
    producer.produces<std::vector<sim::Waveform>>();
    
    if(fStoreDriftedElectronClusters) producer.produces< std::vector<sim::SimDriftedElectronCluster> >();
    
    return;
}

// Set up output data products
void ElectronDrift2D::setupOutputDataProducts()
{
    fSimChannelVec  = std::make_unique<std::vector<sim::SimChannel>>(std::vector<sim::SimChannel>());
    fSimWaveformVec = std::make_unique<std::vector<sim::Waveform>>(std::vector<sim::Waveform>());

    if(fStoreDriftedElectronClusters)
        fSimDriftedElectronClustersVec = std::make_unique<std::vector<sim::SimDriftedElectronCluster>>(std::vector<sim::SimDriftedElectronCluster>());
    
    return;
}

void ElectronDrift2D::driftElectrons(const size_t                 edIndex,
                                     const sim::SimEnergyDeposit& energyDeposit,
                                     CLHEP::RandGauss&            randGauss,
                                     ChannelToSimEnergyMap&       channelToSimEnergyMap)
{
    // First call the tool to do the drifting
    if (fDiffusionTool->getDiffusionVec(energyDeposit, randGauss))
    {
        // "xyz" is the position of the energy deposit in world
        // coordinates. Note that the units of distance in
        // sim::SimEnergyDeposit are supposed to be cm.
        auto const   mp     = energyDeposit.MidPoint();
        double const xyz[3] = { mp.X(), mp.Y(), mp.Z() };
        
        // Recover drift distance and drift time
        double driftDistance = fDiffusionTool->getDriftDistance();
        double driftTime     = fDiffusionTool->getDriftTime();
        
        // We should be able to determine the first index into the lookup tables before we start the loops over planes/etc.
        // Note that the look up tables start some distance (fDriftPlaneOffset) before the first plane and so we can check
        // if the current position is inside that range or outside...
        size_t startLookupIdx(0);
        double startDriftTime = fDriftPlaneOffset * fRecipDriftVel[0];
        
        if (fDriftPlaneOffset - driftDistance > 0.)
            startLookupIdx = std::round((fDriftPlaneOffset - driftDistance) * fRecipDriftVel[0] * fTicksPerNanosecond);
        
        // Note that "TDrift" is the time to drift electrons from the charge deposit position to the first sense plane
        // The look up tables start at fDriftPlaneOffset distance before that plane so to get the right time we need
        // to reset the total drift back to this position.
        driftTime -= startDriftTime;
        
        // Recover the coordinate indices
        int driftCoordinate       = fDiffusionTool->getDriftCoordinate();
        int transverseCoordinate1 = fDiffusionTool->getTransverseCoordinate1();
        int transverseCoordinate2 = fDiffusionTool->getTransverseCoordinate2();
        
        // Get a pointer to the TPC geometry
        const geo::TPCGeo* tpcGeo = fDiffusionTool->getTPCGeo();

        // Now have everything we need to start filling SimChannel output
        // We need to loop over planes since we'll then need to find the doca to the wire in the plane for each cluster
        for(const auto& planeItr : fPlaneDistWireResponseVecMap)
        {
            int    thisPlane(planeItr.first);
            Eigen::Vector3d driftClusterPos;
            
            driftClusterPos[driftCoordinate] = tpcGeo->PlaneLocation(thisPlane)[driftCoordinate];
            
            // Recover offset map for this plane
            const DistWireOffsetMap& distWireOffsetMap = fPlaneDistWireOffsetMap.at(thisPlane);
            
            // Now loop over each of the clusters we have created
            for (const auto &cluster : fDiffusionTool->getDiffusionValsVec())
            {
                // Define function to handle adding the charge to the waveform
                float numElectrons   = cluster.clusterNumElectrons;
//                auto  sumWaveformVec = [&numElectrons](const float& response, const float& waveform){return waveform + numElectrons * response;};
                
                // Correct drift time for longitudinal diffusion
                double TDiff = driftTime + cluster.longitudinalDiffusion * fRecipDriftVel[0];
                
                // And correct position for transverse...
                driftClusterPos[transverseCoordinate1] = cluster.transverseDiffusion_1;
                driftClusterPos[transverseCoordinate2] = cluster.transverseDiffusion_2;
                
                /// \todo check on what happens if we allow the tdc value to be
                /// \todo beyond the end of the expected number of ticks
                // Add potential decay/capture/etc delay effect, simTime.
                auto const   simTime = energyDeposit.Time();
                unsigned int tdc     = fClock.Ticks(fTimeService->G4ToElecTime(TDiff + simTime));
                
                // Skip this cluster if timing is already past the length of the buffer...
                if (tdc > fMaxWaveformTicks) continue;

                // grab the nearest channel to the fDriftClusterPos position
                geo::WireID      wireID;
                raw::ChannelID_t channel;
                
                try
                {
                    const geo::TPCID& tpcID = tpcGeo->ID();
                    
                    wireID  = fGeometry->NearestWireID(driftClusterPos.data(), thisPlane, tpcID.TPC, tpcID.Cryostat);
                    channel = fGeometry->PlaneWireToChannel(wireID);
                }
                catch(cet::exception &e)
                {
                    mf::LogWarning("ElectronDrift2D") << "unable to drift electrons from point ("
                    << xyz[0] << "," << xyz[1] << "," << xyz[2]
                    << ") with exception " << e;
                    continue;
                } // end try to determine channel
                
                // Find the transverse doca to this wire...
                int distWireToPosIdx = getTransverseBin(wireID, driftClusterPos, transverseCoordinate1, transverseCoordinate2);

                // Now we can get the look up tables we need to complete the task
                const DistWireResponseVecMap& distWireResponseVecMap = planeItr.second;
                
                // Recover the wire information for this plane
                int maxWire     = tpcGeo->Plane(wireID.Plane).Nwires() - 1;
                int centralWire = wireID.Wire;
                
                // The WireResponcseVecMap gives the mapping of wires to the positive side of the central wire...
                // We also need to deal with the wires to the negative side so we do an outer loop to get both
                for(int side = -1; side < 2; side += 2)
                {
                    int transverseIdx = side > 0 ? distWireToPosIdx : fDistDecoderVec.size() - distWireToPosIdx - 1;
                    
                    // Use this to look up the response
                    const WireResponseVecMap& wireResponseVecMap = distWireResponseVecMap.at(transverseIdx);
                    
                    // Recover the map of offsets by wire number
                    const WireOffsetMap& wireOffsetMap = distWireOffsetMap.at(transverseIdx);

                    // Initiate the internal loop over wires
                    for(const auto& wireResponseVec : wireResponseVecMap)
                    {
                        // Don't double count the cental wire
                        if (side < 0 && wireResponseVec.first == 0) continue;
                        
                        // Check for valid wire
                        int curWire = centralWire + side * wireResponseVec.first;
                        
                        if (curWire < 0 || curWire > maxWire) continue;
                        
                        // Skip looping over all wires for the moment... take the central wire to start
                        const ResponseVec& responseVec = wireResponseVec.second;
                        
                        // Recover the SimChannel and handle the bookeeping
                        size_t locChannel = channel + side * wireResponseVec.first;
                        
                        ChannelToWaveformVecMap::iterator waveformVecItr = fChannelToWaveformVecMap.find(locChannel);
                        
                        if(waveformVecItr == fChannelToWaveformVecMap.end())
                            waveformVecItr = fChannelToWaveformVecMap.insert(std::pair<raw::ChannelID_t, WaveformVec>(locChannel,WaveformVec(fMaxWaveformTicks,0.))).first;
//                            waveformVecItr = fChannelToWaveformVecMap.insert(std::pair<raw::ChannelID_t, WaveformVec>(locChannel,lar::sparse_vector<float>())).first;

                        WaveformVec& waveformVec = waveformVecItr->second;
                        
                        // The offset map contains the starting indices in the look up table for the first non-zero (measureable) response so we
                        // want the starting index to be at least that... or later if the track start is later
                        size_t responseFirstIdx = wireOffsetMap.at(wireResponseVec.first) > startLookupIdx ? wireOffsetMap.at(wireResponseVec.first) : startLookupIdx;
                        size_t lastLookupIdx    = responseVec.size();
                        
                        if (tdc + lastLookupIdx >= fMaxWaveformTicks) lastLookupIdx = fMaxWaveformTicks - tdc;
    
                        // It is a bit of a mystery to me why std::transform is slow than a for loop here...
//                        std::transform(waveformVec.begin() + tdc + responseFirstIdx,
//                                       waveformVec.begin() + tdc + lastLookupIdx,
//                                       responseVec.begin() + responseFirstIdx,
//                                       waveformVec.begin() + tdc + responseFirstIdx,
//                                       sumWaveformVec);

                        // Update the contents of the waveform vector
                        for(size_t tick = responseFirstIdx; tick < lastLookupIdx; tick++)
                            waveformVec[tdc + tick] += numElectrons * responseVec[tick];
                        
//                        std::vector<float> tempVec(lastLookupIdx - responseFirstIdx);
//
//                        for(size_t tick = 0; tick < lastLookupIdx - responseFirstIdx; tick++) tempVec[tick] = cluster.clusterNumElectrons * responseVec[responseFirstIdx + tick];
//
//                        std::vector<std::pair<size_t,size_t>> pairVector;
//
//                        bool rangeOverlap = !waveformVec.empty();    // If the waveformVec is empty then treat as overlapped
//
//                        // Check for range overlaps
//                        if (!rangeOverlap)
//                        {
//                            for(size_t rangeIdx = 0; rangeIdx < waveformVec.n_ranges(); rangeIdx++)
//                            {
//                                const lar::sparse_vector<float>::datarange_t& range = waveformVec.range(rangeIdx);
//
//                                if (tdc + lastLookupIdx >= range.offset && tdc + responseFirstIdx <= range.last)
//                                {
//                                    rangeOverlap = true;
//                                    break;
//                                }
//                            }
//                        }
//
//                        if (!rangeOverlap) waveformVec.add_range(tdc + responseFirstIdx, tempVec);
//                        else               waveformVec.combine_range(tdc + responseFirstIdx, tempVec, std::plus<float>());
                    }
                }
                
                // Recoer the SimChannel to fill
                sim::SimChannel* channelPtr = getSimChannel(channelToSimEnergyMap, *fSimChannelVec, edIndex, channel);
                
                // And fill it...
                channelPtr->AddIonizationElectrons(energyDeposit.TrackID(),
                                                   tdc,
                                                   cluster.clusterNumElectrons,
                                                   xyz,
                                                   cluster.clusterEnergy);

                if(fStoreDriftedElectronClusters)
                    fSimDriftedElectronClustersVec->push_back(sim::SimDriftedElectronCluster(cluster.clusterNumElectrons,
                                                                                             TDiff + simTime,        // timing
                                                                                             {mp.X(),mp.Y(),mp.Z()}, // mean position of the deposited energy
                                                                                             {driftClusterPos[0],driftClusterPos[1],driftClusterPos[2]}, // final position of thedrifted  cluster
                                                                                             {fDiffusionTool->getLongitudinalDiffusionSig(),
                                                                                              fDiffusionTool->getTransverseDiffusionSig(),
                                                                                              fDiffusionTool->getTransverseDiffusionSig()}, // Longitudinal (X) and transverse (Y,Z) diffusion
                                                                                             cluster.clusterEnergy, //deposited energy that originated this cluster
                                                                                             energyDeposit.TrackID()) );
            } // end loop over clusters
        } // end loop over planes
    }

    return;
}
    
int ElectronDrift2D::getTransverseBin(const geo::WireID& wireID, const Eigen::Vector3d& driftClusterPos, int transverseCoordinate1, int transverseCoordinate2) const
{
    // Find the transverse doca to this wire...
    Eigen::Vector3d wireEndPoint1;
    Eigen::Vector3d wireEndPoint2;
    
    fGeometry->WireEndPoints(wireID, &wireEndPoint1[0], &wireEndPoint2[0]);
    
    // Want the distance to the wire in the plane transverse to the drift
    Eigen::Vector2f wireEndPoint2D1(wireEndPoint1[transverseCoordinate1],wireEndPoint1[transverseCoordinate2]);
    Eigen::Vector2f wireEndPoint2D2(wireEndPoint2[transverseCoordinate1],wireEndPoint2[transverseCoordinate2]);
    Eigen::Vector2f driftPos2D(driftClusterPos[transverseCoordinate1],driftClusterPos[transverseCoordinate2]);
    
    Eigen::Vector2f wireDirVec  = wireEndPoint2D2 - wireEndPoint2D1;
    Eigen::Vector2f driftPosVec = driftPos2D - wireEndPoint2D1;
    
    float wireLength = wireDirVec.norm();
    
    wireDirVec = wireDirVec.normalized();
    
    float arcLenToDoca = wireDirVec.dot(driftPosVec);
    
    if (arcLenToDoca > wireLength || arcLenToDoca < 0.) arcLenToDoca = std::min(wireLength,std::max(float(0.),arcLenToDoca));
    
    Eigen::Vector2f docaPos = wireEndPoint2D1 + arcLenToDoca * wireDirVec;
    Eigen::Vector2f docaVec = driftPos2D - docaPos;
    
    float distWireToPos = docaVec.norm();
    
    docaVec.normalize();
    
    if (wireDirVec[0] * docaVec[1] - wireDirVec[1] * docaVec[0] < 0.) distWireToPos = -distWireToPos;
    
    int distWireToPosIdx = int(std::round((distWireToPos - fDistDecoderVec[0]) / fDriftPosBinSize)) - 1;
    
    return std::max(0,std::min(distWireToPosIdx,int(fDistDecoderVec.size()-1)));
}

sim::SimChannel* ElectronDrift2D::getSimChannel(ChannelToSimEnergyMap& channelToSimEnergyMap, std::vector<sim::SimChannel>& simChannelVec, const size_t edIndex, const raw::ChannelID_t channel) const
{
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

    return &(simChannelVec.at(channelIndex));
}
    
void ElectronDrift2D::put(art::Event& event)
{
    // Want to retrieve our waveforms but want them in order of channel number...
    std::map<raw::ChannelID_t, WaveformVec> orderedChanWaveMap(fChannelToWaveformVecMap.begin(),fChannelToWaveformVecMap.end());

    for(auto& chanWaveItr : orderedChanWaveMap)
    {
        raw::ChannelID_t channel = chanWaveItr.first;
        
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        
        WaveformVec& waveformVec = chanWaveItr.second;

        // Want to drop zeroes off the vector... so search for the first non-zero bin
        WaveformVec::iterator nonZeroElemItr = std::find_if(waveformVec.begin(),waveformVec.end(),[](const auto& val){return std::abs(val) > std::numeric_limits<float>::epsilon();});

        // If the entire vector is zero (unlikely) then skip
        if (nonZeroElemItr != waveformVec.end())
        {
            // Create an empty sparse vector
            lar::sparse_vector<float> sparseWaveformVec = lar::sparse_vector<float>();

            // Get the first offset
            size_t firstNonZeroElem = std::distance(waveformVec.begin(),nonZeroElemItr);

            // Make a current "last" element iterator
            WaveformVec::iterator zeroElemItr = std::find_if(nonZeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) < std::numeric_limits<float>::epsilon();});

            // Look for small gaps
            WaveformVec::iterator nextNonZeroElemItr = std::find_if(zeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) > std::numeric_limits<float>::epsilon();});

            // Want to loop over all possible valid ranges
            while(nonZeroElemItr != waveformVec.end())
            {
                // First we advance
                while(nextNonZeroElemItr != waveformVec.end() && std::distance(zeroElemItr,nextNonZeroElemItr) < 6)
                {
                    zeroElemItr        = std::find_if(nextNonZeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) < std::numeric_limits<float>::epsilon();});
                    nextNonZeroElemItr = std::find_if(zeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) > std::numeric_limits<float>::epsilon();});
                }

                // Add the range
                sparseWaveformVec.add_range(firstNonZeroElem,nonZeroElemItr,zeroElemItr);

                // Set up the next round
                nonZeroElemItr   = nextNonZeroElemItr;
                firstNonZeroElem = std::distance(waveformVec.begin(),nonZeroElemItr);

                zeroElemItr        = std::find_if(nonZeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) < std::numeric_limits<float>::epsilon();});
                nextNonZeroElemItr = std::find_if(zeroElemItr, waveformVec.end(),[](const auto& val){return std::abs(val) > std::numeric_limits<float>::epsilon();});
            }

            //fSimWaveformVec->emplace_back(waveformVec, channel, fGeometry->View(wids[0]));
            fSimWaveformVec->emplace_back(sparseWaveformVec, channel, fGeometry->View(wids[0]));
        }
    }
    
    event.put(std::move(fSimChannelVec));
    event.put(std::move(fSimWaveformVec));
    
    if(fStoreDriftedElectronClusters)
        event.put(std::move(fSimDriftedElectronClustersVec));
    
    // Clear for the next event
    fChannelToWaveformVecMap.clear();
    
    return;
}

    
DEFINE_ART_CLASS_TOOL(ElectronDrift2D)
}
