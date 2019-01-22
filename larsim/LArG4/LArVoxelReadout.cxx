////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.cxx
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

// C/C++ standard library
#include <cstdio> // std::sscanf()
#include <cmath> // std::ceil()
#include <string>
#include <map>
#include <utility> // std::move()
#include <cassert>

// GEANT
#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4TouchableHistory.hh"
#include "Geant4/G4TouchableHandle.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4VPhysicalVolume.hh"

// framework libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft code
#include "larsim/LArG4/LArVoxelReadout.h"
#include "larsim/LArG4/ParticleListAction.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// CLHEP
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"

namespace larg4 {
  
  
  //---------------------------------------------------------------------------------------
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout(std::string const& name)
    : G4VSensitiveDetector(name)
  {
    // Initialize values for the electron-cluster calculation.
    ClearSimChannels();

    const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock = ts->TPCClock();

    // the standard name contains cryostat and TPC;
    // if we don't find it, we will detect the TPC at each Geant hit
    unsigned int cryostat, tpc;
    if (std::sscanf(name.c_str(),"%*19s%1u%*4s%u",&cryostat, &tpc) == 2)
      SetSingleTPC(cryostat, tpc);
    else
      SetDiscoverTPC();
    
  } // LArVoxelReadout::LArVoxelReadout()

  //---------------------------------------------------------------------------------------
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout
    (std::string const& name, unsigned int cryostat, unsigned int tpc)
    : LArVoxelReadout(name)
    { SetSingleTPC(cryostat, tpc); }

  
  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::Setup(Setup_t const& setupData) {
    
    SetOffPlaneChargeRecoveryMargin(setupData.offPlaneMargin);
    SetRandomEngines(setupData.propGen);
    
  } // LArVoxelReadout::Setup()
  
  
  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::SetSingleTPC(unsigned int cryostat, unsigned int tpc) {
    bSingleTPC = true;
    fCstat = cryostat;
    fTPC = tpc;
    MF_LOG_DEBUG("LArVoxelReadout")
      << GetName() << "covers C=" << fCstat << " T=" << fTPC;
  } // LArVoxelReadout::SetSingleTPC()

  void LArVoxelReadout::SetDiscoverTPC() {
    bSingleTPC = false;
    fCstat = 0;
    fTPC = 0;
    MF_LOG_DEBUG("LArVoxelReadout") << GetName() << " autodetects TPC";
  } // LArVoxelReadout::SetDiscoverTPC()
  
  
  //---------------------------------------------------------------------------------------
  LArVoxelReadout::~LArVoxelReadout() {}

  //---------------------------------------------------------------------------------------
  // Called at the start of each event.
  void LArVoxelReadout::Initialize(G4HCofThisEvent*)
  {
    // for c2: larp is unused
    //auto const * larp = lar::providerFrom<detinfo::LArPropertiesService>();
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fElectronLifetime      = detprop->ElectronLifetime();
    for (int i = 0; i<3; ++i)
      fDriftVelocity[i]    = detprop->DriftVelocity(detprop->Efield(i),
						    detprop->Temperature())/1000.;
    
    fElectronClusterSize   = fLgpHandle->ElectronClusterSize();
    fMinNumberOfElCluster  = fLgpHandle->MinNumberOfElCluster();
    fLongitudinalDiffusion = fLgpHandle->LongitudinalDiffusion();
    fTransverseDiffusion   = fLgpHandle->TransverseDiffusion();
    fDontDriftThem         = fLgpHandle->DisableWireplanes();
    fSkipWireSignalInTPCs  = fLgpHandle->SkipWireSignalInTPCs();

    MF_LOG_DEBUG("LArVoxelReadout")  << " e lifetime: "        << fElectronLifetime
                                  << "\n Temperature: "     << detprop->Temperature()
                                  << "\n Drift velocity: "  << fDriftVelocity[0]
                                  <<" "<<fDriftVelocity[1]<<" "<<fDriftVelocity[2];

    fDontDriftThem = (fDontDriftThem || fLgpHandle->NoElectronPropagation() );

    fNSteps=0;
  }

  //---------------------------------------------------------------------------------------
  // Called at the end of each event.
  void LArVoxelReadout::EndOfEvent(G4HCofThisEvent*)
  {
    mf::LogWarning("LArVoxelReadout") << "Total number of steps was " << fNSteps << std::endl;

  } // LArVoxelReadout::EndOfEvent()

  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::clear()
  {
  }

  //--------------------------------------------------------------------------------------
  void LArVoxelReadout::ClearSimChannels() {
    fChannelMaps.resize(fGeoHandle->Ncryostats());
    size_t cryo = 0;
    for (auto& cryoData: fChannelMaps) { // each, a vector of maps
      cryoData.resize(fGeoHandle->NTPC(cryo++));
      for (auto& channelsMap: cryoData) channelsMap.clear(); // each, a map
    } // for cryostats
  } // LArVoxelReadout::ClearSimChannels()

  
  const LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap() const
  {
    if (bSingleTPC) return GetSimChannelMap(fCstat, fTPC);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannelMap() const
  
  LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap() {
    if (bSingleTPC) return GetSimChannelMap(fCstat, fTPC);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannelMap()
  
  
  const LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap
    (unsigned short cryo, unsigned short tpc) const
    { return fChannelMaps.at(cryo).at(tpc); }
  
  LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap
    (unsigned short cryo, unsigned short tpc)
    { return fChannelMaps.at(cryo).at(tpc); }
  
  
  std::vector<sim::SimChannel> LArVoxelReadout::GetSimChannels() const {
    if (bSingleTPC) return GetSimChannels(fCstat, fTPC);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannels()
  
  std::vector<sim::SimChannel> LArVoxelReadout::GetSimChannels
    (unsigned short cryo, unsigned short tpc) const
  { 
    std::vector<sim::SimChannel> channels;
    const ChannelMap_t& chmap = fChannelMaps.at(cryo).at(tpc);
    channels.reserve(chmap.size());
    for(const auto& chpair: chmap) channels.push_back(chpair.second);
    return channels;
  } // LArVoxelReadout::GetSimChannels(short, short)
  
  
  //---------------------------------------------------------------------------------------
  // Called for each step.
  G4bool LArVoxelReadout::ProcessHits( G4Step* step, G4TouchableHistory* pHistory)
  {
    // All work done for the "parallel world" "box of voxels" in
    // LArVoxelReadoutGeometry makes this a fairly simple routine.
    // First, the usual check for non-zero energy:

    // Only process the hit if the step is inside the active volume and
    // it has deposited energy.  The hit being inside the active volume
    // is virtually sure to happen because the LArVoxelReadoutGeometry
    // that this class makes use of only has voxels for inside the TPC.

    // The step can be no bigger than the size of the voxel,
    // because of the geometry set up in LArVoxelGeometry and the
    // transportation set up in PhysicsList.  Find the mid-point
    // of the step.

    if ( step->GetTotalEnergyDeposit() > 0 ){
      
      // Make sure we have the IonizationAndScintillation singleton
      // reset to this step
      larg4::IonizationAndScintillation::Instance()->Reset(step);
      fNSteps++;
      if( !fDontDriftThem ){

        G4ThreeVector midPoint = 0.5*( step->GetPreStepPoint()->GetPosition()
                                       + step->GetPostStepPoint()->GetPosition() );
        double g4time = step->GetPreStepPoint()->GetGlobalTime();
        
        // Find the Geant4 track ID for the particle responsible for depositing the
        // energy.  if we are only storing primary EM shower particles, and this energy
        // is from a secondary etc EM shower particle, the ID returned is the primary
        const int trackID = ParticleListAction::GetCurrentTrackID();
        
        // Find out which TPC we are in.
        // If this readout object covers just one, we already know it.
        // Otherwise, we have to ask Geant where we are.
        unsigned short int cryostat = 0, tpc = 0;
        if (bSingleTPC) {
          cryostat = fCstat;
          tpc = fTPC;
        }
        else {
          // detect the TPC we are in
          const G4VTouchable* pTouchable = step->GetPreStepPoint()->GetTouchable();
          if (!pTouchable) {
            throw cet::exception
              ("LArG4") << "Untouchable step in LArVoxelReadout::ProcessHits()";
          }
          
          // one of the ancestors of the touched volume is supposed to be
          // actually a G4PVPlacementInTPC that knows which TPC it covers;
          // currently, it's the 4th in the ladder:
          // [0] voxel [1] voxel tower [2] voxel plane [3] the full box;
          G4int depth = 0;
          while (depth < pTouchable->GetHistoryDepth()) {
            const G4PVPlacementInTPC* pPVinTPC = 
              dynamic_cast<const G4PVPlacementInTPC*>
              (pTouchable->GetVolume(depth++));
            if (!pPVinTPC) continue;
            cryostat = pPVinTPC->ID.Cryostat;
            tpc = pPVinTPC->ID.TPC;
	    if (Has(fSkipWireSignalInTPCs, tpc))
	    {
	    	return true;
	    }
            break;
          } // while
          if (depth < pTouchable->GetHistoryDepth()) {
            // this is a fundamental error where the step does not happen in
            // any TPC; this should not happen in the readout geometry!
            throw cet::exception
              ("LArG4") << "No TPC ID found in LArVoxelReadout::ProcessHits()";
          } // if
          MF_LOG_DEBUG("LArVoxelReadoutHit") << " hit in C=" << cryostat << " T=" << tpc;
        } // if more than one TPC
        
        // Note that if there is no particle ID for this energy deposit, the
        // trackID will be sim::NoParticleId.
        
        DriftIonizationElectrons(midPoint, g4time, trackID, cryostat, tpc);
      } // end we are drifting
    } // end there is non-zero energy deposition

    return true;
  }
  
  
  //----------------------------------------------------------------------------
  void LArVoxelReadout::SetRandomEngines(CLHEP::HepRandomEngine* pPropGen) {
    assert(pPropGen); // random engine must be present
    fPropGen = pPropGen;
  } // LArVoxelReadout::SetRandomEngines()
  
  
  //----------------------------------------------------------------------------
  geo::Point_t LArVoxelReadout::RecoverOffPlaneDeposit
    (geo::Point_t const& pos, geo::PlaneGeo const& plane) const
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
    
    return plane.ComposePoint<geo::Point_t>(distance, landingPos + offPlane);
    
  } // LArVoxelReadout::RecoverOffPlaneDeposit()
  
  
  //----------------------------------------------------------------------------
  // energy is passed in with units of MeV, dx has units of cm
  void LArVoxelReadout::DriftIonizationElectrons(G4ThreeVector stepMidPoint,
                                                 const double simTime,
                                                 int trackID,
                                                 unsigned short int cryostat, unsigned short int tpc)
  {
    auto const * ts = lar::providerFrom<detinfo::DetectorClocksService>();

    // this must be always true, unless caller has been sloppy
    assert(fPropGen); // No propagation random generator provided?!
    
    CLHEP::RandGauss PropRand(*fPropGen);

    // This routine gets called frequently, once per every particle
    // traveling through every voxel. Use whatever tricks we can to
    // increase its execution speed.

    static double LifetimeCorr_const = -1000. * fElectronLifetime;
    static double LDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    static double TDiff_const        = std::sqrt(2.*fTransverseDiffusion);
    static double fRecipDriftVel[3]   = {1./fDriftVelocity[0], 
                                        1./fDriftVelocity[1], 
                                        1./fDriftVelocity[2]};

    struct Deposit_t {
      double energy = 0.;
      double electrons = 0.;
      
      void add(double more_energy, double more_electrons)
        { energy += more_energy; electrons += more_electrons; }
    }; // Deposit_t
    
    // Map of electrons to store - catalogued by map[channel][tdc]
    std::map<raw::ChannelID_t, std::map<unsigned int, Deposit_t>> DepositsToStore;

    double fDriftClusterPos[3] = {0.};

    double const xyz[3] = {stepMidPoint.x() / CLHEP::cm,
                           stepMidPoint.y() / CLHEP::cm,
                           stepMidPoint.z() / CLHEP::cm};

    // Already know which TPC we're in because we have been told

    try{
      const geo::TPCGeo& tpcGeo = fGeoHandle->TPC(tpc, cryostat);

      // X drift distance - the drift direction can be either in
      // the positive or negative direction, so use std::abs


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
      if(driftsign == 1 && tpcGeo.PlaneLocation(0)[driftcoordinate] < xyz[driftcoordinate] )
	return;
      if(driftsign == -1 && tpcGeo.PlaneLocation(0)[driftcoordinate] > xyz[driftcoordinate] )
	return;

      /// \todo think about effects of drift between planes.
      // Center of plane is also returned in cm units
      double DriftDistance = std::abs(xyz[driftcoordinate] - tpcGeo.PlaneLocation(0)[driftcoordinate]);

      // Get SCE {x,y,z} offsets for particular location in TPC      
      geo::Vector_t posOffsets;
      double posOffsetxyz[3] = {0.0,0.0,0.0}; //need this array for the driftcoordinate and transversecoordinates
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      if (SCE->EnableSimSpatialSCE() == true)
      {
        posOffsets = SCE->GetPosOffsets({ xyz[0], xyz[1], xyz[2] });
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
          
      const double lifetimecorrection = TMath::Exp(TDrift / LifetimeCorr_const);
      const int    nIonizedElectrons  = larg4::IonizationAndScintillation::Instance()->NumberIonizationElectrons();
      const double energy             = larg4::IonizationAndScintillation::Instance()->EnergyDeposit();
      
      // if we have no electrons (too small energy or too large recombination)
      // we are done already here
      if (nIonizedElectrons <= 0) {
        MF_LOG_DEBUG("LArVoxelReadout")
          << "No electrons drifted to readout, " << energy << " MeV lost.";
        return;
      }
      // includes the effect of lifetime
      const double nElectrons   = nIonizedElectrons * lifetimecorrection;

      // Longitudinal & transverse diffusion sigma (cm)
      double SqrtT    = std::sqrt(TDrift);
      double LDiffSig = SqrtT * LDiff_const;
      double TDiffSig = SqrtT * TDiff_const;
      double electronclsize = fElectronClusterSize;
      
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
      
      // Compute arrays of values as quickly as possible.
      std::vector< double > fLongDiff(nClus);
      std::vector< double > fTransDiff1(nClus);
      std::vector< double > fTransDiff2(nClus);
      std::vector< double > fnElDiff(nClus, electronclsize);
      std::vector< double > fnEnDiff(nClus);

      // fix the number of electrons in the last cluster, that has smaller size
      fnElDiff.back() = nElectrons - (nClus-1)*electronclsize;
      
      for(size_t xx = 0; xx < fnElDiff.size(); ++xx){
        if(nElectrons > 0) fnEnDiff[xx] = energy/nElectrons*fnElDiff[xx];
        else               fnEnDiff[xx] = 0.;
      }
      
      double const avegageYtransversePos
        = (stepMidPoint.y()/CLHEP::cm) + posOffsets.Y();
      double const avegageZtransversePos
        = (stepMidPoint.z()/CLHEP::cm) + posOffsets.Z();
      
      // Smear drift times by longitudinal diffusion
      if (LDiffSig > 0.0)
	PropRand.fireArray( nClus, &fLongDiff[0], 0., LDiffSig);
      else
	fLongDiff.assign(nClus, 0.0);
      
      if (TDiffSig > 0.0) {
	// Smear the coordinates in plane perpendicular to drift direction by the transverse diffusion
	PropRand.fireArray( nClus, &fTransDiff1[0], avegagetransversePos1, TDiffSig);
	PropRand.fireArray( nClus, &fTransDiff2[0], avegagetransversePos2, TDiffSig);
      }
      else {
	fTransDiff1.assign(nClus, avegagetransversePos1);
	fTransDiff2.assign(nClus, avegagetransversePos2);
      }

      // make a collection of electrons for each plane
      for(size_t p = 0; p < tpcGeo.Nplanes(); ++p){
        
        geo::PlaneGeo const& plane = tpcGeo.Plane(p);

	fDriftClusterPos[driftcoordinate] = tpcGeo.PlaneLocation(p)[driftcoordinate];
        
        // Drift nClus electron clusters to the induction plane
        for(int k = 0; k < nClus; ++k){

          // Correct drift time for longitudinal diffusion and plane
	  double TDiff = TDrift + fLongDiff[k] * fRecipDriftVel[0];

	  // Take into account different Efields between planes
	  // Also take into account special case for ArgoNeuT (Nplanes = 2 and drift direction = x): plane 0 is the second wire plane
	  for (size_t ip = 0; ip<p; ++ip){
	    TDiff += (tpcGeo.PlaneLocation(ip+1)[driftcoordinate] - tpcGeo.PlaneLocation(ip)[driftcoordinate]) * fRecipDriftVel[(tpcGeo.Nplanes() == 2 && driftcoordinate == 0)?ip+2:ip+1];
	  }

	  fDriftClusterPos[transversecoordinate1] = fTransDiff1[k];
	  fDriftClusterPos[transversecoordinate2] = fTransDiff2[k];  
          
          /// \todo think about effects of drift between planes
          
          // grab the nearest channel to the xyz1 position
          try{
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
            uint32_t channel = fGeoHandle->NearestChannel(fDriftClusterPos, p, tpc, cryostat);
            
            /// \todo check on what happens if we allow the tdc value to be
            /// \todo beyond the end of the expected number of ticks
            // Add potential decay/capture/etc delay effect, simTime.
            unsigned int tdc = fClock.Ticks(ts->G4ToElecTime(TDiff + simTime));

            // Add electrons produced by each cluster to the map
            DepositsToStore[channel][tdc].add(fnEnDiff[k], fnElDiff[k]);
          }
          catch(cet::exception &e){
            mf::LogWarning("LArVoxelReadout") << "unable to drift electrons from point ("
                                              << xyz[0] << "," << xyz[1] << "," << xyz[2]
                                              << ") with exception " << e;
          }
        } // end loop over clusters
      } // end loop over planes
      
      // Now store them in SimChannels
      ChannelMap_t& ChannelDataMap = fChannelMaps[cryostat][tpc];
      
      // browse deposited data on each channel: (channel; deposit data in time)
      for(auto const& deposit_per_channel: DepositsToStore){
        
        raw::ChannelID_t channel = deposit_per_channel.first;
        
        // find whether we already have this channel
        auto iChannelData = ChannelDataMap.find(channel);
        
        // channelData is the SimChannel these deposits are going to be added to
        // If there is such a channel already, use it (first beanch).
        // If it's a new channel, the inner assignment creates a new SimChannel
        // in the map, and we save its reference in channelData
        sim::SimChannel& channelData
          = (iChannelData == ChannelDataMap.end())
          ? (ChannelDataMap[channel] = sim::SimChannel(channel))
          : iChannelData->second;
        
        // go through all deposits, one for each TDC: (TDC, deposit data)
        for(auto const& deposit_per_tdc: deposit_per_channel.second) {
          channelData.AddIonizationElectrons(trackID,
                                             deposit_per_tdc.first,
                                             deposit_per_tdc.second.electrons,
                                             xyz,
                                             deposit_per_tdc.second.energy);
          
        } // for deposit on TDCs
      } // for deposit on channels

    } // end try intended to catch points where TPC can't be found
    catch(cet::exception &e){
      mf::LogWarning("LArVoxelReadout") << "step cannot be found in a TPC\n"
                                        << e;
    }

    return;
  }

  //---------------------------------------------------------------------------------------
  // Never used but still have to be defined for G4
  void LArVoxelReadout::DrawAll()  {}
  void LArVoxelReadout::PrintAll() {}

} // namespace larg4
