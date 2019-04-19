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

// LArSoft code
#include "larsim/LArG4/LArVoxelReadout.h"
#include "larsim/LArG4/ParticleListAction.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/zip.h"
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
    , fGeom(*lar::providerFrom<geo::Geometry>())
    , fLgp(*art::ServiceHandle<sim::LArG4Parameters>())
    , fChannelMaps(fGeom.makeTPCData<ChannelMap_t>())
  {
    // Initialize values for the electron-cluster calculation.
    ClearSimChannels();

    const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock = ts->TPCClock();

    // the standard name contains cryostat and TPC;
    // if we don't find it, we will detect the TPC at each Geant hit
    unsigned int cryostat, tpc;
    if (std::sscanf(name.c_str(),"%*19s%1u%*4s%u",&cryostat, &tpc) == 2)
      SetSingleTPC({ cryostat, tpc });
    else
      SetDiscoverTPC();

  } // LArVoxelReadout::LArVoxelReadout()

  //---------------------------------------------------------------------------------------
  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout(std::string const& name, geo::TPCID const& tpcid)
    : LArVoxelReadout(name)
    { SetSingleTPC(tpcid); }


  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::Setup(Setup_t const& setupData) {

    SetOffPlaneChargeRecoveryMargin(setupData.offPlaneMargin);
    SetRandomEngines(setupData.propGen);

  } // LArVoxelReadout::Setup()


  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::SetSingleTPC(geo::TPCID const& tpcid) {
    fTPCID = tpcid;
    MF_LOG_DEBUG("LArVoxelReadout") << GetName() << "covers " << fTPCID;
  } // LArVoxelReadout::SetSingleTPC()

  void LArVoxelReadout::SetDiscoverTPC() {
    fTPCID = geo::TPCID{};
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

    fElectronClusterSize   = fLgp.ElectronClusterSize();
    fMinNumberOfElCluster  = fLgp.MinNumberOfElCluster();
    fLongitudinalDiffusion = fLgp.LongitudinalDiffusion();
    fTransverseDiffusion   = fLgp.TransverseDiffusion();
    fDontDriftThem         = fLgp.DisableWireplanes();
    fSkipWireSignalInTPCs  = fLgp.SkipWireSignalInTPCs();

    MF_LOG_DEBUG("LArVoxelReadout")  << " e lifetime: "        << fElectronLifetime
                                  << "\n Temperature: "     << detprop->Temperature()
                                  << "\n Drift velocity: "  << fDriftVelocity[0]
                                  <<" "<<fDriftVelocity[1]<<" "<<fDriftVelocity[2];

    fDontDriftThem = (fDontDriftThem || fLgp.NoElectronPropagation() );

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
    fChannelMaps.clear();
  } // LArVoxelReadout::ClearSimChannels()


  const LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap() const
  {
    if (isSingleTPC()) return GetSimChannelMap(fTPCID);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannelMap() const

  LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap() {
    if (isSingleTPC()) return GetSimChannelMap(fTPCID);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannelMap()


  LArVoxelReadout::ChannelMap_t const& LArVoxelReadout::GetSimChannelMap
    (geo::TPCID const& tpcid) const
    { return fChannelMaps.at(tpcid); }

  LArVoxelReadout::ChannelMap_t& LArVoxelReadout::GetSimChannelMap
    (geo::TPCID const& tpcid)
    { return fChannelMaps.at(tpcid); }


  std::vector<sim::SimChannel> LArVoxelReadout::GetSimChannels() const {
    if (isSingleTPC()) return GetSimChannels(fTPCID);
    throw cet::exception("LArVoxelReadout") << "TPC not specified";
  } // LArVoxelReadout::GetSimChannels()

  std::vector<sim::SimChannel> LArVoxelReadout::GetSimChannels
    (geo::TPCID const& tpcid) const
  {
    std::vector<sim::SimChannel> channels;
    const ChannelMap_t& chmap = fChannelMaps.at(tpcid);
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
        geo::TPCID tpcid;
        if (isSingleTPC()) {
          tpcid = fTPCID;
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
            tpcid = pPVinTPC->ID;
            if (Has(fSkipWireSignalInTPCs, tpcid.TPC))
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
          MF_LOG_DEBUG("LArVoxelReadoutHit") << " hit in " << tpcid;
        } // if more than one TPC

        // Note that if there is no particle ID for this energy deposit, the
        // trackID will be sim::NoParticleId.

        DriftIonizationElectrons(midPoint, g4time, trackID, tpcid);
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
                                                 geo::TPCID const& tpcid)
  {
    auto const * ts = lar::providerFrom<detinfo::DetectorClocksService>();
    auto const& SCE = *(lar::providerFrom<spacecharge::SpaceChargeService>());
    // this must be always true, unless caller has been sloppy
    assert(fPropGen); // No propagation random generator provided?!

    CLHEP::RandGauss PropRand(*fPropGen);

    // This routine gets called frequently, once per every particle
    // traveling through every voxel. Use whatever tricks we can to
    // increase its execution speed.

    static double LifetimeCorr_const = -1000. * fElectronLifetime;
    static double LDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    static double TDiff_const        = std::sqrt(2.*fTransverseDiffusion);
    static double RecipDriftVel[3]   = {1./fDriftVelocity[0],
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

    geo::Point_t const ionizationCenter = {
      stepMidPoint.x() / CLHEP::cm,
      stepMidPoint.y() / CLHEP::cm,
      stepMidPoint.z() / CLHEP::cm
    };
    
    // Already know which TPC we're in because we have been told

    try{
      const geo::TPCGeo &tpcg = fGeom.TPC(tpcid);
      
      geo::PlaneGeo const& firstPlane = tpcg.FirstPlane();

      double XDrift = firstPlane.DistanceFromPlane(ionizationCenter);
      if(XDrift < 0.) return;

      // Get SCE {x,y,z} offsets for particular location in TPC
      geo::Vector_t posOffsets;
      if (SCE.EnableSimSpatialSCE() == true)
      {
        posOffsets = SCE.GetPosOffsets(ionizationCenter);
      }
      posOffsets.SetX(-posOffsets.X());

      // Drift time (nano-sec)
      double TDrift;
      XDrift += posOffsets.X();

      // Space charge distortion could push the energy deposit beyond the wire
      // plane (see issue #15131). Given that we don't have any subtlety in the
      // simulation of this region, bringing the deposit exactly on the plane
      // should be enough for the time being.
      if (XDrift < 0.) XDrift = 0.;

      TDrift = XDrift * RecipDriftVel[0];
      if (tpcg.Nplanes() == 2){// special case for ArgoNeuT (plane 0 is the second wire plane)
        TDrift = ((XDrift - tpcg.PlanePitch(0,1)) * RecipDriftVel[0]
                  + tpcg.PlanePitch(0,1) * RecipDriftVel[1]);
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
      std::vector< double > XDiff(nClus);
      std::vector< double > YDiff(nClus);
      std::vector< double > ZDiff(nClus);
      std::vector< double > nElDiff(nClus, electronclsize);
      std::vector< double > nEnDiff(nClus);

      // fix the number of electrons in the last cluster, that has smaller size
      nElDiff.back() = nElectrons - (nClus-1)*electronclsize;

      for(size_t xx = 0; xx < nElDiff.size(); ++xx){
        if(nElectrons > 0) nEnDiff[xx] = energy/nElectrons*nElDiff[xx];
        else               nEnDiff[xx] = 0.;
      }

      double const averageYtransversePos
        = ionizationCenter.Y() + posOffsets.Y();
      double const averageZtransversePos
        = ionizationCenter.Z() + posOffsets.Z();

      // Smear drift times by x position and drift time
      if (LDiffSig > 0.0)
        PropRand.fireArray( nClus, &XDiff[0], 0., LDiffSig);
      else
        XDiff.assign(nClus, 0.0);

      if (TDiffSig > 0.0) {
        // Smear the Y,Z position by the transverse diffusion
        PropRand.fireArray( nClus, &YDiff[0], averageYtransversePos, TDiffSig);
        PropRand.fireArray( nClus, &ZDiff[0], averageZtransversePos, TDiffSig);
      }
      else {
        YDiff.assign(nClus, averageYtransversePos);
        ZDiff.assign(nClus, averageZtransversePos);
      }

      // make a collection of electrons for each plane
      for(size_t p = 0; p < tpcg.Nplanes(); ++p){
        
        geo::PlaneID const planeid(tpcid, p);
        
        geo::PlaneGeo const& plane = tpcg.Plane(planeid);

        double Plane0Pitch = tpcg.Plane0Pitch(p);

        // "-" sign is because Plane0Pitch output is positive. Andrzej
        double const xOnPlane = firstPlane.GetCenter().X() - Plane0Pitch;

        // Drift nClus electron clusters to the induction plane
        for(int k = 0; k < nClus; ++k){
          // Correct drift time for longitudinal diffusion and plane
          double TDiff = TDrift + XDiff[k] * RecipDriftVel[0];
          // Take into account different Efields between planes
          // Also take into account special case for ArgoNeuT where Nplanes = 2.
          for (size_t ip = 0; ip<p; ++ip){
            TDiff += tpcg.PlanePitch(ip,ip+1) * RecipDriftVel[tpcg.Nplanes()==3?ip+1:ip+2];
          }
          geo::Point_t clusterCenter { xOnPlane, YDiff[k], ZDiff[k] };

          /// \todo think about effects of drift between planes

          // grab the nearest channel to the clusterCenter position
          try{
            if (fOffPlaneMargin != 0) {
              // get the effective position where to consider the charge landed;
              //
              // Some optimisations are possible; in particular, this method
              // could be extended to inform us if the point was too far.
              // Currently, if that is the case the code will proceed, find the
              // point is off plane, emit a warning and skip the deposition.
              //
              clusterCenter = RecoverOffPlaneDeposit(clusterCenter, plane);

            } // if charge lands off plane
            raw::ChannelID_t const channel = fGeom.NearestChannel(clusterCenter, planeid);

            /// \todo check on what happens if we allow the tdc value to be
            ///       beyond the end of the expected number of ticks
            // Add potential decay/capture/etc delay effect, simTime.
            unsigned int const tdc = fClock.Ticks(ts->G4ToElecTime(TDiff + simTime));

            // Add electrons produced by each cluster to the map
            DepositsToStore[channel][tdc].add(nEnDiff[k], nElDiff[k]);
          }
          catch(cet::exception &e){
            mf::LogWarning("LArVoxelReadout") << "unable to drift electrons from point "
                                              << ionizationCenter << " with exception " << e;
          }
        } // end loop over clusters
      } // end loop over planes

      // Now store them in SimChannels
      ChannelMap_t& ChannelDataMap = fChannelMaps[tpcid];

      // browse deposited data on each channel: (channel; deposit data in time)
      for(auto const& deposit_per_channel: DepositsToStore){

        raw::ChannelID_t const channel = deposit_per_channel.first;

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
          std::array<double, 3U> xyz;
          geo::vect::fillCoords(xyz, ionizationCenter);
          channelData.AddIonizationElectrons(trackID,
                                             deposit_per_tdc.first,
                                             deposit_per_tdc.second.electrons,
                                             xyz.data(),
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
