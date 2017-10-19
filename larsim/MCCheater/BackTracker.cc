////////////////////////////////////////////////////////////////////
//
//
// \file BackTracker.cc
// \brief The functions needed for the BackTracker class needed by the BackTracker service in order to connect truth information with reconstruction.
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
//
///////////////////////////////////////////////////////////////////

//Includes
#include "BackTracker.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/SimListUtils.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"


namespace cheat{

  //-----------------------------------------------------------------------
  BackTracker::BackTracker(const fhiclConfig& config, const cheat::ParticleInventory* partInv, const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock )
    :fPartInv(partInv),fGeom(geom),fDetClocks(detClock),fG4ModuleLabel(config.G4ModuleLabel()),
     fHitLabel(config.DefaultHitModuleLabel()),fMinHitEnergyFraction(config.MinHitEnergyFraction())
  {
  }

  //-----------------------------------------------------------------------
  BackTracker::BackTracker(const fhicl::ParameterSet& pSet, const cheat::ParticleInventory* partInv, const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock)
    :fPartInv(partInv),fGeom(geom),fDetClocks(detClock),
     fG4ModuleLabel       (pSet.get<art::InputTag>("G4ModuleLabel", "largeant")),
     fHitLabel            (pSet.get<art::InputTag>("DefaultHitModuleLabel", "hitfd")),
     fMinHitEnergyFraction(pSet.get<double>       ("MinHitEnergyFraction", 0.010))
  {
  }

  //-----------------------------------------------------------------------
  BackTracker::~BackTracker()
  {
  }


  //-----------------------------------------------------------------------
  void BackTracker::ClearEvent(){
    fSimChannels.clear();
    //Do not clear the inventory here. This is something the service will do, or the use themselves.
  }


  //-----------------------------------------------------------------------
  const std::vector< const sim::IDE* > BackTracker::TrackIdToSimIDEs_Ps( int const& id ) const{
    std::vector< const sim::IDE* > ideps;
    for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
      const auto & tdcidemap = fSimChannels[sc]->TDCIDEMap(); //This returns a reference.
      // loop over the IDEMAP
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
        // loop over the vector of IDE objects.
        const std::vector<sim::IDE>& idevec = (*mapitr).second; //note, mapitr.second returns the actual data from the map, not a copy
        for(size_t iv = 0; iv < idevec.size(); ++iv){
          const sim::IDE* const idep = &idevec[iv];
          //if( abs(idevec[iv].trackID) == id) continue;
          //ideps.push_back(idep);
          if( abs(idevec[iv].trackID) == id) ideps.push_back(&(idevec[iv]));
        }//end for index in idevec
      } // end loop over map from sim::SimChannel
    } // end loop over sim::SimChannels
    return ideps;
  }


  //-----------------------------------------------------------------------
  const std::vector<const sim::IDE* >   BackTracker::TrackIdToSimIDEs_Ps (int const& id, const geo::View_t view) const
  {  
    std::vector<const sim::IDE*> ide_Ps;
    for(const sim::SimChannel* sc : fSimChannels){
      if (fGeom->View(sc->Channel()) != view) continue;

      // loop over the IDEMAP
      for(const auto & item : sc->TDCIDEMap()){

        // loop over the vector of IDE objects.
        for(const sim::IDE & ide : item.second){
          if (abs(ide.trackID) == id) ide_Ps.push_back(&ide);
        }
      } // end loop over map from sim::SimChannel
    } // end loop over sim::SimChannels

    return ide_Ps;
  }


  //-----------------------------------------------------------------------
  const sim::SimChannel* BackTracker::FindSimChannel(raw::ChannelID_t channel) const{
    const sim::SimChannel* chan = 0;
    auto ilb = std::lower_bound(fSimChannels.begin(),fSimChannels.end(),channel,[](const sim::SimChannel *a, raw::    ChannelID_t channel) {return(a->Channel()<channel);});
    if (ilb != fSimChannels.end())
      if ( (*ilb)->Channel() == channel) {chan = *ilb;}
    if(!chan)
      throw cet::exception("BackTracker") << "No sim::SimChannel corresponding "
        << "to channel: " << channel << "\n";
    return chan;
  }


  //-----------------------------------------------------------------------
  const std::vector< sim::TrackIDE > BackTracker::ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const{
    std::vector< sim::TrackIDE > trackIDEs;
    double totalE=0.;
    try{
      const sim::SimChannel* schannel = this->FindSimChannel(channel);

      // loop over the electrons in the channel and grab those that are in time
      // with the identified hit start and stop times
      int start_tdc = fDetClocks->TPCTick2TDC( hit_start_time );
      int end_tdc   = fDetClocks->TPCTick2TDC( hit_end_time   );
      if(start_tdc<0) start_tdc = 0;
      if(end_tdc<0) end_tdc = 0;
      std::vector<sim::IDE> simides = schannel->TrackIDsAndEnergies(start_tdc, end_tdc);

      // first get the total energy represented by all track ids for
      // this channel and range of tdc values
      for(size_t e = 0; e < simides.size(); ++e)
        totalE += simides[e].energy;


      // protect against a divide by zero below
      if(totalE < 1.e-5) totalE = 1.;

      // loop over the entries in the map and fill the input vectors

      for(size_t e = 0; e < simides.size(); ++e){

        if(simides[e].trackID == sim::NoParticleId) continue;

        sim::TrackIDE info;
        info.trackID    = simides[e].trackID;
        info.energyFrac = simides[e].energy/totalE;
        info.energy     = simides[e].energy;

        trackIDEs.push_back(info);

      }
    }// end try
    catch(cet::exception e){
      mf::LogWarning("BackTracker") << "caught exception \n"
        << e;
    }

    return trackIDEs;

  }


  //-----------------------------------------------------------------------
  const std::vector< sim::TrackIDE > BackTracker::HitToTrackIDEs( recob::Hit const& hit) const {
    std::vector<  sim::TrackIDE > trackIDEs;
    const double start = hit.PeakTimeMinusRMS();
    const double end   = hit.PeakTimePlusRMS();
    trackIDEs = this->ChannelToTrackIDEs(hit.Channel(), start, end);
    return trackIDEs;
  }


  //-----------------------------------------------------------------------
  const std::vector< int > BackTracker::HitToTrackIds(recob::Hit const& hit) const {
    std::vector< int > retVec;
    for(auto const trackIDE : this->HitToTrackIDEs( hit ) ){
      retVec.push_back( trackIDE.trackID );
    }
    return retVec;
  }

  //-----------------------------------------------------------------------
  std::vector < art::Ptr< recob::Hit > > BackTracker::TrackIdToHits_Ps( const int& tkId, std::vector<art::Ptr<recob::Hit>> const& hitsIn ) const{
    // returns a subset of the hits in the allhits collection that are matched
    // to the given track

    // temporary vector of TrackIDs and Ptrs to hits so only one
    // loop through the (possibly large) allhits collection is needed
    std::vector< art::Ptr<recob::Hit>> hitList;
    std::vector<sim::TrackIDE> trackIDE;
    for(auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      trackIDE.clear();
      art::Ptr<recob::Hit> const& hit = *itr;
      trackIDE = this->ChannelToTrackIDEs( hit->Channel(), hit->PeakTimeMinusRMS(), hit->PeakTimePlusRMS());
      for(auto itr_trakIDE = trackIDE.begin(); itr_trakIDE != trackIDE.end(); ++itr_trakIDE) {
        if(itr_trakIDE->trackID == tkId && itr_trakIDE->energyFrac > fMinHitEnergyFraction)
          hitList.push_back( hit);
      } // itr_trakIDE
    } // itr
    return hitList;
  }


  //-----------------------------------------------------------------------
  //This function could clearly be made by calling TrackIdToHits for each trackId, but that would be significantly slower because we would loop through all hits many times.
  std::vector< std::vector< art::Ptr<recob::Hit> > > BackTracker::TrackIdsToHits_Ps( std::vector<int> const& tkIds, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const{
    // returns a subset of the hits in the allhits collection that are matched
    // to MC particles listed in tkIds

    // temporary vector of TrackIDs and Ptrs to hits so only one
    // loop through the (possibly large) allhits collection is needed
    std::vector<std::pair<int, art::Ptr<recob::Hit>>> hitList;
    std::vector<sim::TrackIDE> tids;
    for(auto itr = hitsIn.begin(); itr != hitsIn.end(); ++itr) {
      tids.clear();
      art::Ptr<recob::Hit> const& hit = *itr;
      tids=this->ChannelToTrackIDEs( hit->Channel(), hit->PeakTimeMinusRMS(), hit->PeakTimePlusRMS());
      for(auto itid = tids.begin(); itid != tids.end(); ++itid) {
        for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
          if(itid->trackID == *itkid) {
            if(itid->energyFrac > fMinHitEnergyFraction)
              hitList.push_back(std::make_pair(*itkid, hit));
          }
        } // itkid
      } // itid
    } // itr

    // now build the truHits vector that will be returned to the caller
    std::vector<std::vector<art::Ptr<recob::Hit>>> truHits;
    // temporary vector containing hits assigned to one MC particle
    std::vector<art::Ptr<recob::Hit>> tmpHits;
    for(auto itkid = tkIds.begin(); itkid != tkIds.end(); ++itkid) {
      tmpHits.clear();
      for(auto itr = hitList.begin(); itr != hitList.end(); ++itr) {
        if(*itkid == (*itr).first) tmpHits.push_back((*itr).second);
      }
      truHits.push_back(tmpHits);
    }
    return truHits;
  }

  //-----------------------------------------------------------------------
  //Cannot be returned as a pointer, as these IDEs do not exist in the event. They are constructed on the fly.
  const std::vector<  sim::IDE > BackTracker::HitToAvgSimIDEs (recob::Hit const& hit) const{
    // Get services.

    int start_tdc = fDetClocks->TPCTick2TDC( hit.PeakTimeMinusRMS() );
    int end_tdc   = fDetClocks->TPCTick2TDC( hit.PeakTimePlusRMS()   );
    if(start_tdc<0) start_tdc = 0;
    if(end_tdc<0) end_tdc = 0;

    return (this->FindSimChannel(hit.Channel()))->TrackIDsAndEnergies(start_tdc, end_tdc);
  }

  //-----------------------------------------------------------------------
  const std::vector< const sim::IDE* > BackTracker::HitToSimIDEs_Ps (recob::Hit const& hit) const{
    std::vector< const sim::IDE* > retVec;
    const auto start_tdc = hit.PeakTimeMinusRMS();
    const auto end_tdc = hit.PeakTimePlusRMS();
    if(start_tdc > end_tdc){throw;}
    std::vector< std::pair<unsigned short, std::vector<sim::IDE>> > tdcIDEMap = (this->FindSimChannel(hit.Channel()))->TDCIDEMap(); //This in fact does not return a map. It returns a vector... with no guarantee that it is sorted...

    auto pairSort = [](auto& a, auto& b) { return a.first < b.first ; } ;
    if( !std::is_sorted( tdcIDEMap.begin(), tdcIDEMap.end(), pairSort)) {
      std::sort (tdcIDEMap.begin(), tdcIDEMap.end(), pairSort);
    }

    std::vector<sim::IDE> dummyVec; //I need something to stick in a pair to compare pair<tdcVal, IDE>. This is an otherwise useless "hack".
    std::pair<double, std::vector<sim::IDE>> start_tdcPair = std::make_pair(start_tdc,dummyVec); //This pair is a "hack" to make my comparison work for lower and upper bound.
    std::pair<double, std::vector<sim::IDE>> end_tdcPair = std::make_pair(end_tdc,dummyVec);
    std::vector<std::pair<unsigned short, std::vector<sim::IDE> > >::iterator  //iterator to the first interesting IDE
      mapFirst = std::lower_bound(tdcIDEMap.begin(), tdcIDEMap.end(), start_tdcPair, pairSort);
    std::vector<std::pair<unsigned short, std::vector<sim::IDE> > >::iterator  //iterator to just after the last interesting IDE
      mapLast  = std::upper_bound(tdcIDEMap.begin(), tdcIDEMap.end(), end_tdcPair, pairSort);
    for( auto mapitr = mapFirst; mapitr != mapLast; ++mapitr ){
      for( auto ide : mapitr->second){ retVec.push_back(&ide);} //Add all interesting IDEs to the retVec
    }
    return retVec;
  }

  //------------------------------------------------------------------------------
}//End namespace cheat
