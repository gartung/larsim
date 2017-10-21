/////////////////////////////////////////////////////////////////////////////////
// \file BackTrackerService.h
// \A service for backtracking reconstruction information to its truth information
//
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov
//
/////////////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_BACKTRACKERSERVICESERVICE_H
#define CHEAT_BACKTRACKERSERVICESERVICE_H

#include "BackTracker.h"

namespace cheat{
  class BackTrackerService
  {
    public:

      struct fhiclConfig{
        fhicl::Table < BackTracker::fhiclConfig > BackTrackerTable{
          fhicl::Name("providerBackTracker"),
            fhicl::Comment("This is the fhicl configuration of the BackTracker service provider.") };
      };


      using provider_type = BackTracker;
      const BackTracker* provider() const {return fBackTracker;}

      BackTrackerService(fhicl::ParameterSet const& pSet, art::ActivityRegistry& reg);
      BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg);
      ~BackTrackerService();

      //Move this function into the .cpp file and give it the appropriate prep and check functions
      const std::vector < art::Ptr < sim::SimChannel > >& const { return fBackTracker.SimChannels(); }

      const std::vector < const sim::IDE* > TrackIdToSimIDEs_Ps(int const& id) const;
      const std::vector < const sim::IDE* > TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const;

      art::Ptr < sim::SimChannel > FindSimChannel( raw::ChannelID_t channel ) const;

      const std::vector < sim::TrackIDE > ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const;

      const std::vector < sim::TrackIDE > HitToTrackIDEs(recob::Hit const& hit) const;
      const std::vector < sim::TrackIDE > HitToTrackIDEs(art::Ptr < recob::Hit > const& hit) const;

      const std::vector < int > HitToTrackIds(recob::Hit const& hit) const ;

      const std::vector < sim::TrackIDE > HitToEveTrackIDEs(recob::Hit const& hit) const;
      const std::vector < sim::TrackIDE > HitToEveTrackIDEs(art::Ptr < recob::Hit > const& hit) const;

      const std::vector < art::Ptr < recob::Hit > > TrackIdToHits_Ps( const int& tkId, std::vector < art::Ptr < recob::Hit > > const& hitsIn ) const;
      const std::vector < art::Ptr < recob::Hit > > TrackIdToHits_Ps( const int& tkId ) const;

      const std::vector < std::vector < art::Ptr < recob::Hit > > > TrackIdsToHits_Ps( std::vector < int > const& tkIds, std:: vector < art::Ptr < recob::Hit > > const& hitsIn ) const;
      const std::vector < std::vector < art::Ptr < recob::Hit > > > TrackIdsToHits_Ps( std::vector < int > const& tkIds ) const;

      const std::vector< sim::IDE > HitToAvgSimIDEs ( recob::Hit const& hit) const;
      const std::vector< sim::IDE > HitToAvgSimIDEs ( art::Ptr<recob::Hit> hit) const;

      const std::vector< const sim::IDE* > HitToSimIDEs_Ps (recob::Hit const& hit) const;
      const std::vector< const sim::IDE* > HitToSimIDEs_Ps (art::Ptr< recob::Hit > const& hit) const;

      const std::vector<double> SimIDEsToXYZ( std::vector< sim::IDE > const& ides) const;
      const std::vector<double> SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) const;

      const std::vector<double> HitToXYZ(const recob::Hit& hit) const;
      const std::vector<double> HitToXYZ(art::Ptr<recob::Hit> const& hit) const;

      const double HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits) const;
      const double HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const&     hits) const;

      const double HitCollectionEfficiency( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits, std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const ;

      const double HitChargeCollectionEfficiency( std::set<int> trackIds, std::vector< art::Ptr<recob::Hit> > const&        hits,        std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view) const;

      const std::set<int> GetSetOfTrackIds() const;
      const std::set<int> GetSetOfEveIds() const;

      const std::set<int> GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const;
      const std::set<int> GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits ) const;






    private:
      //FHICL parameters.
      //The BackTracker service has no parameters.

      cheat::BackTracker fBackTracker;

      const art::Event* fEvt;

      //Prep functions go here.
      void priv_PrepEvent ();
      void priv_PrepSimChannels ();
      void priv_PrepAllHitList ();

      bool priv_CanRun (const art::Event& evt);

      bool priv_SimChannelsReady() { return fBackTracker.SimChannelsReady();}
      bool priv_AllHitListReady() { return fBackTracker.AllHitListReady();}




  }; //class BackTrackerService

}


#endif //CHEAT_BACKTRACKERSERVICESERVICE_H
