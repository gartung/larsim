////////////////////////////////////////////////////////////////////////////
//
// \file BackTracker.h
// \brief Functions needed by the BackTracker service in order to connect truth information with reconstruction.
//
// \author jason.stock@mines.sdsmt.edu
//
// Based on the original BackTracker by brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_BACKTRACKER_H
#define CHEAT_BACKTRACKER_H

//Includes
#include <vector>

#include "ParticleInventory.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

/*namespace recob{
    class SpacePoint;
}*/


namespace cheat{
  class BackTracker{
    public:
      BackTracker(const cheat::ParticleInventory* partInv,
                  const geo::GeometryCore* geom,
                  std::string g4label, double minHitFrac);
      BackTracker(const cheat::ParticleInventory* partInv,
                  const geo::GeometryCore* geom);
      ~BackTracker();

      template<typename Evt>
        void PrepEvent ( const Evt& evt );

      template<typename Evt>
        void PrepSimChannels ( const Evt& evt );

      template<typename Evt>
        void PrepRecoHitList ( const Evt& evt);

      template<typename Evt>
        void CheckCanRun( Evt& evt );

      void ClearEvent();

      const std::vector<const sim::SimChannel*>& SimChannels() const { return fSimChannels; }

      const std::vector<const sim::IDE* >      TrackIdToSimIDEs_Ps(int const& id) const;
      //const std::vector<const sim::IDE>       TrackIdToSimIDEs (int const& id) const; //I don't like this. It is poor use of memory. I may not include it in the release.
      const std::vector<const sim::IDE* >      TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const; //This function is removed as it requires the geometry service, not the geometry service provider.
      //std::vector<const sim::IDE>       TrackIdToSimIDEs (int const& id, const geo::View_t view) const; //I don't like this. It is poor use of memory. I may not include it in the release.

      const sim::SimChannel* FindSimChannel( raw::ChannelID_t channel ) const;

      const std::vector< sim::TrackIDE > ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const;

      //Track IDEs cannot be returned as pointers, as they dont exist in the data product, and we will not be storing them.
      const std::vector< sim::TrackIDE> HitToTrackIDE(recob::Hit const& hit) const;
//      std::vector< const sim::TrackIDE> HitToTrackIDE(art::Ptr<recob::Hit> const& hit) { return this->HitToTrackIDE(*hit); }

      const std::vector< int > HitToTrackId(recob::Hit const& hit) const ;
//      std::vector< const int> HitToTrackId(art::Ptr<recob::Hit> const& hit) { return this->HitToTrackId(*hit); }


      //I will not return these by copy,  as that could get very large very quickly.
      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const; 
      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId ) const;
      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds ) const; //planned
      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const;

      std::vector< const sim::IDE* const > HitToSimIDEsPs (recob::Hit const& hit) const;
//      std::vector< const sim::IDE* const > HitToSimIDEsPs (art::Ptr< recob::Hit > const& hit) const { return this->HitToSimIDEsPs (*hit); }
      std::vector< const sim::IDE > HitToSimIDEs (recob::Hit const& hit);
//      std::vector< const sim::IDE > HitToSimIDEs (art::Ptr< recob::Hit > const& hit) { return this->HitToSimIDEsPs (*hit); }

      std::vector<double> SimIDEsToXYZ( std::vector< sim::IDE > const& ides) const;

      std::vector<double> HitToXYZ(art::Ptr<recob::Hit> const& hit) const;

      std::vector< art::Ptr< recob::Hit > > SpacePointToHits(art::Ptr<recob::SpacePoint> const& spt) const;

      std::vector< double > SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt ) const;

      double HitCollectionPurity( std::set<int> const& trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits);

//      double HitCollectionEfficiency( std::set<int> const& trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits,
//          std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view); //This function removed as it depends on view, which causes issues with the geom service provider

      double HitChargeCollectionPurity( std::set<int> const& trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits);

//      double HitChargeCollectionEfficiency( std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits,
//          std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view); //This function removed as it depends on view, which causes issues with the geom service provider

      std::set<int> GetSetOfTrackIds();
      std::set<int> GetSetOfEveIDs();

      std::set<int> GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits );
      std::set<int> GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits );


    private:
      cheat::ParticleInventory const* fPartInv; //The constructor needs to put something in here
      geo::GeometryCore const* fGeom;
      std::string fG4ModuleLabel;
      double      fMinHitEnergyFraction;

      bool fCanRun=0;
      std::vector<const sim::SimChannel*>             fSimChannels;
      std::vector< art::Ptr<recob::Hit> >             fAllHits;

  };//end class BackTracker

}//end namespace cheat

#include "BackTracker.tpp"

#endif //CHEAT_BACKTRACKER_H
