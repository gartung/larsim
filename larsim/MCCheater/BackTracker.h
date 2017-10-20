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

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

/*namespace recob{
  class SpacePoint;
  }*/


namespace cheat{

  class BackTracker{
    public:
      struct fhiclConfig{
        fhicl::Atom<art::InputTag> G4ModuleLabel{fhicl::Name("G4ModuleLabel"), fhicl::Comment("The label of the LArG4   module used to produce the art file we will be using."), "largeant"};
        fhicl::Atom<art::InputTag> DefaultHitModuleLabel{fhicl::Name("DefaultHitModuleLabel"), fhicl::Comment("The label  of the module used to produce the hits in the art file we will default to when no hitlist is provided."), "hitfd"};
        fhicl::Atom<double> MinHitEnergyFraction{fhicl::Name("MinHitEnergyFraction"), fhicl::Comment("The minimum     contribution an energy deposit must make to a Hit to be considered part of that hit."),0.010};
      };
      BackTracker(const fhiclConfig& config, const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock);
      BackTracker(const fhicl::ParameterSet& pSet, const cheat::ParticleInventory* partInv,
          const geo::GeometryCore* geom, const detinfo::DetectorClocks* detClock);
      ~BackTracker();

      template<typename Evt>
        void PrepEvent ( const Evt& evt );

      template<typename Evt>
        void PrepSimChannels ( const Evt& evt );

      template<typename Evt>
        void PrepAllHitList ( const Evt& evt);

      template<typename Evt>
        void CheckCanRun( Evt& evt );
      //These two functions need to be redesigned as templates as they MUST accept the event to use FindManyP.
      /*
         std::vector< art::Ptr< recob::Hit > > SpacePointsToHits_Ps(const recob::SpacePoint& spt) const;

         std::vector< double > SpacePointToXYZ( art::Ptr< recob::SpacePoint > const& spt ) const;
         */

      void ClearEvent();

      const std::vector<const sim::SimChannel*>& SimChannels() const { return fSimChannels; }

      const std::vector<const sim::IDE* >   TrackIdToSimIDEs_Ps(int const& id) const;
      //const std::vector<const sim::IDE>    TrackIdToSimIDEs (int const& id) const; 
      const std::vector<const sim::IDE* >   TrackIdToSimIDEs_Ps(int const& id, const geo::View_t view) const; 
      //std::vector<const sim::IDE>    TrackIdToSimIDEs (int const& id, const geo::View_t view) const; 

      const sim::SimChannel* FindSimChannel( raw::ChannelID_t channel ) const;

      const std::vector< sim::TrackIDE > ChannelToTrackIDEs(raw::ChannelID_t channel, const double hit_start_time, const double hit_end_time) const;


      //Track IDEs cannot be returned as pointers, as they dont exist in the data product, and we will not be storing them.
      const std::vector< sim::TrackIDE> HitToTrackIDEs(recob::Hit const& hit) const;
      const std::vector< sim::TrackIDE> HitToTrackIDEs(art::Ptr<recob::Hit> const& hit) const { return this->HitToTrackIDEs(*hit);}

      const std::vector< int > HitToTrackIds(recob::Hit const& hit) const ;
      //   std::vector< const int> HitToTrackId(art::Ptr<recob::Hit> const& hit) { return this->HitToTrackId(*hit); }

      const std::vector<sim::TrackIDE> HitToEveTrackIDEs(recob::Hit const& hit) const;
      const std::vector<sim::TrackIDE> HitToEveTrackIDEs(art::Ptr<recob::Hit> const& hit) const{ return this->HitToEveTrackIDEs(*hit);}

      //I will not return these by copy, as that could get very large very quickly.
      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const; 
      std::vector< art::Ptr<recob::Hit> > TrackIdToHits_Ps( const int& tkId ) const
      {return this->TrackIdToHits_Ps(tkId, fAllHitList); }

      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds, std::vector< art::Ptr< recob::Hit > > const& hitsIn ) const;
      std::vector< std::vector< art::Ptr<recob::Hit> > > TrackIdsToHits_Ps( std::vector<int> const& tkIds ) const
      {return this->TrackIdsToHits_Ps(tkIds, fAllHitList);}

      const std::vector< sim::IDE > HitToAvgSimIDEs ( recob::Hit const& hit) const;
      const std::vector< sim::IDE > HitToAvgSimIDEs ( art::Ptr<recob::Hit> hit) const{ return this->HitToAvgSimIDEs(*hit);}

      const std::vector< const sim::IDE* > HitToSimIDEs_Ps (recob::Hit const& hit) const;
      const std::vector< const sim::IDE* > HitToSimIDEs_Ps (art::Ptr< recob::Hit > const& hit) const { return this->HitToSimIDEs_Ps (*hit); }

      //const std::vector< sim::IDE > HitToSimIDEs (recob::Hit const& hit);
      //   std::vector< const sim::IDE > HitToSimIDEs (art::Ptr< recob::Hit > const& hit) { return this->HitToSimIDEsPs (*hit); }

      std::vector<double> SimIDEsToXYZ( std::vector< sim::IDE > const& ides) const;
      std::vector<double> SimIDEsToXYZ( std::vector< const sim::IDE* > const& ide_Ps) const;

      std::vector<double> HitToXYZ(const recob::Hit& hit) const;
      std::vector<double> HitToXYZ(art::Ptr<recob::Hit> const& hit) const{ return this->HitToXYZ(*hit);}



      double HitCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits);
      double HitChargeCollectionPurity( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits);

      double HitCollectionEfficiency( std::set<int> const& trackIds, std::vector< art::Ptr<recob::Hit> > const& hits, std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view); 


      double HitChargeCollectionEfficiency( std::set<int> trackIds, std::vector< art::Ptr<recob::Hit> > const& hits,        std::vector< art::Ptr<recob::Hit> > const& allhits, geo::View_t const& view); //This function removed as it depends on view, which causes issues with the geom service provider

      std::set<int> GetSetOfTrackIds(){ return fPartInv->GetSetOfTrackIds();}
      std::set<int> GetSetOfEveIds(){ return fPartInv->GetSetOfEveIds();}

      std::set<int> GetSetOfTrackIds( std::vector< art::Ptr< recob::Hit > > const& hits );
      std::set<int> GetSetOfEveIds( std::vector< art::Ptr< recob::Hit > > const& hits );


    private:
      cheat::ParticleInventory const* fPartInv; //The constructor needs to put something in here
      geo::GeometryCore    const* fGeom;
      detinfo::DetectorClocks const* fDetClocks;
      const art::InputTag      fG4ModuleLabel;
      const art::InputTag      fHitLabel;
      const double          fMinHitEnergyFraction;


      std::vector<const sim::SimChannel*>       fSimChannels;
      std::vector< art::Ptr<recob::Hit> >       fAllHitList;

  };//end class BackTracker

}//end namespace cheat

#include "BackTracker.tpp"

#endif //CHEAT_BACKTRACKER_H
