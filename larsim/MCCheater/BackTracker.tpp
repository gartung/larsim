namespace cheat{

  //--------------------------------------------------------------------
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The BackTrackerService is designed to impliment these methods as cleanly as possible within the art framework. This is intended for gallery users.
    void BackTracker::PrepEvent        (const Evt& evt ){
      if( !( this->CanRun( evt ) ) ){       
        throw cet::exception("ParticleInventory")
          << "Particle Inventory cannot function. "
          << "Is this file real data?";
      }
      fSimChannels.clear();
      fAllHitList.clear();
      this->PrepSimChannels( evt );
      //this->PrepAllHitList ( evt ); //This line temporarily commented out until I figure out how I want PrepAllHitList to work.

    }

  //--------------------------------------------------------------------
  template<typename Evt>
    void BackTracker::PrepSimChannels (const Evt& evt){
      if(this->SimChannelsReady()){ return;}
      //The SimChannels list needs to be built.
      const auto& simChannelIn = *(evt.template getValidHandle<art::Ptr<sim::SimChannel>>(fG4ModuleLabel));
      for( const auto& simChan : simChannelIn){
        fSimChannels.push_back(simChan);
      }
    }

  //--------------------------------------------------------------------
  template<typename Evt>
    void BackTracker::PrepAllHitList( const Evt& evt){
      if(this->AllHitListReady()){return;}
      const auto& allHitsIn = *(evt.template getValidHandle<art::Ptr<recob::Hit>>(fHitLabel));
      for (const auto& hit : allHitsIn ){
        fAllHitList.push_back(hit);
      }
    }

}//end namespace
