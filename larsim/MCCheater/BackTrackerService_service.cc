////////////////////////////////////////////////////////////////////////////////////////
// 
// \file BackTrackerService_service.cc
// \brief A service for backtracking reconstruction information to its truth information
// 
// \author jason.stock@mines.sdsmt.edu
// Based on the original BackTracker by Brian Rebel (brebel@fnal.gov)
//
////////////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "BackTrackerService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cheat{
//  lar::providerFrom<detinfo::DetectorClocksService>();
//  lar::providerFrom<geo::Geometry>();
  //---------------------------------------------------------------------
    BackTrackerService::BackTrackerService( const fhicl::ParameterSet& pSet,  art::ActivityRegistry& reg)
      :fBackTracker(
          pSet.get<fhicl::ParameterSet>("providerBackTracker"), 
          lar::providerFrom<cheat::ParticleInventory>(), 
          lar::providerFrom<geo::Geometry>(),
          lar::providerFrom<detinfo::DetectorClocksService>())
    {
          reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
          //setup all services here.
          //Needed services include geo::Geometry, ParticleInventory, and DetectorClocks.
//          fBackTracker(pSet.get<fhicl::ParameterSet>("providerBackTracker")/*add pointers to other service providers*/)
    }

  //---------------------------------------------------------------------
    BackTrackerService::BackTrackerService(const fhiclConfig& config, art::ActivityRegistry& reg)
     :fBackTracker(
         config.BackTrackerTable(),
         lar::providerFrom<cheat::ParticleInventory>(),
         lar::providerFrom<geo::Geometry>(),
         lar::providerFrom<detinfo::DetectorClocksService>())
    {
          reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
    }

  //---------------------------------------------------------------------
    BackTrackerService::~BackTrackerService(){}

  void BackTrackerService::priv_PrepEvent( const art::Event& evt ){
    fEvt=&evt;
    //fBackTracker.ClearEvent();
  }




}//end namespace cheat
