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

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace cheat{

  //---------------------------------------------------------------------
  BackTrackerService::BackTrackerService( const fhicl::ParameterSet& pSet,  art::ActivityRegistry& reg)
    :fBackTracker(
        pSet.get<fhicl::ParameterSet>("providerBackTracker"), 
        lar::providerFrom<cheat::ParticleInventory>(), 
        lar::providerFrom<geo::Geometry>(),
        lar::providerFrom<detinfo::DetectorClocksService>())
  {
    reg.sPreProcessEvent.watch(this, &BackTrackerService::priv_PrepEvent);
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
