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

namespace cheat{
  BackTracker::BackTracker(std::shared_ptr<cheat::ParticleInventory> partInv, std::string g4label, double minHitFrac)//The normal constructor. The module using this should get g4label and minHitFrac from the fhicl config.
  {
    fPartInv              = partInv;
    fG4ModuleLabel        = g4label; //This should not be hardcoded. Get advice about fhicl configuration of ServiceProviders.
    fMinHitEnergyFraction = minHitFrac;
  }
  BackTracker::BackTracker(std::shared_ptr<cheat::ParticleInventory> partInv)//Constructor with default values for Gallery users that don't want to manually configure parameters.
  {
    fPartInv              = partInv;
    fG4ModuleLabel        = "largeant";
    fMinHitEnergyFraction = 0.001;
  }

  BackTracker::~BackTracker()
  {
  }

}//End namespace cheat
