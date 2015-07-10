////////////////////////////////////////////////////////////////////////
/// \file  G4BadIdeaAction.h
/// \brief this UserAction derived class is to implement catches to known bugs
///        in Geant4 that require grabbing const G4 objects and altering them - 
///        a very bad idea in general.  Please do not add to this class without
///        discussing with the LArSoft Conveners
///
/// \version $Id: ParticleListAction.h,v 1.3 2010/04/29 15:39:33 seligman Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

/// This class implements the LArG4::UserAction interface in order to
/// try to side step known bugs in version 4.9.4.p02 of Geant4
//
/// It uses multiple inheritance: it inherits from G4Base::UserAction,
/// in order to take advantage of Geant4's user hooks

#ifndef LArG4_LARLIGHTENERGYACTION_H
#define LArG4_LARLIGHTENERGYACTION_H

#include "G4Base/UserAction.h"
#include <cstring>
#include "Geant4/globals.hh"
#include <map>

#include "Geometry/Geometry.h"

// Forward declarations.
class G4Step;

namespace larg4 {
double energy_deposit_step;
//double phot_refl_ev;
int eventnumber_fast;
  class LaRLightEnergyAction : public g4b::UserAction
  {
  public:
    // Standard constructors and destructors;
    LaRLightEnergyAction(int, std::string, std::string );
    virtual ~LaRLightEnergyAction();

    // UserActions method that we'll override, to obtain access to
    // Geant4's steps
    virtual void SteppingAction    (const G4Step*);

  private:

    art::ServiceHandle<geo::Geometry> fGeo;  //< handle to geometry service
    int fNoIncomingMuons;

  };

} // namespace LArG4

#endif // LArG4_LARLIGHTENERGYACTION_H
