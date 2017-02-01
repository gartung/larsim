////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.h
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.h,v 1.2 2009/03/31 17:58:39 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// One way to implement voxels in Geant4 is to create a parallel
/// "read-out" geometry along with the real, physical geometry.  The
/// read-out geometry is implemented in LArVoxelReadoutGeometry; this
/// class is the sensitive detector for that geometry.  That is,
/// Geant4 will call this routine every time there is a step within a
/// volume of the read-out geometry; this routine then accumulates
/// information from that step.
///
/// In general, Geant4 expects to have per-event user information
/// attached to the G4Event in some way; their G4VSensitiveDetector
/// class supports this by allowing user-defined hit collections to
/// added to a G4HCOfThisEvent object (a collection of hit
/// collections; yes, it makes my head ache too!) that's part of each
/// G4Event.  
///
/// This class works differently, by accumulating the information in
/// its internal sim::LArVoxelList.  See LArVoxelListAction for how
/// this information is made available to the main LArG4 module.
/// 
/// Why define a parallel geometry?  Here are some reasons:
///
/// - The regular LAr TPC is one large volume of liquid argon.  When
///   Geant4 does its physics modeling, it can be unconstrained in
///   step size by the voxels.  Only for readout would the steps be
///   sub-divided.
///
/// - There may be more than one kind of readout, depending on a
///   detector's instrumentation (e.g., OpDets in addition to the wire
///   planes).  It's possible that the voxelization appropriate for
///   the wire planes may not be an appropriate readout for the other
///   readouts.  Geant4 allows the construction of multiple parallel
///   readouts, so this mechanism is relatively easy to extend for
///   each type of readout.
///
///
///
///  This version has undergone changes from the original LArVoxelReadout to
///  remove post G4 energy deposition steps. Ionization/Scintillation,
///  drifitng, and other elements are now separated.
///
///  Wes, February 2017



#ifndef LArG4_LArVoxelReadout_h
#define LArG4_LArVoxelReadout_h

#include <stdint.h>
#include <vector>

#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/globals.hh"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEDep.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/LArG4/IonizationAndScintillation.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;
namespace CLHEP { class HEPRandomEngine; }

namespace larg4 {

  /// Simple structure holding a TPC and cryostat number
  struct TPCID_t {
    unsigned short int Cryostat, TPC;
    bool operator< (const TPCID_t& than) const
      {
        return (Cryostat < than.Cryostat)
          || ((Cryostat == than.Cryostat) && (TPC < than.TPC));
      } // operator< ()
  }; // TPCID_t
  
  /**
   * @brief A G4PVPlacement with an additional identificator
   * @param IDTYPE type of ID class
   * 
   * This class is a G4PVPlacement with in addition an ID parameter.
   * The ID type is an object which can be default-constructed and copied,
   * better to be a POD.
   * 
   * This being a very stupid utility class, only the constructor that we
   * actually use is available. The others can be implemented in the same way.
   * Also the merry company of copy and move constuctors and operators is left
   * to the good will of the compiler, despite the destructor is specified.
   */
  template <class IDTYPE>
  class G4PVPlacementWithID: public G4PVPlacement {
      public:
    typedef IDTYPE ID_t;
    
    ID_t ID; ///< Physical Volume identificator
    
    /// Constructor
    G4PVPlacementWithID(const G4Transform3D& Transform3D, const G4String &pName,
      G4LogicalVolume* pLogical, G4VPhysicalVolume* pMother,
      G4bool pMany, G4int pCopyNo, G4bool pSurfChk = false,
      ID_t id = ID_t()
      ):
      G4PVPlacement
        (Transform3D, pName, pLogical, pMother, pMany, pCopyNo, pSurfChk),
      ID(id)
      {}
    
    /// Virtual destructor: does nothing more
    virtual ~G4PVPlacementWithID() {}
  }; // G4PVPlacementWithID<>
  
  /// A physical volume with a TPC ID
  typedef G4PVPlacementWithID<TPCID_t> G4PVPlacementInTPC;
  
  
  class LArVoxelReadout : public G4VSensitiveDetector
  {
  public:
    /// Constructor. Can detect which TPC to cover by the name
    LArVoxelReadout(std::string const& name);
    
    /// Constructor. Sets which TPC to work on
    LArVoxelReadout
      (std::string const& name, unsigned int, unsigned int);

    // Destructor
    virtual ~LArVoxelReadout();

    // Required for classes that inherit from G4VSensitiveDetector.
    //
    // Called at start and end of each event.
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);

    // Called to clear any accumulated information.
    virtual void clear();

    // The key method of this class.  It's called by Geant4 for each
    // step within the read-out geometry.  It accumulates the energy
    // in the G4Step in the LArVoxelList.
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory* );

    // Empty methods; they have to be defined, but they're rarely
    // used in Geant4 applications.
    virtual void DrawAll();
    virtual void PrintAll();

    const std::vector<sim::SimEDep>& GetSimEDepCollection()
    { return fSimEDepCol; }
    
  private:

    void ProcessStep(G4Step*);
    
    G4ThreeVector                             fStepMidPoint;
    std::vector<sim::SimEDep>                 fSimEDepCol;

  };

}

#endif // LArG4_LArVoxelReadout_h
