////////////////////////////////////////////////////////////////////////
// \file PhotonVisibilityService.h
//
// \brief Service to report opdet visibility to different points in
//         the system
//
// \author bjpjones@mit.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef PHOTONVISIBILITYSERVICE_H
#define PHOTONVISIBILITYSERVICE_H


#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "larsim/PhotonPropagation/PhotonLibrary.h"

#include "TF1.h"

///General LArSoft Utilities
namespace phot{
  
  class PhotonVisibilityService {
  public:
    
    PhotonVisibilityService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    void reconfigure(fhicl::ParameterSet const& p);
    
    double GetQuenchingFactor(double dQdx) const;
    
    static double DistanceToOpDet(          double const* xyz, unsigned int OpDet );
    static double SolidAngleFactor(         double const* xyz, unsigned int OpDet );
    float GetVisibility(                    double const* xyz, unsigned int OpChannel, bool wantReflected=false ) const;         

    float const* GetAllVisibilities( double const* xyz, bool wantReflected=false ) const;

    float GetLibraryEntry( int VoxID, int OpChannel, bool wantReflected=false ) const;
    float const* GetLibraryEntries( int VoxID, bool wantReflected=false ) const;
    float const* GetReflT0s( double const* xyz ) const;
    float const* GetLibraryReflT0Entries( int VoxID ) const;

    float GetLibraryReflT0Entry( int VoxID, int Channel ) const;
    void SetDirectLightPropFunctions(TF1 const* functions[8], double& d_break, double& d_max, double& tf1_sampling_factor) const;
    void SetReflectedCOLightPropFunctions(TF1 const* functions[5], double& t0_max, double& t0_break_point) const;
    
    bool UseParameterization() const {return fParameterization;}
    bool StoreReflected() const { return fStoreReflected; }
    bool StoreReflT0() const { return fStoreReflT0; }
    bool IncludePropTime() const { return fIncludePropTime; }

    sim::PhotonVoxelDef GetVoxelDef() const {return fVoxelDef; }
    size_t NOpChannels() const;
    
  private:
    void LoadLibrary() const;

    float  fXmin, fXmax;
    float  fYmin, fYmax;
    float  fZmin, fZmax;
    int    fNx, fNy, fNz;

    bool fUseCryoBoundary;
    
    bool                 fDoNotLoadLibrary;
    bool                 fParameterization;
    bool                 fStoreReflected;
    bool                 fStoreReflT0;
    bool                 fIncludePropTime;
    TF1 *fparslogNorm;
    TF1 *fparslogNorm_far;
    TF1 *fparsMPV;
    TF1 *fparsMPV_far;
    TF1 *fparsWidth;
    TF1 *fparsCte;
    TF1 *fparsCte_far;
    TF1 *fparsSlope;
    double fD_break, fD_max, fTF1_sampling_factor;
    TF1 *fparslogNorm_refl;
    TF1 *fparsMPV_refl;
    TF1 *fparsWidth_refl;
    TF1 *fparsCte_refl;
    TF1 *fparsSlope_refl;
    double fT0_max, fT0_break_point;
   
    std::string          fLibraryFile;      
    mutable PhotonLibrary* fTheLibrary;
    sim::PhotonVoxelDef  fVoxelDef;
    
    
  }; // class PhotonVisibilityService
} //namespace phot
DECLARE_ART_SERVICE(phot::PhotonVisibilityService, LEGACY)
#endif // UTIL_DETECTOR_PROPERTIES_H
