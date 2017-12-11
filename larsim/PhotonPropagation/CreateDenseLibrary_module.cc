// \file CreateDenseLibrary_module.cc
// \author Chris Backhouse - c.backhouse@ucl.ac.uk 2017
// Based on Ben Jones' SimPhotonCounter

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PhotonLibrary.h"
#include "larsim/PhotonPropagation/OpDetResponseInterface.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <vector>

namespace opdet
{
  class CreateDenseLibrary: public art::EDAnalyzer
  {
  public:
    CreateDenseLibrary(const fhicl::ParameterSet&);
    virtual ~CreateDenseLibrary();

    void analyze(const art::Event&);

    void endJob();

  protected:
    void HandlePhoton(int vox, int opchan, double wavelength, double time);

    art::ServiceHandle<phot::PhotonVisibilityService> pvs;
    art::ServiceHandle<opdet::OpDetResponseInterface> odresponse;

    std::vector<int> fProdCount; ///< How many photons were sim'd on each voxel

    phot::PhotonLibrary* fLib;

    std::string fGeneratorLabel;
    std::string fGeantLabel;

    bool fUseLitePhotons;
    bool fStoreReflected;
    bool fStoreReflT0;
  };

  //------------------------------------------------------------
  CreateDenseLibrary::CreateDenseLibrary(const fhicl::ParameterSet& pset)
    : EDAnalyzer(pset),
      fLib(0),
      fGeneratorLabel(pset.get<std::string>("GeneratorLabel")),
      fGeantLabel(pset.get<std::string>("GeantLabel")),
      fUseLitePhotons(pset.get<bool>("UseLitePhotons")),
      fStoreReflected(pset.get<bool>("StoreReflected")),
      fStoreReflT0(pset.get<bool>("StoreReflT0"))
  {
    fProdCount.resize(pvs->GetVoxelDef().GetNVoxels());

    fLib = new phot::PhotonLibrary;
    fLib->CreateEmptyLibrary(pvs->GetVoxelDef().GetNVoxels(),
                             pvs->NOpChannels(),
                             fStoreReflected, fStoreReflT0);

    if(fStoreReflT0){
      // Initialize all the t0s to large numbers
      for(int vox = 0; vox < pvs->GetVoxelDef().GetNVoxels(); ++vox){
        for(unsigned int opchan = 0; opchan < pvs->NOpChannels(); ++opchan){
          fLib->SetReflT0(vox, opchan, 999);
        }
      }
    }
  }

  //------------------------------------------------------------
  CreateDenseLibrary::~CreateDenseLibrary()
  {
    delete fLib;
  }

  //------------------------------------------------------------
  void CreateDenseLibrary::endJob()
  {
    // Loop through the library normalizing each entry by the number of photons
    // simulated on that voxel.
    for(int vox = 0; vox < pvs->GetVoxelDef().GetNVoxels(); ++vox){
      const double denom = fProdCount[vox];
      if(denom == 0) continue;

      for(unsigned int opchan = 0; opchan < pvs->NOpChannels(); ++opchan){
        fLib->SetCount(vox, opchan,
                       fLib->GetCount(vox, opchan)/denom);
        if(fStoreReflected)
          fLib->SetReflCount(vox, opchan,
                             fLib->GetReflCount(vox, opchan)/denom);
      } // end for opchan
    } // end for vox

    fLib->StoreLibraryToFile(fStoreReflected, fStoreReflT0);
  }

  //------------------------------------------------------------
  void CreateDenseLibrary::HandlePhoton(int vox, int opchan,
                                        double wavelength, double time)
  {
    if(!odresponse->detectedLite(opchan)) return;

    if(!fStoreReflected || wavelength < 200){
      //only store direct/UV light
      fLib->SetCount(vox, opchan, fLib->GetCount(vox, opchan) + 1);
    }
    else if(fStoreReflected && wavelength > 380){
      //shifted light is in visible range
      fLib->SetReflCount(vox, opchan, fLib->GetReflCount(vox, opchan) + 1);

      // find the first visible arrival time
      if(fStoreReflT0){
        fLib->SetReflT0(vox, opchan,
                        std::min(float(time), fLib->GetReflT0(vox, opchan)));
      }
    } // end if reflected
  }

  //------------------------------------------------------------
  void CreateDenseLibrary::analyze(const art::Event& evt)
  {
    // Figure out which voxel the light in this event (from LightSource) was
    // generated in.
    art::Handle<std::vector<simb::MCTruth>> truthcol;
    evt.getByLabel(fGeneratorLabel, truthcol);
    assert(truthcol->size() == 1);
    const TVector3 pt = (*truthcol)[0].GetParticle(0).Position().Vect();
    const int vox = pvs->GetVoxelDef().GetVoxelID(pt);
    // Increment the appropriate denominator in the visibility
    fProdCount[vox] += (*truthcol)[0].NParticles();


    if(fUseLitePhotons){
      art::Handle<std::vector<sim::SimPhotonsLite>> photcol;
      evt.getByLabel(fGeantLabel, photcol);

      const double wavelength = 128; // nm

      for(const sim::SimPhotonsLite& phots: *photcol){
        for(auto it: phots.DetectedPhotons){
          //Get arrival time from phot
          const double time = it.first;

          for(int i = 0; i < it.second; i++){
            HandlePhoton(vox, phots.OpChannel, wavelength, time);
          } // end for i
        } // end for it
      } // end for phots
    } // end if lite
    else{ // otherwise use the full photons
      sim::SimPhotonsCollection hitcol = sim::SimListUtils::GetSimPhotonsCollection(evt, fGeantLabel);

      for(auto it: hitcol){
        const int opchan = it.first;
        const sim::SimPhotons& hit = it.second;

        for(const sim::OnePhoton& phot: hit){
          // Calculate wavelength in nm
          const double wavelength = odresponse->wavelength(phot.Energy);
          const double time = phot.Time;

          HandlePhoton(vox, opchan, wavelength, time);
        } // end for phot
      } // end for it (itcol)
    } // end if full photons
  }

  DEFINE_ART_MODULE(CreateDenseLibrary)

}//end namespace opdet

