#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/LArG4/OpDetPhotonTable.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimPhotons.h"

namespace cheat{


  //----------------------------------------------------------------
  template<typename Evt>
    const bool PhotonBackTracker::CanRun(Evt const& evt) {
      return ! ( evt.isRealData() ) ;
    }


  //----------------------------------------------------------------
  template<typename Evt>
    void PhotonBackTracker::PrepOpDetBTRs(Evt const& evt)
    {
      if(this->BTRsReady()){ return;}
      larg4::OpDetPhotonTable::Instance()->ClearTable(fGeom->NOpDets());
      priv_OpDetBTRs.clear();
      priv_local_OpDetBTRs.clear();
      priv_ptrs_OpDetBTRs.clear();
      typename Evt::template HandleT<std::vector<sim::OpDetBacktrackerRecord>> btrHandle;
      bool has_backtracker_records = evt.getByLabel(fG4ModuleLabel, btrHandle);
      if (!has_backtracker_records) {
        // then there better be sim photons
        auto const& simphotons_handle = evt.template getValidHandle<std::vector<sim::SimPhotons>>(fG4ModuleLabel);
        const std::vector<sim::SimPhotons> simphotons = *simphotons_handle;
        for (const sim::SimPhotons &photons: simphotons) {
          sim::OpDetBacktrackerRecord record(photons.fOpChannel);
          for (const sim::OnePhoton &photon: photons) {
            double xyz[3];
            photon.InitialPosition.GetXYZ(xyz);
            // std::cout << "Adding photon at time: " << photon.Time << std::endl;
            record.AddScintillationPhotons(photon.MotherTrackID, photon.Time, 1, xyz, photon.Energy);
          }
          larg4::OpDetPhotonTable::Instance()->AddOpDetBacktrackerRecord(record);
        }
        priv_local_OpDetBTRs = larg4::OpDetPhotonTable::Instance()->YieldOpDetBacktrackerRecords();
        for (unsigned i = 0; i < priv_local_OpDetBTRs.size(); i++) {
          priv_OpDetBTRs.push_back(&priv_local_OpDetBTRs[i]);
        }
      }
      else {
      //auto const& btrHandle = evt.template getValidHandle < std::vector < sim::OpDetBacktrackerRecord > > (fG4ModuleLabel);
        //      if(btrHandle.failedToGet()){
        /*  mf::LogWarning("PhotonBackTracker") << "failed to get handle to     simb::MCParticle from "
         *              << fG4ModuleLabel
         *                          << ", return";*/ //This is now silent as it is expected to    happen every generation run. It is also temporary while we wait for
        /*if( 0 ){ return;} //Insert check for DivRecs here, or don't use validHandle below.
          auto const& divrecHandle = evt.template getValidHandle <std::vector<sim::OpDetDivRec>>(fWavLabel);
          if(divrecHandle.failedToGet()){
          return;
          }*/
  
        art::fill_ptr_vector(priv_ptrs_OpDetBTRs, btrHandle);
        for (unsigned i = 0; i < priv_ptrs_OpDetBTRs.size(); i++) {
          priv_OpDetBTRs.push_back(priv_ptrs_OpDetBTRs[i].get());
        }
        //art::fill_ptr_vector(priv_DivRecs, divrecHandle);
      }

      auto compareBTRlambda = [](const sim::OpDetBacktrackerRecord *a, const sim::OpDetBacktrackerRecord *b) {return(a->OpDetNum()<b->OpDetNum());};
      if (!std::is_sorted(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda)) 
        std::sort(priv_OpDetBTRs.begin(),priv_OpDetBTRs.end(),compareBTRlambda);
      //auto compareDivReclambda = [](art::Ptr<sim::OpDetDivRec> a, art::Ptr<sim::OpDetDivRec> b) {return(a->OpDetNum() < b->OpDetNum());};
      /*if (!std::is_sorted(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda)) 
        std::sort(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda);*/
      //art::FindManyP<raw::OpDetWaveform, sim::OpDetDivRec> fp(priv_OpDetBTRs, evt, fWavLabel);// fp;
      //art::FindOneP<raw::OpDetWaveform, sim::OpDetDivRec> fp(priv_OpDetBTRs, evt, fWavLabel);// fp;
      //They come in sorted by BTR. Now make an index matched vector of data_t sorted by BTR. No. I need easy, not efficient. Map of DetNum to data_t. data_t is then channel mapped.
      /*
         if (fp.isValid()){
         for( size_t btr_iter=0; btr_iter<priv_OpDetBTRs.size(); ++btr_iter){
         auto btr=priv_OpDetBTRs.at(btr_iter);
         auto od = btr->OpDetNum();
         auto const& dr = fp.data(btr_iter);
         for(auto& d : dr)
         {
         if(!d) continue;
         priv_od_to_DivRec[od]=*d;//->ref();
         }

         }
         }else{throw cet::exception("PhotonBackTracker")<<"find Waveforms and DivRecs from BTRs failed.";}
         */
      // std::cout << "OPDETBACKTRACKER RECORDS\n";
      // for (auto const &record: priv_OpDetBTRs) {
      //   record->Dump(std::cout, "    ");
      // }

      return;
    }

    //----------------------------------------------------------------
    //ToDo: Figure out why I get OpHit* out of here instead of art::Ptr.
    template<typename Evt>
      void PhotonBackTracker::PrepOpFlashToOpHits( Evt const& evt)
      {
        if(this->OpFlashToOpHitsReady()){ return;}
        // std::vector< art::Handle< std::vector < recob::OpFlash >>> flashHandles;
        std::vector<typename Evt::template HandleT<std::vector<recob::OpFlash>>> flashHandles;

        evt.getManyByType(flashHandles);
        for( const auto& handle : flashHandles)
        {
          std::vector< art::Ptr < recob::OpFlash > > flash_vec;
          if(handle.failedToGet())
          {
            mf::LogWarning("PhotonBackTracker")<<" failed to get handle to recob::OpFlash. Has reco run yet?";
            return;
          }
          art::fill_ptr_vector(flash_vec, handle);
          auto tag = art::InputTag( handle.provenance()->moduleLabel() );
          art::FindManyP<recob::OpHit>  flash_hit_assn(flash_vec, evt, tag);
          //          std::cout<<"flash_hit_assn.size: "<<flash_hit_assn.size()<<"\n";
          for ( size_t i = 0; i < flash_vec.size(); ++i)
          {
            art::Ptr< recob::OpFlash > flashp = flash_vec.at(i);
            std::vector< art::Ptr< recob::OpHit > > ophits = flash_hit_assn.at(i);
            auto check = priv_OpFlashToOpHits.emplace(flashp, ophits);
            if ( check.second == false )
            {
              // loop ophit_ps
              // push_back to vector.
              for ( auto& ophitp : ophits )
              {
                check.first->second.push_back(ophitp);
              }
            }
          }
        }
      }

    //----------------------------------------------------------------
    template<typename Evt>
      void PhotonBackTracker::PrepEvent( Evt const& evt)
      {
        if( !(this->CanRun( evt ) ) ){
          throw cet::exception("PhotonBackTracker")
            <<"PhotonBackTracker cannot function."
            <<"Is this file real data?";
        }
        priv_OpDetBTRs.clear();
        this->PrepOpDetBTRs(evt);
        this->PrepOpFlashToOpHits(evt);
      } 
    }
