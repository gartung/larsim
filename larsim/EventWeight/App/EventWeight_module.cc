////////////////////////////////////////////////////////////////////////
// Class:       EventWeight
// Module Type: producer
// File:        EventWeight_module.cc
//
// Generated at Fri Mar 20 09:36:11 2015 by Zarko Pavlovic using artmod
// from cetpkgsupport v1_08_04.
// 
// Ported from uboonecode to larsim on Feb 14 2018 by Marco Del Tutto
////////////////////////////////////////////////////////////////////////

#include <memory>
#include <iostream>
#include <iomanip>
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "lardataobj/Simulation/sim.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Persistency/Common/Assns.h" 
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/EventWeight/Base/Weight_t.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larsim/EventWeight/Base/WeightCalc.h"
#include "larsim/EventWeight/Base/WeightCalcFactory.h"
#include "larsim/EventWeight/Base/WeightManager.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace evwgh {

  class EventWeight : public art::EDProducer {
  public:
    explicit EventWeight(fhicl::ParameterSet const & p);

    EventWeight(EventWeight const &) = delete;
    EventWeight(EventWeight &&) = delete;
    EventWeight & operator = (EventWeight const &) = delete;
    EventWeight & operator = (EventWeight &&) = delete;
  
  private:
    // Required functions.
    void produce(art::Event & e) override;

    //Optional functions.
    void endJob() override;

    // Run initialization (for deferred loading, after GENIE)
    void Initialize();

    WeightManager _wgt_manager;
    std::string fGenieModuleLabel;
    fhicl::ParameterSet fParameterSet;
    bool fConfigured;
  };

  EventWeight::EventWeight(fhicl::ParameterSet const& p)
      : art::EDProducer{p}, fParameterSet(p), fConfigured(false) {
    auto const n_func = _wgt_manager.Configure(p, *this);

    fGenieModuleLabel = \
      p.get<std::string>("genie_module_label", "generator");

    if (n_func > 0) {
      produces<std::vector<MCEventWeight> >();
    }
  }

  void EventWeight::produce(art::Event & e)
  {
    // Initialize if we haven't yet. This defers setting up GENIE reweighting
    // until after GENIEHelper has configured the tune.
    if (!fConfigured) {
      _wgt_manager.ConfigureFunctions();
    }

    // Implementation of required member function here.
    auto mcwghvec = std::make_unique<std::vector<MCEventWeight>>();
  
    // Get the MC generator information out of the event       
    // these are all handles to mc information.
    std::vector<art::Ptr<simb::MCTruth> > mclist;
  
      // Actually go and get the stuff
    auto const mcTruthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);
      art::fill_ptr_vector(mclist, mcTruthHandle);    
  
    // Loop over all neutrinos in this event
    for (unsigned int inu = 0; inu < mclist.size(); ++inu) {
      auto const mcwgh = _wgt_manager.Run(e, inu);
      mcwghvec->push_back(mcwgh);
    }
  
    e.put(std::move(mcwghvec));
  }

  void EventWeight::endJob()
  {
    // Get the map from sting to Weight_t from the manager
    std::map<std::string, Weight_t*> weightCalcMap = \
      _wgt_manager.GetWeightCalcMap();

    std::stringstream job_summary;
    job_summary  <<  std::setprecision(2);
    for (int i=1; i <= 110 ;i++) job_summary << "=";
    job_summary << std::endl;
    job_summary << std::setw(20) << "WeightCalc"
                << std::setw(15) << "Type"
                << std::setw(15) << "#RW neutrinos"
                << std::setw(15) << "#Multisims"
                << std::setw(15) << "Min"
                << std::setw(15) << "Max"
                << std::setw(15) << "Avg"
                << std::endl;
    for (int i=1; i <= 110; i++) job_summary << "=";
    job_summary << std::endl;
    for (auto it = weightCalcMap.begin(); it!=weightCalcMap.end(); it++) {
      job_summary << std::setw(20) << it->first
                  << std::setw(15) << (it->second->fWeightCalcType)
                  << std::setw(15) << (it->second->fNcalls)
                  << std::setw(15) << (it->second->fNmultisims)
                  << std::setw(15) << (it->second->fMinWeight)
                  << std::setw(15) << (it->second->fMaxWeight)
                  << std::setw(15) << (it->second->fAvgWeight)
                  << std::endl;
    }
    for (int i=1; i<=110; i++) job_summary << "=";
    job_summary << std::endl;
    mf::LogInfo("") << job_summary.str();
  }

}  // namespace evwgh

DEFINE_ART_MODULE(evwgh::EventWeight)

