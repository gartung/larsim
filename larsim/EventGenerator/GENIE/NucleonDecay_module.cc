////////////////////////////////////////////////////////////////////////
// Class:       NucleonDecay
// Module Type: producer
// GENIE nucleon decay generator
//
// Converted from gNucleonDecayEvGen.cxx by
// tjyang@fnal.gov
//
//  Nucleon decay mode ID:
// ---------------------------------------------------------
//  ID |   Decay Mode                     |   Current Limit 
//     |                                  |   (1E+34 yrs)
// ---------------------------------------------------------
//   0 |   p --> e^{+}      + \pi^{0}     |   1.3
//   1 |   p --> \mu^{+}    + \pi^{0}     |   1.1
//   2 |   p --> e^{+}      + \eta^{0}    |   0.42
//   3 |   p --> \mu^{+}    + \eta^{0}    |   0.13
//   4 |   p --> e^{+}      + \rho^{0}    |   0.07
//   5 |   p --> \mu^{+}    + \rho^{0}    |   0.02
//   6 |   p --> e^{+}      + \omega^{0}  |   0.03
//   7 |   p --> \mu^{+}    + \omega^{0}  |   0.08
//   8 |   n --> e^{+}      + \pi^{-}     |   0.2
//   9 |   n --> \mu^{+}    + \pi^{-}     |   0.1
//  10 |   p --> \bar{\nu}} + K^{+}       |   0.4
// ---------------------------------------------------------
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "EVGCore/EventRecord.h"
#include "NucleonDecay/NucleonDecayMode.h"

#include <memory>
#include <string>

namespace evgen {
  class NucleonDecay;
}

class evgen::NucleonDecay : public art::EDProducer {
public:
  explicit NucleonDecay(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NucleonDecay(NucleonDecay const &) = delete;
  virtual ~NucleonDecay();
  NucleonDecay(NucleonDecay &&) = delete;
  NucleonDecay & operator = (NucleonDecay const &) = delete;
  NucleonDecay & operator = (NucleonDecay &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  const genie::EventRecordVisitorI * mcgen;
  genie::NucleonDecayMode_t gOptDecayMode    = genie::kNDNull;             // nucleon decay mode
};


evgen::NucleonDecay::NucleonDecay(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  string sname   = "genie::EventGenerator";
  string sconfig = "NucleonDecay";
  genie::AlgFactory * algf = genie::AlgFactory::Instance();
  mcgen =
    dynamic_cast<const genie::EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
    throw cet::exception("NucleonDecay") << "Couldn't instantiate the nucleon decay generator"; 
  }
  int fDecayMode = p.get<int>("DecayMode");
  gOptDecayMode = (genie::NucleonDecayMode_t) fDecayMode;

}

evgen::NucleonDecay::~NucleonDecay(){
  delete mcgen;
}

void evgen::NucleonDecay::produce(art::Event & e)
{
  // Implementation of required member function here.
  genie::EventRecord * event = new genie::EventRecord;
  int target = 1000180400;  //Only use argon target
  int decay  = (int)gOptDecayMode;
  genie::Interaction * interaction = genie::Interaction::NDecay(target,decay);
  event->AttachSummary(interaction);
  
  // Simulate decay     
  mcgen->ProcessEventRecord(event);
  
  LOG_DEBUG("NucleonDecay")
    << "Generated event: " << *event;
  
  delete event;
  return;

}

void evgen::NucleonDecay::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(evgen::NucleonDecay)
