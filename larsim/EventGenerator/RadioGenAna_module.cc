////////////////////////////////////////////////////////////////////////
// Class:       RadioGenAna
// Plugin Type: analyzer (art v2_11_03)
// File:        RadioGenAna_module.cc
//
// Generated at Tue May 21 17:08:27 2019 by Pierre Lasorak using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TProfile.h"
#include "TProfile2D.h"

class RadioGenAna;


class RadioGenAna : public art::EDAnalyzer {
public:
  explicit RadioGenAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RadioGenAna(RadioGenAna const &) = delete;
  RadioGenAna(RadioGenAna &&) = delete;
  RadioGenAna & operator = (RadioGenAna const &) = delete;
  RadioGenAna & operator = (RadioGenAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void beginJob();
  
private:
  std::string fLabel;
  TProfile2D *XY_Prof;
  TProfile2D *XZ_Prof;
  TProfile   *PDG_Prof;
  TProfile   *Mom_Prof;
  TProfile   *Time_Prof;
  // Declare member data here.

};

void RadioGenAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  XY_Prof   = tfs->make<TProfile2D>();
  XZ_Prof   = tfs->make<TProfile2D>();
  PDG_Prof  = tfs->make<TProfile  >();
  Mom_Prof  = tfs->make<TProfile  >();
  Time_Prof = tfs->make<TProfile  >();
}

RadioGenAna::RadioGenAna(fhicl::ParameterSet const & p):
  EDAnalyzer(p),
  fLabel(p.get<std::string>("MCTruthLabel")){
}

void RadioGenAna::analyze(art::Event const & e){
  
  art::Handle< std::vector<simb::MCTruth> > MCT;
  
  if(evt.getByLabel(fLabel, MCT)) {
    for (int ipart=0; ipart<MCTruth->NParticles(); ++ipart) {
      const simb::MCParticle& particle = MCT->GetParticle(ipart);
      
      XY_Prof  ->Fill(particle.Vx(), particle.Vy());
      XZ_Prof  ->Fill(particle.Vx(), particle.Vz());
      if      (particle.PdgCode() == 1000020040){
        PDG_Prof ->Fill(80);
      } else if (particle.PdgCode() >= 100) {
        PDG_Prof ->Fill(90);
      } else {
        PDG_Prof ->Fill(particle.PdgCode());
      }
      Mom_Prof ->Fill(particle.P());
      Time_Prof->Fill(particle.T());
    }
  }
}

DEFINE_ART_MODULE(RadioGenAna)
