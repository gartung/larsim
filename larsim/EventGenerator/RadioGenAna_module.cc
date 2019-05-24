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
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//#include "lardataobj/Simulation/sim.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TH1D.h"
#include "TH2D.h"

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
  void endJob();

private:
  std::string fLabel;
  int nevent;
  TH2D* pos_xy_TH2D;
  TH2D* pos_xz_TH2D;
  TH1D* dir_x_TH1D ;
  TH1D* dir_y_TH1D ;
  TH1D* dir_z_TH1D ;
  TH1D* pdg_TH1D   ;
  TH1D* mom_TH1D   ;
  TH1D* time_TH1D  ;
// Declare member data here.
  int    nbin_x_pos_xy, nbin_y_pos_xy, nbin_x_pos_xz, nbin_y_pos_xz;
  double min__x_pos_xy, min__y_pos_xy, min__x_pos_xz, min__y_pos_xz;
  double max__x_pos_xy, max__y_pos_xy, max__x_pos_xz, max__y_pos_xz;
  int    nbin_x_dir_x,  nbin_x_dir_y,  nbin_x_dir_z;
  double min__x_dir_x,  min__x_dir_y,  min__x_dir_z;
  double max__x_dir_x,  max__x_dir_y,  max__x_dir_z;
  int    nbin_x_pdg,    nbin_x_mom,    nbin_x_time;
  double min__x_pdg,    min__x_mom,    min__x_time;
  double max__x_pdg,    max__x_mom,    max__x_time;
  
};

void RadioGenAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  nevent = 0;
  pos_xy_TH2D = tfs->make<TH2D>("posXY", ";X [cm];Y [cm]",
                                nbin_x_pos_xy, min__x_pos_xy, max__x_pos_xy,
                                nbin_y_pos_xy, min__y_pos_xy, max__y_pos_xy);
  pos_xz_TH2D = tfs->make<TH2D>("posXZ", ";X [cm];Z [cm]",
                                nbin_x_pos_xz, min__x_pos_xz, max__x_pos_xz,
                                nbin_y_pos_xz, min__y_pos_xz, max__y_pos_xz);
  dir_x_TH1D  = tfs->make<TH1D>("dirX", ";X momentum projection", nbin_x_dir_x, min__x_dir_x, max__x_dir_x);
  dir_y_TH1D  = tfs->make<TH1D>("dirY", ";Y momentum projection", nbin_x_dir_y, min__x_dir_y, max__x_dir_y);
  dir_z_TH1D  = tfs->make<TH1D>("dirZ", ";Z momentum projection", nbin_x_dir_z, min__x_dir_z, max__x_dir_z);
  pdg_TH1D    = tfs->make<TH1D>("PDG", ";PDG;n particles", nbin_x_pdg, min__x_pdg, max__x_pdg);
  mom_TH1D    = tfs->make<TH1D>("Momentum", ";Momentum [MeV];n particles", nbin_x_mom, min__x_mom, max__x_mom);
  time_TH1D   = tfs->make<TH1D>("Time", ";Time[ns];n particles", nbin_x_time, min__x_time, max__x_time);

}

void RadioGenAna::endJob() {
  pos_xy_TH2D->Scale(1./nevent);
  pos_xz_TH2D->Scale(1./nevent);
  dir_x_TH1D ->Scale(1./nevent);
  dir_y_TH1D ->Scale(1./nevent);
  dir_z_TH1D ->Scale(1./nevent);
  pdg_TH1D   ->Scale(1./nevent);
  mom_TH1D   ->Scale(1./nevent);
  time_TH1D  ->Scale(1./nevent);

}

RadioGenAna::RadioGenAna(fhicl::ParameterSet const & p):
  EDAnalyzer(p),
  fLabel(p.get<std::string>("MCTruthLabel")),
  nbin_x_pos_xy (p.get<int   >("nbin_x_pos_xy")),
  nbin_y_pos_xy (p.get<int   >("nbin_y_pos_xy")),
  nbin_x_pos_xz (p.get<int   >("nbin_x_pos_xz")),
  nbin_y_pos_xz (p.get<int   >("nbin_y_pos_xz")),
  min__x_pos_xy (p.get<double>("min__x_pos_xy")),
  min__y_pos_xy (p.get<double>("min__y_pos_xy")),
  min__x_pos_xz (p.get<double>("min__x_pos_xz")),
  min__y_pos_xz (p.get<double>("min__y_pos_xz")),
  max__x_pos_xy (p.get<double>("max__x_pos_xy")),
  max__y_pos_xy (p.get<double>("max__y_pos_xy")),
  max__x_pos_xz (p.get<double>("max__x_pos_xz")),
  max__y_pos_xz (p.get<double>("max__y_pos_xz")),
  nbin_x_dir_x  (p.get<int>   ("nbin_x_dir_x")),
  nbin_x_dir_y  (p.get<int>   ("nbin_x_dir_y")),
  nbin_x_dir_z  (p.get<int>   ("nbin_x_dir_z")),
  min__x_dir_x  (p.get<double>("min__x_dir_x")),
  min__x_dir_y  (p.get<double>("min__x_dir_y")),
  min__x_dir_z  (p.get<double>("min__x_dir_z")),
  max__x_dir_x  (p.get<double>("max__x_dir_x")),
  max__x_dir_y  (p.get<double>("max__x_dir_y")),
  max__x_dir_z  (p.get<double>("max__x_dir_z")),
  nbin_x_pdg    (p.get<int   >("nbin_x_pdg")),
  nbin_x_mom    (p.get<int   >("nbin_x_mom")),
  nbin_x_time   (p.get<int   >("nbin_x_time")),
  min__x_pdg    (p.get<double>("min__x_pdg")),
  min__x_mom    (p.get<double>("min__x_mom")),
  min__x_time   (p.get<double>("min__x_time")),
  max__x_pdg    (p.get<double>("max__x_pdg")),
  max__x_mom    (p.get<double>("max__x_mom")),
  max__x_time   (p.get<double>("max__x_time"))
{
}

void RadioGenAna::analyze(art::Event const & evt){
  ++nevent;
  int Event  = evt.event();
  art::Handle< std::vector<simb::MCTruth> > MCT;
  
  if(evt.getByLabel(fLabel, MCT)) {
    for(size_t i = 0; i < MCT->size(); i++) {
      for (int ipart=0; ipart<MCT->at(i).NParticles(); ++ipart) {
        const simb::MCParticle& particle = MCT->at(i).GetParticle(ipart);
        if (particle.PdgCode() == 1000020040){
          pdg_TH1D->Fill(80);
        }// end "If alpha"
        else{
          pdg_TH1D->Fill(particle.PdgCode());
        }// end All standard cases.
        pos_xy_TH2D->Fill(particle.Vx(), particle.Vy());
        pos_xz_TH2D->Fill(particle.Vx(), particle.Vz());
        mom_TH1D   ->Fill(particle.P()*1000.);
        dir_x_TH1D ->Fill(particle.Px()/particle.P());
        dir_y_TH1D ->Fill(particle.Py()/particle.P());
        dir_z_TH1D ->Fill(particle.Pz()/particle.P());
        time_TH1D  ->Fill(particle.T());
      }
    }
  } else {
    std::cout << "This label wasn't found: " << fLabel << "\n";
  }
}

DEFINE_ART_MODULE(RadioGenAna)
