////////////////////////////////////////////////////////////////////////
/// \file  RadioGen_module.cc
/// \brief Generator for radiological decays
///
/// Module designed to produce a set list of particles for MC to model radiological decays
///
/// \author  trj@fnal.gov
//           Rn222 generation feature added by gleb.sinev@duke.edu
//           (based on a generator by jason.stock@mines.sdsmt.edu)
//           JStock. Added preliminary changes to get ready for Ar42 implimentation. This includes allowing for multiple particles from the same decay.
//           Ar42 generation added by JStock (jason.stock@mines.sdsmt.edu).
//             Ar42 is designed to handle 5 separate different
//             beta decay modes, each with it's own chains of
//             possible dexcitation gammas. Because of the
//             high energies, and the relatively high rate
//             expected in the DUNE FD, these chains are
//             completely simulated instead of relying on the
//             existing machinery in this module. To make the
//             treatment of multiple decay products and
//             dexcitation chains generally available would
//             require a significant redesign of the module
//             and possibly of the .root files data structure
//             for each radiological.
//
//           Dec 01, 2017 JStock
//             Adding the ability to make 8.997 MeV Gammas for
//             Ni59 Calibration sources. This is another "hacky"
//             fix to something that really deserves a more
//             elegant and comprehensive solution.
////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <string>
#include <regex>
#include <cmath>
#include <memory>
#include <iterator>
#include <sys/stat.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>
#include <chrono>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "cetlib/exempt_ptr.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

// root includes

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TVector3.h"

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

namespace simb { class MCTruth; }

namespace evgen {
  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes

  class RadioGen : public art::EDProducer {
  public:
    explicit RadioGen(fhicl::ParameterSet const& pset);

  private:
    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void beginJob();
    void endJob();

    typedef int            ti_PDGID;  // These typedefs may look odd, and unecessary. I chose to use them to make the tuples I use later more readable. ti, type integer :JStock
    typedef double         td_Mass;   // These typedefs may look odd, and unecessary. I chose to use them to make the tuples I use later more readable. td, type double  :JStock

    struct ParticleInfo{
      ti_PDGID pdg;
      td_Mass mass;
      TLorentzVector pos;
      TLorentzVector mom;
    };
    
    void SampleOne(int i, simb::MCTruth &mct);

    TLorentzVector dirCalc(double p, double m);

    void readfile(std::string type, std::string const& filename);
    void samplespectrum(std::string nuclideName, ParticleInfo& part);
    void sample_beta_decay_spectrum(double Q, ParticleInfo& part);// PLasorak: A simple rejection method with beta spectrum

    TVector3 GetGoodPosition(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, std::string material, bool& flag);
    
    void Ar42Gamma2(std::vector<ParticleInfo>& v_prods);
    void Ar42Gamma3(std::vector<ParticleInfo>& v_prods);
    void Ar42Gamma4(std::vector<ParticleInfo>& v_prods);
    void Ar42Gamma5(std::vector<ParticleInfo>& v_prods);

    // recoded so as to use the LArSoft-managed random number generator
    double samplefromth1d(TH1D& hist);

    // itype = pdg code:  1000020040: alpha.  itype=11: beta. -11: positron,  itype=22: gamma.  -1: error
    // t is the kinetic energy in GeV  (= E-m)
    // m = mass in GeV
    // p = momentum in GeV

    //double betaphasespace(double mass, double q); // older parameterization.

    // the generator randomly samples points in a rectangular prism of space and time, and only selects those points in
    // volumes with materials that match the regexes in fMaterial.  One can use wildcards * and ? for broader matches.

    std::vector<std::string> fNuclide;   ///< List of nuclides to simulate.  Example:  "39Ar".
    std::vector<std::string> fMaterial;  ///< List of regexes of materials in which to generate the decays.  Example: "LAr"
    std::vector<double> fBq;             ///< Radioactivity in Becquerels (decay per sec) per cubic cm.
    std::vector<double> fT0;             ///< Beginning of time window to simulate in ns
    std::vector<double> fT1;             ///< End of time window to simulate in ns
    std::vector<double> fX0;             ///< Bottom corner x position (cm) in world coordinates
    std::vector<double> fY0;             ///< Bottom corner y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Bottom corner z position (cm) in world coordinates
    std::vector<double> fX1;             ///< Top corner x position (cm) in world coordinates
    std::vector<double> fY1;             ///< Top corner y position (cm) in world coordinates
    std::vector<double> fZ1;             ///< Top corner z position (cm) in world coordinates

    bool                fIsFirstSignalSpecial;
    int trackidcounter;                  ///< Serial number for the MC track ID


    // leftovers from the phase space generator
    // const double gevperamu = 0.931494061;
    // TGenPhaseSpace rg;  // put this here so we don't constantly construct and destruct it


    
    std::vector<std::string> spectrumname;
    
    std::map<std::string,std::unique_ptr<TH1D>> alphaspectrum;
    std::map<std::string,double> alphaintegral;
    
    std::map<std::string,std::unique_ptr<TH1D>> betaspectrum;
    std::map<std::string,double> betaintegral;
    
    std::map<std::string,std::unique_ptr<TH1D>> gammaspectrum;
    std::map<std::string,double> gammaintegral;
    
    std::map<std::string,std::unique_ptr<TH1D>> neutronspectrum;
    std::map<std::string,double> neutronintegral;
    
    CLHEP::HepRandomEngine& fEngine;
    CLHEP::RandFlat        fRandomFlat;
    CLHEP::RandExponential fRandomExponential;
    CLHEP::RandPoisson     fRandomPoisson;
    art::ServiceHandle<geo::Geometry const> fGeo;
    TGeoManager* fGeoManager;
    int nevent;
    TH2D* pos_xy_TH2D;
    TH2D* pos_xz_TH2D;
    TH1D* dir_x_TH1D ;
    TH1D* dir_y_TH1D ;
    TH1D* dir_z_TH1D ;
    TH1D* pdg_TH1D   ;
    TH1D* mom_TH1D   ;
    TH1D* time_TH1D  ;
  };
}

namespace {
  
  constexpr double m_e = 0.000510998928;  // mass of electron in GeV
  constexpr double m_alpha = 3.727379240; // mass of an alpha particle in GeV
  constexpr double m_neutron = 0.9395654133; // mass of a neutron in GeV
  
  constexpr int kAlphaPDG    = 1000020040;
  constexpr int kElectronPDG = 11;
  constexpr int kGammaPDG    = 22;
  constexpr int kNeutronPDG  = 2112;
}

namespace evgen{

  void RadioGen::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    nevent = 0;

    double xmin = *std::min_element(fX0.begin(), fX0.end());
    double xmax = *std::max_element(fX1.begin(), fX1.end());
    double ymin = *std::min_element(fY0.begin(), fY0.end());
    double ymax = *std::max_element(fY1.begin(), fY1.end());
    double zmin = *std::min_element(fZ0.begin(), fZ0.end());
    double zmax = *std::max_element(fZ1.begin(), fZ1.end());
    double tmin = *std::min_element(fT0.begin(), fT0.end());
    double tmax = *std::max_element(fT1.begin(), fT1.end());
    std::cout << "xmin " << xmin << "\n";
    std::cout << "xmax " << xmax << "\n";
    std::cout << "ymin " << ymin << "\n";
    std::cout << "ymax " << ymax << "\n";
    std::cout << "zmin " << zmin << "\n";
    std::cout << "zmax " << zmax << "\n";

    
    pos_xy_TH2D = tfs->make<TH2D>("posXY", ";X [cm];Y [cm]", 100, xmin, xmax, 100, ymin, ymax);
    pos_xz_TH2D = tfs->make<TH2D>("posXZ", ";X [cm];Z [cm]", 100, xmin, xmax, 100, zmin, zmax);
    dir_x_TH1D  = tfs->make<TH1D>("dirX", ";X momentum projection", 100, -1, 1);
    dir_y_TH1D  = tfs->make<TH1D>("dirY", ";Y momentum projection", 100, -1, 1);
    dir_z_TH1D  = tfs->make<TH1D>("dirZ", ";Z momentum projection", 100, -1, 1);
    pdg_TH1D    = tfs->make<TH1D>("PDG", ";PDG;n particles", 100, 0, 100);
    mom_TH1D    = tfs->make<TH1D>("Momentum", ";Momentum [MeV];n particles", 50, 0, 50);
    time_TH1D   = tfs->make<TH1D>("Time", ";Time[ns];n particles", 100, tmin, tmax);
  }
  
  void RadioGen::endJob(){
    pos_xy_TH2D->Scale(1./nevent);
    pos_xz_TH2D->Scale(1./nevent);
    dir_x_TH1D ->Scale(1./nevent);
    dir_y_TH1D ->Scale(1./nevent);
    dir_z_TH1D ->Scale(1./nevent);
    pdg_TH1D   ->Scale(1./nevent);
    mom_TH1D   ->Scale(1./nevent);
    time_TH1D  ->Scale(1./nevent);
  }

//____________________________________________________________________________
  RadioGen::RadioGen(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fNuclide {pset.get< std::vector<std::string>>("Nuclide")}
    , fMaterial{pset.get< std::vector<std::string>>("Material")}
    , fBq{pset.get< std::vector<double> >("BqPercc")}
    , fT0{pset.get< std::vector<double> >("T0")}
    , fT1{pset.get< std::vector<double> >("T1")}
    , fX0{pset.get< std::vector<double> >("X0")}
    , fY0{pset.get< std::vector<double> >("Y0")}
    , fZ0{pset.get< std::vector<double> >("Z0")}
    , fX1{pset.get< std::vector<double> >("X1")}
    , fY1{pset.get< std::vector<double> >("Y1")}
    , fZ1{pset.get< std::vector<double> >("Z1")}
    , fIsFirstSignalSpecial{pset.get< bool >("IsFirstSignalSpecial", false)}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
    , fRandomFlat(fEngine)
    , fRandomExponential(fEngine)
    , fRandomPoisson(fEngine)
    , fGeo()
    , fGeoManager(fGeo->ROOTGeoManager())
    , nevent(0)
    , pos_xy_TH2D  ()
    , pos_xz_TH2D  ()
    , dir_x_TH1D  ()
    , dir_y_TH1D  ()
    , dir_z_TH1D  ()
    , pdg_TH1D ()
    , mom_TH1D ()
    , time_TH1D()
  {
    
    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
 
    // check for consistency of vector sizes
    unsigned int nsize = fNuclide.size();
    if (  fMaterial.size() != nsize ) throw cet::exception("RadioGen") << "Different size Material vector and Nuclide vector";
    if (  fBq.size() != nsize ) throw cet::exception("RadioGen") << "Different size Bq vector and Nuclide vector";
    if (  fT0.size() != nsize ) throw cet::exception("RadioGen") << "Different size T0 vector and Nuclide vector";
    if (  fT1.size() != nsize ) throw cet::exception("RadioGen") << "Different size T1 vector and Nuclide vector";
    if (  fX0.size() != nsize ) throw cet::exception("RadioGen") << "Different size X0 vector and Nuclide vector";
    if (  fY0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y0 vector and Nuclide vector";
    if (  fZ0.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z0 vector and Nuclide vector";
    if (  fX1.size() != nsize ) throw cet::exception("RadioGen") << "Different size X1 vector and Nuclide vector";
    if (  fY1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Y1 vector and Nuclide vector";
    if (  fZ1.size() != nsize ) throw cet::exception("RadioGen") << "Different size Z1 vector and Nuclide vector";

    for(std::string & nuclideName : fNuclide){
      if     (nuclideName=="39Ar"          ){readfile(nuclideName, "Argon_39.root");}
      else if(nuclideName=="60Co"          ){readfile(nuclideName, "Cobalt_60.root");}
      else if(nuclideName=="85Kr"          ){readfile(nuclideName, "Krypton_85.root");}
      else if(nuclideName=="40K"           ){readfile(nuclideName, "Potassium_40.root");}
      else if(nuclideName=="232Th"         ){readfile(nuclideName, "Thorium_232.root");}
      else if(nuclideName=="238U"          ){readfile(nuclideName, "Uranium_238.root");}
      else if(nuclideName=="222Rn"         ){continue;} //Rn222 is handled separately later
      else if(nuclideName=="220Rn"         ){continue;} //Rn220 is handled separately later
      else if(nuclideName=="BiPo(Rn222)"   ){continue;} //BiPo is handled separately later
      else if(nuclideName=="BiPo(Rn220)"   ){continue;} //BiPo is handled separately later
      else if(nuclideName=="59Ni"          ){continue;} //Rn222 is handled separately later
      else if(nuclideName=="42Ar"          ){
        readfile("42Ar_1", "Argon_42_1.root"); //Each possible beta decay mode of Ar42 is given it's own .root file for now.
        readfile("42Ar_2", "Argon_42_2.root"); //This allows us to know which decay chain to follow for the dexcitation gammas.
        readfile("42Ar_3", "Argon_42_3.root"); //The dexcitation gammas are not included in the root files as we want to
        readfile("42Ar_4", "Argon_42_4.root"); //probabilistically simulate the correct coincident gammas, which we cannot guarantee
        readfile("42Ar_5", "Argon_42_5.root"); //by sampling a histogram.
        continue;
      }else{
        std::string searchName = nuclideName;
        searchName+=".root";
        readfile(nuclideName, searchName);
      }
    }
  }

  //____________________________________________________________________________
  void RadioGen::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
  }

  //____________________________________________________________________________
  void RadioGen::produce(art::Event& evt)
  {
    nevent++;
    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);

    trackidcounter = -1;
    for (size_t i=0; i<fNuclide.size(); ++i) {
      SampleOne(i,truth);
    } 

    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  }



  TVector3 RadioGen::GetGoodPosition(double minX, double minY, double minZ, double maxX, double maxY, double maxZ, std::string material, bool& flag)
  {
    int trial = 0;
    TVector3 pos;
    flag = true;
    std::regex regex_material = (std::regex)material;
    
    pos = TVector3(minX + fRandomFlat.fire()*(maxX - minX),
                   minY + fRandomFlat.fire()*(maxY - minY),
                   minZ + fRandomFlat.fire()*(maxZ - minZ));

    std::string volmaterial = fGeoManager->FindNode(pos.X(),pos.Y(),pos.Z())->GetMedium()->GetMaterial()->GetName();
    
    flag = std::regex_match(volmaterial, regex_material);
      
    return pos;
  }
    
  
  //____________________________________________________________________________
  // Generate radioactive decays per nuclide per volume according to the FCL parameters

  void RadioGen::SampleOne(int i, simb::MCTruth &mct)
  {

    // figure out how many decays to generate, assuming that the entire prism consists of the radioactive material.
    // we will skip over decays in other materials later.
    //auto start = std::chrono::steady_clock::now();

    double rate = fabs( fBq[i] * (fT1[i] - fT0[i]) * (fX1[i] - fX0[i]) * (fY1[i] - fY0[i]) * (fZ1[i] - fZ0[i]) ) / 1.0E9;
    long ndecays = fRandomPoisson.shoot(rate);

    for (unsigned int idecay=0; idecay<ndecays; idecay++) {
      // generate just one particle at a time.  Need a little recoding if a radioactive
      // decay generates multiple daughters that need simulation
      // uniformly distributed in position and time
      //
      // JStock: Leaving this as a single position for the decay products. For now I will assume they all come from the same spot.

      // Moved pdgid into the next statement, so that it is localized.
      // electron=11, photon=22, alpha = 1000020040, neutron = 2112

      //JStock: Allow us to have different particles from the same decay. This requires multiple momenta.
      bool flag = true;
      TVector3 pos = GetGoodPosition(fX0[i],fY0[i],fZ0[i],fX1[i],fY1[i],fZ1[i],fMaterial[i], flag);
      if (!flag)
        continue;
      
      double time = (idecay==0 && fIsFirstSignalSpecial) ? 0 : ( fT0[i] + fRandomFlat.fire()*(fT1[i] - fT0[i]));

      std::vector<ParticleInfo> v_prods; //(First is for PDGID, second is mass, third is Momentum)
      
      if (fNuclide[i] == "222Rn")          // Treat 222Rn separately
      {
        double p=0; double t=0.0055904; td_Mass m=m_alpha; ti_PDGID pdgid=kAlphaPDG;
        double energy = t + m; 
        double p2     = energy*energy - m*m;
        if (p2 > 0) p = TMath::Sqrt(p2);
        else        p = 0;
        ParticleInfo part;
        part.mass = m_alpha;
        part.pdg = kAlphaPDG;
        part.mom = dirCalc(p, m);
        part.pos = TLorentzVector(pos, time);
        v_prods.emplace_back(part);
      }
      
      else if (fNuclide[i] == "220Rn")          // Treat 220Rn separately
      {
        double p=0; double t=0.00640474; td_Mass m=m_alpha; ti_PDGID pdgid=kAlphaPDG;
        double energy = t + m; 
        double p2     = energy*energy - m*m;
        if (p2 > 0) p = TMath::Sqrt(p2);
        else        p = 0;
        ParticleInfo part;
        part.mass = m_alpha;
        part.pdg = kAlphaPDG;
        part.mom = dirCalc(p, m);
        part.pos = TLorentzVector(pos, time);
        v_prods.emplace_back(part);
      }
      
      else if (fNuclide[i] == "BiPo_Rn222")          // Treat 222Rn separately
      {
        ParticleInfo part1;
        sample_beta_decay_spectrum(0.003269, part1);
        part1.pos = TLorentzVector(pos, time);
        v_prods.emplace_back(part1);

        ParticleInfo part2;
        part2.mass = m_alpha;
        part2.pdg = kAlphaPDG;
        double t = 0.00783354;
        double energy = t + part2.mass; 
        double p2     = energy*energy - part2.mass*part2.mass;
        double p = 0;
        if (p2 > 0) p = TMath::Sqrt(p2);
        else        p = 0;
        part2.mom = dirCalc(p, part2.mass);
        part2.pos = TLorentzVector(pos, time+fRandomExponential.fire(164000));
        v_prods.emplace_back(part2);
      }
      
      else if (fNuclide[i] == "BiPo_Rn220")          // Treat 222Rn separately
      {
        ParticleInfo part1;
        sample_beta_decay_spectrum(0.0022515, part1);
        part1.pos = TLorentzVector(pos, time);
        v_prods.emplace_back(part1);

        ParticleInfo part2;
        part2.mass = m_alpha;
        part2.pdg = kAlphaPDG;
        double t = 0.0089542;
        double energy = t + part2.mass; 
        double p2     = energy*energy - part2.mass*part2.mass;
        double p = 0;
        if (p2 > 0) p = TMath::Sqrt(p2);
        else        p = 0;
        part2.mom = dirCalc(p, part2.mass);
        part2.pos = TLorentzVector(pos, time+fRandomExponential.fire(299));
        v_prods.emplace_back(part2);
      }
      
      else if(fNuclide[i] == "59Ni"){ //Treat 59Ni Calibration Source separately (as I haven't made a spectrum for it, and ultimately it should be handeled with multiple particle outputs.
        double p=0.008997; // td_Mas=double. ti_PDFID=int. Assigning p directly, as t=p for gammas.
        ParticleInfo part;
        part.mass = 0;
        part.pdg = kGammaPDG;
        part.mom = dirCalc(p, part.mass);
        part.pos = TLorentzVector(pos, time);
        v_prods.emplace_back(part);
      }//end special case Ni59 calibration source
      
      else if(fNuclide[i] == "42Ar")
      {   // Spot for special treatment of Ar42.
        ParticleInfo part;

        double bSelect = fRandomFlat.fire();   //Make this a random number from 0 to 1.
        if(bSelect<0.819){              //beta channel 1. No Gamma. beta Q value 3525.22 keV
          samplespectrum("42Ar_1", part);
          v_prods.emplace_back(part);
          //No gamma here.
        }else if(bSelect<0.9954){       //beta channel 2. 1 Gamma (1524.6 keV). beta Q value 2000.62
          samplespectrum("42Ar_2", part);
          v_prods.emplace_back(part);
          Ar42Gamma2(v_prods);
        }else if(bSelect<0.9988){       //beta channel 3. 1 Gamma Channel. 312.6 keV + gamma 2. beta Q value 1688.02 keV
          samplespectrum("42Ar_3", part);
          v_prods.emplace_back(part);
          Ar42Gamma3(v_prods);
        }else if(bSelect<0.9993){       //beta channel 4. 2 Gamma Channels. Either 899.7 keV (i 0.052) + gamma 2 or 2424.3 keV (i 0.020). beta Q value 1100.92 keV
          samplespectrum("42Ar_4", part);
          v_prods.emplace_back(part);
          Ar42Gamma4(v_prods);
        }else{                          //beta channel 5. 3 gamma channels. 692.0 keV + 1228.0 keV + Gamma 2 (i 0.0033) ||OR|| 1021.2 keV + gamma 4 (i 0.0201) ||OR|| 1920.8 keV + gamma 2 (i 0.041). beta Q value 79.82 keV
          samplespectrum("42Ar_5", part);
          v_prods.emplace_back(part);
          Ar42Gamma5(v_prods);
        }
        
        for (auto& it: v_prods) { // All the Ar42 are produced at the same place.
          it.pos = TLorentzVector(pos, time);
        }
      }
      
      else
      { //General Case.
        ParticleInfo part;
        samplespectrum(fNuclide[i],part);
        part.pos = TLorentzVector(pos, time);
        v_prods.push_back(part);


      }//end else (not RN or other special case


      //JStock: Modify this to now loop over the v_prods.
      for(auto prodEntry : v_prods){
        // set track id to a negative serial number as these are all primary particles and have id <= 0
        int trackid = trackidcounter;
        ti_PDGID pdgid = prodEntry.pdg;
        td_Mass  m = prodEntry.mass;
        TLorentzVector mom_vec = prodEntry.mom;
        TLorentzVector pos_vec = prodEntry.pos;
        trackidcounter--;
        std::string primary("primary");

        // alpha particles need a little help since they're not in the TDatabasePDG table
        // // so don't rely so heavily on default arguments to the MCParticle constructor
        if (pdgid == 1000020040){
          simb::MCParticle part(trackid, pdgid, primary,-1,m,1);
          part.AddTrajectoryPoint(pos_vec, mom_vec);
          pdg_TH1D->Fill(80);
          mct.Add(part);
        }// end "If alpha"
        else{
          simb::MCParticle part(trackid, pdgid, primary);
          part.AddTrajectoryPoint(pos_vec, mom_vec);
          pdg_TH1D->Fill(pdgid);
          mct.Add(part);
        }// end All standard cases.
        pos_xy_TH2D->Fill(prodEntry.pos.X(), prodEntry.pos.Y());
        pos_xz_TH2D->Fill(prodEntry.pos.X(), prodEntry.pos.Z());
        double mom = sqrt(prodEntry.mom.X() * prodEntry.mom.X() +
                          prodEntry.mom.Y() * prodEntry.mom.Y() +
                          prodEntry.mom.Z() * prodEntry.mom.Z());
        mom_TH1D   ->Fill(mom*1000.);
        dir_x_TH1D ->Fill(prodEntry.mom.X()/mom);
        dir_y_TH1D ->Fill(prodEntry.mom.Y()/mom);
        dir_z_TH1D ->Fill(prodEntry.mom.Z()/mom);
        time_TH1D  ->Fill(prodEntry.pos.T());
      }//End Loop over all particles produces in this single decay.
    
    }
 
  }

//Calculate an arbitrary direction with a given magnitude p
  TLorentzVector RadioGen::dirCalc(double p, double m)
  {
  // isotropic production angle for the decay product
  double costheta = (2.0*fRandomFlat.fire() - 1.0);
  if (costheta < -1.0) costheta = -1.0;
  if (costheta > 1.0) costheta = 1.0;
  double const sintheta = sqrt(1.0-costheta*costheta);
  double const phi = 2.0*M_PI*fRandomFlat.fire();
  return TLorentzVector{p*sintheta*std::cos(phi),
      p*sintheta*std::sin(phi),
      p*costheta,
      std::sqrt(p*p+m*m)};
}

  
  std::unique_ptr<TH1D> ParseTGraph(TGraph* graph, std::string type, double& integral) {
    
    if (graph) {
      int np = graph->GetN();
      double *x = graph->GetX();
      double *y = graph->GetY();
      std::string name = "RadioGen_" + type;
      auto hist = std::make_unique<TH1D>(name.c_str(),(name+" Spectrum").c_str(),np,0,np); // TODO!!
      //auto hist = std::make_unique<TH1D>(name.c_str(),(name+" Spectrum").c_str(),np,x[0],x[np-1]);
      //Change the previous line to this at some point
      for (int i=0; i<np; i++) {
        hist->SetBinContent(i+1,y[i]);
        hist->SetBinError  (i+1,0);
      }
      integral = hist->Integral();
      return std::move(hist);
    }
    integral = 0;
    return nullptr;
    
  }

  // only reads those files that are on the fNuclide list.  Copy information from the TGraphs to TH1D's
  void RadioGen::readfile(std::string type, std::string const& filename)
  {
    // bool found{false};

    // for (size_t i=0; i<fNuclide.size(); i++)
    // {
    //   if (fNuclide[i] == nuclide){ //This check makes sure that the nuclide we are searching for is in fact in our fNuclide list. Ar42 handeled separately.
    //       found = true;
    //     break;
    //   } //End If nuclide is in our list. Next is the special case of Ar42
    //     else if (std::regex_match(nuclide, re_argon) && fNuclide[i]=="42Ar") {
    //       found = true;
    //     break;
    //   }
    // }
    
    // if (!found) return;
    //if (type == kUnknown) throw cet::exception("RadioGen") << "The file \"" << filename << "\" will be saved to unknown, so won't be used.";
      
    Bool_t addStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE); // cloned histograms go in memory, and aren't deleted when files are closed.
    // be sure to restore this state when we're out of the routine.

    //spectrumname.push_back(nuclide);
    cet::search_path sp("FW_SEARCH_PATH");
    std::string fn2 = "Radionuclides/";
    fn2 += filename;
    std::string fullname;
    sp.find_file(fn2, fullname);
    struct stat sb;
    if (fullname.empty() || stat(fullname.c_str(), &sb)!=0)
      throw cet::exception("RadioGen") << "Input spectrum file "
        << fn2
        << " not found in FW_SEARCH_PATH!\n";

    TFile f(fullname.c_str(),"READ");
    TGraph *alphagraph   = (TGraph*) f.Get("Alphas");
    TGraph *betagraph    = (TGraph*) f.Get("Betas");
    TGraph *gammagraph   = (TGraph*) f.Get("Gammas");
    TGraph *neutrongraph = (TGraph*) f.Get("Neutrons");


    double integral = 0;
    alphaspectrum  [type] = ParseTGraph(alphagraph,   type, integral); alphaintegral[type]   = integral;
    betaspectrum   [type] = ParseTGraph(betagraph,    type, integral); betaintegral[type]    = integral;
    gammaspectrum  [type] = ParseTGraph(gammagraph,   type, integral); gammaintegral[type]   = integral;
    neutronspectrum[type] = ParseTGraph(neutrongraph, type, integral); neutronintegral[type] = integral;

    f.Close();
    TH1::AddDirectory(addStatus);

    double total = alphaintegral[type] + betaintegral[type] + gammaintegral[type] + neutronintegral[type];
    if (total>0)
    {
      alphaintegral  [type] /= total;
      betaintegral   [type] /= total;
      gammaintegral  [type] /= total;
      neutronintegral[type] /= total;
    }
  }

  void RadioGen::sample_beta_decay_spectrum(double Q, ParticleInfo& part) {  // PLasorak: A simple rejection method with beta spectrum
    
    
  }

  void RadioGen::samplespectrum(std::string type, ParticleInfo& part)
  {

    double t = 0;
    double rtype = fRandomFlat.fire();

    part.pdg = -1;
    part.mass = 0;
    part.mom = TLorentzVector(0,0,0,0);
    
    for (int itry=0;itry<10;itry++) // maybe a tiny normalization issue with a sum of 0.99999999999 or something, so try a few times.
    {
      if (rtype <= alphaintegral.at(type) && alphaspectrum.at(type) != nullptr)
      {
        part.pdg = 1000020040; // alpha
        part.mass = m_alpha;
        t = samplefromth1d(*alphaspectrum.at(type))/1000000.0;
      }
      else if (rtype <= alphaintegral.at(type)+betaintegral.at(type) && betaspectrum.at(type) != nullptr)
      {
        part.pdg = 11; // beta
        part.mass = m_e;
        t = samplefromth1d(*betaspectrum.at(type))/1000000.0;
      }
      else if ( rtype <= alphaintegral.at(type) + betaintegral.at(type) + gammaintegral.at(type) && gammaspectrum.at(type) != nullptr)
      {
        part.pdg = 22; // gamma
        part.mass = 0;
        t = samplefromth1d(*gammaspectrum.at(type))/1000000.0;
      }
      else if( neutronspectrum.at(type) != nullptr)
      {
        part.pdg = 2112;
        part.mass = m_neutron;
        t = samplefromth1d(*neutronspectrum.at(type))/1000000.0;
      }
      
      if (part.pdg >= 0) break;
    }
    
    if (part.pdg == -1)
    {
      throw cet::exception("RadioGen") << "Normalization problem with nuclide: " << type;
    }
    
    double e = t + part.mass;
    double p = e*e - part.mass* part.mass;
    if (p>=0)
    { p = TMath::Sqrt(p); }
    else
    { p=0; }
    part.mom = dirCalc(p, part.mass);
  }

  // this is just a copy of TH1::GetRandom that uses the art-managed CLHEP random number generator instead of gRandom
  // and a better handling of negative bin contents

  double RadioGen::samplefromth1d(TH1D& hist)
  {
    //return hist.GetRandom(); // TODO Change the folling code to that
    int nbinsx = hist.GetNbinsX();
    std::vector<double> partialsum;
    partialsum.resize(nbinsx+1);
    partialsum[0] = 0;

    for (int i=1;i<=nbinsx;i++)
    {
      double hc = hist.GetBinContent(i);
      if (hc < 0) throw cet::exception("RadioGen") << "Negative bin:  " << i << " " << hist.GetName() << "\n";
      partialsum[i] = partialsum[i-1] + hc;
    }
    double integral = partialsum[nbinsx];
    if (integral == 0) return 0;
    // normalize to unit sum
    for (int i=1;i<=nbinsx;i++) partialsum[i] /= integral;

    double r1 = fRandomFlat.fire();
    int ibin = TMath::BinarySearch(nbinsx,&(partialsum[0]),r1);
    Double_t x = hist.GetBinLowEdge(ibin+1);
    if (r1 > partialsum[ibin]) {
      x += hist.GetBinWidth(ibin+1)*(r1-partialsum[ibin])/(partialsum[ibin+1] - partialsum[ibin]);
    }
    //std::cout << x << std::endl;
    return x;
  }


  //Ar42 uses BNL tables for K-42 from Aug 2017
  //beta channel 1. No Gamma. beta Q value 3525.22 keV
  //beta channel 2. 1 Gamma (1524.6 keV). beta Q value 2000.62
  //beta channel 3. 1 Gamma Channel. 312.6 keV + gamma 2. beta Q value 1688.02 keV
  //beta channel 4. 2 Gamma Channels. Either 899.7 keV (i 0.052) + gamma 2 or 2424.3 keV (i 0.020). beta Q value 1100.92 keV
  //beta channel 5. 3 gamma channels. 692.0 keV + 1228.0 keV + Gamma 2 (i 0.0033) ||OR|| 1021.2 keV + gamma 4 (i 0.0201) ||OR|| 1920.8 keV + gamma 2 (i 0.041). beta Q value 79.82 keV
  //No Ar42Gamma1 as beta channel 1 does not produce a dexcitation gamma.
  void RadioGen::Ar42Gamma2(std::vector<ParticleInfo>& v_prods)
  {
    std::vector<double> vd_p = {.0015246};//Momentum in GeV
    for(auto p : vd_p){
      ParticleInfo part;
      part.pdg = kGammaPDG;
      part.mass = 0;
      part.mom = dirCalc(p, part.mass);
      v_prods.emplace_back(part);
    }
  }

  void RadioGen::Ar42Gamma3(std::vector<ParticleInfo>& v_prods)
  {
    std::vector<double> vd_p = {.0003126};
    for(auto p : vd_p){
      ParticleInfo part;
      part.pdg = kGammaPDG;
      part.mass = 0;
      part.mom = dirCalc(p, part.mass);
      v_prods.emplace_back(part);
    }
    Ar42Gamma2(v_prods);
  }

  void RadioGen::Ar42Gamma4(std::vector<ParticleInfo>& v_prods)
  {

    double chan1 = (0.052 / (0.052+0.020) );
    
    if(fRandomFlat.fire()<chan1){
      std::vector<double> vd_p = {.0008997};//Momentum in GeV
      for(auto p : vd_p){
        ParticleInfo part;
        part.pdg = kGammaPDG;
        part.mass = 0;
        part.mom = dirCalc(p, part.mass);
        v_prods.emplace_back(part);
      }
      Ar42Gamma2(v_prods);
    }

    else{
      std::vector<double> vd_p = {.0024243};//Momentum in GeV
      for(auto p : vd_p){
        ParticleInfo part;
        part.pdg = kGammaPDG;
        part.mass = 0;
        part.mom = dirCalc(p, part.mass);
        v_prods.emplace_back(part);
      }
    }
  }

  void RadioGen::Ar42Gamma5(std::vector<ParticleInfo>& v_prods)
  {

    double chan1 = ( 0.0033 / (0.0033 + 0.0201 + 0.041) );
    double chan2 = ( 0.0201 / (0.0033 + 0.0201 + 0.041) );

    double chanPick = fRandomFlat.fire();

    if(chanPick < chan1){
      std::vector<double> vd_p = {0.000692, 0.001228};//Momentum in GeV
      for(auto p : vd_p){
        ParticleInfo part;
        part.pdg = kGammaPDG;
        part.mass = 0;
        part.mom = dirCalc(p, part.mass);
        v_prods.emplace_back(part);
      }
      Ar42Gamma2(v_prods);
    }

    else if (chanPick<(chan1+chan2)){
      std::vector<double> vd_p = {0.0010212};//Momentum in GeV
      for(auto p : vd_p){
        ParticleInfo part;
        part.pdg = kGammaPDG;
        part.mass = 0;
        part.mom = dirCalc(p, part.mass);
        v_prods.emplace_back(part);
      }
      Ar42Gamma4(v_prods);
    }

    else{
      std::vector<double> vd_p = {0.0019208};//Momentum in GeV
      for(auto p : vd_p){
        ParticleInfo part;
        part.pdg = kGammaPDG;
        part.mass = 0;
        part.mom = dirCalc(p, part.mass);
        v_prods.emplace_back(part);
      }
      Ar42Gamma2(v_prods);
    }
  }

  // phase space generator for beta decay -- keep it as a comment in case we ever want to revive it

  // double RadioGen::betaphasespace(double mass, double q)
  //{
  //  CLHEP::RandFlat     flat(fEngine);
  //  double p = 0;
  //  double mi = mass+q+m_e;
  //  TLorentzVector p0(0,0,0,mi);      // pre-decay four-vector
  //  double masses[3] = {0,m_e,mass};  // neutrino, electron, nucleus
  //  rg.SetDecay(p0,3,masses);
  //  double wmax = rg.GetWtMax();
  //  for (int igen=0;igen<1000;igen++)   // cap the retries at 1000
  //    {
  //      double weight = rg.Generate();  // want to unweight these if possible
  //      TLorentzVector *e4v = rg.GetDecay(1);  // get electron four-vector
  //      p = e4v->P();
  //  if (weight >= wmax * flat.fire()) break;
  //   }
  //return p;
  //}

}//end namespace evgen

DEFINE_ART_MODULE(evgen::RadioGen)
