////////////////////////////////////////////////////////////////////////
/// \file  HSNGen_module.cc
/// \brief Generator for heavy sterile neutrinos based on pre-generated HSN fluxes.
///
/// Generator for heavy sterile neutrino decays inside the MicroBooNE detector.
/// The code is largely based on InFlight generator, an event generator for sterile decays at
/// SBL facilities independent of LArSoft, written by Mark Ross-Lonergan and Peter Ballett.
/// Credit for the code goes to the two original authors.
///
/// Any malfunction, instability or bug is to be attributed solely to
/// Salvatore Davide Porzio, responsible for re-writing the code in its current form
/// and porting it to the LArSoft code.
///
/// \author  salvatore.porzio@postgrad.manchester.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_HSNGen_H
#define EVGEN_HSNGen_H

// ROOT includes
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TString.h"
#include "TSystem.h" //need BaseName and DirName

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcoreobj/SummaryData/RunData.h"

// HSNGen includes
#include "DataObjects/FourMomentum.h"
#include "DataObjects/SterileNeutrino.h"
#include "DataObjects/Flux.h"
#include "DataObjects/Observables.h"
#include "DataObjects/Channel.h"
#include "Helpers/Helper.h"
#include "Helpers/Settings.h"

#include <sqlite3.h> 
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "ifdh.h"  //to handle flux files

namespace hsngen
{
  /// A module to check the results from the Monte Carlo generator
  class HSNGen : public art::EDProducer
  {
  public:
    explicit HSNGen(fhicl::ParameterSet const& pset);
    virtual ~HSNGen();                       
    
    void produce(art::Event& evt);  
    void beginJob();
    void endJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
  private:
    // Fcl settings
    bool fPrintHepEvt;
    double fSterileMass;
    int fDecayChannel;
    std::string fFluxFile;
    double fDistance;
    double fGlobalTimeOffset;
    double fBeamWindow;
    std::vector<double> fBoundariesX;
    std::vector<double> fBoundariesY;
    std::vector<double> fBoundariesZ;

    // Analysis variables
    CLHEP::HepJamesRandom gEngine;
    Settings gSett;
    std::vector<double> gModelParams;
    twoIP_channel *gChan;
    FluxFile gFlux;
    int gFakeRunNumber;

    // Auxiliary functions
    void CompressSettings(Settings &set);
  }; // END class HSNGen

  HSNGen::HSNGen(fhicl::ParameterSet const& p)
  : fPrintHepEvt(p.get<bool>("PrintHepEvt")),
    fSterileMass(p.get<double>("SterileMass")),
    fDecayChannel(p.get<double>("DecayChannel")),
    fFluxFile(p.get<std::string>("FluxFile")),
    fDistance(p.get<double>("Distance")),
    fGlobalTimeOffset(p.get<double>("GlobalTimeOffset")),
    fBeamWindow(p.get<double>("BeamWindow")),
    fBoundariesX(p.get<std::vector<double>>("BoundariesX")),
    fBoundariesY(p.get<std::vector<double>>("BoundariesY")),
    fBoundariesZ(p.get<std::vector<double>>("BoundariesZ"))
  {
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });
    this->reconfigure(p);

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >(); 
    return;
  }
  
  HSNGen::~HSNGen()
  {}

  void HSNGen::reconfigure(fhicl::ParameterSet const& p)
  {

    return;
  }

  void HSNGen::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &gEngine = rng->getEngine("gen");
    CompressSettings(gSett);
    FillModel(gEngine, gChan, gModelParams, gSett);
    gFlux = FluxFile(fFluxFile, fSterileMass);
    gFakeRunNumber = 0;
    return;
  }

  void HSNGen::endJob()
  {
    delete gChan;
    return;
  }


  void HSNGen::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    run.put(std::move(runcol));
    return;
  }

  void HSNGen::produce(art::Event& evt)
  {
    // Declare geometry and MCTruth objects
    art::ServiceHandle<geo::Geometry> geom;
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    truth.SetOrigin(simb::kUnknown);

    // Generate observables characterizing the event
    Observables obs; 
    GenerateObservables(gEngine, gChan, gFlux, gSett, obs);
    if (fPrintHepEvt) obs.PrintHepEvt(gFakeRunNumber);

    // Generate MCParts from observables
    simb::MCParticle p1(0,obs.pdg1,"primary",-1,obs.mass1,1);
    simb::MCParticle p2(1,obs.pdg2,"primary",-1,obs.mass2,1);
    TLorentzVector pos(obs.xPos, obs.yPos, obs.zPos, obs.time);// time needs to be in ns to match GENIE, etc
    TLorentzVector mom1(obs.P1[0],obs.P1[1],obs.P1[2],obs.E1);
    TLorentzVector mom2(obs.P2[0],obs.P2[1],obs.P2[2],obs.E2);
    p1.AddTrajectoryPoint(pos,mom1);
    p2.AddTrajectoryPoint(pos,mom2);
    truth.Add(p1);
    truth.Add(p2);
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));

    gFakeRunNumber++;
    return;
  } // END function produce

  // Compress fcl settings to a struct in order to make it easier
  // to pass them to other functions
  void HSNGen::CompressSettings(Settings &set)
  {
    set.printHepEvt = fPrintHepEvt;
    set.sterileMass = fSterileMass;
    set.decayChannel = fDecayChannel;
    set.fluxFile = fFluxFile;
    set.distance = fDistance;
    set.globalTimeOffset = fGlobalTimeOffset;
    set.beamWindow = fBeamWindow;
    set.boundariesX = fBoundariesX;
    set.boundariesY = fBoundariesY;
    set.boundariesZ = fBoundariesZ;
    return;
  }

  DEFINE_ART_MODULE(HSNGen)

} // END namespace hsngen

#endif
