////////////////////////////////////////////////////////////////////////
/// \file  LArG4Dualphase_service.cc
/// \brief Store parameters for running LArG4 dual-phase simulation
///
/// \author andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Andrea Scarpelli APC August 2019


// Framework includes

#include "larsim/Simulation/LArG4Dualphase.h"

namespace sim {

  //--------------------------------------------------------------------------
  LArG4Dualphase::LArG4Dualphase(fhicl::ParameterSet const& pset, art::ActivityRegistry& /* reg */)
  {
    this->reconfigure(pset);
  }

  //--------------------------------------------------------------------------
  void LArG4Dualphase::reconfigure(fhicl::ParameterSet const& pset)
  {
    fGainCornerLEMs = pset.get<double>("GainCornerLEMs");
    fGainCentralLEMs = pset.get<double>("GainCentralLEMs");
    fExtractionField = pset.get<double>("ExtractionField");
    fInductionField = pset.get<double>("InductionField");
    fMinNumberOfElClusterFast = pset.get<double>("MinNumberOfElClusterFast");
    fMinNumberOfElClusterSlow = pset.get<double>("MinNumberOfElClusterSlow");
    fSlowExtractionTimeConstantShape = pset.get< std::string  >("SlowExtractionTimeConstantShape"    );
    fSlowExtractionTimeConstantParams = pset.get< std::vector<double>  >("SlowExtractionTimeConstantParams"    );
    fSlowExtractionTimeConstantCorrectionFactor = pset.get< double  >("SlowExtractionTimeConstantCorrectionFactor"    );

    return;
  }
}


namespace sim {

  DEFINE_ART_SERVICE(LArG4Dualphase)

}
