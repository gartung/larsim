////////////////////////////////////////////////////////////////////////
/// \file  LArG4Dualphase.h
/// \brief Store parameters for running LArG4 dual-phase simulation
///
/// \author andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>

#ifndef LArG4Dualphase_h
#define LArG4Dualphase_h 1

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

namespace sim {

  class LArG4Dualphase {
  public:
    LArG4Parameters(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~LArG4Parameters() {}

    void reconfigure(fhicl::ParameterSet const& pset);

    double GainCornerLEMs()				      const { return fGainCornerLEMs;  }
    double GainCentralLEMs()				      const { return fGainCentralLEMs;  }
    double ExtractionField()				      const { return fExtractionField;  }
    double InductionField()				      const { return fInductionField;  }
    double MinNumberOfElClusterFast()		      	      const { return fMinNumberOfElClusterFast;  }
    double MinNumberOfElClusterSlow()		      	      const { return fMinNumberOfElClusterSlow;  }
    std::string SlowExtractionTimeConstantShape()       const { return fSlowExtractionTimeConstantShape;    }
    std::vector<double> SlowExtractionTimeConstantParams()    const { return fSlowExtractionTimeConstantParams;    }
    double SlowExtractionTimeConstantCorrectionFactor()       const { return fSlowExtractionTimeConstantCorrectionFactor;    }

  private:
    
    double fGainCornerLEMs;
    double fGainCentralLEMs;
    double fExtractionField;
    double fInductionField;
    double fMinNumberOfElClusterFast;
    double fMinNumberOfElClusterSlow;

    std::string fSlowExtractionTimeConstantShape;
    std::vector<double> fSlowExtractionTimeConstantParams;

    double fSlowExtractionTimeConstantCorrectionFactor;
  };
}


DECLARE_ART_SERVICE(sim::LArG4Dualphase, LEGACY)
#endif
