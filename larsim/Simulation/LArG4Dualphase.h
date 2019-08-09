////////////////////////////////////////////////////////////////////////
/// \file  LArG4Dualphase.h
/// \brief Store parameters for running LArG4 dual-phase simulation
///
/// \author andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <algorithm> // std::max()

#ifndef LArG4Dualphase_h
#define LArG4Dualphase_h 1

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

// ROOT
#include "TF1.h"

namespace sim {

  class LArG4Dualphase {
  public:
    LArG4Dualphase(fhicl::ParameterSet const& pset);
    ~LArG4Dualphase() {}

    void reconfigure(fhicl::ParameterSet const& pset);
    //void initGAr();

    double GainCornerLEMs()				      const { return fGainCornerLEMs;  }
    double GainCentralLEMs()				      const { return fGainCentralLEMs;  }
    double ExtractionField()				      const { return fExtractionField;  }
    double InductionField()				      const { return fInductionField;  }
    double MinNumberOfElClusterFast()		      	      const { return fMinNumberOfElClusterFast;  }
    double MinNumberOfElClusterSlow()		      	      const { return fMinNumberOfElClusterSlow;  }
    std::string SlowExtractionTimeConstantShape()       const { return fSlowExtractionTimeConstantShape;    }
    std::vector<double> SlowExtractionTimeConstantParams()    const { return fSlowExtractionTimeConstantParams;    }
    double SlowExtractionTimeConstantCorrectionFactor()       const { return fSlowExtractionTimeConstantCorrectionFactor;    }

    void simExtraction(int nIonizedElectrons, int &nElectrons, int &nClus, std::vector< double > &nElDiff);
    double getExtractionDelay(const int clusterIndex);

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

    //Extraction paramters
    double fastextractioncomponent=1.;
    double p0fast = 2.91754e-03;
    double p1fast = 7.74228e-04;
    double p2fast = 5.82340e-01;
    double p3fast = -2.42075e-01;
    double p4fast = 2.82637e-02;

    double slowextractioncomponent=1.;
    double amplitudeslow = 1.90787e-01;
    double meanslow = 1.17742e+00;
    double sigmaslow = 4.83614e-01;

    TF1* fSlowExtractionTimeConstantFunction;
    TF1 *fslowextractiondelay;

    double slowextractiontimeconstant=0.; //microseconds
    double aslow = 1.23091e+05;
    double bslow = -0.0612014;

    int nClusSlow;
    int nClusFast;

  };
}


DECLARE_ART_SERVICE(sim::LArG4Dualphase, GLOBAL)
#endif
