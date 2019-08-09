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
// Andrea Scarpelli, APC-Université de Paris, August 2019
////////////////////////////////////////////////////////////////////////

// Framework includes

#include "larsim/Simulation/LArG4Dualphase.h"

// ROOT
//#include "TF1.h"

namespace sim {

  //----------------------------------------------------------------------------
  LArG4Dualphase::LArG4Dualphase(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //----------------------------------------------------------------------------
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

    fastextractioncomponent = std::min(1.,p0fast + p1fast*fExtractionField + p2fast*pow(fExtractionField,2) + p3fast*pow(fExtractionField,3) + p4fast*pow(fExtractionField,4));

    slowextractioncomponent = std::min(1.,amplitudeslow*exp(-0.5*pow((fExtractionField-meanslow)/sigmaslow,2)));
    slowextractioncomponent += fastextractioncomponent;
    fastextractioncomponent = 0.;

    std::cout << "Fast extraction component: " << fastextractioncomponent << std::endl;
    std::cout << "Slow extraction component: " << slowextractioncomponent << std::endl;

    // slow extraction time constant FHiCL parameters
    fSlowExtractionTimeConstantFunction = new TF1("fSlowExtractionTimeConstantFunction", fSlowExtractionTimeConstantShape.c_str());
    for(unsigned int i=0; i<fSlowExtractionTimeConstantParams.size(); ++i)
    {
      fSlowExtractionTimeConstantFunction->SetParameter(i, fSlowExtractionTimeConstantParams[i]);
    }
    fSlowExtractionTimeConstantFunction->SetRange(0,100000);
    slowextractiontimeconstant = 1000*fSlowExtractionTimeConstantFunction->Eval(1000*fExtractionField);

    std::cout << "Slow extraction time constant (nanoseconds): " << slowextractiontimeconstant << std::endl;
    slowextractiontimeconstant *=fSlowExtractionTimeConstantCorrectionFactor;
    std::cout << "Slow extraction time constant after correction (nanoseconds): " << slowextractiontimeconstant << std::endl;

    fslowextractiondelay = new TF1("fslowextractiondelay","exp(-x/[0])",0,100000);
    fslowextractiondelay->SetParameter(0,slowextractiontimeconstant);

    return;
  }

//------------------------------------------------------------------------------
/*
void LArG4Dualphase::initGAr()
{

}
*/

//------------------------------------------------------------------------------

  void LArG4Dualphase::simExtraction(int nIonizedElectrons, int &nElectrons, int &nClus, std::vector< double > &nElDiff)
  {
    nClusFast=0;
    nClusSlow=0;

    double totalextractionefficiency = std::min(1.,fastextractioncomponent+slowextractioncomponent);

    nElectrons   = nIonizedElectrons * totalextractionefficiency;
    const double nElectronsFast = nIonizedElectrons * std::min(1.,fastextractioncomponent);
    const double nElectronsSlow = nIonizedElectrons * std::min(1.,slowextractioncomponent);

    //Cluster for fast component
    double electronclsizefast = nElectronsFast / fMinNumberOfElClusterFast;
    nClusFast = fMinNumberOfElClusterFast;
    if (electronclsizefast < 1.0)
    {
      electronclsizefast = 1.0;
      nClusFast = (int) std::ceil(nElectronsFast / electronclsizefast);
    }

    std::vector< double > nElDiffFast(nClusFast, electronclsizefast);
    // fix the number of electrons in the last cluster, that has smaller size
    if(nClusFast>0) nElDiffFast.back() = nElectronsFast - (nClusFast-1)*electronclsizefast;

    //Cluster for slow component
    double electronclsizeslow = nElectronsSlow / fMinNumberOfElClusterSlow;
    nClusSlow = fMinNumberOfElClusterSlow;
    if (electronclsizeslow < 1.0)
    {
      electronclsizeslow = 1.0;
      nClusSlow = (int) std::ceil(nElectronsSlow / electronclsizeslow);
    }

    // Compute arrays of values as quickly as possible.
    std::vector< double > nElDiffSlow(nClusSlow, electronclsizeslow);
    // fix the number of electrons in the last cluster, that has smaller size
    if(nClusSlow>0) nElDiffSlow.back() = nElectronsSlow - (nClusSlow-1)*electronclsizeslow;

    nClus = nClusFast + nClusSlow;

    nElDiff.reserve( nElDiffFast.size() + nElDiffSlow.size() ); // preallocate memory
    nElDiff.insert( nElDiff.end(), nElDiffFast.begin(), nElDiffFast.end() );
    nElDiff.insert( nElDiff.end(), nElDiffSlow.begin(), nElDiffSlow.end() );

    return;
  }

//------------------------------------------------------------------------------

  double LArG4Dualphase::getExtractionDelay(const int clusterIndex)
  {
    double delay=0;

    if( clusterIndex > (nClusFast+nClusSlow)  )
      delay = 1;//fslowextractiondelay->GetRandom();
    else
      delay=0;

    return delay;
  }

//------------------------------------------------------------------------------

  /*
  void LArG4Dualphase::driftInGAr(  )
  {

  }
  */

}

DEFINE_ART_SERVICE(sim::LArG4Dualphase)
