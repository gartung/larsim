#ifndef _WEIGHTCALC_H_
#define _WEIGHTCALC_H_
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TString.h"
#include "TTree.h"
#include "TMatrixD.h"
#include <string>
#include <map>
#include <vector>

//weight calc base
namespace evwgh {
  typedef std::map<std::string, std::vector<double> > WeightMap_t;

  class WeightCalc
  {    
  public:
    virtual void                Configure(fhicl::ParameterSet const& pset,
                                          CLHEP::HepRandomEngine&) = 0;
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e) = 0; 
    void                        SetName(std::string name) {fName=name;}
    std::string                 GetName() {return fName;}

    void SetupParameterTree(std::string treeName, std::string title="");
    
    /**
     * @brief Applies Gaussian smearing to a set of data
     * @param centralValues the values to be smeared
     * @param inputCovarianceMatrix covariance matrix for smearing
     * @param n_multisims number of sets of smeared values to be produced
     * @return a set of n_multisims value sets smeared from the central value
     *
     * If centralValues is of dimension N, inputCovarianceMatrix needs to be NxN,
     * and each of the returned data sets will be also of dimension N.
     */
    static std::vector<std::vector<double> > MultiGaussianSmearing(
			std::vector<double> const& centralValues,
			std::vector< std::vector<double>> const& inputCovarianceMatrix,
			int n_multisims, CLHEP::RandGaussQ& GaussRandom);

    static std::vector<double>              MultiGaussianSmearing(
			                        std::vector<double> const& centralValue, 
                                                TMatrixD* const& inputCovarianceMatrix, 
                                                std::vector<double> rand);

    static std::vector<double>              MultiGaussianSmearing(
			                        std::vector<double> const& centralValue, 
                                                TMatrixD* const& LowerTriangleCovarianceMatrix, 
						bool isDecomposed,
                                                std::vector<double> rand);

  protected:
    TTree* fParameterTree;
    TString* fParameterNameBr;
    std::vector<float>* fParameterValBr;

  private:
    std::string fName;
  };
  
}

#endif // _WEIGHTCALC_H_
