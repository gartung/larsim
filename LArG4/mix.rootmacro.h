#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"


Double_t mixLaga(Double_t *x, Double_t *par) {

     // Control constants
     Double_t np = 100.0;      // number of convolution steps
     Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

     // Variables
     Double_t fmixGaus;
     Double_t sum = 0.0;
     Double_t xlow,xupp;
     Double_t step;

     Double_t shift=par[2];  // zero position

     // Range of convolution integral
     xlow = x[0] - shift - sc * par[1];
     xupp = x[0] - shift + sc * par[1];

     step = (xupp-xlow) / np;


     if(x[0] >= shift - 3*par[1])
     	fmixGaus = par[0]*TMath::Gaus( x[0],shift, par[1]) + par[3]*TMath::Landau(x[0],par[5],par[4]) ; 
	 else 
	   	fmixGaus= 0;

 	 sum+=fmixGaus ;


     return sum; 

}









