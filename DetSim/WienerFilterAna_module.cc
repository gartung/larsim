////////////////////////////////////////////////////////////////////////
// $Id: WienerFilterAna.cxx,
//
// WienerFilterAna class designed to calculate the optimum filter for an event
//  (based strongly on CalWireAna)
//  andrzej.szelc@yale.edu
//
//
////////////////////////////////////////////////////////////////////////
#ifndef WIENERFILTERANA_H
#define WIENERFILTERANA_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>
#include <string>
#include <cmath>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Simulation/sim.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Wire.h"
#include "Utilities/LArFFT.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/DetectorProperties.h"

// ROOT includes
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>

class TH2F;
class TH1F;
class TH1;

namespace geo { class Geometry; }

namespace detsim {

  /// Base class for creation of raw signals on wires. 
  class WienerFilterAna : public art::EDAnalyzer {
    
  public:
        
    explicit WienerFilterAna(fhicl::ParameterSet const& pset); 
    virtual ~WienerFilterAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void endJob();

  private:

    std::string            fDetSimModuleLabel; //< name of module that produced the digits

    // Members to hold signal information for buiding filters
    // C-style arrays used to make data members geometry agnostic
    unsigned int fNCryostats;
    unsigned int fNTPC;
    unsigned int fNPlanes;
    void                   buildDataStructures( unsigned int ncryo, unsigned int ntpc, unsigned int nplane );

    //TH1F*                  fCnoise[10][10][5];
    //TH1F*                  fCsignal[10][10][5];

    TH1F**		   fCnoise_av;
    TH1F**		   fCsignal_av;
    TH1F**		   fFilter_av;
    TH1F**		   fSmoothedFilter_av;
    TH1F*                  getAveNoiseHist( int icryo, int itpc, int iplane );  ///< preferred method for accessing above member
    void                   setAveNoiseHist( int icryo, int itpc, int iplane, TH1F* hist );  ///< preferred method for accessing above member
    TH1F*                  getAveSignalHist( int icryo, int itpc, int iplane ); ///< preferred method for accessing above member
    void                   setAveSignalHist( int icryo, int itpc, int iplane, TH1F* hist );  ///< preferred method for accessing above member
    TH1F*                  getAveFilterHist( int icryo, int itpc, int iplane );
    void                   setAveFilterHist( int icryo, int itpc, int iplane, TH1F* hist );  ///< preferred method for accessing above member
    TH1F*                  getSmoothedFilterHist( int icryo, int itpc, int iplane );
    void                   setSmoothedFilterHist( int icryo, int itpc, int iplane, TH1F* hist );  ///< preferred method for accessing above member


    TH1F*                  ff;
    TH1F*                  hh;
    int 		   fNBins;	
    int                    fSmoothingHalfwidthBins;

    std::vector< int >     fNoiseWireRange_cathode;
    std::vector< int >     fSignalWireRange_cathode;

    std::vector< int >     fNoiseWireRange_inductionU;
    std::vector< int >     fSignalWireRange_inductionU;

    std::vector< int >     fNoiseWireRange_inductionV;
    std::vector< int >     fSignalWireRange_inductionV;
    
    void removePedestal( TH1F* h, int start, int end );
    void getSamplingRange( int cryo, int tpc, int plane, unsigned int& signal_start, unsigned int& signal_end, unsigned int& noise_start, unsigned int& noise_end );    

  }; // class WienerFilterAna

} // End caldata namespace.

#endif // WIENERFILTERANA_H

namespace detsim{

  //-------------------------------------------------
  WienerFilterAna::WienerFilterAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset),
      fDetSimModuleLabel(pset.get< std::string >("DetSimModuleLabel")),
      fSmoothingHalfwidthBins(  std::abs(pset.get< int >( "SmoothingHalfwidthBins" )) ),
      fNoiseWireRange_cathode( pset.get< std::vector<int> >( "CathNoiseWireRange" ) ),
      fSignalWireRange_cathode( pset.get< std::vector<int> >( "CathSignalWireRange" ) ),
      fNoiseWireRange_inductionU( pset.get< std::vector<int> >( "IndUNoiseWireRange" ) ),
      fSignalWireRange_inductionU( pset.get< std::vector<int> >( "IndUSignalWireRange" ) ),
      fNoiseWireRange_inductionV( pset.get< std::vector<int> >( "IndVNoiseWireRange" ) ),
      fSignalWireRange_inductionV( pset.get< std::vector<int> >( "IndVSignalWireRange" ) )
  {
    // the fcl interface needs work. too microboone specific.
    hh = NULL;
    ff = NULL;
  }

  //-------------------------------------------------
  WienerFilterAna::~WienerFilterAna()
  {
    // if ( hh )
    //   delete hh;
    // if ( ff )
    //   delete ff;
  }

  //-------------------------------------------------
  void WienerFilterAna::beginJob() 
  {
    art::ServiceHandle<geo::Geometry> geo;
    buildDataStructures( geo->Ncryostats(), geo->NTPC(), geo->Nplanes() );
    //ff=tfs->make<TH1>(Form("fftwaveform"),Form("fftwaveform"),4096,0,4096);
    //hh=tfs->make<TH1F>(Form("waveform"),Form("waveform"),fNTicks,0,fNTicks);
    return;

  }

  void WienerFilterAna::buildDataStructures( unsigned int ncryo, unsigned int ntpc, unsigned int nplane ) {

    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<util::LArFFT> fFFT;
    int fNTicks = fFFT->FFTSize();
    fNBins=fNTicks/2+1;
    art::ServiceHandle<util::DetectorProperties> detp;
    double samprate=detp->SamplingRate();
    double sampfreq=1./samprate *1e6; // in kHz                                                                                                                                   
    fNCryostats = ncryo;
    fNTPC = ntpc;
    fNPlanes = nplane;

    fCnoise_av = new TH1F*[ fNPlanes*fNTPC*fNCryostats ];
    fCsignal_av = new TH1F*[ fNPlanes*fNTPC*fNCryostats ];
    fFilter_av = new TH1F*[ fNPlanes*fNTPC*fNCryostats ];
    fSmoothedFilter_av = new TH1F*[ fNPlanes*fNTPC*fNCryostats ];

    memset( fCnoise_av, 0, sizeof(TH1F*)*fNPlanes*fNTPC*fNCryostats );
    memset( fCsignal_av, 0, sizeof(TH1F*)*fNPlanes*fNTPC*fNCryostats );
    memset( fFilter_av, 0, sizeof(TH1F*)*fNPlanes*fNTPC*fNCryostats );
    memset( fSmoothedFilter_av, 0, sizeof(TH1F*)*fNPlanes*fNTPC*fNCryostats );
    
    for(unsigned int icstat=0;icstat<fNCryostats;icstat++){
      for(unsigned int itpc=0;itpc<fNTPC;itpc++){
	for(unsigned int iplane=0;iplane<fNPlanes;iplane++) {
	  // TH1F* noisehist = getAveNoiseHist( icstat, itpc, iplane );
	  // TH1F* signalhist = getAveSignalHist( icstat, itpc, iplane );
	  // TH1F* filterhist = getAveFilterHist( icstat, itpc, iplane );
	  // TH1F* smoothhist = getSmoothedFilterHist( icstat, itpc, iplane );

	  // std::cout << "prealloc noishist=" << noisehist << std::endl;
	    
	  TH1F* noisehist = tfs->make<TH1F>(Form("fft_noise_%d_%d_%d_av",icstat,itpc,iplane),
					    Form("fft of ave. noise plane:%d_%d_%d;freq (Hz);Magnitude",icstat,itpc,iplane),
					    fNBins,0,sampfreq/2);
	  TH1F* signalhist = tfs->make<TH1F>(Form("fft_signal_%d_%d_%d_av",icstat,itpc,iplane),
					     Form("fft of ave. signal plane:%d_%d_%d;freq (Hz);Magnitude",icstat,itpc,iplane),
					     fNBins,0,sampfreq/2);
	  TH1F* filterhist = tfs->make<TH1F>(Form("fft_filter_%d_%d_%d_av",icstat,itpc,iplane),
					     Form("fft of ave. filter plane:%d_%d_%d;freq (Hz);Magnitude",icstat,itpc,iplane),
					     fNBins,0,sampfreq/2);
	  TH1F* smoothhist = tfs->make<TH1F>(Form("fft_smoothed_filter_%d_%d_%d_av",icstat,itpc,iplane),
					     Form("fft of smoothed filter plane:%d_%d_%d;freq (Hz);Magnitude",icstat,itpc,iplane),
					     fNBins,0,sampfreq/2);

	  setAveNoiseHist( icstat, itpc, iplane, noisehist );
	  setAveSignalHist( icstat, itpc, iplane, signalhist );
	  setAveFilterHist( icstat, itpc, iplane, filterhist );
	  setSmoothedFilterHist( icstat, itpc, iplane, smoothhist );

	  // std::cout << "test noise " << icstat << " " << itpc << " " << iplane << ": " << noisehist << " vs " << getAveNoiseHist( icstat, itpc, iplane ) << std::endl;
	  // std::cout << "test signal: " << signalhist << " vs " << getAveSignalHist( icstat, itpc, iplane ) << std::endl;
	  // std::cout << "test filter: " << filterhist << " vs " << getAveFilterHist( icstat, itpc, iplane ) << std::endl;
	  // std::cout << "test smoothedFilter: " << smoothhist << " vs " << getSmoothedFilterHist( icstat, itpc, iplane ) << std::endl;

	}//end of plane loop
      }//end of tpc loop
    } // end of cryo loop

    // for(unsigned int icstat=0;icstat<fNCryostats;icstat++){
    //   for(unsigned int itpc=0;itpc<fNTPC;itpc++){
    //     for(unsigned int iplane=0;iplane<fNPlanes;iplane++) {
    // 	  std::cout << "test noise " << icstat << " " << itpc << " " << iplane << ":  " << getAveNoiseHist( icstat, itpc, iplane ) << std::endl;
    // 	  std::cout << "test signal: " << getAveSignalHist( icstat, itpc, iplane ) << std::endl;
    // 	  std::cout << "test filter: " << getAveFilterHist( icstat, itpc, iplane ) << std::endl;
    // 	  std::cout << "test smoothedFilter: " << getSmoothedFilterHist( icstat, itpc, iplane ) << std::endl;
    // 	}
    //   }
    // }
    // std::cout << "Allocated data structures." << std::endl;
    
  }
  
  TH1F* WienerFilterAna::getAveNoiseHist( int icryo, int itpc, int iplane ) {
    TH1F* ptr =  *(fCnoise_av + fNCryostats*icryo + fNTPC*itpc + iplane); // dereference send us to address of pointer
    return ptr;
  }

  void WienerFilterAna::setAveNoiseHist( int icryo, int itpc, int iplane, TH1F* hist ) {
    *(fCnoise_av + fNCryostats*icryo + fNTPC*itpc + iplane) = hist; // dereference send us to address of pointer
  }

  TH1F* WienerFilterAna::getAveSignalHist( int icryo, int itpc, int iplane ) {
    TH1F* ptr =  *(fCsignal_av + fNCryostats*icryo + fNTPC*itpc + iplane); // dereference send us to address of pointer
    return ptr;
  }

  void WienerFilterAna::setAveSignalHist( int icryo, int itpc, int iplane, TH1F* hist ) {
    *(fCsignal_av + fNCryostats*icryo + fNTPC*itpc + iplane) = hist; // dereference send us to address of pointer
  }


  TH1F* WienerFilterAna::getAveFilterHist( int icryo, int itpc, int iplane ) {
    TH1F* ptr =  *(fFilter_av + fNCryostats*icryo + fNTPC*itpc + iplane); // dereference send us to address of pointer
    return ptr;
  }

  void WienerFilterAna::setAveFilterHist( int icryo, int itpc, int iplane, TH1F* hist ) {
    *(fFilter_av + fNCryostats*icryo + fNTPC*itpc + iplane) = hist; // dereference send us to address of pointer
  }

  TH1F* WienerFilterAna::getSmoothedFilterHist( int icryo, int itpc, int iplane ) {
    TH1F* ptr =  *(fSmoothedFilter_av + fNCryostats*icryo + fNTPC*itpc + iplane); // dereference send us to address of pointer
    return ptr;
  }

  void WienerFilterAna::setSmoothedFilterHist( int icryo, int itpc, int iplane, TH1F* hist ) {
    *(fSmoothedFilter_av + fNCryostats*icryo + fNTPC*itpc + iplane) = hist; // dereference send us to address of pointer
  }


  //-------------------------------------------------
  void WienerFilterAna::endJob() 
  {
    
    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    // calculate filters 
    for(unsigned int icstat=0;icstat<fNCryostats;icstat++){
      for(unsigned int itpc=0;itpc<fNTPC;itpc++){
	for (unsigned int pp=0;pp<fNPlanes;pp++){

	  TH1F* signalhist = getAveSignalHist( icstat, itpc, pp );
	  TH1F* noisehist  = getAveNoiseHist( icstat, itpc, pp );
	  TH1F* filterhist = getAveFilterHist( icstat, itpc, pp );
	  TH1F* smoothhist = getSmoothedFilterHist( icstat, itpc, pp );
	  for(int ii=1;ii<signalhist->GetNbinsX();ii++) {
	    double diff = signalhist->GetBinContent(ii)*signalhist->GetBinContent(ii)-noisehist->GetBinContent(ii)*noisehist->GetBinContent(ii);
	    diff = ( diff>=0 ) ? diff : 0;
	    
	    if(signalhist->GetBinContent(ii)>0)
	      filterhist->SetBinContent(ii,(double)((diff)/(signalhist->GetBinContent(ii)*signalhist->GetBinContent(ii))) );
	    else
	      filterhist->SetBinContent(ii,0);
	  }

	  // smooth filter
	  for(int ii=1;ii<signalhist->GetNbinsX();ii++)
	    {
                int start=ii-fSmoothingHalfwidthBins;
                int end = ii+fSmoothingHalfwidthBins;
                if ( start<1 )
                  start = 1;
                if ( end>signalhist->GetNbinsX() )
                  end = signalhist->GetNbinsX();
                double width = double(end-start);
                for (int jj=start; jj<end; jj++)
                  smoothhist->AddBinContent(ii,filterhist->GetBinContent(jj)/width);
	    }//end of loop over Csignal                                                                                                                                                          
	  // remove high-f pedestal
	  int nbackpedbins = 0;
	  double backped = 0;
	  for (int ii=int(0.9*smoothhist->GetNbinsX()); ii<=smoothhist->GetNbinsX(); ii++) {
	    backped += smoothhist->GetBinContent(ii);
	    nbackpedbins++;
	  }
	  backped /= double(nbackpedbins);
	  for (int ii=1; ii<=smoothhist->GetNbinsX(); ii++) {
	    smoothhist->AddBinContent( ii, -backped );
	    if ( smoothhist->GetBinContent( ii )<0.005 )
	      smoothhist->SetBinContent( ii, 0 );
	  }
	  
	  // normalize filter: peak at 1.0
	  double smooth_max = smoothhist->GetMaximum();
	  smoothhist->Scale(1.0/smooth_max);
	  
	  double raw_max = filterhist->GetMaximum();
	  filterhist->Scale(1.0/raw_max);
	  
	  // for smoothed filter. force 0 frequency to have zero power
	  for ( int ii=1; ii<smoothhist->GetNbinsX(); ii++ ) {
	    double w = smoothhist->GetBinLowEdge(ii);
	    double e = exp( 0.05*w )-1.0;
	    double f = smoothhist->GetBinContent(ii);
	    if ( e<f )
	      smoothhist->SetBinContent(ii,e);
	  }
	  // take no chances
	  smoothhist->SetBinContent(1,0);
	  smooth_max = smoothhist->GetMaximum();
	  smoothhist->Scale(1.0/smooth_max);
	  
	} // end loop on planes
      } // end loop on TPCs
    } // end loop on cryostats
    
  }

  //-------------------------------------------------
  void WienerFilterAna::analyze(const art::Event& evt)
  {

    // loop over the raw digits and get the adc vector for each, then compress it and uncompress it
    
    art::Handle< std::vector<raw::RawDigit> > rdHandle;
    evt.getByLabel(fDetSimModuleLabel,rdHandle);
    mf::LogInfo("WienerFilterMicroBooNE") << " readout  Wiener " << rdHandle->size() <<  std::endl;
    // return;
    if (!rdHandle->size())  return;
    mf::LogInfo("WienerFilterMicroBooNE") << "WienerFilterMicroBooNE:: rdHandle size is " << rdHandle->size();
 
    // Read in the digit List object(s). 

    // Use the handle to get a particular (0th) element of collection.
    art::PtrVector<raw::RawDigit> rdvec;
    for(unsigned int i = 0; i < rdHandle->size(); ++i){
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
     // std::cout << " i, rdvec: "<<i <<" " << r->ADC(0) << " "<< rdvec[i]->ADC(0)<< std::endl;
    }

    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArFFT> fft;
    unsigned int noise_start, noise_end, signal_start, signal_end;
    std::vector<double> adc(fft->FFTSize());
    if ( !hh )
      hh = new TH1F( "hwienerviretemp", "temp", fft->FFTSize(), 0, fft->FFTSize() );

    for(unsigned int rd = 0; rd < rdvec.size(); ++rd){
      
      // make decision on wire before filling adcs
      geo::WireID wireid= geom->ChannelToWire(rdvec[rd]->Channel())[0];
      unsigned int plane = wireid.Plane;   /// \todo  Need to change hardcoded values to an automatic determination of noise vs. signal 
      unsigned int wire = wireid.Wire;
      unsigned int cstat = wireid.Cryostat;
      unsigned int tpc = wireid.TPC;

      getSamplingRange( cstat, tpc, plane, signal_start, signal_end, noise_start, noise_end );

      if ( wire>=noise_start && wire<=noise_end ) {
	
	for(unsigned int t = 1; t <= rdvec[rd]->Samples(); t++){
	  hh->SetBinContent(t,rdvec[rd]->ADC(t-1));
	}	
	removePedestal( hh, 0, int(rdvec[rd]->Samples()*0.05) );
	
	if ( !ff )
	  ff=(TH1F*)hh->FFT(NULL,"MAG M");
	else {
	  ff->Reset();
	  ff = (TH1F*)hh->FFT( ff, "MAG" );
	}
	
	TH1F* noisehist = getAveNoiseHist( cstat, tpc, plane );
	//std::cout << "noisehist " << cstat << " " << tpc << " " << plane << ": " << noisehist << std::endl;
	for(int ii=0;ii<fNBins;ii++)
	  noisehist->AddBinContent(ii,ff->GetBinContent(ii)); 
      }	
            
      if ( wire>=signal_start && wire<=signal_end ) {
	
        for(unsigned int t = 1; t <= rdvec[rd]->Samples(); t++){
          hh->SetBinContent(t,rdvec[rd]->ADC(t-1));
	}
	removePedestal( hh, 0, int(rdvec[rd]->Samples()*0.05) );

        if ( !ff )
          ff=(TH1F*)hh->FFT(NULL,"MAG M");
        else {
          ff->Reset();
          ff = (TH1F*)hh->FFT( ff, "MAG" );
	}

        TH1F* signalhist = getAveSignalHist( cstat, tpc, plane );
        for(int ii=0;ii<fNBins;ii++)
          signalhist->AddBinContent(ii,ff->GetBinContent(ii));
      }
      

    }//end loop over rawDigits
    
    return;
  }//end analyze method

  void WienerFilterAna::removePedestal( TH1F* h, int start, int end ) {
    double tot = 0.;
    for (int i=start; i<end; i++) {
      tot += h->GetBinContent( i+1 );
    }
    double ped = tot/std::fabs( double( end-start ) );
    for (int i=1; i<=h->GetNbinsX(); i++) {
      if ( h->GetBinContent(i)>0 )
	h->AddBinContent( i, -ped );
    }
  }
  
  void WienerFilterAna::getSamplingRange( int cryo, int tpc, int plane, 
					  unsigned int& signal_start, unsigned int& signal_end, unsigned int& noise_start, unsigned int& noise_end ) {
    // cryo and tpc does nothing!
    if ( cryo!=0 || tpc!=0 ) {
      std::cout << __PRETTY_FUNCTION__ << ": WARNING -- cryo and tpc indices do nothing now!" << std::endl;
    }
    std::vector<int>* noiseranges[3] = { &fNoiseWireRange_inductionU,
                                         &fNoiseWireRange_inductionV,
                                         &fNoiseWireRange_cathode };
    std::vector<int>* signalranges[3] = { &fSignalWireRange_inductionU,
                                          &fSignalWireRange_inductionV,
                                          &fSignalWireRange_cathode };
    
    if ( noiseranges[plane]->size()<2 || signalranges[plane]->size()<2 ) {
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << std::cout << "ERROR: Either noise or signal wire ranges incomplete" << std::endl;
    }
    
    noise_start = std::abs( noiseranges[plane]->at(0) );
    noise_end   = std::abs( noiseranges[plane]->at(1) );
    signal_start = std::abs( signalranges[plane]->at(0) );
    signal_end   = std::abs( signalranges[plane]->at(1) );
    if ( noise_start>noise_end ) {
      unsigned int c = noise_start;
      noise_start = noise_end;
      noise_end = c;
    }
    if ( signal_start>signal_end ) {
      unsigned int c = signal_start;
      signal_start = signal_end;
      signal_end = c;
    }
    
  }
  
}//end namespace

namespace detsim{

  DEFINE_ART_MODULE(WienerFilterAna)

}
