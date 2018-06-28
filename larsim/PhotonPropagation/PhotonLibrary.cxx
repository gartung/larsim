
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "larcore/Geometry/Geometry.h"

#include "larsim/PhotonPropagation/PhotonLibrary.h"
#include "larsim/Simulation/PhotonVoxels.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"


namespace phot{
  
  std::string const PhotonLibrary::OpChannelBranchName = "OpChannel";
  
  //------------------------------------------------------------

  PhotonLibrary::PhotonLibrary()
  {
    fLookupTable.clear();
    fReflLookupTable.clear();
    fReflTLookupTable.clear();
    fTimingParLookupTable.clear();
  }

  
  //------------------------------------------------------------  

  PhotonLibrary::~PhotonLibrary()
  {
    fLookupTable.clear();
    fReflLookupTable.clear();
    fReflTLookupTable.clear();
    fTimingParLookupTable.clear();
  }
  
  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile(std::string LibraryFile, bool storeReflected, bool storeReflT0, size_t storeTiming)
  {
    mf::LogInfo("PhotonLibrary") << "Writing photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel          = 0;
    Int_t     OpChannel      = 0;
    Float_t   Visibility     = 0;
    Float_t   ReflVisibility = 0;
    Float_t   ReflTfirst     = 0;
    Float_t   *timing_par = nullptr;

    tt->Branch("Voxel",      &Voxel,      "Voxel/I");
    tt->Branch(OpChannelBranchName.c_str(),  &OpChannel,  (OpChannelBranchName + "/I").c_str());
    tt->Branch("Visibility", &Visibility, "Visibility/F");

    if(storeTiming!=0)
    {
      if (!hasTiming()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflected set true
        throw cet::exception("PhotonLibrary")
          << "StoreLibraryToFile() requested to store the time propagation distribution parameters, which was not simulated.";
      }
      tt->Branch("timing_par", timing_par, Form("timing_par[%i]/F",size_t2int(storeTiming)));
       if (fLookupTable.size() != fTimingParLookupTable.size())
        throw cet::exception(" Photon Library ") << "Time propagation lookup table is different size than Direct table \n"
                                                   << "this should not be happening. ";
    }

    if(storeReflected)
    {
      if (!hasReflected()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflected set true
        throw cet::exception("PhotonLibrary")
          << "StoreLibraryToFile() requested to store reflected light, which was not simulated.";
      }
      tt->Branch("ReflVisibility", &ReflVisibility, "ReflVisibility/F");
      if (fLookupTable.size() != fReflLookupTable.size())
          throw cet::exception(" Photon Library ") << "Reflected light lookup table is different size than Direct table \n"
                                                   << "this should not be happening. ";
    }
    if(storeReflT0) {
      if (!hasReflectedT0()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflectedT0 set true
        throw cet::exception("PhotonLibrary")
          << "StoreLibraryToFile() requested to store reflected light timing, which was not simulated.";
      }
      tt->Branch("ReflTfirst", &ReflTfirst, "ReflTfirst/F");
    }
    for(size_t ivox=0; ivox!= fNVoxels; ++ivox)
    {
      for(size_t ichan=0; ichan!= fNOpChannels; ++ichan)
      {
        Visibility = uncheckedAccess(ivox, ichan);
        if(storeReflected)
          ReflVisibility = uncheckedAccessRefl(ivox, ichan);
        if(storeReflT0)
          ReflTfirst = uncheckedAccessReflT(ivox, ichan);
	if(storeTiming!=0)
	{
	  for (size_t i =0; i<storeTiming; i++) timing_par[i] = uncheckedAccessTimingPar(ivox, ichan, i);
	  
	}
        if (Visibility > 0 || ReflVisibility > 0)
        {
          Voxel      = ivox;
          OpChannel  = ichan;
          // visibility(ies) is(are) already set
          tt->Fill();
        }
      }	
    }
  }


  //------------------------------------------------------------

  void PhotonLibrary::CreateEmptyLibrary(
    size_t NVoxels, size_t NOpChannels,
    bool storeReflected /* = false */,
    bool storeReflT0 /* = false */, 
    size_t storeTiming /* = false */
  ) {
    fLookupTable.clear();
    fReflLookupTable.clear();
    fReflTLookupTable.clear();
    fTimingParLookupTable.clear();

    fNVoxels     = NVoxels;
    fNOpChannels = NOpChannels;

    fLookupTable.resize(LibrarySize(), 0.);
    fHasReflected = storeReflected;
    if (storeReflected) fReflLookupTable.resize(LibrarySize(), 0.);
    fHasReflectedT0 = storeReflT0;
    if (storeReflT0) fReflTLookupTable.resize(LibrarySize(), 0.);
    fHasTiming = storeTiming;
    if (storeTiming!=0)
    {
	fTimingParLookupTable.resize(LibrarySize());
    }	
  }


  //------------------------------------------------------------

  void PhotonLibrary::LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels, bool getReflected, bool getReflT0, size_t getTiming)
  {
    fLookupTable.clear();
    fReflLookupTable.clear();
    fReflTLookupTable.clear();
    fTimingParLookupTable.clear();
    fHasTiming = getTiming;
    mf::LogInfo("PhotonLibrary") << "Reading photon library from input file: " << LibraryFile.c_str()<<std::endl;

    TFile *f = nullptr;
    TTree *tt = nullptr;
      
    try
      {
	f  =  TFile::Open(LibraryFile.c_str());
	tt =  (TTree*)f->Get("PhotonLibraryData");
        if (!tt) { // Library not in the top directory
            TKey *key = f->FindKeyAny("PhotonLibraryData");
            if (key) 
                tt = (TTree*)key->ReadObj();
            else {
                mf::LogError("PhotonLibrary") << "PhotonLibraryData not found in file" <<LibraryFile;
            }
        }
      }
    catch(...)
      {
	mf::LogError("PhotonLibrary") << "Error in ttree load, reading photon library: " << LibraryFile.c_str()<<std::endl;
      }

    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
    Float_t   ReflVisibility;
    Float_t   ReflTfirst;
//    Double_t   *timing_par = nullptr;
    std::vector<Float_t>   timing_par;
//Float_t *timing_par = new Float_t[NPAR];

    tt->SetBranchAddress("Voxel",      &Voxel);
    tt->SetBranchAddress("OpChannel",  &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);

    if(fHasTiming!=0)
    {
      timing_par.resize(getTiming);
      tt->SetBranchAddress("timing_par", timing_par.data());
    }



    fHasReflected = getReflected;
    if(getReflected)
      tt->SetBranchAddress("ReflVisibility", &ReflVisibility);
    fHasReflectedT0 = getReflT0;
    if(getReflT0)
      tt->SetBranchAddress("ReflTfirst", &ReflTfirst);
    
    
    fNVoxels     = NVoxels;
    fNOpChannels = PhotonLibrary::ExtractNOpChannels(tt); // EXPENSIVE!!!
    
    fLookupTable.resize(LibrarySize(), 0.);
    if(fHasTiming!=0)
    {
      fTimingParLookupTable.resize(LibrarySize());
      for(size_t k=0;k<LibrarySize();k++) fTimingParLookupTable[k].resize(getTiming,0);
    }
    if(fHasReflected)
      fReflLookupTable.resize(LibrarySize(), 0.);
    if(fHasReflectedT0)
      fReflTLookupTable.resize(LibrarySize(), 0.);
    
    size_t NEntries = tt->GetEntries();

    for(size_t i=0; i!=NEntries; ++i) {
      tt->GetEntry(i);
      // Set the visibility at this optical channel
      uncheckedAccess(Voxel, OpChannel) = Visibility;

      if(fHasReflected)
	uncheckedAccessRefl(Voxel, OpChannel) = ReflVisibility;
      if(fHasReflectedT0)
	uncheckedAccessReflT(Voxel, OpChannel) = ReflTfirst; 
      if(fHasTiming!=0)
      {
//	tt->Draw("timing_par","","goff",1,i);
//	timing_par=tt->GetV1();
	for (size_t k=0;k<fHasTiming;k++){ uncheckedAccessTimingPar(Voxel, OpChannel,k) = timing_par[k];}
      }
    } // for entries
    
    mf::LogInfo("PhotonLibrary") <<"Photon lookup table size : "<<  NVoxels << " voxels,  " << fNOpChannels<<" channels";

    try
      {
	f->Close();
      }
    catch(...)
      {
	mf::LogError("PhotonLibrary") << "Error in closing file : " << LibraryFile.c_str()<<std::endl;
      }

  }

  //----------------------------------------------------

  float PhotonLibrary::GetCount(size_t Voxel, size_t OpChannel) const
  { 
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;   
    else
      return uncheckedAccess(Voxel, OpChannel); 
  }
  //----------------------------------------------------

  float PhotonLibrary::GetTimingPar(size_t Voxel, size_t OpChannel, size_t parnum) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessTimingPar(Voxel, OpChannel, parnum);
  }  //----------------------------------------------------

  float PhotonLibrary::GetReflCount(size_t Voxel, size_t OpChannel) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessRefl(Voxel, OpChannel);
  }
  //----------------------------------------------------

  float PhotonLibrary::GetReflT0(size_t Voxel, size_t OpChannel) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessReflT(Voxel, OpChannel);
  }

  //----------------------------------------------------

  void PhotonLibrary::SetCount(size_t Voxel, size_t OpChannel, float Count) 
  { 
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range"; 
    else
      uncheckedAccess(Voxel, OpChannel) = Count; 
  }
  //----------------------------------------------------

  void PhotonLibrary::SetTimingPar(size_t Voxel, size_t OpChannel, float Count, size_t parnum) 
  { 
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set timing t0 count in voxel " << Voxel<<" which is out of range"; 
    else
      uncheckedAccessTimingPar(Voxel, OpChannel, parnum) = Count; 
  }
  //----------------------------------------------------

  void PhotonLibrary::SetReflCount(size_t Voxel, size_t OpChannel, float Count)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range";
    else
      uncheckedAccessRefl(Voxel, OpChannel) = Count;
  }
  //----------------------------------------------------

  void PhotonLibrary::SetReflT0(size_t Voxel, size_t OpChannel, float Count)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range";
    else
      uncheckedAccessReflT(Voxel, OpChannel) = Count;
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetCounts(size_t Voxel) const
  { 
    if (Voxel >= fNVoxels) return nullptr;
    else return fLookupTable.data() + uncheckedIndex(Voxel, 0);
  }

  //----------------------------------------------------

  const std::vector<float>* PhotonLibrary::GetTimingPars(size_t Voxel) const
  { 
    if (Voxel >= fNVoxels) return nullptr;
    else return fTimingParLookupTable.data() + uncheckedIndex(Voxel, 0);
  }
  //----------------------------------------------------

  float const* PhotonLibrary::GetReflCounts(size_t Voxel) const
  {
    if (Voxel >= fNVoxels) return nullptr;
    else return fReflLookupTable.data() + uncheckedIndex(Voxel, 0);
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetReflT0s(size_t Voxel) const
  {
    if (Voxel >= fNVoxels) return nullptr;
    else return fReflTLookupTable.data() + uncheckedIndex(Voxel, 0);
  }

  //----------------------------------------------------
  
  size_t PhotonLibrary::ExtractNOpChannels(TTree* tree) {
    TBranch* channelBranch = tree->GetBranch(OpChannelBranchName.c_str());
    if (!channelBranch) {
      throw art::Exception(art::errors::NotFound)
        << "Tree '" << tree->GetName() << "' has no branch 'OpChannel'";
    }
    
    // fix a new local address for the branch
    char* oldAddress = channelBranch->GetAddress();
    Int_t channel;
    channelBranch->SetAddress(&channel);
    Int_t maxChannel = -1;
    
    // read all the channel values and pick the largest one
    Long64_t iEntry = 0;
    while (channelBranch->GetEntry(iEntry++)) {
      if (channel > maxChannel) maxChannel = channel;
    } // while
    
    LOG_DEBUG("PhotonLibrary")
      << "Detected highest channel to be " << maxChannel << " from " << iEntry
      << " tree entries";
    
    // restore the old branch address
    channelBranch->SetAddress(oldAddress);
    
    return size_t(maxChannel + 1);
    
  } // PhotonLibrary::ExtractNOpChannels()



}
