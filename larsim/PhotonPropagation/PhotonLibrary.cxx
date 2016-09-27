
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

#include <stdio.h>
#include <string.h>

namespace phot{
  
  std::string const PhotonLibrary::OpChannelBranchName = "OpChannel";
  
  //------------------------------------------------------------

  PhotonLibrary::PhotonLibrary()
  {
    fLookupTable.clear();
    fReflLookupTable.clear();
  }

  
  //------------------------------------------------------------  

  PhotonLibrary::~PhotonLibrary()
  {
    fLookupTable.clear();
    fReflLookupTable.clear();
  }
  
  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile(std::string LibraryFile,bool storeReflected)
  {
    mf::LogInfo("PhotonLibrary") << "Writing photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
    Float_t   ReflVisibility;
   // std::string gdmlFileName;


    tt->Branch("Voxel",      &Voxel,      "Voxel/I");
    tt->Branch(OpChannelBranchName.c_str(),  &OpChannel,  (OpChannelBranchName + "/I").c_str());
    tt->Branch("Visibility", &Visibility, "Visibility/F");
 //   tt->Branch("", &gdmlFileName, "gdmlFileName");
     if(storeReflected)
    { tt->Branch("ReflVisibility", &ReflVisibility, "ReflVisibility/F");
      if (fLookupTable.size() != fReflLookupTable.size())
          throw cet::exception(" Photon Library ") << "Reflected light lookup table is different size than Direct table \n"
				  << "this should not be happening. ";
    }						


    for(size_t ivox=0; ivox!= fNVoxels; ++ivox)
      {
	for(size_t ichan=0; ichan!= fNOpChannels; ++ichan)
	  {
	    Visibility = uncheckedAccess(ivox, ichan);
	    if (Visibility > 0)
	      {
		Voxel      = ivox;
		OpChannel  = ichan;
		if(storeReflected)  ReflVisibility = fReflLookupTable[ivox][ichan];
		
	      

		tt->Fill();
	      }
	  }	
      }
  }

  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile2(std::string LibraryFile, int Nx, int Ny, int Nz, int Nv, std::string gdmlfile,bool storeReflected)
  {
    mf::LogInfo("PhotonLibrary") << "Writing extended photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
    Float_t   ReflVisibility;
    Int_t     nx;
    Int_t     ny;
    Int_t     nz;
    Int_t     nv;
    //const char* gdmlfilepath=gdmlfile.c_str();
   Char_t filepath[100];
    //const char* gdmlfilepath="file";
   //std::string gdmlFileName;
	strcpy(filepath,gdmlfile.c_str());
//filepath="path";
    tt->Branch("Voxel",      &Voxel,      "Voxel/I");
    tt->Branch("OpChannel",  &OpChannel,  "OpChannel/I");
    tt->Branch("Visibility", &Visibility, "Visibility/F");
    tt->Branch("Voxels_x",      &nx,      "nx/I");
    tt->Branch("Voxels_y",  &ny,"ny/I");
    tt->Branch("Voxels_z",      &nz,"nz/I");
    tt->Branch("Voxels_total",  &nv,"nv/I");
    tt->Branch("file_path", filepath, "filepath/C");
     if(storeReflected)
    { tt->Branch("ReflVisibility", &ReflVisibility, "ReflVisibility/F");
      if (fLookupTable.size() != fReflLookupTable.size())
          throw cet::exception(" Photon Library ") << "Reflected light lookup table is different size than Direct table \n"
				  << "this should not be happening. ";
    }						



    for(size_t ivox=0; ivox!= fNVoxels; ++ivox)

      {

        for(size_t ichan=0; ichan!= fNOpChannels; ++ichan)

          {

            Visibility = uncheckedAccess(ivox, ichan);

            if (Visibility > 0)

              {

                Voxel      = ivox;

                OpChannel  = ichan;

                // visibility is already set
		nx=Nx;
		ny=Ny;
		nz=Nz;
		nv=Nx*Ny*Nz;
		// gdmlfilepath=const Char_t("file");
		if(storeReflected)
		ReflVisibility = fReflLookupTable[ivox][ichan];
                tt->Fill();

              }

          }        

      }

  }





  //------------------------------------------------------------

  void PhotonLibrary::CreateEmptyLibrary( size_t NVoxels, size_t NOpChannels)
  {
    fLookupTable.clear();
    fReflLookupTable.clear();

    fNVoxels     = NVoxels;
    fNOpChannels = NOpChannels;


    fLookupTable.resize(LibrarySize(), 0.);

  
    fReflLookupTable.resize(NVoxels);  
    
    for(size_t ivox=0; ivox!=NVoxels; ivox++)
      {
	fReflLookupTable[ivox].resize(NOpChannels,0);
      }

  }


  //------------------------------------------------------------

  void PhotonLibrary::LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels,bool getReflected=false)
  {
    fLookupTable.clear();
    
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

//	if(tt==NULL) tt =  (TTree*)f->Get("pmtresponse/PhotonLibraryData");

      }
    catch(...)
      {
	mf::LogError("PhotonLibrary") << "Error in ttree load, reading photon library: " << LibraryFile.c_str()<<std::endl;
      }
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
    Float_t   ReflVisibility;
    tt->SetBranchAddress("Voxel",      &Voxel);
    tt->SetBranchAddress("OpChannel",  &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);
    if(getReflected)
        tt->SetBranchAddress("ReflVisibility", &ReflVisibility);


    
    fNVoxels     = NVoxels;
    fNOpChannels = PhotonLibrary::ExtractNOpChannels(tt); // EXPENSIVE!!!

    

    fLookupTable.resize(LibrarySize(), 0.);


    size_t NEntries = tt->GetEntries();

    for(size_t i=0; i!=NEntries; ++i) {
      tt->GetEntry(i);

      // Set the visibility at this optical channel
      uncheckedAccess(Voxel, OpChannel) = Visibility;
      // Set # of optical channels to 1 more than largest one seen
      if (OpChannel >= (int)fNOpChannels)
        fNOpChannels = OpChannel+1;




        if(getReflected){
		if (fReflLookupTable[Voxel].size() < fNOpChannels) fReflLookupTable[Voxel].resize(fNOpChannels,0);
	        fReflLookupTable[Voxel].at(OpChannel) = ReflVisibility;

	}
   } // for entries

    // Go through the table and fill in any missing 0's
    for(size_t ivox=0; ivox!=NVoxels; ivox++)
    {
       if(getReflected){
		if (fReflLookupTable[Voxel].size() < fNOpChannels) fReflLookupTable[Voxel].resize(fNOpChannels,0);
	}

    }
    

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

  float PhotonLibrary::GetReflCount(size_t Voxel, size_t OpChannel) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels)||/*(OpChannel<0)||*/(OpChannel>=fNOpChannels))
      return 0;   
    else
      return fReflLookupTable[Voxel].at(OpChannel); 
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

  void PhotonLibrary::SetReflCount(size_t Voxel, size_t OpChannel, float Count) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range"; 
    else
      fReflLookupTable[Voxel].at(OpChannel) = Count; 
  }
  //----------------------------------------------------

  float const* PhotonLibrary::GetCounts(size_t Voxel) const
  { 
    if (Voxel >= fNVoxels) return nullptr;
    else return fLookupTable.data() + uncheckedIndex(Voxel, 0);
  }
  //----------------------------------------------------


  const std::vector<float>* PhotonLibrary::GetReflCounts(size_t Voxel) const
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      return nullptr; 
    else 
      return &(fReflLookupTable[Voxel]);
  }
  

  
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
