
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "Geometry/Geometry.h"

#include "PhotonPropagation/PhotonLibrary.h"
#include "Simulation/PhotonVoxels.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"


#include <stdio.h>
#include <string.h>

namespace phot{
  
  std::vector<float> PhotonLibrary::EmptyChannelsList; // used for invalid return value
  
  //------------------------------------------------------------

  PhotonLibrary::PhotonLibrary()
  {
    fLookupTable.clear();
  }

  
  //------------------------------------------------------------  

  PhotonLibrary::~PhotonLibrary()
  {
    fLookupTable.clear();
  }
  
  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile(std::string LibraryFile)
  {
    mf::LogInfo("PhotonLibrary") << "Writing photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
   // std::string gdmlFileName;


    tt->Branch("Voxel",      &Voxel,      "Voxel/I");
    tt->Branch("OpChannel",  &OpChannel,  "OpChannel/I");
    tt->Branch("Visibility", &Visibility, "Visibility/F");
 //   tt->Branch("", &gdmlFileName, "gdmlFileName");
 

    for(size_t ivox=0; ivox!=fLookupTable.size(); ++ivox)
      {
	for(size_t ichan=0; ichan!=fLookupTable.at(ivox).size(); ++ichan)
	  {
	    if(fLookupTable[ivox].at(ichan) > 0)
	      {
		Voxel      = ivox;
		OpChannel  = ichan;
		Visibility = fLookupTable[ivox][ichan];
		tt->Fill();
	      }
	  }	
      }
  }

  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile2(std::string LibraryFile, int Nx, int Ny, int Nz, int Nv, std::string gdmlfile)
  {
    mf::LogInfo("PhotonLibrary") << "Writing extended photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;
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
 

    for(size_t ivox=0; ivox!=fLookupTable.size(); ++ivox)
      {
	for(size_t ichan=0; ichan!=fLookupTable.at(ivox).size(); ++ichan)
	  {
	    if(fLookupTable[ivox].at(ichan) > 0)
	      {
		Voxel      = ivox;
		OpChannel  = ichan;
		Visibility = fLookupTable[ivox][ichan];
		nx=Nx;
		ny=Ny;
		nz=Nz;
		nv=Nx*Ny*Nz;
		// gdmlfilepath=const Char_t("file");
		tt->Fill();
	      }
	  }	
      }
  }

  //------------------------------------------------------------

  void PhotonLibrary::CreateEmptyLibrary( size_t NVoxels, size_t NOpChannels)
  {
    fLookupTable.clear();

    fNVoxels     = NVoxels;
    fNOpChannels = NOpChannels;

    fLookupTable.resize(NVoxels);    

    for(size_t ivox=0; ivox!=NVoxels; ivox++)
      {
        fLookupTable[ivox].resize(NOpChannels,0);
      }
  }


  //------------------------------------------------------------

  void PhotonLibrary::LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels, size_t NOpChannels)
  {
    fLookupTable.clear();
    
    mf::LogInfo("PhotonLibrary") << "Reading photon library from input file: " << LibraryFile.c_str()<<std::endl;

    fNVoxels     = NVoxels;
    fNOpChannels = NOpChannels;

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
    
    fLookupTable.resize(NVoxels);    

    mf::LogInfo("PhotonLibrary") <<"Photon lookup table size : "<<  NVoxels << " voxels,  " << NOpChannels<<" channels";


    for(size_t ivox=0; ivox!=NVoxels; ivox++)
      {
	fLookupTable[ivox].resize(NOpChannels,0);
      }
    
    tt->SetBranchAddress("Voxel",      &Voxel);
    tt->SetBranchAddress("OpChannel",  &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);
    
    size_t NEntries = tt->GetEntries();

    for(size_t i=0; i!=NEntries; ++i)
      {
	tt->GetEntry(i);
	if((Voxel<0)||(Voxel>=int(NVoxels))||(OpChannel<0)||(OpChannel>=int(NOpChannels)))
	  {
           //      mf::LogError("PhotonLibrary")<<"Error building photon library, entry out of range: " << Voxel<< " " << OpChannel;
	  }
	else
	  {
	    fLookupTable[Voxel].at(OpChannel) = Visibility;
	  }
      }

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

  float PhotonLibrary::GetCount(size_t Voxel, size_t OpChannel) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels)||/*(OpChannel<0)||*/(OpChannel>=fNOpChannels))
      return 0;   
    else
      return fLookupTable[Voxel].at(OpChannel); 
  }

  //----------------------------------------------------

  void PhotonLibrary::SetCount(size_t Voxel, size_t OpChannel, float Count) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range"; 
    else{
          mf::LogInfo("PhotonLibrary")<<" flookup table size " << fLookupTable.size()<<" "<<OpChannel<<std::endl;
      fLookupTable[Voxel].at(OpChannel) = Count; 
      
      }
  }

  //----------------------------------------------------

  const std::vector<float>* PhotonLibrary::GetCounts(size_t Voxel) const
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      return EmptyList(); // FIXME!!! better to throw an exception!
    else 
      return &(fLookupTable[Voxel]);
  }

  
  
}
