////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012
#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "TTree.h"
#include "Simulation/PhotonVoxels.h"


namespace phot{
  
  class PhotonLibrary
  {
  public:
    PhotonLibrary();
    ~PhotonLibrary();

    TTree * ProduceTTree();
    

    float GetCount(size_t Voxel, size_t OpChannel);
    void   SetCount(size_t Voxel, size_t OpChannel, float Count);
    
    float GetReflCount(size_t Voxel, size_t OpChannel);
    void   SetReflCount(size_t Voxel, size_t OpChannel, float Count);
    
    const std::vector<float>* GetCounts(size_t Voxel) const;

    float GetReflCount(size_t Voxel, size_t OpChannel);
    void   SetReflCount(size_t Voxel, size_t OpChannel, float Count);
    const std::vector<float>* GetReflCounts(size_t Voxel) const;
    

    void LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels,bool getReflected);

   void StoreLibraryToFile(std::string LibraryFile,bool storeReflected);
      void StoreLibraryToFile2(std::string LibraryFile, int Nx, int Ny, int Nz, int Nv, std::string
    gdmlfile,bool storeReflected=false);


    void CreateEmptyLibrary(size_t NVoxels, size_t NChannels);
    

    int NOpChannels() const { return fNOpChannels; }
    int NVoxels() const { return fNVoxels; }
    
  private:
    // fLookupTable[Voxel]->at(OpChannel) = Count
    std::vector<std::vector<float> > fLookupTable;
    std::vector<std::vector<float> > fReflLookupTable;
    size_t fNOpChannels;
    size_t fNVoxels;
    
    static std::vector<float> EmptyChannelsList;
    
    static std::vector<float>* EmptyList()
      { EmptyChannelsList.clear(); return &EmptyChannelsList; }
  };

}

#endif
