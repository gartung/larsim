////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012
#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "TTree.h"
#include "larsim/Simulation/PhotonVoxels.h"


namespace phot{
  
  class PhotonLibrary
  {
  public:
    PhotonLibrary();
    ~PhotonLibrary();

    TTree * ProduceTTree();
    

    float GetCount(size_t Voxel, size_t OpChannel) const;
    void SetCount(size_t Voxel, size_t OpChannel, float Count);
    
    /// Returns a pointer to NOpChannels() visibility values, one per channel
    float const* GetCounts(size_t Voxel) const;

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
    
    /// Returns the number of elements in the library
    size_t LibrarySize() const { return fNVoxels * fNOpChannels; }
    
  private:
    // fLookupTable[unchecked_index(Voxel, OpChannel)] = Count
    // for each voxel, all NChannels() channels are stored in sequence
    std::vector<float> fLookupTable;
    std::vector<std::vector<float> > fReflLookupTable;

    size_t fNOpChannels;
    size_t fNVoxels;
    
    /// Returns the index of visibility of specified voxel and cell
    size_t uncheckedIndex(size_t Voxel, size_t OpChannel) const
      { return Voxel * fNOpChannels + OpChannel; }
    
    /// Unchecked access to a visibility datum
    float const& uncheckedAccess (size_t Voxel, size_t OpChannel) const
      { return fLookupTable[uncheckedIndex(Voxel, OpChannel)]; }
    
    /// Unchecked access to a visibility datum
    float& uncheckedAccess(size_t Voxel, size_t OpChannel)
      { return fLookupTable[uncheckedIndex(Voxel, OpChannel)]; }
    
    /// Name of the optical channel number in the input tree
    static std::string const OpChannelBranchName;
    
    /// Returns the number of optical channels in the specified tree
    static size_t ExtractNOpChannels(TTree* tree);
    
  };

}

#endif
