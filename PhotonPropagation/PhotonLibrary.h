////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012

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
    
    const std::vector<float>* GetCounts(size_t Voxel) const;

    
   void StoreLibraryToFile(std::string LibraryFile);
      void StoreLibraryToFile2(std::string LibraryFile, int Nx, int Ny, int Nz, int Nv, std::string
    gdmlfile);
    void LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels, size_t NChannels);
    void CreateEmptyLibrary(size_t NVoxels, size_t NChannels);
    


    
  private:
    // fLookupTable[Voxel]->at(OpChannel) = Count
    std::vector<std::vector<float> > fLookupTable;
    size_t fNOpChannels;
    size_t fNVoxels;
    
    static std::vector<float> EmptyChannelsList;
    
    static std::vector<float>* EmptyList()
      { EmptyChannelsList.clear(); return &EmptyChannelsList; }
  };

}
