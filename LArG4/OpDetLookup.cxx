////////////////////////////////////////////////////////////////////////
/// \file OpDetLookup.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetLookup class.
//
// See comments in the OpDetLookup.h file.
//
// Ben Jones, MIT, 06/04/2010
//


#include "LArG4/OpDetLookup.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"


namespace larg4 {
  OpDetLookup * TheOpDetLookup;
  
  //--------------------------------------------------
  OpDetLookup::OpDetLookup()
  {
    fTheTopChannel=0;
  }

  //--------------------------------------------------
  OpDetLookup * OpDetLookup::Instance()
  {
    if(!TheOpDetLookup){
      TheOpDetLookup = new OpDetLookup;
    }
    return TheOpDetLookup;  
  }

  //--------------------------------------------------
  int OpDetLookup::GetChannel(std::string TheName)
  {
    return fTheChannelMap[TheName];
  }

  //--------------------------------------------------
  int OpDetLookup::GetChannel(G4VPhysicalVolume* TheVolume)
  {
    std::string TheName = TheVolume->GetName();
    return GetChannel(TheName);
  }


  //--------------------------------------------------
  
  int OpDetLookup::FindClosestOpDet(G4VPhysicalVolume* vol, double& distance)
  {
    art::ServiceHandle<geo::Geometry> geom;
    int    ChannelCount = 0;
    
    double MinDistance = UINT_MAX;
    int    ClosestChannel   = -1;
    
    for(size_t c=0; c!=geom->Ncryostats(); c++)
      {
	for(size_t o=0; o!=geom->NOpDet(c); o++)
	  {
	    double xyz[3];
	    geom->Cryostat(c).OpDet(o).GetCenter(xyz);
	    
	    CLHEP::Hep3Vector DetPos(xyz[0],xyz[1],xyz[2]);
	    CLHEP::Hep3Vector ThisVolPos = vol->GetTranslation();
	    
	    ThisVolPos/=cm;
	    
	    //	    std::cout<<"Det: " << xyz[0]<< " " <<xyz[1]<< " " << xyz[2]<<std::endl;
	    //    std::cout<<"Vol: " << ThisVolPos.x()<< " " <<ThisVolPos.y() << " " <<ThisVolPos.z()<<std::endl;
    
	    double Distance = (DetPos-ThisVolPos).mag();
	    if(Distance < MinDistance)
	      {
		MinDistance = Distance;
		ClosestChannel  = geom->OpDetCryoToOpChannel(o, c);
	      }
	    ChannelCount++;
	  }
      }
    if(ClosestChannel<0) 
      {
	throw cet::exception("OpDetLookup Error") << "No nearby OpDet found!\n"; 
      }
    
    distance = MinDistance;
    return ClosestChannel;
  }

  int OpDetLookup::FindClosestOpDet(G4ThreeVector& pos, double& distance, int& ClosestCryo)
  {
    art::ServiceHandle<geo::Geometry> geom;
    int    ChannelCount = 0;
    
    double MinDistance = UINT_MAX;
    int    ClosestChannel   = -1;
    ClosestCryo      = -1;
    
    for(size_t c=0; c!=geom->Ncryostats(); c++)
      {
	for(size_t o=0; o!=geom->NOpDet(c); o++)
	  {
	    double xyz[3];
	    geom->Cryostat(c).OpDet(o).GetCenter(xyz);
	    
	    CLHEP::Hep3Vector DetPos(xyz[0],xyz[1],xyz[2]);
	    CLHEP::Hep3Vector ThisVolPos = pos;
	    
	    ThisVolPos/=cm;
	    
	    //	    std::cout<<"Det: " << xyz[0]<< " " <<xyz[1]<< " " << xyz[2]<<std::endl;
	    //    std::cout<<"Vol: " << ThisVolPos.x()<< " " <<ThisVolPos.y() << " " <<ThisVolPos.z()<<std::endl;
    
	    double Distance = (DetPos-ThisVolPos).mag();
	    if(Distance < MinDistance)
	      {
		MinDistance = Distance;
		ClosestChannel  = geom->OpDetCryoToOpChannel(o, c);
		ClosestCryo     = c;
	      }
	    ChannelCount++;
	  }
      }
    if(ClosestChannel<0) 
      {
	throw cet::exception("OpDetLookup Error") << "No nearby OpDet found!\n"; 
      }
    
    distance = MinDistance;
    return ClosestChannel;
  }


  //--------------------------------------------------
  void OpDetLookup::AddPhysicalVolume(G4VPhysicalVolume * volume)
  {
  
    // mf::LogInfo("Optical") <<"G4 placing sensitive opdet"<<std::endl;
    
    std::stringstream VolName("");
    double Distance     = 0;

    int NearestDetChannel = FindClosestOpDet(volume, Distance);
    
    VolName.flush();
    VolName << volume->GetName() << "_" << NearestDetChannel;
    volume->SetName(VolName.str().c_str());
    
    fTheChannelMap[VolName.str()] = NearestDetChannel;

    // mf::LogInfo("Optical") << "Found closest volume: " << VolName.str().c_str() << " Channel : " << fTheChannelMap[VolName.str()]<<"  distance : " <<Distance<<std::endl; 
     
  }


  //--------------------------------------------------
  int OpDetLookup::GetN()
  {
    return fTheTopChannel;
  }
  
}
