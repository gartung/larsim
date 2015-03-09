////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently mainly used to set optical properties
// for LAr and other optical components
//



#include "LArG4/MaterialPropertyLoader.h"
#include "Utilities/LArProperties.h"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4LogicalSkinSurface.hh"
#include "Geant4/G4OpticalSurface.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {

  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialProperty(std::string Material,
						   std::string Property, 
						   std::map<double, double> PropertyVector,
						   double Unit)
  {
    std::map<double,double> PropVectorWithUnit;
    for(std::map<double,double>::const_iterator it=PropertyVector.begin();
	it!=PropertyVector.end();
	it++)
      {
	PropVectorWithUnit[it->first*eV]=it->second*Unit;
      }
    fPropertyList[Material][Property]=PropVectorWithUnit;
    mf::LogInfo("MaterialPropertyLoader")<<"Added property " 
					 << Material<< "  " 
					 << Property;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialConstProperty(std::string Material, 
							std::string Property, 
							double PropertyValue,
							double Unit)
  {
    fConstPropertyList[Material][Property]=PropertyValue*Unit;
    mf::LogInfo("MaterialPropertyLoader") << "Added const property " 
					  << Material << "  " 
					  << Property << " = " << PropertyValue;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetBirksConstant(std::string Material, 
						double PropertyValue,
						double Unit)
  {
    fBirksConstants[Material]=PropertyValue*Unit;	
    mf::LogInfo("MaterialPropertyLoader") << "Set Birks constant " 
					  << Material;
  }

  //----------------------------------------------  
  void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
  {
    std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
    std::map<std::string,bool> MaterialsSet;
    
    mf::LogInfo("MaterialPropertyLoader") << "UPDATING GEOMETRY";
    std::cout<< "-----------------------------UPDATING GEOMETRY"<<std::endl;
    // Loop over each material with a property vector and create a new material table for it
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
std::cout<<">>>>>>>>>>>>>>material set "<<Material<<std::endl;
      MaterialsSet[Material]=true;
      MaterialTables[Material]=new G4MaterialPropertiesTable;
    }
    
    // Loop over each material with a const property, 
    // if material table does not exist, create one
    for(std::map<std::string,std::map<std::string,double> >::const_iterator i=fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      if(!MaterialsSet[Material]){
	MaterialsSet[Material]=true;
	MaterialTables[Material]=new G4MaterialPropertiesTable;
      }
    }
    
    // For each property vector, convert to an array of g4doubles and 
    // feed to materials table Lots of firsts and seconds!  See annotation 
    // in MaterialPropertyLoader.h to follow what each element is
    
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,std::map<double,double> >::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	std::vector<G4double> g4MomentumVector;
	std::vector<G4double> g4PropertyVector;
	
	for(std::map<double,double>::const_iterator k=j->second.begin(); k!=j->second.end(); k++){
	  g4MomentumVector.push_back(k->first);
	  g4PropertyVector.push_back(k->second);
	}
	int NoOfElements=g4MomentumVector.size();
	MaterialTables[Material]->AddProperty(Property.c_str(),&g4MomentumVector[0], &g4PropertyVector[0],NoOfElements); 
	mf::LogInfo("MaterialPropertyLoader") << "Added property "
					      <<Property
					      <<" to material table " 
					      << Material;

std::cout<< "Added property " <<Property<<" to material table " << Material<<std::endl;
      }
    }
    
    //Add each const property element
    for(std::map<std::string,std::map<std::string,double > >::const_iterator i = fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,double>::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	G4double PropertyValue=j->second;
	MaterialTables[Material]->AddConstProperty(Property.c_str(), PropertyValue); 
	mf::LogInfo("MaterialPropertyLoader") << "Added const property "
					      <<Property
					      <<" to material table " 
					      << Material;
      }
    }
    
    //Loop through geometry elements and apply relevant material table where materials match
    for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){
      G4LogicalVolume* volume = (*i);
      G4Material* TheMaterial = volume->GetMaterial();
      std::string Material = TheMaterial->GetName();
//std::cout<< "Material !!!!!!!!! "<<Material<<" "<<volume->GetName()<<std::endl;


if(Material=="Copper"){
		std::cout<< "copper foil surface set "<<volume->GetName()<<std::endl;
		const G4int num3 = 12;
  		G4double Ephoton3[num3] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV, 7.*eV,15.*eV};
		 G4double Reflectivity_refl2[num3] ={0.902,0.841,0.464,0.379,0.345,0.299,0.287,0.264,0.337,0.3,0.0,0.0};//measurements at Cracow University of Technology- need to be rescaled

  		 G4MaterialPropertiesTable* reflspt2 = new G4MaterialPropertiesTable(); 
  		
  		reflspt2->AddProperty("REFLECTIVITY", Ephoton3, Reflectivity_refl2, num3);

 		 G4OpticalSurface* refl_opsurfc = new G4OpticalSurface("Surface copper",glisur,ground,dielectric_metal);
 		 refl_opsurfc->SetMaterialPropertiesTable(reflspt2);
   		refl_opsurfc -> SetPolish(0.2);
		new G4LogicalSkinSurface("refl_surfacec",volume, refl_opsurfc);
		}

if(Material=="G10"){
		std::cout<< "G10 surface set "<<volume->GetName()<<std::endl;
		const G4int num4 = 12;
  		G4double Ephoton4[num4] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV,6.2*eV,15.0*eV};
		 G4double Reflectivity_refl3[num4] = {0.393,0.405,0.404,0.352,0.323,.243,0.127,0.065,0.068,0.068, 0.0, 0.0};//measurements at Cracow University of Technology- need to be rescaled
  		 G4MaterialPropertiesTable* reflspt3 = new G4MaterialPropertiesTable(); 
  		reflspt3->AddProperty("REFLECTIVITY", Ephoton4, Reflectivity_refl3, num4);

 		 G4OpticalSurface* refl_opsurfg = new G4OpticalSurface("g10 Surface",glisur,ground,dielectric_metal);
 		 refl_opsurfg->SetMaterialPropertiesTable(reflspt3);
   		refl_opsurfg-> SetPolish(0.1);
		new G4LogicalSkinSurface("refl_surfaceg",volume, refl_opsurfg);
		}

	if(Material=="vm2000"){
		std::cout<< "reflector  "<<volume->GetName()<<std::endl;
		const G4int num2 = 12;
  		G4double Ephoton2[num2] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV,6.2*eV,15.0*eV};
		 G4double Reflectivity_refl[num2] = {0.83, 0.83,0.83,0.83,0.83,0.83,0.089,0.089,0.152,0.268,0.0,0.0};
 //{0.33, 0.33,0.33,0.33,0.33,0.33,0.01,0.01,0.07,0.03};//VM2002 data from DOE,  
  		 G4MaterialPropertiesTable* reflspt = new G4MaterialPropertiesTable(); 
  		reflspt->AddProperty("REFLECTIVITY", Ephoton2, Reflectivity_refl, num2);
 		 G4OpticalSurface* refl_opsurf = new G4OpticalSurface("Reflector Surface",unified,groundfrontpainted,dielectric_dielectric);//groundfrontpainted
 		 refl_opsurf->SetMaterialPropertiesTable(reflspt);
		G4double sigma_alpha = 0.8;
refl_opsurf->SetSigmaAlpha(sigma_alpha);
   		//refl_opsurf -> SetPolish(0.65);
		new G4LogicalSkinSurface("refl_surface",volume, refl_opsurf);
		}
      for(std::map<std::string,G4MaterialPropertiesTable*>::const_iterator j=MaterialTables.begin(); j!=MaterialTables.end(); j++){
	if(Material==j->first){
	  TheMaterial->SetMaterialPropertiesTable(j->second);
	  //Birks Constant, for some reason, must be set separately
	  if(fBirksConstants[Material]!=0)
	    TheMaterial->GetIonisation()->SetBirksConstant(fBirksConstants[Material]);
	  volume->SetMaterial(TheMaterial);



		/*
		if(volume->GetName()=="volOpDetSensitive"){
				std::cout<< "Hamamatsu!!!!!!!!! "<<volume->GetName()<<std::endl;
		const G4int num3 = 10;
  		G4double Ephoton2[num3] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV};
		 G4double EFF_pmt[num3] = {0.0, 0.05, 0.24, 0.32, 0.34, 0.36, 0.35, 0.28, 0.20,0.};//VM2002 data from DOE, 
				  G4double Ephoton0[2] = {7.0*eV, 10.14*eV};
  		G4double photocath_ReR2[2]={1.92,1.92};
 		G4double photocath_ImR2[2]={1.69,1.69};
  		 G4MaterialPropertiesTable* reflspt2 = new G4MaterialPropertiesTable(); 
  		reflspt2->AddProperty("EFFICIENCY", Ephoton2, EFF_pmt, num3);
  		reflspt2->AddProperty("REALRINDEX",Ephoton0,photocath_ReR2,2);
  		reflspt2->AddProperty("IMAGINARYRINDEX",Ephoton0,photocath_ImR2,2);
 		 G4OpticalSurface* refl_opsurf2 = new G4OpticalSurface("PMT Surface",glisur,polished,
			 dielectric_metal);
 		 refl_opsurf2->SetMaterialPropertiesTable(reflspt2);
   		refl_opsurf2 -> SetPolish(1.0);
		new G4LogicalSkinSurface("pmt_surface",volume, refl_opsurf2);
		}

		if(volume->GetName()=="volOpDetSensitive_1"){
				std::cout<< "ETL!!!!!!!!! "<<volume->GetName()<<std::endl;
		const G4int num4 = 10;
  		G4double Ephoton3[num4] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV};
		 G4double EFF_pmt2[num4] = {0.0, 0.0, 0.14, 0.18, 0.2, 0.19, 0.18, 0.14, 0.02,0.};//VM2002 data from DOE, 
		  G4double Ephoton[2] = {7.0*eV, 10.14*eV};
  		G4double photocath_ReR3[2]={1.92,1.92};
 		G4double photocath_ImR3[2]={1.69,1.69};
  		 G4MaterialPropertiesTable* reflspt3 = new G4MaterialPropertiesTable(); 
  		reflspt3->AddProperty("EFFICIENCY", Ephoton3, EFF_pmt2, num4);
  		reflspt3->AddProperty("REALRINDEX",Ephoton,photocath_ReR3,2);
  		reflspt3->AddProperty("IMAGINARYRINDEX",Ephoton,photocath_ImR3,2);
 		 G4OpticalSurface* refl_opsurf3 = new G4OpticalSurface("PMT Surface2",glisur,polished,
			 dielectric_metal);
 		 refl_opsurf3->SetMaterialPropertiesTable(reflspt3);
   		refl_opsurf3 -> SetPolish(1.0);
		new G4LogicalSkinSurface("pmt_surface2",volume, refl_opsurf3);
		}
		if(volume->GetName()=="volOpDetSensitive_2"||volume->GetName()=="volOpDetSensitive_3"){
				std::cout<< "sipm!!!!!!!!! "<<volume->GetName()<<std::endl;
		const G4int num5 = 10;
  		G4double Ephoton4[num5] = {1.77*eV, 2.0675*eV, 2.481*eV, 2.819*eV, 2.953*eV, 3.1807*eV, 3.54*eV, 4.135*eV, 4.962*eV, 5.39*eV};
		 G4double EFF_pmt3[num5] = {0.19, 0.32, 0.49, 0.5, 0.48, 0.44, 0.25, 0.15, 0.02,0.};//VM2002 data from DOE, 
  		 G4MaterialPropertiesTable* reflspt4 = new G4MaterialPropertiesTable(); 
  		reflspt4->AddProperty("EFFICIENCY", Ephoton4, EFF_pmt3, num5);

 		 G4OpticalSurface* refl_opsurf4 = new G4OpticalSurface("PMT Surface3",glisur,polished,
			 dielectric_metal);
 		 refl_opsurf4->SetMaterialPropertiesTable(reflspt4);
   		refl_opsurf4 -> SetPolish(1.0);
		new G4LogicalSkinSurface("pmt_surface3",volume, refl_opsurf4);
		}
		*/
	}
      }
    }
  }


  void MaterialPropertyLoader::SetReflectances(std::string Material, std::map<std::string,std::map<double, double> > Reflectances,  std::map<std::string,std::map<double, double> >  DiffuseFractions)
  {
    std::map<double, double> ReflectanceToStore;
    std::map<double, double> DiffuseToStore;
    
    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
	itMat!=Reflectances.end();
	++itMat)
      {
	std::string ReflectancePropName = std::string("REFLECTANCE_") + itMat->first;
	ReflectanceToStore.clear();
//std::cout<<"setting reflectance "<<ReflectancePropName<<std::endl;
	for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
	    itEn!=itMat->second.end();
	    ++itEn)	  
	  {
	    ReflectanceToStore[itEn->first]=itEn->second;
	  }    
	SetMaterialProperty("LAr", ReflectancePropName, ReflectanceToStore,1);
//SetMaterialProperty("TPB", ReflectancePropName, ReflectanceToStore,1);
      }

 for(std::map<std::string,std::map<double,double> >::const_iterator itMat=DiffuseFractions.begin();
	itMat!=DiffuseFractions.end();
	++itMat)
      {
	std::string DiffusePropName = std::string("DIFFUSE_REFLECTANCE_FRACTION_") + itMat->first;
	DiffuseToStore.clear();
	for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
	    itEn!=itMat->second.end();
	    ++itEn)	  
	  {
	    DiffuseToStore[itEn->first]=itEn->second;
	  }    
	SetMaterialProperty("LAr", DiffusePropName, DiffuseToStore,1);
      }
    
  } 
  

  void MaterialPropertyLoader::GetPropertiesFromServices()
  {
    art::ServiceHandle<util::LArProperties>   LarProp;
    
    // wavelength dependent quantities

    SetMaterialProperty( "LAr", "FASTCOMPONENT", LarProp->FastScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "SLOWCOMPONENT", LarProp->SlowScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "LAr", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), cm );
    SetMaterialProperty( "LAr", "RAYLEIGH",      LarProp->RayleighSpectrum(),  cm );
std::cout<<"extra material properties -------------- "<<LarProp->ExtraMatProperties()<<std::endl;
   if(LarProp->ExtraMatProperties()){

std::cout<<"LOADING TPB PROPERTIES !!!!!"<<std::endl;
 SetMaterialProperty( "TPB", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
 SetMaterialProperty( "TPB", "WLSABSLENGTH",        LarProp->TpbAbs(),    m  );
 SetMaterialProperty( "TPB", "WLSCOMPONENT",        LarProp->TpbEm(),    1  );
 SetMaterialConstProperty( "TPB", "WLSTIMECONSTANT",        LarProp->TpbTimeConstant(),    ns  );
 SetMaterialProperty( "vm2000", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
   // SetReflectances("TPB", LarProp->SurfaceTpbReflectances(), LarProp->SurfaceReflectanceTpbDiffuseFractions());
}

    // scalar properties

    SetMaterialConstProperty("LAr", "SCINTILLATIONYIELD",  LarProp->ScintYield(),           1/MeV );
    SetMaterialConstProperty("LAr", "RESOLUTIONSCALE",     LarProp->ScintResolutionScale(), 1);
    SetMaterialConstProperty("LAr", "FASTTIMECONSTANT",    LarProp->ScintFastTimeConst(),   ns);
    SetMaterialConstProperty("LAr", "SLOWTIMECONSTANT",    LarProp->ScintSlowTimeConst(),   ns);
    SetMaterialConstProperty("LAr", "YIELDRATIO",          LarProp->ScintYieldRatio(),      1);
    SetMaterialConstProperty("LAr", "ELECTRICFIELD",       LarProp->Efield(),               kilovolt/cm);

    SetBirksConstant("LAr",LarProp->ScintBirksConstant(), cm/MeV);
    
    SetReflectances("LAr", LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());


    // If we are using scint by particle type, load these

    if(LarProp->ScintByParticleType())
      {
	SetMaterialConstProperty("LAr", "PROTONSCINTILLATIONYIELD",  LarProp->ProtonScintYield(),        1./MeV );
	SetMaterialConstProperty("LAr", "PROTONYIELDRATIO",          LarProp->ProtonScintYieldRatio(),   1.);
	SetMaterialConstProperty("LAr", "MUONSCINTILLATIONYIELD",    LarProp->MuonScintYield(),          1./MeV );
	SetMaterialConstProperty("LAr", "MUONYIELDRATIO",            LarProp->MuonScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "KAONSCINTILLATIONYIELD",    LarProp->KaonScintYield(),          1./MeV );
	SetMaterialConstProperty("LAr", "KAONYIELDRATIO",            LarProp->KaonScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "PIONSCINTILLATIONYIELD",    LarProp->PionScintYield(),          1./MeV );
	SetMaterialConstProperty("LAr", "PIONYIELDRATIO",            LarProp->PionScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "ELECTRONSCINTILLATIONYIELD",LarProp->ElectronScintYield(),      1./MeV );
	SetMaterialConstProperty("LAr", "ELECTRONYIELDRATIO",        LarProp->ElectronScintYieldRatio(), 1.);
      	SetMaterialConstProperty("LAr", "ALPHASCINTILLATIONYIELD",   LarProp->AlphaScintYield(),         1./MeV );
	SetMaterialConstProperty("LAr", "ALPHAYIELDRATIO",           LarProp->AlphaScintYieldRatio(),    1.);
      }
  }

}
