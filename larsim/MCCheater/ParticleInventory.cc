////////////////////////////////////////////////////////////////////////
//
// ParticleInventory.cc
// Author: JStock
// EMail:  jason.stock@mines.sdsmt.edu
// 2017-09-12
//
////////////////////////////////////////////////////////////////////////

//STL includes
//ROOT includes
//Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft includes
#include "larsim/MCCheater/ParticleInventory.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "lardata/Utilities/AssociationUtil.h"


namespace cheat{

  //Class members for Inventory

  ParticleInventory::InvParticle::TrackId_t ParticleInventory::InvParticle::Mother() const{
    return part->Mother();
  }

  //----------------------------------------------------------------------
//  const ParticleInventory::InvParticle::TrackId_t ParticleInventory::InvParticle::MotherTrackId(){
//    return this->Mother();
//  }

  //----------------------------------------------------------------------
  //const ParticleInventory::InvParticle::TrackId_t ParticleInventory::InvParticle::Mother(){
  //  return part->Mother();
  //}

  //----------------------------------------------------------------------
   ParticleInventory::InvParticle::TrackId_t ParticleInventory::InvParticle::TrackId() const{
    return tId;
  }

  //----------------------------------------------------------------------
   ParticleInventory::InvParticle::TrackId_t ParticleInventory::InvParticle::EveId() const{
    return EveId(); //This does not return a correct value if the EveId has not been set.
    //a value of std::numeric_limits::min<int>; should be interpreted as a nonid.
  }

  //----------------------------------------------------------------------
  art::Ptr<simb::MCParticle> ParticleInventory::InvParticle::Particle_P() const{
    return part;
  }

  //----------------------------------------------------------------------
  simb::MCParticle ParticleInventory::InvParticle::Particle() const{
    return *part;
  }

  //----------------------------------------------------------------------
  void ParticleInventory::InvParticle::SetEveId( int id) const{
    eveId=id;
  }

  //----------------------------------------------------------------------
  ParticleInventory::InvParticle::InvParticle(int td, art::Ptr<simb::MCParticle> &pt)
    :tId(td), part(pt)
  { }

  //----------------------------------------------------------------------
  ParticleInventory::InvParticle::InvParticle(int td)
    :tId(td)
  { }

  //----------------------------------------------------------------------
  ParticleInventory::InvParticle::InvParticle()
  { }

  //----------------------------------------------------------------------
  /* this is in the wrong place
  void ParticleInventory::Inventory::Add(art::Ptr<simb::MCParticle> pt)
  {
    isSorted=false;
    auto test = this->find(pt->TrackId());
    if ( test != this->end() ){
        throw cet::exception("ParticleInventory::Inventory")
          << "Particle already exists in Inventory.\n";
    }else{
      particleList.emplace_back(pt->TrackId(), pt);
    }
  }

  //----------------------------------------------------------------------
  int ParticleInventory::Inventory::EveId(int tid)
  {//could make this function interchangeable
    auto tmp_part = this->find(tid);
    int tmp = tmp_part->EveId();
    if ( tmp == std::numeric_limits<int>::min()){
      mf::LogDebug("Eve not known. Searching");
      tmp = tmp_part->Mother();
      while(tmp > 0){
        auto old_tmp_part = tmp_part;
        int old_tmp = tmp;
        tmp_part = this->find(tmp);
        tmp=tmp_part->Mother();
        if(tmp==0){
          old_tmp_part->SetEveId(tmp_part->TrackId());
          tmp_part->SetEveId(0);
          break;
        }
      }
    }else{
      return tmp;
    }
    return tmp_part->TrackId();
  }*/

  //----------------------------------------------------------------------
  /*bool ParticleInventory::Inventory::HasParticle(int tid)
  {
    if ( this -> find(tid) != this->end() )
    {
      return true;
    }else{
      return false;
    }
  }

  //----------------------------------------------------------------------
  size_t ParticleInventory::Inventory::NumberOfParticles(){
    return particleList.size();
  }*/

  //----------------------------------------------------------------------
  ParticleInventory::Inventory::Inventory():isSorted(false)
  {}

  //----------------------------------------------------------------------
  void ParticleInventory::Inventory::Add(std::vector<art::Ptr<simb::MCParticle>>::iterator pt)
  {
    isSorted=false;
    auto test = this->find((*pt)->TrackId());
    if ( test != this->end() ){
        throw cet::exception("ParticleInventory::Inventory")
          << "Particle already exists in Inventory.\n";
    }else{
      cheat::ParticleInventory::InvParticle tmp((*pt)->TrackId(), *pt);
      //particleList.emplace_back((*pt)->TrackId(), *pt);
      particleList.push_back(std::move(tmp));
    }
  }

  //----------------------------------------------------------------------
  int ParticleInventory::Inventory::EveId(int tid) const
  {//could make this function interchangeable
    auto tmp_part = this->find(tid);
    int tmp = tmp_part->EveId();
    if ( tmp == std::numeric_limits<int>::min()){
      mf::LogDebug("Eve not known. Searching");
      tmp = tmp_part->Mother();
      auto old_tmp_part = tmp_part;
      //int old_tmp = tmp;
      while(tmp > 0){
        //old_tmp_part=tmp_part;
        //old_tmp=tmp;
        tmp_part = this->find(tmp);
        tmp=tmp_part->Mother();
      }
      if(tmp==0){ //double check this. How can this go wrong. Is 0 guaranteed? That would be better. abs(tmp)>0? That is better. That is the same as ==0 though. Going back to ==0
        old_tmp_part->SetEveId(tmp_part->TrackId());
        tmp_part->SetEveId(0);
      }
    }else{
      return tmp;
    }
    return tmp_part->TrackId();
  }

  //----------------------------------------------------------------------
  bool ParticleInventory::Inventory::IsEve(int tid) const
  {
    int eid = EveId(tid);
    return (eid<=0); //Double check this behavior. Does it always behave correctly if it's the first time this tid has been called?
  }

  //----------------------------------------------------------------------
  bool ParticleInventory::Inventory::HasParticle(int tid)
  {
    if ( this -> find(tid) != this->end() )
    {
      return true;
    }else{
      return false;
    }
  }

  //----------------------------------------------------------------------
  size_t ParticleInventory::Inventory::NumberOfParticles()
  {
    return particleList.size();
  }

  //----------------------------------------------------------------------
  std::vector<art::Ptr<simb::MCParticle>> ParticleInventory::Inventory::GetPrimaries()
  {
    mf::LogWarning("ParticleInventory")<<"There is no guarantee that the GetPrimaries function will return all primaries.\n";
    std::vector<art::Ptr<simb::MCParticle>> tmp;
    //for( auto itr : particleList )
    for( auto itr = particleList.begin(); itr != particleList.end(); ++itr )
    {
      if( itr->EveId() == 0 )
      {
        tmp.push_back(itr->Particle_P());
      }
    }
    return tmp;
  }

  //----------------------------------------------------------------------
  void ParticleInventory::Inventory::clear()
  {
    particleList.clear();
  }

  //----------------------------------------------------------------------
  bool ParticleInventory::Inventory::empty()
  {
    return particleList.empty();
  }

  //----------------------------------------------------------------------
  std::vector<ParticleInventory::InvParticle>::const_iterator ParticleInventory::Inventory::find( int tid) const
  {
    this->sort();
    InvParticle tmp(tid);
    const auto comp = [](const InvParticle &lhs, const InvParticle &rhs) { return lhs.TrackId() < rhs.TrackId();};
    auto rv = std::lower_bound(particleList.cbegin(), particleList.cend(), tmp, comp);
    if ( rv != particleList.cend() && rv->TrackId() == tid )
    {
      return rv;
    }else{
      throw cet::exception("ParticleInventory::Inventory")
        <<"Track "<<tid<<" not found in Inventory. Failed Find.\n";
    }
  }

  //----------------------------------------------------------------------
  //I know this looks weird as const, but it is correct as I DO want find and IsEve, etc to be const
  void ParticleInventory::Inventory::sort() const
  {
    isSorted=true;
    const auto comp = [](const auto &lhs, const auto &rhs){ return lhs.TrackId() < rhs.TrackId();};
    std::sort(particleList.begin(), particleList.end(), comp);
  }


  //----------------------------------------------------------------------
  ParticleInventory::ParticleInventory(const ParticleInventoryConfig& config )
    :fG4ModuleLabel(config.G4ModuleLabel()),
    fEveIdCalculator(config.EveIdCalculator())
  { }

  //----------------------------------------------------------------------
  ParticleInventory::ParticleInventory(const fhicl::ParameterSet& pSet )
    :fG4ModuleLabel(pSet.get<art::InputTag>("G4ModuleLabel", "largeant")),
    fEveIdCalculator(pSet.get<std::string>("EveIdCalculator", "EmEveIdCalculator"))
  { }


  //-----------------------------------------------------------------------
  void ParticleInventory::SetEveIdCalculator(sim::EveIdCalculator *ec)
  {
    throw cet::exception("ParticleInventory::SetEveIdCalculator")
      <<"This is not currently enabled for the ParticleInventory "
      <<"internal art::MCParticle backtracking. The only available "
      <<"EveIdCalculator is BackToEve backtracking.\n";
  }

  //-----------------------------------------------------------------------
  void ParticleInventory::ClearEvent(){
    fInventory.clear();
    fMCTObj.fMCTruthList.clear();
    fMCTObj.fTrackIdToMCTruthIndex.clear();
  }

  //deliverables

  //-----------------------------------------------------------------------
  //TrackIdToParticlePtr
  art::Ptr<simb::MCParticle> ParticleInventory::TrackIdToParticle_P(int const& id) const {
    const auto &part_it = fInventory.find(id);
    if(part_it == fInventory.cend()){
      mf::LogWarning("ParticleInventory") << "Particle with TrackId: "
        << id << " not found in inventory. "
        << "Returning null pointer.";
      return art::Ptr<simb::MCParticle>();
    }
    return part_it->Particle_P();
  }//End TrackIdToParticle

  //----------------------------------------------------------------------
  bool ParticleInventory::TrackIdIsEve(int tid) const
  {
    return fInventory.IsEve(tid);
  }

  //-----------------------------------------------------------------------
  art::Ptr<simb::MCParticle> ParticleInventory::TrackIdToMotherParticle_P(int const& id) const
  {
    return this->TrackIdToParticle_P(fInventory.find(id)->Particle_P()->Mother());
  }

  //-----------------------------------------------------------------------
  art::Ptr<simb::MCParticle> ParticleInventory::TrackIdToEveParticle_P(int const& id) const
  {
    return this->TrackIdToParticle_P(fInventory.EveId(abs(id)));
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::TrackIdToMCTruth_P(int const& id) const
  {
    // find the entry in the MCTruth collection for this track id
    auto mctItr = fMCTObj.fTrackIdToMCTruthIndex.find(abs(id));
    if(mctItr!=fMCTObj.fTrackIdToMCTruthIndex.end()){
      int partIndex = mctItr->second;
      return fMCTObj.fMCTruthList.at(partIndex);
    }else{
      throw cet::exception("ParticleInventory") << "Attempt to find MCTruth for TrackId: "
        << id <<" has failed.";
    }
  }

  //-----------------------------------------------------------------------
  const art::Ptr<simb::MCTruth>& ParticleInventory::ParticleToMCTruth_P(art::Ptr<simb::MCParticle> p) const
  {
    return this->TrackIdToMCTruth_P(p->TrackId());
  }

  //-----------------------------------------------------------------------
  const std::vector< art::Ptr<simb::MCTruth> >& ParticleInventory::MCTruthVector_Ps() const {
    return fMCTObj.fMCTruthList;
  }

  //-----------------------------------------------------------------------
  const std::vector<art::Ptr<simb::MCParticle>> ParticleInventory::MCTruthToParticles_Ps(art::Ptr<simb::MCTruth> const& mct) const
  {
    std::vector<art::Ptr<simb::MCParticle>> ret;
    // sim::ParticleList::value_type is a pair (track Id, particle pointer)
    for (const auto &  iPart: fInventory) {
      if( this->TrackIdToMCTruth_P(iPart.TrackId()) == mct )
        ret.push_back(iPart.Particle_P());
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfTrackIds() const{
    std::set<int> ret;
    for( auto partItr=fInventory.begin(); partItr!=fInventory.end(); ++partItr){
      ret.emplace(partItr->TrackId());
    }
    return ret;
  }

  //-----------------------------------------------------------------------
  std::set<int> ParticleInventory::GetSetOfEveIds() const{
    std::set<int> ret;
    std::set<int> tIds=this->GetSetOfTrackIds();
    for(auto tId : tIds){
      ret.emplace(fInventory.EveId(tId));
    }
    return ret;
  }


} //namespace
