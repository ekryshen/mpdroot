// $Id$
// Author: artur   2016/04/07


//_____________________________________________________________________________
//
// Cluster
//_____________________________________________________________________________

#include "Cluster.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(Cluster)

//_____________________________________________________________________________
Cluster::Cluster():TNamed("cluster","")
{
   Id_.Set(3);
   Value_.Set(3);
}

//_____________________________________________________________________________
Cluster::Cluster(Int_t uid):TNamed("cluster","") 
{  
   Id_.Set(3);
   Value_.Set(3);
  
   Id_[0] = uid;
}

//__________________________________________________________________
void Cluster::clear() 
{ 
   Value_.Set(3);
   Value_.Reset(0);
}

//_____________________________________________________________________________
void Cluster::SetId(Int_t i, Int_t id) 
{ 
   if (i < 1) return; 
   if (i < Id_.fN) { Id_[i] = id; return; }
   Id_.Set(i+1);
   Id_[i] = id;
}

//_____________________________________________________________________________
void Cluster::SetValue(Int_t i, Double_t v) 
{ 
   if (i < 0) return; 
   if (i < Value_.fN) { Value_[i] = v; return; }
   Value_.Set(i+1);
   Value_[i] = v;
}

//_____________________________________________________________________________
void Cluster::AddValue(Int_t i, Double_t v) 
{ 
  if (i < 0) return; 
  if (i < Value_.fN) { Value_[i] += v; return; }
  Value_.Set(i+1);
  Value_[i] = v;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________ 
const ClusterElement* Cluster::GetElement(Int_t id) const
{
   if (GetId() < 1) return 0;
   if (id < 1) return 0;
   
   CL_ESET::const_iterator it = Elements_.find(id);
   return (it != Elements_.end()) ? it->second : 0;
}

//_____________________________________________________________________________
Int_t* Cluster::GetElementsID() const 
{ 
  Int_t n = Elements_.size();
  if (n < 1) return 0;
  Int_t* ids = new Int_t[n];
  CL_ESET::const_iterator it = Elements_.begin();
  for (n = 0; it != Elements_.end(); it++) ids[n++] = it->first;
  return ids;
}

//_____________________________________________________________________________
Int_t* Cluster::GetCollectionsID() const 
{ 
  Int_t n = CollectionsID_.size();
  if (n < 1) return 0;
  Int_t* ids = new Int_t[n];
  CL_IDS::const_iterator it = CollectionsID_.begin();
  for (n = 0; it != CollectionsID_.end(); it++) ids[n++] = *it;
  return ids;
}

//_____________________________________________________________________________
Bool_t Cluster::FindCollection(Int_t id)  const
{
   CL_IDS::const_iterator it = CollectionsID_.find(id); 
   return (it != CollectionsID_.end());
}
   
//_____________________________________________________________________________ 
Bool_t Cluster::FindElement(const ClusterElement* elem) const
{
   if (GetId() < 1) return kFALSE;
   if (!elem) return kFALSE;
   
   Int_t id = elem->GetId();
   CL_ESET::const_iterator it = Elements_.find(id);
   return (it != Elements_.end());  
}

//_____________________________________________________________________________ 
Bool_t Cluster::FindElement(Int_t id) const
{
   if (GetId() < 1) return kFALSE;
   if (id < 1) return kFALSE;
   
   CL_ESET::const_iterator it = Elements_.find(id);
   return (it != Elements_.end());  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________ 
Bool_t Cluster::AddElement(ClusterElement* elem) 
{
   if (GetId() < 1) {
       cout << "[W]<Cluster::AddElement> Cluster has no ID (=0)" << endl;
       return kFALSE;
   }
   if (!elem) {
       cout << "[W]<Cluster::AddElement> No element" << endl;
       return kFALSE;
   }
   
   Int_t id = elem->GetId();
   CL_ESET::const_iterator it = Elements_.find(id);
   if (it == Elements_.end()) Elements_.insert(std::pair<Int_t,ClusterElement*>(id,elem));
   
   elem->AddCluster(GetId());
   return kTRUE;  
}

//_____________________________________________________________________________ 
Bool_t Cluster::RemoveElement(ClusterElement* elem) 
{
   if (GetId() < 1) return kFALSE;
   if (!elem) return kFALSE;
   
   Int_t id = elem->GetId();
   CL_ESET::iterator it = Elements_.find(id);
   if (it != Elements_.end()) Elements_.erase(it);
   
   elem->RemoveCluster(GetId());
   return kTRUE;  
}

//_____________________________________________________________________________ 
Bool_t Cluster::RemoveElement(Int_t id) 
{
   if (GetId() < 1) return kFALSE;
   if (id < 1) return kFALSE;
   
   CL_ESET::iterator it = Elements_.find(id);
   if (it != Elements_.end()) {
       it->second->RemoveCluster(GetId());
       Elements_.erase(it);
   }
   return kTRUE;  
}

//_____________________________________________________________________________
void Cluster::RemoveElements() {
   //
   CL_ESET::iterator it = Elements_.begin();
   Int_t id = GetId();
   for ( ; it != Elements_.end(); it++) it->second->RemoveCluster(id);   
   Elements_.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
Bool_t Cluster::AddCollection(Int_t id) 
{
   if (GetId() < 1) return kFALSE;
   if (id < 1) return kFALSE;
   std::pair<CL_IDS::iterator,Bool_t> it = CollectionsID_.insert(id); 
   return it.second;
}
   
//_____________________________________________________________________________   
Bool_t Cluster::RemoveCollection(Int_t id) 
{
   Int_t n = CollectionsID_.erase(id);
   return (n > 0);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void Cluster::print(Int_t opt) const
{
  
}

