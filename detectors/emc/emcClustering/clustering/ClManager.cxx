// $Id$
// Author: artur   2016/04/07


//_____________________________________________________________________________
//
// ClManager
//_____________________________________________________________________________

#include "ClManager.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(ClManager)

//_____________________________________________________________________________
ClManager::ClManager():ElementsN_(0),ClustersN_(0),CollectionsN_(0) 
{
  //cout << "<ClManager::ClManager> " << endl;
}

//_____________________________________________________________________________
ClManager::~ClManager() 
{
  DeleteAll();
}

//__________________________________________________________________
void ClManager::DeleteAll() 
{ 
  if (!Collections_.empty()) {
      CL_SSET::iterator it = Collections_.begin();
      for (; it != Collections_.end(); it++) {
           if (it->second) delete it->second;
      }
      Collections_.clear();
  }
  CollectionsN_ = 0; 
  
  if (!Clusters_.empty()) {
      CL_CSET::iterator it = Clusters_.begin();
      for (; it != Clusters_.end(); it++) {
           if (it->second) delete it->second;
      }
      Clusters_.clear();
  }
  ClustersN_ = 0;
  
  if (!Elements_.empty()) {
      CL_ESET::iterator it = Elements_.begin();
      for (; it != Elements_.end(); it++) {
           if (it->second) delete it->second;
      }
      Elements_.clear();
  }
  ElementsN_ = 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
ClusterElement* ClManager::FindElement(int id) const
{
  if (id < 1) return 0;
  CL_ESET::const_iterator it = Elements_.find(id);
  return (it != Elements_.end()) ? it->second : 0;  
}

//_____________________________________________________________________________
Cluster* ClManager::FindCluster(int id) const
{
  if (id < 1) return 0;
  CL_CSET::const_iterator it = Clusters_.find(id);
  return (it != Clusters_.end()) ? it->second : 0;  
}

//_____________________________________________________________________________
ClusterCollection* ClManager::FindCollection(int id) const
{
  if (id < 1) return 0;
  CL_SSET::const_iterator it = Collections_.find(id);
  return (it != Collections_.end()) ? it->second : 0;   
}

//_____________________________________________________________________________
ClusterCollection* ClManager::FindCollection(TString name) const
{
  if (name.IsWhitespace()) return 0;
  CL_SSET::const_iterator it = Collections_.begin();
  ClusterCollection* coll = 0;
  for ( ; it != Collections_.end(); it++) {
       if (name == it->second->GetName()) {   
           if (!coll) coll = it->second;
           else {
              cout << "<-W- <ClManager::FindCollection> " 
                   << "The same name collection found: "
                   << " name = " << name << " id = " << it->first 
                   << endl;               
           }
      }
  }
  return coll;   
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
ClusterElement* ClManager::CreateNewElement() 
{
  CL_ESET::iterator it = Elements_.find(ElementsN_+1); 
  if (it != Elements_.end()) {
      cout << "[E]<ClManager::CreateNewElement> Failed" << endl;
      return 0;
  }
  ElementsN_++;
  ClusterElement* elem = new ClusterElement(ElementsN_);
  Elements_.insert(it,std::pair<Int_t,ClusterElement*>(ElementsN_,elem));
  return elem;
}

//_____________________________________________________________________________    
Cluster* ClManager::CreateNewCluster() 
{
  CL_CSET::iterator it = Clusters_.find(ClustersN_+1); 
  if (it != Clusters_.end()) {
      cout << "[E]<ClManager::CreateNewCluster> Failed " << endl;
      return 0;
  }
  ClustersN_++;
  Cluster* cluster = new Cluster(ClustersN_);
  Clusters_.insert(it,std::pair<Int_t,Cluster*>(ClustersN_,cluster));
  return cluster;
}

//_____________________________________________________________________________
ClusterCollection* ClManager::CreateNewCollection() 
{
  CL_SSET::iterator it = Collections_.find(CollectionsN_+1); 
  if (it != Collections_.end()) {
      cout << "[E]<ClManager::CreateNewCollection> Failed " << endl;
      return 0;
  }
  CollectionsN_++;
  ClusterCollection* collection = new ClusterCollection(CollectionsN_);
  Collections_.insert(it,std::pair<Int_t,ClusterCollection*>(CollectionsN_,collection));
  return collection;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
Bool_t ClManager::InitCLPars() 
{ 
  return (fCLPars.GetSize() > 0); 
} 

//_____________________________________________________________________________
Int_t ClManager::GetCLParsN() 
{ 
  return fCLPars.GetSize(); 
} 

//_____________________________________________________________________________
Double_t ClManager::GetCLPar(Int_t i)
{
  return (i < 0 || i>= fCLPars.GetSize()) ? 0 : fCLPars.At(i);
}

//_____________________________________________________________________________
void ClManager::SetCLPars(unsigned int ncopy, ... ) 
{
  fCLPars.Set(ncopy);
  fCLPars.Reset(0);
  
  va_list vars;
  va_start(vars,ncopy);

  Double_t par = va_arg(vars,Int_t);
  fCLPars.SetAt(par,0);
 
  for (unsigned int i = 1; i<ncopy; i++) {
       par = va_arg(vars,Double_t);
       fCLPars.SetAt(par,i);
       //cout << i << " " << par << " " <<endl;
  }
  va_end(vars);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
Bool_t ClManager::RemoveElement(int id) 
{  
  if (id < 1) return kFALSE;
  if (Elements_.empty()) return kFALSE;
   
  CL_ESET::iterator ite = Elements_.find(id);
  if (ite == Elements_.end()) return kFALSE;
  
  // - remove element from all clusters (without deleting)
  // - remove clusters ids from the element
  
  Int_t* ids = ite->second->GetClustersID();
  bool del = kFALSE;
 
  if (ids) {
      Int_t nc = ite->second->GetNClusters();
      Cluster* cl;
      for (int i(0); i<nc; i++) {
           cl = FindCluster(ids[i]);
           if (cl->RemoveElement(id)) del = kTRUE;
      }
      delete ids;
  }
  
  return del;
}

//_____________________________________________________________________________  
Bool_t ClManager::RemoveCluster(int id) 
{
  if (id < 1) return kFALSE;
  if (Clusters_.empty()) return kFALSE;
   
  CL_CSET::iterator itc = Clusters_.find(id);
  if (itc == Clusters_.end()) return kFALSE;
  
  // - remove cluster from all collections (without deleting) 
  // - remove collections ids from the cluster  
      
  Int_t* ids = itc->second->GetCollectionsID();
  bool del = kFALSE;
  
  if (ids) {
      Int_t nc = itc->second->GetNCollections();
      ClusterCollection* coll;
      for (int i(0); i<nc; i++) {
           coll = FindCollection(ids[i]);
           if (coll && coll->RemoveCluster(id)) del = kTRUE;
      }
      delete ids;
  }
  
  // - remove all elements from the cluster (without deleting)
  // - remove cluster's id from the elements  */
  
  itc->second->RemoveElements();
  
  return del;
}

//_____________________________________________________________________________  
Bool_t ClManager::RemoveCollection(int id) 
{
  if (id < 1) return kFALSE;
  if (Collections_.empty()) return kFALSE;
   
  CL_SSET::iterator its = Collections_.find(id);
  if (its == Collections_.end()) return kFALSE;
  
  // - remove all clusters from the collection (without deleting)
  // - remove collection's id from clusters
   
  its->second->RemoveClusters();
  return kTRUE;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void ClManager::DeleteElement(int id) 
{
  if (id < 1) return;
  if (Elements_.empty()) return;
  
  RemoveElement(id);
   
  CL_ESET::iterator ite = Elements_.find(id);
  if (ite == Elements_.end()) return;
  
  ClusterElement* elem = ite->second;
  Elements_.erase(ite);
  delete elem;
}

//_____________________________________________________________________________  
void ClManager::DeleteCluster(int id) 
{
  if (id < 1) return;
  if (Clusters_.empty()) return;
  
  RemoveCluster(id);
  
  CL_CSET::iterator itc = Clusters_.find(id);
  if (itc == Clusters_.end()) return;
  
  Cluster* cluster = itc->second;
  Clusters_.erase(itc);
  delete cluster;
}

//_____________________________________________________________________________    
void ClManager::DeleteCollection(int id) 
{  
  if (id < 1) return;
  if (Collections_.empty()) return;
  
  RemoveCollection(id);
  
  CL_SSET::iterator it = Collections_.find(id);
  if (it == Collections_.end()) return;
  
  ClusterCollection* coll = it->second;
  Collections_.erase(it);
  delete coll;
}

//__________________________________________________________________
void ClManager::DeleteCollections() 
{ 
  CollectionsN_ = 0; 
  if (Collections_.empty()) return;
  
  CL_SSET::iterator it = Collections_.begin();
  for (; it != Collections_.end(); it++) {     
       it->second->RemoveClusters();
       delete it->second;
  }
  Collections_.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//__________________________________________________________________
void ClManager::Print(Int_t opt) const
{
   cout << "[I]<ClManager::Print> Elements:    " << Elements_.size() << endl;
   cout << "[I]<ClManager::Print> Clusters:    " << Clusters_.size() << endl;
   cout << "[I]<ClManager::Print> Collections: " << Collections_.size() << endl;
}


