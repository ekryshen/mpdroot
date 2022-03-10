// $Id$
// Author: artur   2016/04/07


//_____________________________________________________________________________
//
// ClusterCollection
//_____________________________________________________________________________

#include "ClusterCollection.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(ClusterCollection)

//_____________________________________________________________________________
ClusterCollection::ClusterCollection():TNamed("collection","")
{
   Id_.Set(3);
}

//_____________________________________________________________________________
ClusterCollection::ClusterCollection(Int_t uid):TNamed("collection","") 
{  
   Id_.Set(3);
   Id_[0] = uid;
}

//_____________________________________________________________________________
void ClusterCollection::SetId(Int_t i, Int_t v) 
{ 
   if (i < 1) return; 
   if (i < Id_.fN) { Id_[i] = v; return; } 
   Id_.Set(i+1);
   Id_[i] = v;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________ 
const Cluster* ClusterCollection::GetCluster(Int_t id) const
{
   if (GetId() < 1) return 0;
   if (id < 1) return 0;
   
   CL_CSET::const_iterator it = Clusters_.find(id);
   return (it != Clusters_.end()) ? it->second : 0;
}

//_____________________________________________________________________________ 
Bool_t ClusterCollection::FindCluster(const Cluster* cl) const
{
   if (GetId() < 1) return kFALSE;
   if (!cl) return kFALSE;
   
   Int_t id = cl->GetId();
   CL_CSET::const_iterator it = Clusters_.find(id);
   return (it != Clusters_.end());  
}

//_____________________________________________________________________________ 
Bool_t ClusterCollection::FindCluster(Int_t id) const
{
   if (GetId() < 1) return kFALSE;
   if (id < 1) return kFALSE;
   
   CL_CSET::const_iterator it = Clusters_.find(id);
   return (it != Clusters_.end());  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
Bool_t ClusterCollection::AddCluster(Cluster* cl) 
{
   if (GetId() < 1) {
       cout << "[W]<ClusterCollection::AddCluster> Collection has no ID (<1)" << endl;
       return kFALSE;
   }
   if (!cl) {
       cout << "[W]<ClusterCollection::AddCluster> No cluster " << endl;
       return kFALSE;
   }
   
   Int_t id = cl->GetId();
   if (id < 1) {
       cout << "[W]<ClusterCollection::AddCluster> Cluster has no ID (<1)" << endl;
       return kFALSE;
   }
   
   CL_CSET::iterator it = Clusters_.find(id);
   if (it == Clusters_.end()) Clusters_.insert(std::pair<Int_t,Cluster*>(id,cl));
   
   cl->AddCollection(GetId());
   return kTRUE;  
}

//_____________________________________________________________________________ 
Bool_t ClusterCollection::RemoveCluster(Cluster* cl) 
{
   if (GetId() < 1) return kFALSE;
   if (!cl) return kFALSE;
   
   Int_t id = cl->GetId();
   CL_CSET::iterator it = Clusters_.find(id);
   if (it != Clusters_.end()) Clusters_.erase(it);
   
   cl->RemoveCollection(GetId());
   return kTRUE;  
}

//_____________________________________________________________________________ 
Bool_t ClusterCollection::RemoveCluster(Int_t id) 
{
   if (GetId() < 1) return kFALSE;
   if (id < 1) return kFALSE;
   
   CL_CSET::iterator it = Clusters_.find(id);
   if (it != Clusters_.end()) {
       it->second->RemoveCollection(GetId());
       Clusters_.erase(it);
   }
   return kTRUE;  
}

//_____________________________________________________________________________
void ClusterCollection::RemoveClusters() {
   //
   CL_CSET::iterator it = Clusters_.begin();
   Int_t id = GetId();
   for ( ; it != Clusters_.end(); it++) it->second->RemoveCollection(id);   
   Clusters_.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void ClusterCollection::print(Int_t opt) const
{
  
}



