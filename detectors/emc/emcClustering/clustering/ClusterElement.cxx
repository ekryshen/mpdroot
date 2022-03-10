// $Id$
// Author: artur   2016/04/07


//_____________________________________________________________________________
//
// ClusterElement
//_____________________________________________________________________________

#include "ClusterElement.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(ClusterElement)

//_____________________________________________________________________________
ClusterElement::ClusterElement():TNamed("element","") 
{  
  Id_.Set(3);
  Value_.Set(3);
}

//_____________________________________________________________________________
ClusterElement::ClusterElement(Int_t uid):TNamed("element","") 
{  
  Id_.Set(3);
  Value_.Set(3);
  
  Id_[0]= uid;
}

//_____________________________________________________________________________
ClusterElement::ClusterElement(Int_t uid, Int_t id1, Int_t id2):TNamed("element","") 
{  
  Id_.Set(3);
  Value_.Set(3);
  
  Id_[0]= uid;
  Id_[1]= id1;
  Id_[2]= id2; 
}

//_____________________________________________________________________________
void ClusterElement::clear() 
{   
  Value_.Set(3); 
  Value_.Reset(0);
}

//_____________________________________________________________________________
void ClusterElement::SetId(Int_t i, Int_t id) 
{ 
  if (i < 1) return; 
  if (i < Id_.fN) { Id_[i] = id; return; }
  Id_.Set(i+1);
  Id_[i] = id;
}

//_____________________________________________________________________________
void ClusterElement::SetValue(Int_t i, Double_t v) 
{ 
  if (i < 0) return; 
  if (i < Value_.fN) { Value_[i] = v; return; }
  Value_.Set(i+1);
  Value_[i] = v;
}

//_____________________________________________________________________________
void ClusterElement::AddValue(Int_t i, Double_t v) 
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
Int_t ClusterElement::GetNClusters() const 
{ 
  return ClustersID_.size(); 
}

//_____________________________________________________________________________
Int_t* ClusterElement::GetClustersID() const
{
  Int_t n = ClustersID_.size();
  if (n < 1) return 0;
  Int_t* ids = new Int_t[n];
  CL_IDS::const_iterator it = ClustersID_.begin();
  for ( n = 0; it != ClustersID_.end(); it++) ids[n++] = *it;
  return ids;
}

//_____________________________________________________________________________
Bool_t ClusterElement::FindCluster(Int_t id) const
{
  CL_IDS::const_iterator it = ClustersID_.find(id); 
  return (it != ClustersID_.end());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
//_____________________________________________________________________________
Bool_t ClusterElement::AddCluster(Int_t id) 
{
  if (GetId() < 1) return kFALSE;
  if (id < 1) return kFALSE;
  std::pair<CL_IDS::iterator,Bool_t> it = ClustersID_.insert(id); 
  return it.second;
}
 
//_____________________________________________________________________________   
Bool_t ClusterElement::RemoveCluster(Int_t id) 
{
  Int_t n = ClustersID_.erase(id);
  return (n > 0);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void ClusterElement::print(Int_t opt) const
{
  
}   
   
   