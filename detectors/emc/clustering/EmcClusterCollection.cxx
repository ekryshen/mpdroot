// $Id$
// Author: artur   2016/04/15


//_____________________________________________________________________________
//
// EmcClusterCollection
//_____________________________________________________________________________


#include "EmcClusterCollection.h"
#include <iostream>

using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TH2D.h>

ClassImp(EmcClusterCollection)

//_____________________________________________________________________________
EmcClusterCollection::EmcClusterCollection():ClusterCollection() 
{
  
}

//_____________________________________________________________________________
EmcClusterCollection::EmcClusterCollection(Int_t uid):ClusterCollection(uid) 
{  

}

//_____________________________________________________________________________
EmcClusterCollection::~EmcClusterCollection() 
{
 
}

//_____________________________________________________________________________
Bool_t EmcClusterCollection::FillClusterInfo(TClonesArray* info, Int_t clflag)
{
   if (!info) return kFALSE;
   
   EmcCluster* cluster;
   EmcClusterInfo* clinfo;
   CL_CSET::iterator it = Clusters_.begin();
   
   Int_t n(0);
   
   for (; it != Clusters_.end(); it++) {
        cluster = (EmcCluster*)it->second;
        clinfo = new((*info)[n++]) EmcClusterInfo();
        cluster->FillClusterInfo(clinfo);
        clinfo->SetID(cluster->GetId(1));
        clinfo->SetSecondID(cluster->GetId(2));
        clinfo->SetFlag(clflag);
   }

   return kTRUE;
}

//_____________________________________________________________________________
void EmcClusterCollection::FillH(TH2D* h) 
{
   if (!h) return;   
   h->Reset();
   h->Sumw2();
   
   EmcCluster* cluster;
   
   CL_CSET::iterator it = Clusters_.begin();
   for (; it != Clusters_.end(); it++) {
        cluster = (EmcCluster*)it->second;
        cluster->FillH(h);
   }
}
