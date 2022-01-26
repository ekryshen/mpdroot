// $Id$
// Author: artur   2016/04/15


//_____________________________________________________________________________
//
// EmcCluster
//_____________________________________________________________________________

#include "EmcCluster.h"
#include <iostream>

using std::cout;
using std::endl;

#include <TH2D.h>

ClassImp(EmcCluster)

//_____________________________________________________________________________
EmcCluster::EmcCluster():Cluster() 
{
 
}
//_____________________________________________________________________________
EmcCluster::EmcCluster(Int_t uid):Cluster(uid)
{  
}

//_____________________________________________________________________________
EmcCluster::~EmcCluster() 
{
 
}

//_____________________________________________________________________________
EmcClusterInfo* EmcCluster::CreateClusterInfo()
{
   if (GetNElements() == 0) return 0;
     
   EmcClusterInfo* clinfo = new EmcClusterInfo();
   EmcClusterElement* elem;
   Int_t iphi, iz;
   Float_t dep;
   
   Float_t phi(0), z(0);
   
   CL_ESET::iterator it = Elements_.begin();
   for (; it != Elements_.end(); it++) {
        elem = (EmcClusterElement*)it->second;
        elem->GetID(iphi,iz);
        dep = elem->GetTotalDep();
        phi += iphi*dep;
        z += iz*dep;
        clinfo->AddPoint(iphi,iz,dep);
   }
   
   dep = clinfo->GetDeposit();
   if (dep > 0) clinfo->SetCenter(Int_t(phi/dep),Int_t(z/dep));
   
   return clinfo;
}

//_____________________________________________________________________________
void EmcCluster::FillClusterInfo(EmcClusterInfo* clinfo)
{
   if (!clinfo) return;
   clinfo->Reset();
   
   if (GetNElements() == 0) return;
  
   EmcClusterElement* elem;
   Int_t iphi, iz;
   Float_t dep;
   
   Float_t phi(0), z(0);
   
   CL_ESET::iterator it = Elements_.begin();
   for (; it != Elements_.end(); it++) {
        elem = (EmcClusterElement*)it->second;
        elem->GetID(iphi,iz);
        dep = elem->GetTotalDep();
        phi += iphi*dep;
        z += iz*dep;
        clinfo->AddPoint(iphi,iz,dep);
   }
   
   dep = clinfo->GetDeposit();
   if (dep > 0) clinfo->SetCenter(Int_t(phi/dep),Int_t(z/dep));
}

//_____________________________________________________________________________
void EmcCluster::FillH(TH2D* h) 
{
   if (!h) return;   
   
   EmcClusterElement* elem;
   Int_t iphi, iz;
   Float_t dep;
 
   CL_ESET::iterator it = Elements_.begin();
   for (; it != Elements_.end(); it++) {
        elem = (EmcClusterElement*)it->second;
        elem->GetID(iphi,iz);
        dep = elem->GetTotalDep();
        h->SetBinContent(iphi+1,iz+1,dep);
   }
}
