// $Id$
// Author: artur   2016/04/15


//_____________________________________________________________________________
//
// EmcClManager
//_____________________________________________________________________________

#include "FairRootManager.h"
#include "EmcClManager.h"
#include "EmcClusterInfo.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TMath.h"
#include <set>

#include <iostream>

using std::cout;
using std::endl;

using namespace TMath;
using namespace std;

ClassImp(EmcClManager)

//_____________________________________________________________________________
EmcClManager::EmcClManager():ClManager(),fCLinfo(new EmcCLinfo),fGeoPar(0),fHitsArray(0),
fMCClustersInfo(0),fRCClustersInfo(0),
fHitsMap(0),fHitsMapSelected(0),fHitsMapClustered(0),
fMCClusters(0),fRCClusters(0)
{  
  fGeoPar = new MpdEmcGeoParWrapper(new MpdEmcGeoPar,kTRUE);
}

//_____________________________________________________________________________
EmcClManager::~EmcClManager() 
{
  if (fCLinfo) delete fCLinfo;
  if (fGeoPar) delete fGeoPar;
  
  //if (fHitsMap) delete fHitsMap;
  
  if (fMCClustersInfo) delete fMCClustersInfo;
  if (fRCClustersInfo) delete fRCClustersInfo;
}

//_____________________________________________________________________________
InitStatus EmcClManager::Init() 
{
  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if (!ioman) {
      cout << "-E- <EmcClManager::Init> "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
  }

  // Get input array
  fHitsArray = (TClonesArray*)ioman->GetObject("MpdEmcHitA");
  if (!fHitsArray) {
      cout << "-W- <EmcClManager::Init> " << "No hits (MpdEmcHitA) array!" << endl;
      return kERROR;
  }
  
  // Create and register output array
  fMCClustersInfo = new TClonesArray("EmcClusterInfo0",100);
  ioman->Register("EmcClusterInfo0", "CL_INFO", fMCClustersInfo, kTRUE);
  
  // Create and register output array
  fRCClustersInfo = new TClonesArray("EmcClusterInfo",100);
  ioman->Register("EmcClusterInfo", "RCL_INFO", fRCClustersInfo, kTRUE);
 
  cout << "-I- <EmcClManager::Init> Intialization successfull" << endl;
  
  fNum = 0;

  return kSUCCESS;
}

//_____________________________________________________________________________    
void EmcClManager::Exec(Option_t* opt)
{
   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   cout << "-I- <EmcClManager::Exec> start Event:  " << fNum+1 << endl;
    
   if (!fHitsArray) Fatal("MpdEmcHitProducerA::Exec", "No array of hits");
   
   cout << "-I- <EmcClManager::Exec> Number of hits in the array: " 
        << fHitsArray->GetEntriesFast() << endl;
   
   fNum++;
   
   SearchHitsForErrors(); 
   
   fMCClustersInfo->Delete();
   fRCClustersInfo->Delete();
   
   Fill();
   
   if (Elements_.size() > 0) {
       if (InitCLPars()) {
           Int_t cl_method = GetCLPar(0); 
           if (cl_method == 1) RecoClusters_1(GetCLPar(1),GetCLPar(2));
           else {
               cout << "-I- <EmcClManager::Exec> Unknown clusterization method: " 
                    << cl_method << endl;     
           }
       }
       else {
           cout << "-I- <EmcClManager::Exec> Set DEFAULT clusterization method: 1; " << endl;   
           RecoClusters_1(0,4);
       }
   }
   
   SaveHistos();
   
   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   
   
}

//_____________________________________________________________________________    
void EmcClManager::SaveHistos() 
{
   if (fHitsMap) fHitsMap->Write();
   if (fHitsMapSelected) fHitsMapSelected->Write();
   if (fHitsMapClustered) fHitsMapClustered->Write();
   if (fMCClusters) fMCClusters->Write(); 
   if (fRCClusters) fRCClusters->Write(); 
}

//_____________________________________________________________________________    
void EmcClManager::Finish()
{
  cout << "<EmcClManager::Finish> ... end " << endl;
  
  if (!fCLinfo) return; 
    
  fCLinfo->SetNEvents(fNum);
  fCLinfo->SetCLMethod(GetCLPar(0));
}

//_____________________________________________________________________________
EmcClusterElement* EmcClManager::CreateNewElement() 
{
  CL_ESET::iterator it = Elements_.find(ElementsN_+1); 
  if (it != Elements_.end()) {
      cout << "[E]<EmcClManager::CreateNewElement> Failed " << endl;
      return 0;
  }
  ElementsN_++;
  EmcClusterElement* elem = new EmcClusterElement(ElementsN_);
  Elements_.insert(it,std::pair<Int_t,ClusterElement*>(ElementsN_,elem));
  return elem;
}

//_____________________________________________________________________________    
EmcCluster* EmcClManager::CreateNewCluster() 
{
  CL_CSET::iterator it = Clusters_.find(ClustersN_+1); 
  if (it != Clusters_.end()) {
      cout << "[E]<EmcClManager::CreateNewCluster> Failed " << endl;
      return 0;
  }
  ClustersN_++;
  EmcCluster* cluster = new EmcCluster(ClustersN_);
  Clusters_.insert(it,std::pair<Int_t,Cluster*>(ClustersN_,cluster));
  return cluster; 
}

//_____________________________________________________________________________
EmcClusterCollection* EmcClManager::CreateNewCollection() 
{
  CL_SSET::iterator it = Collections_.find(CollectionsN_+1); 
  if (it != Collections_.end()) {
      cout << "[E]<EmcClManager::CreateNewCollection> Failed " << endl;
      return 0;
  }
  CollectionsN_++;
  EmcClusterCollection* collection = new EmcClusterCollection(CollectionsN_);
  Collections_.insert(it,std::pair<Int_t,ClusterCollection*>(CollectionsN_,collection));
  return collection;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//_____________________________________________________________________________
Bool_t EmcClManager::SearchHitsForErrors(TClonesArray* digits)
{
   if (!digits) digits = fHitsArray;
   if (!digits) return kTRUE;
    
   MpdEmcHitA* hit; 
   Int_t ierr, nerr(0);
   
   Int_t ne = digits->GetEntriesFast();
   
   for (Int_t i(0); i<ne; i++) {
        hit = (MpdEmcHitA*)digits->At(i);
        ierr = hit->CheckHit(); 
        if (ierr) {
            nerr++;
            //cout << "ERROR NUM " << ierr << endl;
        }
   }
  
   cout << "-I- <EmcClManager::SearchHitsForErrors> Hits: " << ne << " Errors: " << nerr << endl;
   
   return (nerr > 0);
}

//_____________________________________________________________________________
void EmcClManager::Clear() 
{
   DeleteAll(); 
   
   if (fHitsMap) { delete fHitsMap; fHitsMap = 0; }
   if (fHitsMapSelected) { delete fHitsMapSelected; fHitsMapSelected = 0; }
   if (fHitsMapClustered) { delete fHitsMapClustered; fHitsMapClustered = 0; }
   if (fMCClusters) { delete fMCClusters; fMCClusters = 0; }
   if (fRCClusters) { delete fRCClusters; fRCClusters = 0; }
}

//_____________________________________________________________________________
void EmcClManager::SaveCLinfo(TString filename, TString path)
{
   if (!fCLinfo) return;
   if (filename.IsWhitespace()) return;
   if (!path.IsWhitespace()) {
       if (!path.EndsWith("/")) path += "/";
   }
   fCLinfo->Save(path+filename);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void EmcClManager::Fill(TClonesArray* digits) 
{    
   Clear();
   
   // ------------------------------ Create BASE Histos -------------------------------
   fHitsMap = CreateHitMapHisto(Form("Hitmap %d",fNum),Form("Hitmap %d",fNum));
   fHitsMapSelected = CreateHitMapHisto(Form("Hitmap selected %d",fNum),Form("Hitmap selected %d",fNum));
   fHitsMapClustered = CreateHitMapHisto(Form("Hitmap clustered %d",fNum),Form("Hitmap clustered %d",fNum));
   fMCClusters = CreateHitMapHisto(Form("MC-Clusters %d",fNum),Form("MC-Clusters %d",fNum));
   fRCClusters = CreateHitMapHisto(Form("RC-Clusters %d",fNum),Form("RC-Clusters %d",fNum)); 
  
   if (!digits) digits = fHitsArray;
   if (!digits) return;  
   
   MpdEmcHitA* hit;
   TClonesArray* hitcont;
   MpdEmcHitCont* cont;
   Int_t iz, iphi;
   Float_t hdep;
   
   EmcClusterElement* elem;
    
   Int_t nc(0);
   Int_t ng(0);
   Bool_t xfill;
   
   Int_t ne = digits->GetEntriesFast();
    
   for (Int_t i(0); i < ne; i++) {
     
        hit = (MpdEmcHitA*)digits->At(i);
        fGeoPar->zphi(hit->GetUid(),iz,iphi);
        
        hdep = hit->GetDeposit();
        fHitsMap->SetBinContent(iphi+1,iz+1,hdep);
        
        elem = CreateNewElement();
        elem->SetID(iphi,iz);
        elem->SetTotalDep(hdep);
        
        hitcont = hit->GetContent();
        
        if (!hitcont) continue;
        
        nc += hit->GetNCont();
        
        for (Int_t ic(0); ic < hit->GetNCont(); ic++) {
             cont = (MpdEmcHitCont*)hitcont->At(ic);
             xfill = kFALSE;
             if (cont->GetMotherPdg() == 111) 
             {               
                 xfill = kTRUE;
                 
                 if (xfill) {
                     elem->AddCont(cont->GetCopy());
                     //fHitsMapSelected->SetBinContent(iphi+1,iz+1,cont->GetDeposit());
                     ng++;
                 }      
             }
        }
        
        //cout << i << " [" << iz << "," << iphi << "] ";
        //hit->Print();
   }
   
   cout << "-I- <EmcClManager::Fill> Elements: " << ne << ", total content = " << nc 
        << ", selected content= " << ng 
        << endl;
        
   if (nc > 0) MakeMCClusters();    
}

//_____________________________________________________________________________
void EmcClManager::FillClusterInfo(Int_t id, TClonesArray* clinfo, 
                                 Int_t clflag, TH2D* htot, TH2D* hcl) 
{
   EmcClusterCollection* clc = (EmcClusterCollection*)FindCollection(id);
   FillClusterInfo(clc,clinfo,clflag,htot,hcl);
}

//_____________________________________________________________________________
void EmcClManager::FillClusterInfo(TString coll_name, TClonesArray* clinfo, 
                                 Int_t clflag, TH2D* htot, TH2D* hcl) 
{
   EmcClusterCollection* clc = (EmcClusterCollection*)FindCollection(coll_name);
   FillClusterInfo(clc,clinfo,clflag,htot,hcl);
}

//_____________________________________________________________________________
void EmcClManager::FillClusterInfo(EmcClusterCollection* clcoll, TClonesArray* clinfo, 
                                   Int_t clflag, TH2D* hcoll, TH2D* hcl)
{
   if (clcoll && hcoll) clcoll->FillH(hcoll);
   if (!clcoll || !clinfo) return;
   if (!clcoll->FillClusterInfo(clinfo,clflag)) return; 
   FillClusterInfoHisto(clinfo,hcl);
}

//_____________________________________________________________________________
void EmcClManager::FillClusterInfoHisto(TClonesArray* clinfo, TH2D* hcl) 
{
   if (!hcl || !clinfo) return;
   hcl->Reset();
   hcl->Sumw2();
   
   Int_t ne = clinfo->GetEntriesFast();
   EmcClusterInfo* info;
   for (Int_t i(0); i<ne; i++) {
        info = (EmcClusterInfo*)clinfo->At(i);
        if (info) {
            hcl->SetBinContent(info->GetPhi()+1,info->GetZ()+1,info->GetDeposit());
            //cout << i << " phi,z " << info->GetPhi() << "," << info->GetZ() << endl;
        }
   } 
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//_____________________________________________________________________________
void EmcClManager::MakeMCClusters()
{
  cout << "-I- <EmcClManager::MakeMCClusters> Elements: " << Elements_.size() << endl;
  
  EmcClusterCollection* clc = CreateNewCollection();
  clc->SetName("MC-Clusters [pi0]");
  
  std::map <Int_t,EmcCluster*> clusters;
  std::map <Int_t,EmcCluster*>::iterator icl;
  
  EmcCluster* cl;
  EmcClusterElement* elem;
  MpdEmcHitCont* hcont;
  TObjArray* cont;
  
  Float_t thr = 0;
  //thr = 0.005;
  
  CL_ESET::iterator ie = Elements_.begin();
  
  Int_t mid;
  for ( ; ie != Elements_.end(); ie++) {
       elem = (EmcClusterElement*)ie->second;
       if (!elem) continue;
       
       if (elem->GetTotalDep() < thr) continue; //FIXME
          
       cont = elem->GetCont();
       if (!cont) continue;
       for (int i(0); i<cont->GetEntries(); i++) {
            hcont = (MpdEmcHitCont*)cont->At(i);
            if (hcont) {
                mid = hcont->GetMotherMCTrackNum();
                icl = clusters.find(mid);
                if (icl == clusters.end()) {
                    cl = CreateNewCluster();
                    cl->AddElement(elem);
                    cl->SetId(1,hcont->GetMotherPdg());
                    cl->SetId(2,mid);
                    clusters.insert(std::pair<Int_t,EmcCluster*>(mid,cl));
                    clc->AddCluster(cl);                    
                }
                else {             
                    icl->second->AddElement(elem);
                } 
            }
       }
  }
  
  clusters.clear();
  
  EmcClusterCollection* clcx = CreateNewCollection();
  clcx->SetName("MC-Clusters [(gamma,e)<-pi0]");
  
  const CL_CSET& pcls = clc->GetClusters();
  CL_CSET::const_iterator it = pcls.begin();
  
  for (; it != pcls.end(); ++it) {
       
       cl = (EmcCluster*)it->second;
       
       const CL_ESET& elems = cl->GetElements();
       CL_ESET::const_iterator eit = elems.begin();
       
       for (; eit != elems.end(); ++eit) {
            elem = (EmcClusterElement*)eit->second;
            cont = elem->GetCont();
            if (!cont) continue;
            for (int i(0); i<cont->GetEntries(); i++) {
                 hcont = (MpdEmcHitCont*)cont->At(i);
                 if (hcont) {
                     mid = hcont->GetAttMCTrackNum();
                     icl = clusters.find(mid);
                     if (icl == clusters.end()) {
                         cl = CreateNewCluster();
                         cl->AddElement(elem);
                         cl->SetId(1,hcont->GetAttPdg());
                         cl->SetId(2,mid);
                         clusters.insert(std::pair<Int_t,EmcCluster*>(mid,cl));
                         clcx->AddCluster(cl);
                     }                
                     else {
                         icl->second->AddElement(elem);
                     } 
                 }
            }  
       }
  }
 
  //FillClusterInfo("MC-Clusters [pi0]",fMCClustersInfo,1,fHitsMapSelected,fMCClusters); 
  FillClusterInfo("MC-Clusters [(gamma,e)<-pi0]",fMCClustersInfo,1,fHitsMapSelected,fMCClusters);
  
  cout << "-I- <EmcClManager::MakeMCClusters> Collection " << clc->GetName() << ": "   
       << "Clusters = " << clc->GetNClusters() << endl;
  cout << "-I- <EmcClManager::MakeMCClusters> Collection " << clcx->GetName() << ": "   
       << "Clusters = " << clcx->GetNClusters() << endl;
  
  //Print();
}

//_____________________________________________________________________________
void EmcClManager::RecoClusters_1(Float_t thr, Int_t ncut)
{
   cout << "-I- <EmcClManager::RecoClusters_1> Elements: " << Elements_.size() << endl;
   cout << "-I- <EmcClManager::RecoClusters_1> Clusterization method: 1; " 
                << "Parameters: [" << thr << "," << ncut << "] "
                << endl;   
  
   EmcClusterElement* elem;
   Float_t dep;
   
   map<pair<Int_t,Int_t>,EmcClusterElement*> eset;
   pair<Int_t, Int_t> epos;
   
   CL_ESET::iterator ie = Elements_.begin();
   
   for ( ; ie != Elements_.end(); ie++) {
        elem = (EmcClusterElement*)ie->second;
        if (!elem) continue;  
        dep = elem->GetTotalDep(); 
        elem->GetID(epos.first,epos.second);
        if (dep > thr) eset[epos] = elem;
   }
   
   cout << "-I- <EmcClManager::RecoClusters_1> Elements [selected(threshold)/total]: " 
        << eset.size() << "/" << Elements_.size() << endl;
   
   const Int_t NN = 8;
   const Int_t neighb[NN][2] = { {-1,-1}, {-1,0}, {-1, 1}, 
                               { 0,-1}, { 0,1},
                               { 1,-1}, { 1,0}, { 1, 1} };
   
   map<pair<Int_t,Int_t>,EmcClusterElement*>::iterator it;
  
   set<EmcClusterElement*> cl;
   set<EmcClusterElement*>::iterator cit;
   
   //Int_t ncltot = 0, netot = 0;
   Int_t phi, z;
   
   EmcClusterCollection* clc = CreateNewCollection();
   clc->SetName("Reco-Clusters(1) Total");
  
   TString cc_name = Form("Reco-Clusters(1) N(elem) > %d",ncut); 
   EmcClusterCollection* cc = CreateNewCollection();
   cc->SetName(cc_name);
 
   while (!eset.empty()) {
     
     it = eset.begin();
       
     EmcCluster* cluster = CreateNewCluster();  
     clc->AddCluster(cluster);
     
     cl.insert((EmcClusterElement*)it->second);
     eset.erase(it);
     
     while (!cl.empty()) {
       
       cit = cl.begin();
       
       elem = (EmcClusterElement*)(*cit);
       elem->GetID(phi,z);
       
       for (Int_t i(0); i<NN; i++) {
         
            epos.first  = phi + neighb[i][0];
            epos.second = z + neighb[i][1];
            
            if (epos.first > -1 && epos.second > -1) {
                it = eset.find(epos);
                if (it != eset.end()) {
                    cl.insert((EmcClusterElement*)it->second);
                    eset.erase(it);
                }
            }
       }
       
       cluster->AddElement(elem);
       cl.erase(cit);
       
       if (cluster->GetNElements() > ncut) cc->AddCluster(cluster);
     }
     
     //ncltot++;
     //netot += cluster->GetNElements();  
     //cout << "Cluster : " << ncltot << " " << cluster->GetNElements() << endl;
     
   } 
   
   cout << "-I- <EmcClManager::RecoClusters_1> Collection \"" << clc->GetName() << "\" : "
        << " Clusters = " << clc->GetNClusters() << endl;
   cout << "-I- <EmcClManager::RecoClusters_1> Collection \"" << cc->GetName() << "\": " 
        << " Clusters = " << cc->GetNClusters() << endl;
   
   //FillClusterInfo("Reco-Clusters(1) Total",fRCClustersInfo,2,fHitsMapClustered,fRCClusters); 
   
   FillClusterInfo(cc_name,fRCClustersInfo,2,fHitsMapClustered,fRCClusters); 
   
   Print();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//_____________________________________________________________________________
TH2D* EmcClManager::CreateHitMapHisto(TString name, TString title) const
{  
   TH2D* h = new TH2D(name,title,                    
                      fGeoPar->GetNPhi(),0,360,
                      fGeoPar->GetNZ(),fGeoPar->GetZmin(),fGeoPar->GetZmax()); 
   
   h->GetXaxis()->SetTitle("#phi, deg");
   h->GetYaxis()->SetTitle("z, cm");
   
   return h;
}

//_____________________________________________________________________________
void EmcClManager::Print(Int_t opt) const
{
   ClManager::Print();
   cout << "[I]<EmcClManager::Print> Elements:    " << Elements_.size() << endl;
   
   EmcClusterCollection* clc;
   CL_SSET::const_iterator it = Collections_.begin();
   for ( ; it != Collections_.end(); ++it) {
        clc = (EmcClusterCollection*)it->second;
        printf("%-5d  N(clusters) = %-7d  Collection ID: %-5d  NAME: \"%-s\" \n",
               it->first,clc->GetNClusters(),clc->GetId(),clc->GetName());
   }
   
}
