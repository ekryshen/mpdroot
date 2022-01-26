// $Id$
// Author: artur   2016/04/15

//_____________________________________________________________________________
//
// MpdEmcHitProducerAA
//_____________________________________________________________________________

#include "MpdEmcHitProducerA.h"


#include <TMath.h>
#include <TClonesArray.h>

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairMCPoint.h"

#include "MpdEmcGeoPar.h"
#include "MpdEmcPoint.h"

#include "MpdMCTrack.h"

#include <iostream>
using std::cout;
using std::endl;

#define ADD_CONT_MODE

ClassImp(MpdEmcHitProducerA)

//_____________________________________________________________________________
MpdEmcHitProducerA::MpdEmcHitProducerA():FairTask("Ideal EMC hit Producer"),
fGeoPar(0),fPointArray(0),fMcTrackArray(0),fDigiArray(0),
fTested(kFALSE)
{
   fHgtheta0 = 0;
   fHgtheta1 = 0;
}

//_____________________________________________________________________________
MpdEmcHitProducerA::~MpdEmcHitProducerA() 
{
   if (fHgtheta0) delete fHgtheta0;
   if (fHgtheta1) delete fHgtheta1;
}

//_____________________________________________________________________________
InitStatus MpdEmcHitProducerA::Init() 
{   
    // Get RootManager
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) {
        cout << "-E- MpdEmcHitProducerA::Init: "
             << "RootManager not instantiated!" << endl;
        return kFATAL;
    }

    fGeoPar = new MpdEmcGeoParWrapper(new MpdEmcGeoPar,kTRUE);

    // Get input array
    fPointArray = (TClonesArray*)ioman->GetObject("MpdEmcPoint");
    if (!fPointArray) {
        cout << "-W- MpdEmcHitProducerA::Init: " << "No EmcPoint array!" << endl;
        return kERROR;
    }
    
    fMcTrackArray = (TClonesArray*)ioman->GetObject("MCTrack");
    
    if (!fMcTrackArray) {
        cout << "-W- MpdEmcHitProducerA::Init: " << "No MCTrack array!" << endl;
        return kERROR;
    }

    // Create and register output array
    fDigiArray = new TClonesArray("MpdEmcHitA",100);
    ioman->Register("MpdEmcHitA", "EMC", fDigiArray, kTRUE);

    cout << "-I- MpdEmcHitProducerA: Intialization successfull" << endl;

    return kSUCCESS;

}

//________________________________________________________________________________________
void MpdEmcHitProducerA::Finish() 
{
    cout << "-I- MpdEmcHitProducerA: Finish" << endl;
}

//________________________________________________________________________________________
MpdEmcHitA* MpdEmcHitProducerA::GetHit(Float_t x, Float_t y, Float_t z)
{
    Int_t uid = fGeoPar->GetModUid(x,y,z);
    
    if (uid < 0) return 0;
      
    std::map <Int_t,MpdEmcHitA*>::iterator it = fHitsColl.find(uid);
    if (it != fHitsColl.end()) return it->second;
    
    Int_t ne = fDigiArray->GetEntriesFast();
    MpdEmcHitA* hit = new((*fDigiArray)[ne]) MpdEmcHitA(uid);
    fHitsColl.insert(std::pair<Int_t,MpdEmcHitA*>(uid,hit));
    
    fGeoPar->GetModPos(uid,x,y,z);
    hit->SetXYZ(x,y,z);
  
    return hit;
}

//________________________________________________________________________________________
void MpdEmcHitProducerA::GetNGen(MpdMCTrack* trk, Int_t& ng) 
{
   ng++;
   Int_t mId = trk->GetMotherId();
   if (mId != -1) GetNGen((MpdMCTrack*)fMcTrackArray->At(mId),ng);
}
 
//________________________________________________________________________________________
void MpdEmcHitProducerA::GetVertexMother(MpdMCTrack* trk, Int_t& itrk) 
{
   Int_t mtrk = trk->GetMotherId(); 
   if (mtrk != -1) {
       itrk = mtrk;
       GetVertexMother((MpdMCTrack*)fMcTrackArray->At(itrk),itrk);
   }
}

//________________________________________________________________________________________
void MpdEmcHitProducerA::AddContent(MpdEmcHitA* hit, FairMCPoint* pnt) 
{   
    Int_t trkId = pnt->GetTrackID();
    MpdMCTrack* trk = (MpdMCTrack*)fMcTrackArray->At(trkId); 
    Int_t pdg = trk->GetPdgCode();  
    MpdMCTrack* mtrk; 
    
    Int_t mtrkId;
    Int_t mpdg = 0;
    Float_t men = 0;
    
//     mtrkId = trk->GetMotherId();     //-> get real mother
//     if (mtrkId != -1) {   
//         mtrk = (MpdMCTrack*)fMcTrackArray->At(mtrkId);
//         //mtrk = (MpdMCTrack*)fMcTrackArray->At(trkId); 
//         mpdg = mtrk->GetPdgCode();
//         men = mtrk->GetEnergy();
//     }
      
    mtrkId = trkId;                //-> get vertex mother
    GetVertexMother(trk,mtrkId);   //-> get vertex mother
   
    mtrk = (MpdMCTrack*)fMcTrackArray->At(mtrkId);
    
    if (mtrk->GetPdgCode() == 111 && mtrk->GetMotherId() == -1) {
        mpdg = mtrk->GetPdgCode();
        men = mtrk->GetEnergy();
    }
    
    Bool_t is_new;
    
    MpdEmcHitCont* hitc = hit->AddContent(pdg,trkId,mpdg,mtrkId,is_new);
    
    if (is_new) {
        hitc->SetEnergy(trk->GetEnergy());
        hitc->SetMotherEnergy(men);   
    }
      
    hitc->AddDep(pnt->GetEnergyLoss()); 
    
//     static int nn(0), nn1(0);
//     if (mtrk->GetPdgCode() == 111) {
//         if (is_new) nn1++;
//         cout << "ADD content " << ++nn << " " << nn1 << endl;
//     }
        
    Int_t ngen = 0; 
    if (is_new) {    
        GetNGen(trk,ngen);
        hitc->SetGeneration(ngen);   
    }
    
    ngen = hitc->GetGeneration();  
    
    if (ngen < 2) return;
    
   
    if (ngen == 2) {
        mtrkId = trk->GetMotherId();
        mtrk = (MpdMCTrack*)fMcTrackArray->At(mtrkId);    
        if (mtrk->GetPdgCode() == 111) {
            hitc->SetAttachedParticle(pdg,trkId,trk->GetEnergy()); 
        }
        return;
    }
    
    trkId = trk->GetMotherId();
    trk = (MpdMCTrack*)fMcTrackArray->At(trkId);
           
    while (kTRUE) {        
           mtrkId = trk->GetMotherId();
           if (mtrkId == -1) return;
           mtrk = (MpdMCTrack*)fMcTrackArray->At(mtrkId);         
           if (mtrk->GetPdgCode() == 111) {
               hitc->SetAttachedParticle(trk->GetPdgCode(),trkId,trk->GetEnergy()); 
               return;
           }
           else {
              trkId = mtrkId;
              trk = mtrk;
           }
    }   
}

//________________________________________________________________________________________
void MpdEmcHitProducerA::Exec(Option_t* opt) 
{
    cout << "-I- MpdEmcHitProducerA::Exec: started with " 
         << " tracks/points = " 
         << fMcTrackArray->GetEntriesFast() << "/" << fPointArray->GetEntriesFast() << endl;
    
    if (!fDigiArray) Fatal("MpdEmcHitProducerA::Exec", "No array of digits");

#ifdef ADD_CONT_MODE    
    //Test();
#endif
    
    fHitsColl.clear();
    fDigiArray->Delete();
    
    FairMCPoint* pnt;
    MpdEmcHitA* hit;
    
    for (Int_t i(0); i<fPointArray->GetEntriesFast(); i++) {
      
         pnt = (FairMCPoint*)fPointArray->At(i); 
         hit = GetHit(pnt->GetX(),pnt->GetY(),pnt->GetZ());
  
         if (!hit) continue;
             
         hit->AddPoint(i);
         hit->AddDeposit(pnt->GetEnergyLoss());
         
#ifdef ADD_CONT_MODE           
         AddContent(hit,pnt);      
#endif
         
    } 
       
//     Int_t np = 0;
//     for (Int_t i(0); i<fDigiArray->GetEntriesFast(); i++) {
//          hit = (MpdEmcHitA*)fDigiArray->At(i);
//          np += hit->GetNPoints();
//          //cout << i << " " << hit->GetNPoints() << " " << hit->GetNCont() << " " << endl;
//     }
    
     cout << "-I- MpdEmcHitProducerA::Exec: Hits/Points = " << fDigiArray->GetEntriesFast() 
                      << "/" << fPointArray->GetEntriesFast() 
//                      << "/" << np 
                      << endl;
  
//    fGeoPar->Print();
    
}

//________________________________________________________________________________________
void MpdEmcHitProducerA::Test()
{
    fTested = kTRUE;
    
    if (!fHgtheta0) fHgtheta0 = new TH1D("htheta 0","htheta 0",180/5,0,180);
    if (!fHgtheta1) fHgtheta1 = new TH1D("htheta 1","htheta 1",180/5,0,180);
    
    Int_t npi0 = 0;
    Int_t nempty = 0;
    Int_t ndaughts = 0, ndaughts_in = 0;
    Int_t ngammas = 0;
   
    Int_t tt1(0), tt2(0);
    Int_t np1 = 0, np2 = 0;
    Int_t ngp = 0;
  
    const Float_t tt_min = fGeoPar->GetTmin();
    const Float_t tt_max = fGeoPar->GetTmax();
    
    for (UInt_t i = 0; i < (UInt_t) fMcTrackArray->GetEntriesFast(); ++i) {
         
         MpdMCTrack* tr = (MpdMCTrack*)fMcTrackArray->At(i);
      
         np1 += tr->GetNPoints(kECAL); //tr->GetTotalPoints();
         np2 += tr->GetNPoints(kECAL);
         
         if (tr->GetMotherId() == -1 && tr->GetPdgCode() == 111) npi0++;
         
         if (tr->GetNPoints(kECAL) > 0) nempty++;
         //if (tr->GetTotalPoints() > 0) nempty++;
         
         Int_t mid = tr->GetMotherId();
         MpdMCTrack* mtr = 0;
         if (mid != -1) mtr = (MpdMCTrack*)fMcTrackArray->At(mid);
                
         if (mtr && mtr->GetMotherId() == -1 && mtr->GetPdgCode() == 111) {
           
             ndaughts++;  
             
             //cout << tr->GetPdgCode() << endl;
             
             TVector3 mom;
             mtr->GetMomentum(mom);                 
             Float_t tt = mom.Theta()*TMath::RadToDeg();
             
             if (tt > tt_min && tt < tt_max) ndaughts_in++;
             
             if (tr->GetPdgCode() == 22) {    
               
                 ngammas++;
                 //ngp += tr->GetTotalPoints();
                 ngp += tr->GetNPoints(kECAL);
             }   
             
             fHgtheta0->Fill(tt);
             tt1++;                
             if (tt > tt_min && tt < tt_max) {
                 fHgtheta1->Fill(tt);
                 tt2++;
             }                            
      
         }
    }
    
    Int_t npi0part = 0;
    Int_t npi0points = 0;
    
    for (UInt_t i = 0; i < (UInt_t) fMcTrackArray->GetEntriesFast(); ++i) {
         
         MpdMCTrack* tr = (MpdMCTrack*)fMcTrackArray->At(i);
         if (!tr) continue;

         if (tr->GetPdgCode() == 111) continue;
           
         Int_t mi = i;
         GetVertexMother(tr,mi);
         MpdMCTrack* mtr = (MpdMCTrack*)fMcTrackArray->At(mi);
         
         if (mtr && mtr->GetPdgCode() == 111 && mtr->GetMotherId() == -1) {
             npi0part++;
             npi0points += tr->GetNPoints(kECAL);
             //npi0points += tr->GetTotalPoints();
         }
    }
    
    FairMCPoint* pp;
    
    std::set<Int_t> gammas;
    std::set<Int_t> electrons;
    std::set<Int_t> allpart;
     
    for (Int_t i(0); i<fPointArray->GetEntriesFast(); i++) {   
         pp = (FairMCPoint*)fPointArray->At(i);       
         MpdMCTrack* tr = (MpdMCTrack*)fMcTrackArray->At(pp->GetTrackID());
         
         allpart.insert(pp->GetTrackID());
         
         if (TMath::Abs(tr->GetPdgCode()) == 11) { 
             electrons.insert(pp->GetTrackID());
         }
         else if (tr->GetPdgCode() == 22) {
             Int_t mid = tr->GetMotherId();
             if (mid != -1) {
                 MpdMCTrack* mtr = (MpdMCTrack*)fMcTrackArray->At(mid);
                 if (mtr && mtr->GetPdgCode() == 111 && mtr->GetMotherId() == -1) {
                     gammas.insert(pp->GetTrackID());
                 }
             }
         }       
    }
    
    std::set<Int_t>::iterator it = gammas.begin();
    Int_t gg, ngg = 0;
    
    Int_t np = 0;
    
    for ( ; it != gammas.end(); it++) {
         MpdMCTrack* tr = (MpdMCTrack*)fMcTrackArray->At(*it);
         gg = 0;
         GetNGen(tr,gg);
         if (gg == 2) ngg++;
         np += tr->GetNPoints(kECAL);
         //np += tr->GetTotalPoints();
    }
    
    cout << "--------------------------------------------------------" << endl;
     
    cout << "-I- MpdEmcHitProducer::Test: Ntr [empty / total] = "
         << nempty << "/" << fMcTrackArray->GetEntriesFast() << endl;
         
    cout << "-I- MpdEmcHitProducer::Test: Ntr [in(daughts) / total(daughts)] = "
         << tt2 << "/" << tt1 << endl; 
         
    cout << "-I- MpdEmcHitProducer::Test: Ntr [pi0 / total] = " << npi0 
         << "/" << fMcTrackArray->GetEntriesFast() << endl; 
         
    cout << "-I- MpdEmcHitProducer::Test: Ntr [ndaughters / 2*pi0 / ngamma / ngamma_points] = " 
         << ndaughts << "/" << npi0*2 << "/" <<  ngammas << "/" << ngp << endl; 
   
    cout << "-I- MpdEmcHitProducer::Test: Ntr [ndaughters (in) / ndaughters (total)] "
         << ndaughts_in << "/" << ndaughts << endl;
         
    //cout << endl;
    
    cout << "-I- MpdEmcHitProducer::Test: Ntr [mc-tracks from pi0] = " << npi0part << endl;
    cout << "-I- MpdEmcHitProducer::Test: Ntr [mc-points from pi0] = " << npi0points << endl;
    
    //cout << endl;
    
    cout << "-I- MpdEmcHitProducer::Test: Points [total / mc-total / mc-ECAL] = " 
         << fPointArray->GetEntriesFast() << "/" << np1 << "/" << np2 << endl;
    
    //cout << endl;
    
    cout << "-I- MpdEmcHitProducer::Test: Points:NGAMMA [tracks / tracks,gen = 2 / points,tracks] = " 
         << gammas.size() << "/" << ngg << "/" << np << endl;
         
    cout << "-I- MpdEmcHitProducer::Test: Points:NELECT [tracks(electrons) / total] = " 
         << electrons.size() << "/" << allpart.size() << endl;     
         
    cout << "--------------------------------------------------------" << endl;      
          
}