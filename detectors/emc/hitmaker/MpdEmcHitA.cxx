// $Id$
// Author: artur   2016/04/15


//_____________________________________________________________________________
//
// MpdEmcHitA
//_____________________________________________________________________________

#include "MpdEmcHitA.h"
#include <TMath.h>
#include <TClonesArray.h>


#include <iostream>
using std::cout;
using std::endl;

ClassImp(MpdEmcHitA)

//_____________________________________________________________________________
MpdEmcHitA::MpdEmcHitA(Int_t uid):FairHit(),fHitCont(0),fDeposit(0)
{
   fDetectorID = uid;
}

//_____________________________________________________________________________
MpdEmcHitA::~MpdEmcHitA() 
{
  if (fHitCont) delete fHitCont;
}

//_____________________________________________________________________________
MpdEmcHitCont* MpdEmcHitA::AddContent(Int_t pdg, Int_t track_num, 
                                      Int_t mpdg, Int_t mtrack_num, Bool_t& is_new) 
{
   Int_t ne = 0;
   MpdEmcHitCont* hc;
   
   if (fHitCont) ne = fHitCont->GetEntriesFast();
   else fHitCont = new TClonesArray("MpdEmcHitCont");
   
   if (ne == 0) {
       hc = new((*fHitCont)[ne]) MpdEmcHitCont(pdg,track_num,mpdg,mtrack_num);
       is_new = kTRUE;
       return hc;
   }
   
   for (Int_t i(0); i<ne; i++) {
       hc = (MpdEmcHitCont*)fHitCont->At(i);
       if (hc && hc->Check(pdg,track_num,mpdg,mtrack_num)) {
           is_new = kFALSE;
           return hc;
       }
   }
   
   hc = new((*fHitCont)[ne]) MpdEmcHitCont(pdg,track_num,mpdg,mtrack_num);
   
   is_new = kTRUE;
   return hc; 
}

//_____________________________________________________________________________
MpdEmcHitCont* MpdEmcHitA::FindContent(Int_t pdg, Int_t track_num) 
{
   if (!fHitCont) return 0;
   
   Int_t ne = fHitCont->GetEntriesFast();
   if (ne == 0) return 0;  
  
   MpdEmcHitCont* hc;
   for (Int_t i(0); i<ne; i++) {
       hc = (MpdEmcHitCont*)fHitCont->At(i);
       if (hc && hc->Check(pdg,track_num)) return hc;
   }
 
   return 0;
}

//_____________________________________________________________________________
Int_t MpdEmcHitA::GetNCont() const 
{    
   return fHitCont->GetEntriesFast(); 
}

//_____________________________________________________________________________
MpdEmcHitCont* MpdEmcHitA::GetContent(Int_t i) 
{
   if (!fHitCont || i < 0 || !(i < fHitCont->GetEntriesFast())) return 0;
   return (MpdEmcHitCont*)fHitCont->At(i);
}

//_____________________________________________________________________________
MpdEmcHitCont* MpdEmcHitA::GetContentFast(Int_t i) 
{
   return (MpdEmcHitCont*)fHitCont->At(i);
}
    
//_____________________________________________________________________________
void MpdEmcHitA::AddPoint(Int_t num) 
{
   Int_t size = fPoints.GetSize(); 
   fPoints.Set(size+1);
   fPoints.AddAt(num,size);
}

//_____________________________________________________________________________
void MpdEmcHitA::Print(const Option_t* opt) const
{
   cout << "<MpdEmcHitA::Print> " << "Uid: " << GetUid() 
        << " Points: " << GetNPoints() << " Content: " << GetNCont() 
        << " Deposit: " << GetDeposit() << endl;
}

//_____________________________________________________________________________
Int_t MpdEmcHitA::CheckHit()
{
   if (!fHitCont) return -1;
     
   Int_t nc = GetNCont() ;
   Int_t np = GetNPoints();
   if (nc < 1 || np < 1 || nc > np) return 1;
   
   Float_t dep = 0;
   for (Int_t i(0); i<nc; i++) dep += GetContentFast(i)->GetDeposit();
   
   if (TMath::Abs(fDeposit-dep) < 1.e-6) return 0;
   return 2;
}



