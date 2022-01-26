// $Id$
// Author: artur   2016/04/15

#ifndef __MPDEMCHITA_H__
#define __MPDEMCHITA_H__

#include <TObject.h>
#include <TArrayI.h>
#include "MpdEmcHitCont.h"

class TClonesArray;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// MpdEmcHitA                                                                 //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class MpdEmcHitA : public FairHit {

public:

    MpdEmcHitA(Int_t uid = -1);
    virtual ~MpdEmcHitA();
    
    void Print(const Option_t* opt = 0) const;
    
    MpdEmcHitCont* AddContent(Int_t pdg, Int_t track_num, Int_t mpdg, Int_t mtrack_num, Bool_t& is_new);
    MpdEmcHitCont* FindContent(Int_t pdg, Int_t track_num);
   
    void SetUid(Int_t id) { fDetectorID = id; }
    void AddDeposit(Float_t dep) { fDeposit += dep; } 
    void AddPoint(Int_t num);
   
    Int_t   GetUid() const     { return fDetectorID; }
    Float_t GetDeposit() const { return fDeposit; }
    Int_t   GetNPoints() const { return fPoints.GetSize(); }
    Int_t   GetNCont() const;
    
    const Int_t*  GetPoints()  { return fPoints.GetArray(); }
    TClonesArray* GetContent() { return fHitCont; }
    
    MpdEmcHitCont* GetContent(Int_t i);
    MpdEmcHitCont* GetContentFast(Int_t i);
     
    Int_t CheckHit();
    
private:
   
    TClonesArray* fHitCont; // Output array of MpdEmcHitCont
    TArrayI       fPoints;  // MC-points (numbers), which have contribution in the hit
    Float_t       fDeposit; // total deposit    

    ClassDef(MpdEmcHitA,1)
};

#endif  /* __MPDEMCHITA_H__ */

