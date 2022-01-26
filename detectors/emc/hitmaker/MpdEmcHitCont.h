// $Id$
// Author: artur   2016/04/15

#ifndef __MPDEMCHITCONT_H__
#define __MPDEMCHITCONT_H__

#include <TObject.h>
#include "FairHit.h"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// MpdEmcHitCont                                                              //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class MpdEmcHitCont: public TObject {

public:

    MpdEmcHitCont();
    MpdEmcHitCont(Int_t pdg, Int_t track = -1, Int_t mpdg = 0, Int_t mtrack = -1);
      
    virtual ~MpdEmcHitCont();
    
    void Clear(Option_t* opt = "");
    void Print(const Option_t* opt = 0) const;
    
    inline void Set(Int_t pdg, Int_t track);
    inline void SetMother(Int_t pdg, Int_t track);
    
    inline void SetEnergy(Float_t e) { fEnergy = e; }
    inline void SetMotherEnergy(Float_t e) { fMotherEnergy = e; }
     
    inline void SetAttachedParticle(Int_t pdg, Int_t itr, Float_t en = 0) 
    { fAttPdg = pdg; fAttMCTrackNum = itr; fAttEnergy = en; }
    
    inline void SetGeneration(Int_t gn) { fGeneration = gn; }
    
    inline void ResetDep()          { fDeposit = 0;    } 
    inline void SetDep(Float_t dep) { fDeposit = dep;  }
    inline void AddDep(Float_t dep) { fDeposit += dep; }
     
    inline Int_t GetPdg()              const { return fPdg; }             
    inline Int_t GetMCTrackNum()       const { return fMCTrackNum; } 
    inline Int_t GetEnergy()           const { return fEnergy; }
    inline Int_t GetMotherPdg()        const { return fMotherPdg;  }     
    inline Int_t GetMotherEnergy()     const { return fMotherEnergy; }
    inline Int_t GetMotherMCTrackNum() const { return fMotherMCTrackNum; }
    inline Int_t GetAttPdg()           const { return fAttPdg;  }     
    inline Int_t GetAttEnergy()        const { return fAttEnergy; }
    inline Int_t GetAttMCTrackNum()    const { return fAttMCTrackNum; }
    
    inline Int_t GetGeneration() const { return fGeneration; }  
    
    inline Float_t GetDeposit()  const { return fDeposit; }
    
    inline Bool_t Check(Int_t track);
    inline Bool_t Check(Int_t pdg, Int_t track);
    inline Bool_t CheckMother(Int_t pdg, Int_t track);
    inline Bool_t Check(Int_t pdg, Int_t track, Int_t mpdg, Int_t mtrack);
     
    virtual MpdEmcHitCont* GetCopy() const;
    
protected:

    Int_t   fPdg;               // particle pdg-number
    Int_t   fMCTrackNum;        // particle mc-track number
    Float_t fEnergy;            // particle energy at start vertex
      
    Int_t   fMotherPdg;         // mother particle pdg-number
    Int_t   fMotherMCTrackNum;  // mother particle mc-track number
    Float_t fMotherEnergy;      // particle energy at start vertex
    
    Int_t   fAttPdg;            // attached particle pdg-number
    Int_t   fAttMCTrackNum;     // attached particle mc-track number
    Float_t fAttEnergy;         // attached particle energy at start vertex
   
    Int_t   fGeneration;        // generation number (fGeneration = 1 if fMotherPdg = -1)
    
    Float_t fDeposit;           // energy deposit
    
    ClassDef(MpdEmcHitCont,1)
};

//_____________________________________________________________________________
inline void MpdEmcHitCont::Set(Int_t pdg, Int_t track) 
{ 
  fPdg = pdg; 
  fMCTrackNum = track; 
}

//_____________________________________________________________________________
inline void MpdEmcHitCont::SetMother(Int_t pdg, Int_t track) 
{ 
  fMotherPdg = pdg; 
  fMotherMCTrackNum = track; 
}

//_____________________________________________________________________________
inline Bool_t MpdEmcHitCont::Check(Int_t track) 
{
  return (track != fMCTrackNum) ? kFALSE : kTRUE;
}

//_____________________________________________________________________________
inline Bool_t MpdEmcHitCont::Check(Int_t pdg, Int_t track) 
{
  return (pdg != fPdg || track != fMCTrackNum) ? kFALSE : kTRUE;
}

//_____________________________________________________________________________
inline Bool_t MpdEmcHitCont::CheckMother(Int_t pdg, Int_t track)
{
  return (pdg != fMotherPdg || track != fMotherMCTrackNum) ? kFALSE : kTRUE;
}
    
//_____________________________________________________________________________
inline Bool_t MpdEmcHitCont::Check(Int_t pdg, Int_t track, Int_t mpdg, Int_t mtrack) 
{
  return 
  (       pdg != fPdg 
       || track != fMCTrackNum 
       || mpdg != fMotherPdg 
       || mtrack != fMotherMCTrackNum
  ) ? kFALSE : kTRUE;
}

#endif  /* __MPDEMCHITCONT_H__ */

