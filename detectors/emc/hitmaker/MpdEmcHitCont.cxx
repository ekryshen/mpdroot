// $Id$
// Author: artur   2016/04/15

//_____________________________________________________________________________
//
// MpdEmcHitCont
//_____________________________________________________________________________

#include "MpdEmcHitCont.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(MpdEmcHitCont)

//_____________________________________________________________________________
MpdEmcHitCont::MpdEmcHitCont():
fPdg(0),fMCTrackNum(-1),fEnergy(0),
fMotherPdg(0),fMotherMCTrackNum(-1),fMotherEnergy(0),
fAttPdg(0),fAttMCTrackNum(-1),fAttEnergy(0),
fGeneration(0),
fDeposit(0)
{
 
}
//_____________________________________________________________________________
MpdEmcHitCont::MpdEmcHitCont(Int_t pdg, Int_t track, Int_t mpdg, Int_t mtrack): 
fPdg(pdg),fMCTrackNum(track),fEnergy(0),
fMotherPdg(mpdg),fMotherMCTrackNum(mtrack),fMotherEnergy(0),
fAttPdg(0),fAttMCTrackNum(-1),fAttEnergy(0),
fGeneration(0),
fDeposit(0)
{
  
}

//_____________________________________________________________________________
MpdEmcHitCont::~MpdEmcHitCont() 
{
 
}

//_____________________________________________________________________________
void MpdEmcHitCont::Clear(Option_t* opt) 
{
   fPdg = 0; 
   fMCTrackNum = -10;
   fEnergy = 0;
   fMotherPdg = 0;
   fMotherMCTrackNum = -10;
   fMotherEnergy = 0;
   fAttPdg = 0;
   fAttMCTrackNum = -10;
   fAttEnergy = 0;
   fGeneration = 0;
   fDeposit = 0;
}

//_____________________________________________________________________________
MpdEmcHitCont* MpdEmcHitCont::GetCopy() const 
{
   MpdEmcHitCont* hcont = new MpdEmcHitCont(fPdg,fMCTrackNum,
                                            fMotherPdg,fMotherMCTrackNum);
   hcont->fGeneration = fGeneration;
   hcont->fEnergy = fEnergy;
   hcont->fMotherEnergy = fMotherEnergy;
   hcont->fDeposit = fDeposit;
   hcont->fAttPdg = fAttPdg;
   hcont->fAttMCTrackNum = fAttMCTrackNum;
   hcont->fAttEnergy = fAttEnergy;
   
   return hcont;
}
 
//_____________________________________________________________________________
void MpdEmcHitCont::Print(const Option_t* opt) const
{
  cout << "-I- <MpdEmcHitCont::Print> " << endl;
  printf("Particle:     %-6d Track: %-6d   Energy: %-12.6f Generation: %-3d Deposit: %-12.6f\n",
         fPdg,fMCTrackNum,fEnergy,fGeneration,fDeposit);
  printf("Prim. Mother: %-6d Track: %-6d   Energy: %-12.6f\n",
         fMotherPdg,fMotherMCTrackNum,fMotherEnergy);
  printf("Attached:     %-6d Track: %-6d   Energy: %-12.6f\n",
         fAttPdg,fAttMCTrackNum,fAttEnergy);
}


