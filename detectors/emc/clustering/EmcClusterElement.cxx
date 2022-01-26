// $Id$
// Author: artur   2016/04/15


//_____________________________________________________________________________
//
// EmcClusterElement
//_____________________________________________________________________________

#include <TObjArray.h>

#include "EmcClusterElement.h"
#include "MpdEmcHitCont.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(EmcClusterElement)

//_____________________________________________________________________________
EmcClusterElement::EmcClusterElement():ClusterElement(),fHitCont(0)
{
  fHitCont = new TObjArray();
}
//_____________________________________________________________________________
EmcClusterElement::EmcClusterElement(Int_t uid):ClusterElement(uid),fHitCont(0)
{
  fHitCont = new TObjArray();  
}

//_____________________________________________________________________________
EmcClusterElement::EmcClusterElement(Int_t uid, Int_t id1, Int_t id2)
: ClusterElement(uid, id1, id2),fHitCont(0)
{  
  fHitCont = new TObjArray();
}

//_____________________________________________________________________________
EmcClusterElement::~EmcClusterElement() 
{
  if (fHitCont) {
      fHitCont->Delete();
      delete fHitCont;
  }
}

//_____________________________________________________________________________
void EmcClusterElement::clear()
{
  ClusterElement::clear();
  
  if (fHitCont) {
      fHitCont->Delete();
      delete fHitCont;
  }
  fHitCont = new TObjArray();
}

//_____________________________________________________________________________ 
void EmcClusterElement::SetID(Int_t id1, Int_t id2)
{
  SetId(1,id1);
  SetId(2,id2);
}

//_____________________________________________________________________________ 
void EmcClusterElement::SetTotalDep(Double_t dep)
{
  SetValue(dep);
}

//_____________________________________________________________________________   
void EmcClusterElement::AddCont(MpdEmcHitCont* cont)
{
  if (fHitCont) fHitCont->Add(cont);
  AddValue(1,cont->GetDeposit());
}

//_____________________________________________________________________________  
Bool_t EmcClusterElement::HasCont()
{
  return (GetNCont() > 0);
}

//_____________________________________________________________________________ 
void EmcClusterElement::GetID(Int_t& id1, Int_t& id2)
{
  id1 = GetId(1);
  id2 = GetId(2);
}

//_____________________________________________________________________________ 
Double_t EmcClusterElement::GetTotalDep() 
{
  return GetValue(); 
}

//_____________________________________________________________________________ 
Double_t EmcClusterElement::GetContDep() 
{
  return GetValue(1); 
}

//_____________________________________________________________________________
Int_t EmcClusterElement::GetNCont() 
{ 
  return (fHitCont) ? fHitCont->GetEntries() : 0; 
}

//_____________________________________________________________________________
MpdEmcHitCont* EmcClusterElement::GetCont(Int_t i) 
{
  if (!fHitCont || i <0 || i >= fHitCont->GetEntries()) return 0;
  return (MpdEmcHitCont*)fHitCont->At(i);
}

//_____________________________________________________________________________   
void EmcClusterElement::print(Int_t opt) const
{ 
}
