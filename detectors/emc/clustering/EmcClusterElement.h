// $Id$
// Author: artur   2016/04/15

#ifndef __EMCCLUSTERELEMENT_H__
#define __EMCCLUSTERELEMENT_H__

#include "ClusterElement.h"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcClusterElement                                                          //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class TObjArray;
class MpdEmcHitCont;

class EmcClusterElement: public ClusterElement {

public:

    EmcClusterElement();
    EmcClusterElement(Int_t uid);
    EmcClusterElement(Int_t uid, Int_t id1, Int_t id2);
   
    virtual ~EmcClusterElement();
    
    virtual void clear(); 
    virtual void print(Int_t opt = 0) const; 
    
    void  SetID(Int_t id1, Int_t id2);
    void  SetTotalDep(Double_t dep);
    void  AddCont(MpdEmcHitCont* cont);
    
    void  GetID(Int_t& id1, Int_t& id2);

    Int_t GetPhi() { return GetId(1); }
    Int_t GetZ()   { return GetId(2); }
    
    Double_t GetTotalDep();
    Double_t GetContDep();
    Int_t    GetNCont();
    Bool_t   HasCont();
    
    MpdEmcHitCont* GetCont(Int_t i);
    TObjArray*     GetCont() { return fHitCont; }
    
private:

    TObjArray* fHitCont;
    
    ClassDef(EmcClusterElement,1)
};

#endif  /* __EMCCLUSTERELEMENT_H__ */

