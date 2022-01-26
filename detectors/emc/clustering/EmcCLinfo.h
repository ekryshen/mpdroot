// $Id$
// Author: artur   2016/04/25

#ifndef __EMCCLINFO_H__
#define __EMCCLINFO_H__

#include <TObject.h>
#include <TString.h>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcCLinfo                                                                  //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class EmcCLinfo: public TObject {

public:

    EmcCLinfo();
    virtual ~EmcCLinfo();
    
    void Print(const Option_t* opt = 0) const;
    void Save(TString fullname);
    
    void SetNEvents(Int_t ne)  { fNEvents = ne;  }
    void SetCLMethod(Int_t nm) { fCLMethod = nm; }
 
private:
  
    Int_t fNEvents;
    Int_t fCLMethod;
    
    ClassDef(EmcCLinfo,1)
};

#endif  /* __EMCCLINFO_H__ */

