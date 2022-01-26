// $Id$
// Author: artur   2016/04/27

#ifndef __MPDEMCGEOPARWRAPPER_H__
#define __MPDEMCGEOPARWRAPPER_H__

#include <TObject.h>
#include "MpdEmcGeoPar.h"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// MpdEmcGeoParWrapper                                                        //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class MpdEmcGeoParWrapper: public TObject {

public:

    MpdEmcGeoParWrapper(MpdEmcGeoPar* pars = 0, Bool_t is_owner = kFALSE);
    virtual ~MpdEmcGeoParWrapper();
    
    void SetParameters(MpdEmcGeoPar* pars);
    void SetOwner(Bool_t b = kTRUE)  { fIsOwner = kTRUE; }
    
    inline MpdEmcGeoPar* GetParameters() { return fGeoPar; }
    
    inline Int_t   GetNPhi()   const { return fNPhi;   }
    inline Float_t GetPhimin() const { return fPhimin; }
    inline Int_t   GetNZ()     const { return fNZ;     }
    inline Float_t GetZmin()   const { return fZmin;   }
    inline Float_t GetZmax()   const { return fZmax;   }
    inline Float_t GetTmin()   const { return fTmin;   }
    inline Float_t GetTmax()   const { return fTmax;   } 
     
    inline Int_t iz(Int_t uid)  const { return (fIsInit) ? uid/fNPhi : -1; }
    inline Int_t iphi(Int_t uid) const { return (fIsInit) ? uid%fNPhi : -1; }
    
    inline void zphi(Int_t uid, Int_t& i_z, Int_t& i_phi) const;
    
    Float_t GetAngle(Float_t x, Float_t y) const;
    Int_t   GetModUid(Float_t x, Float_t y, Float_t z) const ;
    Bool_t  GetModPos(Int_t uid, Float_t& x, Float_t& y, Float_t& z) const;
    
    void Print(const Option_t* opt = 0) const;
    
private:
  
    MpdEmcGeoPar* fGeoPar;
    Bool_t        fIsOwner;
    Bool_t        fIsInit;
    
    // --- "X" ---
    Int_t   fNPhi;
    Float_t fDPhi;
    Float_t fPhimin; // deg
    
    // --- "Y" ---
    Int_t   fNZ;
    Float_t fDZ;    // cm
    Float_t fZmin;  // cm
    Float_t fZmax;  // cm
    Float_t fTmin;  // deg
    Float_t fTmax;  // deg
    
    Float_t fR; //cm
    
    ClassDef(MpdEmcGeoParWrapper,1)
};

//_____________________________________________________________________________
inline void MpdEmcGeoParWrapper::zphi(Int_t uid, Int_t& i_z, Int_t& i_phi) const
{
  i_z = -1; i_phi = -1;
  if (!fIsInit) return;
  i_z = uid/fNPhi;
  i_phi = uid%fNPhi;
}

#endif  /* __MPDEMCGEOPARWRAPPER_H__ */

