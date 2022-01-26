// $Id$
// Author: artur   2016/04/25

#ifndef __EMCCLUSTERINFO_H__
#define __EMCCLUSTERINFO_H__

#include <TObject.h>
#include <map>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// EmcClusterInfo                                                             //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

typedef std::pair<Int_t,Int_t>    CLPoint;
typedef std::map<CLPoint,Float_t> CLPointsMap;

class EmcClusterInfo: public TObject {

public:

    EmcClusterInfo();
    virtual ~EmcClusterInfo();
    
    void Print(const Option_t* opt = 0) const;
    void Reset();
   
    inline void SetID(Int_t id) { fId = id; }
    inline void SetSecondID(Int_t id) { fId2 = id; }
    inline void SetFlag(Int_t f) { fFlag = f; }
    inline void SetCenter(Int_t iphi, Int_t iz) { fPhi = iphi; fZ = iz; }
    
    void AddPoint(Int_t iz, Int_t iphi, Float_t dep);
    
    inline Int_t   GetID()       const { return fId;  }
    inline Int_t   GetSecondID() const { return fId2;  }
    inline Int_t   GetFlag()     const { return fFlag; }
    inline Int_t   GetNPoints()  const { return fPoints.size(); }
    inline Float_t GetDeposit()  const { return fDeposit; }  
    inline Int_t   GetPhi()      const { return fPhi; }
    inline Int_t   GetZ()        const { return fZ; }
    
    const CLPointsMap& GetPointsMap() const { return fPoints; }
    
protected:
  
    Int_t   fId, fId2;
    Int_t   fFlag;
    Int_t   fZ, fPhi;
    Float_t fDeposit;
   
    CLPointsMap fPoints;

    ClassDef(EmcClusterInfo,1)
};

class EmcClusterInfo0: public EmcClusterInfo {
  
public:
  
    EmcClusterInfo0():EmcClusterInfo() {}
    virtual ~EmcClusterInfo0() {}
    
    ClassDef(EmcClusterInfo0,1)
};

#endif  /* __EMCCLUSTERINFO_H__ */

