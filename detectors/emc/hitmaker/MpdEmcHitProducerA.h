// $Id$
// Author: artur   2016/04/15

#ifndef __MPDEMCHITPRODUCERA_H__
#define __MPDEMCHITPRODUCERA_H__

#include "FairTask.h"
#include "MpdEmcHitA.h"
#include "MpdEmcGeoParWrapper.h"
#include <map>

#include <TH1D.h>
#include <TVector3.h>
#include <set>

class TClonesArray;
class FairMCPoint;
class MpdMCTrack;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// MpdEmcHitProducerA                                                         //
//                                                                            //
// <brief class description>                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class MpdEmcHitProducerA: public FairTask {

public:

    MpdEmcHitProducerA();
    virtual ~MpdEmcHitProducerA();
    
    virtual InitStatus Init();   
    virtual void Exec(Option_t* opt);
    virtual void Finish();
    
public:
  
    Bool_t IsTested() { return fTested; }
    
    TH1D* fHgtheta0;
    TH1D* fHgtheta1;

private:
  
    /* methods */
    MpdEmcHitA*  GetHit(Float_t x, Float_t y, Float_t z);
   
    void GetNGen(MpdMCTrack* track, Int_t& ng);
    void GetVertexMother(MpdMCTrack* trk, Int_t& itrk);
    void AddContent(MpdEmcHitA* hit, FairMCPoint* pnt);
    
    void Test();
  
    /* data memebers */
    MpdEmcGeoParWrapper* fGeoPar;   // Geometry parameters
   
    TClonesArray* fPointArray;   //! Input array of MpdEmcPoints 
    TClonesArray* fMcTrackArray; //! Input array of MCTracks
    TClonesArray* fDigiArray;    //! Output array of MpdEmcHitA
    
    std::map <Int_t,MpdEmcHitA*>  fHitsColl; //! Collection of MpdEmcHitA (pointers)
  
    Bool_t fTested;
    
    ClassDef(MpdEmcHitProducerA,1)
};

#endif  /* __MPDEMCHITPRODUCERA_H__ */

