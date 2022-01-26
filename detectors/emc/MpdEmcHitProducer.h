// -------------------------------------------------------------------------
// -----                 MpdEmcHitproducer header file                 -----
// -------------------------------------------------------------------------

#ifndef MPDEMCHITPRODUCER_H
#define MPDEMCHITPRODUCER_H 1

#include <iostream>
#include "FairTask.h"
#include "MpdEmcHit.h"
#include "MpdEmcGeoPar.h"


class TClonesArray;

class MpdEmcHitProducer : public FairTask {
public:

    /** Default constructor **/
    MpdEmcHitProducer();


    /** Destructor **/
    ~MpdEmcHitProducer();


    /** Virtual method Init **/
    virtual InitStatus Init();


    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);
    void virtual Finish();

private:

    /** Input array of MpdEmcPoints **/
    TClonesArray* fPointArray;
    /** Input array of MCTracks **/
    TClonesArray* fMcTrackArray;

    /** Output array of MpdEmcHit **/
    TClonesArray* fDigiArray;

    UInt_t GetSecId(Float_t x, Float_t y, Float_t z);
    UInt_t GetRowId(Float_t z);
    UInt_t GetSupModId(Float_t x, Float_t y, Float_t z, UInt_t sec);
    UInt_t GetModId(Float_t x, Float_t y, UInt_t supMod, UInt_t sec);


    Float_t CalcZCenter(UInt_t sec, UInt_t row, UInt_t mod);
    Float_t CalcPhiCenter(UInt_t sec, UInt_t supMod, UInt_t mod);

    MpdEmcHit* SearchHit(UInt_t sec, UInt_t supMod, UInt_t row, UInt_t mod);

    MpdEmcGeoPar* fGeoPar;

    ClassDef(MpdEmcHitProducer, 2);

};

#endif
