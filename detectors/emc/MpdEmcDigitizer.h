//--------------------------------------------------------------------
//
// Description:
//      MPD EMC Digitizer - takes EmcPoints and makes digits
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-May-2016
//      Alexander Zinchenko LHEP, JINR, Dubna - 24-June-2018 - adapted to projective geometry
//
//--------------------------------------------------------------------

#ifndef MPDEMCDIGITIZER_H
#define MPDEMCDIGITIZER_H 1

#include "MpdEmcDigit.h"

#include "FairTask.h"
#include <iostream>
#include <map>

class MpdEmcGeoParams;
class TClonesArray;

class MpdEmcDigitizer : public FairTask {
public:

    /** Default constructor **/
    MpdEmcDigitizer();


    /** Destructor **/
    ~MpdEmcDigitizer();


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

    void RedoId(TClonesArray *digis, TClonesArray *mctrs);

    void FindChanPhiZ(Double_t &phi, Double_t &z);
    void FindChanPhiThe(Double_t &phi, Double_t &the);

    MpdEmcDigit* SearchHit(TString tower);

    MpdEmcGeoParams* fGeoPar;
    std::map<TString,MpdEmcDigit*> fHitMap; //! 

    ClassDef(MpdEmcDigitizer, 1);

};

#endif
