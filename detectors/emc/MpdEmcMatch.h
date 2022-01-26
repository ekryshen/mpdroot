//--------------------------------------------------------------------
//
// Description:
//      MPD EMC track-to-rec.point match object
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 2-June-2016
//
//--------------------------------------------------------------------

#ifndef MPDEMCMATCH_H
#define MPDEMCMATCH_H 1

#include <TObject.h>

class MpdEmcMatch : public TObject
{

 public:

    /** Default constructor **/
    MpdEmcMatch();

    /** Standard constructor**/
    MpdEmcMatch(Int_t trackInd, Int_t recpInd, Double_t distT, Double_t distL, Double_t ener, Double_t trLeng);

    /** Destructor **/
    virtual ~MpdEmcMatch();

    /** Accessors **/
    Int_t GetTrackInd() const { return fInds[0]; }
    Int_t GetRecpInd() const { return fInds[1]; }
    Double_t GetChi2() const { return fChi2; }
    Double_t GetChi2pi() const { return fChi2pi; }
    Double_t GetEnergy() const { return fEnergy; }
    Double_t GetTrLeng() const { return fTrLeng; }
    Double_t GetDist(Int_t indx = 0) const { return fDist[indx]; }

    /** Modifiers **/
    void SetTrackInd(Int_t ind) { fInds[0] = ind; }
    void SetRecpInd(Int_t ind) { fInds[1] = ind; }
    void SetDist(Double_t dist, Int_t indx = 0) { fDist[indx] = dist; }
    void SetChi2(Double_t c2) { fChi2 = c2; }
    void SetChi2pi(Double_t c2) { fChi2pi = c2; }
    void SetEnergy(Double_t e) { fEnergy = e; }
    void SetTrLeng(Double_t leng) { fTrLeng = leng; } 

protected:

    Int_t fInds[2]; // track and rec. point indices
    Double32_t fDist[2]; // rec.point-to-track distance (across - along track)
    Double32_t fEnergy; // rec.point energy
    Double32_t fTrLeng; // extrap. track length in EMC
    Double32_t fChi2; // chi2 of match
    Double32_t fChi2pi; // chi2 of match for pi+-
    

    ClassDef(MpdEmcMatch, 1);

};

#endif
