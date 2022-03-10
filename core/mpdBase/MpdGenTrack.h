#include <TNamed.h>
#include <TLorentzVector.h>

#ifndef MPDGENTRACK_H
#define MPDGENTRACK_H

class MpdGenTrack : public TNamed  {
public:
    MpdGenTrack();
    ~MpdGenTrack();
    
    void SetXYZT(Double_t x, Double_t y, Double_t z, Double_t t) {
        fX = x;
        fY = y;
        fZ = z;
        fT = t;
    }
    
    TLorentzVector GetCoordinates() {
        TLorentzVector coords;
        coords.SetXYZT(fX, fY, fZ, fT);
        return coords;
    }
  
    void SetPxyz(Double_t px, Double_t py, Double_t pz) {
        fPx = px;
        fPy = py;
        fPz = pz;
    }
    
    TVector3 GetMomentum() {
        TVector3 mom;
        mom.SetXYZ(fPx, fPy, fPz);
        return mom;
    }
    
    void SetE(Double_t ene) {
        fEne = ene;    
    }
    
    void SetImpact(Double_t imp) {
        fImp = imp;
    }
    
    void SetPdg(Int_t pdg) {
        fPdg = pdg;  
    }

    void SetIsUsed(Bool_t flag) {
      isUsed = flag;
    }

    Bool_t GetIsUsed() {
      return isUsed;
    }
    
private:
    Double_t fImp;
    
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fT;
    
    Double_t fPx;
    Double_t fPy;
    Double_t fPz;
    
    Double_t fEne;
    
    Int_t fPdg;

    Bool_t isUsed;
    ClassDef(MpdGenTrack, 1);
};
#endif
