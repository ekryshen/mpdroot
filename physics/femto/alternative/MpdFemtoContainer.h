#ifndef MPDFEMTOCONTAINER_H
#define MPDFEMTOCONTAINER_H 1

#include <iostream>
#include <TNamed.h>
#include <TLorentzVector.h>

using namespace std;

class MpdFemtoContainer : public TNamed {
public:
    MpdFemtoContainer();
    MpdFemtoContainer(Int_t, TLorentzVector, TLorentzVector, Float_t, Float_t, Int_t);
    
    virtual ~MpdFemtoContainer();

    // Getters
    Int_t GetEventNumber() {
        return fEventNumber;
    }

    TLorentzVector Get4Momentum() {
        return fMo;
    }

    TLorentzVector Get4Coordinate() {
        return fCo;
    }
    
    Float_t GetPhi() {
        return fPhi;
    }
    
    Float_t GetTheta() {
        return fTheta;
    }
    
    Int_t GetTrackID() {
        return fTrackID;
    }

    // Setters
    void SetEventNumber(Int_t number) {
        fEventNumber = number;
    }

    void Set4Coordinate(TLorentzVector coord) {
        fCo = coord;
    }

    void Set4Momentum(TLorentzVector mom) {
        fMo = mom;
    }
    
    void SetPhi(Float_t val) {
        fPhi = val;
    }
        
    void SetTheta(Float_t val) {
        fTheta = val;
    }
    
    void SetTrackID(Int_t val) {
        fTrackID = val;
    }
    
private:

    Int_t fEventNumber;
    TLorentzVector fMo; // 4-momentum of particle
    TLorentzVector fCo; // 4-coordinate of particle
       
    Float_t fPhi;
    Float_t fTheta;
    Int_t   fTrackID;

    ClassDef(MpdFemtoContainer, 1)
};

#endif