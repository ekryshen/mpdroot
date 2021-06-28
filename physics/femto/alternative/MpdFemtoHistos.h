#ifndef MPDFEMTOHISTOS_H
#define MPDFEMTOHISTOS_H 1

#include <iostream>
#include <TNamed.h>
#include <TH1.h>
#include <TF1.h>
#include <TF3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TVector3.h>
#include <TH3.h>
#include <TH2.h>
#include <TFile.h>

using namespace std;

class MpdFemtoHistos : public TNamed {
public:
    MpdFemtoHistos() {}
    MpdFemtoHistos(Float_t, const Char_t*);

    virtual ~MpdFemtoHistos();
    
    void SetQinv(Float_t val) {
        fQinv = val;
    }
    
    // Getters and Setters for 1D-analysis
    TH1F* GetNominator() {
        return _hCFQinvNom;
    }

    TH1F* GetNominatorBase() {
        return _hCFQinvNomBase;
    }

    TH1F* GetDenominator() {
        return _hCFQinvDenom;
    }
    
    TH1F* GetCF() {
        return _hCF;
    }
    
    TH1F* GetCFBase() {
        return _hCFBase;
    }
    
    TH2F* GetDeltaEtaDeltaPhi() {
        return _hDeltaPhiDeltaEta;
    }
   
    TH2F* GetDeltaEtaDeltaPhiNomin() {
        return _hDeltaPhiDeltaEtaNomin;
    }
    
//    TH1D* GetDeltaEtaDeltaPhiNominProjX() {
//        return _hDeltaPhiDeltaEtaNominProjX;
//    }
    
    TH2F* GetDeltaEtaDeltaPhiDenom() {
        return _hDeltaPhiDeltaEtaDenom;
    }
    
//    TH1D* GetDeltaEtaDeltaPhiDenomProjX() {
//        return _hDeltaPhiDeltaEtaDenomProjX;
//    }
    
    TH2F* GetEtaPhiStar() {
        return _hEtaPhiStar;
    }
    
    TH1F* GetQuality() {
        return _hQuality;
    }
    
    TH1F* GetSharing() {
        return _hSharing;
    }
    
    TH2F* GetQualityVsSharing() {
        return _hQualityVsSharing;
    }
    
    TH2F* GetQualityVsNhits() {
        return _hQualityVsNhits;
    }
    
    TH1F* GetPtNoSplit() {
        return _hPtNoSplit;
    }
    
    TH1F* GetPtSplit() {
        return _hPtSplit;
    }
    
    TH1F* GetEff() {
        return _hEff;
    }
    
    TH1I* GetNsplits() {
        return _hNsplits;
    }
    
    void SetNominator(TH1F* h) {
        _hCFQinvNom = h;
    }

    void SetNominatorBase(TH1F* h) {
        _hCFQinvNomBase = h;
    }

    void SetDenominator(TH1F* h) {
        _hCFQinvDenom = h;
    }
    
    void SetCF(TH1F* h) {
    _hCF = h;
    }

    void SetCFBase(TH1F* h) {
    _hCFBase = h;
    }
    
    // Getters and Setters for 3D-analysis
    TH3F* GetNominator3D() {
        return _hCFNom3D;
    }

    TH3F* GetDenominator3D() {
        return _hCFDenom3D;
    }
    
    TH3F* GetCF3D() {
        return _hCF3D;
    }
        
    void SetNominator3D(TH3F* h) {
        _hCFNom3D = h;
    }

    void SetDenominator3D(TH3F* h) {
        _hCFDenom3D = h;
    }
    
    void SetCF3D(TH3F* h) {
    _hCF3D = h;
    }
    
    void SetDeltaEtaDeltaPhi(TH2F* h) {
    _hDeltaPhiDeltaEta = h;
    }
    
    void SetDeltaEtaDeltaPhiNomin(TH2F* h) {
    _hDeltaPhiDeltaEtaNomin = h;
    }
    
//    void SetDeltaEtaDeltaPhiNominProjX(TH1D* h) {
//    _hDeltaPhiDeltaEtaNominProjX = h;
//    }
    
    void SetDeltaEtaDeltaPhiDenom(TH2F* h) {
    _hDeltaPhiDeltaEtaDenom = h;
    }
   
//    void SetDeltaEtaDeltaPhiDenomProjX(TH1D* h) {
//    _hDeltaPhiDeltaEtaDenomProjX = h;
//    }
    
    void SetEtaPhiStar(TH2F* h) {
    _hEtaPhiStar = h;
    }
    
    void SetQuality(TH1F* h) {
    _hQuality = h;
    } 
    
    void SetSharing(TH1F* h) {
    _hSharing = h;
    }
    
    void SetQualityVsSharing(TH2F* h) {
    _hQualityVsSharing = h;
    }
    
    void SetQualityVsNhits(TH2F* h) {
    _hQualityVsNhits = h;
    }
    
    void SetPtNoSplit(TH1F* h) {
    _hPtNoSplit = h;
    }
    
    void SetPtSplit(TH1F* h) {
    _hPtSplit = h;
    }
    
    void SetEff(TH1F* h) {
    _hEff = h;
    }
    
    void SetNsplits(TH1I* h) {
    _hNsplits = h;
    }
        
    Double_t* GetFitParams1D();
    Double_t* GetFitParams3D();
    
    inline void DebugInfo() {
        cout << "Service information: " << endl;
        cout << " _hCFQinvNomBase, #entries = " << _hCFQinvNomBase->GetEntries() << endl;    
        cout << " _hCFQinvNom, #entries = " << _hCFQinvNom->GetEntries() << endl;    
        cout << " _hCFQinvDenom, #entries = " << _hCFQinvDenom->GetEntries() << endl;    
        cout << " _hCF, #entries = " << _hCF->GetEntries() << endl; 
        cout << " _hCFBase, #entries = " << _hCFBase->GetEntries() << endl;
        
        cout << " _hCFNom3D, #entries = " << _hCFNom3D->GetEntries() << endl; 
        cout << " _hCFDenom3D, #entries = " << _hCFDenom3D->GetEntries() << endl; 
        cout << " _hCF3D, #entries = " << _hCF3D->GetEntries() << endl; 
        cout << " fBins = " << fBins << endl;
        cout << " fxUp = " << fxUp << endl;
        cout << " fQinv = " << fQinv << endl;
    }
       
private:
    TFile* fOut;
    
    Int_t fBins;
    Float_t fxUp;
    
    Float_t fQinv;
    
    TH1F* _hCFQinvNomBase;
    TH1F* _hCFQinvNom;
    TH1F* _hCFQinvDenom;
    TH1F* _hCF;
    TH1F* _hCFBase;
    
    TH3F* _hCFNom3D;
    TH3F* _hCFDenom3D;
    TH3F* _hCF3D;
    
    TH2F* _hDeltaPhiDeltaEta;
    TH2F* _hDeltaPhiDeltaEtaNomin;
    TH2F* _hDeltaPhiDeltaEtaDenom;
    
    //TH1D* _hDeltaPhiDeltaEtaNominProjX;
    //TH1D* _hDeltaPhiDeltaEtaDenomProjX;
    TH1D* _hDeltaPhiDeltaEtaProjX;
    TH1D* _hDeltaPhiDeltaEtaProjY;
    
    TH2F* _hEtaPhiStar;
    
    TH1F* _hQuality;
    TH1F* _hSharing;
    TH2F* _hQualityVsSharing;
    
    TH1F* _hPtNoSplit;
    TH1F* _hPtSplit;
    TH1F* _hEff;
    
    TH1I* _hNsplits;
    
    TH2F* _hQualityVsNhits; 
    
    ClassDef(MpdFemtoHistos, 1)
};

#endif