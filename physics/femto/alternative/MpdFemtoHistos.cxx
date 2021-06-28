#include <TStyle.h>
#include <TH1.h>
#include "MpdFemtoHistos.h"

//--------------------------------------------------------------------------

MpdFemtoHistos::MpdFemtoHistos(Float_t qInv, const Char_t* out) {
    fQinv = qInv;
    fOut = new TFile(out, "RECREATE");
    const Int_t nBins = 100;
    _hCFQinvNomBase = new TH1F("_hCFQinvNomBase", "_hCFQinvNomBase", nBins, 0., qInv);
    _hCFQinvNom = new TH1F("_hCFQinvNom", "_hCFQinvNom", nBins, 0., qInv);
    _hCFQinvDenom = new TH1F("_hCFQinvDenom", "_hCFQinvDenom", nBins, 0., qInv);
    _hCF = new TH1F("_hCF", "_hCF", nBins, 0., qInv);
    _hCFBase = new TH1F("_hCFBase", "_hCFBase", nBins, 0., qInv);
    
    _hCFNom3D   = new TH3F("_hCFNom3D", "_hCFNom3D", 10, 0., qInv, 10, 0., qInv, 10, 0., qInv);
    _hCFDenom3D = new TH3F("_hCFDenom3D", "_hCFDenom3D", 10, 0., qInv, 10, 0., qInv, 10, 0., qInv);
    _hCF3D      = new TH3F("_hCF3D", "_hCF3D", 10, 0., qInv, 10, 0., qInv, 10, 0., qInv);
    
    const Int_t nBins1 = 30; 
    const Float_t xLow = -0.5;
    const Float_t xUp = -xLow;
    _hDeltaPhiDeltaEta = new TH2F("_hDeltaPhiDeltaEta", "_hDeltaPhiDeltaEta", nBins1, xLow, xUp, nBins1, xLow, xUp);
    _hDeltaPhiDeltaEtaNomin = new TH2F("_hDeltaPhiDeltaEtaNomin", "_hDeltaPhiDeltaEtaNomin", nBins1, xLow, xUp, nBins1, xLow, xUp);
    _hDeltaPhiDeltaEtaDenom = new TH2F("_hDeltaPhiDeltaEtaDenom", "_hDeltaPhiDeltaEtaDenom", nBins1, xLow, xUp, nBins1, xLow, xUp);
//    _hDeltaPhiDeltaEtaNominProjX = new TH1D("_hDeltaPhiDeltaEtaNominProjX", "_hDeltaPhiDeltaEtaNominProjX", 30, -1., 1.);
//    _hDeltaPhiDeltaEtaDenomProjX = new TH1D("_hDeltaPhiDeltaEtaDenomProjX", "_hDeltaPhiDeltaEtaDenomProjX", 30, -1., 1.);
    _hDeltaPhiDeltaEtaProjX = new TH1D("_hDeltaPhiDeltaEtaProjX", "_hDeltaPhiDeltaEtaProjX", nBins1, xLow, xUp);
    _hDeltaPhiDeltaEtaProjY = new TH1D("_hDeltaPhiDeltaEtaProjY", "_hDeltaPhiDeltaEtaProjY", nBins1, xLow, xUp);
    _hEtaPhiStar = new TH2F("_hEtaPhiStar", "_hEtaPhiStar", nBins, 0., 0., 100, 0., 0.);
    
    _hQuality = new TH1F("_hQuality", "_hQuality", 100, -0.6, 1.);
    _hSharing = new TH1F("_hSharing", "_hSharing", 100, 0.0, 0.3);
    _hQualityVsSharing = new TH2F("_hQualityVsSharing", "_hQualityVsSharing", 100, -0.6, 1., 100, 0.0, 0.3); 
    
    _hPtNoSplit = new TH1F("_hPtNoSplit", "_hPtNoSplit", 200, 0.2, 1.0);
    _hPtSplit = new TH1F("_hPtSplit", "_hPtSplit", 200, 0.2, 1.0);
    
    _hEff = new TH1F("_hEff", "_hEff", 200, 0.2, 1.0);
    
    _hNsplits = new TH1I("_hNsplits", "_hNsplits", 19, 1, 20);
    
    _hQualityVsNhits = new TH2F("_hQualityVsNhits", "_hQualityVsNhits", 100, -0.6, 1., 50, 4, 54);
    
}

//--------------------------------------------------------------------------
MpdFemtoHistos::~MpdFemtoHistos() {
    _hCFQinvNomBase->Write();
    _hCFQinvNom->Write();
    _hCFQinvDenom->Write();
    _hCF->Write();
    _hCFBase->Write();
    // _hCFNom3D->Write();
    // _hCFDenom3D->Write();
    // _hCF3D->Write();
    _hDeltaPhiDeltaEta->Write();
    _hQuality->Write();
    _hSharing->Write();
    _hQualityVsSharing->Write();
    _hPtNoSplit->Write();
    _hPtSplit->Write();
    _hEff->Write();
    _hNsplits->Write();
    _hQualityVsNhits->Write();
    _hDeltaPhiDeltaEtaNomin->Write();
    _hDeltaPhiDeltaEtaDenom->Write();
//    _hDeltaPhiDeltaEtaNominProjX->Write();
//    _hDeltaPhiDeltaEtaDenomProjX->Write();
    _hEtaPhiStar->Write();
    _hDeltaPhiDeltaEtaProjX->Write();
    _hDeltaPhiDeltaEtaProjY->Write();

    delete _hCFQinvNomBase;
    delete _hCFQinvNom;
    delete _hCFQinvDenom;
    delete _hCF;
    delete _hCFBase;
    delete _hCFNom3D;
    delete _hCFDenom3D;
    delete _hCF3D;
    delete _hDeltaPhiDeltaEta;
    delete _hDeltaPhiDeltaEtaNomin;
    delete _hDeltaPhiDeltaEtaDenom;
    delete _hQuality;
    delete _hSharing;
    delete _hQualityVsSharing;
    delete _hPtNoSplit;
    delete _hPtSplit;
    delete _hEff;
    delete _hNsplits;
    delete _hQualityVsNhits;
//    delete _hDeltaPhiDeltaEtaNominProjX;
//    delete _hDeltaPhiDeltaEtaDenomProjX;
    delete _hEtaPhiStar;
    delete _hDeltaPhiDeltaEtaProjX;
    delete _hDeltaPhiDeltaEtaProjY;
    
    delete fOut;
}

//--------------------------------------------------------------------------
Double_t* MpdFemtoHistos::GetFitParams1D() {
    _hDeltaPhiDeltaEtaNomin->Sumw2();
    _hDeltaPhiDeltaEtaDenom->Sumw2();
    _hDeltaPhiDeltaEta->Divide(_hDeltaPhiDeltaEtaNomin, _hDeltaPhiDeltaEtaDenom, 1., 1., "B");
    
    TH1D* nomProjX = _hDeltaPhiDeltaEtaNomin->ProjectionX("pxN",   _hDeltaPhiDeltaEtaNomin->GetXaxis()->GetFirst(), _hDeltaPhiDeltaEtaNomin->GetXaxis()->GetLast());
    TH1D* nomProjY = _hDeltaPhiDeltaEtaNomin->ProjectionY("pyN",   _hDeltaPhiDeltaEtaNomin->GetYaxis()->GetFirst(), _hDeltaPhiDeltaEtaNomin->GetYaxis()->GetLast());
    TH1D* denomProjX = _hDeltaPhiDeltaEtaDenom->ProjectionX("pxD", _hDeltaPhiDeltaEtaDenom->GetXaxis()->GetFirst(), _hDeltaPhiDeltaEtaDenom->GetXaxis()->GetLast());
    TH1D* denomProjY = _hDeltaPhiDeltaEtaDenom->ProjectionY("pyD", _hDeltaPhiDeltaEtaDenom->GetYaxis()->GetFirst(), _hDeltaPhiDeltaEtaDenom->GetYaxis()->GetLast());
    
    _hDeltaPhiDeltaEtaProjX->Divide(nomProjX, denomProjX, 1., 1., "B");
    _hDeltaPhiDeltaEtaProjY->Divide(nomProjY, denomProjY, 1., 1., "B");
      
    _hPtSplit->Sumw2();
    _hPtNoSplit->Sumw2();
    //_hPtSplit->Scale(1. / _hPtSplit->GetEntries());
    //_hPtNoSplit->Scale(1. / _hPtNoSplit->GetEntries());
    _hEff->Divide(_hPtSplit, _hPtNoSplit, 1., 1., "B");
    _hEff->Scale(100.);
    
    _hQuality->Scale(1. / _hQuality->GetEntries());
    _hSharing->Scale(1. / _hSharing->GetEntries());
    _hQualityVsSharing->Scale(1. / _hQualityVsSharing->GetEntries());
    _hQualityVsNhits->Scale(1. / _hQualityVsSharing->GetEntries());
        
     
    // Normalization:
    _hCFQinvNom->Sumw2();
    _hCFQinvNomBase->Sumw2();
    _hCFQinvDenom->Sumw2();

   // _hCF->Sumw2();
    _hCF->Divide(_hCFQinvNom, _hCFQinvDenom, 1., 1., "B");
    // _hCFBase->Sumw2();
    _hCFBase->Divide(_hCFQinvNomBase, _hCFQinvDenom, 1., 1., "B");

    Float_t normDenomFactor = _hCFQinvDenom->Integral(0.25 * _hCFQinvDenom->GetNbinsX(), 0.75 * _hCFQinvDenom->GetNbinsX());
    Float_t normNominFactor = _hCFQinvNom->Integral(0.25 * _hCFQinvNom->GetNbinsX(), 0.75 * _hCFQinvNom->GetNbinsX());
    Float_t normNominBaseFactor = _hCFQinvNomBase->Integral(0.25 * _hCFQinvNomBase->GetNbinsX(), 0.75 * _hCFQinvNomBase->GetNbinsX());

    if (normNominFactor > 0. && normNominBaseFactor > 0.) {
        _hCF->Scale(normDenomFactor / normNominFactor);
        _hCFBase->Scale(normDenomFactor / normNominBaseFactor);
    }
       
    //Fitting: fit QS only: C2 = N * [1 + lambda * exp(-q^2r0^2/h/h)], h=0.1973 GeV/fm
    TF1* fqs = new TF1("fqs", "[2] * (1 + [1] * exp(-x * x * [0] * [0] / 0.1973 / 0.1973))", 0, fQinv);
    fqs->SetParName(0, "r_{0}");
    fqs->SetParName(1, "#lambda");
    fqs->SetParName(2, "Norm");
    fqs->SetParameters(10., 1., 1.);

    _hCF->Fit(fqs, "SRQ", " ", 0., fQinv); 
      
    return fqs->GetParameters();
}

//--------------------------------------------------------------------------

Double_t* MpdFemtoHistos::GetFitParams3D() {
    _hCFNom3D->Sumw2();
    _hCFDenom3D->Sumw2();
    _hCF3D->Sumw2();
    _hCF3D->Divide(_hCFNom3D, _hCFDenom3D, 1., 1.);
 
    TF3* fitc = new TF3("fitc", "1 + [3] * exp(-25.76578 * (x * x * [0] * [0] + y * y * [1] * [1] + z * z * [2] * [2]))");
    fitc->SetParameters(3.0, 3.0, 3.0, 1.0);

    _hCF3D->Fit(fitc, "SRQ", " ", 0., fQinv);

    return fitc->GetParameters();
}
