#ifndef MPD_PID_H
#define MPD_PID_H

#define MASS_MU 0.1057
#define MASS_PI 0.1396
#define MASS_PI2 0.0195
#define MASS_KA 0.4937
#define MASS_PR 0.9383
#define MASS_DE 1.876
#define MASS_TR 2.8094
#define MASS_HE3 1.4047
#define MASS_HE4 1.863

#define PDG_DEUTERON 1000010020
#define PDG_TRITON 1000010030
#define PDG_HE3 1000020030
#define PDG_HE4 1000020040

#include <TClonesArray.h>
#include <TCollection.h>
#include <TFriendElement.h>
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TObject.h>
#include <TROOT.h>
#include <TFitter.h>
#include <TLegend.h>
#include "TString.h"

using namespace std;

class MpdPid : public TObject
{
	public:
	
	MpdPid(){}; /// default ctor
	
	MpdPid(Double_t sigmaTof, Double_t sigmaEloss, Double_t sqrts, /// generators: "URQMD", "LAQGSM" ("QGSM"), "DEFAULT", "NSIG" (for n-sigma method)
		Double_t koef = 1., TString Generator = "DEFAULT", TString Tracking = "CF", /// tracking: "HP" (Hit Producer), "CF" (Cluster Finder)
		TString NSigPart = "pikapr"){Init(Generator, Tracking, NSigPart);}; /// possible expressions: el, mu, pi, ka, pr, de, tr, he3, he4
	
	virtual ~MpdPid(){} /// destructor
	
	/// Returns pos./neg. ratio for protons
	Double_t GetPrRat() {return prrat;}
	/// Set pos./neg. ratio for protons
	void SetPrRat(Double_t PrRat) {prrat=PrRat; parPrNegMom->SetParameter(0, (parPrPosMom->GetParameter(0) / prrat));}
	
	/// Polynomials for asymmetry description
	TF1 *fAsymmetryElLowP; TF1 *fAsymmetryElMidP; TF1 *fAsymmetryElHighP; /// electrons
	TF1 *fAsymmetryMuLowP; TF1 *fAsymmetryMuMidP; TF1 *fAsymmetryMuHighP; /// muons
	TF1 *fAsymmetryPiLowP; TF1 *fAsymmetryPiMidP; TF1 *fAsymmetryPiHighP; /// pions
	TF1 *fAsymmetryKaLowP; TF1 *fAsymmetryKaMidP; TF1 *fAsymmetryKaHighP; /// kaons
	TF1 *fAsymmetryPrLowP; TF1 *fAsymmetryPrMidP; TF1 *fAsymmetryPrHighP; /// protons
	TF1 *fAsymmetryDeLowP; TF1 *fAsymmetryDeHighP; /// deuterons
	TF1 *fAsymmetryTrLowP; TF1 *fAsymmetryTrMidP; TF1 *fAsymmetryTrHighP; /// tritons
	TF1 *fAsymmetryHe3LowP; TF1 *fAsymmetryHe3MidP; TF1 *fAsymmetryHe3HighP; /// he3
	TF1 *fAsymmetryHe4LowP; TF1 *fAsymmetryHe4MidP; TF1 *fAsymmetryHe4HighP; /// he4
	
	/// Bethe-Bloch functions for mean energy deposit description
	TF1 *parElBB; TF1 *parMuBB; TF1 *parPiBB; TF1 *parKaBB; TF1 *parPrBB; TF1 *parPrBBLowP;
	TF1 *parDeBB; TF1 *parTrBB; TF1 *parHe3BB; TF1 *parHe4BB;
	TF1 *parDePol1; TF1 *parDePol2; TF1 *parTrPol1; TF1 *parTrPol2;
	TF1 *parHe3Pol1; TF1 *parHe3Pol2; TF1 *parHe3Pol3; TF1 *parHe4Pol1; TF1 *parHe4Pol2;
	
	/// Functions for multiplicity description
	TF1 *parElNegMom; TF1 *parElPosMom; TF1 *parMuNegMom; TF1 *parMuPosMom;
	TF1 *parPiNegMom; TF1 *parPiPosMom; TF1 *parKaNegMom; TF1 *parKaPosMom;
	TF1 *parPrPosMom; TF1 *parPrNegMom;
	TF1 *parDeMom; TF1 *parTrMom; TF1 *parHe3Mom; TF1 *parHe4Mom;
	
	/// Functions for m2 width description
	TF1 *parElM2; TF1 *parMuM2; TF1 *parPiLowPM2; TF1 *parPiHighPM2;
	TF1 *parKaM2; TF1 *parPrLowPM2; TF1 *parPrHighPM2;
	TF1 *parDeM2; TF1 *parTrM2; TF1 *parHe3M2; TF1 *parHe4M2;
	
	/// Subsidiary functions for multiplicity description
	Double_t MomPi(Double_t *x, Double_t *par);
	Double_t MomPr(Double_t *x, Double_t *par);
	
	/// Fill array of probabilities, otherwise return kFALSE
	Bool_t FillProbs(Double_t, Double_t, Int_t); /// variables: full momentum, dE/dx, charge
	Bool_t FillProbs(Double_t, Double_t, Double_t, Int_t); /// variables: full momentum, dE/dx, mass squared, charge
	
	/// Return probabilities 
	Double_t GetProbPi(void){return fProbPi;}
	Double_t GetProbMu(void){return fProbMu;}
	Double_t GetProbPr(void){return fProbPr;}
	Double_t GetProbKa(void){return fProbKa;}
	Double_t GetProbEl(void){return fProbEl;}
	Double_t GetProbDe(void){return fProbDe;}
	Double_t GetProbTr(void){return fProbTr;}
	Double_t GetProbHe3(void){return fProbHe3;}
	Double_t GetProbHe4(void){return fProbHe4;}
	
	/// Returns the most probable pdg code
	Long_t GetMaxProb();
	
	/// Return expected dE/dx values (variable: full momentum)
	Double_t GetDedxPiParam(Double_t);
	Double_t GetDedxPrParam(Double_t);
	Double_t GetDedxKaParam(Double_t);
	Double_t GetDedxElParam(Double_t);
	Double_t GetDedxMuParam(Double_t);
	Double_t GetDedxDeParam(Double_t);
	Double_t GetDedxTrParam(Double_t);
	Double_t GetDedxHe3Param(Double_t);
	Double_t GetDedxHe4Param(Double_t);
	
	/// Returns expected asymmetry value delta (variables: full momentum, specie)
	Double_t GetTailValue(Double_t, Int_t);
	
	/// Returns expected dE/dx width (variables: full momentum, specie)
	Double_t GetDedxWidthValue(Double_t, Int_t);
	
	protected:
	
	TF1 *fgaus;
	TF2 *fgaus2;
	TF1 *fasymgaus;
	TF2 *fasymgaus2;
	
	/// dE/dx width
	TF1* elSigmaLowP; TF1* elSigmaMidP; TF1* elSigmaHighP;
	TF1* muSigmaLowP; TF1* muSigmaMidP; TF1* muSigmaHighP;
	TF1* piSigmaLowP; TF1* piSigmaMid1P; TF1* piSigmaMid2P; TF1* piSigmaHighP;
	TF1* prSigmaLowP; TF1* prSigmaHighP;
	TF1* kaSigmaLowP; TF1* kaSigmaHighP;
	TF1* deSigmaLowP; TF1* deSigmaHighP;
	TF1* trSigmaLowP; TF1* trSigmaMidP; TF1* trSigmaHighP;
	TF1* he3SigmaLowP; TF1* he3SigmaMidP; TF1* he3SigmaHighP;
	TF1* he4SigmaLowP; TF1* he4SigmaMid1P; TF1* he4SigmaMid2P; TF1* he4SigmaHighP;
	
	Double_t fProbEl;
	Double_t fProbMu;
	Double_t fProbPi;
	Double_t fProbKa;
	Double_t fProbPr;
	Double_t fProbDe;
	Double_t fProbTr;
	Double_t fProbHe3;
	Double_t fProbHe4;
	
	Double_t prrat; /// rat is pos./neg.
	Double_t fSigmaDedx_1, fSigmaDedx_2, fSigmaDedx_3;
	Double_t fKoef; /// scale of dedx
	
	Double_t kSigmaTof;
	Double_t kSigmaEloss;
	
	Bool_t fTracking;
	Int_t fCharge;
	Int_t Multiplicities[14];
	Double_t fEnergy;
	Double_t AsymGaus(Double_t*, Double_t*);
	Double_t AsymGaus2(Double_t*, Double_t*);
	
	void Init(TString, TString, TString);
	
	/// GetDedxProb_asym variables: cut_dedx, p, dedx, n, emean, sige, specie
	Double_t GetDedxProb_asym(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Int_t); 
	
	/// GetCombProb_asym variables: cut_dedx, cut_m2, p, dedx, m2, n, emean, mmean, sige, sigm, specie
	Double_t GetCombProb_asym(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Int_t);
	
	//ClassDef(MpdPid,3);
};

#endif
