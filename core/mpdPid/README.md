2021-06-15-README

MpdPid How To
=============

* * *

The following instruction provides the description how to use MpdPid class in analysis.

Table of content:
-----------------

1. [MpdPid class](#MpdPid)
   *   [1.1 Constructor with arguments](#MpdPid_1)
   *   [1.2 Probability vector](#MpdPid_2)
      *   [1.2.1 How to fill probability vector](#MpdPid_2_1)
      *   [1.2.2 How to get probabilities](#MpdPid_2_2)
      *   [1.2.3 Getting the most probable pdg value](#MpdPid_2_3)
   *   [1.3 dE/dx methods](#MpdPid_3)
   *   [1.4 Mass-squared methods](#MpdPid_4)
   *   [1.5 Particle multiplicity methods](#MpdPid_5)
   *   [1.6 Other methods](#MpdPid_6)
   *   [1.7 Example](#MpdPid_7)
2. [MpdPidQA class](#MpdPidQA)
   *   [2.1 Constructor with arguments](#MpdPidQA_1)
   *   [2.2 Usage](#MpdPidQA_2)
      *   [2.2.1 How to fill QA histograms](#MpdPidQA_2_1)
      *   [2.2.2 How to save results](#MpdPidQA_2_2)
      *   [2.2.3 Example](#MpdPidQA_2_3)
   *   [2.3 Result structure](#MpdPidQA_3)

<a name="MpdPid"></a>
1\. MpdPid class
----------------
<a name="MpdPid_1"></a>
### 1.1 Constructor with arguments

Constructor of MpdPid class:

```
MpdPid(Double_t sigM, Double_t sigE, Double_t E, Double_t C, TString generator, TString tracking, TString nSigPart);
```

1.  `sigM` ‒ non-zero distance from the average mass-squared value (in terms of standard deviations);
2.  `sigE` ‒ non-zero distance from the average dE/dx value (in terms of standard deviations);
3.  `E` ‒ the collision energy;
4.  `C` ‒ scale coefficient of dE/dx, should be used if dE/dx has been multiplied by this value during the reconstruction process (if in doubt, put 1);
5.  `generator` ‒ the model which has been used in simulation, possible expressions are `"LAQGSM"` (`"QGSM"`), `"EPOS"`, `"URQMD"`, `"PHSD"`, `"PHSD_CENT"` (PHSD central events, 0 < b < 3 fm @ 11 GeV), `"PHSD_CSR"`, `"PHSD_NOCSR"` (PHSD + Chiral Symmetry Restoration (CSR) mechanism on/off), `"NSIG"` (model-independent n-sigma method) and `"DEFAULT"` (the “average” value, set by default);
6.  `tracking` ‒ can be `"HP"` (Hit Producer), `"CF"` (Cluster Finder MLEM) and `"CFHM"` (Cluster Finder MLEM + HEED, used by default);
7.  `nSigPart` ‒ string of particles which are used in n-sigma method, possible expressions are `"el"`, `"mu"`, `"pi"`, `"ka"`, `"pr"`, `"de"`, `"tr"`, `"he3"`, `"he4"` or their combinations.

For example:

```
MpdPid *pid = new MpdPid(4.0, 4.0, 11.0, 1.0, "URQMD", "CF", "pikapr");
```

**IMPORTANT NOTE**:
In order to include particle in n-sigma method use the following method
```
void SetParticleEnabledNsig(MpdPidUtils::ePartType iType);
```
where `iType` is a type of particle which should be included (for example `MpdPidUtils::kProton`)

<a name="MpdPid_2"></a>
### 1.2 Probability vector

Probability vector is assigned to the track as PID result. Each element contains the probability to be one of the species.

<a name="MpdPid_2_1"></a>
#### 1.2.1 How to fill probability vector

To fill the probability vector, use `FillProbs` function for every particular track. It can be done with different arguments:

```
Bool_t MpdPid::FillProbs(MpdTrack *track);
Bool_t MpdPid::FillProbs(Double_t p, Double_t dedx, Int_t charge);
Bool_t MpdPid::FillProbs(Double_t p, Double_t dedx, Double_t m2, Int_t charge);

```

Function returns `kTRUE` if probability vector is successfully filled, otherwise `kFALSE`.

There is a TPC-TOF mismatch effect, the situation when track in TPC is associated with TOF hit produced by another particle. This effect is significant in low momenta. In this situation dE/dx value is reconstructed correctly, but m<sup>2</sup> value is far from expected. In that case PID cannot assign species to the track, thus `FillProbs` function return `kFALSE`. Nevertheless there is a possibility to fill probability vector using only dE/dx value and ignoring m<sup>2</sup> (e.g. using `FillProbs(p, dedx, charge)` function). The suggestion is to do it only for low-momenta particles.

<a name="MpdPid_2_2"></a>
#### 1.2.2 How to get probabilities

To get probabilities assigned to the species, use the following functions:
```
Double_t MpdPid::GetProbEl(void);
Double_t MpdPid::GetProbMu(void);
Double_t MpdPid::GetProbPi(void);
Double_t MpdPid::GetProbPr(void);
Double_t MpdPid::GetProbKa(void);
Double_t MpdPid::GetProbDe(void);
Double_t MpdPid::GetProbTr(void);
Double_t MpdPid::GetProbHe3(void);
Double_t MpdPid::GetProbHe4(void);
```
or
```
Double_t GetProb(MpdPidUtils::ePartType iType);
```
where `iType` is a type of particle which probability is returned.

Each function returns the probability normalized to 1.

<a name="MpdPid_2_3"></a>
#### 1.2.3 Getting the most probable pdg value

Function `Long_t MpdPid::GetMaxProb()` returns the most probable PDG-code (including charge information). It may be called when `MpdPid::FillProbs(...)` has already called.

<a name="MpdPid_3"></a>
### 1.3 dE/dx methods

1) In order to receive the expected <dE/dx>, width and asymmetry parameter (so-called *tail*), use the following methods:
```
Double_t GetDedxParam(Double_t p, MpdPidUtils::ePartType iType);
Double_t GetDedxWidthValue(Double_t p, MpdPidUtils::ePartType iType);
Double_t GetTailValue(Double_t p, MpdPidUtils::ePartType iType);
```
where p is the reconstructed full momentum and `iType` is a type of particle.

2) Parameterizations of <dE/dx>, width and tail in MpdPid are stored as vectors (std::vector) of pointers to 1-dimentional function (TF1*). This type is defined in MpdPid.h as **vecTF1ptrs**:
```
typedef vector<TF1*> vecTF1ptrs;
```

In order to receive the vector corresponded to the particular parameterization (<dE/dx>, width or tail) and particle species `iType`, use methods:
```
vecTF1ptrs GetVecdEdxMean(MpdPidUtils::ePartType iType);
vecTF1ptrs GetVecdEdxWidth(MpdPidUtils::ePartType iType);
vecTF1ptrs GetVecdEdxAsym(MpdPidUtils::ePartType iType);
```

<a name="MpdPid_4"></a>
### 1.4 m<sup>2</sup> methods

1) In order to receive the expected m<sup>2</sup> width, use method:
```
Double_t Getm2WidthParam(Double_t p, MpdPidUtils::ePartType iType);
```
where p is the reconstructed full momentum and `iType` is a type of particle.

2) Call the following method to receive the vector of m<sup>2</sup> width parameterization corresponded to the particle species `iType`:
```
vecTF1ptrs GetVecm2Width(MpdPidUtils::ePartType iType);
```

<a name="MpdPid_5"></a>
### 1.5 Particle multiplicity (yield) methods

1) In order to extract Bayes coefficients of particular track with reconstructed full momentum *p*, type *iType* and charge *ech*, use the method:
```
Double_t GetBayesCoefficient(Double_t p, MpdPidUtils::ePartType iType, MpdPidUtils::ePartCharge ech);
```

2) In order to improve PID quality user may set proton/antiproton ratio corresponding to the current data set before analysis. It can be done using the following function:

```
void MpdPid::SetPrRat(Double_t PrRat);
```
The corresponding getter is:
```
Double_t GetPrRat(void);
```

Default values of proton/antiproton ratio are `1000` (E < 7 GeV) and `100.0` (E > 7 GeV).

3) Call the following method to receive the vector of particle multiplicity parameterization corresponded to the particle species `iType`:
```
vecTF1ptrs GetVecYield(MpdPidUtils::ePartType iType);
```

<a name="MpdPid_6"></a>
#### 1.6 Other methods

1) Function `Double_t GetNsigmaToBetheBloch(MpdPidUtils::ePartType iType)` returns the distance from dE/dx value of the track to dE/dx value expected from Bethe-Bloch function in terms of sigmas. The input parameter is a type of particle.

Function `Double_t GetNsigmaToAverageMass2(MpdPidUtils::ePartType iType)` does the same but for m<sup>2</sup> value.

**IMPORTANT NOTE**: distances are calculating while `MpdPid::FillProbs` execution, so use `Double_t GetNsigmaToBetheBloch` and `Double_t GetNsigmaToAverageMass2` **AFTER** that!

2) Method
```
void SetParticleEnabledNsig(MpdPidUtils::ePartType iType);
```
includes particle `iType` in n-sigma method. Use it after MpdPid initialization with `"NSIG"` parameter to include particles consistently.

<a name="MpdPid_7"></a>
#### 1.7 Example

The following code illustrates how `MpdPid` can be used in external macro:
```
#include <MpdPid.h>

Double_t dedx, kfP, m2, pion, electron, kaon, proton, maxprob, fProbCut = 0.6;
Int_t charge, pidpdg;
Bool_t matchingExists, PIDresult;
/// inside of event- and kalmanTrack-loops
dedx = kalmanTrack→GetDedx(); /// CFHM
//dedx = kalmanTrack→GetPartID(); /// HP and CF
kfP = kalmanTrack→Momentum();
charge = kalmanTrack→Charge();
/// After cuts
if (matchingExists)
{
  m2 = mpdTOFMatchingData→GetMass2();
  PIDresult = pid->FillProbs(kfP, dedx, m2, charge);
  if (  (!PIDresult) && (kfP < 0.8) )
  {
    PIDresult = pid->FillProbs(kfP, dedx, charge);
    if (!PIDresult) continue;
  }
  pion = pid->GetProb(MpdPidUtils::kPion);
  electron = pid->GetProb(MpdPidUtils::kElectron);
  proton = pid→GetProb(MpdPidUtils::kProton);
  kaon = pid→GetProb(MpdPidUtils::kKaon);
  Double_t Probs[] = {0, electron, pion, kaon, proton};
  maxprob =  TMath::MaxElement(5, Probs);
  if (maxprob < fProbCut) continue;
  pidpdg = pid→GetMaxProb();
}
else
{
/// Do the same using FillProbs(kfP, dedx, charge)...
}
```

<a name="MpdPidQA"></a>
2\. MpdPidQA class
------------------

In order to tune `MpdPid` to the arbitrary data set, `MpdPidQA` class can be used.

<a name="MpdPidQA_1"></a>
### 2.1 Constructor with arguments

Constructor of `MpdPidQA` class has 7 arguments. They are the same as `MpdPid` constructor's arguments except of `nSigmaPart`.

`nSigPart` – string of particles which are used in n-sigma method and **string of particles, which will be checked for the PID quality.** Constructor

```
Double_t sigM = 4.0, sigE = 4.0, energy = 11.0, koef = 1.;
TString generator = “LAQGSM”, tracking = “CF”, nSigPart = “pikaprde”;
MpdPidQA *pidQA = new MpdPidQA(sigM, sigE, energy, koef, generator, tracking, nSigPart);
```

creates the element of `MpdPidQA` class to check PID quality for pions, kaons, protons and deuterons.

<a name="MpdPidQA_2"></a>
### 2.2 Usage
<a name="MpdPidQA_2_1"></a>
#### 2.2.1 How to fill QA histograms

`MpdPidQA` class consists of the set of histograms. Quality assurance can be achieved after these histograms have been filled. To fill histograms which are responsible to dE/dx QA use:

```
void MpdPidQA::FillDedxHists(Double_t p, Double_t dedx, Int_t realPDG);
```

To m<sup>2</sup> QA use:

```
void MpdPidQA::Fillm2Hists(Double_t p, Double_t m2, Int_t realPDG);
```

To check the particle yields use:

```
void MpdPidQA::FillAmplHists(Double_t p, Int_t realPDG);
```

And finally to fill efficiency and contamination histograms use:

```
Bool_t MpdPidQA::FillEffContHists(MpdTrack *track, Int_t realPDG, Double_t probCut);
Bool_t MpdPidQA::FillEffContHists(Double_t p, Double_t dedx, Int_t charge, Int_t realPDG, Double_t probCut);
Bool_t MpdPidQA::FillEffContHists(Double_t p, Double_t dedx, Double_t m2, Int_t charge, Int_t realPDG, Double_t probCut);
```

where `realPDG` is PDG code from Monte Carlo and `probCut` is the minimal value of probability for the particle identification (`probCut = 0` means that the particle species are defined as maximum probability).

`FillEffContHists` returns the same value as `FillProbs` function with the same parameters. But `FillEffContHists` returns `kFALSE` in case when the maximum probability is less than `probCut`.

These functions should be called in track loop for each particular track in event.

<a name="MpdPidQA_2_2"></a>
#### 2.2.2 How to save results

In order to save `MpdPidQA` results use the following functions:

```
void MpdPidQA::GetDedxQA("/path/to/some/folder/");
void MpdPidQA::Getm2QA("/path/to/some/folder/");
void MpdPidQA::GetAmplQA("/path/to/some/folder/");
void MpdPidQA::GetEffContQA("/path/to/some/folder/");
```

<a name="MpdPidQA_2_3"></a>
#### 2.2.3 Example

The following code illustrates how `MpdPidQA` can be used:

```
#include <MpdMCTrack.h>
#include <MpdPidQA.h>
#include <MpdVertex.h>
#include <MpdTofMatching.h>

MpdHelix MakeHelix(const MpdKalmanTrack *tr) {
	const Double_t F_CUR0 = 0.3 * 0.01 * 5 / 10;
    Double_t r = tr->GetPosNew();
    Double_t phi = tr->GetParam(0) / r;
    Double_t x = r * TMath::Cos(phi);
    Double_t y = r * TMath::Sin(phi);
    Double_t dip = tr->GetParam(3);
    Double_t cur = F_CUR0 * TMath::Abs (tr->GetParam(4));
    TVector3 o(x, y, tr->GetParam(1));
    Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
    MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
    return helix;
}

int main() {
	const Int_t nPart = 10;          /// pi+, pi-, K+, K-, p, anti-p, d, t, he3, he4 ( DON'T CHANGE! )

	/// customization
	Bool_t PIDQA_dEdx     =  kTRUE;  /// dE/dx VS p plots, ratio plots
	Bool_t PIDQA_m2       =  kTRUE;  /// m^2 VS p plots, sigma m^2 VS p plots
	Bool_t PIDQA_Yield    =  kFALSE; /// particle yield plots
	Bool_t PIDQA_EffCont  =  kFALSE; /// Efficiency and Contamination plots

	const Double_t ENERGY = 8.8;     /// [GeV] sqrt(s_NN)
	Int_t nMaxEvents = 0;            /// max number of events, put 0 to execute all events
	const Double_t ImpParMax = -1.0; /// [fm] max impact parameter, put "-1.0" for mb
	const Double_t VzMax = 50.0;     /// [cm] max |Z_vertex| coordinate, put "-1.0" to exclude cut
	const Double_t AbsEtaMax = 1.3;  /// max |eta| value (RECO loop)
	const Int_t RECOmID = 0;         /// use MC mother ID in RECO track loop, put "1" to select primary tracks only, "2" to secondary tracks only, "0" to exclude cut
	const Int_t nHitsMin = 20;       /// min N_hits in TPC track
	const Int_t PIDmode = 1;         /// 0 - dE/dx PID only; 1 - comb PID only, 2 - comb + dE/dx PID
	const Bool_t DCAmode = 1;        /// put "1" to use cut or "0" to exclude cut
	const Double_t fProbCutArr[nPart] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; /// min value of PID probability
	const Double_t DCAMaxArr[nPart]   = { 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0 }; /// [cm] max DCA value in RECO track loop (uses if cut is switched on)

	const Double_t sigM = 3.0;       /// n-sigma band along m^2 for PID selection
	const Double_t sigE = 3.0;       /// n-sigma band along dE/dx for PID selection
	const Double_t coef = 1.0;       /// dE/dx scale factor (normally is 1.0)
	TString Generator = "PHQMD";
	TString Tracking = "CFHM";
	TString inname = "/path/to/data/*.root";
	TString outpath = "/outpath/phqmd-8.8GeV-mb-50cmVz-1.6AbsEta-primMC_DCA-20nH-PID1.raw/";
	/// end of customization

	TChain chain("mpdsim");
	chain.Add(inname);
	Int_t nev = chain.GetEntries();
	nMaxEvents = nMaxEvents == 0 ? nev : nMaxEvents;

	TClonesArray *mcTracks = NULL;
	chain.SetBranchAddress("MCTrack", &mcTracks);
	FairMCEventHeader *fmcHeader = NULL;
	chain.SetBranchAddress("MCEventHeader.", &fmcHeader);

	TClonesArray *tofMatches = NULL;
	chain.SetBranchAddress("TOFMatching", &tofMatches);
	TClonesArray *kalmanTracks = NULL;
	chain.SetBranchAddress("TpcKalmanTrack", &kalmanTracks);
	TClonesArray *vertexes = NULL;
	chain.SetBranchAddress("Vertex", &vertexes);

	MpdMCTrack *pMCTrack = new MpdMCTrack;
	MpdKalmanTrack *pKFTrack = new MpdKalmanTrack;
	MpdVertex *pVertex = new MpdVertex;
	MpdTofMatchingData *match = new MpdTofMatchingData;
	map<Int_t,Int_t> mapTof;

	MpdPid *pid = new MpdPid(sigM, sigE, ENERGY, coef, Generator, Tracking, "pikaprdetrhe3he4");
	MpdPidQA *pidQA = new MpdPidQA(sigM, sigE, ENERGY, coef, Generator, Tracking, "pikaprdetrhe3he4");

	Int_t mcNtracks, nKalmanTracks, nTofMatch, mcID, pdg, nHits, charge, kfID, PIDPDG;
	Double_t kfP, dedx, Theta, AbsEta, m2, fProbCut = 0.0;
	Bool_t ret, pidFlag;
	TVector3 pca, primaryVertex;

	for ( Int_t iEvent = 0; iEvent < nMaxEvents; iEvent++ ) {
		cout << iEvent+1 << "/" << nMaxEvents << " event is processed... \r" << flush;
		chain.GetEntry(iEvent); /// current entry
		if ( ImpParMax != -1.0 ) { if ( fmcHeader->GetB() > ImpParMax ) continue; }
		if ( VzMax != -1.0 ) { if ( TMath::Abs( fmcHeader->GetZ() ) > VzMax ) continue; }

		mcNtracks = mcTracks->GetEntries();
		nKalmanTracks = kalmanTracks->GetEntries();
		nTofMatch = tofMatches->GetEntriesFast();

		pVertex = (MpdVertex*) vertexes->First(); pVertex->Position(primaryVertex);

		for ( Int_t itof = 0; itof < nTofMatch; ++itof ) {
			match = (MpdTofMatchingData*) tofMatches->UncheckedAt(itof);
			mapTof[match->GetKFTrackIndex()] = itof;
		}

		for ( Int_t iTrack = 0; iTrack < nKalmanTracks; iTrack++ ) {
			pKFTrack = (MpdKalmanTrack*) kalmanTracks -> UncheckedAt(iTrack);
			mcID = pKFTrack->GetTrackID();
			pMCTrack = (MpdMCTrack*) mcTracks -> UncheckedAt(mcID);

			Theta = TMath::PiOver2() - pKFTrack->GetParam(3);
			AbsEta = TMath::Abs( -TMath::Log(TMath::Tan(0.5*Theta)) );
			nHits = pKFTrack->GetNofHits();

			if ( RECOmID == 1 && pMCTrack->GetMotherId() != -1 ) continue;
			if ( RECOmID == 2 && pMCTrack->GetMotherId() == -1 ) continue;
			if ( AbsEta > AbsEtaMax ) continue;
			if ( nHits < nHitsMin ) continue;

			kfP = pKFTrack->Momentum();
			pdg = pMCTrack->GetPdgCode();
			dedx = pKFTrack->GetDedx();
			charge = pKFTrack->Charge();

			if ( PIDQA_dEdx )  pidQA->FillDedxHists(kfP, dedx, pdg);
			if ( PIDQA_Yield ) pidQA->FillAmplHists(kfP, pdg);

			ret = kFALSE;
			if (mapTof.count(iTrack) > 0) {
				m2 = ((MpdTofMatchingData*)tofMatches->UncheckedAt(mapTof[iTrack]))->GetMass2();
				pidFlag = kTRUE;
				if ( PIDQA_m2 )      pidQA->Fillm2Hists(kfP, m2, pdg);
			}

			if ( PIDQA_EffCont ) {
				/// Apply PID
				if ( PIDmode == 0 ) ret = pid->FillProbs(kfP, dedx, charge); /// dE/dx PID only
				else if ( PIDmode == 1 ) { /// Comb PID only
					if ( pidFlag ) ret = pid->FillProbs(kfP, dedx, m2, charge);
				} else if ( PIDmode == 2 ) { /// Comb + dE/dx PID
					if ( pidFlag ) {
						ret = pid->FillProbs(kfP, dedx, m2, charge);
						if ( !ret && kfP < 0.8 ) ret = pid->FillProbs(kfP, dedx, charge);
					} else ret = pid->FillProbs(kfP, dedx, charge);
				} else continue; /// Wrong PIDmode expression

				if ( ret ) {
					PIDPDG = pid->GetMaxProb();
					switch ( PIDPDG ) {
						case 211:            fProbCut = fProbCutArr[0]; break;
						case -211:           fProbCut = fProbCutArr[1]; break;
						case 321:            fProbCut = fProbCutArr[2]; break;
						case -321:           fProbCut = fProbCutArr[3]; break;
						case 2212:           fProbCut = fProbCutArr[4]; break;
						case -2212:          fProbCut = fProbCutArr[5]; break;
						case PDG_DEUTERON:   fProbCut = fProbCutArr[6]; break;
						case PDG_TRITON:     fProbCut = fProbCutArr[7]; break;
						case PDG_HE3:        fProbCut = fProbCutArr[8]; break;
						case PDG_HE4:        fProbCut = fProbCutArr[9]; break;
						default:             fProbCut = 0.0; break;
					}
				}
				else fProbCut = 0.0;

				if ( PIDmode == 0 ) pidQA->FillEffContHists(kfP, dedx, charge, pdg, fProbCut); /// dE/dx PID only
				else if ( PIDmode == 1 ) { /// Comb PID only
					if ( pidFlag ) pidQA->FillEffContHists(kfP, dedx, m2, charge, pdg, fProbCut);
					else continue;
				} else if ( PIDmode == 2 ) { /// Comb + dE/dx PID
					if ( pidFlag ) {
						ret = pid->FillProbs(kfP, dedx, m2, charge);
						if ( !ret && kfP < 0.8 ) pidQA->FillEffContHists(kfP, dedx, charge, pdg, fProbCut);
						else pidQA->FillEffContHists(kfP, dedx, m2, charge, pdg, fProbCut);
					} else pidQA->FillEffContHists(kfP, dedx, charge, pdg, fProbCut);
				} else continue; /// Wrong PIDmode expression
			}
		}
		mapTof.clear();
	}
	if ( PIDQA_dEdx )     pidQA->GetDedxQA(outpath);
	if ( PIDQA_m2 )       pidQA->Getm2QA(outpath);
	if ( PIDQA_Yield )    pidQA->GetAmplQA(outpath);
	if ( PIDQA_EffCont )  pidQA->GetEffContQA(outpath);

	cout << endl << "Macro is finished successfully" << endl;
	return(0);
}
```

<a name="MpdPidQA_3"></a>
### 2.3 Result structure

`GetDedxQA` function creates *dEdXHists.root* file in the choosen directory. It consists of:

1) 1-dimentional dE/dx histograms of choosen particle species in several bins of the full momentum ( *Pi_0, Pi_1*... )

2) 2-dimentional dE/dx VS. p histograms of choosen particle species ( *hPiBB_QA, hKaBB_QA*... )

3) 2-dimentional dE/dx VS. p histograms - sum for all choosen particles ( *fSumBetheBlochHists* and *fChBetheBlochHists* where negative charged particles have negative value of *p* )

4) Graphs illustrated the ratio of dE/dx value in asymmetric gaussian peak over dE/dx value expected from Bethe-Bloch function ( used for estimating the parameterization quality: the closer the ratio to one the better the Bethe-Bloch description. Error bars show dE/dx resolution )

`Getm2QA` function creates *m2Hists.root* file in the choosen directory. It consists of:

1) 2-dimentional m<sup>2</sup> VS. p histograms of choosen particle species ( *mass2Pi, mass2Ka*... )

2) 2-dimentional m<sup>2</sup> VS. p histograms for hadrons ( *fm2LightHist* for pions, kaons and protons ) and light nuclei ( *fm2HeavyHist* for deuterons, tritons, he-3 and he-4 )

`GetAmplQA` function creates *ParticleYields.root* file in the choosen directory. It consists of 1-dimentional distributions of the choosen particles yields VS. p ( *hYield_pipos, hYield_pineg*... )

`GetEffContQA` function creates *EffCont.root* file in the choosen directory. It consists of:

1) PID Efficiency VS. p histograms of choosen particle species ( *IdRight_poschar_pi, IdRight_negchar_K*... )

2) PID Contamination VS. p histograms of choosen particle species ( *IdWrong_poschar_p, IdWrong_d*... )

In case when some parameterizations are not quite good for user, it can be corrected by the hand. All parameterizations are defined in `MpdPid::Init` function.
