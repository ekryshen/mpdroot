//
// Track cut monitor for basic analysis
//

// MpdFemtoMaker headers
#include "MpdFemtoTrack.h"
#include "phys_constants.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoFxtTrackCutMonitor.h"

// C++ headers
#include <cstdio>
#include <string>

// ROOT headers
#include "TMath.h"

ClassImp(MpdFemtoFxtTrackCutMonitor);

//_________________
MpdFemtoFxtTrackCutMonitor::MpdFemtoFxtTrackCutMonitor(const char* name,
						       const double ParticleMass) :
MpdFemtoBaseCutMonitor(),
mDCAGlobal(nullptr),
mNhits(nullptr),
mP(nullptr),
mPt(nullptr),
mPtVsNsigmaPion(nullptr),
mPtVsNsigmaKaon(nullptr),
mPtVsNsigmaProton(nullptr),
mPvsDedx(nullptr),
mRapidity(nullptr),
mPseudoRapidity(nullptr),
mPvsMassSqr(nullptr),
mPvsInvBeta(nullptr),
mPtVsEta(nullptr),
mInvBetaDiffPionVsP(nullptr),
mInvBetaDiffKaonVsP(nullptr),
mInvBetaDiffProtonVsP(nullptr),
monMass(ParticleMass) {
  // Constructor
  string s("MpdFemtoFxtTrackCutMonitor");
  string n(name);
  mDCAGlobal = new TH1F((s + n + "mDCAGlobal").c_str(),
			"DCA Global; DCA Global (cm)",
			100, 0., 5.);
  mNhits = new TH1F((s + n + "mNhits").c_str(), "nHits;nHits",
		    80, -0.5, 79.5);
  mP = new TH1F((s + n + "mP").c_str(), "Momentum;p (GeV/c);Entries",
		76, 0., 1.9);
  mPt = new TH1F((s + n + "mPt").c_str(), "Transverse momentum;p_{T} (GeV/c);Entries",
		 76, 0., 1.9);
  mPtVsNsigmaPion = new TH2F((s + n + "mPtVsNsigmaPion").c_str(),
			     "n#sigma(#pi) vs. p_{T};p_{T} (GeV/c);n#sigma(#pi)",
			     76, 0., 1.9, 100, -10., 10.);
  mPtVsNsigmaKaon = new TH2F((s + n + "mPtVsNsigmaKaon").c_str(),
			     "n#sigma(K) vs. p_{T};p_{T} (GeV/c);n#sigma(K)",
			     76, 0., 1.9, 100, -10., 10.);
  mPtVsNsigmaProton = new TH2F((s + n + "mPtVsNsigmaProton").c_str(),
			       "n#sigma(p) vs. p_{T};p_{T} (GeV/c);n#sigma(p)",
			       76, 0., 1.9, 100, -10., 10.);
  mPvsDedx = new TH2F((s + n + "mPvsDedx").c_str(),
		      "dE/dx vs. momentum;p (GeV/c);dE/dx (keV/cm)",
		      76, 0.0, 1.9, 100, 0., 15000.);
  mRapidity = new TH1F((s + n + "mRapidity").c_str(), "Rapidity;y;Entries",
		       56, -1.2, 1.2);
  mPseudoRapidity = new TH1F((s + n + "mPseudoRapidity").c_str(), "Pseudorapidity;#eta;Entries",
			     56, -1.2, 1.2);
  mPvsMassSqr = new TH2F((s + n + "mPvsMassSqr").c_str(),
			 "m^{2} vs. p_{T};p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",
			 76, 0., 1.9, 120, -0.1, 1.1);
  mPvsInvBeta = new TH2F((s + n + "mPvsInvBeta").c_str(),
			 "1/#beta vs. p;p (GeV/c);1/#beta",
			 76, 0., 1.9, 60, 0.8, 2.);
  mPtVsEta = new TH2F((s + n + "mPtVsEta").c_str(), "#eta vs. p_{T};#eta;p_{T} (GeV/c)",
		      44, -1.1, 1.1, 76, 0., 1.9);
  mInvBetaDiffPionVsP = new TH2F((s + n + "mInvBetaDiffPionVsP").c_str(),
				 "1/#beta-1/#beta(#pi) vs. p;p (GeV/c);1/#beta-1/#beta(#pi)",
				 76, 0., 1.9, 200, -0.1, 0.1);
  mInvBetaDiffKaonVsP = new TH2F((s + n + "mInvBetaDiffKaonVsP").c_str(),
				 "1/#beta-1/#beta(K) vs. p;p (GeV/c);1/#beta-1/#beta(K)",
				 76, 0., 1.9, 200, -0.1, 0.1);
  mInvBetaDiffProtonVsP = new TH2F((s + n + "mInvBetaDiffProtonVsP").c_str(),
				   "1/#beta-1/#beta(p) vs. p;p (GeV/c);1/#beta-1/#beta(p)",
				   76, 0., 1.9, 200, -0.1, 0.1);
}

//_________________
MpdFemtoFxtTrackCutMonitor::MpdFemtoFxtTrackCutMonitor(const MpdFemtoFxtTrackCutMonitor& cutMoni) :
MpdFemtoBaseCutMonitor(cutMoni) {
  // Copy constructor
  if (mDCAGlobal) delete mDCAGlobal;
  mDCAGlobal = new TH1F(*(cutMoni.mDCAGlobal));
  if (mNhits) delete mNhits;
  mNhits = new TH1F(*(cutMoni.mNhits));
  if (mP) delete mP;
  mP = new TH1F(*(cutMoni.mP));
  if (mPt) delete mPt;
  mPt = new TH1F(*(cutMoni.mPt));
  if (mPtVsNsigmaPion) delete mPtVsNsigmaPion;
  mPtVsNsigmaPion = new TH2F(*(cutMoni.mPtVsNsigmaPion));
  if (mPtVsNsigmaKaon) delete mPtVsNsigmaKaon;
  mPtVsNsigmaKaon = new TH2F(*(cutMoni.mPtVsNsigmaKaon));
  if (mPtVsNsigmaProton) delete mPtVsNsigmaProton;
  mPtVsNsigmaProton = new TH2F(*(cutMoni.mPtVsNsigmaProton));
  if (mPvsDedx) delete mPvsDedx;
  mPvsDedx = new TH2F(*(cutMoni.mPvsDedx));
  if (mRapidity) delete mRapidity;
  mRapidity = new TH1F(*(cutMoni.mRapidity));
  if (mPseudoRapidity) delete mPseudoRapidity;
  mPseudoRapidity = new TH1F(*(cutMoni.mPseudoRapidity));
  if (mPvsMassSqr) delete mPvsMassSqr;
  mPvsMassSqr = new TH2F(*(cutMoni.mPvsMassSqr));
  if (mPvsInvBeta) delete mPvsInvBeta;
  mPvsInvBeta = new TH2F(*(cutMoni.mPvsInvBeta));
  if (mPtVsEta) delete mPtVsEta;
  mPtVsEta = new TH2F(*(cutMoni.mPtVsEta));
  if (mInvBetaDiffPionVsP) delete mInvBetaDiffPionVsP;
  mInvBetaDiffPionVsP = new TH2F(*(cutMoni.mInvBetaDiffPionVsP));
  if (mInvBetaDiffKaonVsP) delete mInvBetaDiffKaonVsP;
  mInvBetaDiffKaonVsP = new TH2F(*(cutMoni.mInvBetaDiffKaonVsP));
  if (mInvBetaDiffProtonVsP) delete mInvBetaDiffProtonVsP;
  mInvBetaDiffProtonVsP = new TH2F(*(cutMoni.mInvBetaDiffProtonVsP));
  monMass = cutMoni.monMass;
}

//_________________
MpdFemtoFxtTrackCutMonitor MpdFemtoFxtTrackCutMonitor::operator=(const MpdFemtoFxtTrackCutMonitor& cutMoni) {
  // Assignment operator
  if (this != &cutMoni) {
    if (mDCAGlobal) delete mDCAGlobal;
    mDCAGlobal = new TH1F(*(cutMoni.mDCAGlobal));
    if (mNhits) delete mNhits;
    mNhits = new TH1F(*(cutMoni.mNhits));
    if (mP) delete mP;
    mP = new TH1F(*(cutMoni.mP));
    if (mPt) delete mPt;
    mPt = new TH1F(*(cutMoni.mPt));
    if (mPtVsNsigmaPion) delete mPtVsNsigmaPion;
    mPtVsNsigmaPion = new TH2F(*(cutMoni.mPtVsNsigmaPion));
    if (mPtVsNsigmaKaon) delete mPtVsNsigmaKaon;
    mPtVsNsigmaKaon = new TH2F(*(cutMoni.mPtVsNsigmaKaon));
    if (mPtVsNsigmaProton) delete mPtVsNsigmaProton;
    mPtVsNsigmaProton = new TH2F(*(cutMoni.mPtVsNsigmaProton));
    if (mPvsDedx) delete mPvsDedx;
    mPvsDedx = new TH2F(*(cutMoni.mPvsDedx));
    if (mRapidity) delete mRapidity;
    mRapidity = new TH1F(*(cutMoni.mRapidity));
    if (mPseudoRapidity) delete mPseudoRapidity;
    mPseudoRapidity = new TH1F(*(cutMoni.mPseudoRapidity));
    if (mPvsMassSqr) delete mPvsMassSqr;
    mPvsMassSqr = new TH2F(*(cutMoni.mPvsMassSqr));
    if (mPvsInvBeta) delete mPvsInvBeta;
    mPvsInvBeta = new TH2F(*(cutMoni.mPvsInvBeta));
    if (mPtVsEta) delete mPtVsEta;
    mPtVsEta = new TH2F(*(cutMoni.mPtVsEta));
    if (mInvBetaDiffPionVsP) delete mInvBetaDiffPionVsP;
    mInvBetaDiffPionVsP = new TH2F(*(cutMoni.mInvBetaDiffPionVsP));
    if (mInvBetaDiffKaonVsP) delete mInvBetaDiffKaonVsP;
    mInvBetaDiffKaonVsP = new TH2F(*(cutMoni.mInvBetaDiffKaonVsP));
    if (mInvBetaDiffProtonVsP) delete mInvBetaDiffProtonVsP;
    mInvBetaDiffProtonVsP = new TH2F(*(cutMoni.mInvBetaDiffProtonVsP));
    monMass = cutMoni.monMass;
  } // if (this != &c)

  return *this;
}

//_________________
MpdFemtoFxtTrackCutMonitor::~MpdFemtoFxtTrackCutMonitor() {
  if (mDCAGlobal) {
    delete mDCAGlobal;
    mDCAGlobal = nullptr;
  }
  if (mNhits) {
    delete mNhits;
    mNhits = nullptr;
  }
  if (mP) {
    delete mP;
    mP = nullptr;
  }
  if (mPt) {
    delete mPt;
    mPt = nullptr;
  }
  if (mPtVsNsigmaPion) {
    delete mPtVsNsigmaPion;
    mPtVsNsigmaPion = nullptr;
  }
  if (mPtVsNsigmaKaon) {
    delete mPtVsNsigmaKaon;
    mPtVsNsigmaKaon = nullptr;
  }
  if (mPtVsNsigmaProton) {
    delete mPtVsNsigmaProton;
    mPtVsNsigmaProton = nullptr;
  }
  if (mPvsDedx) {
    delete mPvsDedx;
    mPvsDedx = nullptr;
  }
  if (mRapidity) {
    delete mRapidity;
    mRapidity = nullptr;
  }
  if (mPseudoRapidity) {
    delete mPseudoRapidity;
    mPseudoRapidity = nullptr;
  }
  if (mPvsMassSqr) {
    delete mPvsMassSqr;
    mPvsMassSqr = nullptr;
  }
  if (mPvsInvBeta) {
    delete mPvsInvBeta;
    mPvsInvBeta = nullptr;
  }
  if (mPtVsEta) {
    delete mPtVsEta;
    mPtVsEta = nullptr;
  }
  if (mInvBetaDiffPionVsP) {
    delete mInvBetaDiffPionVsP;
    mInvBetaDiffPionVsP = nullptr;
  }
  if (mInvBetaDiffKaonVsP) {
    delete mInvBetaDiffKaonVsP;
    mInvBetaDiffKaonVsP = nullptr;
  }
  if (mInvBetaDiffProtonVsP) {
    delete mInvBetaDiffProtonVsP;
    mInvBetaDiffProtonVsP = nullptr;
  }
}

//_________________
void MpdFemtoFxtTrackCutMonitor::fill(const MpdFemtoTrack* track) {
  /*
    float TPhi = atan2(track->P().y(),track->P().x());
    if(TPhi<0.0) { TPhi += 2.*TMath::Pi(); }
    TPhi *= 180./TMath::Pi();
  */

  mDCAGlobal->Fill(track->gDCA().Mag());
  mNhits->Fill(track->nHits());
  mP->Fill(track->p().Mag());
  mPt->Fill(track->pt());
  mPtVsNsigmaPion->Fill(track->pt(), track->nSigmaPion());
  mPtVsNsigmaKaon->Fill(track->pt(), track->nSigmaKaon());
  mPtVsNsigmaProton->Fill(track->pt(), track->nSigmaProton());
  mPvsDedx->Fill(track->ptot(), track->dEdxInKeV());
  mPseudoRapidity->Fill(track->eta());
  mPvsMassSqr->Fill(track->ptot(), track->massSqr());
  mPvsInvBeta->Fill(track->ptot(), (1. / track->beta()));
  mPtVsEta->Fill(track->eta(), track->pt());
  mInvBetaDiffPionVsP->Fill(track->ptot(),
			    track->invBeta() -
			    TMath::Sqrt(M_PION_PLUS * M_PION_PLUS +
					track->p().Mag2()) / track->ptot());
  mInvBetaDiffKaonVsP->Fill(track->ptot(),
			    track->invBeta() -
			    TMath::Sqrt(M_KAON_PLUS * M_KAON_PLUS +
					track->p().Mag2()) / track->ptot());
  mInvBetaDiffProtonVsP->Fill(track->ptot(),
			      track->invBeta() -
			      TMath::Sqrt(M_PROTON * M_PROTON +
					  track->p().Mag2()) / track->ptot());
}

//_________________
void MpdFemtoFxtTrackCutMonitor::writeOutHistos() {
  // Write all histograms
  mDCAGlobal->Write();
  mNhits->Write();
  mP->Write();
  mPt->Write();
  mPtVsNsigmaPion->Write();
  mPtVsNsigmaKaon->Write();
  mPtVsNsigmaProton->Write();
  mPvsDedx->Write();
  mRapidity->Write();
  mPseudoRapidity->Write();
  mPvsMassSqr->Write();
  mPvsInvBeta->Write();
  mPtVsEta->Write();
  mInvBetaDiffPionVsP->Write();
  mInvBetaDiffKaonVsP->Write();
  mInvBetaDiffProtonVsP->Write();
}

//_________________
TList* MpdFemtoFxtTrackCutMonitor::getOutputList() {
  // Prepare ouput list with histograms
  TList *outputList = new TList();

  outputList->Add(mDCAGlobal);
  outputList->Add(mNhits);
  outputList->Add(mP);
  outputList->Add(mPt);
  outputList->Add(mPtVsNsigmaPion);
  outputList->Add(mPvsDedx);
  outputList->Add(mRapidity);
  outputList->Add(mPseudoRapidity);
  outputList->Add(mPvsMassSqr);
  outputList->Add(mPvsInvBeta);
  outputList->Add(mPtVsEta);
  outputList->Add(mInvBetaDiffPionVsP);
  outputList->Add(mInvBetaDiffKaonVsP);
  outputList->Add(mInvBetaDiffProtonVsP);

  return outputList;
}
