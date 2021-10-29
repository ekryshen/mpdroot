/**
 * \class MpdFemtoFxtTrackCutMonitor
 * \brief Track cut monitor for basic analysis
 *
 * The class provides histograms for monitoring track cuts
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoFxtTrackCutMonitor_h
#define MpdFemtoFxtTrackCutMonitor_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCutMonitor.h"

// ROOT headers
#include "TH1F.h"
#include "TH2F.h"

// Forward declaration
class MpdFemtoTrack;

//_________________
class MpdFemtoFxtTrackCutMonitor : public MpdFemtoBaseCutMonitor {
 public:
  /// Constructor with mass
  MpdFemtoFxtTrackCutMonitor(const char*, const double);
  /// Copy constructor
  MpdFemtoFxtTrackCutMonitor(const MpdFemtoFxtTrackCutMonitor& cutMoni);
  /// Assignment operator
  MpdFemtoFxtTrackCutMonitor operator=(const MpdFemtoFxtTrackCutMonitor& c);
  /// Destructor
  virtual ~MpdFemtoFxtTrackCutMonitor();

  /// Fill histograms with track information
  virtual void fill(const MpdFemtoTrack* track);

  /// Write all histograms
  void writeOutHistos();

  /// Ouput list with histograms
  virtual TList* getOutputList();

  // These dummy Fill() functions were introduced to remove a compiler
  // warning related to overloaded base-class Fill() functions being
  // hidden by a single version of Fill() in this derived class
  void fill(const MpdFemtoEvent* /* d */) {
    ;
  }

  void fill(const MpdFemtoV0* /* d */) {
    ;
  }

  void fill(const MpdFemtoXi* /* xi */) {
    ;
  }

  void fill(const MpdFemtoKink* /* d */) {
    ;
  }

  void fill(const MpdFemtoPair* /* d */) {
    ;
  }

  void fill(const MpdFemtoParticleCollection* /* d */) {
    ;
  }

  void fill(const MpdFemtoEvent* /* d1 */, const MpdFemtoParticleCollection* /* d2 */) {
    ;
  }

  void fill(const MpdFemtoParticleCollection* /* d1 */, const MpdFemtoParticleCollection* /* d2 */) {
    ;
  }

  /// DCA of the track
  TH1F* dcaGlobal() {
    return (mDCAGlobal) ? mDCAGlobal : nullptr;
  }
  /// DCA of the track
  TH1F* DCAGlobal() {
    return dcaGlobal();
  }
  /// Number of hits
  TH1F* nHits() {
    return (mNhits) ? mNhits : nullptr;
  }
  /// Number of hits
  TH1F* NHits() {
    return nHits();
  }
  /// Track momentum
  TH1F* p() {
    return (mP) ? mP : nullptr;
  }
  /// Track momentum
  TH1F* P() {
    return p();
  }
  /// Track transverse momentum
  TH1F* pt() {
    return (mPt) ? mPt : nullptr;
  }
  /// Track transverse momentum
  TH1F* Pt() {
    return pt();
  }
  /// nSigma(pion) vs pT
  TH2F* nSigmaPionVsPt() {
    return (mPtVsNsigmaPion) ? mPtVsNsigmaPion : nullptr;
  }
  /// nSigma(pion) vs pT
  TH2F* NSigmaPionVsPt() {
    return nSigmaPionVsPt();
  }
  /// nSigma(pion) vs pT
  TH2F* PtVsNsigmaPion() {
    return nSigmaPionVsPt();
  }
  /// nSigma(kaon) vs pT
  TH2F* nSigmaKaonVsPt() {
    return (mPtVsNsigmaKaon) ? mPtVsNsigmaKaon : nullptr;
  }
  /// nSigma(kaon) vs pT
  TH2F* NSigmaKaonVsPt() {
    return nSigmaKaonVsPt();
  }
  /// nSigma(kaon) vs pT
  TH2F* PtVsNsigmaKaon() {
    return nSigmaKaonVsPt();
  }
  /// nSigma(proton) vs pT
  TH2F* nSigmaProtonVsPt() {
    return (mPtVsNsigmaProton) ? mPtVsNsigmaProton : nullptr;
  }
  /// nSigma(proton) vs pT
  TH2F* NSigmaProtonVsPt() {
    return nSigmaProtonVsPt();
  }
  /// nSigma(proton) vs pT
  TH2F* PtVsNsigmaProton() {
    return nSigmaProtonVsPt();
  }
  /// dE/dx vs momentum
  TH2F* dEdxVsP() {
    return (mPvsDedx) ? mPvsDedx : nullptr;
  }
  /// dE/dx vs momentum
  TH2F* DedxVsP() {
    return dEdxVsP();
  }
  /// dE/dx vs momentum
  TH2F* PvsDedx() {
    return dEdxVsP();
  }
  /// Pseudorapidity
  TH1F* eta() {
    return (mPseudoRapidity) ? mPseudoRapidity : nullptr;
  }
  /// Pseudorapidity
  TH1F* Eta() {
    return eta();
  }
  /// Pseudorapidity
  TH1F* PseudoRapidity() {
    return eta();
  }
  /// Square of mass vs. momentum
  TH2F* massSqrVsP() {
    return (mPvsMassSqr) ? mPvsMassSqr : nullptr;
  }
  /// Square of mass vs. momentum
  TH2F* MassSqrVsP() {
    return massSqrVsP();
  }
  /// Square of mass vs. momentum
  TH2F* PvsMassSqr() {
    return massSqrVsP();
  }
  /// Inversed beta vs. momentum
  TH2F* inversedBetaVsP() {
    return (mPvsInvBeta) ? mPvsInvBeta : nullptr;
  }
  /// Inversed beta vs. momentum
  TH2F* InversedBetaVsP() {
    return inversedBetaVsP();
  }
  /// Inversed beta vs. momentum
  TH2F* PvsInvBeta() {
    return inversedBetaVsP();
  }
  /// Transverse momenum vs. pseudorapidity
  TH2F* ptVsEta() {
    return (mPtVsEta) ? mPtVsEta : nullptr;
  }
  /// Transverse momenum vs. pseudorapidity
  TH2F* PtVsEta() {
    return ptVsEta();
  }
  /// Inveresed beta minus inversed beta for pion hypothesis vs. momentum
  TH2F* inversedBetaDiffPionVsP() {
    return (mInvBetaDiffPionVsP) ? mInvBetaDiffPionVsP : nullptr;
  }
  /// Inveresed beta minus inversed beta for pion hypothesis vs. momentum
  TH2F* InversedBetaDiffPionVsP() {
    return inversedBetaDiffPionVsP();
  }
  /// Inveresed beta minus inversed beta for kaon hypothesis vs. momentum
  TH2F* inversedBetaDiffKaonVsP() {
    return (mInvBetaDiffKaonVsP) ? mInvBetaDiffKaonVsP : nullptr;
  }
  /// Inveresed beta minus inversed beta for kaon hypothesis vs. momentum
  TH2F* InversedBetaDiffKaonVsP() {
    return inversedBetaDiffKaonVsP();
  }
  /// Inveresed beta minus inversed beta for proton hypothesis vs. momentum
  TH2F* inversedBetaDiffProtonVsP() {
    return (mInvBetaDiffProtonVsP) ? mInvBetaDiffProtonVsP : nullptr;
  }
  /// Inveresed beta minus inversed beta for proton hypothesis vs. momentum
  TH2F* InversedBetaDiffProtonVsP() {
    return inversedBetaDiffProtonVsP();
  }

 private:

  /// DCA of the track to primary vertex
  TH1F* mDCAGlobal;
  /// Number of hits
  TH1F* mNhits;
  /// Momentum
  TH1F* mP;
  /// Transverse momentum
  TH1F* mPt;
  /// nSigma(pi) vs. transverse momentum
  TH2F* mPtVsNsigmaPion;
  /// nSigma(K) vs. transverse momentum
  TH2F* mPtVsNsigmaKaon;
  /// nSigma(p) vs. transverse momentum
  TH2F* mPtVsNsigmaProton;
  /// dE/dx vs. momentum
  TH2F* mPvsDedx;
  /// Rapidity
  TH1F* mRapidity;
  /// Pseudorapidity
  TH1F* mPseudoRapidity;
  /// Squared mass vs. momentum
  TH2F* mPvsMassSqr;
  /// Inversed beta vs. momentum
  TH2F* mPvsInvBeta;
  /// Pseudorapidity vs. transverse momentum
  TH2F* mPtVsEta;
  /// Inversed beta minus inversed beta for pion hypothesis vs. momentum
  TH2F* mInvBetaDiffPionVsP;
  /// Inversed beta minus inversed beta for kaon hypothesis vs. momentum
  TH2F* mInvBetaDiffKaonVsP;
  /// Inversed beta minus inversed beta for proton hypothesis vs. momentum
  TH2F* mInvBetaDiffProtonVsP;
  /// Particle mass
  double monMass;

  ClassDef(MpdFemtoFxtTrackCutMonitor, 3)
};

#endif // #define MpdFemtoFxtTrackCutMonitor_h
