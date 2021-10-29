/**
 * \class MpdFemtoDummyTrackCut
 * \brief The basic cut for tracks.
 *
 * Cuts on particle identification, transverse momentum, rapidity, distance
 * of closest approach to primary vertex and charge.
 *
 * ## Default cut values
 *
 * |  value   |  min | max |
 * |----------|------|-----|
 * | pt (GeV) |    0 | 100 |
 * | σ-pion   | -100 | 100 |
 * | σ-kaon   | -100 | 100 |
 * | σ-proton | -100 | 100 |
 * | \# Hits  |   10 | 180 |
 * | rapidity |   -2 |   2 |
 * | DCA (cm) |   -1 |  20 |
 *
 * Charge defaults to 1 (positive)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoDummyTrackCut_h
#define MpdFemtoDummyTrackCut_h

#include "MpdFemtoBaseTrackCut.h"

//_________________
class MpdFemtoDummyTrackCut : public MpdFemtoBaseTrackCut {
 public:
  /// Default constructor
  MpdFemtoDummyTrackCut();
    
  /// Test the particle and return true if it meets all criteria. false otherwise.
  virtual bool pass(const MpdFemtoTrack* tr);

  virtual MpdFemtoString report();
  virtual TList *listSettings();

  void setNSigmaPion(const float& lo, const float& hi) {
    fNSigmaPion[0] = lo;
    fNSigmaPion[1] = hi;
  }

  void setNSigmaKaon(const float& lo, const float& hi) {
    fNSigmaKaon[0] = lo;
    fNSigmaKaon[1] = hi;
  }

  void setNSigmaProton(const float& lo, const float& hi) {
    fNSigmaProton[0] = lo;
    fNSigmaProton[1] = hi;
  }

  void setNHits(const int& lo, const int& hi) {
    fNHits[0] = lo;
    fNHits[1] = hi;
  }

  void setPt(const float& lo, const float& hi) {
    fPt[0] = lo;
    fPt[1] = hi;
  }

  void setRapidity(const float& lo, const float& hi) {
    fRapidity[0] = lo;
    fRapidity[1] = hi;
  }

  void setDCA(const float& lo, const float& hi) {
    fDCA[0] = lo;
    fDCA[1] = hi;
  }

  void setCharge(const int& ch) {
    fCharge = ch;
  }

 protected:

  /// Charge of the track - if 0 the charge is not checked
  int fCharge;
  /// Bounds for nsigma dEdx from pion band
  float fNSigmaPion[2];
  /// Bounds for nsigma dEdx from kaon band
  float fNSigmaKaon[2];
  /// Bounds for nsigma dEdx from proton band
  float fNSigmaProton[2];
  /// Bounds for number of hits
  int fNHits[2];
  /// Bounds for transverse momentum
  float fPt[2];
  /// Bounds for rapidity
  float fRapidity[2];
  /// Bounds for DCA to primary vertex
  float fDCA[2];

  /// Passed tracks counter
  long fNTracksPassed;
  /// Falied tracks counter
  long fNTracksFailed; 

  ClassDef(MpdFemtoDummyTrackCut, 2)
};

#endif
