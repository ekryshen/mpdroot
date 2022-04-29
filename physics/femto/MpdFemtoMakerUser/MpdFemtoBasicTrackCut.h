/**
 * \class MpdFemtoBasicTrackCut
 * \brief The basic track cut
 *
 * Cuts on particle identification, transverse momentum, rapidity, distance
 * of closest approach to primary vertex and charge. Charge defaults to 1 (positive)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * email: nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBasicTrackCut_h
#define MpdFemtoBasicTrackCut_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseParticleCut.h"
#include "MpdFemtoBaseTrackCut.h"

// ROOT headers
#include "TLorentzVector.h"

// C++ headers
#include <limits>

//_________________
class MpdFemtoBasicTrackCut : public MpdFemtoBaseTrackCut {
public:
   /// Default constructor
   MpdFemtoBasicTrackCut();
   /// Copy constructor
   MpdFemtoBasicTrackCut(const MpdFemtoBasicTrackCut &copy);
   /// Assignment operator
   MpdFemtoBasicTrackCut &operator=(const MpdFemtoBasicTrackCut &copy);
   /// Destructor
   virtual ~MpdFemtoBasicTrackCut();

   /// Test the particle and return true if it meets all criteria. false otherwise.
   virtual bool pass(const MpdFemtoTrack *tr);

   /// Prepare report
   virtual MpdFemtoString report();
   /// Prepare list of settings
   virtual TList *listSettings();

   /// Enumeration of PID
   enum HbtPID { Electron = 1, Pion, Kaon, Proton };

   /*! Select track type:
   \param[in] type <table>
   <tr><td>false (0)</td><td>global</td></tr>
   <tr><td>true (1)</td><td>primary</td></tr>
   </table>
   */
   void selectPrimary(const bool &type) { setType(type); }
   /*! Set track type to select:
   \param[in] type <table>
   <tr><td>false (0)</td><td>global</td></tr>
   <tr><td>true (1)</td><td>primary</td></tr>
   </table>
   */
   void setType(const bool &type) { mType = type; }
   /// Set charge of the track to select
   void setCharge(const short &charge);
   /// Set min and max range of nhits
   void setNHits(const short &lo, const short &hi);
   /// Set min and max range of nHitsFit/nHitsPossible
   void setNHitsFitOverNHitsPoss(const float &lo) { mNHitsRat = lo; }
   /// Set min and max range of nHitsFit/nHitsPossible
   void setAntiSplit(const float &lo) { setNHitsFitOverNHitsPoss(lo); }
   /// Set min and max ranges of transverse momentum
   void setPt(const float &lo, const float &hi)
   {
      mPt[0] = lo;
      mPt[1] = hi;
   }
   /// Set min and max ranges of total momentum
   void setP(const float &lo, const float &hi)
   {
      mP[0] = lo;
      mP[1] = hi;
   }
   /// Set min and max ranges of transverse momentum
   void setTranverseMomentum(const float &lo, const float &hi) { setPt(lo, hi); }
   /// Set min and max ranges of rapidity
   void setRapidity(const float &lo, const float &hi)
   {
      mRapidity[0] = lo;
      mRapidity[1] = hi;
   }
   /// Set min and max ranges of pseudorapidity
   void setEta(const float &lo, const float &hi)
   {
      mEta[0] = lo;
      mEta[1] = hi;
   }
   /// Set min and max ranges of pseudorapidity
   void setPseudoRapidity(const float &lo, const float &hi) { setEta(lo, hi); }
   /// Set min and max range of DCA
   void setDCA(const float &lo, const float &hi)
   {
      mDCA[0] = lo;
      mDCA[1] = hi;
   }

   /*! Set type of detector selection:
   \param[in] selection <table>
   <tr><td>0</td><td>TPC</td></tr>
   <tr><td>1</td><td>TOF</td></tr>
   <tr><td>2</td><td>TPC+TOF</td></tr>
   <tr><td>3</td><td>if(TOF){TPC+TOF} else{TPC}</td></tr>
   <tr><td>4</td><td>if(TOF&&p>pthresh){TPC+TOF} else if (p<pthresh){TPC}</td></tr>
   </table>
   */
   void setDetectorSelection(const int &selection)
   {
      if (selection < 0 || selection > 4) {
         std::cout << "[WARNING] MpdFemtoBasicTrackCut: wrong detector selection: " << selection << " reset to TPC (0)"
                   << std::endl;
         mDetSelection = 0;
      } else {
         mDetSelection = (unsigned char)selection;
      }
   }

   /// Set min and max range of electron nSigma to select
   void setNSigmaElectron(const float &lo, const float &hi)
   {
      mNSigmaElectron[0] = lo;
      mNSigmaElectron[1] = hi;
   }
   /// Set min and max range of pion nSigma to select
   void setNSigmaPion(const float &lo, const float &hi)
   {
      mNSigmaPion[0] = lo;
      mNSigmaPion[1] = hi;
   }
   /// Set min and max range of kaon nSigma to select
   void setNSigmaKaon(const float &lo, const float &hi)
   {
      mNSigmaKaon[0] = lo;
      mNSigmaKaon[1] = hi;
   }
   /// Set min and max range of proton nSigma to select
   void setNSigmaProton(const float &lo, const float &hi)
   {
      mNSigmaProton[0] = lo;
      mNSigmaProton[1] = hi;
   }
   /// Set low and hight values for the exclusion cut
   void setNSigmaOther(const float &lo, const float &hi)
   {
      mNSigmaOther[0] = lo;
      mNSigmaOther[1] = hi;
   }
   /// Set min and max momentum of the track for TPC identification
   void setTpcP(const float &lo, const float &hi)
   {
      mTpcMom[0] = lo;
      mTpcMom[1] = hi;
   }
   /// Set min and max momentum of the track for TPC identification
   void setTpcMomentum(const float &lo, const float &hi) { setTpcP(lo, hi); }
   /// Set min and max momentum of the track for TPC identification
   void setTpcMom(const float &lo, const float &hi) { setTpcP(lo, hi); }

   /// Set min and max mass square estimated using TOF

   void setMassSqr(const float &lo, const float &hi)
   {
      mTofMassSqr[0] = lo;
      mTofMassSqr[1] = hi;
   }
   /// Set min and max momentum for TOF or TPC+TOF identification
   void setTofP(const float &lo, const float &hi)
   {
      mTofMom[0] = lo;
      mTofMom[1] = hi;
   }
   /// Set min and max momentum for TOF or TPC+TOF identification
   void setTofMomentum(const float &lo, const float &hi) { setTofP(lo, hi); }
   /// Set min and max momentum for TOF or TPC+TOF identification
   void setTofMom(const float &lo, const float &hi) { setTofP(lo, hi); }

   /// Set min and max values of electron nSigma for TPC+TOF identification
   void setTnTNSigmaElectron(const float &lo, const float &hi)
   {
      mTnTNSigmaElectron[0] = lo;
      mTnTNSigmaElectron[1] = hi;
   }
   /// Set min and max values of pion nSigma for TPC+TOF identification
   void setTnTNSigmaPion(const float &lo, const float &hi)
   {
      mTnTNSigmaPion[0] = lo;
      mTnTNSigmaPion[1] = hi;
   }
   /// Set min and max values of kaon nSigma for TPC+TOF identification
   void setTnTNSigmaKaon(const float &lo, const float &hi)
   {
      mTnTNSigmaKaon[0] = lo;
      mTnTNSigmaKaon[1] = hi;
   }
   /// Set min and max values of proton nSigma for TPC+TOF identification
   void setTnTNSigmaProton(const float &lo, const float &hi)
   {
      mTnTNSigmaProton[0] = lo;
      mTnTNSigmaProton[1] = hi;
   }

   /// Set PID to select: 1-electron, 2-pion, 3-kaon, 4-proton
   void setHbtPid(const HbtPID &pid) { mPidSelection = pid; }

   /// Is theoretical work
   void setIsTheory(const bool &isTheory) { mIsTheory = isTheory; }
   /// Set PDD ID to select (will be checked if IsTheory set to true)
   void setPdgId(const int &pdgId) { mPdgId = pdgId; }
   /// Set min and max 1/beta-1/beta_expected values
   void setInvBetaDiff(const float &lo, const float &hi)
   {
      mInvBetaDiff[0] = lo;
      mInvBetaDiff[1] = hi;
   }
   /// Set threshold for detector selection =4
   void setPthresh(const float &mom) { mPthresh = mom; }

   /*! Set verbose mode:
   \param[in] isVerbose <table>
   <tr><td>true</td><td>print information for each track</td></tr>
   <tr><td>false</td><td>do not print (default)</td></tr>
   </table>
   */
   void setVerboseMode(const bool &isVerbose) { mVerbose = isVerbose; }

protected:
   /// Track type: 0 - global, 1 - primary
   bool mType;
   /// Charge of the track - if 0 the charge is not checked
   char mCharge;
   /// Bounds for number of hits
   unsigned char mNHits[2];
   /// nHitsFit/nHitsPossible
   float mNHitsRat;
   /// Bounds for transverse momentum
   float mPt[2];
   /// Bounds for full momentum
   float mP[2];
   /// Bounds for rapidity
   float mRapidity[2];
   /// Bounds for pseudorapidity
   float mEta[2];
   /// Bounds for DCA to primary vertex
   float mDCA[2];

   /*! Type of detector selection:
   \param[in] selection <table>
   <tr><td>0</td><td>TPC</td></tr>
   <tr><td>1</td><td>TOF</td></tr>
   <tr><td>2</td><td>TPC+TOF</td></tr>
   <tr><td>3</td><td>if(TOF){TPC+TOF} else{TPC}</td></tr>
   <tr><td>4</td><td>if(TOF&&p>pthresh){TPC+TOF} else if (p<pthresh){TPC}</td></tr>
   </table>
   */
   unsigned char mDetSelection;

   // TPC identification

   /// Bounds for nsigma dEdx from electron band
   float mNSigmaElectron[2];
   /// Bounds for nsigma dEdx from pion band
   float mNSigmaPion[2];
   /// Bounds for nsigma dEdx from kaon band
   float mNSigmaKaon[2];
   /// Bounds for nsigma dEdx from proton band
   float mNSigmaProton[2];
   /// Bounds for exclusion cuts
   float mNSigmaOther[2];
   /// Momentum of the track for pure TPC identification
   float mTpcMom[2];

   // TOF identification

   /// TOF mass square
   float mTofMassSqr[2];
   /// Momentum of the track for the pure TOF identification
   float mTofMom[2];

   // TPC && TOF identification

   /// Bounds for nsigma dEdx from pion band
   float mTnTNSigmaElectron[2];
   /// Bounds for nsigma dEdx from pion band
   float mTnTNSigmaPion[2];
   /// Bounds for nsigma dEdx from kaon band
   float mTnTNSigmaKaon[2];
   /// Bounds for nsigma dEdx from proton band
   float mTnTNSigmaProton[2];

   /*! PID selection:
   \param[in] selection <table>
   <tr><td>1</td><td>electron</td></tr>
   <tr><td>2</td><td>pion</td></tr>
   <tr><td>3</td><td>kaon</td></tr>
   <tr><td>4</td><td>proton</td></tr>
   </table>
   */
   HbtPID mPidSelection;

   /// Passed tracks counter
   unsigned int mNTracksPassed;
   /// Falied tracks counter
   unsigned int mNTracksFailed;

   /// If work with theory
   bool mIsTheory;
   /// PDG pid
   int mPdgId;
   /// Bounds for 1/beta-1/beta(expected)
   float mInvBetaDiff[2];
   /// Momentum threshold for (detector) selection type == 4
   float mPthresh;
   /// Print information for each pair (default is false)
   bool mVerbose;

   //
   // Functions for PID
   //

   /// Is good clean TPC identification
   bool isGoodCleanTpcPid(const MpdFemtoTrack *t, TLorentzVector mom);
   /// Is good clean electron identification in TPC
   bool isTpcCleanElectron(const MpdFemtoTrack *t, TLorentzVector mom);
   /// Is good clean pion identification in TPC
   bool isTpcCleanPion(const MpdFemtoTrack *t, TLorentzVector mom);
   /// Is good clean kaon identification in TPC
   bool isTpcCleanKaon(const MpdFemtoTrack *t, TLorentzVector mom);
   /// Is good clean proton identification in TPC
   bool isTpcCleanProton(const MpdFemtoTrack *t, TLorentzVector mom);

   /// Is good particle identification by TOF
   bool isGoodTofPid(const MpdFemtoTrack *t);

   /// Is TPC electron by nSigma(electron) only (uses TnT sigma cut values)
   bool isTpcElectron(const MpdFemtoTrack *t);
   /// Is TPC electron by nSigma(pion) only (uses TnT sigma cut values)
   bool isTpcPion(const MpdFemtoTrack *t);
   /// Is TPC electron by nSigma(kaon) only (uses TnT sigma cut values)
   bool isTpcKaon(const MpdFemtoTrack *t);
   /// Is TPC electron by nSigma(proton) only (uses TnT sigma cut values)
   bool isTpcProton(const MpdFemtoTrack *t);

   /// Check that passes TOF and TPC identification (TnT = TPC && TOF)
   bool isGoodTnT(const MpdFemtoTrack *t);
   /// Check that passes TPC+TOF identification in case of TOF signal,
   /// and clean TPC in the other case
   bool isGoodToT(const MpdFemtoTrack *t, TLorentzVector vec);
   /// Check that passes TPC+TOF identification in case of TOF signal and
   /// momentum greater or equal to threshold, and clean TPC in the other case
   bool isGoodToTThresh(const MpdFemtoTrack *t, TLorentzVector vec);

   /// Check kinematics
   bool isGoodKine(const MpdFemtoTrack *t, TLorentzVector vec);

   ClassDef(MpdFemtoBasicTrackCut, 2);
};

#endif // MpdFemtoBasicTrackCut_h
