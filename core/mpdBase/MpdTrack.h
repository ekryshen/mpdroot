// Author: Emelyanov D.
// Update: Oleg Rogachevsky 2009-09-17 17:42:25+0400
// Copyright: 2009 (C) Oleg Rogachevsky

#ifndef ROOT_MpdTrack
#define ROOT_MpdTrack
#include "TObject.h"
#include "MpdHelix.h"

// Number of lines in TPC dE/dx PID and in TOF PID
static constexpr int kNPidLines = 8;
enum { kEl, kPi, kK, kP, kDeutron, kTriton, kHe3, kHe4 };

class MpdTrack : public TObject {

public:
   MpdTrack()          = default;
   virtual ~MpdTrack() = default;

   // Basic properties
   int       GetID() const { return fID; }
   int       GetNofHits() const { return fNofHits; }
   int       GetNofHitsPossTpc() const { return fNofHitsPossTpc; }
   int       GetNofHitsFitTpc() const { return fNofHitsFitTpc; }
   float     GetChi2() const { return fChi2; }
   int16_t   GetHelixQ() const { return fHelixQ; }
   bool      GetEdgeCut() { return fEdgeCut; }
   ULong64_t GetLayerHitMap() const { return fHitMap; };
   ULong64_t GetSharedHitMap() const { return fSharedHitMap; };
   int       GetNSharedTpcHits() const;
   MpdHelix  GetHelix() const;

   void SetID(int n) { fID = n; }
   void SetNofHits(int n) { fNofHits = n; }
   void SetNofHitsPossTpc(int n) { fNofHitsPossTpc = n; }
   void SetNofHitsFitTpc(int n) { fNofHitsFitTpc = n; }
   void SetChi2(float n) { fChi2 = n; }
   void SetHelixQ(int16_t n) { fHelixQ = n; }
   void SetEdgeCut(bool n) { fEdgeCut = n; }
   void SetLayerHitMap(ULong64_t map) { fHitMap = map; };
   void SetSharedHitMap(ULong64_t map) { fSharedHitMap = map; }

   // Kinematics
   float GetCharge() const { return (fPt > 0) ? -1. : 1.; }
   float GetPt() const { return TMath::Abs(fPt); } // Note that charge removed here!
   float GetTheta() const { return fTheta; }
   float GetPhi() const { return fPhi; }
   float GetPtError() const { return fPtError; }
   float GetThetaError() const { return fThetaError; }
   float GetPhiError() const { return fPhiError; }
   float GetPx() const { return TMath::Abs(fPt) * TMath::Cos(fPhi); }
   float GetPy() const { return TMath::Abs(fPt) * TMath::Sin(fPhi); }
   float GetPz() const
   {
      double tt = TMath::Tan(fTheta);
      return tt ? TMath::Abs(fPt) / tt : 9999.;
   }
   float GetEta() const { return -TMath::Log(TMath::Tan(0.5 * fTheta)); }
   void  SetPt(float n) { fPt = n; }
   void  SetTheta(float n) { fTheta = n; }
   void  SetPhi(float n) { fPhi = n; }
   void  SetPtError(float n) { fPtError = n; }
   void  SetThetaError(float n) { fThetaError = n; }
   void  SetPhiError(float n) { fPhiError = n; }

   // DCA: Distance of Closest Approach perpendicular and along Z axis
   float GetDCAX() const { return fDCAX; }
   float GetDCAY() const { return fDCAY; }
   float GetDCAZ() const { return fDCAZ; }
   float GetDCAGlobalX() const { return fDCAGlobalX; }
   float GetDCAGlobalY() const { return fDCAGlobalY; }
   float GetDCAGlobalZ() const { return fDCAGlobalZ; }
   float GetFirstPointX() const { return fFirstPointX; }
   float GetFirstPointY() const { return fFirstPointY; }
   float GetFirstPointZ() const { return fFirstPointZ; }
   float GetLastPointX() const { return fLastPointX; }
   float GetLastPointY() const { return fLastPointY; }
   float GetLastPointZ() const { return fLastPointZ; }
   float GetNSigmaDCAx() const { return fDCAxNSigma; }
   float GetNSigmaDCAy() const { return fDCAyNSigma; }
   float GetNSigmaDCAz() const { return fDCAzNSigma; }

   void SetDCAX(float n) { fDCAX = n; }
   void SetDCAY(float n) { fDCAY = n; }
   void SetDCAZ(float n) { fDCAZ = n; }
   void SetDCAGlobalX(float n) { fDCAGlobalX = n; }
   void SetDCAGlobalY(float n) { fDCAGlobalY = n; }
   void SetDCAGlobalZ(float n) { fDCAGlobalZ = n; }
   void SetFirstPointX(float n) { fFirstPointX = n; }
   void SetFirstPointY(float n) { fFirstPointY = n; }
   void SetFirstPointZ(float n) { fFirstPointZ = n; }
   void SetLastPointX(float n) { fLastPointX = n; }
   void SetLastPointY(float n) { fLastPointY = n; }
   void SetLastPointZ(float n) { fLastPointZ = n; }
   void SetNSigmaDCAx(float n) { fDCAxNSigma = n; }
   void SetNSigmaDCAy(float n) { fDCAyNSigma = n; }
   void SetNSigmaDCAz(float n) { fDCAzNSigma = n; }

   // TPC dEdx
   float GetdEdXTPC() const { return fdEdXTPC; }
   float GetTPCNSigma(int icut = kEl) const { return fTPCNSigma[icut]; }
   float GetNSigmaElectron() const { return fTPCNSigma[kEl]; }
   float GetNSigmaPion() const { return fTPCNSigma[kPi]; }
   float GetNSigmaKaon() const { return fTPCNSigma[kK]; }
   float GetNSigmaProton() const { return fTPCNSigma[kP]; }

   void SetdEdXTPC(float n) { fdEdXTPC = n; }
   void SetTPCNSigma(float *cuts)
   {
      for (int i = kNPidLines; i--;) fTPCNSigma[i] = cuts[i];
   }
   void SetNSigmaElectron(float n) { fTPCNSigma[kEl] = n; }
   void SetNSigmaPion(float n) { fTPCNSigma[kPi] = n; }
   void SetNSigmaKaon(float n) { fTPCNSigma[kK] = n; }
   void SetNSigmaProton(float n) { fTPCNSigma[kP] = n; }

   // Bayesian TPC PID
   float GetTPCPidProbElectron() const { return fPidTPCProb[kEl]; }
   float GetTPCPidProbPion() const { return fPidTPCProb[kPi]; }
   float GetTPCPidProbKaon() const { return fPidTPCProb[kK]; }
   float GetTPCPidProbProton() const { return fPidTPCProb[kP]; }
   void  SetTPCpidProb(float n1, float n2, float n3, float n4, int flag)
   {
      fPidTPCProb[kEl] = n1;
      fPidTPCProb[kPi] = n2;
      fPidTPCProb[kK]  = n3;
      fPidTPCProb[kP]  = n4;
      SetTofFlag(flag);
   }

   // TOF PID
   float GetTOFNSigma(int icut = kEl) const { return fTOFNSigma[icut]; }
   float GetTofBeta() const { return fTofBeta; }
   float GetTofMass2() const { return fTofMass2; }
   int   GetTofHitIndex() const { return fTofHitIndex; };
   int   GetTofFlag() const { return fTofFlag; }
   // Track-TOF matching in sigma
   float GetTofDphiSigma() const { return fTofDphiSigma; }
   float GetTofDzSigma() const { return fTofDzSigma; }

   void SetTOFNSigma(float *cuts)
   {
      for (int i = kNPidLines; i--;) fTOFNSigma[i] = cuts[i];
   }
   void SetTofBeta(float n) { fTofBeta = n; }
   void SetTofMass2(float n) { fTofMass2 = n; }
   void SetTofHitIndex(int n) { fTofHitIndex = n; }
   void SetTofFlag(int n) { fTofFlag = fTofFlag | n; }
   void SetTofMatching(float dphi, float dz)
   {
      fTofDphiSigma = dphi;
      fTofDzSigma   = dz;
   }

   // TOF bayesian PID
   float GetTOFPidProbElectron() const { return fPidTOFProb[kEl]; }
   float GetTOFPidProbPion() const { return fPidTOFProb[kPi]; }
   float GetTOFPidProbKaon() const { return fPidTOFProb[kK]; }
   float GetTOFPidProbProton() const { return fPidTOFProb[kP]; }
   void  SetTOFpidProb(float n1, float n2, float n3, float n4, int flag)
   {
      fPidTOFProb[kEl] = n1;
      fPidTOFProb[kPi] = n2;
      fPidTOFProb[kK]  = n3;
      fPidTOFProb[kP]  = n4;
      SetTofFlag(flag);
   }

   // ECAL
   float GetEp() const { return fEp; }
   float GetEt() const { return fEt; }
   float GetEl() const { return fEl; }
   int   GetEi() const { return fEi; }
   // Track-cluster matching in sigmal
   float GetECALDphiSigma() const { return fEcalDphiSigma; }
   float GetECALDzSigma() const { return fEcalDzSigma; }

   void SetEp(float ep) { fEp = ep; }
   void SetEt(float ep) { fEt = ep; }
   void SetEl(float ep) { fEl = ep; }
   void SetEi(int ep) { fEi = ep; }

   void SetEcalMatching(float dphi, float dz)
   {
      fEcalDphiSigma = dphi;
      fEcalDzSigma   = dz;
   }

   // Global Bayesian probability
   float GetPidProbElectron() const { return fPidProb[kEl]; }
   float GetPidProbPion() const { return fPidProb[kPi]; }
   float GetPidProbKaon() const { return fPidProb[kK]; }
   float GetPidProbProton() const { return fPidProb[kP]; }

   void SetPidProbElectron(float n) { fPidProb[kEl] = n; }
   void SetPidProbPion(float n) { fPidProb[kPi] = n; }
   void SetPidProbKaon(float n) { fPidProb[kK] = n; }
   void SetPidProbProton(float n) { fPidProb[kP] = n; }
   void SetCombPidProb(float n1, float n2, float n3, float n4)
   {
      fPidProb[kEl] = n1;
      fPidProb[kPi] = n2;
      fPidProb[kK]  = n3;
      fPidProb[kP]  = n4;
   }

private:
   int   fID             = -1; // Track number from reconstruction
   int   fNofHits        = 0;  // Number of track hits
   int   fNofHitsPossTpc = 0;  //
   int   fNofHitsFitTpc  = 0;  //
   float fChi2           = 0.; // Track Chi2

   // Kinematics
   float fPt         = 0.; // Signed Transverse momentum
   float fTheta      = 0.; // Theta angle (from beam line)
   float fPhi        = 0.; // Phi angle
   float fPtError    = 0.; // Pt error
   float fThetaError = 0.; // Theta error
   float fPhiError   = 0.; // Phi error

   // DCA
   float fDCAX        = 0.; //
   float fDCAY        = 0.; //
   float fDCAZ        = 0.; //
   float fDCAGlobalX  = 0.; //
   float fDCAGlobalY  = 0.; //
   float fDCAGlobalZ  = 0.; //
   float fFirstPointX = 0.; // point closest to beam line
   float fFirstPointY = 0.; // point closest to beam line
   float fFirstPointZ = 0.; // point closest to beam line
   float fLastPointX  = 0.; //
   float fLastPointY  = 0.; //
   float fLastPointZ  = 0.; //
   float fDCAxNSigma  = 0.; //
   float fDCAyNSigma  = 0.; //
   float fDCAzNSigma  = 0.; //

   // TPC sigmas
   float fdEdXTPC               = 0.;  // TPC dEdx
   float fTPCNSigma[kNPidLines] = {0}; // distance to corresponding line in sigmas
   // TPC bayesian probability
   float fPidTPCProb[kNPidLines] = {0}; // TPC probabilities

   // TOF data
   float fTOFNSigma[kNPidLines] = {0};  // distance to corresponding line in sigmas
   float fTofBeta               = 0.;   // beta
   float fTofMass2              = 0.;   // identified m^2, can be >0, <0 and Nan!
   float fTofDphiSigma          = 999.; // track-TOF matching in sigmas
   float fTofDzSigma            = 999.; // track-TOF matching in sigmas
   int   fTofHitIndex           = -1;   //
   int   fTofFlag               = 0;    // it is set if tof identification exists

   // TOF bayesian probability
   float fPidTOFProb[kNPidLines] = {0}; // [tof] probability

   // ECAL
   float fEp            = -999.; // ECAL E
   float fEt            = -999.; // ECAL tof
   float fEl            = -999.; // ECAL track length
   int   fEi            = -999.; // ECAL track index
   float fEcalDphiSigma = -999.;
   float fEcalDzSigma   = -999.;

   // Global bayesian probability
   float fPidProb[kNPidLines] = {0};

   int16_t fHelixQ;          //
   bool    fEdgeCut = false; // kTRUE if number of hits closer to boundaries than 1.5 cm divided by nHits is larger than
                             // 50% (else: kFALSE)
   ULong64_t fHitMap       = 0; // n-th bit is 1 if there is hit in layer, 0 otherwise
   ULong64_t fSharedHitMap = 0; // n-th bit if bit if hit in given layer shared with any other track

   ClassDef(MpdTrack, 5);
};

#endif
