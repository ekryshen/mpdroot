#ifndef MPDTRACKPIDMAKER_H
#define MPDTRACKPIDMAKER_H
#include "MpdAnalysisTask.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcKalmanFilter.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"

class MpdTpcKalmanTrack;

class MpdTrackPidMaker : public MpdAnalysisTask {

public:
   MpdTrackPidMaker() = default;
   MpdTrackPidMaker(const char *name, const char *outputName = "taskName");
   ~MpdTrackPidMaker() = default;

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

protected:
   bool selectEvent(MpdAnalysisEvent &event);

   float dEdx_sigma_El(float dEdx, float mom) const;
   float dEdx_sigma_Pi(float dEdx, float mom) const;
   float dEdx_sigma_K(float dEdx, float mom) const;
   float dEdx_sigma_P(float dEdx, float mom) const;
   float dEdx_sigma_De(float dEdx, float mom) const;
   float dEdx_sigma_Tr(float dEdx, float mom) const;
   float dEdx_sigma_He3(float dEdx, float mom) const;
   float dEdx_sigma_He4(float dEdx, float mom) const;

   float TofdPhiMatch(int charge, float pt, float dphi) const;
   float TofdZedMatch(float pt, float dzed) const;
   float Beta_sigma_El(float beta, float mom) const;
   float Beta_sigma_Pi(float beta, float mom) const;
   float Beta_sigma_K(float beta, float mom) const;
   float Beta_sigma_P(float beta, float mom) const;
   float Beta_sigma_De(float beta, float mom) const;
   float Beta_sigma_Tr(float beta, float mom) const;
   float Beta_sigma_He3(float beta, float mom) const;
   float Beta_sigma_He4(float beta, float mom) const;

   float EmcdPhiMatch(int charge, float pt, float dphi) const;
   float EmcdZedMatch(float pt, float dzed) const;

private:
   float    cen;
   TVector3 mPrimaryVertex;
   TRandom  RND;

   // DCAs
   float                  eta_max    = 1.5;
   float                  eta_min    = -1.5;
   float                  cent_max   = 100.;
   float                  cent_min   = 0.0;
   static constexpr short neta_bins  = 30;
   static constexpr short ncent_bins = 10;
   TF1                   *f_dca_xy[neta_bins][ncent_bins];
   TF1                   *f_dca_z[neta_bins][ncent_bins];

   MpdKalmanFilter *pKF = nullptr;

   TH1F *mhEvents     = nullptr;
   TH1F *mhVertex     = nullptr;
   TH1F *mhCentrality = nullptr;
   TH2F *mhDCAXs      = nullptr;
   TH2F *mhDCAYs      = nullptr;
   TH2F *mhDCAZs      = nullptr;
   TH2F *mhTPCEl      = nullptr;
   TH2F *mhTPCPi      = nullptr;
   TH2F *mhTPCK       = nullptr;
   TH2F *mhTPCP       = nullptr;
   TH2F *mhTPCDe      = nullptr;
   TH2F *mhTPCTr      = nullptr;
   TH2F *mhTPCHe3     = nullptr;
   TH2F *mhTPCHe4     = nullptr;
   TH2F *mhTOFsdPhi   = nullptr;
   TH2F *mhTOFsdZed   = nullptr;
   TH2F *mhTOFEl      = nullptr;
   TH2F *mhTOFPi      = nullptr;
   TH2F *mhTOFK       = nullptr;
   TH2F *mhTOFP       = nullptr;
   TH2F *mhTOFDe      = nullptr;
   TH2F *mhTOFTr      = nullptr;
   TH2F *mhTOFHe3     = nullptr;
   TH2F *mhTOFHe4     = nullptr;
   TH2F *mhECALsdPhi  = nullptr;
   TH2F *mhECALsdZed  = nullptr;
   TH2F *mhEP         = nullptr;
   TH2F *mhEM2        = nullptr;

   TFile *dcaFile;

   ClassDef(MpdTrackPidMaker, 2);
};
#endif
