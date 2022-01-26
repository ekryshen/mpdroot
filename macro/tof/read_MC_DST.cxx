
#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TCutG.h"

// Fair includes
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"
#include "FairTrackParP.h"
#include "FairMCEventHeader.h"
#include "FairMCTrack.h"

// MPD includes
#include "MpdEvent.h"

#include "MpdKalmanFilter.h"
#include "TpcLheHitsMaker.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"

#include "MpdTofPoint.h"
#include "MpdEtofPoint.h"
#include "MpdTofHit.h"
#include "MpdEtofHit.h"

#include "MpdTofMatching.h"
#include "MpdEtofMatching.h"

#include <iostream>
#include <fstream>
#include <map>
#include <set>
using namespace std;
#endif

//------------------------------------------------------------------------------------------------------------------------------

TH1D *Make1D(const char *, const char *, const char *, const char *, TList *, Int_t, Double_t, Double_t);
TH2D *Make2D(const char *, const char *, const char *, const char *, TList *, Int_t, Double_t, Double_t, Int_t,
             Double_t, Double_t);

// make clone with new name
TH1D *Make1D(const char *, const TH1D *, TList *);
TH2D *Make2D(const char *, const TH2D *, TList *);

// ----------- DST data ----------- tracks sorting by  bigger  fNofHits; // Number of track hits
struct ltDST {
   bool operator()(MpdTrack *s1, MpdTrack *s2) const { return s1->GetNofHits() > s2->GetNofHits(); }
};
typedef set<MpdTrack *, ltDST> DstTrackDef;
typedef DstTrackDef::iterator  DstTrackDefIter;

// ----------- MpdKalmanTrack data ----------- tracks sorting by  bigger  fNofHits; // Number of track hits
struct ltKF {
   bool operator()(MpdKalmanTrack *s1, MpdKalmanTrack *s2) const { return s1->GetNofHits() > s2->GetNofHits(); }
};
typedef set<MpdKalmanTrack *, ltKF> KfTrackDef;
typedef KfTrackDef::iterator        KfTrackDefIter;

// ----------- MpdTofPoint data ----------- MC points sorting by  smaller  fTime; // Time since event start [ns]
struct ltTofPoint {
   bool operator()(const MpdTofPoint *s1, MpdTofPoint *s2) const { return s1->GetTime() < s2->GetTime(); }
};
typedef set<MpdTofPoint *, ltTofPoint> TofPointDef;
typedef TofPointDef::iterator          TofPointDefIter;

// ----------- MpdETofPoint data ----------- MC points sorting by  smaller  fTime; // Time since event start [ns]
struct ltETofPoint {
   bool operator()(const MpdEtofPoint *s1, MpdEtofPoint *s2) const { return s1->GetTime() < s2->GetTime(); }
};
typedef set<MpdEtofPoint *, ltETofPoint> ETofPointDef;
typedef ETofPointDef::iterator           ETofPointDefIter;

typedef multimap<Int_t, MpdTofHit *> TofHitDef; // key - MC track index
typedef TofHitDef::iterator          TofHitDefIter;

typedef multimap<Int_t, MpdEtofHit *> ETofHitDef; // key - MC track index
typedef ETofHitDef::iterator          ETofHitDefIter;

typedef multimap<Int_t, MpdTofMatchingData *> TofMatchDef; // key - MC track index
typedef TofMatchDef::iterator                 TofMatchDefIter;

typedef multimap<Int_t, MpdEtofMatchingData *> ETofMatchDef; // key - MC track index
typedef ETofMatchDef::iterator                 ETofMatchDefIter;

// --------------------------------------------  1 track <-link-> a few objects(set)
struct mpdTestTrack {
   DstTrackDef  aDSTTrack;  // DST data	{MpdTrack}
   KfTrackDef   aKFTrack;   // reco data 	{MpdKalmanTrack}
   TofPointDef  aTofPoint;  // Tof data 	{MpdTofPoint}
   ETofPointDef aETofPoint; // ETof data 	{MpdEtofPoint}

   void clear(void)
   {
      aDSTTrack.clear();
      aKFTrack.clear();

      aTofPoint.clear();
      aETofPoint.clear();
   }
};

typedef map<Int_t, mpdTestTrack> EventDef; // key - MC track index
typedef EventDef::iterator       EventDefIter;
//--------------------------------------------------------------------------------------------------------------
void doTest(void)
{

   TList fList, fBarrel, fEndcap;
   TH2D *h_nPoints_P = Make2D("h_nPoints_P", "Number of MpdTofPoints per one mcTrack vs P", "N_{tofPoints}", "P, GeV/c",
                              &fList, 21, -0.5, 20.5, 1000, 0., 5.);
   TH2D *h_nHits_P = Make2D("h_nHits_P", "Number of MpdTofHits per one mcTrack vs P", "N_{tofHits}", "P, GeV/c", &fList,
                            21, -0.5, 20.5, 1000, 0., 5.);
   TH2D *h_nMatchs_P = Make2D("h_nMatchs_P", "Number of MpdTofMatchingDatas per one mcTrack vs P", "N_{tofMatches}",
                              "P, GeV/c", &fList, 21, -0.5, 20.5, 1000, 0., 5.);

   TH2D *h_nEPoints_P = Make2D("h_nEPoints_P", "Number of MpdETofPoints per one mcTrack vs P", "N_{etofPoints}",
                               "P, GeV/c", &fList, 21, -0.5, 20.5, 1000, 0., 5.);
   TH2D *h_nEHits_P   = Make2D("h_nEHits_P", "Number of MpdETofHits per one mcTrack vs P", "N_{etofHits}", "P, GeV/c",
                             &fList, 21, -0.5, 20.5, 1000, 0., 5.);
   TH2D *h_nEMatchs_P = Make2D("h_nEMatchs_P", "Number of MpdETofMatchingDatas per one mcTrack vs P", "N_{etofMatches}",
                               "P, GeV/c", &fList, 21, -0.5, 20.5, 1000, 0., 5.);

   TH2D *h_nKFtrack_P  = Make2D("h_nKFtrack_P", "Number of KFTracks per one mcTrack vs P", "N_{kfTracks}", "P, GeV/c",
                               &fList, 21, -0.5, 20.5, 1000, 0., 5.);
   TH2D *h_nDSTtrack_P = Make2D("h_nDSTtrack_P", "Number of DSTTracks per one mcTrack vs P", "N_{dstTracks}",
                                "P, GeV/c", &fList, 21, -0.5, 20.5, 1000, 0., 5.);

   // barrel M2 vs dEdX
   TH2D *hProb_e  = Make2D("hProb_e", "", "Prob_{M2}", "Prob_{dEdX}", &fBarrel, 1100, -1.1, 1.1, 1100, -1.1, 1.1);
   TH2D *hProb_p  = Make2D("hProb_p", hProb_e, &fBarrel);
   TH2D *hProb_pi = Make2D("hProb_pi", hProb_e, &fBarrel);
   TH2D *hProb_K  = Make2D("hProb_K", hProb_e, &fBarrel);

   // barrel  dEdX VS p
   TH2D *hProbP_e  = Make2D("hPProbdEdX_e", "", "P, GeV/c", "Prob_{dEdX}", &fBarrel, 1000, 0., 5., 1100, -1.1, 1.1);
   TH2D *hProbP_p  = Make2D("hPProbdEdX_p", hProbP_e, &fBarrel);
   TH2D *hProbP_pi = Make2D("hPProbdEdX_pi", hProbP_e, &fBarrel);
   TH2D *hProbP_K  = Make2D("hPProbdEdX_K", hProbP_e, &fBarrel);
   // barrel  M2 VS p
   TH2D *hProbP1_e  = Make2D("hPProbM2_e", "", "P, GeV/c", "Prob_{M2}", &fBarrel, 1000, 0, 5., 1100, -1.1, 1.1);
   TH2D *hProbP1_p  = Make2D("hPProbM2_p", hProbP1_e, &fBarrel);
   TH2D *hProbP1_pi = Make2D("hPProbM2_pi", hProbP1_e, &fBarrel);
   TH2D *hProbP1_K  = Make2D("hPProbM2_K", hProbP1_e, &fBarrel);

   // endcap M2 vs dEdX
   TH2D *hEProb_e  = Make2D("hEProb_e", hProb_e, &fEndcap);
   TH2D *hEProb_p  = Make2D("hEProb_p", hProb_e, &fEndcap);
   TH2D *hEProb_pi = Make2D("hEProb_pi", hProb_e, &fEndcap);
   TH2D *hEProb_K  = Make2D("hEProb_K", hProb_e, &fEndcap);

   // endcap  dEdX VS p
   TH2D *hEProbP_e  = Make2D("hEPProbdEdX_e", hProbP_e, &fEndcap);
   TH2D *hEProbP_p  = Make2D("hEPProbdEdX_p", hProbP_e, &fEndcap);
   TH2D *hEProbP_pi = Make2D("hEPProbdEdX_pi", hProbP_e, &fEndcap);
   TH2D *hEProbP_K  = Make2D("hEPProbdEdX_K", hProbP_e, &fEndcap);
   // endcap  M2 VS p
   TH2D *hEProbP1_e  = Make2D("hEPProbM2_e", hProbP1_e, &fEndcap);
   TH2D *hEProbP1_p  = Make2D("hEPProbM2_p", hProbP1_e, &fEndcap);
   TH2D *hEProbP1_pi = Make2D("hEPProbM2_pi", hProbP1_e, &fEndcap);
   TH2D *hEProbP1_K  = Make2D("hEPProbM2_K", hProbP1_e, &fEndcap);

   // Matching
   TH1D *hEta_MC        = Make1D("Eta_MC", "", "#||{#eta}", "", &fList, 15, -1.5, 1.5);
   TH1D *hEta_KF        = Make1D("Eta_KF", hEta_MC, &fList);
   TH1D *hEta_Matching  = Make1D("Eta_Matching", hEta_MC, &fList);
   TH1D *hEta_Matching1 = Make1D("Eta_Matching1", hEta_MC, &fList);

   TH1D *hP_MC        = Make1D("P_MC", "", "P, GeV/c", "Events", &fList, 30, 0., 3.);
   TH1D *hP_KF        = Make1D("P_KF", hP_MC, &fList);
   TH1D *hP_Matching  = Make1D("P_Matching", hP_MC, &fList);
   TH1D *hP_Matching1 = Make1D("P_Matching1", hP_MC, &fList);

   TH2D *hEta_P_MC            = Make2D("Eta_P_MC", "", "#||{#eta}", "P, GeV/c", &fList, 1000, -1.6, 1.6, 1000, 0., 3.);
   TH2D *hEta_P_KF            = Make2D("Eta_P_KF", hEta_P_MC, &fList);
   TH2D *hEta_P_Matching      = Make2D("Eta_P_Matching", hEta_P_MC, &fList);
   TH2D *hEta_P_Matching_lost = Make2D("Eta_P_Matching_lost", hEta_P_MC, &fList);
   TH2D *hEta_P_MatchingKF_lost = Make2D("Eta_P_MatchingKF_lost", hEta_P_MC, &fList);

   // EMatching
   TH1D *hEEta_MC        = Make1D("E_Eta_MC", "", "#||{#eta}", "", &fList, 30, 0., 6.);
   TH1D *hEEta_KF        = Make1D("E_Eta_KF", hEEta_MC, &fList);
   TH1D *hEEta_Matching  = Make1D("E_Eta_Matching", hEEta_MC, &fList);
   TH1D *hEEta_Matching1 = Make1D("E_Eta_Matching1", hEEta_MC, &fList);

   TH1D *hEP_MC        = Make1D("E_P_MC", "", "P, GeV/c", "Events", &fList, 30, 0., 3.);
   TH1D *hEP_KF        = Make1D("E_P_KF", hEP_MC, &fList);
   TH1D *hEP_Matching  = Make1D("E_P_Matching", hEP_MC, &fList);
   TH1D *hEP_Matching1 = Make1D("E_P_Matching1", hEP_MC, &fList);

   TH2D *hEEta_P_MC            = Make2D("E_Eta_P_MC", "", "#||{#eta}", "P, GeV/c", &fList, 1000, 0., 6.1, 1000, 0., 3.);
   TH2D *hEEta_P_KF            = Make2D("E_Eta_P_KF", hEEta_P_MC, &fList);
   TH2D *hEEta_P_Matching      = Make2D("E_Eta_P_Matching", hEEta_P_MC, &fList);
   TH2D *hEEta_P_Matching_lost = Make2D("E_Eta_P_Matching_lost", hEEta_P_MC, &fList);
   TH2D *hEEta_P_MatchingKF_lost = Make2D("E_Eta_P_MatchingKF_lost", hEEta_P_MC, &fList);
   // -----------------------------------------------------------------------------------------------
   TString OutFlnm = "result.root";

   TChain *treeMC = new TChain("mpdsim");
   treeMC->Add("evetest.root");
   //	treeMC->Add("evetest_1.root");

   TChain *treeDST = new TChain("mpdsim");
   treeDST->Add("mpddst.root");

   Int_t nDSTentries = treeDST->GetEntries();
   Int_t nMCentries  = treeMC->GetEntries();
   cout << "\n   MC events: " << nMCentries << ", reco events: " << nDSTentries;

   // Activate branches  - DST
   TClonesArray *aTofHit = NULL;
   treeDST->SetBranchAddress("TOFHit", &aTofHit);
   TClonesArray *aETofHit = NULL;
   treeDST->SetBranchAddress("ETOFHit", &aETofHit);
   TClonesArray *aTofMatching = NULL;
   treeDST->SetBranchAddress("TOFMatching", &aTofMatching);
   TClonesArray *aETofMatching = NULL;
   treeDST->SetBranchAddress("ETOFMatching", &aETofMatching);
   TClonesArray *aKFTrack = NULL;
   treeDST->SetBranchAddress("TpcKalmanTrack", &aKFTrack);
   TClonesArray *aEKFTrack = NULL;
   treeDST->SetBranchAddress("EctTrack", &aEKFTrack);
   MpdEvent *pMPDevent = NULL;
   treeDST->SetBranchAddress("MPDEvent.", &pMPDevent);

   // Activate branches  - MC
   TClonesArray *aTofPoint = NULL;
   treeMC->SetBranchAddress("TOFPoint", &aTofPoint);
   TClonesArray *aETofPoint = NULL;
   treeMC->SetBranchAddress("ETOFPoint", &aETofPoint);
   TClonesArray *aMcTrack = NULL;
   treeMC->SetBranchAddress("MCTrack", &aMcTrack);

   MpdTofMatchingData * pTofMatching;
   MpdTofHit *          pTofHit;
   MpdTofPoint *        pTofPoint;
   MpdEtofMatchingData *pETofMatching;
   MpdEtofHit *         pETofHit;
   MpdEtofPoint *       pETofPoint;

   MpdKalmanTrack *pKFtrack;
   FairMCTrack *   pMCtrack;
   TClonesArray *  aDSTtrack;
   MpdTrack *      pDSTtrack;

   mpdTestTrack    TRACK;
   DstTrackDefIter ItDST;
   EventDef        mapEVENT;
   EventDefIter    It;

   TofHitDef        mmapHITS;
   TofHitDefIter    ItHit;
   TofMatchDef      mmapMATCHS;
   TofMatchDefIter  ItMatch;
   ETofHitDef       mmapEHITS;
   ETofHitDefIter   ItEHit;
   ETofMatchDef     mmapEMATCHS;
   ETofMatchDefIter ItEMatch;

   Int_t    mcTrackIndex, PDGcode, counter, flag;
   Double_t P, eta, Pm2, Pdedx;
   TVector3 momentum;
   //	nDSTentries = 10;
   for (int event = 0; event < nDSTentries; event++) {
      treeDST->GetEntry(event);
      treeMC->GetEntry(event);

      mapEVENT.clear();

      mmapHITS.clear();
      mmapMATCHS.clear();
      mmapEHITS.clear();
      mmapEMATCHS.clear();

      Int_t nMC = aMcTrack->GetEntries();

      aDSTtrack       = pMPDevent->GetGlobalTracks();
      Int_t nDSTTrack = aDSTtrack->GetEntriesFast();

      Int_t nKF = 0, nEKF = 0;
      if (aKFTrack) nKF = aKFTrack->GetEntriesFast();
      if (aEKFTrack) nEKF = aEKFTrack->GetEntriesFast();

      Int_t nMatch = 0, nPoint = 0, nHit = 0, nEMatch = 0, nEPoint = 0, nEHit = 0;
      if (aTofPoint) nPoint = aTofPoint->GetEntries();
      if (aTofHit) nHit = aTofHit->GetEntries();
      if (aTofMatching) nMatch = aTofMatching->GetEntries();
      if (aETofPoint) nEPoint = aETofPoint->GetEntries();
      if (aETofHit) nEHit = aETofHit->GetEntries();
      if (aETofMatching) nEMatch = aETofMatching->GetEntries();

      cout << "\n Event " << event << " have " << nMC << " MC tracks, " << nKF + nEKF << " KFtracks, " << nDSTTrack
           << " DST tracks, " << nPoint << " Points, " << nHit << " Hits, " << nMatch << " Matching, " << nEPoint
           << " EPoints, " << nEHit << " EHits, " << nEMatch << " EMatching.";

      for (Int_t dstTrack = 0; dstTrack < nDSTTrack; dstTrack++) {
         pDSTtrack    = (MpdTrack *)aDSTtrack->UncheckedAt(dstTrack); // DST
         mcTrackIndex = pDSTtrack->GetID();

         It = mapEVENT.find(mcTrackIndex);
         if (It == mapEVENT.end())                                                   // first time
            It = (mapEVENT.insert(EventDef::value_type(mcTrackIndex, TRACK))).first; // add blank slote, update iter

         It->second.aDSTTrack.insert(DstTrackDef::value_type(pDSTtrack)); // add MpdTrack
      }

      //-------------------------- Kalman tracks --------------------------
      for (Int_t kfTrack = 0; kfTrack < nKF; kfTrack++) // barrel
      {
         pKFtrack     = (MpdKalmanTrack *)aKFTrack->UncheckedAt(kfTrack);
         mcTrackIndex = pKFtrack->GetTrackID();

         It = mapEVENT.find(mcTrackIndex);
         if (It == mapEVENT.end())                                                   // first time
            It = (mapEVENT.insert(EventDef::value_type(mcTrackIndex, TRACK))).first; // add blank slote, update iter

         It->second.aKFTrack.insert(KfTrackDef::value_type(pKFtrack)); // add MpdKalmanTrack
      }

      for (Int_t kfTrack = 0; kfTrack < nEKF; kfTrack++) // endcap
      {
         pKFtrack     = (MpdKalmanTrack *)aEKFTrack->UncheckedAt(kfTrack);
         mcTrackIndex = pKFtrack->GetTrackID();

         It = mapEVENT.find(mcTrackIndex);
         if (It == mapEVENT.end())                                                   // first time
            It = (mapEVENT.insert(EventDef::value_type(mcTrackIndex, TRACK))).first; // add blank slote, update iter

         It->second.aKFTrack.insert(KfTrackDef::value_type(pKFtrack)); // add MpdKalmanTrack
      }
      //------------------------------------------------- MC -------------------------------------------------
      for (Int_t tofPoint = 0; tofPoint < nPoint; tofPoint++) {
         pTofPoint    = (MpdTofPoint *)aTofPoint->UncheckedAt(tofPoint);
         mcTrackIndex = pTofPoint->GetTrackID();

         It = mapEVENT.find(mcTrackIndex);
         if (It == mapEVENT.end())                                                   // first time
            It = (mapEVENT.insert(EventDef::value_type(mcTrackIndex, TRACK))).first; // add blank slote, update iter

         It->second.aTofPoint.insert(TofPointDef::value_type(pTofPoint)); // add MpdTofPoint
      }

      for (Int_t etofPoint = 0; etofPoint < nEPoint; etofPoint++) {
         pETofPoint   = (MpdEtofPoint *)aETofPoint->UncheckedAt(etofPoint);
         mcTrackIndex = pETofPoint->GetTrackID();

         It = mapEVENT.find(mcTrackIndex);
         if (It == mapEVENT.end())                                                   // first time
            It = (mapEVENT.insert(EventDef::value_type(mcTrackIndex, TRACK))).first; // add blank slote, update iter

         It->second.aETofPoint.insert(ETofPointDef::value_type(pETofPoint)); // add MpdEtofPoint
      }

      //------------------------------------------------- Hits -------------------------------------------------
      Int_t    nLinks;
      FairLink link;
      for (Int_t tofHit = 0; tofHit < nHit; tofHit++) // mmapHITS.size >= nHit !!!!!!!!!!!!!!!
      {
         pTofHit = (MpdTofHit *)aTofHit->UncheckedAt(tofHit);

         if (nLinks = pTofHit->GetNLinks())
            for (Int_t i = 0; i < nLinks; i++) // cycle by links
            {
               link = pTofHit->GetLink(i);
               if (link.GetType() == MpdTofUtils::IsMCTrackIndex)
                  mmapHITS.insert(TofHitDef::value_type(link.GetIndex(), pTofHit)); // key - mcTrackIndex, add MpdTofHit
            }
      }

      for (Int_t etofHit = 0; etofHit < nEHit; etofHit++) // mmapHITS.size >= nHit !!!!!!!!!!!!!!!
      {
         pETofHit = (MpdEtofHit *)aETofHit->UncheckedAt(etofHit);

         if (nLinks = pETofHit->GetNLinks())
            for (Int_t i = 0; i < nLinks; i++) // cycle by links
            {
               link = pETofHit->GetLink(i);
               if (link.GetType() == MpdTofUtils::IsMCTrackIndex)
                  mmapEHITS.insert(
                     ETofHitDef::value_type(link.GetIndex(), pETofHit)); // key - mcTrackIndex, add MpdETofHit
            }
      }
      //------------------------------------------------- Matchings -------------------------------------------------
      for (Int_t tofMatching = 0; tofMatching < nMatch; tofMatching++) {
         pTofMatching = (MpdTofMatchingData *)aTofMatching->At(tofMatching);
         pKFtrack     = (MpdKalmanTrack *)aKFTrack->UncheckedAt(pTofMatching->GetKFTrackIndex());
         mcTrackIndex = pKFtrack->GetTrackID();

         mmapMATCHS.insert(TofMatchDef::value_type(mcTrackIndex, pTofMatching)); // add MpdTofMatchingData
      }

      for (Int_t etofMatching = 0; etofMatching < nEMatch; etofMatching++) {
         pETofMatching = (MpdEtofMatchingData *)aETofMatching->At(etofMatching);
         pKFtrack      = (MpdKalmanTrack *)aEKFTrack->UncheckedAt(pETofMatching->GetKFTrackIndex());
         mcTrackIndex  = pKFtrack->GetTrackID();

         mmapEMATCHS.insert(ETofMatchDef::value_type(mcTrackIndex, pETofMatching)); // add MpdTofMatchingData
      }
      //----------------------------------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------------------------------

      // Number of MpdTofHits per one mcTrack vs P
      for (TofHitDefIter it = mmapHITS.begin(); it != mmapHITS.end();) {
         counter = mmapHITS.count(it->first);
         h_nHits_P->Fill(counter, ((FairMCTrack *)aMcTrack->UncheckedAt(it->first))->GetP());
         for (int i = 0; i < counter; i++) ++it;
      }

      // Number of MpdTofMatchingDatas per one mcTrack vs P
      for (TofMatchDefIter it = mmapMATCHS.begin(); it != mmapMATCHS.end();) {
         counter = mmapMATCHS.count(it->first);
         h_nMatchs_P->Fill(counter, ((FairMCTrack *)aMcTrack->UncheckedAt(it->first))->GetP());
         for (int i = 0; i < counter; i++) ++it;
      }

      // Number of MpdETofHits per one mcTrack vs P
      for (ETofHitDefIter it = mmapEHITS.begin(); it != mmapEHITS.end();) {
         counter = mmapEHITS.count(it->first);
         h_nEHits_P->Fill(counter, ((FairMCTrack *)aMcTrack->UncheckedAt(it->first))->GetP());
         for (int i = 0; i < counter; i++) ++it;
      }

      // Number of MpdETofMatchingDatas per one mcTrack vs P
      for (ETofMatchDefIter it = mmapEMATCHS.begin(); it != mmapEMATCHS.end();) {
         counter = mmapEMATCHS.count(it->first);
         h_nEMatchs_P->Fill(counter, ((FairMCTrack *)aMcTrack->UncheckedAt(it->first))->GetP());
         for (int i = 0; i < counter; i++) ++it;
      }

      //----------------------------------------------------------------------------------------------------------
      for (It = mapEVENT.begin(); It != mapEVENT.end(); ++It) // cycle by selected MC tracks
      {
         mcTrackIndex = It->first;
         pMCtrack     = (FairMCTrack *)aMcTrack->At(mcTrackIndex);
         PDGcode      = pMCtrack->GetPdgCode();
         pMCtrack->GetMomentum(momentum);
         P   = pMCtrack->GetP();
         eta = momentum.Eta();

         h_nPoints_P->Fill(It->second.aTofPoint.size(), P);   // Number of TofPoints per one mcTrack vs P
         h_nEPoints_P->Fill(It->second.aETofPoint.size(), P); // Number of EtofPoints per one mcTrack vs P

         h_nKFtrack_P->Fill(It->second.aKFTrack.size(), P);   // Number of KFTracks per one mcTrack vs P
         h_nDSTtrack_P->Fill(It->second.aDSTTrack.size(), P); // Number of DSTTracks per one mcTrack vs P

         counter = It->second.aDSTTrack.size();
         if (counter >= 1) // mc track have dst track(s)
         {
            for (ItDST = It->second.aDSTTrack.begin(); ItDST != It->second.aDSTTrack.end(); ++ItDST) {

               flag = (*ItDST)->GetTofFlag();

               if (PDGcode == 11 || PDGcode == -11) {
                  Pm2   = (*ItDST)->GetPidProbElectron();
                  Pdedx = (*ItDST)->GetTPCPidProbElectron();

                  if (flag & BIT(1)) hProbP1_e->Fill(P, Pm2);  // tof barrel
                  if (flag & BIT(2)) hProbP_e->Fill(P, Pdedx); // tpc barrel
                  if ((flag & BIT(1)) && (flag & BIT(2))) hProb_e->Fill(Pm2, Pdedx);

                  if (flag & BIT(3)) hEProbP1_e->Fill(P, Pm2);  // tof endcap
                  if (flag & BIT(4)) hEProbP_e->Fill(P, Pdedx); // tpc endcap
                  if ((flag & BIT(3)) && (flag & BIT(4))) hEProb_e->Fill(Pm2, Pdedx);
               } else if (PDGcode == 211 || PDGcode == -211) {
                  Pm2   = (*ItDST)->GetPidProbPion();
                  Pdedx = (*ItDST)->GetTPCPidProbPion();

                  if (flag & BIT(1)) hProbP1_pi->Fill(P, Pm2);
                  if (flag & BIT(2)) hProbP_pi->Fill(P, Pdedx);
                  if ((flag & BIT(1)) && (flag & BIT(2))) hProb_pi->Fill(Pm2, Pdedx);

                  if (flag & BIT(3)) hEProbP1_pi->Fill(P, Pm2);
                  if (flag & BIT(4)) hEProbP_pi->Fill(P, Pdedx);
                  if ((flag & BIT(3)) && (flag & BIT(4))) hEProb_pi->Fill(Pm2, Pdedx);
               } else if (PDGcode == 2212 || PDGcode == -2212) {
                  Pm2   = (*ItDST)->GetPidProbProton();
                  Pdedx = (*ItDST)->GetTPCPidProbProton();

                  if (flag & BIT(1)) hProbP1_p->Fill(P, Pm2);
                  if (flag & BIT(2)) hProbP_p->Fill(P, Pdedx);
                  if ((flag & BIT(1)) && (flag & BIT(2))) hProb_p->Fill(Pm2, Pdedx);

                  if (flag & BIT(3)) hEProbP1_p->Fill(P, Pm2);
                  if (flag & BIT(4)) hEProbP_p->Fill(P, Pdedx);
                  if ((flag & BIT(3)) && (flag & BIT(4))) hEProb_p->Fill(Pm2, Pdedx);
               } else if (PDGcode == 321 || PDGcode == -321) {
                  Pm2   = (*ItDST)->GetPidProbKaon();
                  Pdedx = (*ItDST)->GetTPCPidProbKaon();

                  if (flag & BIT(1)) hProbP1_K->Fill(P, Pm2);
                  if (flag & BIT(2)) hProbP_K->Fill(P, Pdedx);
                  if ((flag & BIT(1)) && (flag & BIT(2))) hProb_K->Fill(Pm2, Pdedx);

                  if (flag & BIT(3)) hEProbP1_K->Fill(P, Pm2);
                  if (flag & BIT(4)) hEProbP_K->Fill(P, Pdedx);
                  if ((flag & BIT(3)) && (flag & BIT(4))) hEProb_K->Fill(Pm2, Pdedx);
               }
            }
         }

         // Matching
         counter = It->second.aTofPoint.size();
         if (counter >= 1) // mc track have tof point(s)
         {
            hEta_MC->Fill(momentum.Eta(), counter);
            hP_MC->Fill(P, counter);
            hEta_P_MC->Fill(eta, P, counter);

            counter = mmapMATCHS.count(mcTrackIndex);
            if (counter >= 1) {
               hEta_Matching->Fill(eta, counter);
               hP_Matching->Fill(P, counter);
               hEta_P_Matching->Fill(eta, P, counter);
            } else
               hEta_P_MatchingKF_lost->Fill(eta, P);

            counter = It->second.aKFTrack.size();
            if (counter >= 1) // mc track have tof point(s) & KF track(s)
            {
               hEta_KF->Fill(eta, counter);
               hP_KF->Fill(P, counter);
               hEta_P_KF->Fill(eta, P, counter);

               counter = mmapMATCHS.count(mcTrackIndex);
               if (counter >= 1) // mc track have tof point(s) & KF track(s) & Matching(s)
               {
                  hEta_Matching1->Fill(momentum.Eta(), counter);
                  hP_Matching1->Fill(P, counter);
               } else
                  hEta_P_Matching_lost->Fill(eta, P);
            }
         }

         // EMatching
         counter = It->second.aETofPoint.size();
         if (counter >= 1) // mc track have etof point(s)
         {
            hEEta_MC->Fill(momentum.Eta(), counter);
            hEP_MC->Fill(P, counter);
            hEEta_P_MC->Fill(eta, P, counter);

            counter = mmapEMATCHS.count(mcTrackIndex);
            if (counter >= 1) {
               hEEta_Matching->Fill(eta, counter);
               hEP_Matching->Fill(P, counter);
               hEEta_P_Matching->Fill(eta, P, counter);
            } else
               hEEta_P_MatchingKF_lost->Fill(eta, P);

            counter = It->second.aKFTrack.size();
            if (counter >= 1) // mc track have etof point(s) & KF track(s)
            {
               hEEta_KF->Fill(eta, counter);
               hEP_KF->Fill(P, counter);
               hEEta_P_KF->Fill(eta, P, counter);

               counter = mmapEMATCHS.count(mcTrackIndex);
               if (counter >= 1) // mc track have etof point(s) & KF track(s) & Matching(s)
               {
                  hEEta_Matching1->Fill(momentum.Eta(), counter);
                  hEP_Matching1->Fill(P, counter);
               } else
                  hEEta_P_Matching_lost->Fill(eta, P);
            }
         }

         //	counter = It->second.aKFTrack.size();

         //	counter = mmapMATCHS.count(mcTrackIndex);
      }
      //----------------------------------------------------------------------------------------------------------

   } // event

   //----------------------------------------------------------------------------------------------------------
   //----------------------------------------------------------------------------------------------------------
   //----------------------------------------------------------------------------------------------------------
   // Matching
   TH1D *h1 = (TH1D *)hP_Matching->Clone("P_Matching_Eff");
   h1->Divide(hP_KF);
   fList.Add(h1);
   h1 = (TH1D *)hEta_Matching->Clone("Eta_Matching_Eff");
   h1->Divide(hEta_KF);
   fList.Add(h1);

   h1 = (TH1D *)hP_Matching1->Clone("P_Matching1_Eff");
   h1->Divide(hP_KF);
   fList.Add(h1);
   h1 = (TH1D *)hEta_Matching1->Clone("Eta_Matching1_Eff");
   h1->Divide(hEta_KF);
   fList.Add(h1);

   // EMatching
   h1 = (TH1D *)hEP_Matching->Clone("E_P_Matching_Eff");
   h1->Divide(hEP_KF);
   fList.Add(h1);
   h1 = (TH1D *)hEEta_Matching->Clone("E_Eta_Matching_Eff");
   h1->Divide(hEEta_KF);
   fList.Add(h1);

   h1 = (TH1D *)hEP_Matching1->Clone("E_P_Matching1_Eff");
   h1->Divide(hEP_KF);
   fList.Add(h1);
   h1 = (TH1D *)hEEta_Matching1->Clone("E_Eta_Matching1_Eff");
   h1->Divide(hEEta_KF);
   fList.Add(h1);

   TFile file(OutFlnm.Data(), "RECREATE");
   fList.Write();
   TDirectory *bar = file.mkdir("barrel");
   bar->cd();
   fBarrel.Write();
   TDirectory *end = file.mkdir("endcap");
   end->cd();
   fEndcap.Write();
   file.Close();
}
//--------------------------------------------------------------------------------------------------------------
int main()
{
   TStopwatch timer;
   timer.Start();

   doTest();

   cout << endl;
   timer.Print();
   cout << " All ok." << endl;

   return 0;
}
//--------------------------------------------------------------------------------------------------------------
TH2D *Make2D(const char *name, const char *title, const char *Xtitle, const char *Ytitle, TList *list, Int_t nbinx,
             Double_t xmin, Double_t xmax, Int_t nbiny, Double_t ymin, Double_t ymax)
{
   TH2D *ptr = new TH2D(name, title, nbinx, xmin, xmax, nbiny, ymin, ymax);
   ptr->GetXaxis()->SetTitle(Xtitle);
   ptr->GetYaxis()->SetTitle(Ytitle);
   if (list) {
      list->Add(ptr);
      ptr->SetDirectory(0);
   }
   return ptr;
}
//------------------------------------------------------------------------------------------------------------------------
TH1D *Make1D(const char *name, const char *title, const char *Xtitle, const char *Ytitle, TList *list, Int_t nbinx,
             Double_t xmin, Double_t xmax)
{
   TH1D *ptr = new TH1D(name, title, nbinx, xmin, xmax);
   ptr->GetXaxis()->SetTitle(Xtitle);
   ptr->GetYaxis()->SetTitle(Ytitle);
   if (list) {
      list->Add(ptr);
      ptr->SetDirectory(0);
   }
   return ptr;
}
//------------------------------------------------------------------------------------------------------------------------
TH2D *Make2D(const char *name, const TH2D *h2, TList *list)
{
   TH2D *ptr = (TH2D *)h2->Clone();
   ptr->SetName(name);
   if (list) {
      list->Add(ptr);
      ptr->SetDirectory(0);
   }
   return ptr;
}
//------------------------------------------------------------------------------------------------------------------------
TH1D *Make1D(const char *name, const TH1D *h1, TList *list)
{
   TH1D *ptr = (TH1D *)h1->Clone();
   ptr->SetName(name);
   if (list) {
      list->Add(ptr);
      ptr->SetDirectory(0);
   }
   return ptr;
}
