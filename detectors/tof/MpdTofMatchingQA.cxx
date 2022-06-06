//------------------------------------------------------------------------------------------------------------------------
#include <map>
#include <iostream>
#include <algorithm>

#include <TGraph.h>
#include <TClonesArray.h>
#include <TFile.h>
#include "FairLogger.h"
#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTof.h"
#include "MpdTofMatching.h"

#include "MpdTofMatchingQA.h"
using namespace std;

//------------------------------------------------------------------------------------------------------------------------
void Lsets::Reset()
{
   for_each(fData.begin(), fData.end(), [](std::set<size_t> &s) { s.clear(); });
}
//------------------------------------------------------------------------------------------------------------------------
void Lsets::Insert(size_t index, size_t value)
{
   fData.at(index).insert(value);
}
//------------------------------------------------------------------------------------------------------------------------
void Lsets::Print(const char *comment, ostream &os) const
{
   if (comment != nullptr) os << comment;

   Ts unionSet;
   for (auto it = fData.begin(), itEnd = fData.end(); it != itEnd; it++) unionSet.insert((*it).begin(), (*it).end());

   os << "\n sizes: ";
   size_t n = 0;
   for (auto it = fData.begin(), itEnd = fData.end(); it != itEnd; it++) os << n++ << ":(" << (*it).size() << "), ";
   os << "\tall=(" << unionSet.size() << ")";

   for (auto it = unionSet.begin(), itEnd = unionSet.end(); it != itEnd; it++) {
      os << "\n (" << (*it) << ") ";
      for (size_t i = 0, iEnd = fData.size(); i != iEnd; i++) {
         if (fData.at(i).find(*it) != fData.at(i).end())
            os << "  +";
         else
            os << "  -";
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
MpdTofMatchingQA::MpdTofMatchingQA(const char *flnm, bool isEndcap) : fFlnm(flnm), fIsEndcap(isEndcap)
{
   const double Pmax   = 5.; // [GeV/c]
   const double EtaMax = 5.;

   fList.SetOwner(); //  all objects will be deleted whenever the collection itself is delete.

   Add(pEff1P = new TEfficiency(mangling("Eff1P"), ";P, GeV/c;Efficiency", 100, 0., Pmax));
   Add(pEff1Eta = new TEfficiency(mangling("Eff1Eta"), ";#eta;Efficiency", 100, -EtaMax, EtaMax));
   Add(pEff1EtaP = new TEfficiency(mangling("Eff1EtaP"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pCont1P = new TEfficiency(mangling("Cont1P"), ";P, GeV/c;Contamination", 100, 0., Pmax));
   Add(pCont1Eta = new TEfficiency(mangling("Cont1Eta"), ";#eta;Contamination", 100, -EtaMax, EtaMax));
   Add(pCont1EtaP = new TEfficiency(mangling("Cont1EtaP"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));

   Add(pEff2P = new TEfficiency(mangling("Eff2P"), ";P, GeV/c;Efficiency", 100, 0., Pmax));
   Add(pEff2Eta = new TEfficiency(mangling("Eff2Eta"), ";#eta;Efficiency", 100, -EtaMax, EtaMax));
   Add(pEff2EtaP = new TEfficiency(mangling("Eff2EtaP"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pCont2P = new TEfficiency(mangling("Cont2P"), ";P, GeV/c;Contamination", 100, 0., Pmax));
   Add(pCont2Eta = new TEfficiency(mangling("Cont2Eta"), ";#eta;Contamination", 100, -EtaMax, EtaMax));
   Add(pCont2EtaP = new TEfficiency(mangling("Cont2EtaP"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));

   Add(pEffEtaPIdeal = (TEfficiency *)pEff1P->Clone(mangling("Eff_EtaP_Ideal")));
   pEffEtaPIdeal->SetTitle("ideal matching by using MC data");
   Add(pContEtaPIdeal = (TEfficiency *)pCont1P->Clone(mangling("Cont_EtaP_Ideal")));
   pContEtaPIdeal->SetTitle("ideal matching by using MC data");

   Add(hDeltaPoint = new TH2D(mangling("DeltaPoint"), "MC point - est. point;#Delta_{Z}, cm;#Delta_{#phi}, cm", 1000,
                              -300., 300., 1000, -300., 300.));
   Add(hDeltaPoint_Dev = new TH2D(mangling("DeltaPoint_Dev"), "MC point - est. point;#Delta, cm;P, GeV/c", 1000, 0.,
                                  10., 1000, 0., Pmax));
   Add(hDeltaPoint_dR = new TH2D(mangling("DeltaPoint_dR"), "MC point - est. point;#Delta_{R}, cm;P, GeV/c", 1000, -10.,
                                 10., 1000, 0., Pmax));
   Add(hDeltaPoint_dZ = new TH2D(mangling("DeltaPoint_dZ"), "MC point - est. point;#Delta_{Z}, cm;P, GeV/c", 1000, -10.,
                                 10., 1000, 0., Pmax));
   Add(hDeltaPoint_dPhi = new TH2D(mangling("DeltaPoint_dPhi"), "MC point - est. point;#Delta_{#phi}, cm;P, GeV/c",
                                   1000, -10., 10., 1000, 0., Pmax));

   Add(hDeltaHitTrue = new TH2D(mangling("DeltaHitTrue"), "est. point <-> Hit point;#Delta_{Z}, cm;#Delta_{#phi}, cm",
                                1000, -25., 25., 1000, -25., 25.));
   Add(hDeltaHitMis = (TH2D *)hDeltaHitTrue->Clone(mangling("DeltaHitMis")));

   Add(hNMcKfTracks = new TH2D(mangling("NMcKfTracks"), "count tracks with Tof points only;N mc tacks;N kf tracks",
                               1000, 0.5, 2000.5, 1000, 0.5, 2000.5));
   Add(htKFTrack = new TH2D(mangling("KFTrack"), ";#eta;P, GeV/c", 1000, -EtaMax, EtaMax, 1000, 0., Pmax));
   htKFTrack->SetTitle("[MpdTofMatching::FillWeightMatrix] cycle by TmPt& mPt");
   Add(htKFTrackPoint = (TH2D *)htKFTrack->Clone(mangling("KFTrackPoint")));
   htKFTrackPoint->SetTitle("[MpdTofMatching::FillWeightMatrix] cycle by TmPt& mPt(cuts: mcInfo::TofTouch=true)");
   Add(htMisMatch = (TH2D *)htKFTrack->Clone(mangling("MisMatchings")));
   htMisMatch->SetTitle("[MpdTofMatchingQA::FillMatchingEfficiency] cycle by aTPCkfTracks(cuts: IsMatchingExist && "
                        "MpdTofHit::CheckTid(tid)=false)");
   Add(htTrueMatch = (TH2D *)htKFTrack->Clone(mangling("TrueMatchings")));
   htTrueMatch->SetTitle("[MpdTofMatchingQA::FillMatchingEfficiency] cycle by aTPCkfTracks(cuts: IsMatchingExist && "
                         "MpdTofHit::CheckTid(tid)=true)");

   Add(hPointInsideDetector = new TH2D(mangling("test_IsPointInsideDetectors_function"), "; Z, cm; #Phi, degree", 1000,
                                       -300., 300., 1000, -180., 180.));

   // matching weight parameter tests
   const double maxWeight = 1.5E+3;
   Add(hWeightT = new TH1D(mangling("WeightTrue"), ";weight;Events", 10000, 0., maxWeight));
   Add(hWeightF = (TH1D *)hWeightT->Clone(mangling("WeightFalse")));
   Add(hNormWeightT = new TH1D(mangling("NormWeightTrue"), ";norm. weight;Events", 1000, 0., 1.01));
   Add(hNormWeightF = (TH1D *)hNormWeightT->Clone(mangling("NormWeightFalse")));

   // kf tracks tests
   const double Ptmax = 3., dPmin = -20., dPmax = 10.;
   Add(hPvsPt = new TH2D(mangling("dP_dPt"), "; #Delta P, %; #Delta Pt, %", 1000, dPmin, dPmax, 1000, dPmin, dPmax));
   Add(hPvsNp = new TH2D(mangling("dP_Np"), "; #Delta P, %; NofTrHits", 1000, dPmin, dPmax, 60, 0.5, 60.5));
   Add(hPtvsNp = new TH2D(mangling("Pt_Np"), "; #Delta Pt, %; NofTrHits", 1000, dPmin, dPmax, 60, 0.5, 60.5));
   Add(hPvsP =
          new TH2D(mangling("dP_P"), "; #Delta P, % (mc<est|0|mc>est); P, GeV/c", 1000, dPmin, dPmax, 1000, 0., Ptmax));
   Add(hPtvsPt = new TH2D(mangling("dPt_Pt"), "; #Delta Pt, % (mc<est|0|mc>est); Pt, GeV/c", 1000, dPmin, dPmax, 1000,
                          0., Ptmax));

   // smeared tracks tests
   Add(hSmeared_dZ = new TH2D(mangling("Smeared_dZ"), "; #Delta Z1 - Z2, cm; #Delta Z1 - est. point, cm", 1000, 0., 50.,
                              1000, 0., 25.));
   Add(hSmeared_dPhi = new TH2D(mangling("Smeared_dPhi"), "; #Delta #phi1 - #phi2, cm; #Delta #phi1 - est. point, cm",
                                1000, 0., 50., 1000, 0., 25.));
   Add(hSmeared_dZPhi =
          new TH2D(mangling("Smeared_dZ_dPhi"), "; #Delta Z, cm; #Delta #phi, cm", 1000, 0., 50., 1000, 0., 50.));

   Add(h1 = new TH2D(mangling("h1"), "; #Delta Z, cm; #Delta #phi, cm", 1000, -100., 100., 1000, -100., 100.));
   Add(h2 = new TH2D(mangling("h2"), "; #Delta Z, cm; #Delta #phi, cm", 1000, -100., 100., 1000, -100., 100.));
   Add(h3 = new TH2D(mangling("h3"), "; #Delta Z, cm; #Delta #phi, cm", 1000, -100., 100., 1000, -100., 100.));

   Add(pEffKalman = new TEfficiency(mangling("EffKalman"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pEffMatch = new TEfficiency(mangling("EffMatch"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));

   Add(pEffKalmanP = new TEfficiency(mangling("EffKalman_P"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pEffMatchP = new TEfficiency(mangling("EffMatch_P"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));

   Add(pEffKalmanPi = new TEfficiency(mangling("EffKalman_Pi"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pEffMatchPi = new TEfficiency(mangling("EffMatch_Pi"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));

   Add(pEffKalmanK = new TEfficiency(mangling("EffKalman_K"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
   Add(pEffMatchK = new TEfficiency(mangling("EffMatch_K"), ";#eta;P, GeV/c", 100, -EtaMax, EtaMax, 100, 0., Pmax));
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillParameter(bool isTrueMatching, double weight, double normWeight)
{
   if (isTrueMatching) {
      hWeightT->Fill(weight);
      hNormWeightT->Fill(normWeight);
   } else {
      hWeightF->Fill(weight);
      hNormWeightF->Fill(normWeight);
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillHitDeviation(const MpdTofMatchingData &data)
{
   double delta, deltaZ, deltaR, deltaPhi;
   MpdTof::GetDelta(data.fHitPosition, data.fEstPoint, delta, deltaZ, deltaR, deltaPhi);

   if (data.fIsTrueMatching)
      hDeltaHitTrue->Fill(deltaZ, deltaPhi);
   else
      hDeltaHitMis->Fill(deltaZ, deltaPhi);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillSmearedPoints(const TVector3 &estPosition, const TVector3 *smearedPosition,
                                         const TVector3 *smearedMomentum)
{
   const double dZ   = (smearedPosition[0] - smearedPosition[2]).Mag(),
                dPhi = (smearedPosition[1] - smearedPosition[3]).Mag();

   hSmeared_dZ->Fill(dZ, (estPosition - smearedPosition[0]).Mag());
   hSmeared_dPhi->Fill(dPhi, (estPosition - smearedPosition[1]).Mag());
   hSmeared_dZPhi->Fill(dZ, dPhi);

   double delta, deltaZ1, deltaZ2, deltaR, deltaPhi;
   MpdTof::GetDelta(estPosition, smearedPosition[0], delta, deltaZ1, deltaR, deltaPhi);
   h1->Fill(deltaZ1, deltaPhi);

   MpdTof::GetDelta(estPosition, smearedPosition[1], delta, deltaZ2, deltaR, deltaPhi);
   h2->Fill(deltaZ2, deltaPhi);

   MpdTof::GetDelta(smearedPosition[0], smearedPosition[1], delta, deltaZ1, deltaR, deltaPhi);
   h3->Fill(deltaZ1, deltaPhi);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillMcPointDeviation(const TVector3 &mcPosition, const TVector3 &estPointOnPlate,
                                            const TVector3 &mcMom, const TVector3 &kfMom, int NofTrHits)
{
   double delta, deltaZ, deltaR, deltaPhi;
   MpdTof::GetDelta(mcPosition, estPointOnPlate, delta, deltaZ, deltaR, deltaPhi);

   const double kfP = kfMom.Mag();
   hDeltaPoint_Dev->Fill(delta, kfP);
   hDeltaPoint_dR->Fill(deltaR, kfP);
   hDeltaPoint_dZ->Fill(deltaZ, kfP);
   hDeltaPoint_dPhi->Fill(deltaPhi, kfP);
   hDeltaPoint->Fill(deltaZ, deltaPhi);

   const double dP  = (mcMom.Mag() - kfMom.Mag()) / mcMom.Mag() * 100.;    //[%]
   const double dPt = (mcMom.Perp() - kfMom.Perp()) / mcMom.Perp() * 100.; //[%]

   hPvsPt->Fill(dP, dPt);
   hPvsNp->Fill(dP, NofTrHits);
   hPtvsNp->Fill(dPt, NofTrHits);

   hPvsP->Fill(dP, mcMom.Mag());
   hPtvsPt->Fill(dPt, mcMom.Perp());
}

//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillNtracks(const TClonesArray *aTracks, size_t NkfTracks)
{
   size_t NmcTracks = 0;
   for (Int_t i = 0, iMax = aTracks->GetEntriesFast(); i < iMax; i++)
      if (((MpdMCTrack *)aTracks->UncheckedAt(i))->GetNPoints(kTOF)) NmcTracks++;

   hNMcKfTracks->Fill(NmcTracks, NkfTracks);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillTrackEfficiency(bool kalman, bool matched, int pdgcode, const TVector3 &momentum)
{
   pEffKalman->Fill(kalman, momentum.Eta(), momentum.Mag());
   pEffMatch->Fill(matched, momentum.Eta(), momentum.Mag());

   if (2212 == pdgcode || -2212 == pdgcode) {
      pEffKalmanP->Fill(kalman, momentum.Eta(), momentum.Mag());
      pEffMatchP->Fill(matched, momentum.Eta(), momentum.Mag());
   } else if (211 == pdgcode || -211 == pdgcode) {
      pEffKalmanPi->Fill(kalman, momentum.Eta(), momentum.Mag());
      pEffMatchPi->Fill(matched, momentum.Eta(), momentum.Mag());
   } else if (321 == pdgcode || -321 == pdgcode) {
      pEffKalmanK->Fill(kalman, momentum.Eta(), momentum.Mag());
      pEffMatchK->Fill(matched, momentum.Eta(), momentum.Mag());
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::Finish()
{
   // fill weigh thresh. parameter graph
   double Ntracks = htKFTrackPoint->Integral();
   if (Ntracks > 0.) {
      double Eff[100], Cont[100], matchT, matchF;
      TAxis *xaxis   = hNormWeightT->GetXaxis();
      Int_t  binLast = xaxis->FindBin(1.);

      for (int i = 0; i < 100; i++) {
         double param    = 0.01 * (i + 1); // [0.01;1]
         Int_t  binThres = xaxis->FindBin(param);
         matchT          = hNormWeightT->Integral(binThres, binLast);
         matchF          = hNormWeightF->Integral(binThres, binLast);

         Eff[i]  = matchT / Ntracks;
         Cont[i] = matchF / (matchT + matchF);
         cout << "\n par=" << param << " Eff.= " << Eff[i] << " Cont.= " << Cont[i];
      }

      TGraph *graph = new TGraph(100, Eff, Cont);
      fList.Add(graph);
      graph->SetMarkerStyle(21);
      graph->SetMarkerColor(kRed);
      graph->SetName(mangling("Eff_vs_Cont(thresh)"));
      graph->SetTitle("Efficiency vs Contamination as function of weight threshould.");
      graph->GetHistogram()->SetXTitle("efficiency");
      graph->GetHistogram()->SetYTitle("contamination");
   }

   // write histo to file
   LOG(debug2) << "[MpdTofMatchingQA::Finish] Update  " << fFlnm.Data() << " file. ";
   auto  tmp = gFile;
   TFile file(fFlnm.Data(), "RECREATE");
   fList.Write();
   file.Close();
   gFile = tmp;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillIdealMatching(const TmPt &mPt, const TmmT2H &mmT2H, const TClonesArray *aMCTracks)
{
   TVector3 momentum;

   for (const auto &it : mPt) // cycle by TPC KF tracks
   {
      auto pKfTrack = (const MpdTpcKalmanTrack *)it.first;
      // Int_t KfIndex = it.second;
      Int_t mcTrackIndex = pKfTrack->GetTrackID();

      size_t Nhits = mmT2H.count(mcTrackIndex); // number of hits for current tid
      if (Nhits == 0) continue;

      auto pMCtrack = (MpdMCTrack *)aMCTracks->UncheckedAt(mcTrackIndex);
      pMCtrack->GetMomentum(momentum);
      // double Eta = momentum.Eta(),
      double P = momentum.Mag();

      bool IsTrueMatching = (Nhits == 1);
      pEffEtaPIdeal->Fill(IsTrueMatching, P);   //  = 1 hit per track
      pContEtaPIdeal->Fill(!IsTrueMatching, P); //  > 1 hit per track
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofMatchingQA::FillMatchingEfficiency(const TmPt &tids, const TClonesArray *aTofMatching,
                                              const TClonesArray *aTofHits, const TClonesArray *aMCTracks)
{
   // Sorting & mapping matchings
   map<int, const MpdTofHit *>      mMatchings;   // pair <kfTrackIndex, MpdTofHit*>
   multimap<const MpdTofHit *, int> mmHitCounter; // pair <MpdTofHit*, matchingIndex>

   for (int index = 0, iEnd = aTofMatching->GetEntriesFast(); index != iEnd; index++) // cycle by the  matching N
   {
      auto  pData    = (const MpdTofMatchingData *)aTofMatching->At(index);
      Int_t hitIndex = pData->GetTofHitIndex();

      if (hitIndex >= 0) // matching with the fired strip, i.e. with the hit
      {
         auto hit = (MpdTofHit *)aTofHits->At(hitIndex);

         mMatchings.insert(make_pair(pData->GetKFTrackIndex(), hit));
         mmHitCounter.insert(make_pair(hit, index));
      }
   }

   //--------------------------------------------------------------------------------------------------------------------
   // algorithm efficiency definition:
   //
   // efficiency = N true matchings / N tpc kf tracks having TOF hit
   // contamination = N wrong matchings / ( N true matchings + N wrong matchings)
   //
   // ALICE TDR TOF 2000 definitions:
   //
   // N = N mis + N match
   // N mis = N0 + N2 ( matched to to blank strip or dead regions + matched to fired strip was more than one track)
   // N match = N t + N w ( N true matching + N wrong matching)
   //
   // efficiency = N t / N
   // contamination = N w / N match
   //--------------------------------------------------------------------------------------------------------------------

   TVector3 momentum;
   for (const auto &iter : tids) // cycle by TPC KF tracks(Pt descending)
   {
      auto  pKfTrack = (const MpdTpcKalmanTrack *)iter.first;
      Int_t KfIndex  = iter.second;

      Int_t mcTrackIndex = pKfTrack->GetTrackID();
      auto  pMCtrack     = (MpdMCTrack *)aMCTracks->UncheckedAt(mcTrackIndex);
      bool  mcTofTouch   = pMCtrack->GetNPoints(fIsEndcap ? kETOF : kTOF);

      // if(mcTofTouch) tidsTofTouch.Insert(2, mcTrackIndex);

      pMCtrack->GetMomentum(momentum);
      double Eta = momentum.Eta(), P = momentum.Mag();
      auto   Iter            = mMatchings.find(KfIndex);
      bool   IsMatchingExist = (Iter != mMatchings.end());
      bool   IsTrueMatching  = false;
      bool   IsNmatch        = false;

      if (IsMatchingExist) {
         const MpdTofHit *hit = Iter->second;

         IsTrueMatching = hit->CheckTid(mcTrackIndex);
         IsNmatch       = (1 == mmHitCounter.count(hit)); // hit linked ONLY ONE matching

         if (IsTrueMatching)
            htTrueMatch->Fill(Eta, P);
         else
            htMisMatch->Fill(Eta, P);
      }

      if (mcTofTouch) {
         pEff1P->Fill(IsTrueMatching, P);
         pEff1Eta->Fill(IsTrueMatching, Eta);
         pEff1EtaP->Fill(IsTrueMatching, Eta, P);
      }

      if (-1.35 < Eta && Eta < 1.35) {
         pEff2P->Fill(IsTrueMatching, P);
         pEff2Eta->Fill(IsTrueMatching, Eta);
         pEff2EtaP->Fill(IsTrueMatching, Eta, P);
      }

      if (IsMatchingExist) {
         pCont1P->Fill(!IsTrueMatching, P);
         pCont1Eta->Fill(!IsTrueMatching, Eta);
         pCont1EtaP->Fill(!IsTrueMatching, Eta, P);
      }

      if (IsNmatch) {
         pCont2P->Fill(!IsTrueMatching, P);
         pCont2Eta->Fill(!IsTrueMatching, Eta);
         pCont2EtaP->Fill(!IsTrueMatching, Eta, P);
      }

   } // cycle by TPC KF tracks
}
