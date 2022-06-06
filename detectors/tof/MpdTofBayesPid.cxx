//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofBayesPid
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>

#include "FairLogger.h"
#include "MpdTofMatchingData.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTofHit.h"
#include "MpdMCTrack.h"

#include "MpdTofBayesPid.h"

//#define _PDFbyBeta 1

using namespace std;

ClassImp(MpdBayesPriors);
const char *MpdBayesPriors::fSpeciesNames[fNpdf] = {"Proton", "Pion", "Kaon", "Electron"};
const TH1D *MpdBayesPriors::fPriorOrigin = new TH1D("MpdBayesPriors_fPriorOrigin", ";P, GeV/c; events", 1000, 0., 10.);
//------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------	MpdBayesPriors
//----------
//------------------------------------------------------------------------------------------------------------------------
MpdBayesPriors::MpdBayesPriors(size_t iterNumb, const char *flnm, bool doLoad) : fIterNmb(iterNumb)
{
   fFlnm = _getFullFlnm(iterNumb, flnm);

   // load or create priors
   if (doLoad) {
      if (iterNumb == 1) {
         _createPriors();
         _setDefaultPriors();
         Write();
      } else {
         bool ok = _loadPriors();
         assert(ok);
      }
   } else
      _createPriors();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::_createPriors() // create TH1D prior histos
{
   cout << "\n [MpdBayesPriors:_createPriors] Create priors at file: " << fFlnm;

   for (size_t pid = 0; pid < fNpdf; pid++) {
      TString name = _getTHname(fIterNmb, pid);
      fPriors[pid] = dynamic_cast<TH1D *>(fPriorOrigin->Clone(name.Data()));

      _add(fPriors[pid]);
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::_setDefaultPriors() // set default value = 1. for all bins
{
   cout << "\n [MpdBayesPriors:_setDefaultPriors] Set priors to default at file: " << fFlnm;

   for (size_t pid = 0; pid < fNpdf; pid++) {
      for (size_t bin = 1, nBins = fPriorOrigin->GetXaxis()->GetNbins(); bin <= nBins; bin++)
         fPriors[pid]->SetBinContent(bin, 1.);
   }
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdBayesPriors::_loadPriors()
{
   cout << "\n [MpdBayesPriors:_loadPriors] Load priors from file: " << fFlnm;

   auto  ptr = gFile;
   TFile file(fFlnm.Data(), "READ");
   if (file.IsZombie()) {
      cerr << "\n [MpdBayesPriors:_loadPriors] error open file: " << fFlnm << endl;
      return false;
   }

   bool ok = true;
   for (size_t pid = 0; pid < fNpdf; pid++) {
      TString name = _getTHname(fIterNmb, pid);
      auto    h1   = (TH1D *)file.Get(name.Data());
      _add(fPriors[pid] = h1);
      if (h1 == nullptr) cerr << "\n [MpdBayesPriors:_loadPriors] error load prior: " << name << endl;

      ok &= (bool)h1;
   }

   file.Close();
   gFile = ptr;

   return ok;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::Yields2Priors(TH1D **data)
{
   /*	// Prior_n+1(H) = Y(H)/ Y(pion)
      for(size_t pid = 0; pid < fNpdf; pid++)
      {
         fPriors[pid]->Divide(data[pid], data[1]);
      }
   */
   // Prior_n+1(H) = Y(H)/ SUM(Y(pid))
   auto hsum = (TH1D *)data[0]->Clone("Yields2Priors_Sum");

   for (size_t pid = 1; pid < fNpdf; pid++) {
      hsum->Add(data[pid]);
   }

   for (size_t pid = 0; pid < fNpdf; pid++) {
      fPriors[pid]->Divide(data[pid], hsum);
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::Write()
{
   cout << "[MpdBayesPriors::Write] Update  " << fFlnm.Data() << " file. ";

   auto  ptr = gFile;
   TFile file(fFlnm.Data(), "RECREATE");
   fList.Write();
   file.Close();
   gFile = ptr;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::_add(TH1D *h1)
{
   h1->SetDirectory(nullptr);
   fList.Add(h1);
}
//------------------------------------------------------------------------------------------------------------------------
TString MpdBayesPriors::_getFullFlnm(size_t iterNmb, const char *flnm) const
{
   return TString::Format("%s_iter%d.root", flnm, (int)iterNmb);
}
//------------------------------------------------------------------------------------------------------------------------
TString MpdBayesPriors::_getTHname(size_t iterNmb, size_t pid) const
{
   return TString::Format("MpdBayesPriors_%s_iter%zu", fSpeciesNames[pid], iterNmb);
}
//------------------------------------------------------------------------------------------------------------------------
double MpdBayesPriors::GetPrior(size_t pid, double beta) const
{
   assert(pid < fNpdf);

   return fPriors[pid]->Interpolate(beta);
}
//------------------------------------------------------------------------------------------------------------------------
const char *MpdBayesPriors::GetSpeciesName(size_t index)
{
   assert(index < fNpdf);

   return fSpeciesNames[index];
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::FillPrior(size_t pid, double beta, double weight)
{
   assert(pid < fNpdf);

   fPriors[pid]->Fill(beta, weight);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdBayesPriors::DumpPriors(const char *comment, std::ostream &os) const
{
   auto prec = os.precision(5);
   os << "\n [MpdBayesPriors::DumpPriors]-------------------------------------------------------------------->>> ";

   if (comment != nullptr) os << comment;

   for (size_t pid = 0; pid < fNpdf; pid++)
      os << "\n Integral=" << fPriors[pid]->Integral() << "  H=(" << fSpeciesNames[pid] << ")";

   os << "\n [MpdBayesPriors::DumpPriors]--------------------------------------------------------------------<<< ";
   os.precision(prec);
}
//------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------	PidProbabilitiesMatrix	--
//------------------------------------------------------------------------------------------------------------------------
PidProbabilitiesMatrix::PidProbabilitiesMatrix(double thresh, const char *prefix)
{
   SetThresh(thresh);

   if (prefix) fFlnm = prefix;

   // create histos
   for (size_t pid = 0; pid < fNpdf; pid++) {
      fA_true[pid] = CreatePriorClone(TString::Format("TofPidMatrix_fA_true_%s", MpdBayesPriors::GetSpeciesName(pid)));
      fA_meas[pid] = CreatePriorClone(TString::Format("TofPidMatrix_fA_meas_%s", MpdBayesPriors::GetSpeciesName(pid)));

      for (size_t index = 0; index < fNpdf; index++) {
         fProb[index][pid] = CreatePriorClone(TString::Format("TofPidMatrix_fProb_%zd_%zd", index, pid));
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
PidProbabilitiesMatrix::~PidProbabilitiesMatrix()
{
   fList.Delete();
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::Update(double *prob_H_S, double Pt, Int_t pdgcode)
{
   auto index = MpdTofBayesPid::GetIndex(pdgcode);
   if (index >= 0) {
      fA_true[index]->Fill(Pt, 1.);

      if (fPidThres < 0.) //  max. probability method(1)
      {
         double maxProb = 0.;
         size_t maxPid  = -1;
         for (size_t pid = 0; pid < fNpdf; pid++) {
            double prob = prob_H_S[pid]; //->Interpolate(beta);

            if (prob > maxProb) {
               maxPid  = pid;
               maxProb = prob;
            }
         }

         if (maxProb > 0.) // make id decision
         {
            fProb[index][maxPid]->Fill(Pt, 1.);
            fA_meas[maxPid]->Fill(Pt, 1.);
         }
      } else // fixed threshold method(2)
      {
         for (size_t pid = 0; pid < fNpdf; pid++) {
            double prob = prob_H_S[pid]; //->Interpolate(beta);

            if (prob > fPidThres) // make id decision
            {
               fProb[index][pid]->Fill(Pt, 1.);
               fA_meas[pid]->Fill(Pt, 1.);
            }
         }
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::FinalizeMatrices()
{
   // create&fill Efficiency and Contamination histos
   for (size_t pid = 0; pid < fNpdf; pid++) {
      // Efficiency_ii = N_i identified as i (εPID) / A_true_i
      TH1D *h1 = CreatePriorClone(TString::Format("TofPidMatrix_Efficiency_%s", MpdBayesPriors::GetSpeciesName(pid)));
      h1->Divide(fProb[pid][pid], fA_true[pid]);

      // The contamination of the species j due to a different species i (c ji ) is the number of particles belonging
      // to species i that are wrongly identified as j (Ni identified as j ), divided by the total number of identified
      // j particles (Ameas j) Contamination_ji = N_i identified as j  / A_meas_j

      for (size_t pid2 = 0; pid2 < fNpdf; pid2++) {
         if (pid == pid2) continue;
         h1 =
            CreatePriorClone(TString::Format("TofPidMatrix_Contamination_%s_by_%s", MpdBayesPriors::GetSpeciesName(pid),
                                             MpdBayesPriors::GetSpeciesName(pid2)));
         h1->Divide(fProb[pid2][pid], fA_meas[pid]);
      }
   }

   // Normalize probability matrices
   /*	for(size_t index = 0; index < fNpdf; index++)
      {
         //  true abundances of species i.
         double A_true_Integral = fA_true[index]->Integral();

         if(A_true_Integral > 0.)
         {
            for(size_t pid = 0; pid < fNpdf; pid++)
            {
               fProb[index][pid]->Divide(fA_true[index]);
            }
         }
      }
   */
}
//------------------------------------------------------------------------------------------------------------------------
TH1D *PidProbabilitiesMatrix::CreatePriorClone(TString name)
{
   TH1D *h1 = new TH1D(name.Data(), "; P, GeV/c", 1000, 0., 10.);
   _add(h1);

   return h1;
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::Print(const char *comment, std::ostream &os)
{
   auto prec = os.precision(5);

   if (comment != nullptr) os << comment;

   os << std::scientific << "\n\n\t" << fName;
   os << "\n\t\tProton\t\tPion\t\tKaon\t\tElectron";

   _printLine(fA_true, "\n A_true=\t");
   _printLine(fA_meas, "\n A_meas=\t");

   os << "\n A_calc_meas=\t(";
   for (size_t pid = 0; pid < fNpdf; pid++) os << Calc_A_meas(pid) << ",\t";
   os << ")";

   os << "\n\n EpsPID matrix:";
   _printLine(fProb[0], "\n Proton=  \t");
   _printLine(fProb[1], "\n Pion=    \t");
   _printLine(fProb[2], "\n Kaon=    \t");
   _printLine(fProb[3], "\n Electron=\t");

   os.precision(prec);
}
//------------------------------------------------------------------------------------------------------------------------
double PidProbabilitiesMatrix::Calc_A_meas(size_t index) const
{
   // A_meas =  εPID  * A_true (bin by bin)
   double sum = 0.;
   for (size_t pid = 0; pid < fNpdf; pid++) {
      for (Int_t bin = 1, lastBin = MpdBayesPriors::GetPriorXbins(); bin <= lastBin; bin++) {
         sum += fProb[pid][index]->GetBinContent(bin) * fA_true[pid]->GetBinContent(bin);
      }
   }

   return sum;
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::_printLine(TH1D **data, const char *comment, std::ostream &os)
{
   if (comment != nullptr) os << comment;

   os << "(";
   for (size_t pid = 0; pid < fNpdf; pid++) os << data[pid]->Integral() << ",\t";
   os << ")";
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::SetThresh(double value)
{
   fPidThres = value; // negative or [0., 1.]

   // Set default name
   if (fPidThres < 0.)
      fName = "max. probability method.";
   else
      fName.Form("fixed threshold(%f) method.", fPidThres);
}
//------------------------------------------------------------------------------------------------------------------------
void PidProbabilitiesMatrix::_add(TH1D *h1)
{
   h1->SetDirectory(nullptr);
   fList.Add(h1);
}
//------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------	MpdTofBayesPid
//----------
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofBayesPid);
//------------------------------------------------------------------------------------------------------------------------
MpdTofBayesPid::MpdTofBayesPid(const char *flnm, size_t iterNmb, const char *priorPrefix, const char *nextPriorPrefix,
                               const char *QAflnm, const char *taskName, Int_t verbose, Bool_t useMCdata)
   : FairTask(taskName, verbose), fUseMCData(useMCdata), fDoPDFs(false), fPDFflnm(flnm), fQAflnm(QAflnm)
{
   if (priorPrefix) fPriorPrefix = priorPrefix;
   TString nextPrefix = (nextPriorPrefix != nullptr) ? TString(nextPriorPrefix) : fPriorPrefix;

   fDoTest = (QAflnm != nullptr);

   if (fDoTest) {
      // create QA histos
      for (size_t pid = 0; pid < MpdBayesPriors::GetNpdf(); pid++) {
         TString name = TString::Format("Tof_CalcProb_H_S_%s", MpdBayesPriors::GetSpeciesName(pid));
         Add(fQAlist, fH_S[pid] = new TH2D(name.Data(), "", 1000, 0., 1.1, 1000, 0., 1.1));
      }
   }

   //  create TH1D fProb_H_S and fProb_S_H histos
   for (size_t pid = 0; pid < MpdBayesPriors::GetNpdf(); pid++) {
      TString name = TString::Format("Tof_Prob_Yield_%s", MpdBayesPriors::GetSpeciesName(pid));
      fYield[pid]  = dynamic_cast<TH1D *>(MpdBayesPriors::GetPriorOrigin()->Clone(name.Data()));
   }

   // load current priors
   fPriors = new MpdBayesPriors(iterNmb, fPriorPrefix.Data(), true);

   // create next iteration priors containers
   fNextPrior = new MpdBayesPriors(iterNmb + 1, nextPrefix.Data(), false);

   // load PDFs (fDoPDFs = false)
   bool res = CreateLoadPDFs();

   // create PidProbabilitiesMatrix with thresholds
   for (size_t i = 0; i < fNmatrix; i++) {
      double thresh = (0 == i) ? -1. : i * 0.1; // default threshoulds: -1., 0.1, 0.2, ... ,0.9 (Use
                                                // MpdTofBayesPid::SetProbMatrixThresh() to change.)

      fEffMatrix[i] = new PidProbabilitiesMatrix(thresh, nextPriorPrefix);
   }

   assert(res);
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofBayesPid::MpdTofBayesPid(const char *flnm)
   : FairTask("TOF_pid_createPDFs", 1), fUseMCData(true), fDoTest(false), fDoPDFs(true), fPDFflnm(flnm)
{
   // create PDFs (fDoPDFs = true)
   CreateLoadPDFs();
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofBayesPid::~MpdTofBayesPid()
{
   // delete histos
   fQAlist.Delete();
   fPDFlist.Delete();

   delete fPriors;
   delete fNextPrior;

   for (size_t i = 0; i < fNmatrix; i++) delete fEffMatrix[i];
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdTofBayesPid::Init()
{
   aTPCkfTracks  = (TClonesArray *)FairRootManager::Instance()->GetObject("TpcKalmanTrack");
   aTofMatchings = (TClonesArray *)FairRootManager::Instance()->GetObject("TOFMatching");
   aTofHits      = (TClonesArray *)FairRootManager::Instance()->GetObject("TOFHit");
   assert(aTofMatchings);
   assert(aTPCkfTracks);
   assert(aTofHits);

   if (fUseMCData) {
      aMcTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
      assert(aMcTracks);
   }

   LOG(info) << "[MpdTofBayesPid::Init] Initialization finished succesfully.";

   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::Exec(Option_t *option)
{
   if (fDoPDFs) {
      FillPDFs();
      return;
   }

   size_t TofMatchingNmb = aTofMatchings->GetEntriesFast();
   for (size_t i = 0; i < TofMatchingNmb; i++) // cycle by Tof matchings
   {
      auto         match = (MpdTofMatchingData *)aTofMatchings->At(i);
      const double mass2 = match->GetMass2(), P = match->GetMomentum().Mag(), beta = match->GetBeta();

      auto tofHit  = (MpdTofHit *)aTofHits->UncheckedAt(match->GetTofHitIndex());
      auto kfTrack = (MpdTpcKalmanTrack *)aTPCkfTracks->UncheckedAt(match->GetKFTrackIndex());
      assert(tofHit);
      assert(kfTrack);

      // cuts:
      if (kfTrack->GetNofHits() < 20) continue; // pass Kalman tracks with hits number >  20
      if (fUseMCData) {
         auto mcTrack = (MpdMCTrack *)aMcTracks->UncheckedAt(kfTrack->GetTrackID());
         if (mcTrack->GetMotherId() != -1) continue; // pass primary tracks only
      }

      Int_t rank = 0;
#ifdef _PDFbyBeta
      if (GetProb_S_H(mass2, beta, kfTrack->GetDedx(), rank, fProb_S_H)) // notZeroExist = true
#else
      if (GetProb_S_H(mass2, P, kfTrack->GetDedx(), rank, fProb_S_H)) // notZeroExist = true
#endif
      {
         // calc. P(H_pid|S)
         CalcProb_H_S(P);

         // fill fProb_H_S to MpdTofMatchingData
         std::copy(begin(fProb_H_S), end(fProb_H_S), begin(match->fBayesPid));

         // update yields
         for (size_t pid = 0; pid < fNpdf; pid++) {
            /////////////////				if(fProb_H_S[pid] > 0.) fYield[pid]->Fill(P, fProb_H_S[pid]);
            if (fProb_H_S[pid] > 0.) fNextPrior->FillPrior(pid, P, fProb_H_S[pid]); // update Priors
         }

         if (fUseMCData) {
            auto mcTrack = (MpdMCTrack *)aMcTracks->UncheckedAt(kfTrack->GetTrackID());
            assert(mcTrack);
            Int_t pdg = mcTrack->GetPdgCode();

            // update all PidEfficiencyMatrix
            for (size_t k = 0; k < fNmatrix; k++) {
               fEffMatrix[k]->Update(fProb_H_S, P, pdg);
            }

            if (fDoTest) {
               auto index = MpdTofBayesPid::GetIndex(pdg);
               if (index >= 0) // accepted true pid index
               {
                  for (size_t pid = 0; pid < fNpdf; pid++) {
                     if (pid != index) fH_S[index]->Fill(fProb_H_S[index], fProb_H_S[pid]);
                  }
               }
            }
         }
         // DumpProbs(fProb_H_S, "fProb_H_S", pdg);
      }
   } // cycle by Tof matchings
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::Finish()
{
   auto WriteTList = [](const TString &flnm, const TList &list) {
      LOG(debug2) << "[MpdTofBayesPid::Finish] Update  " << flnm.Data() << " file. ";
      auto  ptr = gFile;
      TFile file(flnm.Data(), "RECREATE");
      list.Write();
      file.Close();
      gFile = ptr;
   };

   if (fDoPDFs) {
      WriteTList(fPDFflnm, fPDFlist);
   } else {
      if (fDoTest) WriteTList(fQAflnm, fQAlist);

      // calc. new iteration priors from yields and write to file
      if (fNextPrior) {
         /////////////// update ALREADY	PRIOR	fNextPrior->Yields2Priors(fYield);
         fNextPrior->Write();
      }

      if (fUseMCData) {
         for (size_t k = 0; k < fNmatrix; k++) {
            // normalize PidProbabilitiesMatrix histos and calc. efficiency and contamination histos
            fEffMatrix[k]->FinalizeMatrices();
            fEffMatrix[k]->Print();

            // write PidProbabilitiesMatrix histos
            TString flnm;
            flnm.Form("PidMatrices_%s_%zd.root", fEffMatrix[k]->GetFlnm(), k);
            WriteTList(flnm, fEffMatrix[k]->GetHistos());
         }
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::CalcProb_H_S(double P)
{
   // calc. denominator(5) = Sum_by_k(P(S|H_k) * C(H_k))
   double norm = 0.;
   for (size_t pid = 0; pid < fNpdf; pid++) {
      norm += fProb_S_H[pid] * fPriors->GetPrior(pid, P);
   }

   // calc. fProb_H_S:  P(H_k|S) = P(S|H_k) * C(H_k) / Sum_by_k(P(S|H_k) * C(H_k))
   for (size_t pid = 0; pid < fNpdf; pid++) {
      fProb_H_S[pid] = (norm > 0.) ? fProb_S_H[pid] * fPriors->GetPrior(pid, P) / norm : 0.;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::DumpProbs(double *data, const char *comment, Int_t pdg, std::ostream &os) const
{
   auto prec = os.precision(5);
   os << "\n [MpdTofBayesPid::DumpProbs]-------------------------------------------------------------------->>> ";

   if (comment != nullptr) os << comment;
   if (pdg != 0) os << " pdg = " << pdg;

   for (size_t pid = 0; pid < fNpdf; pid++)
      os << "\n Prob=" << data[pid] << "  H=(" << MpdBayesPriors::GetSpeciesName(pid) << ")";

   os << "\n [MpdTofBayesPid::DumpProbs]--------------------------------------------------------------------<<< ";
   os.precision(prec);
}
//------------------------------------------------------------------------------------------------------------------------
// increase square rank until the all prob. sums becomes not zero or will reach the limit; return - notZeroExist and
// rank
bool MpdTofBayesPid::GetProb_S_H(double mass2, double P, double dedx, Int_t &rank, double *data)
{
   rank              = 0;
   bool notZeroExist = false;
   do {
      if (10 == rank) break; // square 20x20x20 bins limit

      notZeroExist = _getProbs(rank, mass2, P, dedx, data);
      rank++;
   } while (!notZeroExist);

   return notZeroExist;
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofBayesPid::_getProbs(Int_t rank, double mass2, double P, double dedx, double *data)
{
   static const auto  Xaxis = hPdf[0]->GetXaxis(), Yaxis = hPdf[0]->GetYaxis(), Zaxis = hPdf[0]->GetZaxis();
   static const Int_t XNbins = Xaxis->GetNbins(), YNbins = Yaxis->GetNbins(), ZNbins = Zaxis->GetNbins();

   Int_t Xbin = Xaxis->FindBin(mass2), Ybin = Yaxis->FindBin(P), Zbin = Zaxis->FindBin(dedx);
   Int_t Xfirst = (Xbin - rank >= 1) ? Xbin - rank : 1,
         Xlast  = (Xbin + rank <= XNbins) ? Xbin + rank : XNbins; // clamp(Xbin+-rank, 1, XNbins);
   Int_t Yfirst = (Ybin - rank >= 1) ? Ybin - rank : 1,
         Ylast  = (Ybin + rank <= YNbins) ? Ybin + rank : YNbins; // clamp(Ybin+-rank, 1, YNbins);
   Int_t Zfirst = (Zbin - rank >= 1) ? Zbin - rank : 1,
         Zlast  = (Zbin + rank <= ZNbins) ? Zbin + rank : ZNbins; // clamp(Zbin+-rank, 1, ZNbins);

   // sum. of probability inside square [Xbin-rank, Xbin+rank] x [Ybin-rank, Ybin+rank] x [Zbin-rank, Zbin+rank]
   bool notZeroExist = false;
   for (size_t pid = 0; pid < fNpdf; pid++) {
      data[pid] = 0.;

      for (Int_t ix = Xfirst; ix <= Xlast; ix++) // cycle by square bins
      {
         for (Int_t iy = Yfirst; iy <= Ylast; iy++) {
            for (Int_t iz = Zfirst; iz <= Zlast; iz++) {
               data[pid] += hPdf[pid]->GetBinContent(ix, iy, iz);
            }
         }
      }

      assert(fPdfIntegrals[pid] > 0.);

      data[pid] /= fPdfIntegrals[pid];

      notZeroExist |= (data[pid] > 0.);
   }

   // DumpProbs(data, TString::Format("\n VVVVVVVVVVVVVVVV  rank=%d, notZeroExist=%d", rank, notZeroExist).Data(), 0);

   return notZeroExist;
}
//------------------------------------------------------------------------------------------------------------------------
Int_t MpdTofBayesPid::GetIndex(Int_t pdg)
{
   if (2212 == pdg || -2212 == pdg)
      return 0;
   else if (211 == pdg || -211 == pdg)
      return 1;
   else if (321 == pdg || -321 == pdg)
      return 2;
   else if (11 == pdg || -11 == pdg)
      return 3;

   return -1;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::FillPDFs()
{
   size_t TofMatchingNmb = aTofMatchings->GetEntriesFast();
   for (size_t i = 0; i < TofMatchingNmb; i++) // cycle by Tof matchings
   {
      auto match = (MpdTofMatchingData *)aTofMatchings->At(i);

      auto tofHit  = (MpdTofHit *)aTofHits->UncheckedAt(match->GetTofHitIndex());
      auto kfTrack = (MpdTpcKalmanTrack *)aTPCkfTracks->UncheckedAt(match->GetKFTrackIndex());
      assert(tofHit);
      assert(kfTrack);

      if (MpdTofUtils::HaveTail & tofHit->GetFlag()) continue; // pass only single mcPoint hits

      Int_t tid = kfTrack->GetTrackID();
      if (!tofHit->CheckTid(tid)) continue; // pass only true matching

      auto mcTrack = (MpdMCTrack *)aMcTracks->UncheckedAt(tid);
      assert(mcTrack);

      auto index = GetIndex(mcTrack->GetPdgCode());
#ifdef _PDFbyBeta
      if (index >= 0) hPdf[index]->Fill(match->GetMass2(), match->GetBeta(), kfTrack->GetDedx());
#else
      if (index >= 0) hPdf[index]->Fill(match->GetMass2(), match->GetMomentum().Mag(), kfTrack->GetDedx());
#endif
   } // cycle by Tof matchings
}
//--------------------------------------------------------------------------------------------------------------------------------------
bool MpdTofBayesPid::CreateLoadPDFs()
{
   static const TString pdfPrefix("MpdTOFpid_Pdf_");

   if (fDoPDFs) // create new PDFs
   {
      cout << "\n [MpdTofBayesPid:CreateLoadPDFs] Create new PDFs histos: ";

      for (size_t pid = 0; pid < fNpdf; pid++) {
         TString name = pdfPrefix;
         name += MpdBayesPriors::GetSpeciesName(pid);
         cout << "\n " << name;

//			hPdf[pid] = new TH3D(name.Data(), ";mass^{2}, GeV^{2};P, GeV/c;dE/dx", 500, -0.5, 3., 500, 0., 10., 100,
// 0., 10.);
#ifdef _PDFbyBeta
         hPdf[pid] = new TH3D(name.Data(), ";mass^{2}, GeV^{2};#beta;dE/dx", 500, -1.5, 6., 500, 0., 1.1, 100, 0., 10.);
#else
         hPdf[pid] =
            new TH3D(name.Data(), ";mass^{2}, GeV^{2};P, GeV/c;dE/dx", 500, -1.5, 6., 500, 0., 10., 100, 0., 10.);
#endif
         Add(fPDFlist, hPdf[pid]);
      }

      return true;
   } else // load PDFs from file
   {
      cout << "\n [MpdTofBayesPid:CreateLoadPDFs] Load PDF histos from file: " << fPDFflnm;

      auto  ptr = gFile;
      TFile file(fPDFflnm.Data(), "READ");
      if (file.IsZombie()) {
         cerr << "\n [MpdTOFpid:CreateLoadPDFs] error open file: " << fPDFflnm << endl;
         return false;
      }

      bool retval = true;

      for (size_t n = 0; n < fNpdf; n++) {
         TString name = pdfPrefix;
         name += MpdBayesPriors::GetSpeciesName(n);
         auto h3 = (TH3D *)file.Get(name.Data());
         if (h3) {
            fPdfIntegrals[n] = h3->Integral();
            assert(fPdfIntegrals[n] > 0.);

            Add(fPDFlist, hPdf[n] = h3);
            retval &= true;
         } else {
            cerr << "\n [MpdTofBayesPid:CreateLoadPDFs] error load PDF name: " << name << endl;
            retval = false;
         }
      }

      file.Close();
      gFile = ptr;

      return retval;
   }
}
//--------------------------------------------------------------------------------------------------------------------------------------
template <typename T>
void MpdTofBayesPid::Add(TList &list, T *obj)
{
   assert(obj->InheritsFrom("TH1") || !std::string(obj->ClassName()).compare("TEfficiency"));

   obj->SetDirectory(nullptr);
   list.Add(obj);
}
//--------------------------------------------------------------------------------------------------------------------------------------
void MpdTofBayesPid::SetProbMatrixThresh(size_t index, double value)
{
   assert(index < fNmatrix);
   assert(fEffMatrix[index]);

   fEffMatrix[index]->SetThresh(value);
}
//------------------------------------------------------------------------------------------------------------------------
