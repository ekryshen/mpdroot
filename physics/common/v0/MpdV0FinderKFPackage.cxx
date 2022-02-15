/*
 * MpdV0FinderKFPackage.cxx
 *
 *  Created on: 19 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0FinderKFPackage.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TObjArray.h>
#include <utility>

#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFParticleDef.h"
#include "KFParticleTopoReconstructor.h"
#include "KFPTrackVector.h"
#include "KFPVertex.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniTrackCovMatrix.h"
#include "MpdV0CandidateCut.h"
#include "MpdV0DaughterCut.h"
#include "MpdV0Particle.h"
#include "MpdV0Track.h"
#include "KFPTrackVector.h"
#include "KFPTrack.h"

MpdV0FinderKFPackage::MpdV0FinderKFPackage(Int_t pidMom, Int_t pidFirstDau, Int_t pidSecDau)
   : MpdV0Finder(pidMom, pidFirstDau, pidSecDau), fKFFinder(nullptr), fMiniCovMatrix(nullptr)
{
}

void MpdV0FinderKFPackage::ExecDst(Option_t *option) {}

void MpdV0FinderKFPackage::ExecMiniDst(Option_t *option)
{
   fKFFinder->Clear();
   KFPVertex kfVertex;
   kfVertex.SetXYZ(fEventVertex.X(), fEventVertex.Y(), fEventVertex.Z());
   kfVertex.SetCovarianceMatrix(0, 0, 0, 0, 0, 0);
   kfVertex.SetChi2(-100);
   kfVertex.SetNContributors(fMiniTracks->GetEntriesFast());
   std::vector<int> pvTrackIds;
   KFVertex         pv(kfVertex);
   fKFFinder->AddPV(pv, pvTrackIds);

   // fKFFinder->AddPV(kfVertex);
   std::pair<TObject *, Int_t> data;
   KFPTrackVector              tracksOut;
   int                         bufferedTrack = 0;
   Double_t                    params[6];

   std::vector<int> index1, index2;
   for (int iTrack = 0; iTrack < fMiniTracks->GetEntriesFast(); iTrack++) {
      MpdMiniTrack *track = (MpdMiniTrack *)fMiniTracks->UncheckedAt(iTrack);
      if (fPositiveDaughterCut->PassMiniDstTrack(*track)) {
         index1.push_back(iTrack);
      }
      if (fNegativeDaughterCut->PassMiniDstTrack(*track)) {
         index2.push_back(iTrack);
      }
   }
   fInputTracks.Resize(index1.size() + index2.size());
   for (unsigned int iTrack = 0; iTrack < index1.size(); iTrack++) {
      MpdMiniTrack *track = (MpdMiniTrack *)fMiniTracks->UncheckedAt(index1[iTrack]);
      fInputTracks.SetPDG(fPidDauPos, bufferedTrack);
      fInputTracks.SetQ(track->charge(), bufferedTrack);
      fInputTracks.SetId(index1[iTrack], bufferedTrack);
      if (track->isPrimary()) {
         fInputTracks.SetPVIndex(0, bufferedTrack);
      } else {
         fInputTracks.SetPVIndex(-1, bufferedTrack);
      }

      std::vector<float> covMat = GetCovMatrixMini(index1[iTrack], params);
      for (int i = 0; i < 6; i++) fInputTracks.SetParameter(params[i], i, bufferedTrack);
      for (int i = 0; i < 26; i++) fInputTracks.SetCovariance(covMat[i], i, bufferedTrack);

      bufferedTrack++;
   }

   for (unsigned int iTrack = 0; iTrack < index2.size(); iTrack++) {
      MpdMiniTrack *track = (MpdMiniTrack *)fMiniTracks->UncheckedAt(index2[iTrack]);
      fInputTracks.SetPDG(fPidDauNeg, bufferedTrack);
      fInputTracks.SetQ(track->charge(), bufferedTrack);
      fInputTracks.SetId(index2[iTrack], bufferedTrack);
      if (track->isPrimary()) {
         fInputTracks.SetPVIndex(0, bufferedTrack);
      } else {
         fInputTracks.SetPVIndex(-1, bufferedTrack);
      }

      std::vector<float> covMat = GetCovMatrixMini(index2[iTrack], params);
      for (int i = 0; i < 6; i++) fInputTracks.SetParameter(params[i], i, bufferedTrack);
      for (int i = 0; i < 26; i++) fInputTracks.SetCovariance(covMat[i], i, bufferedTrack);

      bufferedTrack++;
   }
   fKFFinder->Init(fInputTracks, tracksOut);
   fKFFinder->SortTracks();
   fKFFinder->ReconstructParticles();
   WriteCandidates();
}

InitStatus MpdV0FinderKFPackage::Init()
{
   fKFFinder = new KFParticleTopoReconstructor();
   fCandicateCut->SetupKF(fKFFinder);
   if (FairRootManager::Instance()->GetObject("Event") != nullptr) { // minist file
      fMiniCovMatrix = (TClonesArray *)FairRootManager::Instance()->GetObject("TrackCovMatrix");
   }
   return MpdV0Finder ::Init();
}
std::vector<float> MpdV0FinderKFPackage::GetCovMatrixMini(Int_t index, Double_t *newParams)
{
   MpdMiniTrack *         mini    = (MpdMiniTrack *)fMiniTracks->UncheckedAt(index);
   MpdMiniTrackCovMatrix *cMatrix = (MpdMiniTrackCovMatrix *)fMiniCovMatrix->UncheckedAt(index);
   TVector3               vtx;
   MpdTpcKalmanTrack      kalmanIn(0, vtx), kalmanOut;
   kalmanIn.SetPos(cMatrix->R0());
   kalmanIn.SetPosNew(cMatrix->R());
   kalmanIn.SetNofHits(mini->nHits());
   TMatrixD stateVector = cMatrix->stateVector();
   kalmanIn.SetParam(stateVector);

   TMatrixDSym covMatrix = cMatrix->covarianceMatrix();
   kalmanIn.SetCovariance(covMatrix);
   KalmanToInnerTpc(kalmanIn, kalmanOut);
   Double_t phi = kalmanOut.GetParam(0) / kalmanOut.GetPosNew();
   newParams[0] = kalmanOut.GetPosNew() * TMath::Cos(phi);
   newParams[1] = kalmanOut.GetPosNew() * TMath::Sin(phi);
   newParams[2] = kalmanOut.GetZ();

   newParams[3] = kalmanOut.Momentum3().X();
   newParams[4] = kalmanOut.Momentum3().Y();
   newParams[5] = kalmanOut.Momentum3().Z();
   std::vector<float> cov21 =
      DoErrorPropagationToXYZPxPyPz(kalmanOut.GetCovariance(), kalmanOut.GetParam(), &kalmanOut);
   return cov21;
}

void MpdV0FinderKFPackage::KalmanToInnerTpc(MpdTpcKalmanTrack &in, MpdTpcKalmanTrack &out)
{
   MpdTpcKalmanTrack tmp(in);
   tmp.SetDirection(MpdKalmanTrack::kOutward);

   /// Propagate TPC track to 27 cm inward (fPos == 27)
   MpdKalmanHit hitTmp;
   hitTmp.SetType(MpdKalmanHit::kFixedR);
   hitTmp.SetPos(tmp.GetPos()); // fPos ~= 27 cm

   tmp.SetParamNew(*tmp.GetParam());
   tmp.SetPos(tmp.GetPosNew());
   tmp.ReSetWeight();
   TMatrixDSym w = *tmp.GetWeight(); // save current weight matrix

   Bool_t ok = MpdKalmanFilter::Instance()->PropagateToHit(&tmp, &hitTmp, kFALSE, kFALSE);

   if (!ok) {
      out = tmp;
      return;
   }
   TMatrixDSym cov(*tmp.Weight2Cov());
   TMatrixD    param(*tmp.GetParamNew());

   tmp.SetParam(param);
   tmp.SetCovariance(cov);

   out = tmp;
}

std::vector<float> MpdV0FinderKFPackage::DoErrorPropagationToXYZPxPyPz(TMatrixDSym *cov, TMatrixD *params,
                                                                       MpdTpcKalmanTrack *track)
{
   // State track vector in MPD:
   // 0) r\Phi_{T}
   // 1) z
   // 2) \Phi
   // 3) \Lambda (dip angle)
   // 4) -q / Pt

   // J : A --> B
   // A : x, y, z, px, py, pz
   // B : x, y, z, \Phi, \Lambda, -q / Pt

   // Getting Jacobian ...
   Float_t J[6][6];
   for (Int_t i = 0; i < 6; i++)
      for (Int_t j = 0; j < 6; j++) J[i][j] = 0;

   // Getting corresponding track parameters ...
   Int_t q  = track->Charge();
   Int_t pt = track->Pt();

   Double_t Phi    = (*params)(2, 0);
   Double_t Lambda = (*params)(3, 0);

   // Coordinate transformations ...
   J[0][0] = 1.; // dx / dx
   J[1][1] = 1.; // dy / dy
   J[2][2] = 1.; // dz / dz

   // Momentum transformations ...
   J[3][3] = -pt * TMath::Sin(Phi);         // dPx / d\Phi
   J[3][5] = pt * pt * TMath::Cos(Phi) / q; // dPx / d(-q / Pt)

   J[4][3] = pt * TMath::Cos(Phi);          // dPy / d\Phi
   J[4][5] = pt * pt * TMath::Sin(Phi) / q; // dPy / d(-q / Pt)

   J[5][4] = pt / (TMath::Cos(Lambda) * TMath::Cos(Lambda)); // dPz / dLambda
   J[5][5] = TMath::Tan(Lambda) * pt * pt / q;               // dPz / d(-q / pt)

   // Extending track covariance matrix by one row and one column ...
   Float_t CovIn[6][6]; // triangular -> symmetric matrix

   CovIn[0][0] = 1e-4; // dx. start from nowhere

   for (Int_t i = 1; i < 6; i++) {
      CovIn[i][0] = 0;
      CovIn[0][i] = 0;
   }

   for (Int_t i = 1; i < 6; i++) {
      for (Int_t j = 1; j <= i; j++) {
         CovIn[i][j] = (*cov)(i - 1, j - 1);
         CovIn[j][i] = (*cov)(i - 1, j - 1);
      }
   }

   Float_t CovInJt[6][6]; // CovInJt = CovIn * J^t
   for (Int_t i = 0; i < 6; i++)
      for (Int_t j = 0; j < 6; j++) {
         CovInJt[i][j] = 0;
         for (Int_t k = 0; k < 6; k++) CovInJt[i][j] += CovIn[i][k] * J[j][k];
      }

   Float_t CovOut[6][6]; // CovOut = J * CovInJt
   for (Int_t i = 0; i < 6; i++)
      for (Int_t j = 0; j < 6; j++) {
         CovOut[i][j] = 0;
         for (Int_t k = 0; k < 6; k++) CovOut[i][j] += J[i][k] * CovInJt[k][j];
      }

   std::vector<float> KFPCov; // symmetric matrix -> triangular

   for (Int_t i = 0; i < 6; i++)
      for (Int_t j = 0; j <= i; j++) KFPCov.push_back(CovOut[i][j]);

   return KFPCov;
}

void MpdV0FinderKFPackage::WriteCandidates()
{
   //   kfvector_int          ids    = tracks->Id();

   for (const auto &particle : fKFFinder->GetParticles()) {
      MpdV0Track       candidate;
      std::vector<int> dau = particle.DaughterIds();
      KFPTrack         firstDau, secondDau;
      fInputTracks.GetTrack(firstDau, dau[0]);
      fInputTracks.GetTrack(secondDau, dau[1]);

      int id1 = firstDau.Id();
      int id2 = secondDau.Id();
      candidate.SetPositiveDaughterIndex(id1);
      candidate.SetNegativeDaughterIndex(id2);
      candidate.SetMomPositiveDaughter(TVector3(firstDau.GetPx(), firstDau.GetPy(), firstDau.GetPz()));
      candidate.SetMomPositiveDaughter(TVector3(secondDau.GetPx(), secondDau.GetPy(), secondDau.GetPz()));
      candidate.SetChi2(particle.GetChi2());
      candidate.Recalculate(fEventVertex);
      candidate.SetMomentum(TVector3(particle.Px(), particle.Py(), particle.Pz()));
      //   if (fCandicateCut->Pass(candidate)) {
      MpdV0Track *v0 = (MpdV0Track *)fV0s->ConstructedAt(fV0s->GetEntriesFast());
      *v0            = candidate;
      //  }
   }
   LOG(WARNING) << "found" << fV0s->GetEntries();
}

MpdV0FinderKFPackage::~MpdV0FinderKFPackage()
{
   //  if (fKFFinder) delete fKFFinder;
}
