/// \class MpdKalmanFilter
///
/// Kalman filter track reconstructor in the MPD detector
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKalmanFilter.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdCodeTimer.h"
#include "FairField.h"
//#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairTask.h"
#include "MpdConstField.h"
#include "MpdMultiField.h"
#include <TMath.h>
#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

//#include "mkl.h"
//#include "mkl_lapacke.h"

#ifdef _OPENMP
#include "omp.h"
omp_lock_t geoManagerLock;
#endif

#include <iostream>
using std::cout;
using std::endl;

const Int_t      MpdKalmanFilter::fgkTriesMax = 1000;
const Double_t   MpdKalmanFilter::fgkEpsilon  = 0.001;
MpdKalmanFilter *MpdKalmanFilter::fgKF        = 0x0;

//__________________________________________________________________________
MpdKalmanFilter::MpdKalmanFilter(const char *name, const char *title)
   : FairTask("Kalman Filter engine"), fGeoScheme(new MpdKalmanGeoScheme), fNumer(1)
{
   /// Constructor
   fgKF = this;
}

//__________________________________________________________________________
MpdKalmanFilter *MpdKalmanFilter::Instance()
{
   /// Get pointer to the Kalman filter reconstructor singleton object
   if (!fgKF) {
      fgKF = new MpdKalmanFilter;
#ifdef _OPENMP
      omp_init_lock(&geoManagerLock);
#endif

      // automatic destroy is supposed
      std::atexit(DestroyInstance);
   }
   return fgKF;
}

//__________________________________________________________________________
MpdKalmanFilter *MpdKalmanFilter::Instance(const char *name, const char *title)
{
   /// Get pointer to the Kalman filter reconstructor singleton object
   if (!fgKF) {
      fgKF = new MpdKalmanFilter(name, title);
#ifdef _OPENMP
      omp_init_lock(&geoManagerLock);
#endif

      // automatic destroy is supposed
      std::atexit(DestroyInstance);
   }
   return fgKF;
}

//__________________________________________________________________________
MpdKalmanFilter::~MpdKalmanFilter()
{
   /// Destructor
   // FairRootManager *manager = FairRootManager::Instance();
   // manager->Write();
   delete fGeoScheme;
   fgKF = NULL;
}

//__________________________________________________________________________
InitStatus MpdKalmanFilter::Init()
{

   cout << "InitStatus MpdKalmanFilter::Init\n\n";

   // Check if mag. field was set - if not create the default

   fMagField = FairRunAna::Instance()->GetField();

   if (!fMagField || TMath::Abs(fMagField->GetBz(0, 0, 0)) < 0.01) {
      cout << " -I- Using the default constant magnetic field Bz = 5 kG " << endl;
      fMagField            = new MpdMultiField();
      MpdConstField *field = new MpdConstField();
      field->SetField(0., 0., 5.);                            // values are in kG
      field->SetFieldRegion(-230, 230, -230, 230, -375, 375); // cm
      ((MpdMultiField *)fMagField)->AddField(field);
      FairRunAna::Instance()->SetField(fMagField);
      fMagField->Init();
      cout << " -I- The magnetic field at (0,0,0) = (" << fMagField->GetBx(0, 0, 0) << "," << fMagField->GetBy(0, 0, 0)
           << "," << fMagField->GetBz(0, 0, 0) << ") kG" << endl;
   } else {
      cout << " -I- The magnetic field at (0,0,0) = (" << fMagField->GetBx(0, 0, 0) << "," << fMagField->GetBy(0, 0, 0)
           << "," << fMagField->GetBz(0, 0, 0) << ") kG" << endl;
      cout << "OVERHEAD FIELD INFO: " << ((MpdMultiField *)fMagField)->GetFieldList()->UncheckedAt(0) << " "
           << ((FairField *)((MpdMultiField *)fMagField)->GetFieldList()->UncheckedAt(0))->GetType() << endl;
   }

   fJacob.ResizeTo(5, 5);
   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdKalmanFilter::ReInit()
{
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdKalmanFilter::Reset()
{
   ///
   // cout << " MpdKalmanFilter::Reset  " << endl;
}

//__________________________________________________________________________
void MpdKalmanFilter::Register()
{
   ///
}

//__________________________________________________________________________
void MpdKalmanFilter::Finish()
{
   ///
}

//__________________________________________________________________________
void MpdKalmanFilter::Exec(Option_t *option)
{
   ///
}

//__________________________________________________________________________
void MpdKalmanFilter::GoOutward(MpdKalmanTrack *track)
{
   /// Propagate track in the outward direction

   TMatrixDSym weight = *track->GetWeight(); // save starting weight matrix
   // RunKalmanFilter(track);
   // track->SetLastLay();
   BackTrace(track, weight);
}

//__________________________________________________________________________
void MpdKalmanFilter::BackTrace(MpdKalmanTrack *track, TMatrixDSym &weight0)
{
   /// Propagate track using already found hits

   track->StartBack();
   Int_t nHits = track->GetNofHits();
   if (nHits == 1) track->SetWeight(weight0); // restore original weight

   TMatrixD    param(5, 1);
   TMatrixDSym weight(5), pointWeight(5);

   for (Int_t i = 1; i < nHits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)track->GetHits()->UncheckedAt(i);
      PropagateToHit(track, hit);
      Double_t dChi2 = 0.;
      dChi2          = FilterHit(track, hit, pointWeight, param);
      track->SetChi2(track->GetChi2() + dChi2);
      weight = *track->GetWeight();
      weight += pointWeight;
      track->SetWeight(weight);
      track->SetParamNew(param);
   }
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateToHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng, Bool_t local,
                                       Double_t stepBack)
{
   /// Propagate track to given hit

   // Start timing
   // if (FairRun::Instance()->GetTask("Code timer")) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);

   track->SetParam(*track->GetParamNew());
   track->SetPos(track->GetPosNew());
   track->SetNode(track->GetNodeNew());
   // TVector3 pos = fGeoScheme->GlobalPos(hit);
   Double_t sign = 1.;

   if (hit->GetType() == MpdKalmanHit::kFixedR) {
      // if (TMath::Abs(pos.Pt() - track->GetPos()) < fgkEpsilon) return kTRUE; // the "same" pos
      if (track->GetNode() == "" && TMath::Abs(hit->GetPos() - track->GetPos()) < fgkEpsilon) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kTRUE; // the "same" pos
      }
      if (!PropagateParamR(track, hit, calcLeng)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
      Double_t phi0 = track->GetParamNew(0) / track->GetPosNew();
      Double_t dphi = Proxim(phi0, track->GetParamNew(2)) - phi0;
      sign          = TMath::Sign(1., -dphi);
   } else if (hit->GetType() == MpdKalmanHit::kFixedZ) {
      // if (TMath::Abs(pos.Z() - track->GetPos()) < fgkEpsilon) return kTRUE; // the same pos
      if (TMath::Abs(hit->GetPos() - track->GetPos()) < fgkEpsilon) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kTRUE; // the same pos
      }
      if (!PropagateParamZ(track, hit, calcLeng)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
   } else {
      if (!PropagateParamP(track, hit, calcLeng, local, stepBack)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
      if (track->GetNodeNew() != "") {
         Double_t v3[7];
         SetGeantParamB(track, v3, 1);
         Double_t phi0 = TMath::ATan2(v3[1], v3[0]);
         Double_t dphi = Proxim(phi0, track->GetParamNew(2)) - phi0;
         sign          = TMath::Sign(1., -dphi);
      }
   }

   // cout << " PropagateParam ok " << track->GetPos() << " " << track->GetPosNew() << endl;
   // if (!PropagateWeight(track, hit, sign, stepBack)) return kFALSE;
   if (!PropagateWeight(track, hit, sign, stepBack)) {
      // 15-may-12
      track->SetParamNew(*track->GetParam());
      track->SetPosNew(track->GetPos());
      track->SetNodeNew(track->GetNode());
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kFALSE;
   }
   // cout << " PropagateWeight ok " << track->GetPos() << " " << track->GetPosNew() << endl;
   if (track->GetType() != MpdKalmanTrack::kEndcap && hit->GetType() == MpdKalmanHit::kFixedZ)
      track->SetType(MpdKalmanTrack::kEndcap);
   else if (track->GetType() == MpdKalmanTrack::kEndcap && hit->GetType() != MpdKalmanHit::kFixedZ)
      track->SetType(MpdKalmanTrack::kBarrel); // change of track type

   // Stop timing
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateParamToHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng,
                                            Bool_t local, Double_t stepBack)
{
   /// Propagate track parameters to given hit

   // Start timing
   // if (FairRun::Instance()->GetTask("Code timer")) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);

   track->SetParam(*track->GetParamNew());
   track->SetPos(track->GetPosNew());
   track->SetNode(track->GetNodeNew());
   // TVector3 pos = fGeoScheme->GlobalPos(hit);
   Double_t sign = 1.;

   if (hit->GetType() == MpdKalmanHit::kFixedR) {
      // if (TMath::Abs(pos.Pt() - track->GetPos()) < fgkEpsilon) return kTRUE; // the "same" pos
      if (TMath::Abs(hit->GetPos() - track->GetPos()) < fgkEpsilon) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kTRUE; // the "same" pos
      }
      if (!PropagateParamR(track, hit, calcLeng)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
      Double_t phi0 = track->GetParamNew(0) / track->GetPosNew();
      Double_t dphi = Proxim(phi0, track->GetParamNew(2)) - phi0;
      sign          = TMath::Sign(1., -dphi);
   } else if (hit->GetType() == MpdKalmanHit::kFixedZ) {
      // if (TMath::Abs(pos.Z() - track->GetPos()) < fgkEpsilon) return kTRUE; // the same pos
      if (TMath::Abs(hit->GetPos() - track->GetPos()) < fgkEpsilon) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kTRUE; // the same pos
      }
      if (!PropagateParamZ(track, hit, calcLeng)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
   } else {
      if (!PropagateParamP(track, hit, calcLeng, local, stepBack)) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
      if (track->GetNodeNew() != "") {
         Double_t v3[7];
         SetGeantParamB(track, v3, 1);
         Double_t phi0 = TMath::ATan2(v3[1], v3[0]);
         Double_t dphi = Proxim(phi0, track->GetParamNew(2)) - phi0;
         sign          = TMath::Sign(1., -dphi);
      }
   }
   /*
   //cout << " PropagateParam ok " << track->GetPos() << " " << track->GetPosNew() << endl;
   //if (!PropagateWeight(track, hit, sign, stepBack)) return kFALSE;
   if (!PropagateWeight(track, hit, sign, stepBack)) {
     // 15-may-12
     track->SetParamNew(*track->GetParam());
     track->SetPosNew(track->GetPos());
     track->SetNodeNew(track->GetNode());
     if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__);
     return kFALSE;
   }
   //cout << " PropagateWeight ok " << track->GetPos() << " " << track->GetPosNew() << endl;
   if (track->GetType() != MpdKalmanTrack::kEndcap && hit->GetType() == MpdKalmanHit::kFixedZ)
     track->SetType(MpdKalmanTrack::kEndcap);
   else if (track->GetType() == MpdKalmanTrack::kEndcap && hit->GetType() != MpdKalmanHit::kFixedZ)
     track->SetType(MpdKalmanTrack::kBarrel); // change of track type
   */
   // Stop timing
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateParamR(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng, Double_t *vert)
{
   /// Propagate parameter vector to the hit (from one R to another R)

   // TVector3 pos = fGeoScheme->GlobalPos(hit);
   // Double_t hitR = pos.Pt(), hitRphi = hitR * pos.Phi();
   Double_t hitR = hit->GetPos(), hitRphi = hit->GetMeas(0);
   track->SetNodeNew("");
   TString node = track->GetNode();

   Double_t v3[7], v3new[7], dR = 0, dir = 1, v000[3] = {0.};
   if (vert == 0x0) vert = v000;
   Int_t nTries = 0;
   track->SetParamNew(*track->GetParam());
   track->SetPosNew(track->GetPos());
   // track->GetParamNew()->Print();
   if (track->GetType() == MpdKalmanTrack::kBarrel)
      SetGeantParamB(track, v3, dir);
   else
      SetGeantParamE(track, v3, dir);

   // Propagation step estimate
   Double_t r0  = TMath::Sqrt(v3[0] * v3[0] + v3[1] * v3[1]);
   Double_t dR0 = hitR - r0;
   dir          = TMath::Sign(1.0, dR0); // direction
   if (dir < 0.) {
      for (Int_t i = 3; i < 6; ++i) v3[i] = -v3[i];
      if (hitR < 1.e-7) {
         // Going to PCA - check if is not already passed
         TVector3 tr(v3[3], v3[4], 0.0);
         TVector3 rad(v3[0] - vert[0], v3[1] - vert[1], 0.0);
         Double_t sign = TMath::Sign(1., tr * rad);
         if (sign > 0) return kFALSE;
      }
   }

   Double_t step = dR0 / TMath::Cos(track->GetParam(3));

   // sign of charge (sign of Pt if forward motion) must be changed if
   // backward extrapolation
   Double_t charge = -dir * TMath::Sign(1.0, track->GetParam(4));
   // cout << " Start:" << endl; for (Int_t i = 0; i < 7; ++i) cout << v3[i] << " "; cout << track->GetR() << " " <<
   // hit->GetR() << " " << step << endl;

   // Check if overstep or within precision
   Double_t step0 = 0.;
   do {
      step  = TMath::Abs(step);
      step0 = step;
      // Propagate parameters
      // ExtrapOneStepRungekutta(charge,step,v3,v3new);
      ExtrapOneStepHelix(charge, step, v3, v3new);
      Double_t extrapRad = TMath::Sqrt(v3new[0] * v3new[0] + v3new[1] * v3new[1]);
      if (hitR > 1.e-7) {
         dR = hitR - extrapRad;
         step *= dR0 / (extrapRad - r0);
      } else {
         // Find PCA
         // TVector3 tr(v3new[0]-v3[0],v3new[1]-v3[1],v3new[2]-v3[2]); // in 3-D
         // TVector3 rad(-v3new[0],-v3new[1],-v3new[2]);
         TVector3 tr(v3new[3], v3new[4], 0.0); // in 2-D
         TVector3 rad(v3new[0] - vert[0], v3new[1] - vert[1], 0.0);
         Double_t sign = TMath::Sign(1., tr * rad);
         // TVector3 tr(v3new[0]-v3[0],v3new[1]-v3[1]); // in 2-D
         // TVector3 rad(-v3new[0],-v3new[1]);
         // Double_t sign = -TMath::Sign(1.,tr * rad);
         dR = extrapRad * sign;
         // step *= dR0 / (-sign * extrapRad - r0);
         step *= dR0 / (-dR - r0);
      }
      nTries++;
   } while (dR * dir < 0 && TMath::Abs(dR) > fgkEpsilon);

   // for (Int_t i = 0; i < 7; ++i) cout << v3new[i] << " "; cout << nTries << endl;
   GetFromGeantParamB(track, v3new, dir);
   if (calcLeng) track->UpdateLength(step0); // update track length
   track->SetNode("");
   // track->GetParamNew()->Print();
   //  Position adjustment (until within tolerance)
   step  = 1.;
   step0 = -1.;
   while (TMath::Abs(dR) > fgkEpsilon && step > fgkEpsilon) {
      dR0              = hitR - track->GetPosNew();
      Double_t  dir1   = TMath::Sign(1.0, dR0);
      TMatrixD *parNew = track->GetParamNew();
      // Double_t dPhi = hitR > 1.e-7 ?
      // hitRphi / hitR - (*parNew)(2,0) : hitRphi - (*parNew)(2,0);
      Double_t dPhi = 0.; // 24-04-11
      if (hitR > 1.e-7 || step0 < 0)
         step = dR0 / TMath::Cos(dPhi) / TMath::Cos((*parNew)(3, 0));
      else
         step = step0;
      step = TMath::Abs(step);
      SetGeantParamB(track, v3, dir1);
      do {
         // binary search
         // ExtrapOneStepRungekutta(charge,step,v3,v3new);
         ExtrapOneStepHelix(charge, step, v3, v3new);
         Double_t extrapRad = TMath::Sqrt(v3new[0] * v3new[0] + v3new[1] * v3new[1]);
         // !!! Check this if backpropagation is needed
         if (hitR > 1.e-7) {
            if (dir1 * extrapRad < dir1 * track->GetPosNew()) {
               if (track->GetType() == MpdKalmanTrack::kBarrel) return kFALSE; // curling track
               if (dir1 * v3[2] < dir1 * v3new[2]) return kFALSE;              // curling track
            }
            dR = hitR - extrapRad;
         } else {
            // Find PCA
            // TVector3 tr(v3new[0]-v3[0],v3new[1]-v3[1],v3new[2]-v3[2]); // in 3-D
            // TVector3 rad(-v3new[0],-v3new[1],-v3new[2]);
            TVector3 tr(v3new[3], v3new[4], 0.0); // in 2-D
            TVector3 rad(v3new[0] - vert[0], v3new[1] - vert[1], 0.0);
            Double_t sign = TMath::Sign(1., tr * rad);
            dR            = sign;
            // TVector3 tr(v3new[0]-v3[0],v3new[1]-v3[1]); // in 2-D
            // TVector3 rad(-v3new[0],-v3new[1]);
            // Double_t sign = -TMath::Sign(1.,tr * rad);
            // dR = sign;
         }
         // !!!
         step0 = step;
         step /= 2.;
         nTries++;
         if (nTries > fgkTriesMax) {
            cout << " Too many tries: " << nTries << endl;
            // exit(0);
            return kFALSE;
         }
         // exit(0);
         // if (nTries > 4) cout << " *** " << nTries << " " << dPhi << " " << step << " " << extrapRad << endl;
      } while (dR * dir1 < 0);

      GetFromGeantParamB(track, v3new, dir1);
      if (calcLeng) track->UpdateLength(step0); // update track length
   }                                            // while (TMath::Abs(dR) > fgkEpsilon && step > fgkEpsilon
   /*
   if (hitR < 1.e-7) {
     cout << " Tries: " << nTries << " " << dR << " " << step << " " << track->GetPosNew() << endl;
     //for (Int_t i = 0; i < 7; ++i) cout << v3new[i] << " "; cout << endl;
     //if (track->GetParamNew(3) < 0) { cout << track->GetParamNew(3) << endl; exit(0); }
     track->GetParamNew()->Print();
   }
   */
   // if (hitR < 1.e-7) cout << "nTries: " << nTries << endl;

   track->SetNode(node);
   return kTRUE;
}

//__________________________________________________________________________
void MpdKalmanFilter::SetGeantParamB(MpdKalmanTrack *track, Double_t *v3, Double_t dir)
{
   /// Set vector of Geant3 parameters pointed to by "v3"
   /// from track parameters (for "barrel" tracks)

   TMatrixD *parNew = track->GetParamNew();
   Double_t  rad    = track->GetPosNew();
   if (track->GetNode() == "") {
      // Global coordinates
      // Double_t phi = 0;
      Double_t phi = (*parNew)(2, 0);
      if (rad > 1.e-7) phi = (*parNew)(0, 0) / rad;
      v3[0] = rad * TMath::Cos(phi); // X
      v3[1] = rad * TMath::Sin(phi); // Y
      v3[2] = (*parNew)(1, 0);       // Z
   } else {
      // Local coordinates
      // Double_t pos[3] = {(*parNew)(0,0), (*parNew)(1,0), rad};
      Double_t pos[3] = {-(*parNew)(0, 0), (*parNew)(1, 0), rad};
      if (track->GetUniqueID()) {
         // ITS
         pos[0] = -pos[0]; //
         pos[1] = rad;
         pos[2] = (*parNew)(1, 0);
      }
#ifdef _OPENMP
      omp_set_lock(&geoManagerLock);
#endif
      gGeoManager->cd(track->GetNode());
      if (track->GetUniqueID() == 0) gGeoManager->CdUp(); // TPC
      gGeoManager->LocalToMaster(pos, v3);
#ifdef _OPENMP
      omp_unset_lock(&geoManagerLock);
#endif
   }
   Double_t cosTh = TMath::Cos((*parNew)(3, 0));               // cos(Theta)
   v3[3]          = dir * TMath::Cos((*parNew)(2, 0)) * cosTh; // Px/Ptot
   v3[4]          = dir * TMath::Sin((*parNew)(2, 0)) * cosTh; // Py/Ptot
   v3[5]          = dir * TMath::Sin((*parNew)(3, 0));         // Pz/Ptot
   v3[6]          = 1. / TMath::Abs((*parNew)(4, 0)) / cosTh;  // Ptot
}

//__________________________________________________________________________
void MpdKalmanFilter::GetFromGeantParamB(MpdKalmanTrack *track, Double_t *v3, Double_t dir)
{
   /// Get track parameters from vector of Geant3 parameters pointed
   /// to by "v3" (for "barrel" tracks)

   TMatrixD *parNew = track->GetParamNew();
   if (track->GetNodeNew() == "") {
      // Global coordinates
      Double_t rad = TMath::Sqrt(v3[0] * v3[0] + v3[1] * v3[1]);
      Double_t phi = TMath::ATan2(v3[1], v3[0]);
      track->SetPosNew(rad);
      (*parNew)(0, 0) = rad * phi; // RPhi
      (*parNew)(1, 0) = v3[2];     // Z
   } else {
      // Local coordinates
      Double_t pos[3];
#ifdef _OPENMP
      omp_set_lock(&geoManagerLock);
#endif
      gGeoManager->cd(track->GetNodeNew());
      if (track->GetUniqueID() == 0) gGeoManager->CdUp(); // !!! for TPC
      gGeoManager->MasterToLocal(v3, pos);
#ifdef _OPENMP
      omp_unset_lock(&geoManagerLock);
#endif
      //(*parNew)(0,0) = pos[0];
      if (track->GetUniqueID()) {
         // ITS
         (*parNew)(0, 0) = pos[0]; //
         (*parNew)(1, 0) = pos[2];
         track->SetPosNew(pos[1]);
      } else {
         (*parNew)(0, 0) = -pos[0]; //
         (*parNew)(1, 0) = pos[1];
         track->SetPosNew(pos[2]);
      }
   }
   Double_t phi0   = (*parNew)(2, 0);
   (*parNew)(2, 0) = TMath::ATan2(dir * v3[4], dir * v3[3]); // track Phi
   (*parNew)(2, 0) = Proxim(phi0, (*parNew)(2, 0));          // adjust angle for continuity
   (*parNew)(3, 0) = TMath::ASin(dir * v3[5]);               // Theta
   (*parNew)(4, 0) = 1. / v3[6] / TMath::Cos((*parNew)(3, 0)) * TMath::Sign(1., (*parNew)(4, 0)); // q/Pt
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateParamZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng)
{
   /// Propagate parameter vector to the hit at fixed Z

   if (!fNumer && track->GetType() == MpdKalmanTrack::kEndcap && hit->GetType() == MpdKalmanHit::kFixedZ) {
      // Analytical extrapolation
      // TMatrixD tmp = *track->GetParamNew();
      AnalParamZ(track, hit, calcLeng);
      // track->SetParamNew(tmp);
      track->SetPosNew(hit->GetDist());
      return kTRUE;
   }

   Double_t hitZ = hit->GetPos();

   // Propagation step estimate
   Double_t dZ0  = hitZ - track->GetZ();
   Double_t step = dZ0 / TMath::Sin(track->GetParam(3));

   // Track parameter extrapolation to Z using Rungekutta algorithm.
   // dZ0 = TMath::Abs(hitZ) - TMath::Abs(track->GetZ());
   // Double_t dir = TMath::Sign(1.0,dZ0); // direction
   Double_t dir = TMath::Sign(1.0, TMath::Abs(hitZ) - TMath::Abs(track->GetZ())); // direction
   // sign of charge (sign of Pt if forward motion) must be changed if
   // backward extrapolation
   Double_t charge = -dir * TMath::Sign(1.0, track->GetParam(4));
   // if (track->GetDirection() == TpcLheKalmanTrack::kOutward) charge = -charge;
   Double_t v3[7], v3new[7], dZ = 0;
   Int_t    nTries = 0;
   track->SetParamNew(*track->GetParam());
   track->SetPosNew(track->GetPos());
   // track->GetParamNew()->Print();
   if (track->GetType() == MpdKalmanTrack::kBarrel)
      SetGeantParamB(track, v3, dir);
   else
      SetGeantParamE(track, v3, dir);
   // cout << " Start:" << endl; for (Int_t i = 0; i < 7; ++i) cout << v3[i] << " "; cout << track->GetR() << " " <<
   // hit->GetR() << " " << step << endl;

   // Check if overstep or within precision
   Double_t step0 = 0.;
   do {
      step  = TMath::Abs(step);
      step0 = step;
      // Propagate parameters
      // ExtrapOneStepRungekutta(charge,step,v3,v3new);
      ExtrapOneStepHelix(charge, step, v3, v3new);
      Double_t extrapZ = v3new[2];
      // for (Int_t i = 0; i < 7; ++i) cout << v3new[i] << " "; cout << extrapRad << endl;
      // dZ = hitZ - extrapZ;
      dZ = TMath::Abs(hitZ) - TMath::Abs(extrapZ);
      step *= dZ0 / (extrapZ - track->GetZ());
      nTries++;
   } while (dZ * dir < 0 && TMath::Abs(dZ) > fgkEpsilon);

   GetFromGeantParamE(track, v3new, dir);
   if (calcLeng) track->UpdateLength(step0); // update track length
   // track->GetParamNew()->Print();
   //  Position adjustment (until within tolerance)
   while (TMath::Abs(dZ) > fgkEpsilon) {
      dZ0 = hitZ - track->GetPosNew();
      // Double_t dir1 = TMath::Sign(1.0,dZ0);
      Double_t  dir1   = TMath::Sign(1.0, TMath::Abs(hitZ) - TMath::Abs(track->GetPosNew()));
      TMatrixD *parNew = track->GetParamNew();
      step             = dZ0 / TMath::Sin((*parNew)(3, 0));
      step             = TMath::Abs(step);
      SetGeantParamE(track, v3, dir1);
      do {
         // binary search
         // ExtrapOneStepRungekutta(charge,step,v3,v3new);
         ExtrapOneStepHelix(charge, step, v3, v3new);
         Double_t extrapZ = v3new[2];
         // !!! Check if backpropagation is needed
         // if (dir1*extrapRad < dir1*track->GetPosNew()) return kFALSE; // curling track
         // !!!
         // dZ = hitZ - extrapZ;
         dZ    = TMath::Abs(hitZ) - TMath::Abs(extrapZ);
         step0 = step;
         step /= 2.;
         nTries++;
         if (nTries > fgkTriesMax) {
            cout << " Too many tries: " << nTries << endl;
            // exit(0);
            return kFALSE;
         }
         // if (nTries > 4) cout << " *** " << nTries << " " << dPhi << " " << step << " " << extrapRad << endl;
      } while (dZ * dir1 < 0);

      GetFromGeantParamE(track, v3new, dir1);
      if (calcLeng) track->UpdateLength(step0); // update track length
   }
   // cout << " Tries: " << nTries << " " << dR << endl;
   // for (Int_t i = 0; i < 7; ++i) cout << v3new[i] << " "; cout << endl;
   // if (track->GetParamNew(3) < 0) { cout << track->GetParamNew(3) << endl; exit(0); }
   // track->GetParamNew()->Print();

   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::AnalParamZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng)
{
   /// Propagate parameter vector to the hit at fixed Z analytically

   static Int_t    first = 1;
   static Double_t coef  = 0.;
   if (first) {
      first       = 0;
      Double_t bZ = TMath::Abs(fMagField->GetBz(0, 0, 0));
      coef        = 0.3 * 0.01 * bZ / 10;
   }

   Double_t  dz     = hit->GetDist() - track->GetPos();
   Double_t  r      = 1. / track->GetParam(4) / coef;
   Double_t  ph0    = track->GetParam(2);
   Double_t  ph     = ph0 + dz / r / TMath::Tan(track->GetParam(3));
   Double_t  x      = track->GetParam(0) + r * (TMath::Sin(ph) - TMath::Sin(ph0));
   Double_t  y      = track->GetParam(1) - r * (TMath::Cos(ph) - TMath::Cos(ph0));
   TMatrixD &parNew = *track->GetParamNew();
   parNew(0, 0)     = x;
   parNew(1, 0)     = y;
   parNew(2, 0)     = ph;
   parNew(3, 0)     = track->GetParam(3);
   parNew(4, 0)     = track->GetParam(4);
   if (calcLeng) {
      Double_t step = r * (ph - ph0);
      track->UpdateLength(TMath::Abs(step / TMath::Cos(track->GetParam(3)))); // update track length
   }
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::AnalParamX(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng, Bool_t local)
{
   /// Propagate parameter vector to the hit at fixed X analytically

   static Int_t    first = 1;
   static Double_t coef  = 0.;
   if (first) {
      first       = 0;
      Double_t bZ = TMath::Abs(fMagField->GetBz(0, 0, 0));
      coef        = 0.3 * 0.01 * bZ / 10;
   }

   // cout << track->GetNode() << " " << fGeoScheme->Path(hit->GetDetectorID()) << endl;
   //  Coordinate transformation
   gGeoManager->cd(fGeoScheme->Path(hit->GetDetectorID()));
   if (track->GetUniqueID() == 0) gGeoManager->CdUp(); // !!! for TPC
   Double_t v3[7], v31[3], dir = 1.0;
   fGeoScheme->GlobalPos(hit).GetXYZ(v3);
   gGeoManager->MasterToLocal(v3, v31);
   v31[0]      = -v31[0]; //
   Double_t dx = v31[2], posNew = v31[2];
   if (track->GetUniqueID()) dx = posNew = v31[1]; // ITS
   Double_t phSec   = MpdKalmanFilter::Proxim(track->GetParam(2), TMath::ATan2(v3[1], v3[0]));
   TString  pathHit = gGeoManager->GetPath();
   TString  pathTra = track->GetNode();
   if (track->GetUniqueID() == 0) {
      // if (track->GetUniqueID() >= 0) {
      //  TPC
      Int_t ip = pathTra.Last('/');
      if (ip > 0) pathTra.Replace(ip, pathTra.Length() - ip + 1, "");
   } // else cout << pathHit << " " << pathTra << endl;

   Double_t y0        = track->GetParam(0);
   Double_t z0        = track->GetParam(1);
   Int_t    changeVol = 0;
   if (pathTra != pathHit) {
      if (track->GetDirection() == MpdKalmanTrack::kInward) dir = -1.0;
      track->SetParamNew(*track->GetParam());
      track->SetPosNew(track->GetPos());
      SetGeantParamB(track, v3, dir); // global coordinates
      gGeoManager->cd(pathHit);
      gGeoManager->MasterToLocal(v3, v31);
      if (track->GetUniqueID() == 0) {
         // TPC
         v31[0] = -v31[0]; //
         y0     = v31[0];
         z0     = v31[1];
         dx -= v31[2];
      } else {
         v31[0] = -v31[0]; //
         y0     = v31[0];
         z0     = v31[2];
         dx -= v31[1];
      }
      // Change of volume - set flag to use numerical weight propagation
      changeVol = 1;
      //} else dx -= track->GetPos();
   } else {
      dx -= track->GetPos();
      if (track->GetUniqueID()) y0 = -y0; // ITS
   }

   Double_t r = 1. / track->GetParam(4) / coef;
   // Double_t r = -1. / track->GetParam(4) / coef;
   Double_t ph0 = track->GetParam(2) - phSec;
   // Double_t ph0 = -track->GetParam(2) + phSec;
   Double_t arg = TMath::Sin(ph0) + dx / r;
   if (TMath::Abs(arg) > 0.999999) return kFALSE;
   Double_t  ph     = TMath::ASin(arg);
   Double_t  y      = y0 - r * (TMath::Cos(ph) - TMath::Cos(ph0));
   Double_t  step   = r * (ph - ph0);
   Double_t  z      = z0 + step * TMath::Tan(track->GetParam(3));
   TMatrixD &parNew = *track->GetParamNew();
   parNew(0, 0)     = y;
   if (track->GetUniqueID()) parNew(0, 0) = -y; //
   parNew(1, 0) = z;
   parNew(2, 0) = ph + phSec;
   // parNew(2,0) = -ph + phSec;
   parNew(3, 0) = track->GetParam(3);
   parNew(4, 0) = track->GetParam(4);
   track->SetPosNew(posNew);
   if (calcLeng) track->UpdateLength(TMath::Abs(step / TMath::Cos(track->GetParam(3)))); // update track length
   track->SetChi2Vertex(dx);                                                             // just to pass value
   if (changeVol) track->SetChi2Vertex(999999);
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateParamP(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng, Bool_t local,
                                        Double_t stepBack)
{
   /// Propagate parameter vector to the hit at fixed plane
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);

   if (!fNumer && track->GetType() == MpdKalmanTrack::kBarrel && hit->GetType() == MpdKalmanHit::kFixedP) {
      // Analytical extrapolation
      if (local) track->SetNodeNew(fGeoScheme->Path(hit->GetDetectorID())); // in local coordinates
      AnalParamX(track, hit, calcLeng, kTRUE);
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kTRUE;
   }

   TVector3 pos      = fGeoScheme->GlobalPos(hit);
   TVector3 norm     = fGeoScheme->Normal(hit);
   Double_t plane[6] = {0};
   pos.GetXYZ(plane);
   norm.GetXYZ(&plane[3]);
   // cout << fGeoScheme->Path(hit->GetDetectorID()) << endl;
   // ITS if (local) track->SetNodeNew(fGeoScheme->Path(hit->GetDetectorID() % 1000000 | 1)); // in local coordinates
   if (local) track->SetNodeNew(fGeoScheme->Path(hit->GetDetectorID())); // in local coordinates
   if (!PropagateParamP(track, plane, calcLeng, local, stepBack)) {
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kFALSE;
   }

   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateParamP(MpdKalmanTrack *track, const Double_t *plane, Bool_t calcLeng, Bool_t local,
                                        Double_t stepBack)
{
   /// Propagate parameter vector to the plane
   /// plane[0-2] - x, y, z of a point on the plane;
   /// plane[3-5] - direction cosines of the normal to the plane
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   Double_t dir = 1., norm = 1., scalar = 0., normal[3] = {0}, normLeng = 0;

   if (track->Momentum() < 0.001) return kFALSE; // AZ - safety 19.04.2020

   for (Int_t i = 0; i < 3; ++i) {
      normal[i] = plane[i + 3];
      normLeng += normal[i] * normal[i];
      scalar += plane[i] * normal[i];
   }
   normLeng = TMath::Sqrt(normLeng); // length of the normal vector
   // if (scalar < 0) norm = -1; // normal direction
   if (scalar < 0 && (TMath::Abs(plane[0]) > 3. || TMath::Abs(plane[1]) > 3.)) norm = -1; // normal direction
   normLeng *= norm;

   for (Int_t i = 0; i < 3; ++i) normal[i] /= normLeng; // "canonical" normal vector

   Double_t v3[7], v3new[7];
   Int_t    nTries = 0;
   track->SetParamNew(*track->GetParam());
   track->SetPosNew(track->GetPos());

   // Double_t dZ0 = hit->GetZ() - track->GetZ();
   // Double_t dir = TMath::Sign(1.0,dZ0); // direction
   if (track->GetDirection() == MpdKalmanTrack::kInward) dir = -1.;

   Double_t charge = -dir * TMath::Sign(1.0, track->GetParam(4));
   if (track->GetType() == MpdKalmanTrack::kBarrel)
      SetGeantParamB(track, v3, dir);
   else {
      SetGeantParamE(track, v3, dir);
      // Fatal(" PropagateParamP "," !!! Local track parameterization is not implemented for Endcap tracks !!!");
   }

   Double_t d = 0.; // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
   for (Int_t i = 0; i < 3; ++i) d -= plane[i] * normal[i];
   // Distance to the plane
   Double_t dist0 = d;
   for (Int_t i = 0; i < 3; ++i) dist0 += v3[i] * normal[i];
   if (TMath::Abs(dist0) < fgkEpsilon) {
      GetFromGeantParamB(track, v3, dir); // 11-may-12
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kTRUE; // 02-11-11
   }
   if (dist0 * dir > 0) {
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kFALSE; // track already passed the plane
   }
   // if (dist0 < 0) return kFALSE; // track last point behind plane (farther from interaction point)
   //  Angle between momentum and plane normal
   Double_t cosA = 0;
   for (Int_t i = 0; i < 3; ++i) cosA += v3[i + 3] * normal[i];
   if (TMath::Abs(cosA) < 0.5) cosA = 0.5; // to constrain step size

   // Propagation step estimate
   Double_t step = dist0 / cosA;
   // step = TMath::Max (TMath::Abs(step), 0.5*TMath::Abs(TMath::Sqrt(v3[0]*v3[0]+v3[1]*v3[1])-
   //						      TMath::Sqrt(plane[0]*plane[0]+plane[1]*plane[1])));
   if (stepBack > 0.0) step = 0.9 * stepBack;

   // cout << " " << norm << " " << normLeng << " " << step << " " << dist0 << endl;
   // cout << " Start:" << endl; for (Int_t i = 0; i < 7; ++i) cout << v3[i] << " "; cout << track->GetR() << " " <<
   // hit->GetR() << " " << step << endl;

   // Check if overstep or within precision
   Double_t step0 = 0., dist = 0.;
   do {
      step  = TMath::Abs(step);
      step0 = step;
      // Propagate parameters
      // ExtrapOneStepRungekutta(charge,step,v3,v3new);
      ExtrapOneStepHelix(charge, step, v3, v3new);
      dist = d;
      for (Int_t i = 0; i < 3; ++i) dist += v3new[i] * normal[i];
      // for (Int_t i = 0; i < 7; ++i) cout << v3new[i] << " "; cout << endl;
      // step *= dist0 / (TMath::Abs(dist0) - TMath::Abs(dist));
      step *= dist0 / (dist0 - dist);
      nTries++;
      //} while (dist * dir < 0 && TMath::Abs(dist) > fgkEpsilon);
   } while (dist * dir > 0 && TMath::Abs(dist) > fgkEpsilon);
   //} while (dist * dir * norm > 0 && TMath::Abs(dist) > fgkEpsilon);

   GetFromGeantParamE(track, v3new, dir);
   if (calcLeng) track->UpdateLength(step0); // update track length
   // track->GetParamNew()->Print();

   // Position adjustment (until within tolerance)
   Double_t dir1 = dir;
   // dist0 = TMath::Abs(dist);
   // AZ dist0 = dist;
   if (TMath::Abs(dist) < TMath::Abs(dist0)) dist0 = dist; // AZ

   while (TMath::Abs(dist) > fgkEpsilon) {
      dist = d;
      cosA = 0.;
      // dir = 1.;
      for (Int_t i = 0; i < 3; ++i) {
         dist += v3new[i] * normal[i];
         cosA += v3[i + 3] * normal[i];
      }
      step = dist / cosA;
      step = TMath::Abs(step);
      dir1 = TMath::Sign(1., -dist);
      // SetGeantParamE(track,v3,dir);
      SetGeantParamE(track, v3, dir1);
      do {
         // binary search
         // ExtrapOneStepRungekutta(charge,step,v3,v3new);
         ExtrapOneStepHelix(charge, step, v3, v3new);
         step0 = step;
         step /= 2.;
         nTries++;
         if (nTries > fgkTriesMax) {
            cout << " Too many tries: " << nTries << endl;
            // exit(0);
            if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
            return kFALSE;
         }
         dist = d;
         for (Int_t i = 0; i < 3; ++i) dist += v3new[i] * normal[i];
         if (TMath::Abs(dist) < fgkEpsilon) break; // 02-11-11
         // if (TMath::Abs(dist) > dist0) return kFALSE; // track going in different direction
         if (dist > 0 && dist0 > 0 && dist > dist0 || dist < 0 && dist0 < 0 && dist < dist0)
            return kFALSE; // track going in different direction
                           //} while (dist * dir < 0);
      } while (dist * dir > 0);
      //} while (dist * dir1 < 0);
      //} while (dist*dir*norm > 0);

      // GetFromGeantParamE(track,v3new,dir);
      GetFromGeantParamE(track, v3new, dir1);
      if (calcLeng) track->UpdateLength(step0); // update track length
   }

   // if (track->GetType() == MpdKalmanTrack::kBarrel) GetFromGeantParamB(track,v3new,dir);
   // 2-may-12 if (track->GetType() == MpdKalmanTrack::kBarrel) GetFromGeantParamB(track,v3new,dir1);
   GetFromGeantParamB(track, v3new, dir1);

   // cout << " Tries: " << nTries << endl;
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return kTRUE;
}

//__________________________________________________________________________
void MpdKalmanFilter::SetGeantParamE(MpdKalmanTrack *track, Double_t *v3, Double_t dir)
{
   /// Set vector of Geant3 parameters pointed to by "v3"
   /// from track parameters (for "end-cap" tracks)

   TMatrixD *parNew = track->GetParamNew();
   v3[0]            = (*parNew)(0, 0);                           // X
   v3[1]            = (*parNew)(1, 0);                           // Y
   v3[2]            = track->GetPosNew();                        // Z
   Double_t cosTh   = TMath::Cos((*parNew)(3, 0));               // cos(Theta)
   v3[3]            = dir * TMath::Cos((*parNew)(2, 0)) * cosTh; // Px/Ptot
   v3[4]            = dir * TMath::Sin((*parNew)(2, 0)) * cosTh; // Py/Ptot
   v3[5]            = dir * TMath::Sin((*parNew)(3, 0));         // Pz/Ptot
   v3[6]            = 1. / TMath::Abs((*parNew)(4, 0)) / cosTh;  // Ptot
}

//__________________________________________________________________________
void MpdKalmanFilter::GetFromGeantParamE(MpdKalmanTrack *track, Double_t *v3, Double_t dir)
{
   /// Get track parameters from vector of Geant3 parameters pointed
   /// to by "v3" (for "end-cap" tracks)

   TMatrixD *parNew = track->GetParamNew();
   track->SetPosNew(v3[2]); // Z
   (*parNew)(0, 0) = v3[0]; // X
   (*parNew)(1, 0) = v3[1]; // Y
   Double_t phi0   = (*parNew)(2, 0);
   (*parNew)(2, 0) = TMath::ATan2(dir * v3[4], dir * v3[3]); // track Phi
   (*parNew)(2, 0) = Proxim(phi0, (*parNew)(2, 0));          // adjust angle for continuity
   (*parNew)(3, 0) = TMath::ASin(dir * v3[5]);               // Theta
   (*parNew)(4, 0) = 1. / v3[6] / TMath::Cos((*parNew)(3, 0)) * TMath::Sign(1., (*parNew)(4, 0)); // q/Pt
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::FindPca(MpdKalmanTrack *track, Double_t *vert)
{
   /// Propagate parameter vector to PCA w.r.t. point vert

   // Double_t v3[7], v3new[7], dR = 0, dir0 = -1, dir = -1;
   Double_t v3[7], v3new[7], dR = 0, dir0 = 1, dir = 1;
   Int_t    nTries = 0;
   track->SetParamNew(*track->GetParam());
   track->SetPosNew(track->GetPos());
   // track->GetParamNew()->Print();
   if (track->GetType() == MpdKalmanTrack::kBarrel)
      SetGeantParamB(track, v3, dir);
   else
      SetGeantParamE(track, v3, dir);

   Double_t step = 0.1, dist2d = 0;
   Double_t charge = -dir * TMath::Sign(1.0, track->GetParam(4));
   // for (Int_t i = 3; i < 6; ++i) v3[i] = -v3[i];
   // cout << " Start:" << endl; for (Int_t i = 0; i < 7; ++i) cout << v3[i] << " "; cout << track->GetR() << " " <<
   // hit->GetR() << " " << step << endl;

   // Go to the approximate position of PCA
   TVector3 rad(v3[0] - vert[0], v3[1] - vert[1], 0.0);
   TVector3 tr0(v3[3], v3[4], v3[5]);
   dist2d = tr0 * rad / tr0.Pt();
   if (dist2d > 0) {
      // Change direction
      dir = -dir;
      for (Int_t i = 3; i < 6; ++i) v3[i] = -v3[i];
   }
   // dir = TMath::Sign(1.,-dist2d);
   step = TMath::Abs(dist2d / TMath::Sin(tr0.Theta()));

   TVector3 tr;
   for (Int_t i = 0; i < 1; ++i) {
      // ExtrapOneStepRungekutta(charge,step,v3,v3new);
      ExtrapOneStepHelix(dir * charge, step, v3, v3new);
      rad.SetXYZ(v3new[0] - vert[0], v3new[1] - vert[1], 0);
      tr.SetXYZ(v3new[3], v3new[4], v3new[5]);
      dist2d = tr * rad / tr.Pt();
      // if (dist2d * dir > 0) {
      if (dist2d > 0) {
         // Change direction
         dir = -dir;
         for (Int_t j = 3; j < 6; ++j) v3new[j] = -v3new[j];
      }
      // dir = TMath::Sign(1.,-dist2d);
      step = TMath::Abs(dist2d / TMath::Sin(tr.Theta()));
      nTries++;
   }

   // Position adjustment (until within tolerance)
   // step = 1.;
   while (step > fgkEpsilon) {
      for (Int_t i = 0; i < 7; ++i) v3[i] = v3new[i];
      // binary search
      // ExtrapOneStepRungekutta(charge,step,v3,v3new);
      ExtrapOneStepHelix(dir * charge, step, v3, v3new);
      // TVector3 tr(v3new[0]-v3[0],v3new[1]-v3[1]);
      tr.SetXYZ(v3new[3], v3new[4], v3new[5]);
      rad.SetXYZ(v3new[0] - vert[0], v3new[1] - vert[1], 0.);
      dist2d = tr * rad;
      // if (dist2d * dir > 0) {
      if (dist2d > 0) {
         // Change direction
         dir = -dir;
         step /= 2.;
         for (Int_t i = 3; i < 6; ++i) v3new[i] = -v3new[i];
      }
      ++nTries;
      if (nTries > fgkTriesMax) {
         cout << " Too many tries: " << nTries << endl;
         // exit(0);
         return kFALSE;
      }
      // exit(0);
      // if (nTries > 4) cout << " *** " << nTries << " " << dPhi << " " << step << " " << extrapRad << endl;
   }

   // if (track->GetType() == MpdKalmanTrack::kFixedR) GetFromGeantParamR(track,v3new,dir0);
   // else GetFromGeantParamZ(track,v3new,dir0);
   tr.SetXYZ(v3new[3], v3new[4], v3new[5]);
   dir = 1;
   if (tr0 * tr < 0) dir = -1;
   if (track->GetType() == MpdKalmanTrack::kBarrel)
      GetFromGeantParamB(track, v3new, dir);
   else
      GetFromGeantParamE(track, v3new, dir);

   /*cout << " Tries: " << nTries << " " << dist2d << " " << dir << " "
        << track->GetTrackID() << " "
        << TMath::Sqrt((v3new[0]-vert[0])*(v3new[0]-vert[0])+
             (v3new[1]-vert[1])*(v3new[1]-vert[1])) << endl;
   cout << tr0.X() << " " << tr0.Y() << " " << tr0.Z() << " "
   << tr.X() << " " << tr.Y() << " " << tr.Z() << endl;*/
   return kTRUE;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Proxim(Double_t phi0, Double_t phi, Double_t scale)
{
   /// Adjust angle phi to be "around" phi0 - to avoid discontinuity around +- Pi:
   /// scale - scale factor due to stereo angle

   Double_t dPhi = phi0 - phi;
   // if (TMath::Abs(dPhi) > TMath::Pi()) phi += TMath::Pi() * 2 * TMath::Sign(1.,dPhi);
   if (TMath::Abs(dPhi) > TMath::Pi() * scale) phi += TMath::Pi() * scale * TMath::Sign(2., dPhi);
   return phi;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit, Double_t scale)
{
   /// Adjust hit coord. R-Phi to be "around" track R-Phi - to avoid
   /// discontinuity around +- Pi: scale - scale factor due to stereo angle

   // TVector3 pos = fGeoScheme->GlobalPos(hit);
   // Double_t hitPhi = pos.Phi();
   Double_t hitPhi = hit->GetMeas(0) / hit->GetPos();
   Double_t phi0   = track->GetParamNew(0) / track->GetPosNew();
   return hit->GetPos() * Proxim(phi0, hitPhi, scale);
   // return pos.Pt() * Proxim(phi0, hitPhi, scale);
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Proxim(const MpdKalmanHit *hit0, const MpdKalmanHit *hit, Double_t scale)
{
   /// Adjust hit coord. R-Phi to be "around" hit0 R-Phi - to avoid
   /// discontinuity around +- Pi: scale - scale factor due to stereo angle

   // TVector3 pos0 = fGeoScheme->GlobalPos(hit0);
   // TVector3 pos = fGeoScheme->GlobalPos(hit);
   // Double_t phi0 = pos0.Phi();
   // Double_t phi = pos.Phi();
   Double_t phi0 = hit0->GetMeas(0) / hit0->GetPos();
   Double_t phi  = hit->GetMeas(0) / hit->GetPos();
   return hit->GetPos() * Proxim(phi0, phi, scale);
   // return pos.Pt() * Proxim(phi0, phi, scale);
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::PropagateWeight(MpdKalmanTrack *track, const MpdKalmanHit *hit, Double_t sign,
                                        Double_t stepBack)
{
   /// Propagation of the weight matrix
   /// W = DtWD, where D is Jacobian
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);

   // TMatrixD jacob(5,5);
   fJacob          = 0.0;
   TMatrixD &jacob = fJacob;

   if (!fNumer && track->GetType() == MpdKalmanTrack::kEndcap && hit->GetType() == MpdKalmanHit::kFixedZ) {
      // Analytical track propagation (Jacobian computation)
      AnalyticJacob(track, hit, jacob);
      // jacob.Print();
      jacob.Invert();
      // jacob.Print();
      TMatrixD tmp(*track->GetWeight(), TMatrixD::kMult, jacob); // WD
      TMatrixD weight1(jacob, TMatrixD::kTransposeMult, tmp);    // DtWD
      track->SetWeight(weight1);
      // jacob = 0;
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kTRUE;
   }
   if (!fNumer && track->GetType() == MpdKalmanTrack::kBarrel && hit->GetType() == MpdKalmanHit::kFixedP &&
       track->GetChi2Vertex() < 99999) {
      // Analytical track propagation (Jacobian computation) - no ITS at present !!!
      AnalyticJacobX(track, hit, jacob);
      // jacob.Print();
      jacob.Invert();
      // jacob.Print();
      TMatrixD tmp(*track->GetWeight(), TMatrixD::kMult, jacob); // WD
      TMatrixD weight1(jacob, TMatrixD::kTransposeMult, tmp);    // DtWD
      track->SetWeight(weight1);
      // jacob = 0;
      if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
      return kTRUE;
   }

   MpdKalmanTrack tr     = *track;
   MpdKalmanHit   hitTmp = *hit;

   // Save initial and propagated parameters
   TMatrixD param0    = *tr.GetParam();
   TMatrixD paramNew0 = *tr.GetParamNew();
   Double_t pos0      = tr.GetPosNew();
   Bool_t   invert    = kTRUE;
   // if (tr.GetDirection() == MpdKalmanTrack::kInward &&
   //   tr.GetNode() != "" && tr.GetNodeNew() != "") {
   if (tr.GetDirection() == MpdKalmanTrack::kInward && tr.GetNode() != "" && tr.GetNodeNew() != "" &&
       tr.GetNode() != tr.GetNodeNew()) {
      // Change direction
      param0    = *tr.GetParamNew();
      paramNew0 = *tr.GetParam();
      tr.SetPos(pos0);
      invert          = kFALSE;
      TString detName = tr.GetNode();
      if (track->GetUniqueID()) {
         // ITS
         detName = detName(16, detName.Length());
         detName += "#0";
      }
      hitTmp.SetDetectorID(fGeoScheme->DetId(detName));
      tr.SetNode(track->GetNodeNew());
      tr.SetNodeNew(track->GetNode());
      tr.SetDirection(MpdKalmanTrack::kOutward);
   }

   // param0.Print();
   // paramNew0.Print();
   TMatrixDSym *covar = tr.Weight2Cov(); // get covariance matrix
   Double_t     rNew = 0., phi0 = 0.;
   if (hit->GetType() == MpdKalmanHit::kFixedR) {
      // if (track->GetType() == MpdKalmanTrack::kBarrel) {
      rNew = pos0; // barrel
      phi0 = paramNew0(0, 0) / rNew;
   }
   // covar->Print();

   Bool_t local = kFALSE;
   if (tr.GetNodeNew() != "") local = kTRUE;
   Double_t dPar;
   // Loop over parameters to find change of the propagated vs initial ones
   Bool_t ok = kTRUE;
   for (Int_t i = 0; i < 5; ++i) {
      dPar = TMath::Sqrt((*covar)(i, i));
      if (i < 4)
         dPar = TMath::Min(dPar, 0.1);
      else
         dPar = TMath::Min(dPar, 0.1 * TMath::Abs(param0(4, 0)));
      tr.SetParam(param0);
      if (i == 4)
         dPar *= TMath::Sign(1., -param0(4, 0)); // 1/p
      else if (i == 2)
         dPar *= sign;
      else if (i == 3)
         dPar *= TMath::Sign(1., -param0(3, 0)); // dip-angle
      tr.SetParam(i, param0(i, 0) + dPar);
      if (hit->GetType() == MpdKalmanHit::kFixedR)
         ok = PropagateParamR(&tr, &hitTmp, kFALSE);
      else if (hit->GetType() == MpdKalmanHit::kFixedZ)
         ok = PropagateParamZ(&tr, &hitTmp, kFALSE);
      // else ok = PropagateParamP(track, hit, kFALSE, kTRUE, stepBack);
      // else ok = PropagateParamP(&tr, &hitTmp, kFALSE, kTRUE);
      else
         ok = PropagateParamP(&tr, &hitTmp, kFALSE, local);
      if (!ok) {
         if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
         return kFALSE;
      }
      for (Int_t j = 0; j < 5; ++j) {
         if (j == 0 && hit->GetType() == MpdKalmanHit::kFixedR) {
            Double_t dPhi = Proxim(phi0, tr.GetParamNew(j) / rNew) - phi0;
            jacob(j, i)   = dPhi * rNew / dPar;
            /*} else if (j == 0 && hit->GetType() == MpdKalmanHit::kFixedP) {
            Double_t phi = Proxim(phi0, track->GetParamNew(j)/track->GetPosNew());
            jacob(j,i) = (phi * track->GetPosNew() - phi0 * rNew) / dPar;*/
         } else
            jacob(j, i) = (tr.GetParamNew(j) - paramNew0(j, 0)) / dPar;
         // if (hitTmp.GetDetectorID()) cout << i << " " << j << " " << dPar << " " << tr.GetParamNew(j) << " " <<
         // paramNew0(j,0) << endl;
         //    << " " << track->GetPos() << " " << track->GetPosNew() << " " << jacob(j,i) << endl;
      }
   }

   // jacob.Print();
   if (invert) jacob.Invert();

   TMatrixD tmp(*tr.GetWeight(), TMatrixD::kMult, jacob);  // WD
   TMatrixD weight1(jacob, TMatrixD::kTransposeMult, tmp); // DtWD
   // track->SetWeight(TMatrixD(jacob,TMatrixD::kTransposeMult,weight1)); // DtWD
   track->SetWeight(weight1);

   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::AnalyticJacob(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixD &jacob)
{
   /// Analytical computation of Jacobian matrix (matrix of parameter derivatives)

   static Int_t    first = 1;
   static Double_t coef  = 0.;
   if (first) {
      first       = 0;
      Double_t bZ = TMath::Abs(fMagField->GetBz(0, 0, 0));
      coef        = 0.3 * 0.01 * bZ / 10;
   }

   Double_t cosPh1 = TMath::Cos(track->GetParamNew(2));
   Double_t cosPh0 = TMath::Cos(track->GetParam(2));
   Double_t dCosPh = cosPh1 - cosPh0;
   Double_t sinPh1 = TMath::Sin(track->GetParamNew(2));
   Double_t sinPh0 = TMath::Sin(track->GetParam(2));
   Double_t dSinPh = sinPh1 - sinPh0;
   Double_t r      = 1. / track->GetParamNew(4) / coef;
   Double_t dz     = track->GetPosNew() - track->GetPos();
   Double_t sinTh  = TMath::Sin(track->GetParamNew(3));
   Double_t sinTh2 = sinTh * sinTh;
   Double_t tanTh  = TMath::Tan(track->GetParamNew(3));

   // dX1 / dP0
   jacob(0, 0) = 1.;                                      // dX1/dX0
   jacob(0, 2) = r * dCosPh;                              // dX1/dPh0
   jacob(0, 3) = -cosPh1 * dz / sinTh2;                   // dX1/dTh
   jacob(0, 4) = r * (-r * dSinPh + cosPh1 * dz / tanTh); // dX1/d(1/R)
   jacob(0, 4) *= coef;

   // dY1 / dP0
   jacob(1, 1) = 1.;                                     // dY1/dY0
   jacob(1, 2) = r * dSinPh;                             // dY1/dPh0
   jacob(1, 3) = -sinPh1 * dz / sinTh2;                  // dY1/dTh0
   jacob(1, 4) = r * (r * dCosPh + sinPh1 * dz / tanTh); // dY1/d(1/R)
   jacob(1, 4) *= coef;

   // dPh1 / dP0
   jacob(2, 2) = 1.;               // dPh1/dPh0
   jacob(2, 3) = -dz / sinTh2 / r; // dPh1/dTh0
   jacob(2, 4) = dz / tanTh;       // dPh1/d(1/R)
   jacob(2, 4) *= coef;

   jacob(3, 3) = jacob(4, 4) = 1.;
   return kTRUE;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::AnalyticJacobX(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixD &jacob)
{
   /// Analytical computation of Jacobian matrix (matrix of parameter derivatives)

   static Int_t    first = 1;
   static Double_t coef  = 0.;
   if (first) {
      first       = 0;
      Double_t bZ = TMath::Abs(fMagField->GetBz(0, 0, 0));
      coef        = 0.3 * 0.01 * bZ / 10;
   }

   // Sector angle
   gGeoManager->cd(fGeoScheme->Path(hit->GetDetectorID()));
   if (track->GetUniqueID() == 0) gGeoManager->CdUp(); // !!! for TPC
   Double_t v3[7];
   fGeoScheme->GlobalPos(hit).GetXYZ(v3);
   // Double_t phSec = TMath::ATan2(v3[1],v3[0]);
   Double_t phSec = MpdKalmanFilter::Proxim(track->GetParam(2), TMath::ATan2(v3[1], v3[0]));

   Double_t ph1    = track->GetParamNew(2) - phSec;
   Double_t ph0    = track->GetParam(2) - phSec;
   Double_t cosPh1 = TMath::Cos(ph1);
   Double_t cosPh0 = TMath::Cos(ph0);
   Double_t dCosPh = cosPh1 - cosPh0;
   Double_t tanPh1 = TMath::Tan(ph1);
   Double_t r      = 1. / track->GetParamNew(4) / coef;
   // Double_t r = -1. / track->GetParamNew(4) / coef;
   // Double_t dx = track->GetPosNew() - track->GetPos();
   Double_t dx     = track->GetChi2Vertex();
   Double_t cosTh  = TMath::Cos(track->GetParamNew(3));
   Double_t cosTh2 = cosTh * cosTh;
   Double_t tanTh  = TMath::Tan(track->GetParamNew(3));
   Double_t dPh    = ph1 - ph0;

   // dY1 / dP0
   jacob(0, 0) = 1.; // dY1/dY0
   /*
   if (track->GetNode() != "") {
     TString pathTra = track->GetNode();
     TString pathTraNew = track->GetNodeNew();
     Int_t ip = pathTra.Last('/');
     Int_t ipnew = pathTraNew.Last('/');
     if (ip > 0) pathTra.Replace(ip,pathTra.Length()-ip+1,"");
     if (ipnew > 0) pathTraNew.Replace(ipnew,pathTraNew.Length()-ipnew+1,"");
     if (pathTra != pathTraNew)
       { jacob(0,0) = TMath::Cos(TMath::DegToRad()*30.); cout << pathTra << " " << pathTraNew << endl; }
   }
   */
   jacob(0, 2) = r * sin(dPh) / cosPh1;          // dY1/dPh0
   jacob(0, 3) = 0;                              // dY1/dTh
   jacob(0, 4) = r * (r * dCosPh + tanPh1 * dx); // dY1/d(1/R)
   jacob(0, 4) *= coef;

   // dZ1 / dP0
   jacob(1, 1) = 1;                            // dZ1/dZ0
   jacob(1, 2) = -r * tanTh * dCosPh / cosPh1; // dZ1/dPh0
   // AZ jacob(1,3) = -r * dPh / cosTh2; // dZ1/dTh0
   jacob(1, 3) = r * dPh / cosTh2;                     // dZ1/dTh0
   jacob(1, 4) = r * tanTh * (-r * dPh + dx / cosPh1); // dZ1/d(1/R)
   jacob(1, 4) *= coef;

   // dPh1 / dP0
   jacob(2, 2) = cosPh0 / cosPh1; // dPh1/dPh0
   jacob(2, 3) = 0;               // dPh1/dTh0
   jacob(2, 4) = dx / cosPh1;     // dPh1/d(1/R)
   jacob(2, 4) *= coef;

   jacob(3, 3) = jacob(4, 4) = 1.;
   return kTRUE;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                                    TMatrixD &paramTmp)
{
   /// Hit filtering
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);

   Double_t chi2 = 0.0;
   // AZ if (hit->GetType() == MpdKalmanHit::kFixedR) return FilterHitR(track, hit, pointWeight, paramTmp);
   if (hit->GetNofDim() > 1 && hit->GetType() != MpdKalmanHit::kFixedZ)
      chi2 = FilterHitR(track, hit, pointWeight, paramTmp);
   else if (hit->GetType() == MpdKalmanHit::kFixedZ)
      chi2 = FilterHitZ(track, hit, pointWeight, paramTmp);
   else
      chi2 = FilterStrip(track, hit, pointWeight, paramTmp);
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return chi2;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterHitR(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                                     TMatrixD &paramTmp)
{
   /// Compute Chi2 increment if the "barrel" hit is assigned to the track
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   TMatrixD point(5, 1);
   if (hit->GetType() == MpdKalmanHit::kFixedR)
      point(0, 0) = Proxim(track, hit);
   else
      point(0, 0) = hit->GetMeas(0);
   point(1, 0)       = hit->GetMeas(1);
   pointWeight       = 0.0;
   pointWeight(0, 0) = 1. / hit->GetErr(0) / hit->GetErr(0);
   pointWeight(1, 1) = 1. / hit->GetErr(1) / hit->GetErr(1);

   // Solving linear system (W+U)p' = U(m-p) + (W+U)p
   // TMatrixD wu = *fWeight;
   TMatrixDSym wu = *track->GetWeight();
   wu += pointWeight; // W+U
   paramTmp = point;
   // paramTmp -= *fTrackParNew; // m-p
   paramTmp -= *track->GetParamNew();                      // m-p
   TMatrixD right(pointWeight, TMatrixD::kMult, paramTmp); // U(m-p)
   // TMatrixD right1(wu,TMatrixD::kMult,*fTrackParNew); // (W+U)p
   TMatrixD right1(wu, TMatrixD::kMult, *track->GetParamNew()); // (W+U)p
   right += right1;                                             // U(m-p) + (W+U)p

   // wu.Print();
   // wu.Invert();
   Int_t iok = 0;
   fgKF->MnvertLocal(wu.GetMatrixArray(), 5, 5, 5, iok);
   // wu.Print();

   paramTmp = TMatrixD(wu, TMatrixD::kMult, right);

   right1 = paramTmp;
   right1 -= point; // p'-m
   point = paramTmp;
   point -= *track->GetParamNew();                                // p'-p
   right = TMatrixD(*track->GetWeight(), TMatrixD::kMult, point); // W(p'-p)
   TMatrixD value(point, TMatrixD::kTransposeMult, right);        // (p'-p)'W(p'-p)
   Double_t dChi2 = value(0, 0);
   right          = TMatrixD(pointWeight, TMatrixD::kMult, right1);    // U(p'-m)
   value          = TMatrixD(right1, TMatrixD::kTransposeMult, right); // (p'-m)'U(p'-m)
   dChi2 += value(0, 0);

   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   return dChi2;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterHitZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                                     TMatrixD &parFilt)
{
   /// Compute Chi2 increment if the "end-cap" hit is assigned to the track
   /// (assuming 1-d hit measurement)

   TMatrixD meas(2, 1); // measurements
   TMatrixD h(2, 5);    // transformation matrix from hit meas. to track params
   TMatrixD u(2, 2);    // weight matrix of measurements
   // h(0,0) = -TMath::Sin(hit->GetPhi());
   // h(0,1) = TMath::Cos(hit->GetPhi());
   h(0, 0) = -hit->GetCosSin(1);
   h(0, 1) = hit->GetCosSin(0);
   // point(0,0) = Proxim(track,hit);
   meas(0, 0) = hit->GetMeas(0);
   // pointWeight(0,0) = 1. / hit->GetRphiErr() / hit->GetRphiErr();
   u(0, 0) = 1. / hit->GetErr(0) / hit->GetErr(0);
   if (hit->GetNofDim() > 1) {
      meas(1, 0) = hit->GetMeas(1);
      u(1, 1)    = 1. / hit->GetErr(1) / hit->GetErr(1);
      // h(1,0) = TMath::Sin(hit->GetPhi());
      // h(1,1) = -TMath::Cos(hit->GetPhi());
      h(1, 0) = hit->GetCosSin(0);
      h(1, 1) = hit->GetCosSin(1);
   }

   // Solving linear system (W+H'UH)p' = H'U(m-Hp) + (W+H'UH)p
   // TMatrixD uh(pointWeight,TMatrixD::kMult,h); // UH
   TMatrixD    uh(u, TMatrixD::kMult, h);            // UH
   TMatrixD    huh(h, TMatrixD::kTransposeMult, uh); // H'UH
   TMatrixDSym huhSym = pointWeight;
   huhSym.SetMatrixArray(huh.GetMatrixArray());

   TMatrixDSym wu = *track->GetWeight();
   wu += huhSym; // W+H'UH

   TMatrixD hp(h, TMatrixD::kMult, *track->GetParamNew()); // Hp
   TMatrixD tmp = meas;                                    // m
   tmp -= hp;                                              // m-Hp

   // TMatrixD right0(pointWeight,TMatrixD::kMult,paramTmp); // U(m-Hp)
   TMatrixD right0(u, TMatrixD::kMult, tmp);            // U(m-Hp)
   TMatrixD right(h, TMatrixD::kTransposeMult, right0); // H'U(m-Hp)

   TMatrixD right1(wu, TMatrixD::kMult, *track->GetParamNew()); // (W+H'UH)p
   right += right1;                                             // H'U(m-Hp) + (W+H'UH)p

   // wu.Print();
   // wu.Invert();
   Int_t iok = 0;
   fgKF->MnvertLocal(wu.GetMatrixArray(), 5, 5, 5, iok);
   // wu.Print();

   parFilt = TMatrixD(wu, TMatrixD::kMult, right); // p'
   // parFilt.Print();

   // right1 = h;
   // right1 *= parFilt; // Hp'
   // right1 -= meas; // Hp'-m
   tmp = TMatrixD(h, TMatrixD::kMult, parFilt);
   tmp -= meas;

   TMatrixD point = parFilt;
   point -= *track->GetParamNew();                                // p'-p
   right = TMatrixD(*track->GetWeight(), TMatrixD::kMult, point); // W(p'-p)
   TMatrixD value(point, TMatrixD::kTransposeMult, right);        // (p'-p)'W(p'-p)
   Double_t dChi2 = value(0, 0);
   // right = TMatrixD(pointWeight,TMatrixD::kMult,right1); // U(Hp'-m)
   // right = TMatrixD(u,TMatrixD::kMult,right1); // U(Hp'-m)
   // value = TMatrixD(right1,TMatrixD::kTransposeMult,right); // (Hp'-m)'U(Hp'-m)
   meas  = TMatrixD(u, TMatrixD::kMult, tmp);             // U(Hp'-m)
   value = TMatrixD(tmp, TMatrixD::kTransposeMult, meas); // (Hp'-m)'U(Hp'-m)
   dChi2 += value(0, 0);
   pointWeight = huhSym;

   // pointWeight.Print();
   // cout<<"dChi2====="<<dChi2<<'\n';
   // sleep(1);

   return dChi2;
}

/*
//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterHitZ(MpdKalmanTrack *track, const MpdKalmanHit *hit,
                TMatrixDSym &pointWeight, TMatrixD &parFilt)
{
  /// Compute Chi2 increment if the "end-cap" hit is assigned to the track
  /// (assuming 1-d hit measurement)

  //double meas[3][1], h[3][5], u[3][2];
  double* u = (double *)calloc(6, sizeof(double));
  double* h = (double *)calloc(15, sizeof(double));
  double* meas = (double *)calloc(3, sizeof(double));
  //TMatrixD meas(2,1); // measurements
  //TMatrixD h(2,5); // transformation matrix from hit meas. to track params
  //TMatrixD u(2,2); // weight matrix of measurements
  MKL_INT mM = 2, nM = 1, mH = 2, nH = 5, mU = 2, nU = 2, mP = 5, nP = 1, info;

  Double_t HitCosSin1 = hit->GetCosSin(1), HitCosSin0 = hit->GetCosSin(0), HitErr0 = hit->GetErr(0);

  //h[0][0] = -HitCosSin1;
  //h[0][1] = HitCosSin0;
  h[0] = -HitCosSin1;
  h[1] = HitCosSin0;

  //meas[0][0] = hit->GetMeas(0);
  meas[0] = hit->GetMeas(0);

  //u[0][0] = 1. / HitErr0 / HitErr0;
  u[0] = 1. / HitErr0 / HitErr0;

  if (hit->GetNofDim() > 1) {
    //meas[1][0] = hit->GetMeas(1);
    meas[1] = hit->GetMeas(1);

    Double_t HitErr1 = hit->GetErr(1);

    //u[1][1] = 1. / HitErr1 / HitErr1;
    u[3] = 1. / HitErr1 / HitErr1;

    //h[1][0] = HitCosSin0;
    //h[1][1] = HitCosSin1;
    h[5] = HitCosSin0;
    h[6] = HitCosSin1;
  }

  // Solving linear system (W+H'UH)p' = H'U(m-Hp) + (W+H'UH)p

  //TMatrixD uh(u,TMatrixD::kMult,h); // UH
  //double uh[3][5];
  double uh[15], huhSym[30], tmp[3], param[6], tmp1[3], right[6];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mU, nH, nU, 1, u, nU, h, nH, 0, uh, nH);
  //cout<<"UH:
"<<uh[0][0]<<';'<<uh[0][1]<<';'<<uh[0][2]<<';'<<uh[0][3]<<';'<<uh[0][4]<<';'<<uh[1][0]<<';'<<uh[1][1]<<';'<<uh[1][2]<<';'<<uh[1][3]<<';'<<uh[1][4]<<'\n';

  //TMatrixD huh(h,TMatrixD::kTransposeMult,uh); // H'UH
  //double huhSym[6][5];
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nH,  nH,  mH, 1, h, nH, uh, nH, 0, pointWeight.GetMatrixArray(),
nH);
  //cout<<"HUH:
"<<huhSym[0][0]<<';'<<huhSym[0][1]<<';'<<huhSym[0][2]<<';'<<huhSym[0][3]<<';'<<huhSym[0][4]<<';'<<huhSym[1][0]<<';'<<huhSym[1][1]<<';'<<huhSym[1][2]<<';'<<huhSym[1][3]<<';'<<huhSym[1][4]
  //
<<huhSym[2][0]<<';'<<huhSym[2][1]<<';'<<huhSym[2][2]<<';'<<huhSym[2][3]<<';'<<huhSym[2][4]<<';'<<huhSym[3][0]<<';'<<huhSym[3][1]<<';'<<huhSym[3][2]<<';'<<huhSym[3][3]<<';'<<huhSym[3][4]
  //     <<huhSym[4][0]<<';'<<huhSym[4][1]<<';'<<huhSym[4][2]<<';'<<huhSym[4][3]<<';'<<huhSym[4][4]<<'\n';

  //pointWeight = huhSym;
  //memcpy(pointWeight.GetMatrixArray(), huhSym, 25*sizeof(double));
  //pointWeight.Print();

  //TMatrixDSym huhSym = pointWeight;
  //huhSym.SetMatrixArray(huh.GetMatrixArray());
  //TMatrixDSym wu = *track->GetWeight();
  //wu += huhSym; // W+H'UH
  //double wu[6][5];
  memcpy(huhSym, (*track->GetWeight()).GetMatrixArray(), 25*sizeof(double));
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nH, nH, mH, 1, h, nH, uh, nH, 1, huhSym, nH);
  //cout<<"WU:
"<<wu[0][0]<<';'<<wu[0][1]<<';'<<wu[0][2]<<';'<<wu[0][3]<<';'<<wu[0][4]<<';'<<wu[1][0]<<';'<<wu[1][1]<<';'<<wu[1][2]<<';'<<wu[1][3]<<';'<<wu[1][4]
  //
<<wu[2][0]<<';'<<wu[2][1]<<';'<<wu[2][2]<<';'<<wu[2][3]<<';'<<wu[2][4]<<';'<<wu[3][0]<<';'<<wu[3][1]<<';'<<wu[3][2]<<';'<<wu[3][3]<<';'<<wu[3][4]
  //     <<wu[4][0]<<';'<<wu[4][1]<<';'<<wu[4][2]<<';'<<wu[4][3]<<';'<<wu[4][4]<<'\n';

  //TMatrixD hp(h,TMatrixD::kMult,*track->GetParamNew()); // Hp
  //TMatrixD tmp = meas; // m
  //tmp -= hp; // m-Hp
  //double tmp[3][1];
  memcpy(tmp, meas, 2*sizeof(double));
  //double param[6][1];
  memcpy(param, (*track->GetParamNew()).GetMatrixArray(), 5*sizeof(double));
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mH, nP, nH, -1, h, nH, param, nP, 1, tmp, nP);
  //cout<<"m-Hp: "<<tmp[0][0]<<';'<<tmp[1][0]<<'\n';

  //TMatrixD right0(u,TMatrixD::kMult,tmp); // U(m-Hp)
  //double tmp1[3][1];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mU, nP, nU, 1, u, nU, tmp, nP, 0, tmp1, nP);
  //cout<<"U(m-Hp): "<<tmp1[0][0]<<';'<<tmp1[1][0]<<'\n';

  //TMatrixD right(h,TMatrixD::kTransposeMult,right0); // H'U(m-Hp)
  //double right[6][1];
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,  nH, nP, mH, 1, h, nH, tmp1, nP, 0, right, nP);
  //cout<<"H'U(m-Hp):  "<<right[0][0]<<';'<<right[1][0]<<';'<<right[2][0]<<';'<<right[3][0]<<';'<<right[4][0]<<'\n';

  //TMatrixD right1(wu,TMatrixD::kMult,*track->GetParamNew()); // (W+H'UH)p
  cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, nH, nP, 1, huhSym, nH, param, nP, 1, right, nP);
  //cout<<"H'U(m-Hp) + (W+H'UH)p:
"<<right[0][0]<<';'<<right[1][0]<<';'<<right[2][0]<<';'<<right[3][0]<<';'<<right[4][0]<<'\n';

  //Int_t iok = 0;
  //fgKF->MnvertLocal(wu.GetMatrixArray(), 5, 5, 5, iok);
  //parFilt = TMatrixD(wu,TMatrixD::kMult,right); // p'
  int* ipiv = new int[nH];
  //double LU[6][5];
  //memcpy(LU, wu, 25*sizeof(double));;
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, nH, 1, huhSym,  nH, ipiv, right,  nH);
  //if (info > 0){
  //    cout<<"Error\n";
  //    return 0;
  //}
  //info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 5, 5, (double*)LU, 5, &ipiv);
  //info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nH, right_hand, (double*)LU, nH, &ipiv, (double*)right, nH);
  memcpy(parFilt.GetMatrixArray(), right, 5*sizeof(double));
  //parFilt.Print();

  //tmp = TMatrixD(h,TMatrixD::kMult,parFilt);
  //tmp -= meas;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mH, nP, nH, 1, h, nH, right, nP, -1, meas, nP);
  //cout<<"meas: "<<meas[0][0]<<';'<<meas[1][0]<<'\n';

  //TMatrixD point = parFilt;
  //point -= *track->GetParamNew(); // p'-p
  double value[2] = {1, 0};
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mP, nP, nP, 1, right, nP, value, nP, -1, param, nP);
  //cout<<"p'-p: "<<param[0][0]<<';'<<param[1][0]<<';'<<param[2][0]<<';'<<param[3][0]<<';'<<param[4][0]<<'\n';

  //right = TMatrixD(*track->GetWeight(),TMatrixD::kMult,point); // W(p'-p)
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nH, nP, nH, 1, (*track->GetWeight()).GetMatrixArray(), nH,
param, nP, 0, right, nP);
  //cout<<"W(p'-p):  "<<right[0][0]<<';'<<right[1][0]<<';'<<right[2][0]<<';'<<right[3][0]<<';'<<right[4][0]<<'\n';

  //TMatrixD value(point,TMatrixD::kTransposeMult,right); // (p'-p)'W(p'-p)
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nP, nP, mP, 1, param, nP, right, nP, 0, value, nP);
  //cout<<"(p'-p)'W(p'-p): "<<value[0][0]<<'\n';

  //meas = TMatrixD(u,TMatrixD::kMult,tmp); // U(Hp'-m)
  Double_t dChi2 = value[0];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mU, nP, nU, 1, u, nU, meas, nP, 0, tmp, nP);
  //cout<<"U(Hp'-m): "<<tmp[0][0]<<';'<<tmp[1][0]<<'\n';

  //value = TMatrixD(tmp,TMatrixD::kTransposeMult,meas); // (Hp'-m)'U(Hp'-m)
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nP, nP, mU, 1, meas, nP, tmp, nP, 0, value, nP);
  //cout<<"(Hp'-m)'U(Hp'-m): "<<value[0][0]<<'\n';

  dChi2 += value[0];

  //cout<<"dChi2====="<<dChi2<<'\n';

  ///sleep(5);

  return dChi2;
}*/

//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterStrip(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                                      TMatrixD &parFilt)
{
   /// Kalman filtering for "strip" hits

   TMatrixD meas(2, 1); // measurements
   TMatrixD h(2, 5);    // transformation matrix from track params to measurements
   TMatrixD u(2, 2);    // weight matrix of measurements
   h(0, 0) = hit->GetCosSin(0);
   h(0, 1) = hit->GetCosSin(1);
   // h(0,0) = -TMath::Sin(hit->GetPhi());
   // h(0,1) = TMath::Cos(hit->GetPhi());
   meas(0, 0) = hit->GetMeas(0);
   u(0, 0)    = 1. / hit->GetErr(0) / hit->GetErr(0);
   /*
   if (hit->GetNofDim() > 1) {
     meas(1,0) = hit->GetR();
     u(1,1) = 1. / hit->GetRerr() / hit->GetRerr();
     //h(1,0) = TMath::Sin(hit->GetPhi());
     //h(1,1) = -TMath::Cos(hit->GetPhi());
     h(1,0) = TMath::Cos(hit->GetPhi());
     h(1,1) = TMath::Sin(hit->GetPhi());
   }
   */

   // Solving linear system (W+H'UH)p' = H'U(m-Hp) + (W+H'UH)p
   // TMatrixD uh(pointWeight,TMatrixD::kMult,h); // UH
   TMatrixD    uh(u, TMatrixD::kMult, h);            // UH
   TMatrixD    huh(h, TMatrixD::kTransposeMult, uh); // H'UH
   TMatrixDSym huhSym = pointWeight;
   huhSym.SetMatrixArray(huh.GetMatrixArray());

   TMatrixDSym wu = *track->GetWeight();
   wu += huhSym; // W+H'UH

   TMatrixD hp(h, TMatrixD::kMult, *track->GetParamNew()); // Hp
   TMatrixD tmp = meas;                                    // m
   tmp -= hp;                                              // m-Hp
   // tmp.Print();

   // TMatrixD right0(pointWeight,TMatrixD::kMult,paramTmp); // U(m-Hp)
   TMatrixD right0(u, TMatrixD::kMult, tmp);            // U(m-Hp)
   TMatrixD right(h, TMatrixD::kTransposeMult, right0); // H'U(m-Hp)

   TMatrixD right1(wu, TMatrixD::kMult, *track->GetParamNew()); // (W+H'UH)p
   right += right1;                                             // H'U(m-Hp) + (W+H'UH)p

   // wu.Print();
   // wu.Invert();
   Int_t iok = 0;
   fgKF->MnvertLocal(wu.GetMatrixArray(), 5, 5, 5, iok);
   // wu.Print();

   parFilt = TMatrixD(wu, TMatrixD::kMult, right); // p'

   // right1 = h;
   // right1 *= parFilt; // Hp'
   // right1 -= meas; // Hp'-m
   tmp = TMatrixD(h, TMatrixD::kMult, parFilt);
   tmp -= meas;

   TMatrixD point = parFilt;
   point -= *track->GetParamNew();                                // p'-p
   right = TMatrixD(*track->GetWeight(), TMatrixD::kMult, point); // W(p'-p)
   TMatrixD value(point, TMatrixD::kTransposeMult, right);        // (p'-p)'W(p'-p)
   Double_t dChi2 = value(0, 0);
   // right = TMatrixD(pointWeight,TMatrixD::kMult,right1); // U(Hp'-m)
   // right = TMatrixD(u,TMatrixD::kMult,right1); // U(Hp'-m)
   // value = TMatrixD(right1,TMatrixD::kTransposeMult,right); // (Hp'-m)'U(Hp'-m)
   meas  = TMatrixD(u, TMatrixD::kMult, tmp);             // U(Hp'-m)
   value = TMatrixD(tmp, TMatrixD::kTransposeMult, meas); // (Hp'-m)'U(Hp'-m)
   dChi2 += value(0, 0);
   pointWeight = huhSym;
   return dChi2;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::FilterStripLocal(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                                           TMatrixD &parFilt, Double_t &posNew)
{
   /// Kalman filtering for "strip" hits (in local detector coordinates)

   TMatrixD meas(2, 1); // measurements
   TMatrixD h(2, 5);    // transformation matrix from track params to measurements
   TMatrixD u(2, 2);    // weight matrix of measurements
   if (hit->GetNofDim() == 1) {
      h(0, 0)    = hit->GetCosSin(0);
      h(0, 1)    = hit->GetCosSin(1);
      meas(0, 0) = hit->GetMeas(0);
      u(0, 0)    = 1. / hit->GetErr(0) / hit->GetErr(0);
   } else {
      h(0, 0)    = 1;
      h(1, 1)    = 1;
      meas(0, 0) = hit->GetMeas(0);
      meas(1, 0) = hit->GetMeas(1);
      u(0, 0)    = hit->GetErr(0) * hit->GetErr(0);
      u(0, 1) = u(1, 0) = -u(0, 0) / TMath::Tan(hit->GetPhi());
      u(1, 1)           = TMath::Abs(u(1, 0) / TMath::Tan(hit->GetPhi())) +
                u(0, 0) / TMath::Sin(hit->GetPhi()) / TMath::Sin(hit->GetPhi());
      u.Invert();
   }
   /*
   if (hit->GetNofDim() > 1) {
     meas(1,0) = hit->GetR();
     u(1,1) = 1. / hit->GetRerr() / hit->GetRerr();
     //h(1,0) = TMath::Sin(hit->GetPhi());
     //h(1,1) = -TMath::Cos(hit->GetPhi());
     h(1,0) = TMath::Cos(hit->GetPhi());
     h(1,1) = TMath::Sin(hit->GetPhi());
   }
   */

   // Transform track params to local coordinates
   /*
   Double_t r = track->GetPosNew();
   Double_t phi0 = track->GetParamNew(0) / r;
   Double_t v3loc[3], xyz[3] = {r*TMath::Cos(phi0), r*TMath::Sin(phi0), track->GetParamNew(1)};
   TGeoNode *node = gGeoManager->FindNode(xyz[0], xyz[1], xyz[2]);
   TString path = gGeoManager->GetPath();
   if (!path.Contains("sensor")) {
     path += "/sts01sensor";
     Int_t layer = hit->GetLayer() / 2 + 1;
     path += layer;
     path += "_0";
     gGeoManager->cd(path);
   }
   cout << path << endl;
   gGeoManager->MasterToLocal(xyz, v3loc);
   */

   // Solving linear system (W+H'UH)p' = H'U(m-Hp) + (W+H'UH)p
   TMatrixD    uh(u, TMatrixD::kMult, h);            // UH
   TMatrixD    huh(h, TMatrixD::kTransposeMult, uh); // H'UH
   TMatrixDSym huhSym = pointWeight;
   huhSym.SetMatrixArray(huh.GetMatrixArray());

   TMatrixDSym wu = *track->GetWeight();
   wu += huhSym; // W+H'UH

   /*
   TMatrixD trackPar = *track->GetParamNew();
   TMatrixD trackParLoc = *track->GetParamNew();
   trackParLoc(0,0) = v3loc[0]; // Xlocal
   trackParLoc(1,0) = v3loc[2]; // Zlocal
   */

   TMatrixD hp(h, TMatrixD::kMult, *track->GetParamNew()); // Hp
   // TMatrixD hp(h,TMatrixD::kMult,trackParLoc); // Hp
   TMatrixD tmp = meas; // m
   tmp -= hp;           // m-Hp
   // tmp.Print();

   TMatrixD right0(u, TMatrixD::kMult, tmp);            // U(m-Hp)
   TMatrixD right(h, TMatrixD::kTransposeMult, right0); // H'U(m-Hp)

   TMatrixD right1(wu, TMatrixD::kMult, *track->GetParamNew()); // (W+H'UH)p
   // TMatrixD right1(wu,TMatrixD::kMult,trackPar); // (W+H'UH)p
   right += right1; // H'U(m-Hp) + (W+H'UH)p

   // wu.Print();
   // wu.Invert();
   Int_t iok = 0;
   MnvertLocal(wu.GetMatrixArray(), 5, 5, 5, iok);
   // wu.Print();

   parFilt = TMatrixD(wu, TMatrixD::kMult, right); // p'
   /*
   TMatrixD parFiltLoc = parFilt;
   phi0 = parFilt(0,0) / r;
   Double_t z = parFilt(1,0);
   Double_t xyz1[3] = {r*TMath::Cos(phi0), r*TMath::Sin(phi0), z};
   gGeoManager->MasterToLocal(xyz1, v3loc);
   parFiltLoc(0,0) = v3loc[0]; // Xlocal
   parFiltLoc(1,0) = v3loc[2]; // Zlocal
   */

   // right1 = h;
   // right1 *= parFilt; // Hp'
   // right1 -= meas; // Hp'-m
   tmp = TMatrixD(h, TMatrixD::kMult, parFilt);
   // tmp = TMatrixD(h,TMatrixD::kMult,parFiltLoc);
   tmp -= meas;

   TMatrixD point = parFilt;
   point -= *track->GetParamNew(); // p'-p
   // point -= trackPar; // p'-p
   right = TMatrixD(*track->GetWeight(), TMatrixD::kMult, point); // W(p'-p)
   TMatrixD value(point, TMatrixD::kTransposeMult, right);        // (p'-p)'W(p'-p)
   Double_t dChi2 = value(0, 0);
   // right = TMatrixD(pointWeight,TMatrixD::kMult,right1); // U(Hp'-m)
   // right = TMatrixD(u,TMatrixD::kMult,right1); // U(Hp'-m)
   // value = TMatrixD(right1,TMatrixD::kTransposeMult,right); // (Hp'-m)'U(Hp'-m)
   meas  = TMatrixD(u, TMatrixD::kMult, tmp);             // U(Hp'-m)
   value = TMatrixD(tmp, TMatrixD::kTransposeMult, meas); // (Hp'-m)'U(Hp'-m)
   dChi2 += value(0, 0);
   pointWeight = huhSym;

   // Transform track params to global system
   /*
   v3loc[0] = parFilt(0,0);
   v3loc[2] = parFilt(1,0);
   gGeoManager->LocalToMaster(v3loc, xyz);
   parFilt(1,0) = xyz[2];
   Double_t phi = TMath::ATan2 (xyz[1], xyz[0]);
   r = xyz[0] * xyz[0] + xyz[1] * xyz[1];
   r = TMath::Sqrt(r);
   parFilt(0,0) = r * Proxim(phi0,phi);
   */
   // posNew = r;

   return dChi2;
}

//__________________________________________________________________________
Bool_t MpdKalmanFilter::Refit(MpdKalmanTrack *track, Int_t iDir, Bool_t updLeng)
{
   /// Refit track using track hits (from first to last if iDir == -1)

#pragma omp critical
   {
      track->GetHits()->Sort();
   }
   // cout << track->GetChi2() << endl;
   track->SetChi2(0.);
   // fTrackDir = kOutward;
   //*fWeight = *fCovar;
   // Int_t iok;
   // MpdKalmanFilter::Instance()->MnvertLocal(fWeight->GetMatrixArray(), 5, 5, 5, iok);

   // if (GetNofHits() == 1) return;
   TMatrixDSym *w     = track->GetWeight();
   Int_t        nHits = track->GetNofHits(), nHits2 = nHits * nHits;
   ;
   for (Int_t i = 0; i < 5; ++i) {
      for (Int_t j = i; j < 5; ++j) {
         // if (j == i) (*w)(i,j) /= 1000.;
         // if (j == i) (*w)(i,j) /= nHits;
         if (j == i) (*w)(i, j) /= nHits2;
         // if (j == i) (*w)(i,j) /= 1;
         else
            (*w)(i, j) = (*w)(j, i) = 0.;
      }
   }

   TMatrixD    param(5, 1);
   TMatrixDSym weight(5), pointWeight(5);
   // Int_t ibeg = iDir > 0 ? 1 : nHits-2;
   Int_t ibeg = iDir > 0 ? 0 : nHits - 1;
   Int_t iend = iDir > 0 ? nHits - 1 : 0;
   iend += iDir;

   MpdKalmanHit *hit = 0x0;
   // for (Int_t i = ibeg; i >= iend; i+=iDir) {
   for (Int_t i = ibeg; i != iend; i += iDir) {
      hit = (MpdKalmanHit *)track->GetHits()->UncheckedAt(i);
      // if (!MpdKalmanFilter::Instance()->PropagateToHit(track, hit, updLeng)) return kFALSE;
      if (!MpdKalmanFilter::Instance()->PropagateToHit(track, hit, updLeng, kTRUE)) return kFALSE;
      Double_t dChi2 = FilterHit(track, hit, pointWeight, param);
      track->SetChi2(track->GetChi2() + dChi2);
      weight = *track->GetWeight();
      weight += pointWeight;
      track->SetWeight(weight);
      track->SetParamNew(param);
      // cout << i << " " << dChi2 << " " << 1./track->GetParamNew(4) << " " << hit->GetLayer() << " " << hit->GetDist()
      // << endl;
   }
   // if (hit) {
   if (0) {
      TVector3 pos = fGeoScheme->GlobalPos(hit);
      cout << " Hit: " << pos.Pt() * pos.Phi() << " " << pos.Pt() << " " << pos.Z() << endl;
      cout << " Track: " << nHits << " " << track->GetParamNew(0) << " " << track->GetPosNew() << " "
           << track->GetParamNew(1) << " " << track->GetChi2() << endl;
   }
   return kTRUE;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Interp(Int_t nDim, const Double_t *mom, const Double_t *dedxm, Double_t p)
{
   // Energy loss interpolation

   // if (p < mom[0]) return dedxm[0];
   if (p > mom[nDim - 1]) return dedxm[nDim - 1];

   // Binary search
   Int_t i1 = 0, i2 = nDim - 1, i = i2;
   if (p < mom[0]) {
      // Extrapolation
      i1 = 0;
      i2 = 1;
   } else {
      do {
         i = i1 + (i2 - i1) / 2;
         if (p > mom[i])
            i1 = i;
         else
            i2 = i;
      } while (i2 - i1 > 1);
   }
   // Linear interpolation
   return dedxm[i1] + (dedxm[i2] - dedxm[i1]) / (mom[i2] - mom[i1]) * (p - mom[i1]);
}

//__________________________________________________________________________
void MpdKalmanFilter::Convert(const MpdKalmanTrack *track1, MpdKalmanTrack *track2)
{
   // Convert track of diff. types (from "R" to "Z" at the moment)

   if (track1->GetType() == track2->GetType()) return;

   track2->SetPos(track1->GetParamNew(1));                                        // Z
   track2->SetParam(1, track1->GetPosNew());                                      // R
   track2->SetParam(3, TMath::PiOver2() - track1->GetParamNew(3));                // Theta
   track2->SetParam(4, track1->GetParamNew(4) * TMath::Sin(track2->GetParam(3))); // q/P
}

//__________________________________________________________________________
void MpdKalmanFilter::ExtrapOneStepRungekutta(Double_t charge, Double_t step, Double_t *vect, Double_t *vout)
{
   /// <pre>
   ///	******************************************************************
   ///	*								 *
   ///	*  Runge-Kutta method for tracking a particle through a magnetic *
   ///	*  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of	 *
   ///	*  Standards, procedure 25.5.20)				 *
   ///	*								 *
   ///	*  Input parameters						 *
   ///	*	CHARGE    Particle charge				 *
   ///	*	STEP	  Step size					 *
   ///	*	VECT	  Initial co-ords,direction cosines,momentum	 *
   ///	*  Output parameters						 *
   ///	*	VOUT	  Output co-ords,direction cosines,momentum	 *
   ///	*  User routine called  					 *
   ///	*	CALL GUFLD(X,F) 					 *
   ///	*								 *
   ///	*    ==>Called by : <USER>, GUSWIM				 *
   ///	*	Authors    R.Brun, M.Hansroul  *********		 *
   ///	*		   V.Perevoztchikov (CUT STEP implementation)	 *
   ///	*								 *
   ///	*								 *
   ///	******************************************************************
   /// </pre>

   Double_t h2, h4, f[4];
   Double_t xyzt[3], a, b, c, ph, ph2;
   Double_t secxs[4], secys[4], seczs[4], hxp[3];
   Double_t g1, g2, g3, g4, g5, g6, ang2, dxt, dyt, dzt;
   Double_t est, at, bt, ct, cba;
   Double_t f1, f2, f3, f4, rho, tet, hnorm, hp, rho1, sint, cost;

   Double_t x;
   Double_t y;
   Double_t z;

   Double_t xt;
   Double_t yt;
   Double_t zt;

   Double_t maxit  = 1992;
   Double_t maxcut = 11;

   const Double_t kdlt   = 1e-4;
   const Double_t kdlt32 = kdlt / 32.;
   const Double_t kthird = 1. / 3.;
   const Double_t khalf  = 0.5;
   const Double_t kec    = 2.9979251e-4;

   const Double_t kpisqua = 9.86960440109;
   const Int_t    kix     = 0;
   const Int_t    kiy     = 1;
   const Int_t    kiz     = 2;
   const Int_t    kipx    = 3;
   const Int_t    kipy    = 4;
   const Int_t    kipz    = 5;

   // *.
   // *.    ------------------------------------------------------------------
   // *.
   // *             this constant is for units cm,gev/c and kgauss
   // *
   Int_t iter = 0;
   Int_t ncut = 0;
   for (Int_t j = 0; j < 7; j++) vout[j] = vect[j];

   Double_t pinv = kec * charge / vect[6];
   Double_t tl   = 0.;
   Double_t h    = step;
   Double_t rest;

   do {
      rest = step - tl;
      if (TMath::Abs(h) > TMath::Abs(rest)) h = rest;
      // cmodif: call gufld(vout,f) changed into:

      GetField(vout, f);

      // *
      // *             start of integration
      // *
      x = vout[0];
      y = vout[1];
      z = vout[2];
      a = vout[3];
      b = vout[4];
      c = vout[5];

      h2       = khalf * h;
      h4       = khalf * h2;
      ph       = pinv * h;
      ph2      = khalf * ph;
      secxs[0] = (b * f[2] - c * f[1]) * ph2;
      secys[0] = (c * f[0] - a * f[2]) * ph2;
      seczs[0] = (a * f[1] - b * f[0]) * ph2;
      ang2     = (secxs[0] * secxs[0] + secys[0] * secys[0] + seczs[0] * seczs[0]);
      if (ang2 > kpisqua) break;

      dxt = h2 * a + h4 * secxs[0];
      dyt = h2 * b + h4 * secys[0];
      dzt = h2 * c + h4 * seczs[0];
      xt  = x + dxt;
      yt  = y + dyt;
      zt  = z + dzt;
      // *
      // *              second intermediate point
      // *

      est = TMath::Abs(dxt) + TMath::Abs(dyt) + TMath::Abs(dzt);
      if (est > h) {
         if (ncut++ > maxcut) break;
         h *= khalf;
         continue;
      }

      xyzt[0] = xt;
      xyzt[1] = yt;
      xyzt[2] = zt;

      // cmodif: call gufld(xyzt,f) changed into:
      GetField(xyzt, f);

      at = a + secxs[0];
      bt = b + secys[0];
      ct = c + seczs[0];

      secxs[1] = (bt * f[2] - ct * f[1]) * ph2;
      secys[1] = (ct * f[0] - at * f[2]) * ph2;
      seczs[1] = (at * f[1] - bt * f[0]) * ph2;
      at       = a + secxs[1];
      bt       = b + secys[1];
      ct       = c + seczs[1];
      secxs[2] = (bt * f[2] - ct * f[1]) * ph2;
      secys[2] = (ct * f[0] - at * f[2]) * ph2;
      seczs[2] = (at * f[1] - bt * f[0]) * ph2;
      dxt      = h * (a + secxs[2]);
      dyt      = h * (b + secys[2]);
      dzt      = h * (c + seczs[2]);
      xt       = x + dxt;
      yt       = y + dyt;
      zt       = z + dzt;
      at       = a + 2. * secxs[2];
      bt       = b + 2. * secys[2];
      ct       = c + 2. * seczs[2];

      est = TMath::Abs(dxt) + TMath::Abs(dyt) + TMath::Abs(dzt);
      if (est > 2. * TMath::Abs(h)) {
         if (ncut++ > maxcut) break;
         h *= khalf;
         continue;
      }

      xyzt[0] = xt;
      xyzt[1] = yt;
      xyzt[2] = zt;

      // cmodif: call gufld(xyzt,f) changed into:
      GetField(xyzt, f);

      z = z + (c + (seczs[0] + seczs[1] + seczs[2]) * kthird) * h;
      y = y + (b + (secys[0] + secys[1] + secys[2]) * kthird) * h;
      x = x + (a + (secxs[0] + secxs[1] + secxs[2]) * kthird) * h;

      secxs[3] = (bt * f[2] - ct * f[1]) * ph2;
      secys[3] = (ct * f[0] - at * f[2]) * ph2;
      seczs[3] = (at * f[1] - bt * f[0]) * ph2;
      a        = a + (secxs[0] + secxs[3] + 2. * (secxs[1] + secxs[2])) * kthird;
      b        = b + (secys[0] + secys[3] + 2. * (secys[1] + secys[2])) * kthird;
      c        = c + (seczs[0] + seczs[3] + 2. * (seczs[1] + seczs[2])) * kthird;

      est = TMath::Abs(secxs[0] + secxs[3] - (secxs[1] + secxs[2])) +
            TMath::Abs(secys[0] + secys[3] - (secys[1] + secys[2])) +
            TMath::Abs(seczs[0] + seczs[3] - (seczs[1] + seczs[2]));

      if (est > kdlt && TMath::Abs(h) > 1.e-4) {
         if (ncut++ > maxcut) break;
         h *= khalf;
         continue;
      }

      ncut = 0;
      // *               if too many iterations, go to helix
      if (iter++ > maxit) break;

      tl += h;
      if (est < kdlt32) h *= 2.;
      cba     = 1. / TMath::Sqrt(a * a + b * b + c * c);
      vout[0] = x;
      vout[1] = y;
      vout[2] = z;
      vout[3] = cba * a;
      vout[4] = cba * b;
      vout[5] = cba * c;
      rest    = step - tl;
      if (step < 0.) rest = -rest;
      if (rest < 1.e-5 * TMath::Abs(step)) return;

   } while (1);

   // angle too big, use helix

   f1  = f[0];
   f2  = f[1];
   f3  = f[2];
   f4  = TMath::Sqrt(f1 * f1 + f2 * f2 + f3 * f3);
   rho = -f4 * pinv;
   tet = rho * step;

   hnorm = 1. / f4;
   f1    = f1 * hnorm;
   f2    = f2 * hnorm;
   f3    = f3 * hnorm;

   hxp[0] = f2 * vect[kipz] - f3 * vect[kipy];
   hxp[1] = f3 * vect[kipx] - f1 * vect[kipz];
   hxp[2] = f1 * vect[kipy] - f2 * vect[kipx];

   hp = f1 * vect[kipx] + f2 * vect[kipy] + f3 * vect[kipz];

   rho1 = 1. / rho;
   sint = TMath::Sin(tet);
   cost = 2. * TMath::Sin(khalf * tet) * TMath::Sin(khalf * tet);

   g1 = sint * rho1;
   g2 = cost * rho1;
   g3 = (tet - sint) * hp * rho1;
   g4 = -cost;
   g5 = sint;
   g6 = cost * hp;

   vout[kix] = vect[kix] + g1 * vect[kipx] + g2 * hxp[0] + g3 * f1;
   vout[kiy] = vect[kiy] + g1 * vect[kipy] + g2 * hxp[1] + g3 * f2;
   vout[kiz] = vect[kiz] + g1 * vect[kipz] + g2 * hxp[2] + g3 * f3;

   vout[kipx] = vect[kipx] + g4 * vect[kipx] + g5 * hxp[0] + g6 * f1;
   vout[kipy] = vect[kipy] + g4 * vect[kipy] + g5 * hxp[1] + g6 * f2;
   vout[kipz] = vect[kipz] + g4 * vect[kipz] + g5 * hxp[2] + g6 * f3;
}

//__________________________________________________________________________
void MpdKalmanFilter::GetField(Double_t *position, Double_t *field)
{
   /// Just simple interface
   // if (fMagField) fMagField->GetFieldValue(position,field);
   if (fMagField)
      ((FairField *)((MpdMultiField *)fMagField)->GetFieldList()->UncheckedAt(0))->Field(position, field);
   else {
      cout << "F-MpdKalmanFilter::GetField: fMagField = 0x0" << endl;
      exit(-1);
   }
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Scattering(MpdKalmanTrack *track, Double_t dx, TString mass2, Int_t charge)
{
   /// Compute multiple scattering angle. dx - material thickness (in X0 - for normal incidence)

   // const Double_t piMass2 = 0.13957 * 0.13957; // pion hypothesis
   Double_t piMass2 = mass2.Atof();
   Double_t ph      = track->GetParamNew(2);
   Double_t cosTh   = TMath::Cos(track->GetParamNew(3));
   Double_t sinTh   = TMath::Sin(track->GetParamNew(3));
   Double_t ppp     = charge / track->GetParamNew(4) / cosTh;
   ppp              = TMath::Abs(ppp);
   ppp              = TMath::Max(ppp, 0.1);
   ppp              = TMath::Min(ppp, 3.);
   TVector3       p(0, 0, ppp);
   TLorentzVector pe(p, TMath::Sqrt(ppp * ppp + piMass2));
   Double_t       beta = pe.Beta(), phPos = 0, scaleL = 0;
   if (track->GetType() == MpdKalmanTrack::kBarrel) {
      // Path length correction - barrel
      phPos  = track->GetParamNew(0) / track->GetPosNew();
      scaleL = TMath::Cos(phPos - ph) * cosTh;
   } else {
      phPos  = TMath::ATan2(track->GetParamNew(1), track->GetParamNew(0));
      scaleL = sinTh; // path length correction - endcap
   }
   Double_t dxx0  = dx / TMath::Abs(scaleL);
   Double_t angle = 0.0136 / beta / ppp * charge * (1. + 0.038 * TMath::Log(dxx0));
   // cout << angle << " " << beta << " " << scaleL << " " << ppp << endl;
   // exit(0);
   return angle * angle * dxx0;
}

//__________________________________________________________________________
Double_t MpdKalmanFilter::Scattering(MpdKalmanTrack *track, Double_t x0, Double_t step, TString mass2, Int_t charge)
{
   /// Compute multiple scattering angle. x0 - rad. length, step - track path length

   // const Double_t piMass2 = 0.13957 * 0.13957; // pion hypothesis
   Double_t piMass2 = mass2.Atof();
   // Double_t ph = track->GetParamNew(2);
   Double_t cosTh = TMath::Cos(track->GetParamNew(3));
   // Double_t sinTh = TMath::Sin (track->GetParamNew(3));
   Double_t ppp = charge / track->GetParamNew(4) / cosTh;
   ppp          = TMath::Abs(ppp);
   ppp          = TMath::Max(ppp, 0.1);
   ppp          = TMath::Min(ppp, 3.);
   TVector3       p(0, 0, ppp);
   TLorentzVector pe(p, TMath::Sqrt(ppp * ppp + piMass2));
   Double_t       beta = pe.Beta();
   /*
   Double_t beta = pe.Beta(), phPos = 0, scaleL = 0;
   if (track->GetType() == MpdKalmanTrack::kFixedR) {
     // Path length correction - fixed R
     phPos = track->GetParamNew(0) / track->GetPosNew();
     scaleL = TMath::Cos(phPos-ph) * cosTh;
   } else {
     phPos = TMath::ATan2 (track->GetParamNew(1), track->GetParamNew(0));
     scaleL = sinTh; // path length correction - fixed Z
   }
   Double_t dxx0 = dx / TMath::Abs(scaleL);
   */
   Double_t dxx0  = step / x0;
   Double_t angle = 0.0136 / beta / ppp * charge * (1. + 0.038 * TMath::Log(dxx0));
   // cout << angle << " " << beta << " " << scaleL << " " << ppp << endl;
   // exit(0);
   return angle * angle * dxx0;
}

//__________________________________________________________________________
/*
TVector3 TpcLheKalmanFilter::GetTrackMom(const TpcLheHit *hit)
{
  /// Returns track momentum for given hit

  TpcPoint* point = (TpcPoint*) fTpcPoints->UncheckedAt(hit->GetUniqueID());
  TVector3 p(point->GetPx(),point->GetPy(),point->GetPz());
  return p;
}
*/

//__________________________________________________________________________
void MpdKalmanFilter::ExtrapOneStepHelix(Double_t charge, Double_t step, Double_t *vect, Double_t *vout)
{
   //    ******************************************************************
   //    *                                                                *
   //    *  Performs the tracking of one step in a magnetic field         *
   //    *  The trajectory is assumed to be a helix in a constant field   *
   //    *  taken at the mid point of the step.                           *
   //    *  Parameters:                                                   *
   //    *   input                                                        *
   //    *     STEP =arc length of the step asked                         *
   //    *     VECT =input vector (position,direction cos and momentum)   *
   //    *     CHARGE=  electric charge of the particle                   *
   //    *   output                                                       *
   //    *     VOUT = same as VECT after completion of the step           *
   //    *                                                                *
   //    *    ==>Called by : <USER>, GUSWIM                               *
   //    *       Author    m.hansroul  *********                          *
   //    *       modified  s.egli, s.v.levonian                           *
   //    *       modified  v.perevoztchikov
   //    *                                                                *
   //    ******************************************************************
   //

   // modif: everything in double precision

   Double_t xyz[3], h[4], hxp[3];
   Double_t h2xy, hp, rho, tet;
   Double_t sint, sintt, tsint, cos1t;
   Double_t f1, f2, f3, f4, f5, f6;

   const Int_t kix  = 0;
   const Int_t kiy  = 1;
   const Int_t kiz  = 2;
   const Int_t kipx = 3;
   const Int_t kipy = 4;
   const Int_t kipz = 5;
   const Int_t kipp = 6;

   const Double_t kec = 2.9979251e-4;
   //
   //    ------------------------------------------------------------------
   //
   //       units are kgauss,centimeters,gev/c
   //
   vout[kipp] = vect[kipp];
   if (TMath::Abs(charge) < 0.00001) {
      for (Int_t i = 0; i < 3; i++) {
         vout[i]     = vect[i] + step * vect[i + 3];
         vout[i + 3] = vect[i + 3];
      }
      return;
   }
   xyz[0] = vect[kix] + 0.5 * step * vect[kipx];
   xyz[1] = vect[kiy] + 0.5 * step * vect[kipy];
   xyz[2] = vect[kiz] + 0.5 * step * vect[kipz];

   // cmodif: call gufld (xyz, h) changed into:
   GetField(xyz, h);

   h2xy = h[0] * h[0] + h[1] * h[1];
   h[3] = h[2] * h[2] + h2xy;
   if (h[3] < 1.e-12) {
      for (Int_t i = 0; i < 3; i++) {
         vout[i]     = vect[i] + step * vect[i + 3];
         vout[i + 3] = vect[i + 3];
      }
      return;
   }
   if (h2xy < 1.e-12 * h[3]) {
      ExtrapOneStepHelix3(charge * h[2], step, vect, vout);
      return;
   }
   h[3] = TMath::Sqrt(h[3]);
   h[0] /= h[3];
   h[1] /= h[3];
   h[2] /= h[3];
   h[3] *= kec;

   hxp[0] = h[1] * vect[kipz] - h[2] * vect[kipy];
   hxp[1] = h[2] * vect[kipx] - h[0] * vect[kipz];
   hxp[2] = h[0] * vect[kipy] - h[1] * vect[kipx];

   hp = h[0] * vect[kipx] + h[1] * vect[kipy] + h[2] * vect[kipz];

   rho = -charge * h[3] / vect[kipp];
   tet = rho * step;

   if (TMath::Abs(tet) > 0.15) {
      sint  = TMath::Sin(tet);
      sintt = (sint / tet);
      tsint = (tet - sint) / tet;
      cos1t = 2. * (TMath::Sin(0.5 * tet)) * (TMath::Sin(0.5 * tet)) / tet;
   } else {
      tsint = tet * tet / 36.;
      sintt = (1. - tsint);
      sint  = tet * sintt;
      cos1t = 0.5 * tet;
   }

   f1 = step * sintt;
   f2 = step * cos1t;
   f3 = step * tsint * hp;
   f4 = -tet * cos1t;
   f5 = sint;
   f6 = tet * cos1t * hp;

   vout[kix] = vect[kix] + f1 * vect[kipx] + f2 * hxp[0] + f3 * h[0];
   vout[kiy] = vect[kiy] + f1 * vect[kipy] + f2 * hxp[1] + f3 * h[1];
   vout[kiz] = vect[kiz] + f1 * vect[kipz] + f2 * hxp[2] + f3 * h[2];

   vout[kipx] = vect[kipx] + f4 * vect[kipx] + f5 * hxp[0] + f6 * h[0];
   vout[kipy] = vect[kipy] + f4 * vect[kipy] + f5 * hxp[1] + f6 * h[1];
   vout[kipz] = vect[kipz] + f4 * vect[kipz] + f5 * hxp[2] + f6 * h[2];

   return;
}

//__________________________________________________________________________
void MpdKalmanFilter::ExtrapOneStepHelix3(Double_t field, Double_t step, Double_t *vect, Double_t *vout) const
{
   //
   //     ******************************************************************
   //     *                                                                *
   //     *       Tracking routine in a constant field oriented            *
   //     *       along axis 3                                             *
   //     *       Tracking is performed with a conventional                *
   //     *       helix step method                                        *
   //     *                                                                *
   //     *    ==>Called by : <USER>, GUSWIM                               *
   //     *       Authors    R.Brun, M.Hansroul  *********                 *
   //     *       Rewritten  V.Perevoztchikov
   //     *                                                                *
   //     ******************************************************************
   //

   Double_t hxp[3];
   Double_t h4, hp, rho, tet;
   Double_t sint, sintt, tsint, cos1t;
   Double_t f1, f2, f3, f4, f5, f6;

   const Int_t kix  = 0;
   const Int_t kiy  = 1;
   const Int_t kiz  = 2;
   const Int_t kipx = 3;
   const Int_t kipy = 4;
   const Int_t kipz = 5;
   const Int_t kipp = 6;

   const Double_t kec = 2.9979251e-4;

   //
   //     ------------------------------------------------------------------
   //
   //       units are kgauss,centimeters,gev/c
   //
   vout[kipp] = vect[kipp];
   h4         = field * kec;

   hxp[0] = -vect[kipy];
   hxp[1] = +vect[kipx];

   hp = vect[kipz];

   rho = -h4 / vect[kipp];
   tet = rho * step;
   if (TMath::Abs(tet) > 0.15) {
      sint  = TMath::Sin(tet);
      sintt = (sint / tet);
      tsint = (tet - sint) / tet;
      cos1t = 2. * TMath::Sin(0.5 * tet) * TMath::Sin(0.5 * tet) / tet;
   } else {
      tsint = tet * tet / 36.;
      sintt = (1. - tsint);
      sint  = tet * sintt;
      cos1t = 0.5 * tet;
   }

   f1 = step * sintt;
   f2 = step * cos1t;
   f3 = step * tsint * hp;
   f4 = -tet * cos1t;
   f5 = sint;
   f6 = tet * cos1t * hp;

   vout[kix] = vect[kix] + f1 * vect[kipx] + f2 * hxp[0];
   vout[kiy] = vect[kiy] + f1 * vect[kipy] + f2 * hxp[1];
   vout[kiz] = vect[kiz] + f1 * vect[kipz] + f3;

   vout[kipx] = vect[kipx] + f4 * vect[kipx] + f5 * hxp[0];
   vout[kipy] = vect[kipy] + f4 * vect[kipy] + f5 * hxp[1];
   vout[kipz] = vect[kipz] + f4 * vect[kipz] + f6;

   return;
}

//__________________________________________________________________________
void MpdKalmanFilter::MnvertLocal(Double_t *a, Int_t l, Int_t, Int_t n, Int_t &ifail)
{
   ///*-*-*-*-*-*-*-*-*-*-*-*Inverts a symmetric matrix*-*-*-*-*-*-*-*-*-*-*-*-*
   ///*-*                    ==========================
   ///*-*        inverts a symmetric matrix.   matrix is first scaled to
   ///*-*        have all ones on the diagonal (equivalent to change of units)
   ///*-*        but no pivoting is done since matrix is positive-definite.
   ///*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   // taken from TMinuit package of Root (l>=n)
   // fVERTs, fVERTq and fVERTpp changed to localVERTs, localVERTq and localVERTpp
   //  Double_t localVERTs[n], localVERTq[n], localVERTpp[n];
   Double_t *localVERTs  = new Double_t[n];
   Double_t *localVERTq  = new Double_t[n];
   Double_t *localVERTpp = new Double_t[n];
   // fMaxint changed to localMaxint
   Int_t localMaxint = n;

   /* System generated locals */
   Int_t aOffset;

   /* Local variables */
   Double_t si;
   Int_t    i, j, k, kp1, km1;

   /* Parameter adjustments */
   aOffset = l + 1;
   a -= aOffset;

   /* Function Body */
   ifail = 0;
   if (n < 1) goto L100;
   if (n > localMaxint) goto L100;
   //*-*-                  scale matrix by sqrt of diag elements
   for (i = 1; i <= n; ++i) {
      si = a[i + i * l];
      if (si <= 0) goto L100;
      localVERTs[i - 1] = 1 / TMath::Sqrt(si);
   }
   for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
         a[i + j * l] = a[i + j * l] * localVERTs[i - 1] * localVERTs[j - 1];
      }
   }
   //*-*-                                       . . . start main loop . . . .
   for (i = 1; i <= n; ++i) {
      k = i;
      //*-*-                  preparation for elimination step1
      if (a[k + k * l] != 0)
         localVERTq[k - 1] = 1 / a[k + k * l];
      else
         goto L100;
      localVERTpp[k - 1] = 1;
      a[k + k * l]       = 0;
      kp1                = k + 1;
      km1                = k - 1;
      if (km1 < 0)
         goto L100;
      else if (km1 == 0)
         goto L50;
      else
         goto L40;
   L40:
      for (j = 1; j <= km1; ++j) {
         localVERTpp[j - 1] = a[j + k * l];
         localVERTq[j - 1]  = a[j + k * l] * localVERTq[k - 1];
         a[j + k * l]       = 0;
      }
   L50:
      if (k - n < 0)
         goto L51;
      else if (k - n == 0)
         goto L60;
      else
         goto L100;
   L51:
      for (j = kp1; j <= n; ++j) {
         localVERTpp[j - 1] = a[k + j * l];
         localVERTq[j - 1]  = -a[k + j * l] * localVERTq[k - 1];
         a[k + j * l]       = 0;
      }
   //*-*-                  elimination proper
   L60:
      for (j = 1; j <= n; ++j) {
         for (k = j; k <= n; ++k) {
            a[j + k * l] += localVERTpp[j - 1] * localVERTq[k - 1];
         }
      }
   }
   //*-*-                  elements of left diagonal and unscaling
   for (j = 1; j <= n; ++j) {
      for (k = 1; k <= j; ++k) {
         a[k + j * l] = a[k + j * l] * localVERTs[k - 1] * localVERTs[j - 1];
         a[j + k * l] = a[k + j * l];
      }
   }
   delete[] localVERTs;
   delete[] localVERTq;
   delete[] localVERTpp;
   return;
//*-*-                  failure return
L100:
   delete[] localVERTs;
   delete[] localVERTq;
   delete[] localVERTpp;
   ifail = 1;
} /* mnvertLocal */
//__________________________________________________________________________

Double_t *MpdKalmanFilter::ExtrapOneStep(MpdKalmanTrack *track, Double_t step, Int_t flag)
{
   // Extrapolate one step: if flag!=0 then use the cached parameters

   Double_t vect[7], charge;

   if (flag == 0) {
      if (track->GetType() == MpdKalmanTrack::kBarrel)
         SetGeantParamB(track, vect, 1);
      else
         SetGeantParamE(track, vect, 1);
      charge = fVectorG[7] = -TMath::Sign(1.0, track->GetParam(4));
   } else {
      for (Int_t j = 0; j < 7; ++j) vect[j] = fVectorG[j];
      charge = fVectorG[7];
   }

   ExtrapOneStepRungekutta(charge, step, vect, fVectorG);
   return fVectorG;
}
//__________________________________________________________________________

ClassImp(MpdKalmanFilter)
