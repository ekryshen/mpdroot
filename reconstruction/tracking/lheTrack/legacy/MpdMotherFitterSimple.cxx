/// \class MpdMotherFitterSimple
///
/// Simple mother fitter for the MPD detector (using helices and straight lines)
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdMotherFitterSimple.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdHelix.h"

#include "FairField.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"
#include "FairTask.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TMinuit.h>
#include <TVector3.h>

#include <iostream>
using std::cout;
using std::endl;

MpdMotherFitterSimple *MpdMotherFitterSimple::fgMF = 0x0;

//__________________________________________________________________________
MpdMotherFitterSimple::MpdMotherFitterSimple() : FairTask()
{
   /// Default constructor
   for (Int_t i = 0; i < 3; ++i) fTrC[i] = NULL;
}

//__________________________________________________________________________
MpdMotherFitterSimple::MpdMotherFitterSimple(const char *name, const char *title) : FairTask(name)
{
   /// Constructor
   for (Int_t i = 0; i < 3; ++i) fTrC[i] = NULL;
   fgMF = this;
}

//__________________________________________________________________________
MpdMotherFitterSimple *MpdMotherFitterSimple::Instance()
{
   /// Get pointer to the mother fitter singleton object
   if (!fgMF) {
      fgMF = new MpdMotherFitterSimple;
      fgMF->Init();
      // automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgMF;
}

//__________________________________________________________________________
MpdMotherFitterSimple *MpdMotherFitterSimple::Instance(const char *name, const char *title)
{
   /// Get pointer to the mother fitter singleton object
   if (!fgMF) {
      fgMF = new MpdMotherFitterSimple(name, title);
      fgMF->Init();
      // automatic destruction is forced
      std::atexit(DestroyInstance);
   }
   return fgMF;
}

//__________________________________________________________________________
MpdMotherFitterSimple::~MpdMotherFitterSimple()
{
   /// Destructor
   // FairRootManager *manager = FairRootManager::Instance();
   // manager->Write();
   for (Int_t i = 0; i < 3; ++i) delete fTrC[i];
}

//__________________________________________________________________________
InitStatus MpdMotherFitterSimple::Init()
{

   cout << "InitStatus MpdMotherFitterSimple::Init\n\n";
   Double_t   bZ       = 5.;
   FairField *magField = FairRunAna::Instance()->GetField();
   if (!magField || TMath::Abs(magField->GetBz(0, 0, 0)) < 0.01) {
      cout << " -I- MpdMotherFitterSimple::Init - Using the default constant magnetic field Bz = 5 kG " << endl;
   } else
      bZ = TMath::Abs(magField->GetBz(0, 0, 0));
   fieldConst = 0.3 * 0.01 * bZ / 10;
   for (Int_t i = 0; i < 3; ++i) fTrC[i] = new MpdHelix(0.0, 0.0, 0.0, TVector3(0, 0, 0));

   if (!gMinuit) new TMinuit(2);
   return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdMotherFitterSimple::ReInit()
{
   // fMagField = FairRunAna::Instance()->GetField(); // !!! interim solution !!!
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdMotherFitterSimple::Reset()
{
   ///
   // cout << " MpdMotherFitterSimple::Reset  " << endl;
}

//__________________________________________________________________________
void MpdMotherFitterSimple::Register()
{
   ///
   // FairRootManager::Instance()->
   // Register("TpcLheKalmanTrack","TpcLheKalmanTrack", fTracks, kFALSE);
}

//__________________________________________________________________________
void MpdMotherFitterSimple::Finish()
{
   ///
}

//__________________________________________________________________________
void MpdMotherFitterSimple::Exec(Option_t *option)
{
   ///
}

//__________________________________________________________________________
TLorentzVector MpdMotherFitterSimple::BuildMother(vector<Double_t> &mass, vector<MpdKalmanTrack *> &vDaughtTr,
                                                  vector<TLorentzVector *> &vDaughtLv)
{
   // Build mother from tracks and/or Lorentz vectors

   Int_t           nDaughtTr = vDaughtTr.size();
   Int_t           nDaughtLv = vDaughtLv.size();
   TVector3        vtx;
   TLorentzVector  mother;
   MpdKalmanTrack *tr = NULL;

   if (nDaughtLv == 0) {
      // Only tracks
      FindVertex(vDaughtTr);
      for (Int_t i = 0; i < nDaughtTr; ++i) {
         // tr = vDaughtTr[i];
         Double_t       pt  = Pt(fTrC[i]);
         Double_t       phi = Phi(fTrC[i]);
         Double_t       dip = fTrC[i]->dipAngle();
         TLorentzVector part;
         part.SetXYZM(pt * TMath::Cos(phi), pt * TMath::Sin(phi), pt * TMath::Tan(dip), mass[i]);
         mother += part;
      }
      return mother;
   }

   // At present this part works only for one neutral (LorentzVector) and one charged track
   // (to reconstruct cascades)
   // Create MpdHelix
   tr           = vDaughtTr[0];
   Double_t r   = tr->GetPosNew();
   Double_t phi = tr->GetParam(0) / r;
   Double_t x   = r * TMath::Cos(phi);
   Double_t y   = r * TMath::Sin(phi);
   Double_t dip = tr->GetParam(3);
   Double_t cur = TMath::Abs(tr->GetParam(4)) * fieldConst;
   TVector3 o(x, y, tr->GetParam(1));
   Int_t    h = (Int_t)TMath::Sign(1.1, tr->GetParam(4));
   fTrC[0]->setParameters(cur, dip, tr->GetParam(2) - TMath::PiOver2() * h, o, h);

   /*
   Double_t s = trC.pathLength(primVert);
   TVector3 pca = trC.at(s);
   pca -= primVert;
   // !!! Cut on K-
   if (tr->GetChi2Vertex() < gChi2K || pca.Mag() < gDcaK) continue;
   */

   // Find decay vertex of a neutral and charged track
   TLorentzVector lambda = *vDaughtLv.front();
   fMomN                 = lambda.Vect();
   fMomN *= (1. / fMomN.Mag());
   // vtxN = vtxL0; ///
   // void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t))
   gMinuit->SetFCN(FcnDist);

   Double_t arglist[10];
   Int_t    ierflg = 0;
   arglist[0]      = 1;
   gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = -1; // 1; //-1;
   gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

   Double_t        vstart[2] = {-0.1, -0.1};
   static Double_t step[2]   = {0.1, 0.1};
   gMinuit->mnparm(0, "lengN", vstart[0], step[0], 0, 0, ierflg);
   gMinuit->mnparm(1, "lengC", vstart[1], step[1], 0, 0, ierflg);

   // Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

   // Print results
   //*
   Double_t amin, edm, errdef;
   Int_t    nvpar, nparx, icstat;
   gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
   fDist = amin;
   // gMinuit->mnprin(3,amin);
   //*/
   //  !!! Cut on distance between daughters
   // if (amin > gDistLK) continue;

   Double_t lengC, lengN, parErr;
   gMinuit->GetParameter(0, lengN, parErr);
   gMinuit->GetParameter(1, lengC, parErr);
   TVector3 posC  = fTrC[0]->at(lengC);
   TVector3 posN  = TVector3(fVtx) + fMomN * lengN;
   TVector3 posOm = (posN + posC) * 0.5;
   posOm.GetXYZ(fVtx);

   Double_t pt = TMath::Abs(1. / tr->GetParam(4));
   Double_t s  = fTrC[0]->pathLength(posC); // PCA with decay vertex
   fTrC[0]->moveOrigin(s);
   phi = Phi(fTrC[0]);
   dip = tr->GetParam(3);
   TLorentzVector part;
   part.SetXYZM(pt * TMath::Cos(phi), pt * TMath::Sin(phi), pt * TMath::Tan(dip), mass[0]);

   mother = lambda;
   mother += part;
   return mother;
}

//____________________________________________________________________________
void FcnDist(Int_t &npar, Double_t *gin, Double_t &rf, Double_t *par, Int_t iflag)
{
   // Compute distance between straight line and helix

   TVector3 mom = MpdMotherFitterSimple::Instance()->GetMomN();
   mom *= par[0];
   TVector3 posN = TVector3(MpdMotherFitterSimple::Instance()->GetVertex());
   posN += mom;

   TVector3 posC = MpdMotherFitterSimple::Instance()->GetHelix(0)->at(par[1]);
   posC -= posN;
   rf = posC.Mag();
   // cout << par[0] << " " << par[1] << " " << f << endl;
}

//__________________________________________________________________________
void MpdMotherFitterSimple::FindVertex(vector<MpdKalmanTrack *> vDaught)
{
   /// Find decay vertex position as PCA of helices

   // Create MpdHelix's
   Int_t nDaught = vDaught.size();
   if (nDaught > 2) Fatal("MpdMotherFitterSimple::FindVertex", " Too many tracks !!!");
   TVector3 sum, p1, p2;
   for (Int_t i = 0; i < nDaught; ++i) {
      Double_t r   = vDaught[i]->GetPosNew();
      Double_t phi = vDaught[i]->GetParam(0) / r;
      Double_t x   = r * TMath::Cos(phi);
      Double_t y   = r * TMath::Sin(phi);

      Double_t dip = vDaught[i]->GetParam(3);
      Double_t cur = TMath::Abs(vDaught[i]->GetParam(4)) * fieldConst;
      TVector3 o(x, y, vDaught[i]->GetParam(1));

      Int_t h = (Int_t)TMath::Sign(1.1, vDaught[i]->GetParam(4));
      // helix[itr++] = MpdHelix(cur, dip, phi, o, h);
      // setParameters(double c, double dip, double phase, const TVector3 o, int h);
      fTrC[i]->setParameters(cur, dip, vDaught[i]->GetParam(2) - TMath::PiOver2() * h, o, h);
      // cout << helix[i]->xcenter() << " " << helix[i]->ycenter() << endl;
      // sum += o;
   }

   Int_t    ncombs = 0, nD1 = nDaught - 1;
   Double_t pathMax = 0.0;
   for (Int_t i = 0; i < nD1; ++i) {
      for (Int_t j = 1; j < nDaught; ++j) {
         pair<Double_t, Double_t> paths = fTrC[i]->pathLengths(*fTrC[j]);
         // cout << " Intersection: " << paths.first << " " << paths.second << endl;
         p1 = fTrC[i]->at(paths.first);
         p2 = fTrC[j]->at(paths.second);
         sum += (p1 + p2);
         ncombs += 2;
         pathMax = TMath::Max(pathMax, TMath::Abs(paths.first));
         pathMax = TMath::Max(pathMax, TMath::Abs(paths.second));
      }
   }

   if (pathMax > 80.)
      sum.SetXYZ(0, 0, 0); // check constrains
   else
      sum *= (1. / ncombs);
   sum.GetXYZ(fVtx);
   fDist = (p1 - p2).Mag();
   // Move helices to the vertex
   for (Int_t i = 0; i < nDaught; ++i) {
      Double_t s = fTrC[i]->pathLength(sum); // PCA with decay vertex
      fTrC[i]->moveOrigin(s);
   }
}

//__________________________________________________________________________
Double_t MpdMotherFitterSimple::Phi(MpdHelix *helix) const
{
   // Compute Phi of the track

   return helix->phase() + TMath::PiOver2() * helix->h(); // check !!!
}

//__________________________________________________________________________

ClassImp(MpdMotherFitterSimple);
