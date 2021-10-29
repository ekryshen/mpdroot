#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include "TpcPoint.h"
#include "MpdHelix.h"
#include "MpdItsKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdKfV0Fitter.h"
#include "MpdTrackFinderIts.h"
#include "MpdVertex.h"

// CBM includes
#include "FairMCPoint.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"

// ROOT includes
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
//#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <Riostream.h>

#include <set>
#include <map>
#include <vector>

#endif

Int_t pdgCodePos= 2212; // proton
Int_t pdgCodeNeg= -211; // pi-
Int_t pdgCodeL0 = 3122; // Lambda0 (1.11568)
Int_t pdgCodeK0 = 310; // K0s
Int_t pdgCodeOm = 3334; // Omega- (1.67245 GeV/c)
Int_t pdgCodeXi = 3312; // Xi- (1.3213 GeV/c)
Int_t pdgCodeKm = -321; // K-
Int_t nITS = 0, idPlus = -1;
Int_t *lays = 0x0;
MpdHelix trC;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks;
TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
FILE *lun = fopen("event.dat","w");
FILE *lun1 = 0x0; //fopen("ids.dat","w");
const Double_t gChi2pos = 26.;
const Double_t gChi2neg = 26.;
const Double_t gImpPos = 0.1;
const Double_t gImpNeg = 0.2;
const Double_t gDistV0 = 0.12;
const Double_t gDecayV0 = 1.5; //1.0;
const Double_t gChi2V0 = 10.0;
const Double_t gDcaL0 = 0.15; //0.1;
const Double_t gDcaK = 0.3;
const Double_t gChi2K = 100; //50;
const Double_t gDecayOm = 1.0;
const Double_t gDistLK = 1.15; //0.2;
const Double_t gDcaOm = 0.1; //0.15;

typedef MpdItsKalmanTrack MpdTrack;

//__________________________________________________________________________
void CalcMass(Int_t decayInd[2], Double_t rads[2], Double_t mass[2], Int_t flag)
{
  // Compute invariant mass of V0

  TLorentzVector lorVec[2];
  for (Int_t i = 0; i < 2; ++i) {
    MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(decayInd[i]);
    MpdTrack track = *tr;
    track.SetPos(track.GetPosNew());
    track.SetParamNew(*track.GetParam());
    MpdKalmanHit hit;
    hit.SetType(MpdKalmanHit::kFixedR);
    hit.SetDist(rads[i]);
    MpdKalmanFilter::Instance()->PropagateParamR(&track, &hit, kFALSE);
    Double_t pt = TMath::Abs (1. / track.GetParamNew(4));
    Double_t phi = track.GetParamNew(2);
    Double_t dip = track.GetParamNew(3);
    lorVec[i].SetXYZM(pt*TMath::Cos(phi),pt*TMath::Sin(phi),pt*TMath::Tan(dip),mass[i]);
    if (mass[i] < 0.2) {
      FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(tr->GetTrackID());
      TVector3 mom;
      mcTr->GetMomentum(mom);
      Double_t dphi = MpdKalmanFilter::Instance()->Proxim(mom.Phi(),lorVec[i].Phi()) - mom.Phi();
    }
  }
  TLorentzVector lorInv = lorVec[0] + lorVec[1];
  if (flag) {
    // True vertex
    cout << " V0 mass = " << lorInv.Mag() << endl;
    ((TH1D*)gROOT->FindObjectAny("hInvMass0"))->Fill(lorInv.Mag());
    //if (lun) fprintf(lun,"%f ",lorInv.Mag());
  } else ((TH1D*)gROOT->FindObjectAny("hInvMassBkg"))->Fill(lorInv.Mag());
}

//__________________________________________________________________________
TLorentzVector CalcMass(MpdTrack *trs[2], Double_t mass[2], TVector3 vtx, MpdTrack trVtx[2])
{
  // Compute invariant mass of V0

  TLorentzVector lorVec[2];
  for (Int_t i = 0; i < 2; ++i) {
    MpdTrack *tr = trs[i];
    Double_t pt = tr->Pt();
    Double_t phi = tr->Phi();
    Double_t dip = tr->GetParam(3);
    lorVec[i].SetXYZM(pt*TMath::Cos(phi),pt*TMath::Sin(phi),pt*TMath::Tan(dip),mass[i]);
  }
  TLorentzVector lorInv = lorVec[0] + lorVec[1];
  //cout << " V0 mass = " << lorInv.Mag() << endl;
  ((TH1D*)gROOT->FindObjectAny("hInvMass1"))->Fill(lorInv.Mag());
  //if (lun) fprintf(lun,"%f ",lorInv.Mag());

  // Just propagate to DCA w.r.t. secondary vertex
  Double_t v[3];
  vtx.GetXYZ(v);
  for (Int_t i = 0; i < 2; ++i) {
    MpdTrack track = trVtx[i];
    track.SetPos(track.GetPosNew());
    track.SetParamNew(*track.GetParam());
    MpdKalmanFilter::Instance()->FindPca(&track,v);
    Double_t pt = TMath::Abs (1. / track.GetParamNew(4));
    Double_t phi = track.GetParamNew(2);
    Double_t dip = track.GetParamNew(3);
    lorVec[i].SetXYZM(pt*TMath::Cos(phi),pt*TMath::Sin(phi),pt*TMath::Tan(dip),mass[i]);
  }
  lorInv = lorVec[0] + lorVec[1];
  //cout << " V0 mass = " << lorInv.Mag() << endl;
  //((TH1D*)gROOT->FindObjectAny("hInvMass1"))->Fill(lorInv.Mag());
  //if (lun) fprintf(lun,"%f ",lorInv.Mag());
  return lorInv;
}

//__________________________________________________________________________
TLorentzVector CalcMass3(MpdTrack *tr, Double_t radC, Double_t mass, TLorentzVector lambda, Int_t signal)
{
  // Compute invariant mass of Xi- or Omega-

  TLorentzVector lorVec;

  MpdTrack track = *tr;
  track.SetPos(track.GetPosNew());
  track.SetParamNew(*track.GetParam());
  MpdKalmanHit hit;
  hit.SetType(MpdKalmanHit::kFixedR);
  hit.SetDist(radC);
  MpdKalmanFilter::Instance()->PropagateParamR(&track, &hit, kFALSE);
  Double_t pt = TMath::Abs (1. / track.GetParamNew(4));
  Double_t phi = track.GetParamNew(2);
  Double_t dip = track.GetParamNew(3);
  lorVec.SetXYZM(pt*TMath::Cos(phi),pt*TMath::Sin(phi),pt*TMath::Tan(dip),mass);

  TLorentzVector lorInv = lorVec + lambda;
  if (signal) {
    cout << " Xi- (Omega-) mass = " << lorInv.Mag() << " " << endl;
    if (mass < 0.200) ((TH1D*)gROOT->FindObjectAny("hInvMassXi"))->Fill(lorInv.Mag());
    else {
      ((TH1D*)gROOT->FindObjectAny("hInvMassOm"))->Fill(lorInv.Mag());
    }
  } else ((TH1D*)gROOT->FindObjectAny("hInvMassXiBkg"))->Fill(lorInv.Mag());

  //if (lun) fprintf(lun,"%f ",lorInv.Mag());
  return lorInv;
}

//____________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Compute distance between straight line and helix

  TVector3 mom = momN;
  mom *= par[0];
  TVector3 posN = vtxN;
  posN += mom;
  
  TVector3 posC = trC.at(par[1]);
  posC -= posN;
  f = posC.Mag();
  //cout << par[0] << " " << par[1] << " " << f << endl;
}

//__________________________________________________________________________
void CheckOmega(Int_t mothID, TLorentzVector lambda, TVector3 vtxL0, Int_t iev, Int_t jPi = -1)
{
  // Process Omega- or Xi-

  Double_t massLamb = lambda.Mag();
  //if (massLamb < 1.10713 || massLamb > 1.12423) return; // lambda mass +- 5*1.71 MeV
  if (massLamb < 1.11055 || massLamb > 1.12081) return; // lambda mass +- 3*1.71 MeV
 
  TVector3 dist = primVert - vtxL0;
  TVector3 momL0 = lambda.Vect();
  Double_t cosa = dist.Dot(momL0) / momL0.Mag() / dist.Mag();
  Double_t impL0 = dist.Mag() * TMath::Sqrt (1.-cosa*cosa);
  // !!! Cut on lambda
  if (impL0 < gDcaL0) return;

  Int_t gmothID = -1; 
  Int_t pdgCode = 0;
  if (mothID > -1) {
    gmothID = ((FairMCTrack*)mcTracks->UncheckedAt(mothID))->GetMotherId(); // grandmother
    if (gmothID == -1) {
      ((TH1D*)gROOT->FindObjectAny("hImpLamb"))->Fill(impL0);
      //return;
    } else {
      FairMCTrack *gmoth = (FairMCTrack*) mcTracks->UncheckedAt(gmothID); // grandmother
      pdgCode = gmoth->GetPdgCode();
    }

    if (pdgCode != pdgCodeOm && pdgCode != pdgCodeXi) {
      ((TH1D*)gROOT->FindObjectAny("hImpLamb"))->Fill(impL0);
      //return;
    }
  }

  for (Int_t j = 0; j < nITS; ++j) {
    if (j == jPi) continue; // the same pion as in Lambda combination
    MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
    //if (tr->GetNofIts() < 1) continue;
    if (tr->GetChi2() < -8) continue;
    Int_t id = tr->GetTrackID();
    //if (lays[id] < 0) continue; // primary
    //if (TMath::Abs(lays[id]) < 49) continue;
    if (TMath::Abs(lays[id]) < 41) continue;
    Int_t signal = 0;
    if (gmothID > -1 && ((FairMCTrack*) mcTracks->UncheckedAt(id))->GetMotherId() == gmothID) {
      cout << TDatabasePDG::Instance()->GetParticle(((FairMCTrack*) mcTracks->UncheckedAt(gmothID))->GetPdgCode())->GetName() << endl;
      signal = 1;
      //exit(0);
    }
    if (((FairMCTrack*) mcTracks->UncheckedAt(id))->GetPdgCode() != pdgCodeNeg) continue; // for Xi-

    // Create MpdHelix
    Double_t r = tr->GetPosNew();
    Double_t phi = tr->GetParam(0) / r;
    Double_t x = r * TMath::Cos(phi);
    Double_t y = r * TMath::Sin(phi);
    //cout << r << " " << x << " " << y << " " << tr->GetParam(2) << " " << phi << " " << tr->GetParam(4) << endl;
    Double_t dip = tr->GetParam(3);
    Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
    cur *= TMath::Abs (tr->GetParam(4));
    TVector3 o(x, y, tr->GetParam(1));
    Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
    //helix[itr++] = MpdHelix(cur, dip, phi, o, h);
    trC = MpdHelix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);

    Double_t s = trC.pathLength(primVert);
    TVector3 pca = trC.at(s);
    pca -= primVert;
    // !!! Cut on K-
    if (tr->GetChi2Vertex() < gChi2K || pca.Mag() < gDcaK) continue;

    // Find decay vertex of a neutral and charged track 
    momN = lambda.Vect();
    momN *= (1. / momN.Mag());
    //momN.Dump();
    vtxN = vtxL0;
    //vtxN.Dump();
    gMinuit->SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    arglist[0] = -1; //1; //-1;
    gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

    Double_t vstart[2] = {-0.1,-0.1};
    static Double_t step[2] = {0.1, 0.1};
    gMinuit->mnparm(0, "lengN", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "lengC", vstart[1], step[1], 0,0,ierflg);
    
    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    
    // Print results
    //*
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    //gMinuit->mnprin(3,amin);
    //*/
    // !!! Cut on distance between daughters
    if (amin > gDistLK) continue;

    Double_t lengC, lengN, parErr;
    gMinuit->GetParameter(0,lengN,parErr);
    gMinuit->GetParameter(1,lengC,parErr);
    TVector3 posC = trC.at(lengC);
    Double_t radC = posC.Pt();
    TVector3 posN = vtxN + momN * lengN;
    TVector3 posOm = (posN + posC) * 0.5;
    // !!! Cut on Omega decay length
    if (posOm.Mag() < gDecayOm) continue;

    Double_t massC = TDatabasePDG::Instance()->GetParticle(((FairMCTrack*) mcTracks->UncheckedAt(id))->GetPdgCode())->Mass();
    TLorentzVector omega = CalcMass3(tr, radC, massC, lambda, signal);

    dist = primVert - posOm;
    TVector3 momOm = omega.Vect();
    cosa = dist.Dot(momOm) / momOm.Mag() / dist.Mag();
    Double_t impOm = dist.Mag() * TMath::Sqrt (1.-cosa*cosa);
    // !!! Cut on omega impact
    if (impOm > gDcaOm) { 
      if (signal) ((TH1D*)gROOT->FindObjectAny("hInvMassXi"))->Fill(omega.Mag(),-1.);
      else ((TH1D*)gROOT->FindObjectAny("hInvMassXiBkg"))->Fill(omega.Mag(),-1.);
      //return;
      continue;
    }
    if (massC > -0.2) {
      //if (lun) fprintf(lun,"omega: %d\n",iev);
      if (signal) {
	((TH1D*)gROOT->FindObjectAny("hImpOmK3"))->Fill(pca.Mag());
	((TH1D*)gROOT->FindObjectAny("hChi2VtxOmK"))->Fill(tr->GetChi2Vertex());
	((TH1D*)gROOT->FindObjectAny("hDistOmD"))->Fill(amin);
	((TH1D*)gROOT->FindObjectAny("hImpOmLamb"))->Fill(impL0);
	((TH1D*)gROOT->FindObjectAny("hRadOm"))->Fill(posOm.Mag());
	((TH1D*)gROOT->FindObjectAny("hImpOm"))->Fill(impOm);
      } else {
	((TH1D*)gROOT->FindObjectAny("hImpK3"))->Fill(pca.Mag());
	((TH1D*)gROOT->FindObjectAny("hChi2VtxK"))->Fill(tr->GetChi2Vertex());
	((TH1D*)gROOT->FindObjectAny("hDistOmDBkg"))->Fill(amin);
	((TH1D*)gROOT->FindObjectAny("hRadOmBkg"))->Fill(posOm.Mag());
	((TH1D*)gROOT->FindObjectAny("hImpOmBkg"))->Fill(impOm);
	if (lun) {
	  Int_t idPi = ((MpdKalmanTrack*) itsTracks->UncheckedAt(jPi))->GetTrackID();
	  TH1D *hh = (TH1D*) gROOT->FindObjectAny("hInvMassXiBkg");
	  Double_t mOm = omega.Mag();
	  if (mOm > hh->GetXaxis()->GetXmin() && mOm < hh->GetXaxis()->GetXmax()) fprintf(lun,"%d %d %d %d %d %d %d %d %f\n",iev,mothID,id,((FairMCTrack*) mcTracks->UncheckedAt(id))->GetMotherId(),idPi,((FairMCTrack*) mcTracks->UncheckedAt(idPi))->GetMotherId(),idPlus,((FairMCTrack*) mcTracks->UncheckedAt(idPlus))->GetMotherId(),mOm);
	}
      }
    }
  }
}

//__________________________________________________________________________
void SecondVertex(vector<int>& indx, Int_t iev)
{
  // Reconstruct secondary vertex

  // Get track mothers
  Int_t nTracks = indx.size();
  multimap<int,int> v0;
  TVector3 pos;
  for (Int_t i = 0; i < nTracks; ++i) {
    Int_t j = indx[i];
    MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
    Int_t id = tr->GetTrackID();
    FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
    Int_t pdg = mcTr->GetPdgCode();
    Int_t iQ = Int_t (TDatabasePDG::Instance()->GetParticle(pdg)->Charge());
    mcTr->GetStartVertex(pos);
    Int_t mothID = mcTr->GetMotherId();
    v0.insert(pair<int,int>(mothID,j));
    Int_t mothPDG = ((FairMCTrack*)mcTracks->UncheckedAt(mothID))->GetPdgCode(); // mother PDG
    Int_t gmothID = ((FairMCTrack*)mcTracks->UncheckedAt(mothID))->GetMotherId(); // grandmother
    cout << id << " " << pdg << " " << mothID << " " << mothPDG << " " << gmothID << " " << pos.Pt() << endl;
    if (gmothID > -1) cout << TDatabasePDG::Instance()->GetParticle(((FairMCTrack*) mcTracks->UncheckedAt(gmothID))->GetPdgCode())->GetName() << endl;
  }
  // Check V0's
  Int_t nV0 = v0.size();
  cout << nTracks << " " << nV0 << endl;
  multimap<int,int>::iterator it, it1;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;

  for (it = v0.begin(); it != v0.end(); it++) {
    Int_t mothID = (*it).first;
    if ((*it).second < 0) continue; // already analyzed vertex
    if (v0.count(mothID) != 2) continue; // only one decay track found
    cout << mothID << " " << v0.count(mothID) << endl;
    ret = v0.equal_range(mothID);
    
    Int_t itr = 0, decayInd[2] = {0}; 
    Double_t mass[2] = {0};
    MpdHelix helix[2];
    MpdTrack *trs[2];
    Int_t jPi = -1;
    for (it1 = ret.first; it1 != ret.second; ++it1) {
      Int_t j = (*it1).second;
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      trs[itr] = tr;
      decayInd[itr] = j;
      cout << " " << (*it1).second << " " << tr->GetTrackID() << " " << tr->GetNofIts() << " " << tr->GetChi2() << endl;
      // Create MpdHelix
      Double_t r = tr->GetPosNew();
      Double_t phi = tr->GetParam(0) / r;
      Double_t x = r * TMath::Cos(phi);
      Double_t y = r * TMath::Sin(phi);
      cout << r << " " << x << " " << y << " " << tr->GetParam(2) << " " << phi << " " << tr->GetParam(4) << endl;
      Double_t dip = tr->GetParam(3);
      Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
      cur *= TMath::Abs (tr->GetParam(4));
      TVector3 o(x, y, tr->GetParam(1));
      
      Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
      //helix[itr++] = MpdHelix(cur, dip, phi, o, h);
      helix[itr] = MpdHelix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
      cout << helix[itr].xcenter() << " " << helix[itr].ycenter() << endl;
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(tr->GetTrackID());
      mcTr->GetStartVertex(pos);
      mass[itr++] = TDatabasePDG::Instance()->GetParticle(mcTr->GetPdgCode())->Mass();
      (*it1).second = -1;
      if (((FairMCTrack*)mcTracks->UncheckedAt(tr->GetTrackID()))->GetPdgCode() == pdgCodeNeg) jPi = j;
      if (((FairMCTrack*)mcTracks->UncheckedAt(tr->GetTrackID()))->GetPdgCode() == pdgCodePos) idPlus = tr->GetTrackID();
    }
    //v0.erase(ret.first);
    //v0.erase(ret.second);
    pair<Double_t,Double_t> paths = helix[0].pathLengths(helix[1]);
    cout << " Intersection: " << paths.first << " " << paths.second << endl;
    Double_t rads[2] = {helix[0].at(paths.first).Pt(), helix[1].at(paths.second).Pt()};
    TVector3 pos1 = helix[0].at(paths.first);
    TVector3 pos2 = helix[1].at(paths.second);
    TVector3 d12 = pos2 - pos1;
    TVector3 sum = pos2 + pos1;
    sum *= 0.5;

    if (d12.Mag() > gDistV0) continue; // selection !!!
    if (sum.Mag() < gDecayV0) continue;

    //((TH1D*)gROOT->FindObjectAny("hdcaV0"))->Fill(d12.Mag());
    //((TH1D*)gROOT->FindObjectAny("hdcaV0t"))->Fill(d12.Pt());
    
    cout <<  rads[0] << " " << rads[1] << endl; 
    paths = helix[0].pathLength(pos.Pt(),0.,0.);
    cout << " Moth. rad " << pos.Pt() << " " << paths.first << " " << paths.second << endl;
    paths = helix[1].pathLength(pos.Pt(),0.,0.);
    cout << " Moth. rad " << pos.Pt() << " " << paths.first << " " << paths.second << endl;

    paths = helix[0].pathLength(34.19,0.,0.);
    pos = helix[0].at(paths.first);
    cout << paths.first << " " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
    pos = helix[0].at(paths.second);
    cout << paths.second << " " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;

    paths = helix[1].pathLength(34.19,0.,0.);
    pos = helix[1].at(paths.first);
    cout << paths.first << " " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
    pos = helix[1].at(paths.second);
    cout << paths.second << " " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
    //if (lun) fprintf(lun," %d\n",iev);
    //*
    TVector3 vtx;
    MpdTrack trVtx[2] = {*trs[0], *trs[1]};
    Double_t chi2d = MpdKfV0Fitter::Instance()->FindVertex(trs[0], trs[1], vtx);
    cout << chi2d << " " << trs[0]->GetPosNew() << " " << trs[1]->GetPosNew() << endl;

    if (chi2d > gChi2V0) continue; // selection !!!
    CalcMass(decayInd, rads, mass, 1);

    ((TH1D*)gROOT->FindObjectAny("hChi2V0L0"))->Fill(chi2d);
    ((TH1D*)gROOT->FindObjectAny("hDistV0L0"))->Fill(d12.Mag());
    ((TH1D*)gROOT->FindObjectAny("hRadV0L0"))->Fill(sum.Mag());
    TLorentzVector lambda = CalcMass(trs, mass, vtx, trVtx);
    //if (lun) fprintf(lun,"%f\n",chi2d);

    // Check if Lambda
    if (((FairMCTrack*)mcTracks->UncheckedAt(mothID))->GetPdgCode() == pdgCodeL0) 
      CheckOmega(mothID, lambda, vtx, iev, jPi);

  } // for (it = v0.begin();
}

//__________________________________________________________________________
void FakeVertex(vector<int>& indx, Int_t iev)
{
  // Reconstruct fake vertex

  // Get track mothers
  Int_t nTracks = indx.size();
  multimap<int,int> v0;
  TVector3 pos;
  for (Int_t i = 0; i < nTracks; ++i) {
    Int_t j = indx[i];
    v0.insert(pair<int,int>(0,j));
  }
  // Check V0's
  multimap<int,int>::iterator it, it1;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;

  for (it = v0.begin(); it != v0.end(); it++) {
    Int_t mothID = (*it).first;
    if ((*it).second < 0) continue; // already analyzed vertex
    if (v0.count(mothID) != 2) continue; // only one decay track found
    //cout << mothID << " " << v0.count(mothID) << endl;
    ret = v0.equal_range(mothID);
    
    Int_t itr = 0, decayInd[2] = {0}; 
    Double_t mass[2] = {0};
    MpdHelix helix[2];
    MpdTrack *trs[2];
    Int_t jPi = -1;
    for (it1 = ret.first; it1 != ret.second; ++it1) {
      Int_t j = (*it1).second;
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      trs[itr] = tr;
      decayInd[itr] = j;
      //cout << " " << (*it1).second << " " << tr->GetTrackID() << " " << tr->GetNofIts() << " " << tr->GetChi2() << endl;
      // Create MpdHelix
      Double_t r = tr->GetPosNew();
      Double_t phi = tr->GetParam(0) / r;
      Double_t x = r * TMath::Cos(phi);
      Double_t y = r * TMath::Sin(phi);
      Double_t dip = tr->GetParam(3);
      Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
      cur *= TMath::Abs (tr->GetParam(4));
      TVector3 o(x, y, tr->GetParam(1));
      
      Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
      //helix[itr++] = MpdHelix(cur, dip, phi, o, h);
      helix[itr] = MpdHelix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
      //cout << helix[itr].xcenter() << " " << helix[itr].ycenter() << endl;
      Double_t massPart = TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
      if (tr->GetParam(4) < 0) massPart = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
      mass[itr++] = massPart;
      (*it1).second = -1;
      if (((FairMCTrack*)mcTracks->UncheckedAt(tr->GetTrackID()))->GetPdgCode() == pdgCodeNeg) jPi = j;
      if (((FairMCTrack*)mcTracks->UncheckedAt(tr->GetTrackID()))->GetPdgCode() == pdgCodePos) idPlus = tr->GetTrackID();
    }
    pair<Double_t,Double_t> paths = helix[0].pathLengths(helix[1]);
    TVector3 pos1 = helix[0].at(paths.first);
    TVector3 pos2 = helix[1].at(paths.second);
    TVector3 d12 = pos2 - pos1;
    TVector3 sum = pos2 + pos1;
    sum *= 0.5;

    if (d12.Mag() > gDistV0) continue; // selection !!!
    if (sum.Mag() < gDecayV0) continue;
    
    Double_t rads[2] = {helix[0].at(paths.first).Pt(), helix[1].at(paths.second).Pt()};
    //*
    TVector3 vtx;
    MpdTrack trVtx[2] = {*trs[0], *trs[1]};
    //cout << " here " << d12.Mag() << " " << rads[0] << " " << rads[1] << " " << trs[0]->GetPosNew() << " " << trs[1]->GetPosNew() << " " << trs[0]->GetParam(4) << " " << trs[1]->GetParam(4) << " " << trs[0]->GetPos() << " " << trs[1]->GetPos() << " " << pos1.Z() << " " << trs[0]->GetParam(1) << endl;
    //Double_t chi2d = MpdKfV0Fitter::Instance()->FindVertex(trs[0], trs[1], vtx);
    Double_t chi2d = MpdKfV0Fitter::Instance()->FindVertex(&(trVtx[0]), &(trVtx[1]), vtx);

    if (chi2d > gChi2V0) continue; // selection !!!
    CalcMass(decayInd, rads, mass, 0);

    //cout << " here1 " << chi2d << endl;
    ((TH1D*)gROOT->FindObjectAny("hChi2V0"))->Fill(chi2d);
    ((TH1D*)gROOT->FindObjectAny("hDistV0"))->Fill(d12.Mag());
    ((TH1D*)gROOT->FindObjectAny("hRadV0"))->Fill(sum.Mag());

    TLorentzVector lambda = CalcMass(trs, mass, vtx, trVtx);

    // Find background for omega
    CheckOmega(-1, lambda, vtx, iev, jPi);

  } // for (it = v0.begin();
}

//__________________________________________________________________________
void FindV0(Int_t indx, Int_t iev)
{
  // Find V0

  //return;
  //Int_t nITS = itsTracks->GetEntriesFast();
  MpdTrack *tr0 = (MpdTrack*) itsTracks->UncheckedAt(indx);
  Int_t iQ = Int_t (TMath::Sign(1.,-tr0->GetParam(4)));
  Int_t id0 = tr0->GetTrackID();
  Int_t mothID0 = ((FairMCTrack*)mcTracks->UncheckedAt(id0))->GetMotherId();
  Int_t mothPDG = -1;
  if (mothID0 > -1) {
    FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothID0);
    mothPDG = moth->GetPdgCode();
  }
  vector<int> decay;
  decay.push_back(indx);
  decay.push_back(indx);
    
  for (Int_t j = indx+1; j < nITS; ++j) {
    MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
    //if (tr->GetNofIts() < 1) continue;
    if (tr->GetChi2() < -8) continue; 

    Int_t charge = tr->Charge();
    if (iQ * charge > 0) continue; // the same charge
    Int_t id = tr->GetTrackID();
    Int_t mothID = ((FairMCTrack*)mcTracks->UncheckedAt(id))->GetMotherId();
    if (mothID == mothID0 && mothPDG == pdgCodeL0) continue; // the same mother - skip, because
    //            combinations of true lambdas with other tracks are treated separately 
    Int_t pdg = ((FairMCTrack*)mcTracks->UncheckedAt(id))->GetPdgCode();
    if (pdg != pdgCodeNeg && pdg != pdgCodePos) continue; // only for pi- or proton

    decay[1] = j;
    FakeVertex(decay,iev);
  }
}
    
//__________________________________________________________________________
void AnalOm(Int_t n1 = 0, Int_t n2 = 0)
{
  // Analyze ITS reco data - reconstruct lambda0

  // Add particles to PDG database
  gROOT->ProcessLine(".x $VMCWORKDIR/macro/mpd/AddToPdg.C");
  
  // Load basic libraries
  gROOT->ProcessLine(".x ../loadlibs.C");
  //FairRunAna ana;
  //MpdKalmanFilter::Instance("KF")->Init();

  //gROOT->ProcessLine(".x ./Chain.C(100,\"./urqmd_7gev_its_1.root\")");
  gROOT->ProcessLine(".x ./Chain.C(100,\"./mc_1.root\")");
  TChain *simMC = (TChain*) gROOT->FindObject("cbmsim");
  simMC->SetName("cbmsim1");
  //cout << simMC->GetEntries() << endl;
  //gROOT->ProcessLine(".x ./Chain.C(100,\"./urqmd_7gev_its_1.reco.root\")");
  gROOT->ProcessLine(".x ./Chain.C(100,\"./reco_1.root\")");
  TChain *simITS = (TChain*) gROOT->FindObject("cbmsim");
  //simITS->ls();
  //cout << simITS->GetEntries() << endl;

  //TTree *simITS = (TTree*) fileITS.Get("cbmsim");
  //TTree *simMC = (TTree*) fileMC.Get("cbmsim");

  TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
  TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
  fileMC.Get("FairBaseParSet");
  TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
  simITS->SetBranchAddress("Vertex",&vtxs);
  TBranch *vtxB = simITS->GetBranch("Vertex");
  itsTracks = (TClonesArray*) fileITS.FindObjectAny("ItsTrack");
  simITS->SetBranchAddress("ItsTrack",&itsTracks);
  TBranch *itsRecoB = simITS->GetBranch("ItsTrack");
  if (itsTracks == 0x0) {
    itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
    simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
    itsRecoB = simITS->GetBranch("TpcKalmanTrack");
  }
  //TClonesArray *tpcTracks = (TClonesArray*) fileTPC.FindObjectAny("TpcKalmanTrack");
  //simTPC->SetBranchAddress("TpcKalmanTrack",&tpcTracks);
  //TBranch *tpcRecoB = simTPC->GetBranch("TpcKalmanTrack");
  TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
  simMC->SetBranchAddress("TpcPoint",&tpcPoints);
  TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
  TClonesArray *itsPoints = (TClonesArray*) fileMC.FindObjectAny("StsPoint");
  simMC->SetBranchAddress("StsPoint",&itsPoints);
  TBranch *itsSimB = simMC->GetBranch("StsPoint");
  mcTracks = (TClonesArray*) fileMC.FindObjectAny("MCTrack");
  simMC->SetBranchAddress("MCTrack",&mcTracks);
  TBranch *mcBranch = simMC->GetBranch("MCTrack");

  TFile out("xi_full.histo.root","recreate");

  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();
  MpdTrackFinderIts* recoIts = new MpdTrackFinderIts("ITS track finder");
  recoIts->FillGeoScheme();

  // Book histos

  TH1D *hImp = new TH1D("hImp","DCA (primary in last layer)",100,-1.,1.);
  TH1D *hImpP = new TH1D("hImpP","DCA for pos. (primary in last layer)",100,0.,0.4);
  TH1D *hImpN = new TH1D("hImpN","DCA for neg. (primary in last layer)",100,0.,0.4);
  TH1D *hImpP3 = new TH1D("hImpP3","DCA for pos. (primary in last layer)",100,0.,0.5);
  TH1D *hImpN3 = new TH1D("hImpN3","DCA for neg. (primary in last layer)",100,0.,0.5);
  TH1D *hImpK3 = new TH1D("hImpK3","DCA for kaons (primary in last layer)",100,0.,1.);
  TH1D *hImp2 = new TH1D("hImp2","DCA (secondary in last layer)",100,-1.,1.);
  TH1D *hImpLamb = new TH1D("hImpLamb","DCA of lambda",100,0,1.);
  TH1D *hChi2Vtx = new TH1D("hChi2Vtx","Chi2-distance to prim. vertex (primary in last layer)",100,0.,200.);
  TH1D *hChi2VtxPos = new TH1D("hChi2VtxPos","Chi2-distance to prim. vertex for pos. (primary in last layer)",100,0.,200.);
  TH1D *hChi2VtxNeg = new TH1D("hChi2VtxNeg","Chi2-distance to prim. vertex for neg. (primary in last layer)",100,0.,200.);
  TH1D *hChi2VtxK = new TH1D("hChi2VtxK","Chi2-distance to prim. vertex for bkg. K- (in last layer)",100,0.,200.);
  TH1D *hChi2VtxP = new TH1D("hChi2VtxP","Chi2-distance to prim. vertex (from vertex in last layer)",100,0.,100.);
  TH1D *hImpK0 = new TH1D("hImpK0","DCA (K0 in last layer)",100,-1.,1.);
  TH1D *hImpL0 = new TH1D("hImpL0","DCA (L0 in last layer)",100,-1.,1.);
  TH1D *hImpL0P = new TH1D("hImpL0P","DCA for pos. (L0 in last layer)",100,0.,0.4);
  TH1D *hImpL0N = new TH1D("hImpL0N","DCA for neg. (L0 in last layer)",100,0.,0.4);
  TH1D *hImpL0P3 = new TH1D("hImpL0P3","DCA for pos. (L0 in last layer)",100,0.,0.5);
  TH1D *hImpL0N3 = new TH1D("hImpL0N3","DCA for neg. (L0 in last layer)",100,0.,0.5);
  TH1D *hImpOmK3 = new TH1D("hImpOmK3","DCA for kaons (omega in last layer)",100,0.,1.);
  TH1D *hImpOmLamb = new TH1D("hImpOmLamb","DCA of lambda from omega",100,0,1.);
  TH1D *hImpOm = new TH1D("hImpOm","DCA of omega",100,0,1.);
  TH1D *hImpOmBkg = new TH1D("hImpOmBkg","DCA of omega bkg.",100,0,1.);
  TH1D *hChi2VtxL0 = new TH1D("hChi2VtxL0","Chi2-distance to prim. vertex (L0 in last layer)",100,0.,200.);
  TH1D *hChi2VtxL0P = new TH1D("hChi2VtxL0P","Chi2-distance to prim. vertex for pos. (L0 in last layer)",100,0.,200.);
  TH1D *hChi2VtxL0N = new TH1D("hChi2VtxL0N","Chi2-distance to prim. vertex for neg. (L0 in last layer)",100,0.,200.);
  TH1D *hChi2VtxOmK = new TH1D("hChi2VtxOmK","Chi2-distance to prim. vertex for kaons (omega in last layer)",100,0.,200.);
  TH1D *hChi2V0 = new TH1D("hChi2V0","Chi2 of V0 reco. (in last layer)",100,0.,100.);
  TH1D *hDistV0 = new TH1D("hDistV0","Dist. between daughters in V0 (in last layer)",100,0.,4.);
  TH1D *hRadV0 = new TH1D("hRadV0","Decay length of V0 (in last layer)",100,0.,10.);
  TH1D *hChi2V0L0 = new TH1D("hChi2V0L0","Chi2 of V0 reco. (L0 in last layer)",100,0.,100.);
  TH1D *hDistV0L0 = new TH1D("hDistV0L0","Dist. between daughters in V0 (L0 in last layer)",100,0.,4.);
  TH1D *hDistOmD = new TH1D("hDistOmD","Dist. between daughters in omega vertex (in last layer)",100,0.,1.);
  TH1D *hDistOmDBkg = new TH1D("hDistOmDBkg","Dist. between daughters in omega bkg. vertex (in last layer)",100,0.,1.);
  TH1D *hRadV0L0 = new TH1D("hRadV0L0","Decay length of V0 (L0 in last layer)",100,0.,10.);
  TH1D *hRadOm = new TH1D("hRadOm","Decay length of omega (in last layer)",100,0.,10.);
  TH1D *hRadOmBkg = new TH1D("hRadOmBkg","Decay length of omega bkg. (in last layer)",100,0.,10.);
  //TH1D *hInvMass0 = new TH1D("hInvMass0","Invariant mass",50,1.100,1.130);
  TH1D *hInvMass0 = new TH1D("hInvMass0","Invariant mass",50,1.070,1.170);
  //TH1D *hInvMass0 = new TH1D("hInvMass0","Invariant mass",50,0.400,0.600);
  //TH1D *hInvMass1 = new TH1D("hInvMass1","Invariant mass",50,1.100,1.130);
  TH1D *hInvMass1 = new TH1D("hInvMass1","Invariant mass",50,1.070,1.170);
  //TH1D *hInvMass1 = new TH1D("hInvMass1","Invariant mass",50,0.400,0.600);
  //TH1D *hInvMassBkg = new TH1D("hInvMassBkg","Invariant mass",50,1.100,1.130);
  TH1D *hInvMassBkg = new TH1D("hInvMassBkg","Invariant mass",50,1.070,1.170);
  TH1D *hInvMassXi = new TH1D("hInvMassXi","Invariant mass",50,1.260,1.360);
  //TH1D *hInvMassOm = new TH1D("hInvMassOm","Invariant mass",50,1.650,1.690);
  TH1D *hInvMassOm = new TH1D("hInvMassOm","Invariant mass",50,1.630,1.730);
  TH1D *hInvMassXiBkg = new TH1D("hInvMassXiBkg","Invariant mass",50,1.260,1.360);

  Double_t pmom, eta1, dpp, rorig, ptt;
  Int_t prim, evNo, idtr, np, moth, pdg;
  TTree *tree = new TTree("tracks","Tracks");
  tree->Branch("mom",&pmom,"pmom/D");
  tree->Branch("ptt",&ptt,"ptt/D");
  tree->Branch("eta",&eta1,"eta1/D");
  tree->Branch("dpp",&dpp,"dpp/D");
  tree->Branch("rorig",&rorig,"rorig/D");
  tree->Branch("prim",&prim,"prim/I");
  tree->Branch("id",&idtr,"idtr/I");
  tree->Branch("evNo",&evNo,"evNo/I");
  tree->Branch("np",&np,"np/I");
  tree->Branch("moth",&moth,"moth/I");
  tree->Branch("pdg",&pdg,"pdg/I");

  Int_t events = simITS->GetEntries();
  if (n2 != 0) events = TMath::Min (events, n2);
  cout << " Number of events = " << events << endl;

  for (Int_t i = 0; i < events; ++i) {
    if (i < n1) continue;
    simMC->GetEntry(i);
    simITS->GetEntry(i);
    evNo = i + 1;

    // For ITS points
    TVector3 mom; 
    set<Int_t> idxs;
    if (itsPoints) {
      Int_t nITSp = itsPoints->GetEntriesFast();
      for (Int_t j = 0; j < nITSp; ++j) {
	FairMCPoint *p1 = (FairMCPoint*) itsPoints->UncheckedAt(j);
	Int_t id = p1->GetTrackID();
	if (idxs.find(id) != idxs.end()) continue;
	idxs.insert(id);
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	if (mcTr->GetMotherId() == -1) continue;
	mcTr->GetMomentum(mom);
	//if (TMath::Abs(mom.Eta()) < 1.2) {
	if (TMath::Abs(mom.Eta()) < 1.3) {
	  pdg = mcTr->GetPdgCode();
	  moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	  //if (TMath::Abs(pdg) == 211 && moth == 3122) hsecR[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212 && moth == 3122) hsecR[1]->Fill(mom.Pt()); // proton
	  idxs.insert(id);
	}
      }
    }
    // For UrQMD tracks
    Int_t nMC = mcTracks->GetEntriesFast();
    Int_t skip = 0;
    for (Int_t j = 0; j < nMC; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      //if (mcTr->GetPdgCode() == pdgCodeXi) { skip = 1; break; }
      if (mcTr->GetMotherId() == -1) continue;
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      mcTr->GetMomentum(mom);
      //if (TMath::Abs(mom.Eta()) < 1.2) {
      if (TMath::Abs(mom.Eta()) < 1.3) {
	pdg = mcTr->GetPdgCode();
	moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	if (moth == 3122) {
	  //if (TMath::Abs(pdg) == 211) hsecS[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212) hsecS[1]->Fill(mom.Pt()); // proton
	  if (lun1 && idxs.find(j) == idxs.end()) fprintf(lun1,"%d %d %f %f\n",pdg,j,mom.Pt(),mom.Eta());
	}
      }
    }
    if (skip) continue;

    Int_t nTPC = 0; //tpcTracks->GetEntriesFast();
    nITS = itsTracks->GetEntriesFast();
    Int_t nVert = vtxs->GetEntriesFast();

    cout << " *** Event No: " << i << ", reco tracks in TPC, ITS: " << nTPC << " " << nITS 
	 << ", vertices: " << nVert << endl;
    MpdVertex *vtx = (MpdVertex*) vtxs->First();
    vtx->Position(primVert);
    TArrayI *indxs = vtx->GetIndices();
    set<int> indxVert;
    for (Int_t k = 0; k < indxs->GetSize(); ++k) indxVert.insert((*indxs)[k]);
    cout << indxVert.size() << endl;
    
    // Find TPC track IDs 
    Int_t nPoints = tpcPoints->GetEntriesFast(), idMax = 0;
    for (Int_t j = 0; j < nPoints; ++j) {
      FairMCPoint *p1 = (FairMCPoint*) tpcPoints->UncheckedAt(j);
      idMax = TMath::Max(idMax,p1->GetTrackID());
    }
    cout << " Max ID: " << idMax << endl;
    Int_t *ids = new Int_t [idMax+1];
    lays = new Int_t [idMax+1];
    Int_t *moths = new Int_t [idMax+1];
    Int_t *pdgs = new Int_t [idMax+1];
    Double_t *pt = new Double_t [idMax+1];
    Double_t *th = new Double_t [idMax+1];
    Double_t *rad = new Double_t [idMax+1];
    FairMCPoint **point = new FairMCPoint* [idMax+1];
    MpdTrack **track = new MpdTrack* [idMax+1];

    for (Int_t j = 0; j <= idMax; ++j) { 
      ids[j] = lays[j] = 0; 
      point[j] = 0x0;
      track[j] = 0x0;
    }

    // Get max. reached layer No.
    // !!!!!!!!! Geometry dependent values
    Double_t rMin = 34.19, rMax = 99.81, dR = (rMax-rMin)/50; // 1.3124; // new version (with dead material)
    // !!!!!!!!!
    for (Int_t j = 0; j < nPoints; ++j) {
      FairMCPoint *p1 = (FairMCPoint*) tpcPoints->UncheckedAt(j);
      Int_t id = p1->GetTrackID();
      Double_t radius = p1->GetX() * p1->GetX() + p1->GetY() * p1->GetY();
      radius = TMath::Sqrt(radius);
      Int_t lay = (Int_t) ((radius - rMin) / dR);
      lays[id] = TMath::Max (lay, lays[id]);
    }

    for (Int_t j = 0; j < nITS; ++j) {
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      ids[tr->GetTrackID()]++;
      if (ids[tr->GetTrackID()] > 1) cout << " More than 1 reco track ID: " << tr->GetTrackID() << " " 
					 << ids[tr->GetTrackID()] << endl;
    }
    //cout << " *** " << endl;

    // Get sim. track Pt and Theta - take the last point of track with correct ID
    for (Int_t j = 0; j < nITS; ++j) {
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      Int_t id = tr->GetTrackID();
      if (track[id] == 0x0) track[id] = tr;
      // Get ITS info
      TClonesArray *hits = tr->GetTrHits();
      Int_t nHits = hits->GetEntriesFast();
      FairMCPoint *p1 = 0x0;
      for (Int_t ih = nHits-1; ih >= 0; --ih) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(ih);
	//if (hit->GetDist() < rMin) continue; // ITS hit
	if (hit->GetUniqueID()) continue; // ITS hit
        //TpcPoint *p1 = (TpcPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
        p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
	if (p1->GetTrackID() != id) continue;
	if (ids[id] > 1 && point[id]) {
	  // More than 1 reco track with the same ID
	  //cout << " Time: " << id << " " << p1 << " " << point[id] << " " << p1->GetTime() << " " << point[id]->GetTime() << endl;
	  if (p1 == point[id]) {
	    // The same 1st hit - take "better" track
	    if (tr->GetNofTrHits() - tr->GetNofWrong() < track[id]->GetNofTrHits() - track[id]->GetNofWrong()) {
	      tr->SetChi2(-9.); // exclude this track from further consideration
	      break;
	    } else {
	      // Exclude previous track from further consideration
	      track[id]->SetChi2(-9.);
	      track[id] = tr;
	    }
	  } else if (p1->GetTime() > point[id]->GetTime()) {
	    tr->SetChi2(-9.); // exclude this track from further consideration
	    break;
	  } else {
	    // Exclude previous track from further consideration
	    track[id]->SetChi2(-9.);
	    track[id] = tr;
	  }
	} 
	point[id] = p1;
	p1->Momentum(mom);
	//pt[id] = mom.Pt();
	//th[id] = mom.Theta();
	//Double_t ptRec = TMath::Abs (1./tr->GetParam(4));
	//Double_t thRec = TMath::PiOver2() - tr->GetParam(3);
	//hdpt[0]->Fill(ptRec-pt[id]);
	//hdth[0]->Fill(thRec-th[id]);
	// MC track
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	mcTr->GetMomentum(mom);
	pt[id] = mom.Pt();
	th[id] = mom.Theta();
	break;
      }
    } // for (Int_t j = 0; j < nITS;

    // Flag primary tracks
    multimap<Int_t,Int_t> mapLamb;
    for (Int_t j = 0; j <= idMax; ++j) {
      //if (lays[j] == 0) continue;
      //TpcLheKalmanTrack *tr = (TpcLheKalmanTrack*) tpcTracks->UncheckedAt(j);
      //Int_t id = tr->GetTrackID();
      //if (tr->GetChi2() < -8) { cout << " Loop: " << id << endl; continue; }
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      Int_t mothID = mcTr->GetMotherId();
      if (mothID < 0 && lays[j] != 0) {
	lays[j] = -lays[j];
	//href->Fill(mom.Eta());
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      rad[j] = pos.Pt();
      moths[j] = 0;
      pdgs[j] = mcTr->GetPdgCode();
      if (mothID >= 0) {
	moths[j] = ((FairMCTrack*) mcTracks->UncheckedAt(mothID))->GetPdgCode();
	if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<Int_t,Int_t>(mothID,j));
      }
      //if ((pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg) && moths[j] == pdgCodeL0) 
      //hPtVsEtaS->Fill(TMath::Abs(mom.Eta()),mom.Pt());
    }

    multimap<int,int>::iterator mit, mit1;
    pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;

    for (mit = mapLamb.begin(); mit != mapLamb.end(); ++mit) {
      Int_t mothID = mit->first;
      if (mapLamb.count(mothID) != 2) continue; // only one decay particle
      cout << mapLamb.count(mothID) << endl;
      ret = mapLamb.equal_range(mothID);
      Int_t nppi[2] = {0}, nok = 0;
      for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(mit1->second);
	if (mcTr->GetPdgCode() == pdgCodePos) nppi[0] = 1;
	else if (mcTr->GetPdgCode() == pdgCodeNeg) nppi[1] = 1;
	cout << mcTr->GetPdgCode() << endl;
	mcTr->GetMomentum(mom);
	if (TMath::Abs(mom.Eta()) < 1.2) ++nok;
      }
      if (nppi[0] != 1 || nppi[1] != 1) { cout << " Wrong decay mode !!! " << endl; exit(0); }
    }

    // Track selection 
    for (Int_t j = 0; j < nITS; ++j) {
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      //if (tr->GetNofIts() > 0) continue;
      Int_t id = tr->GetTrackID();
      Double_t thRec = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      if (TMath::Abs(lays[id]) < 41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      Int_t iQ = tr->Charge();
      if (iQ > 0) {
	//if (pdgs[id] != pdgCodePos || tr->GetChi2Vertex() < gChi2pos) tr->SetChi2(-9.);
	if (tr->GetChi2Vertex() < gChi2pos) tr->SetChi2(-9.);
      } else {
	//if (pdgs[id] != pdgCodeNeg || tr->GetChi2Vertex() < gChi2neg) tr->SetChi2(-9.); 
	if (tr->GetChi2Vertex() < gChi2neg) tr->SetChi2(-9.); 
      }
      if (tr->GetNofHits() - tr->GetNofIts() < 10) tr->SetChi2(-9.);
      if (tr->GetChi2() < -8) continue;
      // Create MpdHelix
      Double_t r = tr->GetPosNew();
      Double_t phi = tr->GetParam(0) / r;
      Double_t x = r * TMath::Cos(phi);
      Double_t y = r * TMath::Sin(phi);
      Double_t dip = tr->GetParam(3);
      Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
      cur *= TMath::Abs (tr->GetParam(4));
      TVector3 o(x, y, tr->GetParam(1));
      Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
      MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;
      if (iQ < 0) {
	if (pdgs[id] != pdgCodeKm && pca.Mag() < gImpNeg) tr->SetChi2(-9.);
      }
      else if (iQ > 0 && pca.Mag() < gImpPos) tr->SetChi2(-9.);
    }

    Int_t nPrim = 0, nPrimOut = 0, nSec = 0, nSecIn = 0;
    vector<int> vecK0;
    vector<int> vecL0;
    // Fill resolution histos
    for (Int_t j = 0; j < nITS; ++j) {
      MpdTrack *tr = (MpdTrack*) itsTracks->UncheckedAt(j);
      //if (tr->GetNofIts() < 1) continue;
      Int_t id = tr->GetTrackID();
      if (tr->GetChi2() < -8) { /*cout << " Loop: " << id << " " << tr->GetNofTrHits() << endl;*/ continue; }

      // Create MpdHelix
      Double_t r = tr->GetPosNew();
      Double_t phi = tr->GetParam(0) / r;
      Double_t x = r * TMath::Cos(phi);
      Double_t y = r * TMath::Sin(phi);
      Double_t dip = tr->GetParam(3);
      Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
      cur *= TMath::Abs (tr->GetParam(4));
      TVector3 o(x, y, tr->GetParam(1));
      Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
      MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      vtx->Position(primVert);
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;

      Double_t ptRec = tr->Pt();
      Double_t thRec = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      Double_t pRec = tr->Momentum();
      Double_t pSim = pt[id] / TMath::Sin(th[id]);
      Double_t dpt = (ptRec - pt[id]) / pt[id];
      dpp = (pRec - pSim) / pSim;
      eta1 = -TMath::Log (TMath::Tan(th[id]/2));
      idtr = id;
      np = tr->GetNofTrHits();
      pmom = pSim;
      MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->First();
      FairMCPoint *p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());

      // DCA w.r.t. reconstructed primary vertex
      MpdTrack track1 = *tr;
      track1.SetPos(track1.GetPosNew());
      track1.SetParamNew(*track1.GetParam());
      Double_t vert[3];
      primVert.GetXYZ(vert);
      MpdKalmanFilter::Instance()->FindPca(&track1,vert);

      Double_t dca2 = track1.GetPosNew();
      phi = 0;
      if (dca2 > 1.e-7) phi = track1.GetParam(0) / dca2;
      TVector3 radV(dca2*TMath::Cos(phi), dca2*TMath::Sin(phi), 0);
      TVector3 dirT(TMath::Cos(track1.GetParam(2)), TMath::Sin(track1.GetParam(2)), 0);
      TVector3 vecP = radV.Cross(dirT);

      if (lays[id] < 0) {
	// Primary
	//if (TMath::Abs(lays[id]) < 49) hdpt[1]->Fill(dpt);
	//if (TMath::Abs(lays[id]) < 41) hdpt[1]->Fill(dpt);
	prim = 1;
      }
      //if (TMath::Abs(lays[id]) >= 49) {
      if (TMath::Abs(lays[id]) >= 41 && TMath::Abs(etaRec) < 1.3) {
	// In last layer
	prim = 0;
	rorig = rad[id];
	moth = moths[id];
	pdg = pdgs[id];
	ptt = pt[id];
	TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(pdgs[id]);
	if (part == 0x0) exit(0);
	//Int_t iQ = Int_t (part->Charge());
	Int_t iQ = Int_t (TMath::Sign(1.,-tr->GetParam(4)));
	//cout << j << " ok " << endl;
	if (pdg == pdgCodeNeg || pdg == pdgCodePos)
	  FindV0(j,i); // reconstruct V0
	//cout << " ok " << endl;

	if (lays[id] < 0) {
	  // Primary
	  hImp->Fill(dca2*TMath::Sign(1.,vecP[2]));
	  hChi2Vtx->Fill(tr->GetChi2Vertex());
	  if (iQ > 0) {
	    hImpP->Fill(dca2);
	    hImpP3->Fill(pca.Mag());
	    hChi2VtxPos->Fill(tr->GetChi2Vertex());
	  } else {	   
	    hImpN->Fill(dca2);
	    hImpN3->Fill(pca.Mag());
	    hChi2VtxNeg->Fill(tr->GetChi2Vertex());
	  }
	  prim = 1;
	  ++nPrim;
	  if (indxVert.find(j) == indxVert.end()) ++nPrimOut;
	} else {
	  // Secondary
	  hImp2->Fill(dca2*TMath::Sign(1.,vecP[2]));
	  if (moths[id] == pdgCodeK0) {
	    // K0s
	    hImpK0->Fill(dca2*TMath::Sign(1.,vecP[2]));
	    vecK0.push_back(j);
	    //if (pdgs[id] == -pdgCodeNeg) hdpL0p->Fill(pRec - pSim); 
	    //if (pdgs[id] == pdgCodeNeg) hdpL0n->Fill(pRec - pSim); 
	  } else if (moths[id] == pdgCodeL0 && pdgs[id] < 9999) {
	    // Lambda0
	    hImpL0->Fill(dca2*TMath::Sign(1.,vecP[2]));
	    vecL0.push_back(j);
	    hChi2VtxL0->Fill(tr->GetChi2Vertex());
	    if (pdgs[id] == pdgCodePos) {
	      hImpL0P->Fill(dca2);
	      hImpL0P3->Fill(pca.Mag());
	      hChi2VtxL0P->Fill(tr->GetChi2Vertex());
	    }
	    if (pdgs[id] == pdgCodeNeg) {
	      hImpL0N->Fill(dca2);
	      hImpL0N3->Fill(pca.Mag());
	      hChi2VtxL0N->Fill(tr->GetChi2Vertex());
	    }
	  }
	  ++nSec;
	  if (indxVert.find(j) != indxVert.end()) ++nSecIn;
	}
	tree->Fill();
	if (indxVert.find(j) != indxVert.end()) hChi2VtxP->Fill(tr->GetChi2Vertex());
      } // if (TMath::Abs(lays[id]) >= 49)
    } // for (Int_t j = 0; j < nITS;

    //cout << " Here " << endl;
    SecondVertex(vecL0,i);
    //SecondVertex(vecK0,i);

    delete [] lays;
    delete [] ids;
    delete [] moths;
    delete [] pdgs;
    delete [] pt;
    delete [] th;
    delete [] point;
    delete [] rad;
    delete [] track;
  } // for (Int_t i = 0; i < events;
  
  //TFile out("tpc.histo.root","recreate");
  /*
  TH1D *hEff1 = (TH1D*) hrec->Clone("hEff");
  hEff1->Sumw2();
  hEff1->Divide(href);
  hLayEff[0]->Divide(hLayEff[1]);
  hLayEff[2]->Divide(hLayEff[3]);
  TH2D *hEffV0 = (TH2D*) hPtVsEtaR->Clone("hEffV0");
  hEffV0->Sumw2();
  hEffV0->Divide(hPtVsEtaS);
  */

  out.Write();
  out.Close();
  if (lun) fclose(lun);
  if (lun1) fclose(lun1);
}
