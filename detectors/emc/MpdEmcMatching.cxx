//--------------------------------------------------------------------
//
// Description:
//      MPD TPC-EMC Matching
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 1-June-2016
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-June-2018 - adapted for projective geometry
//
//--------------------------------------------------------------------

#include "MpdEmcMatching.h"
#include "MpdEmcGeoParams.h"
//#include "MpdEmcDigit.h"
#include "MpdEmcMatch.h"
#include "MpdEmcPoint.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcHit.h"
#include "MpdTpcKalmanTrack.h"

#include "FairMCPoint.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TSpline.h>

#include <iterator>

using namespace std;
using namespace TMath;

//FILE *lun = fopen ("file.txt","w"); //AZ

// -----   Default constructor   -------------------------------------------

MpdEmcMatching::MpdEmcMatching() :
FairTask("TPC-EMC matching") {

}

// -----   Destructor   ----------------------------------------------------

MpdEmcMatching::~MpdEmcMatching() {
}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcMatching::Init() {

    cout << "******************* EMC Matching INIT *********************" << endl;

    // Get RootManager
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) {
        cout << "-E- MpdEmcMatching::Init: "
                << "RootManager not instantiated!" << endl;
        return kFATAL;
    }

    fGeoPar = new MpdEmcGeoParams();

    // Containers for rec. points from all+2 sectors
    Int_t nSec = fGeoPar->GetNsec() / 2; // number of long sectors
    multimap<Double_t,Int_t> aaa;
    fRecPoints.assign(nSec+2,aaa);

    // Fill ending rows of sectors
    const Int_t nSecRows = fGeoPar->GetNrows() / (nSec-1); // rows per wide sector
    fSecRows0.insert(nSecRows);
    for (Int_t isec = 1; isec < nSec; ++isec) {
      if (isec == 2 || isec == 6) fSecRows0.insert(*fSecRows0.rbegin()+nSecRows/2);
      else fSecRows0.insert(*fSecRows0.rbegin()+nSecRows);
    }

    // Get input array
    fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
    if (!fPointArray) {
        cout << "-W- MpdEmcMatching::Init: " << "No EmcPoint array!" << endl;
        return kERROR;
    }
    
    fMcTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMcTrackArray) {
        cout << "-W- MpdEmcMatching::Init: " << "No MCTrack array!" << endl;
        return kERROR;
    }

    fRecPointArray = (TClonesArray*) ioman->GetObject("EmcRecPoint");
    if (!fRecPointArray) {
      cout << "-W- MpdEmcMatching::Init: " << "No RecPoint array!" << endl;
      return kERROR;
    }

    fTpcTracks = (TClonesArray*) ioman->GetObject("TpcKalmanTrack");
    if (!fTpcTracks) {
      cout << "-W- MpdEmcMatching::Init: " << "No TPC track array!" << endl;
      return kERROR;
    }

    // Create and register output array
    //*
    fMatchArray = new TClonesArray("MpdEmcMatch", 100);
    ioman->Register("EmcMatch", "EMC", fMatchArray, kTRUE);
    //*/

    cout << "-I- MpdEmcMatching: Intialization successfull" << endl;

    return kSUCCESS;

}

//__________________________________________________________________________

void MpdEmcMatching::Finish() {
    cout << "-I- MpdEmcMatching: Finish" << endl;
}

//__________________________________________________________________________

void MpdEmcMatching::Exec(Option_t* opt) 
{
  // Main processing engine

  cout << "MpdEmcMatching::Exec started" << endl;

  static const Int_t nSec = fGeoPar->GetNsec() / 2; // number of EMC long sectors

  // Reset output Array
  //if (!fDigiArray) Fatal("MpdEmcMatching::Exec", "No array of digits");
  fMatchArray->Delete();

  // Clear rec. point containers
  for (Int_t i = 0; i < nSec+2; ++i) fRecPoints[i].clear();
  
  // Fill rec. point containers
  Int_t nRecP = fRecPointArray->GetEntriesFast();
  cout << " Total number of rec. points: " << nRecP << endl;

  for (Int_t i = 0; i < nRecP; ++i) {
    MpdTpcHit *recp = (MpdTpcHit*) fRecPointArray->UncheckedAt(i);
    Int_t isec = recp->GetLayer();
    //fRecPoints[isec].insert(pair<Double_t,Int_t>(recp->GetZ(),i)); // sort according to theta
    fRecPoints[isec].insert(pair<Double_t,Int_t>(recp->GetLocalZ(),i)); // sort according to theta bins
  }

  Int_t ntpc = fTpcTracks->GetEntriesFast();
  if (ntpc == 0) return;

  for (Int_t i = 0; i < ntpc; ++i) {
    //if (isec < nSec && (*tr->GetParamAtHit())(1,0) > 20 && (*tr->GetParamAtHit())(3,0) > 0) continue; // different sides
    //else if (isec >= nSec && (*tr->GetParamAtHit())(1,0) < -20 && (*tr->GetParamAtHit())(3,0) < 0) continue; // different sides
    DoMatching(i);
  }
}

//__________________________________________________________________________

void MpdEmcMatching::DoMatching(Int_t itrack) 
{
  // Match track with EMC clusters

  static const Int_t nSec = fGeoPar->GetNsec() / 2; // number of EMC long sectors
  static const Double_t rMin = fGeoPar->GetRmin(); // inner radius !!!
  static const Double_t zMax = fGeoPar->GetLength();

  static Int_t first = 1;
  static TSpline3 *rMaxS;

  if (first) {
    // Get rmax-vs-theta dependence
    first = 0;
    const vector<Double_t> &thes = fGeoPar->GetThetaBox();
    const vector<Double_t> &rhos = fGeoPar->GetRhoCenterBox();
    const vector<Double_t> &zs =   fGeoPar->GetZCenterBox();

    Int_t nthe = thes.size();
    Double_t *the = new Double_t [nthe];
    Double_t *rmax = new Double_t [nthe];
    Double_t height = fGeoPar->GetLengthBox(); // tower half-height
    
    Int_t j1 = -1;
    for (Int_t j = nthe-1; j >= 0; --j) {
      Double_t rho = rhos[j];
      Double_t z = zs[j];
      Double_t theta1 = thes[j];
      if (j < nthe-1 && thes[j] <= thes[j+1]+0.1) theta1 = 180 - theta1;
      Double_t costhe = TMath::Cos(theta1*TMath::DegToRad());
      Double_t sinthe = TMath::Sin(theta1*TMath::DegToRad());
      rho += height * sinthe;
      z += height * costhe;
      the[++j1] = TMath::ATan2(rho,z) * TMath::RadToDeg();
      rmax[j1] = rho;
    }
    rMaxS = new TSpline3("rmax",the,rmax,nthe); // rmax vs theta
    delete [] the;
    delete [] rmax;
  }

  MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(itrack);
  if (TMath::Abs((*tr->GetParamAtHit())(1,0)) > zMax) return;

  MpdTpcKalmanTrack tr1(*tr);
  tr1.SetParam(*tr1.GetParamAtHit());
  tr1.SetParamNew(*tr1.GetParamAtHit());
  tr1.SetWeight(*tr1.GetWeightAtHit());
  tr1.SetPos(tr1.GetPosAtHit());
  tr1.SetPosNew(tr1.GetPos());
  tr1.SetLength(tr1.GetLengAtHit());
  //Double_t eta = TMath::Abs (tr->Momentum3().Eta());
  Double_t theta = tr->Theta() * TMath::RadToDeg();
  Double_t pt = tr->Pt();
  
  // For chi2: parameters obtained for box e (pT = 0.2-1.9, eta = -1.1-1.1)
  //Double_t depth = 2.3 + 8.8 * eta;
  //Double_t sigmaL = 0.82 + 1.50 * eta + 0.91 * eta * eta;
  //Double_t sigmaL2 = sigmaL * sigmaL, sigmaT2 = 0.66 * 0.66;
  Double_t sigTel[3] = {1.004, 0.0, 1.386};
  sigTel[1] = 0.051 / TMath::Exp(TMath::Log(pt*pt+1.494e-4*1.494e-4)*0.616) + 0.115;
  Double_t sigLel[4] = {0.789, 0.0, 1.785, 0.528};
  Double_t depth[4] = {0,0,0,-2.0};
  depth[1] = 0.387 / TMath::Exp(TMath::Log(pt*pt+0.382*0.382)*0.737) + 0.210;
  sigLel[1] = 0.064 / TMath::Exp(TMath::Log(pt*pt+0.648*0.648)*2.426) + 0.146;
  depth[0] = -6.080 / TMath::Exp(TMath::Log(pt*pt+0.813*0.813)*2.588) - 2.830;

  // For chi2pi: parameters obtained for box pi (pT = 0.1-2.0, eta = -1.1-1.1)
  //Double_t sig1 = 7.8 - 10.9 * eta + 8.6 * eta * eta;
  //Double_t sig2 = 7.1 - 8.1 * eta + 5.3 * eta * eta;
  //Double_t et = TMath::Min (eta, 1.0);
  //Double_t sigT = 1.4 - 2.2 * et + 5.6 * et * et - 3.2 * et * et * et;
  Double_t sigTpi[3] = {1.222, 0.350, 1.337};
  Double_t sigLpi[3] = {3.21, 0.0, 3.09};
  Double_t depthPi = 0.611 / TMath::Exp(TMath::Log(pt*pt+0.651*0.651)*1.278) + 0.307;
  sigLpi[1] = 0.197 / TMath::Exp(TMath::Log(pt*pt+0.202*0.202)*0.859) + 0.126;
  
  Int_t secs[2], iok = 0;
  Double_t thmin = 999, thmax = -999, thetr[2], phitr[2];
  pair<multimap<Double_t,Int_t>::iterator,multimap<Double_t,Int_t>::iterator> pits[2]; // up to 2 sectors to check

  // Loop over TPC tracks
  MpdKalmanFilter *pKF = MpdKalmanFilter::Instance("KF","KF");

  MpdKalmanHit hEnd; 
  hEnd.SetType(MpdKalmanHit::kFixedR);
  
  // Propagate to EMC min. radius
  hEnd.SetPos(rMin);
  iok = pKF->PropagateToHit(&tr1,&hEnd,kTRUE);
  if (!iok) return;
  Double_t phi = tr1.GetParamNew(0) / tr1.GetPosNew();
  //if (phi < 0) phi += TMath::TwoPi();
  //else if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
  //phitr[0] = phi;
  TVector3 pos (tr1.GetPosNew()*TMath::Cos(phi), tr1.GetPosNew()*TMath::Sin(phi), tr1.GetParamNew(1));
  GetTowerCoords(pos,0,phitr[0],thetr[0]);
  thmin = TMath::Min (thmin,thetr[0]);
  thmax = TMath::Max (thmax,thetr[0]);

  // Propagate to EMC max. radius by small steps
  Double_t r = rMin, rmax = rMaxS->Eval(theta), dr = 4.0, z = 0;
  while (r < rmax) {
    hEnd.SetPos(r);
    iok = pKF->PropagateToHit(&tr1,&hEnd,kTRUE);
    if (!iok) break;
    //if (r < rmin + 1 && TMath::Abs(theT) > zmax) break;
    r += dr;
  }
  if (iok) {
    z = tr1.GetParamNew(1);
    r = tr1.GetPosNew();
    phi = tr1.GetParamNew(0) / r;
  } else {
    z = tr1.GetParam(1);
    r = tr1.GetPos();
    phi = tr1.GetParam(0) / r;
  }
  //phitr[1] = phi;
  pos.SetXYZ(r*TMath::Cos(phi),r*TMath::Sin(phi),z);
  GetTowerCoords(pos,1,phitr[1],thetr[1]);
  thmin = TMath::Min (thmin,thetr[1]);
  thmax = TMath::Max (thmax,thetr[1]);
  thmin -= 2;
  thmax += 2;
  
  for (Int_t io = 0; io < 2; ++io) {
    set<Int_t>::iterator sit = fSecRows0.upper_bound(phitr[io]);
    secs[io] = std::distance(fSecRows0.begin(),sit);
  }
  Int_t dsecs = secs[1] - secs[0];
  if (TMath::Abs(dsecs) > 1 && TMath::Abs(dsecs) < nSec / 2) secs[1] = secs[0] + TMath::Sign(1,dsecs);
  else if (TMath::Abs(dsecs) > nSec / 2) secs[1] = secs[0] - TMath::Sign(1,dsecs);
  if (secs[1] < 0) secs[1] += nSec;
  /*
  for (Int_t io = 0; io < 2; ++io) {
    Double_t r = (io == 0) ? rMin : rMaxS->Eval(theta);
    hEnd.SetPos(r);
    iok = pKF->PropagateToHit(&tr1,&hEnd,kTRUE);
    if (!iok) break;
    ztr[io] = tr1.GetParamNew(1);
    Double_t phi = tr1.GetParamNew(0) / tr1.GetPosNew();
    if (phi < 0) phi += TMath::TwoPi();
    else if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
    phitr[io] = phi;
    zminmax[0] = TMath::Min (zminmax[0],ztr[io]);
    zminmax[1] = TMath::Max (zminmax[1],ztr[io]);
  }
  if (!iok) return;
  zminmax[0] -= 2 * dzTower;
  zminmax[1] += 2 * dzTower;
  
  secs[0] = phitr[0] / dphiSec;
  secs[1] = phitr[1] / dphiSec;
  Int_t dsecs = secs[1] - secs[0];
  if (TMath::Abs(dsecs) > 1 && TMath::Abs(dsecs) < nSec / 2) secs[1] = secs[0] + TMath::Sign(1,dsecs);
  else if (TMath::Abs(dsecs) > nSec / 2 && TMath::Abs(dsecs) < nSec - 1) secs[1] = secs[0] - TMath::Sign(1,dsecs);
  if (secs[1] < 0) secs[1] += nSec;
  */
  
  // Select rec. points (road) to compute distance to the track
  for (Int_t io = 0; io < 2; ++io) {
    if (fRecPoints[secs[io]].count(thmin) > 1 || fRecPoints[secs[io]].count(thmax) > 1)
      { cout << " !!! The same coordinate !!! " << endl; exit(0); }
    pits[io].first = fRecPoints[secs[io]].lower_bound(thmin);
    pits[io].second = fRecPoints[secs[io]].upper_bound(thmax);
  }
   
  phitr[1] = MpdKalmanFilter::Instance()->Proxim(phitr[0],phitr[1],fGeoPar->GetNrows()/TMath::TwoPi());
  TVector2 track (thetr[1]-thetr[0], phitr[1]-phitr[0]);
  Double_t trLeng = track.Mod();
  Int_t imatch = fMatchArray->GetEntriesFast();
  
  for (Int_t io = 0; io < 2; ++io) {
    if (io && secs[io] == secs[0]) continue; // the same sector
    
    for (multimap<Double_t,Int_t>::iterator it = pits[io].first; it != pits[io].second; ++it) {
      MpdTpcHit *recp = (MpdTpcHit*) fRecPointArray->UncheckedAt(it->second);
      Double_t phiRecp = recp->GetLocalX();
      phiRecp = MpdKalmanFilter::Instance()->Proxim(phitr[0],phiRecp,fGeoPar->GetNrows()/TMath::TwoPi());
      TVector2 vect (it->first-thetr[0], phiRecp-phitr[0]);
      Double_t distT = TMath::Sin(vect.DeltaPhi(track)) * vect.Mod();
      Double_t distL = TMath::Cos(vect.DeltaPhi(track)) * vect.Mod();
      MpdEmcMatch *match = new ((*fMatchArray)[imatch++]) MpdEmcMatch(itrack, it->second, distT, distL, 
								      recp->GetQ(), trLeng); 
      // Electrons
      Int_t ireg = 0;
      if (distL >= -1.0 && distL <= trLeng + 1.0) ireg = 1;
      else if (distL > trLeng) { distL -= (trLeng + 1.0); ireg = 2; }
      Double_t chi2 = distT * distT / sigTel[ireg] / sigTel[ireg];
      if (distL > -3.0 && distL < -1.0) ireg = 3;
      chi2 += (distL - depth[ireg]) * (distL - depth[ireg]) / sigLel[ireg] / sigLel[ireg];
      match->SetChi2(chi2);
      // Pions
      ireg = 0;
      if (distL >= -1.0 && distL <= trLeng + 1.0) ireg = 1;
      else if (distL > trLeng) { distL -= (trLeng + 1.0); ireg = 2; }
      chi2 = distT * distT / sigTpi[ireg] / sigTpi[ireg];
      if (ireg == 1) chi2 += (distL - depthPi) * (distL - depthPi) / sigLpi[ireg] / sigLpi[ireg];
      else chi2 += (distL * distL / sigLpi[ireg] / sigLpi[ireg]);
      match->SetChi2pi(chi2);
    }
  }
}

//__________________________________________________________________________

void MpdEmcMatching::GetTowerCoords(TVector3 &pos, Int_t io, Double_t &phiT, Double_t &theT)
{
  // Transform point coordinates to tower indices

  static Int_t first = 1, offset = 0, nphi = 0;
  static TSpline3 *phiS, *theSmin, *theSmax;

  if (first) {
    // Get phi and theta angles of the tower centers at their inner face                                
    first = 0;
    const vector<Double_t> &phis = fGeoPar->GetPhiRow();
    const vector<Double_t> &thes = fGeoPar->GetThetaBox();
    const vector<Double_t> &rhos = fGeoPar->GetRhoCenterBox();
    const vector<Double_t> &zs =   fGeoPar->GetZCenterBox();

    nphi = phis.size();
    // Offset due to the fact that the 1'st sector starts at phi = -Phi_sec/2; 
    offset = nphi / (fGeoPar->GetNsec()/2 - 1) / 2;
    Double_t *phia = new Double_t [nphi];
    Double_t *ind = new Double_t [nphi];

    for (Int_t j = 0; j < nphi; ++j) {
      phia[j] = phis[j];
      ind[j] = j;
    }
    phiS = new TSpline3("phis",phia,ind,nphi); // ind vs phi                     
    delete [] phia;
    delete [] ind;

    Int_t nthe = thes.size();
    Double_t *theMin = new Double_t [nthe];
    Double_t *theMax = new Double_t [nthe];
    ind = new Double_t [nthe];
    Double_t height = fGeoPar->GetLengthBox(); // tower half-height
    
    Int_t j1 = -1;
    for (Int_t j = nthe-1; j >= 0; --j) {
      Double_t rho = rhos[j];
      Double_t z = zs[j];
      Double_t theta1 = thes[j];
      if (j < nthe-1 && thes[j] <= thes[j+1]+0.1) theta1 = 180 - theta1;
      Double_t costhe = TMath::Cos(theta1*TMath::DegToRad());
      Double_t sinthe = TMath::Sin(theta1*TMath::DegToRad());
      rho -= height * sinthe;
      z -= height * costhe;
      theMin[++j1] = TMath::ATan2(rho,z) * TMath::RadToDeg(); // at inner tower face
      ind[j1] = j;
      rho += 2 * height * sinthe;
      z += 2 * height * costhe;
      theMax[j1] = TMath::ATan2(rho,z) * TMath::RadToDeg(); // at outer face
    }
    theSmin = new TSpline3("theMin",theMin,ind,nthe); // ind vs theta
    theSmax = new TSpline3("theMax",theMax,ind,nthe); // ind vs theta
    delete [] theMin;
    delete [] theMax;
    delete [] ind;
  }

  // Get Theta and Phi tower indices at inner or outer face
  Double_t phi = pos.Phi() * TMath::RadToDeg();
  Double_t theta = pos.Theta() * TMath::RadToDeg();
  if (phi < 0) phi += 360;
  phiT = phiS->Eval(phi) + offset;
  if (phiT > nphi) phiT -= nphi;
  if (io == 0) theT = theSmin->Eval(theta);
  else theT = theSmax->Eval(theta);
}
//__________________________________________________________________________

ClassImp(MpdEmcMatching)