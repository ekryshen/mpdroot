//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class MpdEmcClusterFinderAZ
//      see MpdEmcClusterFinderAZ.h for details
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 18-May-2016
//      Alexander Zinchenko LHEP, JINR, Dubna - 24-June-2018 - adapted for projective geometry
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "MpdEmcClusterFinderAZ.h"

// Collaborating Class Headers --------
#include "MpdTpc2dCluster.h"
#include "MpdEmcDigit.h"
#include "MpdEmcGeoPar.h"
#include "MpdEmcGeoParams.h"
#include "MpdEmcPoint.h"
#include "MpdTpcHit.h"

#include "MpdMCTrack.h"
#include "FairRootManager.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TTree.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <iterator>
#include <math.h>
#include <set>
#include <vector>

using namespace std;

//__________________________________________________________________________

MpdEmcClusterFinderAZ::MpdEmcClusterFinderAZ()
  : FairTask("EMC Cluster finder AZ"), fPersistence(kFALSE), fThresh(0.005), fTimeWindow(0)
{
}

//__________________________________________________________________________

MpdEmcClusterFinderAZ::~MpdEmcClusterFinderAZ()
{
}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::FinishTask()
{
}

//__________________________________________________________________________

InitStatus MpdEmcClusterFinderAZ::Init()
{

  // Create containers for digits
  //fEmcGeo = new MpdEmcGeoPar();
  fEmcGeo = new MpdEmcGeoParams();
  Int_t nSec = fEmcGeo->GetNsec() / 2; // number of long sectors
  set<Int_t> aaa;
  fDigiSet.assign(nSec,aaa);
  //for (Int_t i = 0; i < fgkNsec2; ++i) fDigiSet[i] = new set<Int_t> [nRows];

  // Containers for all+2 sectors
  //const Int_t nPhi = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi() * (nSec + 2);
  //const Int_t nZ = fEmcGeo->GetNrows() * fEmcGeo->GetNModInSuperModByZ() * 2;
  const Int_t nSecRows = fEmcGeo->GetNrows() / (nSec-1); // rows per wide sector
  const Int_t nPhi = fEmcGeo->GetNrows() + nSecRows * 2; // extra 2 wide sectors
  const Int_t nZ = fEmcGeo->GetNmod();

  // Fill starting rows of sectors
  fSecRows0.push_back(0);
  for (Int_t isec = 1; isec <= nSec; ++isec) {
    if (isec == 3 || isec == 7) fSecRows0.push_back(fSecRows0[isec-1]+nSecRows/2);
    else fSecRows0.push_back(fSecRows0[isec-1]+nSecRows);
  }

  // Fill ending rows of sectors
  fSecRows1.insert(nSecRows);
  for (Int_t isec = 1; isec < nSec; ++isec) {
    if (isec == 2 || isec == 6) fSecRows1.insert(*fSecRows1.rbegin()+nSecRows/2);
    else fSecRows1.insert(*fSecRows1.rbegin()+nSecRows);
  }

  vector<Double_t> dVec(nZ);
  vector<Int_t> iVec(nZ);
  fCharges.assign(nPhi,dVec);
  fFlags.assign(nPhi,iVec);
  fDigis.assign(nPhi,iVec);

  //Get ROOT Manager
  FairRootManager* ioman = FairRootManager::Instance();

  if (ioman == 0) {
    Error("MpdEmcClusterFinderAZ::Init","RootManager not instantiated!");
    return kERROR;
  }
  
  // Get input collection
  fDigiArray = (TClonesArray*) ioman->GetObject("EmcDigit");
  
  if (fDigiArray == 0) {
    Error("EmcClusterFinderAZ::Init","Array of digits not found!");
    return kERROR;
  }
  
  fMcTrArray = (TClonesArray*) ioman->GetObject("MCTrack");
  if (!fMcTrArray) {
    Error("EmcClusterFinderAZ::Init","Array of MCTracks not found!");
    return kERROR;
  }

  // Create and register output array
  SetPersistence();
  fClusArray = new TClonesArray("MpdTpc2dCluster"); 
  ioman->Register("EmcCluster", "Emc", fClusArray, fPersistence);
  fHitArray = new TClonesArray("MpdTpcHit"); 
  ioman->Register("EmcRecPoint", "Emc", fHitArray, fPersistence);

  cout << fEmcGeo->GetPhiRow().size() << " " << fEmcGeo->GetThetaBox().size() << endl;
  Double_t theBeg = fEmcGeo->GetThetaBox()[0];
  theBeg = TMath::Max (theBeg,180-theBeg);
  cout << fEmcGeo->GetPhiRow()[0] << " " << *fEmcGeo->GetPhiRow().rbegin() << " " 
       << theBeg << " " << *fEmcGeo->GetThetaBox().rbegin() << endl;
  return kSUCCESS;
}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::Exec(Option_t* opt)
{

  fClusArray->Delete();
  fHitArray->Delete();
  static const Int_t nSec = fEmcGeo->GetNsec() / 2; // number of EMC sectors
  //static const Int_t nTowSec = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi();
  // Clear digi containers
  for (Int_t i = 0; i < nSec; ++i) fDigiSet[i].clear();
  
  // Fill digi containers
  Int_t nDigis = fDigiArray->GetEntriesFast();
  cout << " Total number of digits: " << nDigis << endl;
  for (Int_t i = 0; i < nDigis; ++i) {
    MpdEmcDigit *digi = (MpdEmcDigit*) fDigiArray->UncheckedAt(i);
    //Int_t isec = digi->GetChanPhiId() / nTowSec;
    Int_t isec = digi->Sector();
    fDigiSet[isec].insert(i);
  }

  FillEmcInfo(); // fill info (matrices)

  Int_t nSum = 0;
  // Loop over sectors
  for (Int_t isec = 0; isec < nSec; ++isec) {
    if (fDigiSet[isec].size() == 0) continue;
    //cout << " Sector, digits: " << isec << " " << fDigiSet[isec].size() << endl; 
    ProcessSector(isec);
    nSum += fDigiSet[isec].size();
  }
  cout << " Control sum: " << nSum << endl;

  // Find hits
  FindHits();
  //FindHitsLocMax();
}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::FillEmcInfo()
{
  // Fill EMC info (matrices)

  static const Int_t nPhi = fFlags.size(), nZ = fFlags[0].size(), nZov2 = nZ / 2;
  static const Int_t nSec = fEmcGeo->GetNsec() / 2; // number of EMC sectors
  //static const Int_t nTowSec = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi();
  static const Int_t nSecRows = nPhi / (nSec + 1);
  static const Double_t rmin = fEmcGeo->GetRmin();

  for (Int_t i = 0; i < nPhi; ++i) fFlags[i].assign(nZ,0);

  set<Int_t>::iterator it;

  for (Int_t isec = 0; isec < nSec; ++isec) {
 
    for (it = fDigiSet[isec].begin(); it != fDigiSet[isec].end(); ++it) {
      MpdEmcDigit *digi = (MpdEmcDigit*) fDigiArray->UncheckedAt(*it);

      // Apply threshold
      if (digi->GetE() < fThresh) continue;

      // Apply time window
      //Double_t dist = TMath::Sqrt (digi->GetZcenter() * digi->GetZcenter() + rmin * rmin);
      Double_t dist = rmin / TMath::Sin(digi->GetZcenter());
      Double_t dt = digi->GetTimeStamp() - dist / 30.; // c = 30 cm/ns       
      if (fTimeWindow && (dt < -0.5 || dt > 2.0)) continue;

      Int_t iphi = digi->GetChanPhiId(), iz = digi->GetChanZId();
      Int_t ix = iphi + nSecRows; // offset by one sector
      //if (digi->Side() == 1) iz = -iz - 1;
      if (digi->Side() == 1) iz = -iz - 0; // shift by 1 already accounted for in ChanZId
      iz += nZov2;
      fCharges[ix][iz] = digi->GetE();
      fFlags[ix][iz] = 1;
      fDigis[ix][iz] = *it; // index of digit
      if (isec == 0) {
	// Fill extra sector (for edge effect)
	Int_t offset = fEmcGeo->GetNrows();
	Int_t ix1 = ix + offset;
	fCharges[ix1][iz] = digi->GetE();
	fFlags[ix1][iz] = 1;
	fDigis[ix1][iz] = *it;
      } else if (isec == nSec - 1) {
	Int_t ix1 = (ix + 2 * nSecRows) % nPhi;
	fCharges[ix1][iz] = digi->GetE();
	fFlags[ix1][iz] = 1;
	fDigis[ix1][iz] = *it;
      }
    }
  }
}
  
//__________________________________________________________________________

void MpdEmcClusterFinderAZ::ProcessSector(Int_t isec)
{
  // Process one sector

  const Int_t nSec = fEmcGeo->GetNsec() / 2; // number of EMC sectors
  //const Int_t nTowSec = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi();
  const Int_t nPhi = fFlags.size(), nZ = fFlags[0].size();
  static const Int_t nSecRows = nPhi / (nSec + 1);

  // Find (pre)clusters in the central sector
  Int_t nclus0 = fClusArray->GetEntriesFast(), nclus = nclus0;
  //Int_t phiBeg = (isec + 1) * nTowSec, phiEnd = phiBeg + nTowSec;
  Int_t phiBeg = fSecRows0[isec] + nSecRows, phiEnd = fSecRows0[isec+1] + nSecRows;

  for (Int_t iphi = phiBeg; iphi < phiEnd; ++iphi) {

    for (Int_t iz = 0; iz < nZ; ++iz) {
      if (fFlags[iphi][iz] <= 0) continue;
      // New cluster
      MpdTpc2dCluster* clus = new ((*fClusArray)[nclus++]) MpdTpc2dCluster(0, isec); 
      clus->Insert(fDigis[iphi][iz], 0, iphi-nSecRows, iz, fCharges[iphi][iz]);
      clus->SetID(fDigis[iphi][iz]); // trackID
      clus->SetErrY(fCharges[iphi][iz]); // ADC counts
      fFlags[iphi][iz] = -1;
      //if (isec == 0) fFlags[iphi+nPhi-2*nTowSec][iz] = -1; 
      //else if (isec == nSec-1) fFlags[iphi-nPhi+2*nTowSec][iz] = -1; 
      if (isec == 0) fFlags[iphi+nPhi-2*nSecRows][iz] = -1; 
      else if (isec == nSec-1) fFlags[(iphi+2*nSecRows)%nPhi][iz] = -1; 
      for (Int_t ip = -1; ip < 2; ++ip) {
	for (Int_t it = -1; it < 2; ++it) {
	  //if (it == ip) continue;
	  if (it == 0 && ip == 0) continue;
	  Int_t ip1 = iphi + ip, iz1 = iz + it;
	  if (ip1 < 0 || ip1 >= nPhi) continue;
	  if (iz1 < 0 || iz1 >= nZ) continue;
	  if (fFlags[ip1][iz1] <= 0) continue;
	  NextPixel(clus, ip1, iz1);
	}
      }
    }
  }
  //cout << " Found preclusters: " << nclus - nclus0 << endl;
}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::NextPixel(MpdTpc2dCluster* clus, Int_t iphi, Int_t iz)
{
  // Add next pixel to the cluster

  const Int_t nPhi = fFlags.size(), nZ = fFlags[0].size();
  const Int_t nSec = fEmcGeo->GetNsec() / 2; // number of EMC sectors
  //const Int_t nTowSec = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi();
  const Int_t nSecRows = nPhi / (nSec + 1);

  clus->Insert(fDigis[iphi][iz], clus->Row(), iphi-nSecRows, iz, fCharges[iphi][iz]);
  clus->SetID (TMath::Min(clus->ID(),fDigis[iphi][iz])); // min trackID
  clus->SetErrY (TMath::Max(Double_t(clus->GetErrY()),fCharges[iphi][iz])); // max ADC counts
  fFlags[iphi][iz] = -1;
  // Edges
  if (iphi < nSecRows) fFlags[iphi+nPhi-2*nSecRows][iz] = -1; // last sector
  else if (iphi >= nPhi-nSecRows) fFlags[iphi-nPhi+2*nSecRows][iz] = -1; // first sector
  else if (iphi >= nPhi-2*nSecRows) fFlags[iphi-nPhi+2*nSecRows][iz] = -1; // last sector
  else if (iphi < 2*nSecRows) fFlags[iphi+nPhi-2*nSecRows][iz] = -1; // first sector

  for (Int_t ip = -1; ip < 2; ++ip) {
    for (Int_t it = -1; it < 2; ++it) {
      //if (TMath::Abs(it) == TMath::Abs(ip)) continue;
      Int_t ip1 = iphi + ip, it1 = iz + it;
      if (ip1 < 0 || ip1 >= nPhi) continue;
      if (it1 < 0 || it1 >= nZ) continue;
      if (fFlags[ip1][it1] <= 0) continue;
      NextPixel(clus, ip1, it1);
    }
  }
}

//__________________________________________________________________________

/*
void MpdEmcClusterFinderAZ::FindHits()
{
  // Reconstruct hits (one hit per precluster)

  const Int_t nTowSec = fEmcGeo->GetNsupMod() * fEmcGeo->GetNModInSuperModByPhi();
  TVector3 p3loc, p3glob, p3err(0.05,0.0,0.1);
  Int_t nclus = fClusArray->GetEntriesFast();

  for (Int_t iclus = 0; iclus < nclus; ++iclus) {
    MpdTpc2dCluster* clus = (MpdTpc2dCluster*) fClusArray->UncheckedAt(iclus);
    Int_t nDigis = clus->NDigits();
    Double_t phiMean = 0, zMean = 0;

    for (Int_t idig = 0; idig < nDigis; ++idig) {
      phiMean += clus->Col(idig) * clus->Adc(idig);
      zMean += clus->Bkt(idig) * clus->Adc(idig);
    }
    phiMean /= clus->GetADC();
    phiMean -= nTowSec;
    zMean /= clus->GetADC();
    clus->SetX(phiMean);
    clus->SetZ(zMean);
    phiMean *= fEmcGeo->GetAngleOfModule() * TMath::DegToRad();
    //phiMean += clus->GetSect() * fEmcGeo->GetAngleOfSector() * TMath::DegToRad();
    clus->SetGlobX(phiMean);
    zMean -= (fFlags[0].size() / 2 - 0.5);
    zMean *= fEmcGeo->GetLengthOfModuleByZ();
    clus->SetGlobZ(zMean);

  }
}
*/

//__________________________________________________________________________
 /*
void MpdTpcClusterFinderAZ::FindHitsLocMax()
{
  // Reconstruct hits (find local maxima)

  TVector3 p3loc, p3glob, p3err(0.05,0.0,0.1);
  Int_t nclus = fClusArray->GetEntriesFast(), ihit = 0;

  for (Int_t iclus = 0; iclus < nclus; ++iclus) {
    MpdTpc2dCluster* clus = (MpdTpc2dCluster*) fClusArray->UncheckedAt(iclus);
    Int_t nDigis = clus->NDigits();

    memset(fFlags, 0, sizeof(fFlags[0][0]) * fgkNpads * fgkNtimes);
    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t ipad = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      fCharges[ipad][itime] = clus->Adc(idig);
      fFlags[ipad][itime] = 1;
      fDigis[ipad][itime] = idig;
    }

    // Exclude pads which are not local maxima
    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t ipad = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      for (Int_t ip = -1; ip < 2; ++ip) {
	for (Int_t it = -1; it < 2; ++it) {
	  //if (TMath::Abs(it) == TMath::Abs(ip)) continue; // exclude diagonals
	  if (it == 0 && ip == 0) continue;
	  Int_t ip1 = ipad + ip, it1 = itime + it;
	  if (ip1 < 0 || ip1 >= fgkNpads) continue;
	  if (it1 < 0 || it1 >= fgkNtimes) continue;
	  if (fFlags[ip1][it1] == 0) continue;
	  if (clus->Adc(idig) < fCharges[ip1][it1]) { fFlags[ipad][itime] = -1; break; }
	}
      }
    }
    multimap<Double_t,Int_t> localMax;
    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t ipad = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      if (fFlags[ipad][itime] <= 0) continue;
      localMax.insert(pair<Double_t,Int_t>(clus->Adc(idig),idig));
      cout << clus->Col(idig) << " " << clus->Bkt(idig) << " " << clus->Adc(idig) << endl;
    }
    cout << " Local max: " << localMax.size() << endl;

    multimap<Double_t,Int_t>::reverse_iterator rit = localMax.rbegin(), rit1 = rit;
    vector<Int_t> vecDig;
    for ( ; rit != localMax.rend(); ++rit) {
      Int_t idig = rit->second;
      Int_t ipad = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      if (fFlags[ipad][itime] <= 0) continue; // merged local max
      vecDig.clear();

      Double_t padMean = 0, timeMean = 0, adcTot = 0;
      for (Int_t ip = -1; ip < 2; ++ip) {
	for (Int_t it = -1; it < 2; ++it) {
	  //if (TMath::Abs(it) == TMath::Abs(ip)) continue;
	  Int_t ip1 = ipad + ip, it1 = itime + it;
	  if (ip1 < 0 || ip1 >= fgkNpads) continue;
	  if (it1 < 0 || it1 >= fgkNtimes) continue;
	  if (fFlags[ip1][it1] == 0) continue;
	  if (fFlags[ip1][it1] > 0) fFlags[ip1][it1] = -1; // merge neighbour local max
	  padMean += ip1 * fCharges[ip1][it1];
	  timeMean += it1 * fCharges[ip1][it1];
	  adcTot += fCharges[ip1][it1];
	  vecDig.push_back(fDigis[ip1][it1]);
	}
      }
      padMean /= adcTot;
      timeMean /= adcTot;

      Double_t xloc = fSecGeo->Pad2Xloc(padMean,clus->Row());
      Int_t padID = fSecGeo->PadID(clus->GetSect() % fSecGeo->NofSectors(), clus->Row());
      Double_t yloc = fSecGeo->LocalPadPosition(padID).Y();
      Double_t zloc = fSecGeo->TimeBin2Z(timeMean);
      p3loc.SetXYZ(xloc, yloc, zloc);

      // Apply corrections
      TVector3 p3errCor(p3err);
      CorrectReco(p3loc, p3errCor, clus->GetNumPads(), adcTot);

      fSecGeo->Local2Global(clus->GetSect(), p3loc, p3glob);
      if (clus->GetSect() >= fSecGeo->NofSectors()) p3glob[2] = -p3glob[2];

      MpdTpcHit* hit = new ((*fHitArray)[ihit++]) MpdTpcHit(padID, p3glob, p3errCor, iclus); 
      hit->SetLayer(clus->Row());
      hit->SetLocalPosition(p3loc); // point position
      hit->SetEnergyLoss(adcTot);
      hit->SetStep(0.0);
      hit->SetModular(1); // modular geometry flag
      hit->SetPad(Int_t(padMean));
      hit->SetBin(Int_t(timeMean));

      // !!! Warning: FairLinks are not persistent !!! 
      //hit->AddLink(FairLink(MpdTpcHit::PointIndex, pointIndx));
      Int_t ndig = vecDig.size();
      for (Int_t idig = 0; idig < ndig; ++idig)
	hit->AddLink(FairLink(MpdTpcHit::MCTrackIndex, clus->Sec(vecDig[idig]))); // trackID stored in Sec
      //cout << hit->GetNLinks();

    } // for ( ; rit != localMax.rend();
  } // for (Int_t iclus = 0; iclus < nclus;
}
*/
//__________________________________________________________________________

void MpdEmcClusterFinderAZ::FindHits()
{
  // Reconstruct hits (find local maxima)
  
  static const Int_t nPhi = fFlags.size(), nZ = fFlags[0].size();
  //static const Double_t dphiSec = fEmcGeo->GetAngleOfSector() * TMath::DegToRad();

  TVector3 p3loc, p3glob, p3err(0.05,0.0,0.1);
  Int_t nclus = fClusArray->GetEntriesFast(), ihit = 0;

  for (Int_t iclus = 0; iclus < nclus; ++iclus) {
    MpdTpc2dCluster* clus = (MpdTpc2dCluster*) fClusArray->UncheckedAt(iclus);
    Int_t nDigis = clus->NDigits();
    for (Int_t i = 0; i < nPhi; ++i) fFlags[i].assign(nZ,0);
    Int_t isec = clus->GetSect();
    Int_t ishift = 0;
    //if (isec == 0) ishift = nPhi / 30; // just to shift index up from 0
    if (isec == 0) ishift = 50; // just to shift index up from 0

    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t iphi = clus->Col(idig) + ishift;
      Int_t iz = clus->Bkt(idig);
      fCharges[iphi][iz] = clus->Adc(idig);
      fFlags[iphi][iz] = 1;
      //fDigis[iphi][iz] = idig;
      fDigis[iphi][iz] = clus->Sec(idig);
    }

    // Exclude pads which are not local maxima
    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t iphi = clus->Col(idig) + ishift;
      Int_t iz = clus->Bkt(idig);
      for (Int_t ip = -1; ip < 2; ++ip) {
	for (Int_t it = -1; it < 2; ++it) {
	  //if (TMath::Abs(it) == TMath::Abs(ip)) continue; // exclude diagonals
	  if (it == 0 && ip == 0) continue;
	  Int_t ip1 = iphi + ip, it1 = iz + it;
	  if (ip1 < 0 || ip1 >= nPhi) continue;
	  if (it1 < 0 || it1 >= nZ) continue;
	  if (fFlags[ip1][it1] == 0) continue;
	  if (clus->Adc(idig) < fCharges[ip1][it1]) { fFlags[iphi][iz] = -1; break; }
	}
      }
    }
    multimap<Double_t,Int_t> localMax;
    for (Int_t idig = 0; idig < nDigis; ++idig) {
      Int_t iphi = clus->Col(idig) + ishift;
      Int_t iz = clus->Bkt(idig);
      if (fFlags[iphi][iz] <= 0) continue;
      localMax.insert(pair<Double_t,Int_t>(clus->Adc(idig),idig));
      //cout << clus->Col(idig) << " " << clus->Bkt(idig) << " " << clus->Adc(idig) << endl;
    }
    //cout << " Local max: " << clus->GetSect() << " " << clus->Row() << " " << localMax.size() << endl;

    // Remove local maxima not separated by valleys - Peak and Valley
    //Int_t nLocMax0 = localMax.size();
    if (localMax.size() > 1) PeakAndValley(clus, localMax, ishift);

    multimap<Double_t,Int_t>::reverse_iterator rit = localMax.rbegin();
    //vector<Int_t> vecDig;
    map<Int_t,Double_t> mapIdQ;
    for ( ; rit != localMax.rend(); ++rit) {
      Int_t idig = rit->second;
      Int_t iphi = clus->Col(idig) + ishift;
      Int_t iz = clus->Bkt(idig);
      if (fFlags[iphi][iz] <= 0) continue; // merged local max
      //vecDig.clear();
      mapIdQ.clear();
      set<pair<Int_t,Int_t> > pixels;

      Double_t phiMean = 0, zMean = 0, adcTot = 0, sum2t = 0, sum2p = 0;

      // Process simple cluster (only 1 local max)
      if (localMax.size() == 1) {
      //AZ if (nLocMax0 == 1) {
	for (Int_t idig1 = 0; idig1 < nDigis; ++idig1) {
	  Int_t ip = clus->Col(idig1) + ishift;
	  Int_t it = clus->Bkt(idig1);
	  phiMean += ip * fCharges[ip][it];
	  zMean += it * fCharges[ip][it];
	  adcTot += fCharges[ip][it];
	  //vecDig.push_back(fDigis[ip][it]);
	  //Int_t id = clus->Sec(fDigis[ip][it]);
	  //if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = fCharges[ip][it];
	  //else mapIdQ[id] = mapIdQ[id] + fCharges[ip][it];
	  //else mapIdQ[id] = TMath::Max (mapIdQ[id], fCharges[ip][it]);
	  //if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = 1;
	  //else mapIdQ[id] = mapIdQ[id] + 1;
	  MpdEmcDigit *dig = (MpdEmcDigit*) fDigiArray->UncheckedAt(fDigis[ip][it]);
	  map<Int_t,Float_t> contrib = dig->GetContrib();
	  RedoId(contrib);
	  for (map<Int_t,Float_t>::iterator mit = contrib.begin(); mit != contrib.end(); ++mit) {
	    if (mapIdQ.find(mit->first) == mapIdQ.end()) mapIdQ[mit->first] = mit->second;
	    else mapIdQ[mit->first] = mapIdQ[mit->first] + mit->second;
	  }

	  sum2t += it * it * fCharges[ip][it];
	  sum2p += ip * ip * fCharges[ip][it];
	}
      } else {
	// Process complex cluster - start from maximum and go ouward (up tp 5 steps), 
	// adding pixels with repectively lower charges
	for (Int_t idirp = -1; idirp < 2; idirp += 2) {
	  for (Int_t ip = 0; ip < 5; ++ip) {
	    if (idirp > 0 && ip == 0) continue;
	    Int_t ipsign = ip * idirp;
	    Int_t ip1 = iphi + ipsign;
	    if (ip1 < 0 || ip1 >= nPhi) break;

	    for (Int_t idirt = -1; idirt < 2; idirt += 2) {
	      for (Int_t it = 0; it < 5; ++it) {
		if (idirt > 0 && it == 0) continue;
		Int_t itsign = it * idirt;

		Int_t it1 = iz + itsign;
		if (it1 < 0 || it1 >= nZ) break;
		if (fFlags[ip1][it1] == 0) continue;

		Int_t add = 1;
		if (ip || it) {
		  if (fFlags[ip1][it1] > 0) continue; // case when 2 local max next to each other on 1 diagonal 
		  // Check gradient
		  add = 0;
		  Int_t ipprev = ipsign, itprev = itsign;
		  if (it) {
		    itprev -= idirt;
		    Int_t it10 = iz + itprev;
		    if (it10 >= 0 && it10 < nZ && fFlags[ip1][it10] != 0) {
		      if (pixels.find(pair<Int_t,Int_t>(ip1,it10)) != pixels.end() 
			  && fCharges[ip1][it1] <= fCharges[ip1][it10]) add = 1;
		    }
		  }
		  if (add == 0 && ip) {
		    ipprev -= idirp;
		    Int_t ip10 = iphi + ipprev;
		    if (ip10 >= 0 && ip10 < nPhi && fFlags[ip10][it1] != 0) {
		      if (pixels.find(pair<Int_t,Int_t>(ip10,it1)) != pixels.end() 
			  && fCharges[ip1][it1] <= fCharges[ip10][it1]) add = 1;
		    }
		  }
		}
		
		if (!add) break;
		phiMean += ip1 * fCharges[ip1][it1];
		zMean += it1 * fCharges[ip1][it1];
		adcTot += fCharges[ip1][it1];
		//vecDig.push_back(fDigis[ip1][it1]);
		//Int_t id = clus->Sec(fDigis[ip1][it1]);
		//if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = fCharges[ip1][it1];
		//else mapIdQ[id] = mapIdQ[id] + fCharges[ip1][it1];
		//else mapIdQ[id] = TMath::Max (mapIdQ[id], fCharges[ip1][it1]);
		//if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = 1;
		//else mapIdQ[id] = mapIdQ[id] + 1; // number of digits with the same ID
		MpdEmcDigit *dig = (MpdEmcDigit*) fDigiArray->UncheckedAt(fDigis[ip1][it1]);
		map<Int_t,Float_t> contrib = dig->GetContrib();
		RedoId(contrib);
		for (map<Int_t,Float_t>::iterator mit = contrib.begin(); mit != contrib.end(); ++mit) {
		  if (mapIdQ.find(mit->first) == mapIdQ.end()) mapIdQ[mit->first] = mit->second;
		  else mapIdQ[mit->first] = mapIdQ[mit->first] + mit->second;
		}
		pixels.insert(pair<Int_t,Int_t>(ip1,it1));

		sum2t += it1 * it1 * fCharges[ip1][it1];
		sum2p += ip1 * ip1 * fCharges[ip1][it1];
	      }
	    }
	  }
	}
      }

      phiMean /= adcTot;
      zMean /= adcTot;
      //phiMean += (0.5 - ishift); //
      phiMean += (0.0 - ishift); //
      //zMean += 0.5; //

      Double_t xloc = phiMean; //fSecGeo->Pad2Xloc(phiMean,clus->Row());
      Int_t padID = 0; //fSecGeo->PadID(clus->GetSect() % fSecGeo->NofSectors(), clus->Row());
      Double_t yloc = 0; //fSecGeo->LocalPadPosition(padID).Y();
      Double_t zloc = zMean; //fSecGeo->TimeBin2Z(zMean);
      p3loc.SetXYZ(xloc, yloc, zloc);
      Double_t rmsZ = TMath::Sqrt (sum2t/adcTot - zMean*zMean);
      Double_t rmsX = TMath::Sqrt (sum2p/adcTot - phiMean*phiMean);
      //p3err[1] = fSecGeo->TimeBin2Z(rms); // to pass the value
      //p3err[1] = rms;
      //cout << " Result: " << nLocMax0 << " " << pixels.size() << " " << xloc << " " << zloc << endl;

      // Apply corrections
      TVector3 p3errCor(p3err);
      /*
      CorrectReco(p3loc, p3errCor, clus->GetNumPads(), adcTot);

      if (clus->GetSect() >= fSecGeo->NofSectors()) p3loc[2] = -p3loc[2];
      fSecGeo->Local2Global(clus->GetSect(), p3loc, p3glob);

      // Dip angle correction (for interaction point with Z = 0)
      Double_t dip = TMath::ATan2 (TMath::Abs(p3glob[2]),p3glob.Pt()) * TMath::RadToDeg();
      p3errCor[2] = TMath::Max (p3errCor[2], 0.116 - 0.0046 * dip + 0.00015 * dip * dip);
      p3errCor[2] = TMath::Min (p3errCor[2], 0.5);
      // Correct Z-coordinate for rows = 0 and 27
      if (clus->Row() == 0 || clus->Row() == 27) {
	Double_t zcor = 0;
	if (clus->Row() == 0) zcor = -0.011 + 0.002 * dip;
	else zcor = -0.029 + 0.005 * dip + 3.886e-5 * dip * dip;
	zcor = TMath::Max (zcor, 0.);
	zcor = TMath::Min (zcor, 0.6);
	p3loc[2] -= TMath::Sign (zcor, p3loc[2]);
	p3glob[2] -= TMath::Sign (zcor, p3glob[2]);
      }
      */
      //zMean -= (fFlags[0].size() / 2 - 0.5);
      //AZ zMean -= (fFlags[0].size() / 2 - 0.0);
      //FIXME phiMean *= fEmcGeo->GetAngleOfModule() * TMath::DegToRad();
      //zMean *= fEmcGeo->GetLengthOfModuleByZ(); // FIXME - Z-position calculation
      GetPhiTheta(phiMean, zMean); // get angular coordinates phi and theta at (near) the inner faces of towers 
      p3glob.SetXYZ(phiMean,0,zMean);

      MpdTpcHit* hit = new ((*fHitArray)[ihit++]) MpdTpcHit(padID, p3glob, p3errCor, iclus); 
      //hit->SetLayer(clus->GetSect()); // sector number
      set<Int_t>::iterator sit = fSecRows1.upper_bound(xloc);
      hit->SetLayer(std::distance(fSecRows1.begin(),sit)); // sector number
      hit->SetLocalPosition(p3loc); // point position
      hit->SetEnergyLoss(adcTot);
      //Int_t ireg = 0; //(clus->Row() < fSecGeo->NofRowsReg(0)) ? 0 : 1;
      Double_t step = 1; //fSecGeo->PadHeight(ireg);
      hit->SetStep(step);
      hit->SetModular(1); // modular geometry flag
      hit->SetPad(Int_t(phiMean));
      hit->SetBin(Int_t(zMean));
      hit->SetRMS(rmsX, 0);
      hit->SetRMS(rmsZ, 1);
      //if (nLocMax0 == 1) hit->SetNdigits(nDigis);
      if (pixels.size() == 0) hit->SetNdigits(nDigis);
      else hit->SetNdigits(pixels.size());

      // !!! Warning: FairLinks are not persistent !!! 
      //hit->AddLink(FairLink(MpdTpcHit::PointIndex, pointIndx));
      vector<Int_t> vecDig;
      // Store maximum 5 track IDs with the highest charge contribution
      multimap<Double_t,Int_t> mapDig;
      for (map<Int_t,Double_t>::iterator mit = mapIdQ.begin(); mit != mapIdQ.end(); ++mit) 
	mapDig.insert(pair<Double_t,Int_t>(-mit->second,mit->first));
      multimap<Double_t,Int_t>::iterator mit = mapDig.begin();
      while (mit != mapDig.end() && vecDig.size() < 5) {
      //while (mit != mapDig.end() && vecDig.size() < 1) {
	vecDig.push_back(mit->second);
	++mit;
      }
      Int_t ndig = vecDig.size();
      for (Int_t idig1 = 0; idig1 < ndig; ++idig1) {
	hit->AddLink(FairLink(MpdTpcHit::MCTrackIndex, vecDig[idig1], idig1)); // weight = idig
	hit->AddID(vecDig[idig1]);
      }
      //cout << " EMC IDs: " << hit->GetNLinks() << " " << *(vecDig.begin()) << " " << hit->GetLinksWithType(MpdTpcHit::MCTrackIndex).GetLink(0).GetIndex() << " " << hit->GetLinksWithType(MpdTpcHit::MCTrackIndex).GetLink(hit->GetNLinks()-1).GetIndex() << endl;

    } // for ( ; rit != localMax.rend();
  } // for (Int_t iclus = 0; iclus < nclus;
}
//*/
  
//__________________________________________________________________________

void MpdEmcClusterFinderAZ::GetPhiTheta(Double_t &phi, Double_t &theta)
{
  // Convert COG in units of bins to angles

  static Int_t first = 1, offset = 0;
  static TSpline3 *phiS, *theS;

  if (first) {
    // Get phi and theta angles of the tower centers at their inner face                                
    first = 0;
    const vector<Double_t> &phis = fEmcGeo->GetPhiRow();
    const vector<Double_t> &thes = fEmcGeo->GetThetaBox();
    const vector<Double_t> &rhos = fEmcGeo->GetRhoCenterBox();
    const vector<Double_t> &zs =   fEmcGeo->GetZCenterBox();

    Int_t nphi = phis.size();
    // Offset due to the fact that the 1'st sector starts at phi = -Phi_sec/2; 
    offset = nphi / (fEmcGeo->GetNsec()/2 - 1) / 2;
    Double_t *phia = new Double_t [nphi];
    Double_t *ind = new Double_t [nphi];

    for (Int_t j = 0; j < nphi; ++j) {
      phia[j] = phis[j];
      ind[j] = j;
    }
    phiS = new TSpline3("grs",ind,phia,nphi); // phi vs ind                     
    delete [] phia;
    delete [] ind;

    Int_t nthe = thes.size();
    Double_t *the = new Double_t [nthe];
    Double_t *ind1 = new Double_t [nthe];
    Double_t height = fEmcGeo->GetLengthBox(); // tower half-height
    
    for (Int_t j = nthe-1; j >= 0; --j) {
      Double_t rho = rhos[j];
      Double_t z = zs[j];
      Double_t theta1 = thes[j];
      if (j < nthe-1 && thes[j] <= thes[j+1]+0.1) theta1 = 180 - theta1;
      Double_t costhe = TMath::Cos(theta1*TMath::DegToRad());
      Double_t sinthe = TMath::Sin(theta1*TMath::DegToRad());
      rho -= height * sinthe;
      z -= height * costhe;
      the[j] = TMath::ATan2(rho,z) * TMath::RadToDeg();
      ind1[j] = j; // - nthe/2;
    }
    theS = new TSpline3("grs1",ind1,the,nthe); // theta vs ind                  
    delete [] the;
    delete [] ind1;
  }

  phi = phiS->Eval(phi-offset);
  if (phi > 180) phi -= 360;
  theta = theS->Eval(theta);

}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::PeakAndValley(const MpdTpc2dCluster* clus, multimap<Double_t,Int_t> &localMax, Int_t ishift)
{
  // Apply peak-and-valley cuts to remove some local maxima
  //*
  static const Double_t ratio = 1.3;
  multimap<Double_t,Int_t>::reverse_iterator rit = localMax.rbegin(), rit1;
  ++rit;

  for ( ; rit != localMax.rend(); ++rit) {
    Int_t idig = rit->second;
    Int_t ipad = clus->Col(idig) + ishift;
    Int_t itime = clus->Bkt(idig);
    if (fFlags[ipad][itime] <= 0) continue; // merged local max

    // Loop over higher peaks
    rit1 = localMax.rbegin();
    for ( ; rit1 != rit; ++rit1) {
      Int_t idig0 = rit1->second;
      Int_t ipad0 = clus->Col(idig0) + ishift;
      Int_t itime0 = clus->Bkt(idig0);
      if (fFlags[ipad0][itime0] <= 0) continue; // merged local max

      Int_t dpad = ipad - ipad0, dt = itime - itime0;
      Int_t i0 = itime0, i1 = itime, j0 = ipad0, j1 = ipad, intime = 1;
      if (TMath::Abs(dpad) > TMath::Abs(dt)) i0 = ipad0, i1 = ipad, j0 = itime0, j1 = itime, intime = 0;
      Int_t stepi = TMath::Sign(1,i1-i0), stepj = TMath::Sign(1,j1-j0);
      
      //Int_t valOk = 0;
      Int_t merge = 0;
      //if (TMath::Abs(dpad) <= 1 && TMath::Abs(dt) <= 1) merge = 1;
      //else {
      if (TMath::Abs(dpad) <= 1 && TMath::Abs(dt) <= 1) {
	if ( (fFlags[ipad-dpad][itime] && fCharges[ipad-dpad][itime] > fCharges[ipad][itime] / ratio) ||
	     (fFlags[ipad][itime-dt] && fCharges[ipad][itime-dt] > fCharges[ipad][itime] / ratio) ) merge = 1;
      } else {
	for (Int_t ii = i0 + stepi; ii != i1; ii += stepi) {
	  merge = 0;
	  for (Int_t jj = j0; jj != j1 + stepj; jj += stepj) {
	    if (TMath::Abs(jj-j0) > TMath::Abs(ii-i0)) continue;
	    if (TMath::Abs(jj-j1) > TMath::Abs(ii-i1)) continue;
	    Int_t ip = ii, it = jj;
	    if (intime) ip = jj, it = ii;
	    //if (fCharges[ip][it] < fCharges[ipad][itime] / ratio) { valOk = 1; break; }
	    if (fFlags[ip][it] && fCharges[ip][it] > fCharges[ipad][itime] / ratio) { merge = 1; break; }
	  }
	  //if (valOk) break;
	  if (!merge) break;
	}
      }
      //if (!valOk) fFlags[ipad][itime] = -fFlags[ipad][itime]; // not deep enough valley
      if (merge) fFlags[ipad][itime] = -fFlags[ipad][itime]; // not deep enough valley
    }
  }

  // Remove failed peaks
  multimap<Double_t,Int_t>::iterator it = localMax.begin();
  //for ( ; it != localMax.end(); ++it) cout << it->second << " ";
  //cout << endl;

  it = localMax.begin();
  for ( ; it != localMax.end(); ++it) {
    Int_t idig = it->second;
    Int_t ipad = clus->Col(idig) + ishift;
    Int_t itime = clus->Bkt(idig);
    if (fFlags[ipad][itime] > 0) continue;
    //cout << " Before: " << idig << " " << itime << " " << ipad << " " << it->first << ", ";
    localMax.erase(it);
    //cout << " After: " << it->first << endl;
  }
  //*/
}

//__________________________________________________________________________

void MpdEmcClusterFinderAZ::RedoId(map<Int_t,Float_t>& contrib)
{
  // Redo track ID numbering 
  // Take IDs of particles produced outside EMC

  static const Double_t rmin = fEmcGeo->GetRmin(), rmax = fEmcGeo->GetRmax(); 
  TVector3 vert;

  map<Int_t,Float_t> copy(contrib);
  contrib.clear();

  for (map<Int_t,Float_t>::iterator it = copy.begin(); it != copy.end(); ++it) {
    Int_t id = it->first, idm = -1;
    MpdMCTrack *mctr = (MpdMCTrack*) fMcTrArray->UncheckedAt(id);
    mctr->GetStartVertex(vert);
    Double_t rvert = vert.Pt();
    while (rvert > rmin && rvert < rmax) {
      // Born inside EMC - find ancestor
      idm = mctr->GetMotherId();
      if (idm < 0) break;
      mctr = (MpdMCTrack*) fMcTrArray->UncheckedAt(idm);
      mctr->GetStartVertex(vert);
      rvert = vert.Pt();
    }
    if (idm >= 0) id = idm;
    if (contrib.find(id) == contrib.end()) contrib[id] = it->second;
    else contrib[id] += it->second;
  }
}    

//__________________________________________________________________________
ClassImp(MpdEmcClusterFinderAZ)