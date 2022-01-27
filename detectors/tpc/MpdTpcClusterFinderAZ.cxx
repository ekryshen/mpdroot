//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class MpdTpcClusterFinderAZ
//      see MpdTpcClusterFinderAZ.h for details
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - 14-July-2015
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "MpdTpcClusterFinderAZ.h"

// Collaborating Class Headers --------
#include "MpdTpc2dCluster.h"
#include "MpdTpcDigit.h"
#include "MpdTpcHit.h"
#include "MpdTpcSectorGeo.h"
#include "TpcPoint.h"
//#include "TpcGas.h"
//#include "TpcPrimaryCluster.h"
//#include "LinearInterpolPolicy.h"

#include "FairRootManager.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TTree.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <math.h>
#include <set>
#include <vector>

using namespace std;

//__________________________________________________________________________

MpdTpcClusterFinderAZ::MpdTpcClusterFinderAZ() : FairTask("TPC Cluster finder AZ"), fPersistence(kFALSE)
{
   /*
   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas= new TpcGas(tpcGasFile, 130);
   std::cout<<*fGas<<std::endl;
   */
}

//__________________________________________________________________________

MpdTpcClusterFinderAZ::~MpdTpcClusterFinderAZ()
{
   // delete fGas;
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::FinishTask()
{
   for (Int_t i = 0; i < fgkNsec2; ++i) {
      // fPrimArray[i]->Delete();
      delete[] fDigiSet[i];
   }
}

//__________________________________________________________________________

InitStatus MpdTpcClusterFinderAZ::Init()
{

   // Create containers for digits
   fSecGeo = MpdTpcSectorGeo::Instance();
   // Int_t nRows = MpdTpcSectorGeo::Instance()->NofRows();
   Int_t nRows = fSecGeo->NofRows();
   for (Int_t i = 0; i < fgkNsec2; ++i) fDigiSet[i] = new set<Int_t>[nRows];

   // Get ROOT Manager
   FairRootManager *ioman = FairRootManager::Instance();

   if (ioman == 0) {
      Error("MpdTpcClusterFinderAZ::Init", "RootManager not instantiated!");
      return kERROR;
   }

   // Get input collection
   fDigiArray = (TClonesArray *)ioman->GetObject("MpdTpcDigit");

   if (fDigiArray == 0) {
      Error("TpcClusterFinderAZ::Init", "Array of digits not found!");
      return kERROR;
   }

   // Create and register output array
   // SetPersistence();
   fClusArray = new TClonesArray("MpdTpc2dCluster");
   ioman->Register("TpcCluster", "Tpc", fClusArray, fPersistence);
   fHitArray = new TClonesArray("MpdTpcHit");
   ioman->Register("TpcRecPoint", "Tpc", fHitArray, fPersistence);

   return kSUCCESS;
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::Exec(Option_t *opt)
{

   fClusArray->Delete();
   fHitArray->Delete();
   const Int_t nSec = fgkNsec2 / 2; // number of TPC readout sectors
   // Clear digi containers
   // Int_t nRows = MpdTpcSectorGeo::Instance()->NofRows();
   Int_t nRows = fSecGeo->NofRows();
   for (Int_t i = 0; i < fgkNsec2; ++i) {
      for (Int_t j = 0; j < nRows; ++j) fDigiSet[i][j].clear();
   }

   // Fill digi containers
   Int_t nDigis = fDigiArray->GetEntriesFast();
   cout << " Total number of digits: " << nDigis << endl;
   for (Int_t i = 0; i < nDigis; ++i) {
      MpdTpcDigit *digi = (MpdTpcDigit *)fDigiArray->UncheckedAt(i);
      Int_t        isec = digi->GetSector();
      Int_t        irow = digi->GetRow();
      fDigiSet[isec][irow].insert(i);
   }

   Int_t nSum = 0;
   // Loop over sectors
   for (Int_t isec = 0; isec < fgkNsec2; ++isec) {

      // Loop over padrows
      for (Int_t irow = 0; irow < nRows; ++irow) {
         if (fDigiSet[isec][irow].size() == 0) continue;
         // cout << " Sector, row, digits: " << isec << " " << irow << " " << fDigiSet[isec][irow].size() << endl;
         ProcessPadrow(isec, irow);
         nSum += fDigiSet[isec][irow].size();
      }
   }
   cout << " Control sum: " << nSum << endl;

   // Find hits
   FindHits();
   // FindHitsLocMax();
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::ProcessPadrow(Int_t isec, Int_t irow)
{
   // Process one padrow of a sector

   memset(fFlags, 0, sizeof(fFlags[0][0]) * fgkNpads * fgkNtimes);

   set<Int_t>::iterator it = fDigiSet[isec][irow].begin();

   for (; it != fDigiSet[isec][irow].end(); ++it) {
      MpdTpcDigit *digi = (MpdTpcDigit *)fDigiArray->UncheckedAt(*it);
      Int_t        ipad = digi->GetPad(), itime = digi->GetTimeBin();
      if (ipad >= fgkNpads) Fatal("ProcessPadrow", "Too few pads!!! %i %i", ipad, fgkNpads);
      if (itime >= fgkNtimes) Fatal("ProcessPadrow", "Too few time bins!!! %i %i", itime, fgkNtimes);
      fCharges[ipad][itime] = digi->GetAdc();
      // fDigis[ipad][itime] = *it;
      fDigis[ipad][itime] = digi->GetOrigin(); // store trackID with max charge contribution
      fFlags[ipad][itime] = 1;
   }

   // Find (pre)clusters
   Int_t nclus0 = fClusArray->GetEntriesFast(), nclus = nclus0;
   for (Int_t ipad = 0; ipad < fgkNpads; ++ipad) {
      for (Int_t itime = 0; itime < fgkNtimes; ++itime) {
         if (fFlags[ipad][itime] <= 0) continue;
         // New cluster
         MpdTpc2dCluster *clus = new ((*fClusArray)[nclus++]) MpdTpc2dCluster(irow, isec);
         clus->Insert(fDigis[ipad][itime], irow, ipad, itime, fCharges[ipad][itime]);
         clus->SetID(fDigis[ipad][itime]);     // trackID
         clus->SetErrY(fCharges[ipad][itime]); // ADC counts
         fFlags[ipad][itime] = -1;
         for (Int_t ip = 0; ip < 2; ++ip) {
            for (Int_t itt = 0; itt < 2; ++itt) {
               if (itt == ip) continue;
               Int_t ip1 = ipad + ip, it1 = itime + itt;
               if (ip1 >= fgkNpads) continue;
               if (it1 >= fgkNtimes) continue;
               if (fFlags[ip1][it1] <= 0) continue;
               NextPixel(clus, ip1, it1);
            }
         }
      }
   }
   // cout << " Found preclusters: " << nclus - nclus0 << endl;
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::NextPixel(MpdTpc2dCluster *clus, Int_t ipad, Int_t itime)
{
   // Add next pixel to the cluster

   clus->Insert(fDigis[ipad][itime], clus->Row(), ipad, itime, fCharges[ipad][itime]);
   clus->SetID(TMath::Min(clus->ID(), fDigis[ipad][itime]));                    // min trackID
   clus->SetErrY(TMath::Max(Double_t(clus->GetErrY()), fCharges[ipad][itime])); // max ADC counts
   fFlags[ipad][itime] = -1;
   for (Int_t ip = -1; ip < 2; ++ip) {
      for (Int_t it = -1; it < 2; ++it) {
         if (TMath::Abs(it) == TMath::Abs(ip)) continue;
         Int_t ip1 = ipad + ip, it1 = itime + it;
         if (ip1 < 0 || ip1 >= fgkNpads) continue;
         if (it1 < 0 || it1 >= fgkNtimes) continue;
         if (fFlags[ip1][it1] <= 0) continue;
         NextPixel(clus, ip1, it1);
      }
   }
}

//__________________________________________________________________________
/*

void MpdTpcClusterFinderAZ::FindHits()
{
  // Reconstruct hits (one hit per precluster)

  TVector3 p3loc, p3glob, p3err(0.05,0.0,0.1);
  Int_t nclus = fClusArray->GetEntriesFast();

  for (Int_t iclus = 0; iclus < nclus; ++iclus) {
    MpdTpc2dCluster* clus = (MpdTpc2dCluster*) fClusArray->UncheckedAt(iclus);
    Int_t nDigis = clus->NDigits();
    Double_t padMean = 0, timeMean = 0;

    for (Int_t idig = 0; idig < nDigis; ++idig) {
      padMean += clus->Col(idig) * clus->Adc(idig);
      timeMean += clus->Bkt(idig) * clus->Adc(idig);
    }
    padMean /= clus->GetADC();
    timeMean /= clus->GetADC();

    Double_t xloc = fSecGeo->Pad2Xloc(padMean,clus->Row());
    Int_t padID = fSecGeo->PadID(clus->GetSect() % fSecGeo->NofSectors(), clus->Row());
    Double_t yloc = fSecGeo->LocalPadPosition(padID).Y();
    Double_t zloc = fSecGeo->TimeBin2Z(timeMean);
    p3loc.SetXYZ(xloc, yloc, zloc);

    // Apply corrections
    TVector3 p3errCor(p3err);
    CorrectReco(p3loc, p3errCor, clus->GetNumPads(), clus->GetADC());

    fSecGeo->Local2Global(clus->GetSect(), p3loc, p3glob);
    if (clus->GetSect() >= fSecGeo->NofSectors()) p3glob[2] = -p3glob[2];

    MpdTpcHit* hit = new ((*fHitArray)[iclus]) MpdTpcHit(padID, p3glob, p3errCor, iclus);
    hit->SetLayer(clus->Row());
    hit->SetLocalPosition(p3loc); // point position
    hit->SetEnergyLoss(clus->GetADC());
    hit->SetStep(0.0);
    hit->SetModular(1); // modular geometry flag
    hit->SetPad(Int_t(padMean));
    hit->SetBin(Int_t(timeMean));

    // !!! Warning: FairLinks are not persistent !!!
    //hit->AddLink(FairLink(MpdTpcHit::PointIndex, pointIndx));
    for (Int_t idig = 0; idig < nDigis; ++idig)
      hit->AddLink(FairLink(MpdTpcHit::MCTrackIndex, clus->Sec(idig))); // trackID stored in Sec
    //cout << hit->GetNLinks();
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

void MpdTpcClusterFinderAZ::FindHits()
{
   // Reconstruct hits (find local maxima)

   /*
     fCharges[ipad][itime] = digi->GetAdc();
     //fDigis[ipad][itime] = *it;
     fDigis[ipad][itime] = digi->GetOrigin(); // store trackID with max charge contribution
     fFlags[ipad][itime] = 1;
     */
   TVector3 p3loc, p3glob, p3err(0.05, 0.0, 0.1);
   Int_t    nclus = fClusArray->GetEntriesFast(), ihit = 0;

   for (Int_t iclus = 0; iclus < nclus; ++iclus) {
      MpdTpc2dCluster *clus   = (MpdTpc2dCluster *)fClusArray->UncheckedAt(iclus);
      Int_t            nDigis = clus->NDigits();

      memset(fFlags, 0, sizeof(fFlags[0][0]) * fgkNpads * fgkNtimes);
      for (Int_t idig = 0; idig < nDigis; ++idig) {
         Int_t ipad            = clus->Col(idig);
         Int_t itime           = clus->Bkt(idig);
         fCharges[ipad][itime] = clus->Adc(idig);
         fFlags[ipad][itime]   = 1;
         fDigis[ipad][itime]   = idig;
      }

      // Exclude pads which are not local maxima
      for (Int_t idig = 0; idig < nDigis; ++idig) {
         Int_t ipad  = clus->Col(idig);
         Int_t itime = clus->Bkt(idig);
         for (Int_t ip = -1; ip < 2; ++ip) {
            for (Int_t it = -1; it < 2; ++it) {
               if (TMath::Abs(it) == TMath::Abs(ip)) continue; // exclude diagonals
               if (it == 0 && ip == 0) continue;
               Int_t ip1 = ipad + ip, it1 = itime + it;
               if (ip1 < 0 || ip1 >= fgkNpads) continue;
               if (it1 < 0 || it1 >= fgkNtimes) continue;
               if (fFlags[ip1][it1] == 0) continue;
               if (clus->Adc(idig) < fCharges[ip1][it1]) {
                  fFlags[ipad][itime] = -1;
                  break;
               }
            }
         }
      }
      multimap<Double_t, Int_t> localMax;
      for (Int_t idig = 0; idig < nDigis; ++idig) {
         Int_t ipad  = clus->Col(idig);
         Int_t itime = clus->Bkt(idig);
         if (fFlags[ipad][itime] <= 0) continue;
         localMax.insert(pair<Double_t, Int_t>(clus->Adc(idig), idig));
         // cout << clus->Col(idig) << " " << clus->Bkt(idig) << " " << clus->Adc(idig) << endl;
      }
      // cout << " Local max: " << clus->GetSect() << " " << clus->Row() << " " << localMax.size() << endl;

      // Remove local maxima not separated by valleys - Peak and Valley
      Int_t nLocMax0 = localMax.size();
      if (localMax.size() > 1) PeakAndValley(clus, localMax);

      multimap<Double_t, Int_t>::reverse_iterator rit = localMax.rbegin();
      // vector<Int_t> vecDig;
      map<Int_t, Double_t> mapIdQ;
      for (; rit != localMax.rend(); ++rit) {
         Int_t idig  = rit->second;
         Int_t ipad  = clus->Col(idig);
         Int_t itime = clus->Bkt(idig);
         if (fFlags[ipad][itime] <= 0) continue; // merged local max
         // vecDig.clear();
         mapIdQ.clear();
         set<pair<Int_t, Int_t>> pixels;

         Double_t padMean = 0, timeMean = 0, adcTot = 0, sum2t = 0, sum2p = 0;

         // Process simple cluster (only 1 local max)
         // if (localMax.size() == 1) {
         if (nLocMax0 == 1) {
            for (Int_t idig1 = 0; idig1 < nDigis; ++idig1) {
               Int_t ip = clus->Col(idig1);
               Int_t it = clus->Bkt(idig1);
               padMean += ip * fCharges[ip][it];
               timeMean += it * fCharges[ip][it];
               adcTot += fCharges[ip][it];
               // vecDig.push_back(fDigis[ip][it]);
               Int_t id = clus->Sec(fDigis[ip][it]);
               // if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = fCharges[ip][it];
               // else mapIdQ[id] = mapIdQ[id] + fCharges[ip][it];
               // else mapIdQ[id] = TMath::Max (mapIdQ[id], fCharges[ip][it]);
               if (mapIdQ.find(id) == mapIdQ.end())
                  mapIdQ[id] = 1;
               else
                  mapIdQ[id] = mapIdQ[id] + 1;

               sum2t += it * it * fCharges[ip][it];
               sum2p += ip * ip * fCharges[ip][it];
            }
         } else {
            // Process complex cluster - start from maximum and go outward (up to 5 steps),
            // adding pixels with respectively lower charges
            for (Int_t idirp = -1; idirp < 2; idirp += 2) {
               for (Int_t ip = 0; ip < 5; ++ip) {
                  if (idirp > 0 && ip == 0) continue;
                  Int_t ipsign = ip * idirp;
                  Int_t ip1    = ipad + ipsign;
                  if (ip1 < 0 || ip1 >= fgkNpads) break;

                  for (Int_t idirt = -1; idirt < 2; idirt += 2) {
                     for (Int_t it = 0; it < 5; ++it) {
                        if (idirt > 0 && it == 0) continue;
                        Int_t itsign = it * idirt;

                        Int_t it1 = itime + itsign;
                        if (it1 < 0 || it1 >= fgkNtimes) break;
                        if (fFlags[ip1][it1] == 0) continue;

                        Int_t add = 1;
                        if (ip || it) {
                           if (fFlags[ip1][it1] > 0) continue; // case when 2 local max next to each other on 1 diagonal
                           // Check gradient
                           add          = 0;
                           Int_t ipprev = ipsign, itprev = itsign;
                           if (it) {
                              itprev -= idirt;
                              Int_t it10 = itime + itprev;
                              if (it10 >= 0 && it10 < fgkNtimes && fFlags[ip1][it10] != 0) {
                                 if (pixels.find(pair<Int_t, Int_t>(ip1, it10)) != pixels.end() &&
                                     fCharges[ip1][it1] <= fCharges[ip1][it10])
                                    add = 1;
                              }
                           }
                           if (add == 0 && ip) {
                              ipprev -= idirp;
                              Int_t ip10 = ipad + ipprev;
                              if (ip10 >= 0 && ip10 < fgkNpads && fFlags[ip10][it1] != 0) {
                                 if (pixels.find(pair<Int_t, Int_t>(ip10, it1)) != pixels.end() &&
                                     fCharges[ip1][it1] <= fCharges[ip10][it1])
                                    add = 1;
                              }
                           }
                        }

                        if (!add) break;
                        padMean += ip1 * fCharges[ip1][it1];
                        timeMean += it1 * fCharges[ip1][it1];
                        adcTot += fCharges[ip1][it1];
                        // vecDig.push_back(fDigis[ip1][it1]);
                        Int_t id = clus->Sec(fDigis[ip1][it1]);
                        // if (mapIdQ.find(id) == mapIdQ.end()) mapIdQ[id] = fCharges[ip1][it1];
                        // else mapIdQ[id] = mapIdQ[id] + fCharges[ip1][it1];
                        // else mapIdQ[id] = TMath::Max (mapIdQ[id], fCharges[ip1][it1]);
                        if (mapIdQ.find(id) == mapIdQ.end())
                           mapIdQ[id] = 1;
                        else
                           mapIdQ[id] = mapIdQ[id] + 1; // number of digits with the same ID
                        pixels.insert(pair<Int_t, Int_t>(ip1, it1));

                        sum2t += it1 * it1 * fCharges[ip1][it1];
                        sum2p += ip1 * ip1 * fCharges[ip1][it1];
                     }
                  }
               }
            }
         }

         padMean /= adcTot;
         timeMean /= adcTot;

         Double_t xloc  = fSecGeo->Pad2Xloc(padMean, clus->Row());
         Int_t    padID = fSecGeo->PadID(clus->GetSect() % fSecGeo->NofSectors(), clus->Row());
         Double_t yloc  = fSecGeo->LocalPadPosition(padID).Y();
         Double_t zloc  = fSecGeo->TimeBin2Z(timeMean);
         p3loc.SetXYZ(xloc, yloc, zloc);
         Double_t rmsZ = TMath::Sqrt(sum2t / adcTot - timeMean * timeMean);
         Double_t rmsX = TMath::Sqrt(sum2p / adcTot - padMean * padMean);
         // p3err[1] = fSecGeo->TimeBin2Z(rms); // to pass the value
         // p3err[1] = rms;
         // cout << " Result: " << nLocMax0 << " " << pixels.size() << " " << xloc << " " << zloc << endl;

         // Apply corrections
         TVector3 p3errCor(p3err);
         CorrectReco(p3loc, p3errCor, clus->GetNumPads(), adcTot);

         if (clus->GetSect() >= fSecGeo->NofSectors()) p3loc[2] = -p3loc[2];
         fSecGeo->Local2Global(clus->GetSect(), p3loc, p3glob);

         // Dip angle correction (for interaction point with Z = 0)
         Double_t dip = TMath::ATan2(TMath::Abs(p3glob[2]), p3glob.Pt()) * TMath::RadToDeg();
         p3errCor[2]  = TMath::Max(p3errCor[2], 0.116 - 0.0046 * dip + 0.00015 * dip * dip);
         p3errCor[2]  = TMath::Min(p3errCor[2], 0.5);
         // Correct Z-coordinate for rows = 0 and 27
         // if ((clus->Row() == 0 || clus->Row() == 27) && rmsZ > 1.0) {
         if (clus->Row() == 0 && rmsZ > 1.0) {
            Double_t zcor = 0;
            if (clus->Row() == 0)
               zcor = -0.011 + 0.002 * dip;
            else
               zcor = -0.029 + 0.005 * dip + 3.886e-5 * dip * dip;
            zcor = TMath::Max(zcor, 0.);
            zcor = TMath::Min(zcor, 0.6);
            p3loc[2] -= TMath::Sign(zcor, p3loc[2]);
            p3glob[2] -= TMath::Sign(zcor, p3glob[2]);
         }

         MpdTpcHit *hit = new ((*fHitArray)[ihit++]) MpdTpcHit(padID, p3glob, p3errCor, iclus);
         hit->SetLayer(clus->Row());
         hit->SetLocalPosition(p3loc); // point position
         hit->SetEnergyLoss(adcTot);
         Int_t    ireg = (clus->Row() < fSecGeo->NofRowsReg(0)) ? 0 : 1;
         Double_t step = fSecGeo->PadHeight(ireg);
         hit->SetStep(step);
         hit->SetModular(1); // modular geometry flag
         hit->SetPad(Int_t(padMean));
         hit->SetBin(Int_t(timeMean));
         hit->SetRMS(rmsX, 0);
         hit->SetRMS(rmsZ, 1);
         if (nLocMax0 == 1)
            hit->SetNdigits(nDigis);
         else
            hit->SetNdigits(pixels.size());

         // !!! Warning: FairLinks are not persistent !!!
         // hit->AddLink(FairLink(MpdTpcHit::PointIndex, pointIndx));
         vector<Int_t> vecDig;
         // Store maximum 3 track IDs with the highest charge contribution
         multimap<Double_t, Int_t> mapDig;
         for (map<Int_t, Double_t>::iterator mit = mapIdQ.begin(); mit != mapIdQ.end(); ++mit)
            mapDig.insert(pair<Double_t, Int_t>(-mit->second, mit->first));
         multimap<Double_t, Int_t>::iterator mit = mapDig.begin();
         while (mit != mapDig.end() && vecDig.size() < 9) {
            // while (mit != mapDig.end() && vecDig.size() < 1) {
            vecDig.push_back(mit->second);
            ++mit;
         }
         Int_t ndig = vecDig.size();
         /*
         for (Int_t idig = 0; idig < ndig; ++idig)
      hit->AddLink(FairLink(MpdTpcHit::MCTrackIndex, clus->Sec(vecDig[idig]))); // trackID stored in Sec
         */
         for (Int_t idig1 = 0; idig1 < ndig; ++idig1) {
            hit->AddLink(FairLink(MpdTpcHit::MCTrackIndex, vecDig[idig], idig1)); // weight = idig
            hit->AddID(vecDig[idig1]);
         }
         // cout << hit->GetNLinks() << " " << *(vecDig.begin()) << " " <<
         // hit->GetLinksWithType(MpdTpcHit::MCTrackIndex).GetLink(0).GetIndex() << " " <<
         // hit->GetLinksWithType(MpdTpcHit::MCTrackIndex).GetLink(hit->GetNLinks()-1).GetIndex() << endl;

      } // for ( ; rit != localMax.rend();
   }    // for (Int_t iclus = 0; iclus < nclus;
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::PeakAndValley(const MpdTpc2dCluster *clus, multimap<Double_t, Int_t> &localMax)
{
   // Apply peak-and-valley cuts to remove some local maxima

   const Double_t                              ratio = 1.3;
   multimap<Double_t, Int_t>::reverse_iterator rit   = localMax.rbegin(), rit1;
   ++rit;

   for (; rit != localMax.rend(); ++rit) {
      Int_t idig  = rit->second;
      Int_t ipad  = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      if (fFlags[ipad][itime] <= 0) continue; // merged local max

      // Loop over higher peaks
      rit1 = localMax.rbegin();
      for (; rit1 != rit; ++rit1) {
         Int_t idig0  = rit1->second;
         Int_t ipad0  = clus->Col(idig0);
         Int_t itime0 = clus->Bkt(idig0);
         if (fFlags[ipad0][itime0] <= 0) continue; // merged local max

         Int_t dpad = ipad - ipad0, dt = itime - itime0;
         Int_t i0 = itime0, i1 = itime, j0 = ipad0, j1 = ipad, intime = 1;
         if (TMath::Abs(dpad) > TMath::Abs(dt)) i0 = ipad0, i1 = ipad, j0 = itime0, j1 = itime, intime = 0;
         Int_t stepi = TMath::Sign(1, i1 - i0), stepj = TMath::Sign(1, j1 - j0);

         // Int_t valOk = 0;
         Int_t merge = 0;
         // if (TMath::Abs(dpad) <= 1 && TMath::Abs(dt) <= 1) merge = 1;
         // else {
         if (TMath::Abs(dpad) <= 1 && TMath::Abs(dt) <= 1) {
            if (fCharges[ipad - dpad][itime] > fCharges[ipad][itime] / ratio ||
                fCharges[ipad][itime - dt] > fCharges[ipad][itime] / ratio)
               merge = 1;
         } else {
            for (Int_t ii = i0 + stepi; ii != i1; ii += stepi) {
               merge = 0;
               for (Int_t jj = j0; jj != j1 + stepj; jj += stepj) {
                  if (TMath::Abs(jj - j0) > TMath::Abs(ii - i0)) continue;
                  if (TMath::Abs(jj - j1) > TMath::Abs(ii - i1)) continue;
                  Int_t ip = ii, it = jj;
                  if (intime) ip = jj, it = ii;
                  // if (fCharges[ip][it] < fCharges[ipad][itime] / ratio) { valOk = 1; break; }
                  if (fCharges[ip][it] > fCharges[ipad][itime] / ratio) {
                     merge = 1;
                     break;
                  }
               }
               // if (valOk) break;
               if (!merge) break;
            }
         }
         // if (!valOk) fFlags[ipad][itime] = -fFlags[ipad][itime]; // not deep enough valley
         if (merge) fFlags[ipad][itime] = -fFlags[ipad][itime]; // not deep enough valley
      }
   }

   // Remove failed peaks
   multimap<Double_t, Int_t>::iterator it = localMax.begin();
   // for ( ; it != localMax.end(); ++it) cout << it->second << " ";
   // cout << endl;

   it = localMax.begin();
   for (; it != localMax.end(); ++it) {
      Int_t idig  = it->second;
      Int_t ipad  = clus->Col(idig);
      Int_t itime = clus->Bkt(idig);
      if (fFlags[ipad][itime] > 0) continue;
      // cout << " Before: " << idig << " " << itime << " " << ipad << " " << it->first << ", ";
      localMax.erase(it);
      // cout << " After: " << it->first << endl;
   }
}

//__________________________________________________________________________

void MpdTpcClusterFinderAZ::CorrectReco(TVector3 &p3loc, TVector3 &p3errCor, Int_t nPads, Double_t &adc)
{
   // Correct reconstructed coordinates and charge

   Double_t xsec = (p3loc.Y() + fSecGeo->GetMinY()) * TMath::Tan(fSecGeo->Dphi() / 2);
   Double_t xloc = p3loc.X(), xloc0 = xloc, edge = 0.0; // distance to sector boundary
   if (xloc < 0)
      edge = xloc + xsec;
   else
      edge = xloc - xsec;
   if (TMath::Abs(edge) > 2.0) return; // no corrections
   Double_t adc0 = adc, coef = 1.0;
   // Correct charge
   if (nPads == 1 && adc0 < 5000) {
      coef = -0.0065 + 0.00015 * adc0 - 2.528e-8 * adc0 * adc0 + 1.380e-12 * adc0 * adc0 * adc0;
      coef = TMath::Max(coef, 0.02) / 0.74 * 0.91;
      adc /= coef;
   } else if (nPads == 2 && adc0 < 10000) {
      coef = 0.0139 + 0.000201 * adc0 - 2.351e-8 * adc0 * adc0 + 9.750e-13 * adc0 * adc0 * adc0;
      coef = TMath::Max(coef, 0.05) / 0.74 * 0.91;
      adc /= coef;
   } else if (nPads == 3 && adc0 < 6500) {
      coef = -0.068 + 0.000333 * adc0 - 4.401e-8 * adc0 * adc0 + 1.834e-12 * adc0 * adc0 * adc0;
      coef = TMath::Max(coef, 0.02) / 0.74 * 0.91;
      adc /= coef;
   }

   if (TMath::Abs(edge) > 1.5 || nPads > 3) return; // no corrections
   // Correct coordinates
   if (nPads == 1) {
      xloc -= TMath::Sign(0.520, edge);
      Double_t corr = -0.15;
      if (adc0 < 50000)
         corr = 0.206 - 2.532e-5 * adc0 + 8.715e-10 * adc0 * adc0 - 1.610e-14 * adc0 * adc0 * adc0 +
                1.193e-19 * adc0 * adc0 * adc0 * adc0;
      if (edge < 0) corr = -corr;
      xloc -= corr;
      p3errCor[0] = 0.116;
   } else if (nPads == 2) {
      xloc -= TMath::Sign(0.200, edge);
      Double_t corr = -0.07;
      if (adc0 < 200000)
         corr = 0.374 - 1.046e-5 * adc0 + 9.595e-11 * adc0 * adc0 - 4.098e-16 * adc0 * adc0 * adc0 +
                6.811e-22 * adc0 * adc0 * adc0 * adc0;
      if (edge < 0) corr = -corr;
      xloc -= corr;
      p3errCor[0] = 0.110;
   } else if (nPads == 3) {
      xloc -= TMath::Sign(0.033, edge);
      Double_t corr = 0.01;
      if (adc0 < 200000)
         corr = 0.188 - 4.007e-6 * adc0 + 3.364e-11 * adc0 * adc0 - 1.265e-16 * adc0 * adc0 * adc0 +
                1.771e-22 * adc0 * adc0 * adc0 * adc0;
      if (edge < 0) corr = -corr;
      xloc -= corr;
      p3errCor[0] = 0.073;
   } else
      return;
   if (TMath::Abs(xloc) > xsec) xloc = TMath::Sign(xsec, xloc0);
   p3loc.SetX(xloc);
}

//__________________________________________________________________________

ClassImp(MpdTpcClusterFinderAZ);