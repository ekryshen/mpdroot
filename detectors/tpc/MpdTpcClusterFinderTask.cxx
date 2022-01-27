//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcClusterFinderTask
//
// Environment:
//      Software developed for the MPD at NICA.
//
// Author List:
//      Sergey Merts
//      Artem Basalaev
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "MpdTpcClusterFinderTask.h"
#include "MpdTpcHit.h"
#include "TpcCluster.h"
#include "TpcPoint.h"
#include "MpdTpcPeak.h"
#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "TClonesArray.h"
#include <TSystem.h>
#include <iostream>
#include "TaskHelpers.h"
#include <TGraph.h>
#include "MpdTpcSector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TFile.h"
#include "sys/timeb.h"
#include <sys/time.h>

using namespace std;
using namespace TMath;

ClassImp(MpdTpcClusterFinderTask);

static Float_t kOneOverSqrt12   = 1 / Sqrt(12);
static Float_t kErrorCorrection = 0.115;

static clock_t tStart  = 0;
static clock_t tFinish = 0;
static clock_t tAll    = 0;

MpdTpcClusterFinderTask::MpdTpcClusterFinderTask()
   : fDigits(nullptr), fHitsArray(nullptr), fMCPointArray(nullptr), fMCTracks(nullptr), fDigitsArray(nullptr),
     fHisto(nullptr), fPersistence(kTRUE), fPrintDebugInfo(kFALSE), isHistogramsInitialized(kFALSE), fMakeQA(kFALSE),
     fCalcResiduals(kFALSE), fitGamma(kFALSE), fNumOfPadsInRow(nullptr)
{
   inputBranchName  = "MpdTpcDigit";
   outputBranchName = "MpdTpcHit";

   string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas = new TpcGas(tpcGasFile, 130);
}

MpdTpcClusterFinderTask::~MpdTpcClusterFinderTask()
{
   if (isHistogramsInitialized) delete fHisto;
   delete fHitsArray;
   delete fGas;
}

InitStatus MpdTpcClusterFinderTask::Init()
{

   FairRootManager *ioman = FairRootManager::Instance();

   if (!ioman) {
      cout << "\n-E- [MpdTpcClusterFinderTask::Init]: RootManager not instantiated!" << endl;
      return kFATAL;
   }
   if (fCalcResiduals) {
      fMCPointArray = (TClonesArray *)ioman->GetObject("TpcPoint");
      fMCTracks     = (TClonesArray *)ioman->GetObject("MCTrack");
   }
   fDigits    = (TClonesArray *)ioman->GetObject(inputBranchName);
   fHitsArray = new TClonesArray(outputBranchName);
   ioman->Register(outputBranchName, "TPC", fHitsArray, fPersistence);

   TpcSector *sector = new TpcSector();
   zDrift            = sector->GetLength(); // cm
   nSect             = sector->GetNSectors();
   nRows             = sector->GetNumRows();
   fNumOfPadsInRow   = (UInt_t *)sector->GetArrayPadsInRow();
   nTimeBins         = sector->GetNTimeBins();
   pwIn              = sector->GetInnerPadWidth();
   pwOut             = sector->GetOuterPadWidth();
   phIn              = sector->GetInnerPadHeight();
   phOut             = sector->GetOuterPadHeight();
   nInRows           = sector->GetNumInnerRows();
   nOutRows          = sector->GetNumOuterRows();
   fSectInHeight     = sector->GetSectInnerHeight();
   r_min             = sector->GetRmin();
   delete sector;

   if (!isHistogramsInitialized && fMakeQA) {
      fHisto = new MpdTpcClusterFinderQAHistograms();
      fHisto->Initialize();
      isHistogramsInitialized = kTRUE;
   }

   fNoiseThreshold = 1000.0; // electrons
   fSpread         = 0.196;  // cm  // Value is given by TPC group
   fEventId        = 0;

   OnePeakCntr         = 0;
   MoreThenOnePeakCntr = 0;
   OneDigitHitCntr     = 0;
   AllCntr             = 0;
   ThreePadsCntr       = 0;

   // memory allocating for input array
   fDigitsArray = new DigOrigArray **[nRows];
   for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
      fDigitsArray[iRow] = new DigOrigArray *[fNumOfPadsInRow[iRow] * 2];
      for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
         fDigitsArray[iRow][iPad] = new DigOrigArray[nTimeBins];
         for (UInt_t iTime = 0; iTime < nTimeBins; ++iTime) {
            fDigitsArray[iRow][iPad][iTime].isOverlap = kFALSE;
            fDigitsArray[iRow][iPad][iTime].origin    = 0;
            fDigitsArray[iRow][iPad][iTime].origins.clear();
            fDigitsArray[iRow][iPad][iTime].signal = 0.0;
         }
      }
   }

   cout << "-I- MpdTpcClusterFinderTask: Initialization successful." << endl;
   return kSUCCESS;
}

void MpdTpcClusterFinderTask::Exec(Option_t *opt)
{

   cout << "MpdTpcClusterFinder::Exec started! Event #" << ++fEventId << endl;

   tStart = clock();

   // Reset output Array
   if (!fHitsArray) Fatal("MpdTpcClusterFinder::Exec)", "No FoundClustersArray");
   fHitsArray->Delete();

   // parallel processing for different sectors
#pragma omp parallel // num_threads(8)
   {
#pragma omp for
      for (UInt_t iSec = 0; iSec < nSect; ++iSec) {

         for (Int_t i = 0; i < fDigits->GetEntriesFast(); ++i) {
            MpdTpcDigit *fDigit = (MpdTpcDigit *)fDigits->At(i);
            if ((UInt_t)fDigit->GetSector() == iSec) {
               fDigitsArray[fDigit->GetRow()][fDigit->GetPad()][fDigit->GetTimeBin()].origin = fDigit->GetOrigin();
               fDigitsArray[fDigit->GetRow()][fDigit->GetPad()][fDigit->GetTimeBin()].signal = fDigit->GetAdc();
            }
         }

         for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
            vector<MpdTpcFoundHit *>  tpcHitList;  // output vector of TpcHits
            vector<MpdTpc2dCluster *> extClusters; // vector of extended clusters
            vector<MpdTpcPeak *>      peakList;    // vector of all peaks in cluster (for each cluster))
            Find2DClusters(&extClusters, fDigitsArray[iRow], iRow, iSec);
            for (UInt_t iClust = 0; iClust < extClusters.size(); ++iClust) {
               FindPeaksInCluster(extClusters.at(iClust), &peakList, fDigitsArray[iRow]);

               // mybeg
               // Peaks drawing
               if (fMakeQA) {
                  Int_t numSect = fHisto->NumSector_hist;
                  if (numSect == -1) {
                     if (iRow == fHisto->NumRow_hist) {
                        for (UInt_t iPeak = 0; iPeak < peakList.size(); iPeak++) {
                           MpdTpcPeak *peak = peakList.at(iPeak);
                           fHisto->_hXT_peak_row->Fill(peak->Col(), peak->MaxAdcBkt(), /*peak->SumADC()*/ 0.9);
                        }
                     }
                  } else {
                     if (iRow == fHisto->NumRow_hist && iSec == (UInt_t)numSect) {
                        for (UInt_t iPeak = 0; iPeak < peakList.size(); iPeak++) {
                           MpdTpcPeak *peak = peakList.at(iPeak);
                           fHisto->_hXT_peak_row->Fill(peak->Col(), peak->MaxAdcBkt(), /*peak->SumADC()*/ 0.9);
                        }
                     }
                  }
               }

               Float_t iSizeCell = 0.1; // for collection painting (cell of some coolection has its size)
               // myend
               vector<MpdTpcPeak *> collectedPeakList; // for each cluster
               if (peakList.empty()) {
                  CreateHit(peakList, extClusters.at(iClust), &tpcHitList, fDigitsArray[iRow]);
               } else {
                  while (!peakList.empty()) {
                     CollectPeaks(peakList, extClusters.at(iClust), &collectedPeakList);

                     // mybeg
                     // Collected Peaks drawing
                     if (fMakeQA) {
                        if (fHisto->NumSector_hist == -1) {
                           if (iRow == fHisto->NumRow_hist) {
                              for (UInt_t iPeak = 0; iPeak < collectedPeakList.size(); iPeak++) {
                                 MpdTpcPeak *cpeak = collectedPeakList.at(iPeak);
                                 if (cpeak->Chi2() > 0.0) {
                                    fHisto->_hXT_collected_peak_row->Fill(cpeak->Col(), cpeak->PeakTime(),
                                                                          /*peak->SumADC()*/ iSizeCell);

                                 } else {
                                    fHisto->_hXT_collected_peak_row->Fill(cpeak->Col(), cpeak->Mean() + cpeak->BktOff(),
                                                                          /*peak->SumADC()*/ iSizeCell);
                                 }
                              }
                           }
                        } else {
                           if (iRow == fHisto->NumRow_hist && iSec == (UInt_t)fHisto->NumSector_hist) {
                              for (UInt_t iPeak = 0; iPeak < collectedPeakList.size(); iPeak++) {
                                 MpdTpcPeak *cpeak = collectedPeakList.at(iPeak);
                                 if (cpeak->Chi2() > 0.0) {
                                    fHisto->_hXT_collected_peak_row->Fill(cpeak->Col(), cpeak->PeakTime(),
                                                                          /*peak->SumADC()*/ iSizeCell);
                                 } else {
                                    fHisto->_hXT_collected_peak_row->Fill(cpeak->Col(), cpeak->Mean() + cpeak->BktOff(),
                                                                          /*peak->SumADC()*/ iSizeCell);
                                 }
                              }
                           }
                        }
                        iSizeCell += 0.1; // increse in size for the next collection
                        if (iSizeCell > 0.7) iSizeCell = 0.1;
                     }
                     // myend

                     CreateHit(collectedPeakList, extClusters.at(iClust), &tpcHitList, fDigitsArray[iRow]);
                     for (UInt_t i = 0; i < collectedPeakList.size(); ++i) {
                        vector<MpdTpcPeak *>::iterator pkItr =
                           find(peakList.begin(), peakList.end(), collectedPeakList[i]);
                        delete *pkItr;
                        peakList.erase(pkItr);
                     }
                     /// fixed bug with creating multiple identical clusters when cluster is oneHitCluster AND
                     /// peakList.size()>1
                     if (extClusters.at(iClust)->GetNumPads() == 1 || extClusters.at(iClust)->GetNumTimeBins() == 1)
                        break;
                  }
               }
               collectedPeakList.clear();
               peakList.clear();
            }

            for (UInt_t iHit = 0; iHit < tpcHitList.size(); ++iHit) {
               MpdTpcFoundHit *hit = tpcHitList.at(iHit);

               const UInt_t  Sect    = hit->GetSect();
               const Float_t sectPhi = Sect * TwoPi() / 12;

               const Float_t lX = hit->GetLocalX(); // local X coordinate
               const Float_t lY = hit->GetLocalY(); // local Y coordinate
               const Float_t lZ = hit->GetLocalZ(); // local Z coordinate

               const Float_t gX = (lY + r_min) * Cos(sectPhi) - lX * Sin(sectPhi); // global X coordinate
               const Float_t gY = (lY + r_min) * Sin(sectPhi) + lX * Cos(sectPhi); // global Y coordinate
               const Float_t gZ = (Sect < nSect / 2) ? lZ : lZ * (-1);             // global Z coordinate

               hit->SetGlobalX(gX);
               hit->SetGlobalY(gY);
               hit->SetGlobalZ(gZ);

               const Float_t Q = hit->QADC(); // charge of hit

               if (fCalcResiduals) CalcResiduals(hit, sectPhi);

               Int_t      row   = hit->Cluster()->Row();
               UInt_t     sec12 = (iSec < (nSect / 2)) ? iSec : (iSec - (nSect / 2));
               UInt_t     padID = (sec12 | (row << 5));
               MpdTpcHit *outHit;
#pragma omp critical
               {
                  outHit = new ((*fHitsArray)[fHitsArray->GetEntriesFast()]) MpdTpcHit(
                     padID, TVector3(gX, gY, gZ), TVector3(hit->errX(), hit->errY(), hit->errZ()), hit->GetOrigin());
               }

               outHit->SetQ(Q);
               outHit->SetLocalXYZ(lX, lY, lZ);
               outHit->SetLayer(row);
               outHit->SetModular(1);

               if (fMakeQA) {
#pragma omp critical
                  {
                     fHisto->_hXY->Fill(lX, lY, Q);
                     fHisto->_hYZ_local->Fill(lZ, lY, Q);
                     fHisto->_hRZ_global->Fill(gZ, gY, Q);
                     fHisto->_hX->Fill(lX, Q);
                     fHisto->_hY->Fill(lY, Q);
                     fHisto->_hZ->Fill(lZ, Q);
                     fHisto->_hX_global->Fill(gX, Q);
                     fHisto->_hY_global->Fill(gY, Q);
                     fHisto->_hZ_global->Fill(gZ, Q);
                     fHisto->_hXY_global->Fill(gX, gY, Q);
                     fHisto->_hPeak->Fill(Q);
                     if ((UInt_t)hit->PadRow() < nInRows) {
                        fHisto->_hErrX_inner->Fill(hit->errX());
                        fHisto->_hErrZ_inner->Fill(hit->errZ());
                     } else {
                        fHisto->_hErrX_outer->Fill(hit->errX());
                        fHisto->_hErrZ_outer->Fill(hit->errZ());
                     }
                     fHisto->_hErrY->Fill(hit->errY());
                     fHisto->_h3D->Fill(gX, gY, gZ, Q);
                     if (hit->GetSect() == 3 && lY <= 1.0) fHisto->_hXT_clust_row1->Fill(lX, lZ, Q); // bmy

                     // mybeg
                     // Clusters and hits drawing (together - on one histogram and separately - on different histograms)
                     if (fHisto->NumSector_hist == -1) {
                        if ((UInt_t)hit->Cluster()->Row() == fHisto->NumRow_hist) {
                           for (UInt_t iDig = 0; iDig < hit->Cluster()->GetNumDigits(); ++iDig) {
                              fHisto->_hXT_clust_hit_row->Fill(hit->Cluster()->Col(iDig), hit->Cluster()->Bkt(iDig),
                                                               /*hit->Cluster()->Adc(iDig)*/ 0.001);
                              fHisto->_hXT_clust_row->Fill(hit->Cluster()->Col(iDig), hit->Cluster()->Bkt(iDig),
                                                           hit->Cluster()->Adc(iDig));
                           }
                           fHisto->_hXT_clust_hit_row->Fill(hit->PadCol(), hit->TimeBkt(), /*hit->QADC()*/ 1000.0);
                           fHisto->_hXT_hit_row->Fill(hit->PadCol(), hit->TimeBkt(), hit->QADC());
                        }
                     } else {
                        if ((UInt_t)hit->Cluster()->Row() == fHisto->NumRow_hist &&
                            hit->GetSect() == fHisto->NumSector_hist) {
                           for (UInt_t iDig = 0; iDig < hit->Cluster()->GetNumDigits(); ++iDig) {
                              fHisto->_hXT_clust_hit_row->Fill(hit->Cluster()->Col(iDig), hit->Cluster()->Bkt(iDig),
                                                               /*hit->Cluster()->Adc(iDig)*/ 0.001);
                              fHisto->_hXT_clust_row->Fill(hit->Cluster()->Col(iDig), hit->Cluster()->Bkt(iDig),
                                                           hit->Cluster()->Adc(iDig));
                           }
                           fHisto->_hXT_clust_hit_row->Fill(hit->PadCol(), hit->TimeBkt(), /*hit->QADC()*/ 1000.0);
                           fHisto->_hXT_hit_row->Fill(hit->PadCol(), hit->TimeBkt(), hit->QADC());
                        }
                     }
                     // myend
                     fHisto->_hNumOfDigitsInCluster->Fill(hit->Cluster()->GetNumDigits());
                     fHisto->_hSect->Fill(hit->Cluster()->GetSect(), hit->Cluster()->GetADC());
                     fHisto->_hNumOfPadsInCluster->Fill(hit->Cluster()->GetNumPads());
                     fHisto->_hNumOfTimeBinsInCluster->Fill(hit->Cluster()->GetNumTimeBins());
                  }
               }
            }
            if (fMakeQA && tpcHitList.size() != 0) fHisto->_hHitDistr->Fill(tpcHitList.size());

            for (UInt_t i = 0; i < extClusters.size(); ++i) delete extClusters.at(i);
            extClusters.clear();
            for (UInt_t i = 0; i < tpcHitList.size(); ++i) delete tpcHitList.at(i);
            tpcHitList.clear();

            for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
               for (UInt_t iTime = 0; iTime < nTimeBins; ++iTime) {
                  if (fDigitsArray[iRow][iPad][iTime].signal > 0.0) {
                     fDigitsArray[iRow][iPad][iTime].origin = -1;
                     fDigitsArray[iRow][iPad][iTime].origins.clear();
                     fDigitsArray[iRow][iPad][iTime].signal    = 0.0;
                     fDigitsArray[iRow][iPad][iTime].isOverlap = kFALSE;
                  }
               }
            }
         }
      } // for (UInt_t iSec = 0; iSec < nSect; ++iSec)
   }

   tFinish = clock();
   tAll    = tAll + (tFinish - tStart);
   if (fPrintDebugInfo) cout << "Number of Found Clusters = " << fHitsArray->GetEntriesFast() << endl;
   cout << "MpdTpcClusterFinder::Exec finished" << endl;
}

Bool_t MpdTpcClusterFinderTask::Find2DClusters(vector<MpdTpc2dCluster *> *extClusters, DigOrigArray **f2dArray,
                                               UInt_t row, UInt_t sec)
{

   UInt_t curDigit[3];
   Bool_t result    = kFALSE;
   UInt_t cluscount = 0;

   Bool_t **fADCMarks = new Bool_t *[fNumOfPadsInRow[row] * 2];
   for (UInt_t i = 0; i < fNumOfPadsInRow[row] * 2; ++i) {
      fADCMarks[i] = new Bool_t[nTimeBins];
      for (UInt_t j = 0; j < nTimeBins; ++j) fADCMarks[i][j] = kFALSE;
   }

   for (UInt_t pad = 0; pad < fNumOfPadsInRow[row] * 2; ++pad) {
      for (UInt_t tBin = 0; tBin < nTimeBins; ++tBin) {
         if (f2dArray[pad][tBin].signal < fNoiseThreshold) continue;

         curDigit[0] = row;
         curDigit[1] = pad;
         curDigit[2] = tBin;

         if (!fADCMarks[pad][tBin]) { // if current ADC has no mark yet

            MpdTpc2dCluster *cluster = new MpdTpc2dCluster(row, sec);
            fADCMarks[pad][tBin]     = kTRUE;
            cluster->SetID(cluscount++);

            if (!cluster->Insert(row, pad, tBin, f2dArray[pad][tBin].signal)) {
#pragma omp critical
               {
                  cout << "Failed to insert digit into cluster..." << endl;
               }
               return kFALSE;
            }
            // now that we have the first digit in the 2d cluster, the following
            //  will find the rest of the connected digits
            result = GetNextDigit(curDigit, cluster, f2dArray, fADCMarks);

            if (!result) return kFALSE;
            extClusters->push_back(cluster);
         }
      }
   }

   for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[row] * 2; ++iPad) delete[] fADCMarks[iPad];
   delete[] fADCMarks;
   return kTRUE;
}

Bool_t MpdTpcClusterFinderTask::GetNextDigit(UInt_t *currdig, MpdTpc2dCluster *Clus2d, DigOrigArray **fDigitsArray1,
                                             Bool_t **fADCMarks)
{

   UInt_t thisCol, thisRow, thisBkt;
   UInt_t nextCol, nextBkt;

   thisRow = currdig[0];
   thisCol = currdig[1];
   thisBkt = currdig[2];
   for (int i = -1; i < 2; ++i) {

      // check left-right bounds...
      if (thisCol == 0 && i < 0) continue;
      if (thisCol == fNumOfPadsInRow[thisRow] * 2 - 1 && i >= 0) continue;
      nextCol = thisCol + i; // look right, left, check if we're at the edge

      for (int j = -1; j < 2; ++j) {
         if (i == 0 && j == 0) continue;
         // check up-down bounds...
         if (thisBkt == 0 && j < 0) continue;
         if (thisBkt == nTimeBins - 1 && j >= 0) continue;
         nextBkt = thisBkt + j; // look down, up, check if we're at the edge

         if (fDigitsArray1[nextCol][nextBkt].signal < fNoiseThreshold) continue;
         UInt_t nextdig[3] = {thisRow, nextCol, nextBkt};

         if (!fADCMarks[nextCol][nextBkt]) {
            if (!Clus2d->Insert(thisRow, nextCol, nextBkt, fDigitsArray1[nextCol][nextBkt].signal)) {
#pragma omp critical
               {
                  cout << "Failed to insert digit into cluster..." << endl;
               }
               return kFALSE;
            }
            fADCMarks[nextCol][nextBkt] = kTRUE;
            GetNextDigit(nextdig, Clus2d, fDigitsArray1, fADCMarks);
         }
      }
   }
   return kTRUE;
}

void MpdTpcClusterFinderTask::FindPeaksInCluster(MpdTpc2dCluster *clust, vector<MpdTpcPeak *> *PeakList,
                                                 DigOrigArray **fOriginsArray)
{

   const UInt_t fPeakValleyRatio = 2; // maximum of peak / minimum of peak
   const Int_t  mincol           = clust->MinCol();
   const Int_t  minbkt           = clust->MinBkt();
   const Int_t  nbkt             = clust->GetNumTimeBins();
   const Int_t  ncol             = clust->GetNumPads();
   const Int_t  ndig             = clust->GetNumDigits();

   Float_t xyADC[ncol][nbkt]; // box around 2D cluster

   for (Int_t i = 0; i < ncol; ++i)
      for (Int_t j = 0; j < nbkt; ++j) xyADC[i][j] = 0.0;

   for (Int_t i = 0; i < ndig; ++i) // data extraction block.
      xyADC[clust->Col(i) - mincol][clust->Bkt(i) - minbkt] = clust->Adc(i);

   // now loop through each pad and extrapolate individual peaks

   /// TODO!!! Set proper values for our TPC
   /// minthresh is a minimum ADC required for peak to be created
   const Float_t adcmin    = 1000.0;
   const Float_t minthresh = adcmin + 2 * sqrt(adcmin); //???

   for (Int_t i = 0; i < ncol; ++i) {

      Float_t thresh = 0.;
      Float_t max    = 0.;
      Float_t min    = 1.e9;
      // Int_t imax = -1;
      // Int_t imin = -1;
      Int_t npeak = 0;

      MpdTpcPeak *peak = new MpdTpcPeak();
      peak->SetCol(i + mincol);
      peak->AttachCluster(clust);
      Bool_t FoundPeak = kFALSE;

      for (Int_t j = 0; j < nbkt; ++j) {

         Float_t adc = xyADC[i][j];

         if (adc > minthresh) {
            if (peak->NSamples() == 0) peak->SetBktOff(minbkt + j);
            peak->Insert(adc, i + mincol, j + minbkt);
            if (adc > max) {
               max = adc;
               // imax = i;
               thresh = max - 2 * sqrt(max); // ???
            }
         }

         if (FoundPeak && adc < min) {
            min = adc;
            // imin = i;
         }
         if (adc < thresh) FoundPeak = kTRUE;

         if ((FoundPeak && adc < minthresh) ||
             // ^-- here we're below our min. threshold, so we save
             // the peak info. and create a new peak
             (FoundPeak && adc > min)) {
            //   ^-- here we've found a peak but the slope is rising... could be
            // entering a new peak, so save the peak info. and create a new one

            if (adc > min) // remove current ADC hit from peak
               peak->Remove(peak->NSamples() - 1);

            FoundPeak = kFALSE;
            if (peak->NSamples() > 1 && (peak->Max() / peak->Min() > fPeakValleyRatio)) {
               peak->SetOrigin(fOriginsArray[peak->MaxAdcCol()][peak->MaxAdcBkt()].origin);
               PeakList->push_back(peak);
               ++npeak;
               peak = new MpdTpcPeak();
            } else {
               peak->Clear();
            }

            peak->SetCol(i + mincol);
            peak->SetBktOff(minbkt + j);
            peak->AttachCluster(clust);
            if (adc > minthresh) peak->Insert(adc, i + mincol, j + minbkt);

            thresh = 0.;
            max    = 0.;
            min    = 1.e9;
            // imax = -1;
            // imin = -1;

         } // end check for 2nd peak (or clear end of 1st peak)
      }    // end loop through buckets in this column

      // need to add check that _last_ peak was not already added to list!

      if (FoundPeak && peak->NSamples() > 1 && (peak->Max() / peak->Min() > fPeakValleyRatio)) {
         // check that we haven't already added this peak, just to be safe...
         vector<MpdTpcPeak *>::iterator pkItr = find(PeakList->begin(), PeakList->end(), peak);
         if (pkItr == PeakList->end()) {
            peak->SetOrigin(fOriginsArray[peak->MaxAdcCol()][peak->MaxAdcBkt()].origin);
            PeakList->push_back(peak);
            ++npeak;
         }
      }

      if (npeak == 0) { // no peaks were found, so we'll store the ADC info
         // as a single "peak" anyway...
         peak->Clear();
         peak->SetCol(i + mincol);
         peak->AttachCluster(clust);

         for (Int_t j = 0; j < nbkt; ++j) {
            if (xyADC[i][j] > minthresh) {
               if (peak->NSamples() == 0) peak->SetBktOff(minbkt + j);
               peak->Insert(xyADC[i][j], i + mincol, j + minbkt);
            }
         }

         if (peak->NSamples() > 1) {
            peak->SetOrigin(fOriginsArray[peak->MaxAdcCol()][peak->MaxAdcBkt()].origin);
            PeakList->push_back(peak);
         }
      }

      // finally, check that peak has not been added to the peak list...
      // if not, then delete it!

      vector<MpdTpcPeak *>::iterator pkItr = find(PeakList->begin(), PeakList->end(), peak);
      if (pkItr == PeakList->end()) delete peak;

   } // end loop thru pads in this column
}

void MpdTpcClusterFinderTask::CollectPeaks(vector<MpdTpcPeak *> peakList, MpdTpc2dCluster *clust,
                                           vector<MpdTpcPeak *> *collectedPeakList)
{

   const Int_t mincol = clust->MinCol();
   const Int_t ncol   = clust->GetNumPads();
   Int_t       npeaks[ncol];

   for (Int_t i = 0; i < ncol; ++i) npeaks[i] = 0;

   if (fitGamma)
      FitPeaks(peakList, npeaks, mincol);
   else {
      for (UInt_t i = 0; i < peakList.size(); ++i) {
         MpdTpcPeak *peak = peakList[i];
         peak->SetChi2(0.0); // TMP
      }
   }

   collectedPeakList->clear();
   collectedPeakList->push_back(peakList[0]);
   // collect all peaks that line-up with each other
   Float_t tprev, tnext, tdiff; // time of current and previous peaks and time shift between two peaks
   Int_t   cprev, cnext, cdiff; // pad of current and previous peaks and pad shift between two peaks

   if (peakList[0]->Chi2() > 0.)
      tprev = peakList[0]->PeakTime();
   else {
      // Float_t x1 = peakList[0]->NSamples();
      // Float_t dt1 = par[0]*(1 - exp(-(par[1] - x1)*(par[1] - x1) / par[2]));
      tprev = peakList[0]->Mean() + peakList[0]->BktOff(); // - dt1;
   }

   cprev = peakList[0]->Col();

   for (UInt_t i = 1; i < peakList.size(); ++i) {
      if (peakList[i]->Chi2() > 0.)
         tnext = peakList[i]->PeakTime();
      else {
         // Float_t x1 = peakList[i]->NSamples();
         // Float_t dt1 = par[0]*(1 - exp(-(par[1] - x1)*(par[1] - x1) / par[2]));
         tnext = peakList[i]->Mean() + peakList[i]->BktOff(); // - dt1;
      }
      tdiff = tnext - tprev;
      cnext = peakList[i]->Col();
      cdiff = cnext - cprev;

      if (fabs(tdiff) < 2. && abs(cdiff) == 1) {
         collectedPeakList->push_back(peakList[i]);
         cprev = cnext;
      }
   }
}

void MpdTpcClusterFinderTask::CreateHit(vector<MpdTpcPeak *> collectedPeakList, MpdTpc2dCluster *clust,
                                        vector<MpdTpcFoundHit *> *hitList, DigOrigArray **fOriginsArray)
{

   AllCntr++;
   if (clust->GetNumPads() == 3) ThreePadsCntr++;

   Bool_t oneHitCluster = kFALSE;
   if (clust->GetNumPads() == 1 || clust->GetNumTimeBins() == 1) oneHitCluster = kTRUE;
   if (!oneHitCluster && !collectedPeakList.empty()) {

      //        TF1* fitfcn = nullptr;
      Float_t sigmaX = -1.0;

      MpdTpcFoundHit *hit = new MpdTpcFoundHit(clust);

      // now calculate <t>, sig<t>, <x> and sig<x>:
      Float_t avgt = 0.0, avgx = 0.0;
      Float_t sumq = 0.0;
      Float_t sigt = 0.0, sigx = 0.0;
      Float_t w = 0.0, dw = 0.;
      Bool_t  GoodHit = kFALSE;

      if (collectedPeakList.size() == 1) { // handle the case where a single _clear_
         // peak is found, but no neighboring peaks ...

         OnePeakCntr++;

         MpdTpcPeak *peak  = collectedPeakList[0];
         Float_t     pTime = peak->PeakTime();
         Float_t     pInt  = peak->Integral();
         Int_t       pCol  = peak->Col();
         // if (peak->Chi2() > 0. && peak->Chi2() < 30. && peak->Max() > 100.) {

         hit->SetQFit(pInt);
         hit->SetSigQFit(peak->IntegSig());
         hit->SetType(MpdTpcFoundHit::kFitPeak);

         // now estimate sigx and sigt...
         // draw a 3x3 box around peak
         Float_t xyADC[3][5];
         Float_t xADC[3];
         memset(xyADC, 0, sizeof(xyADC));
         memset(xADC, 0, sizeof(xADC));
         int ncol = 1;

         for (int ih = 0; ih < clust->NDigits(); ++ih) {
            int ic = (Int_t)(clust->Col(ih) - pCol);
            int ib = (Int_t)(clust->Bkt(ih) - pTime);
            if (abs(ic) <= 1 && abs(ib) <= 2) {
               xyADC[ic + 1][ib + 2] = clust->Adc(ih);
               xADC[ic + 1] += xyADC[ic + 1][ib + 2];
            }
         }

         if (xADC[0] > 0) ++ncol;
         if (xADC[2] > 0) ++ncol;

         if (ncol > 1) {
            // first deal with x-position
            avgx = 0.33 * xADC[0] * (pCol - 1) + 0.33 * xADC[2] * (pCol + 1) + pInt * pCol;
            avgx /= (0.33 * (xADC[0] + xADC[2]) + pInt);
            sigx = 0.33 * xADC[0] * (pCol - 1 - avgx) * (pCol - 1 - avgx) +
                   0.33 * xADC[2] * (pCol + 1 - avgx) * (pCol + 1 - avgx) + pInt * (pCol - avgx) * (pCol - avgx);
            sigx /= (0.33 * xADC[0] + pInt + 0.33 * xADC[2]) * (ncol - 1);

            if (sigx < 0.0001) // this is completely unrealistic, so we assume something really bad happened here...
               sigx = kOneOverSqrt12;
            else
               sigx = sqrt(sigx) * kErrorCorrection;

            // now deal with y-position
            Float_t y1 = 0;
            for (int ib = -2; ib <= 2; ++ib) y1 += (Int_t)(pTime + ib) * xyADC[0][ib + 2];
            if (xADC[0] > 0) y1 /= xADC[0];

            Float_t y2 = 0;
            for (int ib = -2; ib <= 2; ++ib) y2 += (Int_t)(pTime + ib) * xyADC[2][ib + 2];
            if (xADC[2] > 0) y2 /= xADC[2];

            avgt = 0.33 * xADC[0] * y1 + 0.33 * xADC[2] * y2 + pTime * pInt;
            avgt /= (0.33 * (xADC[0] + xADC[2]) + pInt);

            if (fabs(pTime - avgt) < 1.) {
               // looks reasonable, so continue

               sigt = 0.33 * xADC[0] * (y1 - avgt) * (y1 - avgt) + 0.33 * xADC[2] * (y2 - avgt) * (y2 - avgt) +
                      pInt * (pTime - avgt) * (pTime - avgt);
               sigt /= (0.33 * xADC[0] + pInt + 0.33 * xADC[2]) * (ncol - 1);

               if (sigt < 0.001)
                  sigt = peak->SigMean() / peak->SumADC();
               else
                  sigt = sqrt(sigt);
               GoodHit = kTRUE;
            }
         }
         //}
      } else { // more than one peak found in cluster...
         MoreThenOnePeakCntr++;
         GoodHit = kTRUE;

         //            TH1F* h = new TH1F("h", "h", 50, 0, 50);

         for (UInt_t i = 0; i < collectedPeakList.size(); ++i) {
            MpdTpcPeak *peak = collectedPeakList[i];
            if (peak->Chi2() > 0.) { // this means the peak was successfully fit.
               w  = peak->Integral();
               dw = peak->IntegSig();
               hit->SetQFit(hit->QFit() + w);
               hit->SetSigQFit(hit->QFitSig() + dw * dw);
               hit->SetType(hit->Type() | MpdTpcFoundHit::kFitPeak);
               Float_t myt = peak->PeakTime();
               /// Correction for TPC MPD may be added here
               // if (fUseTCorr) myt -= fTCorrTable->TCorr(peak->Row(),peak->Col());
               avgt += myt * w;
               sigt += myt * myt * w;
            } else { // fit was not successfully fit
               w = peak->SumADC();
               hit->SetQADC(hit->QADC() + peak->SumADC());
               hit->SetType(hit->Type() | MpdTpcFoundHit::kWMPeak);
               hit->AddOrigin(peak->GetOrigin());
               // Float_t x1 = peak->NSamples();
               // Float_t dt1 = par[0]*(1 - exp(-(par[1] - x1)*(par[1] - x1) / par[2]));
               Float_t myt = peak->Mean() + peak->BktOff(); // - dt1;
               /// Correction for TPC MPD may be added here
               // 	  if (fUseTCorr) myt -= fTCorrTable->TCorr(peak->Row(),peak->Col());
               avgt += myt * w;
               sigt += myt * myt * w;
            }
            //                h->Fill(peak->Col(), w);
            avgx += peak->Col() * w;
            sigx += peak->Col() * peak->Col() * w;
            sumq += w;
         }
         //            if (h->GetEntries() > 3) {
         //                h->Fit("gaus", "wQ");
         //                fitfcn = h->GetFunction("gaus");
         //                if (fitfcn) sigmaX = fitfcn->GetParameter(2);
         //            }
         //            delete h;

         hit->SetSigQFit(Sqrt(hit->QFitSig()));
         avgt /= sumq;
         avgx /= sumq;
         sigt = (sigt / sumq - avgt * avgt);
         sigx = (sigx / sumq - avgx * avgx);

         if (sigt <= 0.0)
            sigt = kOneOverSqrt12; //(hitPeaks.size()) * kOneOverSqrt12;
         else
            sigt = Sqrt(sigt);

         if (sigx <= 0.0)
            sigx = kOneOverSqrt12; // sigx = (hitPeaks.size()) * kOneOverSqrt12;
         else
            sigx = Sqrt(sigx) * kErrorCorrection;
      }

      Float_t pos[3], dpos[3];

      if (GoodHit) {

         Float_t padW, padH;
         UInt_t  Row = (UInt_t)clust->Row();

         if (Row < nInRows) {
            padW   = pwIn;
            padH   = phIn;
            pos[1] = padH * ((Float_t)Row + 0.5); // y-coordinate of pad center
         } else {
            padW   = pwOut;
            padH   = phOut;
            pos[1] = fSectInHeight + ((Float_t)(Row - nInRows) + 0.5) * padH; // y-coordinate of pad center
         }
         pos[0] = padW * ((Float_t)avgx - (Float_t)fNumOfPadsInRow[Row] + 0.5);
         pos[2] = ((Float_t)avgt + 0.5) * (zDrift / nTimeBins);

         hit->SetPadCol(avgx);
         hit->SetTimeBkt(avgt);
         // adding errors
         const Float_t timeBktSize = zDrift / nTimeBins;
         // cout << sigx * padW << " " << sigmaX * padW << endl;
         if (sigmaX > 0.0)
            dpos[0] = sigmaX * padW; // dx
         else
            dpos[0] = sigx * padW;        // dx
         dpos[1] = padH * kOneOverSqrt12; // dy
         dpos[2] = sigt * timeBktSize;    // dz

         hit->SetPos(pos, dpos);

         hitList->push_back(hit);
      } else
         delete hit;

   } else { // oneHitCluster is true  or collectedPeakList is empty

      Int_t   iPad, iRow, iTime;
      Float_t y;

      // Initialize variables to 0
      Float_t avgx = 0.0, avgz = 0.0;
      Float_t rmsx = 0.0, rmsz = 0.0;
      Float_t adc = 0.0, sumadc = 0.0, invsum = 0.0;

      UInt_t ndig = (UInt_t)clust->GetNumDigits();
      Int_t  Origin{0};

      for (UInt_t i = 0; i < ndig; ++i) {

         iPad  = clust->Col(i);
         iTime = clust->Bkt(i);
         iRow  = clust->Row(i);
         adc   = clust->Adc(i);

         Int_t CurrDigOrigin = fOriginsArray[iPad][iTime].origin;
         Origin              = CurrDigOrigin;

         // y-coordinate of pad center
         y = ((UInt_t)iRow < nInRows) ? (phIn * ((Float_t)iRow + 0.5))
                                      : (fSectInHeight + ((Float_t)(iRow - nInRows) + 0.5) * phOut);

         // Compute weighted sums of x, t, x^2, and t^2
         avgx += adc * iPad;
         avgz += adc * iTime;
         rmsx += adc * iPad * iPad;
         rmsz += adc * iTime * iTime;
         sumadc += adc;
      }
      // Make sure that we had at least one digit with charge
      if (sumadc > 0.0) {
         invsum = 1.0 / sumadc;
         avgx *= invsum;
         avgz *= invsum;
         rmsx = rmsx * invsum - avgx * avgx;
         rmsz = rmsz * invsum - avgz * avgz;
      }

      if (clust->GetNumPads() == 1 || rmsx <= 0.0) {
         Float_t xErrOnePad = (y < nInRows * phIn) ? (pwIn * kOneOverSqrt12) : (pwOut * kOneOverSqrt12);
         clust->SetErrX(xErrOnePad);
      } else {
         Float_t xErr = (y < nInRows * phIn) ? (pwIn * Sqrt(rmsx)) : (pwOut * Sqrt(rmsx));
         clust->SetErrX(xErr);
      }

      if (rmsz <= 0.0)
         clust->SetErrZ(zDrift / nTimeBins * kOneOverSqrt12);
      else
         clust->SetErrZ(Sqrt(rmsz) * zDrift / nTimeBins);

      /// mistake fixed!
      // if (clust->GetSect() > nSect / 2) avgz *= -1;

      clust->SetX(avgx);
      clust->SetY(y);
      clust->SetZ(avgz);

      Float_t yErrOneRow = (y < nInRows * phIn) ? (phIn * kOneOverSqrt12) : (phOut * kOneOverSqrt12);
      clust->SetErrY(yErrOneRow);
      clust->SetADC(sumadc);

      Float_t pos[3];
      Float_t hitErr[3];

      const Float_t padW = ((UInt_t)clust->Row() < nInRows) ? pwIn : pwOut;
      pos[0]             = padW * ((Float_t)clust->GetX() - (Float_t)fNumOfPadsInRow[clust->Row()] + 0.5);

      pos[1] = y;
      pos[2] = ((Float_t)clust->GetZ()) * (zDrift / nTimeBins);

      // adding errors
      hitErr[0] = clust->GetErrX(); // dx
      hitErr[1] = clust->GetErrY(); // dy
      hitErr[2] = clust->GetErrZ(); // dz

      MpdTpcFoundHit *tmphit = new MpdTpcFoundHit(pos, hitErr);
      tmphit->AttachCluster(clust);
      tmphit->SetTimeBkt(clust->GetZ());
      tmphit->SetPadCol(clust->GetX());
      tmphit->SetNumHits(1);
      tmphit->SetType(MpdTpcFoundHit::kWMPeak);
      tmphit->SetQADC(clust->GetADC());
      tmphit->AddOrigin(Origin);
      hitList->push_back(tmphit);

      OneDigitHitCntr++;
   }
}

void MpdTpcClusterFinderTask::FitPeaks(vector<MpdTpcPeak *> PeakList, Int_t npeaks[], Int_t mincol)
{

   Float_t par[3], epar[3], chi2;
   Float_t lnq[100], w[100], dt[100], lndt[100];

   for (UInt_t i = 0; i < PeakList.size(); ++i) {
      MpdTpcPeak *peak = PeakList[i];
      Int_t       ny   = peak->NSamples();

      // skip the fitting if this peak has too few bins...

      // Int_t fMinNBktGamma = 5; ///because this is the minimum number of bins in hit, must be re-checked!
      if (ny < 5) continue;

      Float_t yADC[ny];

      npeaks[peak->Col() - mincol] += 1;

      for (Int_t j = 0; j < ny; ++j) {
         yADC[j] = peak->Sample(j);
         lnq[j]  = log(yADC[j]);
         w[j]    = yADC[j] * yADC[j];
      }

      // search for good guess for t0
      Float_t t0, t0_save;
      t0_save         = -100.;
      Float_t minchi2 = 1.e30;

      for (Float_t k = -2.; k < 2.1; k += 0.5) {
         t0 = peak->BktOff() + k;
         for (Int_t j = 0; j < ny; ++j) {
            dt[j] = peak->BktOff() + (Float_t)j - t0;
            if (dt[j] < 1.e-3) dt[j] = 1.e-10;
            lndt[j] = Log(dt[j]);
         }
         LinGammaFit(lnq, w, dt, lndt, ny, par, epar, chi2);

         if (fPrintDebugInfo) {
#pragma omp critical
            {
               cout << "t0 = " << t0 << endl;
               for (int j = 0; j < 3; ++j) cout << "\tpar[" << j << "] = " << par[j] << endl;
               cout << "\t chi2 = " << chi2 << endl;
            }
         }

         if (chi2 > 0 && par[0] == par[0] && par[1] == par[1] && par[2] == par[2] && chi2 < minchi2) {
            minchi2 = chi2;
            t0_save = t0;
         }
      }

      if (fPrintDebugInfo)
#pragma omp critical
      {
         cout << "t0_save = " << t0_save << endl;
      }

      // check that we found a good guess for t0; if not, skip this peak
      if (t0_save < peak->BktOff() - 2.) continue;

      Float_t t0_start     = t0_save - 0.5;
      Float_t t0_stop      = t0_save + 0.6;
      t0_save              = -100.;
      minchi2              = 1.e30;
      Float_t par_save[3]  = {0., 0., 0.};
      Float_t epar_save[3] = {0., 0., 0.};
      Float_t chi2_save    = 0.;

      // now do a finer-binned search for the best t0...
      for (t0 = t0_start; t0 < t0_stop; t0 += 0.1) {
         for (Int_t j = 0; j < ny; ++j) {
            dt[j] = peak->BktOff() + (Float_t)j - t0;
            if (dt[j] < 1.e-3) dt[j] = 1.e-10;
            lndt[j] = log(dt[j]);
         }

         LinGammaFit(lnq, w, dt, lndt, ny, par, epar, chi2);

         if (fPrintDebugInfo) {
#pragma omp critical
            {
               cout << "t0 = " << t0 << endl;
               for (Int_t j = 0; j < 3; ++j) cout << "\tpar[" << j << "] = " << par[j] << endl;
               cout << "\t chi2 = " << chi2 << endl;
            }
         }

         if (chi2 > 0 && par[0] == par[0] && par[1] == par[1] && par[2] == par[2] && chi2 < minchi2) {
            minchi2 = chi2;
            t0_save = t0;
            for (int m = 0; m < 3; ++m) {
               par_save[m]  = par[m];
               epar_save[m] = epar[m];
               chi2_save    = chi2;
            }
         }
      }

      if (fPrintDebugInfo) {
#pragma omp critical
         {
            cout << "Final fit result:" << endl;
            cout << "\t t0 = " << t0_save << endl;
            for (int m = 0; m < 3; ++m)
               cout << "\t par[" << m << "] = " << par_save[m] << " +- " << epar_save[m] << endl;
            cout << "\t chi^2/ndf = " << chi2_save << endl;
         }
      }

      // now estimate the uncertainty on the integrated charge
      Float_t gam1  = Gamma(par_save[2] + 1);
      Float_t dgam1 = (Gamma(par_save[2] + 1 + 1.e-5) - Gamma(par_save[2] + 1 - 1.e-5)) / 2.e-5;
      Float_t Q     = par_save[0] * par_save[1] * gam1;
      Float_t dQ    = epar_save[0] * epar_save[0] / (par_save[0] * par_save[0]) +
                   epar_save[1] * epar_save[1] / (par_save[1] * par_save[1]) +
                   dgam1 * dgam1 * epar_save[2] * epar_save[2] / (gam1 * gam1);
      dQ = sqrt(fabs(dQ));
      dQ *= Q;

      peak->SetPeakTime(t0_save + par_save[1] * par_save[2]);
      //         peak->SetPeakTime(WeightedAverage(peak));

      peak->SetIntegral(Q);
      peak->SetIntegSig(dQ);
      peak->SetChi2(chi2_save);
   }
}

//.......................................................................
// FitCluster finds "peaks" in each pad of a cluster, then tries to
// fit these peaks to a Gamma distribution.  We then collect all peaks
// that are within 2 time buckets of each other to form a hit.  The hit
// time position is taken as the weighted mean of the peak positions, and
// the x-position is taken as the weighted pad number, where the weights
// in both cases are the total charge.  For peaks that were fit to Gamma,
// the charge is the integral of the Gamma Dist.  For peaks that were not
// fit to Gamma, the charge is set to the summed ADC values.
//.......................................................................

void MpdTpcClusterFinderTask::LinGammaFit(Float_t lnq[], Float_t w[], Float_t dt[], Float_t lndt[], Int_t nbin,
                                          Float_t par[], Float_t epar[], Float_t &chi2)
{

   Float_t a[3], b[3];
   Float_t C[3][3]; // C is symmetric
   Float_t CInv[3][3];

   for (Int_t i = 0; i < 3; ++i) {
      a[i] = b[i] = 0.;
      for (Int_t j = 0; j < 3; ++j) C[i][j] = CInv[i][j] = 0.;
   }

   Float_t sumlnq2 = 0;

   for (int i = 0; i < nbin; ++i) {
      sumlnq2 += w[i] * lnq[i] * lnq[i];
      b[0] += w[i] * lnq[i];
      b[1] += w[i] * dt[i] * lnq[i];
      b[2] += w[i] * lndt[i] * lnq[i];

      C[0][0] += w[i];
      C[0][1] += w[i] * dt[i];
      C[0][2] += w[i] * lndt[i];
      C[1][1] += w[i] * dt[i] * dt[i];
      C[1][2] += w[i] * dt[i] * lndt[i];
      C[2][2] += w[i] * lndt[i] * lndt[i];
   }

   C[1][0] = C[0][1];
   C[2][0] = C[0][2];
   C[2][1] = C[1][2];

   Float_t Cdet = C[0][0] * C[1][1] * C[2][2] + C[0][1] * C[1][2] * C[2][0] + C[0][2] * C[1][0] * C[2][1] -
                  C[0][2] * C[1][1] * C[2][0] - C[1][2] * C[2][1] * C[0][0] - C[2][2] * C[0][1] * C[1][0];

   CInv[0][0] = (C[1][1] * C[2][2] - C[1][2] * C[2][1]) / Cdet;
   CInv[0][1] = (C[0][2] * C[2][1] - C[0][1] * C[2][2]) / Cdet;
   CInv[0][2] = (C[0][1] * C[1][2] - C[0][2] * C[1][1]) / Cdet;

   CInv[1][0] = (C[1][2] * C[2][0] - C[1][0] * C[2][2]) / Cdet;
   CInv[1][1] = (C[0][0] * C[2][2] - C[0][2] * C[2][0]) / Cdet;
   CInv[1][2] = (C[0][2] * C[1][0] - C[0][0] * C[1][2]) / Cdet;

   CInv[2][0] = (C[1][0] * C[2][1] - C[1][1] * C[2][0]) / Cdet;
   CInv[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) / Cdet;
   CInv[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) / Cdet;

   for (Int_t i = 0; i < 3; ++i)
      for (Int_t j = 0; j < 3; ++j) a[i] += CInv[i][j] * b[j];

   par[2] = a[2];
   par[1] = -1. / a[1];
   par[0] = pow(par[1], par[2]) * exp(a[0]);

   chi2 = sumlnq2 - 2. * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) + a[0] * a[0] * C[0][0] + a[1] * a[1] * C[1][1] +
          a[2] * a[2] * C[2][2] + 2. * (a[0] * a[1] * C[0][1] + a[0] * a[2] * C[0][2] + a[1] * a[2] * C[1][2]);

   chi2 /= nbin - 2;

   Float_t A    = par[0];
   Float_t K    = a[2];
   Float_t tau  = par[1];
   Float_t tau2 = tau * tau;
   Float_t tau3 = tau * tau2;
   Float_t lnT  = log(par[1]);

   epar[0] = sqrt(fabs(A * A / (C[0][0] * (1 - a[0]) + b[0] - K * C[0][2] + C[0][1] / tau)));
   epar[1] = C[0][0] * ((K / tau2) * (K + a[0])) - C[0][1] * ((2 / (tau3)) * (a[0] + 3 * K / 2)) +
             C[0][2] * K * K / tau2 + C[1][1] / (tau2 * tau2) - 2 * K * C[1][2] / (tau * tau * tau) +
             2 * C[2][2] / (tau2 * tau2) - K * b[0] / tau2 + 2 * b[1] / tau3;

   epar[1] = sqrt(fabs(1. / epar[1]));
   epar[2] = sqrt(fabs(1. / (lnT * lnT * C[0][0] - lnT * C[0][2] + C[1][1])));
}

void MpdTpcClusterFinderTask::CalcResiduals(MpdTpcFoundHit *hit, Float_t sectPhi)
{
   Int_t   mcPntIdxZ = 0;
   Int_t   mcPntIdxX = 0;
   TGraph *gResX     = new TGraph();
   TGraph *gResZ     = new TGraph();

   Int_t         j, mctrackid;
   TpcPoint *    point;
   const Float_t gX         = hit->GetGlobalX();
   const Float_t gY         = hit->GetGlobalY();
   const Float_t gZ         = hit->GetGlobalZ();
   const Float_t lX         = hit->GetLocalX();
   const Float_t lY         = hit->GetLocalY();
   const Float_t lZ         = hit->GetLocalZ();
   const Bool_t  isInnerPad = Sqrt(gX * gX + gY * gY) < fSectInHeight + r_min ? kTRUE : kFALSE;
   const Float_t cosPhi     = Cos(sectPhi);
   const Float_t sinPhi     = Sin(sectPhi);
   for (j = 0; j < fMCPointArray->GetEntriesFast(); ++j) {
      point     = (TpcPoint *)fMCPointArray->At(j);
      mctrackid = point->GetTrackID();
      if (mctrackid == hit->GetOrigin()) {
         gResX->SetPoint(mcPntIdxX++, -point->GetX() * sinPhi + point->GetY() * cosPhi,
                         point->GetX() * cosPhi + point->GetY() * sinPhi - r_min);
         gResZ->SetPoint(mcPntIdxZ++, point->GetZ(), Sqrt(Power(point->GetX(), 2) + Power(point->GetY(), 2)));
      }
   }

   if (gResX->GetN() < 2 || gResZ->GetN() < 2) return;

   gResX->Fit("pol1", "Q");
   TF1 *         fitfcn_x = gResX->GetFunction("pol1");
   const Float_t bx       = fitfcn_x->GetParameter(0);
   const Float_t ax       = fitfcn_x->GetParameter(1);
   const Float_t xRes     = (lY - bx) / ax - lX;

   gResZ->Fit("pol1", "Q");
   TF1 *         fitfcn_z = gResZ->GetFunction("pol1");
   const Float_t bz       = fitfcn_z->GetParameter(0);
   const Float_t az       = fitfcn_z->GetParameter(1);
   const Float_t zRes     = (Sqrt(gX * gX + gY * gY) - bz) / az - gZ;

   fstream f;
   f.open(fNameResFile, fstream::out | fstream::app);
   f << isInnerPad << " " << xRes << " " << zRes << " " << hit->errX() << " " << hit->errZ() << " " << gX << " " << gY
     << " " << gZ << " " << lX << " " << lY << " " << lZ << " " << hit->QADC() << " " << hit->GetOrigin() << endl;
   f.close();

   delete gResX;
   delete gResZ;
}

void MpdTpcClusterFinderTask::CalcNewErrors(MpdTpcFoundHit *hit)
{
   Float_t sig_preamp2 = 0.0 * 0.0; // FIXME!!!
   Float_t sig_prf2    = fSpread * fSpread;
   Float_t Ldrift      = zDrift - hit->GetGlobalZ();
   Float_t Dt2         = fGas->Dt() * fGas->Dt();
   Float_t Dl2         = fGas->Dl() * fGas->Dl();
   Float_t tanA2       = 0.33;      // tg(30) * tg(30) ---> supremum
   Float_t tanB2       = 1.42;      // tg(55) * tg(55) ---> supremum???
   Float_t Lpad2       = 0.4 * 0.4; // TODO: get from geometry

   Float_t SigT = Sqrt(Dt2 * Ldrift + sig_preamp2 + tanA2 * Lpad2 / 12); // sigma time
   Float_t SigP = Sqrt(Dl2 * Ldrift + sig_prf2 + tanB2 * Lpad2 / 12);    // sigma pad

   hit->SetXerr(SigP * 0.4);
   hit->SetZerr(SigT * 0.4);
}

void MpdTpcClusterFinderTask::Finish()
{
   cout << "ClusterFinder work time = " << ((Float_t)tAll) / CLOCKS_PER_SEC << endl;
   if (fMakeQA) {
      toDirectory("QA/TPC");
      fHisto->Write();
      gFile->cd();
   }

   for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
      for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) delete[] fDigitsArray[iRow][iPad];
      delete[] fDigitsArray[iRow];
   }
   delete[] fDigitsArray;
   delete[] fNumOfPadsInRow;

   // mybeg
   /*
       TH2F *_thist = new TH2F("thist","thist_title",300,-150,150,300,-150,150);
       for(int i=0; i<fHitsArray->GetEntriesFast(); i++) {
           MpdTpcHit *fHit = (MpdTpcHit*) fHitsArray->At(i);
           _thist->Fill(fHit->GetX(), fHit->GetY(), fHit->GetQ());
           tfile2 <<i<<") "<<"XYZ="<<fHit->GetLocalX()<<", "<<fHit->GetLocalY()<<", "<<fHit->GetLocalZ()<<";  Layer=
      "<<fHit->GetLayer()<< ";  Pad="<<fHit->GetPad()<<";  Bin"<<fHit->GetBin()<<"\n";
       }

       _thist->Write("THIST", kOverwrite);

   */
   // myend
}