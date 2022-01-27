//-----------------------------------------------------------
//
// Description:
//      Implementation of class MpdTpcDigitizerTask
//      see MpdTpcDigitizerTask.h for details
//
// Author List:
//      Sergey Merts
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "MpdTpcDigitizerTask.h"

// C/C++ Headers ----------------------
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "MpdMultiField.h"
#include "FairRunAna.h"
#include "FairEventHeader.h"
#include "TpcPoint.h"
#include "TLorentzVector.h"
#include "FairRootManager.h"
#include "FairRunSim.h"
#include <TGeoManager.h>
#include <TRefArray.h>
#include "TClonesArray.h"
#include "TpcGas.h"
#include "TRandom.h"
#include "TMath.h"
#include "MpdTpcSector.h"
#include "TSystem.h"
#include "TaskHelpers.h"
#include "MpdTpcDigit.h"
#include "TFile.h"
// Class Member definitions -----------

using namespace std;
using namespace TMath;

static Int_t nOverlapDigit;
static Int_t nAllLightedDigits;

static clock_t tStart  = 0;
static clock_t tFinish = 0;
static clock_t tAll    = 0;

MpdTpcDigitizerTask::MpdTpcDigitizerTask()
   : fPersistence(kTRUE), fResponse(kTRUE), fDistribute(kTRUE), fAttach(kFALSE), fDiffuse(kFALSE), fDistort(kFALSE),
     fPrintDebugInfo(kFALSE), fIsHistogramsInitialized(kFALSE), fMakeQA(kFALSE), fHisto(nullptr), fPRF(nullptr),
     fNumOfPadsInRow(nullptr), fMCPointArray(nullptr), fMCTracksArray(nullptr), fDigits(nullptr), fSector(nullptr),
     fDigits4dArray(nullptr)
{
   fInputBranchName  = "TpcPoint";
   fOutputBranchName = "MpdTpcDigit";

   string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas = new TpcGas(tpcGasFile, 130);
}

MpdTpcDigitizerTask::~MpdTpcDigitizerTask()
{
   if (fIsHistogramsInitialized) {
      delete fHisto;
   }
   delete fGas;
   delete fPRF;
   delete fSector;
}

InitStatus MpdTpcDigitizerTask::Init()
{

   // Get ROOT Manager
   FairRootManager *ioman = FairRootManager::Instance();
   fMagField              = FairRunSim::Instance()->GetField();

   if (!ioman) {
      cout << "\n-E- [MpdTpcDigitizerTask::Init]: RootManager not instantiated!" << endl;
      return kFATAL;
   }
   fMCPointArray  = (TClonesArray *)ioman->GetObject(fInputBranchName);
   fMCTracksArray = (TClonesArray *)ioman->GetObject("MCTrack");

   fSector       = new TpcSector();
   nTimeBackets  = fSector->GetNTimeBins();
   nSectors      = fSector->GetNSectors();
   pwIn          = fSector->GetInnerPadWidth();
   pwOut         = fSector->GetOuterPadWidth();
   phIn          = fSector->GetInnerPadHeight();
   phOut         = fSector->GetOuterPadHeight();
   nRows         = fSector->GetNumRows();
   nInRows       = fSector->GetNumInnerRows();
   nOutRows      = fSector->GetNumOuterRows();
   fSectInHeight = fSector->GetSectInnerHeight();
   fSectHeight   = fSector->GetSectHeight();
   r_min         = fSector->GetRmin();
   zCathode      = fSector->GetLength(); // cm

   fNumOfPadsInRow = fSector->GetArrayPadsInRow();
   if (fPrintDebugInfo) {
      cout << "Number of pads in every rows is ";
      for (UInt_t k = 0; k < nRows; ++k) cout << fNumOfPadsInRow[k] * 2 << " ";
      cout << endl;
   }

   // memory allocating for output array
   fDigits4dArray = new DigOrigArray **[nRows];
   for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
      fDigits4dArray[iRow] = new DigOrigArray *[fNumOfPadsInRow[iRow] * 2];
      for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
         fDigits4dArray[iRow][iPad] = new DigOrigArray[nTimeBackets];
         for (UInt_t iTime = 0; iTime < nTimeBackets; ++iTime) {
            fDigits4dArray[iRow][iPad][iTime].isOverlap = kFALSE;
            fDigits4dArray[iRow][iPad][iTime].origin    = 0;
            fDigits4dArray[iRow][iPad][iTime].origins.clear();
            fDigits4dArray[iRow][iPad][iTime].signal = 0.0;
         }
      }
   }

   fDigits = new TClonesArray(fOutputBranchName);
   ioman->Register(fOutputBranchName, "TPC", fDigits, fPersistence);

   fNoiseThreshold = 1000.0; // electrons
   fGain           = 5000.0; // electrons
   if (fResponse) {
      fSpread = 0.196; // cm  // Value is given by TPC group
      k1      = 1.0 / (Sqrt(TwoPi()) * fSpread);
      k2      = -0.5 / fSpread / fSpread;
   } else {
      fSpread = 0.0; // cm  // FOR TEST ONLY. NO RESPONSE.
      k1 = k2 = 1.0;
   }

   if (!fIsHistogramsInitialized && fMakeQA) {
      fHisto = new MpdTpcDigitizerQAHistograms();
      fHisto->Initialize();
      fIsHistogramsInitialized = kTRUE;
   }
   fPRF              = padResponseFunction();
   nOverlapDigit     = 0;
   nAllLightedDigits = 0;

   cout << "-I- MpdTpcDigitizerTask: Initialization successful." << endl;
   return kSUCCESS;
}

void MpdTpcDigitizerTask::Exec(Option_t *opt)
{

   tStart = clock();

   cout << "MpdTpcDigitizer::Exec started" << endl;
   fDigits->Delete();

   Int_t nPoints = fMCPointArray->GetEntriesFast();
   if (nPoints < 2) {
      Warning("MpdTpcDigitizerTask::Exec", "Not enough Hits in TPC for Digitization (<2)");
      return;
   }

   if (fPrintDebugInfo) cout << "Number of MC points is " << nPoints << endl << endl;

   const Float_t phiStep = TwoPi() / nSectors * 2;

   for (UInt_t iSec = 0; iSec < nSectors; ++iSec) {

      // get memory for array
      TpcPoint *prePoint = (TpcPoint *)fMCPointArray->At(0);
      for (UInt_t i = 1; i < nPoints; i++) {
         TpcPoint *curPoint = (TpcPoint *)fMCPointArray->At(i);
         Float_t   globPhi  = ATan2(curPoint->GetY(), curPoint->GetX()); // angle in global coordinates
         if (globPhi < 0) globPhi += TwoPi();
         UInt_t curSectID = (UInt_t)(globPhi / phiStep + 0.5); // index of current sector
         if (curSectID == nSectors / 2) curSectID = 0;
         if (curPoint->GetZ() < 0.0) curSectID += (nSectors / 2);

         if (curSectID == iSec) { // is current point in sector
            TpcProcessing(prePoint, curPoint, iSec, i, nPoints);
         }
         prePoint = curPoint;
      }

      for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
         for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
            for (UInt_t iTime = 0; iTime < nTimeBackets; ++iTime) {
               if (fDigits4dArray[iRow][iPad][iTime].signal > fNoiseThreshold) {
                  Int_t outSize = fDigits->GetEntriesFast();
                  new ((*fDigits)[outSize]) MpdTpcDigit(CalcOrigin(fDigits4dArray[iRow][iPad][iTime]), iPad, iRow,
                                                        iTime, iSec, fDigits4dArray[iRow][iPad][iTime].signal);
               }
               if (fDigits4dArray[iRow][iPad][iTime].signal > 0.0) {
                  fDigits4dArray[iRow][iPad][iTime].origins.clear();
                  fDigits4dArray[iRow][iPad][iTime].origin    = -1;
                  fDigits4dArray[iRow][iPad][iTime].signal    = 0.0;
                  fDigits4dArray[iRow][iPad][iTime].isOverlap = kFALSE;
               }
            }
         }
      }
   }

   tFinish = clock();
   tAll    = tAll + (tFinish - tStart);
   cout << "MpdTpcDigitizer::Exec finished" << endl;
}

Int_t MpdTpcDigitizerTask::CalcOrigin(const DigOrigArray dig)
{

   Float_t max     = 0.0;
   Int_t   maxOrig = -1;
   if (dig.origins.size() > 1) {
      for (map<Int_t, Float_t>::const_iterator it = dig.origins.begin(); it != dig.origins.end(); ++it) {
         if (it->second > max) {
            maxOrig = it->first;
            max     = it->second;
         }
      }
   } else {
      maxOrig = dig.origins.begin()->first;
   }
   return maxOrig;
}

void MpdTpcDigitizerTask::PadResponse(Float_t x, Float_t y, UInt_t timeID, Int_t origin, DigOrigArray ***arr)
{

   vector<UInt_t>  lightedPads;
   vector<UInt_t>  lightedRows;
   vector<Float_t> amps;

   Float_t avAmp     = 0.0;
   Float_t amplSum   = 0.0;
   Float_t amplitude = 0.0;

   GetArea(x, y, fSpread * 3, lightedPads, lightedRows);

   for (UInt_t i = 0; i < lightedPads.size(); ++i) {
      amplitude = CalculatePadResponse(lightedPads.at(i), lightedRows.at(i), x, y);
      amps.push_back(amplitude);
      amplSum += amplitude;
   }

   if (amplSum > 0.0) {
      map<Int_t, Float_t>::iterator it;
      avAmp = fGain / amplSum; // Normalize amplitudes
      for (UInt_t i = 0; i < amps.size(); ++i) {
         arr[lightedRows.at(i)][lightedPads.at(i)][timeID].signal += (amps.at(i) * avAmp);
         it = arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.find(origin);
         if (it != arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.end()) {
            it->second += (amps.at(i) * avAmp);
         } else {
            arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.insert(
               pair<Int_t, Float_t>(origin, amps.at(i) * avAmp));
         }
      }
   }
}

void MpdTpcDigitizerTask::GetArea(Float_t xEll, Float_t yEll, Float_t radius, vector<UInt_t> &padIDs,
                                  vector<UInt_t> &rowIDs)
{

   Float_t padW = 0.0, padH = 0.0;
   Float_t y = 0.0, x = 0.0;
   UInt_t  pad = 0, row = 0;
   Float_t delta;
   if (fResponse)
      delta = 0.0;
   else
      delta = -1000.0; // for test only!!!

   if (yEll > nInRows * phIn + radius) {
      padW = pwOut;
      padH = phOut;

      x = xEll - radius - padW;
      do {
         x += padW;
         y = yEll - radius - padH;
         do {
            y += padH;
            row = (UInt_t)((y - fSectInHeight) / padH) + nInRows;
            pad = (x > 0.0) ? (fNumOfPadsInRow[row] + floor(x / padW)) : (fNumOfPadsInRow[row] - 1 + ceil(x / padW));
            if (row >= nRows || pad >= fNumOfPadsInRow[row] * 2) continue;
            padIDs.push_back(pad);
            rowIDs.push_back(row);
         } while (y < yEll + radius + delta);
      } while (x < xEll + radius + delta);
   } else if (yEll < nInRows * phIn - radius) {
      padW = pwIn;
      padH = phIn;

      x = xEll - radius - padW;
      do {
         x += padW;
         y = yEll - radius - padH;
         do {
            y += padH;
            row = (UInt_t)(y / padH);
            pad = (x > 0.0) ? (fNumOfPadsInRow[row] + floor(x / padW)) : (fNumOfPadsInRow[row] - 1 + ceil(x / padW));
            if (row >= nRows || pad >= fNumOfPadsInRow[row] * 2) continue;
            padIDs.push_back(pad);
            rowIDs.push_back(row);
         } while (y < yEll + radius + delta);
      } while (x < xEll + radius + delta);
   } else {
      x = xEll - radius - pwIn;
      do {
         x += pwIn;
         row = nInRows - 1;
         pad = (x > 0.0) ? (fNumOfPadsInRow[row] + floor(x / pwIn)) : (fNumOfPadsInRow[row] - 1 + ceil(x / pwIn));
         if (pad >= fNumOfPadsInRow[row] * 2) continue;
         padIDs.push_back(pad);
         rowIDs.push_back(row);
      } while (x < xEll + radius + delta);
      do {
         x += pwOut;
         row = nInRows;
         pad = (x > 0.0) ? (fNumOfPadsInRow[row] + floor(x / pwOut)) : (fNumOfPadsInRow[row] - 1 + ceil(x / pwOut));
         if (pad >= fNumOfPadsInRow[row] * 2) continue;
         padIDs.push_back(pad);
         rowIDs.push_back(row);
      } while (x < xEll + radius + delta);
   }
}

Float_t MpdTpcDigitizerTask::CalculatePadResponse(UInt_t padID, UInt_t rowID, Float_t x, Float_t y)
{

   Float_t padX, padY;
   Float_t padW, padH;
   if (rowID < nInRows) {
      padW = pwIn;
      padH = phIn;
      padX = padW * ((Float_t)padID - (Float_t)fNumOfPadsInRow[rowID] + 0.5); // x-coordinate of pad center
      padY = padH * ((Float_t)rowID + 0.5);                                   // y-coordinate of pad center
   } else {
      padW = pwOut;
      padH = phOut;
      padX = padW * ((Float_t)padID - (Float_t)fNumOfPadsInRow[rowID] + 0.5); // x-coordinate of pad center
      padY = fSectInHeight + ((Float_t)(rowID - nInRows) + 0.5) * padH;       // y-coordinate of pad center
   }

   const Float_t maxX = x - (padX - padW / 2);
   const Float_t minX = x - (padX + padW / 2);
   const Float_t maxY = y - (padY - padH / 2);
   const Float_t minY = y - (padY + padH / 2);

   const Float_t i1 = (Exp(k2 * minX * minX) + Exp(k2 * maxX * maxX)) * k1 * padW / 2;
   const Float_t i2 = (Exp(k2 * minY * minY) + Exp(k2 * maxY * maxY)) * k1 * padH / 2;

   return i1 * i2;
}

TF1 *MpdTpcDigitizerTask::padResponseFunction()
{
   if (fPRF) return fPRF;

   fPRF = new TF1("Gaus PRF", "gaus", -5, 5);
   fPRF->SetParameter(0, 1.0 / (sqrt(2.0 * TMath::Pi()) * fSpread));
   fPRF->SetParameter(1, 0);
   fPRF->SetParameter(2, fSpread);

   return fPRF;
}

Bool_t MpdTpcDigitizerTask::isSubtrackInInwards(const TpcPoint *p1, const TpcPoint *p2)
{ // WHAT AM I DOING???
   const Float_t x1   = p1->GetX();
   const Float_t x2   = p2->GetX();
   const Float_t y1   = p1->GetY();
   const Float_t y2   = p2->GetY();
   const Float_t a    = (y1 - y2) / (x1 - x2);
   const Float_t b    = (y1 * x2 - x1 * y2) / (x2 - x1);
   const Float_t minR = fabs(b) / sqrt(a * a + 1);

   if (minR < r_min) // then check if minimal distance is between our points
   {
      const Float_t x = -a * b / (a * a + 1);
      const Float_t y = b / (a * a + 1);
      if ((x1 - x) * (x2 - x) < 0 && (y1 - y) * (y2 - y) < 0) {
         return kTRUE;
      }
   }
   return kFALSE;
}

void MpdTpcDigitizerTask::TpcProcessing(const TpcPoint *prePoint, const TpcPoint *curPoint, const UInt_t secID,
                                        const UInt_t iPoint, const UInt_t nPoints)
{

   Float_t        dE       = 0.0;             // energy loss
   UInt_t         qTotal   = 0;               // sum of clusters charges (=sum of electrons between two TpcPoints)
   UInt_t         qCluster = 0;               // charge of cluster (= number of electrons)
   TLorentzVector curPointPos;                // coordinates for current TpcPoint
   TLorentzVector prePointPos;                // coordinates for previous TpcPoint
   TLorentzVector diffPointPos;               // steps for clusters creation
   TVector3       diffuse;                    // vector of diffuse for every coordinates
   TVector3       distort;                    // vector of distortion for every coordinates
   TLorentzVector electronPos;                // coordinates for created electrons
   TLorentzVector clustPos;                   // coordinates for created clusters
   Float_t        driftl = 0.0;               // length for drifting
   vector<UInt_t> clustArr;                   // vector of clusters between two TpcPoints
   Float_t        localX = 0.0, localY = 0.0; // local coordinates of electron (sector coordinates)

   if (fPrintDebugInfo && (iPoint % 1000 == 0))
      cout << UInt_t(iPoint * 1.0 / nPoints * 100.0) << " % of TPC points processed" << endl;
   //        curPoint = (TpcPoint*) fMCPointArray->At(i);

   if (fOnlyPrimery == kTRUE) {
      MpdMCTrack *tr = (MpdMCTrack *)fMCTracksArray->At(curPoint->GetTrackID());
      if (tr->GetMotherId() != -1) return;
   }
   // check if hits are on the same track
   if (curPoint->GetTrackID() == prePoint->GetTrackID() && !isSubtrackInInwards(prePoint, curPoint)) {

      dE = curPoint->GetEnergyLoss() * 1E9; // convert from GeV to eV
      if (dE < 0) {
         Error("MpdTpcDigitizerTask::Exec", "Negative Energy loss!");
         return;
      }

      curPointPos.SetXYZT(curPoint->GetX(), curPoint->GetY(), curPoint->GetZ(), curPoint->GetTime());
      prePointPos.SetXYZT(prePoint->GetX(), prePoint->GetY(), prePoint->GetZ(), prePoint->GetTime());
      if ((curPointPos.T() < 0) || (prePointPos.T() < 0)) {
         Error("MpdTpcDigitizerTask::Exec", "Negative Time!");
         return;
      }

      diffPointPos = curPointPos - prePointPos; // differences between two points by coordinates
      diffPointPos *=
         (1 / diffPointPos.Vect().Mag()); // directional cosines //TODO! Do we need this??? Look at line #297

      qTotal = (UInt_t)floor(fabs(dE / fGas->W()));

      // while still charge not used-up distribute charge into next cluster

      if (fDistribute) {
         while (qTotal > 0) {
            // roll dice for next cluster
            qCluster = fGas->GetRandomCSUniform();
            if (qCluster > qTotal) qCluster = qTotal;
            qTotal -= qCluster;
            clustArr.push_back(qCluster);
         } // finish loop for cluster creation
      } else {
         clustArr.push_back(qTotal); // JUST FOR TEST. NO CLUSTER DISTRIBUTION!
         //             clustArr.push_back(1); // JUST FOR TEST. NO CLUSTER DISTRIBUTION ONLY ONE ELECTRON IN CLUSTER!
      }

      diffPointPos *=
         (diffPointPos.Vect().Mag() /
          clustArr.size()); // now here are steps between clusters by coordinates TODO: correct distribution
      clustPos = prePointPos;
      for (UInt_t iClust = 0; iClust < clustArr.size(); ++iClust) {
         clustPos += diffPointPos;
         driftl = zCathode - fabs(clustPos.Z());
         for (UInt_t iEll = 0; iEll < clustArr.at(iClust); ++iEll) {

            // attachment
            if (fAttach)
               if (exp(-driftl * fGas->k()) < gRandom->Uniform()) continue; // FIXME

            // diffusion
            if (fDiffuse) {
               const Float_t sqrtDrift = sqrt(driftl);
               const Float_t sigmat    = fGas->Dt() * sqrtDrift;
               const Float_t sigmal    = fGas->Dl() * sqrtDrift;
               diffuse.SetXYZ(gRandom->Gaus(0, sigmat), gRandom->Gaus(0, sigmat), gRandom->Gaus(0, sigmal));
            }

            if (fDistort) {

               const Float_t dt  = 1E-03;   // time step [s]
               const Float_t mu  = 4.23;    // electron mobility [m^2 / s / V]
               const Float_t mu2 = mu * mu; // just square of mu

               const TVector3 E(0.0, 0.0, fGas->E() * 100); // vector of electric field components (now is constant and
                                                            // parallel to Z axes) // 100 - convert Ez from V/cm to V/m
               TVector3 B(0.0, 0.0, 0.0);                   // vector of magnetic field components
               TVector3 v;                                  // vector of current velocity components
               TVector3 posCur;                             // vector of current position components
               TVector3 posPre = clustPos.Vect();           // vector of previous position components
               TVector3 EBCross(0.0, 0.0, 0.0);             // vector product of E and B vectors

               Bool_t inTpc = kTRUE;
               while (inTpc) {

                  B.SetXYZ(fMagField->GetBx(posPre.X(), posPre.Y(), posPre.Z()) * 0.1,
                           fMagField->GetBy(posPre.X(), posPre.Y(), posPre.Z()) * 0.1,
                           fMagField->GetBz(posPre.X(), posPre.Y(), posPre.Z()) * 0.1);
                  EBCross = E.Cross(B);

                  v      = mu / (1 + mu2 * B.Mag2()) * (E - mu * EBCross + mu2 * B * (E * B));
                  posCur = v * dt + posPre;
                  //                            cout << "X = " << posCur.X() << " Y = " << posCur.Y() << " Z = " <<
                  //                            posCur.Z() << " Vx = " << v.X() << " Vy = " << v.Y() << " Vz = " <<
                  //                            v.Z() << endl;
                  if ((posCur.Perp() > r_min + fSectHeight) || (posCur.Perp() < r_min) || (Abs(posCur.Z()) > zCathode))
                     inTpc = kFALSE;
                  posPre = posCur;
               }
               distort.SetX(posCur.X() - clustPos.Vect().X());
               distort.SetY(posCur.Y() - clustPos.Vect().Y());
               distort.SetZ(0.0); // FIXME
            }

            electronPos.SetVect(clustPos.Vect() + diffuse + distort);
            electronPos.SetT(clustPos.T() +
                             (zCathode - fabs(electronPos.Z())) / fGas->VDrift()); // Do we need to use clustPos.T() ???

            const Float_t phiStep = TwoPi() / nSectors * 2;
            const Float_t sectPhi = secID * phiStep;

            localY = electronPos.X() * Cos(sectPhi) + electronPos.Y() * Sin(sectPhi) -
                     r_min; // converting from global to local coordinates
            localX = -electronPos.X() * Sin(sectPhi) +
                     electronPos.Y() * Cos(sectPhi); // converting from global to local coordinates
            if ((localY < 0.0) || (Abs(localX) > fSector->GetMaxX()) || (localY > fSectHeight)) continue; // FIXME!!!
            const Float_t timeStep  = (zCathode / nTimeBackets) / fGas->VDrift();
            const UInt_t  curTimeID = (UInt_t)((zCathode / fGas->VDrift() - electronPos.T()) / timeStep);
            if (curTimeID >= nTimeBackets) continue;
            if (fMakeQA) {
               fHisto->_hX_local->Fill(localX);
               fHisto->_hY_local->Fill(localY);
               fHisto->_hXY_local->Fill(localX, localY);
               fHisto->_hYZ_local->Fill(fabs(electronPos.Z()), localY);
               fHisto->_hXY_global->Fill(electronPos.X(), electronPos.Y());
               fHisto->_hRZ_global->Fill(electronPos.Z(), electronPos.Y());
               fHisto->_h3D_el->Fill(electronPos.X(), electronPos.Y(), electronPos.Z());
               fHisto->_hZ_local->Fill(fabs(electronPos.Z()));
               fHisto->_hX_global->Fill(electronPos.X());
               fHisto->_hY_global->Fill(electronPos.Y());
               fHisto->_hZ_global->Fill(electronPos.Z());
               fHisto->_hDiffuseXY->Fill(diffuse.X(), diffuse.Y());
               fHisto->_hDistortXY->Fill(distort.X(), distort.Y());
            }

            Int_t origin = prePoint->GetTrackID();
            PadResponse(localX, localY, curTimeID, origin, fDigits4dArray);
         }
      }
   } // end check for same track
   clustArr.clear();
   clustArr.resize(0);
}

void MpdTpcDigitizerTask::Finish()
{

   cout << "Digitizer work time = " << ((Float_t)tAll) / CLOCKS_PER_SEC << endl;

   if (fMakeQA) {
      toDirectory("QA/TPC");
      Float_t digit        = 0.0;
      UInt_t  iPad_shifted = 0; // needed for correct drawing of fDigitsArray

      for (UInt_t iRows = 0; iRows < nRows; ++iRows) {
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[iRows] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[iRows];
            for (UInt_t iTime = 0; iTime < nTimeBackets; ++iTime) {
               digit = fDigits4dArray[iRows][iPads][iTime].signal;
               fHisto->_hXY_dig->Fill(iPad_shifted, iRows, digit);
               //                        fHisto->_hSect_dig->Fill(iSect, digit);
               fHisto->_hX_dig->Fill(iPad_shifted, digit);
               fHisto->_hY_dig->Fill(iRows, digit);
               fHisto->_hZ_dig->Fill(iTime, digit);
               fHisto->_h3D_dig->Fill(iPad_shifted, iRows, iTime, digit);
               if (digit > 1000.0) fHisto->_hADC_dig->Fill(digit);
            }
         }
      }

      for (UInt_t iRows = 0; iRows < nRows; ++iRows) {
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[iRows] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[iRows];
            for (UInt_t iTime = 0; iTime < nTimeBackets; ++iTime) {
               digit = fDigits4dArray[iRows][iPads][iTime].signal;
               // pad activity
               // if (digit > 1000.0) {
               //    fHisto->_hXY_dig->Fill(iPad_shifted, iRows, 1.0);
               //}
               //                    fHisto->_hXY_dig->Fill(iPad_shifted, iRows, digit);
               fHisto->_h3D_dig->Fill(iPad_shifted, iRows, iTime, digit);
            }
         }
      }

      for (UInt_t iTime = 0; iTime < nTimeBackets; ++iTime) {
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[1] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[1];
            digit        = fDigits4dArray[1][iPads][iTime].signal;
            fHisto->_hXT_dig_1->Fill(iPad_shifted, iTime, digit);
         }
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[5] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[5];
            digit        = fDigits4dArray[5][iPads][iTime].signal;
            fHisto->_hXT_dig_5->Fill(iPad_shifted, iTime, digit);
         }
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[10] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[10];
            digit        = fDigits4dArray[10][iPads][iTime].signal;
            fHisto->_hXT_dig_10->Fill(iPad_shifted, iTime, digit);
         }
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[20] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[20];
            digit        = fDigits4dArray[20][iPads][iTime].signal;
            fHisto->_hXT_dig_20->Fill(iPad_shifted, iTime, digit);
         }
         for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[40] * 2; ++iPads) {
            iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[40];
            digit        = fDigits4dArray[40][iPads][iTime].signal;
            fHisto->_hXT_dig_40->Fill(iPad_shifted, iTime, digit);
         }
      }

      fHisto->Write();
      gFile->cd();
   }
}

ClassImp(MpdTpcDigitizerTask);
