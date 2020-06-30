//-----------------------------------------------------------
//
// Description:
//      Implementation of class MpdTpcDigitizerAZ
//      see MpdTpcDigitizerAZ.h for details
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 17-July-2015
//      (modified version of MpdTpcDigitizerTask.cxx)
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "MpdTpcDigitizerAZlt.h"

// MPD Headers ----------------------
#include "MpdTpcDigit.h"
#include "MpdMultiField.h"
#include "MpdTpcSector.h"
#include "MpdTpcSectorGeo.h"
#include "TpcGas.h"
#include "TpcPoint.h"

// FAIR Headers ----------------------
#include "FairRunAna.h"
#include "FairEventHeader.h"
#include "FairRootManager.h"
#include "FairRunSim.h"

// ROOT Headers ----------------------
#include "TClonesArray.h"
#include <TGeoManager.h>
#include "TLorentzVector.h"
#include "TRandom.h"
#include <TRefArray.h>
#include "TMath.h"
#include "TSystem.h"
#include "TaskHelpers.h"
#include <TVirtualFFT.h>

// C/C++ Headers ----------------------
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

// Class Member definitions -----------

using namespace std;
using namespace TMath;

static Int_t nOverlapDigit;
static Int_t nAllLightedDigits;

static clock_t tStart = 0;
static clock_t tFinish = 0;
static clock_t tAll = 0;

//FILE *lunAZ = NULL; //fopen("gasGain.dat","w");
//---------------------------------------------------------------------------

MpdTpcDigitizerAZlt::MpdTpcDigitizerAZlt()
  : FairTask("TPC digitizerAZlt"),
  fOnlyPrimary(kFALSE),
  fPersistence(kTRUE),
  fResponse(kTRUE),
  fDistribute(kTRUE),
  fAttach(kFALSE),
  fDiffuse(kTRUE),
  fDistort(kFALSE),
  fPrintDebugInfo(kFALSE),
  fIsHistogramsInitialized(kFALSE),
  fMakeQA(kFALSE),
  fHisto(NULL),
  fPRF(NULL),
  fNumOfPadsInRow(NULL),
  fMCPointArray(NULL),
  fMCTracksArray(NULL),
  fDigits(NULL),
  fSector(NULL),
  fDigits4dArray(NULL) 
  
{
  fInputBranchName = "TpcPoint";
  fOutputBranchName = "MpdTpcDigit";
  
  string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
  tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
  fGas = new TpcGas(tpcGasFile, 130);
}

//---------------------------------------------------------------------------
MpdTpcDigitizerAZlt::~MpdTpcDigitizerAZlt() 
{
  if (fIsHistogramsInitialized) {
    delete fHisto;
  }
  delete fGas;
  delete fPRF;
  delete fSector;
}

//---------------------------------------------------------------------------
InitStatus MpdTpcDigitizerAZlt::Init() 
{

  //Get ROOT Manager
  FairRootManager* ioman = FairRootManager::Instance();
  if (FairRunSim::Instance() == NULL)
    fMagField = FairRunAna::Instance()->GetField();
  else fMagField = FairRunSim::Instance()->GetField();
  //fMagField = FairRunSim::Instance()->GetField();
  
  if (!ioman) {
    cout << "\n-E- [MpdTpcDigitizerAZlt::Init]: RootManager not instantiated!" << endl;
    return kFATAL;
  }
  fMCPointArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
  fMCTracksArray = (TClonesArray*) ioman->GetObject("MCTrack");
  
  fSector = new TpcSector();
  /*
  nTimeBackets = fSector->GetNTimeBins();
  nSectors = fSector->GetNSectors();
  pwIn = fSector->GetInnerPadWidth();
  pwOut = fSector->GetOuterPadWidth();
  phIn = fSector->GetInnerPadHeight();
  phOut = fSector->GetOuterPadHeight();
  nRows = fSector->GetNumRows();
  nInRows = fSector->GetNumInnerRows();
  nOutRows = fSector->GetNumOuterRows();
  fSectInHeight = fSector->GetSectInnerHeight();
  fSectHeight = fSector->GetSectHeight();
  r_min = fSector->GetRmin();
  zCathode = fSector->GetLength(); //cm  
  */
  fSecGeo = MpdTpcSectorGeo::Instance();
  fNTimeBins = fSecGeo->GetNTimeBins();
  nSectors = fSecGeo->NofSectors() * 2;
  pwIn = fSecGeo->PadWidth(0);
  pwOut = fSecGeo->PadWidth(1);
  phIn = fSecGeo->PadHeight(0);
  phOut = fSecGeo->PadHeight(1);
  nRows = fSecGeo->NofRows();
  nInRows = fSecGeo->NofRowsReg(0);
  nOutRows = fSecGeo->NofRowsReg(1);
  fSectInHeight = fSecGeo->GetRocY(1) - fSecGeo->GetRocY(0);
  fSectHeight = fSecGeo->GetRocY(2) - fSecGeo->GetRocY(0);
  r_min = fSecGeo->GetMinY();
  zCathode = fSecGeo->GetZmax(); //cm  
  
  //fNumOfPadsInRow = fSector->GetArrayPadsInRow();
  fNumOfPadsInRow = fSecGeo->NPadsInRows();
  //if (fPrintDebugInfo) {
  if (1) {
    cout << "Number of pads in every rows is ";
    for (UInt_t k = 0; k < nRows; ++k)
      cout << fNumOfPadsInRow[k] * 2 << " ";
    cout << endl;
  }
  
  //memory allocating for output array
  fDigits4dArray = new DigOrigArray** [nRows];
  for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
    fDigits4dArray[iRow] = new DigOrigArray* [fNumOfPadsInRow[iRow] * 2];
    for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
      fDigits4dArray[iRow][iPad] = new DigOrigArray [fNTimeBins];
      for (UInt_t iTime = 0; iTime < fNTimeBins; ++iTime) {
	fDigits4dArray[iRow][iPad][iTime].isOverlap = kFALSE;
	//AZ fDigits4dArray[iRow][iPad][iTime].origin = 0;
	fDigits4dArray[iRow][iPad][iTime].origin = -1;
	fDigits4dArray[iRow][iPad][iTime].origins.clear();
	fDigits4dArray[iRow][iPad][iTime].signal = 0.0;
      }
    }
  }
  
  fDigits = new TClonesArray(fOutputBranchName);
  ioman->Register(fOutputBranchName, "TPC", fDigits, fPersistence);
  
  //AZ fNoiseThreshold = 1000.0; // electrons    
  fNoiseThreshold = 30.0; // ADC counts
  //AZ fGain = 5000.0; //electrons
  fGain = 1000.0; //electrons
  if (fResponse) {
    fSpread = 0.196; // cm  // Value is given by TPC group
    k1 = 1.0 / (Sqrt(TwoPi()) * fSpread);
    k2 = -0.5 / fSpread / fSpread;
  } else {
    fSpread = 0.0; // cm  // FOR TEST ONLY. NO RESPONSE.
    k1 = k2 = 1.0;
  }
  
  if (!fIsHistogramsInitialized && fMakeQA) {
    fHisto = new MpdTpcDigitizerQAHistograms();
    fHisto->Initialize();
    fIsHistogramsInitialized = kTRUE;
  }
  fPRF = padResponseFunction();
  nOverlapDigit = 0;
  nAllLightedDigits = 0;
  
  cout << "-I- MpdTpcDigitizerAZlt: Initialization successful." << endl;
  return kSUCCESS;
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::Exec(Option_t* opt) 
{

  tStart = clock();
        
  cout << "MpdTpcDigitizerAZlt::Exec started" << endl;
  fDigits->Delete();
  
  Int_t nPoints = fMCPointArray->GetEntriesFast();
  if (nPoints < 2) {
    Warning("MpdTpcDigitizerAZlt::Exec", "Not enough Hits in TPC for Digitization (<2)");
    return;
  }
  
  if (fPrintDebugInfo) cout << "Number of MC points is " << nPoints << endl << endl;
  
  const Float_t phiStep = TwoPi() / nSectors * 2;
  
  multimap<Double_t,Int_t>* pointsM = new multimap<Double_t,Int_t>[nSectors]; 
  for (Int_t ip = 0; ip < nPoints; ++ip) {
    TpcPoint* point = (TpcPoint*) fMCPointArray->UncheckedAt(ip);
    Float_t phi = ATan2(point->GetY(), point->GetX()); //angle in global coordinates
    if (phi < 0) phi += TMath::TwoPi();
    Int_t isec = (Int_t) (phi / phiStep + 0.5); //index of current sector
    if (isec == nSectors / 2) isec = 0;
    if (point->GetZ() < 0.0) isec += (nSectors / 2);
    pointsM[isec].insert(pair<Float_t,Int_t>(point->GetTime(),ip));
  }
  
  TpcPoint* virtPoint = new TpcPoint;
  TpcPoint tmpPoint;
  TpcPoint *ppp = &tmpPoint;

  for (UInt_t iSec = 0; iSec < nSectors; ++iSec) {
    
    if (pointsM[iSec].size() == 0) continue;
    multimap<Int_t,Int_t> pointID;
    multimap<Double_t,Int_t>::iterator mit = pointsM[iSec].begin();
    for ( ; mit != pointsM[iSec].end(); ++mit) {
      TpcPoint* point = (TpcPoint*) fMCPointArray->UncheckedAt(mit->second);
      pointID.insert(pair<Int_t,Int_t>(point->GetTrackID(),mit->second));
    }
    
    multimap<Int_t,Int_t>::iterator it = pointID.begin();
    TpcPoint* prePoint = (TpcPoint*) fMCPointArray->UncheckedAt(it->second);
    //cout << " prePoint: " << prePoint->GetTrackID() << " " << prePoint->GetX() << " " 
    //<< prePoint->GetY() << " " << prePoint->GetZ() << endl;
    Check4Edge(iSec, prePoint, virtPoint); // check for edge-effect
    if (prePoint != virtPoint) ++it;
    
    for ( ; it != pointID.end(); ++it) {
      TpcPoint* curPoint = (TpcPoint*) fMCPointArray->UncheckedAt(it->second);
      //if (curPoint->GetTrackID() != prePoint->GetTrackID()) Check4Edge(iSec, prePoint, virtPoint); // check for edge-effect 
      //if (curPoint->GetTrackID() != prePoint->GetTrackID()) {
      if (curPoint->GetTrackID() != prePoint->GetTrackID() || curPoint->GetLength() - prePoint->GetLength() > 1.1) { // 14.11.2017
	TVector3 pos;
	curPoint->Position(pos);
	tmpPoint.SetPosition(pos);
	curPoint->Momentum(pos);
	tmpPoint.SetMomentum(pos);
	tmpPoint.SetTrackID(curPoint->GetTrackID());
	ppp = &tmpPoint;
	Check4Edge(iSec, ppp, virtPoint); // check for edge-effect 
	prePoint = ppp;
	//cout << " --- " << curPoint->GetTrackID() << " " << curPoint->GetX() << " " << curPoint->GetY() 
	//   << " " << curPoint->GetZ() << endl;
	//cout << " --- " << prePoint->GetTrackID() << " " << prePoint->GetX() << " " << prePoint->GetY() 
	//   << " " << prePoint->GetZ() << endl;
      }
      TpcProcessing(prePoint, curPoint, iSec, it->second, nPoints);
      prePoint = curPoint;
    }
    
    //AZ Electronics response
    SignalShaping();
 
    Int_t maxTimeBin = MpdTpcSectorGeo::Instance()->TimeMax() / MpdTpcSectorGeo::Instance()->TimeBin() + 1;
    for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
      for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
	//AZ for (UInt_t iTime = 0; iTime < fNTimeBins; ++iTime) {
	for (UInt_t iTime = 0; iTime < maxTimeBin; ++iTime) {
	  if (fDigits4dArray[iRow][iPad][iTime].signal > fNoiseThreshold) {
	    Int_t outSize = fDigits->GetEntriesFast();
	    Int_t id = CalcOrigin(fDigits4dArray[iRow][iPad][iTime]);
	    if (id >= 0) new((*fDigits)[outSize]) MpdTpcDigit(id, iPad, iRow, iTime, iSec, fDigits4dArray[iRow][iPad][iTime].signal);
	  }
	  //if (fDigits4dArray[iRow][iPad][iTime].signal > 0.0) {
	  if (CalcOrigin(fDigits4dArray[iRow][iPad][iTime]) >= 0) {
	    fDigits4dArray[iRow][iPad][iTime].origins.clear();
	    fDigits4dArray[iRow][iPad][iTime].origin = -1;
	    fDigits4dArray[iRow][iPad][iTime].signal = 0.0;
	    fDigits4dArray[iRow][iPad][iTime].isOverlap = kFALSE;
	  }
	}
      }
    }
  } // for (UInt_t iSec = 0; iSec < nSectors;
  
  tFinish = clock();
  tAll = tAll + (tFinish - tStart);
  delete [] pointsM;
  delete virtPoint;
  cout << "MpdTpcDigitizerAZ::Exec finished" << endl;
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::Check4Edge(UInt_t iSec, TpcPoint* &prePoint, TpcPoint* virtPoint)
{
  // Check for edge-effect for track entering TPC from inside 
  // (and correct for it if necessary)

  TVector3 posG, posL;
  prePoint->Position(posG);
  Int_t row0 = MpdTpcSectorGeo::Instance()->Global2Local(posG, posL, iSec%(nSectors/2));
  row0 = MpdTpcSectorGeo::Instance()->PadRow(row0);
  //cout << " Row: " << row0 << " " << iSec << " " << posL[1] << endl;
  if (row0) return;

  // For padrow == 0:  create virtual point to correct for edge effect
  TVector3 mom, posL1;
  prePoint->Momentum(mom);
  if (mom.Pt() < 0.02) return; // do not adjust for very low-Pt tracks 
  if (posL[1] < 0.01) return; // do not adjust - almost at the entrance

  posG += mom;
  MpdTpcSectorGeo::Instance()->Global2Local(posG, posL1, iSec%(nSectors/2));
  mom = posL1;
  mom -= posL; // momentum in sector frame
  if (mom[1] < 0.02) return; // do not adjust - going inward or parallel to sector lower edge

  Double_t scale = mom[1] / posL[1];
  mom.SetMag(mom.Mag() / scale);
  posL -= mom;
  //cout << posL[0] << " " << posL[1] << " " << posL[2] << endl;
  MpdTpcSectorGeo::Instance()->Local2Global(iSec%(nSectors/2), posL, posG);
  virtPoint->SetPosition(posG);
  virtPoint->SetTrackID(prePoint->GetTrackID());
  prePoint->SetEnergyLoss(prePoint->GetEnergyLoss()*1.3); // 29.10.16 - correct for edge-effect
  prePoint = virtPoint;
}

//---------------------------------------------------------------------------

//AZ Int_t MpdTpcDigitizerAZ::CalcOrigin(const DigOrigArray dig) {
Int_t MpdTpcDigitizerAZlt::CalcOrigin(DigOrigArray& dig) 
{

  if (dig.origin >= 0) return dig.origin; // already done before
  if (dig.origins.size() == 0) return -1;
 
  Float_t max = 0.0;
  Int_t maxOrig = -1;
  if (dig.origins.size() > 1) {
    for (map<Int_t, Float_t>::const_iterator it = dig.origins.begin(); it != dig.origins.end(); ++it) {
      if (it->second > max) {
	maxOrig = it->first;
	max = it->second;
      }
    }
  } else {
    maxOrig = dig.origins.begin()->first;
  }
  dig.origin = maxOrig; //AZ
  return maxOrig;
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::PadResponse(Float_t x, Float_t y, UInt_t timeID, Int_t origin, DigOrigArray ***arr) 
{
    
  vector<UInt_t> lightedPads;
  vector<UInt_t> lightedRows;
  vector<Float_t> amps;

  Float_t avAmp = 0.0;
  Float_t amplSum = 0.0;
  Float_t amplitude = 0.0;

  GetArea(x, y, fSpread * 3, lightedPads, lightedRows);
  Double_t gain = fGain * Polya();
  //fprintf(lunAZ,"%f\n",gain/fGain);
  Int_t nPads = lightedPads.size();
    
  for (Int_t i = 0; i < nPads; ++i) {
    //AZ amplitude = CalculatePadResponse(lightedPads.at(i), lightedRows.at(i), x, y);
    amplitude = gain * CalculatePadResponse(i, nPads, lightedPads.at(i), lightedRows.at(i), x, y);
    amps.push_back(amplitude);
    amplSum += amplitude;
  }
  
  /*AZ
  if (amplSum > 0.0) {
    map<Int_t, Float_t>::iterator it;
    avAmp = fGain / amplSum; // Normalize amplitudes
    for (UInt_t i = 0; i < amps.size(); ++i) {
      arr[lightedRows.at(i)][lightedPads.at(i)][timeID].signal += (amps.at(i) * avAmp);
      it = arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.find(origin);
      if (it != arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.end()) {
	it->second += (amps.at(i) * avAmp);
      } else {
	arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.insert(pair<Int_t, Float_t>(origin, amps.at(i) * avAmp));
      }
    }
  }
  */
  if (amplSum > 0.0) {
    map<Int_t, Float_t>::iterator it;
    for (UInt_t i = 0; i < amps.size(); ++i) {
      arr[lightedRows.at(i)][lightedPads.at(i)][timeID].signal += amps.at(i);
      it = arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.find(origin);
      if (it != arr[lightedRows.at(i)][lightedPads.at(i)][timeID].origins.end()) {
	it->second += amps.at(i);
      } else {
	arr[lightedRows[i]][lightedPads[i]][timeID].origins.insert(pair<Int_t, Float_t>(origin, amps[i]));
      }
    }
  }
}

//---------------------------------------------------------------------------

Double_t MpdTpcDigitizerAZlt::Polya()
{
  // Gas gain according to Polya-distribution with parameter \theta = 1.5

  static Int_t first = 0;
  static TH1D *hPolya;

  //long rseed;
  Double_t step = 0.01, shift = 0.005, param = 1.5, prob, lambda;

  if (first == 0) {
    hPolya = new TH1D("hPolya","Polya distribution",1000,0,10);
    Double_t param1 = 1 + param;
    Double_t coef = TMath::Exp(TMath::Log(param1)*param1) / TMath::Gamma(param1);
    for (Int_t i = 0; i < 1000; ++i) {
      lambda = i * step + shift;
      //prob = param / 0.8862 * TMath::Sqrt(param*lambda)*TMath::Exp(-param*lambda); // 0.8862 = Gamma(1.5)
      prob = coef * TMath::Exp(TMath::Log(lambda)*param) * TMath::Exp(-param1*lambda);
      hPolya->Fill(lambda,prob);
    }
    first = 1;
  }

  return hPolya->GetRandom();
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::SignalShaping()
{
  // Apply electronics response function

  static Int_t first = 0, nbins = 0, icent = 0, n2 = 0;
  static Double_t *reFilt = NULL, *imFilt = NULL;
  static TVirtualFFT *fft[2] = {NULL,NULL};
  const Double_t sigma = 190./2/TMath::Sqrt(2*TMath::Log(2)), sigma2 = sigma * sigma; // FWHM = 190 ns

  if (first == 0) {
    first = 1;
    nbins = MpdTpcSectorGeo::Instance()->GetNTimeBins();
    if (nbins % 2 == 0) --nbins;
    n2 = nbins / 2 + 1;
    icent = nbins / 2;
    reFilt = new Double_t [nbins];
    imFilt = new Double_t [nbins];
    for (Int_t i = 0; i < nbins; ++i) {
      Double_t t = (i - icent) * MpdTpcSectorGeo::Instance()->TimeBin();
      Double_t ampl = TMath::Exp (-t*t/2/sigma2);
      if (TMath::Abs(t) > 5*sigma) ampl = 0;
      reFilt[i] = ampl;
    }
    fft[0] = TVirtualFFT::FFT(1, &nbins, "R2C ES K");
    //fft[0] = TVirtualFFT::FFT(1, &nbins, "R2C EX K");
    fft[0]->SetPoints(reFilt);
    fft[0]->Transform();
    fft[0]->GetPointsComplex(reFilt,imFilt);
  }

  Double_t *reSig = new Double_t [nbins];
  Double_t *imSig = new Double_t [nbins];
  Double_t *reTot = new Double_t [nbins];
  Double_t *imTot = new Double_t [nbins];

  //AZ Int_t nRows = MpdTpcSectorGeo::Instance()->NofRows();
  for (UInt_t iRow = 0; iRow < nRows; ++iRow) {
    for (UInt_t iPad = 0; iPad < fNumOfPadsInRow[iRow] * 2; ++iPad) {
      Int_t fired = 0;
      for (UInt_t iTime = 0; iTime < nbins; ++iTime) {
	//if (fDigits4dArray[iRow][iPad][iTime].signal > 0) {
	if (CalcOrigin(fDigits4dArray[iRow][iPad][iTime]) >= 0) {
	  // Fired channel
	  fired = 1;
	  reSig[iTime] = fDigits4dArray[iRow][iPad][iTime].signal;
	} else reSig[iTime] = 0;
      }
      if (!fired) continue;
      // Fourier transform
      fft[0]->SetPoints(reSig);
      fft[0]->Transform();
      fft[0]->GetPointsComplex(reSig,imSig);
      // Convolution
      //for (Int_t i = 0; i < nbins; ++i) {
      for (Int_t i = 0; i < n2; ++i) {
	Double_t re = reSig[i] * reFilt[i] - imSig[i] * imFilt[i];
	Double_t im = reSig[i] * imFilt[i] + imSig[i] * reFilt[i];
	reTot[i] = re / nbins;
	imTot[i] = im / nbins;
      }
      // Inverse Fourier transform
      if (!fft[1]) fft[1] = TVirtualFFT::FFT(1, &nbins, "C2R ES K");
      //if (!fft[1]) fft[1] = TVirtualFFT::FFT(1, &nbins, "C2R EX K");
      fft[1]->SetPointsComplex(reTot,imTot);
      fft[1]->Transform();
      fft[1]->GetPoints(reTot);

      for (Int_t i = 0; i < nbins; ++i) {
	if (fDigits4dArray[iRow][iPad][i].origin < 0) continue; // !!! do not add extra time bins due to shaping !!!
	Int_t i1 = i;
	if (i1 <= icent) i1 += icent;
	else i1 -= (icent + 1);
	Double_t ampl = reTot[i1];
	//
	// Scale factor to adjust ADC counts
	ampl /= 30.0;
	if (ampl > 4095.1) ampl = 4095.1; // dynamic range 12 bits
	//
	fDigits4dArray[iRow][iPad][i].signal = ampl;
      }
    }
  }
  delete [] reSig;
  delete [] imSig;
  delete [] reTot;
  delete [] imTot;
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::GetArea(Float_t xEll, Float_t yEll, Float_t radius, vector<UInt_t> &padIDs, vector<UInt_t> &rowIDs) 
{

  Float_t padW = 0.0, padH = 0.0;
  Float_t y = 0.0, x = 0.0;
  UInt_t pad = 0, row = 0;
  Float_t delta, yNext;
  if (fResponse) delta = 0.0;
  else delta = -1000.0; //for test only!!!
  
  fSecGeo->PadID(xEll, yEll, row, pad, yNext);

  for (Int_t ip = -1; ip < 2; ++ip) {
    Int_t pad1 = pad + ip;
    if (pad1 < 0) continue;
    if (pad1 >= MpdTpcSectorGeo::Instance()->NPadsInRows()[row] * 2) break;
    padIDs.push_back(pad1);
    rowIDs.push_back(row);
  }    
  // Add extra row
  if (TMath::Abs(yNext) < radius) {
    Int_t row1 = row;
    if (yNext < 0) {
      --row1;
      if (row1 < 0) return;
      if (fSecGeo->NPadsInRows()[row1] != fSecGeo->NPadsInRows()[row]) --pad; // different number of pads in rows
    } else if (yNext > 0) {
      ++row1;
      if (row1 > nRows - 1) return;
      if (fSecGeo->NPadsInRows()[row1] != fSecGeo->NPadsInRows()[row]) ++pad; // different number of pads in rows
    }
    for (Int_t ip = -1; ip < 2; ++ip) {
      Int_t pad1 = pad + ip;
      if (pad1 < 0) continue;
      if (pad1 >= MpdTpcSectorGeo::Instance()->NPadsInRows()[row1] * 2) break;
      padIDs.push_back(pad1);
      rowIDs.push_back(row1);
    }    
  }
}

//---------------------------------------------------------------------------

Float_t MpdTpcDigitizerAZlt::CalculatePadResponse(Int_t iloop, Int_t nLoop, UInt_t padID, UInt_t rowID, Float_t x, Float_t y) 
{
  // Calculate pad response using lookup table (assumes the same pad width in both readout plane regions)

  const Int_t npads = 5, nposX = 100, nposX2 = nposX * 2; // 5 pads, 100 positions along padrow (step 25 um)
  const Int_t nposY = 100, nposY2 = nposY * 2;
  static Int_t first = 1, padID0, rowID0, istep0, idist0;
  static Double_t chargeX[npads][nposX], stepX = pwIn / nposX2;
  static Double_t chargeY[2][nposY]; // charges on padrow (inner and outer regions)
  static Double_t cy;
  Double_t padW = pwIn, padH = 0, distx, disty, maxX, minX, cx1, cx2, stepY, minY, maxY;

  Double_t sigma = fSpread;

  if (first) {
    // Compute lookup tables
    first = 0;

    // Along X
    for (Int_t i = 0; i < nposX; ++i) {
      distx = i * stepX; // distance to pad center
      maxX = padW / 2 - distx;

      for (Int_t j = 0; j < npads; ++j) {
	if (j == 0) maxX -= 2 * padW;
	minX = maxX - padW;
	if (j == 0) cx2 = TMath::Erf(minX/sigma);
	else cx2 = cx1;
	cx1 = TMath::Erf(maxX/sigma);
	chargeX[j][i] = TMath::Abs(cx1-cx2) / 2.;
	maxX += padW;
      }
    }

    // Along Y
    // 2 pad heights
    for (Int_t i = 0; i < 2; ++i) {
      padH = (i == 0) ? phIn : phOut;
      stepY = padH / nposY / 2;

      for (Int_t j = 0; j < nposY; ++j) {
	disty = j * stepY; // distance to pad center
	maxY = padH / 2 - disty;
	minY = maxY - padH;
	chargeY[i][j] = TMath::Abs (TMath::Erf(maxY/sigma) - TMath::Erf(minY/sigma)) / 2.;
      }
    }

  }
  //exit(0);

  Int_t izone = 0, idist;

  // Different padrow
  if (iloop == 0 || rowID != rowID0) {
    Double_t padY;
    rowID0 = rowID;
    if (rowID < nInRows) {
      //padW = pwIn;
      padH = phIn;
      padY = padH * ((Double_t)rowID + 0.5); // y-coordinate of pad center
    } else {
      //padW = pwOut;
      padH = phOut;
      padY = fSectInHeight + (((Double_t)rowID - nInRows) + 0.5) * padH; // y-coordinate of pad center        
      izone = 1;
    }
    disty = y - padY;
    Int_t istep = TMath::Abs (disty / (padH / nposY2));
    if (istep < nposY) cy = chargeY[izone][istep];
    else {
      istep = nposY2 - istep;
      cy = 1.0 - chargeY[izone][istep];
    }

    Double_t padX = padW * ((Double_t)padID - fNumOfPadsInRow[rowID] + 0.5); // x-coordinate of pad center
    padID0 = padID;
    distx = x - padX;
    idist = TMath::Nint (distx / stepX);
    idist0 = idist;
    istep0 = TMath::Abs (idist);
    istep0 %= nposX2;
    if (istep0 > nposX) istep0 = TMath::Abs(istep0-nposX2);
    istep0 = TMath::Min (istep0, nposX-1);
  } else idist = idist0 - (padID-padID0) * nposX2; 

  Int_t ipad = TMath::Abs (idist / nposX); 
  if (ipad == 0) ipad = npads / 2; // central pad
  //else if (idist > 0) {
  else if (1) {
    if (ipad % 2 == 0) ipad = npads / 2 - ipad / 2;
    else ipad = npads / 2 + (ipad+1) / 2;
  } else {
    if (ipad % 2 != 0) ipad = npads / 2 + (ipad+1) / 2;
    else ipad = npads / 2 - ipad / 2;
  }
  cx1 = chargeX[ipad][istep0];

  Double_t ctot = TMath::Abs(cy) * cx1;
  //if (rowID == 1) cout << " Row, pad: " << rowID << " " << padID << " " << x << " " << istep0 << " " << idist << " " << ipad << " " << cx1 << " " << ctot << endl; 
  return ctot;
}

//---------------------------------------------------------------------------
TF1* MpdTpcDigitizerAZlt::padResponseFunction() {
    if (fPRF)
        return fPRF;

    fPRF = new TF1("Gaus PRF", "gaus", -5, 5);
    fPRF->SetParameter(0, 1.0 / (sqrt(2.0 * TMath::Pi()) * fSpread));
    fPRF->SetParameter(1, 0);
    fPRF->SetParameter(2, fSpread);

    return fPRF;
}

//---------------------------------------------------------------------------
Bool_t MpdTpcDigitizerAZlt::isSubtrackInInwards(const TpcPoint *p1, const TpcPoint *p2) { //WHAT AM I DOING???
    const Float_t x1 = p1->GetX();
    const Float_t x2 = p2->GetX();
    const Float_t y1 = p1->GetY();
    const Float_t y2 = p2->GetY();
    const Float_t a = (y1 - y2) / (x1 - x2);
    const Float_t b = (y1 * x2 - x1 * y2) / (x2 - x1);
    const Float_t minR = fabs(b) / sqrt(a * a + 1);

    if (minR < r_min) //then check if minimal distance is between our points
    {
        const Float_t x = -a * b / (a * a + 1);
        const Float_t y = b / (a * a + 1);
        if ((x1 - x) * (x2 - x) < 0 && (y1 - y) * (y2 - y) < 0) {
            return kTRUE;
        }
    }
    return kFALSE;
}

//---------------------------------------------------------------------------

void MpdTpcDigitizerAZlt::TpcProcessing(const TpcPoint* prePoint, const TpcPoint* curPoint, 
				      const UInt_t secID, const UInt_t iPoint, const UInt_t nPoints) 
{
    
  Float_t dE = 0.0; //energy loss
  UInt_t qTotal = 0; //sum of clusters charges (=sum of electrons between two TpcPoints)
  UInt_t qCluster = 0; //charge of cluster (= number of electrons)
  TLorentzVector curPointPos; // coordinates for current TpcPoint
  TLorentzVector prePointPos; // coordinates for previous TpcPoint
  TLorentzVector diffPointPos; // steps for clusters creation
  TVector3 diffuse; // vector of diffuse for every coordinates
  TVector3 distort; // vector of distortion for every coordinates
  TLorentzVector electronPos; // coordinates for created electrons
  TLorentzVector clustPos; // coordinates for created clusters
  Float_t driftl = 0.0; // length for drifting
  vector<UInt_t> clustArr; // vector of clusters between two TpcPoints
  Float_t localX = 0.0, localY = 0.0; //local coordinates of electron (sector coordinates)
  MpdTpcSectorGeo* secGeo = MpdTpcSectorGeo::Instance();
    
  if ( fPrintDebugInfo && (iPoint % 1000 == 0) ) cout << UInt_t(iPoint * 1.0 / nPoints * 100.0) << " % of TPC points processed" << endl;
  //        curPoint = (TpcPoint*) fMCPointArray->At(i);

  if (fOnlyPrimary == kTRUE) {
    MpdMCTrack* tr = (MpdMCTrack*) fMCTracksArray->At(curPoint->GetTrackID());
    if (tr->GetMotherId() != -1) return;
  }
  //check if hits are on the same track
  if (curPoint->GetTrackID() == prePoint->GetTrackID() && !isSubtrackInInwards(prePoint, curPoint)) {

    dE = curPoint->GetEnergyLoss() * 1E9; //convert from GeV to eV
    if (dE < 0) {
      Error("MpdTpcDigitizerTask::Exec", "Negative Energy loss!");
      return;
    }

    qTotal = (UInt_t) floor(fabs(dE / fGas->W()));
    if (qTotal == 0) return;

    curPointPos.SetXYZT(curPoint->GetX(), curPoint->GetY(), curPoint->GetZ(), curPoint->GetTime());
    prePointPos.SetXYZT(prePoint->GetX(), prePoint->GetY(), prePoint->GetZ(), prePoint->GetTime());
    if ((curPointPos.T() < 0) || (prePointPos.T() < 0)) {
      Error("MpdTpcDigitizerTask::Exec", "Negative Time!");
      return;
    }

    diffPointPos = curPointPos - prePointPos; //differences between two points by coordinates
    //AZ diffPointPos *= (1 / diffPointPos.Vect().Mag()); //directional cosines //TODO! Do we need this??? Look at line #297

    //while still charge not used-up distribute charge into next cluster

    if (fDistribute) {
      while (qTotal > 0) {
	//roll dice for next cluster
	qCluster = fGas->GetRandomCSUniform();
	if (qCluster > qTotal) qCluster = qTotal;
	qTotal -= qCluster;
	clustArr.push_back(qCluster);
      }// finish loop for cluster creation
    } else {
      clustArr.push_back(qTotal); // JUST FOR TEST. NO CLUSTER DISTRIBUTION!
      //             clustArr.push_back(1); // JUST FOR TEST. NO CLUSTER DISTRIBUTION ONLY ONE ELECTRON IN CLUSTER!
    }

    //AZ diffPointPos *= (diffPointPos.Vect().Mag() / clustArr.size()); //now here are steps between clusters by coordinates TODO: correct distribution
    diffPointPos *= (1./ clustArr.size());
    clustPos = prePointPos;

    for (UInt_t iClust = 0; iClust < clustArr.size(); ++iClust) {
      clustPos += diffPointPos;
      driftl = zCathode - fabs(clustPos.Z());
      if (driftl <= 0) continue;

      for (UInt_t iEll = 0; iEll < clustArr.at(iClust); ++iEll) {
	
	//attachment
	if (fAttach)
	  if (exp(-driftl * fGas->k()) < gRandom->Uniform()) continue; // FIXME

	//diffusion
	if (fDiffuse) {
	  const Float_t sqrtDrift = sqrt(driftl);
	  const Float_t sigmat = fGas->Dt() * sqrtDrift;
	  const Float_t sigmal = fGas->Dl() * sqrtDrift;
	  diffuse.SetXYZ(gRandom->Gaus(0, sigmat), gRandom->Gaus(0, sigmat), gRandom->Gaus(0, sigmal));
	}

	if (fDistort) {                   
	  
	  const Float_t dt = 1E-03;                            //time step [s]
	  const Float_t mu = 4.23;                             //electron mobility [m^2 / s / V]
	  const Float_t mu2 = mu * mu;                         //just square of mu

	  const TVector3 E(0.0, 0.0, fGas->E() * 100);         // vector of electric field components (now is constant and parallel to Z axes) // 100 - convert Ez from V/cm to V/m
	  TVector3 B(0.0, 0.0, 0.0);                           // vector of magnetic field components
	  TVector3 v;                                          // vector of current velocity components
	  TVector3 posCur;                                     // vector of current position components
	  TVector3 posPre = clustPos.Vect();                   // vector of previous position components
	  TVector3 EBCross(0.0, 0.0, 0.0);                     // vector product of E and B vectors
                        
	  Bool_t inTpc = kTRUE;
	  while (inTpc) {
                            
	    B.SetXYZ(fMagField->GetBx(posPre.X(), posPre.Y(), posPre.Z()) * 0.1, fMagField->GetBy(posPre.X(), posPre.Y(), posPre.Z()) * 0.1, fMagField->GetBz(posPre.X(), posPre.Y(), posPre.Z()) * 0.1);
	    EBCross = E.Cross(B);
                            
	    v = mu / (1 + mu2 * B.Mag2()) * (E - mu * EBCross + mu2 * B * (E * B));
	    posCur = v * dt + posPre;
	    //                            cout << "X = " << posCur.X() << " Y = " << posCur.Y() << " Z = " << posCur.Z() << " Vx = " << v.X() << " Vy = " << v.Y() << " Vz = " << v.Z() << endl;
	    if ((posCur.Perp() > r_min + fSectHeight) || (posCur.Perp() < r_min) || (Abs(posCur.Z()) > zCathode)) inTpc = kFALSE; 
	    posPre = posCur;
	  }
	  distort.SetX(posCur.X() - clustPos.Vect().X());
	  distort.SetY(posCur.Y() - clustPos.Vect().Y());
	  distort.SetZ(0.0);  //FIXME
	}

	electronPos.SetVect(clustPos.Vect() + diffuse + distort);
	electronPos.SetT(clustPos.T() + (zCathode - fabs(electronPos.Z())) / fGas->VDrift()); // Do we need to use clustPos.T() ???
                    
	/*AZ
	const Float_t phiStep = TwoPi() / nSectors * 2;
	const Float_t sectPhi = secID * phiStep;
	localY =  electronPos.X() * Cos(sectPhi) + electronPos.Y() * Sin(sectPhi) - r_min; //converting from global to local coordinates
	localX = -electronPos.X() * Sin(sectPhi) + electronPos.Y() * Cos(sectPhi);        //converting from global to local coordinates
	if ((localY < 0.0) || (Abs(localX) > fSector->GetMaxX()) || (localY > fSectHeight)) continue; //FIXME!!!
	const Float_t timeStep = (zCathode / fNTimeBins) / fGas->VDrift();
	const UInt_t curTimeID = (UInt_t) ((zCathode / fGas->VDrift() - electronPos.T()) / timeStep);
	*/

	//AZ
	TVector3 xyzLoc;
	if (secGeo->Global2Local(electronPos.Vect(), xyzLoc, secID % secGeo->NofSectors()) < 0) continue;
	localX = xyzLoc[0];
	localY = xyzLoc[1];
	if (!TMath::Finite(localX) || !TMath::Finite(localY)) {
	  // Debug
	  cout << " !!! Not finite " << secID % secGeo->NofSectors() << endl;
	  electronPos.Vect().Print();
	  diffuse.Print();
	}
	UInt_t curTimeID = UInt_t (secGeo->T2TimeBin(electronPos.T()));
	//cout << localX << " " << xyzLoc[0] << " " << localY << " " << xyzLoc[1] << " " << curTimeID << " " << timeBin << endl;
	//AZ
	if (curTimeID >= fNTimeBins) continue;
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
      } // for (UInt_t iEll = 0; iEll < clustArr.at(iClust);
    } // for (UInt_t iClust = 0; iClust < clustArr.size();
  } // if (curPoint->GetTrackID() == prePoint->GetTrackID()

  clustArr.clear();
  clustArr.resize(0);
}

//---------------------------------------------------------------------------
void MpdTpcDigitizerAZlt::Finish() {
    
    cout << "Digitizer work time = " << ((Float_t)tAll) / CLOCKS_PER_SEC << endl;
    
    if (fMakeQA) {
        toDirectory("QA/TPC");
        Float_t digit = 0.0;
        UInt_t iPad_shifted = 0; //needed for correct drawing of fDigitsArray

            for (UInt_t iRows = 0; iRows < nRows; ++iRows) {  
                for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[iRows] * 2; ++iPads) {
                    iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[iRows];
                    for (UInt_t iTime = 0; iTime < fNTimeBins; ++iTime) {
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
                for (UInt_t iTime = 0; iTime < fNTimeBins; ++iTime) {
                    digit = fDigits4dArray[iRows][iPads][iTime].signal;
                    //pad activity
                    //if (digit > 1000.0) {
                    //    fHisto->_hXY_dig->Fill(iPad_shifted, iRows, 1.0);
                    //}
                    //                    fHisto->_hXY_dig->Fill(iPad_shifted, iRows, digit);
                    fHisto->_h3D_dig->Fill(iPad_shifted, iRows, iTime, digit);
                }
            }
        }

        for (UInt_t iTime = 0; iTime < fNTimeBins; ++iTime) {
            for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[1] * 2; ++iPads) {
                iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[1];
                digit = fDigits4dArray[1][iPads][iTime].signal;
                fHisto->_hXT_dig_1->Fill(iPad_shifted, iTime, digit);
            }
            for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[5] * 2; ++iPads) {
                iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[5];
                digit = fDigits4dArray[5][iPads][iTime].signal;
                fHisto->_hXT_dig_5->Fill(iPad_shifted, iTime, digit);
            }
            for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[10] * 2; ++iPads) {
                iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[10];
                digit = fDigits4dArray[10][iPads][iTime].signal;
                fHisto->_hXT_dig_10->Fill(iPad_shifted, iTime, digit);
            }
            for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[20] * 2; ++iPads) {
                iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[20];
                digit = fDigits4dArray[20][iPads][iTime].signal;
                fHisto->_hXT_dig_20->Fill(iPad_shifted, iTime, digit);
            }
            for (UInt_t iPads = 0; iPads < fNumOfPadsInRow[40] * 2; ++iPads) {
                iPad_shifted = iPads + fNumOfPadsInRow[nRows - 1] - fNumOfPadsInRow[40];
                digit = fDigits4dArray[40][iPads][iTime].signal;
                fHisto->_hXT_dig_40->Fill(iPad_shifted, iTime, digit);
            }
        }

        fHisto->Write();
        gFile->cd();
    }
}

//---------------------------------------------------------------------------

ClassImp(MpdTpcDigitizerAZlt)
