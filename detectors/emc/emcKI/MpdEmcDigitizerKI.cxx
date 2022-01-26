//--------------------------------------------------------------------
//
// Description:
//      MPD EMC Digitizer - takes EmcPoints and makes digits
//
//
// Author List:
//                D.Peresunko, KI, 2019
//--------------------------------------------------------------------

#include "MpdEmcDigitizerKI.h"
#include "MpdMCTrack.h"
#include "MpdEmcDigitKI.h"
#include "MpdEmcGeoUtils.h"
#include "MpdEmcPointKI.h"
#include "MpdEmcSimParams.h"

#include "FairEventHeader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRun.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>

using namespace std;
using namespace TMath;

// -----   Default constructor   -------------------------------------------

MpdEmcDigitizerKI::MpdEmcDigitizerKI()
  : FairTask("EMC digitizer"),
    fPointArray(nullptr),
    fDigitsArray(nullptr),
    fSimParams(nullptr),
    fGeom(nullptr),
    fNDigits(0)
{
}

// -----   Destructor   ----------------------------------------------------

MpdEmcDigitizerKI::~MpdEmcDigitizerKI() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcDigitizerKI::Init()
{
  LOG(INFO) << "******************* EMC INIT *********************" << endl;

  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if (!ioman) {
    LOG(FATAL) << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  // Get input array
  fPointArray = (TClonesArray*)ioman->GetObject("EmcPoint");
  if (!fPointArray) {
    LOG(ERROR) << "No EmcPoint array!" << endl;

    return kERROR;
  }
  fMcTrArray = (TClonesArray*)ioman->GetObject("MCTrack");

  // Create and register output array
  fDigitsArray = new TClonesArray(
    "MpdEmcDigitKI", 100); // TODO: Check, can digits with different length of primary list be in TClonesArray???
  ioman->Register("EmcDigit", "EMC", fDigitsArray, kTRUE);

  LOG(INFO) << "Intialization successfull" << endl;

  return kSUCCESS;
}

//__________________________________________________________________________

void MpdEmcDigitizerKI::Finish() { LOG(INFO) << "Finish" << endl; }

//__________________________________________________________________________

void MpdEmcDigitizerKI::Exec(Option_t* opt)
{
  // Add all EMCPoints with energy deposition in same tower
  // Simulate electronic noise
  // Add non-linearity, digitization of energy
  // Remove digits in bad map and below threshold

  LOG(INFO) << " Event No. " << FairRun::Instance()->GetEventHeader()->GetMCEntryNumber();

  // Reset output Array
  if (!fDigitsArray)
    Fatal("MpdEmcDigitizerKI::Exec", "No array of digits");

  // Reset output Array

  if (!fDigitsArray) {
    LOG(FATAL) << "No array of digits";
  }
  fDigitsArray->Clear();
  fNDigits = 0;

  LOG(INFO) << "ECAL: number of points in the event: " << fPointArray->GetEntriesFast();

  // Class with list of parameters
  if (!fSimParams) {
    fSimParams = MpdEmcSimParams::GetInstance();
  }
  if (!fGeom) {
    fGeom = MpdEmcGeoUtils::GetInstance();
  }

  // Find the first cell with signal
  int nPoints = fPointArray->GetEntriesFast();
  int pointIndex = 0;
  int nextSigId = 999999; // in case no ginal in
  MpdEmcPointKI* point = nullptr;
  if (nPoints > 0) {
    point = static_cast<MpdEmcPointKI*>(fPointArray->At(pointIndex));
    nextSigId = point->GetDetectorID();
  }

  int nTotCells = fGeom->GetTotalNCells();
  // go through all cells, either add noisy digit, or signal+noise
  for (Int_t cellId = 0; cellId < nTotCells; cellId++) {
    // If signal exist in this cell, add noise to it, otherwise just create noise digit
    if (cellId == nextSigId) {
      // Create new Digit
      MpdEmcDigitKI* digit = new ((*fDigitsArray)[fNDigits++]) MpdEmcDigitKI(point);
      pointIndex++;
      while (pointIndex < nPoints) {
        point = static_cast<MpdEmcPointKI*>(fPointArray->At(pointIndex));

        if (digit->CanAdd(point)) { // Point from the same Tower
          digit->AddPoint(point);
          pointIndex++;
        } else { // Points are sorted according to cellID. If no more points left, finish
          nextSigId = point->GetDetectorID();
          break;
        }
      }
      // Add Electroinc noise, apply non-linearity, digitize, de-calibrate, time resolution
      double energy = digit->GetE();
      // Emulate Poissonian light collection
      if (fSimParams->SmearLightCollection()) {
        energy = SimulateLightCollection(energy);
      }

      if (fSimParams->SimulateNoise()) {
        // Simulate electronic noise
        energy += SimulateNoiseEnergy();
      }

      if (fSimParams->ApplyNonLinearity()) {
        energy = NonLinearity(energy);
      }
      if (fSimParams->ApplyDigitization()) {
        energy = DigitizeEnergy(energy);
      }

      //Convert energy from GeV to ADC units - as in real data
      energy/=fSimParams->ADCwidth();//V
      digit->SetE(energy);

      if (fSimParams->ApplyTimeResolution()) {
        digit->SetTime(TimeResolution(digit->GetTime(), energy));
      }
    } else { // No signal in this cell,
      // Simulate noise
      if (fSimParams->SimulateNoise()) {
        double energy = SimulateNoiseEnergy();
        double time = SimulateNoiseTime();
        if (energy > fSimParams->ZSthreshold()) {
          energy/=fSimParams->ADCwidth();
          new ((*fDigitsArray)[fNDigits++])
            MpdEmcDigitKI(cellId, energy, time, -1); // current sellId, energy, random time, no primary
        }
      }
    }
  }

  LOG(INFO) << "Emc points done:" << fDigitsArray->GetEntriesFast();
  LOG(INFO) << system("date");
}

//_______________________________________________________________________
double MpdEmcDigitizerKI::SimulateNoiseEnergy()
{
  // Simulation of noise of electronics
  // Here we assume, that noise is independent on signal
  // and just Gaus with fixed width
  return gRandom->Gaus(0., fSimParams->ElectronicNoise());
}
//_______________________________________________________________________
double MpdEmcDigitizerKI::NonLinearity(const double e)
{
  double a, b, c;
  fSimParams->CellNonLineaityParams(a, b, c);
  return e * c * (1. + a * exp(-e / b));
}
//_______________________________________________________________________
double MpdEmcDigitizerKI::DigitizeEnergy(const double e)
{
  // distretize energy if necessary
  double w = fSimParams->ADCwidth();
  return w * ceil(e / w);
}
//_______________________________________________________________________
double MpdEmcDigitizerKI::TimeResolution(const double time, const double e)
{
  // apply time resolution
  if (e <= 0)
    return 0.;
  double a, b;
  fSimParams->TimeResolutionParams(a, b);
  double timeResolution = a + b / e;
  return gRandom->Gaus(time, timeResolution);
}
//_______________________________________________________________________
double MpdEmcDigitizerKI::SimulateNoiseTime()
{
  // Evaluate (random) time of noise digits as uniform in some ranges
  double a, b;
  fSimParams->NoiseTimeRange(a, b);
  return gRandom->Uniform(a, b);
}
//_______________________________________________________________________
double MpdEmcDigitizerKI::SimulateLightCollection(const double lostenergy)
{
  // Emulate production of scintillator light and its collection
  double coef = fSimParams->GetLYPerGeV();
  double lightYield = gRandom->Poisson(coef * lostenergy);
  return lightYield / coef;
}

ClassImp(MpdEmcDigitizerKI)
