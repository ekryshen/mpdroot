//--------------------------------------------------------------------
//
// Description:
//      MPD EMC Digitizer - takes EmcPoints and makes digits
//
//
//--------------------------------------------------------------------

#ifndef MPDEMCDIGITIZERKI_H
#define MPDEMCDIGITIZERKI_H 1

#include "MpdEmcDigitKI.h"

#include "FairTask.h"

class MpdEmcGeoParams;
class TClonesArray;
class MpdEmcSimParams;
class MpdEmcGeoUtils;

class MpdEmcDigitizerKI : public FairTask
{
 public:
  /** Default constructor **/
  MpdEmcDigitizerKI();

  /** Destructor **/
  ~MpdEmcDigitizerKI();

  /** Virtual method Init **/
  virtual InitStatus Init();

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);
  void virtual Finish();

 private:
  double SimulateNoiseEnergy();                             // Simulation of noise of electronics
  double NonLinearity(const double e);                      // simulate non-lineraity
  double DigitizeEnergy(const double e);                    // Account final width of ADC
  double TimeResolution(const double time, const double e); // Apply final time resolution
  double SimulateNoiseTime();                               // calculate time in noise digit
  double SimulateLightCollection(const double lostenergy);  // Simulate Poissonian light production and collection

 private:
  /** Input array of MpdEmcPoints **/
  TClonesArray* fPointArray;   //!

  TClonesArray* fMcTrArray;    //!

  /** Output array of MpdEmcHit **/
  TClonesArray* fDigitsArray;  //!

  MpdEmcSimParams* fSimParams; //! Class with all simulation parameters

  MpdEmcGeoUtils* fGeom;       //! Geometry parameters and methods
  int fNDigits;                //! Number of digits created
  ClassDef(MpdEmcDigitizerKI, 1);
};

#endif
