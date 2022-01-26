//-----------------------------------------------------------
//
// Author List:
//          D.Peresunko, KI, 2019
//-----------------------------------------------------------

#ifndef MPDEMCCLUSTERIZERKI_H
#define MPDEMCCLUSTERIZERKI_H

// Base Class Headers ----------------
#include <map>
#include "FairLogger.h"
#include "FairTask.h"

// Collaborating Class Declarations --
class TClonesArray;
class TObjArray;
class MpdEmcSimParams;
class MpdEmcGeoUtils;
class MpdEmcClusterKI;
class MpdEmcCalibParams;

class MpdEmcClusterizerKI : public FairTask
{
 public:
  // Constructors/Destructors ---------
  MpdEmcClusterizerKI();
  ~MpdEmcClusterizerKI() {}

  virtual InitStatus Init();
  void FinishTask() {}

  virtual void Exec(Option_t* opt);

  static double ShowerShape(double dx, double dz); // Parameterization of EM shower

 protected:
  void PrepareDigits();  // Calibrate, Allpy BadMap, clean...
  void MakeClusters();   // Do the job
  void MakeUnfoldings(); // Find and unfold clusters with few local maxima
  void UnfoldOneCluster(MpdEmcClusterKI* iniClu, Int_t nMax, int* digitId, float* maxAtEnergy);
  // Performs the unfolding of a cluster with nMax overlapping showers
  // Parameters: iniClu cluster to be unfolded
  //             nMax number of local maxima found (this is the number of new clusters)
  //             digitId: index of digits, corresponding to local maxima
  //             maxAtEnergy: energies of digits, corresponding to local maxima
  void EvalClusters();

 private:
  int fNumberOfClusters; // Number of created clusters

  TClonesArray* fDigitsArray;  //! Input digits array
  TObjArray* fClustersArray;   //! output clusters array
  MpdEmcSimParams* fSimParams; //! Configuration parameters
  MpdEmcGeoUtils* fGeom;       //! Geometry class
  MpdEmcCalibParams * fCalibData;  //! Calibration parameters
  ClassDef(MpdEmcClusterizerKI, 1)
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
