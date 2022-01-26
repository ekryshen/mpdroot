#ifndef MPDEMCDETSIMPARAMS_H
#define MPDEMCDETSIMPARAMS_H

#include "TObject.h"

class MpdEmcSimParams : public TObject
{
 public:
  /// Default constructor.
  /// It must be kept public for root persistency purposes,
  /// but should never be called by the outside world
  MpdEmcSimParams() {}
  ///
  /// Destructor.
  ///
  ~MpdEmcSimParams() {}

  ///
  /// \return the pointer of the _existing_ unique instance of the geometry
  /// It should have been set before with GetInstance(name) method
  ///
  static MpdEmcSimParams* GetInstance()
  {
    if (!sParams) {
      sParams = new MpdEmcSimParams();
    }
    return sParams;
  }

  // Parameters to configure Digitization
  bool SimulateNoise() const { return fSimulateNoise; }
  bool CorrectZ() const { return fCorrectZ; }
  bool ApplyNonLinearity() const { return fApplyNonlinearity; }
  bool ApplyDigitization() const { return fApplyDigitization; }
  bool ApplyTimeResolution() const { return fApplyTimeResolution; }

  double ElectronicNoise() const { return fElectronicNoiseWidth; }
  void CellNonLineaityParams(double& a, double& b, double& c) const
  {
    a = fCelNonlinParamA;
    b = fCelNonlinParamC;
    b = fCelNonlinParamC;
  }
  double ADCwidth() const { return fADCWidth; }
  double ZSthreshold() const { return fZSthreshold; }
  void TimeResolutionParams(double& a, double& b)
  {
    a = fTimeResolitionParamA;
    b = fTimeResolitionParamB;
  }
  void NoiseTimeRange(double& a, double& b) const
  {
    a = fNoiseTimeMin;
    b = fNoiseTimeMax;
  }
  bool SmearLightCollection() const { return fSmearLightCollection; }
  double GetLYPerGeV() const { return fEdepToLightYield; }

  // Parameters to configure Clusterization
  bool   AllowMultiSectorClusters() const { return fMultiSectorClusters; }
  double DigitMinEnergy() const { return fDigitMinEnergy; }
  double ClusteringThreshold() const { return fClusteringThreshold; }
  double LogWeight() const { return fLogWeight; }
  double LocalMaximumCut() const { return fLocalMaximumCut; }
  double ClusteringTimeGate() const { return fClusteringTimeGate; }
  bool UnfoldClusters() const { return fUnfoldClusters; }
  double UnfogingEAccuracy() const { return fUnfogingEAccuracy; }
  double UnfogingXZAccuracy() const { return fUnfogingXZAccuracy; }
  double EcoreCut1() const { return fEcoreCut1; }
  double EcoreCut2() const { return fEcoreCut2; }
  double Chi2radiusCut() const { return fChi2radiusCut; }
  int NMaxIterations() const { return fNMaxIterations; }
  int NLMMax() const { return fNLMMax; }
  int NPrimMax() const { return fNPrimMax; }
  double ENonLinCorrection(int index)const{return fNonLinCorrection[index];}

  // Parameters for cluster position correction (See MpdEmcCluster::CorrectVertex)
  double ZcorrSinA(int indx) const { return fZcorrSinA[indx]; }
  double ZcorrSinW(int indx) const { return fZcorrSinW[indx]; }
  double ZcorrA(int indx) const { return fZcorrA[indx]; }
  double ZcorrB(int indx) const { return fZcorrB[indx]; }

 private:
  static MpdEmcSimParams* sParams; // Pointer to the unique instance of the singleton

 public: // One can set parameters in configure macros, but use for access parameters only dedicated methods
  bool fCorrectZ = false; // Apply theta/zed correction
  bool fSmearLightCollection = false; // Emulate smearing and collecting light using fEdepToLightYield photoelectrons per GeV
  bool fSimulateNoise = false; // Simulate electronic noise in HitCreation
  bool fApplyNonlinearity = false;  // Apply energy non-linearity in HitCreation
  bool fApplyDigitization = false;   // Apply digitization of energy in HitCreation
  bool fApplyTimeResolution = false; // Apply time resolution in HitCreation

  double fElectronicNoiseWidth = 0.005;  // Width of Gaussian electronic noise in GeV
  double fCelNonlinParamA = -0.02;       // Cell energy non-linearity parameterization
  double fCelNonlinParamB = 0.5;         // in the form
  double fCelNonlinParamC = 1.0;         // e=e*c(1-a*exp(-e/b))
  double fADCWidth = 0.005;              // Widht of one ADC count in GeV (used in energy digitization in HitCreation)
  double fZSthreshold = 0.0015;          // ZeroSuppression threshold (remove digits below) in GeV
  double fTimeResolitionParamA = 5.e-10; // Parameters used for time resolution simulation
  double fTimeResolitionParamB = 2.e-11; // in the form width = a + b/e (in seconds)
  double fNoiseTimeMin = -100.e-9;       // simulate noise signal
  double fNoiseTimeMax = 100.e-9;        // in this range (in seconds)
  double fEdepToLightYield = 200000.;    // Number of photoelectrons per GeV

  // Clusterization
  bool   fMultiSectorClusters = true; // allow clusters with digits in different sectors
  double fLogWeight = 3.1;              // cutoff used in log. weight calculation
  double fDigitMinEnergy = 0.0015;      // Minimal energy of digits to be used in cluster (GeV)
  double fClusteringThreshold = 0.010; // Minimal energy of digit to start clustering (GeV)
  double fLocalMaximumCut = 0.0015;     // minimal height of local maximum over neighbours
  double fClusteringTimeGate = 1e9;  // maximal time difference between digits to be accepted to clusters (in ns)
  bool fUnfoldClusters = true;         // to perform cluster unfolding
  double fUnfogingEAccuracy = 1.e-4;   // accuracy of energy calculation in unfoding prosedure (GeV)
  double fUnfogingXZAccuracy = 1.e-2;  // accuracy of position calculation in unfolding procedure (cm)
  double fEcoreCut1 = 0.01;           // threshold for Ecore calculation E_p1
  double fEcoreCut2 = 0.02;           // threshold for Ecore calculation E_p2
  double fChi2radiusCut = 0.0001;      // cut in dispersion Chi2
  int fNMaxIterations = 100;             // maximal number of iterations in unfolding procedure
  int fNLMMax = 30;                    // Maximal number of local maxima in unfolding
  int fNPrimMax = 5;                   // Maximal number of primaries in list (sorted with deposited energy)
  float fNonLinCorrection[2] = {0.0, 1.0} ; //Parameters for Nonlinearity correction: Ecorr=[0]+[1]*E
  float fZcorrSinA[2] = {-0.5588, -0.042};  // Parameters for cluster position correction (See MpdEmcCluster::CorrectVertex)
  float fZcorrSinW[2] = {52.7, 0.4};        // Parameters for cluster position correction
  float fZcorrA[2] = {0.07887, 0.0101};     // Parameters for cluster position correction
  float fZcorrB[2] = {0.00256, 0.00032};    // Parameters for cluster position correction

  ClassDef(MpdEmcSimParams, 1)
};

#endif
