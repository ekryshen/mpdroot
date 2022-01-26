#ifndef MPDEMCCLUSTERKI_H
#define MPDEMCCLUSTERKI_H 1

#include <map>
#include "TMath.h"
#include "TObject.h"
class MpdEmcDigitKI;
class TLorentzVector;
class TVector3;

using namespace std;

class MpdEmcClusterKI : public TObject
{
 public:
  /** Default constructor **/
  MpdEmcClusterKI();

  /** Constructor with hit parameters **/
  MpdEmcClusterKI(const MpdEmcDigitKI* digit);

  /** Destructor **/
  virtual ~MpdEmcClusterKI();

  // Returm momentum of photon assuming it came from the provided vertex
  void GetMomentum(TLorentzVector& p, const TVector3* vertex) const;

  void AddDigit(const MpdEmcDigitKI* digit, Double_t edep = 0);
  void EvalAll(); // Evaluate cluster parameters

  void Print(const Option_t* opt = 0) const;

  void Purify(); // Remove digits below threshold

  void CorrectVertex(double vZ);

  int GetNumberOfLocalMax(int* maxAt, float* maxAtEnergy) const; // Finds local maxima

  int GetDigitTowerId(int i) const { return fDigitIDEnergy.at(i).first; } // detectorId of i-th digit

  // Get tower ID and energy of i-th cotributing digit
  void GetDigitParams(int i, int& detId, float& eDigit) const
  {
    if (i >= 0 && i < fNDigits) {
      detId = fDitisId[i];
      eDigit = fDigitsE[i];
    } else {
      detId = -1;
      eDigit = 0;
    }
  }
  // Method used in unfolding and should not be used by analysers
  void GetTransientDigitParams(int i, int& detId, float& eDigit) const
  {
    std::pair<int, float> p = fDigitIDEnergy.at(i);
    detId = p.first;
    eDigit = p.second;
  }

  // Number of MC tracks
  int GetNumberOfTracks() const { return fNPrimaries; }
  void GetMCTrack(Int_t i, Int_t& trackId, Float_t& trackEdep) const
  {
    if (i >= 0 && i < fNPrimaries) {
      trackId = fPrimId[i];
      trackEdep = fPrimE[i];
    } else {
      trackId = -1;
      trackEdep = 0;
    }
  }
  Float_t GetE() const { return fE; };

  Float_t GetEcore() const { return fEcore; };

  Float_t GetEcore_1p() const { return fEcore1p; };

  Float_t GetEcore_2p() const { return fEcore2p; };

  Float_t GetChi2() const { return fChi2; };

  Float_t GetTime() const { return fTime; };

  Float_t GetPhi() const { return TMath::ATan2(fY, fX); };

  Float_t GetRho() const { return TMath::Sqrt(fX * fX + fY * fY); };

  Float_t GetX() const { return fX; };

  Float_t GetY() const { return fY; };

  Float_t GetZ() const { return fZ; };

  Float_t GetDPhi() const { return fdPhi; };

  Float_t GetDZ() const { return fdZ; };

  Float_t GetRad() const { return TMath::Sqrt(fX * fX + fY * fY + fZ * fZ); };

  Int_t GetMultiplicity() const
  {
    return TMath::Max((int)fDigitIDEnergy.size(), fNDigits); // either from std::vector or in final cluster - from array
  };

  Float_t GetTrackIndex() const { return fTrackId; };

  Int_t GetNLM() const { return fNExLM; }

  void SetNLM(int n) { fNExLM = n; }

  /*   void SetFlag(Int_t flag)  {
         fFlag = flag;
     };
 */

  void SetTime(Float_t time) { fTime = time; };

  void SetTrackIndex(Int_t ind) { fTrackId = ind; }
  // Distance to track extrapolation
  void SetTrackDxDz(Float_t dphi, Float_t dz)
  {
    fdPhi = dphi;
    fdZ = dz;
  }

  // Dispersion parameters
  void GetLambdas(Float_t& l1, Float_t& l2)
  {
    l1 = fLambda1;
    l2 = fLambda2;
  }

 protected:
  void FillArrays(); // Fill arrays to store to disk

 protected:
  // Parameters used for re-calibration if necessary
  std::vector<std::pair<int, float>> fDigitIDEnergy; //! transient list of contributed digits with energies
  std::map<int, float> fMCTracks;                    //! transient trackID and energy deposit

  Int_t fNDigits;    // Digit multiplicity
  Int_t* fDitisId;   //[fNDigits] cellId
  Float_t* fDigitsE; //[fNDigits] deposited energy

  Int_t fNPrimaries; // Number of primaries
  Int_t* fPrimId;    //[fNPrimaries] cellId
  Float_t* fPrimE;   //[fNPrimaries] deposited energy per cell

  Float_t fE;       // cluster energy
  Float_t fEcore;   // cluster energy core
  Float_t fEcore1p; // cluster energy core
  Float_t fEcore2p; // cluster energy core

  Float_t fTime; // cluster time

  // cluster coordinates in global system
  Float_t fX; // x-coordinate of cluster
  Float_t fY; // y-coordinate of cluster
  Float_t fZ; // z-coordinate of cluster

  // Distance to closest track
  Float_t fdPhi;  // dPhi
  Float_t fdZ;    // dZed
  Int_t fTrackId; // Index of matched (closest) track if any

  // Dispersion and shower shape parameters
  Float_t fDisp;    // Dispersion
  Float_t fChi2;    // Chi2 of a fit with EM shape
  Float_t fLambda1; // smaller Disp axis
  Float_t fLambda2; // larger Disp axis

  Int_t fNExLM; // Number of local maxima or NLM in parent cluster before unfolding

  ClassDef(MpdEmcClusterKI, 1)
};

#endif
