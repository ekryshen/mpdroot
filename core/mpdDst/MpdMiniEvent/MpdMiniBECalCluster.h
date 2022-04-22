/**
 * \class MpdMiniBECalCluster
 * \brief Holds information about ECal cluster
 *
 * The class holds information about the cluster reconstructed
 * in the barrel Electromagnetic Calorimeter (ECal)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI), Pavel Batyuk (JINR)
 * email: nigmatkulov@gmail.com; ganigmatkulov@mephi.ru; pavel.batyuk@jinr.ru
 * \date June 1, 2020
 */

#ifndef MpdMiniBECalCluster_h
#define MpdMiniBECalCluster_h

// ROOT headers
#include "TObject.h"
#include "TMath.h"

// C++ headers
#include <vector>

//_________________
class MpdMiniBECalCluster : public TObject {

public:
   /// Default constructor
   MpdMiniBECalCluster();
   /// Copy constructor
   MpdMiniBECalCluster(const MpdMiniBECalCluster &hit);
   /// Destructor
   virtual ~MpdMiniBECalCluster();
   /// Print cluster information
   virtual void Print(const Char_t *option = "") const;

   //
   // Getters
   //

   /// Return vector with digit IDs
   std::vector<UShort_t> digitIds() const { return fDigitIds; }
   /// Return vector with digit energies
   std::vector<Float_t> digitEDeps() const { return fDigitEDeps; }
   /// Return i-th digit parameters (DigitId and energy deposition)
   void digitParams(UShort_t i, Int_t &digId, Float_t &eneDep);
   /// Return number of digits in the cluster
   UInt_t nDigits() const { return fDigitIds.size(); }

   /// Return vector with MC track IDs
   std::vector<UShort_t> mcTrackIds() const { return fMcTrackIds; }
   /// Return vector with MC track energy deposited in ECal cell
   std::vector<Float_t> mcTrackEDeps() const { return fMcTrackEDeps; }
   /// Return i-th MC track parameters (mcTrackId, energy deposition)
   void mcTrackParams(UShort_t i, Int_t &mcTrkId, Float_t &eneDep);
   /// Return number of tracks that deposited energy in the cluster
   UInt_t nMcTracks() const { return fMcTrackIds.size(); }

   /// Return energy deposited in the cluster (GeV)
   Float_t energy() const { return fEnergy; }
   /// Return cluster energy core (GeV)
   Float_t eCore() const { return fEcore; }
   /// Return cluster energy core (GeV)
   Float_t eCore1p() const { return fEcore1p; }
   /// Return cluster energy core (GeV)
   Float_t eCore2p() const { return fEcore2p; }

   /// Return cluster average time (ns)
   Float_t time() const { return fTime; }

   /// Return sector index
   UShort_t sector(Int_t digitId) const { return digitId / 768; }
   /// Return chamber index
   UShort_t chamber(Int_t digitId) const { return digitId / 19200; }
   /// Return crate index
   UShort_t crate(Int_t digitId) const { return digitId / 128; }

   /// Return X coordinate of the cluster
   Float_t x() const { return fX; }
   /// Return Y coordinate of the cluster
   Float_t y() const { return fY; }
   /// Return Z coordinate of the cluster
   Float_t z() const { return fZ; }
   /// Return azimuthal angle of the cluster
   Float_t phi() const { return TMath::ATan2(fY, fX); }
   /// Return rho of the cluster
   Float_t rho() const { return TMath::Sqrt(fX * fX + fY * fY); }
   /// Return radius vector length of the cluster
   Float_t rad() const { return TMath::Sqrt(fX * fX + fY * fY + fZ * fZ); }

   /// Return distance of closest approach (in azimuth)
   /// from cluster center to track
   Float_t dPhi() const { return fdPhi; }
   /// Return distance of closest approach (along z axis)
   /// from cluster center to track
   Float_t dZ() const { return fdZ; }

   /// Return index of the track that produced the cluster
   Int_t trackId() const { return (Int_t)fTrackId; }
   /// Return matching flag (true - matched, false - not matched)
   Bool_t becalMatchFlag() const { return (trackId() < 0) ? kFALSE : kTRUE; }

   // Return dispertion of the shower
   // Float_t dispertion() const           { return fDisp; }
   /// Return smaller dispertion axis
   Float_t lambda1() const { return fLambda1; }
   /// Return larger dispertion axis
   Float_t lambda2() const { return fLambda2; }
   /// Return chi2 of fit with EM shape
   Float_t chi2() const { return (Float_t)fChi2 / 100; }
   /// Return number of local maxima in parent cluster before unfolding
   Int_t nlm() const { return (Int_t)fNExLM; }

   //
   // Setters
   //

   /// Set digit IDs
   void setDigitIds(std::vector<Int_t> ids);
   /// Set digit energies
   void setDigitEDeps(std::vector<Float_t> eDeps);
   /// Add digit with the give ID and energy
   void addDigit(Int_t id, Float_t eDep);

   /// Set MC track IDs
   void setMcTrackIds(std::vector<Int_t> ids);
   /// Set MC track energy deposition in ECal cells
   void setMcTrackEDeps(std::vector<Float_t> eDeps);
   /// Add MC track ID and energy deposition in cell
   void addMcTrack(Int_t id, Float_t eDep);
   /// Replace MCTrack index with miniMcTrack index
   void changeMcTrackIdx(Int_t mcTrkIdx, Int_t mcMiniTrkIdx);

   /// Set deposited energy in the cluster (GeV)
   void setEnergy(Float_t energy) { fEnergy = energy; }
   /// Set core energy (GeV)
   void setECore(Float_t energy) { fEcore = energy; }
   /// Set core energy 1p (GeV)
   void setECore1p(Float_t energy) { fEcore1p = energy; }
   /// Set core energy 2p (GeV)
   void setECore2p(Float_t energy) { fEcore2p = energy; }
   /// Set cluster time
   void setTime(Float_t t) { fTime = t; }
   /// Set cluster coordinates
   void setXYZ(Float_t x, Float_t y, Float_t z)
   {
      fX = x;
      fY = y;
      fZ = z;
   }
   /// Set distance of closest approach in phi
   void setDPhi(Float_t dPhi) { fdPhi = dPhi; }
   /// Set distance of closest approach in z
   void setDz(Float_t dz) { fdZ = dz; }
   /// Set index of track that is closest to center (-1 if not match)
   void setTrackId(Int_t id);

   // Set dispertion of a shower
   // void setDispertion(Float_t disp)        { fDisp = disp; }

   /// Set smaller disperion axis
   void setLambda1(Float_t lam) { fLambda1 = lam; }
   /// Set larger dispertion axis
   void setLambda2(Float_t lam) { fLambda2 = lam; }
   /// Set dispersion axes
   void setLambdas(Float_t lam1, Float_t lam2)
   {
      fLambda1 = lam1;
      fLambda2 = lam2;
   }
   /// Set chi2 of a fit with EM shape
   void setChi2(Float_t chi2);
   /// Set number of local maxima in parent cluster before unfolding
   void setNLM(Int_t nlm);

protected:
   /// Digit ID (cellId)
   std::vector<UShort_t> fDigitIds;
   /// Deposited energy per digit
   std::vector<Float16_t> fDigitEDeps;

   /// Index of the MC track
   std::vector<UShort_t> fMcTrackIds;
   /// Energy deposited by MC track
   std::vector<Float16_t> fMcTrackEDeps;

   /// Energy deposited in the cluster (GeV)
   Float16_t fEnergy;
   /// Cluster energy core (from MpdEmcClusterKI)
   Float16_t fEcore;
   /// Cluster energy core (from MpdEmcClusterKI)
   Float16_t fEcore1p;
   /// Cluster energy core (from MpdEmcClusterKI)
   Float16_t fEcore2p;

   /// Cluster average time (ns)
   Float16_t fTime;

   /// X coordinate of the cluster in global system
   Float16_t fX;
   /// Y coordinate of the cluster in global system
   Float16_t fY;
   /// Z coordinate of the cluster in global system
   Float16_t fZ;

   /// Distance to closest track in phi
   Float16_t fdPhi;
   /// Distance to closest track in z
   Float16_t fdZ;
   /// Index of closest matched track if any (-1 means no match)
   Short_t fTrackId;

   // Dispersion of a shower
   // Float16_t fDisp;
   /// Smaller dispersion axis
   Float16_t fLambda1;
   /// Larger dispersion axis
   Float16_t fLambda2;
   /// Chi2 of a fit with EM shape (chi2 * 100)
   UShort_t fChi2;

   /// Number of local maxima in parent cluster before unfolding
   UChar_t fNExLM;

   ClassDef(MpdMiniBECalCluster, 1);
};

#endif // #define MpdMiniBECalCluster_h
