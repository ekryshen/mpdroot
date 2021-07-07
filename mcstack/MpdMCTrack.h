// -------------------------------------------------------------------------
// -----                      FairMCTrack header file                  -----
// -----                  Created 03/08/04  by V. Friese               -----
// -----              adopted for NICA/MPD 29/03/10  (litvin)          -----
// -----              adopted for NICA/MPD 20/12/19  (ABychkov)        -----
// -------------------------------------------------------------------------


/** FairMCTrack.h
 *@author V.Friese <v.friese@gsi.de>
 **
 ** Data class for storing Monte Carlo tracks processed by the MpdStack (FairStack).
 ** A MCTrack can be a primary track put into the simulation or a
 ** secondary one produced by the transport through decay or interaction.
 **
 ** Redesign 13/06/07 by V. Friese
 **/


#ifndef MPDMCTRACK_H
#define MPDMCTRACK_H 1


#include "MpdDetectorList.h"

#include "TLorentzVector.h"
#include "TObject.h"
#include "TParticle.h"
#include "TVector3.h"

#ifndef ROOT_TParticlePDG
 #include "TParticlePDG.h"
#endif
#ifndef ROOT_TDatabasePDG
 #include "TDatabasePDG.h"
#endif


class MpdMCTrack : public TObject
{

 public:


  /**  Default constructor  **/
  MpdMCTrack();


  /**  Standard constructor  **/
  MpdMCTrack(Int_t pdgCode, Int_t motherID, Double_t px, Double_t py,
	     Double_t pz, Double_t x, Double_t y, Double_t z, 
	     Double_t t, Int_t nPoints);

  /**  Copy constructor  **/
  MpdMCTrack(const MpdMCTrack& track);


  /**  Constructor from TParticle  **/
  MpdMCTrack(TParticle* particle);


  /**  Destructor  **/
  virtual ~MpdMCTrack();


  /**  Output to screen  **/
  void Print(Int_t iTrack=0) const;


  /**  Accessors  **/
  Int_t    GetPdgCode()  const { return fPdgCode; }
  Int_t    GetMotherId() const { return fMotherId; }
  Double_t GetPx()       const { return fPx; }
  Double_t GetPy()       const { return fPy; }
  Double_t GetPz()       const { return fPz; }
  Double_t GetStartX()   const { return fStartX; }
  Double_t GetStartY()   const { return fStartY; }
  Double_t GetStartZ()   const { return fStartZ; }
  Double_t GetStartT()   const { return fStartT; }
  Double_t GetMass()     const;
  Double_t GetEnergy()   const;
  Double_t GetPt()       const { return TMath::Sqrt(fPx*fPx+fPy*fPy); }
  Double_t GetP() const { return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz); }
  Double_t GetRapidity() const;
  void GetMomentum(TVector3& momentum); 
  void Get4Momentum(TLorentzVector& momentum);
  void GetStartVertex(TVector3& vertex);
  Double_t GetPolar(Int_t indx) { return fPolar[indx]; }
  Double_t GetWeight() { return fWeight; }

  /** Accessors to the number of MCPoints in the detectors **/
  Int_t GetNPoints(DetectorIdMPD detId)  const;


  /**  Modifiers  **/
  void SetMotherId(Int_t id) { fMotherId = id; }
  void SetNPoints(Int_t iDet, Int_t np);
  void SetPolar(Int_t indx, Double_t polar) { fPolar[indx] = polar; }
  void SetWeight(Double_t weight) { fWeight = weight; }

private:

  /**  PDG particle code  **/
  Int_t  fPdgCode;

  /**  Index of mother track. -1 for primary particles.  **/
  Int_t  fMotherId;

  /** Momentum components at start vertex [GeV]  **/
  Double32_t fPx, fPy, fPz;

  /** Coordinates of start vertex [cm, ns]  **/
  Double32_t fStartX, fStartY, fStartZ, fStartT;

  /**  Bitvector representing the number of MCPoints for this track in 
   **  each subdetector. The detectors are represented by
   **  kSTS:         Bit  0      (1 bit,  max. value  1)
   **  kTPC:         Bit  1 
   **  kTOF:         Bit  2 
   **  kETOF:        Bit  3 
   **  kFD:          Bit  4 
   **  kECT:         Bit  5 
   **  kECAL:        Bit  6
   **  kNDET:        Bit  7 
   **  kCPC:         Bit  8 
   **  kBBC:         Bit  9
   **  kZDC:         Bit  10 
   **  kFSA:         Bit  11
   **  The respective point numbers can be accessed and modified 
   **  with the inline functions. 
   **  Bits 12-31 are spare for potential additional detectors.
   **/
  Int_t fNPoints;
  /**AZ Polarization  **/
  Double32_t fPolar[3]; 

  /** Accounts for the length of the polarization vector **/
  Double32_t fWeight;

  ClassDef(MpdMCTrack,5);

};



// ==========   Inline functions   ========================================

inline Double_t MpdMCTrack::GetEnergy() const {
  Double_t mass = GetMass();
  return TMath::Sqrt(mass*mass + fPx*fPx + fPy*fPy + fPz*fPz );
}


inline void MpdMCTrack::GetMomentum(TVector3& momentum) { 
  momentum.SetXYZ(fPx,fPy,fPz); 
}


inline void MpdMCTrack::Get4Momentum(TLorentzVector& momentum) {
  momentum.SetXYZT(fPx,fPy,fPz,GetEnergy()); 
}


inline void MpdMCTrack::GetStartVertex(TVector3& vertex) { 
  vertex.SetXYZ(fStartX,fStartY,fStartZ); 
}





#endif
