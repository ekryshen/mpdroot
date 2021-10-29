#ifndef NICATPCSECTORGEO_H
#define NICATPCSECTORGEO_H

/// \ingroup rec
/// \class NicaTpcSectorGeo
/// \brief Geometry configuration of MPD TPC sector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include "NicaHelix.h"
#include <TObject.h>
#include <TVector3.h>

class NicaTpcSectorGeo : public TObject {
public:
  enum Shifts { kSectorS = 0, kPadrowS = 5, kPadS = 13, kPadSignS = 30 };
  enum Masks { kSectorM = 31, kPadrowM = 255, kPadM = 255, kPadSignM = 1 };

public:
  static NicaTpcSectorGeo* Instance(TString geo = "full");  ///< get singleton instance
  Int_t Global2Local(const Double_t* xyzGlob,
                     Double_t* xyzLoc,
                     Int_t iSec = -1);  ///< transform global coordinates to
                                        ///< local (sector) - returns padID
  Int_t Global2Local(const TVector3& xyzGlob,
                     TVector3& xyzLoc,
                     Int_t iSec = -1);  ///< transform global coordinates to
                                        ///< local (sector) - returns padID
  void Local2Global(Int_t iSec,
                    const Double_t* xyzLoc,
                    Double_t* xyzGlob);  ///< transform local coordinates of
                                         ///< sector iSec to global
  void Local2Global(Int_t iSec,
                    const TVector3& xyzLoc,
                    TVector3& xyzGlob);                                   ///< transform local coordinates of
                                                                          ///< sector iSec to global
  Int_t PadRow(Int_t padID) { return (padID >> kPadrowS) & kPadrowM; }    ///< pad rows from padID
  Int_t Sector(Int_t padID) { return (padID >> kSectorS) & kSectorM; }    ///< sector No. from padID
  Int_t PadID(Int_t sec, Int_t row) { return sec |= (row << kPadrowS); }  ///< padID from sector and padrow numbers
  void PadID(Float_t xloc, Float_t yloc, UInt_t& row, UInt_t& pad,
             Float_t& yNext);              ///< get row, pad and distance to nearest row
  TVector2 LocalPadPosition(Int_t padID);  ///< get local pad position for padID
  Int_t NofSectors() { return fgkNsect; }  ///< get number of TPC r/out sectors
  // Int_t NofRows() { return fgkNrows; } ///< get number of pad rows
  Int_t NofRows() { return fNrows[0] + fNrows[1]; }                          ///< get number of pad rows
  Int_t NofRowsReg(Int_t ireg) { return fNrows[ireg]; }                      ///< get number of pad rows in ROC region ireg
  Double_t Dphi() { return fDphi; }                                          ///< get sector angle
  Double_t GetMinY() { return fYsec[0]; }                                    ///< get minimum sector Y
  Double_t GetMaxY() { return fYsec[2]; }                                    ///< get maximum sector Y
  Double_t GetRocY(Int_t i) { return fYsec[i]; }                             ///< get ROC region Y-coordinates
  Double_t PadHeight(Int_t ireg = 0) { return fPadH[ireg]; }                 ///< get pad height in ROC region ireg
  Double_t PadWidth(Int_t ireg = 0) { return fPadW[ireg]; }                  ///< get pad width in ROC region ireg
  Double_t SectorAngle(Int_t iSec) { return fPhi0 + (iSec + 0.5) * fDphi; }  ///< get azimuthal angle of sector axis
  Int_t GetZmin() { return fZmin; }                                          ///< sensitive volume Zmin
  Int_t GetZmax() { return fZmax; }                                          ///< sensitive volume Zmax
  const Int_t* NPadsInRows() const { return fNPadsInRows; }                  ///< numbers of pads in rows
  Double_t Pad2Xloc(Double_t pad,
                    Int_t row);  ///< Pad number-to-Xlocal conversion
  void CalculatePads(const NicaHelix& helix, Float_t* Paths, Short_t* Pads, Bool_t shift = kFALSE);
  void SetNofRows(Int_t nRows, Int_t ireg = 0) { fNrows[ireg] = nRows; }
  void SetPadHeight(Double_t height, Int_t ireg = 0) { fPadH[ireg] = height; }
  void SetMinY(Double_t rmin) { fYsec[0] = rmin; }

protected:
  // virtual void Finish();
  NicaTpcSectorGeo(TString geo = "full") : fSectAngle(TMath::Pi() / 12.0) { Init(geo); }  ///< Default ctor
  Bool_t GetPad(const NicaHelix& helix, Double_t R, Double_t& S, Int_t& Sect);
  Bool_t CheckSector(const NicaHelix& helix, Double_t R, Double_t& S, Int_t Sect);
  Bool_t AdditionalCheck(const NicaHelix& helix, Double_t S, Double_t R, Double_t sect) const;
  virtual ~NicaTpcSectorGeo() { ; }  ///< Destructor
  Int_t FindFirstPad(const NicaHelix& helix, Float_t* Paths, Int_t& sect);
  Int_t FindLastPad(const NicaHelix& helix, Float_t* Paths, Int_t first_pad, Int_t& sect);

private:
  void Init(TString geo);
  // Automatic deletion when application exits
  static void DestroyInstance() {
    if (fgTpcSec) delete fgTpcSec;
  }

  static NicaTpcSectorGeo* fgTpcSec;  //! pointer to Singleton instance
  static const Int_t fgkNsect = 12;   // number of TPC sectors
  // static const Int_t fgkNrows = 50; // number of padrows
  Int_t fNrows[2];      // number of padrows 2 ROC regions
  Double_t fPhi0;       // phi0 of the first sector
  Double_t fDphi;       // sector angle
  Double_t fYsec[3];    // coordinates of ROC regions (height direction): minimum,
                        // boundary, maximum
  Double_t fPadH[2];    // pad heights for 2 ROC regions
  Double_t fPadW[2];    // pad widths for 2 ROC regions
  Double_t fZmin;       // sensitive volume Zmin
  Double_t fZmax;       // sensitive volume Zmax
  Int_t* fNPadsInRows;  // numbers of pads in rows
                        /** for pad calcuation **/
  Int_t fLastSector;
  Double_t fSecCH;
  Double_t fSinPhase;
  Double_t fCosPhase;
  Double_t fCHCosDipAngle;
  Double_t fSinDipAngle;
  Double_t fHelixR;
  Double_t fPeriod;
  const Double_t fSectAngle;
  ClassDef(NicaTpcSectorGeo, 1);
};

//__________________________________________________________________________
inline Double_t NicaTpcSectorGeo::Pad2Xloc(Double_t pad, Int_t row) {
  // Pad number-to-Xlocal conversion

  Double_t padW = (row < NofRowsReg(0)) ? fPadW[0] : fPadW[1];
  return padW * (pad - fNPadsInRows[row] + 0.5);
}

//__________________________________________________________________________

#endif
