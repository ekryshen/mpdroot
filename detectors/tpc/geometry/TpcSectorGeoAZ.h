#ifndef TPCSECTORGEOAZ_H
#define TPCSECTORGEOAZ_H

/// \ingroup rec
/// \class MpdTpcSectorGeo
/// \brief Geometry configuration of MPD TPC sector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna
///         
///        singleton removal, BaseTpcGeo inheritance port:
///        Slavomir Hnatic, LIT JINR Dubna, September 2022

#include <BaseTpcGeo.h>
#include <TVector3.h>

class TpcGas;

class TpcSectorGeoAZ : public BaseTpcGeo {

public:
   enum Shifts { kSectorS = 0, kPadrowS = 5, kPadS = 13, kPadSignS = 30 };
   enum Masks { kSectorM = 31, kPadrowM = 255, kPadM = 255, kPadSignM = 1 };

public:
   TpcSectorGeoAZ();              ///< Default ctor
   virtual ~TpcSectorGeoAZ(){ ; } ///< Destructor   

   Int_t    Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc,
                         Int_t iSec = -1); ///< transform global coordinates to local (sector) - returns padID
   Int_t    Global2Local(const TVector3 &xyzGlob, TVector3 &xyzLoc,
                         Int_t iSec = -1); ///< transform global coordinates to local (sector) - returns padID
   void     Local2Global(Int_t iSec, const Double_t *xyzLoc,
                         Double_t *xyzGlob); ///< transform local coordinates of sector iSec to global
   void     Local2Global(Int_t iSec, const TVector3 &xyzLoc,
                         TVector3 &xyzGlob); ///< transform local coordinates of sector iSec to global
   Int_t    PadRow(Int_t padID) { return (padID >> kPadrowS) & kPadrowM; }   ///< pad rows from padID
   Int_t    Sector(Int_t padID) { return (padID >> kSectorS) & kSectorM; }   ///< sector No. from padID
   Int_t    PadID(Int_t sec, Int_t row) { return sec |= (row << kPadrowS); } ///< padID from sector and padrow numbers
   void     PadID(Float_t xloc, Float_t yloc, UInt_t &row, UInt_t &pad,
                  Float_t &yNext);            ///< get row, pad and distance to nearest row
   TVector2 LocalPadPosition(Int_t padID);    ///< get local pad position for padID
   Int_t    NofSectors() { return fgkNsect; } ///< get number of TPC r/out sectors
   // Int_t NofRows() { return fgkNrows; } ///< get number of pad rows
   Int_t        NofRows() { return fNrows[0] + fNrows[1]; }       ///< get number of pad rows
   Int_t        NofRowsReg(Int_t ireg) { return fNrows[ireg]; }   ///< get number of pad rows in ROC region ireg
   Double_t     Dphi() { return fDphi; }                          ///< get sector angle
   Double_t     GetMinY() { return fYsec[0]; }                    ///< get minimum sector Y
   Double_t     GetMaxY() { return fYsec[2]; }                    ///< get maximum sector Y
   Double_t     GetRocY(Int_t i) { return fYsec[i]; }             ///< get ROC region Y-coordinates
   Double_t     PadHeight(Int_t ireg = 0) { return fPadH[ireg]; } ///< get pad height in ROC region ireg
   Double_t     PadWidth(Int_t ireg = 0) { return fPadW[ireg]; }  ///< get pad width in ROC region ireg
   Double_t     SectorAngle(Int_t iSec) { return fPhi0 + (iSec + 0.5) * fDphi; } ///< get azimuthal angle of sector axis
   Int_t        GetNTimeBins() { return fNTimeBins; }                            ///< number of time bins
   Int_t        GetZmin() { return fZmin; }                                      ///< sensitive volume Zmin
   Int_t        GetZmax() { return fZmax; }                                      ///< sensitive volume Zmax
   const Int_t *NPadsInRows() const { return fNPadsInRows; }                     ///< numbers of pads in rows
   Double_t     Z2TimeBin(Double_t z);                                           ///< Z-to-Time bin conversion
   Double_t     T2TimeBin(Double_t time);                                        ///< time-to-Time bin conversion
   Double_t     TimeBin2Z(Double_t timeBin);                                     ///< Time bin-to-Z conversion
   Double_t     Pad2Xloc(Double_t pad, Int_t row);                               ///< Pad number-to-Xlocal conversion
   Double_t     TimeMax() const { return fTimeMax; }                             ///< max drift time
   Double_t     TimeBin() const { return fTimeBin; }                             ///< time bin length

   void SetNofRows(Int_t nRows, Int_t ireg = 0) { fNrows[ireg] = nRows; }
   void SetPadHeight(Double_t height, Int_t ireg = 0) { fPadH[ireg] = height; }
   void SetMinY(Double_t rmin) { fYsec[0] = rmin; }

protected:
   

private:
      
   static const Int_t      fgkNsect = 12; // number of TPC sectors
   // static const Int_t fgkNrows = 50; // number of padrows
   Int_t    fNrows[2];    // number of padrows 2 ROC regions
   Double_t fPhi0;        // phi0 of the first sector
   Double_t fDphi;        // sector angle
   Double_t fYsec[3];     // coordinates of ROC regions (height direction): minimum,  boundary, maximum
   Double_t fPadH[2];     // pad heights for 2 ROC regions
   Double_t fPadW[2];     // pad widths for 2 ROC regions
   Double_t fZmin;        // sensitive volume Zmin
   Double_t fZmax;        // sensitive volume Zmax
   Double_t fNTimeBins;   // number of time bins
   Int_t   *fNPadsInRows; // numbers of pads in rows
   Double_t fZ2TimeBin;   // Z-to-Time bin conversion coefficient
   TpcGas  *fGas;         // pointer to gas system
   Double_t fTimeMax;     // max drift time
   Double_t fTimeBin;     // time bin length
   // Int_t fTimeBinMax; // max time bin
   Double_t fTimeBinMax; // max time bin

   ClassDef(TpcSectorGeoAZ, 1);
};

//__________________________________________________________________________
inline Double_t TpcSectorGeoAZ::Z2TimeBin(Double_t z)
{
   // Z-to-Time bin conversion

   return (fZmax - TMath::Abs(z)) * fZ2TimeBin;
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAZ::T2TimeBin(Double_t time)
{
   // Time-to-Time bin conversion

   return time / fTimeBin;
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAZ::TimeBin2Z(Double_t timeBin)
{
   // Time bin-to-Z conversion

   // return (fTimeBinMax - timeBin - 0.5) / fZ2TimeBin;
   return (fTimeBinMax - timeBin - 0.5 + 0.037) / fZ2TimeBin; // extra correction
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAZ::Pad2Xloc(Double_t pad, Int_t row)
{
   // Pad number-to-Xlocal conversion

   Double_t padW = (row < NofRowsReg(0)) ? fPadW[0] : fPadW[1];
   return padW * (pad - fNPadsInRows[row] + 0.5);
}
//__________________________________________________________________________

#endif
