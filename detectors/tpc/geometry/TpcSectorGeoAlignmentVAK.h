#ifndef TPCSECTORGEOALIGNMENTVAK_H
#define TPCSECTORGEOALIGNMENTVAK_H
/// \ingroup rec
/// \class TpcSectorGeoAlignmentVAK
/// \brief Geometry configuration of MPD TPC sector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna
/// \corrected by Valentin Kuzmin, SINP of Moscow State University

#include <BaseTpcSectorGeo.h>
#include <TVector3.h>
#include <TRotation.h>

class TpcSectorGeoAlignmentVAK : public BaseTpcSectorGeo
/*
  The TpcSectorGeoAlignmentVK class describes the geometry of the TPC sectors.
It provides transformation from the global MPD TPC coordinate system to
sectors local coordinate systems, defines the PadID structure and the
numbering sectors, pad rows and pads in rows. The transformation take
in account the TPC sectors allignment.
  Sector numbering for the camera is counterclockwise (looking at the end
face outside the camera). For a camera with Z>0, the numbering starts at 0
from the X axis (TPC coordinate system). For a camera with Z<0, the
numbering starts at 12 from the negative x axis (TPC coordinate system).
  The coordinate axes of the local coordinate system of the sector form the
right three. The Y axis is directed from the center of the chamber cylinder
along it radius. The Z-axis of the sector for the positive part of the
camera coincides with the direction of the z-axis of the global system,
in the negative half they are opposite. From the point of view outside the
camera on its end face, the direction of the local x axis is such that its
rotation towards the local Y axis occurs counterclockwise for both cameras.
 The local coordinate center is on the camera cilindr axis (xGlob=yGlob=0)
*/
{

public:
   enum Shifts { kSectorS = 0, kPadrowS = 5, kPadS = 13, kPadSignS = 30 };
   enum Masks { kSectorM = 31, kPadrowM = 255, kPadM = 255, kPadSignM = 1 };

   TpcSectorGeoAlignmentVAK()
   { ///< Default constructor
      for (Int_t i = 0; i < 24; i++) {
         aR0_lsc[i] = TVector3(0., 0., 0.);
         alpha_E[i] = 0.;
         beta_E[i]  = 0.;
         gamma_E[i] = 0.;
      }
      Init();
      PutSectorAlignment();
   };

   virtual ~TpcSectorGeoAlignmentVAK() { ; } ///< Destructor

   // *************** obsolate **********************
   Int_t Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc,
                      Int_t iSec); ///< transform global coordinates
                                   /// to local (sector) - returns padID
   Int_t Global2Local(const TVector3 &xyzGlob, TVector3 &xyzLoc,
                      Int_t iSec); ///< transform global coordinates to
                                   ///  local (sector) - returns padID

   Int_t Global2Local(TVector3 xyzGlob, TVector3 &xyzLoc);
   void  PadID(Double_t xloc, Double_t yloc, UInt_t &row, UInt_t &pad, Float_t &yNext);

   ///< Get the allignment from the file
   void TpcSectorGeoA(TString *inFile)
   {
      GetSectorAlignment(inFile);
      PutSectorAlignment();
   }
   /// Get the sector allignment as a shift and Euler's angles
   void GetAlignment(Int_t iSec, TVector3 &r0_loc, Double_t &alpha, Double_t &beta, Double_t &gamma);
   ///< to use a special allignment
   void TpcSectorGeoA(Double_t r0_lsc[3][24], Double_t alpha[24], Double_t beta[24], Double_t gamma[24]);

   /// LCS->GCS geolocation of the sector (the shift & the transformation)
   void SectorTransformation(const Int_t iSec, TVector3 &SecCoor_shift_toglo, TRotation &fromLocSecCoor);
   /// GCS->LCS geolocation of the sector (the shift & the transformation)
   void SectorBackTransformation(const Int_t iSec, TVector3 &SecCoor_shift, TRotation &fromGloSecCoor);
   /// TLCS->GCS geolocation of the sector (the shift & the transformation)
   void TLCS2GCStransformation(const Int_t iSec, TVector3 &TLCS2GCSTshift, TRotation &TLCS2GCST);
   /// GCS->TLCS geolocation of the sector (the shift & the transformation)
   void GCS2TLCStransformation(const Int_t iSec, TVector3 &GCS2TLCSshift, TRotation &GCS2TLCS);
   /// TLCS->LCS geolocation of the sector (the shift & the transformation)
   void TLCS2LCStransformation(const Int_t iSec, TVector3 &TLCS2LCSshift, TRotation &TLCS2LCS);
   /// LCS->TLCS geolocation of the sector (the shift & the transformation)
   void LCS2TLCStransformation(const Int_t iSec, TVector3 &LCS2TLCSshift, TRotation &LCS2TLCS);
   ///< transform local designed coordinates of sector iSec to global
   void DesignSectorTransformation(const Int_t iSec, TRotation &fromDesSecCoor) { fromDesSecCoor = T_tl2glo[iSec]; }

   /// Inishalization of current allignment parameters of one sector
   void PutSectorAlignment(Int_t iSec, TVector3 R0, Double_t alpha, Double_t beta, Double_t gamma);
   /// Inishalization of current allignment parameters by given values
   void PutSectorAlignment();
   /// Sectors Transformations due to Alignment(the shift & the transformation)
   void PrintAlignmentTransformation();
   ///< transform global coordinates to local (sector) & returns padID
   Int_t Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc, Double_t &dist);
   ///< transform global coordinates to local (sector) - returns padID
   Int_t Global2Local(const TVector3 xyzGlob, TVector3 &xyzLoc, Double_t &dist);
   ///< transform local coordinates of sector iSec to global
   void Local2Global(Int_t iSec, const Double_t *xyzLoc, Double_t *xyzGlob);
   ///< transform local coordinates of sector iSec to global
   void Local2Global(Int_t iSec, const TVector3 xyzLoc, TVector3 &xyzGlob);

   ///< pad rows from padID
   Int_t PadRow(Int_t padID) { return (padID >> kPadrowS) & kPadrowM; }
   ///< pad number in the row from padID
   Int_t Pad(Int_t padID) { return (padID >> kPadS) & kPadM; }
   ///< sector No. from padID
   Int_t Sector(Int_t padID) { return (padID >> kSectorS) & kSectorM; }
   ///< padID from sector and padrow numbers
   Int_t PadID(Int_t sec, Int_t row) { return sec | (row << kPadrowS); }
   ///< padID from sector, padrow and pad numbers
   Int_t PadID(Int_t sec, Int_t row, Int_t pad) { return sec | ((row << kPadrowS) | (pad << kPadS)); }
   ///< pads plane z coordinate
   TVector2 zPads() { return TVector2(fPadZ[0], fPadZ[1]); }

   ///< get row, pad and distance to nearest row
   void PadID(Double_t xloc, Double_t yloc, UInt_t &row, UInt_t &pad, Double_t &yNext);
   ///< check is the local point inside the sensitive sector area  (ROC)
   Bool_t IsItInside(Double_t xloc, Double_t yloc);
   ///< check is the local point inside a pad & get its local pad position
   Bool_t IsItPad(Double_t xloc, Double_t yloc, Int_t &irow, Int_t &ipad, Double_t &x1, Double_t &x2, Double_t &y1,
                  Double_t &y2);
   ///< get local pad position for padID
   TVector2 LocalPadPosition(Int_t padID);
   ///< get local pad position of the pad
   TVector2 LocalPadPosition(Int_t irow, Int_t ipad);
   ///< get number of ROC
   Int_t NofRoc(Int_t irow)
   {
      if (irow < fNrows[0])
         return 0;
      else
         return 1;
   }
   ///< get number of row
   Int_t NofRow(Double_t yloc);
   ///< get number of pad
   Int_t NofPad(Int_t row, Double_t xloc);
   ///< get number of TPC r/out sectors
   Int_t NofSectors() { return fgkNsect; }
   ///< get number of pad rows
   Int_t NofRows() { return fFNRows; }
   ///< get number of pad rows in ROC region ireg
   Int_t NofRowsReg(Int_t ireg) { return fNrows[ireg]; }
   ///< get number of pads in row
   Int_t NofPadsInRow(Int_t irow) { return fNPadsInRow[irow]; }
   ///< get sector angle
   Double_t Dphi() { return fDphi; }
   ///< get minimum sector Y
   Double_t GetMinY() { return fYsec[0]; }
   ///< get maximum sector Y
   Double_t GetMaxY() { return fYsec[2]; }
   ///< get local minimum sector Y
   Double_t GetMinYl() { return flYsec[0]; }
   ///< get local maximum sector Y
   Double_t GetMaxYl() { return flYsec[2]; }
   ///< get the half the lowest pads row
   Double_t GetMinX() { return fXsec[0]; }
   ///< get the half the highest pads row
   Double_t GetMaxX() { return fXsec[1]; }
   ///< get ROC region Y-coordinates
   Double_t GetRocY(Int_t i) { return fYsec[i]; }
   ///< get min&max pad coordinates
   void MinMaxPadPosition(Int_t irow, Int_t ipad, Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax);
   ///< get pad height of the row
   Double_t PadHeight(Int_t irow = 0)
   {
      if (irow < fNrows[0])
         return fPadH[0];
      else
         return fPadH[1];
   }
   ///< get low pad row edge of the row
   Double_t LowRowEdge(Int_t irow = 0)
   {
      if (irow < fNrows[0])
         return fPadH[0] * irow;
      else
         return flYsec[1] + fPadH[1] * (irow - fNrows[0]);
   }
   ///< get uper pad row edge of the row
   Double_t UpRowEdge(Int_t irow = 0)
   {
      if (irow < fNrows[0])
         return fPadH[0] * (irow + 1);
      else
         return flYsec[1] + fPadH[1] * (irow - fNrows[0] + 1);
   }
   ///< get pad squarE_trans_lsce for the ROC region
   Double_t GetPadSforROC(Int_t iroc) { return fPadS[iroc]; }
   ///< get pad square
   Double_t GetPadS(Int_t irow = 0)
   {
      if (irow < fNrows[0])
         return fPadS[0];
      else
         return fPadS[1];
   }
   ///< get pad width in ROC region ireg
   Double_t PadWidth(Int_t ireg = 0) { return fPadW[ireg]; }
   ///< get azimuthal angle of sector boundary
   Double_t SectorAngle(Int_t iSec)
   {
      Double_t phi = fPhi0 + (iSec - iSec / 12 * 12) * fDphi;
      if ((phi - TMath::TwoPi()) > 0) phi -= TMath::TwoPi();
      return phi;
   }
   ///< get azimuthal angle of sector boundary
   Double_t SectorAxisAngle(Int_t iSec) { return (iSec - iSec / 12 * 12) * fDphi; }
   ///< numbers of pads in rows below
   const Int_t *NPadsInRows() const { return fNPadsInRows; }
   ///< numbers of pads in a sector
   const Int_t NPadsInSec() const { return fNPads; }

   ///< number of time bins
   Int_t GetNTimeBins() { return fNTimeBins; }
   ///< sensitive volume Zmin
   Int_t GetZmin() { return fZmin; }
   ///< sensitive volume Zmax
   Int_t GetZmax() { return fZmax; }
   ///< Z-to-Time bin conversion
   Double_t T2TimeBin(Double_t time);            ///< time-to-Time bin conversionм
   Double_t TimeBin2Z(Double_t timeBin);         ///< Time bin-to-Z conversion
   Double_t Pad2Xloc(Double_t pad, Int_t row);   ///< Pad number-to-Xlocal conversion
   Double_t TimeMax() const { return fTimeMax; } ///< max drift time
   Double_t TimeBin() const { return fTimeBin; } ///< time bin length

   void SetNofRows(Int_t nRows, Int_t ireg = 0) { fNrows[ireg] = nRows; }
   void SetPadHeight(Double_t height, Int_t ireg = 0) { fPadH[ireg] = height; }
   void SetMinY(Double_t rmin) { fYsec[0] = rmin; }
   void PrintSectorParam();
   void PrintSectorAlignment(Int_t isec); // prits the alignment of the sector isec
                                          // if isec<0 prints for all sectors

private:
   // Initialization of current allignment parameters from MPD DB
   void GetSectorAlignment(TString *inFile);
   // Inishalization of the sector parameters
   void Init();

   static const Int_t fgkNsect = 12; // number of TPC sectors in one chamber

   TVector3 R0_lsc[24];  // the shift of the local coordinate system (LCS)
   TVector3 G0_gsc[24];  // the shift of the global coordinate system
                         // inside the local coordinate system
   TVector3 R0_tlsc[24]; // the shift of the Theoretical LCS
   TVector3 G0_tlsc[24]; // the shift of the global coordinate system
                         // inside the Theoretical LCS
   TVector3 aR0_lsc[24]; // the shift of the LCS in TLCS by alligment
   TVector3 aG0_lsc[24]; // the shift of the inside LCS

   Int_t fFNRows; // Full number of rows

   // GCS - detector global coordinate system
   // TLCS - theoretical sector coordinate system
   // LCS - sector local coordinate system due to alignment
   Double_t alpha_E[24]; // Euler angle alpha(phi) - X rotation
   Double_t beta_E[24];  // Euler angle alpha(theta) - Y' rotation
   Double_t gamma_E[24]; // Euler angle gamma(psi) - Z'' rotation

   TRotation T_tl2glo[24]; // Sector rotation TLCS -> GCS
   TRotation T_glo2tl[24]; // Sector rotation TLCS -> GCS
   TRotation A_l2tl[24];   // Euler matrix LCS -> TLCS
   TRotation A_tl2l[24];   // Back transformation TLCS -> LCS

   TRotation R_loc2glo[24]; // sector rotation to global TPC coordinates
   TRotation R_glo2loc[24]; // Inverse global sector rotation to TPC coordinates

   Double_t fPadZ[2];
   ;                               // Z coordinate of the pads plane
   Int_t    fNrows[2];             // number of padrows 2 ROC regions
   Double_t fPhi0;                 // phi0 of the first sector
   Double_t fDphi;                 // sector angle
   Double_t fYsec[3];              // "global" coordinates of ROC regions (heights):
                                   // minimum, boundary, maximum
   Double_t flYsec[3];             // local coordinates ROC regions (heights)
   Double_t fXsec[2];              // sector bottom half width, sector top half width
   Double_t fPadH[2];              // pad heights for 2 ROC regions
   Double_t fPadW[2];              // pad widths for 2 ROC regions
   Double_t fPadS[2];              // pad widths for 2 ROC regions
   Int_t   *fNPadsInRow;           // numbers of pads in the row
   Int_t   *fNPadsInRows;          // numbers of pads in rows below
   Int_t    N_InnerPads_min;       // number of pads in the first inner row
   Int_t    fNPads;                // number of pads in the sector
   Double_t Z2TimeBin(Double_t z); ///< Z-to-Time bin conversion
   Double_t fZmin;                 // sensitive volume Zmin
   Double_t fZmax; // sensitive volume ZmaxНайдётся всё. человек который изменил все
   Double_t fNTimeBins;  // number of time bins
   Double_t fZ2TimeBin;  // Z-to-Time bin conversion coefficient
   Double_t fTimeMax;    // max drift time
   Double_t fTimeBin;    // time bin length
   Double_t fTimeBinMax; // max time bin

   ClassDef(TpcSectorGeoAlignmentVAK, 1);
};

//__________________________________________________________________________
inline Double_t TpcSectorGeoAlignmentVAK::Z2TimeBin(Double_t z)
{
   // Z-to-Time bin conversion

   return (fZmax - TMath::Abs(z)) * fZ2TimeBin;
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAlignmentVAK::T2TimeBin(Double_t time)
{
   // Time-to-Time bin conversion

   return time / fTimeBin;
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAlignmentVAK::TimeBin2Z(Double_t timeBin)
{
   // Time bin-to-Z conversion

   // return (fTimeBinMax - timeBin - 0.5) / fZ2TimeBin;
   return (fTimeBinMax - timeBin - 0.5 + 0.037) / fZ2TimeBin; // extra correction
}

//__________________________________________________________________________
inline Double_t TpcSectorGeoAlignmentVAK::Pad2Xloc(Double_t pad, Int_t row)
{
   // Pad number-to-Xlocal conversion

   Double_t padW = (row < NofRowsReg(0)) ? fPadW[0] : fPadW[1];
   return padW * (pad - fNPadsInRows[row] + 0.5);
}
//__________________________________________________________________________

#endif
