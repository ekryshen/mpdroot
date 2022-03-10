// -------------------------------------------------------------------------
//                            MpdFieldMapData header file              -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldMapData (PNDROOT)             -----
// -------------------------------------------------------------------------

#ifndef MPDMAGFIELDMAPDATA_H
#define MPDMAGFIELDMAPDATA_H

#include "TNamed.h"

class TArrayF;

class MpdFieldMap;

class MpdFieldMapData : public TNamed {
  
public:
  
  /** Default constructor **/
  MpdFieldMapData();

  /** Standard constructor **/
  MpdFieldMapData(const char* name);

  /** Constructor from an existing MpdFieldMap **/
  MpdFieldMapData(const char* name, const MpdFieldMap& map);
  
  /** Destructor **/
  virtual ~MpdFieldMapData();
  
  /** Accessors to field parameters in local coordinate system **/
  Int_t    GetType()  const { return fType; }
  Double_t GetXmin()  const { return fXmin; } 
  Double_t GetYmin()  const { return fYmin; }
  Double_t GetZmin()  const { return fZmin; }
  Double_t GetXmax()  const { return fXmax; }  
  Double_t GetYmax()  const { return fYmax; }
  Double_t GetZmax()  const { return fZmax; }
  Int_t    GetNx()    const { return fNx; }
  Int_t    GetNy()    const { return fNy; }
  Int_t    GetNz()    const { return fNz; }
  
  /** Accessors to the field value arrays **/
  TArrayF* GetBx() const { return fBx; }
  TArrayF* GetBy() const { return fBy; }
  TArrayF* GetBz() const { return fBz; }
	
 private:
  
  /** Type of map. 1 = MpdFieldMap, 2 = Sym2, 3 = Sym3 **/
  Int_t fType;
  
  /** Field limits in local coordinate system **/
  Double_t fXmin, fXmax;
  Double_t fYmin, fYmax;
  Double_t fZmin, fZmax;
  
  /**Original units of the map */
  Double_t fUnit; 

  /** Number of grid points  **/
  Int_t fNx, fNy, fNz;

  /** Arrays with the field values in T **/
  TArrayF* fBx;
  TArrayF* fBy;
  TArrayF* fBz;

  ClassDef(MpdFieldMapData,1) 
    };

#endif
