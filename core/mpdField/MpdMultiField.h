// -------------------------------------------------------------------------
//                            MpdMultiField header file                -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdMultiField    (PNDROOT)            -----
// -------------------------------------------------------------------------

#ifndef MPDMULTIMAP_H
#define MPDMULTIMAP_H 1

#include "FairField.h"
#include "TObjArray.h"
#include <map>

class MpdRegion;
class MpdMultiFieldPar;
class MpdMultiField : public FairField {
  
 public:
  
  /** Default constructor **/
  MpdMultiField();
  
  /** Constructor from MpdFieldPar **/
  MpdMultiField(MpdMultiFieldPar* fieldPar);
  
  /** Destructor **/
  virtual ~MpdMultiField();
  
  /** Initialisation (read map from file) **/
  void Init();
  
  /**Adding a field to the collection*/

  void AddField(FairField *field);
  TObjArray *GetFieldList(){return fMaps; }

  /** Get the field components at a certain point 
   ** @param x,y,z     Point coordinates (global) [cm]
   ** @value Bx,By,Bz  Field components [kG]
   **/
  
  virtual Double_t GetBx(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBy(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBz(Double_t x, Double_t y, Double_t z);
  
  void GetFieldValue(const Double_t point[3], Double_t* bField);
  
  /** Screen output **/
  virtual void Print();
  ClassDef(MpdMultiField,1) 

    protected:
  
  TObjArray *fMaps;
  Int_t     fNoOfMaps;
  std::map <MpdRegion*, FairField* > fFieldMaps;//!
};

#endif



