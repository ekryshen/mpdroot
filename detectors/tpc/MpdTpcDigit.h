//--------------------------------------------------------------------
//
// Description:
//      Class for storing TPC Digits
//
//
// Author List:
//      Sergey Merts
//      7.12.2019 - Alexander Zinchenko: change of design
//
//--------------------------------------------------------------------

#ifndef MPDTPCDIGIT_HH
#define MPDTPCDIGIT_HH

#include <TObject.h>
#include <map>

class MpdTpcDigit : public TObject {

public:
   MpdTpcDigit();
   MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc);
   virtual ~MpdTpcDigit();

   Int_t   GetOrigin() const { return fOrigin; }
   Int_t   GetPad() const { return fPad; }
   Int_t   GetRow() const { return fRow; }
   Int_t   GetTimeBin() const { return fTimeBin; }
   Float_t GetAdc() const { return fAdc; }
   Int_t   GetSector() const { return fSector; }

   void SetPad(Int_t pad) { fPad = pad; }
   void SetRow(Int_t row) { fRow = row; }
   void SetTimeBin(Int_t bin) { fTimeBin = bin; }
   void SetAdc(Float_t adc) { fAdc = adc; }
   void SetSector(Int_t sec) { fSector = sec; }

private:
   Int_t   fPad;
   Int_t   fRow;
   Int_t   fTimeBin;
   Int_t   fSector;
   Float_t fAdc;
   Int_t   fOrigin;

   ClassDef(MpdTpcDigit, 3);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
