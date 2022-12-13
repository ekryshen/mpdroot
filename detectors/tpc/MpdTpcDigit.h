//--------------------------------------------------------------------
//
// Description:
//      Class for storing TPC Digits
//
//
// Author List:
//      Sergey Merts
//      7.12.2019 - Alexander Zinchenko: change of design
//        12.2022 - Slavomir Hnatic: BaseTpcDigit port
//--------------------------------------------------------------------

#ifndef MPDTPCDIGIT_HH
#define MPDTPCDIGIT_HH

#include "BaseTpcDigit.h"

class MpdTpcDigit : public BaseTpcDigit {

public:
   MpdTpcDigit();
   MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc);
   virtual ~MpdTpcDigit();

   Int_t   GetOrigin() const { return fOrigin; }
   Int_t   GetPad() const { return pad; }
   Int_t   GetRow() const { return row; }
   Int_t   GetTimeBin() const { return timeBin; }
   Float_t GetAdc() const { return signal; }
   Int_t   GetSector() const { return sector; }

   void SetPad(Int_t newPad) { pad = newPad; }
   void SetRow(Int_t newRow) { row = newRow; }
   void SetTimeBin(Int_t newTimeBin) { timeBin = newTimeBin; }
   void SetAdc(Float_t newAdc) { signal = newAdc; }
   void SetSector(Int_t newSector) { sector = newSector; }

private:
   Int_t fOrigin;

   ClassDef(MpdTpcDigit, 4);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
