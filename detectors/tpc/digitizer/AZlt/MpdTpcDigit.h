//--------------------------------------------------------------------
//
// Description:
//      Class for storing TPC Digits
//
//
// Author List:
//      Sergey Merts
//      7.12.2019 - Alexander Zinchenko: change of design
//        12.2022 - Slavomir Hnatic: AbstractTpcDigit port
//--------------------------------------------------------------------

#ifndef MPDTPCDIGIT_HH
#define MPDTPCDIGIT_HH

#include "AbstractTpcDigit.h"

class MpdTpcDigit : public AbstractTpcDigit {

public:
   MpdTpcDigit();
   MpdTpcDigit(Int_t ori, Int_t pad, Int_t row, Int_t bin, Int_t sec, Float_t adc);
   virtual ~MpdTpcDigit();

   Int_t   GetOrigin() const { return fOrigin; }
   Int_t   GetPad() const { return fPad; }
   Int_t   GetRow() const { return fRow; }
   Int_t   GetTimeBin() const { return fTimeBin; }
   Float_t GetSignal() const { return fAdc; }
   Int_t   GetSector() const { return fSector; }

   std::map<int, float> GetTrackSignals() const
   {
      std::map<int, float> trackID{{fOrigin, -1.}};
      return trackID;
   }
   // port to interface: function alias GetAdc = GetSignal
   inline Float_t GetAdc() const { return GetSignal(); }

   void SetPad(Int_t pad) { fPad = pad; }
   void SetRow(Int_t row) { fRow = row; }
   void SetTimeBin(Int_t timeBin) { fTimeBin = timeBin; }
   void SetAdc(Float_t adc) { fAdc = adc; }
   void SetSector(Int_t sector) { fSector = sector; }

private:
   Int_t   fOrigin;
   Int_t   fPad;
   Int_t   fRow;
   Int_t   fTimeBin;
   Int_t   fSector;
   Float_t fAdc;

   ClassDef(MpdTpcDigit, 5);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
