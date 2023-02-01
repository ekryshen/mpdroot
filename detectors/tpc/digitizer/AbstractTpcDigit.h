//--------------------------------------------------------------------
// Description:
//      Test-driven interface (Abstract base class) for TPC Digits
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, December, 2022
//--------------------------------------------------------------------

#ifndef ABSTRACTTPCDIGIT_HH
#define ABSTRACTTPCDIGIT_HH

#include <TObject.h>
#include <map>

class AbstractTpcDigit : public TObject {

public:
   AbstractTpcDigit() {}
   virtual ~AbstractTpcDigit() {}

   virtual int   GetSector() const  = 0;
   virtual int   GetRow() const     = 0;
   virtual int   GetPad() const     = 0;
   virtual int   GetTimeBin() const = 0;
   virtual float GetSignal() const  = 0;

   // for debug purposes
   // amount of signal coming from different MC tracks <MC track ID, trackSignal>
   virtual std::map<int, float> GetTrackSignals() const = 0;

private:
   ClassDef(AbstractTpcDigit, 1);
};

#endif
