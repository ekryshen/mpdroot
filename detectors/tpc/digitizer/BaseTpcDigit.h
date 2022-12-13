//--------------------------------------------------------------------
// Description:
//      Base class for storing TPC Digits
//
//      Any custom TPC Digit class must inherit from this class
//
//      Note: - this class is a data structure by its nature
//            - has no getter/setter mechanism
//            - contains debug information in a map,
//              initialized only when calling debug constructor
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, December, 2022
//--------------------------------------------------------------------

#ifndef BASETPCDIGIT_HH
#define BASETPCDIGIT_HH

#include <TObject.h>
#include <map>

class BaseTpcDigit : public TObject {

public:
   BaseTpcDigit() : sector(-1), row(-1), pad(-1), timeBin(-1), signal(-1.0) {}
   BaseTpcDigit(int sector, int row, int pad, int timeBin, float signal)
      : sector(sector), row(row), pad(pad), timeBin(timeBin), signal(signal)
   {
   }

   // debug constructor
   BaseTpcDigit(int sector, int row, int pad, int timeBin, float signal, std::map<int, float> trackSignals)
      : sector(sector), row(row), pad(pad), timeBin(timeBin), signal(signal), trackSignals(trackSignals)
   {
   }

   virtual ~BaseTpcDigit() {}

   int   sector;
   int   pad;
   int   row;
   int   timeBin;
   float signal;

   // for debug purposes
   // amount of signal coming from different MC tracks <MC track ID, trackSignal>
   std::map<int, float> trackSignals;

private:
   ClassDef(BaseTpcDigit, 1);
};

#endif
