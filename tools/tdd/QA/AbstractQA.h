//--------------------------------------------------------------------
// Description:
//      Abstract base class for generating QA plots
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, March, 2023
//--------------------------------------------------------------------

#ifndef ABSTRACTQA_HH
#define ABSTRACTQA_HH

#include <TObject.h>
#include <TString.h>

enum EQAMode {
   OFF = 0, // QA turned off by default
   TPCCLUSTERHITFINDER,
   ALL
};

class AbstractQA : public TObject {

public:
   AbstractQA() {}
   virtual ~AbstractQA() {}

   static EQAMode qaEngineMode;

private:
   ClassDef(AbstractQA, 1);
};

#endif
