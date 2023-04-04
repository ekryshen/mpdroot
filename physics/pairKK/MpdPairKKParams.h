#ifndef MPDPAIRKKPARAMS_H
#define MPDPAIRKKPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdPairKKParams : public TObject {

public:
   //
   // Event selection cuts
   float mZvtxCut = 140.; // (V) event selection cut (cm)

   // PID cuts
   float mPIDsigTPC = 2.0; // (V)
   float mPIDsigTOF = 2.0; // (V)

   int   mNofHitsCut = 10;   // (V) minimal number of hits
   float mEtaCut     = 1.0;  // (V) maximal pseudorapidity 
   float mPtminCut   = 0.05; // (V) minimal pt
   float mDCACut     = 2.0 ; // (V) maximum DCA

   float mYCut       = 0.5 ; // (V) pair rapidity

   void ReadFromFile(std::string fname = "ConvDef");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdPairKKParams, 1);
};
#endif // MPDPAIRKKPARAMS_H
