#ifndef MPDPAIRPIKSPARAMS_H
#define MPDPAIRPIKSPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdPairPiKsParams : public TObject {

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
   float mDCACut     = 2.0;  // (V) maximum DCA

   float mYCut = 0.5; // (V) pair rapidity

   float mChi2PionKs = 7.0;    // (V) Chi2-to-PV for pion from Ks
   float mKsEtaCut   = 1.2;    // maximum pseudorapidity for Ks
   float mChi2Ks     = 3.0;    // maximum Chi2 for Ks secondary vertex
   float mPAKs       = 0.1;    // maximum pointing angle for Ks
   float mDecayKs    = 0.5;    // minimum decay distance for Ks
   float mDistKs     = 1.0;    // maximum distance between p-pi in the secondary vertex from Ks
   float mNSigmaKs   = 3.0;    // n-sigma selection for Ks candidates
   float mWidthKs    = 5.0e-3; // width of Ks peak

   void ReadFromFile(std::string fname = "ConvDef");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdPairPiKsParams, 1);
};
#endif // MPDPAIRPIKSPARAMS_H
