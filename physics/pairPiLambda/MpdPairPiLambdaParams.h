#ifndef MPDPAIRPILAMBDAPARAMS_H
#define MPDPAIRPILAMBDAPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdPairPiLambdaParams : public TObject {

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

   float mChi2PionLam = 7.0;    // (V) Chi2-to-PV for pion from Lambda
   float mChi2ProtLam = 3.0;    // (V) Chi2-to-PV for proton from Lambda
   float mLamEtaCut   = 1.2;    // maximum pseudorapidity for Lambda
   float mChi2Lam     = 3.0;    // maximum Chi2 for Lambda secondary vertex
   float mPALam       = 0.1;    // maximum pointing angle for Lambda
   float mDecayLam    = 0.5;    // minimum decay distance for Lambda
   float mDistLam     = 1.0;    // maximum distance between p-pi in the secondary vertex from Lambda
   float mNSigmaLam   = 3.0;    // n-sigma selection for Lambda candidates
   float mWidthLam    = 2.0e-3; // width of Lambda peak

   void ReadFromFile(std::string fname = "ConvDef");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdPairPiLambdaParams, 1);
};
#endif // MPDPAIRPILAMBDAPARAMS_H
