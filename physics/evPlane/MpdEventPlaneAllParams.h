#ifndef MPDEVENTPLANEALLPARAMS_H
#define MPDEVENTPLANEALLPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdEventPlaneAllParams : public TObject {

public:
   //
   // Event selection cuts
   float mZvtxCut = 130.; //(V) event selection cut (cm)

   // Track selection cuts for EP from TPC
   int   mNofHitsCut = 16;  //(V) minimal number of hits to accept track
   float mEtaCut     = 1.5; //(V) maximal pseudorapidity accepted
   float mEtaGapCut =
      0.1; //(V) pseudorapidity gap between 2 TPC sub events (deltaEtaGap). Default 0.1 -> from -0.05 to 0.05.
   float mPtminCut = 0.1; //(V) minimal pt used in analysis
   float mPtmaxCut = 2.0; //(V) maximal pt used in analysis
   float mDcaCut   = 2.0; //(V) maximal DCA accepted

   std::string mInFileEpCorr = "ANY"; //(V) input file with QA histograms and EP corrections profiles

   void ReadFromFile(std::string fname = "EpQa");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdEventPlaneAllParams, 1);
};
#endif // MPDEVENTPLANEALLPARAMS_H