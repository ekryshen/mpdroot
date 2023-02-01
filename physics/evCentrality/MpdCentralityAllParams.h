#ifndef MPDCENTRALITYALLPARAMS_H
#define MPDCENTRALITYALLPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdCentralityAllParams : public TObject {

public:
   //
   // Event selection cuts
   float mZvtxCut = 140.; //(V) event selection cut (cm)

   // PID cuts
   int   mNofHitsCut = 10;   //(V) minimal number of hits to accept track
   float mEtaCut     = 0.5;  //(V) maximal pseudorapidity accepted
   float mPtminCut   = 0.05; //(V) minimal pt used in analysis

   void ReadFromFile(std::string fname = "ConvDef");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdCentralityAllParams, 1);
};
#endif // MPDCENTRALITYALLPARAMS_H
