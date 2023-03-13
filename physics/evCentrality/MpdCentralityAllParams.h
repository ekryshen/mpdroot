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
   float mZvtxCut = 130.; //(V) event selection cut (cm)

   // Track selection cuts
   int   mNofHitsCut = 10;  //(V) minimal number of hits to accept track
   float mEtaCut     = 0.5; //(V) maximal pseudorapidity accepted
   float mPtminCut   = 0.1; //(V) minimal pt used in analysis
   float mDcaCut     = 2.0; //(V) maximal DCA accepted

   std::string mProdGenerator = "ANY"; //(V) production and event generator
   std::string mInFileConvert = "ANY"; //(V) input file name with track-to-centrality converter
   std::string mInFileTrEff   = "ANY"; //(V) input file name with track reconstruction efficiecnies

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
