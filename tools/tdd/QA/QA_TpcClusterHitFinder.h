#ifndef QA_TPCCLUSTERHITFINDER_HH
#define QA_TPCCLUSTERHITFINDER_HH

#include "BaseQA.h"

#include "AbstractTpcHit.h"

class QA_TpcClusterHitFinder : public BaseQA {

public:
   QA_TpcClusterHitFinder() {}
   virtual ~QA_TpcClusterHitFinder() {}

   void WriteToFile();
   void ReadFromFile(TString moduleSuffix, TString directory = TString(""));

   std::vector<TClonesArray *> eventClusArray;
   std::vector<TClonesArray *> eventHitArray;

   std::map<int, std::vector<int>> MCTracksFromTpcTracks(int event);

private:
   void WriteTpcClusterHitFinderQA();
   void ReadTpcClusterHitFinderQA(TString moduleSuffix, TString directory = TString(""));

   ClassDef(QA_TpcClusterHitFinder, 1);
};

#endif
