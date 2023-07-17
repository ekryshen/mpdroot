#ifndef QA_TPCCLUSTERHITFINDER_HH
#define QA_TPCCLUSTERHITFINDER_HH

#include "BaseQA.h"

#include "AbstractTpcHit.h"

class QA_TpcClusterHitFinder : public BaseQA {

public:
   QA_TpcClusterHitFinder() {}
   virtual ~QA_TpcClusterHitFinder() {}

   void WriteToFile();

   std::vector<TClonesArray *> eventClusArray;
   std::vector<TClonesArray *> eventHitArray;

private:
   void WriteTpcClusterHitFinderQA();

   ClassDef(QA_TpcClusterHitFinder, 1);
};

#endif
