//--------------------------------------------------------------------
// Description:
//      QA for ClusterHitFinder modules
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, March, 2023
//--------------------------------------------------------------------

#ifndef QA_TPCCLUSTERHITFINDER_HH
#define QA_TPCCLUSTERHITFINDER_HH

#include "AbstractQA.h"

class QA_TpcClusterHitFinder : public AbstractQA {

public:
   QA_TpcClusterHitFinder() {}
   virtual ~QA_TpcClusterHitFinder() {}

private:
   ClassDef(QA_TpcClusterHitFinder, 1);
};

#endif
