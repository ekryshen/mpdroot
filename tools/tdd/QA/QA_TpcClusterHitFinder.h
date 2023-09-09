#ifndef QA_TPCCLUSTERHITFINDER_HH
#define QA_TPCCLUSTERHITFINDER_HH

#include "BaseQA.h"

#include "AbstractTpcHit.h"
#include "TVector3.h"

class QA_TpcClusterHitFinder : public BaseQA {

public:
   QA_TpcClusterHitFinder() {}
   virtual ~QA_TpcClusterHitFinder() {}

   void WriteToFile();
   void ReadFromFile(TString moduleSuffix, TString directory = TString(""));

   std::vector<TClonesArray *> eventClusArray;
   std::vector<TClonesArray *> eventHitArray;

   std::map<int, std::vector<int>> MCTracksFromTpcTracks(int event);

   TH2F                 *GenerateSingleCluster(int event, int clusterIndex);
   std::vector<TVector3> GetClusterHits(int event, int clusterIndex);
   std::vector<int>      GetTrackHitsIndices(int event, int tpcTrackIndex);

private:
   void WriteTpcClusterHitFinderQA();
   void ReadTpcClusterHitFinderQA(TString moduleSuffix, TString directory = TString(""));

   ClassDef(QA_TpcClusterHitFinder, 1);
};

#endif
