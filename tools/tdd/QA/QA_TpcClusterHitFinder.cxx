#include "QA_TpcClusterHitFinder.h"
#include "FairLogger.h"

ClassImp(QA_TpcClusterHitFinder);

void QA_TpcClusterHitFinder::WriteToFile()
{
   WriteBaseQA(moduleNameSuffix);
   WriteTpcClusterHitFinderQA();
}

void QA_TpcClusterHitFinder::WriteTpcClusterHitFinderQA()
{
   TString outputFilename = TString("QA_TpcClusterHitFinder_") + moduleNameSuffix + TString(".root");
   LOG(info) << " Writing QA file: " << outputFilename;

   TFile outputFile(outputFilename, "recreate");
   TTree tree("QA", "QA data tree");

   TClonesArray *clusters     = (moduleNameSuffix.EqualTo(TString("Mlem"))) ? new TClonesArray("MpdTpc2dCluster")
                                                                            : new TClonesArray("Tpc2dClusterFast");
   TClonesArray &clusters_ref = *clusters;
   tree.Branch("Clusters", &clusters, 256000, 1);

   TClonesArray *hits     = new TClonesArray("MpdTpcHit");
   TClonesArray &hits_ref = *hits;
   tree.Branch("Hits", &hits, 256000, 1);

   int eventCount = eventNumber.size();

   for (int event = 0; event < eventCount; ++event) {
      clusters_ref.Clear();
      hits_ref.Clear();
      clusters_ref = *(eventClusArray.at(event));
      hits_ref     = *(eventHitArray.at(event));
      tree.Fill();
   }

   tree.Write();
   outputFile.Close();
}
