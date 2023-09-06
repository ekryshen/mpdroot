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

void QA_TpcClusterHitFinder::ReadFromFile(TString suffix, TString directory)
{
   ReadBaseQA(suffix, directory);
   ReadTpcClusterHitFinderQA(suffix, directory);
}

void QA_TpcClusterHitFinder::ReadTpcClusterHitFinderQA(TString suffix, TString directory)
{
   TString inputFilename("");
   if (directory.Length() > 0) inputFilename = directory + TString("/");
   inputFilename += TString("QA_TpcClusterHitFinder_") + suffix + TString(".root");
   LOG(info) << " Reading QA file: " << inputFilename;

   TFile *f    = new TFile(inputFilename);
   TTree *tree = (TTree *)f->Get("QA");

   TClonesArray *clusters =
      (suffix.EqualTo(TString("Mlem"))) ? new TClonesArray("MpdTpc2dCluster") : new TClonesArray("Tpc2dClusterFast");
   tree->GetBranch("Clusters")->SetAutoDelete(kFALSE);
   tree->SetBranchAddress("Clusters", &clusters);

   TClonesArray *hits = new TClonesArray("MpdTpcHit");
   tree->GetBranch("Hits")->SetAutoDelete(kFALSE);
   tree->SetBranchAddress("Hits", &hits);

   int nEvents = tree->GetEntries();
   for (int i = 0; i < nEvents; ++i) {
      clusters->Clear();
      hits->Clear();

      tree->GetEntry(i);
      eventClusArray.push_back(new TClonesArray());
      eventClusArray.back() = (TClonesArray *)clusters->Clone();
      eventHitArray.push_back(new TClonesArray());
      eventHitArray.back() = (TClonesArray *)hits->Clone();
   }
}
