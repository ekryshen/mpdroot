#include "QA_TpcClusterHitFinder.h"
#include "MpdTpcHit.h"
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
   moduleNameSuffix = suffix;
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

//-------------------------------------------------------------------------------------------------
// Method to identify disconnected tracks
// For given event generates map of Monte-Carlo ID, vector of TPC track ID's
// NOTE: MC track with major contribution is taken.
//       There can be contributions from multiple MC tracks
//       - then such method must return std::unordered_map with key
//       std::vector<std::pair<MC_ID, total charge contribution from MC_ID>>

std::map<int, std::vector<int>> QA_TpcClusterHitFinder::MCTracksFromTpcTracks(int event)
{
   std::map<int, std::vector<int>> idMap;

   int nTpcTracks = eventTpcTracksArray[event]->GetEntriesFast();
   // LOG(info) << "Event: " << event << "    Number of TPC tracks: " << nTpcTracks;

   for (int iTrack = 0; iTrack < nTpcTracks; ++iTrack) {

      MpdTpcKalmanTrack *tpcTrack   = (MpdTpcKalmanTrack *)eventTpcTracksArray[event]->UncheckedAt(iTrack);
      int                nTrackHits = tpcTrack->GetNofTrHits();
      TClonesArray      *trackHits  = tpcTrack->GetTrHits();
      // LOG(info) << "   Track: " << iTrack << "    Number of Hits: " << nTrackHits;

      std::map<int, int> mcIDs; // key: Monte Carlo ID, value: its count from digits

      for (int h = 0; h < nTrackHits; ++h) {
         MpdKalmanHit *trackHit     = (MpdKalmanHit *)trackHits->UncheckedAt(h);
         int           hitIndex     = trackHit->GetIndex();
         MpdTpcHit    *hit          = (MpdTpcHit *)eventHitArray[event]->UncheckedAt(hitIndex);
         int           clusterIndex = hit->GetClusterID();

         AbstractTpc2dCluster *cluster;
         if (moduleNameSuffix.EqualTo("Mlem"))
            cluster = (MpdTpc2dCluster *)eventClusArray[event]->UncheckedAt(clusterIndex);
         if (moduleNameSuffix.EqualTo("Fast"))
            cluster = (Tpc2dClusterFast *)eventClusArray[event]->UncheckedAt(clusterIndex);
         std::vector<AbstractTpcDigit *> clusterDigits = cluster->GetClusterDigits();
         int                             digitsNumber  = clusterDigits.size();

         for (int d = 0; d < digitsNumber; ++d) {
            std::map<int, float> signalMap = clusterDigits.at(d)->GetTrackSignals();
            int mainSignal = signalMap.begin()->first; // currently, MC track with major contribution to digit is stored
            mcIDs[mainSignal]++;
         }
      }

      auto mcMaxSignal =
         std::max_element(mcIDs.begin(), mcIDs.end(), [](const auto &x, const auto &y) { return x.second < y.second; });
      // LOG(info) << "Track max MC signal ID: " << mcMaxSignal->first;

      std::map<int, std::vector<int>>::iterator it = idMap.find(mcMaxSignal->first);
      if (it != idMap.end())
         it->second.push_back(iTrack);
      else
         idMap.insert(make_pair(mcMaxSignal->first, std::vector<int>{iTrack}));
   }

   return idMap;
}

//-------------------------------------------------------------------------------------------------