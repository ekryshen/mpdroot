//--------------------------------------------------------------------
// Description:
//      Base class for generating QA plots
//      - contains global reconstruction QA's
//          Track distribution vs number of hits
//          Track distribution vs number of true hits (true - MC)
//          Matching between MC tracks & reconstructed tracks
//          Tracking efficiency
//          Pools for all track parameters
//
//      - specific module related QA's must be implemented in derived classes
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, September, 2023
//--------------------------------------------------------------------

#ifndef BASEQA_HH
#define BASEQA_HH

#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "TClonesArray.h"
#include "TH2F.h"

#include <vector>

#include "MpdTpcKalmanTrack.h"
#include "Tpc2dClusterFast.h"
#include "MpdTpc2dCluster.h"

enum EQAMode {
   OFF = 0, // QA turned off by default
   BASIC,   // QA's from this class only
   TPCCLUSTERHITFINDER
};

class BaseQA : public TObject {

public:
   BaseQA() : moduleNameSuffix(TString("")) {}
   virtual ~BaseQA() {}

   virtual void WriteToFile() { return WriteBaseQA(); }
   virtual void ReadFromFile(TString directory = TString("")) { return ReadBaseQA(TString(""), directory); }

   std::vector<int> eventNumber;
   TString          moduleNameSuffix;

   std::vector<TClonesArray *> eventMCTracksArray;
   std::vector<TClonesArray *> eventTpcTracksArray;

protected:
   void WriteBaseQA(TString suffix = TString(""));
   void ReadBaseQA(TString suffix = TString(""), TString dir = TString(""));

private:
   ClassDef(BaseQA, 1);
};

#endif
