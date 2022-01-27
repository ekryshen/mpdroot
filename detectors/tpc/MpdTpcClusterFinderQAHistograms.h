//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Histograms for MpdTpcClusterFinderTaskQA
//
//
// Environment:
//      Software developed for the MPD at NICA.
//
// Author List:
//      Roman Salmin            (original author)
//
//
//-----------------------------------------------------------

#ifndef MPDTPCCLUSTERFINDERQAHistograms_H
#define MPDTPCCLUSTERFINDERQAHistograms_H

#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "TaskHelpers.h"

class MpdTpcClusterFinderQAHistograms : public TObject {
private:
   std::string _prefix;

public:
   TH1F *_hHitDistr;

   TH2F *_hRZ_global;
   TH2F *_hYZ_local;

   TH2F *_hXY;
   TH3F *_h3D;
   TH1F *_hSect;

   TH1F *_hErrX_inner;
   TH1F *_hErrX_outer;
   TH1F *_hErrY;
   TH1F *_hErrZ_inner;
   TH1F *_hErrZ_outer;

   TH1F *_hX;
   TH1F *_hY;
   TH1F *_hZ;

   TH2F *_hXY_global;
   TH1F *_hX_global;
   TH1F *_hY_global;
   TH1F *_hZ_global;

   TH1F *_hN;
   TH1F *_hUniqueN;
   TH1F *_hUniqueTracks;

   TH1F *_hPointNotFoundError;

   TH1F *_hGlobalDeltaX;
   TH1F *_hGlobalDeltaY;
   TH1F *_hGlobalDeltaZ;
   TH2F *_hGlobalDeltaXY;

   TH1F *_hDeltaXLocal;
   TH1F *_hDeltaYLocal;
   TH1F *_hDeltaZLocal;
   TH2F *_hDeltaXYLocal;

   TH1F *_hNumOfPadsInCluster;
   TH1F *_hNumOfTimeBinsInCluster;
   TH1F *_hNumOfDigitsInCluster;

   TH2F *_hXT_clust_row1;

   TH2F *_hDeltaXLocalVsPadRowDistance;
   TH2F *_hDeltaYLocalVsPadRowDistance;
   TH2F *_hDeltaZLocalVsPadRowDistance;

   TH2F *_hDeltaXLocalVsXInterval;
   TH2F *_hDeltaYLocalVsXInterval;
   TH2F *_hDeltaZLocalVsXInterval;

   TH1F *_hDeltaRXY;
   TH1F *_hDeltaRXYZ;

   TH2F *_hDeltaXLocalVsSector;
   TH2F *_hDeltaYLocalVsSector;
   TH2F *_hDeltaZLocalVsSector;

   TH1F *_hPeak;
   TH1F *_hPeakValeyRatio;

   // mybeg
   TH2F *_hXT_clust_row;
   TH2F *_hXT_hit_row;
   TH2F *_hXT_clust_hit_row;
   TH2F *_hXT_peak_row;
   TH2F *_hXT_collected_peak_row;

   UInt_t NumRow_hist;    // Number of PadRow that we wanna draw on the histogram
   Int_t  NumSector_hist; // Number of Sector to draw
   // myend

   MpdTpcClusterFinderQAHistograms();
   MpdTpcClusterFinderQAHistograms(const std::string &suffix);

   virtual ~MpdTpcClusterFinderQAHistograms();

   void Initialize();
   void Write();

   ClassDef(MpdTpcClusterFinderQAHistograms, 8);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------