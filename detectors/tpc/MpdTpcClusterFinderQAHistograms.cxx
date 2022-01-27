//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcClusterFinderQAHistograms
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

// Panda Headers ----------------------

// This Class' Header ------------------
#include "MpdTpcClusterFinderQAHistograms.h"
#include "MpdTpcSector.h"
#include <sstream>

ClassImp(MpdTpcClusterFinderQAHistograms);

MpdTpcClusterFinderQAHistograms::MpdTpcClusterFinderQAHistograms()
   : _hHitDistr(0), _hRZ_global(0), _hYZ_local(0), _hXY(0), _h3D(0), _hSect(0), _hErrX_inner(0), _hErrX_outer(0),
     _hErrY(0), _hErrZ_inner(0), _hErrZ_outer(0), _hX(0), _hY(0), _hZ(0), _hXY_global(0), _hX_global(0), _hY_global(0),
     _hZ_global(0), _hN(0), _hUniqueN(0), _hUniqueTracks(0), _hPointNotFoundError(0), _hGlobalDeltaX(0),
     _hGlobalDeltaY(0), _hGlobalDeltaZ(0), _hGlobalDeltaXY(0), _hDeltaXLocal(0), _hDeltaYLocal(0), _hDeltaZLocal(0),
     _hDeltaXYLocal(0), _hNumOfPadsInCluster(0), _hNumOfTimeBinsInCluster(0), _hNumOfDigitsInCluster(0),
     _hXT_clust_row1(0), _hDeltaXLocalVsPadRowDistance(0), _hDeltaYLocalVsPadRowDistance(0),
     _hDeltaZLocalVsPadRowDistance(0), _hDeltaXLocalVsXInterval(0), _hDeltaYLocalVsXInterval(0),
     _hDeltaZLocalVsXInterval(0), _hDeltaRXY(0), _hDeltaRXYZ(0), _hDeltaXLocalVsSector(0), _hDeltaYLocalVsSector(0),
     _hDeltaZLocalVsSector(0), _hPeak(0), _hPeakValeyRatio(0), _hXT_clust_row(0), _hXT_hit_row(0),
     _hXT_clust_hit_row(0), _hXT_peak_row(0), _hXT_collected_peak_row(0)
{
}

MpdTpcClusterFinderQAHistograms::MpdTpcClusterFinderQAHistograms(const std::string &prefix)
   : _prefix(prefix), _hHitDistr(0), _hRZ_global(0), _hYZ_local(0), _hXY(0), _h3D(0), _hZ(0), _hN(0), _hUniqueN(0),
     _hUniqueTracks(0), _hPointNotFoundError(0), _hGlobalDeltaX(0), _hGlobalDeltaY(0), _hGlobalDeltaZ(0),
     _hGlobalDeltaXY(0), _hDeltaXLocal(0), _hDeltaYLocal(0), _hDeltaZLocal(0), _hDeltaXYLocal(0),
     _hDeltaXLocalVsPadRowDistance(0), _hDeltaYLocalVsPadRowDistance(0), _hDeltaZLocalVsPadRowDistance(0),
     _hDeltaXLocalVsXInterval(0), _hDeltaYLocalVsXInterval(0), _hDeltaZLocalVsXInterval(0), _hDeltaRXY(0),
     _hDeltaRXYZ(0), _hDeltaXLocalVsSector(0), _hDeltaYLocalVsSector(0), _hDeltaZLocalVsSector(0), _hPeak(0),
     _hPeakValeyRatio(0)
{
}

MpdTpcClusterFinderQAHistograms::~MpdTpcClusterFinderQAHistograms()
{

   delete _hRZ_global;
   delete _hYZ_local;
   delete _hXY;
   delete _hXY_global;
   delete _h3D;
   delete _hErrX_inner;
   delete _hErrX_outer;
   delete _hErrY;
   delete _hErrZ_inner;
   delete _hErrZ_outer;
   delete _hX;
   delete _hY;
   delete _hZ;
   delete _hX_global;
   delete _hY_global;
   delete _hZ_global;
   delete _hSect;
   delete _hNumOfPadsInCluster;
   delete _hNumOfDigitsInCluster;
   delete _hNumOfTimeBinsInCluster;
   delete _hHitDistr;

   delete _hN;
   delete _hUniqueN;
   delete _hUniqueTracks;
   delete _hPointNotFoundError;
   delete _hGlobalDeltaX;
   delete _hGlobalDeltaY;
   delete _hGlobalDeltaZ;
   delete _hGlobalDeltaXY;
   delete _hDeltaXLocal;
   delete _hDeltaYLocal;
   delete _hDeltaZLocal;
   delete _hDeltaXLocalVsSector;
   delete _hDeltaYLocalVsSector;
   delete _hDeltaZLocalVsSector;
   delete _hDeltaXLocalVsPadRowDistance;
   delete _hDeltaYLocalVsPadRowDistance;
   delete _hDeltaZLocalVsPadRowDistance;
   delete _hDeltaXLocalVsXInterval;
   delete _hDeltaYLocalVsXInterval;
   delete _hDeltaZLocalVsXInterval;
   delete _hDeltaXYLocal;

   delete _hXT_clust_row1;

   delete _hDeltaRXY;
   delete _hDeltaRXYZ;

   delete _hPeak;
   delete _hPeakValeyRatio;

   // mybeg
   delete _hXT_clust_hit_row;
   delete _hXT_clust_row;
   delete _hXT_hit_row;
   delete _hXT_peak_row;
   delete _hXT_collected_peak_row;
   // myend
}

void MpdTpcClusterFinderQAHistograms::Initialize()
{

   // TODO: GET PARAMETERS FROM DATA BASE
   Float_t zMax  = 170; // cm
   UInt_t  nSect = 12;
   Float_t maxX  = 30.73; // cm
   Float_t maxY  = 79.2;  // cm
   Float_t Rmax  = 121.0;

   // mybeg
   NumRow_hist    = 1;
   NumSector_hist = 3; // if "NumSector_hist = -1" : Show all sectors

   UInt_t     nTimeBin = 512;
   TpcSector *sector   = new TpcSector();
   // UInt_t nSectors = sector->GetNSectors();
   UInt_t nRows             = sector->GetNumRows();
   UInt_t nPads             = sector->GetNumPadsInRow(nRows - 1) * 2;
   UInt_t nPadsInCurrentRow = sector->GetNumPadsInRow(NumRow_hist) * 2;
   delete sector;
   // myend

   _hYZ_local  = CreateHistogram2("local ZY distribution of Found Hits", _prefix, "Z, cm", "Y, cm", 500, 0.0, zMax, 500,
                                  0.0, maxY);
   _hRZ_global = CreateHistogram2("global ZR distribution of Found Hits", _prefix, "Z, cm", "R, cm", 800, -zMax, zMax,
                                  800, -Rmax, Rmax);

   _h3D = CreateHistogram3("Found Clusters 3D", _prefix, "X (cm)", "Y (cm)", "Z (cm)", 200, -150, 150, 200, -150, 150,
                           200, -zMax, zMax);

   _hHitDistr = CreateHistogram1("Hits distribution", _prefix, "Number of found hits per PadRow", 10, 0, 10);

   _hErrZ_inner =
      CreateHistogram1("Distribution of Found Hits Z Errors in inner sector", _prefix, "Z error (cm)", 300, 0, 1);
   _hErrZ_outer =
      CreateHistogram1("Distribution of Found Hits Z Errors in outer sector", _prefix, "Z error (cm)", 300, 0, 1);
   _hErrY = CreateHistogram1("Distribution of Found Hits Y Errors", _prefix, "Y error (cm)", 300, 0, 1);
   _hErrX_outer =
      CreateHistogram1("Distribution of Found Hits X Errors in outer sector", _prefix, "X error (cm)", 300, 0, 1);
   _hErrX_inner =
      CreateHistogram1("Distribution of Found Hits X Errors in inner sector", _prefix, "X error (cm)", 300, 0, 1);

   _hXT_clust_row1 =
      CreateHistogram2("XT distribution of TPC hits on PadRow #1", "X, cm", "Z, cm", 200, -30, 30, 200, 0, zMax);

   _hXY =
      CreateHistogram2("XY distribution of Found Hits", _prefix, "X (cm)", "Y (cm)", 800, -maxX, maxX, 800, 0.0, maxY);
   _hZ = CreateHistogram1("local Z distribution of Found Hits", _prefix, "Z (cm)", 300, 0.0, zMax);
   _hY = CreateHistogram1("local Y distribution of Found Hits", _prefix, "Y (cm)", 100, 0, maxY);
   _hX = CreateHistogram1("local X distribution of Found Hits", _prefix, "X (cm)", 100, -maxX, maxX);

   _hXY_global = CreateHistogram2("global XY distribution of Found Hits", _prefix, "X (cm)", "Y (cm)", 800, -150, 150,
                                  800, -150, 150);
   _hZ_global  = CreateHistogram1("global Z distribution of Found Hits", _prefix, "Z (cm)", 200, -zMax, zMax);
   _hY_global  = CreateHistogram1("global Y distribution of Found Hits", _prefix, "Y (cm)", 200, -150, 150);
   _hX_global  = CreateHistogram1("global X distribution of Found Hits", _prefix, "X (cm)", 200, -150, 150);

   _hSect = CreateHistogram1("Sector distribution of Found Hits", _prefix, "#Sect", nSect, 0, nSect);

   _hNumOfPadsInCluster     = CreateHistogram1("Pads distribution of clusters", _prefix, "Pads", 20, 0, 20);
   _hNumOfTimeBinsInCluster = CreateHistogram1("Time bins distribution of clusters", _prefix, "Time bins", 50, 0, 50);
   _hNumOfDigitsInCluster   = CreateHistogram1("Cluster size distribution", _prefix, "Digits", 100, 0, 100);

   _hPeak = CreateHistogram1("Charge distribution of found Hits", _prefix, "ADC Value", 200, 0, 2000000);

   // mybeg
   _hXT_clust_hit_row = CreateHistogram2("XT distribution of 2D-clusters and hits on PadRow #", "X, pads",
                                         "T, time bins", nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_clust_row = CreateHistogram2("_XT distribution of 2D-clusters on PadRow #", "X, pads", "T, time bins", nPads, 0,
                                     nPads, nTimeBin, 0, nTimeBin);
   _hXT_hit_row   = CreateHistogram2("_XT distribution of hits on PadRow #", "X, pads", "T, time bins", nPads, 0, nPads,
                                     nTimeBin, 0, nTimeBin);
   _hXT_peak_row  = CreateHistogram2("(tmp)_XT distribution of Peaks on PadRow #", "X, pads", "T, time bins", nPads, 0,
                                     nPads, nTimeBin, 0, nTimeBin);
   _hXT_collected_peak_row = CreateHistogram2("(tmp)_XT distribution of Collected Peaks on PadRow #", "X, pads",
                                              "T, time bins", nPads, 0, nPads, nTimeBin, 0, nTimeBin);

   // substring to set hist title for current row and sector
   TString hSubTitle_hist;
   hSubTitle_hist = " ";
   hSubTitle_hist += NumRow_hist;
   hSubTitle_hist += " (Sector#";
   if (NumSector_hist == -1)
      hSubTitle_hist += "All";
   else
      hSubTitle_hist += NumSector_hist;
   hSubTitle_hist += ")";

   _hXT_clust_hit_row->SetTitle(_hXT_clust_hit_row->GetTitle() + hSubTitle_hist);
   _hXT_clust_row->SetTitle(_hXT_clust_row->GetTitle() + hSubTitle_hist);
   _hXT_hit_row->SetTitle(_hXT_hit_row->GetTitle() + hSubTitle_hist);
   _hXT_peak_row->SetTitle(_hXT_peak_row->GetTitle() + hSubTitle_hist);
   _hXT_collected_peak_row->SetTitle(_hXT_collected_peak_row->GetTitle() + hSubTitle_hist);

   _hXT_clust_hit_row->SetBins(nPadsInCurrentRow, 0, nPadsInCurrentRow, nTimeBin, 0, nTimeBin);
   _hXT_clust_row->SetBins(nPadsInCurrentRow, 0, nPadsInCurrentRow, nTimeBin, 0, nTimeBin);
   _hXT_hit_row->SetBins(nPadsInCurrentRow, 0, nPadsInCurrentRow, nTimeBin, 0, nTimeBin);
   _hXT_peak_row->SetBins(nPadsInCurrentRow, 0, nPadsInCurrentRow, nTimeBin, 0, nTimeBin);
   _hXT_collected_peak_row->SetBins(nPadsInCurrentRow, 0, nPadsInCurrentRow, nTimeBin, 0, nTimeBin);
   ;

   _hXT_clust_row->SetOption("BOX");
   _hXT_hit_row->SetOption("BOX");
   _hXT_peak_row->SetOption("BOX");
   _hXT_collected_peak_row->SetOption("BOX");

   _hXT_clust_row->SetMaximum(1);
   _hXT_hit_row->SetMaximum(1);
   _hXT_peak_row->SetMaximum(1);
   _hXT_collected_peak_row->SetMinimum(1);

   _hXT_clust_row->SetFillColor(kBlue);
   _hXT_hit_row->SetFillColor(kRed);
   _hXT_peak_row->SetFillColor(kGreen);
   _hXT_collected_peak_row->SetFillColor(kOrange);

   _hXT_clust_hit_row->SetOption("COL");
   _hXT_clust_hit_row->SetContour(5);
   _hXT_clust_hit_row->SetContourLevel(1, 100);
   _hXT_clust_hit_row->SetContourLevel(2, 101);
   _hXT_clust_hit_row->SetContourLevel(3, 102);
   _hXT_clust_hit_row->SetContourLevel(4, 103);
   // myend
}

void MpdTpcClusterFinderQAHistograms::Write()
{
   _hRZ_global->Write(0, kOverwrite);
   _hYZ_local->Write(0, kOverwrite);
   _hHitDistr->Write(0, kOverwrite);
   _hXY->Write(0, kOverwrite);
   _hXY_global->Write(0, kOverwrite);
   _h3D->Write(0, kOverwrite);
   _hX->Write(0, kOverwrite);
   _hY->Write(0, kOverwrite);
   _hZ->Write(0, kOverwrite);
   _hX_global->Write(0, kOverwrite);
   _hY_global->Write(0, kOverwrite);
   _hZ_global->Write(0, kOverwrite);
   _hSect->Write(0, kOverwrite);
   _hErrX_inner->Write(0, kOverwrite);
   _hErrX_outer->Write(0, kOverwrite);
   _hErrY->Write(0, kOverwrite);
   _hErrZ_inner->Write(0, kOverwrite);
   _hErrZ_outer->Write(0, kOverwrite);
   _hNumOfPadsInCluster->Write(0, kOverwrite);
   _hNumOfTimeBinsInCluster->Write(0, kOverwrite);
   _hNumOfDigitsInCluster->Write(0, kOverwrite);
   _hPeak->Write(0, kOverwrite);
   _hXT_clust_row1->Write(0, kOverwrite);

   // _hN->Write(0, kOverwrite);
   // _hUniqueN->Write(0, kOverwrite);
   // _hUniqueTracks->Write(0, kOverwrite);
   //_hPointNotFoundError->Write(0, kOverwrite);

   //_hGlobalDeltaX->Write(0, kOverwrite);
   //_hGlobalDeltaY->Write(0, kOverwrite);
   //_hGlobalDeltaZ->Write(0, kOverwrite);
   //_hGlobalDeltaXY->Write(0, kOverwrite);

   //_hDeltaRXY->Write(0, kOverwrite);
   //_hDeltaRXYZ->Write(0, kOverwrite);

   //_hDeltaXLocal->Write(0, kOverwrite);
   //_hDeltaYLocal->Write(0, kOverwrite);
   //_hDeltaZLocal->Write(0, kOverwrite);
   //_hDeltaXLocalVsSector->Write(0, kOverwrite);
   //_hDeltaYLocalVsSector->Write(0, kOverwrite);
   //_hDeltaZLocalVsSector->Write(0, kOverwrite);
   //_hDeltaXYLocal->Write(0, kOverwrite);

   //_hDeltaXLocalVsPadRowDistance->Write(0, kOverwrite);
   //_hDeltaYLocalVsPadRowDistance->Write(0, kOverwrite);
   //_hDeltaZLocalVsPadRowDistance->Write(0, kOverwrite);

   //_hDeltaXLocalVsXInterval->Write(0, kOverwrite);
   //_hDeltaYLocalVsXInterval->Write(0, kOverwrite);
   //_hDeltaZLocalVsXInterval->Write(0, kOverwrite);

   // _hPeakValeyRatio->Write(0, kOverwrite);

   // mybeg
   _hXT_clust_hit_row->Write(0, kOverwrite);
   _hXT_clust_row->Write(0, kOverwrite);
   _hXT_hit_row->Write(0, kOverwrite);
   _hXT_peak_row->Write(0, kOverwrite);
   _hXT_collected_peak_row->Write(0, kOverwrite);
   // myend
}