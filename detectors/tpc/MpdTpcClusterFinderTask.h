//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Takes ADC signals and find primary cluster
//
//
// Environment:
//      Software developed for the MPD at NICA.
//
// Author List:
//      Sergey Merts
//      Dmitry Sosnov
//      Artem Basalaev
//
//
//-----------------------------------------------------------

#ifndef MPDTPCCLUSTERFINDERTASK_HH
#define MPDTPCCLUSTERFINDERTASK_HH

#include "TMath.h"
#include "FairTask.h"
#include "MpdTpcClusterFinderQAHistograms.h"
#include "MpdTpcDigitizerTask.h"
#include "MpdTpc2dCluster.h"
#include "TRefArray.h"
#include "MpdTpcPeak.h"
#include "MpdTpcFoundHit.h"
#include "MpdTpcDigit.h"
#include "TpcGas.h"

using namespace std;

class TpcGas;
class TClonesArray;

class MpdTpcClusterFinderTask : public FairTask {
public:
   MpdTpcClusterFinderTask();
   virtual ~MpdTpcClusterFinderTask();

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }
   void SetDebug(Bool_t opt = kTRUE) { fPrintDebugInfo = opt; }
   void SetMakeQA(Bool_t opt = kTRUE) { fMakeQA = opt; }
   void SetCalcResiduals(Bool_t opt = kTRUE) { fCalcResiduals = opt; }
   void SetResFileName(TString opt) { fNameResFile = opt; }

   // Operations ----------------------
   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();

   Bool_t Find2DClusters(vector<MpdTpc2dCluster *> *extClust, DigOrigArray **digArr, UInt_t row, UInt_t sec);
   void   FindPeaksInCluster(MpdTpc2dCluster *clust, vector<MpdTpcPeak *> *peakList, DigOrigArray **fOriginsArray);
   // void FindPeaksInClusterNew(MpdTpc2dCluster* clust, vector<MpdTpcPeak*> *peakList, DigOrigArray** fOriginsArray);
   void   CollectPeaks(vector<MpdTpcPeak *> peakList, MpdTpc2dCluster *clust, vector<MpdTpcPeak *> *collectedPeakList);
   void   CreateHit(vector<MpdTpcPeak *> collectedPeakList, MpdTpc2dCluster *clust, vector<MpdTpcFoundHit *> *hitList,
                    DigOrigArray **fOriginsArray);
   Bool_t GetNextDigit(UInt_t *currDig, MpdTpc2dCluster *Clus2d, DigOrigArray **digArr, Bool_t **fADCMarks);
   void   FitPeaks(vector<MpdTpcPeak *> PeakList, Int_t npeaks[], Int_t mincol);
   void LinGammaFit(Float_t lnq[], Float_t w[], Float_t dt[], Float_t lndt[], Int_t nbin, Float_t par[], Float_t epar[],
                    Float_t &chi2);
   void CalcResiduals(MpdTpcFoundHit *hit, Float_t phi);
   void CalcNewErrors(MpdTpcFoundHit *hit);

private:
   // Private Data Members ------------
   TString inputBranchName;
   TString outputBranchName;

   Int_t OnePeakCntr, MoreThenOnePeakCntr, OneDigitHitCntr, AllCntr, ThreePadsCntr; // TMP

   TpcGas *        fGas;          // pointer to gas system
   TClonesArray *  fDigits;       // input TClonesArray of TPC digits
   TClonesArray *  fHitsArray;    // output TClonesArray of TPC hits
   TClonesArray *  fMCPointArray; // input TClonesArray of MC points. It's needed only for residuals calculating
   TClonesArray *  fMCTracks;     // input TClonesArray of MC tracks. It's needed only for residuals calculating
   DigOrigArray ***fDigitsArray;  // input 4D-array of TPC digits
   MpdTpcClusterFinderQAHistograms *fHisto; // pointer to object needed for QA creating

   // set of boolean flags for manage of work process
   Bool_t fPersistence;            // print or not output array into tree
   Bool_t fPrintDebugInfo;         // print or not additional information in output
   Bool_t isHistogramsInitialized; // is QA histograms initialized or not
   Bool_t fMakeQA;                 // create or not in output tree branch with QA histograms
   Bool_t fCalcResiduals;          // calculate residuals of TPC hits or not
   Bool_t fitGamma;                // not used now

   Float_t fSpread;         // sigma for pad response function
   Float_t zDrift;          // length of TPC
   Float_t r_min;           // minimal radius of TPC
   UInt_t  fEventId;        // event number
   UInt_t  nSect;           // number of sectors in TPC
   UInt_t  nTimeBins;       // number of time samples in TPc
   Float_t pwIn, phIn;      // inner pad width and height
   Float_t pwOut, phOut;    // outer pad width and height
   UInt_t  nRows;           // number of rows in TPC sector
   UInt_t  nInRows;         // number of rows in inner part of TPC sector
   UInt_t  nOutRows;        // number of rows in outer part of TPC sector
   UInt_t *fNumOfPadsInRow; // array of pads quantity in every TPC rows
   Float_t fSectInHeight;   // height of inner part of TPC sector
   Float_t fNoiseThreshold; // threshold for signal separation
   TString fNameResFile;    // name of output residuals file

public:
   ClassDef(MpdTpcClusterFinderTask, 6);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------