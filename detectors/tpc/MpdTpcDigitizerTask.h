//--------------------------------------------------------------------
//
// Description:
//      Tpc Digitizer reads array of MC points and produces TpcDigits
//
//
// Author List:
//      Sergey Merts
//
//--------------------------------------------------------------------

#ifndef MPDTPCDIGITIZERTASK_HH
#define MPDTPCDIGITIZERTASK_HH

// Base Class Headers ----------------

#include "FairTask.h"
#include "MpdMCTrack.h"
#include "FairField.h"
#include "TpcPoint.h"
#include "MpdTpcDigitizerQAHistograms.h"

class TpcGas;

using namespace std;

struct DigOrigArray {
   Float_t             signal;
   map<Int_t, Float_t> origins; // array of pairs <MC-track ID, its signal> for calculating origin
   Int_t               origin;  // MC-track index for testing
   Bool_t              isOverlap;
};

class MpdTpcDigitizerTask : public FairTask {
public:
   // Constructors/Destructors ---------
   MpdTpcDigitizerTask();
   virtual ~MpdTpcDigitizerTask();

   Bool_t isSubtrackInInwards(const TpcPoint *p1, const TpcPoint *p2);

   void SetPrimBranchName(const TString &name) { fInputBranchName = name; }
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }
   void SetAttach(Bool_t opt = kTRUE) { fAttach = opt; }
   void SetDiffuse(Bool_t opt = kTRUE) { fDiffuse = opt; }
   void SetDistort(Bool_t opt = kTRUE) { fDistort = opt; }
   void SetDebug(Bool_t opt = kTRUE) { fPrintDebugInfo = opt; }
   void SetMakeQA(Bool_t opt = kFALSE) { fMakeQA = opt; }
   void SetDistribute(Bool_t opt = kFALSE) { fDistribute = opt; }
   void SetResponse(Bool_t opt = kFALSE) { fResponse = opt; }
   void SetOnlyPrimary(Bool_t opt = kFALSE) { fOnlyPrimery = opt; }

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();

   void    PadResponse(Float_t x, Float_t y, UInt_t timeID, Int_t origin, DigOrigArray ***arr);
   TF1 *   padResponseFunction();
   Float_t CalculatePadResponse(UInt_t padID, UInt_t rowID, Float_t x, Float_t y);
   void    GetArea(Float_t xEll, Float_t yEll, Float_t radius, vector<UInt_t> &padIDs, vector<UInt_t> &rowIDs);
   Int_t   CalcOrigin(const DigOrigArray dig);
   void    TpcProcessing(const TpcPoint *prePoint, const TpcPoint *curPoint, const UInt_t secID, const UInt_t iPoint,
                         const UInt_t nPoints);

private:
   // Private Data Members ------------
   TString fInputBranchName;
   TString fOutputBranchName;

   TClonesArray *               fMCPointArray;  // input array of MC points
   TClonesArray *               fMCTracksArray; // input array of MC tracks
   TClonesArray *               fDigits;        // output array of TPC digits (write into output tree)
   TpcGas *                     fGas;           // pointer to gas system
   DigOrigArray ***             fDigits4dArray; // output array of digital signals (don't write into output tree)
   TpcSector *                  fSector;        // object for getting geometrical parameters of TPC sector
   MpdTpcDigitizerQAHistograms *fHisto;         // pointer to object needed only for QA creating
   TF1 *                        fPRF;           // pad response function
   FairField *                  fMagField;      // magnetic field

   Float_t fGain;           // coefficient for avalanches calculating
   Float_t zCathode;        // length of TPC
   Float_t fNoiseThreshold; // threshold for signal separation
   Int_t * fNumOfPadsInRow; // array of pads quantity in every TPC rows
   UInt_t  nSectors;        // number of TPC sectors
   UInt_t  nTimeBackets;    // number of time samples in TPC
   UInt_t  nRows;           // number of rows in TPC sector
   UInt_t  nInRows;         // number of rows in inner part of TPC sector
   UInt_t  nOutRows;        // number of rows in outer part of TPC sector
   Float_t fSectInHeight;   // height of inner part of TPC sector
   Float_t fSectHeight;     // height of TPC sector
   Float_t r_min;           // minimal radius of TPC
   Float_t fSpread;         // sigma for pad response function
   Float_t k1, k2;          // coefficients for padRespose calculating
   Float_t pwIn;            // inner pad width
   Float_t phIn;            // inner pad height
   Float_t pwOut;           // outer pad width
   Float_t phOut;           // outer pad height

   // set of boolean flags for manage of work process
   Bool_t fIsHistogramsInitialized; // is QA histograms initialized or not
   Bool_t fMakeQA;                  // create or not in output tree branch with QA histograms
   Bool_t fOnlyPrimery;             // Take into account only primery particle or not
   Bool_t fPersistence;             // print or not output array into tree
   Bool_t fAttach;                  // attach electrons in gas or not
   Bool_t fDiffuse;                 // diffuse electrons in TPC  or not
   Bool_t fDistort;                 // not implemented yet
   Bool_t fResponse;                // to do pad response or not
   Bool_t fDistribute;              // distribute electrons between two MC points or not
   Bool_t fPrintDebugInfo;          // print or not additional information in output

public:
   ClassDef(MpdTpcDigitizerTask, 5);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
