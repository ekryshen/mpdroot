#ifndef MINITPCALIGNMENT_H
#define MINITPCALIGNMENT_H

/// \ingroup rec
/// \class miniTpcAlignment
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

using namespace std;
#include <iostream>
#include <TObject.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Mille.h>

// Collaborating Class Headers -------

#include "TpcSectorGeoAlignmentVAK.h"
#include "DerLMpdTpc.h"
#include "RecoMiniMC.h"
#include "DerHMpdTpc.h"

// Collaborating Class Declarations --

class DerHMpdTpc;

class miniTpcAlignment : public TObject
/*
  The MpdTpcMiniMC class reads generated events with the track hits
  from the root-file and finds new Allignment which provides for the
  sample the best chi2 betwen jits and a model (line or helix.
*/
{

public:
   ///< Constructor
   miniTpcAlignment() { ; }
   miniTpcAlignment(BaseTpcSectorGeo &tpcGeo, const char *InDataFile, const char *aFile)
   {
      Int_t Sectors[24] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      minNhits          = 20;
      Int_t Event1      = 1;
      Int_t Event2      = 0;
      Init(tpcGeo, Event1, Event2, minNhits, Sectors, InDataFile, aFile);
   }
   ///< Constructor 2
   miniTpcAlignment(BaseTpcSectorGeo &tpcGeo, Int_t Event1, Int_t Event2, Int_t minHits, Int_t *Sectors,
                    const char *InDataFile, const char *aFile)
   {
      Init(tpcGeo, Event1, Event2, minNhits, Sectors, InDataFile, aFile);
   }
   ///< Destructor
   virtual ~miniTpcAlignment() { ; } ///< Destructor

   void  Init(BaseTpcSectorGeo &tpcGeo, Int_t Event1, Int_t Event2, Int_t minHits, Int_t *Sectors,
              const char *InDataFile, const char *aFile);
   void  alignment();
   void  Alignment(Double_t &F, Double_t *gradF);
   void  alignmentD(Double_t *der);
   Int_t getNLV() { return NLV; }
   Int_t getNvar() { return N; }
   void  getDFDP(Double_t *der);
   void  results(Double_t *R0x, Double_t *R0y, Double_t *R0z, Double_t *alpha, Double_t *beta, Double_t *gamma);
   void  SetPrintLevel(Int_t prt) { printL = prt; }
   void  SetOffset(Bool_t OS)
   {
      anchor = OS;
      if (OS)
         Nanchor = 3;
      else
         Nanchor = 0;
   }
   void  SetLaserMode(Int_t LM) { laserM = LM; }
   Int_t GetEvent() { return CurrentEvent; }
   //  void finish();

   //______________________________________________________________________

private:
   static const Int_t maxH = 1000;
   static const Int_t NLV  = 0;
   static const Int_t NSEC = 12;
   static const Int_t NGV  = 6 * NSEC;
   static const Int_t NV   = NLV + NGV;
   Int_t              CurrentEvent;
   // anchored the 0th sector
   Bool_t anchor = false;
   // laser rays mode special mode
   Int_t    laserM = 0;
   Double_t lW[24] = {1,         0.8242372, 0.4108678, 1,         0.8242372, 0.4108678, 1, 0.8242372,
                      0.4108678, 0.4108678, 1,         0.8242372, 1,         1,         1, 1,
                      1,         1,         1,         1,         1,         1,         1, 1};
   // print level
   Int_t Nanchor = 0;
   // sectors to reconstuct
   Int_t nS, sectors[24];
   // comment
   const char *comment;
   // Maximum number of pads in one cluster
   Bool_t tLine;
   // print level
   Int_t printL = 0;
   // Maximum number of pads in one cluster
   Int_t minNhits = 0;
   // first&last events for analysis
   Int_t event1, event2;
   // File name with alignment input data
   TString aInFile;
   // File name with input data
   TString InHitsFile;
   // output directory name with input data
   TString outDir;
   // File with alignment data
   TFile *alignmentF;
   // File with input data
   TFile *InputMCfile;
   // File with input data
   TFile *OutputReco;
   // trees
   TTree *tin_alignment;
   TTree *tin_data;
   TTree *tin_parameters;
   // full hits number
   Long64_t fNhits = 0;
   // full chi2 sum of all hits
   Double_t fSumChi2 = 0;
   // sector hits number
   Long64_t fNhitSec[24];
   // full chi2 sum of all hits
   Double_t fSecChi2[24];
   // proceeded
   TTree *t_fullChi2;
   // sector position in the quetion system
   Int_t secpos[24], possec[24];

   // Sectors LSC shifts due to the alignment
   Double_t R0_A[3][24], old_R0_A[3][24];
   // alpha angle GSC-LSC due to the alignment
   Double_t alpha_A[24], old_alpha_A[24];
   // beta angle GSC-LSC due to the alignment
   Double_t beta_A[24], old_beta_A[24];
   // gamma angle GSC-LSC due to the alignment
   Double_t gamma_A[24], old_gamma_A[24];

   // Magnetic field value (default), T
   Double_t b_mf[3];

   // Branch "data" if the input file
   Double_t trackPin[6], trackPout[6], trackDiff[6];
   // initial & reconstructed track parameters
   //  [0],[1],[2] - vertex coordinates
   //  [3],[4] - theta & phi angles (line)
   //  [3] - radius with the sign (helix)
   //  [4] - alpha_0 starting angle (helix)
   //  [5] -lambda - angle track-MF (helix)
   Double_t Chi2Fit;
   Int_t    RecoMiniMC_Nhits;
   Double_t RecoMiniMC_XHits[maxH], RecoMiniMC_YHits[maxH],
      RecoMiniMC_ZHits[maxH], // clusters coordinates
      RecoMiniMC_WHits[maxH],
      RecoMiniMC_chi2s[maxH];   // clusters weights&chi2
   Int_t RecoMiniMC_iSec[maxH]; // sector number

   // proceeded events number
   Int_t proEvents;

   TVector3 oldG0_gsc[24];    // the shift of the old minhitsglobal coordinate system
                              // inside the local coordinate system
   TRotation old_glo2loc[24]; // Rotation matrix GSC->LSC for the current sector
                              // inside the local coordinate system
   TVector3 newR0_shift[24];  // the shift of the new locfal coordinate system
                              // inside the global coordinate system
   TVector3 newG0_gsc[24];    // the shift of the new global coordinate system
                              // inside the local coordinate system
   TRotation new_loc2glo[24]; // Rotation matrix LSC->GSC for the current sector
                              // inside the local coordinate system
   TRotation new_glo2loc[24]; // Rotation matrix LSC->GSC for the current sector
                              // inside the local coordinate system
   TVector3  oG0_gsc;         // Shift from GSC to LSC in GSC for the current sector
   TRotation oglo2loc;        // Rotation matrix from LSC to GSC for the current sector
   TVector3  R0shift;         // Shift from GSC to LSC in GSC for the current sector
   TRotation loc2glo;         // Rotation matrix from LSC to GSC for the current sector
   TVector3  G0_gsc;          // Shift from GSC to LSC in GSC for the current sector
   TRotation glo2loc;         // Rotation matrix from LSC to GSC for the current sector
   Double_t  dFdq[6], dFdp[6], d2Fdqdq[6][6], d2Fdqdp[6][6], d2Fdpdp[6][6], d3Fdqdp2[6][6][6], d3Fdp3[6][6][6];
   Double_t  A[NV][NV], B[NV], C[NV], X[NV], dx[NV], DFDQ[NLV], DFDP[NGV], R2[NV];
   Int_t     ldp[NV];
   Int_t     N; // number variables (NLV + NGV of involving sectors)

   // -----  -----

   TpcSectorGeoAlignmentVAK *fTpcSecGeo; // TPC geometry
   Mille                    *fMille;
   DerLMpdTpc               *fLParDer; // Line model partial derivatives
   DerHMpdTpc               *fHParDer; // Helix mjdel partial derivatives

   // -----  ----- methods --------------
   Int_t linsys(Int_t n, Double_t *x, Double_t a[NV][NV], Double_t *b);
   // -----  ----- fits --------------

   ClassDef(miniTpcAlignment, 1);
};

#endif
