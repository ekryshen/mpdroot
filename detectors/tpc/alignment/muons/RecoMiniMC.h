#ifndef RECOMINIMC_H
#define RECOMINIMC_H

/// \ingroup rec
/// \class MpdTpcMiniMC
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

using namespace std;
#include <iostream>
#include <TObject.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TRandom.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <TMath.h>

// Collaborating Class Headers -------

#include "TpcSectorGeoAlignmentVAK.h"

static const Int_t RecoMiniMC_maxH = 1000;
static Int_t       RecoMiniMC_Nhits;
static Double_t    RecoMiniMC_XHits[RecoMiniMC_maxH], RecoMiniMC_YHits[RecoMiniMC_maxH],
   RecoMiniMC_ZHits[RecoMiniMC_maxH], // clusters coordinates

   RecoMiniMC_XBHits[RecoMiniMC_maxH], // clusters coordinates
   RecoMiniMC_YBHits[RecoMiniMC_maxH], // in BCS
   RecoMiniMC_ZBHits[RecoMiniMC_maxH],

   RecoMiniMC_WHits[RecoMiniMC_maxH], // clusters weights&chi2
   RecoMiniMC_chi2s[RecoMiniMC_maxH], RecoMiniMC_tHelix[RecoMiniMC_maxH];
static Int_t RecoMiniMC_iSec[RecoMiniMC_maxH], // sector number
   RecoMiniMC_qHelix;                          //    sign of the particle (+-1)

static Int_t    sectors[24]; // sectors to reconstuct
static Double_t trackParOut[6];

class RecoMiniMC : public TObject
/*
  The MpdTpcMiniMC class reads generated events with the track hits
  from the root-file and reconstructs track parameters using hits.
  The TPC gometry can be taken as from the file root-file with hits
  as from the separated root-file with the TPC alignment. The output
  file with hits and new reconstrucntion wil be wtitten.

*/
{

public:
   Int_t RecoMiniMC_IEVT;

   ///< Constructor
   RecoMiniMC() { ; }
   RecoMiniMC(BaseTpcSectorGeo &secGeo, const char *InDataFile, const char *AInFile = " ")
   {
      Int_t Event1 = 1, Event2 = 0;
      Int_t Sectors[24];
      for (Int_t i = 0; i < 24; i++) Sectors[i] = 1;
      Int_t       MinHits     = 0;
      const char *OutDir      = ".";
      const char *OutDataFile = ".";
      const char *Comment     = "";
      Init(secGeo, Event1, Event2, Sectors, MinHits, InDataFile, AInFile, OutDir, Comment, OutDataFile);
   }

   ///< Constructor 2 uses a event range
   RecoMiniMC(BaseTpcSectorGeo &secGeo, Int_t Event1, Int_t Event2, const char *InDataFile, const char *AInFile = " ")
   {
      Int_t Sectors[24];
      for (Int_t i = 0; i < 24; i++) Sectors[i] = 1;
      Int_t       MinHits     = 0;
      const char *OutDir      = ".";
      const char *OutDataFile = ".";
      const char *Comment     = "";
      Init(secGeo, Event1, Event2, Sectors, MinHits, InDataFile, AInFile, OutDir, Comment, OutDataFile);
   }

   ///< Constructor 3 uses histograms parameters
   RecoMiniMC(BaseTpcSectorGeo &secGeo, Int_t Event1, Int_t Event2, Int_t *Sectors, Int_t MinHits,
              const char *InDataFile, // full path
              const char *AInFile,    // full path
              const char *OutDir,     // full pat
              const char *Comment, const char *OutDataFile)
   { // only name
      Init(secGeo, Event1, Event2, Sectors, MinHits, InDataFile, AInFile, OutDir, Comment, OutDataFile);
   }

   ///< Destructor
   virtual ~RecoMiniMC() { ; } ///< Destructor

   void Init(BaseTpcSectorGeo &secGeo, Int_t Event1, Int_t Event2, Int_t *Sectors, Int_t MinHits,
             const char *InDataFile, // full path
             const char *AInFile,    // full path
             const char *OutDir,     // full pat
             const char *Comment,
             const char *OutDataFile); // only name
   void reco();
   void SetPrintLevel(Int_t prt) { printL = prt; }
   void getChi2(Long64_t *NhitSec, Double_t *SecChi2)
   {
      for (Int_t i = 0; i < 24; i++) {
         NhitSec[i] = fNhitSec[i];
         SecChi2[i] = fSecChi2[i];
      }
   }
   void SetEvalCut(Double_t cut) { EvalpHit = cut; }

   //______________________________________________________________________

private:
   // recalculate local parameters
   Bool_t newfit;
   // Magnetic field value (default), T
   Double_t b_mf[3];
   // Electric field direction
   Double_t e_ef[3];
   // Track start point on the beam, (cm)
   Double_t r0_beam[3];
   // Track momentum value at the start point, (GeV/c)
   Double_t p_las[3];
   // Track width on the pads plane
   Double_t d_track;
   // Minimum momentum value, (GeV/c)
   Double_t p_min;
   // Maximum momentum value, (GeV/c)
   Double_t p_max;
   // Minimum pad sensivity threshould, (cm^2)
   Double_t min_padS[2];
   // Maximum pad sensivity threshould, (cm^2)
   Double_t max_padS[2];
   // Sigma of the pad signal distribution, (%/100)
   Double_t sig_padS[2];
   // Sigma of the time signal distribution, (%/100)
   Double_t sig_timeBins;
   // Sigma of the z pozition , (%/100)
   Double_t sig_Z_pos;
   // first used in the simulation
   Int_t sector0;
   // last used in the simulation
   Int_t sector1;
   // print level
   Int_t printL = 0;
   // comment
   const char *comment;
   // calculation mode of the elevation angle of cosmic muons
   Int_t k_par; // 0 - const; 1 - "-x"; 2 - "cos{x)"; 3 - beam
   // Maximum number of pads in one cluster
   Int_t minHits;
   // Maximum number of pads in one cluster
   Int_t maxPcluster;
   // Initial seed for random generator
   UInt_t fSeed;
   // Production Mode  (0 - uniform in the space, 1 - beam)
   Bool_t beam_point;
   // Magnetic field mode  (0 - no field, 1 - field)
   Bool_t MFieldMode;
   // File name with alignment input data
   TString aInFile;
   // File name with input data
   TString InHitsFile;
   // File name with input data
   TString OutHitsFile;
   // output directory name with input data
   TString outDir;
   // File with input data
   TFile *InputMCfile;
   // File with alignment data
   TFile *alignmentF;
   // File with input data
   TFile *OutputReco;
   // trees
   TTree *old_alignment; // copy of the input Alignment for output file
   TTree *tin_alignment; // input Alignment
   TTree *tin_parameters;
   TTree *tin_data;

   TTree *tout_alignment; // "new" Alignment used at the reconstruction
   TTree *tout_parameters;
   TTree *tout_data;
   TTree *tout_fullChi2;

   // Sectors LSC shifts due to the alignment
   Double_t R0_A[3][24], old_R0_A[3][24];
   // alpha angle GSC-LSC due to the alignment
   Double_t alpha_A[24], old_alpha_A[24];
   // beta angle GSC-LSC due to the alignment
   Double_t beta_A[24], old_beta_A[24];
   // gamma angle GSC-LSC due to the alignment
   Double_t gamma_A[24], old_gamma_A[24];

   Int_t Nreconstructed;
   // full hits number
   Long64_t fNhits = 0;
   // full chi2 sum of all hits
   Double_t fSumChi2 = 0;
   // sector hits number
   Long64_t fNhitSec[24];
   // full chi2 sum of all hits
   Double_t fSecChi2[24];
   // proceeded events number
   Int_t proEvents = 0;

   TVector3 oldR0_shift[24];  // the shift of the old locfal coordinate system
                              // inside the global coordinate system
   TVector3 oldG0_gsc[24];    // the shift of the old global coordinate system
                              // inside the local coordinate system
   TRotation old_loc2glo[24]; // Rotation matrix GSC->LSC for the current sector
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

   TRotation B2G, // Rotation matrix BSC->GSC for the current sector
      G2B;        // Rotation matrix GSC->BSC for the current sector

   Int_t event1, event2 = 0; // first&last events for reconstruction

   // -----  -----          constants             --------------
   const Double_t pi        = TMath::Pi();    // sector angle
   const Double_t twopi     = TMath::TwoPi(); // sector angle
   const Double_t r2d       = 180 / pi;       // sector angle
   const Double_t phi_shift = TMath::Pi() / 12;
   const Double_t phi_sec   = TMath::Pi() / 6; // sector angle
   const Double_t epsMin    = 1e-12;           // minimal number
                                               // -----  -----

   TVector3 bMf;
   TVector3 eEf;
   TVector3 r0Beam;
   TVector3 pLas;

   Double_t trackPin[6], trackPout[6], trackDiff[6];
   // initial & reconstructed track parameters
   //  [0],[1],[2] - vertex coordinates
   //  [3],[4] - theta & phi angles (line)
   //  [3] - radius with the sign (helix)
   //  [4] - alpha_0 starting angle (helix)
   //  [5] -lambda - angle track-MF (helix)

   Double_t Chi2Fit;
   Double_t EvalpHit = 10;

   TpcSectorGeoAlignmentVAK *fTpcSecGeo;

   // -----  ----- methods --------------
   //  void Init();

   // -----  ----- fits --------------
   void  convert();
   Int_t fitLine();
   Int_t fitHelix();
   void  putChi2L();
   void  putChi2H();
   void  CircleCenter(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3, Double_t &xc,
                      Double_t &yc);

   ClassDef(RecoMiniMC, 1);
};

#endif
