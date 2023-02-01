#ifndef TPCMINIMC_H
#define TPCMINIMC_H

/// \ingroup rec
/// \class TpcMiniMC
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
#include <TRandom.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <TMath.h>
// Collaborating Class Headers -------
#include "TpcSectorGeoAlignmentVAK.h"

static const Int_t MpdTpcMiniMC_maxH = 1000;
static Int_t       MpdTpcMiniMC_Nhits;
static Double_t    MpdTpcMiniMC_XHits[MpdTpcMiniMC_maxH],
   MpdTpcMiniMC_YHits[MpdTpcMiniMC_maxH], // clusters coordinates
   MpdTpcMiniMC_ZHits[MpdTpcMiniMC_maxH],

   MpdTpcMiniMC_XBHits[MpdTpcMiniMC_maxH], // clusters coordinates
   MpdTpcMiniMC_YBHits[MpdTpcMiniMC_maxH], // in BCS
   MpdTpcMiniMC_ZBHits[MpdTpcMiniMC_maxH],

   MpdTpcMiniMC_WHits[MpdTpcMiniMC_maxH],          // clusters weights
   MpdTpcMiniMC_chi2s[MpdTpcMiniMC_maxH];          // clusters chi2
static Int_t MpdTpcMiniMC_iSec[MpdTpcMiniMC_maxH], // sector number
   MpdTpcMiniMC_qHelix;                            //    sign of the particle (+-1)

static const Int_t mxp = 124, mxr = 53;
// indexPad shows which pad band croses the track
static Int_t indexPad[mxp][mxr][3];
// 1000 - left, 0100 - upper, 0010 - right, 0001 - lower
static Double_t SPads[mxp], xc1[mxp][mxr][2], yc1[mxp][mxr][2], xc2[mxp][mxr][2], yc2[mxp][mxr][2], xminPad[mxp][mxr],
   yminPad[mxp][mxr], xmaxPad[mxp][mxr], ymaxPad[mxp][mxr];

class TpcMiniMC : public TObject
/*
  The TpcMiniMC class generates single tracks in TPC with and
  without magnetic field, which can start in any point inside the
  TPC or crosses it. The result of the simulation is an event with
  track hits.
*/
{

public:
   ///< Constructor 1 uses the default allignment and parameters
   TpcMiniMC(BaseTpcSectorGeo &tpcGeo, Int_t mode, Int_t Chambers)
   // mode://  ROOT::Minuit2::Minuit2Minimizer minFit;

   //  0  -  magnetic field is on, particle
   //        production in the interaction point
   //        using default parameters mode == 0 || mode == 2
   //  1  -  magnetic field is on, uniform particle
   //        production inside the chamber
   //        using default parameters
   //  2  -  magnetic field is off, particle
   //        production in the interaction point
   //        using default parameters
   //  3  -  magnetic field is off, uniform particle
   //        production inside the chamber
   //        using default parameters
   // chambers:
   //  0  -  both chambers 0
   //  1  -  only chamber 0
   //  2  -  only chamber 1
   {
      MFieldMode = mode < 2;               // true if the Magnetic field is on
      beam_point = mode == 0 || mode == 2; // true if a beam
      chambers   = Chambers;
      if (Chambers == 0) {
         sector0 = 0;
         sector1 = 24;
      } else if (Chambers == 1) {
         sector0 = 0;
         sector1 = 12;
      } else if (Chambers == 2) {
         sector0 = 0;
         sector1 = 24;
      } else
         cout << " MiniMC: no chambers mode " << Chambers << endl;
      Init(tpcGeo);
   }

   //______________________________________________________________________

   ///< Constructor 2 difines all parameters
   TpcMiniMC(BaseTpcSectorGeo &tpcGeo, Double_t *B_mf, Bool_t Beam_point, Double_t *E_ef, Double_t *P_las,
             Double_t D_track, Double_t *Min_padS, Double_t *Max_padS, Double_t *Sig_padS, Double_t Sig_timeBins,
             Double_t Sig_Z_pos, Double_t *R0_beam, Double_t K_par, Double_t P_min, Double_t P_max, Int_t Chambers,
             Int_t MaxPcluster, Int_t MinHits, UInt_t Seed, TString *aFileIn, TString *OutDir, TString *OutFile);

   //______________________________________________________________________

   ///< Constructor 3 difines all parameters
   TpcMiniMC(BaseTpcSectorGeo &tpcGeo, UInt_t Seed, TString *ParFile, TString *aFileIn, TString *OutDir,
             TString *OutFile);

   ///< Destructor
   virtual ~TpcMiniMC() { ; } ///< Destructor

   void  Init(BaseTpcSectorGeo &tpcGeo); // write header and open data files
   void  simulate(Int_t nEvents);
   Int_t GetNhits() { return MpdTpcMiniMC_Nhits; } // number of hits
   void  GetHit(Int_t i, double &x, double &y, double &z, double &w);
   void  getChi2(Long64_t *NhitSec, Double_t *SecChi2)
   {
      for (Int_t i = 0; i < 24; i++) {
         NhitSec[i] = fNhitSec[i];
         SecChi2[i] = fSecChi2[i];
      }
   }
   void SetPrintLevel(Int_t prt) { printL = prt; }
   void SetEvalCut(Double_t cut) { EvalpHit = cut; }
   //______________________________________________________________________

private:
   // ----- default parameters -----
   // current number of reconstructed events
   Int_t Nreconstructed;
   // Print level
   Int_t printL = 0;
   // Number of events to simulate
   Int_t fNevents = 10;
   // Number of events to simulate
   Int_t fminHits = 20;
   // Magnetic field value (default), T
   Double_t b_mf[3] = {0., 0., 3.0};
   // Electric field direction
   Double_t e_ef[3] = {0., 0., 1.};
   // Track width on the pads plane
   Double_t d_track = 0.03;
   // Minimum pad sensivity threshould, (cm^2)
   Double_t min_padS[2] = {0.0005, 0.0005};
   // Maximum pad sensivity threshould, (cm^2)
   Double_t max_padS[2] = {0.0054, 0.0081};
   // Sigma of the pad signal distribution, (%/100)
   Double_t sig_padS[2] = {0.05, 0.05};
   // Sigma of the time signal distribution, (%/100)
   Double_t sig_timeBins = 0.02;
   // Sigma of the z pozition , (%/100)
   Double_t sig_Z_pos = 0.02;
   // Minimum momentum value, (GeV/c)
   Double_t p_min = 0.2;
   // Maximum momentum value, (GeV/c)
   Double_t p_max = 2.0;
   // Track start point on the beam, (cm)
   Double_t r0_beam[3] = {0, 0, 0};
   // Track momentum value at the start point, (GeV/c)
   Double_t p_las[3] = {1., 0., 1.};
   // calculation mode of thetha&phi angles of  muons
   Int_t k_par = 0; // 0 - const;
                    // 1 - theha "x*dx" phi is uniform (cosmic muon)
                    // 2 - theha "cos{x)" phi is uniform (cosmic muon)
                    // 3 - beam theha&phi are uniform
                    // 4 - TPC laser system rays
   // Chamber Mode
   Int_t chambers = 0;
   // Initial seed for random generator
   UInt_t fSeed = 12345678;
   //  Generator
   TRandom *fRandom = new TRandom();
   // Maximum number of attemps to simulate event
   Int_t maxSim = 1000;
   // Maximum number of pads in one cluster
   Int_t maxPcluster = 3;
   // File name with parameters
   TString *parFile;
   // File name with alignment data
   TString *aFile;
   // output directory name
   TString *outDir;
   // output file name
   TString    *outFile;
   const char *outfile;
   // Sectors LSC shifts due to the alignment
   Double_t R0_A[3][24] = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
   // alpha angle GSC-LSC due to the alignment
   Double_t alpha_A[24];
   // beta angle GSC-LSC due to the alignment
   Double_t beta_A[24];
   // gamma angle GSC-LSC due to the alignment
   Double_t gamma_A[24];
   // Radius of the TPC sensitive volume
   Double_t fr_svol = 133.;
   // Length of the drift volume
   Double_t fL_svol, fZmin, fZmax;

   // Length of the drift volume
   Double_t sZ_min;
   // -----  -----

   Bool_t   beam_point; // Production Mode  (0 - uniform in the space, 1 - beam)
   Bool_t   MFieldMode; // Magnetic field mode  (0 - no field, 1 - field)
   Int_t    sector0;    // first used in the simulation
   Int_t    sector1;    // last used in the simulation
   Int_t    iSector;    // current setor
   TVector3 eGy = TVector3(0, 1, 0);
   TVector3 bMf;
   TVector3 eEf;
   TVector3 r0Beam;
   TVector3 pLas;
   TVector3 fP11, fP12, // start & finish points of the track in the 1st chamber
      fP21, fP22;       // start & finish points of the track in the 2nd chamber

   TVector3 fR1_c, fR2_c; // points on the particle track for the current chamber

   TVector2 fR0sec, // start point of the band center line
      fR0sec1,      // end point of the band center line
      fR1sec,       // start point of the left band line
      fR2sec,       // start point of the right band line
      fR1sec1,      // start point of the left band line
      fR2sec1,      // start point of the right band line
      fEsec,        // band direction (fR0sec1-fR0sec)
      uEsec,        // band direction (fR0sec1-fR0sec)/|fR0sec1-fR0sec|
      Uniper;
   // --- global sector parameters ----
   Double_t dh_track, fWpad, fHpad[2]; // pad width & hight
   Int_t    fNrows[3];
   // "null"  alignment
   TVector3 R0;       // Shift from GSC to LSC in GSC for the current sector
                      //  TRotation l2g;  // Rotation matrix LSC->GSC for the current sector
                      //  TRotation g2l;  // Rotation matrix GSC->LSC for the current sector
   TVector3 e3;       // Unit vector on the LSC Z-axis in GSC
   TVector3 G0;       // Shift from LSC to GSC in LSC for the current sector
                      // current "known" alignment
   TVector3  R0shift; // Shift from GSC to LSC in GSC for the current sector
   TRotation loc2glo; // Rotation matrix LSC->GSC for the current sector
   TRotation glo2loc; // Rotation matrix GSC->LSC for the current sector
   TVector3  e_3;     // Unit vector on the LSC Z-axis in GSC
   TVector3  G0shift; // Shift from LSC to GSC in LSC for the current sector

   /*   Int_t fNhits;
      Double_t fXHits[1000],fYHits[1000],fZHits[1000] // clusters coordinates
              ,fWHits[1000];                                // clusters weights
   */
   // --- cluster room ----
   Double_t xPad[10], yPad[10], zPad[10], sPad[10];

   Double_t trackPin[6], trackPout[6], trackDiff[6];
   // initial & reconstructed track parameters
   //  [0],[1],[2] - vertex coordinates
   //  [3],[4] - theta & phi angles (line)
   //  [3] - radius with the sign (helix)
   //  [4] - alpha_0 starting angle (helix)
   //  [5] -lambda - angle track-MF (helix)
   Double_t EvalpHit = 2;
   Double_t Chi2Fit;

   // ----- laser system variables  -----
   Bool_t      ls1, ls2; // true if to get parameters
   Int_t       i_ls1, j_ls1, a_ls1, i_ls2, j_ls2, a_ls2;
   Int_t       j1_ls = -4, j2_ls = 4;
   const Int_t i1_ls = 0, i2_ls = 4, a11_ls = -1, a12_ls = 1, a21_ls = -2, a22_ls = 2;
   // -----  -----

   TFile                    *alignmentF;
   TFile                    *outmMCfile;
   TTree                    *t_alignment;
   TTree                    *t_parameters;
   TTree                    *t_data;
   TTree                    *tout_fullChi2;
   TpcSectorGeoAlignmentVAK *fTpcSecGeo;
   // -----  -----
   Int_t miniEvent();              // generate an event
   Int_t laser_ray();              // laser ray simulation
   Int_t lineClusters(Int_t isec); // find clusters of direct line track
   Int_t ClustersInRow(Int_t iRow, Int_t iPad1, Int_t iPad2, Int_t iRoc,
                       Double_t *sPads);                // find clusters in direct line track
                                                        //  Int_t cosmic_noMF(); // cosmic ray simulation without MF
   Double_t theta_cosm();                               // simulate theta angle for cosmic muon
   Double_t theta_beam();                               // simulate theta angle for produced muon
   Double_t RandomPadSignal(Int_t iRoc, Double_t sPad); // simulate pad signal
   Bool_t   newCluster(Double_t hpad, Int_t nPads);     // create new cluster
   Double_t simLineZpos(TVector2 LCxy);                 // simulate z pozition
   Double_t simHelixZpos(TVector2 LCxy);                // simulate z pozition
   Double_t f1(Double_t x, Double_t y);                 // left track bound line function
   Double_t f2(Double_t x, Double_t y);                 // right track bound line function
   Double_t padANDband(Double_t s1, Double_t s2, Double_t b1, Double_t b2); // pad
                                                                            // square covered by the track band
   Double_t theta_a[24] = {24 * 0.};
   Double_t beta_a[24]  = {24 * 0.};
   Double_t alpha_a[24] = {24 * 0.};
   Double_t sTriangle(Int_t line, Double_t s, Double_t b);
   Double_t sTrapezoid(Int_t line, Int_t axes, Double_t s, Double_t b1, Double_t b2);
   Double_t sParallelogram(Int_t axes, Double_t s1, Double_t s2, Double_t b1, Double_t b2);
   // full hits number
   Long64_t fNhits = 0;
   // full chi2 sum of all hits
   Double_t fSumChi2 = 0;
   // sector hits number
   Long64_t fNhitSec[24];
   // full chi2 sum of all hits
   Double_t fSecChi2[24];

   // -----  -----          constants             --------------
   const Double_t pi        = TMath::Pi();    // sector angle
   const Double_t twopi     = TMath::TwoPi(); // sector angle
   const Double_t r2d       = 180 / pi;
   const Double_t d2r       = pi / 180;
   const Double_t phi_shift = TMath::Pi() / 12;
   const Double_t phi_sec   = TMath::Pi() / 6; // sector angle
   const Double_t epsMin    = 1e-12;           // minimal number
   const Double_t epsMax    = 1e-5;            // non zero number

   // -----  ----- magnetic field case parameters --------------
   Int_t q_muon, // particle charge
      q_sign,    // charge/|charge|
      signRotB,  // rotation in BCS: -1(counterclockwise)
                 //            1(clockwise)
      signRotL   // rotation in LCS: -1(counterclockwise)
                 //            1(clockwise)
      ;
   Double_t minMuonTheta, // minimum muon angle with Z-axis
      xcB, ycB,           // coordinates of the cilinder center in BCS
      xcL, ycL,           // coordinates of the ellipse  center in LCS
      rH, rHa[3],         // radiuses of the helixes
      phi_0B,             // start angle of the muon to X-axis in BCS
      phi_E,              // rotation angle in LCS to ECS
      phi_0E,             // start angle of the muon to X-axis in ECS
      Ae, Be,             // large and small semi-axes of the ellipse
      pTMuonB,            // muon momentum in BCS
      helix_Pz,           // Pz for the helix formula
      e11, e12, e21, e22,
      bMF_mod, // module of magnetic field vector
      t_end    // helix patameter t at the chamber crossing
      ;
   TVector3 r0BeamB,      // production point in BCS
      r0B2G,              // shift GCS in BGS
      eBx, eBy, eBz,      // unit coordinates vectors in BCS
      rc0_G,              // coordinates of the helix center in GCS
      rc0_B,              // coordinates of the helix center in BCS (magnetic field)
      rc0l_G,             // coordinates of the ellipse center in GCS
      rc0_L,              // coordinates of the ellipse center in LCS (sector)
      bMf_L,              // magnetic field vector in LCS (sector)
      rGt, rBt, rLt, rLt0 // helix point coordinates in GCS,BCS,LCS
      ;
   TRotation B2G, // Rotation matrix BSC->GSC for the current sector
      G2B,        // Rotation matrix GSC->BSC for the current sector
      L2E         // Rotation matrix BSC->GSC for the current sector
                  //            E2L;   // Rotation matrix GSC->BSC for the current sector
      ;
   TVector3 pMuonB;
   // -----  ----- methods --------------  TVector3 pLas;

   Bool_t        crossCheck(Int_t iSec, TVector3 r0, TVector3 r1);
   void          crossLinePoints();
   void          GetCosmicMuon();
   void          GetLaserRay();
   void          setMuon(Double_t x, Double_t y, Double_t z, Double_t theta, Double_t phi);
   void          GetBeamMuon();
   void          GetPointMuon();
   Int_t         passMuoninMF();
   void          CheckPadsinMF();
   Int_t         coordinates(Double_t t, Double_t R);
   inline double PhiMinusPhi0(double phi0, double phi, Int_t q)
   {
      if (q * (phi0 - phi) > 0)
         return abs(phi0 - phi);
      else
         return twopi - abs(phi0 - phi);
   }
   // -----  ----- fits --------------
   Int_t fitLine();
   Int_t fitHelix();
   void  putChi2L();
   void  putChi2H();
   Int_t fitCenter(Double_t x0, Double_t y0, Double_t &xc, Double_t &yc);
   ClassDef(TpcMiniMC, 1);
};

//__________________________________________________________________________
inline Double_t TpcMiniMC::f1(Double_t x, Double_t y)
{
   // left track bound line function
   return fEsec.Y() * (x - fR1sec.X()) - fEsec.X() * (y - fR1sec.Y());
}

//__________________________________________________________________________
inline Double_t TpcMiniMC::f2(Double_t x, Double_t y)
{
   // right track bound line function
   return fEsec.Y() * (x - fR2sec.X()) - fEsec.X() * (y - fR2sec.Y());
}

#endif
