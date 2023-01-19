void missMC(const char *aReal = "0r5_0g5.root", const char *aUsed = "0r0_0g0.root", const char *outDir = "miniDST",
            const char *outFile = "lrc0_0r5_0g5_cl0_dt08evt10000")
{

   Int_t    Chambers    = 1;
   Double_t B_mf[3]     = {0, 0, 0}; // magnetic field, Tesla
   Int_t    MaxPcluster = 333;       // maximum pad number in cluster
   Int_t    minHits     = 20;        // minimum hits number in an event
   Int_t    events      = 10000;     // number of events to be simulated

   Bool_t   Beam_point; // single muon production point
   Double_t D_track;    // track width at the ROC level
   Double_t K_par;      // cosmic or laser rays

   Beam_point = true;
   D_track    = 0.8;
   K_par      = 3; // laser rays
                   //  Beam_point=false;D_track=0.8;K_par=2;                        // cosmic rays

   // file name with the alignment to be simulated
   TString *aRealIn = new TString(aReal);
   // file name with the alignment used at reconstruction
   TString *aUsedIn = new TString(aUsed);
   // directory to load results
   TString *OutDir = new TString(outDir);
   // output root file (the name without ".root" extension
   TString *OutFile = new TString(outFile);

   Double_t E_ef[3]      = {0, 0, 3};        // electrical field direction (arbitrary units)
   Double_t R0_beam[3]   = {0., 0., 0.};     // muon production point coordinates
   Double_t P_las[3]     = {1., 0., 1.};     // muon momentum direction (single production)
   Double_t P_min        = 1.;               // minimum muon momentum, GeV
   Double_t P_max        = 1.;               // maximum muon momentum, GeV
   Double_t Sig_timeBins = 0.02;             // not used sigma of the time digitizer
   Double_t Sig_padS[2]  = {0.05, 0.05};     // pad signal sigma, %
   Double_t Min_padS[2]  = {0.0005, 0.0005}; // maximum pad sgnal, cm^2
   Double_t Max_padS[2]  = {0.0054, 0.0081}; // minimum pad sgnal, cm^2
   Double_t Sig_Z_pos    = 0.01;             // z-coordinate (gaz velocity) sigma, %
   UInt_t   seed         = 12345678;         // initial seed number for random generator
   // Null alignment
   BaseTpcSectorGeo *fTpcSecGeo = new TpcSectorGeoAlignmentVAK();
   TpcMissAlignment *mc = new TpcMissAlignment(*fTpcSecGeo, B_mf, Beam_point, E_ef, P_las, D_track, Min_padS, Max_padS,
                                               Sig_padS, Sig_timeBins, Sig_Z_pos, R0_beam, K_par, P_min, P_max,
                                               Chambers, MaxPcluster, minHits, seed, aRealIn, aUsedIn, OutDir, OutFile);

   Long64_t *NhitSec = new Long64_t[24];
   Double_t *SecChi2 = new Double_t[24];

   mc->SetPrintLevel(1);          // print informatin per number of events
   mc->SetEvalCut(4);             // set maximum Xi^2 value at the event reconstruction
   mc->simulate(events);          // simulate the "events" number
   mc->getChi2(NhitSec, SecChi2); // get sector hits numbers and Xi^2

   Double_t sum = 0;
   Long64_t nhs = 0;
   for (Int_t is = 0; is < 24; is++)
      if (NhitSec[is] > 0) {
         sum += SecChi2[is];
         nhs += NhitSec[is];
         Double_t ws = (NhitSec[is] > 0) ? SecChi2[is] / NhitSec[is] : 0;
         cout << "sector=" << is << "  chi2/ndf=" << ws << " nhits=" << NhitSec[is] << endl;
      }
   sum = (nhs > 0) ? sum / nhs : 0;
   cout << "sum for chi2/ndf=" << sum << "  at nhits=" << nhs << endl;
}
