// -------------------------------------------------------------------------
// -----                MpdDCMSMMGenerator source file                 -----
// -----                Created 27-AUG-2019  by Igor Rufanov           -----
// -------------------------------------------------------------------------
#include "MpdDCMSMMGenerator.h"

#include "FairPrimaryGenerator.h"
#include "FairMCEventHeader.h"

#include "TRandom.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "FairIon.h"
#include "FairParticle.h"
#include "FairRunSim.h"

using namespace std;
using namespace TMath;

// -----   Default constructor   ------------------------------------------

MpdDCMSMMGenerator::MpdDCMSMMGenerator() : FairGenerator(), fInputFile(NULL), fFileName(NULL) {}
// ------------------------------------------------------------------------

// -----   Standard constructor   -----------------------------------------

MpdDCMSMMGenerator::MpdDCMSMMGenerator(const char *fileName)
   : FairGenerator(), fInputFile(NULL), fFileName(fileName), fSpectatorsON(kFALSE), fFixedTarget(kFALSE), fGammaCM(1.),
     fBetaCM(0.)
{
   cout << "-I MpdDCMSMMGenerator: Opening input file " << fFileName << endl;
#ifdef GZIP_SUPPORT
   fInputFile = gzopen(fFileName, "rb");
#else
   fInputFile = fopen(fFileName, "r");
#endif
   if (!fInputFile) {
      Fatal("MpdDCMSMMGenerator", "Cannot open input file.");
      exit(1);
   }
   // read file header
   char     read[200];
   Int_t    A1, Z1, A2, Z2;
   Double_t T0, sqS;
   for (Int_t i = 0; i < 3; i++) {
#ifdef GZIP_SUPPORT
      char *ch = gzgets(fInputFile, read, 200);
#else
      char *ch = fgets(read, 200, fInputFile);
#endif
      cout << "-I MpdDCMSMMGenerator:" << read;
      if (i == 0) {
         string str0(read);
         A1 = stoi(str0.substr(str0.find("A1=") + 3, 3));
         Z1 = stoi(str0.substr(str0.find("Z1=") + 3, 3));
         A2 = stoi(str0.substr(str0.find("A2=") + 3, 3));
         Z2 = stoi(str0.substr(str0.find("Z2=") + 3, 3));
         // Int_t A1= 0;
      } else if (i == 1) {
         string str0(read);
         T0  = stof(str0.substr(str0.find("T0=") + 3, 8));
         sqS = stof(str0.substr(str0.find("sqrt(s)=") + 8, 8));
         // Double_t mp=0.938272;
         Double_t mp   = 0.940; // to obtain equivalence of read "sqrt(s)=" and calculated below mCMS
         Double_t e    = T0 + mp;
         Double_t p    = pow(e * e - mp * mp, 0.5);
         Double_t mCMS = pow(2. * mp * (e + mp), 0.5);
         cout << "mCMS=" << mCMS << endl;
         fGammaCM = (e + mp) / mCMS;
         fBetaCM  = pow(1. - 1. / (fGammaCM * fGammaCM), 0.5);
      }
   }
   cout << "-I MpdDCMSMMGenerator: A1=" << A1 << " Z1=" << Z1 << " A2=" << A2 << " Z2=" << Z2 << " T0=" << T0
        << " sqS=" << sqS << endl;

   if (fSpectatorsON) {
      Int_t n = RegisterIons();
   }
}
// ------------------------------------------------------------------------

// -----   Destructor   ---------------------------------------------------

MpdDCMSMMGenerator::~MpdDCMSMMGenerator()
{
   if (fInputFile) {
#ifdef GZIP_SUPPORT
      gzclose(fInputFile);
#else
      fclose(fInputFile);
#endif
      fInputFile = NULL;
   }
   cout << "Leave Destructor of MpdDCMSMMGenerator" << endl;
}
// ------------------------------------------------------------------------

// -----   Public method ReadEvent   --------------------------------------

Bool_t MpdDCMSMMGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{

   // ---> Check for input file
   // cout<<"MpdDCMSMMGenerator::ReadEvent -------------------------"<<endl;
   if (!fInputFile) {
      cout << "-E MpdDCMSMMGenerator: Input file not open! " << endl;
      return kFALSE;
   }

   // ---> Check for primary generator
   if (!primGen) {
      cout << "-E- MpdDCMSMMGenerator::ReadEvent: "
           << "No PrimaryGenerator!" << endl;
      return kFALSE;
   }

   char read[80];
#ifdef GZIP_SUPPORT
   char *ch = gzgets(fInputFile, read, 80);
#else
   char *ch = fgets(read, 80, fInputFile);
#endif
   Int_t   evnr = 0;
   Float_t b, bimpX, bimpY;
   sscanf(read, "%d %f %f %f", &evnr, &b, &bimpX, &bimpY);

#ifdef GZIP_SUPPORT
   if (gzeof(fInputFile)) {
#else
   if (feof(fInputFile)) {
#endif
      cout << "-I MpdDCMSMMGenerator : End of input file reached." << endl;
      const Bool_t ZeroSizeEvents = kFALSE;
      if (ZeroSizeEvents) {
         return kTRUE; // gives zero multiplicity events after EOF
      } else {         // gives geant crash after EOF and one empty event in .root file
#ifdef GZIP_SUPPORT
         gzclose(fInputFile);
#else
         fclose(fInputFile);
#endif
         fInputFile = NULL;
         return kFALSE;
      }
   }

   Float_t phi = atan2(bimpY, bimpX);

   // Set event id and impact parameter in MCEvent if not yet done
   FairMCEventHeader *event = primGen->GetEvent();
   if (event && (!event->IsSet())) {
      event->SetEventID(evnr);
      event->SetB(b);
      event->MarkSet(kTRUE);
      event->SetRotZ(phi);
   }

   TDatabasePDG *dbPDG = TDatabasePDG::Instance();

   Float_t px, py, pz;
   for (Int_t ibeam = 0; ibeam < 3; ibeam++) { // spectators pz+, spectators pz-, participants.
      Int_t np = 0;
#ifdef GZIP_SUPPORT
      ch = gzgets(fInputFile, read, 80);
#else
      ch = fgets(read, 80, fInputFile);
#endif
      sscanf(read, "%d", &np); // cout<<np<<" "<<endl;

      for (Int_t i = 0; i < np; i++) {
#ifdef GZIP_SUPPORT
         ch = gzgets(fInputFile, read, 80);
#else
         ch = fgets(read, 80, fInputFile);
#endif
         Int_t   iN, iQ, iS = 0;
         Float_t xxx = 0., mass;
         if (ibeam < 2)
            sscanf(read, "%d %d %f %f %f %f", &iN, &iQ, &xxx, &px, &py, &pz);
         else
            sscanf(read, "%d %d %d %f %f %f %f", &iN, &iQ, &iS, &px, &py, &pz, &mass);

         Int_t pdgID = 0;
         if (ibeam == 2) { // participants
            Double_t massFactor = 1.;
            pdgID               = FindPDGCodeParticipant(iN, iS, iQ, mass, massFactor);
            if (massFactor != 1.) {
               px *= massFactor;
               py *= massFactor;
               pz *= massFactor;
            }
            if (pdgID > 1000030000) cout << pdgID << " " << iN << " " << iS << " " << iQ << " " << mass << endl;
         } else { // spectators
            if (fSpectatorsON) {
               Int_t dN = -999; // difference of number of nucleons iN-NTab between DCM-DCMSMM and registered ions
               pdgID    = FindPDGCodeSpectator(iN, iQ, dN);
            }
         }

         if (pdgID) {
            if (fFixedTarget) {
               TParticlePDG *particle = dbPDG->GetParticle(pdgID);
               Double_t      m        = particle->Mass();
               Double_t      e        = pow(px * px + py * py + pz * pz + m * m, 0.5);
               Double_t      pzF      = fGammaCM * (pz + fBetaCM * e); // Lorentz transformation to fixed target system
               Double_t      eF       = pow(px * px + py * py + pzF * pzF + m * m, 0.5);
               // cout<<fGammaCM<<" "<<fBetaCM<<" "<<pdgID<<" "<<m<<" "<<pz<<" "<<pzF<<" "<<eF-m<<endl;
               pz = pzF;
            }
            // Int_t Geant3ID= dbPDG->ConvertPdgToGeant3(pdgID);
            primGen->AddTrack(pdgID, px, py, pz, 0., 0., 0.);
         } else {
            if (ibeam == 2 || (ibeam == 2 && fSpectatorsON))
               cout << "-I MpdDCMSMMGenerator : unknown particle N=" << iN << " Q=" << iQ << endl;
         }
      }
   }

   return kTRUE;
}

// ------------------------------------------------------------------------

Bool_t MpdDCMSMMGenerator::SkipEvents(Int_t count)
{
   if (count <= 0) return kTRUE;
   for (Int_t ii = 0; ii < count; ii++) {
      // ---> Check for input file
      if (!fInputFile) {
         cout << "-E MpdDCMSMMGenerator: Input file not open! " << endl;
         return kFALSE;
      }
      char read[80];
#ifdef GZIP_SUPPORT
      char *ch = gzgets(fInputFile, read, 80);
#else
      char *ch = fgets(read, 80, fInputFile);
#endif
      Int_t   evnr = 0;
      Float_t b, bimpX, bimpY;
      sscanf(read, "%d %f %f %f", &evnr, &b, &bimpX, &bimpY);
      for (Int_t ibeam = 0; ibeam < 3; ibeam++) {
         Int_t np = 0;
#ifdef GZIP_SUPPORT
         ch = gzgets(fInputFile, read, 80);
#else
         ch = fgets(read, 80, fInputFile);
#endif
         sscanf(read, "%d", &np);
         for (Int_t i = 0; i < np; i++) {
#ifdef GZIP_SUPPORT
            ch = gzgets(fInputFile, read, 80);
#else
            ch = fgets(read, 80, fInputFile);
#endif
            Int_t   iN, iQ, iS = 0;
            Float_t xxx = 0., mass, px, py, pz;
            if (ibeam < 2)
               sscanf(read, "%d %d %f %f %f %f", &iN, &iQ, &xxx, &px, &py, &pz); // cout<<np<<" "<<endl;
            else
               sscanf(read, "%d %d %d %f %f %f %f", &iN, &iQ, &iS, &px, &py, &pz, &mass);
         }
      }
   }
   return kTRUE;
}

// ------------------------------------------------------------------------

Int_t MpdDCMSMMGenerator::FindPDGCodeParticipant(Int_t A, Int_t S, Int_t Z, Float_t mass, Double_t &massFactor)
{
   Int_t k  = 0;
   Int_t aA = abs(A), aS = abs(S), aZ = abs(Z);

   if (aA == 0) {                        // mesons and gamma =====================================
      if (S == 0) {                      // not strange
         if (aZ == 0) {                  // neutral
            const Int_t n000        = 7; // gamma, pi0,    eta,   rho, omega, etaPri,  phi
            Float_t     M000M[n000] = {0.0, 0.13497, 0.54745, 0.7690, 0.78265, 0.95778, 1.019455};
            Float_t     M000L[n000] = {0.0, 0.13400, 0.54700, 0.7680, 0.78200, 0.95700, 1.019390};
            Float_t     M000U[n000] = {0.0, 0.13510, 0.54810, 0.7711, 0.78310, 0.95810, 1.019610};
            Int_t       C000[n000]  = {22, 111, 211, 113, 223, 331, 333};
            for (Int_t i = 0; i < n000; i++) {
               if (mass < M000L[i] || mass > M000U[i]) continue;
               k = C000[i];
               break;
            }

            if (k == 0 && mass > 0.4979 && mass < 0.4981) k = 22; // rare particle in pairs in begining of list
            // invariant mass of these pairs (considering them as gamma) is close to mass of pion
         } else { // non-stranges charged mesons
            if (mass > 0.1389 && mass < 0.1401)
               k = 211 * Z; // pi+/-
            else if (mass > 0.134 && mass < 0.136)
               k = 211 * Z; // also pi+/-
            else if (mass == 0.776)
               k = 213 * Z; // rho+/-
            else if (mass > 0.0004 && mass < 0.0011)
               k = -11 * Z; // e+/e-
            else if (mass > 0.1055 && mass < 0.1065)
               k = -13 * Z; // mu+/mu-
         }
      } else {                                   // A==0 && S!=0 strange mesons
         if (aZ == 0) {                          // neutral K0, K*0
            if (mass > 0.491 && mass < 0.4981) { // K0 S=1 311, AK0 S=-1 -311
               Double_t KSL = gRandom->Uniform(0., 1.);
               if (KSL < 0.5)
                  k = 130; // KL
               else
                  k = 310; // KS
            } else if (mass == 0.8922)
               k = 313 * S; // K*0
         }
         if (aZ == 1) {
            if (mass > 0.491 && mass < 0.4951)
               k = 321 * Z; // K+
            else if (mass > 0.8815 && mass > 0.8895)
               k = 323 * Z; // K*+ (not occured)
         }
      }
   } else if (aA == 1) { // barions ============================================
      if (S == 0) {      // Light baryons: n,p,Delta
         if (Z == 0) {
            if ((mass > 0.9389 && mass < 0.940) || mass < -1.)
               k = 2112 * A; // n
            else if (mass == 1.232)
               k = 2114 * A; // Delta0 (not observed)
         } else if (aZ == 1) {
            if ((mass > 0.9379 && mass < 0.9401) || mass < -1.)
               k = 2212 * A; // p
            else if (mass == 1.232) {
               if (A * Z == 1)
                  k = 2214 * A; // Delta+ (not observed)
               else if (A * Z == -1)
                  k = 1114 * A; // Delta- (not observed)
            }
         } else if (aZ == 2) {
            if (mass == 1.232) k = 2224 * A; // Delta++ (not observed)
         }
      } else if (aS == 1) { // Lambda, Sigma, Sigma*
         if (Z == 0) {
            if (mass > 1.1149 && mass < 1.1161)
               k = 3122 * A; // Lambda
            else if (mass > 1.189 && mass < 1.1921)
               k = 3212 * A; // Sigma0
            else if (mass == 1.382)
               k = 3214 * A;   // Sigma*0  (not observed)
         } else if (aZ == 1) { // Sigma-, Sigma+, Sigma*-, Sigma*+
            if (mass > 1.1889 && mass < 1.1921)
               k = 3222 * A; // Sigma+
            else if (mass > 1.196 && mass < 1.1971)
               k = 3112 * A; // Sigma-
            else if (mass == 1.3823)
               k = 3224 * A; // Sigma*+ (not observed)
            else if (mass == 1.3875)
               k = 3114 * A; // Sigma*- (not observed)
         }
      } else if (aS == 2) { // Xi, Xi*
         if (Z == 0) {      // Xi0, Xi*0
            if (mass > 1.3139 && mass < 1.3161)
               k = 3322 * A; // Xi0
            else if (mass == 1.5318)
               k = 3324 * A;   // Xi*0
         } else if (aZ == 1) { // Xi-, Xi*-
            if (mass > 1.3139 && mass < 1.3221)
               k = 3312 * A; // Xi-
            else if (mass == 1.535)
               k = 3314 * A; // Xi*-
         }
      } else if (aS == 3) { // Omega-
         if (aZ == 1) {
            if (mass > 1.6 && mass < 1.801) k = 3334 * A; // Omega- 1.67245
         }
      }
   } else if (A == 2) { // d is only valid ion (particle). All over are unstable. ==

      if (S == 0) {
         if (Z == 1 && (mass < -1 || (mass > 1.875 && mass < 1.877))) k = 1000010020; // d
      } else if (S == -1) { // unstable strange pairs => take the most heavy nucleon
         if (Z == -1) {
            k          = 3112;
            massFactor = 1.197436 / mass;
         } // Sigma- n 2.13699 GeV (not occured)
         else if (Z == 0 && mass > 2.0545 && mass < 2.0555) {
            k          = 3122;
            massFactor = 1.11568 / mass;
         } // Lambda n
         else if (Z == 1 && mass > 2.0535 && mass < 2.0545) {
            k          = 3122;
            massFactor = 1.11568 / mass;
         } // Lambda p
         else if (Z == 2) {
            k          = 3222;
            massFactor = 1.18937 / mass;
         } // Sigma+ p 2.12764 GeV (not occured)
      } else if (S == -2) {
         if (Z == -1 && mass > 2.2545 && mass < 2.2615) {
            k          = 3312;
            massFactor = 1.32132 / mass;
         } // Xi- n M=2.26088 ( also ghost at 2.255 GeV)
         else if (Z == 0) {
            if (mass > 2.2305 && mass < 2.2315) {
               k          = 3122;
               massFactor = 1.11568 / mass;
            } // Lambda Lambda M=2.23136
            else if (mass > 2.2525 && mass < 2.2555) {
               k          = 3322;
               massFactor = 1.19255 / mass;
            } // Xi0 n M=2.25446
            else if (mass > 2.2595 && mass < 2.2605) {
               k          = 3312;
               massFactor = 1.197436 / mass;
            } // Xi- p M=2.25959
         } else if (Z == 1 && mass > 2.2525 && mass < 2.2535) {
            k          = 3322;
            massFactor = 1.3149 / mass;
         } // Xi0 p M=2.25317 (132)
      } else if (S == -3) {
         if (Z == -1 && mass > 2.4305 && mass < 2.4375) {
            k          = 3312;
            massFactor = 1.197436 / mass;
         } // Xi- Lambda M=2.43700 ( also ghost at 2.431 GeV) (37)
         else if (Z == 0 && mass > 2.4301 && mass < 2.4315) {
            k          = 3322;
            massFactor = 1.3149 / mass;
         } // Xi0 Lambda M=2.43058 (36)
      } else if (S == -4 && Z == -1 && mass > 2.6295 && mass < 2.6365) {
         k          = 3312;
         massFactor = 1.32132 / mass;
      } // Xi- Xi0 M=2.63622 (1)
   } else if (A == 3) {
      //
      if (S == 0) {
         if (Z == 1 && mass > 2.8085 && mass < 2.8095)
            k = 1000010030; // t
         else if (Z == 2 && mass > 2.8085 && mass < 2.8095)
            k = 1000020030; // He3
      } else if (S == -1) {
         if (Z == 1 && mass > 2.9915 && mass < 2.9925) k = 1010010030; // H3L
      }
   } else if (A == 4) {
      if (S == 0) {
         if (Z == 2 && mass > 3.7275 && mass < 3.7285) k = 1000020040; // Alpha
         // else if( Z==1 && mass>3.7275 && mass<3.7285) k=1000020040; // H4
      } else if (S == -1) {
         if (Z == 1 && mass > 3.9245 && mass < 3.9255) k = 1010010040; // H4L
         if (Z == 2 && mass > 3.9245 && mass < 3.9255) k = 1010020040; // He4L
      } else if (S == -2) {
         if (Z == 1 && mass > 4.1065 && mass < 4.1075) k = 1020010040; // H4LL (6 per 100k MPD events)
      }
   } else if (A == 5) {
      if (S == -2 && mass > 5.0405 && mass < 5.0415) k = 1020020050; // He5LL (0 per 100k MPD events)
   }

   // if( aA<=2 && k==0) {
   if (mass > -10 && k == 0) {
      // if( A==0 && S==0 && Z==0 && fabs( mass)-0.498 < 0.001) {}
      // else
      cout << "A=" << A << " S=" << S << " Z=" << Z << "  mass= " << mass << "   Pdg=" << k << endl;
   }
   // k= 1000030060;
   return k;
}

// ------------------------------------------------------------------------

Int_t MpdDCMSMMGenerator::FindPDGCodeSpectator(Int_t N, Int_t Z, Int_t &dN)
{
   Int_t k = 0;
   dN      = 0;
   if (Z == 0 && N == 1)
      k = 2112; // n
   else if (Z == 1 && N == 1)
      k = 2212; // p
   else if (Z >= 1 && Z <= fZMax) {
      Int_t NTab;
      if (N < fN1[Z]) {
         NTab = fN1[Z];
         dN   = N - fN1[Z];
      } else if (N <= fN2[Z]) {
         NTab = N;
         dN   = 0;
      } else {
         NTab = fN2[Z];
         dN   = N - fN2[Z];
      }
      k = 1000000000 + Z * 10000 + NTab * 10;
   }
   return k;
}

// ------------------------------------------------------------------------

Int_t MpdDCMSMMGenerator::RegisterIons(void)
{
   // see ./fairsoft/transport/geant3/TGeant3/TGeant3.cxx::AddParticlesToPdgDataBase()
   // or ./fairsoft/transport/geant4_vmc/source/physics/src/TG4ParticlesManager::DefineParticles()
   // Deuteron, Triton, Alpha, HE3, HyperTriton and their ANTI-particle are added there
   // also Special particles: "Cherenkov", "FeedbackPhoton".
   // He4L, H4L are added in mpdroot/gconfig/UserDecay.C

   struct gpions {
      Int_t    N, Z;
      Char_t   Name[6];
      Double_t Mass;
   };
   const Int_t         NSI       = 96;
   const struct gpions ions[NSI] = {{2, 1, "d", 0.0},
                                    {3, 1, "t", 0.0},
                                    {3, 2, "He3", 0.0},
                                    {4, 2, "He4", 0.0},
                                    {5, 2, "He5q", 4.6688534},
                                    {6, 2, "He6", 5.6065596},
                                    {6, 3, "Li6", 5.60305},
                                    {7, 3, "Li7", 6.53536},
                                    {8, 3, "Li8", 7.4728996},
                                    {9, 3, "Li9", 8.40840118},
                                    {7, 4, "Be7", 6.53622},
                                    {8, 4, "Be8q", 7.4568945},
                                    {9, 4, "Be9", 8.39479},
                                    {10, 5, "B10", 9.32699},
                                    {11, 5, "B11", 10.25510},
                                    {10, 6, "C10", 9.3306397},
                                    {11, 6, "C11", 10.257085},
                                    {12, 6, "C12", 11.17793},
                                    {13, 6, "C13", 12.112548},
                                    {14, 6, "C14", 13.043937},
                                    {14, 7, "N14", 13.04378},
                                    {16, 8, "O16", 14.89917},
                                    {19, 9, "F19", 17.69690},
                                    {20, 10, "Ne20", 18.62284},
                                    {23, 11, "Na23", 21.41483},
                                    {24, 12, "Mg24", 22.34193},
                                    {27, 13, "Al27", 25.13314},
                                    {28, 14, "Si28", 26.06034},
                                    {31, 15, "P31", 28.85188},
                                    {32, 16, "S32", 29.78180},
                                    {35, 17, "Cl35", 32.57328},
                                    {36, 18, "Ar36", 33.50356},
                                    {39, 19, "K39", 36.29447},
                                    {40, 20, "Ca40", 37.22492},
                                    {45, 21, "Sc45", 41.87617},
                                    {48, 22, "Ti48", 44.66324},
                                    {51, 23, "V51", 47.45401},
                                    {52, 24, "Cr52", 48.38228},
                                    {55, 25, "Mn55", 51.17447},
                                    {56, 26, "Fe56", 52.10307},
                                    {59, 27, "Co59", 54.89593},
                                    {58, 28, "Ni58", 53.96644},
                                    {63, 29, "Cu63", 58.61856},
                                    {64, 30, "Zn64", 59.54963},
                                    {69, 31, "Ga69", 64.2037653},
                                    {74, 32, "Ge74", 68.85715},
                                    {75, 33, "As75", 69.78902528},
                                    {80, 34, "Se80", 74.44178},
                                    {81, 35, "Br81", 75.373047},
                                    {84, 36, "Kr84", 78.16309},
                                    {85, 37, "Rb85", 79.09483},
                                    {88, 38, "Sr88", 81.88358},
                                    {89, 39, "Y89", 82.81527},
                                    {90, 40, "Zr90", 83.74571},
                                    {93, 41, "Nb93", 86.541743},
                                    {98, 42, "Mo98", 91.19832},
                                    {99, 43, "Tc99", 92.13059},
                                    {100, 44, "Ru100", 93.060191},
                                    {103, 45, "Rh103", 93.9934966},
                                    {106, 46, "Pd106", 98.64997},
                                    {109, 47, "Ag109", 101.444134},
                                    {114, 48, "Cd114", 106.10997},
                                    {117, 49, "In117", 108.89586},
                                    {120, 50, "Sn120", 111.68821},
                                    {123, 51, "Sb123", 114.48455},
                                    {125, 52, "Te125", 116.3477},
                                    {127, 53, "I127", 118.210768},
                                    {132, 54, "Xe132", 122.86796},
                                    {135, 55, "Cs135", 125.66412195},
                                    {138, 56, "Ba138", 128.45793},
                                    {139, 57, "La139", 129.3904488},
                                    {140, 58, "Ce140", 130.32111},
                                    {141, 59, "Pr141", 131.2546475},
                                    {145, 60, "Nd145", 134.9852},
                                    {149, 61, "Pm149", 138.716549},
                                    {152, 62, "Sm152", 141.51236},
                                    {153, 63, "Eu153", 142.44522},
                                    {157, 64, "Gd157", 146.17374},
                                    {161, 65, "Tb161", 149.90308},
                                    {164, 66, "Dy164", 152.69909},
                                    {165, 67, "Ho165", 153.63162},
                                    {168, 68, "Er168", 156.42801},
                                    {171, 69, "Tm171", 159.2262758},
                                    {174, 70, "Yb174", 162.02245},
                                    {175, 71, "Lu175", 162.956297},
                                    {178, 72, "Hf178", 165.7535059},
                                    {181, 73, "Ta181", 168.55199},
                                    {184, 74, "W184", 171.34924},
                                    {185, 75, "Re185", 172.28258},
                                    {188, 76, "Os188", 175.079754},
                                    {191, 77, "Ir191", 177.878667},
                                    {194, 78, "Pt194", 180.67513},
                                    {197, 79, "Au197", 183.47324},
                                    {202, 80, "Hg202", 188.13451},
                                    {205, 81, "Tl205", 190.93247},
                                    {208, 82, "Pb208", 193.72907}};

   TDatabasePDG *dbPDG = TDatabasePDG::Instance();
   for (Int_t i = 0; i < NSI; i++) {
      if (ions[i].Mass < 0.1) continue; // expected that it is already registered
      Int_t ID = 1000000000 + ions[i].Z * 10000 + ions[i].N * 10;
      if (dbPDG->GetParticle(ID)) continue;
      FairIon *ion = new FairIon(ions[i].Name, ions[i].Z, ions[i].N, ions[i].Z);
      FairRunSim::Instance()->AddNewIon(ion);
   }
   for (Int_t i = 0; i < NSI; i++) {
      Int_t z = ions[i].Z;
      if (z < 0 && z > fZMax) continue;
      fN2[z] = ions[i].N;
   }
   for (Int_t i = NSI - 1; i >= 0; i--) {
      Int_t z = ions[i].Z;
      if (z < 0 && z > fZMax) continue;
      fN1[z] = ions[i].N;
   }
   return NSI;
}
ClassImp(MpdDCMSMMGenerator);
