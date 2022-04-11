// -------------------------------------------------------------------------
// -----                    MpdDecayer source file                     -----
// -----                      Created 22/04/2020                       -----
// -----                  R. Akhat,  A. Zinchenko                      -----
// -----                  External decayer for MPD                     -----
// -------------------------------------------------------------------------

#include "MpdDecayer.h"

#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include <TDecayChannel.h>
#include "TF1.h"
//#include "TGeant3.h"
#include "TH1.h"
#include <THashList.h>
#include "TLorentzVector.h"
#include "TParticle.h"
#include <TParticlePDG.h>
#include "TPDGCode.h"
#include "TPythia6.h"
#include "TVirtualMC.h"

//#include "G4ParticleDefinition.hh"
//#include "G4ParticleTable.hh"
//#include "G4DecayTable.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::set;

ClassImp(MpdDecayer);

MpdDecayer *MpdDecayer::fgInstance = 0;

////////////////////////////////////////////////////////////////////////////////
/// Get the singleton object.

MpdDecayer *MpdDecayer::Instance()
{
   if (!fgInstance) fgInstance = new MpdDecayer;
   return fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

MpdDecayer::MpdDecayer()
   //    : fDecay(kMaxDecay),
   : fDecay(501), fBraPart(501)
{
   fBraPart.Reset(1);
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the decayer

void MpdDecayer::Init()
{
   static Bool_t init = kFALSE;
   if (init) return;
   init = kTRUE;
   // ForceDecay();
   fParticles = new TClonesArray("TParticle");

   // Get branching of Lmabda to p + \pi-
   TPythia6 *pythia = TPythia6::Instance();
   Int_t     ccLamb = pythia->Pycomp(3122); // compressed code for Lambda
   Int_t     idcb   = pythia->GetMDCY(ccLamb, 2);
   Int_t     idce   = idcb + pythia->GetMDCY(ccLamb, 3);
   ;
   // for (Int_t idc = idcb; idc < idce; ++idc) cout << pythia->GetBRAT(idc) << " ";
   // cout << endl;
   fBranch = pythia->GetBRAT(idcb);

   // Define additional particles
   DefineParticles();

   // Since the pi0 lifetime is less than 1.E-15 sec, if pi0 is produced as a decay
   // product in pythia6, e.g. KL0 -> pi0 pi+ pi-, the pi0 will be immediately
   // decayed by pythia6 to 2 gammas, and the KL0 decay product list passed
   // back to the transport mechanism will be "gamma gamma pi+ pi-", i.e.
   // the pi0 will not be visible in the list of secondaries passed back to
   // the transport mechanism and will not be pushed to the stack for possible
   // storage to the stdhep output array.
   // To avoid this, the pi0 is set to stable in pythia6, and its decay
   // will be handled by Geant.
   pythia->SetMDCY(pythia->Pycomp(111), 1, 0);

   if (TString(gMC->GetName()).Contains("TGeant4")) return;
   /*
   TGeant3 *g3 = (TGeant3*) gMC;
   // Get decay modes
   Int_t ipart = g3->IdFromPDG(3122);
   Int_t jpart = g3->Gclink()->jpart;
   Int_t jpa = g3->Lq()[jpart-ipart];
   Int_t jpa1 = g3->Lq()[jpa-1];
   Int_t jpa2 = g3->Lq()[jpa-2];
   for (Int_t i = 1; i <= 6; ++i) {
     fBratio[i-1] = g3->Q()[jpa1+i];
     fMode[i-1] = g3->Iq()[jpa2+i];
     std::cout << i << " " << fBratio[i-1] << " " << fMode[i-1] << std::endl;
   }
   */
}

////////////////////////////////////////////////////////////////////////////////
/// Decay a particle of type IDPART (PDG code) and momentum P.

void MpdDecayer::Decay(Int_t idpart, TLorentzVector *p)
{
   // Decay function

   //{ cout << " !!! Out !!! " << endl; exit(0); }
   // Reset internal particle container
   fParticles->Delete();
   fSourceFlag = kPythia;
   fMother.SetXYZT(0, 0, 0, 0);

   if (!p) return;
   TPythia6 *pythia = TPythia6::Instance();

   TParticle *part = gMC->GetStack()->GetCurrentTrack();
   // part->Print();

   if (fMothersPdg.find(idpart) == fMothersPdg.end()) {
      // Particle is not defined
      new ((*fParticles)[0]) TParticle(*part); // store mother particle
      fSourceFlag = kCustom;
      return;
   }

   TVector3 polar;
   part->GetPolarisation(polar);
   // polar.SetX(1.0); // just for test
   // polar.Print();
   // p->Print();
   // exit(0);
   if (TMath::Abs(polar.X()) < 0.0001 && TMath::Abs(polar.Z()) < 0.0001) {
      // if (polar.X() < -0.0001 && polar.Z() < 0.0001) {
      //  No polarisation - check this !!!
      pythia->Py1ent(0, idpart, p->Energy(), p->Theta(), p->Phi());
      pythia->GetPrimaries();
      /*
      Int_t kcLamb = pythia->Pycomp(3122); // compressed code for Lambda
      cout << " ----------------- " << pythia->GetMDCY(kcLamb,1) << " " << pythia->GetMDCY(kcLamb,2) << " "
      << pythia->GetMDCY(kcLamb,3) << endl;
      Int_t idcb = pythia->GetMDCY(kcLamb,2);
      Int_t idce = idcb + pythia->GetMDCY(kcLamb,3);;
      for (Int_t idc = idcb; idc < idce; ++idc) cout << pythia->GetBRAT(idc) << " ";
      cout << endl;
      */
   } else {
      // Polarized lambda - simulate anysotropic decay only to p + \pi-
      if (gRandom->Rndm() < fBranch) {
         // Anysotropic decay
         // cout << " ----------------- " << fBranch << endl;
         new ((*fParticles)[0]) TParticle(*part); // store mother particle
         Gdecay(idpart, p);
      } else {
         // Force the other channels
         Int_t idcb = pythia->GetMDCY(pythia->Pycomp(3122), 2);
         pythia->SetMDME(idcb, 1, 0); // disable decay to p + \pi-
         pythia->Py1ent(0, idpart, p->Energy(), p->Theta(), p->Phi());
         pythia->GetPrimaries();
         pythia->SetMDME(idcb, 1, 1); // enable decay to p + \pi-
      }
   }
   /*
   // Restore decay modes - to do the decay in Geant3
   TGeant3 *g3 = (TGeant3*) gMC;
   g3->Gsdk(g3->IdFromPDG(3122), fBratio, fMode);
   //g3decay_ ();
   g3decayaz_();
   g3->SetUserDecay(3122);
   */
}

////////////////////////////////////////////////////////////////////////////////
/// Get the decay products into the passed PARTICLES TClonesArray of
/// TParticles

Int_t MpdDecayer::ImportParticles(TClonesArray *particles)
{
   Int_t npart = 0;

   if (fSourceFlag == kPythia) {
      fParticles->Delete();
      // npart = TPythia6::Instance()->ImportParticles(particles,"All");
      npart = TPythia6::Instance()->ImportParticles(fParticles, "All");
      TLorentzVector lvDaugh;

      for (Int_t j = 0; j < npart; ++j) {
         TParticle *part = (TParticle *)fParticles->UncheckedAt(j);
         part->ProductionVertex(lvDaugh);
         lvDaugh += fMother;
         part->SetProductionVertex(lvDaugh);
         new ((*particles)[j]) TParticle(*part);
      }
   } else {
      npart = fParticles->GetEntriesFast();
      for (Int_t j = 0; j < npart; ++j) new ((*particles)[j]) TParticle(*((TParticle *)fParticles->UncheckedAt(j)));
   }
   /*
   cout << " Number of particles: " << npart << endl;
   for (Int_t j = 0; j < npart; ++j) {
     TParticle *part = (TParticle*) particles->UncheckedAt(j);
     part->Print();
     cout << " " << part->GetFirstMother() << " " << part->GetSecondMother() << " "
     << part->GetFirstDaughter() << " " << part->GetLastDaughter() << " " << part->GetNDaughters()
     << " " << part->GetStatusCode() << endl;
   }
   */
   return npart;
}

////////////////////////////////////////////////////////////////////////////////
/// Force a particular decay type

void MpdDecayer::SetForceDecay(Int_t type)
{
   // if (type > kMaxDecay) {
   if (type > 501) {
      Warning("SetForceDecay", "Invalid decay mode: %d", type);
      return;
   }
   // fDecay = EDecayType(type);
}

////////////////////////////////////////////////////////////////////////////////
/// Force a particle decay mode

void MpdDecayer::ForceDecay() {}

////////////////////////////////////////////////////////////////////////////////
/// Get the partial branching ratio for a particle of type IPART (a
/// PDG code).

Float_t MpdDecayer::GetPartialBranchingRatio(Int_t ipart)
{
   /*
     Int_t kc = TPythia6::Instance()->Pycomp(TMath::Abs(ipart));
     // return TPythia6::Instance()->GetBRAT(kc);
     return fBraPart[kc];
   */
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the life-time of a particle of type KF (a PDG code).

Float_t MpdDecayer::GetLifetime(Int_t kf)
{
   /*
     Int_t kc=TPythia6::Instance()->Pycomp(TMath::Abs(kf));
     return TPythia6::Instance()->GetPMAS(kc,4) * 3.3333e-12;
   */
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Define additional particles

void MpdDecayer::DefineParticles()
{

   SetDecayTableFile("pythia6.dat");
   std::ofstream fout(fDecayTableFile);
   // Lambda(1520)
   // PrintPDG(TDatabasePDG::Instance()->GetParticle(3124),fout);
   MakeDecayList(3124, fout);
   fout.close();
   ReadDecayTable();

   TPythia6 *pythia = TPythia6::Instance();
   Int_t     ccLamb = pythia->Pycomp(2214); // compressed code for phi
   Int_t     idcb   = pythia->GetMDCY(ccLamb, 2);
   Int_t     idce   = idcb + pythia->GetMDCY(ccLamb, 3);
   ;
   for (Int_t idc = idcb; idc < idce; ++idc) cout << pythia->GetBRAT(idc) << " ";
   cout << endl;
   for (Int_t idc = idcb; idc < idce; ++idc) cout << pythia->GetMDME(idc, 2) << " ";
   cout << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Read in particle data from an ASCII file.   The file name must
/// previously have been set using the member function
/// SetDecayTableFile.

void MpdDecayer::ReadDecayTable()
{
   if (fDecayTableFile.IsNull()) {
      Warning("ReadDecayTable", "No file set");
      return;
   }
   Int_t lun = 15;
   TPythia6::Instance()->OpenFortranFile(lun, const_cast<char *>(fDecayTableFile.Data()));
   TPythia6::Instance()->Pyupda(3, lun);
   TPythia6::Instance()->CloseFortranFile(lun);
}

// ===================================================================
// BEGIN COMMENT
//
// It would be better if the particle and decay information could be
// read from the current TDatabasePDG instance.
//
// However, it seems to me that some information is missing.  In
// particular
//
//   - The broadning cut-off,
//   - Resonance width
//   - Color charge
//   - MWID (?)
//
// Further more, it's not clear to me at least, what all the
// parameters Pythia needs are.
//
// Code like the below could be used to make a temporary file that
// Pythia could then read in.   Ofcourse, one could also manipulate
// the data structures directly, but that's propably more dangerous.
//
#if 1
void MpdDecayer::PrintPDG(TParticlePDG *pdg, std::ofstream &fout)
{
   TParticlePDG *anti     = pdg->AntiParticle();
   const char   *antiName = (anti ? anti->GetName() : "");
   Int_t         color    = 0;
   switch (TMath::Abs(pdg->PdgCode())) {
   case 1:
   case 2:
   case 3:
   case 4:
   case 5:
   case 6:
   case 7:
   case 8: // Quarks
      color = 1;
      break;
   case 21: // Gluon
      color = 2;
      break;
   case 1103:
   case 2101:
   case 2103:
   case 2203:
   case 3101:
   case 3103:
   case 3201:
   case 3203:
   case 3303:
   case 4101:
   case 4103:
   case 4201:
   case 4203:
   case 4301:
   case 4303:
   case 4403:
   case 5101:
   case 5103:
   case 5201:
   case 5203:
   case 5301:
   case 5303:
   case 5401:
   case 5403:
   case 5503:
      // Quark combinations
      color = -1;
      break;
   case 1000001:
   case 1000002:
   case 1000003:
   case 1000004:
   case 1000005:
   case 1000006: // super symmetric partners to quars
      color = 1;
      break;
   case 1000021: // ~g
      color = 2;
      break;
   case 2000001:
   case 2000002:
   case 2000003:
   case 2000004:
   case 2000005:
   case 2000006: // R hadrons
      color = 1;
      break;
   case 3000331:
   case 3100021:
   case 3200111:
   case 3100113:
   case 3200113:
   case 3300113:
   case 3400113:
      // Technicolor
      color = 2;
      break;
   case 4000001:
   case 4000002: color = 1; break;
   case 9900443:
   case 9900441:
   case 9910441:
   case 9900553:
   case 9900551:
   case 9910551: color = 2; break;
   }
   // std::cout << std::right
   fout << std::right << " " << std::setw(9) << pdg->PdgCode() << "  " << std::left << std::setw(16) << pdg->GetName()
        << "  " << std::setw(16) << antiName << std::right << std::setw(3) << Int_t(pdg->Charge()) << std::setw(3)
        << color << std::setw(3) << (anti ? 1 : 0) << std::fixed << std::setprecision(5) << std::setw(12) << pdg->Mass()
        << std::setw(12)
        << pdg->Width()
        //<< std::setw(12) << 0.0 // Broad
        << std::setw(12) << 4 * pdg->Width() // Broad
                                             //<< std::setw(12) << 8*pdg->Width() // Broad
        << std::scientific
        //<< " " << std::setw(13) << pdg->Lifetime()
        << std::setw(13) << pdg->Lifetime() << std::setw(3) << 0 // MWID
        << std::setw(3) << !pdg->Stable() << std::endl;
}

// ===================================================================

void MpdDecayer::MakeDecayList(Int_t pdgCode, std::ofstream &fout)
{
   TDatabasePDG *pdgDB = TDatabasePDG::Instance();
   if (!pdgDB->ParticleList()) pdgDB->ReadPDGTable();
   const THashList *pdgs = pdgDB->ParticleList();
   TParticlePDG    *pdg  = 0;
   TIter            nextPDG(pdgs);
   while ((pdg = static_cast<TParticlePDG *>(nextPDG()))) {
      // std::cout << "Processing " << pdg->GetName() << std::endl;
      if (pdg->PdgCode() != pdgCode) continue;
      PrintPDG(pdg, fout);

      TObjArray     *decays = pdg->DecayList();
      TDecayChannel *decay  = 0;
      TIter          nextDecay(decays);
      while ((decay = static_cast<TDecayChannel *>(nextDecay()))) {
         // std::cout << "Processing decay number " << decay->Number() << std::endl;
         Int_t ndecay = decay->NDaughters();
         fout << std::right << std::setw(15) << 1 << std::setw(5) << decay->MatrixElementCode() << std::fixed
              << std::setw(12) << std::setprecision(6) << decay->BranchingRatio();

         for (Int_t id = 0; id < 5; ++id) {
            if (id < ndecay)
               fout << std::setw(10) << decay->DaughterPdgCode(id);
            else
               fout << std::setw(10) << 0;
         }
         fout << "\n";
         //	"          %5d%5d%12.5f%10d%10d%10d%10d%10d\n"
      }
   }
}
#endif
// END COMMENT
// ===================================================================

////////////////////////////////////////////////////////////////////////////////
/// write particle data to an ASCII file.   The file name must
/// previously have been set using the member function
/// SetDecayTableFile.
///
/// Users can use this function to make an initial decay list file,
/// which then can be edited by hand, and re-loaded into the decayer
/// using ReadDecayTable.
///
/// The file syntax is
///
///     particle_list : partcle_data
///                   | particle_list particle_data
///                  ;
///     particle_data : particle_info
///                   | particle_info '\n' decay_list
///                   ;
///     particle_info : See below
///                   ;
///     decay_list    : decay_entry
///                   | decay_list decay_entry
///                   ;
///     decay_entry   : See below
///
/// The particle_info consists of 13 fields:
///
///     PDG code             int
///     Name                 string
///     Anti-particle name   string  if there's no anti-particle,
///                                  then this field must be the
///                                  empty string
///     Electic charge       int     in units of |e|/3
///     Color charge         int     in units of quark color charges
///     Have anti-particle   int     1 of there's an anti-particle
///                                  to this particle, or 0
///                                  otherwise
///     Mass                 float   in units of GeV
///     Resonance width      float
///     Max broadning        float
///     Lifetime             float
///     MWID                 int     ??? (some sort of flag)
///     Decay                int     1 if it decays. 0 otherwise
///
/// The format to write these entries in are
///
///     " %9  %-16s  %-16s%3d%3d%3d%12.5f%12.5f%12.5f%13.gf%3d%d\n"
///
/// The decay_entry consists of 8 fields:
///
///     On/Off               int     1 for on, -1 for off
///     Matrix element type  int
///     Branching ratio      float
///     Product 1            int     PDG code of decay product 1
///     Product 2            int     PDG code of decay product 2
///     Product 3            int     PDG code of decay product 3
///     Product 4            int     PDG code of decay product 4
///     Product 5            int     PDG code of decay product 5
///
/// The format for these lines are
///
///    "          %5d%5d%12.5f%10d%10d%10d%10d%10d\n"
///

void MpdDecayer::WriteDecayTable()
{
   if (fDecayTableFile.IsNull()) {
      Warning("ReadDecayTable", "No file set");
      return;
   }
   Int_t lun = 15;
   TPythia6::Instance()->OpenFortranFile(lun, const_cast<char *>(fDecayTableFile.Data()));
   TPythia6::Instance()->Pyupda(1, lun);
   TPythia6::Instance()->CloseFortranFile(lun);
}

////////////////////////////////////////////////////////////////////////////////
/// Count number of decay products

Int_t MpdDecayer::CountProducts(Int_t channel, Int_t particle)
{
   /*
     Int_t np = 0;
     for (Int_t i = 1; i <= 5; i++)
     if (TMath::Abs(TPythia6::Instance()->GetKFDP(channel,i)) == particle) np++;
     return np;
   */
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Simulate two-body decay process of Lambda-hyperon to p+\pi-

void MpdDecayer::Gdecay(Int_t idpart, TLorentzVector *p)
{

   Double_t pcm[2][4] = {{0}, {0}};
   Double_t xm0 = p->M(), xm1 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
   Double_t xm2 = TDatabasePDG::Instance()->GetParticle(-211)->Mass();

   Gdeca2(xm0, xm1, xm2, pcm);

   // First, second decay products.

   TLorentzVector daughter1, daughter2;
   daughter1.SetPxPyPzE(pcm[0][0], pcm[0][1], pcm[0][2], pcm[0][3]);
   daughter2.SetPxPyPzE(pcm[1][0], pcm[1][1], pcm[1][2], pcm[1][3]);

   // Boost to lab frame
   // TVector3 boost = p->BoostVector(); //AZ - this is wrong !!!!
   TVector3 boost(0, 0, p->P() / p->E());
   daughter1.Boost(boost);
   daughter2.Boost(boost);
   daughter1.Transform(fRotation);
   daughter2.Transform(fRotation);

   TLorentzVector pos;
   gMC->TrackPosition(pos);
   // pos.Print();

   Int_t npart = fParticles->GetEntriesFast();
   new ((*fParticles)[npart++]) TParticle(2212, 1, 1, -1, 0, 0, daughter1, pos);
   new ((*fParticles)[npart]) TParticle(-211, 1, 1, -1, 0, 0, daughter2, pos);

   fSourceFlag = kCustom;

   // Modify mother particle parameters
   TParticle *mother = (TParticle *)fParticles->UncheckedAt(0);
   mother->SetFirstMother(1);
   mother->SetLastMother(-1);
   mother->SetFirstDaughter(2);
   mother->SetLastDaughter(3);
   mother->SetStatusCode(11);
}

////////////////////////////////////////////////////////////////////////////////
/// Simulate two-body decay process with anisotropic angular distribution in CMS.

void MpdDecayer::Gdeca2(Double_t xm0, Double_t xm1, Double_t xm2, Double_t pcm[2][4])
{

   Double_t random[2], pvert[3], costh = 0.0, sinth = 0.0, phi = 0.0;

   Double_t e1 = (xm0 * xm0 + xm1 * xm1 - xm2 * xm2) / (2. * xm0);
   Double_t p1 = TMath::Sqrt(TMath::Abs((e1 - xm1) * (e1 + xm1)));

   static TRandom *myRandom = NULL;
   if (!myRandom) {
      myRandom = new TRandom();
      myRandom->SetSeed(0);
   }
   TRandom *oldRandom = gRandom;
   gRandom            = myRandom;

   gRandom->RndmArray(2, random);

   // Sanity check - should not happen (legacy code)
   TParticle *part = gMC->GetStack()->GetCurrentTrack();
   if (part->GetPdgCode() != 3122) {
      cout << " ??? Not Lambda - exit " << part->GetPdgCode() << endl;
      exit(0);
   }

   TLorentzVector lorV;
   part->Momentum(lorV);
   lorV.Vect().GetXYZ(pvert);

   TVector3 polar;
   part->GetPolarisation(polar);

   if (TMath::Abs(polar.X()) < 0.0001 && TMath::Abs(polar.Z()) < 0.0001) {
      costh = 2. * random[0] - 1.0;
      phi   = TMath::TwoPi() * random[1];
   } else {
      // Anisotropic decay angular distribution (0.5*(1+alpha*P*cos)).
      Anisotropy(pvert, random, polar.X(), phi, costh);
   }

   if (TMath::Abs(costh) >= 1.0)
      costh = TMath::Sign(1.0, costh);
   else
      sinth = TMath::Sqrt((1. - costh) * (1. + costh));

   // Polar co-ordinates to momentum components.

   pcm[0][0] = p1 * sinth * TMath::Cos(phi);
   pcm[0][1] = p1 * sinth * TMath::Sin(phi);
   pcm[0][2] = p1 * costh;
   pcm[0][3] = e1;

   // Generate second decay product.

   pcm[1][0] = -pcm[0][0];
   pcm[1][1] = -pcm[0][1];
   pcm[1][2] = -pcm[0][2];
   pcm[1][3] = TMath::Sqrt(pcm[1][0] * pcm[1][0] + pcm[1][1] * pcm[1][1] + pcm[1][2] * pcm[1][2] + xm2 * xm2);

   gRandom = oldRandom;
}

////////////////////////////////////////////////////////////////////////////////
/// Simulate anisotropic (due to polarization) decay of lambda

void MpdDecayer::Anisotropy(Double_t *pvert, Double_t *rndm, Double_t polar, Double_t &phi, Double_t &costh)
{
   // Simulate anisotropic (due to polarization) decay of lambda

   // exit(0);
   // std::ofstream outfile ("costh.txt");
   // freopen ("costh.txt", "w", stdout);

   // const Double_t alpha = 0.642;
   static TF1 f("f", "1+0.642*[0]*x", -1, 1);
   f.SetParameter(0, polar);
   f.SetNpx(1000);

   Double_t costhe = f.GetRandom();

   Double_t sinth = 0.0;
   if (TMath::Abs(costhe) >= 1.0)
      costhe = TMath::Sign(1.0, costhe);
   else
      sinth = TMath::Sqrt(1.0 - costhe * costhe);
   phi = TMath::TwoPi() * rndm[1];
   // std::cout << " Cos: " << costhe << " " << phi << " " << pvert[0] << " " << pvert[1] << " " << pvert[2] <<
   // std::endl;

   // Compute normal vector
   TVector3 beam(0.0, 0.0, 1.0);
   TVector3 lambda(pvert[0], pvert[1], pvert[2]);
   TVector3 norm = beam.Cross(lambda).Unit();

   // Unit vector with theta and phi w.r.t. normal
   TVector3 unit(sinth * TMath::Cos(phi), sinth * TMath::Sin(phi), costhe);

   // Coordinate system transformation
   // (from lambda to lab)
   TVector3 lambU = lambda.Unit(), lambZ(0, 0, 1), lambX(1, 0, 0), lambY(0, 1, 0);
   lambZ.RotateUz(lambU);
   lambX.RotateUz(lambU);
   lambY.RotateUz(lambU);
   TRotation rotL;
   rotL.RotateAxes(lambX, lambY, lambZ); // transformation to lab. system
   fRotation = rotL;
   rotL.Invert();        // from lab. to lambda
   unit.RotateUz(norm);  // to lab. system
   unit.Transform(rotL); // to lambda system
   /*
   TVector3 normZ(0,0,1), normX(1,0,0), normY(0,1,0);
   normZ.RotateUz(norm);
   normX.RotateUz(norm);
   normY.RotateUz(norm);
   TRotation rotN;
   rotN.RotateAxes(normX,normY,normZ); // transformation to lab. system
   unit.Transform(rotN); // to lab. system
   unit.Transform(rotL); // to lambda system
   */
   costh = TMath::Cos(unit.Theta());
   phi   = unit.Phi();
   // std::cout << costh << " " << phi << std::endl;
   // outfile << i << " " << costh << endl;
   // outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
/// Add PDG of particle to be decayed by this package

void MpdDecayer::AddMotherPdg(Int_t pdg)
{

   if (fMothersPdg.find(pdg) != fMothersPdg.end()) return;
   fMothersPdg.insert(pdg);
   /*
   G4DecayTable *decayTbl = NULL;
   if (TString(gMC->GetName()) == "TGeant4") {
     // Get decay table to feed it in again - bypass TGeant4 feature of zeroing decay table
     G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
     G4ParticleDefinition* particleDefinition = particleTable->FindParticle(pdg);
     decayTbl = particleDefinition->GetDecayTable();
   }
   */
   gMC->SetUserDecay(pdg);
   /*
   if (TString(gMC->GetName()) == "TGeant4") {
     G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
     G4ParticleDefinition* particleDefinition = particleTable->FindParticle(pdg);
     particleDefinition->SetDecayTable(decayTbl);
   }
   */
}

////////////////////////////////////////////////////////////////////////////////
/// Decay a TParticle

void MpdDecayer::Decay(TParticle *p)
{
   // Decay function

   //{ cout << " !!! Out !!! " << endl; exit(0); }
   // Reset internal particle container
   fParticles->Delete();
   fSourceFlag = kPythia;

   if (!p) return;
   p->ProductionVertex(fMother);
   TPythia6 *pythia = TPythia6::Instance();

   pythia->Py1ent(0, p->GetPdgCode(), p->Energy(), p->Theta(), p->Phi());
   pythia->GetPrimaries();
}
