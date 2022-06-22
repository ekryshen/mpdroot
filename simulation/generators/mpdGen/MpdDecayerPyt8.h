// -------------------------------------------------------------------------
// -----                 MpdDecayerPyt8 header file                    -----
// -----                      Created 23/07/2020                       -----
// -----                          A. Zinchenko                         -----
// -----             External decayer for MPD using Pythia8            -----
// -----                (adapted from TPythia8Decayer)                 -----
// -------------------------------------------------------------------------

#ifndef MPDDECAYERPYT8_H
#define MPDDECAYERPYT8_H

#include <TLorentzVector.h>
#include <TVirtualMCDecayer.h>
#include <TRotation.h>

#include <Pythia8/Pythia.h>

class TClonesArrray;
class TObjArrray;
class TParticle;
class TPythia8;
class TRandom;

#include <set>
#if PYTHIA_VERSION_INTEGER < 8304
typedef Pythia8::ParticleDataEntry *ParticleDataEntryPtr;
#else
typedef Pythia8::ParticleDataEntryPtr ParticleDataEntryPtr;
#endif

class MpdDecayerPyt8 : public TVirtualMCDecayer {

public:
   enum SourceFlag { kPythia, kCustom }; // particle container source

   virtual ~MpdDecayerPyt8() { ; }

   virtual void    Init();
   virtual void    Decay(Int_t pdg, TLorentzVector *p);
   virtual Int_t   ImportParticles(TClonesArray *particles);
   virtual void    SetForceDecay(Int_t type);
   virtual void    ForceDecay();
   virtual Float_t GetPartialBranchingRatio(Int_t ipart);
   virtual Float_t GetLifetime(Int_t kf);
   virtual void    ReadDecayTable();

   virtual void SetDebugLevel(Int_t debug) { fDebug = debug; }

   static MpdDecayerPyt8 *Instance();
   void                   Decay(TParticle *part);
   void                   AddMotherPdg(Int_t pdg);
   // using   TPythia8Decayer::Decay; // makes base class function visible

protected:
   MpdDecayerPyt8();
   void AppendParticle(Int_t pdg, TLorentzVector *p);
   void ClearEvent();

private:
   void Gdecay(Int_t idpart, TLorentzVector *p);
   void Gdeca2(Double_t xm0, Double_t xm1, Double_t xm2, Double_t pcm[2][4]);
   void Anisotropy(Double_t *rndm, TVector3 &polar, Double_t &phi, Double_t &costh);
   void ChangeBranchings();
   void ChangeParticleBr(TObjArray *tokens);

   static MpdDecayerPyt8 *fgInstance; // singleton instance

   TPythia8 *fPythia8;     // Pointer to pythia8
   Int_t     fDebug;       // Debug level

   Float_t              fBranch; // branching of lambda to p + \pi-
   ParticleDataEntryPtr fLambda; // pointer to Lambda in Pythia8 database
   TLorentzVector       fMother; // mother decay point position
   TClonesArray        *fParticles;
   SourceFlag           fSourceFlag;
   std::set<Int_t>      fMothersPdg;
   TRandom             *fRandom;

   ClassDef(MpdDecayerPyt8, 0); // MPD Particle Decayer using Pythia8
};
#endif
