// -------------------------------------------------------------------------
// -----          TShield source file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

#include "TShield.h"

// class TShield implementation
TShield *TShield::fInstance = 0;

//______________________________________________________________________________
TShield::TShieldClean::~TShieldClean()
{
   // delete the unique copy of TShield
   // if exists
   if (TShield::fInstance) delete TShield::fInstance;
}
//______________________________________________________________________________

ClassImp(TShield)

   //______________________________________________________________________________
   TShield::TShield(bool clean, int FlyOutMax, int AbsorbedMax)
   : fGeometry(nullptr), fParticlesFlyOut(nullptr), fParticlesAbsorbed(nullptr), ParticlesNumFlyOutMax(FlyOutMax),
     ParticlesNumAbsorbedMax(AbsorbedMax), clearAtStart(clean), fTree(nullptr), fCurrentParticle(nullptr)
{
   // Check if another instance of TShield already exists
   if (fInstance) {
      Fatal("TShield", "FATAL ERROR: Another instance of TShield class already exists");
   }
   fInstance = this;
   // create static cleaner
   static TShieldClean cleaner;

   SetAutoSeed();
   // Create arrays for particles storage
   fParticlesFlyOut   = new TClonesArray("TParticle", ParticlesNumFlyOutMax);
   fParticlesAbsorbed = new TClonesArray("TParticle", ParticlesNumAbsorbedMax);
   fGeometry          = 0;
   // init the random generator
   srand(time(0));
   // init callbacks for tree generation and geometry
   InitCallbacks();

   void *t = 0;
   fTree   = new TTree("TShieldTree", "Particles from TShield");
   //     fTree->Branch("particle", "TParticle", &fCurrentParticle); //Don't work, undefined symbol: _ZTI9TParticle
   TBranch *branch = fTree->Branch("Particle", "TParticle", &(t), 32000, 0);
   branch->SetAddress(&fCurrentParticle);
   fCurrentParticle = new TParticle();

   // Setting all default
   shield_set_defaults();
   ClearArrays();
}

TShield *TShield::Instance()
{
   return fInstance ? fInstance : new TShield();
}

//______________________________________________________________________________
TShield::~TShield()
{
   // delete TClonesArray's
   delete fParticlesFlyOut;
   delete fParticlesAbsorbed;
   delete fTree;
   delete fCurrentParticle;
   // if deleted, clear the instance
   TShield::fInstance = 0;
   // useful if needed to create a new instance of TShield in future
   // actually TShield object should be automatically cleaned
   // but user can also delete it manually
}

//______________________________________________________________________________
void TShield::GenerateEvent()
{
   // Clean particles arrays
   if (clearAtStart) {
      ClearArrays();
   }
   shield_clean_geometry();
   shield_run();
}

//______________________________________________________________________________
void TShield::ClearArrays()
{
   // Clean particles arrays
   fParticlesFlyOut->Clear();
   fParticlesAbsorbed->Clear();
   ParticlesNumFlyOut   = 0;
   ParticlesNumAbsorbed = 0;
   fTree->Clear();
}

//______________________________________________________________________________
void TShield::SetGeometry(TGeoManager *geom)
{
   fGeometry = geom;
}

//______________________________________________________________________________
void TShield::PrintGeometry()
{
   printGeometry(fGeometry);
}

//______________________________________________________________________________
TGeoManager *TShield::GetGeometry()
{
   return fGeometry;
}