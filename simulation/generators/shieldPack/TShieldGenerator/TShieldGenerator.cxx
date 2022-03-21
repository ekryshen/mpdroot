// -------------------------------------------------------------------------
// -----          TShieldGenerator source file                         -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGenerator.h"

TShieldGenerator::TShieldGenerator()
   : TGenerator(), fEnergy(1E3), fStartX(0), fStartY(0), fStartZ(0), fPdgCode(1), fNumStat(0)
{
   fShield = new TShield();
}
TShieldGenerator::~TShieldGenerator()
{
   delete fShield;
}

void TShieldGenerator::SetGeometry(TGeoManager *geom)
{
   fGeoManager = geom;
}
void TShieldGenerator::SetEnergy(Float_t energy)
{
   fEnergy = energy;
}
void TShieldGenerator::SetStartPoint(Float_t x, Float_t y, Float_t z)
{
   fStartX = x;
   fStartY = y;
   fStartZ = z;
}
void TShieldGenerator::SetParticleFromPdgCode(Int_t pdg_code, Float_t A, Float_t Z)
{
   fPdgCode = pdg_code;
   fA       = A;
   fZ       = Z;
}
void TShieldGenerator::SetStatistics(Int_t nstat)
{
   fNumStat = nstat;
}

void TShieldGenerator::GenerateEvent()
{
   srand(time(0));
   printf(":::TShieldGenerator:GenerateEvent\n");
   THadgen *currGen = new THadgen();
   currGen->SetAutoSeed();
   currGen->SetRandomSeed(clock() + time(0));
   currGen->SetParticleFromPdgCode(fPdgCode, fA, fZ);
   currGen->SetEnergy(fEnergy);
   currGen->SetNumStat(fNumStat);
   currGen->GenerateEvent();
   TClonesArray *p = new TClonesArray("TParticles");
   currGen->ImportParticles(p, "");

   fShield->SetGeometry(fGeoManager);
   fShield->SetStartPoint(fStartX, fStartY, fStartZ);

   TParticle *par        = 0;
   int        numEntries = ((TClonesArray *)p)->GetEntriesFast();
   printf(":::TShieldGenerator:numEntries = %i\n", numEntries);
   for (Int_t k = 0; k < numEntries; k++) {
      par = (TParticle *)(((TClonesArray *)p)->At(k));
      fShield->SetRandomSeed(clock() + time(0));
      fShield->SetEnergy(par->Energy());
      fShield->SetParticleFromPdgCode(par->GetPdgCode());
      fShield->SetDirectionVector(par->Px(), par->Py(), par->Pz());

      fShield->GenerateEvent();
   }
}
TObjArray *TShieldGenerator::ImportParticles(Option_t *option)
{
   ImportParticles((TClonesArray *)fParticles, option);
   return fParticles;
}
Int_t TShieldGenerator::ImportParticles(TClonesArray *particles, Option_t *option)
{
   //     if (!strcmp(option,"") || !strcmp(option,"Final")) { //Only stable particles;
   //     } else if (!strcmp(option,"All")) { //All particles
   TClonesArray &a = *particles;
   a.Clear();

   // Int_t pdg;
   TTree     *tree = (TTree *)fShield->GetTree();
   TParticle *p    = 0;
   int        numEntries;
   tree->SetBranchAddress("Particle", &p);
   numEntries = ((TTree *)tree)->GetEntriesFast();
   for (Int_t k = 0; k < numEntries; k++) {
      tree->GetEntry(k);
      //             if(p->GetStatusCode()!=3)continue;
      // pdg = p->GetPdgCode();
      new (a[k]) TParticle(*p);
   }
   return numEntries;
}

ClassImp(TShieldGenerator)