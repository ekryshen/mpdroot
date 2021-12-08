// -------------------------------------------------------------------------
// -----          TShield source file                                  -----
// -----          Created by A. Timofeev                               -----
// -------------------------------------------------------------------------

#include "TShield.h"

double TShield::GeoDistCallback(double x, double y, double z, double vx, double vy, double vz)
{
   TShield *fShield = TShield::Instance();
   if (fShield->fGeometry != 0)
      return fShield->FindBoundary(x, y, z, vx, vy, vz);
   else {
      printf("TShield FATAL ERROR: TGeoManager not set\n");
      return -1;
   }
}

int TShield::GeoNextCallback(double x, double y, double z, double vx, double vy, double vz)
{
   TShield *fShield = TShield::Instance();
   if (fShield->fGeometry != 0)
      return fShield->GetNextVolume(x, y, z, vx, vy, vz);
   else {
      printf("TShield FATAL ERROR: TGeoManager not set\n");
      return -1;
   }
}

void TShield::TreeCallback(shield_tree_node *node)
{
   // anyway, all particles are always registered in the tree
   TShield::Instance()->AddTreeParticle(node); // Update fCurrentParticle and fill fTree
   switch (node->event) {
   case EVENT_FLYOUT: TShield::Instance()->AddFliedOutParticle(node); break;
   case EVENT_ABSORPTION: TShield::Instance()->AddAbsorbedParticle(node); break;
   default: break;
   }
}

void TShield::InitCallbacks()
{
   shield_set_tree_callback(&TShield::TreeCallback);
   //     shield_set_gcurzl_callback(&TShield::GeoDistCallback); //Disabled now
   shield_set_gnextz_callback(&TShield::GeoNextCallback);
}

void TShield::AddFliedOutParticle(shield_tree_node *node)
{
   if (ParticlesNumFlyOut <= ParticlesNumFlyOutMax) {
      TClonesArray &a = *fParticlesFlyOut;
      new (a[ParticlesNumFlyOut]) TParticle(*fCurrentParticle);
      ParticlesNumFlyOut++;
   } else {
      printf("TShield ERROR: TClonesArray fParticlesFlyOut overfull\n");
   }
}

void TShield::AddAbsorbedParticle(shield_tree_node *node)
{
   if (ParticlesNumAbsorbed <= ParticlesNumAbsorbedMax) {
      TClonesArray &a = *fParticlesAbsorbed;
      new (a[ParticlesNumAbsorbed]) TParticle(*fCurrentParticle);
      ParticlesNumAbsorbed++;
   } else {
      printf("TShield ERROR: TClonesArray fParticlesAbsorbed overfull\n");
   }
}

void TShield::AddTreeParticle(shield_tree_node *node)
{
   new (fCurrentParticle) TParticle(node->pdg, node->event, 0, 0, 0, 0, node->px_z, node->py_z, node->pz_z,
                                    node->tz + node->weight, node->xz, node->yz, node->zz, 0);
   fTree->Fill();
}

double TShield::FindBoundary(double x, double y, double z, double vx, double vy, double vz)
{
   // TGeoNode *node = fGeometry->FindNode(x, y, z);
   return 0; // line 45: Disabled now
}

int TShield::GetNextVolume(double x, double y, double z, double vx, double vy, double vz)
{
   shield_clean_geometry();
   Geant4ShieldData geometryData = getDataAtPoint(fGeometry, x, y, z, vx, vy, vz);
#ifdef _DEBUG
   printf("TShield::GetNextVolume(%f,%f,%f,%f,%f,%f);\n", x, y, z, vx, vy, vz);
   tgeanttoshield::printGeant4ShieldData(geometryData);
#endif
   bool isOnlyOuterVacuum = true;
   for (UInt_t i = 0; i < geometryData.size(); i++) {
      int countBodies = geometryData.at(i).bodyVector.size();
      // at "23.2.4 Class template vector" says, that elements of a vector are stored contiguously.
      int      *zones  = &(geometryData.at(i).zoneVector[0]);
      SGeoBody *bodies = &(geometryData.at(i).bodyVector[0]);
      if (shield_add_zone(countBodies, zones, bodies, geometryData.at(i).medium) == -1) {
         printf("An error was accured");
      }
      isOnlyOuterVacuum = isOnlyOuterVacuum && (geometryData.at(i).medium.nType == 0);
   }
   shield_init_medium();
   shield_init_geometry();
   return (isOnlyOuterVacuum) ? 0 : 1;
}