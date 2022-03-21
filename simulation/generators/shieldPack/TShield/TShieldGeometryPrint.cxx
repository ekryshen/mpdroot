// -------------------------------------------------------------------------
// -----          TShieldGeometry source file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGeometry.h"

void printGeometry(TGeoManager *geoMan)
{
   using namespace tgeanttoshield;
   TGeoNavigator *nav = geoMan->GetCurrentNavigator();
   nav->DoBackupState();
   nav->CdTop();
   printGeometryBranch(nav);
   nav->DoRestoreState();
}
namespace tgeanttoshield {
void printElement(Element element)
{
   printf("Element: { %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f}\n", element.Nuclid, element.Conc, element.Density, element.Z,
          element.A, element.PureDensity, element.ionEv);
}
void printMediumData(MediumData medium)
{
   printf("Medium:\n");
   printf("\tType: %i\n", medium.nType);
   printf("\tRho: %f\n", medium.Rho);
   printf("\tCount of elements: %i\n", medium.nChemEl);
   for (int i = 0; i < medium.nChemEl; i++) {
      printf("\tElement %i: ", i);
      printElement(medium.Elements[i]);
   }
   printf("\n");
}
void printSGeoBody(SGeoBody body)
{
   printf("Body: %i, { ", body.type);
   for (int i = 0; i < 36; i++) {
      printf("%f ", body.parameters[i]);
   }
   printf("}\n");
}
void printSGeoZone(SGeoZone zone)
{
   printf("Zone:\n");
   printf("\tNumber of medium: %i\n", zone.mediaNumber);
   printf("\tCount of elements: %i\n", zone.countELements);
   printf("\t");
   for (int i = 0; i < zone.countELements; i++) {
      printf("%i ", zone.definition[2 * i] * zone.definition[2 * i + 1]);
   }
   printf("\n");
}
void printGeant4ShieldElement(Geant4ShieldElement element)
{
   printf("Zones: ");
   for (unsigned int k = 0; k < element.zoneVector.size(); ++k) {
      printf("%i ", element.zoneVector.at(k));
   }
   printf("\n  Bodies:\n");
   for (unsigned int k = 0; k < element.bodyVector.size(); ++k) {
      printf("\t");
      printSGeoBody(element.bodyVector.at(k));
   }
   printf("  ");
   printMediumData(element.medium);
}
void printGeant4ShieldData(Geant4ShieldData data)
{
   for (unsigned int i = 0; i < data.size(); ++i) {
      printGeant4ShieldElement(data.at(i));
   }
}
void printGeometryBranch(TGeoNavigator *nav)
{
   printGeant4ShieldData(zoneDataToShield(getCurrentZone(nav)));
   int      countOfDaughters = nav->GetCurrentNode()->GetNdaughters();
   zoneList daughterOuterZones;
   for (int i = 0; i < countOfDaughters; ++i) {
      nav->CdDown(i);
      printGeometryBranch(nav);
      nav->CdUp();
   }
}
} // namespace tgeanttoshield
