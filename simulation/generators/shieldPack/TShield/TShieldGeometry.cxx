// -------------------------------------------------------------------------
// -----          TShieldGeometry source file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGeometry.h"

bool operator<(const tgeanttoshield::zoneElement z1, const tgeanttoshield::zoneElement z2)
{
   if (z1.first != z2.first) {
      return (z1.first > z2.first);
   }
   if (z1.second.type != z2.second.type) {
      return (z1.second.type < z2.second.type);
   }
   for (int k = 0; k < 36; ++k) {
      if (tgeanttoshield::doubleNE(z1.second.parameters[k], z2.second.parameters[k])) {
         return (z1.second.parameters[k] < z2.second.parameters[k]);
      }
   }
   return false;
}

Geant4ShieldData getDataAtPoint(TGeoManager *geoMan, Double_t x, Double_t y, Double_t z, Double_t vx, Double_t vy,
                                Double_t vz)
{
   using namespace tgeanttoshield;
   TGeoNavigator *nav = geoMan->GetCurrentNavigator();
   nav->DoBackupState();
   nav->CdTop();
   nav->SetCurrentDirection(vx, vy, vz);
   zoneData      outputData;
   static double s = TGeoShape::Tolerance();
   zoneData      zone;
   zone = getZonesAtPoint(nav, x + vx * s, y + vy * s, z + vz * s);
   outputData.insert(outputData.end(), zone.begin(), zone.end());
   zone = getZonesAtPoint(nav, x - vx * s, y - vy * s, z - vz * s);
   outputData.insert(outputData.end(), zone.begin(), zone.end());
   nav->DoRestoreState();
   return zoneDataToShield(outputData);
}

namespace tgeanttoshield {
zoneData getZonesAtPoint(TGeoNavigator *nav, Double_t x, Double_t y, Double_t z)
{
   TGeoNode *node = nav->FindNode(x, y, z);
   // Now, current position is deepest volume with point (x,y,z)
   nav->Safety(); // Set correct fBoundary and fOutside values
   zoneData outputData;
   if (!nav->IsOutside() && node) {
      outputData = getCurrentZone(nav);
   } else {
      // Is outer vacuum needed or return empty data?
      nav->CdTop();
      outputData.push_back(addOuterVacuum(nav->GetCurrentNode()));
      //         outputData.push_back(addOuterVacuum(geoMan->GetTopNode()));
   }
   return outputData;
}
} // namespace tgeanttoshield

namespace tgeanttoshield {
zoneData getCurrentZone(TGeoNavigator *nav)
{
   zoneData    outputData;
   int         countOfDaughters   = nav->GetCurrentNode()->GetNdaughters();
   TGeoVolume *volume             = nav->GetCurrentVolume();
   TGeoHMatrix currTransformation = *(nav->GetCurrentMatrix()); // getGlobalMatrix(nav);
   zoneList    internalVolume     = (getZoneFromShape(volume->GetShape(), currTransformation)).first;
   TGeoNode   *daughter;
   zoneList    daughterOuterZones;
   for (int i = 0; i < countOfDaughters; ++i) {
      daughter           = volume->GetNode(i);
      daughterOuterZones = getZoneFromBody(daughter, currTransformation).second;
      internalVolume     = andZone(internalVolume, notZone(daughterOuterZones)); // Only outel shell of body.
   }
   outputData.push_back(std::make_pair(internalVolume, getMedium(nav->GetCurrentNode())));
   return outputData;
}

Geant4ShieldData zoneDataToShield(zoneData data)
{
   Geant4ShieldData                   out;
   Geant4ShieldElement                tmpElement;
   std::pair<zoneList, MediumData>    currentPair;
   std::vector<SGeoBody>              bodyVector;
   std::vector<int>                   zoneVector;
   std::set<std::vector<zoneElement>> geoData; // Converting to set for removing duplicates.
   std::set<zoneElement>              curSet;  // Set have only unique keys.
   for (unsigned int kk = 0; kk < data.size(); ++kk) {
      bodyVector.clear();
      zoneVector.clear();
      currentPair = data.at(kk); // This is data for one volume
      geoData     = std::set<std::vector<zoneElement>>(currentPair.first.begin(), currentPair.first.end());
      for (std::set<std::vector<zoneElement>>::const_iterator k = geoData.begin(); k != geoData.end(); ++k) {
         curSet = std::set<zoneElement>(
            (std::vector<zoneElement>::const_iterator)k->begin(),
            (std::vector<zoneElement>::const_iterator)k->end()); // As fact, removing duplicates of zoneElements
         for (std::set<zoneElement>::const_iterator i = curSet.begin(); i != curSet.end(); ++i) {
            zoneVector.push_back(i->first);
            bodyVector.push_back(i->second);
         }
         zoneVector.push_back(0);
      }
      zoneVector.pop_back(); // Removing of tail '0' element.
      tmpElement.medium     = currentPair.second;
      tmpElement.zoneVector = zoneVector;
      tmpElement.bodyVector = bodyVector;
      out.push_back(tmpElement);
   }
   return out;
}

std::pair<zoneList, MediumData> addOuterVacuum(TGeoNode *world)
{
   zoneList   out        = getZoneFromBody(world, TGeoHMatrix(), 1.5 * getDefaultScale()).second;
   zoneList   outerWorld = getZoneFromBody(world, TGeoHMatrix()).second;
   MediumData medium     = {0, 0, 0, {}};
   return std::make_pair(andZone(out, notZone(outerWorld)), medium);
}
} // namespace tgeanttoshield

// Obsolete versions:
Geant4ShieldData convertTGeantToShield(TGeoManager *geoMan)
{
   printf("TShieldGeometry: using obsolete convertTGeantToShield function!\n");
   using namespace tgeanttoshield;
   TGeoNode *world     = geoMan->GetTopNode();
   zoneData  zoneWorld = convertBodyRecursive(world);
   //     zoneWorld.push_back(addOuterVacuum(world));
   return zoneDataToShield(zoneWorld);
}
namespace tgeanttoshield {
zoneData convertBodyRecursive(TGeoNode *node)
{
   return convertBodyRecursive(node, TGeoHMatrix()).first;
}
std::pair<zoneData, zoneList> convertBodyRecursive(TGeoNode *node, TGeoHMatrix oldTransformation)
{
   zoneData                      out;
   std::pair<zoneList, zoneList> tmp               = getZoneFromBody(node, oldTransformation);
   TGeoHMatrix                   newTransformation = *(node->GetMatrix());
   TGeoHMatrix                   currMatrix        = oldTransformation * newTransformation;
   zoneList                      internalVolume    = tmp.first;
   zoneList                      outerVolume       = tmp.second;

   TGeoVolume                   *volume           = node->GetVolume();
   int                           countOfDaughters = volume->GetNdaughters();
   TGeoNode                     *daughter;
   std::pair<zoneData, zoneList> daughterZones;
   for (int i = 0; i < countOfDaughters; ++i) {
      daughter      = volume->GetNode(i);
      daughterZones = convertBodyRecursive(daughter, currMatrix);
      out.insert(out.end(), daughterZones.first.begin(),
                 daughterZones.first.end());                                   // Insert daughter zones to output zones
      internalVolume = andZone(internalVolume, notZone(daughterZones.second)); // Only outel shell of body.
   }
   //     out.push_back(std::make_pair(internalVolume, getMedium(node)));
   out.insert(out.begin(),
              std::make_pair(internalVolume, getMedium(node))); // Add internalVolume to first place at output.
   return std::make_pair(out, outerVolume);
}
} // namespace tgeanttoshield
