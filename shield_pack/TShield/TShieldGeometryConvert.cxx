// -------------------------------------------------------------------------
// -----          TShieldGeometry source file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGeometry.h"
namespace tgeanttoshield {
std::pair<zoneList, zoneList> getZoneFromBody(TGeoNode *node, TGeoHMatrix oldTransformation, double scale)
{
   TGeoHMatrix newTransformation  = *(node->GetMatrix());
   TGeoHMatrix currTransformation = oldTransformation * newTransformation;
   TGeoShape  *shape              = node->GetVolume()->GetShape();
   return getZoneFromShape(shape, currTransformation, scale);
}

std::pair<zoneList, zoneList> getZoneFromShape(TGeoShape *shape, TGeoHMatrix currTransformation, double scale)
{
   std::pair<zoneList, zoneList> zonePair;
   if (dynamic_cast<TGeoTubeSeg *>(shape) != nullptr) {
      zonePair = tubeSegToZones((TGeoTubeSeg *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoTube *>(shape) != nullptr) {
      zonePair = tubeToZones((TGeoTube *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoCtub *>(shape) != nullptr) { // TODO
      printf("TGeoCtub not implemented yet\n");
   } else if (dynamic_cast<TGeoConeSeg *>(shape) != nullptr) {
      zonePair = coneSegToZones((TGeoConeSeg *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoCone *>(shape) != nullptr) {
      zonePair = coneToZones((TGeoCone *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoPara *>(shape) != nullptr) { // TODO
      printf("TGeoPara not implemented yet\n");
   } else if (dynamic_cast<TGeoTrd1 *>(shape) != nullptr) {
      zonePair = trd1ToZones((TGeoTrd1 *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoTrd2 *>(shape) != nullptr) {
      zonePair = trd2ToZones((TGeoTrd2 *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoTrap *>(shape) != nullptr) {
      zonePair = trapToZones((TGeoTrap *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoSphere *>(shape) != nullptr) {
      zonePair = sphereToZones((TGeoSphere *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoPgon *>(shape) != nullptr) {
      zonePair = polyhedraToZones((TGeoPgon *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoPcon *>(shape) != nullptr) {
      zonePair = polyconeToZones((TGeoPcon *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoEltu *>(shape) != nullptr) {
      zonePair = ellipticalTubeToZones((TGeoEltu *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoArb8 *>(shape) != nullptr) {
      zonePair = arb8ToZones((TGeoArb8 *)shape, currTransformation, scale);
   } else if (dynamic_cast<TGeoGtra *>(shape) != nullptr) { // TODO
      printf("TGeoGtra not implemented yet\n");
   } else if (dynamic_cast<TGeoHype *>(shape) != nullptr) { // TODO
      printf("TGeoHype not implemented yet\n");
   } else if (dynamic_cast<TGeoTorus *>(shape) != nullptr) { // TODO
      printf("TGeoTorus not implemented yet\n");
   } else if (dynamic_cast<TGeoParaboloid *>(shape) != nullptr) { // TODO
      printf("TGeoParaboloid not implemented yet\n");
   } else if (dynamic_cast<TGeoXtru *>(shape) != nullptr) { // TODO
      printf("TGeoXtru not implemented yet\n");
   } else if (dynamic_cast<TGeoHalfSpace *>(shape) != nullptr) { // TODO
      printf("TGeoHalfSpace not implemented yet\n");
   } else if (dynamic_cast<TGeoCompositeShape *>(shape) != nullptr) { // TODO
      printf("TGeoCompositeShape not implemented yet\n");
   } else if (dynamic_cast<TGeoBBox *>(shape) != nullptr) { // All of this shapes in subclasses of TGeoBBox
      zonePair = boxToZones((TGeoBBox *)shape, currTransformation, scale);
   } else {
      printf("ELSE");
   }
   return zonePair;
}

SGeoBody createCone(TGeoTranslation startPoint, TGeoTranslation vectorToEnd, double rStart, double rEnd)
{
   SGeoBody outBody;
   if (doubleEQ(rStart, rEnd)) {
      outBody.type          = 5;
      outBody.parameters[0] = 0;
      outBody.parameters[1] = 0;
      outBody.parameters[2] = 0;
      outBody.parameters[3] = 0;
      outBody.parameters[4] = 0;
      outBody.parameters[5] = 0;
      outBody.parameters[6] = rStart;
      addVectorToElement(outBody, 0, startPoint);
      addVectorToElement(outBody, 3, vectorToEnd);
   } else {
      if (rStart <= rEnd) startPoint = startPoint + vectorToEnd;
      if (rStart <= rEnd) vectorToEnd = vectorToEnd.Inverse();
      double rMax           = (rStart > rEnd) ? rStart : rEnd;
      double rMin           = (rStart > rEnd) ? rEnd : rStart;
      outBody.type          = 7;
      outBody.parameters[0] = 0;
      outBody.parameters[1] = 0;
      outBody.parameters[2] = 0;
      outBody.parameters[3] = 0;
      outBody.parameters[4] = 0;
      outBody.parameters[5] = 0;
      outBody.parameters[6] = rMax;
      outBody.parameters[7] = rMin;
      addVectorToElement(outBody, 0, startPoint);
      addVectorToElement(outBody, 3, vectorToEnd);
   }
   return outBody;
}
zoneList cutByPhi(TGeoHMatrix currTransformation, double sPhi, double dPhi, double halfX, double halfY, double halfZ)
{
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   TGeoTranslation currTranslation = TGeoTranslation(currTransformation);
   double          ePhi            = sPhi + dPhi;
   if (doubleGE(dPhi, 360)) return zoneList();
   TGeoHMatrix tmp = TGeoHMatrix();
   tmp.RotateZ(sPhi);
   TGeoHMatrix     innerRot      = currRotation * tmp;
   TGeoTranslation startFirstBox = currTranslation + TGeoTranslation(innerRot * TGeoTranslation(-halfX, 0, -halfZ));
   TGeoTranslation vec11         = TGeoTranslation(innerRot * TGeoTranslation(2 * halfX, 0, 0));
   TGeoTranslation vec12         = TGeoTranslation(innerRot * TGeoTranslation(0, -halfY, 0));
   TGeoTranslation vec13         = TGeoTranslation(innerRot * TGeoTranslation(0, 0, 2 * halfZ));
   SGeoBody        box1          = {3, {}};
   addVectorToElement(box1, 0, startFirstBox);
   addVectorToElement(box1, 3, vec11);
   addVectorToElement(box1, 6, vec12);
   addVectorToElement(box1, 9, vec13);
   tmp = TGeoHMatrix();
   tmp.RotateZ(ePhi);
   TGeoHMatrix     outerRot       = currRotation * tmp;
   TGeoTranslation startSecondBox = currTranslation + TGeoTranslation(outerRot * TGeoTranslation(halfX, 0, -halfZ));
   TGeoTranslation vec21          = TGeoTranslation(outerRot * TGeoTranslation(-2 * halfX, 0, 0));
   TGeoTranslation vec22          = TGeoTranslation(outerRot * TGeoTranslation(0, halfY, 0));
   TGeoTranslation vec23          = TGeoTranslation(outerRot * TGeoTranslation(0, 0, 2 * halfZ));
   SGeoBody        box2           = {3, {}};
   addVectorToElement(box2, 0, startSecondBox);
   addVectorToElement(box2, 3, vec21);
   addVectorToElement(box2, 6, vec22);
   addVectorToElement(box2, 9, vec23);
   zoneList outList;
   if (doubleLE(dPhi, 180)) {
      std::vector<zoneElement> vec_el;
      vec_el.push_back(std::make_pair(-1, box1));
      vec_el.push_back(std::make_pair(-1, box2));
      outList.push_back(vec_el);
   } else {
      outList.push_back(std::vector<zoneElement>(1, std::make_pair(-1, box1)));
      outList.push_back(std::vector<zoneElement>(1, std::make_pair(-1, box2)));
   }
   return outList;
}
inline std::pair<zoneList, zoneList> boxToZones(TGeoBBox *box, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double          x           = (double)box->GetDX() * scale;
   double          y           = (double)box->GetDY() * scale;
   double          z           = (double)box->GetDZ() * scale;
   TGeoTranslation startVector = TGeoTranslation(currTransformation * TGeoTranslation(-x, -y, -z));
   TGeoTranslation vec1        = TGeoTranslation(currRotation * TGeoTranslation(2 * x, 0, 0));
   TGeoTranslation vec2        = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * y, 0));
   TGeoTranslation vec3        = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
   SGeoBody        tmp         = {3, {}};
   addVectorToElement(tmp, 0, startVector);
   addVectorToElement(tmp, 3, vec1);
   addVectorToElement(tmp, 6, vec2);
   addVectorToElement(tmp, 9, vec3);
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, tmp)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> tubeSegToZones(TGeoTubeSeg *tube, TGeoHMatrix currTransformation, double scale)
{
   zoneList                      out, outerShell; // outer shell for conjunction at generating mother volume
   std::pair<zoneList, zoneList> tubeOut = tubeToZones(tube, currTransformation, scale);
   double                        sPhi    = (double)tube->GetPhi1(); // As I understand, in gedrees
   double                        dPhi    = (double)tube->GetPhi2() - sPhi;
   double                        rOut    = (double)tube->GetRmax() * scale;
   double                        z       = (double)tube->GetDz() * scale;
   zoneList                      phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, z);
   out                                   = andZone(tubeOut.first, phiZone);
   outerShell                            = andZone(tubeOut.second, phiZone);
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> tubeToZones(TGeoTube *tube, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double          rIn        = (double)tube->GetRmin() * scale;
   double          rOut       = (double)tube->GetRmax() * scale;
   double          z          = (double)tube->GetDz() * scale;
   TGeoTranslation startTube  = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -z));
   TGeoTranslation endTubeVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
   SGeoBody        innerTube  = {5, {}};
   addVectorToElement(innerTube, 0, startTube);
   addVectorToElement(innerTube, 3, endTubeVec);
   innerTube.parameters[6] = rIn;
   SGeoBody outerTube      = {5, {}};
   addVectorToElement(outerTube, 0, startTube);
   addVectorToElement(outerTube, 3, endTubeVec);
   outerTube.parameters[6] = rOut;
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, outerTube)));
   outerShell = out;
   if (doubleNE(rIn, 0)) {
      out = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, innerTube))));
   }
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> coneSegToZones(TGeoConeSeg *cone, TGeoHMatrix currTransformation, double scale)
{
   zoneList out, outerShell; // outer shell for conjunction at generating mother volume
   double   rMax1 = (double)cone->GetRmax1() * scale;
   double   rMax2 = (double)cone->GetRmax2() * scale;
   double   z     = (double)cone->GetDz() * scale;
   double   sPhi  = (double)cone->GetPhi1(); // As I understand, in gedrees
   double   dPhi  = (double)cone->GetPhi2() - sPhi;

   std::pair<zoneList, zoneList> coneOut = coneToZones(cone, currTransformation, scale);

   double   rOut    = (rMax1 > rMax2) ? rMax1 : rMax2;
   zoneList phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, z);
   out              = andZone(coneOut.first, phiZone);
   outerShell       = andZone(coneOut.second, phiZone);
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> coneToZones(TGeoCone *cone, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double rMin1 = (double)cone->GetRmin1() * scale;
   double rMin2 = (double)cone->GetRmin2() * scale;
   double rMax1 = (double)cone->GetRmax1() * scale;
   double rMax2 = (double)cone->GetRmax2() * scale;
   double z     = (double)cone->GetDz() * scale;

   TGeoTranslation startTubeVector = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -z));
   TGeoTranslation endTubeVec      = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
   SGeoBody        outerCone       = createCone(startTubeVector, endTubeVec, rMax1, rMax2);
   SGeoBody        innerCone       = createCone(startTubeVector, endTubeVec, rMin1, rMin2);

   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, outerCone)));
   outerShell = out;
   if (doubleNE(rMin1, 0) || doubleNE(rMin2, 0)) {
      out = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, innerCone))));
   }

   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trd1ToZones(TGeoTrd1 *trd, TGeoHMatrix currTransformation, double scale)
{
   zoneList        out, outerShell; // outer shell for conjunction at generating mother volume
   double          dz  = (double)trd->GetDz() * scale;
   double          dx1 = (double)trd->GetDx1() * scale;
   double          dy1 = (double)trd->GetDy() * scale;
   double          dx2 = (double)trd->GetDx2() * scale;
   double          dy2 = (double)trd->GetDy() * scale;
   TGeoTranslation v1  = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, -dy1, -dz));
   TGeoTranslation v2  = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, dy1, -dz));
   TGeoTranslation v3  = TGeoTranslation(currTransformation * TGeoTranslation(dx1, dy1, -dz));
   TGeoTranslation v4  = TGeoTranslation(currTransformation * TGeoTranslation(dx1, -dy1, -dz));
   TGeoTranslation v5  = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, -dy2, dz));
   TGeoTranslation v6  = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, dy2, dz));
   TGeoTranslation v7  = TGeoTranslation(currTransformation * TGeoTranslation(dx2, dy2, dz));
   TGeoTranslation v8  = TGeoTranslation(currTransformation * TGeoTranslation(dx2, -dy2, dz));
   SGeoBody        tmp = {2, {}};
   addVectorToElement(tmp, 0, v1);
   addVectorToElement(tmp, 3, v2);
   addVectorToElement(tmp, 6, v3);
   addVectorToElement(tmp, 9, v4);
   addVectorToElement(tmp, 12, v5);
   addVectorToElement(tmp, 15, v6);
   addVectorToElement(tmp, 18, v7);
   addVectorToElement(tmp, 21, v8);
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, tmp)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trd2ToZones(TGeoTrd2 *trd, TGeoHMatrix currTransformation, double scale)
{
   zoneList        out, outerShell; // outer shell for conjunction at generating mother volume
   double          dz  = (double)trd->GetDz() * scale;
   double          dx1 = (double)trd->GetDx1() * scale;
   double          dy1 = (double)trd->GetDy1() * scale;
   double          dx2 = (double)trd->GetDx2() * scale;
   double          dy2 = (double)trd->GetDy2() * scale;
   TGeoTranslation v1  = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, -dy1, -dz));
   TGeoTranslation v2  = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, dy1, -dz));
   TGeoTranslation v3  = TGeoTranslation(currTransformation * TGeoTranslation(dx1, dy1, -dz));
   TGeoTranslation v4  = TGeoTranslation(currTransformation * TGeoTranslation(dx1, -dy1, -dz));
   TGeoTranslation v5  = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, -dy2, dz));
   TGeoTranslation v6  = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, dy2, dz));
   TGeoTranslation v7  = TGeoTranslation(currTransformation * TGeoTranslation(dx2, dy2, dz));
   TGeoTranslation v8  = TGeoTranslation(currTransformation * TGeoTranslation(dx2, -dy2, dz));
   SGeoBody        tmp = {2, {}};
   addVectorToElement(tmp, 0, v1);
   addVectorToElement(tmp, 3, v2);
   addVectorToElement(tmp, 6, v3);
   addVectorToElement(tmp, 9, v4);
   addVectorToElement(tmp, 12, v5);
   addVectorToElement(tmp, 15, v6);
   addVectorToElement(tmp, 18, v7);
   addVectorToElement(tmp, 21, v8);
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, tmp)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trapToZones(TGeoTrap *trap, TGeoHMatrix currTransformation, double scale)
{
   zoneList out, outerShell; // outer shell for conjunction at generating mother volume
   double   dz     = (double)trap->GetDz() * scale;
   double   dy1    = (double)trap->GetH1() * scale;
   double   dy2    = (double)trap->GetH2() * scale;
   double   bl1    = (double)trap->GetBl1() * scale;
   double   tl1    = (double)trap->GetTl1() * scale;
   double   bl2    = (double)trap->GetBl2() * scale;
   double   tl2    = (double)trap->GetTl2() * scale;
   double   theta  = (double)trap->GetTheta() / (180 / M_PI);
   double   phi    = (double)trap->GetPhi() / (180 / M_PI);
   double   alpha1 = (double)trap->GetAlpha1() / (180 / M_PI);
   double   alpha2 = (double)trap->GetAlpha2() / (180 / M_PI);

   // double rCenter = dz / cos(theta);
   double xC = dz * tan(theta) * cos(phi); // rCenter * sin(theta) * cos(phi);
   double yC = dz * tan(theta) * sin(phi); // rCenter * sin(theta) * sin(phi);
   double x1 = dy1 * tan(alpha1);
   double x2 = dy2 * tan(alpha2);

   TGeoTranslation v1  = TGeoTranslation(currTransformation * TGeoTranslation(-xC + x1 - tl1, -yC + dy1, -dz));
   TGeoTranslation v2  = TGeoTranslation(currTransformation * TGeoTranslation(-xC + x1 + tl1, -yC + dy1, -dz));
   TGeoTranslation v3  = TGeoTranslation(currTransformation * TGeoTranslation(-xC - x1 + bl1, -yC - dy1, -dz));
   TGeoTranslation v4  = TGeoTranslation(currTransformation * TGeoTranslation(-xC - x1 - bl1, -yC - dy1, -dz));
   TGeoTranslation v5  = TGeoTranslation(currTransformation * TGeoTranslation(xC + x2 - tl2, yC + dy2, dz));
   TGeoTranslation v6  = TGeoTranslation(currTransformation * TGeoTranslation(xC + x2 + tl2, yC + dy2, dz));
   TGeoTranslation v7  = TGeoTranslation(currTransformation * TGeoTranslation(xC - x2 + bl2, yC - dy2, dz));
   TGeoTranslation v8  = TGeoTranslation(currTransformation * TGeoTranslation(xC - x2 - bl2, yC - dy2, dz));
   SGeoBody        tmp = {2, {}};
   addVectorToElement(tmp, 0, v1);
   addVectorToElement(tmp, 3, v2);
   addVectorToElement(tmp, 6, v3);
   addVectorToElement(tmp, 9, v4);
   addVectorToElement(tmp, 12, v5);
   addVectorToElement(tmp, 15, v6);
   addVectorToElement(tmp, 18, v7);
   addVectorToElement(tmp, 21, v8);
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, tmp)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> sphereToZones(TGeoSphere *sphere, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);

   double rIn  = (double)sphere->GetRmin() * scale;
   double rOut = (double)sphere->GetRmax() * scale;
   double sPhi = (double)sphere->GetPhi1();
   double dPhi = (double)sphere->GetPhi2() - sPhi;
   double sTheta =
      (double)sphere->GetTheta1() / (180 / M_PI); // convert to radians becouse cmath functions work with them
   double   dTheta           = (double)sphere->GetTheta2() / (180 / M_PI) - sTheta;
   double   eTheta           = sTheta + dTheta;
   SGeoBody innerSphere      = {0, {}};
   innerSphere.parameters[3] = rIn;
   SGeoBody outerSphere      = {0, {}};
   outerSphere.parameters[3] = rOut;
   {
      TGeoTranslation currTranslation = TGeoTranslation(currTransformation);
      addVectorToElement(innerSphere, 0, currTranslation);
      addVectorToElement(outerSphere, 0, currTranslation);
   }
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, outerSphere)));
   outerShell = out;
   if (doubleNE(rIn, 0)) {
      out = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, innerSphere))));
   }
   zoneList phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, rOut);
   out              = andZone(out, phiZone);
   outerShell       = andZone(outerShell, phiZone);

   if (doubleNE(sTheta, 0) && doubleNE(sTheta, M_PI) && doubleNE(sTheta, M_PI / 2.0)) {
      TGeoTranslation trc1Start = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -rOut * cos(sTheta)));
      TGeoTranslation trc1Vec   = TGeoTranslation(currRotation * TGeoTranslation(0, 0, (rOut - rIn) * cos(sTheta)));
      SGeoBody        trc1      = {7, {}};
      addVectorToElement(trc1, 0, trc1Start);
      addVectorToElement(trc1, 3, trc1Start);
      trc1.parameters[6] = (double)(rOut * tan(sTheta));
      trc1.parameters[7] = (double)(rIn * sin(sTheta));
      int mul            = (sTheta < (M_PI / 2.0)) ? -1 : 1;
      out                = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(mul, trc1))));
      outerShell         = andZone(outerShell, zoneList(1, std::vector<zoneElement>(1, std::make_pair(mul, trc1))));
   } else if (sTheta == (M_PI / 2.0)) {
      TGeoTranslation db1s  = TGeoTranslation(currTransformation * TGeoTranslation(-rOut, -rOut, 0));
      TGeoTranslation db1v1 = TGeoTranslation(currRotation * TGeoTranslation(2 * rOut, 0, 0));
      TGeoTranslation db1v2 = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * rOut, 0));
      TGeoTranslation db1v3 = TGeoTranslation(currRotation * TGeoTranslation(0, 0, -rOut));
      SGeoBody        db1   = {3, {}};
      addVectorToElement(db1, 0, db1s);
      addVectorToElement(db1, 3, db1v1);
      addVectorToElement(db1, 6, db1v2);
      addVectorToElement(db1, 9, db1v3);
      out        = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, db1))));
      outerShell = andZone(outerShell, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, db1))));
   }

   if (doubleNE(eTheta, 0) && doubleNE(eTheta, M_PI) && doubleNE(eTheta, M_PI / 2.0)) {
      TGeoTranslation trc2Start = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -rOut * cos(eTheta)));
      TGeoTranslation trc2Vec   = TGeoTranslation(currRotation * TGeoTranslation(0, 0, (rOut - rIn) * cos(eTheta)));
      SGeoBody        trc2      = {7, {}};
      addVectorToElement(trc2, 0, trc2Start);
      addVectorToElement(trc2, 3, trc2Start);
      trc2.parameters[6] = (double)(rOut * tan(eTheta));
      trc2.parameters[7] = (double)(rIn * sin(eTheta));
      int mul            = (eTheta > (M_PI / 2.0)) ? -1 : 1;
      out                = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(mul, trc2))));
      outerShell         = andZone(outerShell, zoneList(1, std::vector<zoneElement>(1, std::make_pair(mul, trc2))));
   } else if (eTheta == (M_PI / 2.0)) {
      TGeoTranslation db2s  = TGeoTranslation(currTransformation * TGeoTranslation(-rOut, -rOut, 0));
      TGeoTranslation db2v1 = TGeoTranslation(currRotation * TGeoTranslation(2 * rOut, 0, 0));
      TGeoTranslation db2v2 = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * rOut, 0));
      TGeoTranslation db2v3 = TGeoTranslation(currRotation * TGeoTranslation(0, 0, rOut));
      SGeoBody        db2   = {3, {}};
      addVectorToElement(db2, 0, db2s);
      addVectorToElement(db2, 3, db2v1);
      addVectorToElement(db2, 6, db2v2);
      addVectorToElement(db2, 9, db2v3);
      out        = andZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, db2))));
      outerShell = andZone(outerShell, zoneList(1, std::vector<zoneElement>(1, std::make_pair(-1, db2))));
   }
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyconeToZones(TGeoPcon *polycone, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double sPhi = (double)polycone->GetPhi1(); // In Pcon angle converted to degrees, and TMatrixes work with degrees.
                                              // (?)
   double          dPhi = (double)polycone->GetDphi();
   int             nz   = polycone->GetNz();
   double          zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev;
   double         *zValues    = polycone->GetZ();
   double         *rMinValues = polycone->GetRmin();
   double         *rMaxValues = polycone->GetRmax();
   double          zMin = zValues[0] * scale, zMax = zValues[nz - 1] * scale, rMaxMax = rMaxValues[0] * scale;
   TGeoTranslation startPoint, endVec;
   zoneElement     currentOuterCone, currentInnerCone;
   double          zCenter = (zMax + zMin) / 2.0;
   for (int i = 1; i < nz; ++i) {
      zPrev    = zCenter + (zValues[i - 1] - zCenter) * scale;
      zCurrent = zCenter + (zValues[i] - zCenter) * scale;
      if (zPrev == zCurrent) continue;
      rMaxPrev    = rMaxValues[i - 1] * scale;
      rMaxCurrent = rMaxValues[i] * scale;
      rMinPrev    = rMinValues[i - 1] * scale;
      rMinCurrent = rMinValues[i] * scale;
      startPoint  = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, zPrev));
      endVec      = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zCurrent - zPrev));
      rMaxMax     = (rMaxCurrent > rMaxMax) ? rMaxCurrent : rMaxMax;
      if (doubleNE(rMaxPrev, 0) || doubleNE(rMaxCurrent, 0)) {
         currentOuterCone = std::make_pair(1, createCone(startPoint, endVec, rMaxPrev, rMaxCurrent));
         outerShell       = orZone(outerShell, zoneList(1, std::vector<zoneElement>(1, currentOuterCone)));
      }
      if (doubleNE(rMinPrev, 0) || doubleNE(rMinCurrent, 0)) {
         currentInnerCone = std::make_pair(1, createCone(startPoint, endVec, rMinPrev, rMinCurrent));
         out              = orZone(out, zoneList(1, std::vector<zoneElement>(1, currentInnerCone)));
      }
   }
   out                           = andZone(outerShell, notZone(out));
   TGeoHMatrix phiTransformation = currTransformation * TGeoTranslation(0, 0, zMin + (zMax - zMin) / 2.0);
   zoneList    phiZone           = cutByPhi(phiTransformation, sPhi, dPhi, rMaxMax, rMaxMax, (zMax - zMin) / 2.0);
   //     startPoint = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zMin + (zMax - zMin) / 2.0));
   //     zoneList phiZone = cutByPhi(startPoint, currRotation, sPhi, dPhi, rMaxMax, rMaxMax, (zMax - zMin) / 2.0);
   out        = andZone(out, phiZone);
   outerShell = andZone(outerShell, phiZone);
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyhedraToZones(TGeoPgon *polyhedra, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoRotation(currTransformation);
   double  sPhi = (double)polyhedra->GetPhi1(); // In Pcon angle converted to degrees, and TMatrixes work with degrees.
   double  dPhi = (double)polyhedra->GetDphi();
   int     numSides = polyhedra->GetNedges();
   int     nz       = polyhedra->GetNz();
   double  da       = dPhi / numSides;
   double  zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev;
   double *zValues    = polyhedra->GetZ();
   double *rMinValues = polyhedra->GetRmin();
   double *rMaxValues = polyhedra->GetRmax();
   double  zMin = zValues[0] * scale, zMax = zValues[nz - 1] * scale;
   TGeoTranslation startZPoint, endZPoint, vert0, vert1, vert2, vert3, vecR, vecL;
   TGeoHMatrix     phiRotation;
   SGeoBody        arb;
   double          zCenter = (zMax + zMin) / 2.0;
   for (int i = 1; i < nz; ++i) {
      zPrev    = zCenter + (zValues[i - 1] - zCenter) * scale;
      zCurrent = zCenter + (zValues[i] - zCenter) * scale;
      if (doubleEQ(zPrev, zCurrent)) continue;
      rMaxPrev    = rMaxValues[i - 1] * scale / cos(da / 2 * M_PI / 180); // Convert angle to radians
      rMaxCurrent = rMaxValues[i] * scale / cos(da / 2 * M_PI / 180);
      rMinPrev    = rMinValues[i - 1] * scale / cos(da / 2 * M_PI / 180);
      rMinCurrent = rMinValues[i] * scale / cos(da / 2 * M_PI / 180);
      startZPoint = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, zPrev));
      endZPoint   = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, zCurrent));
      for (int k = 0; k < numSides; ++k) {
         phiRotation = TGeoHMatrix();
         phiRotation.RotateZ(sPhi + da * k);
         phiRotation = currRotation * phiRotation;
         vecR        = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); // Vector to right vertex of secctor
         phiRotation = TGeoHMatrix();
         phiRotation.RotateZ(sPhi + da * (k + 1));
         phiRotation = currRotation * phiRotation;
         vecL        = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); // Vector to left vertex of secctor
         if (doubleNE(rMaxPrev, 0) || doubleNE(rMaxCurrent, 0)) {               // Work with external polyhedra
            vert0    = startZPoint + vecR * rMaxPrev;                           // Right bottom vertex
            vert1    = startZPoint + vecL * rMaxPrev;                           // Left bottom vertex
            vert2    = endZPoint + vecL * rMaxCurrent;                          // Left upper vertex
            vert3    = endZPoint + vecR * rMaxCurrent;                          // Right upper vertex
            arb.type = 2;                                                       // arb is ARB
            addVectorToElement(arb, 0, vert0);
            addVectorToElement(arb, 3, vert1);
            addVectorToElement(arb, 6, vert2);
            addVectorToElement(arb, 9, vert3);
            addVectorToElement(arb, 12, startZPoint);
            addVectorToElement(arb, 15, startZPoint);
            addVectorToElement(arb, 18, endZPoint);
            addVectorToElement(arb, 21, endZPoint);
            outerShell = orZone(outerShell, zoneList(1, std::vector<zoneElement>(1, std::make_pair(1, arb))));
         }
         if (doubleNE(rMinPrev, 0) || doubleNE(rMinCurrent, 0)) { // Work with internal polyhedra
            vert0    = startZPoint + vecR * rMinPrev;             // Right bottom vertex
            vert1    = startZPoint + vecL * rMinPrev;             // Left bottom vertex
            vert2    = endZPoint + vecL * rMinCurrent;            // Left upper vertex
            vert3    = endZPoint + vecR * rMinCurrent;            // Right upper vertex
            arb.type = 2;                                         // wed is ARB
            addVectorToElement(arb, 0, vert0);
            addVectorToElement(arb, 3, vert1);
            addVectorToElement(arb, 6, vert2);
            addVectorToElement(arb, 9, vert3);
            addVectorToElement(arb, 12, startZPoint);
            addVectorToElement(arb, 15, startZPoint);
            addVectorToElement(arb, 18, endZPoint);
            addVectorToElement(arb, 21, endZPoint);
            out = orZone(out, zoneList(1, std::vector<zoneElement>(1, std::make_pair(1, arb))));
         }
      }
   }
   out = andZone(outerShell, notZone(out));
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> ellipticalTubeToZones(TGeoEltu *tube, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double          dx         = (double)tube->GetA() * scale;
   double          dy         = (double)tube->GetB() * scale;
   double          dz         = (double)tube->GetDz() * scale;
   TGeoTranslation startPoint = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -dz));
   TGeoTranslation endVec     = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * dz));
   TGeoTranslation xAxis      = TGeoTranslation(currRotation * TGeoTranslation(dx, 0, 0));
   TGeoTranslation yAxis      = TGeoTranslation(currRotation * TGeoTranslation(0, dy, 0));
   SGeoBody        rec        = {6, {}};
   addVectorToElement(rec, 0, startPoint);
   addVectorToElement(rec, 3, endVec);
   addVectorToElement(rec, 6, xAxis);
   addVectorToElement(rec, 9, yAxis);
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, rec)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> arb8ToZones(TGeoArb8 *gtrap, TGeoHMatrix currTransformation, double scale)
{
   zoneList    out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoHMatrix currRotation = TGeoHMatrix(currTransformation);
   currRotation.SetTranslation(kNullVector);
   double    dz = (double)gtrap->GetDz() * scale;
   Double_t  vertexes[8][2];
   Double_t *tmpPointer = gtrap->GetVertices();
   memcpy(&vertexes[0][0], tmpPointer, 8 * 2 * sizeof(Double_t));
   TGeoTranslation tmp;
   double          cx, cy, cz;
   SGeoBody        trap = {2, {}};
   for (int i = 0; i < 8; i++) {
      cx  = vertexes[i][0] * scale;
      cy  = vertexes[i][1] * scale;
      cz  = (i < 4) ? -dz : dz;
      tmp = TGeoTranslation(currTransformation * TGeoTranslation(cx, cy, cz));
      addVectorToElement(trap, i * 3, tmp);
   }
   out.push_back(std::vector<zoneElement>(1, std::make_pair(1, trap)));
   outerShell = out;
   return std::make_pair(out, outerShell);
}

std::pair<zoneList, zoneList> compositeShapeToZones(TGeoCompositeShape *shape, TGeoHMatrix currTransformation,
                                                    double scale)
{
   zoneList                      out, outerShell; // outer shell for conjunction at generating mother volume
   TGeoBoolNode                 *boolNode = shape->GetBoolNode();
   TGeoMatrix                   *lm       = boolNode->GetLeftMatrix();
   TGeoShape                    *ls       = boolNode->GetLeftShape();
   TGeoHMatrix                   ln       = (*lm) * currTransformation;
   std::pair<zoneList, zoneList> lo       = getZoneFromShape(ls, ln, scale);
   TGeoMatrix                   *rm       = boolNode->GetRightMatrix();
   TGeoShape                    *rs       = boolNode->GetRightShape();
   TGeoHMatrix                   rn       = (*rm) * currTransformation;
   std::pair<zoneList, zoneList> ro       = getZoneFromShape(rs, rn, scale);
   switch (boolNode->GetBooleanOperator()) {
   case TGeoBoolNode::kGeoUnion:
      return std::pair<zoneList, zoneList>(orZone(lo.first, ro.first), orZone(lo.second, ro.second));
      break;
   case TGeoBoolNode::kGeoIntersection:
      return std::pair<zoneList, zoneList>(andZone(lo.first, ro.first), orZone(lo.second, ro.second));
      break;
   case TGeoBoolNode::kGeoSubtraction:
      return std::pair<zoneList, zoneList>(andZone(lo.first, ro.first), notZone(orZone(lo.second, ro.second)));
      break;
   default:
      printf("Something is wrong with Node. Returning subtraction");
      return std::pair<zoneList, zoneList>(andZone(lo.first, ro.first), notZone(orZone(lo.second, ro.second)));
   }
}

MediumData getMedium(TGeoNode *node)
{
   double        g      = 1;
   double        cm3    = getDefaultScale() * getDefaultScale() * getDefaultScale();
   MediumData    out    = {0, 0, 0, {}};
   TGeoMaterial *medium = node->GetMedium()->GetMaterial();
   out.nChemEl          = medium->GetNelements();
   out.nType            = (out.nChemEl == 1) ? 1 : 2;
   out.Rho              = medium->GetDensity() / (g / cm3);
   if (out.Rho == 0) {
      printf("Is medium a Vacuum?");
      out.nType = 1000;
      return out;
   }
   for (int i = 0; i < out.nChemEl; ++i) {
      out.Elements[i].Nuclid  = medium->GetElement(i)->Z();
      out.Elements[i].A       = medium->GetElement(i)->N();
      out.Elements[i].Z       = medium->GetElement(i)->Z();
      out.Elements[i].Density = (medium->IsMixture()) ? out.Rho * (((TGeoMixture *)medium)->GetWmixt())[i] : out.Rho;
      out.Elements[i].Conc    = out.Elements[i].Density / medium->GetElement(i)->A() * TMath::Na() * 1E-27;
      //         out.Elements[i].Conc = mat->GetVecNbOfAtomsPerVolume()[i] * cm3 * 1E-27; //Correct
      //         out.Elements[i].ionEv = mat->GetElement(i)->GetIonisation()->GetMeanExcitationEnergy() / eV; //As I
      //         unserstand, it is, but it can fill automatically.
      for (int k = 0; k < i; ++k) {
         if (out.nType != 3 && out.Elements[i].Z == out.Elements[k].Z) out.nType = 3;
      }
   }
   return out;
}
} // namespace tgeanttoshield
