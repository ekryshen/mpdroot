//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>
#include <iostream>
#include <fstream>

#include <TGeoManager.h>
#include <TRandom2.h>
#include <TRotation.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TGeoNode.h>

#include "FairRootManager.h"
#include "FairLogger.h"

#include "MpdTof.h"
#include "MpdTofHitProducerQA.h"
#include "MpdTofPoint.h"
#include "MpdTofGeoUtils.h"

using namespace std;

MpdTofGeoUtils *MpdTofGeoUtils::instance = nullptr;
//------------------------------------------------------------------------------------------------------------------------
void MpdTofGeoUtils::FindNeighborStrips(Double_t thresh, MpdTofHitProducerQA *pQA, bool forced)
{
   static bool found = false;
   if (!forced && found) return; // neighbor strips already found and filled

   size_t NR = 0, NL = 0;
   for (auto &it1 : mStrips) // cycle1 by strips
   {
      LStrip *strip1 = &(it1.second);

      for (const auto &it2 : mStrips) // cycle2 by strips
      {
         const LStrip *strip2 = &(it2.second);

         if (strip1 == strip2) continue; // skip calculate yourself

         // CATION: Ckeck  only left and right sides(one row strips NOW)
         double distance = strip1->Distance(LStrip::kRight, *strip2);

         if (pQA) pQA->FillDistanceToStrip(distance);
         if (distance < thresh && MpdTofPoint::IsSameGap(strip1->volumeUID, strip2->volumeUID)) {
            strip1->neighboring[LStrip::kRight] = strip2->volumeUID;
            NR++;
            if (pQA) pQA->FillNeighborPairs(strip2->stripID, strip1->stripID);
         }

         distance = strip1->Distance(LStrip::kLeft, *strip2);

         if (pQA) pQA->FillDistanceToStrip(distance * -1.);
         if (distance < thresh && MpdTofPoint::IsSameGap(strip1->volumeUID, strip2->volumeUID)) {
            strip1->neighboring[LStrip::kLeft] = strip2->volumeUID;
            NL++;
            if (pQA) pQA->FillNeighborPairs(strip2->stripID, strip1->stripID);
         }

      } // cycle2 by strips

   } // cycle1 by strips

   found = true;

   LOG(debug1) << "[MpdTofGeoUtils::FindNeighborStrips] Neighbor strips: left = " << NL << ", right = " << NR << ".";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofGeoUtils::ParseTGeoManager(MpdTofHitProducerQA *pQA, bool forced, const char *flnm)
{
   constexpr double degree180 = TMath::Pi();
   assert(gGeoManager);

   if (!forced && !mStrips.empty()) return; // already parsed and filled

   mStrips.clear();

   ofstream ofs;
   if (flnm) ofs.open(flnm);

   const TString pathTOF = "/cave_1/tof1_0";
   gGeoManager->cd(pathTOF);

   Double_t *local = new Double_t[3](), master[3] = {0}, dX, dY, dZ;
   int       nSectors = 0, nDetectors = 0, nStrips = 0;

   TObjArray *array = gGeoManager->GetCurrentVolume()->GetNodes();
   auto       it1   = array->MakeIterator();

   vector<Tinterval> vDetectorsZ, vDetectorsPhi;
   vector<Tinterval> vStripsZ, vStripsPhi;

   while (TGeoNode *sectorNode =
             (TGeoNode *)it1->Next()) // SECTORS vSector_#, DetectorBoxV_# +, StripBoxV_#, StripActiveGasV_# +
   {
      if (!TString(sectorNode->GetName()).Contains("tof1Sector")) continue;

      TString PATH1    = pathTOF + "/" + sectorNode->GetName();
      Int_t   sectorID = sectorNode->GetNumber(); // sector [1,...,14]
      nSectors++;

      if (ofs.is_open())
         ofs << "\n SECTOR: " << sectorNode->GetName() << ", copy# " << sectorID << " path= " << PATH1.Data();

      auto it3 = sectorNode->GetNodes()->MakeIterator();
      while (TGeoNode *detectorNode = (TGeoNode *)it3->Next()) // DETECTORS
      {
         TString detName = detectorNode->GetName();
         if (!detName.Contains("DetectorBox")) continue;

         TString PATH3      = PATH1 + "/" + detName;
         Int_t   detectorID = detectorNode->GetNumber(); // detector [1,...,20]
         nDetectors++;

         if (ofs.is_open())
            ofs << "\n\t\t DETECTOR: " << detName.Data() << ", copy# " << detectorID << " path= " << PATH3.Data();

         gGeoManager->cd(PATH3); // cd to detector box node

         TGeoMatrix *matrix = gGeoManager->GetCurrentMatrix();

         TGeoBBox *box = (TGeoBBox *)gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
         dX            = box->GetDX();
         dY            = box->GetDY();
         dZ            = box->GetDZ();

         TVector3 A, B, C, D; // detector edges
         auto     Z_range = make_pair(1.e+10, -1.e+10), Phi_range(Z_range);

         local[0] = -dX;
         local[1] = -dY + 1.3 - 0.26;
         local[2] = +dZ; // shift by Y to 1th gas sensitive layer of strips (for rough track propagate to plate)
         localToMaster(matrix, local, A, Z_range, Phi_range);

         local[0] = +dX;
         local[1] = -dY + 1.3 - 0.26;
         local[2] = +dZ;
         localToMaster(matrix, local, B, Z_range, Phi_range);

         local[0] = +dX;
         local[1] = -dY + 1.3 - 0.26;
         local[2] = -dZ;
         localToMaster(matrix, local, C, Z_range, Phi_range);

         local[0] = -dX;
         local[1] = -dY + 1.3 - 0.26;
         local[2] = -dZ;
         localToMaster(matrix, local, D, Z_range, Phi_range);

         if (pQA) pQA->FillDetectorsMap(A, B, C, D);

         // Add Detector location interval to intervalTree
         Int_t       duid = MpdTofPoint::ClearGap(MpdTofPoint::GetSuid72(
                  sectorID, detectorID, 0)); // suid for 0th strip of the selected detector, gap reset to 0
         LRectangle *det  = new LRectangle(duid, A, B, C, D);
         det->InitCenterPerp();

         vDetectorsZ.push_back(Tinterval(Z_range.first, Z_range.second, det));

         if (Phi_range.second - Phi_range.first > degree180) // segment overlap the boundary [-180.,+180[
         {
            vDetectorsPhi.push_back(Tinterval(-degree180, Phi_range.first, det));
            vDetectorsPhi.push_back(Tinterval(Phi_range.second, degree180, det));
            cout << "\n OVERLAPED   phimin=" << Phi_range.first << ", phimax=" << Phi_range.second;
         } else
            vDetectorsPhi.push_back(Tinterval(Phi_range.first, Phi_range.second, det));

         auto it4 = detectorNode->GetNodes()->MakeIterator();
         while (TGeoNode *stripBoxNode = (TGeoNode *)it4->Next()) // STRIPS
         {
            TString stripBoxName = stripBoxNode->GetName();
            if (!stripBoxName.Contains("StripActiveGas")) continue;

            TString PATH4   = PATH3 + "/" + stripBoxName;
            Int_t   stripID = stripBoxNode->GetNumber(); // strip[1,...,72]
            nStrips++;

            if (ofs.is_open())
               ofs << "\n\t\t\t STRIP: " << stripBoxName.Data() << ", copy# " << stripID << " path= " << PATH4.Data();

            gGeoManager->cd(PATH4); // cd to strip sensitive gas node

            TGeoMatrix *matrix4 = gGeoManager->GetCurrentMatrix();

            TGeoBBox *box4 = (TGeoBBox *)gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
            dX             = box4->GetDX();
            dY             = box4->GetDY();
            dZ             = box4->GetDZ();

            Int_t suid = MpdTofPoint::GetSuid72(sectorID, detectorID, stripID);

            if (ofs.is_open()) {
               MpdTofPoint::PrintSuid(suid, "\n\t\t\t\t Add strip:", ofs);
               ofs << " center= " << master[0] << ", " << master[1] << ", " << master[2] << ", dXYZ= " << dX << ", "
                   << dY << ", " << dZ;
            }

            LStrip strip(suid, sectorID, detectorID, stripID);
            strip.fMatrix = (*matrix4);

            auto sZ_range = make_pair(1.e+10, -1.e+10), sPhi_range(sZ_range); // init by big values
                                                                              //
                                                                              //   ^ +Z
            //   |     A             right side               B
            //   |       |----------------------------------|
            //   |       |-------------strip----------------|
            //   |     D             left side                C
            //   |-------------------------------------------------> +X
            //   | -Z

            // edges on the front plate of the strips. perp along Y.
            local[0] = -dX;
            local[1] = -dY;
            local[2] = +dZ;
            localToMaster(matrix4, local, strip.A, sZ_range, sPhi_range);

            local[0] = +dX;
            local[1] = -dY;
            local[2] = +dZ;
            localToMaster(matrix4, local, strip.B, sZ_range, sPhi_range);

            local[0] = +dX;
            local[1] = -dY;
            local[2] = -dZ;
            localToMaster(matrix4, local, strip.C, sZ_range, sPhi_range);

            local[0] = -dX;
            local[1] = -dY;
            local[2] = -dZ;
            localToMaster(matrix4, local, strip.D, sZ_range, sPhi_range);

            strip.InitCenterPerp(); // calc. strip center and perp. AFTER setup A,B,C,D

            if (pQA) pQA->FillStripsMap(strip.A, strip.B, strip.C, strip.D);

            auto pStrip = new LStrip(strip);

            vStripsZ.push_back(Tinterval(sZ_range.first, sZ_range.second, pStrip));

            if (sPhi_range.second - sPhi_range.first > degree180) // segment overlap the boundary [-180.,+180[
            {
               vStripsPhi.push_back(Tinterval(-degree180, sPhi_range.first, pStrip));
               vStripsPhi.push_back(Tinterval(sPhi_range.second, degree180, pStrip));
               if (ofs.is_open())
                  ofs << "\n OVERLAPED strip  phimin=" << sPhi_range.first << ", phimax=" << sPhi_range.second;
            } else
               vStripsPhi.push_back(Tinterval(sPhi_range.first, sPhi_range.second, pStrip));
            /// strip.Dump("\n strip:");
            bool IsUniqueUID = mStrips.insert(make_pair(suid, strip)).second;
            assert(IsUniqueUID);
         } // STRIPS

      } // DETECTORS
   }    // SECTORS

   mStripsZ   = TintervalTree(vStripsZ);
   mStripsPhi = TintervalTree(vStripsPhi);

   mDetectorsZ   = TintervalTree(vDetectorsZ);
   mDetectorsPhi = TintervalTree(vDetectorsPhi);

   if (ofs.is_open()) ofs.close();

   delete[] local;

   LOG(debug1) << "[MpdTof::ParseTGeoManager] Sectors=" << nSectors << ", detectors=" << nDetectors
               << ", strips=" << nStrips << ".";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofGeoUtils::localToMaster(const TGeoMatrix *matrix, Double_t *local, TVector3 &position,
                                   pair<double, double> &Z, pair<double, double> &Phi) const
{
   Double_t master[3];
   matrix->LocalToMaster(local, master);
   position.SetXYZ(master[0], master[1], master[2]);

   Z.first  = std::min(position.Z(), Z.first);
   Z.second = std::max(position.Z(), Z.second);

   Phi.first  = std::min(position.Phi(), Phi.first);
   Phi.second = std::max(position.Phi(), Phi.second);
}
//------------------------------------------------------------------------------------------------------------------------
const LStrip *MpdTofGeoUtils::FindStrip(Int_t suid) const
{
   auto cit = mStrips.find(suid);
   assert(
      !(cit == mStrips.end() && MpdTofPoint::PrintSuid(suid, "\n [MpdTofGeoUtils::FindStrip] : unknown suid. ", cerr)));

   return &(cit->second);
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofGeoUtils::IsPointInsideStrips(const TVector3 &position, vector<Tinterval> &intersect, double Zerror,
                                         double PhiError)
{
   const double Z   = position.Z();
   const double Phi = position.Phi();

   vector<Tinterval> insideZ, insidePhi;
   mStripsZ.findOverlapping(Z - Zerror, Z + Zerror, insideZ);
   mStripsPhi.findOverlapping(Phi - PhiError, Phi + PhiError, insidePhi);

   if (!insideZ.empty() && !insidePhi.empty()) // have overlaped segments both Z and Phi
   {
      intersect.clear();

      // calc. intersection
      sort(insideZ.begin(), insideZ.end(), IntervalValueSorter<LRectangle *>());
      sort(insidePhi.begin(), insidePhi.end(), IntervalValueSorter<LRectangle *>());
      set_intersection(insideZ.begin(), insideZ.end(), insidePhi.begin(), insidePhi.end(), back_inserter(intersect),
                       IntervalValueSorter<LRectangle *>());

      return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdTofGeoUtils::IsPointInsideDetectors(const TVector3 &position, vector<Tinterval> &intersect, double Zerror,
                                            double PhiError)
{
   const double Z   = position.Z();
   const double Phi = position.Phi();

   vector<Tinterval> insideZ, insidePhi;
   mDetectorsZ.findOverlapping(Z - Zerror, Z + Zerror, insideZ);
   mDetectorsPhi.findOverlapping(Phi - PhiError, Phi + PhiError, insidePhi);

   if (!insideZ.empty() && !insidePhi.empty()) // have overlaped segments both Z and Phi
   {
      intersect.clear();

      // calc. intersection
      sort(insideZ.begin(), insideZ.end(), IntervalValueSorter<LRectangle *>());
      sort(insidePhi.begin(), insidePhi.end(), IntervalValueSorter<LRectangle *>());
      set_intersection(insideZ.begin(), insideZ.end(), insidePhi.begin(), insidePhi.end(), back_inserter(intersect),
                       IntervalValueSorter<LRectangle *>());

      return true;
   }

   return false;
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
LRectangle::LRectangle(Int_t uid, const TVector3 &a, const TVector3 &b, const TVector3 &c, const TVector3 &d,
                       bool check)
   : IsInvalid(false), volumeUID(uid), A(a), B(b), C(c), D(d)
{
   if (check) CheckInValid();
}
//------------------------------------------------------------------------------------------------------------------------
void LRectangle::GetRPhiRanges(Double_t &Rmin, Double_t &Rmax, Double_t &PhiMin, Double_t &PhiMax)
{
   Rmin = PhiMin = 1.e+10;
   Rmax = PhiMax = -1.e+10;

   getRPhi(A, Rmin, Rmax, PhiMin, PhiMax); // edge
   getRPhi(B, Rmin, Rmax, PhiMin, PhiMax);
   getRPhi(C, Rmin, Rmax, PhiMin, PhiMax);
   getRPhi(D, Rmin, Rmax, PhiMin, PhiMax);

   getRPhi((A + B) * 0.5, Rmin, Rmax, PhiMin, PhiMax); // middle of side
   getRPhi((B + C) * 0.5, Rmin, Rmax, PhiMin, PhiMax);
   getRPhi((C + D) * 0.5, Rmin, Rmax, PhiMin, PhiMax);
   getRPhi((D + A) * 0.5, Rmin, Rmax, PhiMin, PhiMax);
}
//------------------------------------------------------------------------------------------------------------------------
Double_t LRectangle::DistanceFromPointToLine(const TVector3 *pos, const TVector3 &P1, const TVector3 &P2) const
{
   assert(P1 != P2);

   return ((*pos - P1).Cross(*pos - P2)).Mag() / (P2 - P1).Mag();
}
//------------------------------------------------------------------------------------------------------------------------
void LRectangle::Rotate(const TRotation &rot)
{
   A *= rot;
   B *= rot;
   C *= rot;
   D *= rot;
   center *= rot;
   perp *= rot;
}
//------------------------------------------------------------------------------------------------------------------------
Double_t LRectangle::DistanceFromPointToLineSegment(const TVector3 *pos, const TVector3 &P1, const TVector3 &P2) const
{
   assert(P1 != P2);

   TVector3 v = P2 - P1;
   TVector3 w = (*pos) - P1;

   double c1 = w.Dot(v);
   if (c1 <= 0) return w.Mag();

   double c2 = v.Dot(v);
   if (c2 <= c1) return ((*pos) - P2).Mag();

   TVector3 Pb = P1 + (c1 / c2) * v;
   return ((*pos) - Pb).Mag();
}
//------------------------------------------------------------------------------------------------------------------------
Double_t LRectangle::MinDistanceToEdge(const TVector3 *pos, Side_t &side) const
{
   double right = DistanceFromPointToLineSegment(pos, A, B);
   double left  = DistanceFromPointToLineSegment(pos, C, D);

   // sorting & return minimal value
   if (right <= left) {
      side = LStrip::kRight;
      return right;
   }

   side = LStrip::kLeft;
   return left;
}
//------------------------------------------------------------------------------------------------------------------------
void LRectangle::Dump(const char *comment, ostream &out) const
{
   if (comment) out << comment;
   out << " uid=" << volumeUID << " IsInvalid=" << IsInvalid;
   MpdTof::Print(A, " A:", out);
   MpdTof::Print(B, " B:", out);
   MpdTof::Print(C, " C:", out);
   MpdTof::Print(D, " D:", out);
   MpdTof::Print(perp, " perp:", out);
}
//------------------------------------------------------------------------------------------------------------------------
void LRectangle::Test(const char *comment, ostream &out)
{
   if (comment) out << comment;
   out << "\n [LRectangle::Test]--------------------------->>>";
   // const double  epsilon = 1.E-3;

   for (double xshift = -10; xshift < 10; xshift++)
      for (double yshift = -10; yshift < 10; yshift++) {
         LRectangle A(1, TVector3(0, 0, 0), TVector3(0, 10, 0), TVector3(10, 10, 0), TVector3(10, 0, 0));
         A.InitCenterPerp(); // 10x10
         LRectangle B(2, TVector3(0, 0, 1), TVector3(0, 5, 1), TVector3(5, 5, 1), TVector3(5, 0, 1));
         B.InitCenterPerp(); // 5x5

         B.Shift({xshift, yshift, 0.});

         double square = B.GetIntersectionSquare(A); // smaller intersect biger
         out << "\n shift: (" << xshift << "," << yshift << ") square=" << square;

         // Randomize rectangle orientation
         TRotation rot;
         TVector3  axis, shift;
         for (auto i = 0; i < 100; i++) {
            axis.SetXYZ(gRandom->Uniform(), gRandom->Uniform(), gRandom->Uniform());
            shift.SetXYZ(gRandom->Uniform(-10, 10), gRandom->Uniform(-10, 10), gRandom->Uniform(-10, 10));
            rot.Rotate(gRandom->Uniform(-10., 10.), axis); // random rotate around axis

            A.Rotate(rot);
            B.Rotate(rot);
            A.Shift(shift);
            B.Shift(shift);

            double newsquare = B.GetIntersectionSquare(A);
            if (fabs(newsquare - square) > 1.)
               out << "\n\tperp.: theta=" << A.perp.Theta() * TMath::RadToDeg()
                   << " phi=" << A.perp.Phi() * TMath::RadToDeg() << ", int. square=" << newsquare;
         }
      }
}
//------------------------------------------------------------------------------------------------------------------------
Double_t LRectangle::GetIntersectionSquare(const LRectangle &v1) const
{
   LRectangle v(v1);
   if (!IsSamePlane(v1).first) {
      v = ProjectionToPlane(v1);
   }

   if (!IsSamePlane(v).first) return 0.;

   TVector3     p, AB = B - A, AD = D - A;
   const size_t Nprobes  = 10000;
   size_t       Ncounter = 0;

   for (size_t i = 0; i < Nprobes; i++) {
      p = A + AB * gRandom->Uniform() + AD * gRandom->Uniform();
      if (v.IsPointInside(p)) Ncounter++;
   }

   if (Ncounter == 0) return 0.;

   return Square() * Ncounter / ((double)Nprobes);
}
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
LStrip::LStrip() : LRectangle(), sectorID(kInvalid), detectorID(kInvalid), stripID(kInvalid)
{
   neighboring[kRight] = kInvalid;
   neighboring[kLeft]  = kInvalid;
}
//------------------------------------------------------------------------------------------------------------------------
LStrip::LStrip(Int_t uid, Int_t sector, Int_t detector, Int_t strip)
   : LRectangle(), sectorID(sector), detectorID(detector), stripID(strip)
{
   volumeUID           = uid;
   neighboring[kRight] = kInvalid;
   neighboring[kLeft]  = kInvalid;
}
//------------------------------------------------------------------------------------------------------------------------
void LStrip::Dump(const char *comment, ostream &out) const
{
   if (comment) out << comment;
   out << "  ids: " << sectorID << ", " << detectorID << ", " << stripID;

   LRectangle::Dump(nullptr, out);
}
//------------------------------------------------------------------------------------------------------------------------
Double_t LStrip::Distance(Side_t side, const LStrip &strip) const
{
   Double_t value, min1 = 1.e+10, min2 = 1.e+10; // big value

   if ((*this) == strip) return min1;       // same strip, return big value
   if (!IsSameDetector(strip)) return min1; // different detector, return big value

   const TVector3 *p1{nullptr}, *p2{nullptr}; // TODO test
   switch (side) {
   case kRight:
      p1 = &A;
      p2 = &B;
      break;
   case kLeft:
      p1 = &C;
      p2 = &D;
      break;
   case kInvalid: return min1;
   }

   value = (*p1 - strip.A).Mag();
   min1  = std::min(value, min1);
   value = (*p1 - strip.B).Mag();
   min1  = std::min(value, min1);
   value = (*p1 - strip.C).Mag();
   min1  = std::min(value, min1);
   value = (*p1 - strip.D).Mag();
   min1  = std::min(value, min1);

   value = (*p2 - strip.A).Mag();
   min2  = std::min(value, min2);
   value = (*p2 - strip.B).Mag();
   min2  = std::min(value, min2);
   value = (*p2 - strip.C).Mag();
   min2  = std::min(value, min2);
   value = (*p2 - strip.D).Mag();
   min2  = std::min(value, min2);

   return std::min(min1, min2);
}
//------------------------------------------------------------------------------------------------------------------------
