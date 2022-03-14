// -------------------------------------------------------------------------
// -----                                    source file                -----
// -----                 Created 13/06/2014  by A. Zinchenko           -----
// -------------------------------------------------------------------------

/**  MpdVectorFinder.cxx
 *@author D.Zinchenko, A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD Inner Tracking System (ITS) using Vector Finder approach
 **/

#include "MpdStsGeoPar.h"
#include "MpdVectorFinder.h"
#include "MpdVector.h"
#include "MpdCodeTimer.h"
#include "MpdConstField.h"
#include "MpdItsHit5spd.h"
#include "MpdItsKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdMCTrack.h"
#include "MpdStsPoint.h"
//#include "MpdTpcKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"

#include "FairGeoNode.h"
#include "FairMCPoint.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoTube.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include <TRandom.h>
#include <TArrayI.h>

#include <iostream>
#include <map>
#include <set>

#include <string>

using std::cout;
using std::endl;
// using std::map;

const Double_t MpdVectorFinder::fgkChi2Cut = 10; // 20; //100;

static clock_t tStart  = 0;
static clock_t tFinish = 0;
static clock_t tAll    = 0;

FILE *lunErr = nullptr; // fopen("error1.dat","w");
// FILE* lunErr= fopen ("file.txt","w");

std::ofstream fout_v("log.txt");

//__________________________________________________________________________
MpdVectorFinder::MpdVectorFinder(const char *name, Bool_t sa, Int_t iVerbose)
   : FairTask("MpdVectorFinder", iVerbose), fNPass(6), fExact(0), //(0),
     fNTracks(0), fTrackExact(nullptr), fGeo(0), fSa(sa),         // by default its + tpc tracking is done
     fTpcKF(nullptr)
{
   fKHits1 = new TClonesArray("MpdKalmanHit", 100);
   for (Int_t i = 0; i < 5; ++i) {
      fKHits[i]  = new TClonesArray("MpdVector", 100);
      f2DHits[i] = new TClonesArray("MpdVector", 100);
   }

   // fKHits[0] = f2DHits[0]; //work
   fTracks    = new TClonesArray("MpdItsKalmanTrack", 100);
   fTrackCand = new TClonesArray("MpdItsKalmanTrack", 100);
   if (fExact) fTrackExact = new TClonesArray("MpdItsKalmanTrack", 100);
   fHistoDir    = 0x0;
   fhLays       = new TH1F("hLaysITS", "ITS layers", 10, 0, 10);
   fLayPointers = 0x0;
}

//__________________________________________________________________________
MpdVectorFinder::~MpdVectorFinder()
{
   // delete fKHits;
   // delete fTracks;
   // delete fTrackCand;
   delete[] fLayPointers;
   delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdVectorFinder::Init()
{

   //----- Get pipe radius
   fPipeR = ((TGeoTube *)gGeoManager->GetVolume("pipe1")->GetShape())->GetRmax();
   cout << "************** Pipe radius: " << fPipeR << " cm " << endl;

   // return ReInit();
   if (ReInit() != kSUCCESS) return kERROR;

   // Read database
   FairRuntimeDb *rtdb   = FairRun::Instance()->GetRuntimeDb();
   MpdStsGeoPar  *geoPar = (MpdStsGeoPar *)rtdb->getContainer("MpdStsGeoPar");
   // cout << geoPar << endl;
   TString    volName = "sts01 ", path = "";
   TObjArray *sensNodes = geoPar->GetGeoSensitiveNodes();
   // cout << sensNodes->GetEntriesFast() << " " << geoPar->GetGeoPassiveNodes()->GetEntriesFast() << endl;

   //----- Process modular geometry

   fGeo       = 1;
   Int_t nLay = 5;
   // Double_t safety = 0.005;

   for (Int_t i = 0; i < nLay; ++i) {
      fNladders[i] = fNsectors[i] = 0;
      volName                     = "sts01ladder";
      volName += (i + 1);
      path = "/cave_1/sts01_0/" + volName;
      path += "_";
      TString path0 = path;

      //----- Loop over all ladders to find the one with the smallest radius
      // fRad[i] = 999.9;
      fRad[i]     = 0.;
      Int_t nladd = 0;

      for (Int_t jlad = 1; jlad < 99; ++jlad) {
         path = path0;
         path += jlad;
         if (!gGeoManager->CheckPath(path)) break;
         gGeoManager->cd(path);
         TGeoVolume *ladd    = gGeoManager->GetVolume(volName);
         TGeoBBox   *box     = (TGeoBBox *)ladd->GetShape();
         Double_t    xyzL[3] = {0}, xyzM[3];
         gGeoManager->LocalToMaster(xyzL, xyzM);
         Double_t rad = TMath::Sqrt(xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
         fRad[i] += rad;
         ++nladd;
         /*AZ - 3.12.2017
      fRad[i] = TMath::Min (fRad[i],rad);
      xyzL[0] = box->GetDX();
      gGeoManager->LocalToMaster(xyzL,xyzM);
      rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
      fRad[i] = TMath::Min (fRad[i],rad);
      xyzL[0] = -box->GetDX();
      gGeoManager->LocalToMaster(xyzL,xyzM);
      rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
      fRad[i] = TMath::Min (fRad[i],rad);
         */
      }
      cout << " *** Layer # " << i << " min_R= " << fRad[i] << endl;
      TGeoVolume *ladd = gGeoManager->GetVolume(volName);
      if (!ladd) {
         nLay = i;
         break;
      }
      fRad[i] /= nladd; // mean radius
      /*
       TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
       //safety = -box->GetDY();
       //safety = box->GetDY();
       safety = 2 * box->GetDY(); // new
       fRad[i] -= safety;
      */

      cout << " +++ Layer # " << i << " safety min_R= " << fRad[i] << endl;
   }
   FillGeoScheme();

   //----- Get cables
   TObjArray *vols  = gGeoManager->GetListOfVolumes();
   Int_t      nvols = vols->GetEntriesFast(), ncables = 0;
   for (Int_t i = 0; i < nvols; ++i) {
      TGeoVolume *vol   = (TGeoVolume *)vols->UncheckedAt(i);
      TString     cable = TString(vol->GetName());
      if (!(cable.Contains("sts") && cable.Contains("cable"))) continue;
      //  cout << cable << endl;
      ++ncables;
      TGeoBBox *box = (TGeoBBox *)vol->GetShape();
      TString   lad = cable;
      Int_t     ip  = lad.Index("cable");
      lad.Replace(ip, lad.Length() - ip + 1, "");
      lad.Replace(lad.Length() - 2, 1, "");
      lad += "_1/";
      path = "/cave_1/sts01_0/" + lad + cable;
      path += "_1";
      gGeoManager->cd(path);
      Double_t xyzL[3] = {0}, xyzM[3];
      gGeoManager->LocalToMaster(xyzL, xyzM);
      //  cout << xyzM[2] - box->GetDZ() << " " << xyzM[2] + box->GetDZ() << " " << box->GetDZ() << endl;
      if (xyzM[2] - box->GetDZ() > 0) {
         Int_t lay = TString(lad(lad.Length() - 4, 1)).Atoi();
         fCables[lay - 1].insert(pair<Double_t, Double_t>(xyzM[2] - box->GetDZ(), xyzM[2] + box->GetDZ()));
      }
   }
   cout << "-I- MpdTrackFinderIts5spd: Intialization finished successfully" << endl;
   fout_v << "-Init - MpdTrackFinderIts5spd: Intialization finished successfully" << endl;

   return kSUCCESS;
}

//__________________________________________________________________________

void MpdVectorFinder::FillGeoScheme()
{
   /// Fill Kalman filter geometry manager info

   MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

   TGeoVolume *vSts = gGeoManager->GetVolume("sts01");

   for (Int_t layer = 1; layer < 999; ++layer) {
      // Loop over layers
      TString sladder = "sts01ladder";
      sladder += layer;
      TGeoVolume *vLay = gGeoManager->GetVolume(sladder);
      if (vLay == 0x0 && layer == 1) continue; // no first layer
      if (vLay == 0x0) break;
      sladder += "_";

      //----- Loop over ladders
      for (Int_t ladder = 1; ladder < 999; ++ladder) {
         TString sladder1 = sladder;
         sladder1 += ladder;
         TGeoNode *node = vSts->FindNode(sladder1);
         if (node == 0x0) break;
         ++fNladders[layer - 1];
         TGeoVolume *vLad = node->GetVolume();
         // cout << vLad->GetNodes()->GetEntriesFast() << " " << vLad->GetNdaughters() << endl;
         Int_t      nDaught = vLad->GetNdaughters(), detID = -1, detIDsts = -1;
         TObjArray *daught = vLad->GetNodes();

         //----- Loop over ladder daughters
         Int_t iZ = 0;
         for (Int_t det = 0; det < nDaught; ++det) {
            TString sdet1 = ((TGeoNode *)(daught->UncheckedAt(det)))->GetName();
            /// if (!sdet1.Contains("sensor") && !sdet1.Contains("sector")) continue;
            if (ladder == 1) {
               ++fNsectors[layer - 1];
            }
            Int_t det1    = TString(sdet1(sdet1.Index("_") + 1, 2)).Atoi();
            Int_t secType = -1;
            if (sdet1.Contains("sector")) secType = TString(sdet1(sdet1.Index("_") - 2, 1)).Atoi();

            ++iZ;

            Int_t side = 0.;

            detIDsts = fHitSts.SetDetId(layer, ladder, det1);
            detID    = fHitSts.SetDetId(layer, ladder, iZ);
            fId2Id[layer - 1].insert(pair<Int_t, Int_t>(detIDsts, detID)); // detIDsts == detID for pixel detectors (AZ)
            /// cout << "Detector ids " << detIDsts << " " << detID << endl;
            detID += 1000000 * (layer - 1);

            TString detName = sladder1 + "/" + sdet1 + "#";
            detName += side;
            geo->SetDetId(detName, detID);

            TString path = "/cave_1/sts01_0/" + detName(0, detName.Length() - 2);
            // cout << detName << " " << path << endl;
            gGeoManager->cd(path);
            geo->SetPath(detID, path);
            node = gGeoManager->GetCurrentNode();
            // cout << node->GetName() << " " << detID << endl;

            TGeoVolume *vol = node->GetVolume();
            TGeoBBox   *box = (TGeoBBox *)vol->GetShape();
            TVector2    size(box->GetDX(), box->GetDZ());
            geo->SetSize(detID, size);

            Double_t xyzL[3] = {0}, xyzM[3], vecM[3];
            xyzL[1]          = 1;
            gGeoManager->LocalToMasterVect(xyzL, vecM);

            xyzL[1] = 0;
            gGeoManager->LocalToMaster(xyzL, xyzM);
            Double_t sign = TMath::Sign(1., xyzM[0] * vecM[0] + xyzM[1] * vecM[1]);

            TVector3 pos(xyzM[0], xyzM[1], xyzM[2]);
            geo->SetGlobalPos(detID, pos);

            TVector3 norm(sign * vecM[0], sign * vecM[1], sign * vecM[2]);
            geo->SetNormal(detID, norm);
         }
      }
   }
   fout_v << "- FillGeoScheme - done" << endl;
}

//__________________________________________________________________________
InitStatus MpdVectorFinder::ReInit()
{
   FairRootManager *manager = FairRootManager::Instance();

   fItsPoints = (TClonesArray *)manager->GetObject("StsPoint");
   fItsHits   = (TClonesArray *)manager->GetObject("StsHit");
   if (fItsPoints == 0x0 || fItsHits == 0x0) return kERROR;
   fTpcTracks = (TClonesArray *)manager->GetObject("TpcKalmanTrack");
   fMCTracks  = (TClonesArray *)manager->GetObject("MCTrack");

   manager->Register("ItsTrack", "Its", fTracks, kTRUE);
   manager->Register("ItsTrackCand", "ItsCand", fTrackCand, kTRUE);
   if (fTrackExact) manager->Register("ItsTrackExact", "ItsCand", fTrackExact, kTRUE);

   manager->Register("ItsKHits", "ItsKalHits", fKHits1, kFALSE);

   manager->Register("CellTrack0", "Its0", fKHits[0], kFALSE);
   manager->Register("CellTrack1", "Its1", fKHits[1], kFALSE);
   manager->Register("CellTrack2", "Its2", fKHits[2], kFALSE);
   manager->Register("CellTrack3", "Its3", fKHits[3], kFALSE);
   manager->Register("CellTrack4", "Its4", fKHits[4], kFALSE);

   manager->Register("CellTrack_0", "Its0_2D", f2DHits[0], kFALSE);
   manager->Register("CellTrack_1", "Its1_2D", f2DHits[1], kFALSE);
   manager->Register("CellTrack_2", "Its2_2D", f2DHits[2], kFALSE);
   manager->Register("CellTrack_3", "Its3_2D", f2DHits[3], kFALSE);
   manager->Register("CellTrack_4", "Its4_2D", f2DHits[4], kFALSE);

   // fNPass = 2; // 2 prochoda
   // fNPass = 3;
   fNPass = 6;
   fout_v << "- ReInit - done" << endl;
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdVectorFinder::Reset()
{

   cout << " MpdVectorFinder::Reset  " << endl;
   // fKHits->Clear();
   /// changed list of function calls to a loop
   for (Int_t i = 0; i < 5; ++i) { /// was 4
      fKHits[i]->Delete();
      f2DHits[i]->Delete();
   }
   fKHits1->Delete();
   fTracks->Delete();
   fTrackCand->Delete();
   if (fTrackExact) fTrackExact->Delete();
   delete[] fLayPointers;
   fLayPointers = nullptr;
   fCellMap.clear();
   fout_v << "- Reset - done" << endl;
}

//__________________________________________________________________________
void MpdVectorFinder::SetParContainers()
{
   // Get run and runtime database
   FairRunAna *run = FairRunAna::Instance();
   if (!run) Fatal("SetParContainers", "No analysis run");

   FairRuntimeDb *db = run->GetRuntimeDb();
   if (!db) Fatal("SetParContainers", "No runtime database");

   // Get STS geometry parameter container
   db->getContainer("MpdStsGeoPar");
}

//__________________________________________________________________________
void MpdVectorFinder::Finish()
{
   // Write();
   if (lunErr) fclose(lunErr);
}

//__________________________________________________________________________
void MpdVectorFinder::Exec(Option_t *option)
{
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   static int eventCounter = 0;
   cout << " - - - - \n ItsRec event " << ++eventCounter << endl;

   Reset();

   tStart = clock();
   // Create Kalman hits
   if (fItsHits->GetEntriesFast() == 0) return;
   Build2DHits(); /// MakeKalmanHits stuff was added to Build2DHits

   /// if MpdKalmanHit->GetFlag == -1 then hit was used earlier (see ExcludeHits())

   /// fNPass is set in Reinit function
   /// i is for pass number:
   /// 0 - primary tracks, 5-layer - tight cuts
   //  1 - primary tracks, 5-layer - loose cuts
   /// 2 - primary tracks, 4-layer
   /// 3 - secondary tracks, 5-layer - tight cuts
   /// 4 - secondary tracks, 5-layer - loose cuts
   /// 5 - secondary tracks, 4-layer
   if (fTrackExact) fTrackExact->Delete();

   for (Int_t i = 0; i < fNPass; ++i) {

      /// resetting after each pass updated 19.2.2020
      for (Int_t j = 0; j < 5; ++j) {
         fKHits[j]->Delete();
      }
      fTrackCand->Delete(); // AZ-29.02.20
      fCellMap.clear();
      fVectorCands.clear();
      ///

      MakeTrackCandidates(i); // 21.02
      if (i < 3)
         ExtendCellTracks(i); // 21.02
      else
         ExtendSecondaryTracks(i);

      GetTrackSeeds(i); // attention!!
      DoTracking(i);    // attention!!!
      // cout << " AZ before: " << fTrackCand->GetEntriesFast() << endl;
      RemoveDoubles();
      // cout << " AZ after: " << fTrackCand->GetEntriesFast() << endl;
      StoreTracks();
      fout_v << "  Total number of found tracks: " << fTrackCand->GetEntriesFast() << endl;

      ExcludeHits(); // exclude used hits  test

      Int_t usedHits = 0, unusedHits = 0;
      for (int j = 0; j < fKHits1->GetEntriesFast(); j++) {
         if (((MpdKalmanHit *)fKHits1->UncheckedAt(j))->GetFlag() < 0) {
            usedHits++;
         } else {
            unusedHits++;
         }
      }
      cout << "\n used/unused hits before making track candidates: " << usedHits << " " << unusedHits << endl;
      fout_v << "\n used/unused hits before making track candidates: " << usedHits << " " << unusedHits << endl;
   }

   GetShortTracks();
   AddHits(); // add hit objects to tracks

   tFinish = clock();
   tAll    = tFinish - tStart;

   cout << "  Total number of found tracks: " << fTrackCand->GetEntriesFast() << endl;

   cout << " Exec time: " << ((Float_t)tAll) / CLOCKS_PER_SEC << endl;

   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   fout_v.close();
}

//_________________________________________________________________________
void MpdVectorFinder::Build2DHits()
{
   /// Build track candidates from ITS 2-D hits.

   fItsHits->Sort();
   Int_t    nHits = fItsHits->GetEntriesFast(); //, layMax = 0, nKH = 0;
   Double_t errZ = 0.001, errX = 0.001;         /// 10um for Z and X// 250um in Z, 23um in R-Phi (local X)
   Double_t xloc;                               //, z;
   Double_t r;
   // Double_t dX = 0, dZ = 0;
   MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

   Int_t nKalmanHit = 0, nLayerHit = 0;

   for (Int_t ih = 0; ih < nHits; ++ih) {
      MpdItsHit5spd *h = (MpdItsHit5spd *)fItsHits->UncheckedAt(ih);
      if (h->GetFlag() < 0) continue; /// TODO check if necessary. Build2DHits is done only once
      Int_t lay = h->Layer() - 1;
      // AZ - for debug
      // Int_t id = ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();
      // if (id != 2346 && id != 2305 && id != 2876) continue;

      // AZ Get local Z
      // TString path = geo->Path(h->GetDetectorID()+1000000*lay);
      TString path = geo->Path(fId2Id[lay][h->GetDetectorID()] + 1000000 * lay);
      gGeoManager->cd(path);
      Double_t posLoc[3] = {0}, posMas[3] = {h->GetX(), h->GetY(), h->GetZ()};
      /// transform global coords to local
      gGeoManager->MasterToLocal(posMas, posLoc);
      // AZ

      r    = TMath::Sqrt(h->GetX() * h->GetX() + h->GetY() * h->GetY());
      xloc = h->GetLocalX();
      // z = h->GetZ();
      // AZ Double_t meas[2] = {xloc + h->GetDx() * errX, z + dZ * errZ};
      Double_t meas[2] = {xloc + h->GetDx() * errX, posLoc[2] + h->GetDz() * errZ};
      // AZ Double_t err[2] = {errX, errZ};
      // AZ-161020 Double_t err[2] = {errX*1.2, errZ*1.2}; // effective measure
      Double_t err[2] = {errX * 1.4, errZ * 1.4}; // effective measure
      if (lay > 2) {
         err[0] *= 2;
         err[1] *= 2;
      }                                                                              // AZ-161020 - thick layers
      Double_t       cossin[2] = {TMath::Cos(fStereoA[0]), TMath::Sin(fStereoA[0])}; /// wtf is this
      MpdVector     *hit       = 0x0;                                                // ms
      MpdItsHit5spd *hsts      = h;

      Double_t coord[3] = {h->GetX(), h->GetY(), h->GetZ()};
      TVector3 point(coord);

      /// fout_v << "ih = " << ih << " layer = " << h->Layer()-1 << " detID = " << hsts->GetDetectorID() << " " <<
      /// fId2Id[lay][hsts->GetDetectorID()] << endl;

      MpdKalmanHit *khit = 0x0;
      /// fId2Id[lay][hsts->GetDetectorID()] - encoded hit identifier (# of sensor, layer, ladder)
      khit = new ((*fKHits1)[nKalmanHit++])
         MpdKalmanHit(lay * 1000000 + fId2Id[lay][hsts->GetDetectorID()], // detID
                      2,                                                  // nDim
                      MpdKalmanHit::kFixedP, // HitType // different modules (planes) in local coordinates
                      meas,                  // 2D measurement for kalman filter (x and z)
                      err,                   // err
                      cossin,                //
                      0.,                    // signal
                      r,                     // dist
                      ih);                   // index
      // AZ khit->SetUniqueID(0);
      khit->SetUniqueID(nKalmanHit); // AZ

      nLayerHit = f2DHits[lay]->GetEntriesFast();
      hit       = new ((*f2DHits[lay])[nLayerHit])
         MpdVector(lay * 1000000 + fId2Id[lay][hsts->GetDetectorID()], // detID
                   2,                                                  // nDim
                   MpdVector::kFixedP, // HitType // different modules (planes) in local coordinates
                   point,              // TVector3 measurement
                   err,                // err
                   cossin,             //
                   0.,                 // signal
                   r,                  // dist
                   ih,                 // index
                   0,                  // index1
                   -1,                 // trackNo - deprecated?
                   nullptr,            // TrackPointer
                   khit);              // KalmanHit pointer
      hit->SetUniqueID(nLayerHit);
      /// cout << hit->GetKalmanHit() << " ";
   } // for (Int_t ih = 0; ih < nHits; ++ih)

   fout_v << "khits " << nKalmanHit << endl;
   cout << "- Build2DHits - done" << endl;
}

//_________________________________________________________________________
void MpdVectorFinder::MakeTrackCandidates(Int_t iPass)
{
   // Read 2D Hits and create trackCandidates from 2D Hits at first(last) layer
   // Int_t n = 5 - iPass; // last layer
   // n = TMath::Min (n, 4);
   const Int_t layNo[9] = {4, 4, 3, 4, 4, 3};

   Int_t n   = layNo[iPass];
   Int_t nKH = 0, nHits = f2DHits[n]->GetEntriesFast();

   if (iPass < 3) {
      fout_v << "layer " << n << " hits " << nHits << endl;
      for (Int_t i = 0; i < nHits; ++i) {

         MpdVector *cellTr = (MpdVector *)f2DHits[n]->UncheckedAt(i);
         /// MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits-> UncheckedAt(cellTr->GetIndex(0));

         /// fout_v << h->GetTrackID() << " " << cellTr->GetMeas()[2] << endl;
         cellTr = new ((*fKHits[n])[nKH++]) MpdVector(*cellTr);
         /// fout_v << cellTr->GetKalmanHit() << " ";
         cellTr->SetCode(i);

         cellTr->SetCosSin(
            0, TMath::ATan2(cellTr->GetMeas()[1], cellTr->GetMeas()[0])); /// transverse angle for first layer hits
         Float_t r = TMath::Sqrt(cellTr->GetMeas()[0] * cellTr->GetMeas()[0] +
                                 cellTr->GetMeas()[1] * cellTr->GetMeas()[1]); /// distance from interaction point
         cellTr->SetCosSin(1, TMath::ATan2(r, cellTr->GetMeas()[2]));          /// longitudinal angle

         /// added 17.1.2019
         cellTr->SetUniqueID(nKH - 1);
         /// fout_v << nKH - 1 << endl;
         /// fout_v << cellTr->GetKalmanHit() << " ";
         // AZ
         // MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits->UncheckedAt(cellTr->GetIndex(0));
         // Int_t trID = ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();
         // cout << " AZ-TrackID: " << trID << endl;
      }     // for (Int_t i = 0; ih < nHits; ++i)
   } else { /// added 29.1.2020
      // 4th and 5th pass are done for secondary particles
      // n = 5 - iPass + 2;
      nHits = f2DHits[n]->GetEntriesFast();
      /// First hit of track candidate is a hit from 1st layer
      for (Int_t i = 0; i < nHits; ++i) {
         MpdVector *cellTr = (MpdVector *)f2DHits[n]->UncheckedAt(i);
         if ((cellTr->GetKalmanHit())->GetFlag() < 0) {
            continue;
         }
         cellTr = new ((*fKHits[n])[nKH++]) MpdVector(*cellTr);
         cellTr->SetCode(i);
         cellTr->SetCosSin(
            0, TMath::ATan2(cellTr->GetMeas()[1], cellTr->GetMeas()[0])); /// transverse angle for first layer hits
         Float_t r = TMath::Sqrt(cellTr->GetMeas()[0] * cellTr->GetMeas()[0] +
                                 cellTr->GetMeas()[1] * cellTr->GetMeas()[1]); /// distance from interaction point
         cellTr->SetCosSin(1, TMath::ATan2(r, cellTr->GetMeas()[2])); // cellTr->GetMeas()[2] /// z coordinate instead
                                                                      // of longitudinal angle? /// TODO 21.2.20
         cellTr->SetUniqueID(nKH - 1);
         // AZ
         cellTr->SetDeltaZ(-1 / TMath::Tan(cellTr->GetCosSin(1)));
         // MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits->UncheckedAt(cellTr->GetIndex(0));
         // Int_t trID = ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();
         // cout << " AZ1-TrackID: " << trID << " " <<  cellTr->GetMeas()[2] << endl;
      }
      cout << "Track Candidates on pass " << iPass + 1 << ": " << nKH << endl;
   }

   fout_v << "- MakeTrackCandidates - done" << endl;
}

//_________________________________________________________________________
/// added 30.1.2020
void MpdVectorFinder::ExtendSecondaryTracks(Int_t iPass)
{
   // const Int_t nSigm = 3;
   int count = 0;
   // Double_t  errZ = 0.001, errX = 0.001;/// 10um for Z and X // 250um in Z, 23um in R-Phi (local X)
   multimap<Float_t, MpdVector *> zMultimapHits[5], tMultimapHits[5];

   fout_v << "(transverse angle, longitudinal angle, z | MCtrackID)" << endl;

   // Int_t n = 4 - iPass + 2;
   const Int_t layNo[9] = {0, 0, 0, 3, 3, 2};

   // for (Int_t i234 = n; i234 >= 0; i234--) { /// loop over layers 3-0
   for (Int_t i234 = layNo[iPass]; i234 >= 0; i234--) {                             /// loop over layers 3-0
      Int_t nKH = 0;                                                                //, iprint = 0;
      cout << " AZsec: " << i234 << " " << f2DHits[i234]->GetEntriesFast() << endl; // AZ
      Double_t rmean = 0;                                                           // AZ

      for (Int_t i = 0; i < f2DHits[i234]->GetEntriesFast(); i++) {
         MpdVector *hit = (MpdVector *)f2DHits[i234]->UncheckedAt(i);
         if (hit->GetKalmanHit()->GetFlag() < 0) continue; // AZ

         TVector3 meas = hit->GetMeas();
         Double_t phi  = TMath::ATan2(meas[1], meas[0]);
         zMultimapHits[i234].insert(pair<Float_t, MpdVector *>(meas[2], hit)); /// z coordinate
         tMultimapHits[i234].insert(pair<Float_t, MpdVector *>(phi, hit));     /// transverse, phi

         /// if transverse angle is close to pi, add it to other side of map with angle - 2pi
         if (TMath::Abs(TMath::Pi() - TMath::Abs(phi)) < TMath::Pi() / 9.0) {
            tMultimapHits[i234].insert(pair<Float_t, MpdVector *>(phi - TMath::Sign(1, phi) * 2 * TMath::Pi(),
                                                                  hit)); /// transverse, phi, duplicate hit
         }
         rmean += hit->GetDist(); // AZ
      }
      rmean /= zMultimapHits[i234].size(); // AZ

      int nTracks = fKHits[i234 + 1]->GetEntriesFast(); /// was -1
      cout << "secondary nTracks for layer " << i234 + 1 << ": " << nTracks << " " << rmean << endl;

      for (Int_t i = 0; i < nTracks; i++) {
         MpdVector *cand = (MpdVector *)fKHits[i234 + 1]->UncheckedAt(i); /// was -1
         TVector3   meas = cand->GetMeas();

         MpdItsHit5spd *h     = (MpdItsHit5spd *)fItsHits->UncheckedAt(cand->GetIndex(0));
         Int_t          trID0 = ((MpdStsPoint *)fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();

         Double_t pt = 0, epst_lower, epst_upper;
         // AZ Double_t epsz_lower = 5e-2, epsz_upper = 5e-2;
         // Double_t epsz_lower = 3, epsz_upper = 3;
         Double_t                         epsz_lower = 5, epsz_upper = 5, dphi = 0;
         MpdVector                       *prev = cand->GetPrevTrackPointer(), *preprev = nullptr;
         vector<pair<Double_t, Double_t>> epst(2,
                                               pair<Double_t, Double_t>(999, 999)); // for fine structure (2 branches)

         if (prev != nullptr) {
            preprev            = prev->GetPrevTrackPointer();
            Double_t circle[3] = {0};

            if (preprev == nullptr) {
               // if (1) {
               //  Use primary vertex
               /// fout_v << cand->GetPrevTrackPointer() << endl;
               MpdVector zeroV;
               pt = EvalPt(prev, cand, &zeroV, circle);
               if (TMath::Abs(pt) < 0.02) continue;
               if (!PropagateHelix(rmean, meas, circle, dphi)) continue;
               Double_t ptabs = TMath::Abs(pt);
               if (i234 == 2) {
                  // First region
                  Double_t ptmin1 = TMath::Max(ptabs, 0.1);
                  ptmin1          = TMath::Min(ptmin1, 0.2);
                  Double_t eps1   = 0.1 / TMath::Power(ptmin1 - 0.01, 0.3);
                  Double_t ptmin2 = TMath::Max(ptabs, 0.1);
                  ptmin2          = TMath::Min(ptmin2, 0.3);
                  Double_t eps2   = 0.04 / TMath::Power(ptmin2 - 0.01, 0.8);
                  epst_lower      = (pt > 0) ? eps1 : eps2;
                  epst_upper      = (pt > 0) ? eps2 : eps1;
                  epst[0].first   = -epst_lower;
                  epst[0].second  = epst_upper;
                  // Second region
                  ptmin1          = TMath::Max(ptabs, 0.02);
                  Double_t phicor = 0.0569913 / TMath::Power(ptmin1 - 0.01 + 0.0627805, 0.894263);
                  if (pt < 0) phicor = -phicor;
                  eps1           = 0.01;
                  eps2           = 0.06;
                  epst_lower     = (pt > 0) ? eps2 : eps1;
                  epst_upper     = (pt > 0) ? eps1 : eps2;
                  epst[1].first  = -epst_lower + phicor;
                  epst[1].second = epst_upper + phicor;
               } else if (i234 == 1) {
                  // First region
                  Double_t eps1  = 0.004 / TMath::Power(ptabs - 0.01, 1.2);
                  eps1           = TMath::Max(eps1, 0.03);
                  Double_t eps2  = 0.1;
                  epst_lower     = (pt > 0) ? eps2 : eps1;
                  epst_upper     = (pt > 0) ? eps1 : eps2;
                  epst[0].first  = -epst_lower;
                  epst[0].second = epst_upper;
                  // Second region
                  Double_t ptmin1 = TMath::Max(ptabs, 0.02);
                  Double_t phicor = 0.0276577 / TMath::Power(ptmin1 - 0.01, 0.827315);
                  if (pt < 0) phicor = -phicor;
                  Double_t ptmax1 = TMath::Min(ptmin1, 0.3);
                  eps1            = 0.04 / TMath::Power(ptmax1 - 0.01, 0.5);
                  eps1            = TMath::Min(eps1, 0.1);
                  eps2            = 0.017 / TMath::Power(ptmin1 - 0.01, 0.4);
                  eps2            = TMath::Min(eps2, 0.08);
                  epst_lower      = (pt > 0) ? eps2 : eps1;
                  epst_upper      = (pt > 0) ? eps1 : eps2;
                  epst[1].first   = -epst_lower + phicor;
                  epst[1].second  = epst_upper + phicor;
               }
               // epst_lower = epst_upper = 5 * epst_upper;
               // if (i234 == 1) epst_lower = epst_upper = 0.02 / TMath::Power (ptt-0.01, 2./3);
               // epsz_lower = epsz_upper = 1.;
               Double_t ptt1 = TMath::Min(0.1, ptabs);
               epsz_lower = epsz_upper = 0.7 / TMath::Power(ptt1 - 0.015, 0.5);
               if (i234 == 1) epsz_lower = epsz_upper = 1.5;
               // epsz_lower = epsz_upper = 5 * epsz_upper;
            } else {
               // Do not use primary vertex
               // Double_t circle[3] = {0};
               pt = EvalPt(preprev, prev, cand, circle);
               if (TMath::Abs(pt) < 0.02) continue;
               if (!PropagateHelix(rmean, meas, circle, dphi)) continue;
               Double_t ptt = TMath::Abs(pt);
               ptt          = TMath::Min(0.1, ptt);
               // epst_lower = epst_upper = 0.1;
               epst_lower = epst_upper = 0.04 / TMath::Power(ptt - 0.01, 0.2);
               // if (i234 == 0) epst_lower = epst_upper = 0.02 / TMath::Power (ptt-0.01, 2./3);
               epsz_lower = epsz_upper = 1.;
               if (iPass == 3) {
                  epst_lower /= 2;
                  epst_upper /= 2;
                  epsz_lower /= 2;
                  epsz_upper /= 2;
               }
               if (i234 == 0) {
                  /*
                  //epst_lower = epst_upper = 0.2;
                  Double_t ptt1 = TMath::Abs(pt);
                  ptt1 = TMath::Min (0.15, ptt1);
                  ptt1 = TMath::Max (0.05, ptt1);
                  //epst_lower = 0.025 / TMath::Power (ptt-0.01, 2./3);
                  //epst_upper = 0.02 / TMath::Power (ptt-0.01, 2./3);
                  epst_lower = 0.02 / TMath::Power (ptt-0.01, 2./3);
                  epst_upper = 0.025 / TMath::Power (ptt-0.01, 2./3);
                  */
                  epst_lower = 0.1;
                  epst_upper = 0.15;
                  epsz_lower = epsz_upper = 1.5;
               }
               epst[0].first  = -epst_lower;
               epst[0].second = epst_upper;
            }

         } else {
            // First window
            epst_lower = epst_upper = 0.3; // 5e-1;
            if (i234 == 2) epst_lower = epst_upper = 0.4;
            // epst_lower = epst_upper = epst_upper*5;
            // epsz_lower = epsz_upper = 4.0; //1.;
            // if (i234 == 2) epsz_lower = epsz_upper = 6.0;
            epsz_lower = epsz_upper = 10.0; // 1.;
            if (i234 == 2) epsz_lower = epsz_upper = 20.0;
            if (iPass == 3) {
               epst_lower /= 2;
               epst_upper /= 2;
               epsz_lower /= 2;
               epsz_upper /= 2;
            }
            // epsz_lower = epsz_upper = epsz_upper*5;
            epst[0].first  = -epst_lower;
            epst[0].second = epst_upper;
         }

         Double_t lastTrackZ = meas[2];
         Double_t drad       = TMath::Abs(cand->GetDist() - rmean); // AZ
         Double_t deltaZ     = cand->GetDeltaZ() * drad;            // AZ cand->GetDeltaZ() is tan(theta)

         // multimap<Float_t, MpdVector*>::iterator itlow = zMultimapHits[i234].lower_bound(lastTrackZ - epsz);
         // multimap<Float_t, MpdVector*>::iterator ittop = zMultimapHits[i234].upper_bound(lastTrackZ + epsz);

         // AZ multimap<Float_t, MpdVector*>::iterator itlow = zMultimapHits[i234].lower_bound(lastTrackZ - epsz_lower);
         // AZ multimap<Float_t, MpdVector*>::iterator ittop = zMultimapHits[i234].upper_bound(lastTrackZ + epsz_upper);
         multimap<Float_t, MpdVector *>::iterator itlow =
            zMultimapHits[i234].lower_bound(lastTrackZ - epsz_lower + deltaZ);
         multimap<Float_t, MpdVector *>::iterator ittop =
            zMultimapHits[i234].upper_bound(lastTrackZ + epsz_upper + deltaZ);

         Double_t lastTrackT = 0;
         // AZ if (preprev) lastTrackT = dphi; // circle from 3 hits
         if (prev)
            lastTrackT = dphi; // AZ-141020 - circle from 2 hits and (0,0) or from 3 hits
         else {
            lastTrackT = TMath::ATan2(meas[1], meas[0]);
            lastTrackT -= dphi;
         }

         // multimap<Float_t, MpdVector*>::iterator ittlow = tMultimapHits[i234].lower_bound(lastTrackT - epst_lower);
         // multimap<Float_t, MpdVector*>::iterator itttop = tMultimapHits[i234].upper_bound(lastTrackT + epst_upper);
         multimap<Float_t, MpdVector *>::iterator ittlow = tMultimapHits[i234].lower_bound(lastTrackT + epst[0].first);
         multimap<Float_t, MpdVector *>::iterator itttop = tMultimapHits[i234].upper_bound(lastTrackT + epst[0].second);

         set<MpdVector *> z, t, intersect;
         for (multimap<Float_t, MpdVector *>::iterator itr = itlow; itr != ittop; ++itr) {
            z.insert((*itr).second);
         }
         for (multimap<Float_t, MpdVector *>::iterator itr = ittlow; itr != itttop; ++itr) {
            t.insert((*itr).second);
         }

         if (epst[1].first < 900) {
            // Fine structure - second region
            // lastTrackT = epst[1][2];
            ittlow = tMultimapHits[i234].lower_bound(lastTrackT + epst[1].first);
            itttop = tMultimapHits[i234].upper_bound(lastTrackT + epst[1].second);

            for (multimap<Float_t, MpdVector *>::iterator itr = ittlow; itr != itttop; ++itr) t.insert((*itr).second);
         }

         set_intersection(z.begin(), z.end(), t.begin(), t.end(), std::inserter(intersect, intersect.begin()));

         for (set<MpdVector *>::iterator itr = intersect.begin(); itr != intersect.end(); ++itr) {
            MpdVector *track = (*itr);
            if (track->GetKalmanHit()->GetFlag() < 0) continue; // AZ
            MpdItsHit5spd *hh   = (MpdItsHit5spd *)fItsHits->UncheckedAt(track->GetIndex(0));
            Int_t          trID = ((MpdStsPoint *)fItsPoints->UncheckedAt(hh->GetTrackID()))->GetTrackID();
            if (fExact && trID != trID0) continue; // exact ID match

            MpdVector *trNew = nullptr;
            trNew            = new ((*fKHits[i234])[nKH++]) MpdVector(*track); /// was h->Layer()
            trNew->SetUniqueID(nKH - 1);
            trNew->SetNofDim(2);
            trNew->SetCosSin(0, 0);
            trNew->SetCosSin(1, 0);
            trNew->SetPrevTrackPointer(cand);
            trNew->SetCode(cand->GetCode());
            trNew->SetCode(i);

            /// 21.2.20
            trNew->SetDeltaZ(trNew->GetMeas()[2] - trNew->GetPrevTrackPointer()->GetMeas()[2]);
            /// cout << "trNew " << trNew->GetMeas()[2] << " " << trNew->GetPrevTrackPointer()->GetMeas()[2] << endl;
            ///  25.2.20
            // AZ
            // trNew->SetDeltaZ(trNew->GetDeltaZ()/TMath::Abs(trNew->GetDist()-trNew->GetPrevTrackPointer()->GetDist()));
            trNew->SetDeltaZ(trNew->GetDeltaZ() / (trNew->GetMeas() - trNew->GetPrevTrackPointer()->GetMeas()).Pt());
            // cout << " AZ1-TrackID: " << trID << " " << i234 << endl;
            ///

            Int_t trackID   = 0;
            Int_t MCtrackID = 0;
            if (i234 == 0) {
               Int_t f = 0;

               // AZ Int_t l = iPass - 3;

               // AZ for (int j = 0; j < 4 - l; j++) {
               for (int j = 0; j <= layNo[iPass]; j++) {
                  trackID           = ((MpdItsHit5spd *)fItsHits->UncheckedAt(trNew->GetIndex(0)))->GetTrackID();
                  Int_t tempTrackID = MCtrackID;
                  MCtrackID         = ((MpdStsPoint *)fItsPoints->UncheckedAt(trackID))->GetTrackID();

                  if ((tempTrackID != MCtrackID) && (j != 0)) {
                     f = 1;
                  }

                  fout_v << "("
                         << TMath::ATan2(trNew->GetMeas()[1], trNew->GetMeas()[0])
                         /*<< " " << TMath::ATan2(v, u)*/
                         << " "
                         << TMath::ATan2(TMath::Sqrt(trNew->GetMeas()[0] * trNew->GetMeas()[0] +
                                                     trNew->GetMeas()[1] * trNew->GetMeas()[1]),
                                         trNew->GetMeas()[2])
                         << " " << trNew->GetMeas()[2]
                         << "|" //<< trNew->GetMeas()[1] << "|" << trNew->GetMeas()[2] << " " */
                                //<< " " << trNew->GetKalmanHit()->GetMeas(0)
                         << /*trackID <<*/ " " << MCtrackID << ")";
                  // MpdVector *trOld = trNew;

                  trNew = trNew->GetPrevTrackPointer();
               }

               trackID           = ((MpdItsHit5spd *)fItsHits->UncheckedAt(trNew->GetIndex(0)))->GetTrackID();
               Int_t tempTrackID = MCtrackID;
               MCtrackID         = ((MpdStsPoint *)fItsPoints->UncheckedAt(trackID))->GetTrackID();

               fout_v << "(" << TMath::ATan2(trNew->GetMeas()[1], trNew->GetMeas()[0])
                      << " " /*<< TMath::ATan2(v, u) << " "*/
                      << TMath::ATan2(TMath::Sqrt(trNew->GetMeas()[0] * trNew->GetMeas()[0] +
                                                  trNew->GetMeas()[1] * trNew->GetMeas()[1]),
                                      trNew->GetMeas()[2])
                      << " " << trNew->GetMeas()[2]
                      << "|" //<< trNew->GetMeas()[1] << "|" << trNew->GetMeas()[2] << /*" " << trackID*/
                             //<< " " << trNew->GetKalmanHit()->GetMeas(0)
                      << " " << MCtrackID << ")" << endl;
               if (tempTrackID != MCtrackID) {
                  f = 1;
               }

               if (f == 1) {
                  ///  fout_v << "!!!" << endl;
               }
               count += f;
            }
         } // for (set<MpdVector*>::iterator itr = intersect.begin();
      }    // for (Int_t i = 0; i < nTracks;
   }       // for (Int_t i234 = n; i234 >= 0;

   fout_v << "Count of false reconstructed tracks " << count << endl;
   fout_v << "- ExtendSecondaryTracks - done" << endl;
}

//_________________________________________________________________________

void MpdVectorFinder::ExtendCellTracks(Int_t iPass)
{
   /// Extend cell tracks to layers 4-1
   /// pass 2 - start from the last but one layer

   // const Int_t nSigm = 3;
   int count = 0;
   // Double_t  errZ = 0.001, errX = 0.001;/// 10um for Z and X // 250um in Z, 23um in R-Phi (local X)
   multimap<Float_t, MpdVector *> lMultimapHits[5], tMultimapHits[5];

   fout_v << "(transverse angle, longitudinal angle, z | MCtrackID)" << endl;

   Int_t n = 6 - 2 - iPass; // layer numbers 0,1,2,3,4
   n       = TMath::Min(3, n);

   for (Int_t i234 = n; i234 >= 0; i234--) {                                         /// loop over layers 3-0
      Int_t nKH = 0;                                                                 //, iprint = 0;
      cout << " AZprim: " << i234 << " " << f2DHits[i234]->GetEntriesFast() << endl; // AZ
      /// nKH = fKHits[i234]->GetEntriesFast();
      // transverse and longitudinal hit maps

      for (Int_t i = 0; i < f2DHits[i234]->GetEntriesFast(); i++) {
         MpdVector *hit = (MpdVector *)f2DHits[i234]->UncheckedAt(i);

         TVector3 meas = hit->GetMeas();
         Double_t phi  = TMath::ATan2(meas[1], meas[0]);
         /// cout << meas[0] << " " << meas[1] << " " << meas[2] << endl;
         lMultimapHits[i234].insert(pair<Float_t, MpdVector *>(
            TMath::ATan2(TMath::Sqrt(meas[0] * meas[0] + meas[1] * meas[1]), meas[2]), hit)); /// longitudinal, theta
         tMultimapHits[i234].insert(pair<Float_t, MpdVector *>(phi, hit));                    /// transverse, phi

         /// possible solution: if transverse angle is close to pi, add it to other side of map with angle - 2pi
         if (TMath::Abs(TMath::Pi() - TMath::Abs(phi)) < TMath::Pi() / 9.0) {
            tMultimapHits[i234].insert(pair<Float_t, MpdVector *>(phi - TMath::Sign(1, phi) * 2 * TMath::Pi(),
                                                                  hit)); /// transverse, phi, duplicate hit
         }
      }

      int nTracks = fKHits[i234 + 1]->GetEntriesFast(); /// was -1
      cout << "primary nTracks for layer " << i234 + 1 << ": " << nTracks << endl;

      for (Int_t i = 0; i < nTracks; i++) {
         MpdVector     *cand  = (MpdVector *)fKHits[i234 + 1]->UncheckedAt(i); /// was -1
         MpdItsHit5spd *h     = (MpdItsHit5spd *)fItsHits->UncheckedAt(cand->GetIndex(0));
         Int_t          trID0 = ((MpdStsPoint *)fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();

         if (i234 == 3) {
            /// fout_v << h->GetTrackID() << endl;
            /// fout_v << cand->GetUniqueID() << endl;
         }

         Double_t epsl, epst, pt = 0, dphi = 0;

         if (cand->GetPrevTrackPointer() != 0) {
            /// fout_v << cand->GetPrevTrackPointer() << endl;
            pt = EvalPt(cand->GetPrevTrackPointer(), cand); /// 0-4-5 instead of 0-5-4
            if (TMath::Abs(pt) < 0.02) continue;
            // pt = TMath::Sign (TMath::Max(0.02,TMath::Abs(pt)), pt);
            /// fout_v << "pt " << pt << endl;
            // Double_t par0 = 0.00371767;
            Double_t par0 = 0.00601267; // AZ-131020: new geometry
            // if (i234 < 2) par0 = 0.00243389;
            if (i234 < 2) par0 = 0.0013629; // AZ-131020: new geometry
            Double_t ptt = TMath::Abs(pt);
            dphi         = TMath::Sign(par0, pt) / TMath::Max(ptt, 0.01);
            // if (iPass > 0) ptt = TMath::Min (0.1, ptt);
            if (iPass > 0) ptt = TMath::Min(0.2, ptt); // AZ-131020: new geometry

            // epst = 0.002 / TMath::Power (ptt-0.01, 2./3);
            epst = 0.002 / TMath::Power(ptt - 0.015, 7. / 8); // AZ-131020: new geometry
            if (iPass == 0) epst = TMath::Min(epst, 0.01);
            // else if (ptt < 0.1) epst *= 2;
            // AZ epsl = 1e-2;//5e-3;//4, 15 /// 5e-3
            epsl = 0.003 / TMath::Power(ptt - 0.01, 2. / 3);
            if (iPass == 0) epsl = TMath::Min(epsl, 0.01);
            // else if (ptt < 0.1) epsl *= 2;
            /// epst = 5e-2;
            /// epsl = 5e-2;
         } else {
            epst = 0.05; // 0.2; // 5e-2;
            epsl = 0.01; // 0.05; //2e-2;
            if (iPass > 0) {
               epst = 0.2;                 // 0.2;
               if (iPass == 2) epst = 0.3; // AZ - 131020: new geo - larger distance to previous layer
               epsl = 0.05;                // 0.05;
            }
         }
         /// fout_v << "Cut parameters: transverse " << epst << "; longitudinal " << epsl << ";" << endl;

         TVector3 meas = cand->GetMeas();

         Double_t lastTrackL = TMath::ATan2(TMath::Sqrt(meas[0] * meas[0] + meas[1] * meas[1]), meas[2]);
         multimap<Float_t, MpdVector *>::iterator itlow = lMultimapHits[i234].lower_bound(lastTrackL - epsl);
         multimap<Float_t, MpdVector *>::iterator ittop = lMultimapHits[i234].upper_bound(lastTrackL + epsl);

         Double_t lastTrackT = TMath::ATan2(meas[1], meas[0]);
         lastTrackT -= dphi;

         multimap<Float_t, MpdVector *>::iterator ittlow = tMultimapHits[i234].lower_bound(lastTrackT - epst);
         multimap<Float_t, MpdVector *>::iterator itttop = tMultimapHits[i234].upper_bound(lastTrackT + epst);

         set<MpdVector *> l, t, intersect;
         for (multimap<Float_t, MpdVector *>::iterator itr = itlow; itr != ittop; ++itr) {
            l.insert((*itr).second);
         }
         for (multimap<Float_t, MpdVector *>::iterator itr = ittlow; itr != itttop; ++itr) {
            t.insert((*itr).second);
         }

         set_intersection(l.begin(), l.end(), t.begin(), t.end(), std::inserter(intersect, intersect.begin()));

         /// fout_v << "layer " << i234 << " long in window " << l.size() << " transverse in window " << t.size() << "
         /// intersection " << intersect.size() << " pt " << pt << endl;
         for (set<MpdVector *>::iterator itr = intersect.begin(); itr != intersect.end(); ++itr) {
            MpdVector     *track = (*itr);
            MpdItsHit5spd *hh    = (MpdItsHit5spd *)fItsHits->UncheckedAt(track->GetIndex(0));
            Int_t          trID  = ((MpdStsPoint *)fItsPoints->UncheckedAt(hh->GetTrackID()))->GetTrackID();
            if (fExact && trID != trID0) continue; // exact ID match

            MpdVector *trNew = nullptr;
            trNew            = new ((*fKHits[i234])[nKH++]) MpdVector(*track); /// was h->Layer()
            trNew->SetUniqueID(nKH - 1);
            trNew->SetNofDim(2);
            trNew->SetCosSin(0, 0);
            trNew->SetCosSin(1, 0);
            trNew->SetPrevTrackPointer(cand);
            trNew->SetCode(cand->GetCode());
            trNew->SetCode(i);

            Int_t trackID   = 0;
            Int_t MCtrackID = 0;
            if (i234 == 0) {
               Int_t f = 0;

               Int_t ll = iPass;
               if (iPass > 0) ll -= 1;

               for (int j = 0; j < 4 - ll; j++) {
                  trackID           = ((MpdItsHit5spd *)fItsHits->UncheckedAt(trNew->GetIndex(0)))->GetTrackID();
                  Int_t tempTrackID = MCtrackID;
                  MCtrackID         = ((MpdStsPoint *)fItsPoints->UncheckedAt(trackID))->GetTrackID();

                  if ((tempTrackID != MCtrackID) && (j != 0)) {
                     f = 1;
                  }

                  fout_v << "("
                         << TMath::ATan2(trNew->GetMeas()[1], trNew->GetMeas()[0])
                         /*<< " " << TMath::ATan2(v, u)*/
                         << " "
                         << TMath::ATan2(TMath::Sqrt(trNew->GetMeas()[0] * trNew->GetMeas()[0] +
                                                     trNew->GetMeas()[1] * trNew->GetMeas()[1]),
                                         trNew->GetMeas()[2])
                         << " " << trNew->GetMeas()[2]
                         << "|" //<< trNew->GetMeas()[1] << "|" << trNew->GetMeas()[2] << " " */
                                //<< " " << trNew->GetKalmanHit()->GetMeas(0)
                         << /*trackID <<*/ " " << MCtrackID << ")";
                  MpdVector *trOld = trNew;

                  trNew = trNew->GetPrevTrackPointer();
               }

               trackID           = ((MpdItsHit5spd *)fItsHits->UncheckedAt(trNew->GetIndex(0)))->GetTrackID();
               Int_t tempTrackID = MCtrackID;
               MCtrackID         = ((MpdStsPoint *)fItsPoints->UncheckedAt(trackID))->GetTrackID();

               fout_v << "(" << TMath::ATan2(trNew->GetMeas()[1], trNew->GetMeas()[0])
                      << " " /*<< TMath::ATan2(v, u) << " "*/
                      << TMath::ATan2(TMath::Sqrt(trNew->GetMeas()[0] * trNew->GetMeas()[0] +
                                                  trNew->GetMeas()[1] * trNew->GetMeas()[1]),
                                      trNew->GetMeas()[2])
                      << " " << trNew->GetMeas()[2]
                      << "|" //<< trNew->GetMeas()[1] << "|" << trNew->GetMeas()[2] << /*" " << trackID*/
                             //<< " " << trNew->GetKalmanHit()->GetMeas(0)
                      << " " << MCtrackID << ")" << endl;
               if (tempTrackID != MCtrackID) {
                  f = 1;
               }

               if (f == 1) {
                  ///  fout_v << "!!!" << endl;
               }
               count += f;
            }
         }
      }
   } // for (Int_t i234 = 1; i234 < 4;)

   fout_v << "Count of false reconstructed tracks " << count << endl;
   fout_v << "- ExtendCellTracks - done" << endl;
}

//__________________________________________________________________________
void MpdVectorFinder::GetTrackSeeds(Int_t iPass)
{
   /// Build ITS track seeds from Cell tracks
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   Int_t        nCand = 0;
   TVector3     vert(0., 0., 0.), pmom;
   MpdKalmanHit hit;
   hit.SetType(MpdKalmanHit::kFixedR);

   Int_t layEnd = 0;
   if (iPass == fNPass - 2) layEnd = 1;
   set<MpdVector *> setVecs;

   // Loop over last layers

   for (Int_t lay = 0; lay <= layEnd; ++lay) {
      Int_t nTracks = fKHits[lay]->GetEntriesFast();

      fout_v << "seed ITS tracks:" << nTracks << " fCellMap size " << fCellMap.size() << endl;

      for (Int_t itr = 0; itr < nTracks; ++itr) {
         map<Int_t, MpdVector *> mapTr;
         MpdVector              *vec = (MpdVector *)fKHits[lay]->UncheckedAt(itr);
         if (lay && setVecs.count(vec)) continue; // hit has been already considered
         mapTr[vec->GetLayer()] = vec;

         // Check if this track has not been checked already during earlier passes

         /// if (/*iPass &&*/ fCellMap.find(track5->GetCode()) != fCellMap.end()) {
         ///   cout << " Found: " << track5->GetCode() << " " << fCellMap[track5->GetCode()] << endl;
         ///   continue;
         /// }
         fCellMap.insert(pair<TString, Int_t>(vec->GetCode(), 1));
         /// fout_v << "track " << itr << "; track1 code " << track1->GetCode() << endl; /// TODO not useful

         for (Int_t jl = 0; jl < 4; ++jl) {
            vec = vec->GetPrevTrackPointer();
            if (vec == nullptr) break;
            mapTr[vec->GetLayer()] = vec;
            // Save vector from layer 1 not to double-use it
            if (layEnd > 0 && vec->GetLayer() == 1) setVecs.insert(vec);
         }
         // MpdVector *track2 = track1->GetPrevTrackPointer();
         // MpdVector *track3 = track2->GetPrevTrackPointer();
         // MpdVector *track4 = track3->GetPrevTrackPointer();
         // MpdVector *track5 = track4->GetPrevTrackPointer();

         /// fout_v << "x measure first - last: " << track1->GetMeas().X() << " " << track5->GetMeas().X() << endl;

         // AZ MpdItsKalmanTrack *track = new ((*fTrackCand)[nCand++]) MpdItsKalmanTrack(*track1, vert, iPass); ///TODO
         // vert is not used?
         /// fout_v << "trcand noftrhits " << track->GetNofTrHits() << endl; /// track->GetNofTrHits mostly = 5

         // AZ Double_t pt = EvalPt(track2,track1); /// TODO maybe some other hits should be used for pt evaluation
         Double_t pt = 0;
         // if (track5) pt = EvalPt(track5,track4,track1); /// TODO maybe some other hits should be used for pt
         // evaluation - AZ else pt = EvalPt(track4,track3,track1); //AZ AZ-241020 if (mapTr.count(4)) pt =
         // EvalPt(mapTr[4],mapTr[3],mapTr.begin()->second); /// TODO maybe some other hits should be used for pt
         // evaluation - AZ AZ-241020 else pt = EvalPt(mapTr[3],mapTr[2],mapTr.begin()->second); //AZ
         MpdVector                                *vecs[3]; // uninitialized pointer array?
         Int_t                                     iv   = 2;
         map<Int_t, MpdVector *>::reverse_iterator rmit = mapTr.rbegin();
         if (mapTr.size() == 5) ++rmit; // skip one hit
         for (; rmit != mapTr.rend(); ++rmit) {
            vecs[iv--] = rmit->second;
            if (iv < 0) break;
         }
         pt = EvalPt(vecs[2], vecs[1], vecs[0]); // AZ

         if (TMath::Abs(pt) < 0.025) continue; // AZ

         // MpdItsKalmanTrack *track = new ((*fTrackCand)[nCand++]) MpdItsKalmanTrack(*track1, vert); ///TODO vert is
         // not used?
         MpdItsKalmanTrack *track =
            new ((*fTrackCand)[nCand++]) MpdItsKalmanTrack(*(mapTr.begin()->second), vert); /// TODO vert is not used?
         fVectorCands.push_back(mapTr);

         fout_v << "GetTrackSeeds::pt " << pt << endl;

         // TVector3 v1 = track1->GetMeas();
         // TVector3 v2 = track2->GetMeas();
         // TVector3 v3 = track3->GetMeas();
         // TVector3 v4 = track4->GetMeas();

         Double_t                                  phiOut, phiIn, rOut, rIn;
         map<Int_t, MpdVector *>::reverse_iterator mit      = mapTr.rbegin();
         MpdVector                                *vecStart = mit->second;
         ++mit;
         MpdVector *vecNext = mit->second;

         TVector3 vstart = vecStart->GetMeas();
         TVector3 vnext  = vecNext->GetMeas();
         TVector3 dv     = vstart - vnext;

         phiOut = vnext.Phi();
         phiIn  = vstart.Phi();
         rOut   = vnext.Pt();
         rIn    = vstart.Pt();
         /// fout_v << v1.Pt() << " " << v2.Pt() << endl; /// 17.1.2019

         // AZ track->SetUniqueID(itr+1);
         track->SetUniqueID(nCand);
         /// cout << itr + 1 << endl;
         ///<  parametrs track  in TPC:
         ///< 0: RPhi - coordinate in R-Phi direction
         ///< 1: Z - longitudinal coordinate
         ///< 2: Phi - local azimuthal angle
         ///< 3: Theta - local dip angle (angle w.r.t. the transverse plane)
         ///< 4: q/Pt - signed inverse Pt
         track->SetPos(rIn);          //<- ispravleno!
         track->SetParam(4, 1. / pt); // q/Pt

         // TODO don't quite understand this

         // AZ track->SetParam (3, TMath::PiOver2() - track5->GetCosSin(1)); //longitudinal angL
         // if (iPass < 3) track->SetParam (3, TMath::PiOver2() - track5->GetCosSin(1)); //AZ
         if (iPass < 3) track->SetParam(3, TMath::PiOver2() - mapTr.rbegin()->second->GetCosSin(1)); // AZ
         // if (0); //AZ-241020 (iPass < 3) track->SetParam (3, TMath::PiOver2() - vecStart->GetCosSin(1)); //AZ
         else {
            // if (-track1->GetDeltaZ() >= 0) track->SetParam (3, TMath::PiOver2() - TMath::ATan(-track1->GetDeltaZ()));
            // //AZ else track->SetParam (3, TMath::PiOver2() - (TMath::Pi()+TMath::ATan(-track1->GetDeltaZ()))); //AZ
            // AZ track->SetParam (3, TMath::PiOver2() - (v2-v1).Theta());
            // track->SetParam (3, TMath::PiOver2() - (v5-v4).Theta());
            track->SetParam(3, TMath::PiOver2() - dv.Theta());
         }
         track->SetParam(2, dv.Phi()); // Phi - rough estimate

         // Adjust Phi
         Double_t bz     = FairRunAna::Instance()->GetField()->GetBz(0., 0., 0.);
         Double_t factor = 0.003 * bz / 10.; // 0.3 * 0.01 * 5kG / 10
         Double_t rCirc  = TMath::Abs(pt / factor);
         Double_t ph     = TMath::ASin(dv.Pt() / 2 / rCirc); /// v4 - v5 /// was v3 - v4 // v1 - v2
         // AZ-241020 track->SetParam (2, track->GetParam(2) + TMath::Sign(ph,pt));// - TMath::Sign(ph,pt)
         track->SetParam(2, track->GetParam(2) + TMath::Sign(ph, pt)); // - TMath::Sign(ph,pt)
         track->SetParam(1, vstart.Z());                               // Z - coordinate
         track->SetParam(0, phiIn * rIn);                              /// phiIn*rIn // length of circle arc angle phiIn
         fout_v << " trackcand num: " << itr + 1 << " track parameters: " << track->GetParam(0) << " "
                << track->GetParam(1) << " " << track->GetParam(2) << " " << track->GetParam(3) << " "
                << track->GetParam(4) << endl;

         Double_t parOut[4] = {rOut, phiOut, 0., 0.};
         Double_t parIn[4]  = {rIn, phiIn, 0., 0.};
         // EvalCovar(parOut,parIn,track,track5); /// was track4 /// asuming evalcovar works correctly
         EvalCovar(parOut, parIn, track, vecStart); /// was track4 /// asuming evalcovar works correctly
         track->SetPosNew(track->GetPos());
         track->SetParamNew(*track->GetParam());
         track->SetDirection(MpdKalmanTrack::kOutward);

         hit.SetPos(rIn + 5.0); /// track5->GetMeas().Pt() /// rIn - 0.5 // point further than 5th layer
         MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kFALSE);
         /// fout_v << "GetTrackSeeds " << hit.GetPos() << " " << track->GetPos() << endl;
         /// fout_v << nCand-1 << " " << track->GetTrackID() << endl;
         /// fout_v << track->GetHits()->GetEntriesFast() << " " << track->GetTrHits()->GetEntriesFast() << endl;

         /// TODO 11.7.2019 where does Chi2 come from here

         track->GetHits()->Clear();
         /// track->GetChi2() equals 0 here
         track->SetChi2Its(track->GetChi2()); // temporary storage
         track->SetChi2(0.);
         track->SetDirection(MpdKalmanTrack::kInward); /// was kOutward

         Int_t trackID = ((MpdItsHit5spd *)fItsHits->UncheckedAt(vecStart->GetIndex(0)))->GetTrackID();
         track->SetTrackID(((MpdStsPoint *)fItsPoints->UncheckedAt(trackID))->GetTrackID());

         if (fTrackExact) {
            // For QA
            MpdItsKalmanTrack *trE = new ((*fTrackExact)[fTrackExact->GetEntriesFast()]) MpdItsKalmanTrack(*track);
            // Int_t nh = (Bool_t)track1 + (Bool_t)track2 + (Bool_t)track3 + (Bool_t)track4 + (Bool_t)track5;
            Int_t nh = mapTr.size();
            trE->SetNofHits(nh);
            trE->SetNofIts(iPass);
            trE->SetLastLay(mapTr.begin()->first);
         }
      } // for (Int_t itr = 0; itr < nTracks;
   }    // for (Int_t lay = 0; lay <= layEnd;

   fout_v << " Number of ITS track candidates: " << nCand << endl;
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
   fout_v << "- GetTrackSeeds - done" << endl;
}

//__________________________________________________________________________
// AZ Double_t MpdVectorFinder::EvalPt(const MpdVector *track1, const MpdVector *track2)
Double_t MpdVectorFinder::EvalPt(const MpdVector *track1, const MpdVector *track2, const MpdVector *track0,
                                 Double_t *circle)
{
   /// Evaluate signed track Pt (curvature) using 3 points

   TVector3 v1 = track1->GetMeas(); /// attention!!!
   TVector3 v2 = track2->GetMeas();
   TVector3 v0;                        // AZ
   if (track0) v0 = track0->GetMeas(); // AZ - for secondaries

   TVector2 vec1(v1.X() - v0.X(), v1.Y() - v0.Y());
   TVector2 vec2(v2.X() - v0.X(), v2.Y() - v0.Y());
   TVector2 vec21    = vec1 - vec2;
   Double_t cosAlpha = vec2 * vec21 / vec2.Mod() / vec21.Mod();
   Double_t rad      = vec1.Mod() / 2. / TMath::Sin(TMath::ACos(cosAlpha));
   Double_t bz       = FairRunAna::Instance()->GetField()->GetBz(0., 0., 0.);
   Double_t factor   = 0.003 * bz / 10.; // 0.3 * 0.01 * 5kG / 10
   Double_t phi1     = vec1.Phi();
   Double_t phi2     = vec2.Phi();
   // AZ Double_t charge = phi1 - MpdKalmanFilter::Instance()->Proxim(phi1,phi2);
   Double_t charge = 0;
   if (track0)
      charge = vec21.DeltaPhi(vec2); // AZ
   else
      charge = phi1 - MpdKalmanFilter::Instance()->Proxim(phi1, phi2);
   if (track1->GetLayer() > track2->GetLayer()) charge = -charge;
   charge = TMath::Sign(1, charge);

   if (track0) {
      // cout << rad * TMath::Sin(TMath::ACos(cosAlpha)) << " " << vec1.Mod()/2 << endl;
      Double_t perp = rad * cosAlpha;
      TVector3 v31  = v1;
      v31 += v0;
      v31 *= 0.5;
      Double_t phi31 = vec1.Phi();
      phi31 -= (charge * TMath::PiOver2());
      TVector2 perpV(TMath::Cos(phi31), TMath::Sin(phi31));
      perpV *= perp;
      TVector2 cent(v31.X(), v31.Y());
      cent += perpV;
      /*
      TVector2 vv3(v1.X(),v1.Y()), vv2(v2.X(),v2.Y());
      vv3 -= cent;
      vv2 -= cent;
      cout << rad << " " << vv3.Mod() << " " << vv2.Mod() << " " << " " << cent.Mod() << " " << iq << endl;
      */
      if (circle) {
         circle[0] = cent.X();
         circle[1] = cent.Y();
         circle[2] = TMath::Abs(rad);
      }
   }
   return -factor * TMath::Abs(rad) * charge;
}

///-----------------------------------------------------------------------

Bool_t MpdVectorFinder::PropagateHelix(Double_t rlay, TVector3 pos1, Double_t *circle, Double_t &phi)
{
   // Propagate track as a helix to a circular layer

   // Transverse view - intersection of 2 circles (track and layer)
   Double_t d2 = circle[0] * circle[0] + circle[1] * circle[1]; // dist**2 from (0,0) to track center
   Double_t d  = TMath::Sqrt(d2);
   Double_t a  = rlay * rlay - circle[2] * circle[2] + d2;
   a /= (2 * d); // dist from (0,0) to cross
   Double_t perp2 = rlay * rlay - a * a;
   if (perp2 < 0) return kFALSE;
   TVector3 crossp(circle[0], circle[1], 0);
   crossp *= (a / d);
   Double_t phic = crossp.Phi();
   if (crossp.Cross(pos1).Z() < 0)
      phic += TMath::PiOver2();
   else
      phic -= TMath::PiOver2();
   TVector3 perpV(TMath::Cos(phic), TMath::Sin(phic), 0);
   perpV *= TMath::Sqrt(perp2);
   TVector3 point(crossp);
   point -= perpV;
   phi = point.Phi();
   return kTRUE;
}

///-----------------------------------------------------------------------

Double_t MpdVectorFinder::EvalPhiPt(const MpdVector *track1, const MpdVector *track2)
{
   TVector3 v1 = track1->GetMeas(); /// attention!!!
   TVector3 v2 = track2->GetMeas();
   TVector2 vec1(v1.X(), v1.Y());
   TVector2 vec2(v2.X(), v2.Y());
   TVector2 vec21    = vec1 - vec2;
   Double_t cosAlpha = vec2 * vec21 / vec2.Mod() / vec21.Mod();
   Double_t rad      = vec1.Mod() / 2. / TMath::Sin(TMath::ACos(cosAlpha));

   TVector2 temp(v1.Y(), -v1.X());
   TVector2 center = vec1 / 2. + temp * rad * cosAlpha / temp.Mod();
   temp            = vec2 - center;
   TVector2 res(temp.Y(), -temp.X());
   return TMath::ATan2(res.Y(), res.X());
}

//__________________________________________________________________________
void MpdVectorFinder::EvalCovar(Double_t *parOut, Double_t *parIn, MpdItsKalmanTrack *track, const MpdVector *track1)
{
   /// Evaluate covariance matrix for track seed

   Double_t rIn = parIn[0], phiIn = parIn[1];
   Double_t rOut = parOut[0], phiOut = parOut[1];
   Double_t xIn  = rIn * TMath::Cos(phiIn);
   Double_t yIn  = rIn * TMath::Sin(phiIn);
   Double_t xOut = rOut * TMath::Cos(phiOut);
   Double_t yOut = rOut * TMath::Sin(phiOut);
   parOut[2]     = xOut;
   parOut[3]     = yOut;
   parIn[2]      = xIn;
   parIn[3]      = yIn;

   TMatrixD ww(5, 5);
   ww(0, 0) = track1->GetErr(0) * track1->GetErr(0); // <RphiRphi> //26.11 error x
   ww(0, 0) *= 4.;                                   // extra factor of 4

   ww(1, 1) = track1->GetErr(1) * track1->GetErr(1); // <zz> //error z
   ww(1, 1) *= 4.;                                   // AZ extra factor of 8

   Double_t dx = parOut[2] - parIn[2], dy = parOut[3] - parIn[3];
   Double_t dist2  = dx * dx + dy * dy;
   Double_t sinPhi = TMath::Sin(track->GetParam(2));
   Double_t cosPhi = TMath::Cos(track->GetParam(2));
   Double_t pOut   = TMath::Cos(phiOut) * cosPhi + TMath::Sin(phiOut) * sinPhi;
   Double_t pIn    = TMath::Cos(phiIn) * cosPhi + TMath::Sin(phiIn) * sinPhi;
   ww(2, 2)        = (pOut * pOut + pIn * pIn) / dist2 * ww(0, 0); // <PhiPhi>
   // cout << ww(2,2) << " " << ww(0,0)*2/dist2 << endl;
   ww(2, 2) = ww(0, 0) * 2 / dist2; // AZ-251020
   // ww(2,2) *= 2.; // extra factor of 2
   ww(2, 2) *= 50.; // AZ

   Double_t tanThe = TMath::Tan(track->GetParam(3));
   Double_t dRad   = parOut[0] - parIn[0];
   Double_t denom  = dRad * (1. + tanThe * tanThe);
   ww(3, 3)        = ww(1, 1) * 2. / denom / denom; // <TheThe>
   ww(3, 3) *= 50;                                  // AZ

   // AZ ww(1,1) *= 8.; //AZ extra factor of 8

   ww(4, 4) = (track->GetParam(4) * 1.) * (track->GetParam(4) * 1.); // error 100%
   // ww(4,4) = 100; //AZ
   Int_t    iok   = 0;
   TMatrixD wwTmp = ww;
   MpdKalmanFilter::Instance()->MnvertLocal(ww.GetMatrixArray(), 5, 5, 5, iok);
   track->SetWeight(ww);

   TVector3       pos, posOut, pmom;
   MpdItsHit5spd *hit = (MpdItsHit5spd *)fItsHits->UncheckedAt(track1->GetIndex());
   MpdStsPoint   *p   = (MpdStsPoint *)fItsPoints->UncheckedAt(hit->GetRefIndex());
   p->Position(pos);
   p->PositionOut(posOut);
   p->Momentum(pmom);
   Double_t phi = pmom.Phi();
   // Double_t th = pmom.Theta();
   Double_t r    = (pos.Pt() + posOut.Pt()) / 2; // radius
   Double_t phi1 = pos.Phi();
   Double_t phi2 = posOut.Phi();
   phi1 += MpdKalmanFilter::Instance()->Proxim(phi1, phi2);
   phi1 /= 2;
   Double_t zzz = (pos.Z() + posOut.Z()) / 2;
   // track->SetParam(2,phi); // !!!!! exact value - just for test
   // track->SetParam(3,TMath::PiOver2()-th); // !!!!! exact value - just for test

   if (lunErr && p->GetTrackID() > 2000)
      fprintf(lunErr, "%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n",
              r * (MpdKalmanFilter::Instance()->Proxim(phi1, track->GetParam(0) / r) - phi1), track->GetParam(1) - zzz,
              MpdKalmanFilter::Instance()->Proxim(phi, track->GetParam(2)) - phi,
              TMath::PiOver2() - pmom.Theta() - track->GetParam(3), 1 / track->GetParam(4), TMath::Sqrt(wwTmp(0, 0)),
              TMath::Sqrt(wwTmp(1, 1)), TMath::Sqrt(wwTmp(2, 2)), TMath::Sqrt(wwTmp(3, 3)));
}

//__________________________________________________________________________
void MpdVectorFinder::DoTracking(Int_t iPass)
{
   /// Run Kalman tracking
   fout_v << "- DoTracking -" << endl;
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   Int_t nCand = fTrackCand->GetEntriesFast(), iok = 0; // new 05.03.14
   fout_v << "DoTracking candidates on pass " << iPass << ": " << nCand << endl;
   cout << "DoTracking candidates on pass " << iPass << ": " << nCand << endl;
   Int_t lay0 = ((MpdKalmanHit *)fKHits1->Last())
                   ->GetLayer(); // ms 06.05 /// TODO should be last layer for mpdvector? was First()

   for (Int_t i = 0; i < nCand; ++i) {
      MpdItsKalmanTrack *track = (MpdItsKalmanTrack *)fTrackCand->UncheckedAt(i); // new 05.03.14

      // if (iPass == 2) fout_v << "fTrackCand GetNofTrHits and GetNofHits " << track->GetNofTrHits() << " " <<
      // track->GetNofHits() << endl; /// 17.1.2019

      /// fout_v << " Track seed No. " << i << ", ID: " << track->GetTrackID() << ", Hits: " << track->GetNofTrHits() <<
      /// endl;
      // for (Int_t j = 0; j < track->GetNofTrHits(); ++j) {
      // MpdKalmanHit *h = (MpdKalmanHit* )track->GetTrHits()->UncheckedAt(j);
      // MpdStsHit *hh = (MpdStsHit*) fItsHits->UncheckedAt(h->GetIndex());
      // Int_t id = ((FairMCPoint*) fItsPoints->UncheckedAt(hh->GetRefIndex()))->GetTrackID();
      /// fout_v << j << " " << h->GetDist() << " " << h->GetLayer() << endl;
      //}

      /// RunKalmanFilterCell makes track->GetNofHits() from 0 to 3-5, adding them if they pass chi2 cut
      iok = RunKalmanFilterCell(track, iPass); /// removed if (fGeo) or smth like that// from cell track

      /// 26.2.20 was only (track->GetNofHits() < 3), without else
      // AZ
      // AZ if (iok == -1) {
      if ((iok == -1) || (track->GetNofHits() < 3)) {
         fTrackCand->RemoveAt(i);
         continue;
         // AZ-161020 } else if ((iPass < 3) && ((track->GetNofHits() < 4) || (track->GetNofHits() == 4) && (track->Pt()
         // < 0.1))) { } else if ((iPass < 4) && ((track->GetNofHits() < 4) || (track->GetNofHits() == 4) &&
         // (track->Pt() < 0.1))) {
      } else if ((iPass == fNPass - 3) && (track->GetNofHits() < 5)) {
         fTrackCand->RemoveAt(i);
         continue;
      } else if ((iPass < fNPass - 1) &&
                 ((track->GetNofHits() < 4) || ((track->GetNofHits() == 4) && (track->Pt() < 0.1)))) {
         fTrackCand->RemoveAt(i);
         continue;
      }

      track->SetParam(*track->GetParamNew()); // AZ
      track->SetPos(track->GetPosNew());      // AZ

      /// added 20.2.20
      /* temporary output for found secondary track hits with their layers, distance from IP and original MC track id
      if (iPass >= 3) {
        cout << "GetTrHits!";
        /// get trackCand hits id for iPass 3 and 4
        for (Int_t j = 0; j < track->GetNofTrHits(); ++j) { /// TODO maybe GetNofHits will be more representative
     MpdKalmanHit *h = (MpdKalmanHit*) track->GetTrHits()->UncheckedAt(j);
     Int_t trackID = ((MpdItsHit5spd*) fItsHits->UncheckedAt(h->GetIndex(0)))->GetTrackID();
     Int_t MCtrackID = ((MpdStsPoint*) fItsPoints->UncheckedAt(trackID))->GetTrackID();
     cout << "(" << h->GetLayer() << ", " << h->GetPos() << ", " << MCtrackID << ") ";
        }
        cout << " " << track->GetChi2() << endl;
        // cout << "GetHits!";
        // for (Int_t j = 0; j < track->GetNofHits(); ++j) { /// TODO maybe GetNofHits will be more representative
        // 	MpdKalmanHit *h = (MpdKalmanHit*) track->GetHits()->UncheckedAt(j);
        // 	Int_t trackID = ((MpdItsHit5spd*) fItsHits->UncheckedAt(h->GetIndex(0)))->GetTrackID();
        // 	Int_t MCtrackID = ((MpdStsPoint*) fItsPoints->UncheckedAt(trackID))->GetTrackID();
        // 	cout << "(" << h->GetLayer() << ", " << h->GetPos() << ", " << MCtrackID << ") ";
        // }
        // cout << endl << endl;
      }
      */

      // Propagate track to the beam line
      /// fSa stands for Stand Alone Flag. If it's value == 1 then no matching is done.
      Bool_t okk = kTRUE;
      if (fSa) okk = GoToBeamLine(track);
      if (!okk) fTrackCand->RemoveAt(i);
   } // for (Int_t i = 0; i < nCand;

   fTrackCand->Compress();
   fout_v << "DoTracking candidates after removal " << fTrackCand->GetEntriesFast() << endl;
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
}

//__________________________________________________________________________

Bool_t MpdVectorFinder::GoToBeamLine(MpdItsKalmanTrack *track)
{
   // Go to beam line

   Double_t            vert[3]   = {0.0, 0.0, 0.0};
   TString             mass2     = "0.0194797849"; // pion mass squared
   const Double_t      thick[9]  = {0.005, 0.005, 0.005, 0.07, 0.07};
   MpdKalmanGeoScheme *geoScheme = MpdKalmanFilter::Instance()->GetGeo();
   MpdKalmanHit       *hitp;
   TMatrixD            param(5, 1);
   TMatrixDSym         weight(5), pointWeight(5);

   if (track->GetDirection() == MpdKalmanTrack::kOutward) {
      // Going outward - refit inward
      track->SetDirection(MpdKalmanTrack::kInward);
      TMatrixDSym *w     = track->GetWeight();
      Int_t        nHits = track->GetNofHits(), nHits2 = nHits * nHits;
      track->SetParam(*track->GetParamNew());
      track->SetPos(track->GetPosNew());
      track->SetNode(track->GetNodeNew());
      track->SetChi2(0);

      for (Int_t i = 0; i < 5; ++i) {
         for (Int_t j = i; j < 5; ++j) {
            if (j == i)
               (*w)(i, j) /= nHits2;
            else
               (*w)(i, j) = (*w)(j, i) = 0.;
         }
      }
      //(*w)(4,4) = 1 / track->GetParamNew(4) / track->GetParamNew(4);
      Int_t ibeg = track->GetHits()->GetEntriesFast() - 1, iend = -1, iDir = -1;

      for (Int_t i = ibeg; i != iend; i += iDir) {
         hitp = (MpdKalmanHit *)track->GetHits()->UncheckedAt(i);
         if (!MpdKalmanFilter::Instance()->PropagateToHit(track, hitp, kFALSE, kTRUE)) return kFALSE;
         if (1. / TMath::Abs(track->GetParamNew(4)) < 0.025) return kFALSE;
         // MpdKalmanFilter::Instance()->PropagateToHit(track, hitp, kFALSE, kTRUE);
         Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track, hitp, pointWeight, param);
         track->SetChi2(track->GetChi2() + dChi2);
         weight = *track->GetWeight();
         weight += pointWeight;
         track->SetWeight(weight);
         track->SetParamNew(param);

         // Add scattering in Si
         Int_t lay = hitp->GetLayer();
         cout << " aaaa " << hitp->GetPos() << " " << lay << " " << dChi2 << " " << 1 / param(4, 0) << endl;
         Double_t     x0    = 9.36; // rad. length
         TMatrixDSym *cov   = track->Weight2Cov();
         Double_t     th    = track->GetParamNew(3);
         Double_t     cosTh = TMath::Cos(th);
         TVector3     norm  = geoScheme->Normal(hitp);
         // norm[2] = 0.0;
         norm.SetMag(1.0);
         MpdKalmanTrack trTmp(*track);
         trTmp.SetParam(*trTmp.GetParamNew());
         TVector3 mom3 = trTmp.Momentum3();
         // mom3[2] = 0.0;
         mom3.SetMag(1.0);
         Double_t cosPhi = norm.Dot(mom3);
         // Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, thick[lay]*3/x0/cosPhi, mass2);
         Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, 2 * thick[lay] / cosPhi, mass2);
         (*cov)(2, 2) += (angle2 / cosTh / cosTh);
         (*cov)(3, 3) += angle2;
         Int_t iok = 0;
         MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
         track->SetWeight(*cov);
      }
   } // if (track->GetDirection() == MpdKalmanTrack::kOutward)

   // AZ track->SetParam(*track->GetParamNew());
   // AZ track->SetPos(track->GetPosNew());

   Double_t    pos     = track->GetPosNew();
   TMatrixD    par     = *track->GetParamNew();
   TMatrixDSym cov     = *track->Weight2Cov();
   Double_t    leng    = track->GetLength();
   TString     nodeNew = track->GetNodeNew();

   // Go to beam pipe
   MpdKalmanHit hit;
   hit.SetType(MpdKalmanHit::kFixedR);
   hit.SetPos(fPipeR);
   Int_t iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
   // if (1./TMath::Abs(track->GetParamNew(4)) < 0.025) return kFALSE;
   if (iok != 1) {
      // Restore track
      track->SetParam(par);
      track->SetParamNew(par);
      track->SetCovariance(cov);
      track->ReSetWeight();
      track->SetPos(pos);
      track->SetPosNew(pos);
      track->SetLength(leng);
      track->SetNodeNew(nodeNew);
   } else {
      // Add multiple scattering
      // Double_t dX = 0.05 / 8.9; // 0.5 mm of Al
      Double_t     dX     = 0.1 / 35.28; // 1. mm of Be
      TMatrixDSym *pcov   = track->Weight2Cov();
      Double_t     th     = track->GetParamNew(3);
      Double_t     cosTh  = TMath::Cos(th);
      Double_t     angle2 = MpdKalmanFilter::Instance()->Scattering(track, dX);
      (*pcov)(2, 2) += (angle2 / cosTh / cosTh);
      (*pcov)(3, 3) += angle2;
      Int_t ok = 0;
      MpdKalmanFilter::Instance()->MnvertLocal(pcov->GetMatrixArray(), 5, 5, 5, ok);
      track->SetWeight(*pcov);
   }
   cov = *track->Weight2Cov();

   hit.SetPos(0.);
   hit.SetMeas(0, track->GetParam(2)); // track Phi
   // AZ-030121 pos = track->GetPosNew(); //AZ-231020
   // fout_v << i << " " << track->GetTrackID() << " " << track->GetLength() << " " <<
   // ((MpdKalmanHitR*)track->GetHits()->First())->GetLength() << endl; AZ-030121 iok = 0; //AZ-231020
   // MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
   iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
   if (iok != 1) {
      track->SetNodeNew(""); // for FindPca
      MpdKalmanFilter::Instance()->FindPca(track, vert);
   }
   // AZ-030121 track->SetNode(track->GetNodeNew()); //AZ-231020
   // AZ-030121 track->SetPos(pos); //AZ-231020
   track->SetParam(*track->GetParamNew()); // !!! track params at PCA
   return kTRUE;
}

//__________________________________________________________________________

Int_t MpdVectorFinder::RunKalmanFilterCell(MpdItsKalmanTrack *track, Int_t iPass)
{
   /// Run Kalman filter (fitter) for the hits from the cell track
   /// (might not work when propagating outward!!!)

   TString             mass2     = "0.0194797849"; // pion mass squared
   const Double_t      thick[9]  = {0.005, 0.005, 0.005, 0.07, 0.07};
   MpdKalmanGeoScheme *geoScheme = MpdKalmanFilter::Instance()->GetGeo();

   Int_t layMax = ((MpdKalmanHit *)fKHits1->First())->GetLayer();

   MpdKalmanHit            *hitOK = 0x0;
   MpdKalmanHit             hitTmp;
   MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
   // Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1; /// old coment
   // Int_t layEnd = -1, dLay = -1, layOK = -1;
   // if (trackDir == MpdKalmanTrack::kInward) { /// was kOutward
   //   layEnd = layMax + 1;
   //   dLay = 1;
   // }

   TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
   TMatrixD    param(5, 1), paramTmp(5, 1);

   Double_t posNew = 0.0; // saveZ = 0.0, saveLeng = 0.0, dChi2Min = 0.0,
   Int_t    ok     = 0;

   /// TODO what is happening here?
   /// fout_v << track->GetUniqueID() << endl;
   /*AZ-211020
   MpdVector *cellTr = nullptr;

   MpdVector *tracks[5];

   tracks[0] = (MpdVector*) fKHits[0]->UncheckedAt(track->GetUniqueID() - 1); // layer 1
   ///fout_v << track->GetUniqueID() - 1;

   /// added 19.2.2020, iPass changed to l for pass #2 purposes
   Int_t l = iPass;
   //AZ if (iPass >= 2) l = iPass - 2;
   if (iPass > 0) --l;
   if (iPass >= 3) l = iPass - 3;
   ///

   for (Int_t i = 1; i <= 4 - l; i++) { // for layers 2-5
     tracks[i] = (MpdVector*) fKHits[i]->UncheckedAt(tracks[i - 1]->GetPrevTrackPointer()->GetUniqueID());
   }
   */

   // for (Int_t lay = 4 - l; lay >= 0; lay--) {
   //  Get CellTrack from the inner layer
   // cellTr = (MpdVector*) fKHits[lay]->UncheckedAt(tracks[lay]->GetUniqueID()); // for layers 1-5

   map<Int_t, MpdVector *> &mapTr = fVectorCands[track->GetUniqueID() - 1];

   for (map<Int_t, MpdVector *>::reverse_iterator mit = mapTr.rbegin(); mit != mapTr.rend(); ++mit) {
      // for (map<Int_t,MpdVector*>::iterator mit = mapTr.begin(); mit != mapTr.end(); ++mit) {

      /// if (lay == 0)
      ///   fout_v << "pt " << EvalPt(cellTr,cellTr->GetPrevTrackPointer()) << endl;

      /// MpdKalmanHit *hit = (MpdKalmanHit*) fKHits1->UncheckedAt(cellTr->GetIndex(0));
      // AZ MpdKalmanHit *hit = cellTr->GetKalmanHit();
      MpdKalmanHit *hit = mit->second->GetKalmanHit();
      /// fout_v << hit << " " << hit->GetUniqueID() << endl;
      /// fout_v << track->GetPos() << " " << TMath::Sqrt(cellTr->GetMeas()[0] * cellTr->GetMeas()[0] +
      /// cellTr->GetMeas()[1] * cellTr->GetMeas()[1] + cellTr->GetMeas()[2] * cellTr->GetMeas()[2]) << endl;
      // Propagate to hit (if it is not very close to the track)
      if (TMath::Abs(hit->GetPos() - track->GetPosNew()) > 1.e-4) {
         Double_t leng    = track->GetLength();
         posNew           = track->GetPosNew();
         TMatrixD parNew  = *track->GetParamNew();
         TString  nodeNew = track->GetNodeNew();
         TString  curPath = track->GetNode();
         /// fout_v << track->GetLength() << " " << track << " " << hit << endl;
         ///  PropagateToHit updates track->GetLength()
         if (!MpdKalmanFilter::Instance()->PropagateToHit(track, hit, kTRUE, kTRUE)) {
            // Restore initial parameters for the failed track
            /// cout << "restore" << endl; /// TODO this is happening too often
            track->SetPosNew(posNew);
            track->SetParamNew(parNew);
            track->SetLength(leng);
            track->SetNodeNew(nodeNew);
            track->SetNode(curPath);
            ok = -1;
            break;
         }

         // Double_t step = track->GetLength() - leng;

         /// fout_v << "!!!kalmanfiltercell" << " " << hit->GetLayer() << " " << track->GetLength() << " " << leng << "
         /// " << step << " " << " " << TMath::Sqrt(hit->GetMeas(0) * hit->GetMeas(0) + hit->GetMeas(1) *
         /// hit->GetMeas(1)) << " " << track->GetPosNew() << " " << track->GetUniqueID() << endl;
      } // if (TMath::Abs(hit->GetPos()-track->GetPosNew()) > 1.e-4)

      // Exclude used hits!
      if (hit->GetFlag() < 0) continue; // new 04.03 work 15.04
      // !!! Exact ID match
      if (fExact && TrackID(hit) != track->GetTrackID()) continue;

      Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track, hit, pointWeight, param);
      // fout_v << "dChi2 " << dChi2 << endl;

      // AZ if ((TMath::Abs(dChi2) < fgkChi2Cut) || (iPass == 2)) {
      if (TMath::Abs(dChi2) < fgkChi2Cut) { /// AZ
         // if (iPass == 2) fout_v << "chi2 cut passed" << endl;
         track->GetHits()->Add(hit); /// this changes GetNofHits() result
         track->SetChi2(track->GetChi2() + dChi2);
         TMatrixDSym w = *track->GetWeight();
         w += pointWeight;
         track->SetWeight(w);
         track->SetParamNew(param);
         // Save track params at last hit
         track->SetLengAtHit(track->GetLength());
         track->SetParamAtHit(param);
         track->SetWeightAtHit(*track->GetWeight());
         track->SetPosAtHit(track->GetPosNew());
      }

      // Add scattering in Si
      Int_t lay = hit->GetLayer();
      if (1) {
         Double_t x0 = 9.36; // rad. length
         // Double_t x0 = 30420; // rad. length of Air
         TMatrixDSym *cov   = track->Weight2Cov();
         Double_t     th    = track->GetParamNew(3);
         Double_t     cosTh = TMath::Cos(th);
         TVector3     norm  = geoScheme->Normal(hit);
         // norm[2] = 0.0;
         norm.SetMag(1.0);
         MpdKalmanTrack trTmp(*track);
         trTmp.SetParam(*trTmp.GetParamNew());
         TVector3 mom3 = trTmp.Momentum3();
         // mom3[2] = 0.0;
         mom3.SetMag(1.0);
         Double_t cosPhi = norm.Dot(mom3);
         // AZ Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, 0.005*2/x0, mass2);
         // Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, thick[lay]*2/x0, mass2); //AZ-161020
         // Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, thick[lay]*3/x0, mass2); //AZ-231020
         // Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, thick[lay]*3/x0/cosPhi, mass2); //AZ-231020
         Double_t angle2 =
            MpdKalmanFilter::Instance()->Scattering(track, x0, 3 * thick[lay] / cosPhi, mass2); // AZ-231020
         (*cov)(2, 2) += (angle2 / cosTh / cosTh);
         (*cov)(3, 3) += angle2;
         Int_t iok = 0;
         MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
         track->SetWeight(*cov);
      }

      /*
      Double_t nCables = 0, x0 = 0.0116 * 2 / 9.36; // in rad. length - 116um cable per side
      // Find number of cables crossed
      TString path = gGeoManager->GetPath();

      ///fout_v << path << endl;

      if (!path.Contains("sensor") && !path.Contains("sector")) {
        /// fout_v << " !!! MpdVector::RunKalmanFilter - Outside detector !!! " << endl;
        exit(0);
      }
      Double_t v7[3] = {track->GetParamNew(0), track->GetPosNew(), track->GetParamNew(1)}, v77[3];
      gGeoManager->LocalToMaster(v7,v77);
      Double_t zTr = TMath::Abs (v77[2]); // global Z
      map<Double_t,Double_t>::iterator it;
      for (it = fCables[lay/2].begin(); it != fCables[lay/2].end(); ++it) {
        if (zTr < it->first || zTr > it->second) continue;
        ++nCables;
      }
      if (nCables) {
        x0 *= nCables;
        TMatrixDSym *cov = track->Weight2Cov();
        Double_t th = track->GetParamNew(3);
        Double_t cosTh = TMath::Cos(th);
        Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0);
        //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
        (*cov)(2,2) += (angle2 / cosTh / cosTh);
        (*cov)(3,3) += angle2;
        Int_t iok = 0;
        MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
        track->SetWeight(*cov);
      }
      */
      //} // for (Int_t lay = 0; lay < 5; ++ lay)
   } // for (map<Int_t,MpdVector*>::reverse_iterator mit = mapTr.rbegin();
   // AZ return 0;
   return ok;
}

//__________________________________________________________________________
void MpdVectorFinder::RemoveDoubles()
{
   /// Remove double tracks (keep the ones with better quality)

   Int_t ntracks = fTrackCand->GetEntriesFast(); // new 05.03.14
   fout_v << "- RemoveDoubles -" << endl;
   fout_v << " Total tracks: " << ntracks << endl;

   /// added new version 21.1.2019
   // Fill map with track quality

   /// MpdVector

   multimap<Double_t, MpdItsKalmanTrack *>           qualMap;
   multimap<Double_t, MpdItsKalmanTrack *>::iterator mit1, mit2;

   for (Int_t i = 0; i < ntracks; i++) {
      MpdItsKalmanTrack *tr = (MpdItsKalmanTrack *)fTrackCand->UncheckedAt(i);
      // AZ Double_t quality = tr->GetNofHits() + TMath::Min(tr->GetChi2(),999.0) / 1000.0;
      Double_t quality = tr->GetNofHits() + (1 - TMath::Min(tr->GetChi2(), 999.0) / 1000.0);
      qualMap.insert(pair<Double_t, MpdItsKalmanTrack *>(-quality, tr));
   }

   MpdItsKalmanTrack *tr1, *tr2;

   for (mit1 = qualMap.begin(); mit1 != qualMap.end(); ++mit1) {
      /// fout_v << mit1->first << " " << mit1->second << endl;
      tr1 = mit1->second;
      if (tr1 == nullptr) continue; // killed track
      mit2 = mit1;
      ++mit2;

      for (; mit2 != qualMap.end(); ++mit2) {
         tr2 = mit2->second;
         if (tr2 == nullptr) continue; // killed track
         bool x = AreTracksDoubles(tr1, tr2);
         /// fout_v << x << endl;
         if (!x) continue;

         /// fout_v << "tracks doubles" << endl;

         if (mit1->first < mit2->first) {
            fTrackCand->Remove(tr2);
            mit2->second = nullptr;
         } else {
            if (tr1->GetChi2() <= tr2->GetChi2()) {
               fTrackCand->Remove(tr2);
               mit2->second = nullptr;
            } else {
               fTrackCand->Remove(tr1);
               mit1->second = nullptr;
               break;
            }
         }
      }
   }

   fTrackCand->Compress();
   fNTracks = fTrackCand->GetEntriesFast();
}

//__________________________________________________________________________
Bool_t MpdVectorFinder::AreTracksDoubles(MpdItsKalmanTrack *tr1, MpdItsKalmanTrack *tr2)
{
   /// Searching common hits in 2 tracks to determine doubles
   // track1 contains fewer hits than track2
   /// MpdVector?
   MpdItsKalmanTrack *track1 = tr1, *track2 = tr2;
   if (tr1->GetNofHits() > tr2->GetNofHits()) // 12.02
      track1 = tr2, track2 = tr1;

   Int_t limCommonPoint = (track1->GetNofHits() + 1) / 2; // at least so many common hits should be found
   limCommonPoint       = TMath::Min(limCommonPoint, 2);  // AZ 2 common hits at max
   TObjArray *hits1 = track1->GetHits(), *hits2 = track2->GetHits();
   Int_t      nh1 = hits1->GetEntriesFast(), nh2 = hits2->GetEntriesFast(), nHitsCommon = 0, j = nh2 - 1;

   for (Int_t i = nh1 - 1; i >= 0; i--) {
      MpdKalmanHit *hit1 = (MpdKalmanHit *)hits1->UncheckedAt(i);

      for (; j >= 0; j--) {
         MpdKalmanHit *hit2 = (MpdKalmanHit *)hits2->UncheckedAt(j);

         // Is the hit common for two tracks?
         if (hit1 == hit2) {
            nHitsCommon++;
            if (nHitsCommon == limCommonPoint) return kTRUE; // already enough common hits
            break;
         }

         if (hit2->GetLayer() > hit1->GetLayer()) break; // already closer to beam /// was < /// what does this do?
         // if (hit2->GetLayer() < hit1->GetLayer()) break; // already closer to beam /// was < /// what does this do?
      }

      if ((nh1 - i - 1) + limCommonPoint - nHitsCommon > nh1)
         return kFALSE; // there'll be not enough common hits already
   }

   if (nHitsCommon < limCommonPoint) return kFALSE; // not too many common hits

   return kTRUE;
}

//__________________________________________________________________________
Int_t MpdVectorFinder::TrackID(MpdKalmanHit *hit, Int_t indx)
{
   /// Return track ID of the hit

   // AZ-141020 FairHit *h = (FairHit*) fItsHits->UncheckedAt(hit->GetIndex(indx));
   // AZ-141020 return ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
   MpdItsHit5spd *h = (MpdItsHit5spd *)fItsHits->UncheckedAt(hit->GetIndex(indx));
   return ((MpdStsPoint *)fItsPoints->UncheckedAt(h->GetTrackID()))->GetTrackID();
}

//__________________________________________________________________________
TVector2 MpdVectorFinder::GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit)
{
   /// Compute distance between track and hit

   Int_t lay  = hit->GetLayer();
   Int_t lay2 = lay / 2;
   // Int_t side = lay % 2;
   /// Int_t module0 = ((MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex()))->Module();
   Int_t module0 = ((MpdItsHit5spd *)fItsHits->UncheckedAt(hit->GetIndex()))->Module();

   Double_t zTr   = track->GetParamNew(1);
   Double_t zloc  = zTr + fDz[lay2];
   Int_t    modTr = Int_t(zloc / fZmod[lay2]);
   Int_t    dMod  = modTr - module0;

   Double_t dZ = 0;
   if (dMod != 0) {
      // Not in the same module0 - check Z-distance to module0 edge
      if (dMod > 0)
         dZ = zloc - (module0 + 1) * fZmod[lay2];
      else
         dZ = zloc - module0 * fZmod[lay2];
      if (TMath::Abs(dMod) > 2) return TVector2(TMath::Abs(dZ), 999.); // not in neighbour modules
   }

   // Translate transverse track position to local system
   Double_t xloc = track->GetParamNew(0) * hit->GetCosSin(0) + zTr * hit->GetCosSin(1);
   Double_t phTr = xloc / track->GetPosNew();

   Double_t phHit = hit->GetMeas(0) / fRad[lay];
   // TVector2 dist(TMath::Abs(dZ), TMath::Abs(MpdKalmanFilter::Instance()->Proxim(phTr,phHit)-phTr));
   TVector2 dist(TMath::Abs(dZ),
                 TMath::Abs(MpdKalmanFilter::Instance()->Proxim(phTr, phHit, hit->GetCosSin(0)) - phTr));
   // TVector2 dist(TMath::Abs(zTr-hit->GetZ()),
   // TMath::Abs(MpdKalmanFilter::Instance()->Proxim(track,hit)/hit->GetR()-ph));
   return dist;
}

//__________________________________________________________________________
void MpdVectorFinder::Write()
{
   /// Write

   TFile histoFile("ItsRec.root", "RECREATE");
   Writedir2current(fHistoDir);
   histoFile.Close();
}

//__________________________________________________________________________
void MpdVectorFinder::Writedir2current(TObject *obj)
{
   /// Write

   if (!obj->IsFolder())
      obj->Write();
   else {
      TDirectory *cur = gDirectory;
      TDirectory *sub = cur->mkdir(obj->GetName());
      sub->cd();
      TList *listSub = ((TDirectory *)obj)->GetList();
      TIter  it(listSub);
      while (TObject *obj1 = it()) Writedir2current(obj1);
      cur->cd();
   }
}

//__________________________________________________________________________
void MpdVectorFinder::StoreTracks()
{
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   /// Transfer tracks from fTrackCand to fTracks

   Int_t nFound = fTracks->GetEntriesFast();
   /// cout << "fTrackCand entries " << fTrackCand->GetEntriesFast() << endl;
   for (Int_t i = 0; i < fNTracks; ++i) {
      MpdItsKalmanTrack *track = (MpdItsKalmanTrack *)fTrackCand->UncheckedAt(i);
      track->Weight2Cov();

      // AZ-211020 MpdVector *cellTr = (MpdVector*) fKHits[0]->UncheckedAt(track->GetUniqueID()-1);
      MpdVector *cellTr = (MpdVector *)fVectorCands[track->GetUniqueID() - 1].begin()->second;

      /// added 16/9/2019, changing uniqueid
      // AZ track->SetUniqueID(i + 1); /// uniqueid cannot be 0
      track->SetUniqueID(nFound + 1); // AZ - 17.04.2020
      new ((*fTracks)[nFound++]) MpdItsKalmanTrack(*track);

      cout << track->GetNofTrHits() << " " << track->GetNofHits() << ", ";

      TString code   = cellTr->GetCode();
      fCellMap[code] = -1; // used hit combination
                           // fTrackCand->RemoveAt(i);
   }
   fNTracks = 0;
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
}

//__________________________________________________________________________
void MpdVectorFinder::AddHits()
{
   /// Add hit objects to tracks and compute number of wrongly assigned hits
   /// (hits with ID different from ID of starting TPC track)

   Int_t wrongCount = 0;

   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
   Int_t nFound = fTracks->GetEntriesFast(); // work 07
   // Int_t nFound = fTrackCand->GetEntriesFast();

   fout_v << "nFound " << nFound << endl;

   for (Int_t i = 0; i < nFound; ++i) {
      MpdItsKalmanTrack *track = (MpdItsKalmanTrack *)fTracks->UncheckedAt(i);
      // MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
      /// fout_v << track->GetNode() << " " << track->GetNodeNew() << endl;
      Double_t c2 = track->GetChi2();
      track->SetChi2(track->GetChi2Its());
      track->SetChi2Its(c2);
      Int_t nHits = track->GetNofHits();
      /// if(nHits < 4) continue; // candidate not write track /// was 7
      // if (nHits == 0) { fTracks->RemoveAt(i); continue; }
      TClonesArray &trHits = *track->GetTrHits();
      trHits.Delete(); // AZ GetNofTrHits() will be 0 after this
      TObjArray *hits = track->GetHits();
      SetTrackID(track); // set track ID as ID of majority of its hits
      /// fout_v << nHits << " " << trHits.GetEntriesFast() << " " << track->GetTrackID() << " " <<
      /// track->GetNofTrHits() << endl;

      Int_t nWrong = 0, motherID = track->GetTrackID(); // nMirr = 0,
      // Get track mother ID

      /// cout << "motherID " << motherID << endl;
      MpdMCTrack *mctrack = (MpdMCTrack *)fMCTracks->UncheckedAt(motherID);
      while (mctrack->GetMotherId() >= 0) {
         motherID = mctrack->GetMotherId();
         mctrack  = (MpdMCTrack *)fMCTracks->UncheckedAt(mctrack->GetMotherId());
      }

      // Int_t lastIndx = trHits.GetEntriesFast();
      for (Int_t j = 0; j < nHits; ++j) {
         MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(j);
         // AZ hit->SetUniqueID(1); // flag ITS hits
         new (trHits[j]) MpdKalmanHit(*hit);
         /// fout_v << " " << hit->GetLayer();
         Int_t nOver = hit->Index()->GetSize();
         for (Int_t iov = 0; iov < nOver; ++iov) {
            MpdItsHit5spd *h         = (MpdItsHit5spd *)fItsHits->UncheckedAt(hit->GetIndex(iov));
            Int_t          motherID1 = ((FairMCPoint *)fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
            /// fout_v << "-" << motherID1;
            // Get point mother ID
            mctrack = (MpdMCTrack *)fMCTracks->UncheckedAt(motherID1);
            while (mctrack->GetMotherId() >= 0) {
               motherID1 = mctrack->GetMotherId();
               mctrack   = (MpdMCTrack *)fMCTracks->UncheckedAt(mctrack->GetMotherId());
            }
            if (motherID1 != motherID) ++nWrong;
         }
      }
      /// if (nHits) cout << "\n Wrongs: " << nWrong << endl;

      if (nWrong > 0) {
         wrongCount++;
         ///  cout << "!!!" << endl;
      }

      track->SetNofWrong(nWrong);
      track->SetNofHits(track->GetNofTrHits()); // TPC and ITS hits
      track->SetNofIts(nHits);
      /// fout_v << nHits << " " << track->GetNofTrHits() << " " << track->GetTrackID() << " "
      ///      << track->GetChi2Its() << " " << track->GetChi2() << endl;
   }
   /// fout_v << "wrongCount: " << wrongCount << endl;
   fTracks->Compress();
   if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
}

//__________________________________________________________________________
void MpdVectorFinder::SetTrackID(MpdItsKalmanTrack *track)
{
   /// Set track ID as ID of majority of its hits

   const Int_t   idMax = 9999;
   vector<Int_t> ids(idMax);
   Int_t         nHits = track->GetNofHits(), locMax = 0, size = idMax, id0 = -1;
   TObjArray    *hits = track->GetHits();

   for (Int_t i = 0; i < nHits; ++i) {
      MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(i);
      Int_t         id  = GetHitID(hit);
      if (i == 0) id0 = id; // first hit
      if (id >= size) {
         size = id + 1;
         ids.resize(size);
      }
      ++ids[id];
      if (ids[id] > ids[locMax]) locMax = id;
   }
   if (ids[locMax] > 1)
      track->SetTrackID(locMax);
   else
      track->SetTrackID(id0);
}

//__________________________________________________________________________
Int_t MpdVectorFinder::GetHitID(MpdKalmanHit *hit)
{
   /// Get hit ID from MCPoint ID

   Int_t nOver = hit->Index()->GetSize();
   for (Int_t iov = 0; iov < nOver; ++iov) {
      MpdItsHit5spd *h         = (MpdItsHit5spd *)fItsHits->UncheckedAt(hit->GetIndex(iov));
      Int_t          motherID1 = ((FairMCPoint *)fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
      return motherID1;
   }
   return 0; // 0 points in hit?
}

//__________________________________________________________________________
void MpdVectorFinder::ExcludeHits()
{
   /// Exclude hits, already used for tracking, from consideration during the next passes
   fout_v << "- ExcludeHits -" << endl;
   Int_t nReco = fTracks->GetEntriesFast(); // work 07
   // Int_t nReco = fTrackCand->GetEntriesFast();
   fout_v << " nReco: " << nReco << endl;
   for (Int_t i = 0; i < nReco; ++i) {
      MpdItsKalmanTrack *track = (MpdItsKalmanTrack *)fTracks->UncheckedAt(i); // work 07
      // MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
      Int_t      nhitsKF = track->GetNofHits();
      TObjArray *hits    = track->GetHits();
      for (Int_t j = 0; j < nhitsKF; ++j) {
         MpdKalmanHit  *hit = (MpdKalmanHit *)hits->UncheckedAt(j);
         MpdItsHit5spd *h   = (MpdItsHit5spd *)fItsHits->UncheckedAt(hit->GetIndex());
         hit->SetFlag(-1); // work 15.04
         h->SetFlag(-1);   // work 15.04
      }
   }
}

//__________________________________________________________________________
Double_t MpdVectorFinder::Interp(Double_t angt, Int_t choice, Int_t lay)
{
   // Perform parabolic interpolation for Angl (choice = 0) and dAngt (choice = 1)

   const Int_t    np     = 7;
   const Double_t xp[np] = {0.1 / 0.18, 0.1 / 0.069, 0.1 / 0.043, 0.1 / 0.034, 0.1 / 0.017, 0.1 / 0.0086, 0.1 / 0.0044};
   // const Double_t xp[np] =   {0.05,0.125,0.2, 0.25, 0.50, 1.00, 2.00}; //3.09

   // Double_t dpp05[np] = {0.18, 0.069, 0.043, 0.034, 0.017, 0.0086, 0.0044}; //mean angt
   // Double_t err05[np] = {0.0018, 0.00032, 0.00018, 0.00014, 0.000092, 0.00007, 0.00006}; // error mean angt

   // const Double_t angs[2][np] = {{0.042, 0.019, 0.017, 0.016, 0.017, 0.019, 0.018}, //sigma angl
   //                               {0.054, 0.013, 0.0069, 0.006, 0.0036, 0.0034, 0.003}}; // sigma dangt Box(250)
   // const Double_t dangs[2][np] = {{0.0015, 0.00067, 0.00069, 0.00065, 0.00073, 0.00065, 0.00068}, // error sigma angl
   //                                {0.0027, 0.00045, 0.00027, 0.00018, 0.00014, 0.0001, 0.0001}}; //error sigma dangt
   /*
   const Double_t angls[3][np] = {{0.047, 0.027, 0.025, 0.023, 0.027, 0.026, 0.028},
              {0.04, 0.022, 0.021, 0.02, 0.021, 0.021, 0.021},
              {0.041, 0.014, 0.013, 0.012, 0.012, 0.012, 0.012}}; // sigma deltaz 2,3,4 layers
   const Double_t angts[2][np] = {{0.062, 0.013, 0.0075, 0.0058, 0.0041, 0.0039, 0.0036},
              {0.059, 0.013, 0.0074, 0.006, 0.0035, 0.0029, 0.0026}}; //sigma dangt 3,4 layers
   */

   const Double_t angls[3][np] = {{0.045, 0.026, 0.025, 0.023, 0.027, 0.027, 0.028},
                                  {0.02, 0.011, 0.01, 0.01, 0.011, 0.011, 0.01},
                                  {0.016, 0.0045, 0.0037, 0.0036, 0.0033, 0.0033, 0.0031}}; // sigma deltaz 2,3,4 layers
   const Double_t angts[2][np] = {{0.031, 0.0068, 0.0037, 0.0029, 0.0021, 0.0019, 0.0018},
                                  {0.019, 0.0039, 0.0021, 0.0016, 0.00091, 0.00059, 0.0005}}; // sigma dangt 3,4 layers

   Double_t x0 = 0.1 / angt;
   // Double_t x0 = angt;
   Int_t ok = 0, inds[3] = {np - 3, np - 2, np - 1};
   for (Int_t i = 0; i < np; ++i) {
      if (xp[i] < x0) continue;
      inds[1] = i;
      inds[0] = i - 1;
      inds[2] = i + 1;
      ok      = 1;
      break;
   }

   if (ok) {
      if (inds[0] < 0)
         for (Int_t i = 0; i < 3; ++i) ++inds[i];
      if (inds[2] == np)
         for (Int_t i = 0; i < 3; ++i) --inds[i];
   }
   Double_t a{0.}, b{0.}, c{0.}, dy01{0.}, dy02{0.};

   // Parabolic interpolation
   Double_t dx01 = xp[inds[1]] - xp[inds[0]], dx02 = xp[inds[2]] - xp[inds[0]];
   Double_t x01 = xp[inds[1]] + xp[inds[0]], x02 = xp[inds[2]] + xp[inds[0]];
   // Double_t dy01 = angs[choice][inds[1]] - angs[choice][inds[0]];
   // Double_t dy02 = angs[choice][inds[2]] - angs[choice][inds[0]];

   if (choice == 0) {
      dy01 = angls[lay][inds[1]] - angls[lay][inds[0]];
      dy02 = angls[lay][inds[2]] - angls[lay][inds[0]];
   }
   if (choice == 1) {
      dy01 = angts[lay][inds[1]] - angts[lay][inds[0]];
      dy02 = angts[lay][inds[2]] - angts[lay][inds[0]];
   }

   Double_t slope = dy01 / dx01;
   a              = dy02 - slope * dx02;
   a /= (dx02 * x02 - dx02 * x01);

   b = slope - a * x01;
   if (choice == 0) {
      // c = angs[choice][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
      c = angls[lay][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
   }
   if (choice == 1) {
      // c = angs[choice][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
      c = angts[lay][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
   }
   // for (Int_t i = 0; i < 3; ++i) cout << xp[inds[i]] << " " << angs[choice][inds[i]] << " ";
   // cout << endl << x0 << " " << a * x0 * x0 + b * x0 + c << endl;
   return a * x0 * x0 + b * x0 + c;
}

//__________________________________________________________________________
void MpdVectorFinder::GetShortTracks()
{
   /// Collect remaining tracks with small number of hits.

   map<TString, Int_t>::iterator it;
   const Int_t                   nLays = 5;
   Int_t                         inds[nLays];
   for (it = fCellMap.begin(); it != fCellMap.end(); ++it) {
      // cout << it->first << " " << it->second <<endl;
      if (it->second < 0) continue;
      TString code = it->first, cc;
      Int_t   ibeg = 0, leng = -1;
      // cout << code << endl;
      for (Int_t i = 0; i < nLays; ++i) {
         leng    = code.Index("-", leng);
         cc      = code(ibeg, leng - ibeg);
         inds[i] = cc.Atoi();
         ibeg    = leng + 1;
         ++leng;
      }
   }

   // Int_t left = 0;
   /// fout_v << " Leftovers: " << endl;
   for (it = fCellMap.begin(); it != fCellMap.end(); ++it) {
      /// if (it->second > 0) fout_v << ++left << " " << it->first << endl;
   }
}

//__________________________________________________________________________

Bool_t MpdVectorFinder::Refit(MpdItsKalmanTrack *track, Double_t mass, Int_t charge)
{
   /// Refit track in ITS+TPC using track hits (toward beam line) for some
   /// particle mass and charge hypothesis

   if (fTpcKF == nullptr) fTpcKF = (MpdTpcKalmanFilter *)FairRun::Instance()->GetTask("TPC Kalman filter");
   if (fTpcKF == nullptr) {
      cout << " !!! TPC Kalman Filter was not activated !!! " << endl;
      exit(1);
   }
   MpdKalmanGeoScheme *geoScheme = MpdKalmanFilter::Instance()->GetGeo();
   MpdKalmanHit        hTmp;

   if (track->GetNofTrHits() > track->GetNofIts()) {
      // Refit inside TPC
      // MpdTpcKalmanTrack *tpc = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(track->GetUniqueID()-1);
      track->SetUniqueID(0); // reset ITS flag
      // tr.GetWeightAtHit()->Print();
      // tr.GetParamAtHit()->Print();
      fTpcKF->Refit(track, mass, charge); // this mixes up hit order (ITS and TPC) !!!
      // tr.GetWeight()->Print();
      // tr.GetParam()->Print();
      if (track->GetNofIts() == 0) return kTRUE;

      hTmp.SetType(MpdKalmanHit::kFixedR);
      hTmp.SetPos(track->GetPos());
      track->SetParamNew(*track->GetParam());
      track->SetPos(track->GetPosNew());
      track->ReSetWeight();
      // TMatrixDSym w = *tr.GetWeight(); // save current weight matrix
      // MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,kFALSE);
      MpdKalmanFilter::Instance()->PropagateParamToHit(track, &hTmp, kFALSE);
      // tr.SetWeight(w); // restore original weight matrix (near TPC inner shell)
      track->SetDirection(MpdKalmanTrack::kInward);
      track->SetUniqueID(1);
   } else {
      track->SetDirection(MpdKalmanTrack::kOutward);
      hTmp.SetType(MpdKalmanHit::kFixedR);
      hTmp.SetPos(track->GetPos());
      track->SetPos(track->GetPosNew());
      track->SetParamNew(*track->GetParam());
      Bool_t ok = MpdKalmanFilter::Instance()->PropagateParamToHit(track, &hTmp, kFALSE);
      track->ReSetWeight();
      ok = MpdKalmanFilter::Instance()->Refit(track, -1, 0);
      if (!ok) return ok;
      hTmp.SetPos(((MpdKalmanHit *)track->GetTrHits()->First())->GetDist() + 1.0);
      MpdKalmanFilter::Instance()->PropagateParamToHit(track, &hTmp, kFALSE);
      track->SetDirection(MpdKalmanTrack::kInward);
      track->SetPos(track->GetPosNew());
      track->SetParam(*track->GetParamNew());

      Int_t nHits = track->GetNofHits();
      nHits *= nHits;
      TMatrixDSym *www = track->GetWeight();

      for (Int_t i = 0; i < 5; ++i) {
         for (Int_t j = i; j < 5; ++j) {
            if (j == i)
               (*www)(i, j) /= nHits;
            else
               (*www)(i, j) = (*www)(j, i) = 0.;
         }
      }
   }

   Int_t ntot = track->GetNofTrHits(), nhits = TMath::Min(15, ntot);
   Int_t indx0 = ntot - nhits;

   TMatrixD    param(5, 1);
   TMatrixDSym weight(5), pointWeight(5);
   TString     mass2 = "";
   mass2 += mass * mass;
   Double_t       vert[3]  = {0.0, 0.0, 0.0};
   const Double_t thick[9] = {0.005, 0.005, 0.005, 0.07, 0.07};

   for (Int_t ih = indx0; ih < ntot; ++ih) {
      MpdKalmanHit *hit = (MpdKalmanHit *)track->GetTrHits()->UncheckedAt(ih);
      if (hit->GetUniqueID() == 0) continue; // skip TPC hits
      MpdKalmanFilter::Instance()->PropagateToHit(track, hit, kFALSE,
                                                  kTRUE); // do not adjust track length for the moment
      Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track, hit, pointWeight, param);
      // cout << ih << " " << hit->GetPos() << " " << hit->GetUniqueID() << " " << dChi2 << "\n";

      track->SetChi2(track->GetChi2() + dChi2);
      weight = *track->GetWeight();
      weight += pointWeight;
      track->SetWeight(weight);
      track->SetParamNew(param);

      // Add multiple scattering in the sensor
      //*
      Double_t x0   = 9.36; //, step = 0.005 * 4.0; // rad. length
      TVector3 mom3 = track->Momentum3();
      mom3.SetMag(1.0);
      TVector3 norm   = geoScheme->Normal(hit);
      Double_t cosPhi = norm.Dot(mom3);

      TMatrixDSym *cov   = track->Weight2Cov();
      Double_t     th    = track->GetParamNew(3);
      Double_t     cosTh = TMath::Cos(th);
      // Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, thick[lay]/cosPhistep, mass2);
      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, 4 * thick[hit->GetLayer()] / cosPhi, mass2);
      (*cov)(2, 2) += (angle2 / cosTh / cosTh);
      (*cov)(3, 3) += angle2;
      Int_t iok = 0;
      MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
      track->SetWeight(*cov);
      //*/
   }

   //------ Propagate track to the beam line
   track->SetParam(*track->GetParamNew());
   track->SetPos(track->GetPosNew());
   Double_t    pos     = track->GetPos();
   TMatrixD    par     = *track->GetParam();
   TMatrixDSym cov     = *track->Weight2Cov();
   Double_t    leng    = track->GetLength();
   TString     nodeNew = track->GetNodeNew();
   // cout << " 1: " << nodeNew << ", " << track->GetNode() << " " << track->GetHits()->GetEntriesFast() << endl;

   // Go to beam pipe
   hTmp.SetPos(fPipeR);
   // hTmp.SetPos(2.9); // !!! FIXME
   Int_t iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hTmp, kFALSE);
   if (iok != 1) {
      // Restore track
      track->SetParam(par);
      track->SetParamNew(par);
      track->SetCovariance(cov);
      track->ReSetWeight();
      track->SetPos(pos);
      track->SetPosNew(pos);
      track->SetLength(leng);
      // track->SetNode(node);
      track->SetNodeNew(nodeNew);
   } else {
      // Add multiple scattering
      // Double_t dX = 0.05 / 8.9; // 0.5 mm of Al
      Double_t     dX     = 0.1 / 35.28; // 1. mm of Be
      TMatrixDSym *pcov   = track->Weight2Cov();
      Double_t     th     = track->GetParamNew(3);
      Double_t     cosTh  = TMath::Cos(th);
      Double_t     angle2 = MpdKalmanFilter::Instance()->Scattering(track, dX);
      (*pcov)(2, 2) += (angle2 / cosTh / cosTh);
      (*pcov)(3, 3) += angle2;
      Int_t ok = 0;
      MpdKalmanFilter::Instance()->MnvertLocal(pcov->GetMatrixArray(), 5, 5, 5, ok);
      track->SetWeight(*pcov);
   }
   track->Weight2Cov(); // rebuild covar matrix (it was overwritten by weight above)

   hTmp.SetPos(0.);
   hTmp.SetMeas(0, track->GetParam(2)); // track Phi
   iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hTmp, kFALSE);
   if (iok != 1) {
      track->SetNodeNew(""); // for FindPca
      MpdKalmanFilter::Instance()->FindPca(track, vert);
   }
   track->SetParam(*track->GetParamNew()); // !!! track params at PCA
   return kTRUE;
}

//__________________________________________________________________________

ClassImp(MpdVectorFinder);