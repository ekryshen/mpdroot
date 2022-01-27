/////////////////////////////////////////////////////////////
//
//  MpdStsHitProducerV1
//
//  Filler of MpdStsHit
//
//
//
///////////////////////////////////////////////////////////////

#include "TClonesArray.h"

#include "MpdStsGeoPar.h"
#include "MpdStsHitProducerV1.h"
#include "MpdStsHit.h"
#include "MpdStsPoint.h"
#include "FairDetector.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TVector3.h"

// -----   Default constructor   -------------------------------------------
MpdStsHitProducerV1::MpdStsHitProducerV1(const char *fileGeo) : FairTask("Ideal STS hit Producer")
{
   fFileGeo = fileGeo;
   eneThr   = 0.001; // Energy threshold for STS
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdStsHitProducerV1::~MpdStsHitProducerV1() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------
InitStatus MpdStsHitProducerV1::Init()
{

   cout << " INITIALIZATION *********************" << endl;

   // FairDetector::Initialize();
   // FairRun* sim = FairRun::Instance();
   // FairRuntimeDb* rtdb=sim->GetRuntimeDb();

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdStsHitProducerV1::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array
   fPointArray = (TClonesArray *)ioman->GetObject("StsPoint");
   if (!fPointArray) {
      cout << "-W- MpdStsHitProducerV1::Init: "
           << "No StsPoint array!" << endl;
      return kERROR;
   }

   // Create and register output array
   fDigiArray = new TClonesArray("MpdStsHit");

   ioman->Register("StsHit", "Sts", fDigiArray, kTRUE);
   // ioman->Register("StsHit","Sts",fDigiArray,kFALSE);
   CreateStructure();

   cout << "-I- MpdStsHitProducerV1: Intialization successfull" << endl;

   return kSUCCESS;
}
// -------------------------------------------------------------------------

// -----   Private method SetParContainers   -------------------------------
void MpdStsHitProducerV1::SetParContainers()
{

   // Get run and runtime database
   FairRunAna *run = FairRunAna::Instance();
   if (!run) Fatal("SetParContainers", "No analysis run");

   FairRuntimeDb *db = run->GetRuntimeDb();
   if (!db) Fatal("SetParContainers", "No runtime database");

   // Get STS geometry parameter container
   db->getContainer("MpdStsGeoPar");
}

// -----   Public method Exec   --------------------------------------------
void MpdStsHitProducerV1::Exec(Option_t *opt)
{

   // cout << " DIGI EXECUTION *********************" << endl;
   static Int_t eventCounter = 0, first = 1;

   // Detector parameters
   // static Double_t sterCos[2] = {-0.0075, 0.0275}, sterSin[2] = {0}; // stereo-angles in rad
   static Double_t sterCos[2] = {-0.131, 0.131}, sterSin[2] = {0}; // stereo-angles in rad
   // static Double_t sterCos[2] = {-0., 0.}, sterSin[2] = {0}; // stereo-angles in rad
   const Double_t length = 7.3;    // cm // length of sensor
   const Double_t height = 4.0;    // cm // height of sensor
   const Double_t pitch  = 0.0095; // cm //
   const Int_t    nstrip = 768;    // number of strips in sensor

   if (first) {
      first = 0;
      for (Int_t i = 0; i < 2; ++i) {
         sterSin[i] = TMath::Sin(sterCos[i]); // sin
         sterCos[i] = TMath::Cos(sterCos[i]); // cos
      }
   }

   cout << " StsHitProducerV1 event " << ++eventCounter << endl;

   // Reset output array
   fDigiArray->Clear();

   // Declare some variables
   MpdStsPoint *       point = NULL;
   map<Int_t, Float_t> fTrackEnergy;
   fTrackEnergy.clear();
   map<Int_t, Float_t>::const_iterator p;

   // Loop over MCPoints
   Int_t nPoints = fPointArray->GetEntriesFast();
   cout << nPoints << endl;
   for (Int_t iPoint = 0; iPoint < nPoints; ++iPoint) {
      point = (MpdStsPoint *)fPointArray->At(iPoint);

      Int_t detId = point->GetDetectorID();
      Int_t lay   = ((detId >> 1) & 15);
      Int_t lad   = ((detId >> 5) & 31);
      Int_t det   = (detId >> 10);

      TVector3 pos30;
      for (Int_t j = 0; j < 2; ++j) {

         Double_t pos[3] = {0};
         TVector3 pos3;
         if (j == 0) {
            // One detector side
            point->Position(pos3);
            pos30 = pos3;
         } else {
            // The other side
            point->PositionOut(pos3);
            TVector3 dPos = pos3 - pos30;
            Double_t mag  = dPos.Mag();
            if (mag > 1.e-7) dPos *= (1. / dPos.Mag());
            pos3 -= 0.001 * dPos; // small step back
         }
         pos3.GetXYZ(pos); // global coordinate

         TGeoNode *node = gGeoManager->FindNode(pos[0], pos[1], pos[2]);
         if (node == 0x0) {
            cout << " !!! Missing node " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
            exit(0);
         }
         cout << gGeoManager->GetPath() << " " << detId + j << endl;
         Double_t loc[3] = {0.};
         gGeoManager->MasterToLocal(pos, loc);

         // Two hits per point (for two sides)
         MpdStsHit *hit = AddHit(iPoint, detId + j, j);

         Double_t x = loc[0], y = loc[1], z = loc[2];     // local coordinate
         Double_t xloc = x * sterCos[j] + z * sterSin[j]; // rotated coordinate

         Double_t xpitch = xloc / pitch;
         Int_t    strip  = Int_t(xpitch + (nstrip / 2) + 1);
         if (strip > nstrip) strip = nstrip;
         if (strip < 1) strip = 1;
         hit->SetLocalX(xloc);
         hit->SetStrip(strip);

         //*********************************************************************

         cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << endl;
         cout << loc[0] << " " << loc[1] << " " << loc[2] << " " << endl;

         cout << "xloc=" << x << " xrot=" << xloc << " xpitch=" << xpitch << " strip=" << strip << endl;

         //*********************************************************************
      }
      // cout << point->GetDetectorID() << " " << lay << " " << point->GetZ() << " " << module << " " << detId << endl;
      // fTrackEnergy[point->GetDetectorID()] += point->GetEnergyLoss();
   }
}
// -------------------------------------------------------------------------

// -----   Public method Create Structure   --------------------------------
void MpdStsHitProducerV1::CreateStructure() {}

// -----   Private method AddDigi   --------------------------------------------
MpdStsHit *MpdStsHitProducerV1::AddHit(Int_t trackID, Int_t detID, Int_t side)
{
   // It fills the MpdStsHit category

   // cout << "MpdStsHitProducer: track " << trackID << " evt " << eventID << " sec " << sec << " plane " << pla << "
   // strip " << strip << "box " << box << " tube " << tub << endl;
   MpdStsPoint * point = (MpdStsPoint *)fPointArray->UncheckedAt(trackID);
   TClonesArray &clref = *fDigiArray;
   Int_t         size  = clref.GetEntriesFast();
   TVector3      pos;
   if (side == 0)
      point->Position(pos); // one detector side
   else
      point->PositionOut(pos); // the other side
   return new (clref[size]) MpdStsHit(detID, pos, TVector3(0., 0., 0.), point->GetEnergyLoss(), trackID, 0);
}
// ----

ClassImp(MpdStsHitProducerV1);
