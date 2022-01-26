// -------------------------------------------------------------------------
// -----             MpdStsHitProducerNew2 source file                 -----
// -----                      29.04.2012 VK                            -----
// -------------------------------------------------------------------------

//  Independant strip readout for each sector
//  Geometry: its_sec_cables.geo ( Three types of sectors consist of one, two or three sensors)

#include "TClonesArray.h"

#include "MpdStsGeoPar.h"
#include "MpdStsHitProducerNew2.h"
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
MpdStsHitProducerNew2::MpdStsHitProducerNew2(const char *fileGeo) : FairTask("Ideal STS hit Producer")
{
   fFileGeo = fileGeo;
   eneThr   = 0.001; // Energy threshold for STS
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdStsHitProducerNew2::~MpdStsHitProducerNew2() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------
InitStatus MpdStsHitProducerNew2::Init()
{

   cout << " INITIALIZATION *********************" << endl;

   // FairDetector::Initialize();
   // FairRun* sim = FairRun::Instance();
   // FairRuntimeDb* rtdb=sim->GetRuntimeDb();

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdStsHitProducerNew2::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array
   fPointArray = (TClonesArray *)ioman->GetObject("StsPoint");
   if (!fPointArray) {
      cout << "-W- MpdStsHitProducerNew2::Init: "
           << "No StsPoint array!" << endl;
      return kERROR;
   }

   // Create and register output array
   fDigiArray = new TClonesArray("MpdStsHit");

   ioman->Register("StsHit", "Sts", fDigiArray, kTRUE);
   // ioman->Register("StsHit","Sts",fDigiArray,kFALSE);
   CreateStructure();

   cout << "-I- MpdStsHitProducerNew2: Intialization successfull" << endl;

   return kSUCCESS;
}
// -------------------------------------------------------------------------

// -----   Private method SetParContainers   -------------------------------
void MpdStsHitProducerNew2::SetParContainers()
{

   // Get run and runtime database
   FairRunAna *run = FairRunAna::Instance();
   if (!run) Fatal("SetParContainers", "No analysis run");

   FairRuntimeDb *db = run->GetRuntimeDb();
   if (!db) Fatal("SetParContainers", "No runtime database");

   // Get STS geometry parameter container
   db->getContainer("MpdStsGeoPar");
}
// -------------------------------------------------------------------------

// -----   Public method Exec   --------------------------------------------
void MpdStsHitProducerNew2::Exec(Option_t *opt)
{

   cout << "-I- MpdStsHitProducerNew2: Execution *********************" << endl;
   fDigiArray->Delete();

   static Int_t eventCounter = 0, first = 1;

   // Detector parameters

   // static Double_t sterCos[2] = {-0.131, 0.131}; // stereo-angles in rad (7.5, -7.5)
   // static Double_t sterCos[2] = {0.131, -0.131}; // stereo-angles in rad (7.5, -7.5)
   static Double_t sterCos[2] = {0.131, 0.}; // stereo-angles in rad (7.5, 0)
   static Double_t sterSin[2] = {0., 0.};
   static Double_t sterTan[2] = {0., 0.};
   static Double_t par1[100] = {0.}, par2[100] = {0.}, par3[100] = {0.};
   const Double_t  length_1  = 6.2;     // cm, length of sector#1
   const Double_t  length_2  = 12.4;    // cm, length of sensor#2
   const Double_t  length_3  = 18.6;    // cm, length of sensor#3
   const Double_t  width     = 6.2;     // cm, height of sensor
   const Double_t  pitch     = 0.00605; // cm, readout pitch
   const Int_t     nChannels = 1024;    // number of readout channels in sensor

   Int_t    nStrips, nShortStrips;
   Double_t xShift;

   cout << " StsHitProducerNew2 event " << ++eventCounter << endl;

   if (first) {
      first = 0;
      for (Int_t i = 0; i < 2; ++i) {
         sterSin[i] = TMath::Sin(sterCos[i]);  // sin
         sterCos[i] = TMath::Cos(sterCos[i]);  // cos
         sterTan[i] = sterSin[i] / sterCos[i]; // tan
      }
      nStrips =
         Int_t(nChannels + length_1 / pitch * sterTan[0]); // Total number of strips in sector#1 including short ones
      nShortStrips = Int_t(length_1 / pitch *
                           sterTan[0]); // Number of short strips in sector #1 not crossing readout side (extra ships)

      // Define sensitive volumes
      FairRuntimeDb *rtdb      = FairRun::Instance()->GetRuntimeDb();
      MpdStsGeoPar * geoPar    = (MpdStsGeoPar *)rtdb->getContainer("MpdStsGeoPar");
      TObjArray *    sensNodes = geoPar->GetGeoSensitiveNodes();
      cout << "Number of sensors: " << sensNodes->GetEntriesFast() << endl;
      Int_t nSensors = sensNodes->GetEntriesFast();
      for (Int_t i = 0; i < nSensors; ++i) {
         FairGeoNode *sensVol = (FairGeoNode *)(sensNodes->At(i));
         TArrayD *    params  = sensVol->getParameters();
         par1[i]              = params->At(0); // half length of sensor
         par2[i]              = params->At(1); // half tickness of sensor
         par3[i]              = params->At(2); // half width of sensor

         cout << " *** STS sensitive volume: " << sensVol->GetName() << " " << params->At(0) << " " << params->At(1)
              << " " << params->At(2) << endl;
      }
   }

   // Reset output array
   fDigiArray->Clear();

   // Declare some variables
   MpdStsPoint *point = NULL;

   // Loop over MCPoints
   Int_t nPoints = fPointArray->GetEntriesFast();
   cout << "********* Number of STS points: " << nPoints << endl;
   for (Int_t iPoint = 0; iPoint < nPoints; ++iPoint) {
      point = (MpdStsPoint *)fPointArray->UncheckedAt(iPoint);

      // cout << "# of current STS point: " << iPoint+1 << endl;

      Int_t detId      = point->GetDetectorID();
      Int_t sector     = ((detId >> 1) & 31);
      Int_t sectorType = ((detId >> 14) & 3);
      Int_t ladder     = ((detId >> 6) & 31);
      Int_t layer      = ((detId >> 11) & 7);
      Int_t side       = (detId & 1);

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
            if (dPos.Mag() > 0.) dPos.SetMag(1.); // 1 cm
            pos3 -= 0.001 * dPos;                 // small step back
         }
         pos3.GetXYZ(pos); // global coordinate

         TGeoNode *node = gGeoManager->FindNode(pos[0], pos[1], pos[2]);
         if (node == 0x0) {
            cout << " !!! Missing node " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
            exit(0);
         }

         // cout << gGeoManager->GetPath() << " DetectorID = " << detId+j << endl;
         // cout << " Layer: "<< layer << " Ladder: " << ladder << " Sector: " << sector
         //         << " Sector type: " << sectorType << " Side: " << side+j << " " <<endl;

         Double_t loc[3] = {0.};
         gGeoManager->MasterToLocal(pos, loc);

         // Two hits per point (for two sides)
         MpdStsHit *hit = AddHit(iPoint, detId + j, j);
         hit->SetUniqueID(2); // AZ - flag for sector layout
         hit->SetTrackID(point->GetTrackID());

         Double_t x = loc[0], y = loc[1], z = loc[2]; // local coordinate
         Double_t xReadOutSide;                       // X of projection point along strip to readout edge
         Int_t    iChannel, iStripSector;

         if (j == 0) { // Front side

            if (sectorType == 1) xReadOutSide = x + (0.5 * length_1 - z) * sterTan[0] + 0.5 * width;
            if (sectorType == 2) xReadOutSide = x + (0.5 * length_2 - z) * sterTan[0] + 0.5 * width;
            if (sectorType == 3) xReadOutSide = x + (0.5 * length_3 - z) * sterTan[0] + 0.5 * width;
         }

         else { // Back side

            if (sectorType == 1)
               xReadOutSide = x + (0.5 * length_1 - z) * sterTan[1] + 0.5 * width - length_1 * sterTan[1];
            if (sectorType == 2)
               xReadOutSide = x + (0.5 * length_2 - z) * sterTan[1] + 0.5 * width - length_2 * sterTan[1];
            if (sectorType == 3)
               xReadOutSide = x + (0.5 * length_3 - z) * sterTan[1] + 0.5 * width - length_3 * sterTan[1];
         }

         if ((sectorType == 1 && (xReadOutSide < 0. || xReadOutSide > width + length_1 * sterTan[0])) ||
             (sectorType == 2 && (xReadOutSide < 0. || xReadOutSide > width + length_2 * sterTan[0])) ||
             (sectorType == 3 && (xReadOutSide < 0. || xReadOutSide > width + length_3 * sterTan[0]))) {

            cout << "************************ Point is outside sector: xReadOutSide = " << xReadOutSide << endl;
            cout << "Local coordinates of point in sensor: " << x << " " << y << " " << z << " " << endl;
            cout << gGeoManager->GetPath() << endl;
            xReadOutSide = -pitch;
         }

         iStripSector = Int_t(xReadOutSide / pitch + 1);
         iChannel     = iStripSector % nChannels;

         Double_t xlocSensor = x * sterCos[j] + z * sterSin[j]; // rotated coordinate

         // hit->SetLocalX(xReadOutSide);
         hit->SetLocalX(xlocSensor); // AZ
         hit->SetStrip(iStripSector);

         // hit->SetStrip(iChannel);

         //*********************************************************************

         // cout << "Global coordinates of point in sensor: " << pos[0] << " " << pos[1] << " " << pos[2] << " " <<
         // endl; cout << "Local coordinates of point in sensor: " << loc[0] << " " << loc[1] << " " << loc[2] << " " <<
         // endl; cout << "X local along readout side: " << xReadOutSide << endl;
         // cout << "Strip's number in sector: " << iStripSector << " Channel's number in sector: " << iChannel<< endl;

         //*********************************************************************
      }
      // cout << point->GetDetectorID() << " " << lay << " " << point->GetZ() << " " << sector << " " << detId << endl;
      // fTrackEnergy[point->GetDetectorID()] += point->GetEnergyLoss();
   }
}
// -----------------------------------------------------------------------------

// -----   Public method Create Structure   ------------------------------------
void MpdStsHitProducerNew2::CreateStructure() {}
// -----------------------------------------------------------------------------

// -----   Private method AddDigi   --------------------------------------------
MpdStsHit *MpdStsHitProducerNew2::AddHit(Int_t trackID, Int_t detID, Int_t side)
{
   // It fills the MpdStsHit category
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
// -----------------------------------------------------------------------------

ClassImp(MpdStsHitProducerNew2)
