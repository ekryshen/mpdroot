// -------------------------------------------------------------------------
// -----             MpdItsHitProducer5spd source file                     -----
// -----                      13.12.2016 VK                            -----
// -------------------------------------------------------------------------

//  Independant pixel readout for each sector
//  Geometry: its_5spd.geo
//  spd: one type of sector consisting of one sensor (15*30 mm)

#include "MpdItsHitProducer5spd.h"
#include "MpdStsGeoPar.h"
#include "MpdItsHit5spd.h"
#include "MpdStsPoint.h"

#include "FairDetector.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include <TRandom.h>
#include "TVector3.h"

// -----   Default constructor   -------------------------------------------
MpdItsHitProducer5spd::MpdItsHitProducer5spd()
   : FairTask("Ideal ITS hit Producer5spd"), fPointArray(nullptr), fDigiArray(nullptr)
{
}

// -----   Default constructor   -------------------------------------------
/*
MpdItsHitProducer5spd::MpdItsHitProducer5spd(const char* fileGeo) :
  FairTask("Ideal ITS hit Producer_5spd") {
  fFileGeo=fileGeo;
  eneThr = 0.001; // Energy threshold for ITS
  fDigiArray = nullptr;
}
*/

// -----   Destructor   ----------------------------------------------------
// MpdItsHitProducer5spd::~MpdItsHitProducer5spd() { delete fDigiArray; fDigiArray = nullptr; }
MpdItsHitProducer5spd::~MpdItsHitProducer5spd() {}

// -----   Public method Init   --------------------------------------------
InitStatus MpdItsHitProducer5spd::Init()
{

   cout << "-I- MpdItsHitProducer5spd: Intialization started" << endl;

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdItsHitProducer5spd::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array
   fPointArray = (TClonesArray *)ioman->GetObject("StsPoint");
   if (!fPointArray) {
      cout << "-W- MpdItsHitProducer5spd::Init: "
           << "No ItsPoint array!" << endl;
      return kERROR;
   }

   // Create and register output array
   fDigiArray = new TClonesArray("MpdItsHit5spd");

   // ioman->Register("StsHit","Sts",fDigiArray,kTRUE);
   ioman->Register("StsHit", "Sts", fDigiArray, kFALSE);
   CreateStructure();

   cout << "-I- MpdItsHitProducer5spd: Intialization finished successfully" << endl;

   return kSUCCESS;
}

// -----   Private method SetParContainers   -------------------------------
void MpdItsHitProducer5spd::SetParContainers()
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
void MpdItsHitProducer5spd::Exec(Option_t *opt)
{

   cout << "-I- MpdItsHitProducer5spd: Execution started" << endl;
   fDigiArray->Delete();

   static Int_t eventCounter = 0, first = 1;
   // static Double_t par1[9] = {0}, par2[9] = {0}, par3[9] = {0},
   static Double_t spd_length = 0, spd_width = 0; // Detector parameters
   Int_t           detId = 0;                     //, sector = 0, ladder = 0, layer = 0;

   //------ spd parameters
   // const Double_t spd_length = 3.0;       // cm, length of spd sensor
   // const Double_t spd_width = 1.5;        // cm, width of spdsensor
   const Int_t nCol = 1024; // number of colomns in spd matrix
   const Int_t nRow = 512;  // number of rows in spd matrix

   cout << " *** ItsHitProducer5spd event " << ++eventCounter << endl;

   /*
   if (first) {
     first = 0;
     //------ Define sensitive volumes
     FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
     MpdStsGeoPar *geoPar = (MpdStsGeoPar*) rtdb->getContainer("MpdStsGeoPar");
     TObjArray* sensNodes = geoPar->GetGeoSensitiveNodes();
     cout << "Number of sensors: "<< sensNodes->GetEntriesFast() << endl;
     Int_t nSensors = sensNodes->GetEntriesFast();
     // !!! AZ - Crashes at root exit
     for (Int_t i = 0; i < nSensors; ++i) {
       FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->At(i));
       TArrayD* params = sensVol->getParameters();
       par1[i] = params->At(0); // half length of sensor
       par2[i] = params->At(1); // half tickness of sensor
       par3[i] = params->At(2); // half width of sensor

       //   cout << " *** STS sensitive volume: " << sensVol->GetName() << " " << params->At(0) << " "
       //	   << params->At(1) << " " << params->At(2) << endl;
     }
   }
   */

   //------ Reset output array
   fDigiArray->Clear();

   //------ Declare some variables
   MpdStsPoint *point = nullptr;

   //------ Loop over MCPoints
   Int_t nPoints = fPointArray->GetEntriesFast();
   cout << "********* Number of STS points: " << nPoints << endl;

   for (Int_t iPoint = 0; iPoint < nPoints; ++iPoint) {
      point = (MpdStsPoint *)fPointArray->At(iPoint);

      //    cout << "# of current STS point: " << iPoint+1 << endl;

      detId = point->GetDetectorID();
      // sector = (detId & 127);
      // ladder = ((detId >> 7) & 63);
      // layer = ((detId >> 13) & 7);

      Double_t pos[3] = {0};
      // AZ TVector3 pos3;
      TVector3 pos3, pos3out, mom3; // AZ-161020
      point->Position(pos3);
      // AZ-161020 pos3.GetXYZ(pos); // global coordinate

      // AZ-161020 TGeoNode *node = gGeoManager->FindNode(pos[0],pos[1],pos[2]);
      TGeoNode *node = gGeoManager->FindNode(pos3[0], pos3[1], pos3[2]);
      if (node == 0x0) {
         cout << " !!! Missing node " << pos3[0] << " " << pos3[1] << " " << pos3[2] << endl;
         exit(0);
      }
      // cout << "Node: " << node->GetName() << "  Point position: " << pos[0] << " " << pos[1] << " " << pos[2] <<
      // endl;
      if (first) {
         first      = 0;
         spd_length = ((TGeoBBox *)node->GetVolume()->GetShape())->GetDZ() * 2;
         spd_width  = ((TGeoBBox *)node->GetVolume()->GetShape())->GetDX() * 2;
      }

      // cout << gGeoManager->GetPath() << " DetectorID = " << detId << endl;
      // cout << " Layer: "<< layer << " Ladder: " << ladder << " Sector: " << sector << endl;

      // AZ-161020
      point->PositionOut(pos3out);
      point->Momentum(mom3);
      mom3.SetMag(0.0002); // step 2 um
      mom3 *= -1;
      pos3out += mom3; // step back
      pos3 += pos3out;
      pos3 *= 0.5;
      pos3.GetXYZ(pos); // global coordinate
      // AZ

      Double_t loc[3] = {0.};
      gGeoManager->MasterToLocal(pos, loc);

      Double_t x = loc[0], y = loc[1], z = loc[2]; // local coordinate
      Double_t xlocSensor;
      Int_t    iCol, iRow;

      if (x < -0.5 * spd_width || x > 0.5 * spd_width || z < -0.5 * spd_length || z > 0.5 * spd_length) {

         cout << "***** Point is outside SPD sensor! " << endl;
         cout << "***** Point is in volume: " << gGeoManager->GetPath() << endl;
         cout << "***** Local coordinates of point are: " << x << " " << y << " " << z << " " << endl;
         // AZ x = -spd_width*(0.5+1./nRow);  // iRow=-1
         z = -spd_length * (0.5 + 1. / nCol); // iCol=-1
      }

      iCol = nCol * (z / spd_length + 0.5);
      iRow = nRow * (x / spd_width + 0.5);

      xlocSensor = x;

      //       cout << "Global coordinates of point in sensor: " << pos[0] << " " << pos[1] << " " << pos[2] << " " <<
      //       endl; cout << "Local coordinates of point in sensor: " << loc[0] << " " << loc[1] << " " << loc[2] << " "
      //       << endl; cout << " Col= " << iCol << " Row= " << iRow << " Layer: " << layer << endl;

      MpdItsHit5spd *hit = AddHit(iPoint, detId); // One hit per point
      hit->SetUniqueID(2);
      hit->SetTrackID(iPoint);
      hit->SetLocalX(xlocSensor);
      hit->SetCol(iCol);
      hit->SetRow(iRow);

   } // points
}

//------   Public method Create Structure   ------------------------------------
void MpdItsHitProducer5spd::CreateStructure() {}

//------   Private method AddDigi   --------------------------------------------
MpdItsHit5spd *MpdItsHitProducer5spd::AddHit(Int_t trackID, Int_t detID)
{
   // It fills the MpdStsHit category
   MpdStsPoint * point = (MpdStsPoint *)fPointArray->UncheckedAt(trackID);
   TClonesArray &clref = *fDigiArray;
   Int_t         size  = clref.GetEntriesFast();
   TVector3      pos, posOut;
   point->Position(pos);       // one detector side
   point->PositionOut(posOut); // the other side
   pos += posOut;
   pos *= 0.5;
   Double_t errx, errz;
   gRandom->Rannor(errx, errz);
   return new (clref[size]) MpdItsHit5spd(detID, pos, TVector3(errx, 0., errz), point->GetEnergyLoss(), trackID, 0);
}

//__________________________________________________________________________

void MpdItsHitProducer5spd::Finish()
{
   ///

   delete fDigiArray;
   fDigiArray = nullptr;
}

ClassImp(MpdItsHitProducer5spd);