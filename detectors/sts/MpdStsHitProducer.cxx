/////////////////////////////////////////////////////////////
//
//  MpdStsHitProducer
//
//  Filler of MpdStsHit
//
//
//
///////////////////////////////////////////////////////////////

#include "TClonesArray.h"

#include "MpdStsGeoPar.h"
#include "MpdStsHitProducer.h"
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
MpdStsHitProducer::MpdStsHitProducer(const char *fileGeo) : FairTask("Ideal STS hit Producer")
{
   fFileGeo = fileGeo;
   eneThr   = 0.001; // Energy threshold for STS
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdStsHitProducer::~MpdStsHitProducer() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------
InitStatus MpdStsHitProducer::Init()
{

   cout << " INITIALIZATION *********************" << endl;

   // FairDetector::Initialize();
   // FairRun* sim = FairRun::Instance();
   // FairRuntimeDb* rtdb=sim->GetRuntimeDb();

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdStsHitProducer::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array
   fPointArray = (TClonesArray *)ioman->GetObject("StsPoint");
   if (!fPointArray) {
      cout << "-W- MpdStsHitProducer::Init: "
           << "No StsPoint array!" << endl;
      return kERROR;
   }

   // Create and register output array
   fDigiArray = new TClonesArray("MpdStsHit");

   ioman->Register("StsHit", "Sts", fDigiArray, kTRUE);
   // ioman->Register("StsHit","Sts",fDigiArray,kFALSE);
   CreateStructure();

   cout << "-I- MpdStsHitProducer: Intialization successfull" << endl;

   return kSUCCESS;
}
// -------------------------------------------------------------------------

// -----   Private method SetParContainers   -------------------------------
void MpdStsHitProducer::SetParContainers()
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
void MpdStsHitProducer::Exec(Option_t *opt)
{

   // cout << " DIGI EXECUTION *********************" << endl;
   static Int_t eventCounter = 0, first = 1;
   // Layer Rmin, Rmax, module Z-size, layer Z-half-size
   static Double_t rMin[4] = {0.}, rMax[4] = {0.}, zMod[4] = {0.}, dZ[4] = {0.};
   static Double_t sterCos[2] = {-7.5, 7.5}, sterSin[2] = {0}; // stereo-angles
   const Double_t  size = 6.2, pitch = 0.01;                   // 100 um

   if (first) {
      if (!fDigiArray) Fatal("MpdStsHitProducer::Exec", "No DigiArray");
      // Get detector parameters
      first               = 0;
      FairRuntimeDb *rtdb = FairRun::Instance()->GetRuntimeDb();
      // rtdb->printParamContexts();
      // cout << rtdb->findContainer("TpcGeoPar") << " " << rtdb->findContainer("TpcContFact:") << endl;
      MpdStsGeoPar *geoPar = (MpdStsGeoPar *)rtdb->getContainer("MpdStsGeoPar");
      // cout << geoPar << endl;
      TString    volName   = "sts01 ";
      TObjArray *sensNodes = geoPar->GetGeoSensitiveNodes();
      // cout << sensNodes->GetEntriesFast() << " " << geoPar->GetGeoPassiveNodes()->GetEntriesFast() << endl;
      Int_t nLay = 4;
      for (Int_t i = 0; i < nLay; ++i) {
         volName[5]           = 97 + i; // 'a', 'b', ..
         FairGeoNode *sensVol = (FairGeoNode *)(sensNodes->FindObject(volName));
         // FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->At(0));
         TArrayD *params = sensVol->getParameters();
         rMin[i]         = params->At(0);
         rMax[i]         = params->At(1);
         // dR = (rMax-rMin) / 50;
         dZ[i]       = params->At(2);
         Int_t nMods = Int_t(dZ[i] * 2. / size + 0.1);
         zMod[i]     = dZ[i] * 2. / nMods;
         cout << " *** STS sensitive volume: " << sensVol->GetName() << " " << rMin[i] << " " << dZ[i] << " " << zMod[i]
              << " " << nMods << endl;
      }
      for (Int_t i = 0; i < 2; ++i) {
         sterCos[i] *= TMath::DegToRad();     // angle in rads
         sterSin[i] = TMath::Sin(sterCos[i]); // sin
         sterCos[i] = TMath::Cos(sterCos[i]); // cos
      }
   }

   cout << " StsHitProducer event " << ++eventCounter << endl;

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
      point        = (MpdStsPoint *)fPointArray->At(iPoint);
      Int_t detId  = point->GetDetectorID();
      Int_t lay    = (detId >> 1) - 1;
      Int_t module = Int_t((point->GetZ() + dZ[lay]) / zMod[lay]);
      detId        = detId | (module << 10);
      for (Int_t j = 0; j < 2; ++j) {
         // Two hits per point (for two sides)
         MpdStsHit *hit  = AddHit(iPoint, detId + j, j);
         Double_t   zloc = hit->GetZ() + dZ[lay];
         Double_t   phi  = TMath::ATan2(hit->GetY(), hit->GetX());
         Double_t   rad  = (j == 0) ? rMin[lay] : rMax[lay];
         // Double_t x = rad * phi, xMax = TMath::TwoPi() * rad;
         Double_t x = rad * phi, xMax = TMath::TwoPi() * rad * sterCos[j];
         // Double_t xloc = x * sterCos[j] + zloc * sterSin[j]; // rotated coordinate
         Double_t xloc = x * sterCos[j] + hit->GetZ() * sterSin[j]; // rotated coordinate
         if (xloc < 0) xloc += xMax;
         // xloc0 = xloc0 - TMath::Floor(xloc0/xlocMx) * xlocMx; // modulo
         xloc = fmod(xloc, xMax);
         // cout << xloc0 << " " << fmod(xloc0,xlocMx) << endl;
         // if (TMath::Abs(xloc0-fmod(xloc0,xlocMx)) > 1.e-7) exit(0);
         Int_t strip = Int_t(xloc / pitch);
         hit->SetLocalX(xloc);
         hit->SetStrip(strip);
         // hit->SetDetID(Int_t sectorType, Int_t layer, Int_t ladder, Int_t det, Int_t side);
         hit->SetDetId(0, lay + 1, 0, module, j);
         hit->SetUniqueID(1);
      }
      // cout << point->GetDetectorID() << " " << lay << " " << point->GetZ() << " " << module << " " << detId << endl;
      // fTrackEnergy[point->GetDetectorID()] += point->GetEnergyLoss();
   }

   // Loop to register ZdcHit
   /*
   for(p=fTrackEnergy.begin(); p!=fTrackEnergy.end(); ++p) {
     if ((*p).second>eneThr)
       AddHit(1, (*p).first, (*p).second);
   }
   */
}
// -------------------------------------------------------------------------

// -----   Public method Create Structure   --------------------------------
void MpdStsHitProducer::CreateStructure()
{

   /*
   TString work = getenv("VMCWORKDIR");
   work = work + "/geometry/" + fFileGeo;
   cout << "-I- <MpdStsHitProducer::CreateStructure> Zdc geometry loaded from: "
        << work << endl;

   Int_t detId = -1;
   MpdStsReader read(work);

   for(Int_t module=1; module<=read.GetMaxModules(); module++)
     for(Int_t row=1; row<=read.GetMaxRows(module); row++)
       for(Int_t crystal=1; crystal<=read.GetMaxCrystals(module,row); crystal++) {
    DataG4 data = read.GetData(module,row,crystal);
    for(Int_t copy=1; copy<=20; copy++) {
      detId =  module*100000000 + row*1000000 + copy*10000 + crystal;
      emcX[detId] = data.posX; emcY[detId] = data.posY; emcZ[detId] = data.posZ;
      emcTheta[detId] = data.theta; emcTau[detId] = data.tau;
      if (module==3)
        emcPhi[detId] = fmod(data.phi+90.*(copy-1),360);
      else
        emcPhi[detId] = fmod(data.phi+22.5*(copy-1),360);
    }
       }
   */
}

// -----   Private method AddDigi   --------------------------------------------
MpdStsHit *MpdStsHitProducer::AddHit(Int_t trackID, Int_t detID, Int_t side)
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

ClassImp(MpdStsHitProducer)
