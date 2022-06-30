/*************************************************************************************
 *
 *            MpdZdcDigiProducer
 *    Class to create digital data taken from MpdZdc detector
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  18-Apr-2008
 *  Modified March 2021  by A.Strijak
 *  Modified June 2022  by S.Morozovs
 *
 ************************************************************************************/

#include <iostream>
#include "TClonesArray.h"

#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "MpdZdcDigiProducer.h"
#include "MpdZdcDigi.h"
#include "MpdZdcPoint.h"

#include "TROOT.h"
#include "TVectorT.h"

// -----   Default constructor   -------------------------------------------
MpdZdcDigiProducer::MpdZdcDigiProducer(const char *name) : FairTask(name)
{
   fPointArray = 0;
   fDigiArray  = 0;
   fGeoPar     = 0;
   //  fELossZdc1Value = NULL, fELossZdc2Value = NULLL;
   fPix2Mip     = 15;    // 15 MPPC pixels per MIP
   fMIPEnergy   = 0.005; // 5 MeV
   fMIPNoise    = 0.3;   // 0.2 MIP noise level
   fMIP2GeV     = 0.005;
   fMappingFile = "nofile";
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdZdcDigiProducer::~MpdZdcDigiProducer() {}
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
void MpdZdcDigiProducer::SetParContainers()
{
   cout << "-I- MpdZdcDigiProducer: SetParContainers started..." << endl;

   //   Get run and runtime database
   FairRunAna *run = FairRunAna::Instance();
   if (!run) Fatal("FairMuchDigitize::SetParContainers", "No analysis run");

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   if (!rtdb) Fatal("FairMuchDigitize::SetParContainers", "No runtime database");

   cout << "-I- MpdZdcDigiProducer: SetParContainers continued..." << endl;

   rtdb->activateParIo(rtdb->getFirstInput());
   //   fGeoPar=( MpdZdcGeoPar*) rtdb->getContainer("MpdZdcGeoPar");
   fGeoPar = (MpdZdcGeoPar *)gROOT->FindObject("MpdZdcGeoPar");
   fGeoPar->print();

   cout << "-I- MpdZdcDigiProducer: SetParContainers finished." << endl;
}

// -------------------------------------------------------------------------
// -----   Public method Init   --------------------------------------------
InitStatus MpdZdcDigiProducer::Init()
{

   cout << "-I- MpdZdcDigiProducer: Init started..." << endl;

   fRandom3 = new TRandom3();

   // Get FHCal module mapping (if mapping file is set)

   for (int i = 0; i < 91; i++) {
      fModuleX[i] = -100;
      fModuleY[i] = -100;
   }

   if (!(fMappingFile == "nofile")) {

      TString  dummy;
      ifstream in;

      TString dir  = getenv("VMCWORKDIR");
      TString path = dir + "/input/";
      in.open((path + fMappingFile).Data());
      if (!in.is_open()) {
         printf("Loading FHCal Map from file: %s - file open error!\n", fMappingFile.Data());
         return kFATAL;
      }
      printf("Loading FHCal Map from file: %s\n", fMappingFile.Data());
      in >> dummy >> dummy >> dummy;
      while (!in.eof()) {
         int   id;
         float x, y;
         in >> id >> x >> y;
         if (!in.good()) break;
         printf("%d %f %f\n", id, x, y);
         fModuleX[id] = x;
         fModuleY[id] = y;
      }

      in.close();

   } // if (!(fMappingFile == "nofile"))

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdZdcDigiProducer::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array
   fPointArray = (TClonesArray *)ioman->GetObject("ZdcPoint");
   if (!fPointArray) {
      cout << "-W- MpdZdcDigiProducer::Init: "
           << "No ZdcPoint array!" << endl;
      return kERROR;
   }

   // Create and register output array
   fDigiArray = new TClonesArray("MpdZdcDigi");
   ioman->Register("ZdcDigi", "Zdc", fDigiArray, kTRUE);

   /*
     fELossZdc1Value = new TClonesArray("TParameter<double>");
     ioman->Register("ELossZdc1Value","Zdc",fELossZdc1Value,kTRUE);

     fELossZdc2Value = new TClonesArray("TParameter<double>");
     ioman->Register("ELossZdc2Value","Zdc",fELossZdc2Value,kTRUE);
   */
   MpdZdcDigiScheme *fDigiScheme = MpdZdcDigiScheme::Instance();
   fDigiScheme->Init(fGeoPar, 0, kTRUE, 2);

   cout << "-I- MpdZdcDigiProducer: Intialization successfull" << endl;

   return kSUCCESS;
}
// -----   Public method Exec   --------------------------------------------
void MpdZdcDigiProducer::Exec(Option_t *opt)
{

   //#define EDEBUG
#ifdef EDEBUG
   static Int_t lEDEBUGcounter = 0;
   cout << "EDEBUG-- MpdZdcDigiProducer::Exec() started... " << endl;
   ;
#endif

   if (!fDigiArray) Fatal("Exec", "No DigiArray");

   fDigiArray->Clear();

   MpdZdcDigiScheme *pDigiScheme = MpdZdcDigiScheme::Instance();

   if (!pDigiScheme) Fatal("MpdZdcDigiProducer::Exec", "No DigiScheme");

   Int_t          detID, modID, chanID;
   MpdZdcDigiId_t digiID;

   // marina
   Double_t dEdepSectEv[90][10];

   for (Int_t i = 0; i < 90; i++) {       // mod
      for (Int_t ii = 0; ii < 10; ii++) { // section
         dEdepSectEv[i][ii] = 0.;
      }
   }
   // end marina

   MpdZdcPoint *point = NULL;

   map<MpdZdcDigiId_t, Float_t> fDigiIdEnergy;
   fDigiIdEnergy.clear();
   map<MpdZdcDigiId_t, Float_t>::const_iterator p;

   Int_t    nPoints = fPointArray->GetEntriesFast();
   Double_t e1 = 0, e2 = 0;

   // cout <<"marina " <<nPoints <<endl;
   for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {
      // marina
      point = (MpdZdcPoint *)fPointArray->At(iPoint);
      // cout <<"marina 1 " <<point->GetCopyZdc() <<" " <<point->GetCopyMother() <<" " <<point->GetCopy() <<endl;

      detID = point->GetCopyZdc(); //==1 (z>0), ==2 (z<0)
      // modID  = point->GetCopyMother(); // modules 1-45
      if (detID == 1)
         modID = point->GetCopyMother(); // modules 1-45
      else
         modID = point->GetCopyMother() + 45;       // modules 46-90
      chanID = (Int_t)((point->GetCopy() - 1) / 6); // sections 0-9

      //    cout <<"marina 1 " <<detID <<" " <<modID <<" " <<point->GetCopy() <<' ' <<chanID <<endl;

      dEdepSectEv[modID - 1][chanID] += point->GetEnergyLoss();
      // end marina

      /*
          if (detID == 1) {
            e1 += point->GetEnergyLoss();
          }
          else e2 += point->GetEnergyLoss();
      */
      Int_t pMMcopy = 1 * (point->GetZ() > 0) + 2 * (point->GetZ() < 0);
      digiID = pDigiScheme->GetDigiIdFromVolumeData(point->GetDetectorID(), point->GetCopy(), point->GetCopyMother(),
                                                    pMMcopy);

      if ((digiID[0] != -1) && (digiID[1] != -1)) {

         if (fDigiIdEnergy.find(digiID) == fDigiIdEnergy.end())
            fDigiIdEnergy[digiID] = point->GetEnergyLoss();
         else
            fDigiIdEnergy[digiID] += point->GetEnergyLoss();

         if (pMMcopy == 1) {
            e1 += point->GetEnergyLoss();
         } else {
            e2 += point->GetEnergyLoss();
         }

      } // if ((digiID[0]!=-1)&&(digiID[1]!=-1))

#ifdef EDEBUG
      else {
         if (lEDEBUGcounter < 100) {
            cout << "EDEBUG-- MpdZdcDigiProducer::Exec:  Boundary point? : ";
            point->Print("");
            lEDEBUGcounter++;
         }
      }
#endif

   } // for (Int_t iPoint=0; iPoint<nPoints; iPoint++)

   // cout <<"marina 2 " <<endl;
   /*
     TClonesArray& clref1 = *fELossZdc1Value;
     new(clref1[0]) TParameter<double>("ELossZdc1",e1);
     TClonesArray& clref2 = *fELossZdc2Value;
     new(clref2[0]) TParameter<double>("ELossZdc2",e2);
   */
   e1 = 0;
   e2 = 0;

   // cout <<"marina 3 " <<endl;
   for (Int_t i = 0; i < 90; i++) {       // mod
      for (Int_t ii = 0; ii < 10; ii++) { // section
         // cout <<"dEdepSectEv " <<i <<" " <<ii <<" " <<dEdepSectEv[i][ii] <<endl;
         if (dEdepSectEv[i][ii] > 0) {

            if (i <= 44)
               detID = 1;
            else
               detID = 2;

            Int_t hitID = i * 10 + ii;

            Double_t recEnergy = RecoEnergy(dEdepSectEv[i][ii]);

            if (detID == 1) {
               MpdZdcDigi *digi = AddHit(detID, i + 1, ii + 1, dEdepSectEv[i][ii]);
               digi->ConvertSim();
               digi->SetELossReco(recEnergy);
               e1 += recEnergy;
               fHitMap[hitID] = digi;
               digi->SetModuleX(fModuleX[i + 1]);
               digi->SetModuleY(fModuleY[i + 1]);
            } else {
               MpdZdcDigi *digi = AddHit(detID, i - 45 + 1, ii + 1, dEdepSectEv[i][ii]);
               digi->ConvertSim();
               digi->SetELossReco(recEnergy);
               e2 += recEnergy;
               fHitMap[hitID] = digi;
               digi->SetModuleX(fModuleX[i + 1]);
               digi->SetModuleY(fModuleY[i + 1]);
            }
         }
      }
   } // for(Int_t i=0; i<90; i++)

   //    TClonesArray& clref1 = *fELossZdc1Value;
   //    new(clref1[0]) TParameter<double>("ELossZdc1",e1);
   //    TClonesArray& clref2 = *fELossZdc2Value;
   //    new(clref2[0]) TParameter<double>("ELossZdc2",e2);

   /*
 for(p=fDigiIdEnergy.begin(); p!=fDigiIdEnergy.end(); ++p) {

   pDigiScheme->SplitDigiID((*p).first, detID, modID, chanID);


   if ((detID!=-1)&&(chanID!=-1)) {
     MpdZdcDigi* digi = AddHit(detID, modID, chanID, (*p).second);
     digi->ConvertSim();
#ifdef EDEBUG
     if (lEDEBUGcounter<20) {
  cout << "EDEBUG-- MpdZdcDigiProducer::Exec: "<< detID<< " " << chanID << "   " <<
    (*p).second << "     " << lEDEBUGcounter << endl;
  lEDEBUGcounter++;
     }
#endif
   }

 }
   */

   map<Float_t, Float_t> contribTMP;

   for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {

      point = (MpdZdcPoint *)fPointArray->At(iPoint);

      detID = point->GetCopyZdc(); //==1 (z>0), ==2 (z<0)
      if (detID == 1)
         modID = point->GetCopyMother(); // modules 1-45
      else
         modID = point->GetCopyMother() + 45;       // modules 46-90
      chanID = (Int_t)((point->GetCopy() - 1) / 6); // sections 0-9

      if (dEdepSectEv[modID - 1][chanID] == 0) continue;

      contribTMP.clear();

      Int_t hitID = (modID - 1) * 10 + chanID;

      MpdZdcDigi *hit = SearchHitVR(hitID);

      if (hit == NULL) {
         cout << "Lost ZDC digit??? : " << detID << " " << modID << " " << chanID << " " << hitID << " "
              << dEdepSectEv[modID][chanID] << " " << endl;
      } else {

         contribTMP = point->GetContribTimes();

         for (map<Float_t, Float_t>::iterator it = contribTMP.begin(); it != contribTMP.end(); ++it) {
            hit->IncreaseTimesCorrVR(it->second, (it->first) + (405.6 - fabs(point->GetZ())) / (30 / 1.6));
         }
      }
   } // iPoint

#undef EDEBUG
}
// -------------------------------------------------------------------------

// -----   Private method AddDigi   --------------------------------------------
MpdZdcDigi *MpdZdcDigiProducer::AddHit(Int_t detID, Int_t modID, Int_t chanID, Float_t energy)
{
   TClonesArray &clref = *fDigiArray;
   Int_t         size  = clref.GetEntriesFast();
   // cout <<"size " <<size <<endl;
   MpdZdcDigi *result = new (clref[size]) MpdZdcDigi(detID, modID, chanID, energy);
   // cout <<"result " <<result <<endl;
   return result;
}
// ----

Double_t MpdZdcDigiProducer::RecoEnergy(Double_t pfELoss)
{
   Double_t energyMIP        = pfELoss / fMIPEnergy;
   Double_t energyPix        = fRandom3->Poisson(energyMIP * fPix2Mip);
   Double_t energyMIPSmeared = energyPix / fPix2Mip;
   Double_t noise            = fRandom3->Gaus(0, fMIPNoise);
   energyMIPSmeared += noise;
   return energyMIPSmeared * fMIP2GeV;
}

MpdZdcDigi *MpdZdcDigiProducer::SearchHitVR(Int_t sec)
{
   if (fHitMap.find(sec) == fHitMap.end()) return NULL;
   return fHitMap[sec];
}

ClassImp(MpdZdcDigiProducer);
