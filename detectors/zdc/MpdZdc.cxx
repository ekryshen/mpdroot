/*************************************************************************************
 *
 *         Class MpdZdc
 *
 *  Adopted for MPD by:   Elena Litvinenko (EL)
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *  Last update:  22-Feb-2012 (EL)
 *
 ************************************************************************************/

#include <iostream>

#include "TClonesArray.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include "TGeoArb8.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "MpdZdcGeo.h"
#include "FairGeoRootBuilder.h"
#include "MpdStack.h"
#include "MpdZdc.h"
#include "MpdZdcPoint.h"

#include "FairRootManager.h"
#include "FairVolume.h"
// add on for debug
//#include "FairGeoG3Builder.h"
#include "FairRuntimeDb.h"
#include "TObjArray.h"
#include "FairRun.h"

#include "TParticlePDG.h"

// -----   Default constructor   -------------------------------------------
MpdZdc::MpdZdc()
{
   fZdcCollection = new TClonesArray("MpdZdcPoint");
   volDetector    = 0;
   fPosIndex      = 0;
   // fpreflag = 0;
   // fpostflag = 0;
   fEventID      = -1;
   fHitNb        = 0;
   fVSCVolId     = 0;
   fVSCHVolId    = 0;
   fVerboseLevel = 1;
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
MpdZdc::MpdZdc(const char *name, Bool_t active) : FairDetector(name, active)
{
   fZdcCollection = new TClonesArray("MpdZdcPoint");
   fPosIndex      = 0;
   volDetector    = 0;
   // fpreflag = 0;
   // fpostflag = 0;
   fEventID      = -1;
   fHitNb        = 0;
   fVSCVolId     = 0;
   fVSCHVolId    = 0;
   fVerboseLevel = 1;
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdZdc::~MpdZdc()
{
   if (fZdcCollection) {
      fZdcCollection->Delete();
      delete fZdcCollection;
   }
}
// -------------------------------------------------------------------------

// -----   Public method Intialize   ---------------------------------------
void MpdZdc::Initialize()
{
   // Init function

   FairDetector::Initialize();
   // FairRun* sim = FairRun::Instance();
   // FairRuntimeDb* rtdb=sim->GetRuntimeDb();

   fVSCVolId  = gMC->VolId("zdc01s");
   fVSCHVolId = gMC->VolId("ScH");
}
// -------------------------------------------------------------------------
void MpdZdc::BeginEvent()
{
   // Begin of the event
}

//_____________________________________________________________________________
MpdZdcPoint *MpdZdc::GetHit(Int_t i) const
{
   // Returns the hit for the specified layer.
   // ---

   // return (Ex03CalorHit*)fCalCollection->At(i);
   return (MpdZdcPoint *)fZdcCollection->At(i);
}

//_____________________________________________________________________________
MpdZdcPoint *MpdZdc::GetHit(Int_t vsc, Int_t mod, Int_t zdc) const
{
   // Returns the hit for the specified vsc and module0.
   // ---

   MpdZdcPoint *hit;
   Int_t        nofHits = fZdcCollection->GetEntriesFast();
   // cout <<"NAVetoHit: nofHits " <<nofHits << " vsc " << vsc << " mod " << mod <<endl;

   // for (Int_t i=0; i<nofHits; i++) {
   //  hit =  GetHit(i);
   //  cout <<"GetHit method-> " <<i <<" " <<hit->GetCopy() <<" " <<hit->GetCopyMother()  <<" " <<hit->GetCopyZdc()<<" "
   //  <<hit->GetEnergyLoss() <<" " <<hit->GetZ() <<endl;
   //}

   for (Int_t i = 0; i < nofHits; i++) {
      hit = GetHit(i);
      // cout <<"fEdep " <<i <<" " <<hit->GetEdep() <<endl;
      // int iVSCId = hit->GetVSCId();
      // int iMODId = hit->GetMODId();
      // if(iVSCId == vsc && iMODId == mod)
      if (hit->GetCopy() == vsc && hit->GetCopyMother() == mod && hit->GetCopyZdc() == zdc) return hit;
   }

   return 0;
}
/*
//_____________________________________________________________________________
MpdZdcPoint* MpdZdc::GetHitH(Int_t vsch, Int_t modh) const
{
// Returns the hit for the specified vsc and module0.
// ---

  MpdZdcPoint *hit;
  Int_t nofHits = fZdcCollection->GetEntriesFast();
  //cout <<"NAVetoHit: nofHits " <<nofHits << " vsc " << vsc << " mod " << mod <<endl;

  //for (Int_t i=0; i<nofHits; i++) {
  //hit =  GetHit(i);
  //if(hit->GetCopyMotherH()<0 || hit->GetCopyH()<0) cout <<"GetHitH -59 -> " <<i <<" " <<hit->GetCopyMotherH() <<" "
<<hit->GetEnergyLoss() <<" " <<hit->GetCopyH() <<endl;
  //}

  for (Int_t i=0; i<nofHits; i++) {
    hit =  GetHit(i);
    //cout <<"fEdep " <<i <<" " <<hit->GetEdep() <<endl;
    //int iVSCId = hit->GetVSCId();
    //int iMODId = hit->GetMODId();
    //if(iVSCId == vsc && iMODId == mod)
    //cout <<"GetHitH " <<i <<" " <<hit->GetCopyMotherH() <<" " <<hit->GetEnergyLoss() <<" " <<hit->GetCopyH() <<endl;
    //if(hit->GetCopyH() == vsch && hit->GetCopyMotherH() == modh)
    if(hit->GetCopy() == vsch && hit->GetCopyMother() == modh)
      return hit;
  }

  return 0;
}
*/

// -----   Public method ProcessHits  --------------------------------------
Bool_t MpdZdc::ProcessHits(FairVolume *vol)
{
   // if (TMath::Abs(gMC->TrackCharge()) <= 0) return kFALSE;

   // copyNoVTYVEC,copyNoVMOD,copyNoVZDC;
   //,copyNoVTYVECH,copyNoVMODH,copyNoVZDCH; copyNoVSCCom
   Int_t copyNoVSC{0}, copyNoVSCH{0}, copyNoVTYVECCom{0}, copyNoVMODCom{0}, copyNoVZDCCom{0};

   Int_t          ivol;
   TLorentzVector tPos1, tMom1;
   TLorentzVector tPos, tMom;
   // Int_t copyNo;
   // Int_t ivol1;
   Int_t    module0; //, module_h; //module0, module0 with hole
   Int_t    slice;   //, slice_h;
   Int_t    zdc;     //,zdc_h;
   Double_t time   = 0;
   Double_t length = 0;
   // TParticle* part;
   // Double_t charge;
   Double_t QCF       = 1;    // quenching for Birk
   Double_t BirkConst = 12.6; // 0.126 mm/MeV for polystyrene
   // 0.126 *(0.1/0.001) = 12.6 cm/GeV
   //(0.126 mm/MeV - from Wikipedia, 0.07943mm/MeV ô in Geant4,G4EmSaturation.cc (push))

   //#define EDEBUG
#ifdef EDEBUG
   static Int_t lEDEBUGcounter = 0;
   if (lEDEBUGcounter < 1) std::cout << "EDEBUG-- MpdZdc::ProcessHits: entered" << gMC->CurrentVolPath() << endl;
#endif

   if (gMC->CurrentVolID(copyNoVSC) != fVSCVolId && gMC->CurrentVolID(copyNoVSCH) != fVSCHVolId) {
      return kFALSE;
   }

   /*
   if (gMC->CurrentVolID(copyNoVSC) == fVSCVolId) {
     gMC->CurrentVolOffID(1, slice);
     gMC->CurrentVolOffID(2, module0);
     gMC->CurrentVolOffID(3, zdc);
     //copyNo = (copyNoVMOD-1)*60 + copyNoVTYVEC; //1-6480
     copyNoVTYVECCom = slice; copyNoVMODCom = module0; copyNoVZDCCom = zdc;
   }
   if (gMC->CurrentVolID(copyNoVSCH) == fVSCHVolId) {
     gMC->CurrentVolOffID(1, slice);
     gMC->CurrentVolOffID(2, module0);
     gMC->CurrentVolOffID(3, zdc);
     copyNoVTYVECCom = slice; copyNoVMODCom = module0; copyNoVZDCCom = zdc;
   }
   */

   ivol = vol->getMCid();

   if (gMC->CurrentVolID(copyNoVSC) == fVSCVolId || gMC->CurrentVolID(copyNoVSCH) == fVSCHVolId) {
      gMC->CurrentVolOffID(1, slice);
      gMC->CurrentVolOffID(2, module0);
      gMC->CurrentVolOffID(3, zdc);
      // copyNo = (copyNoVMOD-1)*60 + copyNoVTYVEC; //1-6480
      copyNoVTYVECCom = slice;
      copyNoVMODCom   = module0;
      copyNoVZDCCom   = zdc;
   }

   // cout <<"BEGIN copyNoVZDCCom " <<copyNoVZDCCom <<endl;

   if (gMC->IsTrackEntering()) {

      ResetParameters();
      fELoss = 0.;
      time   = 0.;
      length = 0.;
      gMC->TrackPosition(tPos);
      gMC->TrackMomentum(tMom);

#ifdef EDEBUG
      gMC->TrackPosition(tPos1);
      gMC->TrackMomentum(tMom1);
#endif
      /*
      //Int_t copyNo;
      ivol1 = gMC->CurrentVolID(copyNo);
      gMC->CurrentVolOffID(2, module0);
      gMC->CurrentVolOffID(1, slice);
      */
   } // if (gMC->IsTrackEntering())

   if (gMC->IsTrackInside()) {
      // if (gMC->Edep() <= 0. ) return kFALSE;
      // if (gMC->Edep() > 0. ) fELoss += gMC->Edep();
      ////if (gMC->Edep() <= 1.0e-20 ) return kFALSE;
      ////if (gMC->Edep() > 1.0e-20 ) fELoss += gMC->Edep();

      gMC->TrackPosition(tPos);
      gMC->TrackMomentum(tMom);
      length += gMC->TrackStep();

      // fELoss +=gMC->Edep();
      // Birk corrections
      QCF = 1. + (BirkConst / gMC->TrackStep()) * gMC->Edep();
      fELoss += (gMC->Edep()) / QCF;

      time += gMC->TrackTime() * 1.0e09;

      if (gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {

         // part = gMC->GetStack()->GetCurrentTrack();
         // charge = part->GetPDG()->Charge() / 3. ;

         // Create MpdZdcPoint
         fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
         // time    = gMC->TrackTime() * 1.0e09;
         // length  = gMC->TrackLength();
         // gMC->TrackPosition(tPos);
         // gMC->TrackMomentum(tMom);
         // gMC->CurrentVolOffID(2, module0);
         // gMC->CurrentVolOffID(1, slice);

         if (fELoss > 0) {

            // std::cout << "INSIDE MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  "
            // <<fELoss <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  <<tMom.Pz() <<" " << ivol <<" "
            // <<gMC->CurrentVolOffName(2) << " " <<gMC->CurrentVolOffName(1) << " " << gMC->CurrentVolOffName(0)
            // <<std::endl;

            // cout <<"CHECK INSIDE ivol,copyNo,iCell " <<part->GetPdgCode() <<" " << fTrackID <<" " <<ivol <<" "
            // <<copyNoVZDCCom <<" " <<copyNoVMODCom <<" " <<copyNoVTYVECCom <<" " <<tPos.Z() <<" " <<tMom.Pz() <<" "
            // <<fELoss <<endl;

            /*
            if(copyNoVTYVECCom==slice && copyNoVMODCom==module0 && copyNoVZDCCom==zdc) {//module0
              if ( !GetHit(slice,module0,zdc) ) {
                //cout <<"INSIDE NO GETHIT 1-> " <<slice <<" " <<module0 <<" " <<zdc <<endl;
                //new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, -59, -59,
            TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
                new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(),
            tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
                //cout <<"INSIDE NO GETHIT 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
              }
              else {
                //cout <<"INSIDE GETHIT 1-> " <<slice <<" " <<module0  <<" " <<zdc <<endl;
                GetHit(slice,module0,zdc)->AddVSC(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(),
            tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
                //cout <<"INSIDE GETHIT 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
              }
            }//if(copyNoVTYVECCom==copyNoVTYVEC
            */

            if (copyNoVTYVECCom == slice && copyNoVMODCom == module0 && copyNoVZDCCom == zdc) { // module0
               if (!GetHit(slice, module0, zdc)) {
                  // cout <<"INSIDE NO GETHIT 1-> " <<slice <<" " <<module0 <<" " <<zdc <<endl;
                  // new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, -59, -59,
                  // TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length,
                  // fELoss); new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, zdc,
                  // TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length,
                  // fELoss);
                  AddHit(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),
                         TVector3(tMom.Px(), tMom.Py(), tMom.Pz()), time, length, fELoss);
                  // cout <<"INSIDE NO GETHIT 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
               } else {
                  // cout <<"INSIDE GETHIT 1-> " <<slice <<" " <<module0  <<" " <<zdc <<endl;
                  GetHit(slice, module0, zdc)
                     ->AddVSC(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),
                              TVector3(tMom.Px(), tMom.Py(), tMom.Pz()), time, length, fELoss);
                  // cout <<"INSIDE GETHIT 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
               }
            } // if(copyNoVTYVECCom==copyNoVTYVEC

            /*
            if(copyNoVTYVECCom==slice_h && copyNoVMODCom==module_h && copyNoVZDCCom==zdc_h) {//module0 with hole
              if ( !GetHitH(slice_h,module_h) ) {
                cout <<"INSIDE NO GetHitH 1-> " <<slice_h <<" " <<module_h <<" " <<zdc_h  <<endl;
                //new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, -59, -59, slice_h, module_h,
            TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
                new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice_h, module_h, TVector3(tPos.X(),
            tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"INSIDE NO
            GetHitH 2-> " <<slice_h <<" "<<module_h <<" " <<zdc_h <<endl;
              }
              else {
                cout <<"INSIDE GetHitH 1-> " <<slice_h <<" "<<module_h <<" " <<zdc_h <<endl;
                GetHitH(slice_h,module_h)->AddVSCH(fTrackID, ivol, slice_h, module_h, TVector3(tPos.X(), tPos.Y(),
            tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"INSIDE GetHitH 2-> "
            <<slice_h <<" "<<module_h <<" " <<zdc_h <<endl;
              }
            }//if(copyNoVTYVECCom==copyNoVTYVEC
            */

            /*
              std::cout << "INSIDE MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  "
              <<fELoss <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  << ivol << "=="<<
              gMC->CurrentVolID(copyNo) << ","<< copyNo <<"   "  <<gMC->CurrentVolOffName(2) << " "
              <<gMC->CurrentVolOffName(1) << " " << gMC->CurrentVolOffName(0) <<std::endl;

              cout <<"CHECK INSIDE ivol,copyNo,iCell " <<part->GetPdgCode() <<" " << fTrackID <<" " <<ivol <<" "
              <<module0 <<" " <<slice <<" " <<tPos.Z() <<" " <<tMom.Pz() <<" " <<fELoss <<endl; AddHit(fTrackID, ivol,
              slice, module0, TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time,
              length, fELoss);
            */
            /*
              std::cout << "INSIDE MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  "
              <<fELoss <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  << ivol << "=="<<
              gMC->CurrentVolID(copyNo) << ","<< copyNo <<"   "  <<gMC->CurrentVolOffName(2) << " "
              <<gMC->CurrentVolOffName(1) << " " << gMC->CurrentVolOffName(0) <<std::endl;

              cout <<"CHECK INSIDE ivol,copyNo,iCell " <<part->GetPdgCode() <<" " << fTrackID <<" " <<ivol <<" "
              <<module0 <<" " <<slice <<" " <<tPos.Z() <<" " <<tMom.Pz() <<" " <<fELoss <<endl; AddHit(fTrackID, ivol,
              slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(),
              tMom.Pz()),time, length, fELoss);
            */

         } // if(fELoss>0)
      }    // if ( gMC->IsTrackStop() || gMC->IsTrackDisappeared() )
   }       // if ( gMC->IsTrackInside())

   if (gMC->IsTrackExiting()) {
      // if (gMC->Edep() <= 0. ) return kFALSE;
      // if (gMC->Edep() > 0. ) fELoss += gMC->Edep();
      ////if (gMC->Edep() <= 1.0e-20 ) return kFALSE;
      ////if (gMC->Edep() > 1.0e-20 ) fELoss += gMC->Edep();

      // part    = gMC->GetStack()->GetCurrentTrack();

      // Create MpdZdcPoint
      fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
      time += gMC->TrackTime() * 1.0e09;
      length += gMC->TrackLength();

      // fELoss +=gMC->Edep();
      // Birk corrections
      QCF = 1. + (BirkConst / gMC->TrackStep()) * gMC->Edep();
      fELoss += (gMC->Edep()) / QCF;

      gMC->TrackPosition(tPos);
      gMC->TrackMomentum(tMom);
      // gMC->CurrentVolOffID(2, module0);
      // gMC->CurrentVolOffID(1, slice);

      if (fELoss > 0) {

         // std::cout << "EXIT MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  " <<fELoss
         // <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  <<tMom.Pz() <<" " << ivol <<" "
         // <<gMC->CurrentVolOffName(2) << " " <<gMC->CurrentVolOffName(1) << " " << gMC->CurrentVolOffName(0)
         // <<std::endl;

         // cout <<"CHECK EXIT ivol,copyNo,iCell " <<part->GetPdgCode() <<" " << fTrackID <<" " <<ivol <<" "
         // <<copyNoVZDCCom <<" " <<copyNoVMODCom <<" " <<copyNoVTYVECCom <<" " <<slice <<" " <<module0 <<" " <<zdc <<"
         // "
         // <<tPos.Z() <<" " <<tMom.Pz() <<" " <<fELoss <<endl;

         /*
         if(copyNoVTYVECCom==slice && copyNoVMODCom==module0 && copyNoVZDCCom==zdc) {
           cout <<"EXIT slice module0 " <<slice <<" " <<module0 <<" "  <<zdc <<endl;
           if ( !GetHit(slice,module0,zdc) ) {
               cout <<"EXIT NO GETHIT 1-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
               //new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, -59, -59,
         TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
               new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(),
         tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"EXIT NO GETHIT 2->
         " <<slice <<" "<<module0 <<" " <<zdc <<endl;
           }
           else {
             cout <<"EXIT GetHit 1-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
             GetHit(slice,module0,zdc)->AddVSC(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(),
         tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"EXIT GetHit 2-> " <<slice
         <<" "<<module0 <<" " <<zdc <<endl;
           }
         }//if(copyNoVTYVECCom==copyNoVTYVEC
         */

         if (copyNoVTYVECCom == slice && copyNoVMODCom == module0 && copyNoVZDCCom == zdc) {
            // cout <<"EXIT slice module0 " <<slice <<" " <<module0 <<" "  <<zdc <<endl;
            if (!GetHit(slice, module0, zdc)) {
               // cout <<"EXIT NO GETHIT 1-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
               // new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, -59, -59,
               // TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length,
               // fELoss); new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice, module0, zdc,
               // TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length,
               // fELoss);
               AddHit(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),
                      TVector3(tMom.Px(), tMom.Py(), tMom.Pz()), time, length, fELoss);
               // cout <<"EXIT NO GETHIT 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
            } else {
               // cout <<"EXIT GetHit 1-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
               GetHit(slice, module0, zdc)
                  ->AddVSC(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),
                           TVector3(tMom.Px(), tMom.Py(), tMom.Pz()), time, length, fELoss);
               // cout <<"EXIT GetHit 2-> " <<slice <<" "<<module0 <<" " <<zdc <<endl;
            }
         } // if(copyNoVTYVECCom==copyNoVTYVEC

         /*
         if(copyNoVTYVECCom==slice_h && copyNoVMODCom==module_h && copyNoVZDCCom==zdc_h) {
           cout <<"EXIT slice_h module_h " <<slice <<" " <<module0 <<" "  <<slice_h <<" " <<module_h  <<" " <<zdc <<" "
         <<zdc_h <<endl; if ( !GetHitH(slice_h,module_h) ) { cout <<"EXIT NO GetHitH cxx 1 -> " <<slice_h <<" "
         <<module_h <<" " <<zdc_h <<endl;
             //new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, -59, -59, slice_h, module_h,
         TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss);
             new((*fZdcCollection)[fHitNb++]) MpdZdcPoint(fTrackID, ivol, slice_h, module_h, TVector3(tPos.X(),
         tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"EXIT NO GetHitH
         cxx 2 -> " <<slice_h <<" " <<module_h <<" " <<zdc_h <<endl;
           }
           else {
             cout <<"EXIT GetHitH cxx 1 -> " <<slice_h <<" " <<module_h  <<" " <<zdc_h <<endl;
             GetHitH(slice_h,module_h)->AddVSCH(fTrackID, ivol, slice_h, module_h, TVector3(tPos.X(), tPos.Y(),
         tPos.Z()),TVector3(tMom.Px(), tMom.Py(), tMom.Pz()),time, length, fELoss); cout <<"EXIT GetHitH cxx 2 -> "
         <<slice_h <<" " <<module_h <<" " <<zdc_h <<endl;
           }
         }//if(copyNoVTYVECCom==copyNoVTYVEC
         */

         /*
            std::cout << "EXIT MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  " <<fELoss
         <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  << ivol << "=="<< gMC->CurrentVolID(copyNo) << ","<<
         copyNo <<"   "  <<gMC->CurrentVolOffName(2) << " " <<gMC->CurrentVolOffName(1) << " " <<
         gMC->CurrentVolOffName(0) <<std::endl; cout <<"CHECK EXIT ivol,copyNo,iCell " <<part->GetPdgCode() <<" " <<
         fTrackID <<" " <<ivol <<" " <<module0 <<" " <<slice <<" " <<tPos.X() <<" " <<tPos.Y() <<" " <<tPos.Z() <<" "
         <<tMom.Pz() <<" " <<fELoss <<endl;

         AddHit(fTrackID, ivol, slice, module0, TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(), tMom.Py(),
         tMom.Pz()),time, length, fELoss);
         */
         /*
            std::cout << "EXIT MpdZdc::ProcessHits: TrackID:" <<part->GetPdgCode() <<" " << fTrackID << "  " <<fELoss
            <<" " << gMC->CurrentVolPath() << " " << tPos.Z() <<"   "  << ivol << "=="<< gMC->CurrentVolID(copyNo) <<
            ","<< copyNo <<"   "  <<gMC->CurrentVolOffName(2) << " " <<gMC->CurrentVolOffName(1) << " " <<
            gMC->CurrentVolOffName(0) <<std::endl; cout <<"CHECK EXIT ivol,copyNo,iCell " <<part->GetPdgCode() <<" " <<
            fTrackID <<" " <<ivol <<" " <<module0 <<" " <<slice <<" " <<tPos.X() <<" " <<tPos.Y() <<" " <<tPos.Z() <<" "
            <<tMom.Pz() <<" " <<fELoss <<endl;

            AddHit(fTrackID, ivol, slice, module0, zdc, TVector3(tPos.X(), tPos.Y(), tPos.Z()),TVector3(tMom.Px(),
            tMom.Py(), tMom.Pz()),time, length, fELoss);
         */

      } // if(fELoss>0)
   }    // if ( gMC->IsTrackExiting()) {

   Int_t points = gMC->GetStack()->GetCurrentTrack()->GetMother(1);
   //       Int_t nZdcPoints = (points & (1<<30)) >> 30;
   //       nZdcPoints ++;
   //       if (nZdcPoints > 1) nZdcPoints = 1;
   //      points = ( points & ( ~ (1<<30) ) ) | (nZdcPoints << 30);

   points = (points & (~(1 << 30))) | (1 << 30);

   gMC->GetStack()->GetCurrentTrack()->SetMother(1, points);

   ((MpdStack *)gMC->GetStack())->AddPoint(kZDC);

   //}
   //}
   //     Int_t copyNo;
   //     gMC->CurrentVolID(copyNo);
   //     TString nam = gMC->GetMC()->GetName();
   //    cout<<"name "<<gMC->GetMC()->GetName()<<endl;
   //    ResetParameters();

   return kTRUE;
}

// ----------------------------------------------------------------------------

// -----   Public method EndOfEvent   -----------------------------------------
void MpdZdc::EndOfEvent()
{
   if (fVerboseLevel) Print();
   Reset();
}

// -----   Public method Register   -------------------------------------------
void MpdZdc::Register()
{
   FairRootManager::Instance()->Register("ZdcPoint", "Zdc", fZdcCollection, kTRUE);
}

// -----   Public method GetCollection   --------------------------------------
TClonesArray *MpdZdc::GetCollection(Int_t iColl) const
{
   if (iColl == 0) return fZdcCollection;

   return NULL;
}

// -----   Public method Print   ----------------------------------------------
void MpdZdc::Print() const
{
   Int_t nHits = fZdcCollection->GetEntriesFast();
   cout << "-I- MpdZdc: " << nHits << " points registered in this event." << endl;

   if (fVerboseLevel > 1)
      for (Int_t i = 0; i < nHits; i++) (*fZdcCollection)[i]->Print();
}

// -----   Public method Reset   ----------------------------------------------
void MpdZdc::Reset()
{
   fZdcCollection->Delete();

   fPosIndex = 0;
}

// guarda in FairRootManager::CopyClones
// -----   Public method CopyClones   -----------------------------------------
void MpdZdc::CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset)
{
   Int_t nEntries = cl1->GetEntriesFast();
   // cout << "-I- MpdZdc: " << nEntries << " entries to add." << endl;
   TClonesArray &clref    = *cl2;
   MpdZdcPoint  *oldpoint = NULL;
   for (Int_t i = 0; i < nEntries; i++) {
      oldpoint    = (MpdZdcPoint *)cl1->At(i);
      Int_t index = oldpoint->GetTrackID() + offset;
      oldpoint->SetTrackID(index);
      new (clref[fPosIndex]) MpdZdcPoint(*oldpoint);
      fPosIndex++;
   }
   cout << " -I- MpdZdc: " << cl2->GetEntriesFast() << " merged entries." << endl;
}

// -----   Public method ConstructGeometry   ----------------------------------
void MpdZdc::ConstructGeometry()
{

   TString fileName = GetGeometryFileName();
   if (fileName.EndsWith(".root")) {
      LOG(info) << "Constructing ZDC geometry from ROOT file " << fileName.Data();
      ConstructRootGeometry();
   }
   /*
   else if ( fileName.EndsWith(".geo") )
   {
      FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "Constructing ZDC geometry from ASCII file %s", fileName.Data());
      ConstructAsciiGeometry();
   }
   else	FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "Geometry format of ZDC file %S not supported.",
   fileName.Data());
   */

   FairGeoLoader    *geoLoad = FairGeoLoader::Instance();
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   MpdZdcGeo        *zdcGeo  = new MpdZdcGeo();
   zdcGeo->setGeomFile(GetGeometryFileName());
   geoFace->addGeoModule(zdcGeo);

   Bool_t rc = geoFace->readSet(zdcGeo);
   if (rc) zdcGeo->create(geoLoad->getGeoBuilder());
   TList *volList = zdcGeo->getListOfVolumes();

   // store geo parameter
   FairRun       *fRun       = FairRun::Instance();
   FairRuntimeDb *rtdb       = FairRun::Instance()->GetRuntimeDb();
   MpdZdcGeoPar  *par        = (MpdZdcGeoPar *)(rtdb->getContainer("MpdZdcGeoPar"));
   TObjArray     *fSensNodes = par->GetGeoSensitiveNodes();
   TObjArray     *fPassNodes = par->GetGeoPassiveNodes();

   TListIter      iter(volList);
   FairGeoNode   *node = NULL;
   FairGeoVolume *aVol = NULL;

   while ((node = (FairGeoNode *)iter.Next())) {
      aVol = dynamic_cast<FairGeoVolume *>(node);

      if (node->isSensitive()) {
         fSensNodes->AddLast(aVol);
      } else {
         fPassNodes->AddLast(aVol);
      }
   }
   par->setChanged();
   par->setInputVersion(fRun->GetRunId(), 1);

   ProcessNodes(volList);
}

// Check if Sensitive-----------------------------------------------------------
Bool_t MpdZdc::CheckIfSensitive(std::string name)
{
   TString tsname = name;
   if (tsname.Contains("zdc01s") || tsname.Contains("ScH")) {
      return kTRUE;
   }
   return kFALSE;
}
//-----------------------------------------------------------------------------

// -----   Private method AddHit   --------------------------------------------
MpdZdcPoint *MpdZdc::AddHit(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoZdc, TVector3 pos,
                            TVector3 mom, Double_t time, Double_t length, Double_t eLoss)
{
   TClonesArray &clref = *fZdcCollection;
   Int_t         size  = clref.GetEntriesFast();

   // cout <<"AddHit " <<trackID <<" " <<detID <<" " <<copyNoMother <<" " <<copyNo <<" " <<pos.Z() <<" " <<eLoss <<endl;
   // cout <<"-------------------------- " <<endl;

   // return new(clref[size]) MpdZdcPoint(trackID, detID, copyNo, copyNoMother,copyNo_h, copyNoMother_h,pos, mom, time,
   // length, eLoss);
   return new (clref[size]) MpdZdcPoint(trackID, detID, copyNo, copyNoMother, copyNoZdc, pos, mom, time, length, eLoss);
}
// ----

ClassImp(MpdZdc);
