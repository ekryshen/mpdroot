/////////////////////////////////////////////////////////////
//
//  MpdEmcHitProducer
//
//  Filler of MpdEmcHit
//
/////////////////////////////////////////////////////////////// 


#include "TClonesArray.h"
#include "FairRootManager.h"
#include "TMath.h"
#include "FairRunAna.h"

#include "MpdEmcHitProducer.h"
#include "MpdEmcGeoPar.h"
#include "MpdEmcHit.h"
#include "MpdEmcPoint.h"
#include "FairMCPoint.h"
#include "MpdMCTrack.h"

using namespace std;
using namespace TMath;

// -----   Default constructor   -------------------------------------------

MpdEmcHitProducer::MpdEmcHitProducer() :
FairTask("Ideal EMC hit Producer") {

}


// -----   Destructor   ----------------------------------------------------

MpdEmcHitProducer::~MpdEmcHitProducer() {
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcHitProducer::Init() {

    cout << "******************* EMC INIT *********************" << endl;

    // Get RootManager
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) {
        cout << "-E- MpdEmcHitProducer::Init: "
                << "RootManager not instantiated!" << endl;
        return kFATAL;
    }

    fGeoPar = new MpdEmcGeoPar();

    // Get input array
    fPointArray = (TClonesArray*) ioman->GetObject("MpdEmcPoint");
    if (!fPointArray) {
        cout << "-W- MpdEmcHitProducer::Init: " << "No EmcPoint array!" << endl;
        return kERROR;
    }
    
    fMcTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMcTrackArray) {
        cout << "-W- MpdEmcHitProducer::Init: " << "No MCTrack array!" << endl;
        return kERROR;
    }

    // Create and register output array
    fDigiArray = new TClonesArray("MpdEmcHit", 100);
    ioman->Register("MpdEmcHit", "EMC", fDigiArray, kTRUE);

    cout << "-I- MpdEmcHitProducer: Intialization successfull" << endl;

    return kSUCCESS;

}

void MpdEmcHitProducer::Finish() {
    cout << "-I- MpdEmcHitProducer: Finish" << endl;
}

void MpdEmcHitProducer::Exec(Option_t* opt) {
    cout << "MpdEmcHitProducer::Exec started" << endl;

    // Reset output Array
    if (!fDigiArray) Fatal("MpdEmcHitProducer::Exec)", "No array of digits");
    fDigiArray->Delete();

    for (UInt_t iPnt = 0; iPnt < (UInt_t) fPointArray->GetEntriesFast(); ++iPnt) {
        FairMCPoint* pnt = (FairMCPoint*) fPointArray->At(iPnt);
        Int_t trId = pnt->GetTrackID();
        MpdMCTrack* tr = (MpdMCTrack*) fMcTrackArray->At(trId);
//        if (tr->GetMotherId() != -1) continue;
        Int_t pdg = tr->GetPdgCode();
        Float_t x = pnt->GetX();
        Float_t y = pnt->GetY();
        Float_t z = pnt->GetZ();
        UInt_t sec = GetSecId(x, y, z);
        UInt_t row = GetRowId(z);
        UInt_t supMod = GetSupModId(x, y, z, sec);
        UInt_t mod = GetModId(x, y, supMod, sec);
        Float_t e = pnt->GetEnergyLoss();

        MpdEmcHit* hit = SearchHit(sec, supMod, row, mod);
        if (hit == NULL) {
            hit = new((*fDigiArray)[fDigiArray->GetEntriesFast()]) MpdEmcHit(sec, row, supMod, mod, e);
        } else {
            hit->IncreaseEnergy(e);
        }
        hit->SetNumTracks(hit->GetNumTracks() + 1);
        hit->SetTrackId(trId);
        hit->SetPdg(pdg);
        hit->SetZCenter(CalcZCenter(sec, row, mod));
        hit->SetPhiCenter(CalcPhiCenter(sec, supMod, mod));
    }
    
    for (UInt_t iHit = 0; iHit < (UInt_t) fDigiArray->GetEntriesFast(); ++iHit) {
        MpdEmcHit* hit = (MpdEmcHit*) fDigiArray->At(iHit);
        if (hit->GetNumTracks() > 1) {
            hit->SetPdg(0);
            hit->SetTrackId(-1);
        }
        //hit->Print();
    }    
}

UInt_t MpdEmcHitProducer::GetSecId(Float_t x, Float_t y, Float_t z) {
    Float_t ang = ATan2(y, x);
    if (ang < 0) ang += TwoPi();
    ang *= RadToDeg();
    Int_t nSec = fGeoPar->GetNsec();
    if (z > 0.0)
        return UInt_t(ang / 360 * nSec);
    else
        return nSec + UInt_t(ang / 360 * nSec);
}

UInt_t MpdEmcHitProducer::GetSupModId(Float_t x, Float_t y, Float_t z, UInt_t sec) {
    Float_t ang = ATan2(y, x);
    if (ang < 0) ang += TwoPi();
    ang *= RadToDeg();
    Float_t nDegreesInOneSector = fGeoPar->GetAngleOfSector();
    if (z < 0.0)
        sec -= fGeoPar->GetNsec();

    Float_t secStartAng = sec * nDegreesInOneSector;
    //Float_t secFinishAng = secStartAng + nDegreesInOneSector;
    Float_t localAng = ang - secStartAng;

    return UInt_t(localAng * fGeoPar->GetNsupMod() / nDegreesInOneSector);
}

UInt_t MpdEmcHitProducer::GetModId(Float_t x, Float_t y, UInt_t supMod, UInt_t sec) {
    Float_t ang = ATan2(y, x);
    if (ang < 0) ang += TwoPi();
    ang *= RadToDeg();
        
    Float_t secStartAng = (sec % fGeoPar->GetNsec()) * fGeoPar->GetAngleOfSector();
    Float_t supModStartAng = secStartAng + supMod * fGeoPar->GetAngleOfSuperModule();
    Float_t localAng = ang - supModStartAng;
    
    return UInt_t(localAng * fGeoPar->GetNmod() / fGeoPar->GetAngleOfSuperModule());
}

UInt_t MpdEmcHitProducer::GetRowId(Float_t z) {
    return UInt_t(Abs(z) / fGeoPar->GetLength() * fGeoPar->GetNrows());
}

Float_t MpdEmcHitProducer::CalcZCenter(UInt_t sec, UInt_t row, UInt_t mod) {
    Float_t lengthOfModuleByZ = fGeoPar->GetLengthOfModuleByZ();
    Float_t lengthOfSuperModuleByZ = fGeoPar->GetLengthOfSuperModuleByZ();
    Float_t halfLengthOfModuleByZ = lengthOfModuleByZ / 2.0;
    Float_t leftEdgeOfModuleByZ = row * lengthOfSuperModuleByZ + (mod % fGeoPar->GetNModInSuperModByZ()) * lengthOfModuleByZ;
    Float_t z = leftEdgeOfModuleByZ + halfLengthOfModuleByZ;
    return (sec < fGeoPar->GetNsec()) ? z : -z;
}

Float_t MpdEmcHitProducer::CalcPhiCenter(UInt_t sec, UInt_t supMod, UInt_t mod) {
    Float_t sectorAngle = fGeoPar->GetAngleOfSector();
    Float_t supModAngle = fGeoPar->GetAngleOfSuperModule();
    Float_t modAngle = fGeoPar->GetAngleOfModule();
    
    Int_t modIdInSupModByPhi = mod / fGeoPar->GetNModInSuperModByPhi(); // 0, 1, 2
    Float_t sectorAngleEdge = sectorAngle * (sec % fGeoPar->GetNsec());
    Float_t supModAngleEdge = supModAngle * supMod;
    Float_t modAngleEdge = modAngle * modIdInSupModByPhi;
    return sectorAngleEdge + supModAngleEdge + modAngleEdge + modAngle / 2.0;
}

MpdEmcHit* MpdEmcHitProducer::SearchHit(UInt_t sec, UInt_t supMod, UInt_t row, UInt_t mod) {
    MpdEmcHit* foundHit = NULL;
    for (Int_t i = 0; i < fDigiArray->GetEntriesFast(); ++i) {
        MpdEmcHit* hit = (MpdEmcHit*) fDigiArray->At(i);
        if ((UInt_t)hit->GetSec() != sec) continue;
        if ((UInt_t)hit->GetSupMod() != supMod) continue;
        if ((UInt_t)hit->GetRow() != row) continue;
        if ((UInt_t)hit->GetMod() != mod) continue;
        foundHit = hit;
    }
    return foundHit;
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcHitProducer)