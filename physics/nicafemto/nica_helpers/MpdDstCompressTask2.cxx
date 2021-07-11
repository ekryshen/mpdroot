/*
 * MpdDstCompressTask2.cxx
 *
 *  Created on: 31 mar 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdDstCompressTask2.h"

#include <FairRootManager.h>
#include <FairLogger.h>
#include <TCollection.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TVector3.h>

#include "MpdKalmanTrack.h"
#include "MpdTrack.h"
#include "NicaHelix.h"

MpdDstCompressTask2::MpdDstCompressTask2() :
fUseMC(kFALSE),
fUseTpcKalmans(kFALSE),fUseTpcHits(kFALSE),fUseHeader(kFALSE),
fMCCompression(kFALSE),fFixMom(kFALSE),
fMpdEvent(nullptr),fMCTracks(nullptr),
fTpcKalmans(nullptr),fTpcHits(nullptr),fEventHeader(nullptr),
fMCMapSize(0),fMCIndexMap(nullptr),fPID(nullptr){
    // TODO Auto-generated constructor stub

}

InitStatus MpdDstCompressTask2::CheckBranches() {
    FairRootManager *mngr = FairRootManager::Instance();
    fMpdEvent =(MpdEvent*)mngr->GetObject("MPDEvent.");
    if(fMpdEvent==nullptr) return kFATAL;
    mngr->Register("MPDEvent.", "MPD", fMpdEvent, kTRUE);
    if(fUseMC){
        fMCTracks = (TClonesArray*)mngr->GetObject("MCTrack");
        if(fMCTracks==nullptr){
            LOG(WARNING)<<"MC tracks requested but not found!"<<FairLogger::endl;
            fUseMC=NULL;
        }else{
            mngr->Register("MCTrack", "MC",fMCTracks, kTRUE);
        }
    }
    if(fUseTpcKalmans){
        fTpcKalmans = (TClonesArray*)mngr->GetObject("TpcKalmanTrack");
        if(fTpcKalmans==nullptr){
            LOG(WARNING)<<"Kalman TPC tracks requested but not found!"<<FairLogger::endl;
            fUseTpcKalmans=nullptr;
        }else{
            mngr->Register("TpcKalmanTrack",  "TPC", fTpcKalmans,kTRUE);
        }
    }else{
        fTpcKalmans = (TClonesArray*)mngr->GetObject("TpcKalmanTrack");
    }
    if(fTpcKalmans==nullptr)
    fFixMom = kFALSE;
    if(fUseTpcHits){
        fTpcHits = (TClonesArray*)mngr->GetObject("TpcHit");
        if(fTpcHits==nullptr){
            LOG(WARNING)<<"TPC hits requested but not found!"<<FairLogger::endl;
            fUseTpcHits=nullptr;
        }else{
            mngr->Register("TpcHi.",  "TPC", fTpcHits,kTRUE);
        }
    }
    if(!fUseHeader){
    /*   TODO when FairROOT will be ugpraded
        FairRunAna::Instance()->SetEventHeaderPersistence(kFALSE);
        */
    }
    return kSUCCESS;
}

InitStatus MpdDstCompressTask2::Init() {
    NicaHelix::SetMagField(0.5);
    if(CheckBranches()==kFATAL){
        LOG(FATAL)<<"End of macro MPDEvent not found"<<FairLogger::endl;
        return kFATAL;
    }
    fPID = new MpdPid(0, 0, 0, 1, "DEFAULT", "CF", "pi+ka+el+pr");
    fMCMapSize = 1000;
    fMCIndexMap = new Int_t[fMCMapSize];
    return kSUCCESS;
}

void MpdDstCompressTask2::Exec(Option_t *option) {
    Fix();
    Rewrite();
}

void MpdDstCompressTask2::Rewrite() {
    if(fMCCompression){
        TClonesArray *glob_tracks = fMpdEvent->GetGlobalTracks();
        TClonesArray *prim_tracks = fMpdEvent->GetPrimaryTracks();
        if(fMCTracks->GetEntriesFast()>fMCMapSize){
            delete []fMCIndexMap;
            fMCMapSize = fMCTracks->GetSize()*2;
            fMCIndexMap = new Int_t[fMCMapSize];
        }
        for(int iMCTrack=0;iMCTrack<fMCTracks->GetEntriesFast();iMCTrack++){
            fMCIndexMap[iMCTrack] = -1;
        }
        for(int iTrack=0;iTrack<glob_tracks->GetEntriesFast();iTrack++){
            MpdTrack *track = (MpdTrack*)glob_tracks->UncheckedAt(iTrack);
            Int_t matched = track->GetID();
            if(matched>=0)
                fMCIndexMap[matched]=0;
        }
        for(int iTrack=0;iTrack<prim_tracks->GetEntriesFast();iTrack++){
            MpdTrack *track = (MpdTrack*)prim_tracks->UncheckedAt(iTrack);
            Int_t matched = track->GetID();
            if(matched>=0)
                fMCIndexMap[matched]=0;
        }
        Int_t index = 0;
        for(int iMCTrack=0;iMCTrack<fMCTracks->GetEntriesFast();iMCTrack++){
            if(fMCIndexMap[iMCTrack]>-1){//this track is matched
                fMCIndexMap[iMCTrack] = index++;
            }else{
                fMCTracks->RemoveAt(iMCTrack);
            }
        }
        fMCTracks->Compress();
        //set new matching ID's
        for(int iTrack=0;iTrack<glob_tracks->GetEntriesFast();iTrack++){
            MpdTrack *track = (MpdTrack*)glob_tracks->UncheckedAt(iTrack);
            Int_t matched = track->GetID();
            if(matched>-1)
                track->SetID(fMCIndexMap[matched]);
        }
        for(int iTrack=0;iTrack<prim_tracks->GetEntriesFast();iTrack++){
            MpdTrack *track = (MpdTrack*)prim_tracks->UncheckedAt(iTrack);
            Int_t matched = track->GetID();
            if(matched>-1)
                track->SetID(fMCIndexMap[matched]);
        }
    }
}

void MpdDstCompressTask2::Fix() {
    TClonesArray *glob_tracks = fMpdEvent->GetGlobalTracks();
    TClonesArray *prim_tracks = fMpdEvent->GetPrimaryTracks();
    TVector3 vertex(fMpdEvent->GetPrimaryVerticesX(),
            fMpdEvent->GetPrimaryVerticesY(),
            fMpdEvent->GetPrimaryVerticesZ());

    for(int i=0;i<glob_tracks->GetEntriesFast();i++){
        MpdTrack *track = (MpdTrack*)glob_tracks->UncheckedAt(i);
        FixPID(track);
        FixMom(track, i, vertex);
    }
}

void MpdDstCompressTask2::FixPID(MpdTrack *track) {
    TVector3 mom(track->GetPx(),track->GetPy(),track->GetPz());
    Double_t p = mom.Mag();
    Double_t dedx = track->GetdEdXTPC();
    Double_t dedx_el = fPID->GetDedxElParam(p);
    Double_t dedx_pi = fPID->GetDedxPiParam(p);
    Double_t dedx_ka = fPID->GetDedxKaParam(p);
    Double_t dedx_pr = fPID->GetDedxPrParam(p);
    Double_t sigma_el = fPID->GetDedxWidthValue(p, 4)*dedx_el;
    Double_t sigma_pi = fPID->GetDedxWidthValue(p, 1)*dedx_pi;
    Double_t sigma_ka = fPID->GetDedxWidthValue(p, 2)*dedx_ka;
    Double_t sigma_pr = fPID->GetDedxWidthValue(p, 3)*dedx_pr;
    sigma_el = (dedx-dedx_el)/(sigma_el);
    sigma_pi = (dedx-dedx_pi)/(sigma_pi);
    sigma_ka = (dedx-dedx_ka)/(sigma_ka);
    sigma_pr = (dedx-dedx_pr)/(sigma_pr);
    if(TMath::IsNaN(sigma_el))
        sigma_el = -1E+2;
    if(TMath::IsNaN(sigma_pi))
        sigma_pi = -1E+2;
    if(TMath::IsNaN(sigma_ka))
        sigma_ka = -1E+2;
    if(TMath::IsNaN(sigma_pr))
        sigma_pr = -1E+2;
//  std::cout<<sigma_pi<<" "<<dedx<<" "<<dedx_pi<<" "<<sigma_pi<<std::endl;
    track->SetNSigmaElectron(sigma_el);
    track->SetNSigmaKaon(sigma_ka);
    track->SetNSigmaPion(sigma_pi);
    track->SetNSigmaProton(sigma_pr);
    if(track->GetTofMass2()==0){
        track->SetTofMass2(-5);
    }

}

void MpdDstCompressTask2::FixMom(MpdTrack *track,Int_t index, TVector3 &vertex) {
    MpdKalmanTrack *kftrack = (MpdKalmanTrack*) fTpcKalmans->UncheckedAt(index);

    if (kftrack->GetParam(4) == 0.) track->SetPt(TMath::Sqrt(-1)); /*NaN*/
    else track->SetPt(1. / kftrack->GetParam(4)); /*signed Pt*/

    track->SetTheta(TMath::PiOver2() - kftrack->GetParam(3)); // Theta: angle from beam line
    track->SetPhi(kftrack->GetParam(2)); // Phi
    Double_t phi = kftrack->GetParam(0) / kftrack->GetPosNew();
    track->SetFirstPointX(kftrack->GetPosNew() * TMath::Cos(phi)); // closest to beam line
    track->SetFirstPointY(kftrack->GetPosNew() * TMath::Sin(phi));
    track->SetFirstPointZ(kftrack->GetParam(1));
    TVector3 poz(kftrack->GetPosNew() * TMath::Cos(phi),
            kftrack->GetPosNew() * TMath::Sin(phi),
    kftrack->GetParam(1)
    );
    TVector3 mom(track->GetPx(),track->GetPy(),track->GetPz());
    NicaHelixZ helix(poz,mom,track->GetCharge());
    Double_t s = helix.PathLength(vertex, 1);
    TVector3 newMom = helix.EvalMom(s);
    TVector3 newPos = helix.EvalPos(s);
    track->SetFirstPointX(newPos.X());
    track->SetFirstPointY(newPos.Y());
    track->SetFirstPointZ(newPos.Z());
    if(track->GetCharge()>0){
        track->SetPt(-newMom.Pt());
    }else{
        track->SetPt(newMom.Pt());
    }
    track->SetPhi(newMom.Phi());
    track->SetTheta(newMom.Theta());
// recalculate DCA
    track->SetDCAX(newPos.X()-vertex.X());
    track->SetDCAY(newPos.Y()-vertex.Y());
    track->SetDCAZ(newPos.Z()-vertex.Z());
}

MpdDstCompressTask2::~MpdDstCompressTask2() {
    // TODO Auto-generated destructor stub
}

