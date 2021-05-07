/*
 * MpdPIDOnTheFly.cxx
 *
 *  Created on: 19 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPIDOnTheFly.h"
#include "FairRootManager.h"
#include "MpdTrack.h"
#include "TVector3.h"
#include "FairRunAna.h"
#include "FairField.h"

#include <iostream>

MpdPIDOnTheFly::MpdPIDOnTheFly() :fEvent(nullptr),fPID(nullptr),fMiniTrack(nullptr){
	fEventVector = new TVector3();
	fMCVector = new TVector3();

}

InitStatus MpdPIDOnTheFly::Init() {
	FairRootManager *manager = FairRootManager::Instance();
	fEvent = (MpdEvent*)manager->GetObject("MPDEvent.");
	if(fEvent==nullptr){
	    fMiniTrack = (TClonesArray*)manager->GetObject("Track");
	}

	fPID = new MpdPid(0, 0, 0, 1, "DEFAULT", "CF", "pi+ka+el+pr");
	return kSUCCESS;
}

void MpdPIDOnTheFly::Exec(Option_t* opt) {
    if(fEvent){
        TClonesArray *tracks = fEvent->GetGlobalTracks();
        fEventVector->SetXYZ(fEvent->GetPrimaryVerticesX(),
        fEvent->GetPrimaryVerticesY(),
	    fEvent->GetPrimaryVerticesZ());
        for(int i=0;i<tracks->GetEntriesFast();i++){
            MpdTrack *track  = (MpdTrack*)tracks->UncheckedAt(i);
            FillTrackPID(track);
            FillTrackDCA(track, fEventVector,fMCVector);
        }
    }else{// this is minidst
        for(int i=0;i<fMiniTrack->GetEntriesFast();i++){
            MpdMiniTrack *track = (MpdMiniTrack*)fMiniTrack->UncheckedAt(i);
            FillTrackPID(track);
        }
    }
}

void MpdPIDOnTheFly::FillTrackPID(MpdTrack* track) {
	TVector3 mom(track->GetPx(),track->GetPy(),track->GetPz());
	Double_t p = mom.Mag();
	Double_t dedx = track->GetdEdXTPC();
	Double_t dedx_el = fPID->GetDedxElParam(p);
	Double_t dedx_pi = fPID->GetDedxPiParam(p);
	Double_t dedx_ka = fPID->GetDedxKaParam(p);
	Double_t dedx_pr = fPID->GetDedxPrParam(p);
	Double_t sigma_el = fPID->GetDedxWidthValue(p, MpdPidUtils::kElectron)*dedx_el;
	Double_t sigma_pi = fPID->GetDedxWidthValue(p, MpdPidUtils::kPion)*dedx_pi;
	Double_t sigma_ka = fPID->GetDedxWidthValue(p, MpdPidUtils::kKaon)*dedx_ka;
	Double_t sigma_pr = fPID->GetDedxWidthValue(p, MpdPidUtils::kProton)*dedx_pr;
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
//	std::cout<<sigma_pi<<" "<<dedx<<" "<<dedx_pi<<" "<<sigma_pi<<std::endl;
	track->SetNSigmaElectron(sigma_el);
	track->SetNSigmaKaon(sigma_ka);
	track->SetNSigmaPion(sigma_pi);
	track->SetNSigmaProton(sigma_pr);
	if(track->GetTofMass2()==0){
		track->SetTofMass2(-5);
	}
}

void MpdPIDOnTheFly::FillTrackDCA(MpdTrack* track, TVector3* recoVertex,
		TVector3* mcVertex) {
    MpdHelix helix = track->GetHelix();
    Double_t path_at_mcVertex;
    Double_t path_at_recoVertex;
    path_at_mcVertex = helix.pathLength(*mcVertex);
    path_at_recoVertex = helix.pathLength(*recoVertex);
    TVector3 DCA_MC = helix.at(path_at_mcVertex);
    TVector3 DCA_RECO = helix.at(path_at_recoVertex);
    // set dca global as dca to MC vertex DW
    track->SetDCAGlobalX(DCA_MC.X());
    track->SetDCAGlobalY(DCA_MC.Y());
    track->SetDCAGlobalZ(DCA_MC.Z());
    // set dca as dca to reconstructed vertex DW
    track->SetDCAX(DCA_RECO.X());
    track->SetDCAY(DCA_RECO.Y());
    track->SetDCAZ(DCA_RECO.Z());
}

MpdPIDOnTheFly::~MpdPIDOnTheFly() {
	// TODO Auto-generated destructor stub
}

void MpdPIDOnTheFly::FillTrackPID(MpdMiniTrack *track) {
    Double_t p = track->gPtot();
    Double_t dedx = track->dEdx();
    track->setDedx(track->dEdx());
    Double_t dedx_el = fPID->GetDedxElParam(p);
    Double_t dedx_pi = fPID->GetDedxPiParam(p);
    Double_t dedx_ka = fPID->GetDedxKaParam(p);
    Double_t dedx_pr = fPID->GetDedxPrParam(p);
    Double_t sigma_el = fPID->GetDedxWidthValue(p, MpdPidUtils::kElectron)*dedx_el;
    Double_t sigma_pi = fPID->GetDedxWidthValue(p, MpdPidUtils::kPion)*dedx_pi;
    Double_t sigma_ka = fPID->GetDedxWidthValue(p, MpdPidUtils::kKaon)*dedx_ka;
    Double_t sigma_pr = fPID->GetDedxWidthValue(p, MpdPidUtils::kProton)*dedx_pr;
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
    track->setNSigmaElectron(sigma_el);
    track->setNSigmaKaon(sigma_ka);
    track->setNSigmaPion(sigma_pi);
    track->setNSigmaProton(sigma_pr);
}

