/*
 * MpdCosmicRayGenerator.cxx
 *
 *  Created on: 17 lis 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdCosmicRayGenerator.h"

#include <FairLogger.h>
#include <FairMCEventHeader.h>
#include <FairPrimaryGenerator.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <iostream>

MpdCosmicRayGenerator::MpdCosmicRayGenerator(TString filename, Int_t multi) :
  fFile(nullptr),
  fTree(nullptr),
  fMultiplicity(multi),
  fFileName(filename),
  fPID(0),
  fPosX(0),
  fPosY(0),
  fPosZ(0),
  fPx(0),
  fPy(0),
  fPz(0),
  fShift(0),
  fTime(-1),
  fGeneratedTracks(0) {}

Bool_t MpdCosmicRayGenerator::ReadEvent(FairPrimaryGenerator* primGen) {
  Int_t start               = fGeneratedTracks;
  FairMCEventHeader* header = primGen->GetEvent();
  TVector3 pos;
  header->GetVertex(pos);
  TVector3 rotVect(header->GetRotX(), header->GetRotY(), 0);
  Double_t rotZ = header->GetRotZ();
  for (int i = start; i < start + fMultiplicity; i++) {
    fTree->GetEntry(i);
    if (i > fTree->GetEntriesFast()) return kFALSE;
    fPID = 13;
    TVector3 position(fPosX * 100 - pos.X(), fPosY * 100 - fShift - pos.Y(), fPosZ * 100 - pos.Z());
    TVector3 mom(fPx, fPy, fPz);
    if (rotZ != 0) { pos.RotateZ(-rotZ); }
    if (rotVect.X() != 0) {
      mom.RotateUz(-rotVect);
      pos.RotateUz(-rotVect);
    }

    primGen->AddTrack(fPID, mom.Px(), mom.Py(), mom.Pz(), position.X(), position.Y(), position.Z());
    ++fGeneratedTracks;
  }
  return kTRUE;
}

Bool_t MpdCosmicRayGenerator::Init() {
  fFile = new TFile(fFileName);
  fTree = (TTree*) fFile->Get("CosmicCube");
  if (fTree == nullptr) return kFALSE;
  TBranch* b_PID;     //!
  TBranch* b_wallID;  //!
  TBranch* b_posX;    //!
  TBranch* b_posY;    //!
  TBranch* b_posZ;    //!
  TBranch* b_px;      //!
  TBranch* b_py;      //!
  TBranch* b_pz;      //!
  TBranch* b_pr;      //!
  TBranch* b_ptheta;  //!
  TBranch* b_pphi;    //!
  fTree->SetBranchAddress("PID", &fPID, &b_PID);
  fTree->SetBranchAddress("posX", &fPosX, &b_posX);
  fTree->SetBranchAddress("posY", &fPosY, &b_posY);
  fTree->SetBranchAddress("posZ", &fPosZ, &b_posZ);
  fTree->SetBranchAddress("px", &fPx, &b_px);
  fTree->SetBranchAddress("py", &fPy, &b_py);
  fTree->SetBranchAddress("pz", &fPz, &b_pz);

  TH1D* h          = (TH1D*) fFile->Get("SimulationParameters");
  Double_t simTime = h->GetBinContent(1);
  Double_t edge    = h->GetBinContent(2);
  fShift           = edge * 50.0;
  Double_t m       = fMultiplicity;
  Double_t tracks  = fTree->GetEntriesFast();

  if (fTime > 0) { fMultiplicity = fTime * tracks / simTime; }

  m                    = fMultiplicity;
  Double_t simPerEvent = m / tracks * simTime;
  LOG(info) << "MpdCosmicRayGenerator: simulation time per event: " << simPerEvent << " s";
  LOG(info) << "MpdCosmicRayGenerator: number of tracks per event: " << fMultiplicity;
  LOG(info) << "MpdCosmicRayGenerator: edge size: " << edge << " m";

  return kTRUE;
}

MpdCosmicRayGenerator::~MpdCosmicRayGenerator() {
  if (fFile) delete fFile;
}
