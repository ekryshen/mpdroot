/*
 * MpdCosmicRayGenerator.h
 *
 *  Created on: 17 lis 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_GENERATORS_MPDCOSMICRAYGENERATOR_H_
#define MPDROOT_GENERATORS_MPDCOSMICRAYGENERATOR_H_

#include "FairGenerator.h"

#include <TString.h>
#include <vector>

class TFile;
class TTree;

class MpdCosmicRayGenerator : public FairGenerator {
  TFile* fFile;  //!
  TTree* fTree;  //!
  Int_t fMultiplicity;
  Int_t fGeneratedTracks;
  Int_t fPID;
  Double_t fPosX;
  Double_t fPosY;
  Double_t fPosZ;
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fShift;
  Double_t fTime;
  TString fFileName;

public:
  MpdCosmicRayGenerator(TString filename = "", Int_t multi = 1000);
  /**
  set simulation time in seconds, ovewrites multiplicity
  @param t time in sec
  **/
  void SetTime(Double_t t) { fTime = t; }
  virtual Bool_t Init();
  virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);
  virtual ~MpdCosmicRayGenerator();
  ClassDef(MpdCosmicRayGenerator, 1)
};


#endif /* MPDROOT_GENERATORS_MPDCOSMICRAYGENERATOR_H_ */
