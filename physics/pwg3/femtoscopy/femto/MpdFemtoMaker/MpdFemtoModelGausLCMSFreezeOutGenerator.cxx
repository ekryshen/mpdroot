/**
 * MpdFemtoModelGausLCMSFreezeOutGenerator - freeze-out
 * coordinates generator, generating a 3D gaussian ellipsoid in LCMS
 */

// C++ headers
#include <math.h>
#include <iostream>

// MpdFemtoMaker headers
#include "MpdFemtoModelGausLCMSFreezeOutGenerator.h"
#include "MpdFemtoModelHiddenInfo.h"

// ROOT headers
#include "TMath.h"
#include "TLorentzVector.h"

ClassImp(MpdFemtoModelGausLCMSFreezeOutGenerator);

//_________________
MpdFemtoModelGausLCMSFreezeOutGenerator::MpdFemtoModelGausLCMSFreezeOutGenerator() :
mSizeOut(0), mSizeSide(0), mSizeLong(0) {
  // Default constructor
  mRandom = new TRandom3();
}

//_________________
MpdFemtoModelGausLCMSFreezeOutGenerator::MpdFemtoModelGausLCMSFreezeOutGenerator(const MpdFemtoModelGausLCMSFreezeOutGenerator &aModel) :
MpdFemtoBaseModelFreezeOutGenerator(aModel),
mSizeOut(0), mSizeSide(0), mSizeLong(0) {
  // Copy constructor
  mRandom = new TRandom3();
  setSizeOut(aModel.sizeOut());
  setSizeSide(aModel.sizeSide());
  setSizeLong(aModel.sizeLong());
}

//_________________
MpdFemtoModelGausLCMSFreezeOutGenerator& MpdFemtoModelGausLCMSFreezeOutGenerator::operator=(const MpdFemtoModelGausLCMSFreezeOutGenerator &aModel) {
  // Assignment operator
  if (this != &aModel) {
    mRandom = new TRandom3();
    setSizeOut(aModel.sizeOut());
    setSizeSide(aModel.sizeSide());
    setSizeLong(aModel.sizeLong());
  }

  return *this;
}

//_________________
MpdFemtoModelGausLCMSFreezeOutGenerator::~MpdFemtoModelGausLCMSFreezeOutGenerator() {
  if (mRandom) delete mRandom;
}

//_________________
void MpdFemtoModelGausLCMSFreezeOutGenerator::generateFreezeOut(MpdFemtoPair *aPair) {
  // Generate two-particle emission points with respect
  // to their pair momentum
  // The source is the 3D Gaussian ellipsoid in the LCMS frame

  //MpdFemtoModelHiddenInfo *inf1 = (MpdFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //MpdFemtoModelHiddenInfo *inf2 = (MpdFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();

  // Mike Lisa
  MpdFemtoTrack *inf1 = (MpdFemtoTrack *) aPair->track1()->track();
  MpdFemtoTrack *inf2 = (MpdFemtoTrack *) aPair->track2()->track();

  if ((!inf1) || (!inf2)) {
    std::cout << "Hidden info not created! " << std::endl;
    exit(kFALSE);
  }

  //std::cout<<" we are in Freeze-out Generator inf1 inf2  "<<inf1<<"  "<<inf2<<std::endl;
  //std::cout<<" inf1 GetMass "<<((MpdFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid()<<std::endl;
  //std::cout<<" true mom    " <<((MpdFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->x()<<std::endl;

  double tPx = (((MpdFemtoModelHiddenInfo*) inf1->hiddenInfo())->trueMomentum().X() +
		((MpdFemtoModelHiddenInfo*) inf2->hiddenInfo())->trueMomentum().X());
  double tPy = (((MpdFemtoModelHiddenInfo*) inf1->hiddenInfo())->trueMomentum().Y() +
		((MpdFemtoModelHiddenInfo*) inf2->hiddenInfo())->trueMomentum().Y());
  double tPz = (((MpdFemtoModelHiddenInfo*) inf1->hiddenInfo())->trueMomentum().Z() +
		((MpdFemtoModelHiddenInfo*) inf2->hiddenInfo())->trueMomentum().Z());

  // double tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
  // double tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();

  //std::cout<<" tPx tPy tPz"<<tPx<<" "<<tPy<<" "<<tPz<<std::endl;
  if (!(tPx == 0 && tPy == 0 && tPz == 0)) {

    double tM1 = ((MpdFemtoModelHiddenInfo*) inf1->hiddenInfo())->mass();
    double tM2 = ((MpdFemtoModelHiddenInfo*) inf2->hiddenInfo())->mass();
    double tE1 = TMath::Sqrt(tM1 * tM1 +
			     ((MpdFemtoModelHiddenInfo*) inf1->hiddenInfo())->trueMomentum().Mag2());
    double tE2 = TMath::Sqrt(tM2 * tM2 +
			     ((MpdFemtoModelHiddenInfo*) inf2->hiddenInfo())->trueMomentum().Mag2());
    double tEs = tE1 + tE2;

    //std::cout<<" tM1 tM2 tE1 tE2"<<tM1<<" "<<tM2<<" "<<tE2<<std::endl;
    double tPt = TMath::Sqrt(tPx * tPx + tPy * tPy);

    double tRout = mRandom->Gaus(0.0, mSizeOut);
    double tRside = mRandom->Gaus(0.0, mSizeSide);
    double tRlong = mRandom->Gaus(0.0, mSizeLong);

    double tXout = (tPx * tRout + tPy * tRside) / tPt;
    double tXside = (tPy * tRout - tPx * tRside) / tPt;

    double tBetaz = tPz / tEs;
    double tGammaz = 1.0 / TMath::Sqrt(1 - tBetaz * tBetaz);

    double tXlong = tGammaz * (tRlong + tBetaz * 0);
    double tXtime = tGammaz * (0 + tBetaz * tRlong);

    //std::cout<<" tXout tXside before hidden infor "<<tXout<<" "<<tXside<<std::endl;

    inf1->hiddenInfo()->setEmissionPoint(0., 0., 0., 0.);
    inf2->hiddenInfo()->setEmissionPoint(tXout, tXside, tXlong, tXtime);
    //std::cout<<" after all tXout tXside "<<tXout<<" "<<tXside<<std::endl;

  } // if ( !(tPx==0 && tPy==0 && tPz==0 ) ) {
}
