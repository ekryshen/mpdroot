//
// One-dimensional correlation function for the model estimations
//

// MpdFemtoMaker headers
#include "MpdFemtoModelGausLCMSFreezeOutGenerator.h"
#include "MpdFemtoModelHiddenInfo.h"
#include "MpdFemtoModelManager.h"
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoModelQinvCorrFctn.h"

// ROOT headers
#include "TH1D.h"
#include "TString.h"

ClassImp(MpdFemtoModelQinvCorrFctn)

//_________________
MpdFemtoModelQinvCorrFctn::MpdFemtoModelQinvCorrFctn() :
MpdFemtoBaseCorrFctn(),
mManager(nullptr),
mNumeratorWeighted(nullptr),
mNumerator(nullptr),
mDenominatorWeighted(nullptr),
mDenominator(nullptr),
mIsUseDenominator(true) {
  // Default constructor
  mNumeratorWeighted = new TH1D("ModelNumWeighted", "Weighted numerator;q_{inv} (GeV/c);dN/dq_{inv}", 80, 0., 0.8);
  mNumerator = new TH1D("ModelNum", "Unweighted numerator;q_{inv} (GeV/c);dN/dq_{inv}", 80, 0., 0.8);

  mDenominatorWeighted = new TH1D("ModelDenWeighted", "ModelNumTrueIdeal", 80, 0., 0.8);
  mDenominator = new TH1D("ModelDen", "ModelNumFakeIdeal", 80, 0., 0.8);

  mNumeratorWeighted->Sumw2();
  mNumerator->Sumw2();

  mDenominatorWeighted->Sumw2();
  mDenominator->Sumw2();
}

//_________________
MpdFemtoModelQinvCorrFctn::MpdFemtoModelQinvCorrFctn(const char *title, const int& aNbins,
						     const double& aQinvLo, const double& aQinvHi,
						     const bool useDenominator) :
MpdFemtoBaseCorrFctn(),
mManager(nullptr),
mNumeratorWeighted(nullptr),
mNumerator(nullptr),
mDenominatorWeighted(nullptr),
mDenominator(nullptr),
mIsUseDenominator(useDenominator) {
  // Parametrized constructor
  TString baseName;
  TString baseTitle;
  TString appendix = ";q_{inv} (GeV/c);Entries";

  baseName = title;
  baseName += "_num_wei";
  baseTitle = baseName;
  baseTitle += appendix;
  mNumeratorWeighted = new TH1D(baseName.Data(), baseTitle.Data(),
				aNbins, aQinvLo, aQinvHi);
  mNumeratorWeighted->Sumw2();

  baseName = title;
  baseName += "_num";
  baseTitle = baseName;
  baseTitle += appendix;
  mNumerator = new TH1D(baseName.Data(), baseTitle.Data(),
			aNbins, aQinvLo, aQinvHi);
  mNumerator->Sumw2();

  if ( mIsUseDenominator ) {
    baseName = title;
    baseName += "_den_wei";
    baseTitle = baseName;
    baseTitle += appendix;
    mDenominatorWeighted = new TH1D(baseName.Data(), baseTitle.Data(),
				    aNbins, aQinvLo, aQinvHi);
    mDenominatorWeighted->Sumw2();

    baseName = title;
    baseName += "_den";
    baseTitle = baseName;
    baseTitle += appendix;
    mDenominator = new TH1D(baseName.Data(), baseTitle.Data(),
			    aNbins, aQinvLo, aQinvHi);
    mDenominator->Sumw2();
  }
}

//_________________
MpdFemtoModelQinvCorrFctn::MpdFemtoModelQinvCorrFctn(const MpdFemtoModelQinvCorrFctn& corrFctn) :
MpdFemtoBaseCorrFctn(),
mManager(corrFctn.mManager),
mNumeratorWeighted(nullptr),
mNumerator(nullptr),
mDenominatorWeighted(nullptr),
mDenominator(nullptr),
mIsUseDenominator(corrFctn.mIsUseDenominator) {
  // Copy constructor
  if (corrFctn.mNumeratorWeighted) {
    mNumeratorWeighted = new TH1D(*(corrFctn.mNumeratorWeighted));
  }
  if (corrFctn.mNumerator) {
    mNumerator = new TH1D(*(corrFctn.mNumerator));
  }

  if ( mIsUseDenominator ) {
    if (corrFctn.mDenominatorWeighted) {
      mDenominatorWeighted = new TH1D(*(corrFctn.mDenominatorWeighted));
    }
    if (corrFctn.mDenominator) {
      mDenominator = new TH1D(*(corrFctn.mDenominator));
    }
  }
}

//_________________
MpdFemtoModelQinvCorrFctn::~MpdFemtoModelQinvCorrFctn() {
  // Destructor
  if (mNumeratorWeighted) {
    delete mNumeratorWeighted;
    mNumeratorWeighted = nullptr;
  }
  if (mNumerator) {
    delete mNumerator;
    mNumerator = nullptr;
  }

  if ( mIsUseDenominator ) {
    if (mDenominatorWeighted) {
      delete mDenominatorWeighted;
      mDenominatorWeighted = nullptr;
    }
    if (mDenominator) {
      delete mDenominator;
      mDenominator = nullptr;
    }
  }
  delete mManager;
}

//_________________
MpdFemtoModelQinvCorrFctn& MpdFemtoModelQinvCorrFctn::operator=(const MpdFemtoModelQinvCorrFctn& corrFctn) {
  // Assignment operator
  if (this == &corrFctn) {
    return *this;
  }

  if (mNumeratorWeighted) delete mNumeratorWeighted;
  mNumeratorWeighted = ((corrFctn.mNumeratorWeighted) ? new TH1D(*corrFctn.mNumeratorWeighted) : nullptr);

  if (mNumerator) delete mNumerator;
  mNumerator = ((corrFctn.mNumerator) ? new TH1D(*corrFctn.mNumerator) : nullptr);

  if ( mIsUseDenominator ) {
    if (mDenominatorWeighted) delete mDenominatorWeighted;
    mDenominatorWeighted = ((corrFctn.mDenominatorWeighted) ? new TH1D(*corrFctn.mDenominatorWeighted) : nullptr);

    if (mDenominator) delete mDenominator;
    mDenominator = ((corrFctn.mDenominator) ? new TH1D(*corrFctn.mDenominator) : nullptr);
  }

  mManager = corrFctn.mManager;

  return *this;
}

//_________________
void MpdFemtoModelQinvCorrFctn::connectToManager(MpdFemtoModelManager *manager) {
  mManager = manager;
}

//_________________
MpdFemtoString MpdFemtoModelQinvCorrFctn::report() {
  // Make a report
  MpdFemtoString tStr = "MpdFemtoModelQinvCorrFctn report";

  tStr += TString::Format("Number of entries in numerator (w/ weight)   :\t%E\n",
			  mNumeratorWeighted->GetEntries());
  tStr += TString::Format("Number of entries in numerator (w/o weight)  :\t%E\n",
			  mNumerator->GetEntries());
  
  if ( mIsUseDenominator ) {
    tStr += TString::Format("Number of entries in denominator (w/ weight) :\t%E\n",
			    mDenominatorWeighted->GetEntries());
    tStr += TString::Format("Number of entries in denominator (w/o weight):\t%E\n",
			    mDenominator->GetEntries());
  }

  if (mPairCut) {
    tStr += "Here is the PairCut specific to this CorrFctn\n";
    tStr += mPairCut->report();
  } else {
    tStr += "No PairCut specific to this CorrFctn\n";
  }

  return tStr;
}

//_________________
void MpdFemtoModelQinvCorrFctn::addRealPair(MpdFemtoPair* pair) {
  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  // Retrieve femtoscopic weight
  double weight = mManager->weight(pair);

  mNumeratorWeighted->Fill( fabs(pair->qInv()), weight);
  mNumerator->Fill( fabs(pair->qInv()), 1.);
}

//_________________
void MpdFemtoModelQinvCorrFctn::addMixedPair(MpdFemtoPair* pair) {

  if ( !mIsUseDenominator ) return;
  
  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  double weight = mManager->weight(pair);
  mDenominatorWeighted->Fill( fabs(pair->qInv()), weight);
  mDenominator->Fill( fabs(pair->qInv()), 1.);
}

//_________________
void MpdFemtoModelQinvCorrFctn::eventBegin(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctn::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctn::finish() {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctn::writeOutHistos() {
  // Write out data histos
  mNumeratorWeighted->Write();
  mNumerator->Write();
  if ( mIsUseDenominator ) {
    mDenominatorWeighted->Write();
    mDenominator->Write();
  }
}

//_________________
TList* MpdFemtoModelQinvCorrFctn::getOutputList() {

  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(mNumeratorWeighted);
  tOutputList->Add(mNumerator);
  if ( mIsUseDenominator ) {
    tOutputList->Add(mDenominatorWeighted);
    tOutputList->Add(mDenominator);
  }

  return tOutputList;
}
