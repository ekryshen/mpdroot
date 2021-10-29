//
// Three-dimensional Bertsch-Pratt correlation function in LCMS for the model estimations
//

// MpdFemtoMaker headers
#include "MpdFemtoModelGausLCMSFreezeOutGenerator.h"
#include "MpdFemtoModelHiddenInfo.h"
#include "MpdFemtoModelManager.h"
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoModelBPLCMS3DCorrFctn.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"
#include "TH3F.h"

// C++ headers
#include <cmath>

ClassImp(MpdFemtoModelBPLCMS3DCorrFctnKt)

//_________________
MpdFemtoModelBPLCMS3DCorrFctnKt::MpdFemtoModelBPLCMS3DCorrFctnKt(const char* title, const int& nBins,
								 const double& qLo, const double& qHi,
								 const bool isUseDenominator) :
MpdFemtoBaseCorrFctn(),
mManager(nullptr) {
  // Parametrized constructor

  // Set general parameters
  setUseDenominator( isUseDenominator );

  // Define string parameters
  TString baseName, baseTitle;
  TString appendix = ";q_{out} (GeV/c);q_{side} (GeV/c);q_{long} (GeV/c)";

  // Numerator
  baseName = "hBPLCMSModel_num_";
  baseName += title;
  baseTitle = baseName;
  baseTitle += appendix;
  mNumerator = new TH3F( baseName.Data(), baseTitle.Data(),
			 nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
  mNumerator->Sumw2();

  // Weighted numerator
  baseName = "hBPLCMSModel_numw_";
  baseName += title;
  baseTitle = baseName;
  baseTitle += appendix;
  mNumeratorWeighted = new TH3F( baseName.Data(), baseTitle.Data(),
				 nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
  mNumeratorWeighted->Sumw2();

  // Qinv weighted numerator
  baseName = "hBPLCMSModel_qinvw_";
  baseName += title;
  baseTitle = baseName;
  baseTitle += appendix;
  mNumeratorQinvWeighted = new TH3F( baseName.Data(), baseTitle.Data(),
				     nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
  mNumeratorQinvWeighted->Sumw2();

  if ( mIsUseDenominator ) {

    // Denominator
    baseName = "hBPLCMSModel_den_";
    baseName += title;
    baseTitle = baseName;
    baseTitle += appendix;
    mDenominator = new TH3F( baseName.Data(), baseTitle.Data(),
			     nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
    mDenominator->Sumw2();
    
    // Weighted denominator
    baseName = "hBPLCMSModel_denw_";
    baseName += title;
    baseTitle = baseName;
    baseTitle += appendix;
    mDenominatorWeighted = new TH3F( baseName.Data(), baseTitle.Data(),
				     nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
    mDenominatorWeighted->Sumw2();
    
    // Qinv weighted denominator
    baseName = "hBPLCMSModel_qinvw_";
    baseName += title;
    baseTitle = baseName;
    baseTitle += appendix;
    mDenominatorQinvWeighted = new TH3F( baseName.Data(), baseTitle.Data(),
					 nBins, qLo, qHi, nBins, qLo, qHi, nBins, qLo, qHi );
    mDenominatorQinvWeighted->Sumw2();
  } // if ( mIsUseDenominator )
}

//_________________
MpdFemtoModelBPLCMS3DCorrFctnKt::MpdFemtoModelBPLCMS3DCorrFctnKt(const MpdFemtoModelBPLCMS3DCorrFctnKt& copy) :
  MpdFemtoBaseCorrFctn( copy ),
  mManager( copy.mManager ),
  mIsUseDenominator( copy.mIsUseDenominator ) {
  // Copy constructor
  mNumerator = new TH3F( *copy.mNumerator );
  mNumeratorWeighted = new TH3F( *copy.mNumeratorWeighted );
  mNumeratorQinvWeighted = new TH3F( *copy.mNumeratorQinvWeighted );
  if ( mIsUseDenominator ) {
    mDenominator = new TH3F( *copy.mDenominator );
    mDenominatorWeighted = new TH3F( *copy.mDenominatorWeighted );
    mDenominatorQinvWeighted = new TH3F( *copy.mDenominatorQinvWeighted );
  }
}

//____________________________
MpdFemtoModelBPLCMS3DCorrFctnKt::~MpdFemtoModelBPLCMS3DCorrFctnKt() {
  // Destructor
  if ( mNumerator ) delete mNumerator;
  if ( mNumeratorWeighted ) delete mNumeratorWeighted;
  if ( mNumeratorQinvWeighted ) delete mNumeratorQinvWeighted;
  if ( mDenominator ) delete mDenominator;
  if ( mDenominatorWeighted ) delete mDenominatorWeighted;
  if ( mDenominatorQinvWeighted ) delete mDenominatorQinvWeighted;
}

//_________________________
MpdFemtoModelBPLCMS3DCorrFctnKt& MpdFemtoModelBPLCMS3DCorrFctnKt::operator=(const MpdFemtoModelBPLCMS3DCorrFctnKt& copy) {
  // Assignment operator
  if (this != &copy) {
    mManager = copy.mManager;
    mIsUseDenominator = copy.mIsUseDenominator;
    
    mNumerator = new TH3F( *copy.mNumerator );
    mNumeratorWeighted = new TH3F( *copy.mNumeratorWeighted );
    mNumeratorQinvWeighted = new TH3F( *copy.mNumeratorQinvWeighted );
    if ( mIsUseDenominator ) {
      mDenominator = new TH3F( *copy.mDenominator );
      mDenominatorWeighted = new TH3F( *copy.mDenominatorWeighted );
      mDenominatorQinvWeighted = new TH3F( *copy.mDenominatorQinvWeighted );
    }
  } // if (this != &copy)

  return *this;
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::connectToManager(MpdFemtoModelManager *manager) {
  mManager = manager;
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::eventBegin(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::writeOutHistos() {
  // Write out all histograms to file

  if ( mNumerator ) mNumerator->Write();
  if ( mNumeratorWeighted ) mNumeratorWeighted->Write();
  if ( mNumeratorQinvWeighted ) mNumeratorQinvWeighted->Write();
  if ( mIsUseDenominator ) {
    if ( mDenominator ) mDenominator->Write();
    if ( mDenominatorWeighted ) mDenominatorWeighted->Write();
    if ( mDenominatorQinvWeighted ) mDenominatorQinvWeighted->Write();
  } // if ( mIsUseDenominator )
}

//______________________________
TList* MpdFemtoModelBPLCMS3DCorrFctnKt::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  tOutputList->Add(mNumerator);
  tOutputList->Add(mNumeratorWeighted);
  tOutputList->Add(mNumeratorQinvWeighted);
  if ( mIsUseDenominator ) {
    tOutputList->Add(mDenominator);
    tOutputList->Add(mDenominatorWeighted);
    tOutputList->Add(mDenominatorQinvWeighted);
  }
  return tOutputList;
}

//_________________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::finish() {
  /* empty */
}

//____________________________
MpdFemtoString MpdFemtoModelBPLCMS3DCorrFctnKt::report() {
  // Construct the report
  TString report("LCMS Frame Bertsch-Pratt 3D Correlation Function Report:\n");
  report += TString::Format("Number of entries in numerator:\t%E\n", mNumerator->GetEntries());
  if ( mIsUseDenominator ) {
    report += TString::Format("Number of entries in denominator:\t%E\n", mDenominator->GetEntries());
  }
  
  if (mPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n";
    report += mPairCut->report();
  }
  else {
    report += "No PairCut specific to this CorrFctn\n";
  }
  
  return MpdFemtoString((const char *)report);
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::addRealPair(MpdFemtoPair* pair) {
  // Perform operations on real pairs

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  mNumerator->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), 1.);
  mNumeratorWeighted->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), mManager->weight(pair));
  mNumeratorQinvWeighted->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), pair->qInv());
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::addMixedPair(MpdFemtoPair* pair) {

  if ( !mIsUseDenominator ) return;

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  mDenominator->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), 1.);
  mDenominatorWeighted->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), mManager->weight(pair));
  mDenominatorQinvWeighted->Fill(pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), pair->qInv());
}
