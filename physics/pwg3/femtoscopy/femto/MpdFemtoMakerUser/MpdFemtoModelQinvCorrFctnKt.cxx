//
// One-dimensional correlation function with kT-binning
//

// MpdFemtoMaker headers
#include "MpdFemtoModelGausLCMSFreezeOutGenerator.h"
#include "MpdFemtoModelHiddenInfo.h"
#include "MpdFemtoModelManager.h"
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoModelQinvCorrFctnKt.h"

// ROOT headers
#include "TString.h"
#include "TH1D.h"

// C++ headers
#include <iostream>

ClassImp(MpdFemtoModelQinvCorrFctnKt)

//_________________
MpdFemtoModelQinvCorrFctnKt::MpdFemtoModelQinvCorrFctnKt(const char* title,
							 const int &nQbins, const double &qLo, const double &qHi,
							 const int& nKtBins, const double& ktLo, const double& ktHi,
							 const bool useDenominator) :
MpdFemtoBaseCorrFctn(),
mManager(nullptr) {
  // Constructor

  // Set general parameters
  setKtRange(nKtBins, ktLo, ktHi);
  setUseDenominator( useDenominator );

  TString baseName(title);
  TString baseTitle(title);
  TString appendix = ";q_{inv} (GeV/c);Entries";

  // Histogram loop
  for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {

    // Numerator
    baseName = title;
    baseName += "_num_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH1D *histoNum = new TH1D(baseName.Data(), baseTitle.Data(),
			      nQbins, qLo, qHi);
    histoNum->Sumw2();
    mNumerator.push_back(histoNum);

    // Numerator weighted
    baseName = title;
    baseName += "_num_wei_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH1D *histoNumW = new TH1D(baseName.Data(), baseTitle.Data(),
			       nQbins, qLo, qHi);
    histoNumW->Sumw2();
    mNumeratorWeighted.push_back(histoNumW);

    if ( mIsUseDenominator ) {
      
      // Denominator
      baseName = title;
      baseName += "_den_";
      baseName += iKtBin;
      baseTitle = baseName;
      baseTitle += " ";
      baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
      baseTitle += "#leq k_{T} (GeV/c) #leq";
      baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
      baseTitle += appendix;
      TH1D *histoDen = new TH1D(baseName.Data(), baseTitle.Data(),
				nQbins, qLo, qHi);
      histoDen->Sumw2();
      mDenominator.push_back(histoDen);

      // Denominator weighted
      baseName = title;
      baseName += "_den_wei_";
      baseName += iKtBin;
      baseTitle = baseName;
      baseTitle += " ";
      baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
      baseTitle += "#leq k_{T} (GeV/c) #leq";
      baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
      baseTitle += appendix;
      TH1D *histoDenW = new TH1D(baseName.Data(), baseTitle.Data(),
				 nQbins, qLo, qHi);
      histoDenW->Sumw2();
      mDenominatorWeighted.push_back(histoDen);
      
    } // if ( mIsUseDenominator )
  } // for(int i=0; i<mNumberKt; i++)
}

//_________________
MpdFemtoModelQinvCorrFctnKt::MpdFemtoModelQinvCorrFctnKt(const MpdFemtoModelQinvCorrFctnKt& c) :
MpdFemtoBaseCorrFctn(c),
mManager(c.mManager),
mNKtBins(c.mNKtBins),
mKtStep(c.mKtStep),
mIsUseDenominator(c.mIsUseDenominator) {
  // Copy constructor
  mKtRange[0] = c.mKtRange[0];
  mKtRange[1] = c.mKtRange[1];

  for (int iCF = 0; iCF < mNKtBins; iCF++) {
    mNumerator.push_back( (c.mNumerator.at(iCF)) ?
			  new TH1D(*c.mNumerator.at(iCF)) : nullptr);
    mNumeratorWeighted.push_back( (c.mNumeratorWeighted.at(iCF)) ?
				  new TH1D(*c.mNumeratorWeighted.at(iCF)) : nullptr);
    if ( mIsUseDenominator ) {
      mDenominator.push_back( (c.mDenominator.at(iCF)) ?
			      new TH1D(*c.mDenominator.at(iCF)) : nullptr);
      mDenominatorWeighted.push_back( (c.mDenominatorWeighted.at(iCF)) ?
				      new TH1D(*c.mDenominatorWeighted.at(iCF)) : nullptr);
    }
  } // for ( int iBin=0; iBin<mNumberKt; iBin++ )
}

//_________________
MpdFemtoModelQinvCorrFctnKt& MpdFemtoModelQinvCorrFctnKt::operator=(const MpdFemtoModelQinvCorrFctnKt& c) {
  // Assignment operator
  if (this != &c) {
    mManager = c.mManager;
    mNKtBins = c.mNKtBins;
    mKtRange[0] = c.mKtRange[0];
    mKtRange[1] = c.mKtRange[1];
    mKtStep = c.mKtStep;
    mIsUseDenominator = c.mIsUseDenominator;
    // Loop over kT bins
    for (int iCF = 0; iCF < mNKtBins; iCF++) {
      mNumerator.push_back( (c.mNumerator.at(iCF)) ?
			    new TH1D(*c.mNumerator.at(iCF)) : nullptr);
      mNumeratorWeighted.push_back( (c.mNumeratorWeighted.at(iCF)) ?
				    new TH1D(*c.mNumeratorWeighted.at(iCF)) : nullptr);
      if ( mIsUseDenominator ) {
	mDenominator.push_back( (c.mDenominator.at(iCF)) ?
				new TH1D(*c.mDenominator.at(iCF)) : nullptr);
	mDenominatorWeighted.push_back( (c.mDenominatorWeighted.at(iCF)) ?
					new TH1D(*c.mDenominatorWeighted.at(iCF)) : nullptr);
      }
    } // for ( int iBin=0; iBin<mNumberKt; iBin++ )
  } // if (this != &c)

  return *this;
}

//_________________
MpdFemtoModelQinvCorrFctnKt::~MpdFemtoModelQinvCorrFctnKt() {
  // Destructor
  if (!mNumerator.empty()) {
    mNumerator.clear();
  }
  if (!mNumeratorWeighted.empty()) {
    mNumeratorWeighted.clear();
  }
  if (!mDenominator.empty()) {
    mDenominator.clear();
  }
  if (!mDenominatorWeighted.empty()) {
    mDenominatorWeighted.clear();
  }
  delete mManager;
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::writeOutHistos() {
  // Write histograms to the file
  if (!mNumerator.empty()) {
    for (unsigned int i = 0; i < mNumerator.size(); i++) {
      mNumerator.at(i)->Write();
    }
  }
  if (!mNumeratorWeighted.empty()) {
    for (unsigned int i = 0; i < mNumeratorWeighted.size(); i++) {
      mNumeratorWeighted.at(i)->Write();
    }
  }
  if ( mIsUseDenominator ) {
    if (!mDenominator.empty()) {
      for (unsigned int i = 0; i < mDenominator.size(); i++) {
	mDenominator.at(i)->Write();
      }
    }
    if (!mDenominatorWeighted.empty()) {
      for (unsigned int i = 0; i < mDenominatorWeighted.size(); i++) {
	mDenominatorWeighted.at(i)->Write();
      }
    }
  }
}

//_________________
TList* MpdFemtoModelQinvCorrFctnKt::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *outputList = new TList();
  for (int i = 0; i < mNKtBins; i++) {
    outputList->Add(mNumerator.at(i));
    outputList->Add(mNumeratorWeighted.at(i));
    if ( mIsUseDenominator ) {
      outputList->Add(mDenominator.at(i));
      outputList->Add(mDenominatorWeighted.at(i));
    }
  }
  return outputList;
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::eventBegin(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::finish() {
  /* empty */
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::setKtRange(const int& nbins, const double& ktLo, const double& ktHi) {
  mNKtBins = nbins;
  mKtRange[0] = ktLo;
  mKtRange[1] = ktHi;
  mKtStep = (mKtRange[1] - mKtRange[0]) / mNKtBins;
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::addRealPair(MpdFemtoPair* pair) {

  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  int mIndexKt = (int) ( (pair->kT() - mKtRange[0]) / mKtStep );

  if ((mIndexKt >= 0) && (mIndexKt < mNKtBins)) {
    if ( mNumerator.at(mIndexKt) ) {
      mNumerator.at(mIndexKt)->Fill( fabs(pair->qInv()), 1.);
    } else {
      std::cout << "[ERROR] void MpdFemtoModelQinvCorrFctnKt::addRealPair(MpdFemtoPair* pair) - "
		<< "No correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
    if ( mNumeratorWeighted.at(mIndexKt) ) {
      double weight = mManager->weight(pair);
      
      mNumeratorWeighted.at(mIndexKt)->Fill( fabs(pair->qInv()), weight);
    } else {
      std::cout << "[ERROR] void MpdFemtoModelQinvCorrFctnKt::addRealPair(MpdFemtoPair* pair) - "
		<< "No weighted correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
  } // if ( ( mIndexKt>=0 ) && ( mIndexKt<mNumberKt ) )
}

//_________________
void MpdFemtoModelQinvCorrFctnKt::addMixedPair(MpdFemtoPair* pair) {

  if ( !mIsUseDenominator ) return;
  
  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  int mIndexKt = (int) ( (fabs(pair->kT()) - mKtRange[0]) / mKtStep );

  if ( (mIndexKt >= 0) && (mIndexKt < mNKtBins) ) {
    if (mDenominator.at(mIndexKt)) {
      mDenominator.at(mIndexKt)->Fill( fabs(pair->qInv()) );
    } else {
      std::cout << "[ERROR] void MpdFemtoModelQinvCorrFctnKt::addMixedPair(MpdFemtoPair* pair) - "
		<< "No correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
    if (mDenominatorWeighted.at(mIndexKt)) {
      double weight = mManager->weight(pair);
      mDenominatorWeighted.at(mIndexKt)->Fill( fabs(pair->qInv()), weight);
    } else {
      std::cout << "[ERROR] void MpdFemtoModelQinvCorrFctnKt::addMixedPair(MpdFemtoPair* pair) - "
		<< "No weighted correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
  } // if ( ( mIndexKt>=0 ) && ( mIndexKt<mNumberKt ) )
}

//_________________
MpdFemtoString MpdFemtoModelQinvCorrFctnKt::report() {

  // Make a report
  MpdFemtoString tStr = "MpdFemtoModelQinvCorrFctnKt report";

  int mNumeratorEntries = 0;
  int mNumeratorWEntries = 0;
  int mDenominatorEntries = 0;
  int mDenominatorWEntries = 0;

  for (int i = 0; i < mNKtBins; i++) {
    mNumeratorEntries += (int) mNumerator.at(i)->GetEntries();
    mNumeratorWEntries += (int) mNumeratorWeighted.at(i)->GetEntries();
    if ( mIsUseDenominator ) {
      mDenominatorEntries += (int) mDenominator.at(i)->GetEntries();
      mDenominatorWEntries += (int) mDenominatorWeighted.at(i)->GetEntries();
    }
  } // for (int i=0; i<mNumberKt; i++)

  tStr += TString::Format("\nTotal number of pairs in the numerator         :\t%d\n", mNumeratorEntries);
  tStr += TString::Format("\nTotal number of pairs in the weighted numerator:\t%d\n", mNumeratorEntries);
  if ( mIsUseDenominator ) {
    tStr += TString::Format("Total number of pairs in the denominator         :\t%d\n", mDenominatorEntries);
    tStr += TString::Format("Total number of pairs in the weighted denominator:\t%d\n", mDenominatorEntries);
  }
  for (int i = 0; i < mNKtBins; i++) {
    tStr += TString::Format("Total number of pairs in %d-th numerator         :\t%E\n",
			    i, mNumerator.at(i)->GetEntries());
    tStr += TString::Format("Total number of pairs in %d-th weighted numerator:\t%E\n",
			    i, mNumeratorWeighted.at(i)->GetEntries());
    if ( mIsUseDenominator ) {
      tStr += TString::Format("Total number of pairs in %d-th denominator         :\t%E\n",
			      i, mDenominator.at(i)->GetEntries());
      tStr += TString::Format("Total number of pairs in %d-th weighted denominator:\t%E\n",
			      i, mDenominatorWeighted.at(i)->GetEntries());
    }
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
void MpdFemtoModelQinvCorrFctnKt::connectToManager(MpdFemtoModelManager *manager) {
  mManager = manager;
}
