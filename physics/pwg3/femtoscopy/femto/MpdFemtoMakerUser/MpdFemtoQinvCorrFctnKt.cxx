//
// One-dimensional correlation function with kT-binning
//

// MpdFemtoMakerUser headers
#include "MpdFemtoQinvCorrFctnKt.h"

// ROOT headers
#include "TString.h"
#include "TH1D.h"

// C++ headers
#include <iostream>

ClassImp(MpdFemtoQinvCorrFctnKt)

//_________________
MpdFemtoQinvCorrFctnKt::MpdFemtoQinvCorrFctnKt(const char* title,
					       const int &nQbins, const double &qLo, const double &qHi,
					       const int& nKtBins, const double& ktLo, const double& ktHi) :
MpdFemtoBaseCorrFctn() {
  // Constructor

  // Define number of kT bins and kT step
  setKtRange(nKtBins, ktLo, ktHi);

  TString baseName(title);
  TString baseTitle(title);
  TString appendix = ";q_{inv} (GeV/c);Entries";

  // Histogram loop
  for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {

    baseName = title;
    baseName += "_num_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH1D *histoNum = new TH1D(baseName.Data(), baseTitle.Data(),
			      nQbins, qLo, qHi);
    histoNum->Sumw2();
    mNumerator.push_back(histoNum);


    baseName = title;
    baseName += "_den_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH1D *histoDen = new TH1D(baseName.Data(), baseTitle.Data(),
			      nQbins, qLo, qHi);
    histoDen->Sumw2();
    mDenominator.push_back(histoDen);
  } // for(int i=0; i<mNumberKt; i++)
}

//_________________
MpdFemtoQinvCorrFctnKt::MpdFemtoQinvCorrFctnKt(const MpdFemtoQinvCorrFctnKt& c) :
MpdFemtoBaseCorrFctn(c),
mNKtBins(c.mNKtBins),
mKtStep(c.mKtStep) {
  // Copy constructor

  mKtRange[0] = c.mKtRange[0];
  mKtRange[1] = c.mKtRange[1];

  for (int iCF = 0; iCF < mNKtBins; iCF++) {
    mNumerator.push_back((c.mNumerator.at(iCF)) ? new TH1D(*c.mNumerator.at(iCF)) : nullptr);
    mDenominator.push_back((c.mDenominator.at(iCF)) ? new TH1D(*c.mDenominator.at(iCF)) : nullptr);
  } // for ( int iBin=0; iBin<mNumberKt; iBin++ )
}

//_________________
MpdFemtoQinvCorrFctnKt& MpdFemtoQinvCorrFctnKt::operator=(const MpdFemtoQinvCorrFctnKt& c) {
  // Assignment operator
  if (this != &c) {
    mNKtBins = c.mNKtBins;
    mKtRange[0] = c.mKtRange[0];
    mKtRange[1] = c.mKtRange[1];
    mKtStep = c.mKtStep;
    for (int iCF = 0; iCF < mNKtBins; iCF++) {
      mNumerator.push_back((c.mNumerator.at(iCF)) ? new TH1D(*c.mNumerator.at(iCF)) : nullptr);
      mDenominator.push_back((c.mDenominator.at(iCF)) ? new TH1D(*c.mDenominator.at(iCF)) : nullptr);
    } // for ( int iBin=0; iBin<mNumberKt; iBin++ )
  } // if (this != &c)

  return *this;
}

//_________________
MpdFemtoQinvCorrFctnKt::~MpdFemtoQinvCorrFctnKt() {
  // Destructor
  if (!mNumerator.empty()) {
    mNumerator.clear();
  }
  if (!mDenominator.empty()) {
    mDenominator.clear();
  }
}

//_________________
void MpdFemtoQinvCorrFctnKt::writeOutHistos() {
  // Write histograms to the file
  if (!mNumerator.empty()) {
    for (unsigned int i = 0; i < mNumerator.size(); i++) {
      mNumerator.at(i)->Write();
    }
  }
  if (!mDenominator.empty()) {
    for (unsigned int i = 0; i < mDenominator.size(); i++) {
      mDenominator.at(i)->Write();
    }
  }
}

//_________________
TList* MpdFemtoQinvCorrFctnKt::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *outputList = new TList();
  for (int i = 0; i < mNKtBins; i++) {
    outputList->Add(mNumerator.at(i));
    outputList->Add(mDenominator.at(i));
  }
  return outputList;
}

//_________________
void MpdFemtoQinvCorrFctnKt::eventBegin(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoQinvCorrFctnKt::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________
void MpdFemtoQinvCorrFctnKt::finish() {
  /* empty */
}

//_________________
void MpdFemtoQinvCorrFctnKt::setKtRange(const int& nbins, const double& ktLo, const double& ktHi) {
  mNKtBins = nbins;
  mKtRange[0] = ktLo;
  mKtRange[1] = ktHi;
  mKtStep = (mKtRange[1] - mKtRange[0]) / mNKtBins;
}

//_________________
void MpdFemtoQinvCorrFctnKt::addRealPair(MpdFemtoPair* pair) {

  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  int mIndexKt = (int) (((pair->kT() - mKtRange[0]) / mKtStep));

  if ((mIndexKt >= 0) && (mIndexKt < mNKtBins)) {
    if (mNumerator.at(mIndexKt)) {
      mNumerator.at(mIndexKt)->Fill(fabs(pair->qInv()));
    } else {
      std::cout << "[ERROR] void MpdFemtoQinvCorrFctnKt::addRealPair(MpdFemtoPair* pair) - "
		<< "No correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
  } // if ( ( mIndexKt>=0 ) && ( mIndexKt<mNumberKt ) )
}

//_________________
void MpdFemtoQinvCorrFctnKt::addMixedPair(MpdFemtoPair* pair) {

  // Check if pair passes front-loaded cut if exists
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  int mIndexKt = (int) (((fabs(pair->kT()) - mKtRange[0]) / mKtStep));

  if ((mIndexKt >= 0) && (mIndexKt < mNKtBins)) {
    if (mDenominator.at(mIndexKt)) {
      mDenominator.at(mIndexKt)->Fill(fabs(pair->qInv()));
    } else {
      std::cout << "[ERROR] void MpdFemtoQinvCorrFctnKt::addMixedPair(MpdFemtoPair* pair) - "
		<< "No correlation function with index: " << mIndexKt << " was found"
		<< std::endl;
    }
  } // if ( ( mIndexKt>=0 ) && ( mIndexKt<mNumberKt ) )
}

//_________________
MpdFemtoString MpdFemtoQinvCorrFctnKt::report() {

  // Make a report
  MpdFemtoString tStr = "MpdFemtoQinvCorrFctnKt report";

  int mNumeratorEntries = 0;
  int mDenominatorEntries = 0;

  for (int i = 0; i < mNKtBins; i++) {
    mNumeratorEntries += (int) mNumerator.at(i)->GetEntries();
    mDenominatorEntries += (int) mDenominator.at(i)->GetEntries();
  } // for (int i=0; i<mNumberKt; i++)

  tStr += TString::Format("\nTotal number of pairs in the numerator  :\t%d\n", mNumeratorEntries);
  tStr += TString::Format("Total number of pairs in the denominator:\t%d\n", mDenominatorEntries);
  for (int i = 0; i < mNKtBins; i++) {
    tStr += TString::Format("Total number of pairs in %d-th numerator :\t%E\n",
			    i, mNumerator.at(i)->GetEntries());
    tStr += TString::Format("Total number of pairs in %d-th denominator  :\t%E\n",
			    i, mDenominator.at(i)->GetEntries());
  }

  if (mPairCut) {
    tStr += "Here is the PairCut specific to this CorrFctn\n";
    tStr += mPairCut->report();
  } else {
    tStr += "No PairCut specific to this CorrFctn\n";
  }

  return tStr;
}
