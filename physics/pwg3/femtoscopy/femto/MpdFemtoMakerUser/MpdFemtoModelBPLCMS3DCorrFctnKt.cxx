//
// Three-dimensional Bertsch-Pratt correlation function in LCMS for the model estimations
//

// MpdFemtoMaker headers
#include "MpdFemtoModelGausLCMSFreezeOutGenerator.h"
#include "MpdFemtoModelHiddenInfo.h"
#include "MpdFemtoModelManager.h"
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoModelBPLCMS3DCorrFctnKt.h"

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
								 const int& ktBins, const double& ktLo,
								 const double& ktHi, const bool isUseDenominator) :
MpdFemtoBaseCorrFctn(),
mManager(nullptr) {
  // Parametrized constructor

  // Set general parameters
  setKtRange( ktBins, ktLo, ktHi );
  setUseDenominator( isUseDenominator );

  // Define string parameters
  TString baseName(title);
  TString baseTitle(title);
  TString appendix = ";q_{out} (GeV/c);q_{side} (GeV/c);q_{long} (GeV/c)";

  // Loop over kTBins
  for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {

    // Set numerator
    baseName = title;
    baseName += "_num_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
    baseTitle += "#leq k_{T} (GeV/c) #leq";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
    baseTitle += appendix;
    TH3F *histoNum = new TH3F(baseName.Data(), baseTitle.Data(),
			      nBins, qLo, qHi,
			      nBins, qLo, qHi,
			      nBins, qLo, qHi);
    histoNum->Sumw2();
    mNumerator.push_back( histoNum );

    // Weighted numerator
    baseName = title;
    baseName += "_num_wei_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
    baseTitle += "#leq kT (GeV/c) #leq";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
    baseTitle += appendix;
    TH3F *histoNumW = new TH3F(baseName.Data(), baseTitle.Data(),
			       nBins, qLo, qHi,
			       nBins, qLo, qHi,
			       nBins, qLo, qHi);
    histoNumW->Sumw2();
    mNumeratorWeighted.push_back( histoNumW );

    // Qinv weighted numerator
    baseName = title;
    baseName += "_num_qinv_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
    baseTitle += "#leq kT (GeV/c) #leq";
    baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
    baseTitle += appendix;
    TH3F *histoNumQinvW = new TH3F(baseName.Data(), baseTitle.Data(),
				   nBins, qLo, qHi,
				   nBins, qLo, qHi,
				   nBins, qLo, qHi);
    histoNumQinvW->Sumw2();
    mNumeratorQinvWeighted.push_back( histoNumQinvW );

    if ( mIsUseDenominator ) {

      // Denominator
      baseName = title;
      baseName += "_den_";
      baseName += iKtBin;
      baseTitle = baseName;
      baseTitle += " ";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
      baseTitle += "#leq kT (GeV/c) #leq";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
      baseTitle += appendix;
      TH3F *histoDen = new TH3F(baseName.Data(), baseTitle.Data(),
				nBins, qLo, qHi,
				nBins, qLo, qHi,
				nBins, qLo, qHi);
      histoDen->Sumw2();
      mDenominator.push_back( histoDen );
      
      // Weighted denominator
      baseName = title;
      baseName += "_den_wei_";
      baseName += iKtBin;
      baseTitle = baseName;
      baseTitle += " ";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
      baseTitle += "#leq kT (GeV/c) #leq";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
      baseTitle += appendix;
      TH3F *histoDenW = new TH3F(baseName.Data(), baseTitle.Data(),
				 nBins, qLo, qHi,
				 nBins, qLo, qHi,
				 nBins, qLo, qHi);
      histoDenW->Sumw2();
      mDenominatorWeighted.push_back( histoDenW );

      // Qinv weighted denominator
      baseName = title;
      baseName += "_den_qinv_";
      baseName += iKtBin;
      baseTitle = baseName;
      baseTitle += " ";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * iKtBin) );
      baseTitle += "#leq kT (GeV/c) #leq";
      baseTitle += Form("%3.2f",( mKtRange[0] + mKtStep * (iKtBin+1) ) );
      baseTitle += appendix;
      TH3F *histoDenQinvW = new TH3F(baseName.Data(), baseTitle.Data(),
				     nBins, qLo, qHi,
				     nBins, qLo, qHi,
				     nBins, qLo, qHi);
      histoDenQinvW->Sumw2();
      mDenominatorQinvWeighted.push_back( histoDenQinvW );
    } // if ( mIsUseDenominator )
    
  } // for (int iKtBin=0; iKtBin<mNKtBins; iKtBin++ )
}

//_________________
MpdFemtoModelBPLCMS3DCorrFctnKt::MpdFemtoModelBPLCMS3DCorrFctnKt(const MpdFemtoModelBPLCMS3DCorrFctnKt& copy) :
MpdFemtoBaseCorrFctn(copy),
mManager(copy.mManager),
mNKtBins(copy.mNKtBins),
mKtStep(copy.mKtStep),
mIsUseDenominator(copy.mIsUseDenominator) {
  // Copy constructor
  for (int i = 0; i < 2; i++) {
    mKtRange[i] = copy.mKtRange[i];
  }

  // Loop over kT bins and make histogram copies
  for ( int iCF=0; iCF<mNKtBins; iCF++ ) {
    mNumerator.push_back( ( copy.mNumerator.at(iCF) ) ?
			  new TH3F( *copy.mNumerator.at(iCF) ) : nullptr );
    mNumeratorWeighted.push_back( ( copy.mNumeratorWeighted.at(iCF) ) ?
				  new TH3F( *copy.mNumeratorWeighted.at(iCF) ) : nullptr );
    mNumeratorQinvWeighted.push_back( ( copy.mNumeratorQinvWeighted.at(iCF) ) ?
				  new TH3F( *copy.mNumeratorQinvWeighted.at(iCF) ) : nullptr );
    
    if ( mIsUseDenominator ) {
      mDenominator.push_back( ( copy.mDenominator.at(iCF) ) ?
			      new TH3F( *copy.mDenominator.at(iCF) ) : nullptr );
      mDenominatorWeighted.push_back( ( copy.mDenominatorWeighted.at(iCF) ) ?
				      new TH3F( *copy.mDenominatorWeighted.at(iCF) ) : nullptr );
      mDenominatorQinvWeighted.push_back( ( copy.mDenominatorQinvWeighted.at(iCF) ) ?
				      new TH3F( *copy.mDenominatorQinvWeighted.at(iCF) ) : nullptr );
    } // if ( mIsUseDenominator )
  } // for ( int iCF=0; iCF<mNKtBins; iCF++ )
}

//____________________________
MpdFemtoModelBPLCMS3DCorrFctnKt::~MpdFemtoModelBPLCMS3DCorrFctnKt() {
  // Destructor
  if ( !mNumerator.empty() ) {
    mNumerator.clear();
  }
  if ( !mNumeratorWeighted.empty() ) {
    mNumeratorWeighted.clear();
  }
  if ( !mNumeratorQinvWeighted.empty() ) {
    mNumeratorQinvWeighted.clear();
  }
  if ( !mDenominator.empty() ) {
    mDenominator.clear();
  }
  if ( !mDenominatorWeighted.empty() ) {
    mDenominatorWeighted.clear();
  }
  if ( !mDenominatorQinvWeighted.empty() ) {
    mDenominatorQinvWeighted.clear();
  }
}

//_________________________
MpdFemtoModelBPLCMS3DCorrFctnKt& MpdFemtoModelBPLCMS3DCorrFctnKt::operator=(const MpdFemtoModelBPLCMS3DCorrFctnKt& copy) {
  // Assignment operator
  if (this != &copy) {
    mManager = copy.mManager;
    mNKtBins = copy.mNKtBins;
    mKtStep = copy.mKtStep;
    mKtRange[0] = copy.mKtRange[0];
    mKtRange[1] = copy.mKtRange[1];
    mIsUseDenominator = copy.mIsUseDenominator;

    // Loop over kT bins and make histogram copies
    for ( int iCF=0; iCF<mNKtBins; iCF++ ) {
      mNumerator.push_back( ( copy.mNumerator.at(iCF) ) ?
			    new TH3F( *copy.mNumerator.at(iCF) ) : nullptr );
      mNumeratorWeighted.push_back( ( copy.mNumeratorWeighted.at(iCF) ) ?
				    new TH3F( *copy.mNumeratorWeighted.at(iCF) ) : nullptr );
      mNumeratorQinvWeighted.push_back( ( copy.mNumeratorQinvWeighted.at(iCF) ) ?
					new TH3F( *copy.mNumeratorQinvWeighted.at(iCF) ) : nullptr );
    
      if ( mIsUseDenominator ) {
	mDenominator.push_back( ( copy.mDenominator.at(iCF) ) ?
				new TH3F( *copy.mDenominator.at(iCF) ) : nullptr );
	mDenominatorWeighted.push_back( ( copy.mDenominatorWeighted.at(iCF) ) ?
					new TH3F( *copy.mDenominatorWeighted.at(iCF) ) : nullptr );
	mDenominatorQinvWeighted.push_back( ( copy.mDenominatorQinvWeighted.at(iCF) ) ?
					    new TH3F( *copy.mDenominatorQinvWeighted.at(iCF) ) : nullptr );
      } // if ( mIsUseDenominator )
    } // for ( int iCF=0; iCF<mNKtBins; iCF++ )
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
  if ( !mNumerator.empty() ) {
    for (unsigned int i=0; i<mNumerator.size(); i++ ) {
      mNumerator.at(i)->Write();
    }
  }
  if ( !mNumeratorWeighted.empty() ) {
    for (unsigned int i=0; i<mNumeratorWeighted.size(); i++ ) {
      mNumeratorWeighted.at(i)->Write();
    }
  }
  if ( !mNumeratorQinvWeighted.empty() ) {
    for (unsigned int i=0; i<mNumeratorQinvWeighted.size(); i++ ) {
      mNumeratorQinvWeighted.at(i)->Write();
    }
  }
  if ( !mDenominator.empty() ) {
    for (unsigned int i=0; i<mDenominator.size(); i++ ) {
      mDenominator.at(i)->Write();
    }
  }
  if ( !mDenominatorWeighted.empty() ) {
    for (unsigned int i=0; i<mDenominatorWeighted.size(); i++ ) {
      mDenominatorWeighted.at(i)->Write();
    }
  }
  if ( !mDenominatorQinvWeighted.empty() ) {
    for (unsigned int i=0; i<mDenominatorQinvWeighted.size(); i++ ) {
      mDenominatorQinvWeighted.at(i)->Write();
    }
  }
}

//______________________________
TList* MpdFemtoModelBPLCMS3DCorrFctnKt::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  for (int i = 0; i < mNKtBins; i++) {
    tOutputList->Add(mNumerator.at(i));
    tOutputList->Add(mNumeratorWeighted.at(i));
    tOutputList->Add(mNumeratorQinvWeighted.at(i));
    if ( mIsUseDenominator ) {
      tOutputList->Add(mDenominator.at(i));
      tOutputList->Add(mDenominatorWeighted.at(i));
      tOutputList->Add(mDenominatorQinvWeighted.at(i));
    }
  } // for (int i=0; i<mNKtBins; i++)

  return tOutputList;
}

//_________________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::finish() {
  /* empty */
}

//____________________________
MpdFemtoString MpdFemtoModelBPLCMS3DCorrFctnKt::report() {
  // Construct the report
  TString report("LCMS frame Bertsch-Pratt 3D correlation functions with kT binning report with model weights:\n");
  for (int iKt = 0; iKt < mNKtBins; iKt++) {
    report += TString::Format("Number of entries in %d-th numerator                :\t%E\n",
			      iKt, mNumerator.at(iKt)->GetEntries());
    report += TString::Format("Number of entries in %d-th weighted numerator       :\t%E\n",
			      iKt, mNumeratorWeighted.at(iKt)->GetEntries());
    report += TString::Format("Number of entries in %d-th qInv weighted numerator  :\t%E\n",
			      iKt, mNumeratorWeighted.at(iKt)->GetEntries());
    if ( mIsUseDenominator ) {
      report += TString::Format("Number of entries in %d-th denominator              :\t%E\n",
				iKt, mDenominator.at(iKt)->GetEntries());
      report += TString::Format("Number of entries in %d-th weighted denominator     :\t%E\n",
				iKt, mDenominatorWeighted.at(iKt)->GetEntries());
      report += TString::Format("Number of entries in %d-th qInv weighted denominator:\t%E\n",
				iKt, mDenominatorWeighted.at(iKt)->GetEntries());
    }
  } // for (int iKt=0; iKt<mNKtBins; iKt++)

  if (mPairCut) {
    report += "Report from the front-loaded PairCut specific to this CorrFctn\n";
    report += mPairCut->report();
  } else {
    report += "No front-loaded PairCut specific to this CorrFctn\n";
  }

  return MpdFemtoString((const char *) report);
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::addRealPair(MpdFemtoPair* pair) {
  // Perform operations on real pairs

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  // Find index for the current kT value
  int mIndexKt = (int) ( ( pair->kT() - mKtRange[0] ) / mKtStep );

  // Check that index is in the requested kT range
  if ( ( mIndexKt>=0 ) && ( mIndexKt<mNKtBins ) ) {

    // Check that histrogram exists
    if ( mNumerator.at(mIndexKt) ) {
      mNumerator.at(mIndexKt)->Fill(pair->qOutCMS(),
				    pair->qSideCMS(),
				    pair->qLongCMS(),
				    1.);
    } // if ( mNumerator[mIndexKt] )

    // Check that histrogram exists
    if ( mNumeratorWeighted.at(mIndexKt) ) {
      double weight = mManager->weight(pair);
      mNumeratorWeighted.at(mIndexKt)->Fill(pair->qOutCMS(),
					    pair->qSideCMS(),
					    pair->qLongCMS(),
					    weight);
    } // if ( mNumeratorWeighted[mIndexKt] )

    // Check that histogram exists
    if ( mNumeratorQinvWeighted.at(mIndexKt) ) {
      mNumeratorQinvWeighted.at(mIndexKt)->Fill(pair->qOutCMS(),
						pair->qSideCMS(),
						pair->qLongCMS(),
						pair->qInv());
    }
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::addMixedPair(MpdFemtoPair* pair) {

  if ( !mIsUseDenominator ) return;

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  // Find index for the current kT value
  int mIndexKt = (int) ( (pair->kT() - mKtRange[0]) / mKtStep );

  // Check that index is in the requested kT range
  if ( ( mIndexKt>=0 ) && ( mIndexKt<mNKtBins ) ) {

    // Check that histrogram exists
    if ( mDenominator.at(mIndexKt) ) {
      mDenominator.at(mIndexKt)->Fill(pair->qOutCMS(),
				      pair->qSideCMS(),
				      pair->qLongCMS(),
				      1.);
    } // if ( mDenominator[mIndexKt] )

    // Check that histrogram exists
    if ( mDenominatorWeighted.at(mIndexKt) ) {
      double weight = mManager->weight(pair);
      mDenominatorWeighted.at(mIndexKt)->Fill(pair->qOutCMS(),
					      pair->qSideCMS(),
					      pair->qLongCMS(),
					      weight);
    } // if ( mDenominatorWeighted[mIndexKt] )

    // Check that histrogram exists
    if ( mDenominatorQinvWeighted.at(mIndexKt) ) {
      mDenominatorQinvWeighted.at(mIndexKt)->Fill(pair->qOutCMS(),
						  pair->qSideCMS(),
						  pair->qLongCMS(),
						  pair->qInv());
    } // if ( mDenominatorQinvWeighted[mIndexKt] )
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

//_________________
void MpdFemtoModelBPLCMS3DCorrFctnKt::setKtRange(const int& nKtBins, const float& kTLow, const float& kTHi) {
  if (nKtBins >= 1) {
    mNKtBins = nKtBins;
  } else {
    std::cout << "[WARNING} void MpdFemtoModelBPLCMS3DCorrFctnKt::setKtRange - "
	      << "Number of kT bins must be positive" << std::endl;
    mNKtBins = 1;
  }

  if (kTLow < kTHi) {
    mKtRange[0] = kTLow;
    mKtRange[1] = kTHi;
  } else {
    std::cout << "[WARNING} void MpdFemtoModelBPLCMS3DCorrFctnKt::setKtRange - "
	      << "kTLow < kTHi. Resetting to the integral" << std::endl;
    mKtRange[0] = 0.15;
    mKtRange[1] = 1.05;
  }
  mKtStep = (mKtRange[1] - mKtRange[0]) / mNKtBins;
}
