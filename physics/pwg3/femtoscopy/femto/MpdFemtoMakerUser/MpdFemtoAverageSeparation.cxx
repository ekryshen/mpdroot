//
// One-dimensional distr. of Average Separation
//

// MpdFemtoMaker headers
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoAverageSeparation.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

ClassImp(MpdFemtoAverageSeparation)

//_________________
MpdFemtoAverageSeparation::MpdFemtoAverageSeparation(const char* title,
								     const int& nBinsAvSep, const double& AvSepLo,
								     const double& AvSepHi) :
MpdFemtoBaseCorrFctn() {
  // Parametrized constructor
  // Set initial parameters	
  setHistoParameters(nBinsAvSep, AvSepLo, AvSepHi);
 

  // Define string parameters
  TString baseName(title);
  TString baseTitle(title);
  TString appendix = "; AvSep#Delta#phi^{*}_{min} (rad);#Delta#eta";

    // Set numerator
    baseName = title;
    baseName += "_num_";
    baseTitle += "_num_";

    TH1F *hNumer = new TH1F(baseName.Data(), baseTitle.Data(), nBinsAvSep,
			    AvSepLo, AvSepHi);

    hNumer->Sumw2();
    mNumerator.push_back(hNumer);

    // Denominator
    baseName = title;
    baseName += "_den_";
    baseTitle += "_den_";


    TH1F *hDenom = new TH1F(baseName.Data(), baseTitle.Data(),nBinsAvSep,
			    AvSepLo, AvSepHi);


    hDenom->Sumw2();
    mDenominator.push_back(hDenom);
  
}


//____________________________
MpdFemtoAverageSeparation::~MpdFemtoAverageSeparation() {
  // Destructor
  if (!mNumerator.empty()) {
    mNumerator.clear();
  }
  if (!mDenominator.empty()) {
    mDenominator.clear();
  }
}

//_________________
void MpdFemtoAverageSeparation::eventBegin(const MpdFemtoEvent* event) {
  // setBField(event->bField());
}

//_________________
void MpdFemtoAverageSeparation::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________________
void MpdFemtoAverageSeparation::writeOutHistos() {
  // Write out all histograms to file
  for (unsigned int i = 0; i < mNumerator.size(); i++) {
    mNumerator.at(i)->Write();
  }
  for (unsigned int i = 0; i < mDenominator.size(); i++) {
    mDenominator.at(i)->Write();
  }
}

//______________________________
TList* MpdFemtoAverageSeparation::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  //  for (int i = 0; i < mNKtBins; i++) {
    tOutputList->Add(mNumerator.at(0));
    tOutputList->Add(mDenominator.at(0));
    //  } // for (int i=0; i<mNKtBins; i++)

  return tOutputList;
}

//_________________________
void MpdFemtoAverageSeparation::finish() {
  /* empty */
}

//____________________________
MpdFemtoString MpdFemtoAverageSeparation::report() {
  // Construct the report
  TString report("Average Separation report:\n");
  // for (int iKt = 0; iKt < mNKtBins; iKt++) {
    int iKt =0;
    report += TString::Format("Number of entries in %d-th numerator                :\t%E\n",
			      iKt, mNumerator.at(iKt)->GetEntries());
    report += TString::Format("Number of entries in %d-th denominator              :\t%E\n",
			      iKt, mDenominator.at(iKt)->GetEntries());
    // } // for (int iKt=0; iKt<mNKtBins; iKt++)

  if (mPairCut) {
    report += "Report from the front-loaded PairCut specific to this CorrFctn\n";
    report += mPairCut->report();
  } else {
    report += "No front-loaded PairCut specific to this CorrFctn\n";
  }

  return MpdFemtoString((const char *) report);
}

//-------------------------------------------

void MpdFemtoAverageSeparation::addRealPair(MpdFemtoPair* pair) {
  // Perform operations on real pairs

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  if (pair->qInv() < 0.15) {
    double avsep = pair->nominalTpcAverageSeparation(); 

    // Check that histrogram exists
    if (mNumerator.at(0)) {
      mNumerator.at(0)->Fill(avsep, 1.);
    } // if ( mNumerator.at(mIndexKt) )
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

//_________________
void MpdFemtoAverageSeparation::addMixedPair(MpdFemtoPair* pair) {
  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

 
  if (pair->qInv() < 0.15) {
     double avsep = pair->nominalTpcAverageSeparation(); 

    // Check that histrogram exists
    if (mDenominator.at(0)) {
      mDenominator.at(0)->Fill(avsep, 1.);
    } // if ( mDenominator.at(mIndexKt) )
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

void MpdFemtoAverageSeparation::setHistoParameters(const int& nBinsAvSep, const double& AvSepLo, const double& AvSepHi) {

  
  if (AvSepLo > AvSepHi) {
    std::cout << "[WARNING] void MpdFemtoAverageSeparation::void setHistoParameters - "
	      << " AvSepLo > AvSepHi. Resetting to defaults." << std::endl;
    mAvSepRange[0] = 0;
    mAvSepRange[1] = 100.;
  } 
  else {
    mAvSepRange[0] = AvSepLo;
    mAvSepRange[1] = AvSepHi;
   }

  if (nBinsAvSep < 1) {
    std::cout << "[WARNING] void MpdFemtoAverageSeparation::void setHistoParameters - "
	      << " nBinsAvSep<1. Resetting to defaults." << std::endl;
    mAvSepBins = 200;
     } else {
      mAvSepBins = nBinsAvSep;
     }

}
