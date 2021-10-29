//
// Two-dimensional $\Delta\eta$ vs. $\Delta\phi^{*}_{min}$
//

// MpdFemtoMaker headers
#include "MpdFemtoPair.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoDeltaEtaDeltaPhiStarMinKt.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

ClassImp(MpdFemtoDeltaEtaDeltaPhiStarMinKt)

//_________________
MpdFemtoDeltaEtaDeltaPhiStarMinKt::MpdFemtoDeltaEtaDeltaPhiStarMinKt(const char* title,
								     const int& nBinsEta, const double& etaLo,
								     const double& etaHi,
								     const int& nBinsPhi, const double& phiLo,
								     const double& phiHi,
								     const int& ktBins, const double& ktLo,
								     const double& ktHi) :
MpdFemtoBaseCorrFctn() {
  // Parametrized constructor
  // Set initial parameters	
  setCharges(1, 1);
  setBField(-0.5);
  mRadiusStep = 0.35;
  mRadii[0] = 0.6;
  mRadii[1] = 2.;
  setHistoParameters(nBinsEta, etaLo, etaHi, nBinsPhi, phiLo, phiHi);
  setKtRange(ktBins, ktLo, ktHi);

  // Define string parameters
  TString baseName(title);
  TString baseTitle(title);
  TString appendix = ";#Delta#phi^{*}_{min} (rad);#Delta#eta";

  // Loop over kTBins
  for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {

    // Set numerator
    baseName = title;
    baseName += "_num_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH2F *hNumer = new TH2F(baseName.Data(), baseTitle.Data(),
			    mPhiBins, mPhiRange[0], mPhiRange[1],
			    mEtaBins, mEtaRange[0], mEtaRange[1]);
    hNumer->Sumw2();
    mNumerator.push_back(hNumer);

    // Denominator
    baseName = title;
    baseName += "_den_";
    baseName += iKtBin;
    baseTitle = baseName;
    baseTitle += " ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * iKtBin));
    baseTitle += "#leq k_{T} (GeV/c) #leq ";
    baseTitle += Form("%3.2f", (mKtRange[0] + mKtStep * (iKtBin + 1)));
    baseTitle += appendix;
    TH2F *hDenom = new TH2F(baseName.Data(), baseTitle.Data(),
			    mPhiBins, mPhiRange[0], mPhiRange[1],
			    mEtaBins, mEtaRange[0], mEtaRange[1]);
    hDenom->Sumw2();
    mDenominator.push_back(hDenom);
  } // for (int iKtBin=0; iKtBin<mNKtBins; iKtBin++ )
}

//_________________
MpdFemtoDeltaEtaDeltaPhiStarMinKt::MpdFemtoDeltaEtaDeltaPhiStarMinKt(const MpdFemtoDeltaEtaDeltaPhiStarMinKt& copy) :
MpdFemtoBaseCorrFctn(copy) {
  // Copy constructor
  setHistoParameters(copy.mEtaBins, copy.mEtaRange[0], copy.mEtaRange[1],
		     copy.mPhiBins, copy.mPhiRange[0], copy.mPhiRange[1]);
  setKtRange(copy.mNKtBins, copy.mKtRange[0], copy.mKtRange[1]);
  setCharges(copy.mQ1, copy.mQ2);
  setBField(copy.mBField);
  mRadiusStep = copy.mRadiusStep;
  mRadii[0] = copy.mRadii[0];
  mRadii[1] = copy.mRadii[1];

  // Loop over kT bins and make histogram copies
  for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {
    mNumerator.push_back(copy.mNumerator.at(iKtBin));
    mNumerator.back()->Sumw2();
    mDenominator.push_back(copy.mDenominator.at(iKtBin));
    mDenominator.back()->Sumw2();
  } // for (int iKtBin=0; iKtBin<mNKtBins; iKtBin++)
}

//____________________________
MpdFemtoDeltaEtaDeltaPhiStarMinKt::~MpdFemtoDeltaEtaDeltaPhiStarMinKt() {
  // Destructor
  if (!mNumerator.empty()) {
    mNumerator.clear();
  }
  if (!mDenominator.empty()) {
    mDenominator.clear();
  }
}

//_________________________
MpdFemtoDeltaEtaDeltaPhiStarMinKt& MpdFemtoDeltaEtaDeltaPhiStarMinKt::operator=(const MpdFemtoDeltaEtaDeltaPhiStarMinKt& copy) {
  // Assignment operator
  if (this != &copy) {
    mPhiBins = copy.mPhiBins;
    mEtaBins = copy.mEtaBins;
    mNKtBins = copy.mNKtBins;
    mKtStep = copy.mKtStep;

    for (int i = 0; i < 2; i++) {
      mEtaRange[i] = copy.mEtaRange[i];
      mPhiRange[i] = copy.mPhiRange[i];
      mKtRange[i] = copy.mKtRange[i];
    }

    // Loop over kT bins and make histogram copies
    for (int iKtBin = 0; iKtBin < mNKtBins; iKtBin++) {
      mNumerator.push_back(copy.mNumerator.at(iKtBin));
      mNumerator.back()->Sumw2();
      mDenominator.push_back(copy.mDenominator.at(iKtBin));
      mDenominator.back()->Sumw2();
    } // for (int iKtBin=0; iKtBin<mNKtBins; iKtBin++)
  } // if (this != &copy)

  return *this;
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::eventBegin(const MpdFemtoEvent* event) {
  setBField(event->bField());
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::eventEnd(const MpdFemtoEvent* /* event */) {
  /* empty */
}

//_________________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::writeOutHistos() {
  // Write out all histograms to file
  for (unsigned int i = 0; i < mNumerator.size(); i++) {
    mNumerator.at(i)->Write();
  }
  for (unsigned int i = 0; i < mDenominator.size(); i++) {
    mDenominator.at(i)->Write();
  }
}

//______________________________
TList* MpdFemtoDeltaEtaDeltaPhiStarMinKt::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  for (int i = 0; i < mNKtBins; i++) {
    tOutputList->Add(mNumerator.at(i));
    tOutputList->Add(mDenominator.at(i));
  } // for (int i=0; i<mNKtBins; i++)

  return tOutputList;
}

//_________________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::finish() {
  /* empty */
}

//____________________________
MpdFemtoString MpdFemtoDeltaEtaDeltaPhiStarMinKt::report() {
  // Construct the report
  TString report("DeltaEta vs. DeltaPhiStarMin with kT binning report:\n");
  for (int iKt = 0; iKt < mNKtBins; iKt++) {
    report += TString::Format("Number of entries in %d-th numerator                :\t%E\n",
			      iKt, mNumerator.at(iKt)->GetEntries());
    report += TString::Format("Number of entries in %d-th denominator              :\t%E\n",
			      iKt, mDenominator.at(iKt)->GetEntries());
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
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setRadiiParameters(const int& npoints,
							   const double& radiusLo,
							   const double& radiusHi) {
  if (radiusLo < radiusHi) {
    mRadii[0] = radiusLo;
    mRadii[1] = radiusHi;
  } else {
    mRadii[0] = 0.6;
    mRadii[1] = 2.;
  }

  if (npoints > 0) {
    mRadiusStep = (mRadii[1] - mRadii[0]) / npoints;
  } else {
    mRadiusStep = (mRadii[1] - mRadii[0]) / 2;
  }
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setCharges(const int& charge1, const int& charge2) {
  mQ1 = charge1;
  mQ2 = charge2;
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setBField(const double& bfield) {
  mBField = bfield;
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::addRealPair(MpdFemtoPair* pair) {
  // Perform operations on real pairs

  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  // Find index for the current kT value
  int mIndexKt = (int) ((pair->kT() - mKtRange[0]) / mKtStep);

  // Check that index is in the requested kT range
  if ((mIndexKt >= 0) && (mIndexKt < mNKtBins) && pair->qInv() < 0.15) {

    double dPhiStarMin = -1000.;
    mMomentum1 = pair->track1()->momentum();
    mMomentum2 = pair->track2()->momentum();

    dPhiStarMin = MpdFemtoPair::calculateDPhiStarMin(mMomentum1, mQ1,
						     mMomentum2, mQ2,
						     mRadiusStep, mRadii[0], mRadii[1],
						     mBField);

    // Check that histrogram exists
    if (mNumerator.at(mIndexKt)) {
      mNumerator.at(mIndexKt)->Fill(dPhiStarMin, pair->deltaEta(), 1.);
    } // if ( mNumerator.at(mIndexKt) )
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::addMixedPair(MpdFemtoPair* pair) {
  // Check front-loaded pair cut
  if (mPairCut && !mPairCut->pass(pair)) {
    return;
  }

  // Find index for the current kT value
  int mIndexKt = (int) ((pair->kT() - mKtRange[0]) / mKtStep);

  // Check that index is in the requested kT range
  if ((mIndexKt >= 0) && (mIndexKt < mNKtBins) && pair->qInv() < 0.15) {

    double dPhiStarMin = -1000.;
    mMomentum1 = pair->track1()->momentum();
    mMomentum2 = pair->track2()->momentum();

    dPhiStarMin = MpdFemtoPair::calculateDPhiStarMin(mMomentum1, mQ1,
						     mMomentum2, mQ2,
						     mRadiusStep, mRadii[0], mRadii[1],
						     mBField);

    // Check that histrogram exists
    if (mDenominator.at(mIndexKt)) {
      mDenominator.at(mIndexKt)->Fill(dPhiStarMin, pair->deltaEta(), 1.);
    } // if ( mDenominator.at(mIndexKt) )
  } // if ( (mIndexKt>=0) && (mIndexKt<mNKtBins) )
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setKtRange(const int& nKtBins, const double& kTLow, const double& kTHi) {
  if (nKtBins >= 1) {
    mNKtBins = nKtBins;
  } else {
    std::cout << "[WARNING} void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setKtRange - "
	      << "Number of kT bins must be positive" << std::endl;
    mNKtBins = 1;
  }

  if (kTLow < kTHi) {
    mKtRange[0] = kTLow;
    mKtRange[1] = kTHi;
  } else {
    std::cout << "[WARNING} void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setKtRange - "
	      << "kTLow < kTHi. Resetting to the integral" << std::endl;
    mKtRange[0] = 0.;
    mKtRange[1] = 1.5;
  }
  mKtStep = (mKtRange[1] - mKtRange[0]) / mNKtBins;
}

//_________________
void MpdFemtoDeltaEtaDeltaPhiStarMinKt::setHistoParameters(const int& nBinsEta, const double& etaLo, const double& etaHi,
							   const int& nBinsPhi, const double& phiLo, const double& phiHi) {
  if (etaLo > etaHi) {
    std::cout << "[WARNING] void MpdFemtoDeltaEtaDeltaPhiStarMinKt::void setHistoParameters - "
	      << " etaLo > etaHi. Resetting to defaults." << std::endl;
    mEtaRange[0] = -0.5;
    mEtaRange[1] = 0.5;
  } else {
    mEtaRange[0] = etaLo;
    mEtaRange[1] = etaHi;
  }

  if (nBinsEta < 1) {
    std::cout << "[WARNING] void MpdFemtoDeltaEtaDeltaPhiStarMinKt::void setHistoParameters - "
	      << " nBinsEta<1. Resetting to defaults." << std::endl;
    mEtaBins = 200;
  } else {
    mEtaBins = nBinsEta;
  }

  if (phiLo > phiHi) {
    std::cout << "[WARNING] void MpdFemtoDeltaEtaDeltaPhiStarMinKt::void setHistoParameters - "
	      << " phiLo > phiHi. Resetting to defaults." << std::endl;
    mPhiRange[0] = -0.5;
    mPhiRange[1] = 0.5;
  } else {
    mPhiRange[0] = phiLo;
    mPhiRange[1] = phiHi;
  }

  if (nBinsPhi < 1) {
    std::cout << "[WARNING] void MpdFemtoDeltaEtaDeltaPhiStarMinKt::void setHistoParameters - "
	      << " nBinsPhi<1. Resetting to defaults." << std::endl;
    mPhiBins = 200;
  } else {
    mPhiBins = nBinsPhi;
  }
}
