//
// A class that calculates 3D correlation function of identical particles
//

// MpdFemtoMakerUser headers
#include "MpdFemtoBPLCMS3DCorrFctn.h"
#include "MpdFemtoKtPairCut.h"
//#include "MpdFemtoAnalysisReactionPlane.h"

// ROOT headers
#include "TH3F.h"
#include "TList.h"
#include "TString.h"

ClassImp(MpdFemtoBPLCMS3DCorrFctn);

//_________________
MpdFemtoBPLCMS3DCorrFctn::MpdFemtoBPLCMS3DCorrFctn(const char* ctitle, const int& nbins,
						   const double& qLo, const double& qHi):
  MpdFemtoBaseCorrFctn(),
  mNumerator( nullptr ),
  mDenominator( nullptr ),
  mQinvHisto( nullptr ) {
  // Parametrized constructor

  TString name, title;
  TString appendix = ";q_{out} (GeV/c);q_{side} (GeV/c);q_{long} (GeV/c)";

  // Numerator
  name = "hBPLCMS_num_";
  name += ctitle;
  title = name;
  title += appendix;
  mNumerator = new TH3F( name.Data(), title.Data(),
			 nbins, qLo, qHi, nbins, qLo, qHi, nbins, qLo, qHi );

  // Denominator
  name = "hBPLCMS_den_";
  name += ctitle;
  title = name;
  title += apendix;
  mDenominator = new TH3F( name.Data(), title.Data(),
			   nbins, qLo, qHi, nbins, qLo, qHi, nbins, qLo, qHi );

  // qInv weighted denominator
  name = "hBPLCMS_qinvw_";
  name += ctitle;
  title = name;
  title += appendix;
  mQinvHisto = new TH3F( name.Data(), title.Data(),
			 nbins, qLo, qHi, nbins, qLo, qHi, nbins, qLo, qHi );
  
  // Weight histograms
  mNumerator->Sumw2();
  mDenominator->Sumw2();
  mQinvHisto->Sumw2();
}

//_________________
MpdFemtoBPLCMS3DCorrFctn::MpdFemtoBPLCMS3DCorrFctn(const MpdFemtoBPLCMS3DCorrFctn& c):
  MpdFemtoBaseCorrFctn( c ),
  mNumerator( nullptr ),
  mDenominator( nullptr ),
  mQinvHisto( nullptr ) {
  // Copy constructor

  mNumerator = new TH3F( *c.mNumerator );
  mDenominator = new TH3F( *c.mDenominator );
  mQinvHisto = new TH3F( *c.mQinvHisto );
}

//_________________
MpdFemtoBPLCMS3DCorrFctn::~MpdFemtoBPLCMS3DCorrFctn() {
  // Destructor
  if (mNumerator) delete mNumerator;
  if (mDenominator) delete mDenominator;
  if (mQinvHisto) delete mQinvHisto;
}

//_________________
MpdFemtoBPLCMS3DCorrFctn& MpdFemtoBPLCMS3DCorrFctn::operator=(const MpdFemtoBPLCMS3DCorrFctn& c) {
  // Assignment operator
  if (this != &c) {
    if (mNumerator) delete mNumerator;
    mNumerator = new TH3F( *c.mNumerator );
    if (mDenominator) delete mDenominator;
    mDenominator = new TH3F( *c.mDenominator );
    if (mQinvHisto) delete mQinvHisto;
    mQinvHisto = new TH3F( *c.mQinvHisto );
  }

  return *this;
}

//_________________
void MpdFemtoBPLCMS3DCorrFctn::writeOutHistos() {
  // Write out all histograms to file
  if (mNumerator) mNumerator->Write();
  if (mDenominator) mDenominator->Write();
  if (mQinvHisto) mQinvHisto->Write();
}

//_________________
TList* MpdFemtoBPLCMS3DCorrFctn::getOutputList() {
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(mNumerator);
  tOutputList->Add(mDenominator);
  tOutputList->Add(mQinvHisto);
  return tOutputList;
}

//_________________________
void MpdFemtoBPLCMS3DCorrFctn::finish() {
  /* empty */
}

//_________________
MpdFemtoString MpdFemtoBPLCMS3DCorrFctn::report() {
  // Construct the report
  TString report("LCMS Frame Bertsch-Pratt 3D Correlation Function Report:\n");
  report += TString::Format("Number of entries in numerator:\t%E\n", mNumerator->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%E\n", mDenominator->GetEntries());
  
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
void MpdFemtoBPLCMS3DCorrFctn::addRealPair( MpdFemtoPair* pair) {
  // Perform operations on real pairs
  if ( mPairCut && !mPairCut->pass(pair) ) {
    return;
  }

  mNumerator->Fill( pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS() );
}

//____________________________
void MpdFemtoBPLCMS3DCorrFctn::addMixedPair( MpdFemtoPair* pair) {
  // Perform operations on mixed pairs
  if ( mPairCut && !mPairCut->pass(pair) ) {
    return;
  }

  mDenominator->Fill( pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), 1.);
  mQinvHisto->Fill( pair->qOutCMS(), pair->qSideCMS(), pair->qLongCMS(), pair->qInv() );
}
