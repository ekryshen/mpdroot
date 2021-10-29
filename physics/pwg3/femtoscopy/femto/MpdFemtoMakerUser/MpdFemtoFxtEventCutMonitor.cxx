//
// Event cut monitor for basic analysis
//

// MpdFemtoMaker headers
#include "MpdFemtoEvent.h"
#include "MpdFemtoTypes.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoFxtEventCutMonitor.h"

// C++ headers
#include <cstdio>
#include <cmath>

ClassImp(MpdFemtoFxtEventCutMonitor)

//_________________
MpdFemtoFxtEventCutMonitor::MpdFemtoFxtEventCutMonitor() :
MpdFemtoBaseCutMonitor(),
mVertexZ(nullptr),
mRefMult(nullptr),
mVertexYvsVertexX(nullptr),
mVpdVzDiff(nullptr),
mBTofTrayMultVsRefMult(nullptr),
mBTofMatchedVsRefMult(nullptr) {
  // Constructor
  mVertexYvsVertexX = new TH2F("VertexYvsVertexX",
			       "VertexYvsVertexX;x (cm);y (cm)",
			       300, -6., 6., 300, -6., 6.);
  mVertexZ = new TH1F("VertexZ", "VertexZ; z (cm)",
		      480, -120., 120.);
  mRefMult = new TH1F("RefMult", "RefMult;RefMult",
		      600, -0.5, 599.5);
  mVpdVzDiff = new TH1F("VpdVzDiff", "mVpdVzDiff",
			80, -20., 20.);
  mBTofTrayMultVsRefMult = new TH2F("mBTofTrayMultVsRefMult",
				    "mBTofTrayMultVsRefMult;RefMult;BTof tray multiplicity",
				    600, -0.5, 599.5, 600, -0.5, 599.5);
  mBTofMatchedVsRefMult = new TH2F("mBTofMatchedVsRefMult",
				   "mBTofMatchedVsRefMult;RefMult;BTof-matched multiplicity",
				   600, -0.5, 599.5, 600, -0.5, 599.5);
}

//_________________
MpdFemtoFxtEventCutMonitor::MpdFemtoFxtEventCutMonitor(const char* title1, const char* title2) :
MpdFemtoBaseCutMonitor(),
mVertexZ(nullptr),
mRefMult(nullptr),
mVertexYvsVertexX(nullptr),
mVpdVzDiff(nullptr),
mBTofTrayMultVsRefMult(nullptr),
mBTofMatchedVsRefMult(nullptr) {
  // Parametrized constructor
  char tit1[100];

  sprintf(tit1, "%s%s_VertexYvsVertexX", title1, title2);
  mVertexYvsVertexX = new TH2F(tit1, "VertexYvsVertexX;x (cm);y (cm)",
			       300, -6., 6., 300, -6., 6.);

  sprintf(tit1, "%s%s_VertexZ", title1, title2);
  mVertexZ = new TH1F(tit1, "VertexZ;z (cm);Entries",
		      480, -120., 120.);

  sprintf(tit1, "%s%s_RefMult", title1, title2);
  mRefMult = new TH1F(tit1, "RefMult;RefMult;Entries",
		      600, -0.5, 599.5);

  sprintf(tit1, "%s%s_VpdVzDiff", title1, title2);
  mVpdVzDiff = new TH1F(tit1, "VpdVzDiff", 40, -20., 20.);

  sprintf(tit1, "%s%s_BTofTrayMultVsRefMult", title1, title2);
  mBTofTrayMultVsRefMult = new TH2F(tit1, "mBTofTrayMultVsRefMult;RefMult;BTof tray multiplicity",
				    600, -0.5, 599.5, 600, -0.5, 599.5);

  sprintf(tit1, "%s%s_BTofMatchedVsRefMult", title1, title2);
  mBTofMatchedVsRefMult = new TH2F(tit1, "mBTofMatchedVsRefMult;RefMult;BTof-matched multiplicity",
				   600, -0.5, 599.5, 600, -0.5, 599.5);
}

//_________________
MpdFemtoFxtEventCutMonitor::MpdFemtoFxtEventCutMonitor(const MpdFemtoFxtEventCutMonitor& c) : MpdFemtoBaseCutMonitor(c) {
  // Copy constructor
  if (mVertexZ) delete mVertexZ;
  mVertexZ = new TH1F(*c.mVertexZ);
  if (mRefMult) delete mRefMult;
  mRefMult = new TH1F(*c.mRefMult);
  if (mVertexYvsVertexX) delete mVertexYvsVertexX;
  mVertexYvsVertexX = new TH2F(*c.mVertexYvsVertexX);
  if (mVpdVzDiff) delete mVpdVzDiff;
  mVpdVzDiff = new TH1F(*c.mVpdVzDiff);
  if (mBTofTrayMultVsRefMult) delete mBTofTrayMultVsRefMult;
  mBTofTrayMultVsRefMult = new TH2F(*c.mBTofTrayMultVsRefMult);
  if (mBTofMatchedVsRefMult) delete mBTofMatchedVsRefMult;
  mBTofMatchedVsRefMult = new TH2F(*c.mBTofMatchedVsRefMult);
}

//_________________
MpdFemtoFxtEventCutMonitor MpdFemtoFxtEventCutMonitor::operator=(const MpdFemtoFxtEventCutMonitor& c) {
  // Assignment operator
  if (this != &c) {
    if (mVertexZ) delete mVertexZ;
    mVertexZ = new TH1F(*(c.mVertexZ));
    if (mRefMult) delete mRefMult;
    mRefMult = new TH1F(*(c.mRefMult));
    if (mVertexYvsVertexX) delete mVertexYvsVertexX;
    mVertexYvsVertexX = new TH2F(*(c.mVertexYvsVertexX));
    if (mVpdVzDiff) delete mVpdVzDiff;
    mVpdVzDiff = new TH1F(*(c.mVpdVzDiff));
    if (mBTofTrayMultVsRefMult) delete mBTofTrayMultVsRefMult;
    mBTofTrayMultVsRefMult = new TH2F(*(c.mBTofTrayMultVsRefMult));
    if (mBTofMatchedVsRefMult) delete mBTofMatchedVsRefMult;
    mBTofMatchedVsRefMult = new TH2F(*(c.mBTofMatchedVsRefMult));
  }

  return *this;
}

//_________________
MpdFemtoFxtEventCutMonitor::~MpdFemtoFxtEventCutMonitor() {
  // Destructor
  if (mVertexYvsVertexX) {
    delete mVertexYvsVertexX;
    mVertexYvsVertexX = nullptr;
  }
  if (mVertexZ) {
    delete mVertexZ;
    mVertexZ = nullptr;
  }
  if (mRefMult) {
    delete mRefMult;
    mRefMult = nullptr;
  }
  if (mVpdVzDiff) {
    delete mVpdVzDiff;
    mVpdVzDiff = nullptr;
  }
  if (mBTofTrayMultVsRefMult) {
    delete mBTofTrayMultVsRefMult;
    mBTofTrayMultVsRefMult = nullptr;
  }
  if (mBTofMatchedVsRefMult) {
    delete mBTofMatchedVsRefMult;
    mBTofMatchedVsRefMult = nullptr;
  }
}

//_________________
void MpdFemtoFxtEventCutMonitor::fill(const MpdFemtoEvent* event) {
  // Fill monitor histograms
  mVertexYvsVertexX->Fill(event->primaryVertex().x(), event->primaryVertex().y());
  mVertexZ->Fill(event->primaryVertex().z());
  mRefMult->Fill(event->refMult());
  mVpdVzDiff->Fill(event->vpdVzDiff());
  mBTofTrayMultVsRefMult->Fill(event->refMult(), event->btofTrayMultiplicity());
  mBTofMatchedVsRefMult->Fill(event->refMult(), event->numberOfBTofMatched());
}

//_________________
void MpdFemtoFxtEventCutMonitor::finish() {
  std::cout << " entries in histogram mVertexYvsVertexX : "
            << mVertexYvsVertexX->Integral() << std::endl;
  std::cout << " entries in histogram mVertexZ : "
            << mVertexZ->Integral() << std::endl;
  std::cout << " entries in histogram mRefMult : "
            << mRefMult->Integral() << std::endl;
  std::cout << " entries in histogram mVpdVzDiff : "
            << mVpdVzDiff->Integral() << std::endl;
  std::cout << " entries in histogram mBTofTrayMultVsRefMult : "
            << mBTofTrayMultVsRefMult->Integral() << std::endl;
  std::cout << " entries in histogram mBTofMatchedVsRefMult : "
            << mBTofMatchedVsRefMult->Integral() << std::endl;
}

//_________________
MpdFemtoString MpdFemtoFxtEventCutMonitor::report() {
  string Stemp;
  char Ctemp[100];
  sprintf(Ctemp, " MpdFemtoFxtEventCutMonitor");
  Stemp = Ctemp;
  MpdFemtoString returnThis = Stemp;
  return returnThis;
}

//_________________
void MpdFemtoFxtEventCutMonitor::writeOutHistos() {
  // Write all histograms
  mVertexYvsVertexX->Write();
  mVertexZ->Write();
  mRefMult->Write();
  mVpdVzDiff->Write();
  mBTofTrayMultVsRefMult->Write();
  mBTofMatchedVsRefMult->Write();
}

//_________________
TList* MpdFemtoFxtEventCutMonitor::getOutputList() {
  // Prepare ouput list with histograms
  TList *outputList = new TList();

  outputList->Add(mVertexYvsVertexX);
  outputList->Add(mVertexZ);
  outputList->Add(mRefMult);
  outputList->Add(mVpdVzDiff);
  outputList->Add(mBTofTrayMultVsRefMult);
  outputList->Add(mBTofMatchedVsRefMult);

  return outputList;
}
