//
//  Simple cut to select pairs within the given kT range
//

// MpdFemtoMakerUser headers
#include "MpdFemtoKtPairCut.h"

// ROOT headers
#include <TMath.h>

// C++ headers
#include <string>

ClassImp(MpdFemtoKtPairCut);

//_________________
MpdFemtoKtPairCut::MpdFemtoKtPairCut() : MpdFemtoBasePairCut() {
  // Default constructor
  fKtRange[0] = 0.f;
  fKtRange[1] = 1e6;
  fPhiRange[0] = -1e6;
  fPhiRange[1] = 1e6;
}

//_________________
MpdFemtoKtPairCut::MpdFemtoKtPairCut(double ktLo, double ktHi,
				     double phiLo, double phiHi): MpdFemtoBasePairCut() {
  // Parametrized constructor
  fKtRange[0] = ktLo;
  fKtRange[1] = ktHi;
  fPhiRange[0] = phiLo;
  fPhiRange[1] = phiHi;
}

//_________________
MpdFemtoKtPairCut::MpdFemtoKtPairCut(const MpdFemtoKtPairCut& c):
  MpdFemtoBasePairCut(c) {
  // Copy constructor
  fKtRange[0] = c.fKtRange[0];
  fKtRange[1] = c.fKtRange[1];
  fPhiRange[0] = c.fPhiRange[0];
  fPhiRange[1] = c.fPhiRange[1];
}

//_________________
MpdFemtoKtPairCut::~MpdFemtoKtPairCut() {
  // Destructor
  /* emtpy */
}

//_________________
MpdFemtoKtPairCut& MpdFemtoKtPairCut::operator=(const MpdFemtoKtPairCut& c) {
  // Assignment operator
  if (this != &c) {
    MpdFemtoBasePairCut::operator=(c);
    fKtRange[0] = c.fKtRange[0];
    fKtRange[1] = c.fKtRange[1];
    fPhiRange[0] = c.fPhiRange[0];
    fPhiRange[1] = c.fPhiRange[1];
  }

  return *this;
}

//_________________
MpdFemtoString MpdFemtoKtPairCut::report() {
  
  // Prepare a report from the execution
  MpdFemtoString report = "MpdFemtoKtPairCut \n";
  report += Form("Accept pair with kT in range %4.2f , %4.2f", fKtRange[0], fKtRange[1] );
  report += Form("Accept pair with angle in range %4.2f , %4.2f", fPhiRange[0], fPhiRange[1] );

  return report;
}

//_________________
TList *MpdFemtoKtPairCut::listSettings() {
  // Return a list of settings in a writable form
  TList *settings = new TList();

  settings->AddVector(
    new TObjString( TString::Format("MpdFemtoKtPairCut.ktmax=%4.2f", fKtRange[0] ) ),
    new TObjString( TString::Format("MpdFemtoKtPairCut.ktmin=%4.2f", fKtRange[1] ) ),
    new TObjString( TString::Format("MpdFemtoKtPairCut.phimax=%4.2f", fPhiRange[0] ) ),
    new TObjString( TString::Format("MpdFemtoKtPairCut.phimin=%4.2f", fPhiRange[1] ) ),
    nullptr
  );

  return settings;
}

//_________________
bool MpdFemtoKtPairCut::pass(const MpdFemtoPair* pair) {

  bool isGoodPair = false;
  
  // Apply kT pair cut ( IMPORTANT: a <= x < b )
  isGoodPair = ( fKtRange[0] <= pair->kT() && pair->kT() < fKtRange[1] ) ? true : false;

  /*

    // To be filled with phi cut
    double rpangle = pair->GetPairAngleEP();

    if (rpangle > 180.0) rpangle -= 180.0;
    if (rpangle < 0.0) rpangle += 180.0;
    
    // note: handle "wrap around" if the minimum phi is negative
    if (fPhiMin < 0.0) {
    return (rpangle < fPhiMax) || (180.0+fPhiMin <= rpangle);
    }
    
    // return whether angle is within phi-range
    return (fPhiMin <= rpangle) && (rpangle < fPhiMax);
  */

  return isGoodPair;
}

//_________________
bool MpdFemtoKtPairCut::pass(const MpdFemtoPair* pair, double aRPAngle) {
  // The same as above, but it is defined with RP Angle as input in
  // all the correlation function classes.
  return ( (aRPAngle > 0.) && pass(pair) ) ;
}
