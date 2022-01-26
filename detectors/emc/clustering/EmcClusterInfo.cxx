// $Id$
// Author: artur   2016/04/25


//_____________________________________________________________________________
//
// EmcClusterInfo
//_____________________________________________________________________________

#include "EmcClusterInfo.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(EmcClusterInfo)
ClassImp(EmcClusterInfo0)

//_____________________________________________________________________________
EmcClusterInfo::EmcClusterInfo():fId(-1),fId2(-1),fFlag(-1),fZ(-1),fPhi(-1),fDeposit(0)
{
 
}

//_____________________________________________________________________________
EmcClusterInfo::~EmcClusterInfo() 
{
}

//_____________________________________________________________________________
void EmcClusterInfo::Reset()
{
  fId = -1;
  fId2 = -1;
  fFlag = -1;
  fZ = -1;
  fPhi = -1;
  fDeposit = 0; 
}

//_____________________________________________________________________________
void EmcClusterInfo::AddPoint(Int_t iz, Int_t iphi, Float_t dep)
{
  CLPointsMap::iterator it = fPoints.find(CLPoint(iz,iphi));
  if (it != fPoints.end()) it->second += dep;    
  else fPoints.insert(std::pair<CLPoint,Float_t>(CLPoint(iz,iphi),dep));
  fDeposit += dep;
}


//_____________________________________________________________________________
void EmcClusterInfo::Print(const Option_t* opt) const
{
  cout << "-I- <EmcClusterInfo::Print>" << endl;
  
  cout << "Id:               " << fId << "," << fId2 << endl;
  cout << "Flag:             " << fFlag << endl;
  cout << "Center(phi,z):    " << fPhi << " " << fZ << endl;
  cout << "Points/DEposit:   " << GetNPoints() << "/" << fDeposit << endl;
  cout << endl;
  
  CLPointsMap::const_iterator it = fPoints.begin();
  int n(0);
  for (; it != fPoints.end(); ++it) {
       printf("%3d [%4d, %4d] %12.6f\n",++n,
             it->first.first,it->first.second,it->second);
  }
  cout << endl;
}

