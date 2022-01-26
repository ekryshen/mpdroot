// $Id$
// Author: artur   2016/04/27


//_____________________________________________________________________________
//
// MpdEmcGeoParWrapper
//_____________________________________________________________________________

#include "MpdEmcGeoParWrapper.h"
#include <TMath.h>

#include <iostream>
using std::cout;
using std::endl;

using namespace TMath;

ClassImp(MpdEmcGeoParWrapper)

//_____________________________________________________________________________
MpdEmcGeoParWrapper::MpdEmcGeoParWrapper(MpdEmcGeoPar* pars, Bool_t is_owner)
:fGeoPar(0),fIsOwner(is_owner),fIsInit(kFALSE),
fNPhi(1),fDPhi(0),fPhimin(0),
fNZ(1),fDZ(0),fZmin(0),fZmax(0),fTmin(0),fTmax(0),
fR(0)
{
  if (pars) SetParameters(pars);
}

//_____________________________________________________________________________
MpdEmcGeoParWrapper::~MpdEmcGeoParWrapper() 
{
   if (fIsOwner) delete fGeoPar;
}

//_____________________________________________________________________________
void MpdEmcGeoParWrapper::SetParameters(MpdEmcGeoPar* pars) 
{ 
   fGeoPar = pars;  
   
   fPhimin = 0; // Check it! -> OK (number of supermodules in sector is even = 4)
   
   //fNPhi = fGeoPar->GetNsec()*fGeoPar->GetNsupMod()*fGeoPar->GetNModInSuperModByPhi();
   fNPhi = fGeoPar->GetNsec()*fGeoPar->GetNsupMod()*3; //FIXME!!! 3 -> fGeoPar->GetNModInSuperModByPhi();
   fDPhi = 360./fNPhi;
   
   fZmax = fGeoPar->GetLength(); // cm
   fZmin = -fZmax;
   fNZ = 2*fGeoPar->GetNrows()*fGeoPar->GetNModInSuperModByZ();
   fDZ = (fZmax - fZmin)/fNZ;
   
   fR = fGeoPar->GetRmin()*0.1; // (!) mm -> cm
   
   fTmin = TMath::ATan2(fR,fZmax)*TMath::RadToDeg();  
   fTmax = 180. - fTmin;
   
   fIsInit = kTRUE;
} 

//_____________________________________________________________________________
Float_t MpdEmcGeoParWrapper::GetAngle(Float_t x, Float_t y) const
{
   Float_t ang = ATan2(y, x);
   if (ang < 0) ang += TwoPi();
   ang *= RadToDeg();
   return ang;
}

//_____________________________________________________________________________
Int_t MpdEmcGeoParWrapper::GetModUid(Float_t x, Float_t y, Float_t z) const
{
   if (!fIsInit) return -1;  
   
   Int_t nz, np;
   
   if (z < fZmin) nz = 0;
   else if (z > fZmax) nz = fNZ-1;
   else nz = Int_t((z-fZmin)/fDZ);
   if (nz >= fNZ) nz = fNZ-1; 
   
   Float_t ang = GetAngle(x,y)-fPhimin;
 
   if (ang < 0) ang += 360;
   np = Int_t(ang/fDPhi);
   if (np >= fNPhi) np = fNPhi-1;
   
   //cout << "uid = " << np + nz*fNPhi << " nz = " << nz << " np = " << np << endl;
   //cout << "<MpdEmcGeoParWrapper::GetModUid> " << " ang = " << ang  
   //     << " x = " << x << " y = " << y << endl;
   
   return np + nz*fNPhi;
}

//_____________________________________________________________________________
Bool_t MpdEmcGeoParWrapper::GetModPos(Int_t uid, Float_t& x, Float_t& y, Float_t& z) const
{
   if (!fIsInit) return kFALSE;
  
   Int_t np = uid%fNPhi;
   Int_t nz = uid/fNPhi;
   
   z = fZmin + (nz + 0.5)*fDZ;
   
   Float_t ang = (fPhimin + (np + 0.5)*fDPhi)*DegToRad();
   
   x = fR*Cos(ang); 
   y = fR*Sin(ang); 
   return kTRUE;
}

//_____________________________________________________________________________
void MpdEmcGeoParWrapper::Print(const Option_t* opt) const
{
   cout << "[I]<MpdEmcGeoParWrapper::Print> " << endl;
   cout << "N(z)   = " << fNZ   << " dZ   = " << fDZ   << " [cm]  " << endl;
   cout << "N(phi) = " << fNPhi << " dPhi = " << fDPhi << " [deg] " 
        << " " << fR*fDPhi*DegToRad() << " [cm] " 
        << " T(min,max) = " << fTmin << ", " << fTmax << " [deg] " << endl;
}

