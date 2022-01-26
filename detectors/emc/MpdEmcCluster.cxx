////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcHitCluster                                          //
//  Cluster production for EMC, v03                           //
//  Author List : Martemianov M., ITEP, 2021                  //
//                                                            //
////////////////////////////////////////////////////////////////

#include <iostream>
#include <numeric>
#include "TMath.h"
#include "MpdEmcCluster.h"
#include "MpdEmcClusterCreation.h"

using namespace std;
using namespace TMath;

// -----   Default constructor   -------------------------------------------

MpdEmcCluster::MpdEmcCluster() : fFlag(-1), fE(0.), fTime(0.), fPhi(0.), fRho(0.), fZ(0.) {}

// -----   Standard constructor   ------------------------------------------

MpdEmcCluster::MpdEmcCluster(Float_t energy, Float_t time, TVector3 pos)
   : fFlag(1), fE(energy), fTime(time), fPhi(pos.X()), fRho(pos.Y()), fZ(pos.Z())
{
}

// -----   Destructor   ----------------------------------------------------

MpdEmcCluster::~MpdEmcCluster() {}

// -----  Print  -----------------------------------------------------------

void MpdEmcCluster::Print(const Option_t *opt) const
{
   cout << "MpdEmcCluster: " << endl;
   cout << "\tDeposited energy: " << fE << "\tMean time: " << fTime << "   Rho cluster: " << fRho
        << "   Phi cluster: " << fPhi << "   Z cluster: " << fZ << endl;

   //    cout << "\tNumber of tracks in module: " << fNumTracks << endl;
   //    if (fNumTracks == 1) cout << "PDG code: " << fPDG << "   Track ID: " << fTrackID << endl;
}
// -------------------------------------------------------------------------

// Compute claster radius

Float_t MpdEmcCluster::ComputeClusterRadius(TVector3 fPos, vector<Float_t> xC, vector<Float_t> yC, vector<Float_t> zC,
                                            vector<Float_t> eC)
{

   Float_t phi, theta, valRad, fEnergy, rhoCenter, thetaCenter;
   fEnergy     = accumulate(eC.begin(), eC.end(), 0.0);
   rhoCenter   = sqrt(fPos.Y() * fPos.Y() + fPos.Z() * fPos.Z());
   thetaCenter = ACos(fPos.Z() / rhoCenter);

   valRad = 0.0;
   for (UInt_t iHit = 0; iHit < eC.size(); iHit++) {
      phi = ATan2(yC[iHit], xC[iHit]) * RadToDeg();
      if (phi < 0) phi += 360;
      theta = ACos(zC[iHit] / sqrt(xC[iHit] * xC[iHit] + yC[iHit] * yC[iHit] + zC[iHit] * zC[iHit]));

      valRad += (pow(rhoCenter * (theta - thetaCenter), 2) +
                 pow(fPos.Y() * MpdEmcMath::GetPhiDiff(phi, fPos.X()) * DegToRad(), 2)) *
                eC[iHit];

      /*
            valRad += (pow(sqrt(xC[iHit]*xC[iHit]+yC[iHit]*yC[iHit]+zC[iHit]*zC[iHit])*(theta-thetaCenter),2)+
                        pow(sqrt(xC[iHit]*xC[iHit]+yC[iHit]*yC[iHit])*MpdEmcMath::GetPhiDiff(phi,fPos.X())*DegToRad(),2))*eC[iHit];

      */
   }
   return (sqrt(valRad / fEnergy));
}

ClassImp(MpdEmcCluster)
