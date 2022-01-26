////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcGeoUtils                                            //
//  Geometrical transformations for EMC                	      //
//  Author List : D.Peresunko, KI, 2019                       //
//                                                            //
////////////////////////////////////////////////////////////////
#include "MpdEmcGeoUtils.h"

#include "FairLogger.h"

#include "TGeoCompositeShape.h"
#include "TGeoArb8.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "MpdEmcSimParams.h"

#include <cstdlib>

// these initialisations are needed for a singleton
MpdEmcGeoUtils* MpdEmcGeoUtils::sGeom = nullptr;

ClassImp(MpdEmcGeoUtils)

  // -----   Default constructor   -------------------------------------------
  MpdEmcGeoUtils::MpdEmcGeoUtils()
  : fNTowsersPerModule(16), 
    fNTowersPerCrate(128), 
    fNTowersPerSector(768),
    fNCratesPerSector(6),
    fNSectorsPerChamber(25),      
    fNChambers(2),           
    fNTowersPerChamber(19200),     
    fEcalRmin(0.),           // minimal ECAL radius
    fEcalRmax(0.),           // maximal ECAL radius
    fEcalZmax(0.)            // ECAL half z size
{

}

int MpdEmcGeoUtils::AreNeighbours(int detId1, int detId2) const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  // neighbours are defined as digits having at least a common side
  // The order of detId1 and detId2 is important: first (d1) should be a digit already in a cluster
  //                                      which is compared to a digit (d2)  not yet in a cluster

  int ch1, sector1, iphi1, iz1;
  DetIdToRelIndex(detId1, ch1, sector1, iphi1, iz1);
  int ch2, sector2, iphi2, iz2;
  DetIdToRelIndex(detId2, ch2, sector2, iphi2, iz2);

  if (ch1 > ch2){
    return -1;
  }
  if (ch1 < ch2){
    return 2;
  }

  if(MpdEmcSimParams::GetInstance()->AllowMultiSectorClusters()){

    if (sector1 == sector2) { // start looking at digits in next sector, can not create inter-sector clusters

      int phidiff = std::abs(iphi1 - iphi2);
      int zdiff = std::abs(iz1 - iz2);

      if (((phidiff <= 1) && (zdiff == 0)) || ((phidiff == 0) && (zdiff <= 1))) { // At least common side
        return 1;
      } else {
        if ((iphi2 > iphi1 + 1)) {
	        if(iphi1!=0 && iphi1!=MaxPhiSector(sector1) && sector2>sector1){ //not at the edge and current sector alredy scanned
            return 2; //  Difference in iphi numbers is too large to look further
	        }
	        else{ //current digit at the edge, continue looking at other sectors
            return 0 ;
	        }
        }
        return 0; // continue searching
      }
    }
    //May be edge of sector
    if(IsPreviousSector(sector1,sector2)){  // second digit in previous sector
      if(iphi1==0 && iphi2==MaxPhiSector(sector2) && iz1==iz2){
        return 1;
      }
      else{
        return 0;
      }
    }
    else{
      if(IsPreviousSector(sector2,sector1)){ //second digit in next sector
        if(iphi1==MaxPhiSector(sector1)){
          if(iphi2==0 && iz1==iz2){
            return 1;
          }
          else{
            return 0;
          }
        }
        else{
          if(iphi1!=0 && sector2>sector1 ){ //not at the edge and current sector scanned
            return 2; //  Difference in iphi numbers is too large to look further
          }
          else{ //current digit at the edge, continue looking at other sectors
            return 0 ;
          }
        }
      }
      else{ //not adjacent sectors
        if(iphi1!=0 && iphi1!=MaxPhiSector(sector1) && sector2>sector1 ){
          return 2; //  Difference in iphi numbers is too large to look further
        }
        else{ //current digit at the edge, continue looking at other sectors
          return 0 ;
        }
      }
    }
  }
  else{ // only within sector
    if (sector1 > sector2) { // start looking at digits in next sector, can not create inter-sector clusters
      return -1;
    }
    if (sector1 < sector2) {
      return 2;
    }

    int phidiff = std::abs(iphi1 - iphi2);
    int zdiff = std::abs(iz1 - iz2);

    if (((phidiff <= 1) && (zdiff == 0)) || ((phidiff == 0) && (zdiff <= 1))) { // At least common side
      return 1;
    } else {
      if ((iphi2 > iphi1 + 1)) {
        return 2; //  Difference in iphi numbers is too large to look further
      }

      return 0; // continue searching
    }
  }

}
int MpdEmcGeoUtils::AreNeighboursVertex(int detId1, int detId2) const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching
  //                                       = 1 are neighbour
  // In contrast to previous method, here neighbours are defined as digits having at least a common vertex
  // The order of detId1 and detId2 is important: first (d1) should be a digit already in a cluster
  //                                      which is compared to a digit (d2)  not yet in a cluster

  int ch1, sector1, iphi1, iz1;
  DetIdToRelIndex(detId1, ch1, sector1, iphi1, iz1);
  int ch2, sector2, iphi2, iz2;
  DetIdToRelIndex(detId2, ch2, sector2, iphi2, iz2);

  if (ch1 != ch2)
    return 0;

  if(MpdEmcSimParams::GetInstance()->AllowMultiSectorClusters()){

    if(sector1 == sector2){
      int phidiff = std::abs(iphi1 - iphi2);
      int zdiff = std::abs(iz1 - iz2);

      if (((phidiff <= 1) && (zdiff <= 1))) { // At least common vertex
        return 1;
      } else {
        return 0; // not a neighbour
      }
    }

    //May be edge of sector
    if(IsPreviousSector(sector1,sector2)){  //previous sector
      if(iphi1==0 && iphi2==MaxPhiSector(sector2) && TMath::Abs(iz1-iz2)<=1){
        return 1;
      }
      else{
        return 0;
      }
    }
    if(IsPreviousSector(sector2,sector1)){
      if(iphi1==MaxPhiSector(sector1) && iphi2==0 && TMath::Abs(iz1-iz2)<=1){
        return 1;
      }
    }
    return 0;
  }
  else{ //each sector separately
    if (sector1 != sector2) { // start looking at digits in next sector, can not create inter-sector clusters
      return 0;
    }

    int phidiff = std::abs(iphi1 - iphi2);
    int zdiff = std::abs(iz1 - iz2);

    if (((phidiff <= 1) && (zdiff <= 1))) { // At least common vertex
      return 1;
    } else {
      return 0; // continue searching
    }
  }
}

int MpdEmcGeoUtils::GeantToDetId(int chamberH, int chamber, int sector, int crate, 
                                 int module, int boxA, int boxB) const // Convert Geant volume indexes to abs ID of a channel
{
  // Convert Geant volume indexes to abs ID of a channel
  int num = boxB + 
            crate  * fNTowersPerCrate +  
            sector * fNTowersPerSector + 
            chamber* fNTowersPerChamber ;
  return num;
}
void MpdEmcGeoUtils::DetIdToRelIndex(int detId, int& chamber, int& sector, int& iphi, int& iz) const
{
  // Convert detId to iphi,iz indexes within one sector
  // We use consequent numbering of sectors

  //int chamberH = 0;
  chamber = detId / fNTowersPerChamber;
  detId -= chamber * fNTowersPerChamber;
  sector = detId / fNTowersPerSector;
  detId-=sector * fNTowersPerSector;
  int crate = detId/ fNTowersPerCrate ;
  detId-= crate * fNTowersPerCrate ;
  iz = detId / 2;
  iphi = detId % 2 + 2 * crate;   
}
void MpdEmcGeoUtils::DetIdToGlobalIphiIz(int detId, int& iphi, int& iz) const
{
  int chamber,sect ;
  DetIdToRelIndex(detId,chamber,sect,iphi,iz) ;
  iphi=fNCratesPerSector*2*sect+iphi;
  iz  =chamber*fNTowersPerCrate/2+iz ;
}
void MpdEmcGeoUtils::DetIdToGlobalPosition(int detId, double& x, double& y, double& z) const
{
  // calculates senter of front surfase of tower with index detId
  if (!gGeoManager) {
    LOG(ERROR) << "Can not run without constructed geometry";
    x = y = z = 0;
    return;
  }


  int chamber = detId / (fNTowersPerChamber);
  detId -= chamber * fNTowersPerChamber;
  int sector = detId / fNTowersPerSector;
  detId-=sector* fNTowersPerSector ;
  int crate = detId / fNTowersPerCrate;
  detId-= crate * fNTowersPerCrate ;
  int module = detId/fNTowsersPerModule ;

  // Construct geant pass from detId:
  TString path(Form("cave_1/emcChamber_0/emcChH_%d/emcSector_%d/emcCrate_%d/emcModule%d_0/emc_box%d_%d/", 
       chamber,sector,crate, module, detId/2 + 1, detId));

  if (gGeoManager->CheckPath(path.Data())) {
    gGeoManager->cd(path.Data());
    TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
    // check sizes of this tower
    TGeoVolume* v = gGeoManager->GetVolume(Form("emc_box%d", detId/2 + 1));
    double towerHalfSizeZ =0 ;
    TGeoCompositeShape* sh = dynamic_cast<TGeoCompositeShape*>(v->GetShape());
    if(sh){
      towerHalfSizeZ = sh->GetDZ();
    }
    else{
      TGeoArb8 * sh2 = dynamic_cast<TGeoArb8*>(v->GetShape());
      towerHalfSizeZ = sh2->GetDZ();
    }

towerHalfSizeZ*=0.2; //V
    
    double local[3] = { 0 }, global[3] = { 0 };
    if (chamber == 0) { // Why this is different in Chamber0 and Chamber1? Reflection?
      local[2] = -towerHalfSizeZ;
    } else {
      local[2] = towerHalfSizeZ;
    }
    m->LocalToMaster(local, global);
    x = global[0];
    y = global[1];
    z = global[2];
  } else {
    LOG(ERROR) << "Can not find volume " << path;
    x = y = z = 0;
  }
}

int MpdEmcGeoUtils::MaxPhiSector(int sector)const{
  //returns number of crates per sector
    return fNCratesPerSector*2-1;
}
bool MpdEmcGeoUtils::IsPreviousSector(int sector1, int sector2)const{
  // check if sector1 goes before sector 2
  if(sector1==0) return sector2==24;
  else return sector2==sector1-1;
}
void MpdEmcGeoUtils::GetECALTubeSize(double& rMin, double& rMax, double& zMax)
{
  if (fEcalRmin == 0) {
    TGeoVolume* v = gGeoManager->GetVolume("emcChH"); //TODO!!!!! check!!!
    TGeoTube* sh = dynamic_cast<TGeoTube*>(v->GetShape());
    fEcalRmin = sh->GetRmin();
    fEcalRmax = sh->GetRmax();
    fEcalZmax = sh->GetDz();
  }
  rMin = fEcalRmin;
  rMax = fEcalRmax;
  zMax = fEcalZmax;
}
double MpdEmcGeoUtils::Rperp(double z){
  //Parameterization of the radius in xy plane ot the surface of EMC clusters
  //Parameterization of singe photon simulation from V Riabov (1.12.2019)
  //Chi2                      =       162512  
  //NDf                       =          271
  //p0                        =      188.186   +/-   6.62242e-05
  //p1                        =   0.00882988   +/-   2.19833e-06
  //p2                        = -0.000191223   +/-   1.9986e-08  
  //p3                        =   2.8476e-07   +/-   5.01312e-11 
  const double p0 = 188.186;
  const double p1 = 0.00882988 ;
  const double p2 =-0.000191223 ;
  const double p3 = 2.8476e-07 ;
  return p0 + p1*TMath::Abs(z) + p2*z*z + p3*TMath::Abs(z*z*z) ;  
}