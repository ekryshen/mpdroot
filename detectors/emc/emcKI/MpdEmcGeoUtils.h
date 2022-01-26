#ifndef MPDEMCDETGEOUTILS_H
#define MPDEMCDETGEOUTILS_H

#include "TObject.h"

class MpdEmcGeoUtils : public TObject
{
 public:
  /// Default constructor.
  /// It must be kept public for root persistency purposes,
  /// but should never be called by the outside world
  MpdEmcGeoUtils();
  ///
  /// Destructor.
  ///
  ~MpdEmcGeoUtils() {}

  ///
  /// \return the pointer of the _existing_ unique instance of the geometry
  /// It should have been set before with GetInstance(name) method
  ///
  static MpdEmcGeoUtils* GetInstance()
  {
    if (!sGeom) {
      sGeom = new MpdEmcGeoUtils();
    }
    return sGeom;
  }

  // Check if two towers have common side (for clustering)
  // \return -1: second from prev. sector, 0: no, 1: yes, 2: towers too far apart, no sense to continue searching for
  // neighboirs
  int AreNeighbours(int detId1, int detId2) const;

  // Check if two towers have common vertex (for unfolding)
  // \return -1: second from prev. sector, 0: no, 1: yes, 2: towers too far apart, no sense to continue searching for
  // neighboirs
  int AreNeighboursVertex(int detId1, int detId2) const;

  int GeantToDetId(int chamberH, int chamber, int sector, int crate,
                   int module, int boxA, int boxB) const; // Convert Geant volume indexes to abs ID of a channel
  void DetIdToRelIndex(int detId, int& chamber, int& sector, int& iphi,
                       int& iz) const; // Convert detId to iphi,iz indexes within one sector
  void DetIdToGlobalIphiIz(int detId, int& iphi, int& iz) const; // Convert detId to iphi=0...299,iz=0...127 indexes global 

  int GetTotalNCells() const { return fNTowersPerChamber * fNChambers; }

  void DetIdToGlobalPosition(int detId, double& x, double& y,
                             double& z) const; // calculates senter of front surfase of tower with index detId

  void GetECALTubeSize(double& rMin, double& rMax, double& zMax); // Radius and length of ECAL tube

  int  MaxPhiSector(int sector) const; // number of crates (phi raws) in sector
  bool IsPreviousSector(int sector1, int sector2)const; // check if sector 2 goes before sector1

  double Rperp(double z) ; //Radius in xy plane of the surface where clusters are reconstructed 


 private:
  static MpdEmcGeoUtils* sGeom; // Pointer to the unique instance of the singleton

  int fNTowsersPerModule;   //Number of towers per module 
  int fNTowersPerCrate;    // Number of towsers per crate
  int fNTowersPerSector;   // Number of towsers per sector
  int fNCratesPerSector;    // Number of crates per sector
  int fNCratesPerSector2;   // Number of crates per small sector
  int fNSectorsPerChamber;  // Number of Sectors per Chamber
  int fNChambers;           // Number of chambers
  int fNTowersPerChamber;   // number of tower per chamber

  float fEcalRmin; // Minimal radius of ECAL
  float fEcalRmax; // Maximal radius of ECAL
  float fEcalZmax; // Maximal lenght of ECAL in z direction

  ClassDef(MpdEmcGeoUtils, 1)
};
#endif
