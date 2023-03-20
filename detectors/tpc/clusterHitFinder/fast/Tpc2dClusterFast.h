//-------------------------------------------------------------------------------------------------
// Description:
//      Tpc2dClusterFast class is an adapter class (a data structure)
//      It is designed so that:
//      - its purpose is accessing clusters to generate QA plots for Fast module
//      - it takes the values from Cluster class in TpcClustering.h
//      - the file TpcClustering.h is not modified
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, March, 2023
//-------------------------------------------------------------------------------------------------

#ifndef TPC2DCLUSTERFAST_HH
#define TPC2DCLUSTERFAST_HH

#include "AbstractTpc2dCluster.h"
#include "TpcClustering.h"

class Tpc2dClusterFast : public AbstractTpc2dCluster {
public:
   // Constructors/Destructors ---------
   Tpc2dClusterFast() {}
   virtual ~Tpc2dClusterFast() {}

   // interface implementation
   int                             GetClusterID() const { return clusterID; }
   int                             GetSector() const { return sector; }
   int                             GetRow() const { return row; }
   int                             GetPadMin() const { return padMin; }
   int                             GetPadMax() const { return padMax; }
   int                             GetTimeBinMin() const { return timeBinMin; }
   int                             GetTimeBinMax() const { return timeBinMax; }
   std::vector<AbstractTpcDigit *> GetClusterDigits() const { return digits; }

   // setters (to comply with TClonesArray clusArray way of dealing with things in TpcClusterHitFinderFast.cxx)
   void SetClusterID(int newClusterID) { clusterID = newClusterID; }
   void SetSector(int newSector) { sector = newSector; }
   void SetRow(int newRow) { row = newRow; }
   void SetPadMin(int newPadMin) { padMin = newPadMin; }
   void SetPadMax(int newPadMax) { padMax = newPadMax; }
   void SetTimeBinMin(int newTimeBinMin) { timeBinMin = newTimeBinMin; }
   void SetTimeBinMax(int newTimeBinMax) { timeBinMax = newTimeBinMax; }
   void AddDigit(AbstractTpcDigit *newDigit) { digits.push_back(newDigit); }

protected:
private:
   int                             clusterID;
   int                             sector;
   int                             row;
   int                             padMin;
   int                             padMax;
   int                             timeBinMin;
   int                             timeBinMax;
   std::vector<AbstractTpcDigit *> digits;

   ClassDef(Tpc2dClusterFast, 1);
};

#endif