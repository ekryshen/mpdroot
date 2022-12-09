//-------------------------------------------------------------------------------------------------
// Description:
//      Test-driven interface (Abstract base class) for 2D Clusters in TPC padarea
//      dimension 1 - padrow; dimension 2 - timebin
//      Descriptive drawing:
//      https://git.jinr.ru/nica/docs/-/blob/main/docs/mpdroot/coding/clusterization/TPC_2dCluster.svg
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, December, 2022
//-------------------------------------------------------------------------------------------------

#ifndef ABSTRACTTPC2DCLUSTER_HH
#define ABSTRACTTPC2DCLUSTER_HH

#include <TObject.h>

#include <vector>
#include <memory>

#include "AbstractTpcDigit.h"

class AbstractTpc2dCluster : public TObject {
public:
   // Constructors/Destructors ---------
   AbstractTpc2dCluster() {}
   virtual ~AbstractTpc2dCluster() {}

   virtual int                                            GetClusterID() const     = 0;
   virtual int                                            GetSector() const        = 0;
   virtual int                                            GetRow() const           = 0;
   virtual int                                            GetPadMin() const        = 0;
   virtual int                                            GetPadMax() const        = 0;
   virtual int                                            GetTimeBinMin() const    = 0;
   virtual int                                            GetTimeBinMax() const    = 0;
   virtual std::vector<std::shared_ptr<AbstractTpcDigit>> GetClusterDigits() const = 0;

protected:
private:
   ClassDef(AbstractTpc2dCluster, 1);
};

#endif