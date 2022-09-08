#ifndef MPDTPCDRIFTVELOCITYDB_H
#define MPDTPCDRIFTVELOCITYDB_H

#include "Rtypes.h"

#include <vector>
#include <string>
#include <utility>

class MpdTpcDriftVelocityDB
{
public:
    MpdTpcDriftVelocityDB(/*some coonection params*/);
    virtual ~MpdTpcDriftVelocityDB();
    
public:
    virtual std::vector<std::vector<std::vector<double>>> GetVelocityMap(size_t runID);
    virtual int GetVelocityMap(size_t runID, std::vector<std::vector<std::vector<double>>> & map); //return ok status
    
private:
    const Int_t fNumOfSectors = 24;
    const Int_t fNumZLayers = 4;
    const Int_t fNumYLayers = 1;
    
    ClassDef(MpdTpcDriftVelocityDB, 0) // must be a last line in class
};

#endif

