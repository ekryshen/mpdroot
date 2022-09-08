#include "MpdTpcDriftVelocityDB.h"

MpdTpcDriftVelocityDB::MpdTpcDriftVelocityDB(/*some coonection params*/)
{
    
}

MpdTpcDriftVelocityDB::~MpdTpcDriftVelocityDB()
{
    
}
	
std::vector<std::vector<std::vector<double>>> MpdTpcDriftVelocityDB::GetVelocityMap(size_t runID)
{
    std::vector<std::vector<std::vector<double>>> velMap;
    velMap.resize(fNumOfSectors, std::vector<std::vector<double>>(fNumZLayers, std::vector<double>(fNumYLayers, 0.)));
    
    return velMap;
}

int MpdTpcDriftVelocityDB::GetVelocityMap(size_t runID, std::vector<std::vector<std::vector<double>>> & map)
{
    return -1;
}

ClassImp(MpdTpcDriftVelocityDB)

