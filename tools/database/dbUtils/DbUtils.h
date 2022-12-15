#ifndef DBUTILS_H
#define DBUTILS_H

#include "Rtypes.h"
#include <TString.h>
#include <vector>
#include <string>
#include <utility>

class DbUtils {
public:
   // DbUtils(/*some coonection params*/);
   // ~DbUtils();

public:
   // std::vector<std::vector<std::vector<double>>> GetVelocityMap(size_t periodID,size_t runID);
   static bool getTpcVelocityMap(size_t periodID, size_t runID, std::vector<std::vector<std::vector<double>>> &vecRes);
   static int  setTpcVelocityMap(size_t periodID, size_t runID,
                                 std::vector<std::vector<std::vector<double>>> &map); // return ok status

private:
   DbUtils() { ; };
   // const Int_t fNumOfSectors = 24;
   // const Int_t fNumZLayers = 4;
   // const Int_t fNumYLayers = 1;
   // TString paramName;
   // TString detectorName;
   // ClassDef(DbUtils, 0) // must be a last line in class
};

#endif
