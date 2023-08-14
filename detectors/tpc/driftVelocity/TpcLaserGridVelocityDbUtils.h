#ifndef TPCLASERGRIDVELOCITYDBUTILS_H
#define TPCLASERGRIDVELOCITYDBUTILS_H

#include "Rtypes.h"
#include <TString.h>
#include <vector>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>

class TpcLaserGridVelocityDbUtils {

public:
   static bool getTpcVelocityMap(size_t periodID, size_t runID, std::vector<std::vector<std::vector<double>>> &vecRes);
   static int  setTpcVelocityMap(size_t periodID, size_t runID,
                                 std::vector<std::vector<std::vector<double>>> &map); // return ok status
private:
   TpcLaserGridVelocityDbUtils() { ; };

   // ClassDef(TpcLaserGridVelocityDbUtils, 0) // must be a last line in class
};

#endif
