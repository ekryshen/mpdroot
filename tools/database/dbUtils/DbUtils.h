#ifndef DBUTILS_H
#define DBUTILS_H

#include "Rtypes.h"
#include <TString.h>
#include <vector>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>

class DbUtils {

public:
   // TPC Velocity
   static bool getTpcVelocityMap(size_t periodID, size_t runID, std::vector<std::vector<std::vector<double>>> &vecRes);
   static int  setTpcVelocityMap(size_t periodID, size_t runID,
                                 std::vector<std::vector<std::vector<double>>> &map); // return ok status
   // TPC Const
   static bool           setTpcConsts(int sector_count, double sector_phi_rad, double z_min, double drift_length,
                                      int timebin_count, double timebin_length, std::vector<int> &row_count,
                                      std::vector<double> &pad_height, std::vector<double> &pad_width,
                                      std::vector<int> &pad_count, double ypadplane_offset, double ypadplane2padarea_offset,
                                      size_t periodID = 0, size_t runID = 0);
   static nlohmann::json getTpcConsts(size_t periodID = 0, size_t runID = 0);

   // common
   static bool parceIntVectorFromJson(std::string jsonDump, std::vector<int> &vec, const char *fieldName = NULL);
   static bool parceDoubleVectorFromJson(std::string jsonDump, std::vector<double> &vec, const char *fieldName = NULL);

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
