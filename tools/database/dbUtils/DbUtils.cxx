#include "DbUtils.h"
#include <nlohmann/json.hpp>
#include "UniDbDetectorParameter.h"
#include "UniDbParameter.h"
#include "UniDbDetector.h"
#include <TString.h>
using json = nlohmann::json;
#define DB_UTILS_fNumOfSectors 24
#define DB_UTILS_fNumZLayers 4
#define DB_UTILS_fNumYLayers 1
static const std::string TPC_VELOCITY_DETECTOR_NAME = "TPC";
static const std::string TPC_VELOCITY_PARAMETR_NAME = "TPC_VELOCITY";
static const std::string TPC_CONST_PARAMETR_NAME    = "TPC_CONST";
bool DbUtils::getTpcVelocityMap(size_t periodID, size_t runID, std::vector<std::vector<std::vector<double>>> &vecRes)
{
   // define some const
   //    const std::string detectorName("TEST");
   //    const std::string paramName("tpc_velocity_json");
   const Int_t fNumOfSectors = DB_UTILS_fNumOfSectors;
   const Int_t fNumZLayers   = DB_UTILS_fNumZLayers;
   const Int_t fNumYLayers   = DB_UTILS_fNumYLayers;

   std::vector<std::vector<std::vector<double>>> velMap;
   velMap.resize(fNumOfSectors, std::vector<std::vector<double>>(fNumZLayers, std::vector<double>(fNumYLayers, 0.)));
   UniDbDetectorParameter *pDetectorParameter = UniDbDetectorParameter::GetDetectorParameter(
      TPC_VELOCITY_DETECTOR_NAME, TPC_VELOCITY_PARAMETR_NAME, periodID, runID);
   if (pDetectorParameter == nullptr) {
      cout << "\nMacro finished with errors" << endl;
      return false;
   }
   JSONValue *vv = (JSONValue *)pDetectorParameter->GetValue();
   json       j_out;

   if (vv->getJSON(j_out)) {
      json arr = j_out["array"];
      for (int i1 = 0; i1 < arr.size(); i1++) {
         json                             arr1 = arr.at(i1);
         std::vector<std::vector<double>> map2;
         for (int i2 = 0; i2 < arr1.size(); i2++) {
            json                arr2 = arr1.at(i2);
            std::vector<double> map3;
            // std::vector<double> map3 = map2[i2];
            for (int i3 = 0; i3 < arr2.size(); i3++) {
               double val = arr2.at(i3);
               map3.push_back(val);
            }
            map2.push_back(map3);
         }
         vecRes.push_back(map2);
      }
   } else {
      cout << "\nMacro finished with errors" << endl;
      return false;
   }
   /*    BoolValue* is_on = (BoolValue*) pDetectorParameter->GetValue();
       if (is_on->value)
           cout<<"Detector DCH1 was turned on in run n.77"<<endl;
       else
           cout<<"Detector DCH1 was turned off in run n.77"<<endl;
    */
   return true;
}

int DbUtils::setTpcVelocityMap(size_t periodID, size_t runID, std::vector<std::vector<std::vector<double>>> &map)
{
   // define some const
   const Int_t fNumOfSectors = DB_UTILS_fNumOfSectors;
   const Int_t fNumZLayers   = DB_UTILS_fNumZLayers;
   const Int_t fNumYLayers   = DB_UTILS_fNumYLayers;
   if (UniDbDetector::CheckDetectorExists(TPC_VELOCITY_DETECTOR_NAME) != 1) {
      TString name = "Unknown";
      UniDbDetector::CreateDetector(TPC_VELOCITY_DETECTOR_NAME, &name);
   }
   if (UniDbParameter::CheckParameterExists(TPC_VELOCITY_PARAMETR_NAME) != 1)
      UniDbParameter::CreateParameter(TPC_VELOCITY_PARAMETR_NAME, 17, true);

   json jsOut       = {}; // = R"({"compact": true, "schema": 2})"_json;
   jsOut["version"] = 1.0;
   json arr(json::value_t::array);
   for (int i1 = 0; i1 < map.size(); i1++) {
      json                             arr1(json::value_t::array);
      std::vector<std::vector<double>> map2 = map[i1];
      for (int i2 = 0; i2 < map2.size(); i2++) {
         json                arr2(json::value_t::array);
         std::vector<double> map3 = map2[i2];
         for (int i3 = 0; i3 < map3.size(); i3++) {
            double val = map3[i3];
            arr2.push_back(val);
         }
         arr1.push_back(arr2);
      }
      arr.push_back(arr1);
   }
   jsOut["array"] = arr;
   JSONValue               bValue(jsOut);
   UniDbDetectorParameter *pDetectorParameter = UniDbDetectorParameter::CreateDetectorParameter(
      TPC_VELOCITY_DETECTOR_NAME, TPC_VELOCITY_PARAMETR_NAME, periodID, runID, periodID, runID, &bValue);
   return -1;
}

json DbUtils::getTpcConsts(size_t periodID, size_t runID)
{
   UniDbDetectorParameter *pDetectorParameter = UniDbDetectorParameter::GetDetectorParameter(
      TPC_VELOCITY_DETECTOR_NAME, TPC_CONST_PARAMETR_NAME, periodID, runID);
   if (pDetectorParameter == nullptr) {
      cout << "\nMacro finished with errors" << endl;
      return false;
   }
   JSONValue *vv = (JSONValue *)pDetectorParameter->GetValue();
   json       j_out;

   if (vv->getJSON(j_out)) {
   }
   return j_out;
}
bool DbUtils::setTpcConsts(int sector_count, double sector_phi_rad, double z_min, double drift_length,
                           int timebin_count, double timebin_length, std::vector<int> &row_count,
                           std::vector<double> &pad_height, std::vector<double> &pad_width, std::vector<int> &pad_count,
                           double ypadplane_offset, double ypadplane2padarea_offset, size_t periodID, size_t runID)
{
   if (UniDbDetector::CheckDetectorExists(TPC_VELOCITY_DETECTOR_NAME) != 1) {
      TString name = "Unknown";
      UniDbDetector::CreateDetector(TPC_VELOCITY_DETECTOR_NAME, &name);
   }
   if (UniDbParameter::CheckParameterExists(TPC_CONST_PARAMETR_NAME) != 1)
      UniDbParameter::CreateParameter(TPC_CONST_PARAMETR_NAME, 17, true);
   json jsOut              = {};
   jsOut["version"]        = 1.0;
   jsOut["SECTOR_COUNT"]   = sector_count;
   jsOut["SECTOR_PHI_RAD"] = sector_phi_rad;
   jsOut["Z_MIN"]          = z_min;
   jsOut["DRIFT_LENGTH"]   = drift_length;
   jsOut["TIMEBIN_COUNT"]  = timebin_count;
   jsOut["TIMEBIN_LENGTH"] = timebin_length;
   json         arrRCount(json::value_t::array);
   unsigned int i;
   for (i = 0; i < row_count.size(); i++) {
      int val = row_count[i];
      arrRCount.push_back(val);
   }
   jsOut["ROW_COUNT"] = arrRCount;
   json arrPadH(json::value_t::array);
   for (i = 0; i < pad_height.size(); i++) {
      double val = pad_height[i];
      arrPadH.push_back(val);
   }
   jsOut["PAD_HEIGHT"] = arrPadH;
   json arrPadW(json::value_t::array);
   for (i = 0; i < pad_width.size(); i++) {
      double val = pad_width[i];
      arrPadW.push_back(val);
   }
   jsOut["PAD_WIDTH"] = arrPadW;
   json arrPadC(json::value_t::array);
   for (i = 0; i < pad_count.size(); i++) {
      int val = pad_count[i];
      arrPadC.push_back(val);
   }
   jsOut["PAD_COUNT"]                = arrPadC;
   jsOut["YPADPLANE_OFFSET"]         = ypadplane_offset;
   jsOut["YPADPLANE2PADAREA_OFFSET"] = ypadplane2padarea_offset;
   JSONValue               bValue(jsOut);
   UniDbDetectorParameter *pDetectorParameter = UniDbDetectorParameter::CreateDetectorParameter(
      TPC_VELOCITY_DETECTOR_NAME, TPC_CONST_PARAMETR_NAME, periodID, runID, periodID, runID, &bValue);
   return true;
}

bool DbUtils::parceIntVectorFromJson(std::string jsonDump, std::vector<int> &vec, const char *fieldName)
{
   bool res = false;
   json js;
   try {
      js = json::parse(jsonDump);
   } catch (json::parse_error &e) {
      std::cerr << "message: " << e.what() << '\n'
                << "exception id: " << e.id << '\n'
                << "byte position of error: " << e.byte << std::endl;
      return res;
   }
   json arr;
   if (fieldName != NULL) {
      if (!js.contains(fieldName)) {
         std::cerr << "DbUtils::parceIntVectorFromJson json has not field " << fieldName << std::endl;
         return res;
      }
      arr = js[fieldName];
   } else
      arr = js;

   if (!arr.is_array()) {
      std::cerr << "DbUtils::parceIntVectorFromJson field " << fieldName << " is not array" << std::endl;
      return res;
   }
   for (int i = 0; i < arr.size(); i++) {
      int val = arr[i];
      vec.push_back(val);
   }
   if (vec.size() > 0) {
      res = true;
   }
   return res;
}
bool DbUtils::parceDoubleVectorFromJson(std::string jsonDump, std::vector<double> &vec, const char *fieldName)
{
   bool res = false;
   json js;
   try {
      js = json::parse(jsonDump);
   } catch (json::parse_error &e) {
      std::cerr << "message: " << e.what() << '\n'
                << "exception id: " << e.id << '\n'
                << "byte position of error: " << e.byte << std::endl;
      return false;
   }
   json arr;
   if (fieldName != NULL) {
      if (!js.contains(fieldName)) {
         std::cerr << "DbUtils::parceDoubleVectorFromJson json has not field " << fieldName << std::endl;
         return res;
      }
      arr = js[fieldName];
   } else
      arr = js;

   if (!arr.is_array()) {
      std::cerr << "DbUtils::parceDoubleVectorFromJson field " << fieldName << " is not array" << std::endl;
      return res;
   }
   for (int i = 0; i < arr.size(); i++) {
      double val = arr[i];
      vec.push_back(val);
   }
   if (vec.size() > 0) {
      res = true;
   }
   return res;
}

// ClassImp(DBUtils)
