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
/*DbUtils::DbUtils():detectorName("TEST"),paramName("tpc_velocity_json")
{
    //json j = R"({"compact": true, "schema": 2})"_json;
    //TString paramName = "tpc_velocity_json";
    if(UniDbParameter::CheckParameterExists(paramName)!=1)
        UniDbParameter::CreateParameter(paramName,17,true);

}

DbUtils::~DbUtils()
{

}
*/
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

// ClassImp(DBUtils)
