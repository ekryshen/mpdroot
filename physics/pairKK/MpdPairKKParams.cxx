#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdPairKKParams.h"

using namespace std;

ClassImp(MpdPairKKParams);

void MpdPairKKParams::ReadFromFile(std::string fname)
{
   // Read from file in format
   //  varname  value
   //  comments start from #

   std::string fInputFileTtx = fname;
   fInputFileTtx             = fInputFileTtx + ".txt";

   cout << "Read input text file: " << fInputFileTtx << endl;

   /*
      if (fname.size() == 0) {
         return;
      }
   */

   std::ifstream ifs(fInputFileTtx);
   if (!ifs) {
      cout << "File " << fInputFileTtx << " can not be opened -> using input parameters from the header file instead"
           << endl;
      return;
   }

   std::string a, b;
   while (ifs.good()) {
      ifs >> a;
      // if comment, skip to the enf of line
      if (a.find_first_of('#') == 0) {
         ifs.ignore(999, '\n');
         continue;
      } else {
         ifs >> b;
         mMap.insert({a, b});
      }
   }
   ifs.close();

   // Parse prepared map
   read("mZvtxCut", mZvtxCut);
   // V0 cuts
   read("mPIDsigTPC", mPIDsigTPC);
   read("mPIDsigTOF", mPIDsigTOF);
   read("mNofHitsCut", mNofHitsCut);
   read("mEtaCut", mEtaCut);
   read("mPtminCut", mPtminCut);
}

void MpdPairKKParams::Print() const
{
   cout << "#-------Parameters used for PairKK analysis------" << endl;
   cout << "# Event selection: " << endl;
   cout << "mZvtxCut      " << mZvtxCut << " // cut on vertex z coordinate" << endl;

   cout << "# PID cuts:   " << endl;
   cout << "mPIDsigTPC:   " << mPIDsigTPC << "  // dEdx PID parameters" << endl;
   cout << "mPIDsigTOF:   " << mPIDsigTOF << "  // Beta PID parameters" << endl;

   cout << "mNofHitsCut:  " << mNofHitsCut << "  // minimal number of hits to accept track" << endl;
   cout << "mEtaCut:      " << mEtaCut << "  // maximal pseudorapidity accepted" << endl;
   cout << "mPtminCut:    " << mPtminCut << "  // minimal pt used in analysis" << endl;
   cout << "------------------------------------------" << endl << endl;
}

void MpdPairKKParams::read(std::string name, bool &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      if (search->second.compare("true") == 0 || search->second.compare("TRUE") == 0) {
         b = true;
      } else {
         b = false;
      }
   }
}
void MpdPairKKParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdPairKKParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdPairKKParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   }
}
