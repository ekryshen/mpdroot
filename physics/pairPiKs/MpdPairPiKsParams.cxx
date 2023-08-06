#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdPairPiKsParams.h"

using namespace std;

ClassImp(MpdPairPiKsParams);

void MpdPairPiKsParams::ReadFromFile(std::string fname)
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
      if (a.find_first_of('#') == 0 || a.find_first_of("//") == 0) {
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
   // Track cuts
   read("mPIDsigTPC", mPIDsigTPC);
   read("mPIDsigTOF", mPIDsigTOF);
   read("mNofHitsCut", mNofHitsCut);
   read("mEtaCut", mEtaCut);
   read("mPtminCut", mPtminCut);
   read("mDCACut", mDCACut);
   read("mYCut", mYCut);
   // Ks cuts
   read("mChi2PionKs", mChi2PionKs);
   read("mKsEtaCut", mKsEtaCut);
   read("mChi2Ks", mChi2Ks);
   read("mPAKs", mPAKs);
   read("mDecayKs", mDecayKs);
   read("mDistKs", mDistKs);
   read("mNSigmaKs", mNSigmaKs);
   read("mWidthKs", mWidthKs);
}

void MpdPairPiKsParams::Print() const
{
   cout << "#-------Parameters used for PairPiKs analysis------" << endl;
   cout << "# Event selection: " << endl;
   cout << "mZvtxCut      " << mZvtxCut << " // cut on vertex z coordinate" << endl;

   cout << "# PID cuts:   " << endl;
   cout << "mPIDsigTPC:   " << mPIDsigTPC << "  // dEdx PID parameters" << endl;
   cout << "mPIDsigTOF:   " << mPIDsigTOF << "  // Beta PID parameters" << endl;

   cout << "# Track cuts:   " << endl;
   cout << "mNofHitsCut:  " << mNofHitsCut << "  // minimal number for a track" << endl;
   cout << "mEtaCut:      " << mEtaCut << "  // maximal pseudorapidity for a track" << endl;
   cout << "mPtminCut:    " << mPtminCut << "  // minimal pt for a track" << endl;
   cout << "mDCACut:      " << mDCACut << "  // maximum DCA for a track" << endl;

   cout << "# Pair cuts:  " << endl;
   cout << "mYCut:        " << mYCut << "  // pair cut" << endl;

   cout << "# Ks cuts:  " << endl;
   cout << "mChi2PionKs:   " << mChi2PionKs << "  // minimum Chi2-to-PV for pion from Ks" << endl;
   cout << "mKsEtaCut:     " << mKsEtaCut << "  // maximum pseudorapidity for Ks" << endl;
   cout << "mChi2Ks:       " << mChi2Ks << "  // maximum Chi2 for Ks secondary vertex" << endl;
   cout << "mPAKs:         " << mPAKs << "  // maximum pointing angle for Ks" << endl;
   cout << "mDecayKs:      " << mDecayKs << "  // minimum decay distance for Ks" << endl;
   cout << "mDistKs:       " << mDistKs << "  // maximum distance between p-pi in the secondary vertex from Ks" << endl;
   cout << "mNSigmaKs:     " << mNSigmaKs << "  // n-sigma selection for Ks candidates" << endl;
   cout << "mWidthKs:      " << mWidthKs << "  // pT-averaged width of Ks peak" << endl;

   cout << "------------------------------------------" << endl;
}

void MpdPairPiKsParams::read(std::string name, bool &b)
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
void MpdPairPiKsParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdPairPiKsParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdPairPiKsParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   }
}
