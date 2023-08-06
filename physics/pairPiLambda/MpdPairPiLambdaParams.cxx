#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdPairPiLambdaParams.h"

using namespace std;

ClassImp(MpdPairPiLambdaParams);

void MpdPairPiLambdaParams::ReadFromFile(std::string fname)
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
   // Lambda cuts
   read("mChi2PionLam", mChi2PionLam);
   read("mChi2ProtLam", mChi2ProtLam);
   read("mLamEtaCut", mLamEtaCut);
   read("mChi2Lam", mChi2Lam);
   read("mPALam", mPALam);
   read("mDecayLam", mDecayLam);
   read("mDistLam", mDistLam);
   read("mNSigmaLam", mNSigmaLam);
   read("mWidthLam", mWidthLam);
}

void MpdPairPiLambdaParams::Print() const
{
   cout << "#-------Parameters used for PairPiLambda analysis------" << endl;
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

   cout << "# Lambda cuts:  " << endl;
   cout << "mChi2PionLam:   " << mChi2PionLam << "  // minimum Chi2-to-PV for pion from Lambda" << endl;
   cout << "mChi2ProtLam:   " << mChi2ProtLam << "  // minimum Chi2-to-PV for proton from Lambda" << endl;
   cout << "mLamEtaCut:     " << mLamEtaCut << "  // maximum pseudorapidity for Lambda" << endl;
   cout << "mChi2Lam:       " << mChi2Lam << "  // maximum Chi2 for Lambda secondary vertex" << endl;
   cout << "mPALam:         " << mPALam << "  // maximum pointing angle for Lambda" << endl;
   cout << "mDecayLam:      " << mDecayLam << "  // minimum decay distance for Lambda" << endl;
   cout << "mDistLam:       " << mDistLam << "  // maximum distance between p-pi in the secondary vertex from Lambda"
        << endl;
   cout << "mNSigmaLam:     " << mNSigmaLam << "  // n-sigma selection for Lambda candidates" << endl;
   cout << "mWidthLam:      " << mWidthLam << "  // pT-averaged width of Lambda peak" << endl;

   cout << "------------------------------------------" << endl;
}

void MpdPairPiLambdaParams::read(std::string name, bool &b)
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
void MpdPairPiLambdaParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdPairPiLambdaParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdPairPiLambdaParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   }
}
