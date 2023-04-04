#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdCentralityAllParams.h"

using namespace std;

ClassImp(MpdCentralityAllParams);

void MpdCentralityAllParams::ReadFromFile(std::string fname)
{
   // Read from file in format
   //  varname  value
   //  comments start from #

   std::string fInputFileTtx = fname;
   fInputFileTtx             = fInputFileTtx + ".txt";

   cout << "Read input text file: " << fInputFileTtx << endl;

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
   read("mNofHitsCut", mNofHitsCut);
   read("mEtaCut", mEtaCut);
   read("mPtminCut", mPtminCut);
   read("mDcaCut", mDcaCut);
   read("mProdGenerator", mProdGenerator);
   read("mInFileConvert", mInFileConvert);
   read("mInFileTrEff", mInFileTrEff);
}

void MpdCentralityAllParams::Print() const
{
   cout << "#-------Parameters used for evCentrality analysis------" << endl;
   cout << "# Event selection: " << endl;
   cout << "mZvtxCut       " << mZvtxCut << " // cut on vertex z coordinate" << endl;

   cout << "# Track cuts:  " << endl;
   cout << "mNofHitsCut    " << mNofHitsCut << "  // minimal number of hits to accept track" << endl;
   cout << "mEtaCut        " << mEtaCut << "  // maximal pseudorapidity accepted" << endl;
   cout << "mPtminCut      " << mPtminCut << "  // minimal pt used in analysis" << endl;
   cout << "mDcaCut        " << mDcaCut << "  // maximal DCA accepted" << endl;

   cout << "# Production selection: " << endl;
   cout << "mProdGenerator " << mProdGenerator << "  // Production-Generator" << endl;
   cout << "mInFileConvert " << mInFileConvert << "  // Input file with track-to-centrality converter" << endl;

   cout << "# Track efficiecny corrections: " << endl;
   cout << "mInFileTrEff   " << mInFileTrEff << "  // Input file with track efficiecny corrections" << endl;

   cout << "------------------------------------------" << endl;
}

void MpdCentralityAllParams::read(std::string name, bool &b)
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
void MpdCentralityAllParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdCentralityAllParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdCentralityAllParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   }
}
