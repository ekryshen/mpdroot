#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <sstream>
#include <map>

#include "MpdPhotonAnalysisParams.h"

using namespace std;

ClassImp(MpdPhotonAnalysisParams);

void MpdPhotonAnalysisParams::ReadFromFile(std::string fname)
{
   // Read from file in format
   //  varname  value
   //  comments start from #

   // V  if (fname.size() == 0) {
   // V    return;
   // V }

   // V  std::ifstream ifs(fname);

   std::string fInputFileTtx = fname;
   fInputFileTtx             = fInputFileTtx + ".txt";
   cout << "Read input text file: " << fInputFileTtx << endl;

   std::ifstream ifs(fInputFileTtx);
   if (!ifs) {
      cout << "File " << fInputFileTtx << " can not be opened -> using input parameters from the header file instead"
           << endl;
      return;
   }

   std::string line, a, b;
   while (getline(ifs, line)) {
      stringstream(line) >> a >> b;
      // if comment, skip the line
      if (a.at(0) != '#') mMap.insert({a, b});
   }
   ifs.close();

   // Parse prepared map
   read("mApplySelection", mApplySelection);
   read("mZvtxCut", mZvtxCut);
   read("mNhitsCut", mNhitsCut);
   // V0 cuts
   read("mUseBDT", mUseBDT);
   if (mUseBDT) {
      read("mUseBDTRegP", mUseBDTRegP);
      read("mBDTCut", mBDTCut);
   } else {
      read("mMinR2Cut", mMinR2Cut);
      read("mMaxR2Cut", mMaxR2Cut);
      read("mPIDsigM", mPIDsigM);
      read("mPIDsigE", mPIDsigE);
      read("mPIDenergy", mPIDenergy);
      read("mPIDkoeff", mPIDkoeff);
      read("mPIDgenerator", mPIDgenerator);
      read("mPIDtracking", mPIDtracking);
      read("mPIDparticles", mPIDparticles);
      read("mNofHitsCut", mNofHitsCut);
      read("mEtaCut", mEtaCut);
      read("mPtminCut", mPtminCut);
      read("mProbElCut", mProbElCut);
      read("mdEdxSigmaCut", mdEdxSigmaCut);
      read("mBetaSigmaCut", mBetaSigmaCut);
      read("mRequireTOFpid", mRequireTOFpid);
      read("mMassCut", mMassCut);
      read("mDistCut", mDistCut);
      read("mCosPsiCut", mCosPsiCut);
      read("mCPACut", mCPACut);
      read("mChi2Cut", mChi2Cut);
   }

   read("mCluEmin", mCluEmin);
   read("mCluMult", mCluMult);
   read("mCluTofMin", mCluTofMin);
   read("mCluTofMax", mCluTofMax);
   read("mCluDisp", mCluDisp);
   read("mCluDispEmin", mCluDispEmin);
   read("mCluCPV", mCluCPV);
}
void MpdPhotonAnalysisParams::Print() const
{
   cout << "#-------Parameters used for analysis------" << endl;
   cout << "# Event selection: " << endl;
   cout << "mApplySelection " << mApplySelection << " // apply event, track, V0 and EMC cluster selection" << endl;
   cout << "mZvtxCut " << mZvtxCut << " // cut on vertex z coordinate" << endl;
   cout << "mNhitsCut " << mNhitsCut << " //  number of hits in TPC tracks used for centrality" << endl;

   cout << "# V0 cuts: " << endl;
   if (!mUseBDT) {
      cout << "mMinR2Cut " << mMinR2Cut << " // (cm) Minimal conversion radius (to exclude Dalitz)" << endl;
      cout << "mMaxR2Cut " << mMaxR2Cut << " // (cm) Maximal conversion radius (to exclude poorly reconstructed tracks)"
           << endl;
      cout << "mPIDsigM   " << mPIDsigM << "  // dEdx PID parameters" << endl;
      cout << "mPIDsigE   " << mPIDsigE << "  // dEdx PID parameters" << endl;
      cout << "mPIDenergy " << mPIDenergy << "  // dEdx PID parameters" << endl;
      cout << "mPIDkoeff  " << mPIDkoeff << "  // dEdx PID parameters" << endl;

      cout << "mPIDgenerator " << mPIDgenerator << "  // dEdx PID parameters" << endl;
      cout << "mPIDtracking  " << mPIDtracking << "  // dEdx PID parameters" << endl;
      cout << "mPIDparticles " << mPIDparticles << "  // dEdx PID parameters" << endl;

      cout << "mNofHitsCut  " << mNofHitsCut << "  // minimal number of hits to accept track" << endl;
      cout << "mEtaCut      " << mEtaCut << "  // maximal pseudorapidity accepted" << endl;
      cout << "mPtminCut   " << mPtminCut << "  // minimal pt used in analysis" << endl;
      cout << "mNofHitsCut " << mNofHitsCut << "  // minimal number of hits to accept track" << endl;
      cout << "mProbElCut  " << mProbElCut << "  // minimal dEdx probability for electrons" << endl;
      cout << "mdEdxSigmaCut  " << mdEdxSigmaCut << "  // dEdx cut in sigmas" << endl;
      cout << "mBetaSigmaCut  " << mBetaSigmaCut << "  // beta cut" << endl;
      cout << "mRequireTOFpid  " << mRequireTOFpid << "  // mRequireTOFpid" << endl;
      cout << "mCPACut  " << mCPACut << "  // cos(PA)" << endl;
      cout << "mMassCut  " << mMassCut << "  // e+e- pair mass cut" << endl;
      cout << "mDistCut  " << mDistCut << "  // maximal closest distance between daughters" << endl;
      cout << "mCosPsiCut " << mCosPsiCut << "  // e+e- pair orientation wrt B-filed" << endl;
      cout << "mChi2Cut  " << mChi2Cut << "  // maximal chi2 in Kalman fit" << endl;
   } else {
      cout << "Using BDT selection with threshold " << mBDTCut << endl;
      if (mUseBDTRegP) {
         cout << "Using V0 BDT momemtum correction" << endl;
      }
   }

   cout << "# Cluster cuts: " << endl;

   cout << "mCluEmin  " << mCluEmin << "  // (GeV) minimal cluster energy" << endl;
   cout << "mCluMult  " << mCluMult << "  // minimal number of cells in cluster" << endl;
   cout << "mCluTofMin   " << mCluTofMin << "  // minimal time wrt photon arrival in sigmas" << endl;
   cout << "mCluTofMax   " << mCluTofMax << "  // maximal time wrt photon arrival in sigmas" << endl;
   cout << "mCluDisp  " << mCluDisp << "  // disp cut" << endl;
   cout << "mCluDispEmin  " << mCluDispEmin << "  // Emin for disp cut" << endl;
   cout << "mCluCPV   " << mCluCPV << "  // (sigma) minimal distance to charged track extrapolation" << endl;
   cout << "------------------------------------------" << endl;
}

void MpdPhotonAnalysisParams::read(std::string name, bool &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      if (search->second.compare("true") == 0 || search->second.compare("TRUE") == 0 ||
          search->second.compare("1") == 0) {
         b = true;
      } else {
         b = false;
      }
   }
}
void MpdPhotonAnalysisParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdPhotonAnalysisParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdPhotonAnalysisParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   }
}
