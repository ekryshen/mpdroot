#include <iostream>
#include <fstream> // std::ifstream

#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdZdcDigi.h"
#include "TRandom.h"

#include "MpdEventPlaneAll.h"
#include "TFile.h"

ClassImp(MpdEventPlaneAll);

MpdEventPlaneAll::MpdEventPlaneAll(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

MpdEventPlaneAll::~MpdEventPlaneAll()
{
   // mCorrQxFHCalFAll.clear();
   // mCorrQyFHCalFAll.clear();
   // mCorrQxFHCalNAll.clear();
   // mCorrQyFHCalNAll.clear();
   // mCorrQxFHCalSAll.clear();
   // mCorrQyFHCalSAll.clear();
   // mCorrQxTPCNAll.clear();
   // mCorrQyTPCNAll.clear();
   // mCorrQxTPCSAll.clear();
   // mCorrQyTPCSAll.clear();

   // mCorrCos1FHCalFAll.clear();
   // mCorrCos2FHCalFAll.clear();
   // mCorrCos3FHCalFAll.clear();
   // mCorrCos4FHCalFAll.clear();
   // mCorrCos5FHCalFAll.clear();
   // mCorrCos6FHCalFAll.clear();
   // mCorrCos7FHCalFAll.clear();
   // mCorrCos8FHCalFAll.clear();
   // mCorrSin1FHCalFAll.clear();
   // mCorrSin2FHCalFAll.clear();
   // mCorrSin3FHCalFAll.clear();
   // mCorrSin4FHCalFAll.clear();
   // mCorrSin5FHCalFAll.clear();
   // mCorrSin6FHCalFAll.clear();
   // mCorrSin7FHCalFAll.clear();
   // mCorrSin8FHCalFAll.clear();

   // mCorrCos1FHCalNAll.clear();
   // mCorrCos2FHCalNAll.clear();
   // mCorrCos3FHCalNAll.clear();
   // mCorrCos4FHCalNAll.clear();
   // mCorrCos5FHCalNAll.clear();
   // mCorrCos6FHCalNAll.clear();
   // mCorrCos7FHCalNAll.clear();
   // mCorrCos8FHCalNAll.clear();
   // mCorrSin1FHCalNAll.clear();
   // mCorrSin2FHCalNAll.clear();
   // mCorrSin3FHCalNAll.clear();
   // mCorrSin4FHCalNAll.clear();
   // mCorrSin5FHCalNAll.clear();
   // mCorrSin6FHCalNAll.clear();
   // mCorrSin7FHCalNAll.clear();
   // mCorrSin8FHCalNAll.clear();

   // mCorrCos1FHCalSAll.clear();
   // mCorrCos2FHCalSAll.clear();
   // mCorrCos3FHCalSAll.clear();
   // mCorrCos4FHCalSAll.clear();
   // mCorrCos5FHCalSAll.clear();
   // mCorrCos6FHCalSAll.clear();
   // mCorrCos7FHCalSAll.clear();
   // mCorrCos8FHCalSAll.clear();
   // mCorrSin1FHCalSAll.clear();
   // mCorrSin2FHCalSAll.clear();
   // mCorrSin3FHCalSAll.clear();
   // mCorrSin4FHCalSAll.clear();
   // mCorrSin5FHCalSAll.clear();
   // mCorrSin6FHCalSAll.clear();
   // mCorrSin7FHCalSAll.clear();
   // mCorrSin8FHCalSAll.clear();

   // mCorrCos1TPCNAll.clear();
   // mCorrCos2TPCNAll.clear();
   // mCorrCos3TPCNAll.clear();
   // mCorrCos4TPCNAll.clear();
   // mCorrSin1TPCNAll.clear();
   // mCorrSin2TPCNAll.clear();
   // mCorrSin3TPCNAll.clear();
   // mCorrSin4TPCNAll.clear();

   // mCorrCos1TPCSAll.clear();
   // mCorrCos2TPCSAll.clear();
   // mCorrCos3TPCSAll.clear();
   // mCorrCos4TPCSAll.clear();
   // mCorrSin1TPCSAll.clear();
   // mCorrSin2TPCSAll.clear();
   // mCorrSin3TPCSAll.clear();
   // mCorrSin4TPCSAll.clear();

   // delete mhEvents;
   // delete mhVertex;
   // delete mhHits;
   // delete mhEta;
   // delete mhDca;
   // delete mhPt;
   // delete mhCorrStep;
   // delete mhQxRawFHCalFAll;
   // delete mhQyRawFHCalFAll;
   // delete mhPhiEPRawFHCalFAll;
   // delete mhQxRawFHCalNAll;
   // delete mhQyRawFHCalNAll;
   // delete mhPhiEPRawFHCalNAll;
   // delete mhQxRawFHCalSAll;
   // delete mhQyRawFHCalSAll;
   // delete mhPhiEPRawFHCalSAll;
   // delete mhQxRawTPCNAll;
   // delete mhQyRawTPCNAll;
   // delete mhPhiEPRawTPCNAll;
   // delete mhQxRawTPCSAll;
   // delete mhQyRawTPCSAll;
   // delete mhPhiEPRawTPCSAll;
   // delete mhQxRecFHCalFAll;
   // delete mhQyRecFHCalFAll;
   // delete mhPhiEPRecFHCalFAll;
   // delete mhQxRecFHCalNAll;
   // delete mhQyRecFHCalNAll;
   // delete mhPhiEPRecFHCalNAll;
   // delete mhQxRecFHCalSAll;
   // delete mhQyRecFHCalSAll;
   // delete mhPhiEPRecFHCalSAll;
   // delete mhQxRecTPCNAll;
   // delete mhQyRecTPCNAll;
   // delete mhPhiEPRecTPCNAll;
   // delete mhQxRecTPCSAll;
   // delete mhQyRecTPCSAll;
   // delete mhPhiEPRecTPCSAll;
   // delete mhPhiEPShfFHCalFAll;
   // delete mhPhiEPShfFHCalNAll;
   // delete mhPhiEPShfFHCalSAll;
   // delete mhPhiEPShfTPCNAll;
   // delete mhPhiEPShfTPCSAll;

   // delete mhCorrQxFHCalFAll;
   // delete mhCorrQyFHCalFAll;
   // delete mhCorrQxFHCalNAll;
   // delete mhCorrQyFHCalNAll;
   // delete mhCorrQxFHCalSAll;
   // delete mhCorrQyFHCalSAll;
   // delete mhCorrQxTPCNAll;
   // delete mhCorrQyTPCNAll;
   // delete mhCorrQxTPCSAll;
   // delete mhCorrQyTPCSAll;

   // delete mhCorrCos1FHCalFAll;
   // delete mhCorrCos2FHCalFAll;
   // delete mhCorrCos3FHCalFAll;
   // delete mhCorrCos4FHCalFAll;
   // delete mhCorrCos5FHCalFAll;
   // delete mhCorrCos6FHCalFAll;
   // delete mhCorrCos7FHCalFAll;
   // delete mhCorrCos8FHCalFAll;
   // delete mhCorrSin1FHCalFAll;
   // delete mhCorrSin2FHCalFAll;
   // delete mhCorrSin3FHCalFAll;
   // delete mhCorrSin4FHCalFAll;
   // delete mhCorrSin5FHCalFAll;
   // delete mhCorrSin6FHCalFAll;
   // delete mhCorrSin7FHCalFAll;
   // delete mhCorrSin8FHCalFAll;

   // delete mhCorrCos1FHCalNAll;
   // delete mhCorrCos2FHCalNAll;
   // delete mhCorrCos3FHCalNAll;
   // delete mhCorrCos4FHCalNAll;
   // delete mhCorrCos5FHCalNAll;
   // delete mhCorrCos6FHCalNAll;
   // delete mhCorrCos7FHCalNAll;
   // delete mhCorrCos8FHCalNAll;
   // delete mhCorrSin1FHCalNAll;
   // delete mhCorrSin2FHCalNAll;
   // delete mhCorrSin3FHCalNAll;
   // delete mhCorrSin4FHCalNAll;
   // delete mhCorrSin5FHCalNAll;
   // delete mhCorrSin6FHCalNAll;
   // delete mhCorrSin7FHCalNAll;
   // delete mhCorrSin8FHCalNAll;

   // delete mhCorrCos1FHCalSAll;
   // delete mhCorrCos2FHCalSAll;
   // delete mhCorrCos3FHCalSAll;
   // delete mhCorrCos4FHCalSAll;
   // delete mhCorrCos5FHCalSAll;
   // delete mhCorrCos6FHCalSAll;
   // delete mhCorrCos7FHCalSAll;
   // delete mhCorrCos8FHCalSAll;
   // delete mhCorrSin1FHCalSAll;
   // delete mhCorrSin2FHCalSAll;
   // delete mhCorrSin3FHCalSAll;
   // delete mhCorrSin4FHCalSAll;
   // delete mhCorrSin5FHCalSAll;
   // delete mhCorrSin6FHCalSAll;
   // delete mhCorrSin7FHCalSAll;
   // delete mhCorrSin8FHCalSAll;

   // delete mhCorrCos1TPCNAll;
   // delete mhCorrCos2TPCNAll;
   // delete mhCorrCos3TPCNAll;
   // delete mhCorrCos4TPCNAll;
   // delete mhCorrSin1TPCNAll;
   // delete mhCorrSin2TPCNAll;
   // delete mhCorrSin3TPCNAll;
   // delete mhCorrSin4TPCNAll;

   // delete mhCorrCos1TPCSAll;
   // delete mhCorrCos2TPCSAll;
   // delete mhCorrCos3TPCSAll;
   // delete mhCorrCos4TPCSAll;
   // delete mhCorrSin1TPCSAll;
   // delete mhCorrSin2TPCSAll;
   // delete mhCorrSin3TPCSAll;
   // delete mhCorrSin4TPCSAll;
}

void MpdEventPlaneAll::UserInit()
{

   mParams.ReadFromFile(mParamConfig);
   mParams.Print();

   // Read correction info from the input file with QA and corrections
   TProfile *tmp  = nullptr;
   TH1F     *htmp = nullptr;

   // Pre-setup EP correction maps for reading
   // Get maps for recentering
   mCorrQxFHCalFAll = SetZeroCorr();
   mCorrQyFHCalFAll = SetZeroCorr();
   mCorrQxFHCalNAll = SetZeroCorr();
   mCorrQyFHCalNAll = SetZeroCorr();
   mCorrQxFHCalSAll = SetZeroCorr();
   mCorrQyFHCalSAll = SetZeroCorr();
   mCorrQxTPCNAll   = SetZeroCorr();
   mCorrQyTPCNAll   = SetZeroCorr();
   mCorrQxTPCSAll   = SetZeroCorr();
   mCorrQyTPCSAll   = SetZeroCorr();

   // Get maps for shift
   mCorrCos1FHCalFAll = SetZeroCorr();
   mCorrCos2FHCalFAll = SetZeroCorr();
   mCorrCos3FHCalFAll = SetZeroCorr();
   mCorrCos4FHCalFAll = SetZeroCorr();
   mCorrCos5FHCalFAll = SetZeroCorr();
   mCorrCos6FHCalFAll = SetZeroCorr();
   mCorrCos7FHCalFAll = SetZeroCorr();
   mCorrCos8FHCalFAll = SetZeroCorr();
   mCorrSin1FHCalFAll = SetZeroCorr();
   mCorrSin2FHCalFAll = SetZeroCorr();
   mCorrSin3FHCalFAll = SetZeroCorr();
   mCorrSin4FHCalFAll = SetZeroCorr();
   mCorrSin5FHCalFAll = SetZeroCorr();
   mCorrSin6FHCalFAll = SetZeroCorr();
   mCorrSin7FHCalFAll = SetZeroCorr();
   mCorrSin8FHCalFAll = SetZeroCorr();

   mCorrCos1FHCalNAll = SetZeroCorr();
   mCorrCos2FHCalNAll = SetZeroCorr();
   mCorrCos3FHCalNAll = SetZeroCorr();
   mCorrCos4FHCalNAll = SetZeroCorr();
   mCorrCos5FHCalNAll = SetZeroCorr();
   mCorrCos6FHCalNAll = SetZeroCorr();
   mCorrCos7FHCalNAll = SetZeroCorr();
   mCorrCos8FHCalNAll = SetZeroCorr();
   mCorrSin1FHCalNAll = SetZeroCorr();
   mCorrSin2FHCalNAll = SetZeroCorr();
   mCorrSin3FHCalNAll = SetZeroCorr();
   mCorrSin4FHCalNAll = SetZeroCorr();
   mCorrSin5FHCalNAll = SetZeroCorr();
   mCorrSin6FHCalNAll = SetZeroCorr();
   mCorrSin7FHCalNAll = SetZeroCorr();
   mCorrSin8FHCalNAll = SetZeroCorr();

   mCorrCos1FHCalSAll = SetZeroCorr();
   mCorrCos2FHCalSAll = SetZeroCorr();
   mCorrCos3FHCalSAll = SetZeroCorr();
   mCorrCos4FHCalSAll = SetZeroCorr();
   mCorrCos5FHCalSAll = SetZeroCorr();
   mCorrCos6FHCalSAll = SetZeroCorr();
   mCorrCos7FHCalSAll = SetZeroCorr();
   mCorrCos8FHCalSAll = SetZeroCorr();
   mCorrSin1FHCalSAll = SetZeroCorr();
   mCorrSin2FHCalSAll = SetZeroCorr();
   mCorrSin3FHCalSAll = SetZeroCorr();
   mCorrSin4FHCalSAll = SetZeroCorr();
   mCorrSin5FHCalSAll = SetZeroCorr();
   mCorrSin6FHCalSAll = SetZeroCorr();
   mCorrSin7FHCalSAll = SetZeroCorr();
   mCorrSin8FHCalSAll = SetZeroCorr();

   mCorrCos1TPCNAll = SetZeroCorr();
   mCorrCos2TPCNAll = SetZeroCorr();
   mCorrCos3TPCNAll = SetZeroCorr();
   mCorrCos4TPCNAll = SetZeroCorr();
   mCorrSin1TPCNAll = SetZeroCorr();
   mCorrSin2TPCNAll = SetZeroCorr();
   mCorrSin3TPCNAll = SetZeroCorr();
   mCorrSin4TPCNAll = SetZeroCorr();

   mCorrCos1TPCSAll = SetZeroCorr();
   mCorrCos2TPCSAll = SetZeroCorr();
   mCorrCos3TPCSAll = SetZeroCorr();
   mCorrCos4TPCSAll = SetZeroCorr();
   mCorrSin1TPCSAll = SetZeroCorr();
   mCorrSin2TPCSAll = SetZeroCorr();
   mCorrSin3TPCSAll = SetZeroCorr();
   mCorrSin4TPCSAll = SetZeroCorr();

   if (mParams.mInFileEpCorr != "ANY") {
      cout << "evPlane: Reading out file with EP corrections " << mParams.mInFileEpCorr << endl;
      TFile *inr = new TFile((mParams.mInFileEpCorr).c_str());

      htmp      = (TH1F *)inr->Get("mhCorrStep");
      mCorrStep = 0;
      if (!htmp)
         mCorrStep = 0; // no mhCorrStep histogram was found in the file
      else {
         if (htmp->GetBinContent(1) > 0) mCorrStep = 1; // ready for recentering
         if (htmp->GetBinContent(2) > 0) mCorrStep = 2; // ready for shift
         if (htmp->GetBinContent(3) > 0) mCorrStep = 2; // ready for shift
      }
      delete htmp;

      // Get maps for recentering
      if (mCorrStep == 1 || mCorrStep == 2) {
         tmp              = (TProfile *)inr->Get("mhCorrQxFHCalFAll");
         mCorrQxFHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrQyFHCalFAll");
         mCorrQyFHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrQxFHCalNAll");
         mCorrQxFHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrQyFHCalNAll");
         mCorrQyFHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrQxFHCalSAll");
         mCorrQxFHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrQyFHCalSAll");
         mCorrQyFHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp            = (TProfile *)inr->Get("mhCorrQxTPCNAll");
         mCorrQxTPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp            = (TProfile *)inr->Get("mhCorrQyTPCNAll");
         mCorrQyTPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp            = (TProfile *)inr->Get("mhCorrQxTPCSAll");
         mCorrQxTPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp            = (TProfile *)inr->Get("mhCorrQyTPCSAll");
         mCorrQyTPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
      }

      // Get maps for shift
      if (mCorrStep == 2) {
         tmp                = (TProfile *)inr->Get("mhCorrCos1FHCalFAll");
         mCorrCos1FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos2FHCalFAll");
         mCorrCos2FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos3FHCalFAll");
         mCorrCos3FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos4FHCalFAll");
         mCorrCos4FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos5FHCalFAll");
         mCorrCos5FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos6FHCalFAll");
         mCorrCos6FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos7FHCalFAll");
         mCorrCos7FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos8FHCalFAll");
         mCorrCos8FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin1FHCalFAll");
         mCorrSin1FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin2FHCalFAll");
         mCorrSin2FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin3FHCalFAll");
         mCorrSin3FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin4FHCalFAll");
         mCorrSin4FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin5FHCalFAll");
         mCorrSin5FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin6FHCalFAll");
         mCorrSin6FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin7FHCalFAll");
         mCorrSin7FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin8FHCalFAll");
         mCorrSin8FHCalFAll = ReadEpCorrProfile(tmp);
         delete tmp;

         tmp                = (TProfile *)inr->Get("mhCorrCos1FHCalNAll");
         mCorrCos1FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos2FHCalNAll");
         mCorrCos2FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos3FHCalNAll");
         mCorrCos3FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos4FHCalNAll");
         mCorrCos4FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos5FHCalNAll");
         mCorrCos5FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos6FHCalNAll");
         mCorrCos6FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos7FHCalNAll");
         mCorrCos7FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos8FHCalNAll");
         mCorrCos8FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin1FHCalNAll");
         mCorrSin1FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin2FHCalNAll");
         mCorrSin2FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin3FHCalNAll");
         mCorrSin3FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin4FHCalNAll");
         mCorrSin4FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin5FHCalNAll");
         mCorrSin5FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin6FHCalNAll");
         mCorrSin6FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin7FHCalNAll");
         mCorrSin7FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin8FHCalNAll");
         mCorrSin8FHCalNAll = ReadEpCorrProfile(tmp);
         delete tmp;

         tmp                = (TProfile *)inr->Get("mhCorrCos1FHCalSAll");
         mCorrCos1FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos2FHCalSAll");
         mCorrCos2FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos3FHCalSAll");
         mCorrCos3FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos4FHCalSAll");
         mCorrCos4FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos5FHCalSAll");
         mCorrCos5FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos6FHCalSAll");
         mCorrCos6FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos7FHCalSAll");
         mCorrCos7FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrCos8FHCalSAll");
         mCorrCos8FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin1FHCalSAll");
         mCorrSin1FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin2FHCalSAll");
         mCorrSin2FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin3FHCalSAll");
         mCorrSin3FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin4FHCalSAll");
         mCorrSin4FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin5FHCalSAll");
         mCorrSin5FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin6FHCalSAll");
         mCorrSin6FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin7FHCalSAll");
         mCorrSin7FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp                = (TProfile *)inr->Get("mhCorrSin8FHCalSAll");
         mCorrSin8FHCalSAll = ReadEpCorrProfile(tmp);
         delete tmp;

         tmp              = (TProfile *)inr->Get("mhCorrCos1TPCNAll");
         mCorrCos1TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos2TPCNAll");
         mCorrCos2TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos3TPCNAll");
         mCorrCos3TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos4TPCNAll");
         mCorrCos4TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin1TPCNAll");
         mCorrSin1TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin2TPCNAll");
         mCorrSin2TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin3TPCNAll");
         mCorrSin3TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin4TPCNAll");
         mCorrSin4TPCNAll = ReadEpCorrProfile(tmp);
         delete tmp;

         tmp              = (TProfile *)inr->Get("mhCorrCos1TPCSAll");
         mCorrCos1TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos2TPCSAll");
         mCorrCos2TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos3TPCSAll");
         mCorrCos3TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrCos4TPCSAll");
         mCorrCos4TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin1TPCSAll");
         mCorrSin1TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin2TPCSAll");
         mCorrSin2TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin3TPCSAll");
         mCorrSin3TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
         tmp              = (TProfile *)inr->Get("mhCorrSin4TPCSAll");
         mCorrSin4TPCSAll = ReadEpCorrProfile(tmp);
         delete tmp;
      }

      inr->Close();
   } else {
   }

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // General QA

   cout << "evPlane: Step (" << mCorrStep + 1 << "/3). ";
   if (mCorrStep == 0) cout << "Recentering: collecting. Shift: waiting.";
   if (mCorrStep == 1) cout << "Recentering: applying.   Shift: collecting.";
   if (mCorrStep == 2) cout << "Recentering: applying.   Shift: applying.";
   cout << endl;
   cout << endl;

   mhEvents = new TH1F("mhEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertex);
   mhHits = new TH1F("mhHits", "Number of TPC hits", 100, -0.5, 99.5);
   fOutputList->Add(mhHits);
   mhEta = new TH1F("mhEta", "Eta", 100, -2., 2.);
   fOutputList->Add(mhEta);
   mhDca = new TH1F("mhDca", "DCA", 100, 0., 5.);
   fOutputList->Add(mhDca);
   mhPt = new TH1F("mhPt", "Pt", 300, 0., 3.);
   fOutputList->Add(mhPt);
   mhCorrStep = new TH1F("mhCorrStep", "Correction step: 0 - raw, 1 - rec, 2 - shift", 3, 0., 3.);
   fOutputList->Add(mhCorrStep);
   mhQxRawFHCalFAll = new TH1F("mhQxRawFHCalFAll", "Q_{x}^{Raw} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQxRawFHCalFAll);
   mhQyRawFHCalFAll = new TH1F("mhQyRawFHCalFAll", "Q_{y}^{Raw} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQyRawFHCalFAll);
   mhPhiEPRawFHCalFAll =
      new TH1F("mhPhiEPRawFHCalFAll", "#Psi_{EP}^{Raw} from FHCal F", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRawFHCalFAll);
   mhQxRawFHCalNAll = new TH1F("mhQxRawFHCalNAll", "Q_{x}^{Raw} from FHCal N", 200, -10., 10.);
   fOutputList->Add(mhQxRawFHCalNAll);
   mhQyRawFHCalNAll = new TH1F("mhQyRawFHCalNAll", "Q_{y}^{Raw} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQyRawFHCalNAll);
   mhPhiEPRawFHCalNAll =
      new TH1F("mhPhiEPRawFHCalNAll", "#Psi_{EP}^{Raw} from FHCal N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRawFHCalNAll);
   mhQxRawFHCalSAll = new TH1F("mhQxRawFHCalSAll", "Q_{x}^{Raw} from FHCal S", 200, -10., 10.);
   fOutputList->Add(mhQxRawFHCalSAll);
   mhQyRawFHCalSAll = new TH1F("mhQyRawFHCalSAll", "Q_{y}^{Raw} from FHCal S", 200, -10., 10.);
   fOutputList->Add(mhQyRawFHCalSAll);
   mhPhiEPRawFHCalSAll =
      new TH1F("mhPhiEPRawFHCalSAll", "#Psi_{EP}^{Raw} from FHCal S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRawFHCalSAll);
   mhQxRawTPCNAll = new TH1F("mhQxRawTPCNAll", "Q_{x}^{Raw} from TPC N", 200, -10., 10.);
   fOutputList->Add(mhQxRawTPCNAll);
   mhQyRawTPCNAll = new TH1F("mhQyRawTPCNAll", "Q_{y}^{Raw} from TPC N", 200, -10., 10.);
   fOutputList->Add(mhQyRawTPCNAll);
   mhPhiEPRawTPCNAll = new TH1F("mhPhiEPRawTPCNAll", "#Psi_{EP}^{Raw} from TPC N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRawTPCNAll);
   mhQxRawTPCSAll = new TH1F("mhQxRawTPCSAll", "Q_{x}^{Raw} from TPC S", 200, -10., 10.);
   fOutputList->Add(mhQxRawTPCSAll);
   mhQyRawTPCSAll = new TH1F("mhQyRawTPCSAll", "Q_{y}^{Raw} from TPC S", 200, -10., 10.);
   fOutputList->Add(mhQyRawTPCSAll);
   mhPhiEPRawTPCSAll = new TH1F("mhPhiEPRawTPCSAll", "#Psi_{EP}^{Raw} from TPC S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRawTPCSAll);
   mhQxRecFHCalFAll = new TH1F("mhQxRecFHCalFAll", "Q_{x}^{Rec} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQxRecFHCalFAll);
   mhQyRecFHCalFAll = new TH1F("mhQyRecFHCalFAll", "Q_{y}^{Rec} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQyRecFHCalFAll);
   mhPhiEPRecFHCalFAll =
      new TH1F("mhPhiEPRecFHCalFAll", "#Psi_{EP}^{Rec} from FHCal F", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRecFHCalFAll);
   mhQxRecFHCalNAll = new TH1F("mhQxRecFHCalNAll", "Q_{x}^{Rec} from FHCal N", 200, -10., 10.);
   fOutputList->Add(mhQxRecFHCalNAll);
   mhQyRecFHCalNAll = new TH1F("mhQyRecFHCalNAll", "Q_{y}^{Rec} from FHCal F", 200, -10., 10.);
   fOutputList->Add(mhQyRecFHCalNAll);
   mhPhiEPRecFHCalNAll =
      new TH1F("mhPhiEPRecFHCalNAll", "#Psi_{EP}^{Rec} from FHCal N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRecFHCalNAll);
   mhQxRecFHCalSAll = new TH1F("mhQxRecFHCalSAll", "Q_{x}^{Rec} from FHCal S", 200, -10., 10.);
   fOutputList->Add(mhQxRecFHCalSAll);
   mhQyRecFHCalSAll = new TH1F("mhQyRecFHCalSAll", "Q_{y}^{Rec} from FHCal S", 200, -10., 10.);
   fOutputList->Add(mhQyRecFHCalSAll);
   mhPhiEPRecFHCalSAll =
      new TH1F("mhPhiEPRecFHCalSAll", "#Psi_{EP}^{Rec} from FHCal S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRecFHCalSAll);
   mhQxRecTPCNAll = new TH1F("mhQxRecTPCNAll", "Q_{x}^{Rec} from TPC N", 200, -10., 10.);
   fOutputList->Add(mhQxRecTPCNAll);
   mhQyRecTPCNAll = new TH1F("mhQyRecTPCNAll", "Q_{y}^{Rec} from TPC N", 200, -10., 10.);
   fOutputList->Add(mhQyRecTPCNAll);
   mhPhiEPRecTPCNAll = new TH1F("mhPhiEPRecTPCNAll", "#Psi_{EP}^{Rec} from TPC N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRecTPCNAll);
   mhQxRecTPCSAll = new TH1F("mhQxRecTPCSAll", "Q_{x}^{Rec} from TPC S", 200, -10., 10.);
   fOutputList->Add(mhQxRecTPCSAll);
   mhQyRecTPCSAll = new TH1F("mhQyRecTPCSAll", "Q_{y}^{Rec} from TPC S", 200, -10., 10.);
   fOutputList->Add(mhQyRecTPCSAll);
   mhPhiEPRecTPCSAll = new TH1F("mhPhiEPRecTPCSAll", "#Psi_{EP}^{Rec} from TPC S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPRecTPCSAll);
   mhPhiEPShfFHCalFAll =
      new TH1F("mhPhiEPShfFHCalFAll", "#Psi_{EP}^{Shift} from FHCal F", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPShfFHCalFAll);
   mhPhiEPShfFHCalNAll =
      new TH1F("mhPhiEPShfFHCalNAll", "#Psi_{EP}^{Shift} from FHCal N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPShfFHCalNAll);
   mhPhiEPShfFHCalSAll =
      new TH1F("mhPhiEPShfFHCalSAll", "#Psi_{EP}^{Shift} from FHCal S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPShfFHCalSAll);
   mhPhiEPShfTPCNAll =
      new TH1F("mhPhiEPShfTPCNAll", "#Psi_{EP}^{Shift} from TPC N", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPShfTPCNAll);
   mhPhiEPShfTPCSAll =
      new TH1F("mhPhiEPShfTPCSAll", "#Psi_{EP}^{Shift} from TPC S", 360, -1. * TMath::Pi(), TMath::Pi());
   fOutputList->Add(mhPhiEPShfTPCSAll);
   mhCosFHCalNFHCalSAll = new TProfile("mhCosFHCalNFHCalSAll",
                                       "<cos(#Psi_{EP}^{FHCal N}-#Psi_{EP}^{FHCal S})> vs. centrality", 10, 0., 100.);
   fOutputList->Add(mhCosFHCalNFHCalSAll);
   mhCosTPCNTPCSAll =
      new TProfile("mhCosTPCNTPCSAll", "<cos(2(#Psi_{EP}^{TPC N}-#Psi_{EP}^{TPC S}))> vs. centrality", 10, 0., 100.);
   fOutputList->Add(mhCosTPCNTPCSAll);

   // For the recentering EP correction: <Qx>, <Qy> vs. centrality
   mhCorrQxFHCalFAll = new TProfile("mhCorrQxFHCalFAll", "mhCorrQxFHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQxFHCalFAll);
   mhCorrQyFHCalFAll = new TProfile("mhCorrQyFHCalFAll", "mhCorrQyFHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQyFHCalFAll);
   mhCorrQxFHCalNAll = new TProfile("mhCorrQxFHCalNAll", "mhCorrQxFHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQxFHCalNAll);
   mhCorrQyFHCalNAll = new TProfile("mhCorrQyFHCalNAll", "mhCorrQyFHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQyFHCalNAll);
   mhCorrQxFHCalSAll = new TProfile("mhCorrQxFHCalSAll", "mhCorrQxFHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQxFHCalSAll);
   mhCorrQyFHCalSAll = new TProfile("mhCorrQyFHCalSAll", "mhCorrQyFHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQyFHCalSAll);
   mhCorrQxTPCNAll = new TProfile("mhCorrQxTPCNAll", "mhCorrQxTPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQxTPCNAll);
   mhCorrQyTPCNAll = new TProfile("mhCorrQyTPCNAll", "mhCorrQyTPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQyTPCNAll);
   mhCorrQxTPCSAll = new TProfile("mhCorrQxTPCSAll", "mhCorrQxTPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQxTPCSAll);
   mhCorrQyTPCSAll = new TProfile("mhCorrQyTPCSAll", "mhCorrQyTPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrQyTPCSAll);

   // For the shift EP correction: <cosNphiEP>, <sinNphiEP> vs. centrality
   mhCorrCos1FHCalFAll = new TProfile("mhCorrCos1FHCalFAll", "mhCorrCos1FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos1FHCalFAll);
   mhCorrCos2FHCalFAll = new TProfile("mhCorrCos2FHCalFAll", "mhCorrCos2FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos2FHCalFAll);
   mhCorrCos3FHCalFAll = new TProfile("mhCorrCos3FHCalFAll", "mhCorrCos3FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos3FHCalFAll);
   mhCorrCos4FHCalFAll = new TProfile("mhCorrCos4FHCalFAll", "mhCorrCos4FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos4FHCalFAll);
   mhCorrCos5FHCalFAll = new TProfile("mhCorrCos5FHCalFAll", "mhCorrCos5FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos5FHCalFAll);
   mhCorrCos6FHCalFAll = new TProfile("mhCorrCos6FHCalFAll", "mhCorrCos6FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos6FHCalFAll);
   mhCorrCos7FHCalFAll = new TProfile("mhCorrCos7FHCalFAll", "mhCorrCos7FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos7FHCalFAll);
   mhCorrCos8FHCalFAll = new TProfile("mhCorrCos8FHCalFAll", "mhCorrCos8FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos8FHCalFAll);
   mhCorrSin1FHCalFAll = new TProfile("mhCorrSin1FHCalFAll", "mhCorrSin1FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin1FHCalFAll);
   mhCorrSin2FHCalFAll = new TProfile("mhCorrSin2FHCalFAll", "mhCorrSin2FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin2FHCalFAll);
   mhCorrSin3FHCalFAll = new TProfile("mhCorrSin3FHCalFAll", "mhCorrSin3FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin3FHCalFAll);
   mhCorrSin4FHCalFAll = new TProfile("mhCorrSin4FHCalFAll", "mhCorrSin4FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin4FHCalFAll);
   mhCorrSin5FHCalFAll = new TProfile("mhCorrSin5FHCalFAll", "mhCorrSin5FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin5FHCalFAll);
   mhCorrSin6FHCalFAll = new TProfile("mhCorrSin6FHCalFAll", "mhCorrSin6FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin6FHCalFAll);
   mhCorrSin7FHCalFAll = new TProfile("mhCorrSin7FHCalFAll", "mhCorrSin7FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin7FHCalFAll);
   mhCorrSin8FHCalFAll = new TProfile("mhCorrSin8FHCalFAll", "mhCorrSin8FHCalFAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin8FHCalFAll);

   mhCorrCos1FHCalNAll = new TProfile("mhCorrCos1FHCalNAll", "mhCorrCos1FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos1FHCalNAll);
   mhCorrCos2FHCalNAll = new TProfile("mhCorrCos2FHCalNAll", "mhCorrCos2FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos2FHCalNAll);
   mhCorrCos3FHCalNAll = new TProfile("mhCorrCos3FHCalNAll", "mhCorrCos3FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos3FHCalNAll);
   mhCorrCos4FHCalNAll = new TProfile("mhCorrCos4FHCalNAll", "mhCorrCos4FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos4FHCalNAll);
   mhCorrCos5FHCalNAll = new TProfile("mhCorrCos5FHCalNAll", "mhCorrCos5FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos5FHCalNAll);
   mhCorrCos6FHCalNAll = new TProfile("mhCorrCos6FHCalNAll", "mhCorrCos6FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos6FHCalNAll);
   mhCorrCos7FHCalNAll = new TProfile("mhCorrCos7FHCalNAll", "mhCorrCos7FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos7FHCalNAll);
   mhCorrCos8FHCalNAll = new TProfile("mhCorrCos8FHCalNAll", "mhCorrCos8FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos8FHCalNAll);
   mhCorrSin1FHCalNAll = new TProfile("mhCorrSin1FHCalNAll", "mhCorrSin1FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin1FHCalNAll);
   mhCorrSin2FHCalNAll = new TProfile("mhCorrSin2FHCalNAll", "mhCorrSin2FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin2FHCalNAll);
   mhCorrSin3FHCalNAll = new TProfile("mhCorrSin3FHCalNAll", "mhCorrSin3FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin3FHCalNAll);
   mhCorrSin4FHCalNAll = new TProfile("mhCorrSin4FHCalNAll", "mhCorrSin4FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin4FHCalNAll);
   mhCorrSin5FHCalNAll = new TProfile("mhCorrSin5FHCalNAll", "mhCorrSin5FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin5FHCalNAll);
   mhCorrSin6FHCalNAll = new TProfile("mhCorrSin6FHCalNAll", "mhCorrSin6FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin6FHCalNAll);
   mhCorrSin7FHCalNAll = new TProfile("mhCorrSin7FHCalNAll", "mhCorrSin7FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin7FHCalNAll);
   mhCorrSin8FHCalNAll = new TProfile("mhCorrSin8FHCalNAll", "mhCorrSin8FHCalNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin8FHCalNAll);

   mhCorrCos1FHCalSAll = new TProfile("mhCorrCos1FHCalSAll", "mhCorrCos1FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos1FHCalSAll);
   mhCorrCos2FHCalSAll = new TProfile("mhCorrCos2FHCalSAll", "mhCorrCos2FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos2FHCalSAll);
   mhCorrCos3FHCalSAll = new TProfile("mhCorrCos3FHCalSAll", "mhCorrCos3FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos3FHCalSAll);
   mhCorrCos4FHCalSAll = new TProfile("mhCorrCos4FHCalSAll", "mhCorrCos4FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos4FHCalSAll);
   mhCorrCos5FHCalSAll = new TProfile("mhCorrCos5FHCalSAll", "mhCorrCos5FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos5FHCalSAll);
   mhCorrCos6FHCalSAll = new TProfile("mhCorrCos6FHCalSAll", "mhCorrCos6FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos6FHCalSAll);
   mhCorrCos7FHCalSAll = new TProfile("mhCorrCos7FHCalSAll", "mhCorrCos7FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos7FHCalSAll);
   mhCorrCos8FHCalSAll = new TProfile("mhCorrCos8FHCalSAll", "mhCorrCos8FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos8FHCalSAll);
   mhCorrSin1FHCalSAll = new TProfile("mhCorrSin1FHCalSAll", "mhCorrSin1FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin1FHCalSAll);
   mhCorrSin2FHCalSAll = new TProfile("mhCorrSin2FHCalSAll", "mhCorrSin2FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin2FHCalSAll);
   mhCorrSin3FHCalSAll = new TProfile("mhCorrSin3FHCalSAll", "mhCorrSin3FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin3FHCalSAll);
   mhCorrSin4FHCalSAll = new TProfile("mhCorrSin4FHCalSAll", "mhCorrSin4FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin4FHCalSAll);
   mhCorrSin5FHCalSAll = new TProfile("mhCorrSin5FHCalSAll", "mhCorrSin5FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin5FHCalSAll);
   mhCorrSin6FHCalSAll = new TProfile("mhCorrSin6FHCalSAll", "mhCorrSin6FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin6FHCalSAll);
   mhCorrSin7FHCalSAll = new TProfile("mhCorrSin7FHCalSAll", "mhCorrSin7FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin7FHCalSAll);
   mhCorrSin8FHCalSAll = new TProfile("mhCorrSin8FHCalSAll", "mhCorrSin8FHCalSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin8FHCalSAll);

   mhCorrCos1TPCNAll = new TProfile("mhCorrCos1TPCNAll", "mhCorrCos1TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos1TPCNAll);
   mhCorrCos2TPCNAll = new TProfile("mhCorrCos2TPCNAll", "mhCorrCos2TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos2TPCNAll);
   mhCorrCos3TPCNAll = new TProfile("mhCorrCos3TPCNAll", "mhCorrCos3TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos3TPCNAll);
   mhCorrCos4TPCNAll = new TProfile("mhCorrCos4TPCNAll", "mhCorrCos4TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos4TPCNAll);
   mhCorrSin1TPCNAll = new TProfile("mhCorrSin1TPCNAll", "mhCorrSin1TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin1TPCNAll);
   mhCorrSin2TPCNAll = new TProfile("mhCorrSin2TPCNAll", "mhCorrSin2TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin2TPCNAll);
   mhCorrSin3TPCNAll = new TProfile("mhCorrSin3TPCNAll", "mhCorrSin3TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin3TPCNAll);
   mhCorrSin4TPCNAll = new TProfile("mhCorrSin4TPCNAll", "mhCorrSin4TPCNAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin4TPCNAll);

   mhCorrCos1TPCSAll = new TProfile("mhCorrCos1TPCSAll", "mhCorrCos1TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos1TPCSAll);
   mhCorrCos2TPCSAll = new TProfile("mhCorrCos2TPCSAll", "mhCorrCos2TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos2TPCSAll);
   mhCorrCos3TPCSAll = new TProfile("mhCorrCos3TPCSAll", "mhCorrCos3TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos3TPCSAll);
   mhCorrCos4TPCSAll = new TProfile("mhCorrCos4TPCSAll", "mhCorrCos4TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrCos4TPCSAll);
   mhCorrSin1TPCSAll = new TProfile("mhCorrSin1TPCSAll", "mhCorrSin1TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin1TPCSAll);
   mhCorrSin2TPCSAll = new TProfile("mhCorrSin2TPCSAll", "mhCorrSin2TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin2TPCSAll);
   mhCorrSin3TPCSAll = new TProfile("mhCorrSin3TPCSAll", "mhCorrSin3TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin3TPCSAll);
   mhCorrSin4TPCSAll = new TProfile("mhCorrSin4TPCSAll", "mhCorrSin4TPCSAll", 10, 0., 100.);
   fOutputList->Add(mhCorrSin4TPCSAll);
}

//--------------------------------------
void MpdEventPlaneAll::ProcessEvent(MpdAnalysisEvent &event)
{
   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) {
      return;
   }

   float cent    = event.getCentrTPC();
   int   centBin = GetCentBin(cent);
   if (centBin == -1) return;

   std::vector<float> fhcalEnergy;
   fhcalEnergy.reserve(mFHCalModuleNum);
   for (int i = 0; i < mFHCalModuleNum; i++) {
      fhcalEnergy.push_back(0.);
   }

   // Get energy loss for FHCal modules from fZDCDigit
   for (long int i = 0; i < event.fZDCDigit->GetEntriesFast(); i++) {
      MpdZdcDigi *fhcalhit  = (MpdZdcDigi *)event.fZDCDigit->UncheckedAt(i);
      Int_t       DetId     = fhcalhit->GetDetectorID();
      Int_t       ModId     = fhcalhit->GetModuleID() - 1;
      Int_t       ModNumber = ModId + (mFHCalModuleNum / 2) * (DetId - 1);
      fhcalEnergy.at(ModNumber) += fhcalhit->GetELoss();
   }

   // Get Q-vectors from FHCal
   float Qx_raw_fhcal_F_all = 0., Qy_raw_fhcal_F_all = 0.;
   float Qx_raw_fhcal_N_all = 0., Qy_raw_fhcal_N_all = 0.;
   float Qx_raw_fhcal_S_all = 0., Qy_raw_fhcal_S_all = 0.;
   for (int i = 0; i < mFHCalModuleNum; i++) {
      // Exclude modules 22 and 67 (their position are at beampipe)
      if (i == 22) continue;
      if (i == 67) continue;
      // Full
      float weight = (i < (int)mFHCalMod1Side) ? -1. : 1.;
      Qx_raw_fhcal_F_all += weight * fhcalEnergy.at(i) * cos(1. * GetFHCalPhi(i));
      Qy_raw_fhcal_F_all += weight * fhcalEnergy.at(i) * sin(1. * GetFHCalPhi(i));

      // North (eta<0)
      if (i < (int)mFHCalMod1Side) {
         Qx_raw_fhcal_N_all += fhcalEnergy.at(i) * cos(1. * GetFHCalPhi(i));
         Qy_raw_fhcal_N_all += fhcalEnergy.at(i) * sin(1. * GetFHCalPhi(i));
      }
      // South (eta>0)
      else {
         Qx_raw_fhcal_S_all += fhcalEnergy.at(i) * cos(1. * GetFHCalPhi(i));
         Qy_raw_fhcal_S_all += fhcalEnergy.at(i) * sin(1. * GetFHCalPhi(i));
      }
   }
   float phiEP_raw_fhcal_F_all, phiEP_raw_fhcal_N_all, phiEP_raw_fhcal_S_all;

   phiEP_raw_fhcal_F_all =
      (Qx_raw_fhcal_F_all != 0. && Qy_raw_fhcal_F_all != 0.) ? atan2(Qy_raw_fhcal_F_all, Qx_raw_fhcal_F_all) : -9999.;
   phiEP_raw_fhcal_N_all =
      (Qx_raw_fhcal_N_all != 0. && Qy_raw_fhcal_N_all != 0.) ? atan2(Qy_raw_fhcal_N_all, Qx_raw_fhcal_N_all) : -9999.;
   phiEP_raw_fhcal_S_all =
      (Qx_raw_fhcal_S_all != 0. && Qy_raw_fhcal_S_all != 0.) ? atan2(Qy_raw_fhcal_S_all, Qx_raw_fhcal_S_all) : -9999.;

   // Get Q-vectors from TPC
   mMpdGlobalTracks       = event.fMPDEvent->GetGlobalTracks();
   float Qx_raw_tpc_N_all = 0., Qy_raw_tpc_N_all = 0.;
   float Qx_raw_tpc_S_all = 0., Qy_raw_tpc_S_all = 0.;
   int   mult_N = 0, mult_S = 0;
   for (long int i = 0; i < mMpdGlobalTracks->GetEntriesFast(); i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }

      float weight = mpdtrack->GetPt(); // also can be changed to 1.

      // North (eta<0)
      if (mpdtrack->GetEta() < 0.) {
         Qx_raw_tpc_N_all += weight * cos(2. * mpdtrack->GetPhi());
         Qy_raw_tpc_N_all += weight * sin(2. * mpdtrack->GetPhi());
         mult_N++;
      }
      // South (eta>0)
      else {
         Qx_raw_tpc_S_all += weight * cos(2. * mpdtrack->GetPhi());
         Qy_raw_tpc_S_all += weight * sin(2. * mpdtrack->GetPhi());
         mult_S++;
      }
   }
   float phiEP_raw_tpc_N_all, phiEP_raw_tpc_S_all;

   phiEP_raw_tpc_N_all =
      (Qx_raw_tpc_N_all != 0. && Qy_raw_tpc_N_all != 0.) ? 0.5 * atan2(Qy_raw_tpc_N_all, Qx_raw_tpc_N_all) : -9999.;
   phiEP_raw_tpc_S_all =
      (Qx_raw_tpc_S_all != 0. && Qy_raw_tpc_S_all != 0.) ? 0.5 * atan2(Qy_raw_tpc_S_all, Qx_raw_tpc_S_all) : -9999.;

   if (mult_N < 2) phiEP_raw_tpc_N_all = -9999.;
   if (mult_S < 2) phiEP_raw_tpc_S_all = -9999.;

   // Do recentering correction Q' = Q - <Q>
   float Qx_rec_fhcal_F_all = 0., Qy_rec_fhcal_F_all = 0.;
   float Qx_rec_fhcal_N_all = 0., Qy_rec_fhcal_N_all = 0.;
   float Qx_rec_fhcal_S_all = 0., Qy_rec_fhcal_S_all = 0.;
   float Qx_rec_tpc_N_all = 0., Qy_rec_tpc_N_all = 0.;
   float Qx_rec_tpc_S_all = 0., Qy_rec_tpc_S_all = 0.;
   float phiEP_rec_fhcal_F_all, phiEP_rec_fhcal_N_all, phiEP_rec_fhcal_S_all;
   float phiEP_rec_tpc_N_all, phiEP_rec_tpc_S_all;

   if (mCorrStep == 1 || mCorrStep == 2) {
      Qx_rec_fhcal_F_all = Qx_raw_fhcal_F_all - mCorrQxFHCalFAll.at(centBin);
      Qy_rec_fhcal_F_all = Qy_raw_fhcal_F_all - mCorrQyFHCalFAll.at(centBin);
      Qx_rec_fhcal_N_all = Qx_raw_fhcal_N_all - mCorrQxFHCalNAll.at(centBin);
      Qy_rec_fhcal_N_all = Qy_raw_fhcal_N_all - mCorrQyFHCalNAll.at(centBin);
      Qx_rec_fhcal_S_all = Qx_raw_fhcal_S_all - mCorrQxFHCalSAll.at(centBin);
      Qy_rec_fhcal_S_all = Qy_raw_fhcal_S_all - mCorrQyFHCalSAll.at(centBin);
      Qx_rec_tpc_N_all   = Qx_raw_tpc_N_all - mCorrQxTPCNAll.at(centBin);
      Qy_rec_tpc_N_all   = Qy_raw_tpc_N_all - mCorrQyTPCNAll.at(centBin);
      Qx_rec_tpc_S_all   = Qx_raw_tpc_S_all - mCorrQxTPCSAll.at(centBin);
      Qy_rec_tpc_S_all   = Qy_raw_tpc_S_all - mCorrQyTPCSAll.at(centBin);

      phiEP_rec_fhcal_F_all = (Qx_rec_fhcal_F_all != 0. && Qy_rec_fhcal_F_all != 0.)
                                 ? atan2(Qy_rec_fhcal_F_all, Qx_rec_fhcal_F_all)
                                 : -9999.;
      phiEP_rec_fhcal_N_all = (Qx_rec_fhcal_N_all != 0. && Qy_rec_fhcal_N_all != 0.)
                                 ? atan2(Qy_rec_fhcal_N_all, Qx_rec_fhcal_N_all)
                                 : -9999.;
      phiEP_rec_fhcal_S_all = (Qx_rec_fhcal_S_all != 0. && Qy_rec_fhcal_S_all != 0.)
                                 ? atan2(Qy_rec_fhcal_S_all, Qx_rec_fhcal_S_all)
                                 : -9999.;
      phiEP_rec_tpc_N_all =
         (Qx_rec_tpc_N_all != 0. && Qy_rec_tpc_N_all != 0.) ? 0.5 * atan2(Qy_rec_tpc_N_all, Qx_rec_tpc_N_all) : -9999.;
      phiEP_rec_tpc_S_all =
         (Qx_rec_tpc_S_all != 0. && Qy_rec_tpc_S_all != 0.) ? 0.5 * atan2(Qy_rec_tpc_S_all, Qx_rec_tpc_S_all) : -9999.;

      if (phiEP_raw_fhcal_F_all == -9999.) phiEP_rec_fhcal_F_all = -9999.;
      if (phiEP_raw_fhcal_N_all == -9999.) phiEP_rec_fhcal_N_all = -9999.;
      if (phiEP_raw_fhcal_S_all == -9999.) phiEP_rec_fhcal_S_all = -9999.;
      if (phiEP_raw_tpc_N_all == -9999.) phiEP_rec_tpc_N_all = -9999.;
      if (phiEP_raw_tpc_S_all == -9999.) phiEP_rec_tpc_S_all = -9999.;
   }

   // Do shift correction Psi' = Psi + dPsi - for dPsi see eq.(6) in arXiv:nucl-ex/9805001
   float phiEP_shf_fhcal_F_all, phiEP_shf_fhcal_N_all, phiEP_shf_fhcal_S_all;
   float phiEP_shf_tpc_N_all, phiEP_shf_tpc_S_all;

   if (mCorrStep == 2) {
      if (phiEP_raw_fhcal_F_all != -9999. && phiEP_rec_fhcal_F_all != -9999.) {
         phiEP_shf_fhcal_F_all = phiEP_rec_fhcal_F_all +
                                 (2. / 1.) * (mCorrCos1FHCalFAll.at(centBin) * sin(1. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin1FHCalFAll.at(centBin) * cos(1. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 2.) * (mCorrCos2FHCalFAll.at(centBin) * sin(2. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin2FHCalFAll.at(centBin) * cos(2. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 3.) * (mCorrCos3FHCalFAll.at(centBin) * sin(3. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin3FHCalFAll.at(centBin) * cos(3. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 4.) * (mCorrCos4FHCalFAll.at(centBin) * sin(4. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin4FHCalFAll.at(centBin) * cos(4. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 5.) * (mCorrCos5FHCalFAll.at(centBin) * sin(5. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin5FHCalFAll.at(centBin) * cos(5. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 6.) * (mCorrCos6FHCalFAll.at(centBin) * sin(6. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin6FHCalFAll.at(centBin) * cos(6. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 7.) * (mCorrCos7FHCalFAll.at(centBin) * sin(7. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin7FHCalFAll.at(centBin) * cos(7. * phiEP_rec_fhcal_F_all)) +
                                 (2. / 8.) * (mCorrCos8FHCalFAll.at(centBin) * sin(8. * phiEP_rec_fhcal_F_all) -
                                              mCorrSin8FHCalFAll.at(centBin) * cos(8. * phiEP_rec_fhcal_F_all));
      } else {
         phiEP_shf_fhcal_F_all = -9999.;
      }
      if (phiEP_raw_fhcal_N_all != -9999. && phiEP_rec_fhcal_N_all != -9999.) {
         phiEP_shf_fhcal_N_all = phiEP_rec_fhcal_N_all +
                                 (2. / 1.) * (mCorrCos1FHCalNAll.at(centBin) * sin(1. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin1FHCalNAll.at(centBin) * cos(1. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 2.) * (mCorrCos2FHCalNAll.at(centBin) * sin(2. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin2FHCalNAll.at(centBin) * cos(2. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 3.) * (mCorrCos3FHCalNAll.at(centBin) * sin(3. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin3FHCalNAll.at(centBin) * cos(3. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 4.) * (mCorrCos4FHCalNAll.at(centBin) * sin(4. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin4FHCalNAll.at(centBin) * cos(4. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 5.) * (mCorrCos5FHCalNAll.at(centBin) * sin(5. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin5FHCalNAll.at(centBin) * cos(5. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 6.) * (mCorrCos6FHCalNAll.at(centBin) * sin(6. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin6FHCalNAll.at(centBin) * cos(6. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 7.) * (mCorrCos7FHCalNAll.at(centBin) * sin(7. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin7FHCalNAll.at(centBin) * cos(7. * phiEP_rec_fhcal_N_all)) +
                                 (2. / 8.) * (mCorrCos8FHCalNAll.at(centBin) * sin(8. * phiEP_rec_fhcal_N_all) -
                                              mCorrSin8FHCalNAll.at(centBin) * cos(8. * phiEP_rec_fhcal_N_all));
      } else {
         phiEP_shf_fhcal_N_all = -9999.;
      }
      if (phiEP_raw_fhcal_S_all != -9999. && phiEP_rec_fhcal_S_all != -9999.) {
         phiEP_shf_fhcal_S_all = phiEP_rec_fhcal_S_all +
                                 (2. / 1.) * (mCorrCos1FHCalSAll.at(centBin) * sin(1. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin1FHCalSAll.at(centBin) * cos(1. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 2.) * (mCorrCos2FHCalSAll.at(centBin) * sin(2. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin2FHCalSAll.at(centBin) * cos(2. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 3.) * (mCorrCos3FHCalSAll.at(centBin) * sin(3. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin3FHCalSAll.at(centBin) * cos(3. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 4.) * (mCorrCos4FHCalSAll.at(centBin) * sin(4. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin4FHCalSAll.at(centBin) * cos(4. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 5.) * (mCorrCos5FHCalSAll.at(centBin) * sin(5. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin5FHCalSAll.at(centBin) * cos(5. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 6.) * (mCorrCos6FHCalSAll.at(centBin) * sin(6. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin6FHCalSAll.at(centBin) * cos(6. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 7.) * (mCorrCos7FHCalSAll.at(centBin) * sin(7. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin7FHCalSAll.at(centBin) * cos(7. * phiEP_rec_fhcal_S_all)) +
                                 (2. / 8.) * (mCorrCos8FHCalSAll.at(centBin) * sin(8. * phiEP_rec_fhcal_S_all) -
                                              mCorrSin8FHCalSAll.at(centBin) * cos(8. * phiEP_rec_fhcal_S_all));
      } else {
         phiEP_shf_fhcal_S_all = -9999.;
      }

      if (phiEP_raw_tpc_N_all != -9999. && phiEP_rec_tpc_N_all != -9999.) {
         phiEP_shf_tpc_N_all = phiEP_rec_tpc_N_all +
                               (1. / 2.) * (2. / 1.) *
                                  (mCorrCos1TPCNAll.at(centBin) * sin(1. * 2. * phiEP_rec_tpc_N_all) -
                                   mCorrSin1TPCNAll.at(centBin) * cos(1. * 2. * phiEP_rec_tpc_N_all)) +
                               (1. / 2.) * (2. / 2.) *
                                  (mCorrCos2TPCNAll.at(centBin) * sin(2. * 2. * phiEP_rec_tpc_N_all) -
                                   mCorrSin2TPCNAll.at(centBin) * cos(2. * 2. * phiEP_rec_tpc_N_all)) +
                               (1. / 2.) * (2. / 3.) *
                                  (mCorrCos3TPCNAll.at(centBin) * sin(3. * 2. * phiEP_rec_tpc_N_all) -
                                   mCorrSin3TPCNAll.at(centBin) * cos(3. * 2. * phiEP_rec_tpc_N_all)) +
                               (1. / 2.) * (2. / 4.) *
                                  (mCorrCos4TPCNAll.at(centBin) * sin(4. * 2. * phiEP_rec_tpc_N_all) -
                                   mCorrSin4TPCNAll.at(centBin) * cos(4. * 2. * phiEP_rec_tpc_N_all));
      } else {
         phiEP_shf_tpc_N_all = -9999.;
      }
      if (phiEP_raw_tpc_S_all != -9999. && phiEP_rec_tpc_S_all != -9999.) {
         phiEP_shf_tpc_S_all = phiEP_rec_tpc_S_all +
                               (1. / 2.) * (2. / 1.) *
                                  (mCorrCos1TPCSAll.at(centBin) * sin(1. * 2. * phiEP_rec_tpc_S_all) -
                                   mCorrSin1TPCSAll.at(centBin) * cos(1. * 2. * phiEP_rec_tpc_S_all)) +
                               (1. / 2.) * (2. / 2.) *
                                  (mCorrCos2TPCSAll.at(centBin) * sin(2. * 2. * phiEP_rec_tpc_S_all) -
                                   mCorrSin2TPCSAll.at(centBin) * cos(2. * 2. * phiEP_rec_tpc_S_all)) +
                               (1. / 2.) * (2. / 3.) *
                                  (mCorrCos3TPCSAll.at(centBin) * sin(3. * 2. * phiEP_rec_tpc_S_all) -
                                   mCorrSin3TPCSAll.at(centBin) * cos(3. * 2. * phiEP_rec_tpc_S_all)) +
                               (1. / 2.) * (2. / 4.) *
                                  (mCorrCos4TPCSAll.at(centBin) * sin(4. * 2. * phiEP_rec_tpc_S_all) -
                                   mCorrSin4TPCSAll.at(centBin) * cos(4. * 2. * phiEP_rec_tpc_S_all));
      } else {
         phiEP_shf_tpc_S_all = -9999.;
      }
   }

   // Fill general QA histograms
   mhCorrStep->Fill(mCorrStep);

   mhQxRawFHCalFAll->Fill(Qx_raw_fhcal_F_all);
   mhQyRawFHCalFAll->Fill(Qy_raw_fhcal_F_all);
   mhPhiEPRawFHCalFAll->Fill(phiEP_raw_fhcal_F_all);
   mhQxRawFHCalNAll->Fill(Qx_raw_fhcal_N_all);
   mhQyRawFHCalNAll->Fill(Qy_raw_fhcal_N_all);
   mhPhiEPRawFHCalNAll->Fill(phiEP_raw_fhcal_N_all);
   mhQxRawFHCalSAll->Fill(Qx_raw_fhcal_S_all);
   mhQyRawFHCalSAll->Fill(Qy_raw_fhcal_S_all);
   mhPhiEPRawFHCalSAll->Fill(phiEP_raw_fhcal_S_all);
   mhQxRawTPCNAll->Fill(Qx_raw_tpc_N_all);
   mhQyRawTPCNAll->Fill(Qy_raw_tpc_N_all);
   mhPhiEPRawTPCNAll->Fill(phiEP_raw_tpc_N_all);
   mhQxRawTPCSAll->Fill(Qx_raw_tpc_S_all);
   mhQyRawTPCSAll->Fill(Qy_raw_tpc_S_all);
   mhPhiEPRawTPCSAll->Fill(phiEP_raw_tpc_S_all);
   if (mCorrStep == 1 || mCorrStep == 2) {
      mhQxRecFHCalFAll->Fill(Qx_rec_fhcal_F_all);
      mhQyRecFHCalFAll->Fill(Qy_rec_fhcal_F_all);
      mhPhiEPRecFHCalFAll->Fill(phiEP_rec_fhcal_F_all);
      mhQxRecFHCalNAll->Fill(Qx_rec_fhcal_N_all);
      mhQyRecFHCalNAll->Fill(Qy_rec_fhcal_N_all);
      mhPhiEPRecFHCalNAll->Fill(phiEP_rec_fhcal_N_all);
      mhQxRecFHCalSAll->Fill(Qx_rec_fhcal_S_all);
      mhQyRecFHCalSAll->Fill(Qy_rec_fhcal_S_all);
      mhPhiEPRecFHCalSAll->Fill(phiEP_rec_fhcal_S_all);
      mhQxRecTPCNAll->Fill(Qx_rec_tpc_N_all);
      mhQyRecTPCNAll->Fill(Qy_rec_tpc_N_all);
      mhPhiEPRecTPCNAll->Fill(phiEP_rec_tpc_N_all);
      mhQxRecTPCSAll->Fill(Qx_rec_tpc_S_all);
      mhQyRecTPCSAll->Fill(Qy_rec_tpc_S_all);
      mhPhiEPRecTPCSAll->Fill(phiEP_rec_tpc_S_all);
   }
   if (mCorrStep == 2) {
      mhPhiEPShfFHCalFAll->Fill(phiEP_shf_fhcal_F_all);
      mhPhiEPShfFHCalNAll->Fill(phiEP_shf_fhcal_N_all);
      mhPhiEPShfFHCalSAll->Fill(phiEP_shf_fhcal_S_all);
      mhPhiEPShfTPCNAll->Fill(phiEP_shf_tpc_N_all);
      mhPhiEPShfTPCSAll->Fill(phiEP_shf_tpc_S_all);
   }

   // Fill correction info profiles
   //
   // Info for recentering <Qx>, <Qy>
   mhCorrQxFHCalFAll->Fill(cent, Qx_raw_fhcal_F_all);
   mhCorrQyFHCalFAll->Fill(cent, Qy_raw_fhcal_F_all);
   mhCorrQxFHCalNAll->Fill(cent, Qx_raw_fhcal_N_all);
   mhCorrQyFHCalNAll->Fill(cent, Qy_raw_fhcal_N_all);
   mhCorrQxFHCalSAll->Fill(cent, Qx_raw_fhcal_S_all);
   mhCorrQyFHCalSAll->Fill(cent, Qy_raw_fhcal_S_all);
   mhCorrQxTPCNAll->Fill(cent, Qx_raw_tpc_N_all);
   mhCorrQyTPCNAll->Fill(cent, Qy_raw_tpc_N_all);
   mhCorrQxTPCSAll->Fill(cent, Qx_raw_tpc_S_all);
   mhCorrQyTPCSAll->Fill(cent, Qy_raw_tpc_S_all);

   // Info for shift <cosNpsi>, <sinNpsi>
   if (mCorrStep == 1 || mCorrStep == 2) {
      if (phiEP_rec_fhcal_F_all != -9999.) {
         mhCorrCos1FHCalFAll->Fill(cent, cos(1. * phiEP_rec_fhcal_F_all));
         mhCorrCos2FHCalFAll->Fill(cent, cos(2. * phiEP_rec_fhcal_F_all));
         mhCorrCos3FHCalFAll->Fill(cent, cos(3. * phiEP_rec_fhcal_F_all));
         mhCorrCos4FHCalFAll->Fill(cent, cos(4. * phiEP_rec_fhcal_F_all));
         mhCorrCos5FHCalFAll->Fill(cent, cos(5. * phiEP_rec_fhcal_F_all));
         mhCorrCos6FHCalFAll->Fill(cent, cos(6. * phiEP_rec_fhcal_F_all));
         mhCorrCos7FHCalFAll->Fill(cent, cos(7. * phiEP_rec_fhcal_F_all));
         mhCorrCos8FHCalFAll->Fill(cent, cos(8. * phiEP_rec_fhcal_F_all));
         mhCorrSin1FHCalFAll->Fill(cent, sin(1. * phiEP_rec_fhcal_F_all));
         mhCorrSin2FHCalFAll->Fill(cent, sin(2. * phiEP_rec_fhcal_F_all));
         mhCorrSin3FHCalFAll->Fill(cent, sin(3. * phiEP_rec_fhcal_F_all));
         mhCorrSin4FHCalFAll->Fill(cent, sin(4. * phiEP_rec_fhcal_F_all));
         mhCorrSin5FHCalFAll->Fill(cent, sin(5. * phiEP_rec_fhcal_F_all));
         mhCorrSin6FHCalFAll->Fill(cent, sin(6. * phiEP_rec_fhcal_F_all));
         mhCorrSin7FHCalFAll->Fill(cent, sin(7. * phiEP_rec_fhcal_F_all));
         mhCorrSin8FHCalFAll->Fill(cent, sin(8. * phiEP_rec_fhcal_F_all));
      }

      if (phiEP_rec_fhcal_N_all != -9999.) {
         mhCorrCos1FHCalNAll->Fill(cent, cos(1. * phiEP_rec_fhcal_N_all));
         mhCorrCos2FHCalNAll->Fill(cent, cos(2. * phiEP_rec_fhcal_N_all));
         mhCorrCos3FHCalNAll->Fill(cent, cos(3. * phiEP_rec_fhcal_N_all));
         mhCorrCos4FHCalNAll->Fill(cent, cos(4. * phiEP_rec_fhcal_N_all));
         mhCorrCos5FHCalNAll->Fill(cent, cos(5. * phiEP_rec_fhcal_N_all));
         mhCorrCos6FHCalNAll->Fill(cent, cos(6. * phiEP_rec_fhcal_N_all));
         mhCorrCos7FHCalNAll->Fill(cent, cos(7. * phiEP_rec_fhcal_N_all));
         mhCorrCos8FHCalNAll->Fill(cent, cos(8. * phiEP_rec_fhcal_N_all));
         mhCorrSin1FHCalNAll->Fill(cent, sin(1. * phiEP_rec_fhcal_N_all));
         mhCorrSin2FHCalNAll->Fill(cent, sin(2. * phiEP_rec_fhcal_N_all));
         mhCorrSin3FHCalNAll->Fill(cent, sin(3. * phiEP_rec_fhcal_N_all));
         mhCorrSin4FHCalNAll->Fill(cent, sin(4. * phiEP_rec_fhcal_N_all));
         mhCorrSin5FHCalNAll->Fill(cent, sin(5. * phiEP_rec_fhcal_N_all));
         mhCorrSin6FHCalNAll->Fill(cent, sin(6. * phiEP_rec_fhcal_N_all));
         mhCorrSin7FHCalNAll->Fill(cent, sin(7. * phiEP_rec_fhcal_N_all));
         mhCorrSin8FHCalNAll->Fill(cent, sin(8. * phiEP_rec_fhcal_N_all));
      }

      if (phiEP_rec_fhcal_S_all != -9999.) {
         mhCorrCos1FHCalSAll->Fill(cent, cos(1. * phiEP_rec_fhcal_S_all));
         mhCorrCos2FHCalSAll->Fill(cent, cos(2. * phiEP_rec_fhcal_S_all));
         mhCorrCos3FHCalSAll->Fill(cent, cos(3. * phiEP_rec_fhcal_S_all));
         mhCorrCos4FHCalSAll->Fill(cent, cos(4. * phiEP_rec_fhcal_S_all));
         mhCorrCos5FHCalSAll->Fill(cent, cos(5. * phiEP_rec_fhcal_S_all));
         mhCorrCos6FHCalSAll->Fill(cent, cos(6. * phiEP_rec_fhcal_S_all));
         mhCorrCos7FHCalSAll->Fill(cent, cos(7. * phiEP_rec_fhcal_S_all));
         mhCorrCos8FHCalSAll->Fill(cent, cos(8. * phiEP_rec_fhcal_S_all));
         mhCorrSin1FHCalSAll->Fill(cent, sin(1. * phiEP_rec_fhcal_S_all));
         mhCorrSin2FHCalSAll->Fill(cent, sin(2. * phiEP_rec_fhcal_S_all));
         mhCorrSin3FHCalSAll->Fill(cent, sin(3. * phiEP_rec_fhcal_S_all));
         mhCorrSin4FHCalSAll->Fill(cent, sin(4. * phiEP_rec_fhcal_S_all));
         mhCorrSin5FHCalSAll->Fill(cent, sin(5. * phiEP_rec_fhcal_S_all));
         mhCorrSin6FHCalSAll->Fill(cent, sin(6. * phiEP_rec_fhcal_S_all));
         mhCorrSin7FHCalSAll->Fill(cent, sin(7. * phiEP_rec_fhcal_S_all));
         mhCorrSin8FHCalSAll->Fill(cent, sin(8. * phiEP_rec_fhcal_S_all));
      }

      if (phiEP_rec_tpc_N_all != -9999.) {
         mhCorrCos1TPCNAll->Fill(cent, cos(1. * 2. * phiEP_rec_tpc_N_all));
         mhCorrCos2TPCNAll->Fill(cent, cos(2. * 2. * phiEP_rec_tpc_N_all));
         mhCorrCos3TPCNAll->Fill(cent, cos(3. * 2. * phiEP_rec_tpc_N_all));
         mhCorrCos4TPCNAll->Fill(cent, cos(4. * 2. * phiEP_rec_tpc_N_all));
         mhCorrSin1TPCNAll->Fill(cent, sin(1. * 2. * phiEP_rec_tpc_N_all));
         mhCorrSin2TPCNAll->Fill(cent, sin(2. * 2. * phiEP_rec_tpc_N_all));
         mhCorrSin3TPCNAll->Fill(cent, sin(3. * 2. * phiEP_rec_tpc_N_all));
         mhCorrSin4TPCNAll->Fill(cent, sin(4. * 2. * phiEP_rec_tpc_N_all));
      }

      if (phiEP_rec_tpc_S_all != -9999.) {
         mhCorrCos1TPCSAll->Fill(cent, cos(1. * 2. * phiEP_rec_tpc_S_all));
         mhCorrCos2TPCSAll->Fill(cent, cos(2. * 2. * phiEP_rec_tpc_S_all));
         mhCorrCos3TPCSAll->Fill(cent, cos(3. * 2. * phiEP_rec_tpc_S_all));
         mhCorrCos4TPCSAll->Fill(cent, cos(4. * 2. * phiEP_rec_tpc_S_all));
         mhCorrSin1TPCSAll->Fill(cent, sin(1. * 2. * phiEP_rec_tpc_S_all));
         mhCorrSin2TPCSAll->Fill(cent, sin(2. * 2. * phiEP_rec_tpc_S_all));
         mhCorrSin3TPCSAll->Fill(cent, sin(3. * 2. * phiEP_rec_tpc_S_all));
         mhCorrSin4TPCSAll->Fill(cent, sin(4. * 2. * phiEP_rec_tpc_S_all));
      }
   }

   // Set result in MpdAnalysisEvent (MpdAnalysisEventPlane)
   float phiEP_result_fhcal_F_all, phiEP_result_fhcal_N_all, phiEP_result_fhcal_S_all;
   float phiEP_result_tpc_N_all, phiEP_result_tpc_S_all;

   if (mCorrStep == 0) {
      phiEP_result_fhcal_F_all = phiEP_raw_fhcal_F_all;
      phiEP_result_fhcal_N_all = phiEP_raw_fhcal_N_all;
      phiEP_result_fhcal_S_all = phiEP_raw_fhcal_S_all;
      phiEP_result_tpc_N_all   = phiEP_raw_tpc_N_all;
      phiEP_result_tpc_S_all   = phiEP_raw_tpc_S_all;
   }
   if (mCorrStep == 1) {
      phiEP_result_fhcal_F_all = phiEP_rec_fhcal_F_all;
      phiEP_result_fhcal_N_all = phiEP_rec_fhcal_N_all;
      phiEP_result_fhcal_S_all = phiEP_rec_fhcal_S_all;
      phiEP_result_tpc_N_all   = phiEP_rec_tpc_N_all;
      phiEP_result_tpc_S_all   = phiEP_rec_tpc_S_all;
   }
   if (mCorrStep == 2) {
      phiEP_result_fhcal_F_all = phiEP_shf_fhcal_F_all;
      phiEP_result_fhcal_N_all = phiEP_shf_fhcal_N_all;
      phiEP_result_fhcal_S_all = phiEP_shf_fhcal_S_all;
      phiEP_result_tpc_N_all   = phiEP_shf_tpc_N_all;
      phiEP_result_tpc_S_all   = phiEP_shf_tpc_S_all;
   }

   mhCosFHCalNFHCalSAll->Fill(cent, cos(1. * (phiEP_result_fhcal_N_all - phiEP_result_fhcal_S_all)));
   mhCosTPCNTPCSAll->Fill(cent, cos(2. * (phiEP_result_tpc_N_all - phiEP_result_tpc_S_all)));

   event.fMpdEP.SetPhiEP_FHCal_F_all(phiEP_result_fhcal_F_all);
   event.fMpdEP.SetPhiEP_FHCal_N_all(phiEP_result_fhcal_N_all);
   event.fMpdEP.SetPhiEP_FHCal_S_all(phiEP_result_fhcal_S_all);
   event.fMpdEP.SetPhiEP_TPC_N_all(phiEP_result_tpc_N_all);
   event.fMpdEP.SetPhiEP_TPC_S_all(phiEP_result_tpc_S_all);
}

void MpdEventPlaneAll::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdEventPlaneAll::selectEvent(MpdAnalysisEvent &event)
{
   mhEvents->Fill(0.5);

   // Reject empty events (UrQMD, PHSD)
   mMCTracks = event.fMCTrack;

   int nTrMc = 0;
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      MpdMCTrack *pr  = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
      float       eta = 0.5 * log((pr->GetP() + pr->GetPz()) / (pr->GetP() - pr->GetPz()));
      if (pr->GetMotherId() == -1 && fabs(eta) <= mParams.mEtaCut && fabs(eta) >= mParams.mEtaGapCut / 2.) {
         nTrMc++;
      }
   } // i

   if (nTrMc <= 2) { // check if there're at least 2 mc particles within eta cuts
      return false;
   }

   mhEvents->Fill(1.5);

   // Reject bad vertex
   if (!event.fVertex) {
      return false;
   }

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());

   if (mPrimaryVertex.Z() == 0) { // not reconstructed (==0)
      return false;
   }

   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) { // beyond the limits
      return false;
   }

   mhEvents->Fill(2.5);

   return true;
}

bool MpdEventPlaneAll::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 16

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false;         // |eta| < 1.5
   if (fabs(mpdtrack->GetEta()) < mParams.mEtaGapCut / 2.) return false; // |eta| > 0.05

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 100 MeV/c
   if (pt > mParams.mPtmaxCut) return false; // pT < 2000 MeV/c

   if (fabs(mpdtrack->GetDCAX()) > mParams.mDcaCut) return false; // |DCAx| < 2.
   if (fabs(mpdtrack->GetDCAY()) > mParams.mDcaCut) return false; // |DCAy| < 2.
   if (fabs(mpdtrack->GetDCAZ()) > mParams.mDcaCut) return false; // |DCAz| < 2.

   mhHits->Fill(mpdtrack->GetNofHits());
   mhPt->Fill(pt);
   mhEta->Fill(mpdtrack->GetEta());
   mhDca->Fill(sqrt(mpdtrack->GetDCAX() * mpdtrack->GetDCAX() + mpdtrack->GetDCAY() * mpdtrack->GetDCAY() +
                    mpdtrack->GetDCAZ() * mpdtrack->GetDCAZ()));

   return true;
}

float MpdEventPlaneAll::GetFHCalPhi(int iModule)
{
   const int Nmodules    = 45;
   int       xAxisSwitch = (iModule < Nmodules) ? 1 : -1;
   int       module      = (iModule < Nmodules) ? iModule : iModule - Nmodules;
   float     x, y, phi = -999.;
   if (module >= 0 && module <= 4) {
      y   = 45.;
      x   = -1. * (module - 2) * 15.;
      phi = TMath::ATan2(y, x * xAxisSwitch);
   } else if ((module >= 5) && (module <= 39)) {
      y   = (3 - (module + 2) / 7) * 15.;
      x   = (3 - (module + 2) % 7) * 15.;
      phi = TMath::ATan2(y, x * xAxisSwitch);
   } else if ((module >= 40) && (module <= 44)) {
      y   = -45.;
      x   = -1. * (module - 42) * 15.;
      phi = TMath::ATan2(y, x * xAxisSwitch);
   }
   return phi;
}

int MpdEventPlaneAll::GetCentBin(float cent)
{
   if (cent > 0 && cent <= 10)
      return 0;
   else if (cent > 10 && cent <= 20)
      return 1;
   else if (cent > 20 && cent <= 30)
      return 2;
   else if (cent > 30 && cent <= 40)
      return 3;
   else if (cent > 40 && cent <= 50)
      return 4;
   else if (cent > 50 && cent <= 60)
      return 5;
   else if (cent > 60 && cent <= 70)
      return 6;
   else if (cent > 70 && cent <= 80)
      return 7;
   else
      return -1;
}

std::map<int, float> MpdEventPlaneAll::SetZeroCorr()
{
   std::map<int, float> result;
   for (int i = 0; i < nCentBins; i++) {
      result.insert({i, 0.});
   }
   return result;
}

std::map<int, float> MpdEventPlaneAll::ReadEpCorrProfile(TProfile *const &prof)
{
   std::map<int, float> result = SetZeroCorr();

   if (!prof) {
      std::cerr << "evPlane: Cannot read EP correction TProfile from the file!" << std::endl;
      return result;
   }

   for (int i = 0; i < prof->GetNbinsX(); i++) {
      int centbin = GetCentBin(prof->GetBinCenter(i + 1));
      if (centbin < 0) continue;
      result.at(centbin) = prof->GetBinContent(i + 1);
   }

   // cout << "MpdEventPlaneAll::ReadEpCorrProfile: Contents of the EP correction profile " << prof->GetName() << ":" <<
   // endl; for (auto &element:result) {
   //    cout << "\tCentrality bin " << element.first << ": " << element.second << endl;
   // }

   return result;
}