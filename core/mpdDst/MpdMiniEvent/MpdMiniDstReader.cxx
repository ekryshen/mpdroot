//
// MpdMiniDstReader allows to read miniDst file or a list of files
//

// C++ headers
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <assert.h>

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniDstReader.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniBTofHit.h"
#include "MpdMiniBTofPidTraits.h"
#include "MpdMiniBECalCluster.h"
#include "MpdMiniTrackCovMatrix.h"
#include "MpdMiniFHCalHit.h"
#include "MpdMiniMcEvent.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniArrays.h"
#include "MpdMiniDst.h"

// ROOT headers
#include "TRegexp.h"

ClassImp(MpdMiniDstReader);

//_________________
MpdMiniDstReader::MpdMiniDstReader(const Char_t *inFileName)
   : fMiniDst(new MpdMiniDst()), fChain(nullptr), fTree(nullptr), fEventCounter(0), fMiniArrays{}, fStatusArrays{}
{

   streamerOff();
   createArrays();
   std::fill_n(fStatusArrays, sizeof(fStatusArrays) / sizeof(fStatusArrays[0]), 1);
   mInputFileName = inFileName;
}

//_________________
MpdMiniDstReader::~MpdMiniDstReader()
{
   if (fChain) {
      delete fChain;
   }
   if (fMiniDst) {
      delete fMiniDst;
   }
}

//_________________
void MpdMiniDstReader::clearArrays()
{
   for (Int_t iArr = 0; iArr < MpdMiniArrays::NAllMiniArrays; iArr++) {
      fMiniArrays[iArr]->Clear();
   }
}

//_________________
void MpdMiniDstReader::SetStatus(const Char_t *branchNameRegex, Int_t enable)
{
   if (strncmp(branchNameRegex, "Mpd", 3) == 0) {
      // Ignore first "Mpd"
      branchNameRegex += 3;
   }

   TRegexp re(branchNameRegex, 1);
   for (Int_t iArr = 0; iArr < MpdMiniArrays::NAllMiniArrays; iArr++) {
      Ssiz_t len;
      if (re.Index(MpdMiniArrays::miniArrayNames[iArr], &len) < 0) continue;
      LOG_INFO << "MpdMiniDstMaker::SetStatus " << enable << " to " << MpdMiniArrays::miniArrayNames[iArr] << endm;
      fStatusArrays[iArr] = enable;
   }

   setBranchAddresses(fChain);
}

//_________________
void MpdMiniDstReader::setBranchAddresses(TChain *chain)
{
   if (!chain) return;
   chain->SetBranchStatus("*", 0);
   TString ts;
   for (Int_t i = 0; i < MpdMiniArrays::NAllMiniArrays; ++i) {
      if (fStatusArrays[i] == 0) continue;
      char const *bname = MpdMiniArrays::miniArrayNames[i];
      TBranch    *tb    = chain->GetBranch(bname);
      if (!tb) {
         LOG_WARN << "setBranchAddress: Branch name " << bname << " does not exist!" << endm;
         continue;
      }
      ts = bname;
      ts += "*";
      chain->SetBranchStatus(ts, 1);
      chain->SetBranchAddress(bname, fMiniArrays + i);
      assert(tb->GetAddress() == (char *)(fMiniArrays + i));
   }
   fTree = fChain->GetTree();
}

//_________________
void MpdMiniDstReader::streamerOff()
{
   MpdMiniEvent::Class()->IgnoreTObjectStreamer();
   MpdMiniTrack::Class()->IgnoreTObjectStreamer();
   MpdMiniBTofHit::Class()->IgnoreTObjectStreamer();
   MpdMiniBTofPidTraits::Class()->IgnoreTObjectStreamer();
   MpdMiniBECalCluster::Class()->IgnoreTObjectStreamer();
   MpdMiniTrackCovMatrix::Class()->IgnoreTObjectStreamer();
   MpdMiniFHCalHit::Class()->IgnoreTObjectStreamer();
   MpdMiniMcEvent::Class()->IgnoreTObjectStreamer();
   MpdMiniMcTrack::Class()->IgnoreTObjectStreamer();
}

//_________________
void MpdMiniDstReader::createArrays()
{
   for (Int_t iArr = 0; iArr < MpdMiniArrays::NAllMiniArrays; iArr++) {
      fMiniArrays[iArr] = new TClonesArray(MpdMiniArrays::miniArrayTypes[iArr], MpdMiniArrays::miniArraySizes[iArr]);
   }
   fMiniDst->set(fMiniArrays);
}

//_________________
void MpdMiniDstReader::Finish()
{
   if (fChain) {
      delete fChain;
   }
   fChain = NULL;
}

//_________________
void MpdMiniDstReader::Init()
{

   if (!fChain) {
      fChain = new TChain("MiniDst");
   }

   std::string const dirFile = mInputFileName.Data();

   if (dirFile.find(".list") != std::string::npos || dirFile.find(".lis") != std::string::npos) {

      std::ifstream inputStream(dirFile.c_str());

      if (!inputStream) {
         LOG_ERROR << "ERROR: Cannot open list file " << dirFile << endm;
      }

      Int_t       nFile = 0;
      std::string file;
      while (getline(inputStream, file)) {
         if (file.find(".miniDst.root") != std::string::npos || file.find(".MiniDst.root") != std::string::npos) {
            TFile *ftmp = TFile::Open(file.c_str());
            if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) {
               LOG_INFO << " Read in miniDst file " << file << endm;
               fChain->Add(file.c_str());
               ++nFile;
            } // if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

            if (ftmp) {
               ftmp->Close();
            } // if (ftmp)
         }    // if(file.find(".miniDst.root") != std::string::npos)
      }       // while (getline(inputStream, file))

      LOG_INFO << " Total " << nFile << " files have been read in. " << endm;
   } // if(dirFile.find(".list") != std::string::npos || dirFile.find(".lis" != string::npos))
   else if (dirFile.find(".miniDst.root") != std::string::npos || dirFile.find(".MiniDst.root") != std::string::npos) {
      fChain->Add(dirFile.c_str());
   } else {
      LOG_WARN << " No good input file to read ... extentions must be: "
               << " fname.lis(t), or fname.miniDst.root (fname.MiniDst.root)" << endm;
   }

   if (fChain) {
      setBranchAddresses(fChain);
      fChain->SetCacheSize(50e6);
      fChain->AddBranchToCache("*");
      fMiniDst->set(fMiniArrays);
   }
}

//_________________
Bool_t MpdMiniDstReader::readMiniEvent(Long64_t iEvent __attribute__((unused)))
{

   Int_t fStatusRead = true; // true - okay, false - nothing to read

   if (!fChain) {
      LOG_WARN << " No input files ... ! EXIT" << endm;
      fStatusRead = false;
      return fStatusRead;
   }

   Int_t bytes   = fChain->GetEntry(fEventCounter++);
   Int_t nCycles = 0;
   while (bytes <= 0) {
      if (fEventCounter >= fChain->GetEntriesFast()) {
      }

      LOG_WARN << "Encountered invalid entry or I/O error while reading event " << fEventCounter << " from \""
               << fChain->GetName() << "\" input tree\n";
      bytes = fChain->GetEntry(fEventCounter++);
      nCycles++;
      LOG_WARN << "Not input has been found for: " << nCycles << " times" << endm;
      if (nCycles >= 10) {
         LOG_ERROR << "Terminating MpdMiniDstReader::readMiniEvent(Long64_t) after " << nCycles << " times!" << endm;
         fStatusRead = false;
         break;
      }
   }
   return fStatusRead;
}
