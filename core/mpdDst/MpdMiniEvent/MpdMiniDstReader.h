/**
 * \class MpdMiniDstReader
 * \brief Allows to read miniDst file(s)
 *
 * This class allows to read miniDst.root file or a list of files
 * that contain miniDst and sets up pointers to the miniDst, and
 * certain TClonesArrays that keep Event, Track, BTofHit, etc...
 * One can also turn on or off certain branches using the
 * SetStatus method.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * email: nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniDstReader_h
#define MpdMiniDstReader_h

// ROOT headers
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"

// MiniDst headers
#include "MpdMiniDst.h"
#include "MpdMiniEvent.h"
#include "MpdMiniArrays.h"

//_________________
class MpdMiniDstReader : public TObject {

public:
   /// Constructor that takes either miniDst file or file that
   /// contains a list of miniDst.root files
   MpdMiniDstReader(const Char_t *inFileName);
   /// Destructor
   ~MpdMiniDstReader();

   /// Return a pointer to miniDst (return NULL if no dst is found)
   MpdMiniDst *miniDst() { return fMiniDst; }
   /// Return pointer to the chain of .miniDst.root files
   TChain *chain() { return fChain; }
   /// Return pointer to the current TTree
   TTree *tree() { return fTree; }

   /// Set enable/disable branch matching when reading miniDst
   void SetStatus(const Char_t *branchNameRegex, Int_t enable);

   /// Calls openRead()
   void Init();
   /// Read next event in the chain
   Bool_t readMiniEvent(Long64_t iEvent);
   /// Close files and finilize
   void Finish();

private:
   /// Name of the inputfile (or of the inputfiles.list)
   TString mInputFileName;

   /// Turn off streamers
   void streamerOff();

   /// Create TClonesArray of mini classes and set them to miniDst
   void createArrays();
   /// Clear all TClonesArrays with miniDst classes
   void clearArrays();
   /// Set adresses of miniArrays and their statuses (enable/disable) to chain
   void setBranchAddresses(TChain *chain);

   /// Pointer to the input/output miniDst structure
   MpdMiniDst *fMiniDst;
   /// Pointer to the chain
   TChain *fChain;
   /// Pointer to the current tree
   TTree *fTree;

   /// Event counter
   Int_t fEventCounter;

   /// Pointers to mini arrays
   TClonesArray *fMiniArrays[MpdMiniArrays::NAllMiniArrays];
   /// Status of mini arrays
   Char_t fStatusArrays[MpdMiniArrays::NAllMiniArrays];

   ClassDef(MpdMiniDstReader, 0);
};

#endif // #define MpdMiniDstReader_h
