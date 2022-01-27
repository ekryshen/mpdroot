//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Tpc Drifter reads in TpcPrimaryClusters and produces drifted electrons
//
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - adapted for MPD from PANDARoot
//
//
//-----------------------------------------------------------

#ifndef MPDTPCDRIFTTASK_HH
#define MPDTPCDRIFTTASK_HH

// Base Class Headers ----------------
#include "FairTask.h"

// Collaborating Class Headers -------

// Collaborating Class Declarations --
class TClonesArray;
class TpcGas;

class MpdTpcDriftTask : public FairTask {
public:
   // Constructors/Destructors ---------
   MpdTpcDriftTask();
   ~MpdTpcDriftTask();

   // Operators

   // Accessors -----------------------

   // Modifiers -----------------------
   void SetPrimBranchName(const TString &name) { fPrimBranchName = name; }
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }
   void SetAttach(Bool_t opt = kTRUE) { fAttach = opt; }
   void SetDiffuse(Bool_t opt = kTRUE) { fDiffuse = opt; }
   void SetDistort(Bool_t opt = kTRUE) { fDistort = opt; }

   // Operations ----------------------

   virtual InitStatus Init();

   virtual void Exec(Option_t *opt);
   virtual void Clear(Option_t *opt);

private:
   // Private Data Members ------------
   const Int_t    fkNsec2;
   TString        fPrimBranchName;
   TClonesArray **fPrimArray;
   TClonesArray **fDriftedArray;

   TpcGas * fGas;
   Double_t fzGem; // get from Geom!

   Bool_t fPersistence;
   Bool_t fAttach;
   Bool_t fDiffuse;
   Bool_t fDistort;
   Int_t  fSec; // current sector number

   // Private Methods -----------------

public:
   ClassDef(MpdTpcDriftTask, 1);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
