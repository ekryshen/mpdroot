//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Tpc Clusterizer reads in MCPoints and produces primary electrons
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

#ifndef MPDTPCCLUSTERIZERTASK_HH
#define MPDTPCCLUSTERIZERTASK_HH

// Base Class Headers ----------------
#include "FairTask.h"

// Collaborating Class Headers -------

// Collaborating Class Declarations --
class TClonesArray;
class TpcGas;

class MpdTpcClusterizerTask : public FairTask {
public:
   static const Int_t fgkNsec2; // number of readout sectors (12 * 2)

   // Constructors/Destructors ---------
   MpdTpcClusterizerTask();
   ~MpdTpcClusterizerTask();

   // Operators

   // Accessors -----------------------

   // Modifiers -----------------------
   void SetPointBranchName(const TString &name) { fPointBranchName = name; }
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

   // Operations ----------------------

   virtual InitStatus Init();
   void               FinishTask();

   virtual void Exec(Option_t *opt);
   // virtual void Clear(Option_t* opt);

private:
   // Private Data Members ------------
   TString        fPointBranchName;
   TClonesArray * fPointArray;
   TClonesArray **fPrimArray;

   TpcGas *fGas;

   Bool_t fPersistence;

   // Private Methods -----------------
   void ClearData(Int_t isec); // save and clear data structures

public:
   ClassDef(MpdTpcClusterizerTask, 1);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
