/*
 * File:   MpdSmashGenerator.h
 * Author: Igor Altsybeev et al.
 *
 * Created on 17/07/2020
 */

#ifndef MPDSMASHGENERATOR_H_H
#define MPDSMASHGENERATOR_H_H

#include "FairGenerator.h"
#include "TTree.h"
#include "TRandom3.h"
#include <fstream>

const UInt_t MAX_N_PART = 1e4;

class TVirtualMCStack;
class FairPrimaryGenerator;
class TFile;
// class TTree;

class MpdSmashGenerator : public FairGenerator {
public:
   /** Default constructor without arguments should not be used. **/
   MpdSmashGenerator();

   /** Standard constructor.
    * @param fileName The input file name
    **/
   MpdSmashGenerator(TString fileName);

   /** Destructor. **/
   ~MpdSmashGenerator();

   void SetRandomRP(bool flag_isRandomRP) { fIsRandomRP = flag_isRandomRP; }

   int GetNeventsInTree() { return fInputTree ? fInputTree->GetEntries() : -1; }

   Bool_t ReadEvent(FairPrimaryGenerator *primGen);

   void SkipEvents(Int_t ev) { fEventNumber = ev; }

private:
   TFile  *fInputFile; //!  Input file
   TString fFileName;  //!  Input file name
   TTree  *fInputTree; //! Input SMASH tree
   //    TChain *fInputTree;

   // Vars to be read from SMASH tree:
   Double_t fImpPar;             //!
   Int_t    fNpart;              //!
   Bool_t   fEmptyEv;            //!
   Double_t fPx[MAX_N_PART];     //!
   Double_t fPy[MAX_N_PART];     //!
   Double_t fPz[MAX_N_PART];     //!
   Int_t    fPID[MAX_N_PART];    //!
   Int_t    fCharge[MAX_N_PART]; //!

   //
   Int_t fEventNumber; //!

   Bool_t   fIsRandomRP; //! // random/fixed reaction plane
   Double_t fPsiRP;      //! // reaction plane angle

   TRandom3 *fRandom; //!

   MpdSmashGenerator(const MpdSmashGenerator &);
   MpdSmashGenerator &operator=(const MpdSmashGenerator &);

   ClassDef(MpdSmashGenerator, 2);
};

#endif /* MPDSMASHGENERATOR_H */
