// -------------------------------------------------------------------------
// -----              MpdTPythia8Generator header file                 -----
// -----           Created 10/04/2016  by A.Zinchenko                   -----
// -------------------------------------------------------------------------

/** MpdTPythia8Generator.h
 ** @author
 **
 ** The MpdTPythia8Generator generates particles from TPythia8 output tree.
 ** Derived from FairGenerator.
 */

#ifndef MPDTPYTHIA8GENERATOR_H
#define MPDTPYTHIA8GENERATOR_H

#include "FairGenerator.h"
//#include "Pythia8/Pythia.h"

using namespace std;
// using namespace Pythia8;

class FairPrimaryGenerator;
class TClonesArray;
class TFile;
class TTree;

class MpdTPythia8Generator : public FairGenerator {
public:
   /** Default constructor. */
   MpdTPythia8Generator();

   /** Constructor with file name, tree name, particle branch name
    **/
   // MpdTPythia8Generator(TString fileName, TString treeName = "T", TString branchName = "entry");
   MpdTPythia8Generator(TString fileName, TString treeName = "T", TString branchName = "particles");

   /** Destructor */
   virtual ~MpdTPythia8Generator();

   /** Creates an event with given type and multiplicity.
    **@param primGen  pointer to the CbmPrimaryGenerator
    */
   virtual Bool_t ReadEvent(FairPrimaryGenerator *primGen);

   void SkipEvents(Int_t nEvents) { fNevents = nEvents; }

private:
   TFile        *fInputFile; //! input file
   TTree        *fTree;      //! input tree
   TClonesArray *fParticles; //! particles
   Int_t         fNevents;   //! event number
   Int_t         fEdit;      //! edit flag (as in Pythia)

   void AddParticlesToPdgDataBase();

   ClassDef(MpdTPythia8Generator, 1);
};

#endif
