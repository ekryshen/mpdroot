// -------------------------------------------------------------------------
// -----              MpdHistoGenerator header file                    -----
// -----           Created 03/06/2014  by M.Ilieva, A.Zinchenko        -----
// -------------------------------------------------------------------------

/** MpdHistoGenerator.h
 ** @author
 **
 ** The MpdHistoGenerator generates particles according
 ** to Eta-Pt histogram.
 ** Derived from FairGenerator.
 */

#ifndef MPDHISTOGENERATOR_H
#define MPDHISTOGENERATOR_H

#include "FairGenerator.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;

class FairPrimaryGenerator;

class MpdHistoGenerator : public FairGenerator {
public:
   /** Default constructor. */
   MpdHistoGenerator();

   /** Constructor with PDG-ID, multiplicity
    **@param pdgid Particle type (PDG encoding)
    **@param mult  Multiplicity (default is 1)
    **@param yield
    **/
   MpdHistoGenerator(Int_t pdgid, Int_t mult = 1, Double_t yield = -1.0);

   /** Destructor */
   virtual ~MpdHistoGenerator();

   /** Creates an event with given type and multiplicity.
    **@param primGen  pointer to the CbmPrimaryGenerator
    */
   virtual Bool_t ReadEvent(FairPrimaryGenerator *primGen);
   Bool_t         Init();
   void           SetFileName(TString fileName) { fFileName = fileName; }

private:
   Int_t    fPdgCode;  //! particle PDG code
   Int_t    fMult;     //! particle multiplicity
   TH2D    *fHist;     //! Pointer to the histogram
   Double_t fYield;    ///< yield (particle / event)
   TString  fFileName; ///< input file name

   ClassDef(MpdHistoGenerator, 0);
};

#endif
