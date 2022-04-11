// -------------------------------------------------------------------------
// -----          MpdGeneralGenerator header file                      -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

/**  MpdGeneralGenerator.h
 *@author Dmitry Sosnov (05.2015) <dmitry.e.sosnov@gmail.com>
 *
 * The MpdGeneralGenerator is interface
 * between TGenerator and FairGenerator classes.
 * Derived from FairGenerator.
 **/

#ifndef MPD_GENERAL_GENERATOR_H
#define MPD_GENERAL_GENERATOR_H

#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"
#include "stdio.h"

// class FairPrimaryGenerator;

class MpdGeneralGenerator : public FairGenerator {
   ClassDef(MpdGeneralGenerator, 1);

public:
   /** Default constructor. **/
   MpdGeneralGenerator();
   //         template <class T> MpdGeneralGenerator();

   MpdGeneralGenerator(TGenerator *gen);
   MpdGeneralGenerator(TGenerator &gen);

   /** Destructor **/
   virtual ~MpdGeneralGenerator() { delete fGenerator; };

   /** Modifiers **/
   void        SetGenerator(TGenerator *gen) { fGenerator = gen; }
   TGenerator *GetGenerator(void) { return fGenerator; }

   /** Initializer **/
   Bool_t Init();
   void   SetDebug(bool d = kTRUE) { fDebug = d; }

   virtual Bool_t ReadEvent(FairPrimaryGenerator *primGen);

private:
   /** Copy constructor */
   MpdGeneralGenerator(const MpdGeneralGenerator &G)
      : FairGenerator() {} // TODO Missed s because cannot build. Find error
   /** Assignment operator */
   MpdGeneralGenerator &operator=(const MpdGeneralGenerator &) { return *this; }

private:
   TGenerator *fGenerator; //!
   bool        fDebug;
};

#endif
