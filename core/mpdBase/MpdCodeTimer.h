#ifndef MPDCODETIMER_H
#define MPDCODETIMER_H

/// \ingroup rec
/// \class MpdCodeTimer
/// \brief MPD code timing utility class.
///
/// \author Mariya Ilieva, LHEP JINR Dubna

#include "FairRun.h"
#include "FairTask.h"
#include <TString.h>
//#include <map>
//#include "TMVA/Timer.h"
#include <TMap.h>
class TStopwatch;

class MpdCodeTimer : public FairTask
{
 public:
  // KG public cpnstructor for TStreamer reading
  MpdCodeTimer(const char *name="MpdCodeTimer", const char *title="MPD Task"); ///< Ctor
  static MpdCodeTimer* Instance(); ///< get singleton instance
  static MpdCodeTimer* Instance(const char *name, const char *title="MPD Task"); ///< get singleton instance
  static MpdCodeTimer* Active() { return (MpdCodeTimer*) FairRun::Instance()->GetTask("Code timer"); } ///< get singleton instance
  virtual void Exec(Option_t * option);

  //void Init(Int_t counter);

  void Reset();
  void Register();
  void Start(const TString classname, const TString method);
  void Stop(const TString classname, const TString method);
  void Print();

 protected:

  virtual InitStatus Init();
  virtual InitStatus ReInit();
  virtual void Finish();
  virtual ~MpdCodeTimer(); ///< Destructor

 private:
  // automatic deleting when application exit
  static void DestroyInstance (){
      if (fgCT)
          delete fgCT;
  }
  void Print(TStopwatch *sw, Option_t *opt); // print info from TStopwatch 

  static MpdCodeTimer* fgCT; //! pointer to Singleton instance, may be ! required to exclude recursion

 private:
  TMap fTimeMap; // class-method map

  ClassDef(MpdCodeTimer,1);

};
#endif
