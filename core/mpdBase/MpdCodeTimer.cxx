/// \class MpdCodeTimer
///
/// brief MPD code timing utility class.
/// \author Maria Ilieva (LHEP, JINR, Dubna)

#include "MpdCodeTimer.h"
#include "FairTask.h"
#include <TStopwatch.h>
#include <TObjString.h>
//#include "mkl.h"

//#include "mkl_lapacke.h"

#ifdef _OPENMP
#include "omp.h"
omp_lock_t geoManagerLock;
#endif
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;
//using std::pair;

MpdCodeTimer* MpdCodeTimer::fgCT = 0x0;

//__________________________________________________________________________
MpdCodeTimer::MpdCodeTimer(const char *name, const char *title)
  : FairTask(name)
{
  /// Constructor
  fgCT = this;
}

//__________________________________________________________________________
MpdCodeTimer* MpdCodeTimer::Instance()
{
  /// Get pointer to the code timer singleton object
  if (!fgCT){
    fgCT = new MpdCodeTimer;
    #ifdef _OPENMP
        omp_init_lock(&geoManagerLock);
    #endif

    // automatic destroy is supposed
    std::atexit(DestroyInstance);
  }
  return fgCT;
}

//__________________________________________________________________________
MpdCodeTimer* MpdCodeTimer::Instance(const char *name, const char *title)
{
  /// Get pointer to the Kalman filter reconstructor singleton object
  if (!fgCT){
    fgCT = new MpdCodeTimer(name, title);
    #ifdef _OPENMP
        omp_init_lock(&geoManagerLock);
    #endif

    // automatic destroy is supposed
    std::atexit(DestroyInstance);
  }
  return fgCT;
}

//__________________________________________________________________________
MpdCodeTimer::~MpdCodeTimer() 
{
  /// Destructor

  fTimeMap.Clear();
  fgCT = NULL;
}

//__________________________________________________________________________
InitStatus MpdCodeTimer::Init() {

  fTimeMap.SetOwner(kTRUE);
  fTimeMap.SetOwnerValue(kTRUE);
  return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdCodeTimer::ReInit() 
{
  return kSUCCESS;
}

//__________________________________________________________________________
void MpdCodeTimer::Reset() 
{
  ///
}

//__________________________________________________________________________
void MpdCodeTimer::Register() 
{
  ///
}

//__________________________________________________________________________
void MpdCodeTimer::Finish() 
{
  ///

  Print();
}

//__________________________________________________________________________
void MpdCodeTimer::Exec(Option_t * option) 
{
  ///
}

//__________________________________________________________________________
void MpdCodeTimer::Start(const TString classname, const TString method) 
{
  /// Start the stopwatch. Create new if it does not exist.

  TString classMeth = classname;
  classMeth += "::";
  classMeth += method;

  if (fTimeMap.FindObject(classMeth) == 0x0) {
    // Create new
    TStopwatch *sw = new TStopwatch;
    fTimeMap.Add(new TObjString(classMeth),sw);
    sw->Start();
  } else {
    (static_cast<TStopwatch*>(fTimeMap(classMeth)))->Start(kFALSE);
  }
}

//__________________________________________________________________________
void MpdCodeTimer::Stop(const TString classname, const TString method) 
{
  /// Stop timer

  if (fgCT == NULL) return;
  TString classMeth = classname;
  classMeth += "::";
  classMeth += method;

  if (fTimeMap.FindObject(classMeth) == 0x0) return;
  (static_cast<TStopwatch*>(fTimeMap(classMeth)))->Stop();
}

//__________________________________________________________________________
void MpdCodeTimer::Print() 
{
  ///
  //TFile = new File("time.log","");
  TIterator* it = fTimeMap.MakeIterator();
  TObjString *classMeth = NULL;
  TObjArray tmp;
  while ( (classMeth = static_cast<TObjString*>(it->Next()) ) ) {
    tmp.Add(classMeth);
  }
  tmp.Sort();
  Int_t nEntries = tmp.GetEntriesFast();

  cout << " \n *** MpdCodeTimer info *** " << endl;
  for (Int_t i = 0; i < nEntries; ++i) {
    classMeth = static_cast<TObjString*> (tmp.UncheckedAt(i));
    cout << classMeth->String() << " ";
    //(static_cast<TStopwatch*>(fTimeMap.GetValue(classMeth)))->Print();
    Print(static_cast<TStopwatch*>(fTimeMap.GetValue(classMeth)),"");
  }
    
}

//__________________________________________________________________________
void MpdCodeTimer::Print(TStopwatch *sw, Option_t *opt) 
{
  /// Print TStopwatch information (adapted from ROOT)

  Double_t  realt = const_cast<TStopwatch*>(sw)->RealTime();
  Double_t  cput  = const_cast<TStopwatch*>(sw)->CpuTime();

  Int_t  hours = Int_t(realt / 3600);
  realt -= hours * 3600;
  Int_t  min   = Int_t(realt / 60);
  realt -= min * 60;
  Int_t  sec   = Int_t(realt);

  if (realt < 0) realt = 0;
  if (cput  < 0) cput  = 0;
  
  if (opt && *opt == 'm') {
    Printf("Real time %d:%02d:%06.3f, CP time %.3f, %d slices", hours, min, realt, cput, sw->Counter());
  } else if (opt && *opt == 'u') {
    Printf("Real time %d:%02d:%09.6f, CP time %.3f, %d slices", hours, min, realt, cput, sw->Counter());
  } else {
    Printf("Real time %d:%02d:%02d, CP time %.3f, %d slices", hours, min, sec, cput, sw->Counter());
  }
}

ClassImp(MpdCodeTimer)
