////////////////////////////////////////////
//    BmnEventMonitor.h
//    BM@N Event Monitoring Task declaration
//    Konstantin Gertsenberger
//    Created: Oct. 11 2016
////////////////////////////////////////////

#ifndef BMNEVENTMONITOR_H
#define BMNEVENTMONITOR_H

#define _GLIBCXX_USE_CXX11_ABI 0

#ifndef __MAKECINT__
#include <emon/EventIterator.h>
#include <owl/timer.h>
#include <cmdl/cmdargs.h>
#include <ipc/core.h>
#endif // __MAKECINT__

#include <FairTask.h>
#include <TObject.h>
#include <TString.h>

#include <vector>

class BmnEventMonitor : public FairTask {
public:
   /** Default constructor **/
   BmnEventMonitor();

   /** Standard constructor
    *@param partition_name  partition to work in
    *@param sampling_type   sampling address - sampler type
    *@param sampling_names  sampling address - sampler names
    *@param max_events      maximum event to receive from Event Sampler
    *@param verbose <table>
     <tr><td>0</td><td>print nothing (default)</td></tr>
     <tr><td>1</td><td>print event number and event size?</td></tr>
     <tr><td>2</td><td>print event number, event size, and event data</td></tr>
     </table>
    **/
   BmnEventMonitor(TString partition_name, TString sampling_type, TString sampling_names, int max_events = 100,
                   Int_t verbose = 0);

   /** Destructor **/
   virtual ~BmnEventMonitor();

   /** Set parameters for event receiving. **/
   void SetPartitionName(TString partition_name) { strPartitionName = partition_name; }
   void SetSamplingType(TString sampling_type) { strSamplingType = sampling_type; }
   void SetSamplingName(TString sampling_name) { strSamplingName = sampling_name; }

   void SetMaxEvents(int max_events) { iMaxEvents = max_events; }

   /** Set verbosity level. For this task and all of the subtasks. **/
   virtual void SetVerbose(Int_t verbose) { iVerbose = verbose; }

   /** Executed task for RunMonitor **/
   virtual void Exec(Option_t *option);

   void Reset() {}

protected:
   virtual InitStatus Init();
   /** Action after each event**/
   virtual void Finish();

   // partition name to work in
   TString strPartitionName;
   // address type
   TString strSamplingType;
   // sampler name
   TString strSamplingName;

   // number of events to retrieve (default 10000), -1 means work forever
   Int_t iMaxEvents;

   // use event dispersion
   Bool_t isDispersion;
   // delay in ms before another connection attempt (default 10000)
   Int_t iWait;
   // number of repetitions for the test (default 1), minimum value is 1
   Int_t iRepetitions;
   // maximum number of events in the buffer (default 1000)), minimum value is 1
   Int_t iBufferSize;
   // 1 if monitor should retrieve events asynchronous (default 0)
   Int_t iAsync;
   // timeout for the nextEvent call in synchronous mode (default 10000 ms)
   Int_t iTimeout;

   // Verbosity level
   Int_t iVerbose;

   // array for gem digits to transfer it to the next tasks
   TClonesArray *fGemDigits; //!
   int           iEventCount;

   ClassDef(BmnEventMonitor, 1);
};

#endif
