#ifndef _TPC_COMMON_H_
#define _TPC_COMMON_H_

#include <Rtypes.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TH2.h>
#include <TH3.h>

#include <iostream>
#include <vector>

/************************************************************************/
/* Constants                                                            */
/************************************************************************/

const UInt_t kTPC_SECTORS_NUMBER  = 24;
const UInt_t kTPC_ROWS_IN_SECTOR  = 50;
const UInt_t kTPC_PADS_IN_ROW     = 200;
const UInt_t kTPC_TIMEBINS_IN_PAD = 512;
// Drift velocity
const Double_t kTPC_DRIFT_SPEED = 5.55 * 10 / 1000.; // mm/ns
                                                     // in samples/ns
const Double_t kTPC_READOUT_FREQUENCY = (12.5 * 1000. * 1000.) / (1000. * 1000. * 1000.);
/************************************************************************/
/* Macros                                                               */
/************************************************************************/

#define COUNTOF(arr) (unsigned long)(sizeof(arr) / sizeof(arr[0]))

/************************************************************************/
/* Inline functions                                                     */
/************************************************************************/

inline double Deg2Rad(double deg)
{
   return deg * TMath::Pi() / 180.;
}

inline double Rad2Deg(double rad)
{
   return rad * 180. / TMath::Pi();
}

#endif // _TPC_COMMON_H_
