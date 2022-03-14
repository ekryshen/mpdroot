/*
 * File:   newfile.h
 * Author: niki
 *
 * Created on September 6, 2013, 5:51 PM
 */
#ifndef MPD_V0_CUTS_H
#define MPD_V0_CUTS_H

#include <TROOT.h>

namespace MpdV0Cuts {

// Mother cuts
extern Double_t V0Chi2Cut;
extern Double_t RminCut;
extern Double_t RmaxCut;
extern Double_t cosaCut;
extern Double_t minmassCut;
extern Double_t maxmassCut;
// extern Double_t IPLambdaCut;
// extern Double_t PtlambdaCut;
extern Double_t PminCut;

// Daughter cuts
extern Double_t NofHitsCut;
extern Double_t ProbCut;
// extern Double_t TrChi2posCut;
// extern Double_t TrChi2negCut;
extern Double_t protonChi2vmaxCut;
extern Double_t protonChi2vminCut;
extern Double_t pionChi2vmaxCut;
extern Double_t pionChi2vminCut;
extern Double_t IPposCut;
extern Double_t IPnegCut;
extern Double_t protonPtminCut;
extern Double_t protonPminCut;
extern Double_t pionPtminCut;
extern Double_t pionPminCut;
extern Double_t protonPtmaxCut;
extern Double_t protonPmaxCut;
extern Double_t pionPtmaxCut;
extern Double_t pionPmaxCut;
extern Double_t protonIPCut;
extern Double_t pionIPCut;
extern Double_t DCAV0Cut;

extern Int_t protonPDG; // proton+
extern Int_t pionPDG;   // pi+
extern Int_t lambdaPDG;
} // namespace MpdV0Cuts
#endif /* MPD_V0CUTS_H */
