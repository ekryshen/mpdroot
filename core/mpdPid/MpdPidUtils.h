#ifndef MPD_PID_UTILS_H
#define MPD_PID_UTILS_H

#include <TH1.h>
#include <TH2.h>

#define MASS_MU 0.1057
#define MASS_PI 0.1396
#define MASS_PI2 0.0195
#define MASS_KA 0.4937
#define MASS_PR 0.9383
#define MASS_DE 1.876
#define MASS_TR 2.8094
#define MASS_HE3 1.4047
#define MASS_HE4 1.863

#define PDG_DEUTERON 1000010020
#define PDG_TRITON 1000010030
#define PDG_HE3 1000020030
#define PDG_HE4 1000020040

namespace MpdPidUtils {

const Int_t kNSpecies    = 9;  /// Number of default particle species recognized by the PID
const Int_t kNPosSpecies = 9;  /// Number of positively charged particle species recognized by the PID
const Int_t kNNegSpecies = 5;  /// Number of negatively charged particle species recognized by the PID
const Int_t nQAHists     = 40; /// Max number of dE/dx QA hists

enum ePartType {
   kElectron = 0,
   kMuon     = 1,

   kPion   = 2,
   kKaon   = 3,
   kProton = 4,

   kDeuteron = 5,
   kTriton   = 6,
   kHe3      = 7,
   kHe4      = 8,
   kUnknown  = 9
};

enum ePartCharge { kPos = 0, kNeg = 1 };

enum eTrackingState { kHP = 0, kCF = 1, kCFHM = 2 };

const char *const cParticleName[kNSpecies + 1] = {"electron", "muon",   "pion",     "kaon",  "proton",
                                                  "deuteron", "triton", "helium-3", "alpha", "unknown"};

const char *const cParticleShortName[kNSpecies + 1] = {"e", "mu", "pi", "K", "p", "d", "t", "he3", "alpha", "unknown"};

struct dEdXStruct_t {
   TH1D *dEdXPart[nQAHists];
   TH2D *BetheBlochHist;
   Int_t ibeg;
   Int_t iend;
};
}; // namespace MpdPidUtils
#endif
