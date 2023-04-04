#ifndef MPDANALYSISEVENTPLANE_H
#define MPDANALYSISEVENTPLANE_H

#include "FairEventHeader.h"
#include "FairMCEventHeader.h"

class MpdAnalysisEventPlane {
   // Class to contain basic information about event plane in MPD
public:
   MpdAnalysisEventPlane();
   virtual ~MpdAnalysisEventPlane(); // Destructor

   // Getters
   float GetPhiEP_FHCal_F_all() { return fPhiEP_FHCal_F_all; }
   float GetPhiEP_FHCal_N_all() { return fPhiEP_FHCal_N_all; }
   float GetPhiEP_FHCal_S_all() { return fPhiEP_FHCal_S_all; }
   float GetPhiEP_TPC_N_all() { return fPhiEP_TPC_N_all; }
   float GetPhiEP_TPC_S_all() { return fPhiEP_TPC_S_all; }

   // Setters
   void SetPhiEP_FHCal_F_all(float ep) { fPhiEP_FHCal_F_all = ep; }
   void SetPhiEP_FHCal_N_all(float ep) { fPhiEP_FHCal_N_all = ep; }
   void SetPhiEP_FHCal_S_all(float ep) { fPhiEP_FHCal_S_all = ep; }
   void SetPhiEP_TPC_N_all(float ep) { fPhiEP_TPC_N_all = ep; }
   void SetPhiEP_TPC_S_all(float ep) { fPhiEP_TPC_S_all = ep; }

private:
   float fPhiEP_FHCal_F_all = -9999.; // V event plane angle from all modules in both FHCals
   float fPhiEP_FHCal_N_all = -9999.; // V event plane angle from all modules in north FHCal (eta<0)
   float fPhiEP_FHCal_S_all = -9999.; // V event plane angle from all modules in south FHCal (eta>0)
   float fPhiEP_TPC_N_all   = -9999.; // V event plane angle from all tracks in north part of TPC (eta<0)
   float fPhiEP_TPC_S_all   = -9999.; // V event plane angle from all tracks in south part of TPC (eta>0)

   ClassDef(MpdAnalysisEventPlane, 0);
};

#endif