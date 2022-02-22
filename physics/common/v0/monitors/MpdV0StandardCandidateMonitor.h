/*
 * MpdV0StandardCandidateMonitor.h
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDCANDIDATEMONITOR_H_
#define MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDCANDIDATEMONITOR_H_

#include <Rtypes.h>
#include <RtypesCore.h>

#include "MpdV0Namespace.h"
#include "MpdV0CandidateMonitor.h"
#include "MpdV0Monitor.h"

class TH2D;

class MpdV0StandardCandidateMonitor : public MpdV0CandidateMonitor {
protected:
   MpdV0::EParticleType fType;

public:
   MpdV0StandardCandidateMonitor(MpdV0::EParticleType type = MpdV0::EParticleType::kPdgHypo);
   MpdV0StandardCandidateMonitor(const MpdV0StandardCandidateMonitor &other) = default;
   void SetMinvAxis(Int_t bins, Double_t min, Double_t max) { SetXaxis1d(0, bins, min, max); }
   void SetCosAxis(Int_t bins) { SetXaxis1d(1, bins, -1, 1); }
   void SetArmenterosAlphaAxis(Int_t bins, Double_t min, Double_t max) { SetXaxis2d(0, bins, min, max); }
   void SetArmenterosPtAxis(Int_t bins, Double_t min, Double_t max) { SetYaxis2d(0, bins, min, max); }
   void SetDecayLenghtMonitor(Int_t bins, Double_t min, Double_t max) { SetXaxis1d(2, bins, min, max); };
   void SetDau1to2(Int_t bins, Double_t min, Double_t max) { SetXaxis1d(3, bins, min, max); };
   void Init();
   void SetType(MpdV0::EParticleType type) { fType = type; };
   void Fill(const MpdV0Particle &particle, Bool_t status);
   MpdV0CandidateMonitor *MakeCopy() const { return new MpdV0StandardCandidateMonitor(*this); }
   virtual ~MpdV0StandardCandidateMonitor();
   ClassDef(MpdV0StandardCandidateMonitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0STANDARDCANDIDATEMONITOR_H_ */
