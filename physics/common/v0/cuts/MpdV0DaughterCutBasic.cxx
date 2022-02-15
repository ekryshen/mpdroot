/*
 * MpdV0DaughterCutBasic.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0DaughterCutBasic.h"

#include <TMathBase.h>
#include <TVector3.h>
#include <iostream>

#include "MpdTrack.h"
#include "MpdMiniTrack.h"

MpdV0DaughterCutBasic::MpdV0DaughterCutBasic()
   : MpdV0DaughterCut(), fCharge(1), fDcaXY2(0), fDcaZ(0), fSigmaLow(-2), fSigmaHigh(2), fNHitsTpcLow(0),
     fNHitsTpcHigh(56), fSigmaType(MpdV0::ESigmaType::kPionSigma)
{
}

void MpdV0DaughterCutBasic::SetNTpcHitsCut(Double_t low, Double_t high)
{
   fNHitsTpcLow  = low;
   fNHitsTpcHigh = high;
}

void MpdV0DaughterCutBasic::SetDcaMinCut(Double_t dcaXY, Double_t dcaZ)
{
   fDcaXY2 = dcaXY * dcaXY;
   fDcaZ   = dcaZ;
}

void MpdV0DaughterCutBasic::SetSigmaCut(Double_t low, Double_t high, MpdV0::ESigmaType type)
{
   fSigmaLow  = low;
   fSigmaHigh = high;
   fSigmaType = type;
}

MpdV0DaughterCutBasic::~MpdV0DaughterCutBasic() {}

Bool_t MpdV0DaughterCutBasic::PassDstTrack(MpdTrack &track) const
{
   if (track.GetCharge() != fCharge) return kFALSE;
   if (track.GetNofHits() < fNHitsTpcLow) return kFALSE;
   if (track.GetNofHits() > fNHitsTpcHigh) return kFALSE;
   switch (fSigmaType) {
   case MpdV0::ESigmaType::kPionSigma: {
      if (track.GetNSigmaPion() < fSigmaLow) return kFALSE;
      if (track.GetNSigmaPion() > fSigmaHigh) return kFALSE;
   } break;
   case MpdV0::ESigmaType::kKaonSigma: {
      if (track.GetNSigmaKaon() < fSigmaLow) return kFALSE;
      if (track.GetNSigmaKaon() > fSigmaHigh) return kFALSE;
   } break;
   case MpdV0::ESigmaType::kProtonSigma: {
      if (track.GetNSigmaProton() < fSigmaLow) return kFALSE;
      if (track.GetNSigmaProton() > fSigmaHigh) return kFALSE;
   } break;
   }
   if (TMath::Abs(track.GetDCAGlobalZ()) < fDcaZ) return kFALSE;
   Double_t dcaX = track.GetDCAGlobalX();
   Double_t dcaY = track.GetDCAGlobalY();
   if (dcaX * dcaX + dcaY * dcaY < fDcaXY2) return kFALSE;
   return kTRUE;
}

Bool_t MpdV0DaughterCutBasic::PassMiniDstTrack(MpdMiniTrack &track) const
{

   if (track.charge() != fCharge) return kFALSE;
   if (track.nHits() < fNHitsTpcLow) return kFALSE;
   if (track.nHits() > fNHitsTpcHigh) return kFALSE;
   switch (fSigmaType) {
   case MpdV0::ESigmaType::kPionSigma: {
      if (track.nSigmaPion() < fSigmaLow) return kFALSE;
      if (track.nSigmaPion() > fSigmaHigh) return kFALSE;
   } break;
   case MpdV0::ESigmaType::kKaonSigma: {
      if (track.nSigmaKaon() < fSigmaLow) return kFALSE;
      if (track.nSigmaKaon() > fSigmaHigh) return kFALSE;
   } break;
   case MpdV0::ESigmaType::kProtonSigma: {
      if (track.nSigmaProton() < fSigmaLow) return kFALSE;
      if (track.nSigmaProton() > fSigmaHigh) return kFALSE;
   } break;
   }
   Double_t dcaZ = track.origin().z() - GetVertex().Z();
   Double_t dcaX = track.origin().x() - GetVertex().X();
   Double_t dcaY = track.origin().y() - GetVertex().Y();

   if (TMath::Abs(dcaZ) < fDcaZ) return kFALSE;
   if (dcaX * dcaX + dcaY * dcaY < fDcaXY2) return kFALSE;
   return kTRUE;
}
