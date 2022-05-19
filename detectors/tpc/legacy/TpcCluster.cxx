#include "TpcCluster.h"
#include "TpcSequence.h"
#include "TpcTimeBin.h"
#include "TpcPad.h"

ClassImp(TpcCluster);

TpcCluster::TpcCluster(void)
   : fPadRow(NULL), fSequences(NULL), fRegion(NULL), fLoPad(0), fHiPad(0), fLoBin(0), fHiBin(0)
{
   fSequences = new TObjArray(1);
   fRegion    = new TObjArray(1);
}

TpcCluster::~TpcCluster(void)
{
   if (fSequences) {
      delete fSequences;
   }
   if (fRegion) delete fRegion;
}

Bool_t TpcCluster::RecalcBounds()
{
   if (!fSequences) {
      Error(__FUNCTION__, "fSequences == NULL");
      return kFALSE;
   }

   if (fSequences->IsEmpty()) {
      Error(__FUNCTION__, "Cluster has no sequences");
      fLoPad = fHiPad = fLoBin = fHiBin = 0;
      return kFALSE;
   }

   fLoPad = fLoBin = kMaxInt;
   fHiPad = fHiBin = kMinInt;
   for (Int_t i = 0; i < fSequences->GetEntriesFast(); ++i) {
      TpcSequence *seq = (TpcSequence *)fSequences->At(i);
      fLoPad           = TMath::Min(fLoPad, seq->GetPadId());
      fHiPad           = TMath::Max(fHiPad, seq->GetPadId() + 1);
      fLoBin           = TMath::Min(fLoBin, seq->GetBegin());
      fHiBin           = TMath::Max(fHiBin, seq->GetEnd());
   }

   return kTRUE;
}

void TpcCluster::BeginCluster(TpcPadRow *aPadRow)
{
   fSequences->Delete();
}

void TpcCluster::AddSequence(TpcSequence *aSequence)
{
   if (!aSequence) {
      Error(__FUNCTION__, "Invalid argument. aSequence == NULL");
      return;
   }
   fSequences->Add(aSequence);
}

void TpcCluster::EndCluster()
{
   if (!RecalcBounds()) return;

   MakeRegion();
}

void TpcCluster::MakeRegion()
{
   fRegion->Clear();
   fRegion->Expand((fHiPad - fLoPad) * (fHiBin - fLoBin));

   for (Int_t i = 0; i < fRegion->GetSize(); ++i) {
      fRegion->AddAt(&fBinStub, i);
   }

   for (Int_t i = 0; i < fSequences->GetEntriesFast(); ++i) {
      TpcSequence *seq = (TpcSequence *)fSequences->At(i);
      TpcPad      *pad = seq->GetPad();
      for (Int_t ibin = seq->GetBegin(); ibin < seq->GetEnd(); ++ibin) {
         TpcTimeBin *bin = pad->GetTimeBin(ibin);
         fRegion->AddAt(bin, (pad->GetPadId() - fLoPad) * (fHiBin - fLoBin) + ibin - fLoBin);
      }
   }
}

const TpcTimeBin *TpcCluster::GetTimeBin(Int_t ipad, Int_t ibin)
{
   if (ipad < fLoPad || ipad >= fHiPad) {
      Error(__FUNCTION__, "Invalid argument. Pad number exceed range");
      return NULL;
   }
   if (ibin < fLoBin || ibin >= fHiBin) {
      Error(__FUNCTION__, "Invalid argument. Bin number exceed range");
      return NULL;
   }

   Int_t i = (ipad - fLoPad) * (fHiBin - fLoBin) + (ibin - fLoBin);
   return (TpcTimeBin *)fRegion->At(i);
}
