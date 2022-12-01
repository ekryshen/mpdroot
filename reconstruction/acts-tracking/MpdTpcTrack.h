// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#ifndef MPDTPCTRACK_H
#define MPDTPCTRACK_H

#include <TMath.h>
#include <TObject.h>
#include <TObjArray.h>

class MpdTpcTrack : public TObject {
public:
   explicit MpdTpcTrack(Int_t size);
   virtual ~MpdTpcTrack();

   TObjArray *GetHits() { return fHits; }

   Bool_t IsSortable() const { return kFALSE; }
   Int_t Compare(const TObject *track) const { return 0; }
   void Print(Option_t *opt) {}

private:
   TObjArray *fHits;
   ClassDef(MpdTpcTrack, 0);
};

#endif // MPDTPCTRACK_H
