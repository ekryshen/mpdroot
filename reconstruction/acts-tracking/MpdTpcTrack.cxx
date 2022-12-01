// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTrack.h"

MpdTpcTrack::MpdTpcTrack(Int_t size):
  fHits(new TObjArray(size)) {
}

MpdTpcTrack::~MpdTpcTrack() {
  delete fHits;
}

ClassImp(MpdTpcTrack);
