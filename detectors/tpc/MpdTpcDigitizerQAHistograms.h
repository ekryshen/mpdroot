//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      QA for TpcВшпшешяукTask
//
//
// Environment:
//      Software developed for the MPD at NICA
//
// Author List:
//      Sergey Merts
//
//
//-----------------------------------------------------------

#ifndef MPDTPCDIGITIZERQAHISTOGRAMS_H
#define MPDTPCDIGITIZERQAHISTOGRAMS_H

// Base Class Headers ----------------
#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include "TaskHelpers.h"
#include "MpdTpcSector.h"

class MpdTpcDigitizerQAHistograms : public TObject {
private:
   std::string _suffix;

public:
   MpdTpcDigitizerQAHistograms();
   MpdTpcDigitizerQAHistograms(const std::string &suffix);
   ~MpdTpcDigitizerQAHistograms();

   void Initialize();
   void Write();

   TH2F *_hRZ_global;
   TH2F *_hYZ_local;
   TH2F *_hDiffuseXY;
   TH2F *_hDistortXY;
   TH1F *_hX_global;
   TH1F *_hY_global;
   TH1F *_hZ_global;
   TH1F *_hX_local;
   TH1F *_hY_local;
   TH1F *_hZ_local;
   TH2F *_hXY_local;
   TH2F *_hXY_global;
   TH1F *_hSect_dig;
   TH1F *_hADC_dig;
   TH1F *_hX_dig;
   TH1F *_hY_dig;
   TH1F *_hZ_dig;
   TH2F *_hXY_dig;
   TH3F *_h3D_dig;
   TH3F *_h3D_el;
   TH2F *_hXT_dig_1;
   TH2F *_hXT_dig_5;
   TH2F *_hXT_dig_10;
   TH2F *_hXT_dig_20;
   TH2F *_hXT_dig_40;
   TH2F *_hXT_dig_60;

public:
   ClassDef(MpdTpcDigitizerQAHistograms, 1)
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------