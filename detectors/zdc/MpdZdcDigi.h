/*************************************************************************************
 *
 *            MpdZdcDigi
 *    Class for digital data taken from MpdZdc detector
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  18-Apr-2008
 *  Modified March 2021  by A.Strijak
 *  Modified June 2022  by S.Morozov
 *
 ************************************************************************************/

#ifndef MPDZDCDIGI_H
#define MPDZDCDIGI_H

#include "MpdZdcPoint.h"
#include <map>

class MpdZdcDigi : public TObject {
public:
   MpdZdcDigi();
   MpdZdcDigi(Int_t pfDetectorID, Int_t pfModuleID, Int_t pfChannelID, Double_t pfELoss = 0, Bool_t pIsPsd = kFALSE);
   MpdZdcDigi(MpdZdcPoint *p, Bool_t pIsPsd = kFALSE);

   virtual ~MpdZdcDigi();
   void Clear();

   void         InitStatic(Bool_t pIsPsd = kFALSE);
   virtual void Print(const Option_t *opt = "");

   inline Int_t GetModuleID() { return fModuleID; }

   inline Int_t GetModuleX() { return fModuleX; }
   inline Int_t GetModuleY() { return fModuleY; }

   inline Int_t    GetDetectorID() { return fDetectorID; }
   inline Int_t    GetChannelID() { return fChannelID; }
   inline Int_t    GetADCBits() { return fADCBits; }
   inline UInt_t   GetELossDigi() { return fELossDigi; }
   inline Double_t GetELoss() { return fELoss; }
   inline Double_t GetADCResolution() { return fADCResolution; }
   inline Double_t GetEnergyDigiThreshold() { return fEnergyDigiThreshold; }
   inline Double_t GetELossReco() { return fELossReco; }

   inline char   GetWasInitialized() { return fWasInitialized; }
   inline Bool_t GetIsPsd() { return fIsPsd; }

   inline Int_t SetModuleID(UInt_t pfModuleID)
   {
      fModuleID = pfModuleID;
      return fModuleID;
   }

   inline Double_t SetModuleX(Double_t pfModuleX)
   {
      fModuleX = pfModuleX;
      return fModuleX;
   }

   inline Double_t SetModuleY(Double_t pfModuleY)
   {
      fModuleY = pfModuleY;
      return fModuleY;
   }

   inline Int_t SetDetectorID(UInt_t pfDetectorID)
   {
      fDetectorID = pfDetectorID;
      return fDetectorID;
   }
   inline Int_t SetChannelID(UInt_t pfChannelID)
   {
      fChannelID = pfChannelID;
      return fChannelID;
   }
   inline Int_t SetADCBits(Int_t pfADCBits)
   {
      fADCBits = pfADCBits;
      return fADCBits;
   }
   inline Double_t SetADCResolution(Int_t pfADCResolution)
   {
      fADCResolution = pfADCResolution;
      return fADCResolution;
   }
   inline Double_t SetEnergyDigiThreshold(Int_t pfEnergyDigiThr)
   {
      fEnergyDigiThreshold = pfEnergyDigiThr;
      return fEnergyDigiThreshold;
   }

   inline Double_t AddEloss(Double_t pfELoss)
   {
      fELoss += pfELoss;
      return fELoss;
   }
   UInt_t AddZdcPoint(MpdZdcPoint *p);

   UInt_t Convert()
   {
      fELossDigi = ADC(fELoss);
      fELoss     = 0;
      return fELossDigi;
   }
   UInt_t ConvertSim()
   {
      fELossDigi = ADC(fELoss);
      return fELossDigi;
   }
   UInt_t   ADC(Double_t pfELoss);
   Double_t SetELossReco(Double_t pfELoss)
   {
      fELossReco = pfELoss;
      return fELossReco;
   }

   void IncreaseTimesCorrVR(Float_t e, Float_t time);

   std::map<Float_t, Float_t> GetContribTimesCorr() { return fContribTimesCorr; }

protected:
   static char
      fWasInitialized;       // 0 - from built-in constants; 1 - was attempt to read; 2 - success from MpdZdcDigiScheme
   static Int_t    fADCBits; // Hardware parameter
   static Double_t fADCResolution; // Hardware parameter
   static Double_t fEnergyDigiThreshold;

   Int_t    fDetectorID; // Detector Id
   Int_t    fModuleID;   // Hardware channel number (Module number)
   Int_t    fChannelID;  // Hardware channel number (Sub-module number)
   UInt_t   fELossDigi;  // Sum of the energy losses as digital response of the detector/channel
   Double_t fELoss;      // Sum of the energy losses as analog signal accumulated prior to digitalization
   Double_t fELossReco;  // Sum of the energy losses as analog simulated response of the detector/channel

   Double_t fModuleX; // module X coordinates
   Double_t fModuleY; // module X coordinates

   Bool_t fIsPsd; // static initialization tried from from MpdZdcDigiPar (0) or MpdZdcPsdDigiPar (1)

   std::map<Float_t, Float_t> fContribTimesCorr; // corrected time vs. deposited energy

   ClassDef(MpdZdcDigi, 3);
};

#endif // MPDZDCDIGI_H
