
/*************************************************************************************
 *
 *            MpdZdcDigi
 *    Class for digital data taken from MpdZdc detector
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  18-Apr-2008
 *  Modified March 2021  by A.Strijak
 *
 ************************************************************************************/

#include "MpdZdcDigi.h"
#include "MpdZdcDigiScheme.h"

#include <iostream>
using std::cout;
using std::endl;

static const Int_t    kADCBits             = 30;   // build-in constant
static const Double_t kADCResolution       = 1e-9; // build-in constant
static const Double_t kEnergyDigiThreshold = 0;    // build-in constant

char MpdZdcDigi::fWasInitialized = 0;

Int_t    MpdZdcDigi::fADCBits             = kADCBits;
Double_t MpdZdcDigi::fADCResolution       = kADCResolution;
Double_t MpdZdcDigi::fEnergyDigiThreshold = kEnergyDigiThreshold;

// -------------------------------------------------------------------------

MpdZdcDigi::MpdZdcDigi()
{
   fDetectorID = -1;
   fModuleID   = -1;
   fChannelID  = -1;
   fELoss      = 0;
   fELossDigi  = 0;
   fELossReco  = 0.;
   fIsPsd      = kFALSE;
}

// -------------------------------------------------------------------------

MpdZdcDigi::MpdZdcDigi(Int_t pfDetectorID, Int_t pfModuleID, Int_t pfChannelID, Double_t pfELoss, Bool_t pIsPsd)
{
   fDetectorID = pfDetectorID;
   fModuleID   = pfModuleID;
   fChannelID  = pfChannelID;
   fELossDigi  = 0;

   InitStatic(pIsPsd);
   fELoss = pfELoss;
}

// -------------------------------------------------------------------------

MpdZdcDigi::MpdZdcDigi(MpdZdcPoint *p, Bool_t pIsPsd)
{
   fDetectorID = -1;
   fModuleID   = -1;
   fChannelID  = -1;
   fELoss      = 0;
   fELossDigi  = 0;

   MpdZdcDigiScheme *pDigiScheme = MpdZdcDigiScheme::Instance();

   if ((pDigiScheme) && (p)) {

      Int_t detID, modID, chanID;
      pDigiScheme->SplitDigiID(pDigiScheme->GetDigiIdFromCoords(p->GetX(), p->GetY(), p->GetZ()), detID, modID, chanID);

      fDetectorID = detID;
      fModuleID   = modID;
      fChannelID  = chanID;

      fELoss = p->GetEnergyLoss();
   }

   InitStatic(pIsPsd);
}

// -------------------------------------------------------------------------

MpdZdcDigi::~MpdZdcDigi() {}

// -------------------------------------------------------------------------

void MpdZdcDigi::Clear()
{
   fWasInitialized      = 0;
   fIsPsd               = kFALSE;
   fADCBits             = kADCBits;
   fADCResolution       = kADCResolution;
   fEnergyDigiThreshold = kEnergyDigiThreshold;

   fELossDigi = 0;
   fELoss     = 0;
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

void MpdZdcDigi::InitStatic(Bool_t pIsPsd)
{
   Bool_t do_init = ((fIsPsd != pIsPsd) || (!fWasInitialized));
   if (do_init) {
      Clear();
      fIsPsd                        = pIsPsd;
      fWasInitialized               = 1;
      MpdZdcDigiScheme *pDigiScheme = MpdZdcDigiScheme::Instance();

      if (pDigiScheme) {

         MpdZdcDigiPar *pZdcDigiPar = 0;

         if (pIsPsd)
            pZdcDigiPar = (MpdZdcDigiPar *)pDigiScheme->GetZdcPsdDigiPar();
         else
            pZdcDigiPar = pDigiScheme->GetZdcDigiPar();

         if (pZdcDigiPar) {
            fADCBits             = pZdcDigiPar->GetADCBits();
            fADCResolution       = pZdcDigiPar->GetADCResolution();
            fEnergyDigiThreshold = pZdcDigiPar->GetEnergyDigiThreshold();
            fWasInitialized      = 2;
         }
      }
   }
}

// -------------------------------------------------------------------------

void MpdZdcDigi::Print(const Option_t *opt)
{
   cout << " MpdZdcDigi  DetID:" << fDetectorID << " ModuleID:" << fModuleID << " ChanID:" << fChannelID
        << " ELossDigi:" << fELossDigi << " ELoss: " << fELoss << "  IsPsd: " << fIsPsd << "  [Bits:" << fADCBits
        << ", Resol:" << fADCResolution << ", Thresh:" << fEnergyDigiThreshold << ", Init:" << fWasInitialized << "]"
        << endl;
}

// -------------------------------------------------------------------------

UInt_t MpdZdcDigi::AddZdcPoint(MpdZdcPoint *p)
{
   MpdZdcDigiScheme *pDigiScheme = MpdZdcDigiScheme::Instance();

   if ((pDigiScheme) && (p)) {

      Int_t detID, modID, chanID;
      pDigiScheme->SplitDigiID(pDigiScheme->GetDigiIdFromCoords(p->GetX(), p->GetY(), p->GetZ()), detID, modID, chanID);

      if ((fDetectorID == detID) && (fModuleID == modID) && (fChannelID == chanID)) fELoss += p->GetEnergyLoss();
   }

   return 0;
}

// -------------------------------------------------------------------------

UInt_t MpdZdcDigi::ADC(Double_t pfELoss)
{
   //  if (pfELoss<=fEnergyDigiThreshold) {
   if (pfELoss > 0) {
      UInt_t tmp = 0;
      tmp        = (UInt_t)(pfELoss / fADCResolution);
      return (tmp & ((2 << (fADCBits - 1)) - 1));
   } else
      return 0;
}

// -------------------------------------------------------------------------

ClassImp(MpdZdcDigi);
