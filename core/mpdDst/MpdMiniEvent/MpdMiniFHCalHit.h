/**
 * \class MpdMiniFHCalHit
 * \brief Holds information about FHCal hits
 *
 * The class holds information about the hit in
 * the Forward Hadronic Calorimeter (FHCal)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \email nigmatkulov@gmail.com; ganigmatkulov@mephi.ru
 * \date June 19, 2020
 */

#ifndef MpdMiniFHCalHit_h
#define MpdMiniFHCalHit_h

// ROOT headers
#include "TObject.h"
#include "TMath.h"

//_________________
class MpdMiniFHCalHit : public TObject {

public:
   /// Default constructor
   MpdMiniFHCalHit();
   /// Copy constructor
   MpdMiniFHCalHit(const MpdMiniFHCalHit &hit);
   /// Destructor
   virtual ~MpdMiniFHCalHit();
   /// Print tower information
   virtual void Print(const Char_t *option = "") const;

   //
   // Getters
   //

   /// Return FHCal detector side (1 - East, 2 - West)
   /// \par iTower Tower (aka hit) index
   Int_t detector() const { return fId / 10000; }

   /// Return FHCal module ID (1-45) for each side
   Int_t module() const { return (fId / 100) % 100; }

   /// Channel ID for the given hit
   Int_t channel() const { return fId % 100; }

   /// Energy deposited in the tower
   Float_t eDep() const { return fEDep; }
   // ADC measured in the tower
   // Int_t adc() const              { return (UInt_t)fAdc; }

   /// If channel is bad/dead
   Bool_t isBad() const { return (fEDep < 0) ? true : false; }

   /// Check if unknown channel
   Bool_t isUnknown() const { return (fId == 0) ? true : false; }

   //
   // Setters
   //

   /// Set hit ID
   /// \par det DetectorID
   /// \par mod ModuleID
   /// \par chan ChannelID
   void setId(Int_t det, Int_t mod, Int_t chan) { fId = det * 10000 + mod * 100 + chan; }

   // Set ADC measured in the channel
   // void setAdc(UInt_t adc)   { fAdc = adc; }

   /// Set energy deposition
   void setEDep(Float_t eDep) { fEDep = eDep; }

private:
   /// Hit id: detectorId * 10000 + moduleId * 100 + channelId
   /// DetectorId: 1,2; ModuleId: 1-45; ChannelId: 1-42

   // Another option: ((detId-1)*45+modId)*100 + chId ->
   // ch = fId%100;
   // det = ( (fId/100) > 45) ? 2 : 1;
   // mod = fId/100 - (det-1)*45
   UShort_t fId;

   /// Energy deposition
   Float16_t fEDep;

   // ADC for the hit
   // Float16_t fAdc;

   ClassDef(MpdMiniFHCalHit, 2)
};

#endif // #define MpdMiniFHCalHit_h
