/**
 * \class MpdMiniEvent
 * \brief Stores global information about the event
 *
 * The MpdMiniEvent class keeps variables that characterize event.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * email: nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniEvent_h
#define MpdMiniEvent_h

// C++ headers
#include <vector>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//_________________
class MpdMiniEvent : public TObject {

public:
   /// Default constructor
   MpdMiniEvent();
   /// Copy constructor
   MpdMiniEvent(const MpdMiniEvent &event);
   /// Destructor
   virtual ~MpdMiniEvent();
   /// Print some event information
   virtual void Print(const Char_t *option = "") const;

   //
   // Getters
   //

   /// Return run ID
   Int_t runId() const { return fRunId; }
   /// Return event ID
   Int_t eventId() const { return fEventId; }
   /// Return fill ID
   Int_t fillId() const { return (Int_t)fFillId; }
   /// Return magnetic field (kG)
   Float_t bField() const { return fBField; }
   /// Return time stamp
   Int_t time() const { return fTime; }

   /// Return primary vertex position
   TVector3 primaryVertex() const { return TVector3(fPrimaryVertexX, fPrimaryVertexY, fPrimaryVertexZ); }
   /// Return primary vertex position error
   TVector3 primaryVertexError() const
   {
      return TVector3(fPrimaryVertexErrorX, fPrimaryVertexErrorY, fPrimaryVertexErrorZ);
   }
   /// Return number of tracks that matched ECal
   UShort_t nBECalMatch() const { return fNBECalMatch; }
   /// Return number of tracks that matched TOF
   UShort_t nBTOFMatch() const { return fNBTOFMatch; }

   /// Return trigger list of the current event
   std::vector<unsigned int> triggerIds() const { return fTriggerIds; }
   /// Check if the trigger is in the list of triggers
   /// that were fired in the current event
   bool isTrigger(unsigned int) const;

   /// Return RefMult estimated via positive tracks (-0.5<eta<0.5)
   Int_t refMultPos() const { return (Int_t)fRefMultPos; }
   /// Return RefMult estimated via negative tracks (-0.5<eta<0.5)
   Int_t refMultNeg() const { return (Int_t)fRefMultNeg; }
   /// Return RefMult (-0.5<eta<0.5)
   Int_t refMult() const { return (Int_t)(fRefMultPos + fRefMultNeg); }
   /// Return RefMult estimated via positive tracks (-1<eta<-0.5)
   Int_t refMult2PosEast() const { return (Int_t)fRefMult2PosEast; }
   /// Return RefMult estimated via negative tracks (-1<eta<-0.5)
   Int_t refMult2NegEast() const { return (Int_t)fRefMult2NegEast; }
   /// Return RefMult estimated via positive tracks (0.5<eta<1)
   Int_t refMult2PosWest() const { return (Int_t)fRefMult2PosWest; }
   /// Return RefMult estimated via negative tracks (0.5<eta<1)
   Int_t refMult2NegWest() const { return (Int_t)fRefMult2NegWest; }
   /// Return RefMult estimated via positive particles in (eta<0)
   Int_t refMultHalfPosEast() const { return (Int_t)fRefMultHalfPosEast; }
   /// Return RefMult estimated via negative particles in (eta<0)
   Int_t refMultHalfNegEast() const { return (Int_t)fRefMultHalfNegEast; }
   /// Return RefMult estimated via positive particles in (eta>0)
   Int_t refMultHalfPosWest() const { return (Int_t)fRefMultHalfPosWest; }
   /// Return RefMult estimated via positive negative in (eta>0)
   Int_t refMultHalfNegWest() const { return (Int_t)fRefMultHalfNegWest; }
   /// Return RefMult2 value estimated by east TPC (-1<eta<0.5)
   Int_t refMult2East() const { return (Int_t)(fRefMult2PosEast + fRefMult2NegEast); }
   /// Return RefMult2 value estimated by west TPC (0.5<eta<1)
   Int_t refMult2West() const { return (Int_t)(fRefMult2PosWest + fRefMult2NegWest); }
   /// Return RefMult measured in 0.5<|RefMult|<1
   Int_t refMult2() const { return (Int_t)(fRefMult2PosEast + fRefMult2NegEast + fRefMult2PosWest + fRefMult2NegWest); }
   /// Return RefMult measured in the east TPC
   Int_t refMultHalfEast() const { return (Int_t)(fRefMultHalfPosEast + fRefMultHalfNegEast); }
   /// Return RefMult measured in the west TPC
   Int_t refMultHalfWest() const { return (Int_t)(fRefMultHalfPosWest + fRefMultHalfNegWest); }

   /// Return gRefMult (RefMult by global tracks in |gRefMult|<0.5)
   Int_t grefMult() const { return (Int_t)(fGRefMult); }
   /// Return total number of global tracks that were reconstructed in the event
   UShort_t numberOfGlobalTracks() const { return fNumberOfGlobalTracks; }
   /// Return number of hits in TOF trays
   UShort_t btofTrayMultiplicity() const { return fBTofTrayMultiplicity; }

   /// Return number of hits in the east Fast Forward Detector (FFD)
   Int_t nFfdHitsEast() const { return (Int_t)fNFfdHitsEast; }
   /// Return number of hits in the west Fast Forward Detector (FFD)
   Int_t nFfdHitsWest() const { return (Int_t)fNFfdHitsWest; }
   /// Return number of TOF tracks used for T0 calibration
   Int_t nTofT0() const { return (Int_t)fNTofT0; }
   /// Return z position of the primary vertex estimated by Ffd
   Float_t vzFfd() const { return fVzFfd; }

   /// Return ZDC coincidence rate
   Float_t fhcalX() const { return (Float_t)fFHCalX; }
   /// Return background rate
   Float_t backgroundRate() const { return fBackgroundRate; }
   /// Return east (eta<0) FHCal rate
   Float_t fhcalEastRate() const { return fFHCalEastRate; }
   /// Return west (eta>0) FHCal rate
   Float_t fhcalWestRate() const { return fFHCalWestRate; }

   /// Energy deposit in the FHCal module (90 modules in total)
   Float_t fhcalEnergyDeposit(const Int_t iModule) const
   {
      return (iModule < 0 || iModule > 89) ? -999. : fFHCalEnergyDeposit[iModule];
   }
   /// Energy deposit in the east (eta<0) FHCal module (45 modules)
   Float_t fhcalEastEnergyDeposit(const Int_t iModule) const
   {
      return (iModule < 0 || iModule > 44) ? -999. : fFHCalEnergyDeposit[iModule];
   }
   /// Energy deposit in the west (eta>0) FHCal module (45 modules)
   Float_t fhcalWestEnergyDeposit(const Int_t iModule) const
   {
      return (iModule < 0 || iModule > 44) ? -999. : fFHCalEnergyDeposit[45 + iModule];
   }
   /// Azimuthal angle of the BHCal module
   Double_t fhcalModuleAziAngle(const Int_t iModule) const;
   /// Azimuthal angle of the east (eta<0) BHCal module (45 modules)
   Double_t fhcalEastModuleAziAngle(const Int_t iModule) const;
   /// Azimuthal angle of the west (eta>0) BHCal module (45 modules)
   Double_t fhcalWestModuleAziAngle(const Int_t iModule) const;

   /// Return year
   Int_t year() const;
   /// Return day number
   Int_t day() const;
   /// Return bunch crossing number
   Int_t bunchId() const { return (Int_t)fBunchCrossId; }

   //
   // Setters
   //

   /// Set run ID
   void setRunId(Int_t id) { fRunId = id; }
   /// Set event ID
   void setEventId(Int_t id) { fEventId = id; }
   /// Set fill ID
   void setFillId(Int_t id) { fFillId = (UShort_t)id; }
   /// Set fill ID
   void setFillId(Float_t id) { fFillId = (id > 0) ? (UShort_t)id : 0; }
   /// Set magnetic field
   void setBField(Double_t bField) { fBField = (Float_t)bField; }
   /// Set magnetic field
   void setMagneticField(Double_t bField) { fBField = (Float_t)bField; }
   /// Set time stamp
   void setTime(Int_t time) { fTime = time; }

   /// Set primary vertex position (x,y,z)
   void setPrimaryVertexPosition(Float_t x, Float_t y, Float_t z)
   {
      fPrimaryVertexX = x;
      fPrimaryVertexY = y;
      fPrimaryVertexZ = z;
   }
   /// Set primary vertex position (3-vector)
   void setPrimaryVertexPosition(TVector3 vtxPos)
   {
      fPrimaryVertexX = vtxPos.X();
      fPrimaryVertexY = vtxPos.Y();
      fPrimaryVertexZ = vtxPos.Z();
   }
   /// Set primary vertex position error (ex,ey,ez)
   void setPrimaryVertexPositionError(Float_t x, Float_t y, Float_t z)
   {
      fPrimaryVertexErrorX = x;
      fPrimaryVertexErrorY = y;
      fPrimaryVertexErrorZ = z;
   }
   /// Set primary vertex position error (3-vector)
   void setPrimaryVertexPositionError(TVector3 vtxPosErr)
   {
      fPrimaryVertexErrorX = vtxPosErr.X();
      fPrimaryVertexErrorY = vtxPosErr.Y();
      fPrimaryVertexErrorZ = vtxPosErr.Z();
   }
   /// Set number of ECal-matched tracks
   void setNumberOfBECalMatch(Int_t n) { fNBECalMatch = (UShort_t)n; }
   //// Set number of TOF-matched tracks
   void setNumberOfBTOFMatch(Int_t n) { fNBTOFMatch = (UShort_t)n; }

   /// Set trigger id
   void setTriggerId(UInt_t id);
   /// Set trigger id (pass STL vector with trigger IDs)
   void setTriggerIds(std::vector<unsigned int> ids);

   /// Set RefMult negative (|eta|<0.5)
   void setRefMultNeg(UShort_t mult) { fRefMultNeg = (UShort_t)mult; }
   /// Set RefMult positive (|eta|<0.5)
   void setRefMultPos(UShort_t mult) { fRefMultPos = (UShort_t)mult; }
   /// Set negative RefMult2 east ( -1<eta<-0.5 )
   void setRefMult2NegEast(UShort_t mult) { fRefMult2NegEast = (UShort_t)mult; }
   /// Set positive RefMult2 east ( -1<eta<-0.5 )
   void setRefMult2PosEast(UShort_t mult) { fRefMult2PosEast = (UShort_t)mult; }
   /// Set negative RefMult2 west ( 0.5<eta<1 )
   void setRefMult2NegWest(UShort_t mult) { fRefMult2NegWest = (UShort_t)mult; }
   /// Set positive RefMult2 west ( 0.5<eta<1 )
   void setRefMult2PosWest(UShort_t mult) { fRefMult2PosWest = (UShort_t)mult; }
   /// TPC refMultHalf neg (eta<0)
   void setRefMultHalfNegEast(UShort_t mult) { fRefMultHalfNegEast = (UShort_t)mult; }
   /// TPC refMultHalf pos (eta<0)
   void setRefMultHalfPosEast(UShort_t mult) { fRefMultHalfPosEast = (UShort_t)mult; }
   /// TPC refMultHalf neg (eta>0)
   void setRefMultHalfNegWest(UShort_t mult) { fRefMultHalfNegWest = (UShort_t)mult; }
   /// TPC refMultHalf pos (eta>0)
   void setRefMultHalfPosWest(UShort_t mult) { fRefMultHalfPosWest = (UShort_t)mult; }

   /// Set RefMult estimated by global tracks
   void setGRefMult(UShort_t mult) { fGRefMult = (UShort_t)mult; }
   /// Set number of global tracks reconstructed in the event
   void setNumberOfGlobalTracks(UShort_t mult) { fNumberOfGlobalTracks = (UShort_t)mult; }
   /// Set total number of hits in TOF trays
   void setbTofTrayMultiplicity(UShort_t mult) { fBTofTrayMultiplicity = (UShort_t)mult; }

   /// Set number of hits in the east FFD
   void setNFfdHitsEast(UShort_t nHits) { fNFfdHitsEast = (UChar_t)nHits; }
   /// Set number of hits in the west FFD
   void setNFfdHitsWest(UShort_t nHits) { fNFfdHitsWest = (UChar_t)nHits; };
   /// Set number of T0 particles in BTOF self calibration
   void setNTofT0(Int_t t0) { fNTofT0 = (UShort_t)t0; }
   /// Set Vz of the primary vertex reconstructed by FFD
   void setVzFfd(Float_t FfdVz) { fVzFfd = FfdVz; }

   /// Set FHCal coincidence rate
   void setFHCalx(Float_t fhcalCoinRate) { fFHCalX = (UInt_t)fhcalCoinRate; }
   /// Set background rate
   void setBackgroundRate(Float_t bckgRate) { fBackgroundRate = (Float_t)bckgRate; }
   /// Set east FHCal rate
   void setFHCalEastRate(Float_t fhcalEastRate) { fFHCalEastRate = (Float_t)fhcalEastRate; }
   /// Set west FHCal rate
   void setFHCalWestRate(Float_t fhcalWestRate) { fFHCalWestRate = (Float_t)fhcalWestRate; }

   /// Set FHCal energy deposit measured in i-th module (0-89)
   void setFHCalEnergyDepositInModule(Int_t iModule, Float_t energyDeposit);
   /// Set east FHCal energy deposit measured in i-th module (eta<0, 0-44)
   void setFHCalEastEnergyDepositInModule(Int_t iModule, Float_t energyDeposit);
   /// Set west FHCal energy deposit measured in i-th module (eta<0, 45-89)
   void setFHCalWestEnergyDepositInModule(Int_t iModule, Float_t energyDeposit);

   /// Set bunch crossing ID
   void setBunchId(Int_t id);

private:
   /// Run number (or runId)
   Int_t fRunId;
   /// Event ID
   Int_t fEventId;
   /// Fill number
   UShort_t fFillId;
   /// Number of bunch crossing
   UChar_t fBunchCrossId;
   /// Magnetic field strength (kG)
   Float_t fBField;
   /// GMT time casted into unsigned int
   Int_t fTime;

   /// Primary vertex position X
   Float_t fPrimaryVertexX;
   /// Primary vertex position Y
   Float_t fPrimaryVertexY;
   /// Primary vertex position Z
   Float_t fPrimaryVertexZ;
   /// Primary vertex position error X
   Float_t fPrimaryVertexErrorX;
   /// Primary vertex position error Y
   Float_t fPrimaryVertexErrorY;
   /// Primary vertex position error Z
   Float_t fPrimaryVertexErrorZ;

   /// Number of ECal-matched tracks
   UShort_t fNBECalMatch;
   /// Number of TOF-matched tracks
   UShort_t fNBTOFMatch;

   /// List of triggers that were fired in the current event
   std::vector<unsigned int> fTriggerIds;

   ///
   /// Reference multiplicity (RefMult) - number of tracks with existing momentum (p>0 GeV/c),
   /// nHits>10 and distance of closest approach (DCA) between track and primary
   /// vertex DCA<3 cm.
   ///
   /// Several RefMult quantities can be used for centrality estimation. The importance of
   /// them is to see the effect of the centrality definition on the fluctuation measurements.
   ///

   /// Reference multiplicity estimated via negative tracks (-0.5<eta<0.5)
   UShort_t fRefMultNeg;
   /// Reference multiplicity estimated via positive tracks (-0.5<eta<0.5)
   UShort_t fRefMultPos;
   /// Reference multiplicity of negative tracks in TPC at large eta ( refMult2 -1<eta<-0.5 )
   UShort_t fRefMult2NegEast;
   /// Reference multiplicity of positive tracks in TPC at large eta ( refMult2 -1<eta<-0.5 )
   UShort_t fRefMult2PosEast;
   /// Reference multiplicity of negative tracks in TPC at large eta ( refMult2 0.5<eta<1 )
   UShort_t fRefMult2NegWest;
   /// Reference multiplicity of positive tracks in TPC at large eta ( refMult2 0.5<eta<1 )
   UShort_t fRefMult2PosWest;
   /// TPC refMultHalf for negative particles (eta<0)
   UShort_t fRefMultHalfNegEast;
   /// TPC refMultHalf for positive particles (eta<0)
   UShort_t fRefMultHalfPosEast;
   // TPC refMultHalf for negative particles (eta>0)
   UShort_t fRefMultHalfNegWest;
   /// TPC refMultHalf for positive particles (eta>0)
   UShort_t fRefMultHalfPosWest;

   /// Reference multiplicity estimated by global tracks (DCA cut is not applied)
   UShort_t fGRefMult;
   /// Total number of global tracks reconstructed in the event
   UShort_t fNumberOfGlobalTracks;

   /// Total hit multiplicity in TOF trays
   UShort_t fBTofTrayMultiplicity;

   /// Number of hits in east Fast Forward Detector (FFD)
   UChar_t fNFfdHitsEast;
   /// Number of hits in west Fast Forward Detector (FFD)
   UChar_t fNFfdHitsWest;
   /// Vz estimated via Fast Forward Detector (FFD)
   Float_t fVzFfd;
   /// Number of T0 particles in BTOF self calibration
   UShort_t fNTofT0;

   /// FHCal coincidence rate
   UInt_t fFHCalX;
   /// Background rate
   Float_t fBackgroundRate;
   /// East FHCal rate
   Float_t fFHCalEastRate;
   /// West FHCal rate
   Float_t fFHCalWestRate;

   /// Energy deposit in each of 90 towers 0-44 for eta<0 and 45-89 for eta>0
   Float16_t fFHCalEnergyDeposit[90];

   ClassDef(MpdMiniEvent, 2);
};

#endif // #define MpdMiniEvent_h
