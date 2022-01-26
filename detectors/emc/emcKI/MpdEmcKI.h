#ifndef MPDEMCKIDET_H
#define MPDEMCKIDET_H

#include "FairDetector.h"

class MpdEmcPointKI;
class FairVolume;
class TClonesArray;
class TVector3;
class MpdEmcGeoUtils;

class MpdEmcKI : public FairDetector
{
 public:
  /**      Name :  Detector Name
   *       Active: kTRUE for active detectors (ProcessHits() will be called)
   *               kFALSE for inactive detectors
   */
  MpdEmcKI(const char* Name, Bool_t Active);

  /**      default constructor    */
  MpdEmcKI();

  /**       destructor     */
  virtual ~MpdEmcKI();

  /**      Initialization of the detector is done here    */
  //    virtual void   Initialize();

  /**       this method is called for each step during simulation
   *       (see FairMCApplication::Stepping())
   */
  virtual Bool_t ProcessHits(FairVolume* v = 0);

  /**       Registers the produced collections in FAIRRootManager.     */
  virtual void Register();

  /** Gets the produced collections */
  virtual TClonesArray* GetCollection(Int_t iColl) const;

  /**      has to be called after each event to reset the containers      */
  virtual void Reset();

  /**
    this is called at the end of an event
  */
  virtual void EndOfEvent();

  virtual void FinishEvent();

  /** Screen output of hit collection.     */
  virtual void Print() const;

  /**      Create the detector geometry        */
  void ConstructGeometry();

  /**      Create the detector geometry        */
  void ConstructAsciiGeometry();

  /**      Check sensitivity        */
  virtual Bool_t CheckIfSensitive(std::string name);

  /**      This method is an example of how to add your own point
   *       of type MpdEmcPoint to the clones array
   */
  MpdEmcPointKI* AddHit(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                        Double_t ELoss);

 private:
  std::map<int, int> fSuperParents;     //! map of current tracks to SuperParents: entered ECAL active volumes particles
  TClonesArray* fMpdEmcPointCollection; //! Collection of ECAL hits
  Int_t fNhits;                         //! Number of created so far hits
  Int_t fCurrentTrackID;                //! current track Id
  Int_t fCurrentCellID;                 //! current cell Id
  Int_t fCurentSuperParent;             //! current SuperParent ID: particle entered PHOS
  MpdEmcPointKI* fCurrentHit;             //! current Hit
  MpdEmcGeoUtils* fGeom;                //! pointer to GeoUtils
  Double_t fEcalRmin;                   //! inner and outer
  Double_t fEcalRmax;                   //! radii of ECAL (used for tracking)

  ClassDef(MpdEmcKI, 2)
};

#endif // MPDEMCDET_H
