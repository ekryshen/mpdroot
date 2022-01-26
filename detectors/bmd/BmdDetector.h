//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Detector class - sensitive detector for VMC
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//     Pedro Gonzalez Zamora (original author) 
//
//
//-----------------------------------------------------------

#ifndef BMDDETECTOR_HH
#define BMDDETECTOR_HH

// Base Class Headers ----------------
#include "FairDetector.h"

// Collaborating Class Headers -------
#include <TLorentzVector.h>
#include <ostream> // remove if you do not need streaming op
#include "TVector3.h"
#include "TGeoMedium.h"
#include <map>

using namespace std;

// Collaborating Class Declarations --
class TClonesArray;
class FairVolume;
class BmdPoint;


class BmdDetector : public FairDetector {
public:
  // Constructors/Destructors ---------
  /**      Name :  Detector Name
	   Active: kTRUE for active detectors
	   (ProcessHits() will be called)
	   kFALSE for inactive detectors     */
  BmdDetector(const char * Name, Bool_t Active);

  /**      default constructor    */
  BmdDetector();
  /**       destructor     */
  virtual ~BmdDetector();

  /**      Initialization of the detector is done here    */
  //virtual void   Initialize();

  /**       this method is called for each step
            during simulation (see FairMCApplication::Stepping())    */
  virtual Bool_t ProcessHits( FairVolume *v=0);

  /**       Registers the produced collections in CBMRootManager.     */
  virtual void   Register();

  virtual void Print() const;

  /** Gets the produced collections */
  virtual TClonesArray* GetCollection(Int_t iColl) const ;

  /**      has to be called after each event
	   to reset the containers        */
  virtual void   Reset();

  /**      Create the detector geometry        */
  virtual void ConstructGeometry();

  // Construct the geometry from an ASCII geometry file
  virtual void ConstructAsciiGeometry();

  // Construct the geometry from a GDML geometry file
  virtual void ConstructGDMLGeometry();
  void ExpandNodeForGdml(TGeoNode* node);
  map<TString, TGeoMedium*> fFixedMedia; // List of media "repaired" after importing GMDL

  // Check whether a volume is sensitive.
  // The decision is based on the volume name. Only used in case
  // of GDML and ROOT geometry.
  // @param name    Volume name
  // @value         kTRUE if volume is sensitive, else kFALSE
  virtual Bool_t CheckIfSensitive(std::string name);

  void EndOfEvent();

  BmdPoint* AddHit(Int_t trackID, Int_t detID, TVector3 pos,
 		   TVector3 mom, Double_t time, Double_t length,
 		   Double_t eLoss);
  
  Int_t GetVolumeID(TString volname);
  // Operators


  // Accessors -----------------------


  // Modifiers -----------------------


  // Operations ----------------------

private:

  // Private Data Members ------------
  TClonesArray  	*fBmdPointCollection;
  Int_t			fTrackID;           //!  track index
  Int_t			fVolumeID;          //!  volume id
  TLorentzVector	fPos;               //!  position
  TLorentzVector	fMom;               //!  momentum
  Double_t		fTime;              //!  time
  Double_t		fLength;            //!  length
  Double_t		fELoss;             //!  energy loss
  Int_t                 currentTrackID;     //!
  Int_t                 currentEvent;       //! current event
  Bool_t                fNewTrack;
  const double	nan;			    //!
  
  void 		ResetParameters();

  // Private Methods -----------------

public:

  ClassDef(BmdDetector,1)

};

//------------------------------------------------------------------------------------------------------------------------
inline void BmdDetector::ResetParameters() 
{
	fTrackID = fVolumeID = 0;
	fPos.SetXYZM(nan, nan, nan, nan);
	fMom.SetXYZM(nan, nan, nan, nan);
	fTime = fLength = fELoss = nan;
	//fPosIndex = 0;
};
//-------------------------


#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
