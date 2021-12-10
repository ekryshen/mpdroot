/*
 * MpdFloor.h
 *
 *  Created on: 20 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef Floor_Floor_MPDFloor_H_
#define Floor_Floor_MPDFloor_H_
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "FairDetector.h"

class MpdFloorHit;
class FairVolume;
class MpdFloorGeo;
class MpdFloorPoint;

class MpdFloor : public FairDetector{
  	Int_t		fTrackID;           //!  track index
  	Int_t		fVolumeID;          //!  volume id
  	TLorentzVector	fPos;               //!  position
  	TLorentzVector	fMom;               //!  momentum
  	Double_t	fTime;              //!  time
  	Double_t	fLength;            //!  length
  	Double_t	fELoss;             //!  energy loss

public:
	MpdFloor();
	MpdFloor(const char* name, Bool_t active);
	virtual ~MpdFloor();
	// Defines the action to be taken when a step is inside the
	// active volume. Creates MpdStsPoints and adds them to the collection.
	// @param vol  Pointer to the active volume
    virtual Bool_t  ProcessHits(FairVolume* vol = 0);

   	// If verbosity level is set, print hit collection at the
   	// end of the event and resets it afterwards.
  	virtual void EndOfEvent();

   	// Registers the hit collection in the ROOT manager.
  	virtual void Register();

  	// Accessor to the hit collection
  	virtual TClonesArray* GetCollection(Int_t iColl) const;

   	// Screen output of hit collection.
	virtual void Print() const;

   	// Clears the hit collection
	virtual void Reset();

	// *@param cl1     Origin
	// *@param cl2     Target
	// *@param offset  Index offset
 	virtual void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset);

	// Constructs the STS geometry
 	virtual void ConstructGeometry();


 	virtual Bool_t CheckIfSensitive(std::string name);
 	MpdFloorPoint* AddHit(Int_t trackID, Int_t detId, TVector3 pos, TVector3 mom, Double_t time,
 			Double_t lenght, Double_t eloss);

 	ClassDef(MpdFloor,1)
private:
 	TClonesArray *fPointCollection;
 	static MpdFloorGeo *fgGeo;
};

#endif /* Floor_Floor_MPDFloor_H_ */

