/*
 * NicaMpdEventInterface.h
 *
 *  Created on: 3 sie 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_FORMAT_NICAMPDEVENTINTERFACE_H_
#define INTERFACES_MPDROOT_FORMAT_NICAMPDEVENTINTERFACE_H_
#include "NicaExpEventInterface.h"
#include "MpdEvent.h"
class NicaMpdDstKalmanEventInterface;
class NicaMpdEventInterface: public NicaExpEventInterface {
friend class NicaMpdDstKalmanEventInterface;
	enum kTrackType{
		kAllTracks,
		kPrimaryTracks,
		kGlobalTracks
	};
	kTrackType fMode;
	MpdEvent *fEvent;
protected:
	virtual void ConnectToTree();
	void Register(Bool_t write);
public:
	NicaMpdEventInterface();
	virtual void CopyData(NicaEventInterface *s);
	void OnlyPrimary(){fMode = kPrimaryTracks;};
	void OnlyGlobal(){fMode = kGlobalTracks;};
	virtual void Compress(Int_t *map, Int_t map_size);
	virtual void CopyAndCompress(NicaEventInterface *s, Int_t *map, Int_t map_size);
	virtual void FillTrackInterface(NicaTrackInterface *track, Int_t index);
	virtual Int_t GetTotalTrackNo() const;
	virtual TObject *GetRawEventPointer()const {return fEvent;};
	virtual NicaTrackInterface *GetTrackInterface()const;
	virtual TObject* GetRawTrackPointer(Int_t index)const;
	/** GETTERS AND SETTERS **/
	virtual void SetRunInfoId(Int_t i);
	virtual Int_t GetRunInfoId()const;;
	virtual void SetMagneticField(TVector3 mag)const;
	virtual TVector3 GetMagneticField()const;
	virtual TLorentzVector GetVertexError()const;
	virtual TLorentzVector GetVertex()const;
	virtual ~NicaMpdEventInterface();
	ClassDef(NicaMpdEventInterface,1)
};

#endif /* INTERFACES_MPDROOT_FORMAT_NICAMPDEVENTINTERFACE_H_ */
