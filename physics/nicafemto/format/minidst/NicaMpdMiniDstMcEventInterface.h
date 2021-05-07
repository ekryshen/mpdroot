/*
 * NicaMpdMiniDstMcEventInterface.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDMINIDSTMCEVENTINTERFACE_H_
#define NICAMPDMINIDSTMCEVENTINTERFACE_H_


#include "MpdMiniMcEvent.h"
#include "NicaMCEventInterface.h"
#include "NicaTrackClones.h"

class NicaTrackClones;
class NicaMpdMiniDstMcEvent;

class NicaMpdMiniDstMcEventInterface : public NicaMCEventInterface {
friend class NicaMpdMiniDstMcEvent;
    NicaTrackClones *fEvent;
    NicaTrackClones *fTracks;
protected:
    virtual void ConnectToTree();
    MpdMiniMcEvent *GetEvent()const{return (MpdMiniMcEvent*)fEvent->UncheckedAt(0);}
public:
    NicaMpdMiniDstMcEventInterface();
    virtual Int_t GetTotalTrackNo() const {return fTracks->GetEntriesFast();};
    virtual void Register(Bool_t write);
    virtual void Clear(Option_t *opt="");
    virtual void Compress(Int_t *map, Int_t map_size);
    virtual void CopyData(NicaEventInterface *s);
    virtual void CopyAndCompress(NicaEventInterface *s, Int_t *map, Int_t map_size);
    virtual void FillTrackInterface(NicaTrackInterface *track, Int_t index);
    virtual TObject *GetRawEventPointer()const{return fEvent;};
    virtual TObject * GetRawTrackPointer(Int_t index)const{return fTracks->UncheckedAt(index);}
    virtual NicaTrackInterface *GetTrackInterface()const;
    /** GETTERS SETTERS */
    virtual void SetVertex(Double_t x, Double_t y, Double_t z){GetEvent()->setPrimaryVertex(x, y, z);}
    virtual void SetPhi(Double_t phi, Double_t phi_error){GetEvent()->setReactionPlaneAngle(phi);};
    virtual void SetB(Double_t b){GetEvent()->setImpactParameter(b);};
    virtual Double_t GetB()const{return GetEvent()->b();};
    virtual Double_t GetPhi()const{return GetEvent()->phi();};
    virtual Double_t GetPhiError()const{return 0;};
    virtual TLorentzVector GetVertex()const;
    virtual ~NicaMpdMiniDstMcEventInterface();
    ClassDef(NicaMpdMiniDstMcEventInterface,1)
};

#endif /* NICAMPDMINIDSTMCEVENTINTERFACE_H_ */
