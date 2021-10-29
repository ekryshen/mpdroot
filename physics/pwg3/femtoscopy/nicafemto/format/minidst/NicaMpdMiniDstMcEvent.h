/*
 * NicaMpdMiniDstMcEvent.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDMINIDSTMCEVENT_H_
#define NICAMPDMINIDSTMCEVENT_H_

#include "NicaMCEvent.h"

class NicaMpdMiniDstMcEvent : public NicaMCEvent{
public:
    NicaMpdMiniDstMcEvent();
    NicaMpdMiniDstMcEvent(const NicaMpdMiniDstMcEvent &other);
    virtual void CreateSource();
    virtual void Update();
    virtual void Clear(Option_t *opt);
    virtual void Print();
    virtual Bool_t ExistInTree()const;
    virtual TString GetFormatName()const;
    virtual ~NicaMpdMiniDstMcEvent();
    ClassDef(NicaMpdMiniDstMcEvent,1)
};

#endif /* NICAMPDMINIDSTMCEVENT_H_ */
