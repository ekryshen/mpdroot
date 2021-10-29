/*
 * NicaMpdMiniDstMcTrack.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDMINIDSTMCTRACK_H_
#define NICAMPDMINIDSTMCTRACK_H_

#include "NicaMCTrack.h"
#include "NicaMCTrackInterface.h"

class NicaMpdMiniDstMcTrack : public NicaMCTrack{
public:
    NicaMpdMiniDstMcTrack();
    virtual ~NicaMpdMiniDstMcTrack();
    ClassDef(NicaMpdMiniDstMcTrack,1)
};

class NicaMpdMiniDstMcTrackInterface: public NicaMCTrackInterface{
public:
    NicaMpdMiniDstMcTrackInterface(){};
    virtual ~NicaMpdMiniDstMcTrackInterface(){};
    ClassDef(NicaMpdMiniDstMcTrackInterface,1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MINIDST_NICAMPDMINIDSTMCTRACK_H_ */
