/*
 * MpdPairTpcEntrySimpleCut.h
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICA_MPD_CUTS_PAIRCUTS_SIMPLE_TPC_MPDPAIRTPCENTRYSIMPLECUT_H_
#define NICA_MPD_CUTS_PAIRCUTS_SIMPLE_TPC_MPDPAIRTPCENTRYSIMPLECUT_H_

#include "MpdPairTpcSimpleCut.h"

class MpdPairTpcEntrySimpleCut : public MpdPairTpcSimpleCut{
	Int_t fRmin;
public:
	MpdPairTpcEntrySimpleCut();
	static Int_t XYZ(){return 0;}
	static Int_t XY(){return 1;}
	static Int_t Z(){return 2;}
	Bool_t Pass(NicaTwoTrack *pair);
	virtual Bool_t  Init(Int_t format_id=0);
	virtual ~MpdPairTpcEntrySimpleCut();
	ClassDef(MpdPairTpcEntrySimpleCut,1)
};

#endif /* 2018_DOKTORAT_KODY_NICA_MPD_CUTS_PAIRCUTS_SIMPLE_TPC_MPDPAIRTPCENTRYSIMPLECUT_H_ */
