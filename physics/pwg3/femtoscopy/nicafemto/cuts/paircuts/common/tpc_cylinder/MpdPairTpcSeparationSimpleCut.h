/*
 * MpdPairTpcSeparatioSimpleCut.h
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICA_MPD_CUTS_PAIRCUTS_MPDPAIRTPCSEPARATIOSIMPLECUT_H_
#define NICA_MPD_CUTS_PAIRCUTS_MPDPAIRTPCSEPARATIOSIMPLECUT_H_
#include "MpdPairTpcSimpleCut.h"
/**
 * class for simple calculation of TPC entrance, exit and average separation distances
 * use approximation of TPC as a tube not modular structure
 */

class MpdPairTpcSeparationSimpleCut : public MpdPairTpcSimpleCut{
protected:
	Float_t fSteps;
	Int_t fNSteps;
	Float_t fRMin;
	Float_t fRMax;
	Float_t fZMax;
	Float_t *fRadii; //[fNSteps]
public:
	MpdPairTpcSeparationSimpleCut(Int_t averaging_steps=10);
	static Int_t Entry(){return 0;};
	static Int_t Exit(){return 1;}
	static Int_t Average(){return 2;};
	static Int_t Min(){return 3;};
	static Int_t Max(){return 4;}
	virtual Bool_t Pass(NicaTwoTrack *pair);
	void SetNSteps(Int_t n){fSteps = n;};
	virtual Bool_t  Init(Int_t format_id=0);
	virtual ~MpdPairTpcSeparationSimpleCut();
	ClassDef(MpdPairTpcSeparationSimpleCut,1)
};

class MpdPairTcpSeparationSimpleNegCut : public MpdPairTpcSeparationSimpleCut{
public:
	MpdPairTcpSeparationSimpleNegCut();
	Bool_t Pass(NicaTwoTrack *pair);
	virtual ~MpdPairTcpSeparationSimpleNegCut();
	ClassDef(MpdPairTcpSeparationSimpleNegCut,1)
};

#endif /* 2018_DOKTORAT_KODY_NICA_MPD_CUTS_PAIRCUTS_MPDPAIRTPCSEPARATIOSIMPLECUT_H_ */
