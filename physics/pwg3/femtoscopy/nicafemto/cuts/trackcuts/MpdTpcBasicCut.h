/*
 * NicaTrackTpcSigmaCut.h
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAROOT_CUTS_TRACKCUTS_NICATRACKTPCBASICCUT_H_
#define NICAROOT_CUTS_TRACKCUTS_NICATRACKTPCBASICCUT_H_


#include <TString.h>

#include "NicaTrackTpcCut.h"

class MpdTcpCalibratedCut : public NicaTrackTpcCut{
private:
	Double_t fKaonBB[6];
	Double_t fPionBB[6];
	Double_t fElectronBB[6];
	Double_t fProtonBB[6];
	Double_t fSigmasPi[72];
	Double_t fSigmasP[72];
	Double_t fSigmasK[72];
	Double_t fSigmasE[72];
	Double_t fPInt[72];
	Bool_t fUseDst;//use values from DST files
	Double_t  BetheBlochFunction(Double_t X, Double_t *par);
public:
	MpdTcpCalibratedCut(TString calib_file="");
	virtual Bool_t Pass(NicaTrack *track);
	virtual ~MpdTcpCalibratedCut(){};
	ClassDef(MpdTcpCalibratedCut,1);
};
#endif /* NICAROOT_CUTS_TRACKCUTS_NICATRACKTPCBASICCUT_H_ */
