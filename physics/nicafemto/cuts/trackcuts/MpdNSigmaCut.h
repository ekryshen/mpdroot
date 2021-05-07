/*
 * MpdNSigmaCut.h
 *
 *  Created on: 28 gru 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDNSIGMACUT_H_
#define INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDNSIGMACUT_H_

#include "NicaTrackCut.h"

class TF1;

class MpdNSigmaCut : public NicaTrackCut{
	TF1 *parElBB; TF1 *parPiBB; TF1 *parKaBB; TF1 *parPrBB;
	TF1 *fAsymmetryElLowP; TF1 *fAsymmetryElMidP; TF1 *fAsymmetryElHighP; /// electrons
	TF1 *fAsymmetryPiLowP; TF1 *fAsymmetryPiMidP; TF1 *fAsymmetryPiHighP; /// pions
	TF1 *fAsymmetryKaLowP; TF1 *fAsymmetryKaMidP; TF1 *fAsymmetryKaHighP; /// kaons
	TF1 *fAsymmetryPrLowP; TF1 *fAsymmetryPrMidP; TF1 *fAsymmetryPrHighP; /// protons
	TF1* elSigmaLowP; TF1* elSigmaMidP; TF1* elSigmaHighP;
	TF1* piSigmaLowP; TF1* piSigmaMidP; TF1* piSigmaHighP;
	TF1* prSigmaLowP; TF1* prSigmaHighP;
	TF1* kaSigmaLowP; TF1* kaSigmaHighP;
	Double_t GetDedxWidthValue(Double_t p, Int_t specie);
public:
	MpdNSigmaCut();
	Bool_t Pass(NicaTrack *track);
	virtual ~MpdNSigmaCut();
	static Int_t SigmaPion(){return 0;};
	static Int_t SigmaKaon(){return 1;};
	static Int_t SigmaProton(){return 2;};
	static Int_t SigmaElectron(){return 3;};
	ClassDef(MpdNSigmaCut,1)
};

#endif /* INTERFACES_MPDROOT_CUTS_TRACKCUTS_MPDNSIGMACUT_H_ */
