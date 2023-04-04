//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class ParticleIdentification
//      see ParticleIdentification.h for details
//
// Environment:
//      Software developed for MPD at NICA.
//
// Author List:
//      Sergey Merts
//
//-----------------------------------------------------------

// C/C++ Headers ----------------------
#include <iostream>
#include <TMath.h>
#include <TF1.h>
// This Class' Header ------------------
#include "MpdParticleIdentification.h"

using namespace std;
using namespace TMath;

// Class Member definitions -----------

MpdParticleIdentification::MpdParticleIdentification() {}

MpdParticleIdentification::~MpdParticleIdentification() {}

// Int_t MpdParticleIdentification::GetTpcProbs(Float_t P, Float_t dedx, Int_t nHits, Float_t& Ppi, Float_t& Pk,
// Float_t& Pp, Float_t& Pe, Int_t method) {
//
//    /*Parameters for fitting*/
//    //AZ Float_t *ProtonPar, *PionPar, *KaonPar, *ElectronPar;
//    /************************/
//
//    dedx *= 1000000; // converting to keV/cm
//    const Int_t Ntypes = 4; //order: pion, kaon, proton, electron
//    Float_t ProtonPar[9], PionPar[9], KaonPar[9], ElectronPar[9], resultProb[Ntypes]; //AZ
//    const Int_t Nintervals = 10;
//    Float_t sigmas[Ntypes][Nintervals] = {//array of sigmas for different particles and different intervals of
//    momentum. They were got from 100000 BOX events
//        {0.113, 0.107, 0.100, 0.107, 0.097, 0.111, 0.122, 0.105, 0.115, 0.118},
//        {0.407, 0.276, 0.260, 0.159, 0.185, 0.154, 0.122, 0.105, 0.115, 0.118},
//        {0.507, 0.488, 0.394, 0.330, 0.268, 0.244, 0.260, 0.172, 0.193, 0.118},
//        {0.163, 0.160, 0.149, 0.159, 0.172, 0.205, 0.260, 0.172, 0.193, 0.156}
//    };
//
//    Float_t sigma = 1.0 / TMath::Sqrt(nHits); //for gaussian
//    Float_t sum = 0.0; // for normalizing
//
//    /*A priori coefficients for Bayes approach */
//    Float_t bayesAprioriCoefficients[Ntypes];
//    /*A "measured" probabilities */
//    Float_t gausProb[Ntypes];
//
//    Float_t sigPi, sigPr, sigKa, sigEl; //sigmas for gaussians
//
//    switch (method) {
//
//        case 0: //bethe-bloch approximation with "standard" sigmas and equal bayesian coefficients
//
//            //parameters were got from Bethe-Bloch approximation for 100000 BOX events
//	  /*AZ
//            ProtonPar = new Float_t[4];
//            PionPar = new Float_t[4];
//            KaonPar = new Float_t[4];
//            ElectronPar = new Float_t[4];
//	  */
//            ProtonPar[0] = -3.957;
//            ProtonPar[1] = 3.525;
//            ProtonPar[2] = 16.468;
//            ProtonPar[3] = 0.815;
//            ProtonPar[4] = 5.247;
//            PionPar[0] = -0.739;
//            PionPar[1] = 7.373;
//            PionPar[2] = 3904.790;
//            PionPar[3] = 0.274;
//            PionPar[4] = 5.497;
//            KaonPar[0] = -2.590;
//            KaonPar[1] = 4.918;
//            KaonPar[2] = 79.722;
//            KaonPar[3] = 0.357;
//            KaonPar[4] = 4.511;
//            ElectronPar[0] = -1.552;
//            ElectronPar[1] = 1.748;
//            ElectronPar[2] = 7.425;
//            ElectronPar[3] = 0.980;
//            ElectronPar[4] = 1.604;
//
//            sigma = 0.07 * dedx;
//            sigPi = sigPr = sigKa = sigEl = sigma;
//            bayesAprioriCoefficients[0] = 1.0;
//            bayesAprioriCoefficients[1] = 1.0;
//            bayesAprioriCoefficients[2] = 1.0;
//            bayesAprioriCoefficients[3] = 1.0;
//
//            gausProb[0] = Gaus(dedx, BetheBlochFunction(P, PionPar), sigPi, kTRUE);
//            gausProb[1] = Gaus(dedx, BetheBlochFunction(P, KaonPar), sigKa, kTRUE);
//            gausProb[2] = Gaus(dedx, BetheBlochFunction(P, ProtonPar), sigPr, kTRUE);
//            gausProb[3] = Gaus(dedx, BetheBlochFunction(P, ElectronPar), sigEl, kTRUE);
//
//            sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
//            if (sum == 0.) return 1;
//            gausProb[0] /= sum;
//            gausProb[1] /= sum;
//            gausProb[2] /= sum;
//            gausProb[3] /= sum;
//            break;
//
//        case 1: //bethe-bloch approximation with special different sigmas and byesian coefficients
//
//            //parameters were got from Bethe-Bloch approximation for 100000 BOX events
//	  /*AZ
//            ProtonPar = new Float_t[4];
//            PionPar = new Float_t[4];
//            KaonPar = new Float_t[4];
//            ElectronPar = new Float_t[4];
//	  */
//            ProtonPar[0] = -3.957;
//            ProtonPar[1] = 3.525;
//            ProtonPar[2] = 16.468;
//            ProtonPar[3] = 0.815;
//            ProtonPar[4] = 5.247;
//            PionPar[0] = -0.739;
//            PionPar[1] = 7.373;
//            PionPar[2] = 3904.790;
//            PionPar[3] = 0.274;
//            PionPar[4] = 5.497;
//            KaonPar[0] = -2.590;
//            KaonPar[1] = 4.918;
//            KaonPar[2] = 79.722;
//            KaonPar[3] = 0.357;
//            KaonPar[4] = 4.511;
//            ElectronPar[0] = -1.552;
//            ElectronPar[1] = 1.748;
//            ElectronPar[2] = 7.425;
//            ElectronPar[3] = 0.980;
//            ElectronPar[4] = 1.604;
//
//            if (P < 0.3) {
//                bayesAprioriCoefficients[0] = 0.28;
//                bayesAprioriCoefficients[1] = 0.25;
//                bayesAprioriCoefficients[2] = 0.26;
//                bayesAprioriCoefficients[3] = 0.21;
//                sigPi = sigmas[0][0];
//                sigKa = sigmas[1][0];
//                sigPr = sigmas[2][0];
//                sigEl = sigmas[3][0];
//            } else if ((P >= 0.3) && (P < 0.4)) {
//                bayesAprioriCoefficients[0] = 0.91;
//                bayesAprioriCoefficients[1] = 0.02;
//                bayesAprioriCoefficients[2] = 0.06;
//                bayesAprioriCoefficients[3] = 0.01;
//                sigPi = sigmas[0][1];
//                sigKa = sigmas[1][1];
//                sigPr = sigmas[2][1];
//                sigEl = sigmas[3][1];
//            } else if ((P >= 0.4) && (P < 0.5)) {
//                bayesAprioriCoefficients[0] = 0.70;
//                bayesAprioriCoefficients[1] = 0.11;
//                bayesAprioriCoefficients[2] = 0.13;
//                bayesAprioriCoefficients[3] = 0.06;
//                sigPi = sigmas[0][2];
//                sigKa = sigmas[1][2];
//                sigPr = sigmas[2][2];
//                sigEl = sigmas[3][2];
//            } else if ((P >= 0.5) && (P < 0.6)) {
//                bayesAprioriCoefficients[0] = 0.39;
//                bayesAprioriCoefficients[1] = 0.19;
//                bayesAprioriCoefficients[2] = 0.32;
//                bayesAprioriCoefficients[3] = 0.10;
//                sigPi = sigmas[0][3];
//                sigKa = sigmas[1][3];
//                sigPr = sigmas[2][3];
//                sigEl = sigmas[3][3];
//            } else if ((P >= 0.6) && (P < 0.7)) {
//                bayesAprioriCoefficients[0] = 0.41;
//                bayesAprioriCoefficients[1] = 0.17;
//                bayesAprioriCoefficients[2] = 0.37;
//                bayesAprioriCoefficients[3] = 0.5;
//                sigPi = sigmas[0][4];
//                sigKa = sigmas[1][4];
//                sigPr = sigmas[2][4];
//                sigEl = sigmas[3][4];
//            } else if ((P >= 0.7) && (P < 0.8)) {
//                bayesAprioriCoefficients[0] = 0.43;
//                bayesAprioriCoefficients[1] = 0.16;
//                bayesAprioriCoefficients[2] = 0.31;
//                bayesAprioriCoefficients[3] = 0.10;
//                sigPi = sigmas[0][5];
//                sigKa = sigmas[1][5];
//                sigPr = sigmas[2][5];
//                sigEl = sigmas[3][5];
//            } else if ((P >= 0.8) && (P < 0.9)) {
//                bayesAprioriCoefficients[0] = 0.44;
//                bayesAprioriCoefficients[1] = 0.11;
//                bayesAprioriCoefficients[2] = 0.43;
//                bayesAprioriCoefficients[3] = 0.02;
//                sigPi = sigmas[0][6];
//                sigKa = sigmas[1][6];
//                sigPr = sigmas[2][6];
//                sigEl = sigmas[3][6];
//            } else if ((P >= 0.9) && (P < 1.0)) {
//                bayesAprioriCoefficients[0] = 0.36;
//                bayesAprioriCoefficients[1] = 0.21;
//                bayesAprioriCoefficients[2] = 0.36;
//                bayesAprioriCoefficients[3] = 0.07;
//                sigPi = sigmas[0][7];
//                sigKa = sigmas[1][7];
//                sigPr = sigmas[2][7];
//                sigEl = sigmas[3][7];
//            } else if ((P >= 1.0) && (P < 1.2)) {
//                bayesAprioriCoefficients[0] = 0.32;
//                bayesAprioriCoefficients[1] = 0.32;
//                bayesAprioriCoefficients[2] = 0.32;
//                bayesAprioriCoefficients[3] = 0.04;
//                sigPi = sigmas[0][8];
//                sigKa = sigmas[1][8];
//                sigPr = sigmas[2][8];
//                sigEl = sigmas[3][8];
//            } else if (P >= 1.2) {
//                bayesAprioriCoefficients[0] = 0.34;
//                bayesAprioriCoefficients[1] = 0.27;
//                bayesAprioriCoefficients[2] = 0.31;
//                bayesAprioriCoefficients[3] = 0.08;
//                sigPi = sigmas[0][9];
//                sigKa = sigmas[1][9];
//                sigPr = sigmas[2][9];
//                sigEl = sigmas[3][9];
//            }
//
//            gausProb[0] = Gaus(dedx, BetheBlochFunction(P, PionPar), sigPi, kTRUE);
//            gausProb[1] = Gaus(dedx, BetheBlochFunction(P, KaonPar), sigKa, kTRUE);
//            gausProb[2] = Gaus(dedx, BetheBlochFunction(P, ProtonPar), sigPr, kTRUE);
//            gausProb[3] = Gaus(dedx, BetheBlochFunction(P, ElectronPar), sigEl, kTRUE);
//
//            sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
//            if (sum == 0.) return 1;
//            gausProb[0] /= sum;
//            gausProb[1] /= sum;
//            gausProb[2] /= sum;
//            gausProb[3] /= sum;
//            break;
//
//        case 2: //bethe-bloch approximation with "standard" sigmas and different byesian coefficients
//
//            //parameters were got from Bethe-Bloch approximation for 100000 BOX events
//	  /*AZ
//            ProtonPar = new Float_t[4];
//            PionPar = new Float_t[4];
//            KaonPar = new Float_t[4];
//            ElectronPar = new Float_t[4];
//	  */
//            ProtonPar[0] = -3.957;
//            ProtonPar[1] = 3.525;
//            ProtonPar[2] = 16.468;
//            ProtonPar[3] = 0.815;
//            ProtonPar[4] = 5.247;
//            PionPar[0] = -0.739;
//            PionPar[1] = 7.373;
//            PionPar[2] = 3904.790;
//            PionPar[3] = 0.274;
//            PionPar[4] = 5.497;
//            KaonPar[0] = -2.590;
//            KaonPar[1] = 4.918;
//            KaonPar[2] = 79.722;
//            KaonPar[3] = 0.357;
//            KaonPar[4] = 4.511;
//            ElectronPar[0] = -1.552;
//            ElectronPar[1] = 1.748;
//            ElectronPar[2] = 7.425;
//            ElectronPar[3] = 0.980;
//            ElectronPar[4] = 1.604;
//
//            sigma = 0.07 * dedx;
//            sigPi = sigPr = sigKa = sigEl = sigma;
//
//            if (P < 0.3) {
//                bayesAprioriCoefficients[0] = 0.28;
//                bayesAprioriCoefficients[1] = 0.25;
//                bayesAprioriCoefficients[2] = 0.26;
//                bayesAprioriCoefficients[3] = 0.21;
//            } else if ((P >= 0.3) && (P < 0.4)) {
//                bayesAprioriCoefficients[0] = 0.91;
//                bayesAprioriCoefficients[1] = 0.02;
//                bayesAprioriCoefficients[2] = 0.06;
//                bayesAprioriCoefficients[3] = 0.01;
//            } else if ((P >= 0.4) && (P < 0.5)) {
//                bayesAprioriCoefficients[0] = 0.70;
//                bayesAprioriCoefficients[1] = 0.11;
//                bayesAprioriCoefficients[2] = 0.13;
//                bayesAprioriCoefficients[3] = 0.06;
//            } else if ((P >= 0.5) && (P < 0.6)) {
//                bayesAprioriCoefficients[0] = 0.39;
//                bayesAprioriCoefficients[1] = 0.19;
//                bayesAprioriCoefficients[2] = 0.32;
//                bayesAprioriCoefficients[3] = 0.10;
//            } else if ((P >= 0.6) && (P < 0.7)) {
//                bayesAprioriCoefficients[0] = 0.41;
//                bayesAprioriCoefficients[1] = 0.17;
//                bayesAprioriCoefficients[2] = 0.37;
//                bayesAprioriCoefficients[3] = 0.5;
//            } else if ((P >= 0.7) && (P < 0.8)) {
//                bayesAprioriCoefficients[0] = 0.43;
//                bayesAprioriCoefficients[1] = 0.16;
//                bayesAprioriCoefficients[2] = 0.31;
//                bayesAprioriCoefficients[3] = 0.10;
//            } else if ((P >= 0.8) && (P < 0.9)) {
//                bayesAprioriCoefficients[0] = 0.44;
//                bayesAprioriCoefficients[1] = 0.11;
//                bayesAprioriCoefficients[2] = 0.43;
//                bayesAprioriCoefficients[3] = 0.02;
//            } else if ((P >= 0.9) && (P < 1.0)) {
//                bayesAprioriCoefficients[0] = 0.36;
//                bayesAprioriCoefficients[1] = 0.21;
//                bayesAprioriCoefficients[2] = 0.36;
//                bayesAprioriCoefficients[3] = 0.07;
//            } else if ((P >= 1.0) && (P < 1.2)) {
//                bayesAprioriCoefficients[0] = 0.32;
//                bayesAprioriCoefficients[1] = 0.32;
//                bayesAprioriCoefficients[2] = 0.32;
//                bayesAprioriCoefficients[3] = 0.04;
//            } else if (P >= 1.2) {
//                bayesAprioriCoefficients[0] = 0.34;
//                bayesAprioriCoefficients[1] = 0.27;
//                bayesAprioriCoefficients[2] = 0.31;
//                bayesAprioriCoefficients[3] = 0.08;
//            }
//
//            gausProb[0] = Gaus(dedx, BetheBlochFunction(P, PionPar), sigPi, kTRUE);
//            gausProb[1] = Gaus(dedx, BetheBlochFunction(P, KaonPar), sigKa, kTRUE);
//            gausProb[2] = Gaus(dedx, BetheBlochFunction(P, ProtonPar), sigPr, kTRUE);
//            gausProb[3] = Gaus(dedx, BetheBlochFunction(P, ElectronPar), sigEl, kTRUE);
//
//            sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
//            if (sum == 0.) return 1;
//            gausProb[0] /= sum;
//            gausProb[1] /= sum;
//            gausProb[2] /= sum;
//            gausProb[3] /= sum;
//            break;
//
//        case 3: //parabolic approximation
//
//            //parameters were got from parabolic approximation function for 2000 UrQMD events
//            ProtonPar[0] = 0.031;
//            ProtonPar[1] = -0.124;
//            ProtonPar[2] = 1.412;
//            PionPar[0] = 0.230;
//            PionPar[1] = 0.088;
//            PionPar[2] = 1.072;
//            KaonPar[0] = 0.676;
//            KaonPar[1] = 0.621;
//            KaonPar[2] = 0.831;
//            ElectronPar[0] = 0.000;
//            ElectronPar[1] = -0.018;
//            ElectronPar[2] = 2.055;
//
//            Float_t invP = 1.0 / P;
//
//            gausProb[0] = Gaus(dedx, ParabolicFunction(invP, PionPar), sigPi, kTRUE);
//            gausProb[1] = Gaus(dedx, ParabolicFunction(invP, KaonPar), sigKa, kTRUE);
//            gausProb[2] = Gaus(dedx, ParabolicFunction(invP, ProtonPar), sigPr, kTRUE);
//            gausProb[3] = Gaus(dedx, ParabolicFunction(invP, ElectronPar), sigEl, kTRUE);
//
//            sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
//            if (sum == 0.) return 1;
//            gausProb[0] /= sum;
//            gausProb[1] /= sum;
//            gausProb[2] /= sum;
//            gausProb[3] /= sum;
//
//            break;
//    }
//    //AZ Float_t *resultProb = new Float_t[Ntypes];
//    BayesFunction(gausProb, bayesAprioriCoefficients, resultProb, Ntypes);
//
//    Ppi = resultProb[0];
//    Pk = resultProb[1];
//    Pp = resultProb[2];
//    Pe = resultProb[3];
//
//    return 0;
//}

/*
Int_t MpdParticleIdentification::GetTofProbs(Float_t P, Float_t beta, Float_t& Ppi, Float_t& Pk, Float_t& Pp, Float_t&
Pe, Int_t method) {

    Float_t ProtonParHyper[2], PionParHyper[2], KaonParHyper[2], ElectronParHyper[2];
    ProtonParHyper[0] = 0.938 * 0.938;
    ProtonParHyper[1] = -0.065; //-0.111791

    PionParHyper[0] = 0.139 * 0.139;
    PionParHyper[1] = -0.042;

    KaonParHyper[0] = 0.494 * 0.494;
    KaonParHyper[1] = -0.0491; //-0.0717

    ElectronParHyper[0] = 0.000511 * 0.000511;
    ElectronParHyper[1] = -0.028;

    const Int_t Ntypes = 4; //order: pion, kaon, proton, electron
    const Int_t Nintervals = 72;
    Double_t P_int[Nintervals + 1] =
{0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.6,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.9,0.93,0.96,0.99,1.02,1.05,1.08,
1.11,1.14,1.17,1.2,1.23,1.26,1.29,1.32,1.35,1.38,1.41,1.44,1.47,1.5,1.53,1.56,1.59,1.62,1.65,1.68,1.71,1.74,1.77,1.8,1.83,1.86,1.89,1.92,1.95,1.98,2.01,2.04,2.07,2.1,
2.13,2.16,2.19,2.22, 3};



    if (Abs(beta) < 1e-6) return 1;
    const Float_t oneOverBeta = 1.0 / beta;

    Float_t sigmas[Ntypes][Nintervals] = {//array of sigmas for different particles and different intervals of momentum.
      {0.0693274, 0.0693274, 0.0544018, 0.0434532, 0.0383988, 0.0374991, 0.0346434, 0.03343, 0.0319856, 0.032,
0.0315147, 0.0310805, 0.0298654, 0.0303749, 0.0301869, 0.0296072, 0.0289299, 0.0285614, 0.0298113, 0.0287559, 0.0269426,
0.0276648, 0.0281308, 0.0289822, 0.0260458, 0.0274102, 0.0276939, 0.0264013, 0.0251029, 0.0270367, 0.0269619, 0.0278648,
0.0247407, 0.0273073, 0.0262839, 0.0284527, 0.024663, 0.0270564, 0.0232594, 0.0256739, 0.0257952, 0.0237649, 0.0239696,
0.0248299, 0.0283596, 0.0264846, 0.022326, 0.0208957, 0.0307023, 0.0284478, 0.0229504, 0.0221302, 0.0251733, 0.0258165,
0.0273607, 0.0246405, 0.0229405, 0.024879, 0.0270787, 0.0238393, 0.026333, 0.0112759, 0.0404941, 0.0184759, 0.0225569,
0.0327113, 0.0221111, 0.0187046, 0.0194951, 0.0189474, 0.0211018, 0.00812245},

      {0.213016, 0.213016, 0.213016, 0.213016, 0.213016, 0.0989433, 0.0834294, 0.0675234, 0.0595375, 0.058213,
0.0486858, 0.0448145, 0.0471928, 0.0466786, 0.0451925, 0.0373988, 0.0451043, 0.0343869, 0.039328, 0.038438, 0.0364276,
0.034408, 0.0344429, 0.03274, 0.0368628, 0.0303875, 0.0312515, 0.0331831, 0.0302307, 0.0310654, 0.0382999, 0.0311076,
0.0396424, 0.0350253, 0.0319718, 0.0367985, 0.0300645, 0.0263978, 0.0261712, 0.0314639, 0.0263756, 0.0338912, 0.0313582,
0.0200379, 0.0265785, 0.0300492, 0.0308032, 0.03746, 0.0220187, 0.0222599, 0.0244639, 0.0284196, 0.0286862, 0.0407272,
0.012562, 0.023188, 0.012329, 0.0297552, 0.0208175, 0.0426947, 0.0267256, 0.0157099, 0.0404566, 0.0162424, 0.0144609,
0.0304202, 0.019353, 0.0527709, 0.013652, 0.187007, 0.28284, 0.0161017},

      {0.220549, 0.220549, 0.220549, 0.220549, 0.220549, 0.220549, 0.220549, 0.220549, 0.220549, 0.094147, 0.0942735,
0.0958616, 0.0738113, 0.0776133, 0.0631406, 0.0543964, 0.0672916, 0.0562169, 0.0625832, 0.0549343, 0.0521649, 0.0535391,
0.0529288, 0.0466422, 0.0523035, 0.0472841, 0.0406543, 0.0460797, 0.0530522, 0.0439447, 0.0426407, 0.0495772, 0.0406714,
0.0341736, 0.0376888, 0.0472375, 0.0441835, 0.0363689, 0.0425563, 0.0421598, 0.0371204, 0.0462698, 0.0360139, 0.0384682,
0.0400049, 0.0387571, 0.0328137, 0.0384007, 0.0351909, 0.0594957, 0.0310895, 0.0357163, 0.0387537, 0.0358798, 0.0385559,
0.0277795, 0.0348187, 0.0269486, 0.0523929, 0.0282686, 0.0339755, 0.0238928, 0.0241012, 0.0337723, 0.0334852, 0.0248415,
0.0534019, 0.0270102, 0.0403495, 0.0284284, 0.0289664, 0.0240482}, {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
0.01, 0.01, 0.01, 0.01}
    };

    Int_t i_p = 0; // momentum interval

    Float_t sum = 0.0; // for normalizing
    Float_t gausProb[Ntypes] = {0.0, 0.0, 0.0, 0.0};

    switch (method) {

        case 0: //Hyperbolic (1/beta vs p)

            if (P > 3) i_p = 71;
            for (Int_t k = 0; k < Nintervals; k++) if (P >= P_int[k] && P < P_int[k + 1]) i_p = k;

            gausProb[0] = Gaus(oneOverBeta, HyperbolicFunction(P, PionParHyper), sigmas[0][i_p], kTRUE);
            gausProb[1] = Gaus(oneOverBeta, HyperbolicFunction(P, KaonParHyper), sigmas[1][i_p], kTRUE);
            gausProb[2] = Gaus(oneOverBeta, HyperbolicFunction(P, ProtonParHyper), sigmas[2][i_p], kTRUE);
            gausProb[3] = Gaus(oneOverBeta, HyperbolicFunction(P, ElectronParHyper), sigmas[3][i_p], kTRUE);

            sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
            if (sum == 0.) return 1;
            gausProb[0] /= sum;
            gausProb[1] /= sum;
            gausProb[2] /= sum;
            gausProb[3] /= sum;
            break;

        case 1:

            break;
    }

    Ppi = gausProb[0];
    Pk = gausProb[1];
    Pp = gausProb[2];
    Pe = gausProb[3];

    return 0;
}
*/

Int_t MpdParticleIdentification::GetCombinedProbs(Float_t *tofProbs, Float_t *tpcProbs, Float_t *resultProbs, Int_t N)
{

   /*tofProbs    - measured probabilities from TOF*/
   /*tpcProbs    - measured probabilities from TPC*/
   /*resultProbs - combined probabilities */
   /*N           - number of types {pion, kaon, proton, electron} */

   Float_t sumBayes = 0.0;
   for (Int_t i = 0; i < N; ++i) {
      sumBayes += tofProbs[i] * tpcProbs[i];
   }

   for (Int_t i = 0; i < N; ++i) {
      if ((tpcProbs[i] != 0.0) && (tofProbs[i] != 0.0)) {
         resultProbs[i] = tofProbs[i] * tpcProbs[i] / sumBayes;
      } else if (tpcProbs[i] == 0.0) {
         resultProbs[i] = tofProbs[i] / sumBayes;
      } else if (tofProbs[i] == 0.0) {
         resultProbs[i] = tpcProbs[i] / sumBayes;
      }
   }

   return 0;
}

// Float_t MpdParticleIdentification::BetheBlochFunction(Float_t x, Float_t *p) {
//    Float_t b = 1 / (x / Sqrt(1 + x * x));
//    return p[0] / Power(b, p[3]) * (p[1] - Power(b, p[3]) - Log(p[2] + Power(1 / x, p[4])));
//}

Float_t MpdParticleIdentification::ParabolicFunction(Float_t x, Float_t *p)
{
   return p[0] * x * x + p[1] * x + p[2];
}

/*
Float_t MpdParticleIdentification::HyperbolicFunction(Float_t x, Float_t *p) {
    return Sqrt(1 + p[0] / (x * x)) + p[1];  //p[0]=m^2, p[1]=additive const
}
*/

/*
// Bayes function for probabilities calculation //
Int_t MpdParticleIdentification::BayesFunction(Float_t *measProb, Float_t *aprioriProb, Float_t *bayesProb, Int_t N) {

    // measProb    - measured probabilities from TPC/TOF/etc //
    // aprioriProb - a priori probabilities from some magic //
    // bayesProb   - result bayes probabilities //
    // N           - number of types {pion, kaon, proton, electron} //

    Float_t sumBayes = 0.0;
    for (Int_t i = 0; i < N; ++i) {
        sumBayes += measProb[i] * aprioriProb[i];
    }

    if (sumBayes == 0.) return 1;

    for (Int_t i = 0; i < N; ++i) {
        bayesProb[i] = measProb[i] * aprioriProb[i] / sumBayes;
    }

    return 0;
}
*/
ClassImp(MpdParticleIdentification);
