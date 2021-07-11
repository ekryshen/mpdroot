/*
 * MpdNSigmaCut.cxx
 *
 *  Created on: 28 gru 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdNSigmaCut.h"
#include "NicaExpTrack.h"
#include "NicaMpdConst.h"
#include "NicaTpcTrack.h"
#include <TF1.h>

MpdNSigmaCut::MpdNSigmaCut() : NicaTrackCut(4) {
  SetUnitName("pion_sigma [#sigma]", SigmaPion());
  SetUnitName("kaon_sigma [#sigma]", SigmaKaon());
  SetUnitName("proton_sigma [#sigma]", SigmaProton());
  SetUnitName("electron_sigma [#sigma]", SigmaElectron());
  for (int i = 0; i < 4; i++) SetMinMax(-1E+9, 1E+9, i);
  Double_t PMIN = 0., PMAX = 3.0;
  parElBB = new TF1("parElBB", "[0]/x+[1]", PMIN, PMAX);
  parPiBB = new TF1("parPiBB",
                    "[0]/pow(x/sqrt(x*x+0.01949),[3])*([1]-pow(x/"
                    "sqrt(x*x+0.01949),[3])-log([2]+pow(1./(x/0.1396),[4])) )",
                    PMIN, PMAX);
  parKaBB = new TF1("parKaBB",
                    "[0]/pow(x/sqrt(x*x+0.2437),[3])*([1]-pow(x/"
                    "sqrt(x*x+0.2437),[3])-log([2]+pow(1./(x/0.4937),[4])) )",
                    PMIN, PMAX);
  parPrBB = new TF1("parPrBB",
                    "[0]/pow(x/sqrt(x*x+0.88),[3])*([1]-pow(x/"
                    "sqrt(x*x+0.88),[3])-log([2]+pow(1./(x/0.9383),[4])) )",
                    PMIN, PMAX);
  elSigmaLowP = new TF1("elSigmaLowP", "pol1(0)", 0., 0.55);
  elSigmaLowP->SetParameters(0.0747217, -0.0308101);
  elSigmaMidP = new TF1("elSigmaMidP", "pol1(0)", 0.55, 2.0);
  elSigmaMidP->SetParameters(0.0653074, -0.00546004);
  elSigmaHighP = new TF1("elSigmaHighP", "pol1(0)", 2.0, 3.0);
  elSigmaHighP->SetParameters(0.0572145, -0.00104922);
  piSigmaLowP = new TF1("piSigmaLowP", "pol1(0)", 0., 0.25);
  piSigmaLowP->SetParameters(0.144972, -0.308612);
  piSigmaMidP = new TF1("piSigmaMidP", "pol1(0)", 0.25, 0.4);
  piSigmaMidP->SetParameters(0.0805, -0.0553436);
  piSigmaHighP = new TF1("piSigmaHighP", "pol1(0)", 0.4, 3.0);
  piSigmaHighP->SetParameters(0.0587877, 0.00143879);
  prSigmaLowP = new TF1("prSigmaLowP", "pol1(0)", 0., 0.2);
  prSigmaLowP->SetParameters(0.157451, -0.481992);
  prSigmaHighP = new TF1("prSigmaHighP", "pol1(0)", 0.2, 3.0);
  prSigmaHighP->SetParameters(0.0623534, 0.00287991);
  kaSigmaLowP = new TF1("kaSigmaLowP", "pol1(0)", 0., 0.2);
  kaSigmaLowP->SetParameters(0.120409, -0.289855);
  kaSigmaHighP = new TF1("kaSigmaHighP", "pol1(0)", 0.2, 3.0);
  kaSigmaHighP->SetParameters(0.0630206, -0.00374603);
  fAsymmetryElLowP = new TF1("fAsymmetryElLowP", "pol1(0)", 0., 0.5);
  fAsymmetryElLowP->SetParameters(-0.340074, 1.04407);
  fAsymmetryElMidP = new TF1("fAsymmetryElMidP", "pol1(0)", 0.5, 0.85);
  fAsymmetryElMidP->SetParameters(0.291559, -0.306138);
  fAsymmetryElHighP = new TF1("fAsymmetryElHighP", "pol1(0)", 0.85, 3.0);
  fAsymmetryElHighP->SetParameters(0.0915437, 0.0467649);
  fAsymmetryPiLowP = new TF1("fAsymmetryPiLowP", "pol1(0)", 0., 0.18);
  fAsymmetryPiLowP->SetParameters(3.50026, -18.1992);
  fAsymmetryPiMidP = new TF1("fAsymmetryPiMidP", "pol1(0)", 0.18, 0.9);
  fAsymmetryPiMidP->SetParameters(0.106035, 0.237438);
  fAsymmetryPiHighP = new TF1("fAsymmetryPiHighP", "pol1(0)", 0.9, 3.);
  fAsymmetryPiHighP->SetParameters(0.256077, -0.0069969);
  fAsymmetryPrLowP = new TF1("fAsymmetryPrLowP", "pol1(0)", 0., 0.25);
  fAsymmetryPrLowP->SetParameters(1.68519, -2.75816);
  fAsymmetryPrMidP = new TF1("fAsymmetryPrMidP", "pol1(0)", 0.25, 0.5);
  fAsymmetryPrMidP->SetParameters(1.68519, -2.75816);
  fAsymmetryPrHighP = new TF1("fAsymmetryPrHighP", "pol1(0)", 0.5, 3.);
  fAsymmetryPrHighP->SetParameters(-0.0761382, 0.186669);
  fAsymmetryKaLowP = new TF1("fAsymmetryKaLowP", "pol1(0)", 0., 0.5);
  fAsymmetryKaLowP->SetParameters(2.3588, -4.87114);
  fAsymmetryKaMidP = new TF1("fAsymmetryKaMidP", "pol1(0)", 0.5, 1.05);
  fAsymmetryKaMidP->SetParameters(0.00526557, 0.239378);
  fAsymmetryKaHighP = new TF1("fAsymmetryKaHighP", "pol1(0)", 1.05, 3.);
  fAsymmetryKaHighP->SetParameters(0.0892569, 0.119403);

  Double_t dedxParam = 0;
  Double_t fKoef = 1;
  dedxParam = fKoef * (-65.7432);
  parElBB->SetParameter(0, dedxParam);
  dedxParam = fKoef * 4007.68;
  parElBB->SetParameter(1, dedxParam);
  /// pions
  dedxParam = fKoef * (-333.199);
  parPiBB->SetParameters(dedxParam, -4.76053, 3.73173, 1.82109, -1.2414);
  /// kaons
  dedxParam = fKoef * (-674.513);
  parKaBB->SetParameters(dedxParam, -4.02511, 0.391596, 0.993733, 4.13295);
  /// protons
  dedxParam = fKoef * (-2328.44);
  parPrBB->SetParameters(dedxParam, -0.054584, 1.15196, 0.819829, 2.74767);
}

Bool_t MpdNSigmaCut::Pass(NicaTrack* track) {
  NicaTpcTrack* tpc =
      (NicaTpcTrack*)((NicaExpTrack*)track)->GetDetTrack(MpdDetectorID::kTPC);
  Double_t P = TMath::Sqrt(track->GetPx() * track->GetPx() +
                           track->GetPy() * track->GetPy() +
                           track->GetPz() * track->GetPz());
  Double_t s_pion = GetDedxWidthValue(P, SigmaPion());
  Double_t s_kaon = GetDedxWidthValue(P, SigmaKaon());
  Double_t s_proton = GetDedxWidthValue(P, SigmaProton());
  Double_t s_electron = GetDedxWidthValue(P, SigmaElectron());
  Double_t pion = parPiBB->Eval(P);
  Double_t kaon = parKaBB->Eval(P);
  Double_t proton = parPrBB->Eval(P);
  Double_t electro = parElBB->Eval(P);
  s_pion *= pion;
  s_kaon *= kaon;
  s_proton *= proton;
  s_electron *= electro;

  Double_t dEdX = tpc->GetDeDx();
  // std::cout<<"SIGMA\t"<<P<<"\t"<<dEdX<<"\t"<<pion<<"\t"<<s_pion<<std::endl;
  s_pion = (dEdX - pion) / s_pion;
  s_kaon = (dEdX - kaon) / s_kaon;
  s_proton = (dEdX - proton) / s_proton;
  s_electron = (dEdX - electro) / s_electron;
  // std::cout<<s_pion<<std::endl;
  tpc->SetSigma(s_pion, s_kaon, s_proton, s_electron);
  SetValue(s_pion, SigmaPion());
  SetValue(s_kaon, SigmaKaon());
  SetValue(s_proton, SigmaProton());
  SetValue(s_electron, SigmaElectron());
  return Validate();
}

Double_t MpdNSigmaCut::GetDedxWidthValue(Double_t p, Int_t specie) {
  Double_t WidthValue = 0.;
  if (specie == SigmaPion()) {
    if (p < piSigmaLowP->GetXmax()) WidthValue = piSigmaLowP->Eval(p);
    if ((p >= piSigmaMidP->GetXmin()) && (p < piSigmaMidP->GetXmax()))
      WidthValue = piSigmaMidP->Eval(p);
    if ((p >= piSigmaHighP->GetXmin()) && (p <= piSigmaHighP->GetXmax()))
      WidthValue = piSigmaHighP->Eval(p);
  }
  if (specie == SigmaKaon()) {
    if (p < kaSigmaLowP->GetXmax()) WidthValue = kaSigmaLowP->Eval(p);
    if ((p >= kaSigmaHighP->GetXmin()) && (p <= kaSigmaHighP->GetXmax()))
      WidthValue = kaSigmaHighP->Eval(p);
  }
  if (specie == SigmaProton()) {
    if (p < prSigmaLowP->GetXmax()) WidthValue = prSigmaLowP->Eval(p);
    if ((p >= prSigmaHighP->GetXmin()) && (p <= prSigmaHighP->GetXmax()))
      WidthValue = prSigmaHighP->Eval(p);
  }
  if (specie == SigmaElectron()) {
    if (p < elSigmaLowP->GetXmax()) WidthValue = elSigmaLowP->Eval(p);
    if ((p >= elSigmaMidP->GetXmin()) && (p <= elSigmaMidP->GetXmax()))
      WidthValue = elSigmaMidP->Eval(p);
    if ((p >= elSigmaHighP->GetXmin()) && (p <= elSigmaHighP->GetXmax()))
      WidthValue = elSigmaHighP->Eval(p);
  }
  return WidthValue;
}

MpdNSigmaCut::~MpdNSigmaCut() {
  // TODO Auto-generated destructor stub
}
