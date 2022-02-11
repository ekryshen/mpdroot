/*
 * MpdV0Monitor.cxx
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Monitor.h"
#include <iostream>

MpdV0Monitor::MpdV0Monitor() : fMpdEvent(nullptr), fMiniEvent(nullptr), fInit(kFALSE) {}

void MpdV0Monitor::MakeHistogram1d(TString name, TString xlabel, Int_t histoId)
{
   TVector3 axisConf = fHistogram1dXaxis[histoId];
   TString  namePass = name + "_Passed";
   TString  nameFail = name + "_Failed";
   fHistograms1dPassed[histoId] =
      new TH1D(namePass, Form("%s;%s;", namePass.Data(), xlabel.Data()), axisConf.X(), axisConf.Y(), axisConf.Z());
   fHistograms1dFailed[histoId] =
      new TH1D(nameFail, Form("%s;%s;", nameFail.Data(), xlabel.Data()), axisConf.X(), axisConf.Y(), axisConf.Z());
}

void MpdV0Monitor::MakeHistogram2d(TString name, TString xlabel, TString ylabel, Int_t histoId)
{
   TVector3 axisConfX = fHistogram2dXaxis[histoId];
   TVector3 axisConfY = fHistogram2dYaxis[histoId];
   TString  namePass  = name + "_Passed";
   TString  nameFail  = name + "_Failed";
   fHistograms2dPassed[histoId] =
      new TH2D(namePass, Form("%s;%s;%s", namePass.Data(), xlabel.Data(), ylabel.Data()), axisConfX.X(), axisConfX.Y(),
               axisConfX.Z(), axisConfY.X(), axisConfY.Y(), axisConfY.Z());
   fHistograms2dFailed[histoId] =
      new TH2D(nameFail, Form("%s;%s;%s", nameFail.Data(), xlabel.Data(), ylabel.Data()), axisConfX.X(), axisConfX.Y(),
               axisConfX.Z(), axisConfY.X(), axisConfY.Y(), axisConfY.Z());
}

MpdV0Monitor::MpdV0Monitor(Int_t d1, Int_t d2) : fMpdEvent(nullptr), fMiniEvent(nullptr), fInit(kFALSE)
{
   fHistogram1dXaxis.resize(d1);
   fHistograms1dPassed.resize(d1);
   fHistograms1dFailed.resize(d1);

   fHistograms2dPassed.resize(d2);
   fHistograms2dFailed.resize(d2);

   fHistogram2dXaxis.resize(d2);
   fHistogram2dYaxis.resize(d2);
   for (int i = 0; i < d1; i++) {
      fHistograms1dPassed[i] = nullptr;
      fHistograms1dFailed[i] = nullptr;
   }
   for (int i = 0; i < d2; i++) {
      fHistograms2dPassed[i] = nullptr;
      fHistograms2dFailed[i] = nullptr;
   }
}

MpdV0Monitor::MpdV0Monitor(const MpdV0Monitor &other) : fInit(other.fInit), fMpdEvent(nullptr), fMiniEvent(nullptr)
{

   auto copyStuff = [](std::vector<TH1 *> &this_vec, const std::vector<TH1 *> &other_vec) {
      for (int i = 0; i < other_vec.size(); i++) {
         auto hist = other_vec.at(i);
         if (hist == nullptr) {
            this_vec.push_back(nullptr);
         } else {
            this_vec.push_back((TH1 *)hist->Clone());
         }
      }
   };
   copyStuff(fHistograms1dPassed, other.fHistograms1dPassed);
   copyStuff(fHistograms1dFailed, other.fHistograms1dFailed);
   copyStuff(fHistograms2dPassed, other.fHistograms2dPassed);
   copyStuff(fHistograms2dFailed, other.fHistograms2dFailed);
   fHistogram1dXaxis = other.fHistogram1dXaxis;
   fHistogram2dXaxis = other.fHistogram2dXaxis;
   fHistogram2dYaxis = other.fHistogram2dYaxis;
}

std::vector<TH1 *> MpdV0Monitor::GetHistograms(TString name) const
{
   std::vector<TH1 *> res;
   auto               copyStuff = [](std::vector<TH1 *> &vec, const std::vector<TH1 *> &in, TString prefix) {
      for (int i = 0; i < in.size(); i++) {
         const TH1 *histo   = in.at(i);
         TString    nameNew = histo->GetName();
         nameNew            = prefix + "_" + nameNew;
         TH1 *copy          = (TH1 *)histo->Clone(nameNew);
         vec.push_back(copy);
      }
   };
   copyStuff(res, fHistograms1dPassed, name);
   copyStuff(res, fHistograms1dFailed, name);
   copyStuff(res, fHistograms2dPassed, name);
   copyStuff(res, fHistograms2dFailed, name);

   return res;
}

MpdV0Monitor::~MpdV0Monitor()
{
   if (fInit) {
      for (int i = 0; i < fHistograms1dPassed.size(); i++) {
         delete fHistograms1dPassed[i];
         delete fHistograms1dFailed[i];
      }
      for (int i = 0; i < fHistograms2dPassed.size(); i++) {
         delete fHistograms2dPassed[i];
         delete fHistograms2dFailed[i];
      }
   }
}
