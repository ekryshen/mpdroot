/*
 * MpdV0Monitor.h
 *
 *  Created on: 11 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0MONITOR_H_
#define MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0MONITOR_H_

#include <RtypesCore.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TVector3.h>
#include <utility>

class MpdEvent;
class MpdMiniEvent;

class MpdV0Monitor : public TObject {
private:
   MpdEvent *            fMpdEvent;  //!
   MpdMiniEvent *        fMiniEvent; //!
   std::vector<TH1 *>    fHistograms1dPassed;
   std::vector<TH1 *>    fHistograms1dFailed;
   std::vector<TH1 *>    fHistograms2dPassed;
   std::vector<TH1 *>    fHistograms2dFailed;
   std::vector<TVector3> fHistogram1dXaxis;
   std::vector<TVector3> fHistogram2dXaxis;
   std::vector<TVector3> fHistogram2dYaxis;

protected:
   Bool_t fInit;
   /**
    * fill 1 dimensiona histogram
    * histoId is a position of histogram in vector of 1d histogram
    * passed specify if this entry should go to passed or failed histograms
    */
   void Fill1D(Int_t histoId, Double_t val, Bool_t passed)
   {
      passed ? fHistograms1dPassed[histoId]->Fill(val) : fHistograms1dFailed[histoId]->Fill(val);
   };
   void Fill2D(Int_t histoId, Double_t valX, Double_t valY, Bool_t passed)
   {
      passed ? fHistograms2dPassed[histoId]->Fill(valX, valY) : fHistograms2dFailed[histoId]->Fill(valX, valY);
   };
   void SetXaxis1d(Int_t histoID, Int_t bins, Double_t min, Double_t max)
   {
      fHistogram1dXaxis[histoID].SetXYZ(bins, min, max);
   };
   void SetXaxis2d(Int_t histoID, Int_t bins, Double_t min, Double_t max)
   {
      fHistogram2dXaxis[histoID].SetXYZ(bins, min, max);
   };
   void SetYaxis2d(Int_t histoID, Int_t bins, Double_t min, Double_t max)
   {
      fHistogram2dYaxis[histoID].SetXYZ(bins, min, max);
   };
   void          MakeHistogram1d(TString name, TString xlabel, Int_t histoId);
   void          MakeHistogram2d(TString name, TString xlabel, TString ylabel, Int_t histoId);
   MpdEvent *    GetEvent() const { return fMpdEvent; };
   MpdMiniEvent *GetMiniEvent() const { return fMiniEvent; };
   void          PushHistogram(std::vector<TH1 *> vec, TH1 *histo, TString name) const;

public:
   MpdV0Monitor();
   MpdV0Monitor(Int_t d1, Int_t d2);
   MpdV0Monitor(const MpdV0Monitor &other);
   virtual void       Init() = 0;
   virtual void       SetMiniEvent(MpdMiniEvent *ev) { fMiniEvent = ev; };
   virtual void       SetEvent(MpdEvent *ev) { fMpdEvent = ev; }
   std::vector<TH1 *> GetHistograms(TString name) const;
   virtual ~MpdV0Monitor();
   ClassDef(MpdV0Monitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0MONITOR_H_ */
