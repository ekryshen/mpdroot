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

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TString.h>
#include <TVector3.h>
#include <vector>

#include "MpdMiniBTofPidTraits.h"

class MpdEvent;
class MpdMiniEvent;
class MpdMiniTrack;

class MpdV0Monitor : public TObject {
private:
   MpdEvent *            fMpdEvent;    //!
   MpdMiniEvent *        fMiniEvent;   //!
   TClonesArray *        fMiniTracks;  //!
   TClonesArray *        fMiniTofData; //!
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
   /**
    * makes 1 dimensional histogram, should be called in Init
    * @param name name of histogram
    * @param xlabel label of X-axis
    * @param histoId - position of histogram in array with 1D histograms
    */
   void MakeHistogram1d(TString name, TString xlabel, Int_t histoId);
   /**
    * makes 2D histogram, should be called in Init
    * @param name name of histogram
    * @param xlabel label of X-axis
    * @param ylabel label of Y-axis
    * @param histoId - position of histogram in array of 2D histograms
    */
   void MakeHistogram2d(TString name, TString xlabel, TString ylabel, Int_t histoId);
   /**
    * @return MpdEvent, works only for DST
    */
   MpdEvent *    GetEvent() const { return fMpdEvent; };
   MpdMiniEvent *GetMiniEvent() const { return fMiniEvent; };
   MpdMiniTrack *GetMiniTracks(Int_t index) const { return (MpdMiniTrack *)fMiniTracks->UncheckedAt(index); };

   MpdMiniBTofPidTraits *GetTofPidTraits(Int_t index) const
   {
      return (MpdMiniBTofPidTraits *)fMiniTofData->UncheckedAt(index);
   }

public:
   MpdV0Monitor();
   MpdV0Monitor(Int_t d1, Int_t d2);
   MpdV0Monitor(const MpdV0Monitor &other);
   virtual void Init() = 0;
   virtual void SetMiniEventData(MpdMiniEvent *ev, TClonesArray *tracks, TClonesArray *tof)
   {
      fMiniEvent   = ev;
      fMiniTracks  = tracks;
      fMiniTofData = tof;
   };
   virtual void       SetEventData(MpdEvent *ev) { fMpdEvent = ev; }
   std::vector<TH1 *> GetHistograms(TString name) const;
   virtual ~MpdV0Monitor();
   ClassDef(MpdV0Monitor, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MONITORS_MPDV0MONITOR_H_ */
