//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfdHitProducer
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __HH_MPDFFDHITPRODUCER_H
#define __HH_MPDFFDHITPRODUCER_H

#include <map>

#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TEfficiency.h>

#include "FairTask.h"
//------------------------------------------------------------------------------------------------------------------------
class TClonesArray;
class TGraph;
class MpdFfdPoint;
class MpdFfdHit;

//------------------------------------------------------------------------------------------------------------------------
class MpdFfdHitProducer : public FairTask {
public:
   struct PointData {
      size_t             pointIndex;
      const MpdFfdPoint *pPoint;
      size_t             nPe;

      PointData(size_t index, const MpdFfdPoint *ptr, size_t npe) : pointIndex(index), pPoint(ptr), nPe(npe){};
      PointData() : pointIndex(-1), pPoint(nullptr), nPe(0){};
   };
   typedef std::map<size_t, PointData> Tintegrals;

private:
   TClonesArray *aMcPoints  = nullptr; //! <--- MC input
   TClonesArray *aMcTracks  = nullptr; //! <--- MC input
   TClonesArray *aExpDigits = nullptr; //! <--- Exp input
   TClonesArray *aFfdHits   = nullptr; //! ---> output

   Bool_t   fUseMCData = true, fUseTimeWindow = false;
   Bool_t   fDoTest    = false;
   size_t   fNpeThresh = 70.;                                 // default value = 80 photo electrons
   Double_t fErrXY = 4. / sqrt(12.), fErrZ = 1.5 / sqrt(12.); // 4x4x1.5 [cm]
   TList    fList;
   TString  fFlnm;

   std::multimap<size_t, const MpdFfdHit *> mmOccup;          //!
   Double_t                                 fTimeWindow = 1.; // [ns] channel time window for pe integration
   static const size_t                      fRoIBins    = 1000;
   const double                             fRoISize    = 10; // Region of Interest = [0,10] ns

   // 160 channels(quartz SV), key = timebin
   std::multimap<size_t, PointData> fData[160];

   TH1D        *hSuids   = nullptr;
   TH2D        *hOpYield = nullptr, *hOpYieldPion = nullptr, *hOpYieldProton = nullptr, *hOpYieldElectron = nullptr;
   TH2D        *hPeYield = nullptr, *hPeYieldPion = nullptr, *hPeYieldProton = nullptr, *hPeYieldElectron = nullptr;
   TH2D        *hOpTrash = nullptr, *hXY = nullptr, *hXmap = nullptr, *hYmap = nullptr, *hZmap = nullptr;
   TH2D        *hXcenter = nullptr, *hYcenter = nullptr, *hCenter = nullptr, *hOccup = nullptr, *hChannelsEW = nullptr;
   TH2D        *hTimeWindows = nullptr, *hTimeSize = nullptr;
   TEfficiency *hPMTeff = nullptr;
   static TGraph *gPMTeff_pdf; // PDF for PMT yield

   virtual void SetParContainers();

   void Add(TEfficiency *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }
   void Add(TH1 *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }
   MpdFfdHit *AddHit(const MpdFfdPoint *point, const TVector3 &pos, const TVector3 &dpos, Int_t refIndex, size_t npe,
                     Int_t flag);
   MpdFfdHit *AddHit(const MpdFfdPoint *point, const TVector3 &pos, const TVector3 &dpos, Int_t refIndex, size_t npe,
                     Int_t flag, const Tintegrals &);

   TVector3 GetPadCenter(Int_t suid); // [1,161]
   void     AddEntry(double time, size_t suid, const MpdFfdPoint *ptr, size_t index);

   inline bool GetTimeBin(double time, size_t &bin) // [ns]
   {
      time -= 1.8;                     // shift to 0
      if (0 < time && time < fRoISize) // time inside RoI
      {
         bin = (size_t)(time * fRoIBins / fRoISize); // bin size = fRoISize/fRoIBins
         return true;
      }
      return false;
   }

   inline double GetTime(size_t bin) // [ns]
   {
      assert(0 <= bin && bin < fRoIBins);

      double time = bin * fRoISize / fRoIBins;
      time += 1.8; // shift from 0
      return time;
   }

public:
   MpdFfdHitProducer(const char *name = "FFD hit Producer", bool useMCdata = true, bool useTimeWindow = false,
                     Int_t verbose = 1, const char *QAflnm = nullptr);
   virtual ~MpdFfdHitProducer();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();

   void SetErrors(double XYerr, double Zerr)
   {
      fErrXY = XYerr;
      fErrZ  = Zerr;
   };
   void        SetPeThresh(size_t n) { fNpeThresh = n; };
   void        SetTimeWindow(Double_t v) { fTimeWindow = v; };
   static bool IsPeCreated(double energy); // [eV]

   ClassDef(MpdFfdHitProducer, 3);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __HH_MPDFFDHITPRODUCER_H
