//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfdPoint
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>

#include <TClonesArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TRandom3.h>

#include "FairRootManager.h"
#include "FairDetector.h"
#include "FairBaseParSet.h"
#include "FairGeoNode.h"
#include "FairRunAna.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "MpdMCTrack.h"
#include "MpdFfdGeoPar.h"
#include "MpdFfdHit.h"
#include "MpdFfdPoint.h"

#include "MpdFfdHitProducer.h"
using namespace std;

ClassImp(MpdFfdHitProducer);

TGraph *MpdFfdHitProducer::gPMTeff_pdf = nullptr;
//------------------------------------------------------------------------------------------------------------------------
MpdFfdHitProducer::MpdFfdHitProducer(const char *name, bool useMCdata, bool timeWindow, Int_t verbose, const char *flnm)
   : FairTask(name, verbose), fUseMCData(useMCdata), fUseTimeWindow(timeWindow), fFlnm(flnm)
{
   fDoTest = (flnm != nullptr);

   if (fDoTest) {
      Add(hSuids = new TH1D("FFD_suids", "suid;entries", 200, -0.5, 199.5));

      Add(hOpYield = new TH2D("FFD_opYield", "Yield of op(per quartz 1 cm) vs #beta;#beta;N op", 1000, 0.1, 1.01, 1000,
                              0.5, 10000.5));
      Add(hOpYieldPion = new TH2D("FFD_opYieldPion", "Yield of op(per quartz 1 cm) vs #beta;#beta;N op", 1000, 0.1,
                                  1.01, 1000, 0.5, 10000.5));
      Add(hOpYieldProton = new TH2D("FFD_opYieldProton", "Yield of op(per quartz 1 cm) vs #beta;#beta;N op", 1000, 0.1,
                                    1.01, 1000, 0.5, 10000.5));
      Add(hOpYieldElectron = new TH2D("FFD_opYieldElectron", "Yield of op(per quartz 1 cm) vs #beta;#beta;N op", 1000,
                                      0.1, 1.01, 1000, 0.5, 10000.5));

      Add(hPeYield = new TH2D("FFD_peYield", "Yield of op(per quartz 1 cm) vs #beta;#beta;N pe", 1000, 0.1, 1.01, 1000,
                              0.5, 1000.5));
      Add(hPeYieldPion = new TH2D("FFD_peYieldPion", "Yield of op(per quartz 1 cm) vs #beta;#beta;N pe", 1000, 0.1,
                                  1.01, 1000, 0.5, 1000.5));
      Add(hPeYieldProton = new TH2D("FFD_peYieldProton", "Yield of op(per quartz 1 cm) vs #beta;#beta;N pe", 1000, 0.1,
                                    1.01, 1000, 0.5, 1000.5));
      Add(hPeYieldElectron = new TH2D("FFD_peYieldElectron", "Yield of op(per quartz 1 cm) vs #beta;#beta;N pe", 1000,
                                      0.1, 1.01, 1000, 0.5, 1000.5));

      Add(hOpTrash =
             new TH2D("FFD_opTrashRatio", "trash op ratio;N op; trash op, %", 1000, 0.5, 1000000.5, 1000, 0., 50.));
      Add(hOccup = new TH2D("FFD_Occup", ";occupancy;suid", 100, -0.5, 99.5, 200, -0.5, 199.5));
      Add(hChannelsEW = new TH2D("FFD_ChannelsEW", ";channels, east;channels, west", 100, -0.5, 99.5, 100, -0.5, 99.5));

      Add(hTimeWindows = new TH2D("FFD_TimeWindows", ";t_{mean}, ns;#Deltat, ns", 1000, 0., 10., 1000, 0., 1.));
      Add(hTimeSize = new TH2D("FFD_TimeSize", ";t_{mean}, ns;n pe", 1000, 0., 100., 1000, 0.5, 1000.5));

      Add(hXY = new TH2D("FFD_XY", "; X, cm; Y, cm", 1000, -20., 20., 1000, -20., 20.));
      Add(hXmap = new TH2D("FFD_Xmap", "; X, cm; suid", 1000, -20., 20., 200, -0.5, 199.5));
      Add(hYmap = new TH2D("FFD_Ymap", "; Y, cm; suid", 1000, -20., 20., 200, -0.5, 199.5));
      Add(hZmap = new TH2D("FFD_Zmap", "; Z, cm; suid", 1000, -150., 150., 200, -0.5, 199.5));

      Add(hXcenter = new TH2D("FFD_Xcenter", "; X, cm; suid", 1000, -20., 20., 200, -0.5, 199.5));
      Add(hYcenter = new TH2D("FFD_Ycenter", "; Y, cm; suid", 1000, -20., 20., 200, -0.5, 199.5));
      Add(hCenter = new TH2D("FFD_XYcenterShift", "; X, cm; Y, cm", 1000, -50., 50., 1000, -50., 50.));

      Add(hPMTeff = new TEfficiency("FFD_PMTeff", ";op energy, eV;PMT efficiency", 100, 1., 9.));
   }
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdHitProducer::~MpdFfdHitProducer()
{
   fList.Delete();
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdFfdHitProducer::Init()
{
   if (fUseMCData) {
      aMcPoints = (TClonesArray *)FairRootManager::Instance()->GetObject("FFDPoint");
      aMcTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
      assert(aMcPoints);
      assert(aMcTracks);
   } else {
      aExpDigits = (TClonesArray *)FairRootManager::Instance()->GetObject("??????"); // FIXME: NOW unknown name
      assert(aExpDigits);
   }

   // Create and register output array
   aFfdHits = new TClonesArray("MpdFfdHit");
   FairRootManager::Instance()->Register("FfdHit", "Ffd", aFfdHits, kTRUE);

   LOG(info) << "[MpdFfdHitProducer::Init] Initialization finished succesfully.";

   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdHitProducer::Exec(Option_t *opt)
{
   // Cleanup containers
   aFfdHits->Clear();
   mmOccup.clear();
   if (fUseTimeWindow) {
      for (size_t i = 0; i < 160; i++) fData[i].clear();
   }

   auto CenterTest = [this](MpdFfdPoint *point, Int_t suid, MpdFfdHit *hit, const TVector3 &padCenter) {
      mmOccup.insert(make_pair(suid, hit));

      TVector3 pos;
      point->Position(pos);
      double padCenterX = padCenter.X(), padCenterY = padCenter.Y();
      hCenter->Fill(padCenterX - pos.X(), padCenterY - pos.Y());
      hXcenter->Fill(padCenterX, suid);
      hYcenter->Fill(padCenterY, suid);
   };

   if (fUseMCData) // input data from MC MpdFfdPoint
   {
      size_t nOpTotal = 0, nOpTrash = 0;

      for (Int_t pointIndex = 0, nFfdPoint = aMcPoints->GetEntriesFast(); pointIndex < nFfdPoint;
           pointIndex++) // cycle by FFD points
      {
         auto   pPoint = (MpdFfdPoint *)aMcPoints->UncheckedAt(pointIndex);
         Int_t  suid   = pPoint->GetDetectorID(); // [1,160]
         Int_t  tid    = pPoint->GetTrackID();
         double beta = pPoint->GetBeta(), nOp = pPoint->GetEntries();
         nOpTotal += nOp;

         if (!pPoint->IsClosed()) // invalid parent track parameters, skip MpdFfdPoint data.
         {
            nOpTrash += nOp;
            continue;
         }

         if (fDoTest) {
            hOpYield->Fill(beta, nOp);
            hSuids->Fill(suid);

            double x = pPoint->GetX(), y = pPoint->GetY();
            hXY->Fill(x, y);
            hXmap->Fill(x, suid);
            hYmap->Fill(y, suid);
            hZmap->Fill(pPoint->GetZ(), suid);
         }

         if (fVerbose > 2) pPoint->Print("");

         size_t nOp2_8 = 0, nPe = 0;
         if (pPoint->GetMode() == MpdFfdPoint::kCherenkovPhoton) {
            for (const auto &entry : pPoint->GetData()) // cycle by op
            {
               double energy = entry.first;
               if (2. < energy && energy < 8.) nOp2_8++;
               ; //  	<<<----------- [2,8] eV  cut

               bool pass = IsPeCreated(energy); // simulate PMT efficiency
               if (pass) {
                  nPe++;
                  if (fUseTimeWindow) AddEntry(entry.second, suid, pPoint, pointIndex);
               }

               hPMTeff->Fill(pass, energy);

            }   // cycle by op
         } else // kPhotoElectron mode, saved only pe
         {
            nPe = pPoint->GetEntries();

            if (fUseTimeWindow) {
               for (const auto &entry : pPoint->GetData())
                  AddEntry(entry.second, suid, pPoint, pointIndex); // cycle by pe
            }
         }

         if (!fUseTimeWindow && nPe >= fNpeThresh) // create FFD hit at fast mode
         {
            TVector3 padCenter = GetPadCenter(suid);
            auto     hit       = AddHit(pPoint, padCenter, TVector3(fErrXY, fErrXY, fErrZ), pointIndex, nPe, 0);

            if (fDoTest) CenterTest(pPoint, suid, hit, padCenter);
         }

         if (fDoTest) {
            hPeYield->Fill(beta, nPe);

            auto  mcTrack = (const MpdMCTrack *)aMcTracks->UncheckedAt(tid);
            Int_t pid     = mcTrack->GetPdgCode();

            if (2212 == pid || -2212 == pid) {
               if (nOp2_8) hOpYieldProton->Fill(beta, nOp2_8);
               hPeYieldProton->Fill(beta, nPe);
            } else if (211 == pid || -211 == pid) {
               if (nOp2_8) hOpYieldPion->Fill(beta, nOp2_8);
               hPeYieldPion->Fill(beta, nPe);
            } else if (11 == pid || -11 == pid) {
               if (nOp2_8) hOpYieldElectron->Fill(beta, nOp2_8);
               hPeYieldElectron->Fill(beta, nPe);
            }

            assert(!pPoint->GetData().empty());
            double Tmin = 1.E10, Tmax = -1., Tmean = 0.;
            for (const auto &entry : pPoint->GetData()) // cycle by pe
            {
               double time = entry.second;
               if (Tmin > time) Tmin = time;
               if (time > Tmax) Tmax = time;
               Tmean += time;
            }

            Tmean /= pPoint->GetData().size();
            hTimeWindows->Fill(Tmean, Tmax - Tmin);
            hTimeSize->Fill(Tmean, nPe);
         }

      } // cycle by FFD points

      if (fUseTimeWindow) {
         Int_t nWindowSize = fTimeWindow * fRoIBins / fRoISize;
         nWindowSize       = ((UInt_t)nWindowSize < fRoIBins) ? nWindowSize : fRoIBins - 1;

         for (size_t suid = 0; suid < 160; suid++) {
            auto data = fData[suid];

            for (auto iter = data.begin(), itEnd = data.end(); iter != itEnd;
                 iter++) // cycle by non-zero deposit time bins
            {
               Tintegrals mIntegrals; // key = MpdFfdPoint index, v

               double timeWindowEnd = GetTime(iter->first) + fTimeWindow;
               size_t bin;
               auto   itWindowEnd = GetTimeBin(timeWindowEnd, bin) ? data.upper_bound(bin) : itEnd;

               for (auto it = iter; it != itWindowEnd; it++) // cycle by window bins
               {
                  auto entry = mIntegrals.find(it->second.pointIndex);
                  if (entry != mIntegrals.end())
                     entry->second.nPe += it->second.nPe;
                  else
                     mIntegrals.insert(make_pair(it->second.pointIndex, it->second));
               }

               size_t Integral = 0, maxDeposit = 0, maxIndex = -1;
               for (const auto &entry : mIntegrals) // cycle by MpdFfdPoint deposit into this time window
               {
                  size_t nPe = entry.second.nPe;
                  Integral += nPe;

                  if (nPe > maxDeposit) {
                     maxDeposit = nPe;
                     maxIndex   = entry.first;
                  }
               }

               if (Integral >= fNpeThresh) // create FFD hit at TimeWindow mode
               {
                  TVector3 padCenter = GetPadCenter(suid + 1); // [0,159] -> [1,160]
                  auto     pPoint    = (MpdFfdPoint *)aMcPoints->UncheckedAt(maxIndex);
                  auto     hit       = AddHit(pPoint, padCenter, TVector3(fErrXY, fErrXY, fErrZ), maxIndex, maxDeposit,
                                              mIntegrals.size(), mIntegrals);

                  if (fDoTest) CenterTest(pPoint, suid, hit, padCenter);

                  iter = itWindowEnd; // integrate pe by only ONCE time
                  if (iter == itEnd) break;
               }

            } // cycle by time bins

         } // cycle by suid

      } // fUseTimeWindow

      if (fDoTest) {
         if (nOpTotal) hOpTrash->Fill(nOpTotal, (100. * nOpTrash) / nOpTotal); // [%]

         size_t nE = 0, nW = 0;
         for (auto iter = mmOccup.begin(); iter != mmOccup.end();
              iter      = mmOccup.upper_bound(iter->first)) // cycle by unique channel hits
         {
            size_t suid = iter->first;
            hOccup->Fill(mmOccup.count(suid), suid);
            if (suid <= 80)
               nE++;
            else
               nW++;
         }

         hChannelsEW->Fill(nE, nW);
      }
   } else // input data from experimental MpdFfdDigit
   {
      // FIXME: now not realized
      // AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
      // Int_t flag)
      assert(false);
   }

   LOG(debug1) << "[MpdFfdHitProducer::Exec] FFD hits = " << aFfdHits->GetEntriesFast();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdHitProducer::AddEntry(double time, size_t suid, const MpdFfdPoint *ptr, size_t index)
{
   assert(suid <= 160);

   suid--; // [1,160] -> [0,159]
   size_t timeBin;
   if (GetTimeBin(time, timeBin)) {
      bool found = false;
      auto range = fData[suid].equal_range(timeBin);
      for (auto it = range.first; it != range.second; it++) {
         if (it->second.pointIndex == index) {
            it->second.nPe++;
            found = true;
            break;
         }
      }

      if (!found)
         fData[suid].insert(std::make_pair(
            timeBin, PointData(index, ptr, 1))); // insert first time deposit for pair(timebin, MpdFfdPoint)
   }
}
//------------------------------------------------------------------------------------------------------------------------
TVector3 MpdFfdHitProducer::GetPadCenter(Int_t suid) // [1,161]
{
   assert(1 <= suid && suid <= 161);

   static double *X = nullptr, *Y = nullptr;

   auto calcCenters = [this](double x, double y, size_t index, const char *comment = nullptr,
                             std::ostream &os = std::cout) {
      TVector3  vect(x, y, 0.);
      TRotation rot;
      rot.RotateZ(45. * TMath::DegToRad());
      vect = rot * vect;

      X[index] = vect.X();
      Y[index] = vect.Y();

      if (fDoTest && fVerbose > 1) {
         if (comment) os << comment;
         os << "\n suid=" << (index + 1) << " (" << X[index] << "," << Y[index] << ")";
      }
   };

   if (X == nullptr) // init by first time.
   {
      const double posX[] = {7.9, 14.5, -7.9, -14.5, 0.,  0.,   0.,   0.,    6.6,  -6.6,
                             6.6, -6.6, 6.6,  -6.6,  6.6, -6.6, 13.2, -13.2, 13.2, -13.2}; //[cm]
      const double posY[] = {0.,   0.,   0.,   0.,   7.9,   14.5,  -7.9, -14.5, 6.6,  6.6,
                             -6.6, -6.6, 13.2, 13.2, -13.2, -13.2, 6.6,  6.6,   -6.6, -6.6}; //[cm]

      X = new double[160];
      Y = new double[160];

      size_t index = -1;                      // [0,159]
      for (size_t side = 0; side < 2; side++) // cycle by sides
      {
         for (size_t k = 0; k < 20; k++) // cycle by detectors
         {
            calcCenters(posX[k] + 1.4, posY[k] + 1.4, ++index, "\n-------------------------");
            calcCenters(posX[k] - 1.4, posY[k] + 1.4, ++index);
            calcCenters(posX[k] - 1.4, posY[k] - 1.4, ++index);
            calcCenters(posX[k] + 1.4, posY[k] - 1.4, ++index);
         }
      }
   } // init by first time

   size_t index = suid - 1; // [1,161] -> [0,159]

   double Z = (index < 80) ? 140.7 + 2.55 : -140.7 - 2.55;

   return TVector3(X[index], Y[index], Z);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdHitProducer::Finish()
{
   if (fFlnm.Length() > 0) {
      LOG(debug2) << "[MpdFfdHitProducer::Finish] Update  " << fFlnm.Data() << " file. ";
      auto  ptr = gFile;
      TFile file(fFlnm.Data(), "RECREATE");
      fList.Write();
      file.Close();
      gFile = ptr;
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdHitProducer::SetParContainers()
{
   // Get run and runtime database
   FairRunAna *run = FairRunAna::Instance();
   if (!run) Fatal("SetParContainers", "No analysis run");

   FairRuntimeDb *db = run->GetRuntimeDb();
   if (!db) Fatal("SetParContainers", "No runtime database");

   // Get √ geometry parameter container
   db->getContainer("MpdFfdGeoPar");
}
//------------------------------------------------------------------------------------------------------------------------
bool MpdFfdHitProducer::IsPeCreated(double energy) // [eV]
{
   if (gPMTeff_pdf == nullptr) // init by first fime
   {
      const size_t n         = 11;
      const double length[n] = {155., 170., 200., 250., 300., 350., 400., 450., 500., 550., 600.}; // op wavelength [nm]
      const double eff[n]    = {0.20, 0.20, 0.18, 0.16,  0.18, 0.21,
                             0.22, 0.16, 0.11, 0.035, 0.015}; // pe creating efficiency

      gPMTeff_pdf = new TGraph(n, length, eff);
   }

   if (energy < 2. || energy > 8.) return false; //  [2,8] eV fast cut

   double wavelength = 1239.8 / energy;
   double efficiency = gPMTeff_pdf->Eval(wavelength);

   return (gRandom->Rndm() < efficiency);
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdHit *MpdFfdHitProducer::AddHit(const MpdFfdPoint *point, const TVector3 &pos, const TVector3 &dpos,
                                     Int_t refIndex, size_t npe, Int_t flag)
{
   MpdFfdHit *pHit = new ((*aFfdHits)[aFfdHits->GetEntriesFast()])
      MpdFfdHit(point->GetDetectorID(), pos, dpos, refIndex, point->GetTime(), npe, flag);

   pHit->AddLink(FairLink(1, point->GetTrackID())); // LS: key value = 1  for mc track index

   return pHit;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdHit *MpdFfdHitProducer::AddHit(const MpdFfdPoint *point, const TVector3 &pos, const TVector3 &dpos,
                                     Int_t refIndex, size_t npe, Int_t flag, const Tintegrals &map)
{
   MpdFfdHit *pHit = new ((*aFfdHits)[aFfdHits->GetEntriesFast()])
      MpdFfdHit(point->GetDetectorID(), pos, dpos, refIndex, point->GetTime(), npe, flag);

   for (const auto &entry : map) {
      pHit->AddLink(FairLink(1, entry.second.pPoint->GetTrackID())); // LS: key value = 1  for mc track index
      pHit->AddLink(FairLink(2, entry.second.pointIndex));           // LS: key value = 2  for mc point index
   }
   return pHit;
}
//------------------------------------------------------------------------------------------------------------------------
