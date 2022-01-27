/*************************************************************************************
 *
 *         Class MpdZdcTstSim
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#include "MpdZdcTstSim.h"
#include "MpdZdcPoint.h"
#include "TpcPoint.h"
#include <TFile.h>

#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairBaseParSet.h"

using std::cout;
using std::endl;

MpdZdcTstSim::MpdZdcTstSim()
// : FairTask ("ZdcTstSim")
{
   fVerbose = 0;
   nevents  = 0;
   fTree    = 0;
   fH3      = 0;
}

MpdZdcTstSim::MpdZdcTstSim(const char *name, const char *title, Int_t verbose) : FairTask(name)
{
   if (verbose) cout << "-Info- MpdZdcTstSim::MpdZdcTstSim begin" << endl;

   fVerbose = verbose;
   nevents  = 0;

   fTree        = 0;
   fTreeSummary = 0;
   fH3          = 0;
   //  CreateMyTree();

   if (fVerbose) cout << "-Info- MpdZdcTstSim::MpdZdcTstSim end" << endl;
}

MpdZdcTstSim::MpdZdcTstSim(const char *name, const char *title, Int_t verbose, Int_t flag)
{
   if (verbose) cout << "-Info- MpdZdcTstSim::MpdZdcTstSim begin" << endl;

   fVerbose = verbose;
   nevents  = 0;

   fTree        = 0;
   fTreeSummary = 0;
   //  CreateMyTree();

   if (fVerbose) cout << "-Info- MpdZdcTstSim::MpdZdcTstSim Init starting" << endl;

   Init();

   if (fVerbose) cout << "-Info- MpdZdcTstSim::MpdZdcTstSim end" << endl;
}

MpdZdcTstSim::~MpdZdcTstSim()
{
   if (fTree) fTree->Write();
   if (fTreeSummary) fTreeSummary->Write();

   if (fH3) fH3->Write();
}

void MpdZdcTstSim::CreateMyTree()
{
   if (!fTree)
      fTree = new TNtuple(
         "zdcpoints", "zdc points",
         "tpdg:tmotherid:tpx:tpy:tpz:tx:ty:tz:tmass:tenergy:tpt:tp:trapidity:tstarttime:ptrackid:peventid:px:py:pz:ppx:"
         "ppy:ppz:ptime:plength:peloss:pmcvolumeid:pcopy:pcopymother:tmotherpdg:tmotherz:tZ:tA",
         1000000);

   if (!fTreeSummary)
      fTreeSummary = new TNtupleD("zdcsummary", "zdc summary",
                                  "elossall:elosszdc1:elosszdc2:tpctracks:b:ekin:ekin1:nucleons1:pions1:protons1:ekin2:"
                                  "nucleons2:pions2:protons2",
                                  //				 "elossall:elosszdc1:elosszdc2:tpctracks:b:ekin:ezdc1[4]:ezdc2[4]",
                                  1000000);

   //    if (fVerbose) {
   //     cout << "-Info- MpdZdcTstSim::CreateMyTree zdcpoints:" << (fTree->GetCurrentFile()->GetName())<< endl;
   //     cout << "-Info- MpdZdcTstSim::CreateMyTree zdcsummary:" << (fTreeSummary->GetCurrentFile()->GetName())<< endl;
   //    }

   if (!fH3) {
      fH3 = new TH3F("H3", "Counts(ELoss,TpcTracks,b)", 40, 0, 40, 120, 0, 1200, 32, 0, 16);
      //     fH3 = new TH3F ("H3","Counts(ELoss,TpcTracks,b)",50,0,10,120,0,1200,40,0,2);
      //    fH3->SetDirectory(0);
   }
}

InitStatus MpdZdcTstSim::Init()
{
   if (fVerbose) cout << "-Info- MpdZdcTstSim::Init begin" << endl;

   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdZdcTstSim::Init: "
           << "RootManager not instantised!" << endl;
      return kERROR;
   }

   fMCTrackArray = (TClonesArray *)ioman->GetObject("MCTrack");
   if (!fMCTrackArray) {
      cout << "-E- MpdZdcTstSim::Init: No MCTrack array!" << endl;
      return kERROR;
   }

   fMCZdcPointArray = (TClonesArray *)ioman->GetObject("ZdcPoint");
   if (!fMCZdcPointArray) {
      cout << "-E- MpdZdcTstSim::Init: No ZdcPoint array!" << endl;
      return kERROR;
   }

   fMCTpcPointArray = (TClonesArray *)ioman->GetObject("TpcPoint");
   if (!fMCTpcPointArray) {
      cout << "-E- MpdZdcTstSim::Init: No TpcPoint array!" << endl;
      //    return kERROR;
   }

   fMCEventHeader = (FairMCEventHeader *)ioman->GetObject("MCEventHeader.");
   if (!fMCEventHeader) {
      cout << "-E- MpdZdcTstSim::Init: No MCEventHeader array!" << endl;
      //    return kERROR;
   }

   if (fVerbose) cout << "-Info- MpdZdcTstSim::Init end" << endl;

   CreateMyTree();
   fTree->SetDirectory(ioman->GetOutFile());
   fTreeSummary->SetDirectory(ioman->GetOutFile());
   if (fH3) fH3->SetDirectory(ioman->GetOutFile());
   ioman->Register("zdcpoints", "ZDC", fTree, kTRUE);
   ioman->Register("zdcsummary", "ZDC", fTreeSummary, kTRUE);
   if (fH3) ioman->Register("H3", "ZDC", fH3, kTRUE);

   if (fVerbose) {
      cout << "-Info- MpdZdcTstSim::Init zdcpoints: " << (fTree->GetCurrentFile()->GetName()) << endl;
      cout << "-Info- MpdZdcTstSim::Init zdcsummary: " << (fTreeSummary->GetCurrentFile()->GetName()) << endl;
   }

   return kSUCCESS;
}

InitStatus MpdZdcTstSim::ReInit()
{
   //   //  SetParContainers();
   //    FairRunAna* ana = FairRunAna::Instance();
   //    FairRuntimeDb* rtdb=ana->GetRuntimeDb();
   //    MpdZdcGeoPar *fPar=(MpdZdcGeoPar*)(rtdb->getContainer("MpdZdcGeoPar"));
   return kSUCCESS;
}

void MpdZdcTstSim::Exec(Option_t *option)
{
   if (fVerbose > 3) cout << "-Info- MpdZdcTstSim::Exec begin" << endl;

   Double_t elossall = 0, elosszdc1 = 0, elosszdc2 = 0, e = 0, ezdc1[4], ezdc2[4];
   Int_t    nPoints, nMCTracks;

   if (fTree) {

      TVector3     position, momentum, vertex;
      MpdZdcPoint *pPoint = NULL;
      MpdMCTrack * pTrack = NULL;
      Int_t        pdg    = 0;

      Float_t fArgs[30];

      for (int iii = 0; iii < 4; iii++) {
         ezdc1[iii] = 0;
         ezdc2[iii] = 0;
      }

      nPoints   = fMCZdcPointArray->GetEntriesFast();
      nMCTracks = fMCTrackArray->GetEntriesFast();

      if (fVerbose) {
         cout << " MpdZdcPoint: " << nPoints << " entries" << endl;
         cout << " MCTrack:     " << nMCTracks << " entries" << endl;
      }

      for (Int_t i = 0; i < nPoints; i++) {

         //       if (fVerbose>2)
         // 	cout << " .";

         pPoint = (MpdZdcPoint *)fMCZdcPointArray->UncheckedAt(i);
         if (!pPoint) {
            if (fVerbose > 2) {
               cout << " -";
            }
            continue;
         }

         if (pPoint->GetTrackID() == -2) {
            if (fVerbose > 2) {
               //	  cout << i << ":" <<(pPoint->GetTrackID()) << "!";
               cout << "MpdZdcTstSim::Exec:  MCTrack not saved: particle number:" << i << " must be too slow !" << endl;
            }
            continue;
         }

         pTrack = (MpdMCTrack *)fMCTrackArray->UncheckedAt(pPoint->GetTrackID());
         if (!pTrack) {
            if (fVerbose > 2) {
               cout << i << ":" << (pPoint->GetTrackID()) << "!";
            }
            continue;
         }

         try {
            // 	  if ((fVerbose>2)&&(nMCTracks==174937)) {
            // 	    cout << i << ":" <<(pPoint->GetTrackID()) << " track addr:" << pTrack<< "! ";
            // 	  }
            pdg      = pTrack->GetPdgCode();
            fArgs[0] = pdg;
         } catch (...) {
            cout << "MpdZdcTstSim::Exec:  Exception 1: particle number:" << i
                 << " track number:" << (pPoint->GetTrackID()) << " track addr:" << pTrack << endl;
         }

         fArgs[1] = pTrack->GetMotherId();

         pTrack->GetMomentum(momentum);
         fArgs[2] = momentum.X();
         fArgs[3] = momentum.Y();
         fArgs[4] = momentum.Z();

         pTrack->GetStartVertex(vertex);
         fArgs[5] = vertex.X();
         fArgs[6] = vertex.Y();
         fArgs[7] = vertex.Z();

         Int_t tZ = -1000, tA = -1000;

         TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle((Int_t)floor(fArgs[0]));
         if (particle) {
            fArgs[8] = particle->Mass();
            tZ       = TMath::Nint(particle->Charge() / 3.);
            if (pdg > 1e9) {
               tA = ((pdg % 10000) / 10);
               if (!tZ) tZ = TMath::Nint((pdg - 1e9) / 10000);
            }
         } else {
            fArgs[8] = 0;
            if (pdg > 1e9) {
               tZ = TMath::Nint((pdg - 1e9) / 10000);
               tA = ((pdg % 10000) / 10);
            }
         }

         fArgs[9] = TMath::Sqrt(fArgs[2] * fArgs[2] + fArgs[3] * fArgs[3] + fArgs[4] * fArgs[4] + fArgs[8] * fArgs[8]);

         fArgs[10] = TMath::Sqrt(fArgs[2] * fArgs[2] + fArgs[3] * fArgs[3]);
         fArgs[11] = TMath::Sqrt(fArgs[2] * fArgs[2] + fArgs[3] * fArgs[3] + fArgs[4] * fArgs[4]);
         fArgs[12] = 0.5 * TMath::Log((fArgs[9] + fArgs[4]) / (fArgs[9] - fArgs[4]));

         fArgs[13] = 1; // not implemented on the level of MpdMCTrack now
         // fArgs[13]=pTrack->GetStartTime();

         fArgs[14] = pPoint->GetTrackID();
         //      fArgs[15]= pPoint->GetEventID();
         fArgs[15] = nevents;
         fArgs[16] = pPoint->GetX();
         fArgs[17] = pPoint->GetY();
         fArgs[18] = pPoint->GetZ();
         fArgs[19] = pPoint->GetPx();
         fArgs[20] = pPoint->GetPy();
         fArgs[21] = pPoint->GetPz();
         fArgs[22] = pPoint->GetTime();
         fArgs[23] = pPoint->GetLength();
         fArgs[24] = pPoint->GetEnergyLoss();
         //       fArgs[25]= pPoint->GetModule();
         //       fArgs[26]= pPoint->GetRow();
         //       fArgs[27]= pPoint->GetCopy();
         fArgs[25] = pPoint->GetDetectorID();
         fArgs[26] = pPoint->GetCopy();
         fArgs[27] = pPoint->GetCopyMother();

         if (fArgs[1] != -1) {
            pTrack = (MpdMCTrack *)fMCTrackArray->UncheckedAt(pTrack->GetMotherId());
            if (pTrack) {
               fArgs[28] = pTrack->GetPdgCode();
               TVector3 vertexTemp;
               pTrack->GetStartVertex(vertexTemp);
               fArgs[29] = vertexTemp.Z();
            } else {
               fArgs[28] = -2000000;
               fArgs[29] = -2000000;
            }
         } else {
            fArgs[28] = -1000000;
            fArgs[29] = -1000000;
         }

         // gertsen COMMENT BECAUSE OF OUT OF RANGE!!!
         // fArgs[30]=tZ;
         // fArgs[31]=tA;

         elossall += pPoint->GetEnergyLoss();

         Double_t ekin = fArgs[9] - fArgs[8];

         //      if ((fArgs[25]<=30)&&(fArgs[1]==-1))

         if (fArgs[1] == -1) { // only primary particles:

            //	if (fArgs[27]==1)                        // zdc1
            if (fArgs[18] > 0) // zdc1
               elosszdc1 += pPoint->GetEnergyLoss();
            else
               elosszdc2 += pPoint->GetEnergyLoss();

            if (fArgs[26] == 1) { // zdc entrances

               e += ekin; // kinetic energy of all primary particles at both zdc entrances

               //	  if (fArgs[27]==1) {               // zdc1
               if (fArgs[18] > 0) { // zdc1
                  ezdc1[0] += ekin; // all
                  switch (pdg) {
                  case 2212:           // protons
                     ezdc1[1] += ekin; // nucleons
                     ezdc1[3] += ekin; // protons
                     break;
                  case 2112:           // neutrons
                     ezdc1[1] += ekin; // nucleons
                     break;
                  case 211: // pions
                  case -211: ezdc1[2] += fArgs[9]; break;
                  default: break;
                  }
               } else { // zdc2
                  ezdc2[0] += ekin;
                  switch (pdg) {
                  case 2212:           // protons
                     ezdc2[1] += ekin; // nucleons
                     ezdc2[3] += ekin; // protons
                     break;
                  case 2112:           // neutrons
                     ezdc2[1] += ekin; // nucleons
                     break;
                  case 211: // pions
                  case -211: ezdc2[2] += fArgs[9]; break;
                  default: break;
                  }
               }

            } // end of zdc entrance

            fTree->Fill(fArgs);

         } // end of primary particles

         //     fTree->Fill(fArgs);
      }
      //  cout << " .finished!" << endl;

      nevents++;
   }
   Double_t b = 0;

   if (fMCEventHeader) {
      b = ((FairMCEventHeader *)fMCEventHeader)->GetB();
   }

   if (fTreeSummary) {

      Int_t nTpcTracks = 0;
      Int_t trackId;

      if (fMCTpcPointArray) {

         std::map<Int_t, Bool_t> tpc_tracks;
         tpc_tracks.clear();

         Int_t i_output = 5;
         for (Int_t i = 0; i < fMCTpcPointArray->GetEntriesFast(); i++) {

            trackId = ((TpcPoint *)fMCTpcPointArray->UncheckedAt(i))->GetTrackID();

            if ((trackId >= 0) && (trackId < nMCTracks)) {

               if (((MpdMCTrack *)fMCTrackArray->UncheckedAt(trackId))->GetMotherId() == -1) {

                  if (tpc_tracks.find(trackId) == tpc_tracks.end()) tpc_tracks[trackId] = kTRUE;
               }
            } else {
               if (i_output > 0) {
                  cout << "-E- Wrong TpcPoint GetTrackID = " << trackId << "  for i=" << i << "   b=" << b << endl;
                  i_output--;
               }
            } // else
         }    // for (Int_t i=0;i<fMCTpcPointArray->GetEntriesFast();i++)
         nTpcTracks = tpc_tracks.size();
      }

      fTreeSummary->Fill(elossall, elosszdc1, elosszdc2, nTpcTracks, b, e, ezdc1[0], ezdc1[1], ezdc1[2], ezdc1[3],
                         ezdc2[0], ezdc2[1], ezdc2[2], ezdc2[3]);
      //    cout << " . " << (fTreeSummary) << endl;

      fH3->Fill(elossall, nTpcTracks, b);
   } else
      cout << " !";

   if (fVerbose > 3) cout << "-Info- MpdZdcTstSim::Exec end" << endl;
}

void MpdZdcTstSim::Finish()
{
   fMCTrackArray->Clear();
   fMCZdcPointArray->Clear();
   if (fMCTpcPointArray) fMCTpcPointArray->Clear();
}

ClassImp(MpdZdcTstSim)
