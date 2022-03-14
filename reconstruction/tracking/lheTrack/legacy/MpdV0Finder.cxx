// Task for finding secondary vertices
#include "MpdV0Finder.h"

//__________________________________________________________________________

MpdV0Finder::MpdV0Finder(const char *name, Int_t iVerbose)
   : FairTask(name, iVerbose), V0(0), fEvent(0), fMpdTracks(0), fKFTracks(0), fTofMatching(0), fMCTracks(0), pi(0),
     pr(0)
{
}

//__________________________________________________________________________

MpdV0Finder::~MpdV0Finder()
{
   //  delete V0; V0 = 0;
   fV0Cont->Clear("C");
   // fMpdTracks->Clear("C");
   pi->Clear("C");
   pr->Clear("C");
}

//__________________________________________________________________________

InitStatus MpdV0Finder::Init()
{
   cout << endl;
   cout << " MpdV0Finder::Init  " << endl;
   npr          = 0;
   npi          = 0;
   PVert        = 0.;
   vtx          = 0.;
   fNMPDtracks  = 0;
   fNtracks     = 0;
   fV0Cont      = 0x0;
   fMpdTracks   = 0x0;
   fKFTracks    = 0x0;
   fTofMatching = 0x0;
   fMCTracks    = 0x0;
   fNMPDtracks  = 0;

   nvert = 0;
   Nreco = 0;
   Nlam  = 0;

   pi = new TClonesArray("MpdParticle"); // put in exec to check for mem leak
   pr = new TClonesArray("MpdParticle");

   FairRootManager *manager = FairRootManager::Instance();
   fV0Cont                  = new TClonesArray("MpdV0", nvert); // container for v0s
   manager->Register("V02", "MpdV0", fV0Cont, kTRUE);           // Register v0s
   // fMpdTracks = new TClonesArray("MpdTrack"); //container for MpdTracks               if used in reco chain this is
   // for probability
   fPrimVertex  = (TClonesArray *)manager->GetObject("Vertex");
   fKFTracks    = (TClonesArray *)manager->GetObject("TpcKalmanTrack");
   fTofMatching = (TClonesArray *)manager->GetObject("TOFMatching");
   // fEvent = (MpdEvent*) manager->GetObject("MPDEvent.");

   // fEvent = (MpdEvent*) manager->ActivateBranch("MPDEvent.");
   fEvent = (MpdEvent *)manager->GetObjectFromInTree("MPDEvent.");
   // fMpdTracks = fEvent->GetGlobalTracks();
   // cout<<" GlobalTracks = "<<fMpdTracks->GetEntriesFast()<<endl;

   fMCTracks = (TClonesArray *)manager->GetObject("MCTrack");
   // fMpdTracks = (TClonesArray*) manager->GetObject("MpdTracks");

   return kSUCCESS;
}

//__________________________________________________________________________

InitStatus MpdV0Finder::ReInit()
{
   return kSUCCESS;
}

//__________________________________________________________________________

void MpdV0Finder::Reset()
{
   cout << " MpdV0Finder::Reset  " << endl;
   // fV0Cont->Delete();
   // fMpdTracks->Delete();
   fV0Cont->Clear("C");
   // fMpdTracks->Clear("C");
   pi->Clear("C");
   pr->Clear("C");
   // fMCTracks->Clear("C");
   // fKalmanTrs->Clear();
   // fKalmanTrs->Expand(0);
   // fMpdTracks->Clear();
   Nok1  = 0;
   Nmc1  = 0;
   Nlam1 = 0;
   nvert = 0;
   npr   = 0;
   npi   = 0;

   // pr.resize(0);
   // pi.resize(0);
}

//__________________________________________________________________________

void MpdV0Finder::Finish()
{

   cout << "***For all events*** " << endl;
   Double_t Nbad = Nreco - Nok; // accumulates
   cout << "  Nreco = " << Nreco << "  Nmc = " << Nmc << "  Nok = " << Nok << "  Nbad = " << Nbad << endl;
   Double_t Cont = (Double_t)Nbad / (Double_t)Nreco; // accumulates
   Double_t Eff  = (Double_t)Nok / (Double_t)Nmc;    // accumulates
   cout << "Efficiency = " << Eff * 100 << "%  Contamination = " << Cont * 100 << "%" << endl;
   cout << "MClambdas = " << Nlam << endl;
}

//__________________________________________________________________________

void MpdV0Finder::Exec(Option_t *option)
{
   static int eventCounter = 0;
   ++eventCounter;
   cout << endl;
   if (eventCounter != 796) return;
   time_t now = time(0);     // current date/time based on current system
   char  *dt  = ctime(&now); // convert now to string form
   cout << "The local date and time is: " << dt << endl;
   cout << "Event Number " << eventCounter << endl;
   cout << endl;
   // pi = new TClonesArray("MpdParticle");        //put in exec to check for mem leak
   // pr = new TClonesArray("MpdParticle");

   Reset();
   // GetPrimaryVertex();
   MpdVertex *PrimVtx = (MpdVertex *)fPrimVertex->UncheckedAt(0);
   PVert.SetXYZ(PrimVtx->GetX(), PrimVtx->GetY(), PrimVtx->GetZ());
   FindMCVertices(); // find mc vertices with decay daughters in tpc
   cout << "MC lambdas = " << Nlam1 << endl;
   cout << "Lambdas with both daughters in tpc = " << Nmc1 << endl;
   // FillMpdTracks(); //like in fill dst just used for probability

   //  MpdEvent *ev = new MpdEvent();
   // ev = (MpdEvent*) fEvent->UncheckedAt(0);
   // cout<<fEvent->GetEntries()<<endl;
   // MpdEvent *ev = (MpdEvent*) fEvent->UncheckedAt(1);
   cout << "fEvent->GetRunInfoRunId() = " << fEvent->GetRunInfoRunId() << endl;
   fMpdTracks  = (TClonesArray *)fEvent->GetGlobalTracks();
   fNMPDtracks = fMpdTracks->GetEntriesFast();
   cout << " GlobalTracks = " << fNMPDtracks << endl;

   // cout << "MpdTracks = " << fNMPDtracks << endl;

   SelectTracks();
   // Int_t nprots = pr.size();
   // Int_t npions = pi.size();
   Int_t nprots = pr->GetEntriesFast();
   Int_t npions = pi->GetEntriesFast();

   cout << "Protons for lambda = " << nprots << endl;
   cout << "Pions for lambda = " << npions << endl;

   vector<MpdParticle *> cont; // container for daughters
   for (Int_t i = 0; i < nprots; i++) {
      // MpdParticle *prot = pr.at(i);
      MpdParticle *prot = (MpdParticle *)pr->UncheckedAt(i);
      Int_t        ch1  = prot->GetCharge();
      // Int_t B = ch1; //baryon number
      if (ch1 > 0)
         prot->SetPdg(protonPDG);
      else
         prot->SetPdg(-protonPDG);
      prot->SetMass();
      for (Int_t j = 0; j < npions; j++) {

         // MpdParticle *pion = pi.at(j);
         MpdParticle *pion = (MpdParticle *)pi->UncheckedAt(j);
         Int_t        ch2  = pion->GetCharge();
         if (ch1 * ch2 > 0) continue;
         if (ch2 < 0)
            pion->SetPdg(-pionPDG);
         else
            pion->SetPdg(pionPDG);
         pion->SetMass();

         cont.clear();
         cont.push_back(prot);
         cont.push_back(pion);
         MpdParticle lambda;
         // cout<<i<<"  "<<j<<endl;
         Double_t chi2 = abs(lambda.BuildMother(cont)); // why abs?
         if (chi2 > V0Chi2Cut) continue;
         Double_t mass = lambda.GetMass();
         if (minmassCut > mass || mass > maxmassCut) continue;

         Double_t DCA = lambda.Dca();
         if (abs(DCA) > DCAV0Cut) continue;

         vtx.SetXYZ(lambda.Getx()(0, 0), lambda.Getx()(1, 0), lambda.Getx()(2, 0));
         Double_t xyz[3];
         vtx.GetXYZ(xyz);

         TVector3 R = vtx - PVert; // decay distance from primvert
         Double_t r[3];
         R.GetXYZ(r);
         if (RminCut > R.Mag() || R.Mag() > RmaxCut) continue;

         TVector3 P = lambda.Momentum3();
         if (P.Mag() < PminCut) continue;
         Double_t p[3];
         P.GetXYZ(p);
         Double_t cosa = (R.Dot(P)) / (P.Mag() * R.Mag()); // pointing angle
         if (cosa < cosaCut)
            continue; // negative cosa and cut or absolute?a lot of true are with positive cosa
                      // TMatrixD covar = lambda.GetC();//MpdKfV0Fitter::Instance()->GetCovariance();

         MpdKalmanTrack *tr1 = (MpdKalmanTrack *)fKFTracks->UncheckedAt(
            prot->GetIndx()); // this is necessary because mcid and reco number should be separated
         MpdKalmanTrack *tr2 = (MpdKalmanTrack *)fKFTracks->UncheckedAt(pion->GetIndx());

         Int_t id0 = CheckMCTracks(tr1, tr2);
         //  cout<<"nvert = "<<nvert<<endl;
         V0 = new ((*fV0Cont)[nvert]) MpdV0();
         V0->SetV0(xyz, chi2, p, DCA, cosa, mass, r, ch1, prot->GetIndx(), pion->GetIndx(),
                   id0); // ndf = 2n-3(is there apriory? ) or 2n so its either 4 or 1??? B is baryon number
         ++nvert;        // gets reset for every event
         ++Nreco;        // accumulates
         cout << ".";
      }
   }
   cout << endl;
   cout << "Reconstructed " << nvert << " lambda particles." << endl;
   Double_t Eff1 = (Double_t)Nok1 / (Double_t)Nmc1;
   // Eff1vsZ->Fill(pvert[2],Eff1);
   Double_t Bad1  = nvert - Nok1;
   Double_t Cont1 = (Double_t)Bad1 / (Double_t)nvert;
   cout << "Eff1 = " << Eff1 << endl; // pseudo-efficiency of lambda for this event
   cout << "Cont1 = " << Cont1 << endl;
}

//__________________________________________________________________________

void MpdV0Finder::SelectTracks()
{

   MpdTpcKalmanFilter *recoTpc = new MpdTpcKalmanFilter(); // for refit where should refit?
   // recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
   // recoTpc->FillGeoScheme();

   for (Int_t i = 0; i < fNMPDtracks; i++) {
      MpdKalmanTrack *tr = (MpdKalmanTrack *)fKFTracks->UncheckedAt(i);
      if (tr->GetNofHits() < NofHitsCut) continue;
      MpdTrack *mpdtr = (MpdTrack *)fMpdTracks->UncheckedAt(i);
      Int_t     n     = GetHighestProb(mpdtr);
      if (n == 0 || n == 1 || n == 4) continue; // no need for electrons or kaons
      MpdKalmanTrack trRefit = *tr;
      Double_t       chi2v   = tr->GetChi2Vertex();
      Double_t       Pt      = tr->Pt();
      Double_t       Mom     = tr->Momentum();
      // protons
      if (n == 2) {
         if (protonPtminCut > Pt || Pt > protonPtmaxCut) continue;
         if (protonPminCut > Mom || Mom > protonPmaxCut) continue;
         if (protonChi2vminCut > chi2v || chi2v > protonChi2vmaxCut) continue;
         if (PCA(tr, PVert) < protonIPCut) continue; // cut on particle or charge???
         recoTpc->Refit(&trRefit, 0.93827, 1);       // refit // Where should I refit???
         MpdParticle *part =
            new ((*pr)[npr]) MpdParticle(*tr, i); // MC ID or reco number??? i or GetTrackID <-track id is bad because
                                                  // later you have to loop and check number of track
         //  MpdParticle *part = new ((*pr)[npr]) MpdParticle(trRefit,i);
         npr++;
         // delete part;
         // pions
      } else {
         if (pionPtminCut > Pt || Pt > pionPtmaxCut) continue;
         if (pionPminCut > Mom || Mom > pionPmaxCut) continue;
         if (pionChi2vminCut > chi2v || chi2v > pionChi2vmaxCut) continue;
         if (PCA(tr, PVert) < pionIPCut) continue;
         // recoTpc->Refit(&trRefit, 0.93827, 1);
         MpdParticle *part = new ((*pi)[npi]) MpdParticle(*tr, i);
         // MpdParticle *part = new ((*pr)[npr]) MpdParticle(&trRefit,i);
         npi++;
         // delete part;
      }
      // delete part;   //check if mem leak is fixed
   } // trackloop

   delete recoTpc; // check if mem leak is fixed
}

//__________________________________________________________________________

Double_t MpdV0Finder::PCA(MpdKalmanTrack *tr, TVector3 point)
{
   // Create MpdHelix for 3D DCA as done in AnalXiFull
   Double_t r   = tr->GetPosNew();
   Double_t phi = tr->GetParam(0) / r;
   Double_t x   = r * TMath::Cos(phi);
   Double_t y   = r * TMath::Sin(phi);
   Double_t dip = tr->GetParam(3);
   Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
   cur *= TMath::Abs(tr->GetParam(4));
   TVector3 o(x, y, tr->GetParam(1));
   Int_t    h = (Int_t)TMath::Sign(1.1, tr->GetParam(4));
   MpdHelix helix(cur, dip, tr->GetParam(2) - TMath::PiOver2() * h, o, h);
   TVector3 pca;
   Double_t s = helix.pathLength(point);
   pca        = helix.at(s);
   pca -= point;
   Double_t dist = pca.Mag();

   return dist;
}

//__________________________________________________________________________

Int_t MpdV0Finder::GetHighestProb(MpdTrack *track)
{
   // returns the most probable particle. a cut can still be applied
   Double_t electron = track->GetPidProbElectron();
   Double_t proton   = track->GetPidProbProton();
   Double_t pion     = track->GetPidProbPion();
   Double_t kaon     = track->GetPidProbKaon();
   Double_t Probs[]  = {0, electron, proton, pion, kaon};
   Double_t maxprob  = *max_element(Probs, Probs + 5);
   Int_t    i        = distance(Probs, max_element(Probs, Probs + 5));

   if (maxprob < ProbCut)
      return 0;
   else
      return i;
}

//__________________________________________________________________________
/*
void MpdV0Finder::FillMpdTracks() {
//only for barrel, done for pid (needs to be a separate task in reco imho)
    Int_t nReco = fKFTracks ? fKFTracks->GetEntriesFast() : 0;
    Int_t nMatching = fTofMatching ? fTofMatching->GetEntriesFast() : 0;
    MpdTofMatchingData *pMatchingData;
    bool matchingDataExist;
    Int_t m = 0;
    for (Int_t i = 0; i < nReco; i++) {
        MpdKalmanTrack *kftrack = (MpdKalmanTrack*) fKFTracks->UncheckedAt(i);
        MpdTrack *track = new ((*fMpdTracks)[m]) MpdTrack(); //instead of AddGlobalTrack()
        MpdParticleIdentification *identificator = new MpdParticleIdentification();
        track->SetID(kftrack->GetTrackID());
        Float_t Ppi, Pk, Pe, Pp;
        if (!identificator->GetTpcProbs(kftrack->Momentum3().Mag(), kftrack->GetPartID(), kftrack->GetNofHits(), Ppi,
Pk, Pp, Pe, 0)) { track->SetTPCpidProb(Pe, Ppi, Pk, Pp, BIT(2));  } matchingDataExist = false; for (Int_t tofIndex = 0;
tofIndex < nMatching; tofIndex++) { pMatchingData = (MpdTofMatchingData*) fTofMatching->UncheckedAt(tofIndex); if
(pMatchingData->GetKFTrackIndex() == i) {matchingDataExist = true; break;} } // first matching if (matchingDataExist) {
            track->SetTofMass2(pMatchingData->GetMass2());
            track->SetTofHitIndex(pMatchingData->GetTofHitIndex());
            if (!identificator->GetTofProbs(pMatchingData->GetMomentum().Mag(),  pMatchingData->GetMass2(),
kftrack->GetNofHits(), Ppi, Pk, Pp, Pe, 0)) { track->SetTOFpidProb(Pe, Ppi, Pk, Pp, BIT(1)); }
          }
        Float_t tpcProbs[4] = {track->GetTPCPidProbPion(), track->GetTPCPidProbKaon(), track->GetTPCPidProbProton(),
track->GetTPCPidProbElectron()}; Float_t tofProbs[4] = {track->GetTOFPidProbPion(), track->GetTOFPidProbKaon(),
track->GetTOFPidProbProton(), track->GetTOFPidProbElectron()}; Float_t combProbs[4]; //probabilities combined from TOF &
TPC identificator->GetCombinedProbs(tofProbs, tpcProbs, combProbs, 4); Ppi = combProbs[0]; Pk = combProbs[1]; Pp =
combProbs[2]; Pe = combProbs[3]; track->SetCombPidProb(Pe, Ppi, Pk, Pp); m++;
    }
}
*/
//__________________________________________________________________________

Int_t MpdV0Finder::CheckMCTracks(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2)
{

   // checks if the found vertex of these two tracks is a MC(true) vertex primary lambda have mum -1
   FairMCTrack *mctr1 = (FairMCTrack *)fMCTracks->UncheckedAt(tr1->GetTrackID());
   FairMCTrack *mctr2 = (FairMCTrack *)fMCTracks->UncheckedAt(tr2->GetTrackID());
   Int_t        mum1  = mctr1->GetMotherId();
   Int_t        mum2  = mctr2->GetMotherId();
   if (mum1 == -1 || mum2 == -1) return -2; //  we have a primary track, so vert is not true
   TVector3 vert1, vert2;
   mctr1->GetStartVertex(vert1);
   mctr2->GetStartVertex(vert2);
   if (vert1 != vert2) return -2; // vertex matches
   FairMCTrack *mcmumTr1 = (FairMCTrack *)fMCTracks->UncheckedAt(mum1);
   FairMCTrack *mcmumTr2 = (FairMCTrack *)fMCTracks->UncheckedAt(mum2);
   if (abs(mcmumTr1->GetPdgCode()) != lambdaPDG) return -2;
   ++Nok;
   ++Nok1;
   return mum1;
}

//__________________________________________________________________________

void MpdV0Finder::FindMCVertices()
{

   Int_t fNMC = fMCTracks->GetEntriesFast();
   cout << "Nmctracks = " << fNMC << endl;
   // find all lambdas first
   for (Int_t i = 0; i < fNMC; ++i) {
      FairMCTrack *mctr0 = (FairMCTrack *)fMCTracks->UncheckedAt(i);
      if (mctr0->GetPdgCode() != lambdaPDG) continue;
      ++Nlam;
      ++Nlam1;
   }
   for (Int_t i = 0; i < fNMC; ++i) { // find protons
      FairMCTrack *mctr1    = (FairMCTrack *)fMCTracks->UncheckedAt(i);
      Int_t        mctrpdg1 = mctr1->GetPdgCode();
      if (abs(mctrpdg1) != protonPDG) continue;
      Int_t mom1 = mctr1->GetMotherId();
      if (mom1 == -1) continue;
      // is the mommy lambda
      if (abs(((FairMCTrack *)fMCTracks->UncheckedAt(mom1))->GetPdgCode()) != lambdaPDG) continue;
      for (Int_t j = 0; j < fNMC; ++j) { // find pions for protons
         FairMCTrack *mctr2    = (FairMCTrack *)fMCTracks->UncheckedAt(j);
         Int_t        mctrpdg2 = mctr2->GetPdgCode();
         if (abs(mctrpdg2) != pionPDG) continue;
         Int_t mom2 = mctr2->GetMotherId();
         if (mom2 == -1) continue;
         if (abs(((FairMCTrack *)fMCTracks->UncheckedAt(mom2))->GetPdgCode()) != lambdaPDG) continue;
         // both daughters have detector points??
         if (mctr1->GetNPoints(kTPC) && mctr2->GetNPoints(kTPC)) {
            // daughters have same vertex
            TVector3 vert1, vert2;
            mctr1->GetStartVertex(vert1);
            mctr2->GetStartVertex(vert2);
            if (vert1 != vert2) continue;

            ++Nmc;
            ++Nmc1;

            // mcPt1->Fill(mctr1->GetPt());
            // mcP1 ->Fill(mctr1->GetP());
            // mcPt2->Fill(mctr2->GetPt());
            // mcP2 ->Fill(mctr2->GetP());
            // TVector3 PVert; //primary vertex
            // PVert.SetXYZ(pvert[0], pvert[1], pvert[2]);
            TVector3 mcR  = PVert - vert1;
            Double_t dist = mcR.Mag();
            // hmcR->Fill(dist);
            TLorentzVector lorVec[2];
            FairMCTrack   *trs[2] = {mctr1, mctr2};
            for (Int_t i = 0; i < 2; ++i) {
               FairMCTrack *tr = trs[i];
               lorVec[i].SetXYZM(tr->GetPx(), tr->GetPy(), tr->GetPz(), tr->GetMass());
            }
            TLorentzVector lorInv = lorVec[0] + lorVec[1];
            TVector3       mcP    = lorInv.Vect();
            Double_t       mccosa = (mcR.Dot(mcP)) / (mcP.Mag() * mcR.Mag());
            // mccos->Fill(mccosa);
         }
      }
   }
}
ClassImp(MpdV0Finder);