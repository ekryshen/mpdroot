#include <iostream>
#include "MpdV0.h"
#include "MpdV0Maker.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdVertex.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "TH1.h"
#include "TH2.h"

ClassImp(MpdV0Maker);

using namespace std;
MpdV0Maker::MpdV0Maker(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   // Create V0 engine
}

void MpdV0Maker::UserInit()
{
   cout << "[MpdV0Maker]: Initialization ... " << endl;

   // mParams.ReadFromFile(mParamConfig);
   // mParams.Print();

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   if (mFillEff) {
      const int   nPtbin = 100;
      const float pTmin  = 0.;
      const float pTmax  = 5.;

      // General QA
      mhChi2 = addHist(new TH2F("hChi2", "#chi^{2};#chi^{2};p_{T} (GeV/#it{c})", 100, 0., 50., nPtbin, pTmin, pTmax));
      mhChi2True = addHist((TH2F *)mhChi2->Clone(Form("%sTrue", mhChi2->GetName())));

      mhAlpha = addHist(
         new TH2F("hAlpha", "#alpha distribution;#alpha (rad);p_{T} (GeV/#it{c})", 100, 0., 1., nPtbin, pTmin, pTmax));
      mhAlphaTrue = addHist((TH2F *)mhAlpha->Clone(Form("%sTrue", mhAlpha->GetName())));

      mhDist = addHist(new TH2F("hDist", "track DCA;DCA (cm);p_{T} (GeV/#it{c})", 100, 0., 10., nPtbin, pTmin, pTmax));
      mhDistTrue = addHist((TH2F *)mhDist->Clone(Form("%sTrue", mhDist->GetName())));

      mhMassEE = addHist(
         new TH2F("mEE", "m_{ee};m_{ee} (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", 100, 0., 0.3, nPtbin, pTmin, pTmax));
      mhMassEETrue = addHist((TH2F *)mhMassEE->Clone(Form("%sTrue", mhMassEE->GetName())));

      mhArmPo     = addHist(new TH2F("Armenteros", "Armenteros", 100, -1, 1, 100, 0, 0.3));
      mhArmPoTrue = addHist((TH2F *)mhArmPo->Clone(Form("%sTrue", mhArmPo->GetName())));

      mhAsym = addHist(new TH2F("Asymetry", "Asymmetry;Asymmetry;p_{T} (GeV/#it{c})", 200, 0, 1, nPtbin, pTmin, pTmax));
      mhAsymTrue = addHist((TH2F *)mhAsym->Clone(Form("%sTrue", mhAsym->GetName())));

      mhCosPsi =
         addHist(new TH2F("cosPsi", "cos(#psi);cos(#psi);p_{T} (GeV/#it{c})", 100, -1., 1., nPtbin, pTmin, pTmax));
      mhCosPsiTrue = addHist((TH2F *)mhCosPsi->Clone(Form("%sTrue", mhCosPsi->GetName())));

      mhCutEff = addHist(new TH2F("hCutEff", "track cut efficiency;Cut ID;p_{T} (GeV/#it{c});Efficiency", 10, 0., 10,
                                  nPtbin, pTmin, pTmax));
   }

   cout << "[MpdV0Maker]: Initialization done " << endl << endl;
}
//--------------------------------------
void MpdV0Maker::ProcessEvent(MpdAnalysisEvent &event)
{
   // conbine tracks and apply selection for V0
   if (!mKF) { // not yet initialized
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
   }

   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);

   mMC = (event.fMCEventHeader != nullptr);
   if (mMC) {
      mMCTracks = event.fMCTrack;
   }

   TClonesArray *mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int           ntr              = mMpdGlobalTracks->GetEntriesFast();
   TClonesArray *mKalmanTracks    = event.fTPCKalmanTrack;

   // Clear or create container for output V0s
   TClonesArray *mV0 = new TClonesArray("MpdV0", ntr / 2);
   MpdV0         v;
   int           iv0 = 0;
   for (long int i = 0; i < ntr - 1; i++) {
      MpdTrack          *mpdtrack1 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr1       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack1)) {
         continue;
      }
      for (long int j = i + 1; j < ntr; j++) {
         MpdTrack          *mpdtrack2 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(j);
         MpdTpcKalmanTrack *tr2       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
         if (!selectTrack(mpdtrack2)) {
            continue;
         }
         if (createSelectV0(mpdtrack1, tr1, mpdtrack2, tr2, v)) {
            v.setTr1(i);
            v.setTr2(j);
            new ((*mV0)[iv0++]) MpdV0(v);
         }
      }
   }
   event.fV0 = mV0;
   // cout << event.fV0->GetEntries() << " V0-s!\n";
}
//===============================================
void MpdV0Maker::Finish()
{
   // Post-scan processing not needed
}
//===============================================
bool MpdV0Maker::selectTrack(MpdTrack *mpdtrack)
{

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (mpdtrack->GetNofHits() < mNofHitsCut) return false; // nhits > 10
   if (fabs(mpdtrack->GetEta()) > mEtaCut) return false;   //|eta| < 1.0
   if (pt < mPtminCut) return false;                       // pT > 50 MeV/c

   // PID selection??
   // int charge;
   // if (mpdtrack->GetPt() < 0)
   //    charge = 1;
   // else
   //    charge = -1;

   // bool isGoodPID;
   // if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) {
   //    isGoodPID =
   //       mPID->FillProbs(pt * TMath::CosH(mpdtrack->GetEta()), mpdtrack->GetdEdXTPC(), mpdtrack->GetTofMass2(),
   //       charge);
   // } else {
   //    isGoodPID = mPID->FillProbs(pt * TMath::CosH(mpdtrack->GetEta()), mpdtrack->GetdEdXTPC(), charge);
   // }

   // bool isElectron = false;
   // if (mMC) { // same for true electron tracks
   //    long int prim1 = mpdtrack->GetID();
   //    if (prim1 >= 0) {
   //       isElectron = abs((static_cast<MpdMCTrack *>(mMCTracks->At(prim1)))->GetPdgCode()) == 11;
   //    }
   // }

   // if (isGoodPID && mPID->GetProbEl() < mParams.mProbElCut) {
   //    return false;
   // }

   // float dEdx = dEdx_sigma(mpdtrack->GetdEdXTPC(), sqrt(pow(pt, 2) + pow(mpdtrack->GetPz(), 2)));
   // if ((fabs(dEdx) < mParams.mdEdxSigmaCut) &&
   //     (fabs(Beta_sigma(mpdtrack->GetTofBeta(), sqrt(pow(pt, 2) + pow(mpdtrack->GetPz(), 2)))) <
   //         mParams.mBetaSigmaCut &&
   //      (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6))) {
   //    return true;
   // } else {
   //    return false;
   // }
   // return false;
   return true;
}
//============================================================================================================
bool MpdV0Maker::createSelectV0(MpdTrack *tr1, MpdTpcKalmanTrack *ktr1, MpdTrack *tr2, MpdTpcKalmanTrack *ktr2,
                                MpdV0 &v)
{
   // Construct and check V0

   // Use opposite charge tracks
   int charge1 = tr1->GetCharge(), charge2 = tr2->GetCharge();

   // reject same sign pairs
   if (charge1 * charge2 > 0) return false;

   // Will be used for extrapolation
   MpdTpcKalmanTrack trCorK1(*ktr1);
   MpdHelix          helix1 = MakeHelix(trCorK1);
   MpdParticle       el1(trCorK1, 0);
   el1.SetPdg(-charge1 * mDefPDG);
   el1.SetMass();

   MpdTpcKalmanTrack trCorK2(*ktr2);
   MpdHelix          helix2 = MakeHelix(trCorK2);
   MpdParticle       el2(trCorK2, 0);
   el2.SetPdg(-charge2 * mDefPDG);
   el2.SetMass();

   // pair
   mPartK.clear();
   mPartK.emplace_back(&el1);
   mPartK.emplace_back(&el2);

   MpdParticle gamEE;
   float       chi2 = TMath::Abs(gamEE.BuildMother(mPartK));
   v.setChi2(chi2);
   float pt = gamEE.Pt();
   if (mFillEff) {
      mhChi2->Fill(chi2, pt);
      mhCutEff->Fill(0., pt);
   }

   if (pt < 0.005) { // to avoid fpe
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(1., pt);
   }

   bool isTrue = false; // is true conv pair?
   if (mMC) {           // same for true electrontracks
      long int matched1       = tr1->GetID();
      long int matched2       = tr2->GetID();
      long int commonParentId = FindCommonParent(matched1, matched2, mMCTracks);
      v.setMatched1(matched1);
      v.setMatched2(matched2);
      v.setCommonParent(commonParentId);
      if (commonParentId >= 0) // there is common parent
      {
         isTrue = true;
      }
   }
   if (isTrue && mFillEff) {
      mhChi2True->Fill(chi2, pt);
   }
   if (chi2 > mChi2Cut) {
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(2., pt);
   }

   TVector3 v0(gamEE.Getx()(0, 0), gamEE.Getx()(1, 0), gamEE.Getx()(2, 0));
   v0 -= mPrimaryVertex;

   float rConv = TMath::Sqrt(pow(gamEE.Getx()(0, 0), 2) + pow(gamEE.Getx()(1, 0), 2));
   if (rConv < mMinR2Cut || rConv > mMaxR2Cut) {
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(3., pt);
   }

   float angle = v0.Angle(gamEE.Momentum3());
   v.setPA(angle);
   if (mFillEff) {
      mhAlpha->Fill(angle, pt);
      if (isTrue) { // same for true electrontracks
         mhAlphaTrue->Fill(angle, pt);
      }
   }
   if (angle > mAlphaCut) {
      return false;
   }

   if (mFillEff) {
      mhCutEff->Fill(4., pt);
   }

   // if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
   //   return kFALSE;  // line cut to exclude regions where we do not reconstruct
   // } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) -
   // fLineCutZValueMin)){
   //   return kFALSE;
   // }

   std::pair<float, float> paths = helix1.pathLengths(helix2);
   TVector3                p1    = helix1.at(paths.first);
   TVector3                p2    = helix2.at(paths.second);
   p1 -= p2;
   float dist = p1.Mag(); // Closest distance between daughters

   v.setDaughterDCA(dist);

   if (mFillEff) {
      mhDist->Fill(dist, pt);
      if (isTrue) { // same for true electrontracks
         mhDistTrue->Fill(dist, pt);
      }
   }
   if (dist > mDistCut) {
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(5., pt);
   }

   v.setMass(gamEE.GetMass());
   if (mFillEff) {
      mhMassEE->Fill(gamEE.GetMass(), pt);
      if (isTrue) { // same for true electrontracks
         mhMassEETrue->Fill(gamEE.GetMass(), pt);
      }
   }
   if (gamEE.GetMass() > mMassCut) {
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(6., pt);
   }
   // Pair_chi_1[n_ks] = el1.Chi2Vertex(vertex);
   // Pair_chi_2[n_ks] = el2.Chi2Vertex(vertex);
   // A-P cut
   //  Gamma selection based on QT from Armenteros
   // propagate trCorK1,trCorK2 to conversion point

   MpdKalmanHit hitTmp;
   hitTmp.SetType(MpdKalmanHit::kFixedR);
   hitTmp.SetPos(trCorK1.GetPos());
   trCorK1.SetParamNew(*trCorK1.GetParam());
   trCorK1.SetPos(trCorK1.GetPosNew());
   trCorK1.ReSetWeight();
   //  TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK1, &mKHit, kFALSE, kFALSE)) {
      return false;
   }
   //   trCorK1.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m1 = trCorK1.Momentum3();

   hitTmp.SetPos(trCorK2.GetPos());
   trCorK2.SetParamNew(*trCorK2.GetParam());
   trCorK2.SetPos(trCorK2.GetPosNew());
   trCorK2.ReSetWeight();
   TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK2, &mKHit, kFALSE, kFALSE)) {
      return false;
   }
   //   trCorK2.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m2 = trCorK2.Momentum3();

   float qt, alpha;
   ArmenterosPodolanski(m1, m2, qt, alpha);
   v.setArmenteros(alpha, qt);

   if (mFillEff) {
      mhArmPo->Fill(alpha, qt);
      if (isTrue) {
         mhArmPoTrue->Fill(alpha, qt);
      }
   }
   if (!ArmenterosQtCut(qt, alpha, gamEE)) {
      // return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(7., pt);
   }

   // Asymmetry cut
   float asym1 = m1.Mag() / gamEE.Momentum();
   float asym2 = m2.Mag() / gamEE.Momentum();
   v.setAsymmetry(asym1, asym2);
   if (mFillEff) {
      mhAsym->Fill(asym1, pt);
      mhAsym->Fill(asym2, pt);
      if (isTrue) {
         mhAsymTrue->Fill(asym1, pt);
         mhAsymTrue->Fill(asym2, pt);
      }
   }
   if (!(AsymmetryCut(asym1, pt) && AsymmetryCut(asym2, pt))) {
      return kFALSE;
   }
   if (mFillEff) {
      mhCutEff->Fill(8., pt);
   }

   float cospsi = CosPsiPair(m1, m2);
   v.setCospsi(cospsi);
   if (mFillEff) {
      mhCosPsi->Fill(cospsi, pt);
      if (isTrue) { // same for true electrontracks
         mhCosPsiTrue->Fill(cospsi, pt);
      }
   }
   if (cospsi < mCosPsiCut) {
      return false;
   }
   if (mFillEff) {
      mhCutEff->Fill(9., pt);
   }

   v.setVertex(gamEE.Getx()(0, 0), gamEE.Getx()(1, 0), gamEE.Getx()(2, 0));

   v.SetXYZT((gamEE.Pt()) * TMath::Cos(gamEE.Phi()), (gamEE.Pt()) * TMath::Sin(gamEE.Phi()),
             TMath::Sign(TMath::Sqrt(gamEE.Momentum() * gamEE.Momentum() - gamEE.Pt() * gamEE.Pt()),
                         TMath::Cos(gamEE.Theta())),
             gamEE.Momentum());

   return true;
}

MpdHelix MpdV0Maker::MakeHelix(const MpdKalmanTrack &tr) const
{
   float r   = tr.GetPosNew();
   float phi = tr.GetParam(0) / r;
   float x   = r * TMath::Cos(phi);
   float y   = r * TMath::Sin(phi);
   float dip = tr.GetParam(3);
   float cur = 0.3 * 0.01 * 5 / 10; // 5 kG
   cur *= TMath::Abs(tr.GetParam(4));
   TVector3 o(x, y, tr.GetParam(1));
   Int_t    h = (Int_t)TMath::Sign(1.1, tr.GetParam(4));
   MpdHelix helix(cur, dip, tr.GetParam(2) - TMath::PiOver2() * h, o, h);
   return helix;
}

void MpdV0Maker::ArmenterosPodolanski(TVector3 &m1, TVector3 &m2, float &qt, float &alpha) const
{

   alpha = 0., qt = 0.;

   TVector3 s = m1 + m2;

   float pn  = m1.Mag();
   float pln = m1.Dot(s);
   float plp = m2.Dot(s);

   if (pn == 0.0) return;
   alpha    = (plp - pln) / (plp + pln);
   float sm = s.Mag();
   if (sm > 0) {
      qt = m1.Cross(s).Mag() / sm;
   }
}

//________________________________________________________________________
bool MpdV0Maker::ArmenterosQtCut(float qt, float alpha, MpdParticle &part) const
{ // Armenteros Qt Cut
   // if(mParams.mDo2DQt){
   //   if(mParams.mDoQtGammaSelection==1){
   //     if (
   //     !(TMath::Power(photon->GetArmenterosAlpha()/mParams.mMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/mParams.mQtMax,2)
   //     < 1) ){
   //       return false;
   //     }
   //   } else if(mParams.mDoQtGammaSelection==2){
   //     float qtMaxPtDep = mParams.mQtPtMax*photon->GetPhotonPt();
   //     if (qtMaxPtDep > mParams.mQtMax)
   //       qtMaxPtDep      = mParams.mQtMax;
   //     if (
   //     !(TMath::Power(photon->GetArmenterosAlpha()/mParams.mMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/qtMaxPtDep,2)
   //     < 1) ){
   //       return false;
   //     }
   //   }
   // } else {
   //   if(mParams.mDoQtGammaSelection==1){
   //     if(photon->GetArmenterosQt()>mParams.mQtMax){
   //       return false;
   //     }
   //   } else if(mParams.mDoQtGammaSelection==2){
   //     Float_t qtMaxPtDep = mParams.mQtPtMax*photon->GetPhotonPt();
   //     if (qtMaxPtDep > mParams.mQtMax)
   //       qtMaxPtDep      = mParams.mQtMax;
   //     if(photon->GetArmenterosQt()>qtMaxPtDep){
   //       return false;
   //     }
   //   }
   // }
   return true;
}

///________________________________________________________________________
bool MpdV0Maker::AsymmetryCut(float asym, float pt) const
{
   // Cut on Energy Asymmetry

   // for(Int_t ii=0;ii<2;ii++){

   //   AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

   //   if(fDoPhotonPDependentAsymCut){
   //     float trackNegAsy=0;
   //     if (photon->GetPhotonP()!=0.){
   //         trackNegAsy= track->P()/photon->GetPhotonP();
   //     }

   //     if( trackNegAsy > fFAsymmetryCut->Eval(photon->GetPhotonP()) || trackNegAsy
   //     < 1.-fFAsymmetryCut->Eval(photon->GetPhotonP()) ){
   //       return kFALSE;
   //     }

   //   } else {
   //     if( track->P() > fMinPPhotonAsymmetryCut ){
   //       float trackNegAsy=0;
   //       if (photon->GetPhotonP()!=0.){
   //         trackNegAsy= track->P()/photon->GetPhotonP();
   //       }

   //       if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
   //         return kFALSE;
   //       }
   //     }
   //   }

   // }
   return true;
}
///________________________________________________________________________
float MpdV0Maker::CosPsiPair(TVector3 &p1, TVector3 &p2) const
{

   // float p1[3] = {tr1->GetPx(),tr1->GetPy(),tr1->GetPz()};
   // float p2[3] = {tr2->GetPx(),tr2->GetPy(),tr2->GetPz()};
   // float u[3] = {p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2]};
   TVector3 u = p1 + p2;

   // float normp1 = sqrt( (p1[0]*p1[0]) + (p1[1]*p1[1]) + (p1[2]*p1[2]) );
   // float normp2 = sqrt( (p2[0]*p2[0]) + (p2[1]*p2[1]) + (p2[2]*p2[2]) );
   // float normu  = sqrt( (u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]) );

   // for(int i=3; i--;){
   //   p1[i] /= normp1;
   //   p2[i] /= normp2;
   //   u[i] /= normu;
   // }

   TVector3 v = p1.Cross(p2);
   TVector3 w = u.Cross(v);
   TVector3 z(0, 0, 1.);
   TVector3 wc = u.Cross(z);
   return cos(wc.Angle(w));
}

long int MpdV0Maker::FindCommonParent(long int prim1, long int prim2, TClonesArray *MCTracks)
{
   // Looks through parents and finds if there was commont pi0 among ancestors

   if (!MCTracks->GetEntriesFast()) return -1; // can not say anything

   while (prim1 != -1) {
      long int pr2 = prim2;

      while (pr2 != -1) {
         if (prim1 == pr2) {
            return prim1;
         }
         pr2 = (static_cast<MpdMCTrack *>(MCTracks->At(pr2)))->GetMotherId();
      }
      prim1 = (static_cast<MpdMCTrack *>(MCTracks->At(prim1)))->GetMotherId();
   }
   return -1;
}
