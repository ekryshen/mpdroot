/*
 Standart methods of ellipic flow.
 Author: V.Voronyuk vadimv@jinr.ru
 
 v2{true} is true v2 if Psi_RP=0
 v2{EP} is standart method with full event plane resolution by equal number random sub events 
 v2{2} is cumulant 2 method
 v2{eta} is method with two sub events with eta>0 and eta<0  

 Usage: root v2std.C"(energy, \"dst_filename*.root\")"
 out:   v2std-[energe]Gev-[MC/RECO].root
 
 result:  TGraphErrors *grv2true,*grv2EP,*grv2cum2,*grv2eta (It is needed "APT" options to draw them)
 
 Line "#define MC" defines MC
 Comment it and you obtain RECO 
*/

#define MC

// #define debug

#include <TDatabasePDG.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using TMath::Cos;
using TMath::Sin;
using TMath::ATan2;
using TMath::Sqrt;

using std::cout;
using std::endl;

void v2std(int eBeam=9, TString inFilesChain="mpddst.root")
{
  gStyle->SetFillColor(0);        // white fill color
  gStyle->SetCanvasBorderMode(0); // turn off canvas borders
  gStyle->SetPadBorderMode(0);    // no pad border
  gStyle->SetFrameBorderMode(0);  // no frame border
  
  /* Load basic libraries */
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE,kFALSE); // only reco libs

  const Float_t Cut_PtMin = 0.2;
  const Float_t Cut_PtMax = 1.5;
  const Float_t Cut_Eta   = 1.; // pseudo rapidity
  const Float_t Cut_Prob  = 0.8; // probability
  
  TChain *DSTTree = new TChain("mpdsim");
  DSTTree->Add(inFilesChain.Data());

  TString outFileName = "v2std-";
  outFileName += eBeam;

#ifdef MC
  outFileName += "Gev-MC.root";
#else
  outFileName += "Gev-RECO.root";
#endif

  TFile *outFile = new TFile(outFileName.Data(), "RECREATE");

// ----------------------------------------------->

  TString hName,hNameAdd,hTitle;

#ifdef MC
  hNameAdd="_MC";
#else
  hNameAdd="_RECO";
#endif

  hTitle="#sqrt{s_{nn}}=";
  hTitle +=eBeam;
  
  // true value if PhiRP=0
  
  outFile->mkdir("true_method");
  outFile->cd("true_method");
  
  hName = "v2true_";
  hName += eBeam;
  hName += hNameAdd;
  
  TH1F *hv2true = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hv2true->GetYaxis()->SetTitle("count");
  hv2true->GetXaxis()->SetTitle("v^{true}_2");
  
  hName = "v2true_Pt_";
  hName += eBeam;
  hName += hNameAdd;
  
  TProfile *hp_v2true_Pt = new TProfile(hName.Data(),hTitle.Data(),30,Cut_PtMin,Cut_PtMax,-1.,1.);
  hp_v2true_Pt->GetYaxis()->SetTitle("v^{true}_{2}");
  hp_v2true_Pt->GetXaxis()->SetTitle("p_t");
  
  hName = "v2true_Eta_";
  hName += eBeam;
  hName += hNameAdd;
  
  TProfile *hp_v2true_Eta = new TProfile(hName.Data(),hTitle.Data(),30,-Cut_Eta,Cut_Eta,-1.,1.);
  hp_v2true_Eta->GetYaxis()->SetTitle("v^{true}_{2}");
  hp_v2true_Eta->GetXaxis()->SetTitle("#eta");
  
  outFile->cd();
  
  /* standart method */
  
  outFile->mkdir("std_method");
  outFile->cd("std_method");
  
  hName = "v2obs_";
  hName += eBeam;
  hName += hNameAdd;
  
  TH1F *hv2obs = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hv2obs->GetYaxis()->SetTitle("count");
  hv2obs->GetXaxis()->SetTitle("v^{obs}_2");
  
  hName = "PhiEP_";
  hName += eBeam;
  hName += hNameAdd;
  
  TH1F *hPhiEP = new TH1F(hName.Data(),hTitle.Data(),100,-TMath::Pi()/2.,TMath::Pi()/2.);
  hPhiEP->GetYaxis()->SetTitle("count");
  hPhiEP->GetXaxis()->SetTitle("#Psi_{2,EP}");
  
  hName = "SubResEP_";
  hName += eBeam;
  hName += hNameAdd;
  
  TH1F *hSubResEP = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hSubResEP->GetYaxis()->SetTitle("count");
  hSubResEP->GetXaxis()->SetTitle("cos[2(#Psi^a-#Psi^b)]");
  
  outFile->cd();
  
  // cummulant method
  
  outFile->mkdir("cum_method");
  outFile->cd("cum_method");
  
  hName = "Cum2_";
  hName += eBeam;
  hName += hNameAdd;
  TH1F *hcum2 = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hcum2->GetYaxis()->SetTitle("count");
  hcum2->GetXaxis()->SetTitle("cum2");
  
  outFile->cd();

  // eta method
  
  outFile->mkdir("eta_method");
  outFile->cd("eta_method");
  
  hName = "v2obsEta_";
  hName += eBeam;
  hName += hNameAdd;
  TH1F *hv2obsEta = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hv2obsEta->GetYaxis()->SetTitle("count");
  hv2obsEta->GetXaxis()->SetTitle("v_{2}^{obs}");
  
  hName = "SubResEPEta_";
  hName += eBeam;
  hName += hNameAdd;
  TH1F *hSubResEPEta = new TH1F(hName.Data(),hTitle.Data(),500,-1.,1.);
  hSubResEPEta->GetYaxis()->SetTitle("count");
  hSubResEPEta->GetXaxis()->SetTitle("cos[2(#phi_{EP,#eta_+}-#phi_{EP,#eta_-})]");

  outFile->cd();
  
// ----------------------------------------------->

  MpdEvent *event;
  DSTTree->SetBranchAddress("MPDEvent.", &event);
  TClonesArray *MCTracks;
  DSTTree->SetBranchAddress("MCTrack", &MCTracks);

  Int_t nevents=DSTTree->GetEntries();
  
  cout << "Number of events in DST: " << nevents << endl;
  
  Int_t n=0;
  
#ifdef debug
  nevents=10;
  cout << "Number of events to read: " << nevents << endl;
#endif

  TClonesArray *MpdTracks;
    
  //Loop over events,  nevent --- the number of events 
  for (Int_t j=0; j < nevents; ++j)
  {
    DSTTree->GetEntry(j);
    
#ifdef MC
    Int_t ntracks = MCTracks->GetEntries();
#else
    MpdTracks = event->GetGlobalTracks();
    Int_t ntracks = MpdTracks->GetEntriesFast();
#endif
    
#ifdef debug
    cout << "Number of tracks= " << ntracks << endl;
#endif

    Double_t Pt  = 0.;
    Double_t Eta = 0.;
    Int_t N_PDG  = 0;
    
    FairMCTrack *mctrack;
    MpdTrack *mpdtrack;
    Float_t prob;
    TLorentzVector lvP;
    
    UInt_t nseltracks=0;
    TLorentzVector *sellv = new TLorentzVector[ntracks];
  
    // Loop over tracks in the current event    
    for (UInt_t i=0; i < ntracks; ++i)
    {
      
#ifdef MC
      mctrack = (FairMCTrack*) MCTracks->At(i);
      mctrack->Get4Momentum(lvP);
      N_PDG = mctrack->GetPdgCode();
#else    
      mpdtrack = (MpdTrack*) MpdTracks->UncheckedAt(i);
      if ( !mpdtrack->GetTofFlag() ) continue;  // no tof identification
      prob = mpdtrack->GetPidProbPion();
      if (prob> Cut_Prob) N_PDG=211; // pi+    
      if (prob<-Cut_Prob) N_PDG=-211; // pi-
      prob = mpdtrack->GetPidProbKaon();
      if (prob>Cut_Prob) N_PDG=321; // K+
      if (prob<-Cut_Prob) N_PDG=-321; // K-
      if (N_PDG==0) continue; // other not needed    
      lvP.SetXYZM( mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(),
                   fPDG->GetParticle(N_PDG)->Mass() );
#endif
      
      if (N_PDG!=211 && N_PDG!=-211) continue; // only charged pi mesons (FIXME)
      
      Pt= lvP.Perp();
      if (Pt<Cut_PtMin || Pt>Cut_PtMax) continue; // cut Pt
      Eta=lvP.Eta();
      if (TMath::Abs(Eta)>Cut_Eta) continue;  // cut Eta

#ifdef debug
      cout << "pdg= " << N_PDG << " Pt= " << Pt << " Eta= " << Eta << endl;
#endif

#ifndef MC    
      /* select primary tracks */
      mctrack = (FairMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
#endif
      if (mctrack->GetMotherId() > -1) continue;  // only primary track (FIXME)

      /* store selected */
      sellv[nseltracks] = lvP;
      nseltracks++;
    } // tracks` loop

    if (nseltracks==0) continue; // goto next event
 

    cout << "Number of selected tracks = " << nseltracks << endl;
   
/* ----- start analisys in given event ----- */
  
    /* true value if Phi_RP=0 */
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      Double_t v2true=TMath::Cos(2.*sellv[i].Phi());
      Double_t pt=sellv[i].Perp();
      Double_t eta=sellv[i].Eta();
      
      hv2true->Fill(v2true);
      hp_v2true_Pt->Fill(pt,v2true,1.);
      hp_v2true_Eta->Fill(TMath::Abs(eta),v2true,1.);
      hp_v2true_Eta->Fill(-TMath::Abs(eta),v2true,1.);
    }
    // v2=hv2true->GetMean()
    
    
    /* event plane angle - standart method with weight=pt */
    Double_t Qcos=0., Qsin=0.; // 
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      Double_t pt=sellv[i].Perp();
      Double_t phi=sellv[i].Phi();
      Qcos += pt*TMath::Cos(2.*phi);
      Qsin += pt*TMath::Sin(2.*phi);
    }
    Double_t PhiEP = 0.5*TMath::ATan2(Qsin,Qcos);
#ifdef debug
    cout << "Event plane angle= " << PhiEP << endl;
#endif
    hPhiEP->Fill(PhiEP);


    /* event plane - eta method */
    /* B. I. ABELEV et al. PHYSICAL REVIEW C 77, 054901 (2008) */
    Double_t QcosEtap=0., QsinEtap=0.; // eta>0.075
    Double_t QcosEtam=0., QsinEtam=0.; // eta<-0.075
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      Double_t pt=sellv[i].Perp();
      Double_t phi=sellv[i].Phi();
      Double_t eta=sellv[i].Eta();
      if (eta>0.075)
      {
        QcosEtap += pt*TMath::Cos(2.*phi);
        QsinEtap += pt*TMath::Sin(2.*phi);
      }
      if (eta<-0.075)
      {
        QcosEtam += pt*TMath::Cos(2.*phi);
        QsinEtam += pt*TMath::Sin(2.*phi);
      }
    }
    Double_t PhiEPEtap = 0.5*TMath::ATan2(QsinEtap,QcosEtap);
    Double_t PhiEPEtam = 0.5*TMath::ATan2(QsinEtam,QcosEtam);
    hSubResEPEta->Fill(TMath::Cos(2.*PhiEPEtap-2.*PhiEPEtam));
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      Double_t phi=sellv[i].Phi();
      Double_t eta=sellv[i].Eta();
      if (eta>0.075)
      {
        hv2obsEta->Fill(TMath::Cos(2.*phi-2.*PhiEPEtam));
      }
      if (eta<-0.075)
      {
        hv2obsEta->Fill(TMath::Cos(2.*phi-2.*PhiEPEtap));
      }
    }
    // v2{eta}= hv2obsEta->GetMean()/TMath::Sqrt(hSubResEPEta->GetMean())
  
    
    /* 2 random equal number subevents */
    Int_t subN=0;
    Bool_t *subsel = new Bool_t[nseltracks];
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      subsel[i]=kFALSE;
    }    
    for (UInt_t i=0; i<nseltracks/2; ++i)
    {
      UInt_t isel= (UInt_t) (gRandom->Rndm()*((Float_t) nseltracks));
      subsel[isel]=kTRUE;
      subN++;
    }
#ifdef debug
    cout << "Sub N1= " << subN << " SubN2=" << nseltracks-subN << endl;
#endif
    // subsel[] subevent flag
    
    
    /* cumulant2 method */
    for (UInt_t i=1; i<nseltracks; ++i)
    {
      if (!subsel[i]) continue; // select 1 subevent
      Double_t subphi1=sellv[i].Phi();
      for (UInt_t k=0; k<i-1; ++k)
      {
        if (subsel[k]) continue; // select 2 subevent
        Double_t subphi2=sellv[k].Phi();
        hcum2->Fill(TMath::Cos(2.*subphi1-2.*subphi2));
      }
    }
    // v2{2}=sqrt(hcum2)
    
    
    /* event plane resolution by 2 random subevents */
    // events must be about the same multiplicity (centrality)
    Double_t Qcos1=0.,Qsin1=0.,Qcos2=0.,Qsin2=0.;
    for (UInt_t i=1; i<nseltracks; ++i)
    {
      Double_t pt=sellv[i].Perp();
      Double_t phi=sellv[i].Phi();
      if (!subsel[i]) // select 1 subevent
      {
        Qcos1 += pt*TMath::Cos(2.*phi);
        Qsin1 += pt*TMath::Sin(2.*phi);
      }
      if (subsel[i]) // select 2 subevent
      {
        Qcos2 += pt*TMath::Cos(2.*phi);
        Qsin2 += pt*TMath::Sin(2.*phi);
      }
    }
    Double_t SubPhiEP1 = 0.5*TMath::ATan2(Qsin1,Qcos1);
    Double_t SubPhiEP2 = 0.5*TMath::ATan2(Qsin2,Qcos2);
    Double_t SubRes=TMath::Cos(2.*SubPhiEP1-2.*SubPhiEP2);
    hSubResEP->Fill(SubRes);
    
    for (UInt_t i=0; i<nseltracks; ++i)
    {
      Double_t pt=sellv[i].Perp();
      Double_t phi=sellv[i].Phi();
      Double_t Qs=Qsin-pt*TMath::Sin(2.*phi);
      Double_t Qc=Qcos-pt*TMath::Cos(2.*phi);
      Double_t phiep=0.5*TMath::ATan2(Qs,Qc);
      
      hv2obs->Fill(TMath::Cos(2.*phi-2.*phiep));
    }
    // v2obs=hv2obs->GetMean
    
/* ----- end analisys in given event ----- */
 
    delete [] sellv;
    delete [] subsel;
  } // end loop over events

  Double_t SubEvRes=TMath::Sqrt(hSubResEP->GetMean());
  Double_t SubEvRes1=TMath::Sqrt(hSubResEP->GetMean()+0.01);
  Double_t ResEP=ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes));
  Double_t ResEP1=ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes1));
  Double_t ResEPError= (ResEP1-ResEP)/0.01*hSubResEP->GetMeanError();
  
  // RESULT:
  // v2{true}=hv2true->GetMean() (if Psi_RP=0!)
  // v2{EP}=v2obs/ResEP
  // v2{2}=sqrt(hcum2->GetMean)
  // v2{eta}= hv2obsEta->GetMean()/TMath::Sqrt(hSubResEPEta->GetMean())

  TGraphErrors *grv2true = new TGraphErrors();
  grv2true->SetTitle("v_{2}^{true} (if #Psi_{RP}=0)");
  grv2true->SetMinimum(0.);
  grv2true->SetMaximum(6.);
  grv2true->SetMarkerColor(4);
  grv2true->SetMarkerStyle(30);
  grv2true->SetPoint(0,eBeam,100.*hv2true->GetMean());
  grv2true->SetPointError(0,0.,100.*hv2true->GetMeanError());
  grv2true->GetXaxis()->SetTitle("#sqrt{S_{NN}}");
  grv2true->GetYaxis()->SetTitle("v_{2} , %");
  grv2true->Draw("APT");
  grv2true->GetXaxis()->SetLimits(3.,15.);
  hName = "grv2true_";
  hName += eBeam;
  hName += hNameAdd;
  grv2true->Write(hName.Data());
  
  TGraphErrors *grv2EP = new TGraphErrors();
  grv2EP->SetTitle("v_{2}{EP}");
  grv2EP->SetMinimum(0.);
  grv2EP->SetMaximum(6.);
  grv2EP->SetMarkerColor(4);
  grv2EP->SetMarkerStyle(30);
  grv2EP->SetPoint(0,eBeam,100.*hv2obs->GetMean()/ResEP);
  Double_t v2EPError=TMath::Sqrt(1/ResEP**2*hv2obs->GetMeanError()**2+
                     (hv2obs->GetMean()/ResEP)**2*ResEPError**2);
  grv2EP->SetPointError(0,0.,100.*v2EPError);
  grv2EP->GetXaxis()->SetTitle("#sqrt{S_{NN}}");
  grv2EP->GetYaxis()->SetTitle("v_{2} , %");
  grv2EP->Draw("APT");
  grv2EP->GetXaxis()->SetLimits(3.,15.);
  hName = "grv2EP_";
  hName += eBeam;
  hName += hNameAdd;
  grv2EP->Write(hName.Data());
  
  if (hcum2->GetMean()>0)
  {
    TGraphErrors *grv2cum2 = new TGraphErrors();
    grv2cum2->SetTitle("v_{2}{2}");
    grv2cum2->SetMinimum(0.);
    grv2cum2->SetMaximum(6.);
    grv2cum2->SetMarkerColor(4);
    grv2cum2->SetMarkerStyle(30);
    grv2cum2->SetPoint(0,eBeam,100.*TMath::Sqrt(hcum2->GetMean()));
    Double_t v2cum2Error=hcum2->GetMeanError()/(2.*TMath::Sqrt(hcum2->GetMean()));
    grv2cum2->SetPointError(0,0.,100.*v2cum2Error);
    grv2cum2->GetXaxis()->SetTitle("#sqrt{S_{NN}}");
    grv2cum2->GetYaxis()->SetTitle("v_{2} , %");
    grv2cum2->Draw("APT");
    grv2cum2->GetXaxis()->SetLimits(3.,15.);
    hName = "grv2cum2_";
    hName += eBeam;
    hName += hNameAdd;
    grv2cum2->Write(hName.Data());
  }
  
  if (hSubResEPEta->GetMean()>0)
  {
    TGraphErrors *grv2eta = new TGraphErrors();
    grv2eta->SetTitle("v_{2}{#eta}");
    grv2eta->SetMinimum(0.);
    grv2eta->SetMaximum(6.);
    grv2eta->SetMarkerColor(4);
    grv2eta->SetMarkerStyle(30);
    grv2eta->SetPoint(0,eBeam,100.*hv2obsEta->GetMean()/TMath::Sqrt(hSubResEPEta->GetMean()));
    Double_t v2etaError=TMath::Sqrt(1./hSubResEPEta->GetMean()*hv2obsEta->GetMeanError()**2+
          hv2obsEta->GetMean()**2/(4.*hSubResEPEta->GetMean()**3)*hSubResEPEta->GetMeanError()**2);
    grv2eta->SetPointError(0,0.,100.*v2etaError);
    grv2eta->GetXaxis()->SetTitle("#sqrt{S_{NN}}");
    grv2eta->GetYaxis()->SetTitle("v_{2} , %");
    grv2eta->Draw("APT");
    grv2eta->GetXaxis()->SetLimits(3.,15.);
    hName = "grv2eta_";
    hName += eBeam;
    hName += hNameAdd;
    grv2eta->Write(hName.Data());
  }
  
  
  outFile->Write();

  cout << endl;
  cout << "sqrt(s_{nn})= " << eBeam << " Gev"<< endl; 
  cout << "v2true= " << hv2true->GetMean() << "  only if PhiRP=0!" << endl;
  cout << "v2obs=  " << hv2obs->GetMean() << endl;
  cout << "v2{EP}= " << hv2obs->GetMean()/ResEP << endl;
  cout << " ResEP= " << ResEP << endl;
  if (hcum2->GetMean()>0)
  {
    cout << "v2{2}=  " << TMath::Sqrt(hcum2->GetMean()) << endl;
  }
  else
  {
    cout << "Cum2=" << hcum2->GetMean() << "   <0" << endl;
  }
  if (hSubResEPEta->GetMean()>0)
  {
    cout << "v2{eta}=" << hv2obsEta->GetMean()/TMath::Sqrt(hSubResEPEta->GetMean()) << endl;
  }
  else
  {
    cout << "Sub Event Res {eta}=" << hSubResEPEta->GetMean() << "   <0" << endl;
  }

}


/* ResEventPlane and Chi are taken from AliRoot. Thanks a lot. */
Double_t ResEventPlane(Double_t chi)
{
  // plane resolution as function of chi
  Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
  Double_t arg = chi * chi / 4.;
  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res ;
}
/* ResEventPlane and Chi are taken from AliRoot. Thanks a lot. */
Double_t Chi(Double_t res)
{
  // chi from the event plane resolution
 
  Double_t chi   = 2.0;
  Double_t delta = 1.0;
  for(int i = 0; i < 15; i++)
  {
   if(ResEventPlane(chi) < res) { chi = chi + delta ; }
   else                         { chi = chi - delta ; }
   delta = delta / 2.;
  }

  return chi ;
}


