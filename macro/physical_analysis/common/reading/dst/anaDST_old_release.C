/* Macro reads DST file produced by macro reco.C */

void anaDST(TString inDir = "/opt/exp_soft/mpd/data4mpd/dst/09GeV/",TString inFile = "auau_09gev_0_3fm_0.root")
{
 TStopwatch timer;
 timer.Start();
 
  /* Load basic libraries */
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE,kFALSE); // only reco lib

  TString outDir = "./dst1/";
  //  TStting outFile = "auau_7GeV_0_3fm_h0.root";
  outFile = outDir + inFile;
  fFile = new TFile(outFile, "RECREATE");
  inFile = inDir+inFile;
  //  TString inFile = "./dst/auau_09gev_0_3fm_0.root";
  TFile fileDST(inFile);

  Double_t Pt, Pt_MC, Rap,Rap_MC, fProbP, fProbPi, fProbK, fProbE;
  Int_t ID, Npart, Np, Npi, Nk;
  Int_t NpMC = 0, NpiMC = 0, NkMC = 0;
  Int_t N_ToF = 0, N_ToF_F = 0;
  Double_t fProbCut = 0.8;

  enum DetectorId {kSTS, kTPC, kTOF, kETOF, kFFD, kECT, kECAL, kNDET, kCPC, kBBC, kZDC, kFSA};

  fFile->cd();

  TH1F *hPt = new TH1F("Pt","",100,0.0,3.0);
   hPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hPt->GetYaxis()->SetTitleOffset(1.2);
   hPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hRap = new TH1F("Rap","",100,-2.0,2.0);
   hRap->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hRap->GetYaxis()->SetTitleOffset(1.2);
   hRap->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpart = new TH1F("hpNpart","",50,0.0,200.0);
   hpNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpNpart->GetYaxis()->SetTitleOffset(1.2);
   hpNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPt = new TH1F("hpPt","",100,0.0,3.0);
   hpPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpPt->GetYaxis()->SetTitleOffset(1.2);
   hpPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpRap = new TH1F("hpRap","",100,-2.0,2.0);
   hpRap->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpRap->GetYaxis()->SetTitleOffset(1.2);
   hpRap->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpartMC = new TH1F("hpNpartMC","",50,0.0,200.0);
   hpNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPtMC = new TH1F("hpPtMC","",100,0.0,3.0);
   hpPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpRapMC = new TH1F("hpRapMC","",100,-2.0,2.0);
   hpRapMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpRapMC->GetYaxis()->SetTitleOffset(1.2);
   hpRapMC->GetXaxis()->SetTitle("#eta");

  TH1F *hpiNpart = new TH1F("hpiNpart","",100,0.0,500.0);
   hpiNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpiNpart->GetYaxis()->SetTitleOffset(1.2);
   hpiNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPt = new TH1F("hpiPt","",100,0.0,3.0);
   hpiPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpiPt->GetYaxis()->SetTitleOffset(1.2);
   hpiPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiRap = new TH1F("hpiRap","",100,-2.0,2.0);
   hpiRap->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpiRap->GetYaxis()->SetTitleOffset(1.2);
   hpiRap->GetXaxis()->SetTitle("#eta");

  TH1F *hpiNpartMC = new TH1F("hpiNpartMC","",100,0.0,500.0);
   hpiNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpiNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpiNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPtMC = new TH1F("hpiPtMC","",100,0.0,3.0);
   hpiPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpiPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpiPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiRapMC = new TH1F("hpiRapMC","",100,-2.0,2.0);
   hpiRapMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpiRapMC->GetYaxis()->SetTitleOffset(1.2);
   hpiRapMC->GetXaxis()->SetTitle("#eta");

  TH1F *hkNpart = new TH1F("hkNpart","",50,0.0,100.0);
   hkNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hkNpart->GetYaxis()->SetTitleOffset(1.2);
   hkNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPt = new TH1F("hkPt","",100,0.0,3.0);
   hkPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hkPt->GetYaxis()->SetTitleOffset(1.2);
   hkPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkRap = new TH1F("hkRap","",100,-2.0,2.0);
   hkRap->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hkRap->GetYaxis()->SetTitleOffset(1.2);
   hkRap->GetXaxis()->SetTitle("#eta");

   TH1F *hkNpartMC = new TH1F("hkNpartMC","", 50,0.0,100.0);
   hkNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hkNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hkNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPtMC = new TH1F("hkPtMC","",100,0.0,3.0);
   hkPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hkPtMC->GetYaxis()->SetTitleOffset(1.2);
   hkPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkRapMC = new TH1F("hkRapMC","",100,-2.0,2.0);
   hkRapMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hkRapMC->GetYaxis()->SetTitleOffset(1.2);
   hkRapMC->GetXaxis()->SetTitle("#eta");

  TTree *simTree = (TTree*) fileDST.Get("cbmsim");

  TClonesArray *mpdEvents = (TClonesArray*) fileDST.FindObjectAny("MpdEvent");
  simTree->SetBranchAddress("MpdEvent", &mpdEvents);
  TBranch *dstBranch = simTree->GetBranch("MpdEvent");

  TClonesArray *fT1 = (TClonesArray) fileDST.FindObjectAny("MCTrack");
  simTree->SetBranchAddress("MCTrack",&fT1);
  TBranch *MCBranch = simTree->GetBranch("MCTrack");

  Int_t events = simTree->GetEntries();
  cout << " Number of events in DST file = " << events << endl;
  //   events = 1;
  cout << "     Number of events in loop = " << events << endl;
  Npart = 0;
  /* Events loop */
  for (Int_t i = 0; i < events-1; i++) 
  {
    Np = 0;
    Npi = 0;
    Nk = 0;
    NpMC = 0;
    NpiMC = 0;
    NkMC = 0;

    cout << " *** Event # " << i+1 << endl;
    dstBranch->GetEntry(i);
    MCBranch->GetEntry(i);
    
    MpdEvent *event = (MpdEvent *) mpdEvents->UncheckedAt(i);
    
    //    event->Dump();
    
    TClonesArray *mpdTracks = event->GetGlobalTracks();
    Int_t fNtracks = mpdTracks->GetEntriesFast();
    //    FairMCTrack* mctrack;
    Int_t fMotherID, fPDGID;
    //    cout << " Number of tracks = " << fNtracks << endl;
    for (Int_t k = 0; k < fNtracks; k++)
    {
      MpdTrack *track = mpdTracks->UncheckedAt(k);

      //      track->Dump();
      ID = track->GetID();
      FairMCTrack *mctrack = fT1->UncheckedAt(ID);

      //      hID->Fill(ID);
      fMotherID = mctrack->GetMotherId();

      fPDGID = mctrack->GetPdgCode();

      if (fMotherID < 0 && track->GetTofFlag())
      	{
	  Npart++;

	  Pt = track->GetPt();
	  Pt = TMath::Abs(Pt);
	  hPt->Fill(Pt);
	  Rap = track->GetEta();
	  hRap->Fill(Rap);

	  fProbP = track->GetPidProbProton();
	  fProbPi = track->GetPidProbPion();
	  fProbK = track->GetPidProbKaon();
	  fProbE = track->GetPidProbElectron();

	  if ( fProbP > fProbCut )
	    {
	      Np++;
	      hpPt->Fill(Pt);
	      hpRap->Fill(Rap);
	    }
	  if ( TMath::Abs(fProbPi) > fProbCut )
	    {
	      Npi++;
	      hpiPt->Fill(Pt);
	      hpiRap->Fill(Rap);
	    }

	  if ( TMath::Abs(fProbK) > fProbCut )
	    {
	      Nk++;
	      hkPt->Fill(Pt);
	      hkRap->Fill(Rap);
	    }
	}
    } // track loop
    hpNpart->Fill(Np);
    hpiNpart->Fill(Npi);
    hkNpart->Fill(Nk);

    Int_t mcNtracks = fT1->GetEntries();

    for (Int_t j = 0; j < mcNtracks; j++)
    {
      TVector3 Pxyz;
      mctrack = (FairMCTrack*)fT1->At(j);

      //      FairMCTrack *mctrack = fT1->UncheckedAt(j);

      fMotherID = mctrack->GetMotherId();

      fPDGID = mctrack->GetPdgCode();
	
      if (fMotherID < 0 && mctrack->GetNPoints(kTOF))
      	{
	  Pt_MC = mctrack->GetPt();
          Pt_MC = TMath::Abs(Pt_MC);
          if (Pt_MC > 0)
            {
              Pxyz.SetXYZ(mctrack->GetPx(),mctrack->GetPy(),mctrack->GetPz());
              Rap_MC = Pxyz.PseudoRapidity();
            }
	  if (fPDGID == 2212)
	    {
	      if (TMath::Abs(Pt_MC) > 0.2 && TMath::Abs(Rap_MC) < 1.2) NpMC++;
	      if (TMath::Abs(Rap_MC) < 1.2) hpPtMC->Fill(Pt_MC);
	      hpRapMC->Fill(Rap_MC);
	    }
	  if (TMath::Abs(fPDGID) == 211)
	    {
	      if (TMath::Abs(Pt_MC) > 0.2 && TMath::Abs(Rap_MC) < 1.2) NpiMC++;
	      if (TMath::Abs(Rap_MC) < 1.2) hpiPtMC->Fill(Pt_MC);
	      hpiRapMC->Fill(Rap_MC);
	    }
	  if (TMath::Abs(fPDGID) == 321)
	    {
	      if (TMath::Abs(Pt_MC) > 0.2 && TMath::Abs(Rap_MC) < 1.2) NkMC++;
	      if (TMath::Abs(Rap_MC) < 1.2) hkPtMC->Fill(Pt_MC);
	      hkRapMC->Fill(Rap_MC);
	    }
	}
    } // MC track loop
    hpNpartMC->Fill(NpMC);
    hpiNpartMC->Fill(NpiMC);
    hkNpartMC->Fill(NkMC);
  } // event loop 

  fFile->Write();

  timer.Print();

  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  //  exit(0);
}
