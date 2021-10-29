/* Macro reads DST file produced by macro reco.C */

void anaDST(TString inFile = "dst_auau_04gev_0_3fm_0.root")
{
 TStopwatch timer;
 timer.Start();
 
  /* Load basic libraries */
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(kTRUE,kFALSE); // only reco lib

  TString outDir = "a";

  outFile = outDir + inFile;
  fFile = new TFile(outFile, "RECREATE");

  TFile fileDST(inFile);

  Float_t E, E_MC, P, Pt, Pz, Pt_MC, Pz, Pz_MC, Px, Px_MC, Py, Py_MC, Mt, Mt_MC, Eta,Eta_MC;
  Float_t M2, dEdx;
  Float_t ProbP, ProbPi, ProbK, ProbE;
  Float_t ProbP_dEdx, ProbPi_dEdx, ProbK_dEdx, ProbE_dEdx;
  Int_t ID, Npart=0, Np=0, Nap=0, Npi=0, Nk=0;
  Int_t NpMC = 0, NpiMC = 0, NkMC = 0;
  Int_t N_ToF = 0, N_ToF_F = 0;
  Float_t ProbCut = 0.8, ProbCutdEdx=0.5;
  Int_t MotherID, PDGID;

  enum DetectorId {kSTS, kTPC, kTOF, kETOF, kFFD, kECT, kECAL, kNDET, kCPC, kBBC, kZDC, kFSA};

  const Double_t kProtonMass = 0.938271998;
  const Double_t kPionMass = 0.13957018;
  const Double_t kKaonMass = 0.493677;

  fFile->cd();

  TH1F *hPt = new TH1F("Pt","",100,0.0,3.0);
   hPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hPt->GetYaxis()->SetTitleOffset(1.2);
   hPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hEta = new TH1F("Eta","",100,-2.0,2.0);
   hEta->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hEta->GetYaxis()->SetTitleOffset(1.2);
   hEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpart = new TH1F("hpNpart","",50,0.0,200.0);
   hpNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpNpart->GetYaxis()->SetTitleOffset(1.2);
   hpNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPt = new TH1F("hpPt","",100,0.0,3.0);
   hpPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpPt->GetYaxis()->SetTitleOffset(1.2);
   hpPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpEta = new TH1F("hpEta","",100,-2.0,2.0);
   hpEta->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpEta->GetYaxis()->SetTitleOffset(1.2);
   hpEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpartMC = new TH1F("hpNpartMC","",50,0.0,200.0);
   hpNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPtMC = new TH1F("hpPtMC","",100,0.0,3.0);
   hpPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpEtaMC = new TH1F("hpEtaMC","",100,-2.0,2.0);
   hpEtaMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hpEtaMC->GetXaxis()->SetTitle("#eta");

  TH1F *hpiNpart = new TH1F("hpiNpart","",100,0.0,500.0);
   hpiNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpiNpart->GetYaxis()->SetTitleOffset(1.2);
   hpiNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPt = new TH1F("hpiPt","",100,0.0,3.0);
   hpiPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpiPt->GetYaxis()->SetTitleOffset(1.2);
   hpiPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiEta = new TH1F("hpiEta","",100,-2.0,2.0);
   hpiEta->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpiEta->GetYaxis()->SetTitleOffset(1.2);
   hpiEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpiNpartMC = new TH1F("hpiNpartMC","",100,0.0,500.0);
   hpiNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hpiNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpiNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPtMC = new TH1F("hpiPtMC","",100,0.0,3.0);
   hpiPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hpiPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpiPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiEtaMC = new TH1F("hpiEtaMC","",100,-2.0,2.0);
   hpiEtaMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hpiEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hpiEtaMC->GetXaxis()->SetTitle("#eta");

  TH1F *hkNpart = new TH1F("hkNpart","",50,0.0,100.0);
   hkNpart->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hkNpart->GetYaxis()->SetTitleOffset(1.2);
   hkNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPt = new TH1F("hkPt","",100,0.0,3.0);
   hkPt->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hkPt->GetYaxis()->SetTitleOffset(1.2);
   hkPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkMtMC = new TH1F("hkMtMC","",100,0.0,3.0);
   hkMtMC->GetYaxis()->SetTitle("dN/M_{t}dM_{t}"); // (1/10 MeV/c)");
   hkMtMC->GetYaxis()->SetTitleOffset(1.2);
   hkMtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkEta = new TH1F("hkEta","",100,-2.0,2.0);
   hkEta->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hkEta->GetYaxis()->SetTitleOffset(1.2);
   hkEta->GetXaxis()->SetTitle("#eta");

   TH1F *hkNpartMC = new TH1F("hkNpartMC","", 50,0.0,100.0);
   hkNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); // (1/5)");
   hkNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hkNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPtMC = new TH1F("hkPtMC","",100,0.0,3.0);
   hkPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); // (1/10 MeV/c)");
   hkPtMC->GetYaxis()->SetTitleOffset(1.2);
   hkPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkEtaMC = new TH1F("hkEtaMC","",100,-2.0,2.0);
   hkEtaMC->GetYaxis()->SetTitle("dN/d#eta"); //, (1/0.05)");
   hkEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hkEtaMC->GetXaxis()->SetTitle("#eta");

  TTree *simTree = (TTree*) fileDST.Get("mpdsim");

  //  TClonesArray *mpdEvents = (TClonesArray*) fileDST.FindObjectAny("MPDEvent.");
  MpdEvent *event;
  simTree->SetBranchAddress("MPDEvent.", &event);
  TBranch *dstBranch = simTree->GetBranch("MPDEvent.");

  TClonesArray *fT1 = (TClonesArray) fileDST.FindObjectAny("MCTrack");
  simTree->SetBranchAddress("MCTrack",&fT1);
  TBranch *MCBranch = simTree->GetBranch("MCTrack");

  Int_t events = simTree->GetEntries();
  cout << " Number of events in DST file = " << events << endl;
  //   events = 1;
  cout << "     Number of events in loop = " << events << endl;

  MpdTrack *track;
  FairMCTrack *mctrack;

  /* Events loop */
  for (Int_t i = 0; i < events-1; i++) 
  {
    Np = 0;
    Npi = 0;
    Nk = 0;
    NpMC = 0;
    NpiMC = 0;
    NkMC = 0;

    dstBranch->GetEntry(i);
    MCBranch->GetEntry(i);
    
    //    MpdEvent *event = (MpdEvent *) mpdEvents->UncheckedAt(i);
       
    TClonesArray *mpdTracks = event->GetGlobalTracks();
    Int_t Ntracks = mpdTracks->GetEntriesFast();

    cout << " *** Event # " << i+1 << " No. of tracks = " << Ntracks << endl;

    for (Int_t k = 0; k < Ntracks; k++)
      {
      track = (MpdTrack*) mpdTracks->UncheckedAt(k);

      ID = track->GetID();

      mctrack = (FairMCTrack*) fT1->UncheckedAt(ID);

      MotherID = mctrack->GetMotherId();

      PDGID = mctrack->GetPdgCode();
      Pt = track->GetPt();
      Pt = TMath::Abs(Pt);

      if (MotherID < 0 && track->GetTofFlag() && Pt > 0.2)
      	{
	  Pz = track->GetPz();
	  Px = track->GetPx();
	  Py = track->GetPy();
	  P = TMath::Sqrt(Pz*Pz + Px*Px +Py*Py);
	  hPt->Fill(Pt);
	  Eta = track->GetEta();
	  hEta->Fill(Eta);
	  if ( Pt <= 2.5 && TMath::Abs(Eta) < 1.5)
	    {
	      ProbP = track->GetPidProbProton();
	      ProbP_dEdx = track->GetTPCPidProbProton();
	      ProbPi = track->GetPidProbPion();
	      ProbPi_dEdx = track->GetTPCPidProbPion();
	      ProbK = track->GetPidProbKaon();
	      ProbK_dEdx = track->GetTPCPidProbKaon();
	      ProbE = track->GetPidProbElectron();
	      ProbE_dEdx = track->GetPidProbElectron();

	      if ( TMath::Abs(ProbP) > ProbCut )
		{
		  Np++;
		  hpPt->Fill(Pt);
		  hpEta->Fill(Eta);
		}
	      if ( TMath::Abs(ProbPi) > ProbCut )
		{
		  Npi++;
		  hpiPt->Fill(Pt);
		  hpiEta->Fill(Eta);
		}
	      if ( TMath::Abs(ProbK) > ProbCut )
		{
		  Nk++;
		  hkPt->Fill(Pt);
		  hkEta->Fill(Eta);
	        }
	    }
	}  //  if (fMotherID < 0 && track->GetTofFlag())
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

      MotherID = mctrack->GetMotherId();

      PDGID = mctrack->GetPdgCode();
	
      if (MotherID < 0 && mctrack->GetNPoints(kTOF))
      	{
	  Pt_MC = mctrack->GetPt();
          Pt_MC = TMath::Abs(Pt_MC);
	  Pz_MC = mctrack->GetPz();
	  E_MC = mctrack->GetEnergy();
          if (Pt_MC > 0)
            {
              Pxyz.SetXYZ(mctrack->GetPx(),mctrack->GetPy(),mctrack->GetPz());
              Eta_MC = Pxyz.PseudoRapidity();
            }
	  if ( Pt_MC <= 1.5 && TMath::Abs(Eta_MC) < 1.5)
	    {
	      Pz_MC = mctrack->GetPz();
	      E_MC = mctrack->GetEnergy();

	      if (TMath::Abs(PDGID) == 2212)
		{
		  NpMC++;
		  hpPtMC->Fill(Pt_MC);
		  hpEtaMC->Fill(Eta_MC);
		}
	      if (TMath::Abs(PDGID) == 211)
		{
		  NpiMC++; 
		  hpiPtMC->Fill(Pt_MC);
		  hpiEtaMC->Fill(Eta_MC);
		}
	      if (TMath::Abs(PDGID) == 321)
		{ 
		  NkMC++;
		  hkPtMC->Fill(Pt_MC);
		  hkEtaMC->Fill(Eta_MC);
		}
	    } //  if ( Pt_MC > 0.2 && TMath::Abs(Eta_MC) < 1.5)
	}  //  if (fMotherID < 0 && mctrack->GetNPoints(kTOF))
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
