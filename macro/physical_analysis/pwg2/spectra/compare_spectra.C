/* Macro reads DST file produced by macro reco.C */

void compare_spectra(TString inFile = "dst_auau_09gev_0_3fm_0.root")
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
  Float_t weight, Np_w, Npi_w, Nk_w;
  Int_t NpMC = 0, NpiMC = 0, NkMC = 0;
  Int_t N_ToF = 0, N_ToF_F = 0;
  Float_t ProbCut = 0.8;
  Float_t ProbCut_dEdx = 0.7;
  Int_t MotherID, PDGID;

  enum DetectorId {kSTS, kTPC, kTOF, kETOF, kFFD, kECT, kECAL, kNDET, kCPC, kBBC, kZDC, kFSA};

  const Double_t kProtonMass = 0.938271998;
  const Double_t kPionMass = 0.13957018;
  const Double_t kKaonMass = 0.493677;

  fFile->cd();

  TH1F *hPt = new TH1F("Pt","",100,0.0,3.0);
   hPt->GetYaxis()->SetTitle("dN/dP_{t}"); 
   hPt->GetYaxis()->SetTitleOffset(1.2);
   hPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hEta = new TH1F("Eta","",100,-2.0,2.0);
   hEta->GetYaxis()->SetTitle("dN/d#eta"); 
   hEta->GetYaxis()->SetTitleOffset(1.2);
   hEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpart = new TH1F("hpNpart","",50,0.0,200.0);
   hpNpart->GetYaxis()->SetTitle("dN/dN_{part}");
   hpNpart->GetYaxis()->SetTitleOffset(1.2);
   hpNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpNpartw = new TH1F("hpNpartw","",50,0.0,200.0);
   hpNpartw->GetYaxis()->SetTitle("dN/dN_{part}"); 
   hpNpartw->GetYaxis()->SetTitleOffset(1.2);
   hpNpartw->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPt = new TH1F("hpPt","",100,0.0,3.0);
   hpPt->GetYaxis()->SetTitle("dN/dP_{t}");
   hpPt->GetYaxis()->SetTitleOffset(1.2);
   hpPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpPtw = new TH1F("hpPtw","",100,0.0,3.0);
   hpPtw->GetYaxis()->SetTitle("dN/dP_{t}");
   hpPtw->GetYaxis()->SetTitleOffset(1.2);
   hpPtw->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpMt = new TH1F("hpMt","",100,0.0,3.0);
   hpMt->GetYaxis()->SetTitle("dN/M_{t}dM_{t}");
   hpMt->GetYaxis()->SetTitleOffset(1.2);
   hpMt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpEta = new TH1F("hpEta","",100,-2.0,2.0);
   hpEta->GetYaxis()->SetTitle("dN/d#eta");
   hpEta->GetYaxis()->SetTitleOffset(1.2);
   hpEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpNpartMC = new TH1F("hpNpartMC","",50,0.0,200.0);
   hpNpartMC->GetYaxis()->SetTitle("dN/dN_{part}");
   hpNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpPtMC = new TH1F("hpPtMC","",100,0.0,3.0);
   hpPtMC->GetYaxis()->SetTitle("dN/dP_{t}");
   hpPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpMtMC = new TH1F("hpMtMC","",100,0.0,3.0);
   hpMtMC->GetYaxis()->SetTitle("dN/M_{t}dM_{t}");
   hpMtMC->GetYaxis()->SetTitleOffset(1.2);
   hpMtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpEtaMC = new TH1F("hpEtaMC","",100,-2.0,2.0);
   hpEtaMC->GetYaxis()->SetTitle("dN/d#eta"); 
   hpEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hpEtaMC->GetXaxis()->SetTitle("#eta");

  TH1F *hpM2 = new TH1F("hpM2","",100,0.0,2.0);
   hpM2->GetYaxis()->SetTitle("dN/dM^{2}");
   hpM2->GetYaxis()->SetTitleOffset(1.2);
   hpM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH1F *hpdEdx = new TH1F("hpdEdx","",100,0.0,0.00001);
   hpdEdx->GetYaxis()->SetTitle("dN/(dE/dx)");
   hpdEdx->GetYaxis()->SetTitleOffset(1.2);
   hpdEdx->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hpdEdxM2 = new TH2F("hpdEdxM2","",100,0.0,2.0,100,0.0,0.00001);
   hpdEdxM2->GetYaxis()->SetTitle("dE/dx, GeV/cm"); 
   hpdEdxM2->GetYaxis()->SetTitleOffset(1.2);
   hpdEdxM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hpdEdxP = new TH2F("hpdEdxP","",100,0.0,0.00001,100,0.0,2.0);
   hpdEdxP->GetYaxis()->SetTitle("P, GeV/c");
   hpdEdxP->GetYaxis()->SetTitleOffset(1.2);
   hpdEdxP->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hpM2P = new TH2F("hpM2P","",100,0.0,2.0,100,0.0,2.0);
   hpM2P->GetYaxis()->SetTitle("P, GeV/c"); 
   hpM2P->GetYaxis()->SetTitleOffset(1.2);
   hpM2P->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hpProob2 = new TH2F("hpProb2","",100,0.5,1.0,100,0.5,1.0);
   hpProb2->GetYaxis()->SetTitle("Probability (dE/dx)"); 
   hpProb2->GetYaxis()->SetTitleOffset(1.2);
   hpProb2->GetXaxis()->SetTitle("Probability (TOF)");

  TH1F *hpiNpart = new TH1F("hpiNpart","",100,0.0,500.0);
   hpiNpart->GetYaxis()->SetTitle("dN/dN_{part}"); 
   hpiNpart->GetYaxis()->SetTitleOffset(1.2);
   hpiNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiNpartw = new TH1F("hpiNpartw","",100,0.0,500.0);
   hpiNpartw->GetYaxis()->SetTitle("dN/dN_{part}");
   hpiNpartw->GetYaxis()->SetTitleOffset(1.2);
   hpiNpartw->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPt = new TH1F("hpiPt","",100,0.0,3.0);
   hpiPt->GetYaxis()->SetTitle("dN/dP_{t}"); 
   hpiPt->GetYaxis()->SetTitleOffset(1.2);
   hpiPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiPtw = new TH1F("hpiPtw","",100,0.0,3.0);
   hpiPtw->GetYaxis()->SetTitle("dN/dP_{t}");
   hpiPtw->GetYaxis()->SetTitleOffset(1.2);
   hpiPtw->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiMt = new TH1F("hpiMt","",100,0.0,3.0);
   hpiMt->GetYaxis()->SetTitle("dN/M_{t}dM_{t}");
   hpiMt->GetYaxis()->SetTitleOffset(1.2);
   hpiMt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiMtMC = new TH1F("hpiMtMC","",100,0.0,3.0);
   hpiMtMC->GetYaxis()->SetTitle("dN/M_{t}dM_{t}");
   hpiMtMC->GetYaxis()->SetTitleOffset(1.2);
   hpiMtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiEta = new TH1F("hpiEta","",100,-2.0,2.0);
   hpiEta->GetYaxis()->SetTitle("dN/d#eta");
   hpiEta->GetYaxis()->SetTitleOffset(1.2);
   hpiEta->GetXaxis()->SetTitle("#eta");

  TH1F *hpiNpartMC = new TH1F("hpiNpartMC","",100,0.0,500.0);
   hpiNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); 
   hpiNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hpiNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hpiPtMC = new TH1F("hpiPtMC","",100,0.0,3.0);
   hpiPtMC->GetYaxis()->SetTitle("dN/dP_{t}"); 
   hpiPtMC->GetYaxis()->SetTitleOffset(1.2);
   hpiPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hpiEtaMC = new TH1F("hpiEtaMC","",100,-2.0,2.0);
   hpiEtaMC->GetYaxis()->SetTitle("dN/d#eta");
   hpiEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hpiEtaMC->GetXaxis()->SetTitle("#eta");

  TH1F *hpiM2 = new TH1F("hpiM2","",100,0.0,2.0);
   hpiM2->GetYaxis()->SetTitle("dN/dM^{2}"); 
   hpiM2->GetYaxis()->SetTitleOffset(1.2);
   hpiM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH1F *hpidEdx = new TH1F("hpidEdx","",100,0.0,0.00001);
   hpidEdx->GetYaxis()->SetTitle("dN/(dE/dx)"); 
   hpidEdx->GetYaxis()->SetTitleOffset(1.2);
   hpidEdx->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hpidEdxM2 = new TH2F("hpidEdxM2","",100,0.0,2.0,100,0.0,0.00001);
   hpidEdxM2->GetYaxis()->SetTitle("dE/dx, GeV/cm"); 
   hpidEdxM2->GetYaxis()->SetTitleOffset(1.2);
   hpidEdxM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hpidEdxP = new TH2F("hpidEdxP","",100,0.0,0.00001,100,0.0,2.0);
   hpidEdxP->GetYaxis()->SetTitle("P, GeV/c");
   hpidEdxP->GetYaxis()->SetTitleOffset(1.2);
   hpidEdxP->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hpiM2P = new TH2F("hpiM2P","",100,0.0,2.0,100,0.0,2.0);
   hpiM2P->GetYaxis()->SetTitle("P, GeV/c"); 
   hpiM2P->GetYaxis()->SetTitleOffset(1.2);
   hpiM2P->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hpiProob2 = new TH2F("hpiProb2","",100,0.5,1.0,100,0.5,1.0);
   hpiProb2->GetYaxis()->SetTitle("Probability (dE/dx)"); 
   hpiProb2->GetYaxis()->SetTitleOffset(1.2);
   hpiProb2->GetXaxis()->SetTitle("Probability (TOF)");

  TH1F *hkNpart = new TH1F("hkNpart","",50,0.0,100.0);
   hkNpart->GetYaxis()->SetTitle("dN/dN_{part}");
   hkNpart->GetYaxis()->SetTitleOffset(1.2);
   hkNpart->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkNpartw = new TH1F("hkNpartw","",50,0.0,100.0);
   hkNpartw->GetYaxis()->SetTitle("dN/dN_{part}"); 
   hkNpartw->GetYaxis()->SetTitleOffset(1.2);
   hkNpartw->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPt = new TH1F("hkPt","",100,0.0,3.0);
   hkPt->GetYaxis()->SetTitle("dN/dP_{t}"); 
   hkPt->GetYaxis()->SetTitleOffset(1.2);
   hkPt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkPtw = new TH1F("hkPtw","",100,0.0,3.0);
   hkPtw->GetYaxis()->SetTitle("dN/dP_{t}"); 
   hkPtw->GetYaxis()->SetTitleOffset(1.2);
   hkPtw->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkMt = new TH1F("hkMt","",100,0.0,3.0);
   hkMt->GetYaxis()->SetTitle("dN/M_{t}dM_{t}"); 
   hkMt->GetYaxis()->SetTitleOffset(1.2);
   hkMt->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkMtMC = new TH1F("hkMtMC","",100,0.0,3.0);
   hkMtMC->GetYaxis()->SetTitle("dN/M_{t}dM_{t}");
   hkMtMC->GetYaxis()->SetTitleOffset(1.2);
   hkMtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkEta = new TH1F("hkEta","",100,-2.0,2.0);
   hkEta->GetYaxis()->SetTitle("dN/d#eta");
   hkEta->GetYaxis()->SetTitleOffset(1.2);
   hkEta->GetXaxis()->SetTitle("#eta");

   TH1F *hkNpartMC = new TH1F("hkNpartMC","", 50,0.0,100.0);
   hkNpartMC->GetYaxis()->SetTitle("dN/dN_{part}"); 
   hkNpartMC->GetYaxis()->SetTitleOffset(1.2);
   hkNpartMC->GetXaxis()->SetTitle("N_{part}");

  TH1F *hkPtMC = new TH1F("hkPtMC","",100,0.0,3.0);
   hkPtMC->GetYaxis()->SetTitle("dN/dP_{t}");
   hkPtMC->GetYaxis()->SetTitleOffset(1.2);
   hkPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");

  TH1F *hkEtaMC = new TH1F("hkEtaMC","",100,-2.0,2.0);
   hkEtaMC->GetYaxis()->SetTitle("dN/d#eta");
   hkEtaMC->GetYaxis()->SetTitleOffset(1.2);
   hkEtaMC->GetXaxis()->SetTitle("#eta");

  TH1F *hkM2 = new TH1F("hkM2","",100,0.0,2.0);
   hkM2->GetYaxis()->SetTitle("dN/dM^{2}");
   hkM2->GetYaxis()->SetTitleOffset(1.2);
   hkM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH1F *hkdEdx = new TH1F("hkdEdx","",100,0.0,0.00001);
   hkdEdx->GetYaxis()->SetTitle("dN/(dE/dx)");
   hkdEdx->GetYaxis()->SetTitleOffset(1.2);
   hkdEdx->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hkdEdxM2 = new TH2F("hkdEdxM2","",100,0.0,2.0,100,0.0,0.00001);
   hkdEdxM2->GetYaxis()->SetTitle("dE/dx, GeV/cm");
   hkdEdxM2->GetYaxis()->SetTitleOffset(1.2);
   hkdEdxM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hkdEdxP = new TH2F("hkdEdxP","",100,0.0,0.00001,100,0.0,2.0);
   hkdEdxP->GetYaxis()->SetTitle("P, GeV/c"); 
   hkdEdxP->GetYaxis()->SetTitleOffset(1.2);
   hkdEdxP->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *hkM2P = new TH2F("hkM2P","",100,0.0,2.0,100,0.0,2.0);
   hkM2P->GetYaxis()->SetTitle("P, GeV/c");
   hkM2P->GetYaxis()->SetTitleOffset(1.2);
   hkM2P->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hkProob2 = new TH2F("hkProb2","",100,0.5,1.0,100,0.5,1.0);
   hkProb2->GetYaxis()->SetTitle("Probability (dE/dx)");
   hkProb2->GetYaxis()->SetTitleOffset(1.2);
   hkProb2->GetXaxis()->SetTitle("Probability (TOF)");

  TH2F *hedEdxM2 = new TH2F("hedEdxM2","",100,0.0,2.0,100,0.0,0.00002);
   hedEdxM2->GetYaxis()->SetTitle("dE/dx, GeV/cm"); 
   hedEdxM2->GetYaxis()->SetTitleOffset(1.2);
   hedEdxM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *hedEdxP = new TH2F("hedEdxP","",100,0.0,0.00001,100,0.0,2.0);
   hedEdxP->GetYaxis()->SetTitle("P, GeV/c");
   hedEdxP->GetYaxis()->SetTitleOffset(1.2);
   hedEdxP->GetXaxis()->SetTitle("dE/dx, GeV/cm");

  TH2F *heM2P = new TH2F("heM2P","",100,0.0,2.0,100,0.0,2.0);
   heM2P->GetYaxis()->SetTitle("P, GeV/c"); 
   heM2P->GetYaxis()->SetTitleOffset(1.2);
   heM2P->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");

  TH2F *heProob2 = new TH2F("heProb2","",100,0.5,1.0,100,0.5,1.0);
   heProb2->GetYaxis()->SetTitle("Probability (dE/dx)");
   heProb2->GetYaxis()->SetTitleOffset(1.2);
   heProb2->GetXaxis()->SetTitle("Probability (TOF)");

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
    Np_w = 0;
    Npi_w = 0;
    Nk_w = 0;
    NpMC = 0;
    NpiMC = 0;
    NkMC = 0;

    dstBranch->GetEntry(i);
    MCBranch->GetEntry(i);
      
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
	  M2 = track->GetTofMass2();
	  dEdx = track->GetdEdXTPC();
	  if ( TMath::Abs(Eta) < 1.5)
	    {
	      ProbP = track->GetPidProbProton();
	      ProbP_dEdx = track->GetTPCPidProbProton();
	      ProbPi = track->GetPidProbPion();
	      ProbPi_dEdx = track->GetTPCPidProbPion();
	      ProbK = track->GetPidProbKaon();
	      ProbK_dEdx = track->GetTPCPidProbKaon();
	      ProbE = track->GetPidProbElectron();
	      ProbE_dEdx = track->GetPidProbElectron();

	      weight = TMath::Abs(ProbPi);
	      Npi_w += weight;
	      hpiPtw->Fill(Pt,weight);

	      weight = TMath::Abs(ProbP);
	      Np_w += weight;
	      hpPtw->Fill(Pt,weight);

	      weight = TMath::Abs(ProbK);
	      Nk_w += weight;
	      hkPtw->Fill(Pt,weight);

	      if ( TMath::Abs(ProbP) > ProbCut )
		{
		  Np++;
		  hpPt->Fill(Pt);
		  hpEta->Fill(Eta);
		  hpM2->Fill(M2);
		  hpdEdx->Fill(dEdx);
		  hpdEdxM2->Fill(M2,dEdx);
		  hpdEdxP->Fill(dEdx,P);
		  hpM2P->Fill(M2,P);
		  hpProb2->Fill(ProbP, ProbP_dEdx);
		  E = TMath::Sqrt(Pz*Pz + Px*Px +Py*Py + kProtonMass*kProtonMass);
		  if (E > Pz) 
		    {
		      Mt = TMath::Sqrt(E*E - Pz*Pz);
		      hpMt->Fill(Mt-kProtonMass,1.0/Mt);
		    }
		}
	      if ( TMath::Abs(ProbPi) > ProbCut )
		{
		  Npi++;
		  hpiPt->Fill(Pt);
		  hpiEta->Fill(Eta);
		  hpiM2->Fill(M2);
		  hpidEdx->Fill(dEdx);
		  hpidEdxM2->Fill(M2,dEdx);
		  hpidEdxP->Fill(dEdx,P);
		  hpiM2P->Fill(M2,P);
		  hpiProb2->Fill(ProbPi, ProbPi_dEdx);
		  E = TMath::Sqrt(Pz*Pz + Px*Px +Py*Py + kPionMass*kPionMass);
		  if (E > Pz) 
		    {
		      Mt = TMath::Sqrt(E*E - Pz*Pz);
		      hpiMt->Fill(Mt-kPionMass,1.0/Mt);
		    }
		}
	      if ( TMath::Abs(ProbK) > ProbCut )
		{
		  Nk++;
		  hkPt->Fill(Pt);
		  hkEta->Fill(Eta);
		  hkM2->Fill(M2);
		  hkdEdx->Fill(dEdx);
		  hkdEdxM2->Fill(M2,dEdx);
		  hkdEdxP->Fill(dEdx,P);
		  hkM2P->Fill(M2,P);
		  hkProb2->Fill(ProbK, ProbK_dEdx);
		  E = TMath::Sqrt(Pz*Pz + Px*Px +Py*Py + kKaonMass*kKaonMass);
		  if (E > Pz) 
		    {
		      Mt = TMath::Sqrt(E*E - Pz*Pz);
		      hkMt->Fill(Mt-kKaonMass,1.0/Mt);
		    }
		}
	      if ( TMath::Abs(ProbE) > ProbCut )
		{
		  hedEdxM2->Fill(M2,dEdx);
		  hedEdxP->Fill(dEdx,P);
		  heM2P->Fill(M2,P);
		  heProb2->Fill(ProbK, ProbK_dEdx);
		}
	    }
	}  //  if (fMotherID < 0 && track->GetTofFlag())
    } // track loop
    hpNpart->Fill(Np);
    hpiNpart->Fill(Npi);
    hkNpart->Fill(Nk);

    hpNpartw->Fill(Np_w);
    hpiNpartw->Fill(Npi_w);
    hkNpartw->Fill(Nk_w);

    Int_t mcNtracks = fT1->GetEntries();

    for (Int_t j = 0; j < mcNtracks; j++)
    {
      TVector3 Pxyz;
      mctrack = (FairMCTrack*)fT1->At(j);

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
	  if ( Pt_MC >= 0.2 && TMath::Abs(Eta_MC) < 1.5)
	    {
	      Pz_MC = mctrack->GetPz();
	      E_MC = mctrack->GetEnergy();
	      Mt_MC = TMath::Sqrt(E_MC*E_MC - Pz_MC*Pz_MC);

	      if (TMath::Abs(PDGID) == 2212)
		{
		  NpMC++;
		  hpPtMC->Fill(Pt_MC);
		  hpEtaMC->Fill(Eta_MC);
		  hpMtMC->Fill(Mt_MC-kProtonMass,1.0/Mt_MC);
		}
	      if (TMath::Abs(PDGID) == 211)
		{
		  NpiMC++; 
		  hpiPtMC->Fill(Pt_MC);
		  hpiEtaMC->Fill(Eta_MC);
		  hpiMtMC->Fill(Mt_MC-kPionMass,1.0/Mt_MC);
		}
	      if (TMath::Abs(PDGID) == 321)
		{ 
		  NkMC++;
		  hkPtMC->Fill(Pt_MC);
		  hkEtaMC->Fill(Eta_MC);
		  hkMtMC->Fill(Mt_MC-kKaonMass,1.0/Mt_MC);
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
  exit(0);
}
