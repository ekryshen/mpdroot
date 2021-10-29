/* Macro reads DST file produced by macro reco.C;
It constructs an invariant-mass spectrum of dileptons for all possible e+e-
combinations and gives a background estimate made of e+e+ and e-e- pairs.
It still cheats by using the particle ID of Monte Carlo and the electron mass
in literature.	W. Scheinast, 22/6/2010 */

void ana_dilep_plusplus
	(char* file="/opt/exp_soft/mpd/data4mpd/pluto/AuAu_25gev_ee_central_",
	Int_t istart=0, Int_t iend=1)

{ TStopwatch timer;
  timer.Start();

  /* Load basic libraries if one of PLUTO's classes is not known yet. */
  if (!gROOT->GetListOfClasses()->FindObject("PParticle"))
  { gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
    mpdloadlibs(kTRUE,kTRUE); // not only reco libs
  }

  TH1F* hf0=new TH1F("hf0","dileptons uncut",100,0.,1.6);
  TH1F* hf1=new TH1F("hf1","dileptons cut  ",100,0.,1.6);
  TH1F* hf2=new TH1F("hf2","backg. e+	   ",100,0.,1.6);
  TH1F* hf3=new TH1F("hf3","backg. e-	   ",100,0.,1.6);
  /*TH1F* hf4=new TH1F("hf4","omega component",100,0.,1.6);
  TH1F* hf5=new TH1F("hf5","phi   component",100,0.,1.6);
  TH1F* hf6=new TH1F("hf6","real background",100,0.,1.6);
  TH1F* hf7=new TH1F("hf7","other sources  ",100,0.,1.6);*/
  TH1F* hf8=new TH1F("hf8","fake bck. for uncut",100,0.,1.6);
  TH1F* hf9=new TH1F("hf9","fake bck. for cut",100,0.,1.6);


  Int_t i,j,k;
  char file_[222];
  if (!*file) { cout << "Error: empty input file name." << endl;
		return 0;
	      }
  for (i=0; file[i]; i++) file_[i]=file[i]; k=i; //file_[i]=0;
  if (file[--i]=='/')
  { cout << "Error: input file name required, not only directory!" << endl;
    return 0;
  }


  TChain* dstTree = new TChain("mpdsim");
  cout << endl << endl << "Chain \"mpdsim\":" << endl;
  for (i=istart; i<=TMath::Max(istart,iend); i++)
  { sprintf(file_+k,"%d.root",i);
    if (dstTree->Add(file_)) cout << "	added ";
    else cout << "	failed to add ";
    cout << file_ << endl;
  }
  if (!dstTree->GetEntries())
  { cout << endl << "No input!" << endl << endl;
    return 0;
  }
  //TTree *dstTree = dynamic_cast<TTree*>(fileDST.Get("mpdsim"));


  // Okay, now we have our chain. Let's look into it:
  MpdEvent* event=NULL;
  TClonesArray* MCTracks=NULL;
  dstTree->SetBranchAddress("MPDEvent.", &event);
  dstTree->SetBranchAddress("MCTrack", &MCTracks);
  TClonesArray* mpdTracks;

  Int_t events = dstTree->GetEntries();
  cout << " Number of events in DST files = " << events << endl << endl;
  //events = 1;
  //cout << "  Number of events in loop = " << events << endl << endl;



  TLorentzVector lv1, lv2, dilepton;
  Int_t nback,nback1,nback2,nfake,nfake1,nfake2;
  Int_t n1,n2,ID,PDG,fNtracks,fMCtracks;
  const Double_t emass=.00051099891;	// electron mass in GeV/c^2
  Float_t mass;



  /* Event loop with cut. */
  for (i=0; i<events; i++)

  { j = dstTree->GetEntry(i);
    k = i%1000;  // Only to have some output as a progress display.
    if (!k) cout << "Event " << i << ": " << j << " byte, ";

    mpdTracks = event->GetGlobalTracks();
    fNtracks = mpdTracks->GetEntriesFast();
    fMCtracks = MCTracks->GetEntriesFast();
    if (!k) cout << fNtracks << " global tracks, " << fMCtracks	<< " MCTracks"
	<< endl;


//------------- Now let's analyse!

    // Arrays to keep positrons ("1") and electrons ("2").
    MpdTrack *MPDlept1[50], *MPDlept2[50];	// 50 is necessary!
    FairMCTrack *MClept1[50], *MClept2[50];

    // First we count e+ and e- and store them in arrays.
    for (n1=n2=k=0; k<fNtracks; k++)

    { // As a first guess, we put the particle in the positron array.
      MPDlept1[n1]=static_cast<MpdTrack*>(mpdTracks->UncheckedAt(k));
      // Sometimes the Monte-Carlo ID of a track is screwed up. We skip this.
      if ((ID=MPDlept1[n1]->GetID()) >= fMCtracks) continue;
      MClept1[n1]=static_cast<FairMCTrack*>(MCTracks->UncheckedAt(ID));

      // Here is the room for your cuts and conditions to select the particles;
      // put them in this if statement:
      /*if (TMath::Power(MClept1[n1]->GetStartX(),2)+
	  TMath::Power(MClept1[n1]->GetStartY(),2)<2e-6 &&
	  TMath::Abs(MClept1[n1]->GetStartZ())<.03)*/
      //if (MClept1[n1]->GetMotherId()==-1)
      //if (MPDlept1[n1]->GetPtError()<1)
      if (MPDlept1[n1]->GetThetaError()<.02)
      { PDG=MClept1[n1]->GetPdgCode();	// Here we cheat. This is not analysis!
	if (PDG==-11) n1++;	// Guess was right; just increase e+ counter.
	// Guess was wrong; copy to e- array and increase e- counter n2:
	else if (PDG==11) { MPDlept2[n2]=MPDlept1[n1];
			    MClept2[n2++]=MClept1[n1];
			  }
	// Particles other than e+ and e- are ignored.
      }
    } // Okay, all e+ and e- in this event have been counted.


    // We can only do something if we have at least one e+ and one e-.
    if (n1>0 && n2>0)

    { // Let's assume we have min(n1,n2) "good" (non-background) dilepton pairs.
      nback  = n1*n2-TMath::Min(n1,n2);	// background dileptons
      nfake1 = n1*(n1-1)/2;		// available e+e+ combinations
      nfake2 = n2*(n2-1)/2;		// available e-e- combinations
      nfake  = TMath::Max(nfake1+nfake2,1);  // avoid division by 0
      nback1 = nback*nfake1/nfake;	// integer division, i.e. rounded down!
      // Sometimes we would like to have rounded up, for justice:
      if (gRandom->Rndm()*nfake<nback*nfake1%nfake) nback1++;
      nback2 = nback-nback1;
      // Now nback1 is the part of nback we want to simulate by e+e+
      // and consequently nback2 by e-e-.

      // First we count through the e+.
      for (j=0; j<n1; j++)

      { lv1.SetXYZM (MPDlept1[j]->GetPx(), MPDlept1[j]->GetPy(),
		MPDlept1[j]->GetPz(), emass);

	// ... and for each e+ we take an e-:
	for (k=0; k<n2; k++)

	{ lv2.SetXYZM (MPDlept2[k]->GetPx(), MPDlept2[k]->GetPy(),
		MPDlept2[k]->GetPz(), emass);
	  // emass is also a fraud here; it's not measured.

	  dilepton = lv1 + lv2; // make a 2-particle system
	  mass=dilepton.M();	// invariant mass
	  hf1->Fill(mass);	// fill dilepton mass spectrum
	}

	// And now we make e+e+ combinations with the remnant e+ in the array,
	// at most until nback1 has been counted down to 0.
	for (k=j+1; k<n1 && nback1; k++,nfake1--)

	  /* If we have more e+e+ combinations available (nfake1) than we requi-
	     re (nback1), we decide by a random number whether we take this one
	     or not. nfake1 and nback1 are adjusted for each step in the loop.*/
	  if (nback1==nfake1 || gRandom->Rndm()*nfake1<nback1)
	  { //lv2 = lepton1[k]->Vect4();
	    lv2.SetXYZM (MPDlept1[k]->GetPx(), MPDlept1[k]->GetPy(),
		MPDlept1[k]->GetPz(), emass);
	    dilepton = lv1 + lv2;
	    mass=dilepton.M();
	    hf9->Fill(mass);	// spectrum for all fake dileptons
	    hf2->Fill(mass);	// spectrum for e+ made fake dileptons
	    nback1--;		// one less remaining
	  }
      }

      // This should never happen.
      if (nback1) printf ("Warning: Background with positrons was not computed "
	"entirely. %d dileptons missing\nin event %d.\n",nback1,i);

      // Now we count through the e- to make e-e- combinations.
      // This loop is shorter because the e+e- combinations needn't be done.
      for (j=0; j<n2; j++)

      { // It would be good if the MpdEvent class had a method GetLorentzVector.
	lv1.SetXYZM (MPDlept2[j]->GetPx(), MPDlept2[j]->GetPy(),
		MPDlept2[j]->GetPz(), emass);

	for (k=j+1; k<n2 && nback2; k++,nfake2--)

	if (nback2==nfake2 || gRandom->Rndm()*nfake2<nback2)
	{ lv2.SetXYZM (MPDlept2[k]->GetPx(), MPDlept2[k]->GetPy(),
		MPDlept2[k]->GetPz(), emass);

	  dilepton = lv1 + lv2;
	  mass=dilepton.M();
	  hf9->Fill(mass);
	  hf3->Fill(mass);
	  nback2--;
	}
      }

      if (nback2) printf ("ERROR: Background with electrons was not computed "
	"entirely. %d dileptons missing\nin event %d.\n",nback2,i);
    }




      /* See mpddata/MpdTrack.h for more methods */

    /*} // track loop
    printf("\n");*/
  } // event loop with cut

  cout << endl;



  // Event loop uncut; we repeat the whole process, but without a cut/condition
  // to select "good" particles.
  for (i=0; i<events; i++)

  { j = dstTree->GetEntry(i);
    k = i%1000;
    if (!k) cout << "Event " << i << ": " << j << " byte, ";

    mpdTracks = event->GetGlobalTracks();
    fNtracks = mpdTracks->GetEntriesFast();
    fMCtracks = MCTracks->GetEntriesFast();
    if (!k) cout << fNtracks << " global tracks, " << fMCtracks	<< " MCTracks"
	<< endl;


    MpdTrack *MPDlept1[50], *MPDlept2[50];
    FairMCTrack *MClept1[50], *MClept2[50];

    for (n1=n2=k=0; k<fNtracks; k++)
    { MPDlept1[n1]=static_cast<MpdTrack*>(mpdTracks->UncheckedAt(k));
      if ((ID=MPDlept1[n1]->GetID()) >= fMCtracks) continue;
      MClept1[n1]=static_cast<FairMCTrack*>(MCTracks->UncheckedAt(ID));
      if (true) // This condition is missing here: no cut!
      { PDG=MClept1[n1]->GetPdgCode();
	if (PDG==-11) n1++;
	else if (PDG==11) { MPDlept2[n2]=MPDlept1[n1];
			    MClept2[n2++]=MClept1[n1];
			  }
      }
    }

    if (n1>0 && n2>0)

    { nback  = n1*n2-TMath::Min(n1,n2);
      nfake1 = n1*(n1-1)/2;
      nfake2 = n2*(n2-1)/2;
      nfake  = TMath::Max(nfake1+nfake2,1);
      nback1 = nback*nfake1/nfake;
      if (gRandom->Rndm()*nfake<nback*nfake1%nfake) nback1++;
      nback2 = nback-nback1;

      for (j=0; j<n1; j++)

      { lv1.SetXYZM (MPDlept1[j]->GetPx(), MPDlept1[j]->GetPy(),
		MPDlept1[j]->GetPz(), emass);

	for (k=0; k<n2; k++)

	{ lv2.SetXYZM (MPDlept2[k]->GetPx(), MPDlept2[k]->GetPy(),
		MPDlept2[k]->GetPz(), emass);

	  dilepton = lv1 + lv2;
	  mass=dilepton.M();
	  hf0->Fill(mass);
	}

	for (k=j+1; k<n1 && nback1; k++,nfake1--)

	  if (nback1==nfake1 || gRandom->Rndm()*nfake1<nback1)
	  { lv2.SetXYZM (MPDlept1[k]->GetPx(), MPDlept1[k]->GetPy(),
		MPDlept1[k]->GetPz(), emass);
	    dilepton = lv1 + lv2;
	    mass=dilepton.M();
	    hf8->Fill(mass);
	    nback1--;
	  }
      }

      if (nback1) printf ("Warning: Background with positrons was not computed "
	"entirely. %d dileptons missing\nin event %d.\n",nback1,i);

      for (j=0; j<n2; j++)

      { lv1.SetXYZM (MPDlept2[j]->GetPx(), MPDlept2[j]->GetPy(),
		MPDlept2[j]->GetPz(), emass);

	for (k=j+1; k<n2 && nback2; k++,nfake2--)

	if (nback2==nfake2 || gRandom->Rndm()*nfake2<nback2)
	{ lv2.SetXYZM (MPDlept2[k]->GetPx(), MPDlept2[k]->GetPy(),
		MPDlept2[k]->GetPz(), emass);

	  dilepton = lv1 + lv2;
	  mass=dilepton.M();
	  hf8->Fill(mass);
	  nback2--;
	}
      }

      if (nback2) printf ("ERROR: Background with electrons was not computed "
	"entirely. %d dileptons missing\nin event %d.\n",nback2,i);
    }
  } // event loop uncut

  cout << endl;

  // I misuse the mass variable here as a scaling factor:
  mass = 1./events;
  hf0->Scale(mass);
  hf1->Scale(mass);
  hf2->Scale(mass);
  hf3->Scale(mass);
  hf8->Scale(mass);
  hf9->Scale(mass);
  TH1F* hf7=hf1->Clone("hf7");
  *hf7 = *hf1-*hf9;


  timer.Print();

  cout << endl << "All done." << endl << endl;




  // If there is a canvas, let's make a picture, otherwise not.
  if (gROOT->GetListOfCanvases()->GetEntries())

  { gPad->SetLogy();
    gPad->SetBorderSize(2); //gPad->SetLineWidth(3);
    gPad->SetFillColor(0);
    hf1->SetStats(kFALSE); //hf1->SetMinimum(1);
    //hf1->SetBit(hf1->kNoTitle);
    hf1->SetTitle("#splitline{Au+Au,  #sqrt{s_{NN}} = 7 GeV}{TPC acceptance}");

    // x axis
    hf1->SetXTitle("M (GeV)");
    hf1->GetXaxis()->CenterTitle(true);
    hf1->SetTitleFont(132,"X");
    hf1->SetTitleSize(.04,"X");
    hf1->SetLabelFont(132,"X");
    hf1->SetLabelSize(.03,"X");
    hf1->SetLabelOffset(.005,"X");

    // y axis
    hf1->SetYTitle("(1/N_{ev}) dn/dM, (1/GeV)");
    hf1->GetYaxis()->CenterTitle(true);
    hf1->SetTitleFont(132,"Y");
    hf1->SetTitleSize(.04,"Y");
    hf1->SetLabelFont(132,"Y");
    hf1->SetLabelSize(.03,"Y");
    //hf0->SetLabelOffset(.005,"Y");

    // draw
    gStyle->SetTitleX(.4); gStyle->SetTitleY(.83);
    hf1->SetLineWidth(3); hf1->SetLineColor(1); hf1->Draw();
    hf7->SetLineWidth(3); hf7->SetLineColor(2); hf7->Draw("SAME");
    hf9->SetLineWidth(3); hf9->SetLineColor(4); hf9->Draw("SAME");
    gPad->Modified(); gPad->Update();
    TPaveText* T=gPad->GetPrimitive("title");
    T->SetTextFont(132); T->SetTextSize(.06);
    T->SetFillStyle(0); T->SetBorderSize(0);

    //TLatex t(.7,500.,"all combinations");
    TLatex t;
    t.SetTextFont(132); t.SetTextSize(.04);
    t.DrawLatex(.5,.014,"all combinations")->SetTextColor(1);
    t.DrawLatex(.48,.00015,"difference")->SetTextColor(2);
    t.DrawLatex(.15,.00096,"background estimate")->SetTextColor(4);
    gPad->Modified();
  }
}
