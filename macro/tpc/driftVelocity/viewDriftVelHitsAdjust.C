void viewDriftVelHitsAdjust(std::string inFile0 = "mpddst_laser_makeMap.root",std::string inFile1 = "mpddst_laser_adjMap.root", long long event = 0)
{
    TGraph * gry = new TGraph();
	gry->SetMarkerStyle(kFullCircle);
	gry->SetMarkerColor(kRed);
	gry->SetMarkerSize(0.25);
	TGraph * gry_adj = new TGraph();
	gry_adj->SetMarkerStyle(kFullCircle);
	gry_adj->SetMarkerColor(kBlack);
	gry_adj->SetMarkerSize(0.25);
	
	TGraph2D * gr2d_adj = new TGraph2D();
	gr2d_adj->SetMarkerStyle(kFullCircle);
	gr2d_adj->SetMarkerColor(kBlack);
	gr2d_adj->SetMarkerSize(0.25);
	
	TFile * file0 = TFile::Open(inFile0.c_str());
	TTree * tree0 = (TTree *)file0->Get("mpdsim");
	TClonesArray * fTpcRecPts = 0;
	tree0->SetBranchAddress("TpcRecPoint", &fTpcRecPts);
	tree0->GetEntry(event);
	int nHits = fTpcRecPts->GetEntriesFast();
	
	for (Int_t h = 0; h < nHits; ++h)
	{
		MpdTpcHit * fTpcHit = (MpdTpcHit *)fTpcRecPts->UncheckedAt(h);
		TVector3 hitPos;
		fTpcHit->Position(hitPos);
		gry->SetPoint(h, hitPos.Z(), hitPos.Y());
	}
	
	TFile * file1 = TFile::Open(inFile1.c_str());
	TTree * tree1 = (TTree *)file1->Get("mpdsim");
	TClonesArray * fTpcRecPtsAdj = 0;
	tree1->SetBranchAddress("TpcRecPoint", &fTpcRecPtsAdj);
	tree1->GetEntry(event);
	int nHitsAdj = fTpcRecPtsAdj->GetEntriesFast();
	gr2d_adj->Set(nHitsAdj);

	for (Int_t h = 0; h < nHitsAdj; ++h)
	{
		MpdTpcHit * fTpcHit = (MpdTpcHit *)fTpcRecPtsAdj->UncheckedAt(h);
		TVector3 hitPos;
		fTpcHit->Position(hitPos);
		gry_adj->SetPoint(h, hitPos.Z(), hitPos.Y());
		gr2d_adj->SetPoint(h, hitPos.X(), hitPos.Y(), hitPos.Z());
	}
	
	TCanvas * c = new TCanvas();
	c->Divide(2,1);
	c->cd(1);
	TMultiGraph * mgy = new TMultiGraph();
	mgy->Add(gry, "P");
	mgy->Add(gry_adj, "P");
	mgy->GetXaxis()->SetTitle("Z");
	mgy->GetYaxis()->SetTitle("Y");
	mgy->Draw("A");
	c->cd(2);
	gr2d_adj->Draw();
}
