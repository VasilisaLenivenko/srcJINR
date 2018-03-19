void land(char *fname) {
	gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
	bmnloadlibs();

	TChain *bmnTree = new TChain("cbmsim");
	if (bmnTree->Add(fname) == 0)
	{
		printf("Can't find BMN digits tree in file %s!\n", fname);
		return;
	}
	Int_t nEvents = bmnTree->GetEntries();

	TClonesArray *LANDDigits;
	bmnTree->SetBranchAddress("LAND", &LANDDigits);

	TH1F *hT0 = new TH1F("hT0", "T0", 500, 0, 2000);
	TH2F *hEvsToF = new TH2F("hEvsToF", "Energy vs ToF", 100, 0, 100, 100, 0, 4000);
	TH2F *hEvsVel = new TH2F("hEvsVel", "Energy vs Velocity", 1000, -300, 300, 100, 0, 4000);
	TH2F *hT0vsT = new TH2F("hT0vsT", "T0 vs T", 500, 0, 1000, 500, 0, 2000);
	TH2F *hCutvsToF = new TH2F("hCutvsToF", "Cut vs ToF", 1000, -300, 300, 100, 0, 4000);
	TH2F *hToFvsPlane = new TH2F("hToFvsPlane", "ToF vs plane", 6, 0, 6, 1000, -300, 300);
	TH2F *hTDiffvsBar = new TH2F("hTDiffvsBar", "TDiff vs bar", 120, 0, 120, 500, -300, 300);
	TH2F *hTvsBar = new TH2F("hTvsBar", "T vs bar", 120, 0, 120, 500, 0, 1000);
	TH2F *hToFvsBar = new TH2F("hToFvsBar", "ToF vs bar", 120, 0, 120, 1000, -300, 300);
	TH2F *hYvsX = new TH2F("hYvsX", "Y vs X", 30, -150, 150, 30, -150, 150);

	for (Int_t iEv = 0; iEv < nEvents; iEv++) {
		bmnTree->GetEntry(iEv);

		if (iEv % 10000 == 0) {
			cout << "Event: " << iEv << "/" << nEvents << endl;
		}
		if (0 == LANDDigits->GetEntriesFast()) {
			continue;
		}
		BmnLANDDigit *t0 = (BmnLANDDigit *)LANDDigits->At(0);
		if (!t0->IsT0()) {
			continue;
		}
		hT0->Fill(t0->GetTDiff(0));
		for (Int_t iDig = 1; iDig < LANDDigits->GetEntriesFast(); ++iDig) {
			BmnLANDDigit *digit = (BmnLANDDigit *)LANDDigits->At(iDig);
			double ToF = digit->GetTime() - t0->GetTDiff(0) + 75;
			double dx = digit->GetX() - 131;
			double dy = digit->GetY();
			double dz = 1425 + digit->GetPlane() * 10;
			double vel = sqrt(dx * dx + dy * dy + dz * dz) / ToF;
			hEvsToF->Fill(ToF, digit->GetEnergy());
			hEvsVel->Fill(vel, digit->GetEnergy());
			hT0vsT->Fill(digit->GetTime(), t0->GetTDiff(0));
			if (10 == digit->GetGlobBar() &&
			    abs(digit->GetPosition()) < 30) {
				hCutvsToF->Fill(ToF, digit->GetEnergy());
			}
			hToFvsPlane->Fill(digit->GetPlane(), ToF);
			hTDiffvsBar->Fill(digit->GetGlobBar(), digit->GetTime(1) - digit->GetTime(0));
			hTvsBar->Fill(digit->GetGlobBar(), digit->GetTime());
			hToFvsBar->Fill(digit->GetGlobBar(), ToF);
			if(digit->GetPlane() == 0)
				hYvsX->Fill(digit->GetX(), digit->GetY());
		}
	}

	hT0->Draw();
	TCanvas *c0 = new TCanvas;
	hEvsToF->Draw("colz");
	TCanvas *c1 = new TCanvas;
	hEvsVel->Draw("colz");
	TCanvas *c2 = new TCanvas;
	hT0vsT->Draw("colz");
	TCanvas *c3 = new TCanvas;
	hCutvsToF->Draw("colz");
	TCanvas *c4 = new TCanvas;
	hToFvsPlane->Draw("colz");
	TCanvas *c5 = new TCanvas;
	hTDiffvsBar->Draw("colz");
	TCanvas *c6 = new TCanvas;
	hTvsBar->Draw("colz");
	TCanvas *c7 = new TCanvas;
	hToFvsBar->Draw("colz");
	TCanvas *c8 = new TCanvas;
	hYvsX->Draw("colz");
}
