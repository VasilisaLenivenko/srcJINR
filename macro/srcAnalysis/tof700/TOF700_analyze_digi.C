// file for all of the other detectors..
#include <time.h>
void TOF700_analyze_digi(int run, int startEvent, int stopEvent){
	gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
	bmnloadlibs();
	int nfiles=1;
	TChain* rootTree = new TChain("cbmsim");
	for(int i=0; i<nfiles; i++)rootTree->Add(Form("../../raw/bmn_run%d_digi.root",run+i));
	TClonesArray * EventHead = new TClonesArray("BmnEventHeader");
	TClonesArray * ToF700 = new TClonesArray("BmnTof2Digit");
	TClonesArray * ToF400 = new TClonesArray("BmnTof1Digit");
	TClonesArray * T0 = new TClonesArray("BmnTrigDigit");
	TClonesArray * BC3 = new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * GEM = new TClonesArray("BmnGemStripDigit");
	
	Int_t event_count = rootTree->GetEntries();
	cout << "Number of events in _digi file = " << event_count << endl;
	if ( stopEvent == -1) stopEvent = event_count;
/*
	rootTree->GetBranch("EventHeader")->SetAutoDelete(kFALSE);
	rootTree->GetBranch("TOF700")->SetAutoDelete(kFALSE);
	rootTree->GetBranch("TOF400")->SetAutoDelete(kFALSE);
	rootTree->GetBranch("BC2")->SetAutoDelete(kFALSE);
	rootTree->GetBranch("BC3")->SetAutoDelete(kFALSE);
*/
	rootTree->SetBranchAddress("EventHeader", &EventHead);
	rootTree->SetBranchAddress("TOF700", &ToF700);
	rootTree->SetBranchAddress("TOF400", &ToF400);
	rootTree->SetBranchAddress("T0", &T0);
	rootTree->SetBranchAddress("TQDC_BC3", &BC3);
	rootTree->SetBranchAddress("GEM", &GEM);

	TChain* rootTreeReco = new TChain("cbmsim");
	for(int i=0; i<nfiles; i++)rootTreeReco->Add(Form("../../run/bmndst%d.root",run+i));
	TClonesArray * Tof400Hit = new TClonesArray("BmnTof1HitProducer");
	rootTreeReco->SetBranchAddress("BmnTof1Hit", &Tof400Hit);

	// For ToF700 plotting:
	const int fNStrTOF700 = 32;
	const int NPlaneTOF700 = 56;
	// For ToF400 plotting:
	const int fNStr = 48;
	const int NPlane = 20;
	const int fNDetectors = 20;
	TClonesArray * aTofHits = new TClonesArray("BmnTofHit");
	//BmnTOF1Detector **pDetector = new BmnTOF1Detector *[fNDetectors];
	//for (Int_t i = 0; i < fNDetectors; i++) 
          //      pDetector[i] = new BmnTOF1Detector(i, 0);
        TFile *file=new TFile(Form("analyzed%d.root",run), "RECREATE");

	double fTimeL[fNDetectors][fNStr], fTimeR[fNDetectors][fNStr];
	int hits_usable;

	/*TH1D *tofL[fNDetectors][fNStr+1];
	TH1D *tofR[fNDetectors][fNStr+1];
	TString fname;
	
	for( int k = 0 ; k < fNDetectors ; k++){
		for( int j = 0 ; j < fNStr+1 ; j++){
			fname = Form("HistL_ToF_Plane_%d_str%d", k, j);
			//cout << fname << "\n";
			tofL[k][j] = new TH1D(fname,"",2000,-10000,10000);
			//fname = Form("HistR_ToF_Plane_%d_str%d", k, j);
			tofR[k][j] = new TH1D(fname,"",2000,-10000,10000);
		}
	}*/
	cout<<"declare usefull variables"<<endl;
	double fX,fTime700Ref1,fTime700Ref2,fTime700,fBC3Amp,fBC3Time,fT0,fT0Amp,fAmp,fDiff;
	int GemIndex,fPlane,fStrip,fModule,fStation,fStripNumber,fLayer,fStripSignal;
	bool cutT0;
        double dist=2339;//mm
	double angle=30*TMath::DegToRad();
	double pitch=0.8;//mm
	TVector3 fV;
	cout<<"declare histo"<<endl;
	
	TH2D*	hTof700_vs_BC3=new TH2D("hTof700_vs_BC3","hTof700_vs_BC3",600,0,3000,4000,-20,20);
	TH2D*	hTof700_vs_Amp=new TH2D("hTof700_vs_Amp","hTof700_vs_Amp",4000,-20,20,4000,0,8000);
	TH1D*   hTof700=new TH1D("hTof700","hTof700",4000,-200,200);
	TH1D*   hTime700=new TH1D("hTime700","hTime700",1000,0,4000);
	TH1D*   hTimeDiff700=new TH1D("hTimeDiff700","hTimeDiff700",400,-20,20);
	TH2D*   hTof_vs_Amp_BC3=new TH2D("hTof_vs_Amp_BC3","hTof_vs_Amp_BC3",600,0,3000,4000,-8000,8000);
	TH2D*   hTof_vs_Amp_T0=new TH2D("hTof_vs_Amp_T0","hTof_vs_Amp_T0",300,0,300,4000,0,4000);
	TH2D*   hTof700_vs_strip[NPlaneTOF700];
	TH2D*   hTof700_vs_strip_T0cut[NPlaneTOF700];
	TH2D*   hTof700Diff_vs_strip[NPlaneTOF700];
	for(int i=0; i<NPlaneTOF700; i++){
		hTof700_vs_strip[i]=new TH2D(Form("hTof700_vs_strip_%d",i),Form("hTof700_vs_strip_%d",i),32,0,32,200,-10,10);
		hTof700_vs_strip_T0cut[i]=new TH2D(Form("hTof700_vs_strip_T0cut_%d",i),Form("hTof700_vs_strip_T0cut_%d",i),32,0,32,200,-10,10);
		hTof700Diff_vs_strip[i]=new TH2D(Form("hTof700Diff_vs_strip_%d",i),Form("hTof700Diff_vs_strip_%d",i),32,0,32,200,-20,20);
	}
	TH1I*hTof700Plane=new TH1I("hTof700Plane","hTof700Plane",60,0,60);
	TH1D*hGemTheta[4];
	for(int i=0; i<4; i++) hGemTheta[i]=new TH1D(Form("hGemTheta%d",i),Form("hGemTheta%d",i),180,0,90);
	cout<<"start event loop"<<endl;
	for (Int_t j = startEvent; j < stopEvent; j++){
		if (j%100==0) cout << "\tWorking on entry " << j << "\n";

		int nSingleHits = 0;
		EventHead->Clear();

		ToF700->Clear();
		ToF400->Clear();
		T0->Clear();
		BC3->Clear();
		GEM->Clear();
		Tof400Hit->Clear();

		rootTree->GetEntry(j);
		rootTreeReco->GetEntry(j);

	
		if(T0==NULL){ cout<<"NULL"<<endl; continue;}
		int nT0Digits = T0->GetEntriesFast();
		int nBC3Digits = BC3->GetEntriesFast();
            	if (nT0Digits == 1&&nBC3Digits>0) { // T0 digit should be
	//		cout<<"T0"<<endl;
			BmnTrigDigit* digT0 = (BmnTrigDigit*) T0->At(0);
			fT0=0;//digT0->GetTime();
			fT0Amp=digT0->GetAmp();
			cutT0=(fabs(fT0Amp-20)<0.5);
	//		cout<<"BC3"<<endl;
			hTof_vs_Amp_T0->Fill(fT0Amp,fT0);
			BmnTrigWaveDigit* digBC3 = (BmnTrigWaveDigit*) BC3->At(0);
			fBC3Amp=digBC3->GetPeak();
			fBC3Time=digBC3->GetTime();
			hTof_vs_Amp_BC3->Fill(fBC3Amp,fBC3Time-fT0);
			for (Int_t i = 0; i < fNDetectors; i++){
                    		//pDetector[i]->Clear();
				for (int w = 0 ; w < fNStr ; w++){
					fTimeL[i][w] = 0.;
					fTimeR[i][w] = 0.;
				}
			}
			//cout<<"GEM"<<endl;
			
			for (Int_t iDig = 0; iDig < GEM->GetEntriesFast(); ++iDig) {
				BmnGemStripDigit* digGem = (BmnGemStripDigit*) GEM->At(iDig);
                        	fStation=digGem->GetStation();
                        	fModule=digGem->GetModule();
                        	fLayer=digGem->GetStripLayer();
                        	fStripNumber=digGem->GetStripNumber();
                        	fStripSignal=digGem->GetStripSignal();
				//pDetector[digTof->GetPlane()]->SetDigit(digTof);
				GemIndex=-1;
				if(fStation==0&&fModule==0&&fLayer==0&&fStripSignal>200)
					{fX=fStripNumber*pitch-825.*pitch/2;fV.SetXYZ(fX,0,dist);fV.RotateY(-angle);GemIndex=0;}
				//cout<<"station "<<GemIndex<<endl;
				if(fStation==0&&fModule==1&&fLayer==0&&fStripSignal>200)
					{fX=fStripNumber*pitch-825.*pitch/2;fV.SetXYZ(fX,0,dist);fV.RotateY(-angle);GemIndex=1;}
				//cout<<"station "<<GemIndex<<endl;
				if(fStation==1&&fModule==0&&fLayer==0&&fStripSignal>200)
					{fX=fStripNumber*pitch-825.*pitch/2;fV.SetXYZ(fX,0,dist);fV.RotateY(angle);GemIndex=2;}
				//cout<<"station "<<GemIndex<<endl;
				if(fStation==1&&fModule==1&&fLayer==0&&fStripSignal>200)
					{fX=fStripNumber*pitch-825.*pitch/2;fV.SetXYZ(fX,0,dist);fV.RotateY(angle);GemIndex=3;}
				//cout<<"station "<<GemIndex<<endl;
				if(GemIndex!=-1)hGemTheta[GemIndex]->Fill(fV.X()/fV.Z()*TMath::RadToDeg());
				
			}
		//cout<<"BC3 "<<fBC3Amp<<endl;	
	//		cout<<"TOF700"<<endl;
			for (Int_t iDig = 0; iDig < ToF700->GetEntriesFast(); ++iDig) {
				BmnTof2Digit* digTof700 = (BmnTof2Digit*) ToF700->At(iDig);
                        	
				//pDetector[digTof->GetPlane()]->SetDigit(digTof);
				
				if(digTof700->GetTime() == 0) continue;
				fTime700=digTof700->GetTime();
				fDiff=digTof700->GetDiff();
				fAmp=digTof700->GetAmplitude();
				fPlane=digTof700->GetPlane();
				fStrip=digTof700->GetStrip();
				hTof700_vs_BC3->Fill(fBC3Amp,fTime700);
				//cout<<"T "<<fTime700<<" "<<fTime700-fT0<<endl;	
				//hTime700->Fill(fTime700);
				hTof700->Fill(fTime700);
				hTof700_vs_strip[fPlane]->Fill(fStrip,fTime700);
				if(cutT0)hTof700_vs_strip_T0cut[fPlane]->Fill(fStrip,fTime700);
				hTof700Diff_vs_strip[fPlane]->Fill(fStrip,fDiff);
				hTof700Plane->Fill(fPlane);
				if(fPlane==19&&fStrip==10){hTof700_vs_Amp->Fill(fTime700,fAmp);fTime700Ref1=fTime700;}
				if(fPlane==6&&fStrip==10){fTime700Ref2=fTime700;}
	//			cout<<digTof700->GetPlane()<<endl;
			}
	//		cout<<"TOF400"<<endl;
			
			hTimeDiff700->Fill(fTime700Ref2-fTime700Ref1);
			
			for (Int_t iDig = 0; iDig < ToF400->GetEntriesFast(); ++iDig) {
				BmnTof1Digit* digTof = (BmnTof1Digit*) ToF400->At(iDig);
                        	
				//pDetector[digTof->GetPlane()]->SetDigit(digTof);
				
				if(digTof->GetTime() == 0) continue;

				if(digTof->GetSide()==0){
					fTimeL[digTof->GetPlane()][digTof->GetStrip()] = digTof->GetTime();
					//tofL[digTof->GetPlane()][digTof->GetStrip()]->Fill(digTof->GetTime());
				}		
				if(digTof->GetSide()==1){
					fTimeR[digTof->GetPlane()][digTof->GetStrip()] = digTof->GetTime();
					//tofR[digTof->GetPlane()][digTof->GetStrip()]->Fill(digTof->GetTime());
				}

				//tof[digTof->GetPlane()][digTof->GetStrip()]->Fill(digTof->GetTime());
				//cout << digTof->GetTime() - digT0->GetTime() << "\n";
			}

			for (Int_t iHit = 0; iHit < Tof400Hit->GetEntriesFast(); ++iHit) {
				BmnTof1HitProducer* hit = (BmnTof1HitProducer*) Tof400Hit->At(iHit);
                        	
			}

			hits_usable++;
/*
                	for (Int_t i = 0; i < fNDetectors; i++) {
                     		for (int w = 0 ; w < fNStr ; w++){
					if( (hTimeL[i][w]!=0) || (hTimeR[i][w]!=0) ){
						//cout << "\t\tcandidate event in strip: " << w << "\n";	
						hits_usable++;
					}
				}
				//nSingleHits += pDetector[i] -> FindHits(digT0, aTofHits);
			}
*/
		}
	}
	/*
	for (Int_t i = 0 ; i < fNDetectors; i++){
		cout << "Writing Plane: " << i << "\n";
		TFile *ptr = gFile;
		TList* fList = (TList*) pDetector[i]->GetList(2);
		
		fList->Write();
	}
	*/
	cout << "Number of total hits usable in file: " << hits_usable << "\n";
        file->Write();	
	//file->Close();
	//TFile *ptr = gFile;
	//gFile->cd();
        //TList* fList = (TList*) pDetector[1]->GetList(2);
	//fList->Write();


}
