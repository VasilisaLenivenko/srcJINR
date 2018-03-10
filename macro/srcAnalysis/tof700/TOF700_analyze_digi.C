// file for all of the other detectors..
#include <time.h>

void TOF700_analyze_digi(int run, int startEvent, int stopEvent){
	gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
	bmnloadlibs();

	TChain* rootTree = new TChain("cbmsim");
	for(int i=0; i<2; i++)rootTree->Add(Form("../../raw/bmn_run%d_digi.root",run+i));

	TClonesArray * EventHead = new TClonesArray("BmnEventHeader");
	TClonesArray * ToF700 = new TClonesArray("BmnTof2Digit");
	TClonesArray * ToF400 = new TClonesArray("BmnTof1Digit");
	TClonesArray * T0 = new TClonesArray("BmnTrigDigit");
	TClonesArray * BC3 = new TClonesArray("BmnTrigWaveDigit");
	
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
	rootTree->SetBranchAddress("BC2", &T0);
	rootTree->SetBranchAddress("TQDC_BC3", &BC3);


	// For ToF700 plotting:
	const int fNStrTOF700 = 32;
	const int NPlaneTOF700 = 24;
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
	double fTime700,fBC3Amp,fBC3Time,fT0;

	TH2D*	hTof700_vs_BC3=new TH2D("hTof700_vs_BC3","hTof700_vs_BC3",600,0,3000,4000,-2000,2000);
	TH1D*   hTof700=new TH1D("hTof700","hTof700",1000,0,4000);
	TH1D*   hTime700=new TH1D("hTime700","hTime700",1000,0,4000);
	TH2D*   hTof_vs_Amp_BC3=new TH2D("hTof_vs_Amp_BC3","hTof_vs_Amp_BC3",600,0,3000,1000,-4000,4000);

	for (Int_t j = startEvent; j < stopEvent; j++){
		if (j%100==0) cout << "\tWorking on entry " << j << "\n";

		int nSingleHits = 0;
		EventHead->Clear();

		ToF700->Clear();
		ToF400->Clear();
		T0->Clear();
		BC3->Clear();

		rootTree->GetEntry(j);

	
		if(T0==NULL){ cout<<"NULL"<<endl; continue;}
		int nT0Digits = T0->GetEntriesFast();
		int nBC3Digits = BC3->GetEntriesFast();
            	if (nT0Digits == 1&&nBC3Digits>0) { // T0 digit should be
	//		cout<<"T0"<<endl;
			BmnTrigDigit* digT0 = (BmnTrigDigit*) T0->At(0);
			fT0=digT0->GetTime();
	//		cout<<"BC3"<<endl;
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
		//cout<<"BC3 "<<fBC3Amp<<endl;	
	//		cout<<"TOF700"<<endl;
			for (Int_t iDig = 0; iDig < ToF700->GetEntriesFast(); ++iDig) {
				BmnTof2Digit* digTof700 = (BmnTof2Digit*) ToF700->At(iDig);
                        	
				//pDetector[digTof->GetPlane()]->SetDigit(digTof);
				
				if(digTof700->GetTime() == 0) continue;
				fTime700=digTof700->GetTime();
				hTof700_vs_BC3->Fill(fBC3Amp,fTime700-fT0);
				//cout<<"T "<<fTime700<<" "<<fTime700-fT0<<endl;	
				hTime700->Fill(fTime700);
				hTof700->Fill(fTime700-fT0);
			}
	//		cout<<"TOF400"<<endl;
			
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
	file->Close();
	//TFile *ptr = gFile;
	//gFile->cd();
        //TList* fList = (TList*) pDetector[1]->GetList(2);
	//fList->Write();


}
