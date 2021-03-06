/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <ctime>
/*
#include <TMath.h>
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2S.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TPad.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "BmnTOF1Detector"
*/
Double_t AmpBD, TimeBD, NHitBD;

Double_t CorrLR[20][48] = {0.};

Bool_t FlagHit = kFALSE;
Int_t XforHit = 0, YforHit = 0;

void
Tof400DigitAnalysis_period6(TString file = "", Int_t nEvForRead = 0, Int_t Periud = 7) {
    TStopwatch timer;
    timer.Start();

    gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
    bmnloadlibs(); // load BmnRoot libraries

    ifstream f_call;
    char line[256];
    /*
    Int_t Pl, St;
    Double_t Temp;
    TString NameCallFile = "Tof400LRcorr.dat";
    f_call.open(NameCallFile.Data());
    f_call.getline(line, 256);
    f_call.getline(line, 256);
    while (!f_call.eof()) {
        f_call >> Pl >> St >> Temp >> CorrLR[Pl][St];
        if (TMath::Abs(Temp - CorrLR[Pl][St]) > 2.) CorrLR[Pl][St] = -11.9766;
    }
    */
    BmnTOF1Detector * Plane[20];

    TString name;
    for (Int_t i = 0; i < 20; i++) {
        name = Form("Plane%d", i);
        Plane[i] = new BmnTOF1Detector(i, 2);
        //Plane[i]->SetCorrLR(CorrLR[i]);
        //Plane[i]->SetCorrLR("Tof400LRcorr.dat");
        Plane[i]->SetCorrSlewing("TOF400_SlewingCorr_Period_6.root");
        Plane[i]->SetGeoFile("../../run/geofile_full.root");
    //    Plane[i]->SetGeoFile("geometry_run6.root");
    }
   cout << "finished geo\n";
    TList *fList = new TList();
    TList *fListTof = new TList();
    TH1I *h_dtBDT0 = new TH1I("h_dtBDT0", "h_dtBDT0", 1024, -24, 24);
    fList->Add(h_dtBDT0);
    TH2F *h_dtBDT0_vs_AmpT0 = new TH2F("h_dtBDT0_vs_AmpT0", "h_dtBDT0_vs_AmpT0", 1024, 0, 48, 1024, -24, 24);
    fList->Add(h_dtBDT0_vs_AmpT0);
    TH2F *h_dtBDT0_vs_AmpBD = new TH2F("h_dtBDT0_vs_AmpBD", "h_dtBDT0_vs_AmpBD", 1024, 0, 48, 1024, -24, 24);
    fList->Add(h_dtBDT0_vs_AmpBD);
    TH2F *h_TimeBD_vs_TimeT0 = new TH2F("h_TimeBD_vs_TimeT0", "h_TimeBD_vs_TimeT0", 680, 250, 420, 680, 250, 420);
    fList->Add(h_TimeBD_vs_TimeT0);
    TH1I *h_dtBDT0_Mod = new TH1I("h_dtBDT0_Mod", "h_dtBDT0_Mod", 1024, -24, 24);
    fList->Add(h_dtBDT0_Mod);
    TH1I *h_NHitBD = new TH1I("h_NHitBD", "h_NHitBD", 40, 0, 40);
    fList->Add(h_NHitBD);

    TH1I *h_dtBC2T0 = new TH1I("h_dtBC2T0", "h_dtBC2T0", 1024, -6., 6.);
    fList->Add(h_dtBC2T0);
    TH2F *h_dtBC2T0_vs_AmpT0 = new TH2F("h_dtBC2T0_vs_AmpT0", "h_dtBC2T0_vs_AmpT0", 1024, 0, 48, 1024, -12., 12.);
    fList->Add(h_dtBC2T0_vs_AmpT0);
    TH2F *h_dtBC2T0_vs_AmpBC2 = new TH2F("h_dtBC2T0_vs_AmpBC2", "h_dtBC2T0_vs_AmpBC2", 1024, 0, 48, 1024, -12., 12.);
    fList->Add(h_dtBC2T0_vs_AmpBC2);
    TH2F *h_TimeBC2_vs_TimeT0 = new TH2F("h_TimeBC2_vs_TimeT0", "h_TimeBC2_vs_TimeT0", 680, 250, 420, 680, 250, 420);
    fList->Add(h_TimeBC2_vs_TimeT0);
    TH1I *h_dtBC2T0_Corr = new TH1I("h_dtBC2T0_Corr", "h_dtBC2T0_Corr", 1024, -6., 6.);
    fList->Add(h_dtBC2T0_Corr);

    TH2I *hCrosing = new TH2I("hCrosing", "hCrosing", 48, 0, 48, 48, 0, 48);
    fList->Add(hCrosing);
    TH1I *hCrosingDx = new TH1I("hCrosingDx", "hCrosingDx", 48, 0., 60.);
    fList->Add(hCrosingDx);
    TH1I *hCrosingDy = new TH1I("hCrosingDy", "hCrosingDy", 120, -30., 30.);
    fList->Add(hCrosingDy);
    TH1I *hCrosingDr = new TH1I("hCrosingDr", "hCrosingDr", 600, 0., 60.);
    fList->Add(hCrosingDr);
    TH1I *hCrosingDt = new TH1I("hCrosingDt", "hCrosingDt", 1024, -12., 12.);
    fList->Add(hCrosingDt);
    TH2I *hCrosingWidth_Width = new TH2I("hCrosingWidth_Width", "hCrosingWidth_Width", 1024, 12, 60, 1024, 12, 60);
    fList->Add(hCrosingWidth_Width);
    TH1I *hCrosingDwidth = new TH1I("hCrosingDwidth", "hCrosingDwidth", 512, -12., 12.);
    fList->Add(hCrosingDwidth);

    TH2F *h_HitRPC = new TH2F("h_HitRPC", "h_HitRPC", 96, -96, 0, 5, 0, 5);
    fListTof->Add(h_HitRPC);

    TH2F *h_XYRPC = new TH2F("h_XYRPC", "h_XYRPC", 240, -150, 150, 120, -75, 75);
    fListTof->Add(h_XYRPC);

    TChain *eveTree = new TChain("cbmsim");
    TString inName = Form("/nica/mpd12/segarrae/data/%s", file.Data());
    //inName = file;
    //TString inName = file;
    cout << "Open file " << inName << endl << endl;
    eveTree->Add(inName);

    TClonesArray *ToF400Digits;
    eveTree->SetBranchAddress("TOF400", &ToF400Digits);

    TClonesArray *T0Digits;
    eveTree->SetBranchAddress("BC2", &T0Digits);

    TClonesArray *VetoDigits;
    eveTree->SetBranchAddress("VC", &VetoDigits);

    //TClonesArray *BDDigits;
    //eveTree->SetBranchAddress("BD", &BDDigits);

    //TClonesArray *BC2Digits;
    //eveTree->SetBranchAddress("BC2", &BC2Digits);

    TClonesArray *EventHeader;
    eveTree->SetBranchAddress("EventHeader", &EventHeader);

    Long64_t nEvents = eveTree->GetEntries();
    if (nEvForRead == 0 || nEvForRead > nEvents) nEvForRead = nEvents;
    cout << "Will be readed " << nEvForRead << " events from " << nEvents << endl;

   cout << "finished setting up of trees and histos\n";
    for (Int_t iEv = 0; iEv < nEvForRead; iEv++) {

        if (iEv % 100 == 0) cout << "EVENT: " << iEv << endl;
        eveTree->GetEntry(iEv);
        if ((T0Digits->GetEntriesFast()) != 0
                ) {
            BmnTrigDigit* digT0 = (BmnTrigDigit*) T0Digits->At(0);
            if (digT0->GetMod()!=0) continue;
	    if (digT0->GetAmp() >= 19.5 && digT0->GetAmp() <= 20.5) {

                BmnTrigDigit* digT0 = (BmnTrigDigit*) T0Digits->At(0);

                //--------------------------- EventHeader --------------------------------------------------
                Int_t iEvDig = EventHeader->GetEntriesFast();
                if (iEvDig != 1) cout << "iEvDig == " << iEvDig << endl;
                BmnEventHeader* digEvent = (BmnEventHeader*) EventHeader->At(0);

                //--------------------------- RPC --------------------------------------------------
	        for (Int_t i = 0; i < 20; i++)
                    Plane[i]->Clear();

                for (Int_t iDig = 0; iDig < ToF400Digits->GetEntriesFast(); ++iDig) {
                    FlagHit = kFALSE;
                    BmnTof1Digit* digTof = (BmnTof1Digit*) ToF400Digits->At(iDig);
                    FlagHit = Plane[digTof->GetPlane()]->SetDigit(digTof);
                    if (FlagHit == kTRUE) {
                        if (digTof->GetPlane() < 5) {
                            XforHit = digTof->GetStrip() * -1. - 48;
                            YforHit = 4 - digTof->GetPlane();
                        } else {
                            XforHit = digTof->GetStrip() * -1.;
                            YforHit = 4 - (digTof->GetPlane() - 5);
                        }
                        h_HitRPC->Fill(XforHit, YforHit);
                    }
                }
                TVector3 XYZ;
                XYZ.SetXYZ(0., 0., 0.);
                Double_t ToF = 0;
                FlagHit = kFALSE;
                for (Int_t i = 0; i < 20; i++) {
                    Plane[i] -> FindHits(digT0);
                    for (Int_t s = 0; s < 47; s++) {
                        XYZ.SetXYZ(0., 0., 0.);
                        ToF = 0;
                        FlagHit = kFALSE;
                        FlagHit = Plane[i]->GetXYZTime(s, &XYZ, &ToF);
                        if (FlagHit == kTRUE) h_XYRPC -> Fill(XYZ.x(), XYZ.y());
                    }
                }
                //CrossingDet(Plane[2], Plane[7], hCrosing, hCrosingDx, hCrosingDy, hCrosingDr, hCrosingDt);

                //*/

            }// end             if (digT0->GetAmp() >= 17.3 && digT0->GetAmp() <= 19.2)
        } // end if ((T0Digits->GetEntriesFast()) == 1 && (VetoDigits->GetEntriesFast()) == 0 && (BDDigits->GetEntriesFast()) >= 2)
    } // end for (Int_t iEv = 0; iEv < nEvForRead; iEv++)

    cout << "End Data Analysis" << endl;

    timer.Stop();
    //cout << "Readed " << nEvForRead << " events" << endl;
    cout << "Time          = " << timer.RealTime() << " s" << endl;
    cout << "Time/Event    = " << timer.RealTime() / (Double_t) nEvForRead * 1000. << " ms/Event" << endl;
    cout << "TimeCPU       = " << timer.CpuTime() << " s" << endl;
    cout << "TimeCPU/Event = " << timer.CpuTime() / (Double_t) nEvForRead * 1000. << " ms/Event" << endl;

    ofstream f_time;
    TString NameTimeFile = "Tof400An_Time.dat";
    //TString NameCallFile = file;
    //Point = NameCallFile.First('.');
    //NameCallFile.Replace(Point, 15, "_LRcorr.dat");
    cout << "Write time to " << NameTimeFile.Data() << endl;
    f_time.open(NameTimeFile.Data(), std::ofstream::out | std::ofstream::app);
    f_time << file << "\t" << timer.RealTime() << "\t" << timer.RealTime() / (Double_t) nEvForRead * 1000. << "\t" << timer.CpuTime() << "\t" << timer.CpuTime() / (Double_t) nEvForRead * 1000. << endl;
    f_time.close();

    timer.Reset();
    timer.Start();
    TString outName = "test.root";
    //Int_t Point = outName.First('.');
    //outName.Replace(Point, 15, "_TofAn6.root");
    //outName = Form("/home/storage/analysis/%s", outName.Data());
    cout << "Save data to " << outName.Data() << endl;
    TFile *fileout = new TFile(outName.Data(), "RECREATE");

    TDirectory *Dir;
    Dir = fileout->mkdir("ToF");
    Dir->cd();
    fListTof->Write();

    TDirectory * DirPlane[20];
    TDirectory * Dir1Plane[20];
    TDirectory * Dir2Plane[20];
    for (Int_t i = 0; i < 20; i++) {
        DirPlane[i] = Dir->mkdir(Plane[i]->GetName());
        DirPlane[i] -> cd();
        Plane[i]->GetList(0)->Write();
        Dir1Plane[i] = DirPlane[i]->mkdir("Detector");
        Dir1Plane[i] -> cd();
        Plane[i]->GetList(1)->Write();
        Dir2Plane[i] = DirPlane[i]->mkdir("ToF");
        Dir2Plane[i] -> cd();
        Plane[i]->GetList(2)->Write();
    }

    fileout->cd();
    fList->Write();
    fileout->Close();

    timer.Stop();

    cout << "Time for write root = " << timer.RealTime() << " s" << endl;

}//end of macros

//----Additional functions ----------------------------------------------------------------

