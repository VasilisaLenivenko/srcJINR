// @(#)bmnroot/mwpc:$Id$


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcHitFinderSRC                                                        //
//                                                                            //
// Implementation of an algorithm developed by                                //
// Vasilisa Lenivenko  and Vladimir Palchik                                   //
// to the BmnRoot software                                                    //
//                                                                            //
// The algorithm serves for searching for hits                                //
// in the Multi Wire Prop. Chambers of the BM@N experiment                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <Rtypes.h>
#include <climits>
#include <vector>
#include "TCanvas.h"

#include "BmnMwpcHitFinderSRC.h"
#include "BmnTrack.h"

static Float_t workTime = 0.0;

using namespace std;
using namespace TMath;

BmnMwpcHitFinderSRC::BmnMwpcHitFinderSRC(Bool_t isExp, Int_t runPeriod) :
  fEventNo(0),
  expData(isExp){
    fRunPeriod        = runPeriod;
    fInputBranchName  = "MWPC";
    fOutputBranchName = "BmnMwpcSegment";
    nInputDigits      = 3;
    nTimeSamples      = 3;
    kBig              = 100;
    kNumPairs         = 1;
    kCh_min           = 2;
    kCh_max           = 4;
    if(fRunPeriod == 7){
      kNumPairs       = 2;
      kCh_min         = 0;
    } 
  }
BmnMwpcHitFinderSRC::~BmnMwpcHitFinderSRC() {
}

InitStatus BmnMwpcHitFinderSRC::Init() {
  if (!expData) return kERROR;
  if (fDebug) cout << " BmnMwpcHitFinderSRC::Init() " << endl;

  FairRootManager* ioman = FairRootManager::Instance();
  fBmnMwpcDigitArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
  if (!fBmnMwpcDigitArray){
    cout<<"BmnMwpcHitFinderSRC::Init(): branch "<<fInputBranchName<<" not found! Task will be deactivated"<<endl;
    SetActive(kFALSE);
    return kERROR;
  }

  fBmnMwpcSegmentsArray = new TClonesArray(fOutputBranchName);
  ioman->Register(fOutputBranchName.Data(), "MWPC", fBmnMwpcSegmentsArray, kTRUE);


  fMwpcGeometrySRC = new BmnMwpcGeometrySRC(fRunPeriod); // period number 7
  kNChambers = fMwpcGeometrySRC -> GetNChambers();
  kNPlanes   = fMwpcGeometrySRC -> GetNPlanes();
  kNWires    = fMwpcGeometrySRC -> GetNWires();
  if (fDebug) printf("C-P-W: %d %d %d\n", kNChambers, kNPlanes, kNWires);

  ChCent= new TVector3[kNChambers];
  ChZ =   new Float_t[kNChambers];
  Zmid = new Float_t[kNChambers];

  for (int i=0; i < kNChambers; ++i){ 
    TH1D *h;
    h = new TH1D(Form("Np_best_Ch%d", i), Form("Np_best_Ch%d", i), 6, 1.0, 7.0); fList.Add(h); hNp_best_Ch.push_back(h);
    h = new TH1D(Form("Nbest_Ch%d",   i), Form("Nbest_Ch%d",   i), 6, 0.0, 6.0); fList.Add(h); hNbest_Ch.push_back(h);
    ChCent[i] = fMwpcGeometrySRC->GetChamberCenter(i);
    ChZ[i]    = ChCent[i].Z();
//    ChCent.push_back(fMwpcGeometrySRC->GetChamberCenter(i));
//    ChZ.push_back(ChCent.at(i).Z());
  }

  for (int i=0; i < 24; ++i){
    TH1D *h;
    h = new TH1D(Form("Time%d", i), Form("Time%d", i), 500, 0., 500.);            fList.Add(h); hTime.push_back(h);
    h = new TH1D(Form("Occupancy%d", i), Form("Occupancy%d", i), 100, 0., 100.);    fList.Add(h); hOccupancy.push_back(h);
  }

// Segment finding cuts
  kMinHits  = 4;// for alignment kMinHits  = 6;
  kChi2_Max = 20.;
  kmaxSeg = 10;
  kChMaxAllWires = 200;
// Constants
  dw        = fMwpcGeometrySRC -> GetWireStep(); //0.25; // [cm] // wires step
  dw_half   = 0.5 * dw;
  sq3       = sqrt(3.);
  sq12      = sqrt(12.);
  sigma     = dw / sq12;
  kMiddlePl = 47.25; // Center of wires plane
// Matrices
  matrA  = new Double_t*[4];
  matrb  = new Double_t*[4];
// Arrays
  kPln             = new Int_t*[kNChambers];
  iw_Ch            = new Int_t*[kNChambers];
  wire_Ch          = new Int_t**[kNChambers];
  xuv_Ch           = new Float_t**[kNChambers];
  Wires_Ch         = new Int_t**[kNChambers];    //Wires_Ch1 = new Int_t*[kNPlanes];
  clust_Ch         = new Int_t**[kNChambers];    // clust_Ch1 = new Int_t*[kNPlanes];
  XVU_Ch           = new Float_t**[kNChambers];  // XVU_Ch1 = new Float_t*[kNPlanes];
  Nhits_Ch         = new Int_t*[kNChambers];     // Nhits_Ch1 = new Int_t[kBig];
  Nseg_Ch          = new Int_t[kNChambers];      // Nseg_Ch1
  Nbest_Ch         = new Int_t[kNChambers];      //Nbest_Ch1 = 0;
  ind_best_Ch      = new Int_t*[kNChambers];
  best_Ch_gl       = new Int_t*[kNChambers];
  Chi2_ndf_Ch      = new Double_t*[kNChambers];
  Chi2_ndf_best_Ch = new Double_t*[kNChambers];
  par_ab_Ch        = new Double_t**[kNChambers]; // par_ab_Ch1 = new Double_t*[4];
  XVU              = new Float_t*[kNChambers];   // XVU1 = new Float_t[kNPlanes];
  XVU_cl           = new Float_t*[kNChambers];   // XVU_cl1 = new Float_t[kNPlanes];
  kZ_loc           = new Float_t*[kNChambers];// kZ1_loc = new Float_t[kNPlanes];
  z_gl             = new Float_t*[kNChambers];// z_gl1 = new Float_t[kNPlanes];
  sigm2            = new Float_t[kNPlanes]; 
  ipl              = new Int_t[kNPlanes];
  z2               = new Float_t[kNPlanes];
  for(Int_t i = 0; i < kNChambers; ++i){// for(Int_t ii = kCh_min; ii < kCh_max; ii++){ 
    if (i== 0 || i== 2) { Zmid[i] = (ChZ[i]     - ChZ[i + 1]) *  0.5;}
    if (i== 1 || i== 3) { Zmid[i] = (ChZ[i - 1] - ChZ[i])     * -0.5;}
    printf("Chamber %d Z: %f Zmid: %f\n", i, ChZ[i], Zmid[i]);
    kPln[i]             = new Int_t[kNPlanes];
    iw_Ch[i]            = new Int_t[kNPlanes]; 
    kZ_loc[i]           = new Float_t[kNPlanes];
    z_gl[i]             = new Float_t[kNPlanes];
    Nhits_Ch[i]         = new Int_t[kBig];
    wire_Ch[i]          = new Int_t*[kNWires];
    xuv_Ch[i]           = new Float_t*[kNWires];
    Wires_Ch[i]         = new Int_t*[kNPlanes];
    clust_Ch[i]         = new Int_t*[kNPlanes];
    XVU_Ch[i]           = new Float_t*[kNPlanes];
    par_ab_Ch[i]        = new Double_t*[4];
    XVU[i]              = new Float_t[kNPlanes];
    XVU_cl[i]           = new Float_t[kNPlanes];
    ind_best_Ch[i]      = new Int_t[kmaxSeg];//ind_best_Ch1 = new Int_t[5];
    best_Ch_gl[i]       = new Int_t[kmaxSeg];
    Chi2_ndf_Ch[i]      = new Double_t[kBig];
    Chi2_ndf_best_Ch[i] = new Double_t[kmaxSeg];
    for(int iWire = 0; iWire < kNWires; iWire++){
      wire_Ch[i][iWire] = new Int_t[kNPlanes];
      xuv_Ch[i][iWire]  = new Float_t[kNPlanes];
    }
    for(int iPla = 0; iPla < kNPlanes; ++iPla){
      Wires_Ch[i][iPla] = new Int_t[kBig];
      clust_Ch[i][iPla] = new Int_t[kBig];
      XVU_Ch[i][iPla]   = new Float_t[kBig];
      if ( i == 2 || i == 3){
        kZ_loc[i][iPla] = -0.5 + iPla;
        if(iPla == 4) { kZ_loc[i][iPla] = -2.5;}
        if(iPla == 5) { kZ_loc[i][iPla] = -1.5;}
      }
    }
    for(int ii = 0; ii < 4; ++ii){                 // 4 parameters: tan(x), tan(y), x ,y 
      par_ab_Ch[i][ii] = new Double_t[kBig];
      matrA[ii]        = new Double_t[4];
      matrb[ii]        = new Double_t[4];         
    }
  }

  if (fRunPeriod == 6){
    fill( &kPln[0][0], &kPln[0][0] + sizeof(kPln[0][5]), 0);
    fill( &kPln[1][0], &kPln[1][0] + sizeof(kPln[1][5]), 0);
    fill( &kZ_loc[0][0], &kZ_loc[0][0] + sizeof(kZ_loc[0][5]), 0);
    fill( &kZ_loc[1][0], &kZ_loc[1][0] + sizeof(kZ_loc[1][5]), 0);
   
  }
  else{              // Run 7 for example
    kPln[0][0] = 5;  kZ_loc[0][0] = -1.5;
    kPln[0][1] = 0;  kZ_loc[0][1] = -0.5;
    kPln[0][2] = 1;  kZ_loc[0][2] =  0.5;
    kPln[0][3] = 2;  kZ_loc[0][3] =  1.5;
    kPln[0][4] = 3;  kZ_loc[0][4] =  2.5;
    kPln[0][5] = 4;  kZ_loc[0][5] = -2.5;

    kPln[1][0] = 1;  kZ_loc[1][0] = -1.5;
    kPln[1][1] = 0;  kZ_loc[1][1] = -2.5;
    kPln[1][2] = 5;  kZ_loc[1][2] =  2.5;
    kPln[1][3] = 4;  kZ_loc[1][3] =  1.5;
    kPln[1][4] = 3;  kZ_loc[1][4] =  0.5;
    kPln[1][5] = 2;  kZ_loc[1][5] = -0.5;
}
  kPln[2][0] = 4;
  kPln[2][1] = 5;
  kPln[2][2] = 0;
  kPln[2][3] = 1;
  kPln[2][4] = 2;
  kPln[2][5] = 3;//{4,5,0,1,2,3,  7,11,6,10,9,8,  0,0,0,0,0,0,  0,0,0,0,0,0};

  kPln[3][0] = 4;//1 // run7 like as run5
  kPln[3][1] = 5;
  kPln[3][2] = 0;
  kPln[3][3] = 1;//4
  kPln[3][4] = 2;
  kPln[3][5] = 3;

  if (fRunPeriod == 6){
    kPln[3][0] = 1;//run6-II
    kPln[3][3] = 4;//
  }

// if (fDebug) printf("Run: %d, kPln[0][5]: %d\n", fRunPeriod, kPln[0][5]);

//  if (fDebug) printf("Chamber  Plane  kZ_loc   z_gl\n");
  for(Int_t i = 0; i < kNChambers; ++i){// for(Int_t i = kCh_min; i < kCh_max; i++){
    for(int ii = 0; ii < kNPlanes; ii++){
      z_gl[i][ii] =  Zmid[i] + kZ_loc[i][ii];
      //if (fDebug) printf("%5d  %5d %8.4f  %8.4f\n", i, ii, kZ_loc[i][ii], z_gl[i][ii]);
    }
  }


  return kSUCCESS;
}



void BmnMwpcHitFinderSRC::Exec(Option_t* opt) {
  if (!IsActive()) return;
  clock_t tStart = clock();
  PrepareArraysToProcessEvent();


  if (fDebug) cout << "\n======================== MWPC hit finder exec started =====================\n" << endl;
  printf("Event number: %d\n", fEventNo++);

  Short_t st, wn, pn, pl;
  UInt_t ts;
  for (Int_t iDigit = 0; iDigit < fBmnMwpcDigitArray -> GetEntries(); iDigit++) {
    BmnMwpcDigit* digit = (BmnMwpcDigit*) fBmnMwpcDigitArray ->At (iDigit);
    st = digit -> GetStation();
    wn = digit -> GetWireNumber();
    pl = digit -> GetPlane();
    ts = digit -> GetTime(); //ns
    if ( wn < 0 || wn > 95) cout<<" st "<<st<<" wn = "<<wn<<", pl = "<<pl<<", ts = "<<ts<<endl;

    Int_t ind = st*6 + pl;
    hTime.at(ind) -> Fill(ts);
    hOccupancy.at(ind) -> Fill(wn);
    if ( ts < 35 || ts > 365 ) continue;
    pn = kPln[st][pl];

//  Loop over repeated wires
///*
     Bool_t repeat = 0; 
     if (iw_Ch[st][pn] > 0) {
       for (Int_t ix = 0; ix < iw_Ch[st][pn]; ix++) {
         if (wn == wire_Ch[st][ ix ][pn]  ) {
          repeat = 1;
          break;
         }
        }//ix
     }
     if (repeat) continue;
     
     //  printf("%d %d %d %d\n", st, iw_Ch[st][pn], pn, wn);
     //  */
     // if (fDebug) printf("%d %d %d %d\n", st, iw_Ch[st][pn], pn, wn);
    if (iw_Ch[st][pn] >= 80) continue;
    //if (fEventNo == 2618) continue;
    wire_Ch[st][iw_Ch[st][pn]][pn] = wn;
    xuv_Ch[st][ iw_Ch[st][pn]  ][pn] = (wn - kMiddlePl) * dw;
    //printf("Digit: %d %d %d %d %8.4f --- %d\n", st, pn, iw_Ch[st][pn], wire_Ch[ st ][ iw_Ch[st][pn] ][ pn ], xuv_Ch[st][ iw_Ch[st][pn]  ][pn], fEventNo);
    //cout<<"iw_Ch["<<st<<"]["<<pn<<"]= "<<iw_Ch[st][pn]<<" wire_Ch "<<wire_Ch[ st ][ iw_Ch[st][pn] ][ pn ]<<" xuv_Ch "<<xuv_Ch[st][ iw_Ch[st][pn]  ][pn]<<endl;
    if (pn == 0 || pn == 1 || pn == 5 ) xuv_Ch[st][iw_Ch[st][pn]][pn]= -xuv_Ch[st][iw_Ch[st][pn]][pn]; 

    iw_Ch[st][pn]++;
    //printf("Digit: %d %d %d %d %d\n", st, wn, pl, ts, pn);

  }// iDigit

  //if (fDebug) printf("PS input: chNum sigma  dw/2 z_loc Minhuts Nseg Nhits wiresCh  Nbest ind_best Chi2ndf Chi2ndfbest parab nPln ipl XVU XVU_cl Chi2Max\n");
  
  for (Int_t iChamber = kCh_min; iChamber < kCh_max; iChamber++) {
    Int_t counter = 0;
    for (Int_t iplane = 0; iplane < kNPlanes; iplane++) {
      counter += iw_Ch[iChamber][iplane];
    }
    if (counter < kMinHits || counter > kChMaxAllWires ) continue;

    for(Int_t iCase= 1; iCase < 9; iCase ++){
    // if ( iChamber == 3)
      SegmentFinder(iChamber, Wires_Ch, clust_Ch, XVU_Ch, Nhits_Ch, iw_Ch, Nseg_Ch, wire_Ch, xuv_Ch, kMinHits, iCase, kBig);
    }

// Print segment hits
/*
    printf("SegmentFinder: %d\n", Nseg_Ch[iChamber]);
    for (int ise = 0; ise < Nseg_Ch[iChamber]; ++ise){
      printf("    %d\n", Nhits_Ch[iChamber][ise]);
    }
*/

    if (Nseg_Ch[iChamber] > 0){
      //if (fDebug)  printf("SegmentFinder: %d\n", Nseg_Ch[iChamber]);
      ProcessSegments(iChamber,    sigma,            dw_half,   kZ_loc, 
                      kMinHits,    Nseg_Ch,          Nhits_Ch,  Wires_Ch,  
                      clust_Ch,    XVU_Ch,           Nbest_Ch,  ind_best_Ch, 
                      Chi2_ndf_Ch, Chi2_ndf_best_Ch, par_ab_Ch, kNPlanes, 
                      ipl,         XVU,              XVU_cl,    kChi2_Max);
      hNbest_Ch.at(iChamber) -> Fill(Nbest_Ch[iChamber]);

       //	  cout<<"ProcessSegments: Nbest_Ch["<<iChamber<<"]= "<<Nbest_Ch[iChamber]<<endl;
    }

  }//iChamber

  for (Int_t iChamber = kCh_min; iChamber < kCh_max; iChamber++) {
    for (Int_t ise = 0; ise < Nbest_Ch[iChamber]; ise++) {
      if (Nhits_Ch[iChamber][ ind_best_Ch[iChamber][ise] ] > 3) {
  //   printf("|%6s|%6s|%8s|%5s|\n","Ch.ID", "Pl.ID", "Coord.", "Wires" );
  //    printf("------------------------------\n");
        for (Int_t i = 0; i < 6; i++){
          if ( Wires_Ch[iChamber][i][ind_best_Ch[iChamber][ise]] == -1 ) XVU_Ch[iChamber][i][ind_best_Ch[iChamber][ise]] = -999.;//WARNING!!!
          //  printf("|%6d|%6d|%8.2f|%5d|\n", iChamber, i, XVU_Ch[iChamber][i][ind_best_Ch[iChamber][ise]], Wires_Ch[iChamber][i][ind_best_Ch[iChamber][ise]]);
          //	        cout<<" Ch= "<<iChamber<<" i "<<i<<" XVU_Ch "<<XVU_Ch[iChamber][i][ind_best_Ch[iChamber][ise]]<<" Wires "<<Wires_Ch[iChamber][i][ind_best_Ch[iChamber][ise]] <<endl;      
        }//i6
	/*
	  if (Nhits_Ch[iChamber][ ind_best_Ch[iChamber][ise] ] > 4 && iChamber == 3 )   {
	  cout<<" Ch= "<<iChamber<<" ise "<<ise<<" ind "<<ind_best_Ch[iChamber][ise]<<" Chi2 "<<Chi2_ndf_best_Ch[iChamber][ise]<<" Ax "<<par_ab_Ch[iChamber][0][ind_best_Ch[iChamber][ise]]<<" bx "<<par_ab_Ch[iChamber][1][ind_best_Ch[iChamber][ise]]<<" Ay "<<par_ab_Ch[iChamber][2][ind_best_Ch[iChamber][ise]]<<" by "<<par_ab_Ch[iChamber][3][ind_best_Ch[iChamber][ise]]<<endl;

	  Double_t Um , Xm,  Vm, Vp, Xp, Up;

	  if ( iChamber == 0 ){
	  Um = XVU_Ch[iChamber][5][ind_best_Ch[iChamber][ise]];
	  Xm = XVU_Ch[iChamber][0][ind_best_Ch[iChamber][ise]];
	  Vm = XVU_Ch[iChamber][1][ind_best_Ch[iChamber][ise]];
	  Vp = XVU_Ch[iChamber][2][ind_best_Ch[iChamber][ise]];
	  Xp = XVU_Ch[iChamber][3][ind_best_Ch[iChamber][ise]];
	  Up = XVU_Ch[iChamber][4][ind_best_Ch[iChamber][ise]];
	  }
	  if ( iChamber == 1 ){
	  Vm = XVU_Ch[iChamber][1][ind_best_Ch[iChamber][ise]];
	  Xm = XVU_Ch[iChamber][0][ind_best_Ch[iChamber][ise]];
	  Um = XVU_Ch[iChamber][5][ind_best_Ch[iChamber][ise]];
	  Vp = XVU_Ch[iChamber][4][ind_best_Ch[iChamber][ise]];
	  Xp = XVU_Ch[iChamber][3][ind_best_Ch[iChamber][ise]];
	  Up = XVU_Ch[iChamber][2][ind_best_Ch[iChamber][ise]];
	  }
	  if ( iChamber == 2 ){
	  Vp = XVU_Ch[iChamber][4][ind_best_Ch[iChamber][ise]];
	  Um = XVU_Ch[iChamber][5][ind_best_Ch[iChamber][ise]];
	  Xm = XVU_Ch[iChamber][0][ind_best_Ch[iChamber][ise]];
	  Vm = XVU_Ch[iChamber][1][ind_best_Ch[iChamber][ise]];
	  Up = XVU_Ch[iChamber][2][ind_best_Ch[iChamber][ise]];
	  Xp = XVU_Ch[iChamber][3][ind_best_Ch[iChamber][ise]];
	  }
	  if ( iChamber == 3 ){
	  Vm = XVU_Ch[iChamber][4][ind_best_Ch[iChamber][ise]];
	  Um = XVU_Ch[iChamber][5][ind_best_Ch[iChamber][ise]];
	  Xm = XVU_Ch[iChamber][0][ind_best_Ch[iChamber][ise]];
	  Vp = XVU_Ch[iChamber][1][ind_best_Ch[iChamber][ise]];
	  Up = XVU_Ch[iChamber][2][ind_best_Ch[iChamber][ise]];
	  Xp = XVU_Ch[iChamber][3][ind_best_Ch[iChamber][ise]];
	  }

	  if ( iChamber == 0 ||  iChamber == 1 )
	  cout<<" Chamber "<<iChamber<<" --Ax_x "<<(Xp - Xm )/3.<<"--Ax_uv "<<(Up + Vp - Um - Vm)/3.<<"--Ax_xuv "<<(Xp + Up + Vp - Xm - Um - Vm )/6.<<endl;
	  if ( iChamber == 2 ||  iChamber == 3 ) 
	  cout<<" Chamber "<<iChamber<<" --Ax_x "<<(Xp - Xm )/3.<<"--Ax_uv "<<(Up + Vm - Um - Vp)/4.<<"--Ax_xuv "<<(Xp - Xm + Up + Vm - Um - Vp)/3.5<<endl;
	  cout<<endl;
	 
	  for (Int_t i = 0; i < 6; i++){
	  cout<<" "<<i<<" "<< XVU_Ch[iChamber][i][ind_best_Ch[iChamber][ise]]<<" Z_lock "<<kZ_loc[iChamber][i]<<endl;
	  }
	  cout<<endl;
	}
*/
  //  cout<<" Nhits_Ch "<<Nhits_Ch[iChamber][ ind_best_Ch[iChamber][ise] ]<<endl;
  //  printf("hNp_best_Ch filling; Chamber: %d \n", iChamber);
        hNp_best_Ch.at(iChamber) -> Fill(Nhits_Ch[iChamber][ ind_best_Ch[iChamber][ise] ]);

        BmnTrack *pSeg = new ((*fBmnMwpcSegmentsArray)[fBmnMwpcSegmentsArray->GetEntriesFast()]) BmnTrack();
        pSeg->SetChi2(Chi2_ndf_best_Ch[iChamber][ise]);
        pSeg->SetNHits(Nhits_Ch[iChamber][ ind_best_Ch[iChamber][ise] ]);
              pSeg->SetFlag(ise);
        FairTrackParam pSegParams;
        pSegParams.SetPosition(TVector3(par_ab_Ch[iChamber][1][ind_best_Ch[iChamber][ise]], par_ab_Ch[iChamber][3][ind_best_Ch[iChamber][ise]],ChZ[iChamber])); // pSegParams.SetPosition(TVector3(par_ab_Ch[iChamber][1][ise], par_ab_Ch[iChamber][3][ise],ZCh[iChamber]));
        pSegParams.SetTx(par_ab_Ch[iChamber][0][ind_best_Ch[iChamber][ise]]);//pSegParams.SetTx(par_ab_Ch[iChamber][0][ise]);
        pSegParams.SetTy(par_ab_Ch[iChamber][2][ind_best_Ch[iChamber][ise]]);// pSegParams.SetTy(par_ab_Ch[iChamber][2][ise]);
        pSeg->SetParamFirst(pSegParams);
        FairTrackParam pSegParamsForXUV;
        pSegParamsForXUV.SetX(XVU_Ch[iChamber][0][ind_best_Ch[iChamber][ise]]);
        pSegParamsForXUV.SetY(XVU_Ch[iChamber][1][ind_best_Ch[iChamber][ise]]);
        pSegParamsForXUV.SetZ(XVU_Ch[iChamber][2][ind_best_Ch[iChamber][ise]]);
        pSegParamsForXUV.SetTx(XVU_Ch[iChamber][3][ind_best_Ch[iChamber][ise]]);
        pSegParamsForXUV.SetTy(XVU_Ch[iChamber][4][ind_best_Ch[iChamber][ise]]);
        pSegParamsForXUV.SetQp(XVU_Ch[iChamber][5][ind_best_Ch[iChamber][ise]]);
        pSeg->SetParamLast(pSegParamsForXUV);
      }//ise
    }//if(Nseg_Ch > 0)
  }//[iChamber]

   if (fDebug)  cout << "\n======================== MWPC hit finder exec finished ====================" << endl;
  clock_t tFinish = clock();
  workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
 }




void BmnMwpcHitFinderSRC::SegmentFinder(Int_t chNum, Int_t*** wires_Ch, Int_t ***clust_Ch_, Float_t ***XVU_Ch_, Int_t **Nhits_Ch_, Int_t **iw_Ch_, Int_t *Nseg, Int_t ***wires_glob, Float_t ***xuv_glob, Int_t minHits, Short_t code, Int_t kBig_ ) {
  //xuv_glob     - coordinates of all hits
  //wires_glob   - wires of all hits
  Float_t delta      = 1.125;//Float_t delta = (chNum == 3) ? 1.125 : 0.75; //0.5; //condition for the nearest "lighted" wires
  Int_t cntr_sumWW_x =  95 ;//center wires Wx1+Wx2 due to x-slope
  Int_t min_sumWW_x  =  3;//Int_t min_sumWW_x =  (chNum == 3) ? 3 : 2; //min sum of wires to 95 //tmp
  Int_t min_sumWW    =   3;// Int_t min_sumWW =  (chNum == 3) ? 3 : 2; //min sum of wires to 95  //wide UV  //tmp
  Int_t cntr_sumWW   =  95 ;//center wires W1+W2 due to u,v-slope
  //Int_t min_sumWW  =  2 ; //min sum of wires to 95 //narrow UV
  Int_t minHits4_5   = minHits;

  // code :
  // 1  {X-, V-, U+}
  // 2  {X-, V+, U+}
  // 3  {X-, V-, U-}
  // 7  {X+, V-, U+}
  // 5  {X+, V+, U+}
  // 6  {X+, V-, U-}
  // 4  {X-, V+, U-}
  // 8  {X+, V+, U-}

  Int_t x = 0, v = 1, u = 2 , x1 = 3, v1 = 4, u1 = 5;//MK

  switch (code) {
  case 1:
    x = 0;
    v = 1;
    u = 2;
    x1 = 3;
    v1 = 4;
    u1 = 5;
    break;
  case 2:
    x = 0;
    v = 4;
    u = 2;
    x1 = 3;
    v1 = 1;
    u1 = 5;
    break;
  case 3:
    x = 0;
    v = 1;
    u = 5;
    x1 = 3;
    v1 = 4;
    u1 = 2;
    break;
  case 7:
    x = 3;
    v = 1;
    u = 2;
    x1 = 0;
    v1 = 4;
    u1 = 5;
    break;
  case 5:
    x = 3;
    v = 4;
    u = 2;
    x1 = 0;
    v1 = 1;
    u1 = 5;
    break;
  case 6:
    x = 3;
    v = 1;
    u = 5;
    x1 = 0;
    v1 = 4;
    u1 = 2;
    break;
  case 4:
    x = 0;
    v = 4;
    u = 5;
    x1 = 3;
    v1 = 1;
    u1 = 2;
    break;
  case 8:
    x = 3;
    v = 4;
    u = 5;
    x1 = 0;
    v1 = 1;
    u1 = 2;
    break;
  }

  if (Nseg[chNum] > kBig_ - 2) return;// MP
  if (iw_Ch_[chNum][x] > 0) {
    for (Int_t ix = 0; ix < iw_Ch_[chNum][x]; ix++) {
      if (iw_Ch_[chNum][v] > 0){
        for (Int_t iv = 0; iv < iw_Ch_[chNum][v]; iv++) {
          if (iw_Ch_[chNum][u] > 0) {
            for (Int_t iu = 0; iu < iw_Ch_[chNum][u]; iu++) {
              if (Nseg[chNum] > kBig_ - 2) return;
              Bool_t it_was = 0;
              if (Nseg[chNum] > 0) {
                for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
                  Bool_t it_was_x = 0;
                  Bool_t it_was_v = 0;
                  Bool_t it_was_u = 0;
  // ::out1<<" ch "<<chNum<<" case "<<code<<" iseg "<<iseg<<" ix "<<x<<" wr "<<wires_Ch[chNum][x][iseg]<<" iu "<<u<<" wr "<<wires_Ch[chNum][u][iseg]<<" iv "<<v<<" wr "<<wires_Ch[chNum][v][iseg]<<endl;
                  if (wires_Ch[chNum][x][iseg]   == wires_glob[chNum][ix][x]) it_was_x = 1;
                  if ((clust_Ch_[chNum][x][iseg] == -1 || clust_Ch_[chNum][x][iseg] == 3) && (wires_Ch[chNum][x][iseg] -1 == wires_glob[chNum][ix][x])) it_was_x = 1;
                  if ((clust_Ch_[chNum][x][iseg] ==  1 || clust_Ch_[chNum][x][iseg] == 3) && (wires_Ch[chNum][x][iseg] +1 == wires_glob[chNum][ix][x])) it_was_x = 1;
                  if (wires_Ch[chNum][v][iseg]   == wires_glob[chNum][iv][v]) it_was_v = 1;
                  if ((clust_Ch_[chNum][v][iseg] == -1 || clust_Ch_[chNum][v][iseg] == 3) && (wires_Ch[chNum][v][iseg] -1 == wires_glob[chNum][iv][v])) it_was_v = 1;
                  if ((clust_Ch_[chNum][v][iseg] ==  1 || clust_Ch_[chNum][v][iseg] == 3) && (wires_Ch[chNum][v][iseg] +1 == wires_glob[chNum][iv][v])) it_was_v = 1;
                  if (wires_Ch[chNum][u][iseg]   == wires_glob[chNum][iu][u]) it_was_u = 1;
                  if ((clust_Ch_[chNum][u][iseg] == -1 || clust_Ch_[chNum][u][iseg] == 3) && (wires_Ch[chNum][u][iseg] -1 == wires_glob[chNum][iu][u])) it_was_u = 1;
                  if ((clust_Ch_[chNum][u][iseg] ==  1 || clust_Ch_[chNum][u][iseg] == 3) && (wires_Ch[chNum][u][iseg] +1 == wires_glob[chNum][iu][u])) it_was_u = 1;
  //   bef.clustering    if (wires_Ch[chNum][x][iseg] == wires_glob[chNum][ix][x] && wires_Ch[chNum][v][iseg] == wires_glob[chNum][iv][v] && wires_Ch[chNum][u][iseg] == wires_glob[chNum][iu][u]) {
                  if (it_was_x * it_was_v * it_was_u)  { it_was = 1; break; }
                }  // iseg
              }  //  Nseg[chNum] > 0
              if (it_was) continue;
              if (fabs(xuv_glob[chNum][iu][u] + xuv_glob[chNum][iv][v] - xuv_glob[chNum][ix][x]) < delta) {
                //  3p-candidate new Nseg
                XVU_Ch_[chNum][x][Nseg[chNum]] = xuv_glob[chNum][ix][x];
                XVU_Ch_[chNum][v][Nseg[chNum]] = xuv_glob[chNum][iv][v];
                XVU_Ch_[chNum][u][Nseg[chNum]] = xuv_glob[chNum][iu][u];
                wires_Ch[chNum][x1][Nseg[chNum]] = -1;
                wires_Ch[chNum][v1][Nseg[chNum]] = -1;
                wires_Ch[chNum][u1][Nseg[chNum]] = -1;
                wires_Ch[chNum][x][Nseg[chNum]] = wires_glob[chNum][ix][x];
                wires_Ch[chNum][v][Nseg[chNum]] = wires_glob[chNum][iv][v];
                wires_Ch[chNum][u][Nseg[chNum]] = wires_glob[chNum][iu][u];
                Nhits_Ch_[chNum][Nseg[chNum]] = 3;
                //clustering
                clust_Ch_[chNum][x][Nseg[chNum]]=0;
                clust_Ch_[chNum][v][Nseg[chNum]]=0;
                clust_Ch_[chNum][u][Nseg[chNum]]=0;

                if ( Nseg[chNum] < 20 ) {
                  Int_t below = 0;
                  Int_t above = 0;

                  for (Int_t icl = 0; icl < iw_Ch_[chNum][u]; icl++) {
                    if ( wires_glob[chNum][icl][u] ==  wires_Ch[chNum][u][Nseg[chNum]] -1 ) below = 1;
                    if ( wires_glob[chNum][icl][u] ==  wires_Ch[chNum][u][Nseg[chNum]] +1 ) above = 1;
                  }
                  if (below == 1 && above == 0) 	clust_Ch_[chNum][u][Nseg[chNum]]= -1;
                  if (below == 0 && above == 1) 	clust_Ch_[chNum][u][Nseg[chNum]]= +1;
                  if (below == 1 && above == 1) 	clust_Ch_[chNum][u][Nseg[chNum]]=  3;

                  below = 0;
                  above = 0;
                  for (Int_t icl = 0; icl < iw_Ch_[chNum][x]; icl++) {
                    if ( wires_glob[chNum][icl][x] ==  wires_Ch[chNum][x][Nseg[chNum]] -1 ) below = 1;
                    if ( wires_glob[chNum][icl][x] ==  wires_Ch[chNum][x][Nseg[chNum]] +1 ) above = 1;
                  }
                  if (below == 1 && above == 0) 	clust_Ch_[chNum][x][Nseg[chNum]]= -1;
                  if (below == 0 && above == 1) 	clust_Ch_[chNum][x][Nseg[chNum]]= +1;
                  if (below == 1 && above == 1) 	clust_Ch_[chNum][x][Nseg[chNum]]=  3;

                  below = 0;
                  above = 0;
                  for (Int_t icl = 0; icl < iw_Ch_[chNum][v]; icl++) {
                    if ( wires_glob[chNum][icl][v] ==  wires_Ch[chNum][v][Nseg[chNum]] -1 ) below = 1;
                    if ( wires_glob[chNum][icl][v] ==  wires_Ch[chNum][v][Nseg[chNum]] +1 ) above = 1;
                  }
                  if (below == 1 && above == 0) 	clust_Ch_[chNum][v][Nseg[chNum]]= -1;
                  if (below == 0 && above == 1) 	clust_Ch_[chNum][v][Nseg[chNum]]= +1;
                  if (below == 1 && above == 1) 	clust_Ch_[chNum][v][Nseg[chNum]]=  3;
                }	//	if ( Nseg[chNum] < 20 ) 

                Int_t double_u2 = 0;

                if (iw_Ch_[chNum][u1] > 0) {
                  for (Int_t iu2 = 0; iu2 < iw_Ch_[chNum][u1]; iu2++) {
                    Bool_t ifl = 0;
                    if (abs( wires_Ch[chNum][u][Nseg[chNum]] + wires_glob[chNum][iu2][u1] - cntr_sumWW) <  min_sumWW ) {
                      XVU_Ch_[chNum][u1][Nseg[chNum]] = xuv_glob[chNum][iu2][u1];
                      wires_Ch[chNum][u1][Nseg[chNum]] = wires_glob[chNum][iu2][u1];
                      Nhits_Ch_[chNum][Nseg[chNum]] = 4;
                      clust_Ch_[chNum][u1][Nseg[chNum]]=0;
                      if ( Nseg[chNum] < 20 ) {
                        Int_t below = 0;
                        Int_t above = 0;
                        for (Int_t icl = 0; icl < iw_Ch_[chNum][u1]; icl++) {
                          if ( wires_glob[chNum][icl][u1] ==  wires_Ch[chNum][u1][Nseg[chNum]] -1 ) below = 1;
                          if ( wires_glob[chNum][icl][u1] ==  wires_Ch[chNum][u1][Nseg[chNum]] +1 ) above = 1;
                          if ( wires_glob[chNum][icl][u1] ==  wires_Ch[chNum][u1][Nseg[chNum]] +2 ) {
                            if (abs( wires_Ch[chNum][u][Nseg[chNum]] + wires_glob[chNum][icl][u1] - cntr_sumWW) <  min_sumWW ) double_u2 = icl;			    
                          }
                        }
                        if (below == 1 && above == 0) 	clust_Ch_[chNum][u1][Nseg[chNum]]= -1;
                        if (below == 0 && above == 1) 	clust_Ch_[chNum][u1][Nseg[chNum]]= +1;
                        if (below == 1 && above == 1) 	clust_Ch_[chNum][u1][Nseg[chNum]]=  3;
                      }	//	if ( Nseg[chNum] < 20 ) 
                      ifl = 1;
                    }
                    if ( ifl ) break;
                  }//iu2
                }//u1

                Int_t double_v2 = 0;

                if (iw_Ch_[chNum][v1] > 0) {
                  for (Int_t iv2 = 0; iv2 < iw_Ch_[chNum][v1]; iv2++) {
                    Bool_t ifl = 0;
                    if (abs( wires_Ch[chNum][v][Nseg[chNum]] + wires_glob[chNum][iv2][v1] - cntr_sumWW) <  min_sumWW ) {
                      XVU_Ch_[chNum][v1][Nseg[chNum]] = xuv_glob[chNum][iv2][v1];
                      wires_Ch[chNum][v1][Nseg[chNum]] = wires_glob[chNum][iv2][v1];

                      if ( wires_Ch[chNum][v1][Nseg[chNum]] > -1 ) {
                        Nhits_Ch_[chNum][Nseg[chNum]] = Nhits_Ch_[chNum][Nseg[chNum]] + 1;//5 points
                        clust_Ch_[chNum][v1][Nseg[chNum]]=0;
                        if ( Nseg[chNum] < 20 ) {
                          Int_t below = 0;
                          Int_t above = 0;
                          for (Int_t icl = 0; icl < iw_Ch_[chNum][v1]; icl++) {
                            if ( wires_glob[chNum][icl][v1] ==  wires_Ch[chNum][v1][Nseg[chNum]] -1 ) below = 1;
                            if ( wires_glob[chNum][icl][v1] ==  wires_Ch[chNum][v1][Nseg[chNum]] +1 ) above = 1;
                            if ( wires_glob[chNum][icl][v1] ==  wires_Ch[chNum][v1][Nseg[chNum]] +2 ) {
                              if (abs( wires_Ch[chNum][v][Nseg[chNum]] + wires_glob[chNum][icl][v1] - cntr_sumWW) <  min_sumWW ) double_v2 = icl;			    
                            }
                          }
                          if (below == 1 && above == 0) 	clust_Ch_[chNum][v1][Nseg[chNum]]= -1;
                          if (below == 0 && above == 1) 	clust_Ch_[chNum][v1][Nseg[chNum]]= +1;
                          if (below == 1 && above == 1) 	clust_Ch_[chNum][v1][Nseg[chNum]]=  3;
                        }	//	if ( Nseg[chNum] < 20 ) 
                        ifl = 1;
                      }
                    }//abs( wires_Ch[chNum][v][Nseg[chNum]] +
                    if ( ifl ) break;
                  }//iv2
                }//v1

                Int_t double_x2 = 0;
                if (iw_Ch_[chNum][x1] > 0) {
                  for (Int_t ix2 = 0; ix2 < iw_Ch_[chNum][x1]; ix2++) {
                    Bool_t ifl = 0;
                    if(abs( wires_Ch[chNum][x][Nseg[chNum]] + wires_glob[chNum][ix2][x1] - cntr_sumWW_x) <  min_sumWW_x ) {
                      XVU_Ch_[chNum][x1][Nseg[chNum]] = xuv_glob[chNum][ix2][x1];
                      wires_Ch[chNum][x1][Nseg[chNum]] = wires_glob[chNum][ix2][x1];
                                                    
                      if ( wires_Ch[chNum][x1][Nseg[chNum]] > -1 ) {
                        Nhits_Ch_[chNum][Nseg[chNum]] = Nhits_Ch_[chNum][Nseg[chNum]] + 1;//6 points
                        clust_Ch_[chNum][x1][Nseg[chNum]]=0;

                        if ( Nseg[chNum] < 20 ) {
                          Int_t below = 0;
                          Int_t above = 0;
                          for (Int_t icl = 0; icl < iw_Ch_[chNum][x1]; icl++) {
                            if ( wires_glob[chNum][icl][x1] ==  wires_Ch[chNum][x1][Nseg[chNum]] -1 ) below = 1;
                            if ( wires_glob[chNum][icl][x1] ==  wires_Ch[chNum][x1][Nseg[chNum]] +1 ) above = 1;
                            if ( wires_glob[chNum][icl][x1] ==  wires_Ch[chNum][x1][Nseg[chNum]] +2 ) {
                              if (abs( wires_Ch[chNum][x][Nseg[chNum]] + wires_glob[chNum][icl][x1] - cntr_sumWW_x) <  min_sumWW_x ) double_x2 = icl;			    
                            }
                          }
                          if (below == 1 && above == 0) 	clust_Ch_[chNum][x1][Nseg[chNum]]= -1;
                          if (below == 0 && above == 1) 	clust_Ch_[chNum][x1][Nseg[chNum]]= +1;
                          if (below == 1 && above == 1) 	clust_Ch_[chNum][x1][Nseg[chNum]]=  3;
                        }	//	if ( Nseg[chNum] < 20 ) 
                        ifl = 1;
                      }//wires_Ch[chNum][x1][Nseg[chNum]] > -1
                    }//abs( wires_Ch[chNum][x][Nseg[chNum]] +
                    if ( ifl ) break;
                  }//ix2
                }//x1

                if (Nseg[chNum] > 15) minHits4_5=5;
                if (Nseg[chNum] > 30) minHits4_5=6;

                if (Nhits_Ch_[chNum][Nseg[chNum]] >= minHits4_5) {
                  //		  ::out1<<"find iseg "<<Nseg[chNum]<<" wires_Ch[chNum] "<< wires_Ch[chNum][0][Nseg[chNum]]<<" "<< wires_Ch[chNum][1][Nseg[chNum]]<<" "<<wires_Ch[chNum][2][Nseg[chNum]]<<" "<<wires_Ch[chNum][3][Nseg[chNum]]<<" "<<wires_Ch[chNum][4][Nseg[chNum]]<<" "<<wires_Ch[chNum][5][Nseg[chNum]]<<endl;
                  Nseg[chNum]++;
                  if ( double_u2 > 0 && Nseg[chNum] < 99 ) {
                    //   ::out1<<"   u2 double "<<double_u2<<" Wire "<< wires_glob[chNum][double_u2][u1] <<endl;
                    Nhits_Ch_[chNum][Nseg[chNum]] = Nhits_Ch_[chNum][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x][Nseg[chNum]] = XVU_Ch_[chNum][x][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x1][Nseg[chNum]] = XVU_Ch_[chNum][x1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u][Nseg[chNum]] = XVU_Ch_[chNum][u][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v][Nseg[chNum]] = XVU_Ch_[chNum][v][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v1][Nseg[chNum]] = XVU_Ch_[chNum][v1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u1][Nseg[chNum]] = xuv_glob[chNum][double_u2][u1];
                    wires_Ch[chNum][x][Nseg[chNum]] = wires_Ch[chNum][x][Nseg[chNum]-1];
                    wires_Ch[chNum][x1][Nseg[chNum]] = wires_Ch[chNum][x1][Nseg[chNum]-1];
                    wires_Ch[chNum][u][Nseg[chNum]] = wires_Ch[chNum][u][Nseg[chNum]-1];
                    wires_Ch[chNum][v][Nseg[chNum]] = wires_Ch[chNum][v][Nseg[chNum]-1];
                    wires_Ch[chNum][v1][Nseg[chNum]] = wires_Ch[chNum][v1][Nseg[chNum]-1];
                    wires_Ch[chNum][u1][Nseg[chNum]] = wires_glob[chNum][double_u2][u1];
                    //clustering
                    clust_Ch_[chNum][x][Nseg[chNum]] = clust_Ch_[chNum][x][Nseg[chNum]-1];
                    clust_Ch_[chNum][x1][Nseg[chNum]] = clust_Ch_[chNum][x1][Nseg[chNum]-1];
                    clust_Ch_[chNum][u][Nseg[chNum]] = clust_Ch_[chNum][u][Nseg[chNum]-1];
                    clust_Ch_[chNum][v][Nseg[chNum]] = clust_Ch_[chNum][v][Nseg[chNum]-1];
                    clust_Ch_[chNum][v1][Nseg[chNum]] = clust_Ch_[chNum][v1][Nseg[chNum]-1];
                    clust_Ch_[chNum][u1][Nseg[chNum]]=0;
                    //		        ::out1<<"find doub_u2 iseg "<<Nseg[chNum]<<" wires_Ch[chNum] "<< wires_Ch[chNum][0][Nseg[chNum]]<<" "<< wires_Ch[chNum][1][Nseg[chNum]]<<" "<<wires_Ch[chNum][2][Nseg[chNum]]<<" "<<wires_Ch[chNum][3][Nseg[chNum]]<<" "<<wires_Ch[chNum][4][Nseg[chNum]]<<" "<<wires_Ch[chNum][5][Nseg[chNum]]<<endl;
                    Nseg[chNum]++;
                    if (Nseg[chNum] > kBig_ - 2) break;
                  }

                  if ( double_v2 > 0 && Nseg[chNum] < 99 ) { 
                    // ::out1<<"   v2 double "<<double_v2<<" Wire "<< wires_glob[chNum][double_v2][v1] <<endl;
                    Nhits_Ch_[chNum][Nseg[chNum]] = Nhits_Ch_[chNum][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x][Nseg[chNum]] = XVU_Ch_[chNum][x][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x1][Nseg[chNum]] = XVU_Ch_[chNum][x1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u][Nseg[chNum]] = XVU_Ch_[chNum][u][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v][Nseg[chNum]] = XVU_Ch_[chNum][v][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u1][Nseg[chNum]] = XVU_Ch_[chNum][u1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v1][Nseg[chNum]] = xuv_glob[chNum][double_v2][v1];
                    wires_Ch[chNum][x][Nseg[chNum]] = wires_Ch[chNum][x][Nseg[chNum]-1];
                    wires_Ch[chNum][x1][Nseg[chNum]] = wires_Ch[chNum][x1][Nseg[chNum]-1];
                    wires_Ch[chNum][u][Nseg[chNum]] = wires_Ch[chNum][u][Nseg[chNum]-1];
                    wires_Ch[chNum][v][Nseg[chNum]] = wires_Ch[chNum][v][Nseg[chNum]-1];
                    wires_Ch[chNum][u1][Nseg[chNum]] = wires_Ch[chNum][u1][Nseg[chNum]-1];
                    wires_Ch[chNum][v1][Nseg[chNum]] = wires_glob[chNum][double_v2][v1];
                    //clustering
                    clust_Ch_[chNum][x][Nseg[chNum]] = clust_Ch_[chNum][x][Nseg[chNum]-1];
                    clust_Ch_[chNum][x1][Nseg[chNum]] = clust_Ch_[chNum][x1][Nseg[chNum]-1];
                    clust_Ch_[chNum][u][Nseg[chNum]] = clust_Ch_[chNum][u][Nseg[chNum]-1];
                    clust_Ch_[chNum][v][Nseg[chNum]] = clust_Ch_[chNum][v][Nseg[chNum]-1];
                    clust_Ch_[chNum][u1][Nseg[chNum]] = clust_Ch_[chNum][u1][Nseg[chNum]-1];
                    clust_Ch_[chNum][v1][Nseg[chNum]]=0;
                    //		        ::out1<<"find doub_v2 iseg "<<Nseg[chNum]<<" wires_Ch[chNum] "<< wires_Ch[chNum][0][Nseg[chNum]]<<" "<< wires_Ch[chNum][1][Nseg[chNum]]<<" "<<wires_Ch[chNum][2][Nseg[chNum]]<<" "<<wires_Ch[chNum][3][Nseg[chNum]]<<" "<<wires_Ch[chNum][4][Nseg[chNum]]<<" "<<wires_Ch[chNum][5][Nseg[chNum]]<<endl;
                    Nseg[chNum]++;
                    if (Nseg[chNum] > kBig_ - 2) break;
                  }

                  if ( double_x2 > 0 && Nseg[chNum] < 99 ) { 
                    // ::out1<<"   x2 double "<<double_x2<<" Wire "<< wires_glob[chNum][double_x2][x1] <<endl;
                    Nhits_Ch_[chNum][Nseg[chNum]] = Nhits_Ch_[chNum][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x][Nseg[chNum]] = XVU_Ch_[chNum][x][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u1][Nseg[chNum]] = XVU_Ch_[chNum][u1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][u][Nseg[chNum]] = XVU_Ch_[chNum][u][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v][Nseg[chNum]] = XVU_Ch_[chNum][v][Nseg[chNum]-1];
                    XVU_Ch_[chNum][v1][Nseg[chNum]] = XVU_Ch_[chNum][v1][Nseg[chNum]-1];
                    XVU_Ch_[chNum][x1][Nseg[chNum]] = xuv_glob[chNum][double_x2][x1];
                    wires_Ch[chNum][x][Nseg[chNum]] = wires_Ch[chNum][x][Nseg[chNum]-1];
                    wires_Ch[chNum][u1][Nseg[chNum]] = wires_Ch[chNum][u1][Nseg[chNum]-1];
                    wires_Ch[chNum][u][Nseg[chNum]] = wires_Ch[chNum][u][Nseg[chNum]-1];
                    wires_Ch[chNum][v][Nseg[chNum]] = wires_Ch[chNum][v][Nseg[chNum]-1];
                    wires_Ch[chNum][v1][Nseg[chNum]] = wires_Ch[chNum][v1][Nseg[chNum]-1];
                    wires_Ch[chNum][x1][Nseg[chNum]] = wires_glob[chNum][double_x2][x1];
                    //clustering
                    clust_Ch_[chNum][x][Nseg[chNum]] = clust_Ch_[chNum][x][Nseg[chNum]-1];
                    clust_Ch_[chNum][u1][Nseg[chNum]] = clust_Ch_[chNum][u1][Nseg[chNum]-1];
                    clust_Ch_[chNum][u][Nseg[chNum]] = clust_Ch_[chNum][u][Nseg[chNum]-1];
                    clust_Ch_[chNum][v][Nseg[chNum]] = clust_Ch_[chNum][v][Nseg[chNum]-1];
                    clust_Ch_[chNum][v1][Nseg[chNum]] = clust_Ch_[chNum][v1][Nseg[chNum]-1];
                    clust_Ch_[chNum][x1][Nseg[chNum]]=0;
                    //
                    //		      ::out1<<"find doub_x2 iseg "<<Nseg[chNum]<<" wires_Ch[chNum] "<< wires_Ch[chNum][0][Nseg[chNum]]<<" "<< wires_Ch[chNum][1][Nseg[chNum]]<<" "<<wires_Ch[chNum][2][Nseg[chNum]]<<" "<<wires_Ch[chNum][3][Nseg[chNum]]<<" "<<wires_Ch[chNum][4][Nseg[chNum]]<<" "<<wires_Ch[chNum][5][Nseg[chNum]]<<endl;
                    Nseg[chNum]++;
                    if (Nseg[chNum] > kBig_ - 2) break;
                  }

                  if (Nseg[chNum] > 15) minHits4_5=5;
                  if (Nseg[chNum] > 30) minHits4_5=6;
                } else {
                  Nhits_Ch_[chNum][Nseg[chNum]] = 0;
                  wires_Ch[chNum][x1][Nseg[chNum]] = -1;
                  wires_Ch[chNum][v1][Nseg[chNum]] = -1;
                  wires_Ch[chNum][u1][Nseg[chNum]] = -1;
                  wires_Ch[chNum][x][Nseg[chNum]] = -1;
                  wires_Ch[chNum][v][Nseg[chNum]] = -1;
                  wires_Ch[chNum][u][Nseg[chNum]] = -1;
                }//else
                if (Nseg[chNum] > kBig_ - 2) break;
              }// x v u < delta
              if (Nseg[chNum] > kBig_ - 2) break;
            }//iu
            if (Nseg[chNum] > kBig_ - 2) break;
          }//iw_Ch_[chNum][u]
          if (Nseg[chNum] > kBig_ - 2) break;
        }//iv
        if (Nseg[chNum] > kBig_ - 2) break;
      }//(iw_Ch_[chNum][v]
      if (Nseg[chNum] > kBig_ - 2) break;
    }//ix
    if (Nseg[chNum] > kBig_ - 2)return;//MP
  }//iw_Ch_[chNum][x]

}//SegmentFinder



void BmnMwpcHitFinderSRC::ProcessSegments( Int_t chNum, Double_t sigma_, Float_t dw_half_, Float_t **z_loc, Int_t Min_hits, Int_t *Nseg,   Int_t **Nhits, Int_t ***wires_Ch, Int_t ***clust_Ch_, Float_t ***XVU_Ch_, Int_t *Nbest, Int_t **ind_best, Double_t **Chi2_ndf, Double_t **Chi2_ndf_best, Double_t ***par_ab,  Int_t nPlanes, Int_t* ipl_, Float_t** XVU_, Float_t** XVU_cl_,  Double_t kChi2_Max_ ) {
  cout<<"ProcessSegments: Nseg["<<chNum<<"] "<<Nseg[chNum]<<endl;
  if (Nseg[chNum] > kBig - 2) return; // For further investigation. 16.05.2018 VL
  /*
   if (fDebug) printf("PS input: %d %8.4f %8.4f %8.4f %3d %3d %3d %3d %3d %3d %8.4f %8.4f %8.4f %3d %3d %8.4f %8.4f %8.4f\n", 
        chNum,               sigma_,             dw_half,            z_loc[chNum][0],
        Min_hits,            Nseg[chNum],        Nhits[chNum][0],    wires_Ch[chNum][0][0], 
        Nbest[chNum],        ind_best[chNum][0], Chi2_ndf[chNum][0], Chi2_ndf_best[chNum][0], 
        par_ab[chNum][0][0], nPlanes,            ipl_[0],            XVU_[chNum][0], 
        XVU_cl_[chNum][0],   kChi2_Max_);
  */
  
  Float_t delta =  0.75; //0.5;
  //if ( fRunPeriod == 7 ) delta = (chNum == 2 && chNum == 3 ) ? 1.125 : 0.75; //0.5;   

  Float_t sigma2 = sigma_ * sigma_;
  //  Float_t sigm2[nPlanes] = {sigma2, sigma2, sigma2, sigma2, sigma2, sigma2};
  //  Int_t h6[nPlanes] = {1, 1, 1, 1, 1, 1};   
  Int_t h6[6] = {1, 1, 1, 1, 1, 1};   
  Int_t Min_hits6 = Min_hits;
  for (Int_t Nhitm = nPlanes; Nhitm > Min_hits - 1; Nhitm--) {
    Int_t ifNhitm = 0;

    if (Nhitm < Min_hits6) break;

    for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {

      //  if ( Chi2_ndf[chNum][iseg] > 0 )  cout<<" chNum "<<chNum<<" -- iseg "<<iseg<<" Nhits "<<Nhits[chNum][iseg]<<" Chi2_ndf "<<Chi2_ndf[chNum][iseg]<<endl;

      if (Nhits[chNum][iseg] != Nhitm) continue;
      ifNhitm = 1;
  
      Int_t h[6] = {0, 0, 0, 0, 0, 0};
      //Int_t h[nPlanes] = {0, 0, 0, 0, 0, 0};
      Int_t max_i[7] = {0, 0, 0, 0, 0, 0,   0};
      Int_t min_i[7] = {0, 0, 0, 0, 0, 0,   0};
      //      Int_t ipl_[nPlanes] = {6, 6, 6, 6, 6, 6};
      Int_t ifirst = -1;
      Int_t ilast =  6;
      if ( Nhits[chNum][iseg] > 4) {
	for (Int_t i = 0; i < nPlanes; i++){
	  if (wires_Ch[chNum][i][iseg] > -1) { h[i] = 1;
	    if ( clust_Ch_[chNum][i][iseg] ==  3) {min_i[i]= -2; 	max_i[i]= 2; ilast--; ipl_[ilast] = i;}
	    if ( clust_Ch_[chNum][i][iseg] ==  0) {min_i[i]=  0; 	max_i[i]= 0; ifirst++; ipl_[ifirst] = i;}
	    if (i>1 && i<5 ){
	      if ( clust_Ch_[chNum][i][iseg] == -1) {min_i[i]= -2;   max_i[i]= 0;  }
	      if ( clust_Ch_[chNum][i][iseg] ==  1) {min_i[i]=  0;   max_i[i]= 2;  }
	    } else {
	      if ( clust_Ch_[chNum][i][iseg] ==  1) {min_i[i]= -2;   max_i[i]= 0;  }
	      if ( clust_Ch_[chNum][i][iseg] == -1) {min_i[i]=  0;   max_i[i]= 2;  }
	    }
	  }
	}
      
	for (Int_t i = 0; i < nPlanes; i++){
	  if (wires_Ch[chNum][i][iseg] < 0 ) continue;
	  if ( abs(clust_Ch_[chNum][i][iseg]) == 1) { ifirst++; ipl_[ifirst] = i;}
	}
      }// Nhits[chNum][iseg] > 4 
      else { for (Int_t i = 0; i < nPlanes; i++){ 
    if (wires_Ch[chNum][i][iseg] > -1)  h[i] = 1;
    ipl_[i]=i;
  }
      }
     
    //  cout<<endl;
    //  cout<<" ::ProcessSegments "<<" iseg "<<iseg<<endl;
     
    //  cout<<endl;
      /*
      if ( chNum == 1){
       cout<<" iseg "<<iseg<<" clust_Ch_[chNum] "<<clust_Ch_[chNum][0][iseg] <<" "<<clust_Ch_[chNum][1][iseg] <<" "<<clust_Ch_[chNum][2][iseg] <<" "<<clust_Ch_[chNum][3][iseg] <<" "<<clust_Ch_[chNum][4][iseg]<<" "<<clust_Ch_[chNum][5][iseg]
    	<<" wires_Ch "<<wires_Ch[chNum][0][iseg]<<" "<<wires_Ch[chNum][1][iseg]<<" "<<wires_Ch[chNum][2][iseg]<<" "<<wires_Ch[chNum][3][iseg]<<" "<<wires_Ch[chNum][4][iseg]<<" "<<wires_Ch[chNum][5][iseg]
      	  <<" XVU_Ch_ "<<XVU_Ch_[chNum][0][iseg]<<" "<<XVU_Ch_[chNum][1][iseg]<<" "<<XVU_Ch_[chNum][2][iseg]<<" "<<XVU_Ch_[chNum][3][iseg]<<" "<<XVU_Ch_[chNum][4][iseg]<<" "<<XVU_Ch_[chNum][5][iseg]
	   <<" ipl "<< ipl_[0] <<" "<< ipl_[1] <<" "<< ipl_[2] <<" "<< ipl_[3] <<" "<< ipl_[4] <<" "<< ipl_[5] <<endl;
      }
      */
      //linear fit
      //     Double_t A_[4][4]; //coef matrix
      // Double_t F[4]; //free coef 
      
      Float_t dX[nPlanes];
      //      Float_t XVU[nPlanes];
      //      Float_t XVU_cl[nPlanes];
      Double_t Chi2_min_iseg = 999;
      Double_t Chi2_Super_min = 0.1;
      Double_t par_ab_curr[4];
      Double_t par_ab_cl[4];
      
      for (Int_t i10 = min_i[ipl_[0]]; i10 <= max_i[ipl_[0]] ; i10++){
	 if ( Chi2_min_iseg < Chi2_Super_min ) break;
  //-?	if ( ipl_[0] == 6 ) continue;
  if ( h[ipl_[0]] == 1){
    XVU_[chNum][ipl_[0]] = XVU_Ch_[chNum][ipl_[0]][iseg] + dw_half_*i10 ;
    sigm2[ipl_[0]] = sigma2 ;
    if ( i10 == -1 || i10 == 1 )  sigm2[ipl_[0]] = 0.5*sigma2 ;
  }
  for (Int_t i2 = min_i[ipl_[1]]; i2 <= max_i[ipl_[1]] ; i2++){ 
     if ( Chi2_min_iseg < Chi2_Super_min ) break;
    //-?	 if ( ipl_[1] == 6 ) continue;
    //-?	 XVU_[chNum][ipl_[1]] = XVU_Ch_[chNum][ipl_[1]][iseg];
   if ( h[ipl_[1]] == 1){
      XVU_[chNum][ipl_[1]] = XVU_Ch_[chNum][ipl_[1]][iseg] + dw_half_*i2 ;
      if ( ( (ipl_[0] < 3 && ipl_[1] > 2) || (ipl_[0] > 2 && ipl_[1] < 3) ) &&  abs(ipl_[0] - ipl_[1]) == 3 ){//conjugated coord
        if ( ipl_[0] + ipl_[1] >3   ) {
    //			::out1<<" ipl_[0] "<<ipl_[0]<<" ipl_[1] "<<ipl_[1]<<endl;
    if ( fabs(XVU_[chNum][ipl_[0]] - XVU_[chNum][ipl_[1]]) > 3*dw_half_ ) continue;
    //			::out1<<"     i2 "<<i2<<" XVU_-i2 "<<XVU_[chNum][ipl_[1]]<<" XVU_-conj "<<XVU_Ch_[chNum][ipl_[0]]<<endl;
        }
      }//conjugated coord
      sigm2[ipl_[1]] = sigma2 ;
      if ( i2 == -1 || i2 == 1 )  sigm2[ipl_[1]] = 0.5*sigma2 ;
   }
    for (Int_t i3 = min_i[ipl_[2]]; i3 <= max_i[ipl_[2]] ; i3++){ 
       if ( Chi2_min_iseg < Chi2_Super_min ) break;
      //-?	    if ( ipl_[2] == 6 ) continue;
      //-?	    XVU_[chNum][ipl_[2]] = XVU_Ch_[chNum][ipl_[2]][iseg];
      if ( h[ipl_[2]] == 1){
        XVU_[chNum][ipl_[2]] = XVU_Ch_[chNum][ipl_[2]][iseg] + dw_half_*i3 ;
        Bool_t conj_bad = 0;//conjugated coord
        for (Int_t ic=0; ic<2 ; ic++){
    if ( ( (ipl_[ic] < 3 && ipl_[2] > 2) || (ipl_[ic] > 2 && ipl_[2] < 3) ) &&  abs(ipl_[ic] - ipl_[2]) == 3 ){
      if ( ipl_[ic] + ipl_[2] >3  ) {
        //		      ::out1<<" ipl_[ic] "<<ipl_[ic]<<" ipl_[2] "<<ipl_[2]<<endl;
        if ( fabs(XVU_[chNum][ipl_[ic]] - XVU_[chNum][ipl_[2]]) > 3*dw_half_ ) { conj_bad = 1; break;
          //		        ::out1<<"     i3 "<<i3<<" XVU_-i3 "<<XVU_[chNum][ipl_[2]]<<" XVU_-conj "<<XVU_Ch_[chNum][ipl_[ic]]<<endl;
        }
      }
     
    }
        }
        if (  conj_bad ) continue; //conjugated coord
        sigm2[ipl_[2]] = sigma2 ;
        if ( i3 == -1 || i3 == 1 )  sigm2[ipl_[2]] = 0.5*sigma2 ;
      }
      for (Int_t i4 = min_i[ipl_[3]]; i4 <= max_i[ipl_[3]] ; i4++){ 
	 if ( Chi2_min_iseg < Chi2_Super_min ) break;
        //-?	      if ( ipl_[3] == 6 ) continue;
        //-?	      XVU_[chNum][ipl_[3]] = XVU_Ch_[chNum][ipl_[3]][iseg];
        if ( h[ipl_[3]] == 1){
    XVU_[chNum][ipl_[3]] = XVU_Ch_[chNum][ipl_[3]][iseg] + dw_half_*i4 ;
    Bool_t conj_bad = 0;//conjugated coord
    for (Int_t ic=0; ic < 3 ; ic++){
      if ( ( (ipl_[ic] < 3 && ipl_[3] > 2) || (ipl_[ic] > 2 && ipl_[3] < 3) ) &&  abs(ipl_[ic] - ipl_[3]) == 3 ){
        if ( ipl_[ic] + ipl_[3] >3  ) {
          //		         ::out1<<" ipl_[ic] "<<ipl_[ic]<<" ipl_[3] "<<ipl_[3]<<endl;
          if ( fabs(XVU_[chNum][ipl_[ic]] - XVU_[chNum][ipl_[3]]) > 3*dw_half_ )  { conj_bad = 1; break;
      //			  ::out1<<"     i4 "<<i4<<" XVU_-i4 "<<XVU_[chNum][ipl_[3]]<<" XVU_-conj "<<XVU_Ch_[chNum][ipl_[ic]]<<endl;
          }
        }
        
      }
    }
    if (  conj_bad ) continue; //conjugated coord
    sigm2[ipl_[3]] = sigma2 ;
    if ( i4 == -1 || i4 == 1 )  sigm2[ipl_[3]] = 0.5*sigma2 ;
        }
        for (Int_t i5 = min_i[ipl_[4]]; i5 <= max_i[ipl_[4]] ; i5++){ 
	   if ( Chi2_min_iseg < Chi2_Super_min ) break;
    //-?		if ( ipl_[4] == 6 ) continue;
    //-?		XVU_[chNum][ipl_[4]] = XVU_Ch_[chNum][ipl_[4]][iseg];
    if ( h[ipl_[4]] == 1){
      XVU_[chNum][ipl_[4]] = XVU_Ch_[chNum][ipl_[4]][iseg] + dw_half_*i5 ;
      Bool_t conj_bad = 0;//conjugated coord
      for (Int_t ic=0; ic < 4 ; ic++){
        if ( ( (ipl_[ic] < 3 && ipl_[4] > 2) || (ipl_[ic] > 2 && ipl_[4] < 3) ) &&  abs(ipl_[ic] - ipl_[4]) == 3 ){
          if ( ipl_[ic] + ipl_[4] >3  ) {
      //				::out1<<" ipl_[ic] "<<ipl_[ic]<<" ipl_[4] "<<ipl_[4]<<endl;
      if ( fabs(XVU_[chNum][ipl_[ic]] - XVU_[chNum][ipl_[4]]) > 3*dw_half_ )  { conj_bad = 1; break;
        //			  	  ::out1<<"     i5 "<<i5<<" XVU_-i5 "<<XVU_[chNum][ipl_[4]]<<" XVU_-conj "<<XVU_Ch_[chNum][ipl_[ic]]<<endl;
      }
          }
         
        }
      }
      if (  conj_bad ) continue; //conjugated coord
      sigm2[ipl_[4]] = sigma2 ;
      if ( i5 == -1 || i5 == 1 )  sigm2[ipl_[4]] = 0.5*sigma2 ;
    }
    //??	Double_t Chi2_i6 = 999;
    for (Int_t i6 = min_i[ipl_[5]]; i6 <= max_i[ipl_[5]] ; i6++){
      if ( Chi2_min_iseg < Chi2_Super_min ) break;
      //-?		  if ( ipl_[5] == 6 ) continue;
      //-?		  XVU_[chNum][ipl_[5]] = XVU_Ch_[chNum][ipl_[5]][iseg];
      if ( h[ipl_[5]] == 1){
        XVU_[chNum][ipl_[5]] = XVU_Ch_[chNum][ipl_[5]][iseg] + dw_half_*i6 ;
        Bool_t conj_bad = 0;//conjugated coord
        Float_t conj =0;
        for (Int_t ic=0; ic < 5 ; ic++){
          if ( ( (ipl_[ic] < 3 && ipl_[5] > 2) || (ipl_[ic] > 2 && ipl_[5] < 3) ) &&  abs(ipl_[ic] - ipl_[5]) == 3 ){
      if ( ipl_[ic] + ipl_[5] >3  ) {
        //			  	  ::out1<<" ipl_[ic] "<<ipl_[ic]<<" ipl_[5] "<<ipl_[5]<<endl;
        if ( fabs(XVU_[chNum][ipl_[ic]] - XVU_[chNum][ipl_[5]]) > 3*dw_half_ ) { conj=(XVU_[chNum][ipl_[ic]] - XVU_[chNum][ipl_[5]]); conj_bad = 1; break;
          //			    	     ::out1<<"     i6 "<<i6<<" XVU_-i6 "<<XVU_[chNum][ipl_[5]]<<" XVU_-conj "<<XVU_Ch_[chNum][ipl_[ic]]<<endl;
        }
      }
    
          }
        }
        if (  conj_bad ) continue; //conjugated coord
        sigm2[ipl_[5]] = sigma2 ;
        if ( i6 == -1 || i6 == 1 )  sigm2[ipl_[5]] = 0.5*sigma2 ;
      }
    
      //		   ::out1<<" XVU_ "<<XVU_[chNum][0]<< " "<<XVU_[chNum][1]<<" "<<XVU_[chNum][2]<<" "<<XVU_[chNum][3]<<" "<<XVU_[chNum][4]<<" "<<XVU_[chNum][5]<<" sigm2 "<<sigm2[0]<<" "<<sigm2[1]<<" "<<sigm2[2]<<" "<<sigm2[3]<<" "<<sigm2[4]<<" "<<sigm2[5]<<" ii "<<i10<<i2<<i3<<i4<<i5<<i6<<endl;

      Float_t xx = 0;
      Int_t ii = 0;
      if ( h[0] == 1 ){ ii = ii +1; xx = xx + XVU_[chNum][0]; }
      if ( h[3] == 1 ){ ii = ii +1; xx = xx + XVU_[chNum][3]; xx = xx / ii; }
      Float_t uu = 0;
      ii = 0;
      if ( h[2] == 1 ){ ii = ii +1; uu = uu + XVU_[chNum][2]; }
      if ( h[5] == 1 ){ ii = ii +1; uu = uu + XVU_[chNum][5]; uu = uu / ii; }
      Float_t vv = 0;
      ii = 0;
      if ( h[1] == 1 ){ ii = ii +1; vv = vv + XVU_[chNum][1]; }
      if ( h[4] == 1 ){ ii = ii +1; vv = vv + XVU_[chNum][4]; vv = vv / ii; }
     
      if ( fabs (xx -uu -vv) > delta ) continue;//??
     

      for(Int_t mm=0; mm<4; mm++){
        for(Int_t jj=0; jj<4; jj++){
          matrA[mm][jj] = 0.;
          matrb[mm][jj] = 0.;
        }
      }

      // Double_t A_[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};//coef matrix
      Double_t F1[4] = {0,0,0,0};//free coef 

      if (Nhits[chNum][iseg] == nPlanes) //case 6-point segment
        FillFitMatrix(chNum, matrA, z_loc, sigm2, h6);// FillFitMatrix(A_, z_loc, sigm2, h6);

      else
        FillFitMatrix(chNum, matrA, z_loc, sigm2, h );
    
        

      //  FillFreeCoefVector(chNum, F1, XVU_Ch_, iseg,  z_loc, sigm2, h, kNPlanes); // FillFreeCoefVector(F1, XVU_, z_loc, sigm2, h, kNPlanes);
      FillFreeCoefVectorXUV(chNum, F1, XVU_,  z_loc, sigm2, h);


      Double_t A0[4][4];
      for (Int_t i1 = 0; i1 < 4; i1++){
        for (Int_t j1 = 0; j1 < 4; j1++){
           A0[i1][j1] = matrA[i1][j1];// A0[i1][j1] = A_[i1][j1];
        }
      }
      //	  Double_t b[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
      InverseMatrix(matrA, matrb);// InverseMatrix(A_, b);

      Double_t sum;
      Double_t A_int[4][4];
    //MK
      for (Int_t i1 = 0; i1 < 4; ++i1) {
        for (Int_t j1 = 0; j1 < 4; ++j1) {
          sum = 0;
          for (Int_t k1 = 0; k1 < 4; ++k1) {
	    Double_t a0 = A0[i1][k1];
	    Double_t b0 = matrb[k1][j1];
	    sum += a0*b0;
          }
          A_int[i1][j1] = sum;
        }
      }

      for (Int_t i1 = 0; i1 < 4; i1++) {
        par_ab_curr[i1] = 0;
        for (Int_t j1 = 0; j1 < 4; j1++) {
          par_ab_curr[i1] += matrb[i1][j1] * F1[j1];
        }
      } 

      dX_i = new  Float_t[nPlanes];
       for (Int_t i1 = 0; i1 < nPlanes; i1++) {
          dX_i[i1] = 0.;
       }
      //	  Float_t dX_i[nPlanes] = {0,0,0, 0,0,0};
      Double_t Chi2_curr_iseg = 0;

      //Att
      //   cout<<" iseg "<<iseg<<endl;
      for (Int_t i1 = 0; i1 < nPlanes; i1++) {
        //   dX_i[i1] = 0.;
        if (wires_Ch[chNum][i1][iseg]>-1) {
          if (i1 == 0 || i1 == 3) dX_i[i1] = XVU_[chNum][i1] - par_ab_curr[0] * z_loc[chNum][i1] - par_ab_curr[1];
          if (i1 == 2 || i1 == 5) dX_i[i1] = XVU_[chNum][i1] - 0.5 * (par_ab_curr[0] + sq3 * par_ab_curr[2]) * z_loc[chNum][i1] - 0.5 * (par_ab_curr[1] + sq3 * par_ab_curr[3]);
          if (i1 == 1 || i1 == 4) dX_i[i1] = XVU_[chNum][i1] - 0.5 * (par_ab_curr[0] - sq3 * par_ab_curr[2]) * z_loc[chNum][i1] - 0.5 * (par_ab_curr[1] - sq3 * par_ab_curr[3]);
          Chi2_curr_iseg =  Chi2_curr_iseg + dX_i[i1] * dX_i[i1] / sigm2[i1];//??
	  //  cout<<" i1 "<<i1<<" iseg "<<iseg<<" wires_Ch "<<wires_Ch[chNum][i1][iseg]<<" XVU "<<XVU_[chNum][i1]<<" dX_i "<<dX_i[i1]<<endl;
        } 
      }//i1
           
      // cout<<" iseg "<<iseg<<"  --- Chi2_curr_iseg "<<Chi2_curr_iseg<<" Chi2_min_iseg "<<Chi2_min_iseg<<" i1-6 "<<i10<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<" "<<i6<<endl;

      //??	  if ( Chi2_curr_iseg > Chi2_i6) break;              // optimization
      if ( Chi2_curr_iseg > Chi2_min_iseg * 10. ){ Nhits[chNum][iseg] = 0;  //20.04.2018
	break; }// optimization
      //??	  Chi2_i6 =  Chi2_curr_iseg;

      if ( Chi2_curr_iseg < Chi2_min_iseg) {
        Chi2_min_iseg = Chi2_curr_iseg;
        for (Int_t i1 = 0; i1 < nPlanes; i1++) {
          dX[i1] = dX_i[i1];
          XVU_cl_[chNum][i1] = XVU_[chNum][i1];
	  if (i1 > 3) continue;
          par_ab_cl[i1] =  par_ab_curr[i1];
	  //????
        }
      }

    }//i6
        }//i5
      }//i4
    }//i3
  }//i2
      }//i10

      Chi2_ndf[chNum][iseg] = Chi2_min_iseg; 

      // cout<<"1 iseg "<<iseg<<" Chi2_ndf "<<Chi2_ndf[chNum][iseg]<<endl;

      if (Nhits[chNum][iseg] > 4)  Chi2_ndf[chNum][iseg] = Chi2_ndf[chNum][iseg] / (Nhits[chNum][iseg] - 4);

      if (Chi2_ndf[chNum][iseg] > kChi2_Max_) {
	
	if (Nhits[chNum][iseg] <= Min_hits) { 
	  Nhits[chNum][iseg] = 0; 
	  continue;
	}
	else  { //reject most distant point
	  Float_t Max_dev = 0;
	  Int_t irej = -1;
	  for (Int_t i1 = 0; i1 < nPlanes; i1++){
	    if (wires_Ch[chNum][i1][iseg]>-1){
	      if ( fabs(dX[i1]) > delta ) {
		wires_Ch[chNum][i1][iseg] = -1;
		Nhits[chNum][iseg]--;// reject bad point
		continue;
	      }
	      if (fabs(dX[i1]) > Max_dev) {
		//	cout<<" i1 "<<i1<<" iseg "<<iseg<<" dX[i1] "<<dX[i1]<<" Max_dev "<<Max_dev<<" wires_Ch "<<wires_Ch[chNum][i1][iseg]<<endl;
		irej = i1;
		Max_dev = fabs(dX[i1]);
	      }
	    }//if (wires_Ch[chNum]
	  }//i1
	  //  cout<<" chNum "<<chNum<<" irej "<<irej<<" iseg "<<iseg <<" wires "<<wires_Ch[chNum][irej][iseg]<<" Nhits "<<Nhits[chNum][iseg]<<endl;
	  wires_Ch[chNum][irej][iseg] = -1;
	  Nhits[chNum][iseg]--;// reject most distant point
	  //  continue;
	  //	}
	  if (Nhits[chNum][iseg] <= Min_hits) { 
	    Nhits[chNum][iseg] = 0; 
	  }
	  continue;
	}//	else  { //reject most distant point
      }// if (Chi2_ndf[chNum][iseg]
      else {
  //	cout<<"2 Ch= "<<chNum<<" iseg "<<iseg<<" Nhits(aft.i10) "<<Nhits[chNum][iseg]<<" Chi2/n "<<Chi2_ndf[chNum][iseg]<<" XVU_Ch_[chNum] "<<XVU_Ch_[chNum][0][iseg]<<" "<<XVU_Ch_[chNum][1][iseg]<<" "<<XVU_Ch_[chNum][2][iseg]<<" "<<XVU_Ch_[chNum][3][iseg]<<" "<<XVU_Ch_[chNum][4][iseg]<<" "<<XVU_Ch_[chNum][5][iseg]<<endl;
	for (Int_t i1 = 0; i1 < nPlanes; i1++) {
	  XVU_Ch_[chNum][i1][iseg] = XVU_cl_[chNum][i1];
	  if (i1 > 3) continue;
	  par_ab[chNum][i1][iseg] =  par_ab_cl[i1];
	}
      }
   
    } // iseg

    if (!ifNhitm) continue;

    if (Nbest[chNum] == 0) {
      Double_t Chi2_best = 9999.0;
      Int_t iseg_best = -1;
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
        if (Nhits[chNum][iseg] != Nhitm) continue;
        //cout << "ch2ndf "<< Chi2_ndf[chNum][iseg] <<endl; 
        if (Chi2_ndf[chNum][iseg] >= Chi2_best) continue;
        Chi2_best = Chi2_ndf[chNum][iseg];
        iseg_best = iseg;
      } // iseg

      if (iseg_best == -1) continue;
     
      //  cout<<"3 Ch= "<< chNum <<" iseg_best1 "<< iseg_best <<" Chi2_best "<< Chi2_best<<endl;
  //  <<" wires_Ch "<< wires_Ch[chNum][0][iseg_best]<<" "<< wires_Ch[chNum][1][iseg_best]<<" "<< wires_Ch[chNum][2][iseg_best]<<" "<< wires_Ch[chNum][3][iseg_best]<<" "<< wires_Ch[chNum][4][iseg_best]<<" "<< wires_Ch[chNum][5][iseg_best]
  //   <<" XVU_Ch_[chNum] "<<XVU_Ch_[chNum][0][iseg_best]<<" "<<XVU_Ch_[chNum][1][iseg_best]<<" "<<XVU_Ch_[chNum][2][iseg_best]<<" "<<XVU_Ch_[chNum][3][iseg_best]<<" "<<XVU_Ch_[chNum][4][iseg_best]<<" "<<XVU_Ch_[chNum][5][iseg_best]<<endl;
      

      ind_best[chNum][0] = iseg_best;
      Chi2_ndf_best[chNum][0] = Chi2_best;
      Nbest[chNum] = 1;
   

      //  cout<<"4 ind_best[chNum] "<<ind_best[chNum][0]<< " Chi2_ndf_best "<<Chi2_ndf_best[chNum][0]<<" Nbest[chNum] "<<Nbest[chNum]<<endl;
      
      // if ( chNum == 1) {
  Min_hits6 = Min_hits ;//??
  //	if ( Chi2_best < 0.6) Min_hits6 = Min_hits +1 ;// //	if ( Chi2_best < 0.6) Min_hits6 = Nhitm;// 
  // }

  //	cout<<" Ch= "<<chNum<<" Nhits (bef.rej) "; for (int iseg=0; iseg< Nseg[chNum]; iseg++){  cout<<Nhits[chNum][iseg]<<" "; } cout<<endl;

      //reject(common points)
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
  if (iseg == iseg_best)continue;

  for (Int_t i1 = 0; i1 < nPlanes; i1++) {
    if (wires_Ch[chNum][i1][iseg]>-1) {
     
      if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best]) < 3*dw_half_ ) 
        Nhits[chNum][iseg] = 0; // if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best]) < 3*dw_half_ ) Nhits[chNum][iseg] = 0;
    }
  }
      }// iseg

    }// Nbest[chNum] == 0


    // cout<<" Ch= "<<chNum<<" Nhits_(One best) ";  cout<<" Nseg "<<Nseg[chNum]<<endl; for (int iseg=0; iseg< Nseg[chNum]; iseg++){  cout<<Nhits[chNum][iseg]<<" ";  cout<<Chi2_ndf_best[chNum][iseg]<<" "; } cout<<endl;
  
    if (Nbest[chNum] == 1) {//if (Nbest[chNum] == 1) {
      Double_t Chi2_best = 9999;
      Int_t iseg_best2 = -1;

      //	cout<<"  if (Nbest[chNum] == 1) ind_best[chNum] "<<ind_best[chNum][Nbest[chNum]]<< " Chi2_ndf_best "<<Chi2_ndf_best[chNum][ind_best[chNum][Nbest[chNum]]]<<" Nbest[chNum] "<<Nbest[chNum]<<endl;
  
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
  //	if (Nhits[chNum][iseg] == 0) continue;	//	if (Nhits[chNum][iseg] != Nhitm) continue;//??
  //	cout<<" Nhitm "<<Nhitm<<" Nhits[chNum] "<<Nhits[chNum][iseg]<<" Chi2_ndf "<<Chi2_ndf[chNum][iseg]<<endl;
  if (iseg == ind_best[chNum][0])continue;
  if (Nhits[chNum][iseg] != Nhitm) continue;//
  if (Chi2_ndf[chNum][iseg] > Chi2_best) continue;
  Chi2_best = Chi2_ndf[chNum][iseg];
  iseg_best2 = iseg;
      } // iseg

      if (iseg_best2>-1) {
  //	cout<<" if (iseg_best2>-1) ind_best[chNum] "<<ind_best[chNum][Nbest[chNum]]<< " Chi2_ndf_best "<<Chi2_ndf_best[chNum][ind_best[chNum][Nbest[chNum]]]<<" Nbest[chNum] "<<Nbest[chNum]<<endl;

  ind_best[chNum][Nbest[chNum]] = iseg_best2;
  Chi2_ndf_best[chNum][Nbest[chNum]] = Chi2_best;
  Nbest[chNum]++;

  //	cout<<" Ch= "<< chNum <<" iseg_best2 "<< iseg_best2 <<" Chi2_best "<< Chi2_best<<" wires_Ch "<< wires_Ch[chNum][0][iseg_best2]<<" "
  //   <<wires_Ch[chNum][1][iseg_best2]<<" "<<wires_Ch[chNum][2][iseg_best2]<<" "<<wires_Ch[chNum][3][iseg_best2]<<" "<<wires_Ch[chNum][4][iseg_best2]<<" "<<wires_Ch[chNum][5][iseg_best2]<<endl;
    
      }
    }//Nbest[chNum] == 1

    // cout<<" ind_best[chNum] "<<ind_best[chNum][Nbest[chNum]]<< " Chi2_ndf_best "<<Chi2_ndf_best[chNum][ind_best[chNum][Nbest[chNum]]]<<" Nbest[chNum] "<<Nbest[chNum]<<endl;
    
    
    if (Nbest[chNum] == 2) {

      Double_t Chi2_best = 20;
      Int_t iseg_best3 = -1;
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
  if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1])continue;
  if (Nhits[chNum][iseg] != Nhitm) continue;
  if (Chi2_ndf[chNum][iseg] > Chi2_best) continue;
  Chi2_best = Chi2_ndf[chNum][iseg];
  iseg_best3 = iseg;
      } // iseg

      if (iseg_best3 >-1) {

  ind_best[chNum][Nbest[chNum]] = iseg_best3;
  Chi2_ndf_best[chNum][Nbest[chNum]] = Chi2_best;
  Nbest[chNum]++;

  //reject(common points)
  for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
    if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1] || iseg == iseg_best3)continue;
    for (Int_t i1 = 0; i1 < nPlanes; i1++) {
      if (wires_Ch[chNum][i1][iseg]>-1) {
        if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best3]) < 3*dw_half ) Nhits[chNum][iseg] = 0;
      }
    }
  }
      }

    }
    

    if (Nbest[chNum] == 3) {

      Double_t Chi2_best = 20;
      Int_t iseg_best3 = -1;
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
  if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1] || iseg == ind_best[chNum][2])continue;
  if (Nhits[chNum][iseg] != Nhitm) continue;
  if (Chi2_ndf[chNum][iseg] > Chi2_best) continue;
  Chi2_best = Chi2_ndf[chNum][iseg];
  iseg_best3 = iseg;
      } // iseg

      if (iseg_best3 >-1) {

  ind_best[chNum][Nbest[chNum]] = iseg_best3;
  Chi2_ndf_best[chNum][Nbest[chNum]] = Chi2_best;
  Nbest[chNum]++;

  //reject(common points)
  for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
    if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1] || iseg == iseg_best3)continue;
    for (Int_t i1 = 0; i1 < nPlanes; i1++) {
      if (wires_Ch[chNum][i1][iseg]>-1) {
        if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best3]) < 3*dw_half ) Nhits[chNum][iseg] = 0;
      }
    }
  }
      }

    }
    
    if (Nbest[chNum] == 4) {

      Double_t Chi2_best = 20;
      Int_t iseg_best3 = -1;
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
	if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1] || iseg == ind_best[chNum][2] || iseg == ind_best[chNum][3])continue;
	if (Nhits[chNum][iseg] != Nhitm) continue;
	if (Chi2_ndf[chNum][iseg] > Chi2_best) continue;
	Chi2_best = Chi2_ndf[chNum][iseg];
	iseg_best3 = iseg;
      } // iseg

      if (iseg_best3 >-1) {

	ind_best[chNum][Nbest[chNum]] = iseg_best3;
	Chi2_ndf_best[chNum][Nbest[chNum]] = Chi2_best;
	Nbest[chNum]++;

	//reject(common points)
	for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
	  if (iseg == ind_best[chNum][0] || iseg == ind_best[chNum][1] || iseg == iseg_best3 )continue;
	  for (Int_t i1 = 0; i1 < nPlanes; i1++) {
	    if (wires_Ch[chNum][i1][iseg]>-1) {
	      if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best3]) < 3*dw_half ) Nhits[chNum][iseg] = 0;
	    }
	  }
	}
      }

    }

    // cout<<"Ch= "<<chNum<<" Nhits[chNum](Two-Five bests) ";  
    //	for (int iseg=0; iseg< Nseg[chNum]; iseg++){ 
    //	  cout<<Nhits[chNum][iseg]<<" "; 
    //  	} cout<<endl;
    //	cout<<"Ch= "<<chNum<<" Nbest[chNum] "<<Nbest[chNum]<<endl;
    /*
    if (fDebug){
      for (int iseg=0; iseg< Nbest[chNum]; iseg++){
        //cout<<" Ch= "<<chNum<< " ind "<<ind_best[chNum][iseg]<<" Chi2_ndf_best "<<Chi2_ndf_best[chNum][iseg]<<" par_ ab "<<par_ab[chNum][0][ind_best[chNum][iseg]]<<" "<<par_ab[chNum][1][ind_best[chNum][iseg]]<<" "<<par_ab[chNum][2][ind_best[chNum][iseg]]<<" "<<par_ab[chNum][3][ind_best[chNum][iseg]]<<endl;
        printf("PS out:   %d %8.4f %8.4f %8.4f %3d %3d %3d %3d %3d %3d %8.4f %e %8.4f %3d %3d %8.4f %8.4f %8.4f\n", 
        chNum,               sigma_,             dw_half,            z_loc[chNum][0],
        Min_hits,            Nseg[chNum],        Nhits[chNum][0],    wires_Ch[chNum][0][0], 
        Nbest[chNum],        ind_best[chNum][0], Chi2_ndf[chNum][iseg], Chi2_ndf_best[chNum][ind_best[chNum][iseg]], 
        par_ab[chNum][0][ind_best[chNum][iseg]], nPlanes, ipl_[0],   XVU_[chNum][iseg],     XVU_cl_[chNum][iseg], 
        kChi2_Max_);
		
      }
    }
    */
  } // Nhitm

}// ProcessSegments


void BmnMwpcHitFinderSRC::PrepareArraysToProcessEvent(){
  fBmnMwpcSegmentsArray->Clear();
// Clean and initialize arrays:

 for(Int_t iCh = 0; iCh < kNChambers; iCh++){
  Nseg_Ch[iCh]  = 0;
  Nbest_Ch[iCh] = 0;
  for(Int_t iPl = 0; iPl < kNPlanes; iPl++){
    iw_Ch[iCh][iPl]  = 0;
    XVU[iCh][iPl]    = 0;
    XVU_cl[iCh][iPl] = 0;
    for(Int_t iBig=0; iBig<kBig; iBig++){
      Wires_Ch[iCh][iPl][iBig] = -1;
      clust_Ch[iCh][iPl][iBig] = -1;
      XVU_Ch[iCh][iPl][iBig]   = -999.;    
    }
    for(Int_t iWire=0; iWire<kNWires; iWire++){
      wire_Ch[iCh][iWire][iPl] = 0;//  wire_Ch1[iWire][iPlane] = 0;
      xuv_Ch[iCh][iWire][iPl]  = 0.;// xuv_Ch1[iWire][iPlane] = 0;	
    }
  }

    for(Int_t ii = 0; ii < 4; ii++){
      for(Int_t jj=0; jj < kBig; jj++){
        par_ab_Ch[iCh][ii][jj] = 999.;
      }
    }

    for(Int_t iBig=0; iBig<kBig; iBig++){     
      Nhits_Ch[iCh][iBig]    = 0;
      Chi2_ndf_Ch[iCh][iBig] = 0;
    }

    for(Int_t i=0; i < kmaxSeg; i++){
      ind_best_Ch[iCh][i]      = 0;     // ind_best_Ch1[i] = 0;
      best_Ch_gl[iCh][i]       = -1;    // best_Ch1_gl[i] = -1;
      Chi2_ndf_best_Ch[iCh][i] = -999.; // Chi2_ndf_best_Ch1[i] = -999.;
    }
  }//iCh

  for(Int_t iPl=0; iPl<kNPlanes; iPl++){     
    sigm2[iPl] = sigma*sigma;  
    ipl[iPl] = 6;
    z2[iPl] = 0;
  }

  for(Int_t ii  = 0; ii < 4; ii++){
    for(Int_t jj  = 0; jj < 4; jj++){
      matrA[ii][jj] = 0.;
      matrb[ii][jj] = 0.;
    }
  }
}//PrepareArraysToProcessEvent


void BmnMwpcHitFinderSRC::FillFitMatrix(Int_t chN, Double_t** AA, Float_t** z, Float_t* sigm2_, Int_t* h_) {

  //out1<<" in FillFitMatrix "<<endl;

  // AA - matrix to be filledlayers)
  // sigm2 - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)
  // Float_t z2_[nPlanes];
  Float_t z2_[6] = {z[chN][0] * z[chN][0], z[chN][1] * z[chN][1], z[chN][2] * z[chN][2], z[chN][3] * z[chN][3], z[chN][4] * z[chN][4], z[chN][5] * z[chN][5]}; //cm

  AA[0][0] += 2 * z2_[0] * h_[0] / sigm2_[0] 
           +      z2_[2] * h_[2] / (2 * sigm2_[2]) 
           +      z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z2_[3] * h_[3] / sigm2_[3] 
           +      z2_[5] * h_[5] / (2 * sigm2_[5]) 
           +      z2_[4] * h_[4] / (2 * sigm2_[4]); //Ax

  AA[0][1] += 2 * z[chN][0] * h_[0] / sigm2_[0] 
           +      z[chN][2] * h_[2] / (2 * sigm2_[2]) 
           +      z[chN][1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z[chN][3] * h_[3] / sigm2_[3] 
           +      z[chN][5] * h_[5] / (2 * sigm2_[5]) 
           +      z[chN][4] * h_[4] / (2 * sigm2_[4]); //Bx

  AA[0][2] += sq3 * (z2_[2] * h_[2] / (2 * sigm2_[2]) 
           -         z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +         z2_[5] * h_[5] / (2 * sigm2_[5]) 
           -         z2_[4] * h_[4] / (2 * sigm2_[4])); //Ay

  AA[0][3] += sq3 * (z[chN][2] * h_[2] / (2 * sigm2_[2]) 
           -         z[chN][1] * h_[1] / (2 * sigm2_[1]) 
           +         z[chN][5] * h_[5] / (2 * sigm2_[5]) 
           -         z[chN][4] * h_[4] / (2 * sigm2_[4])); //By

  AA[1][0] = AA[0][1];

  AA[1][1] +=   2 * h_[0] / sigm2_[0] 
           +  0.5 * h_[2] / sigm2_[2] + 0.5 * h_[1] / sigm2_[1] 
           +    2 * h_[3] / sigm2_[3] + 0.5 * h_[5] / sigm2_[5] 
           +  0.5 * h_[4] / sigm2_[4];

  AA[1][2] += sq3 * (z[chN][2] * h_[2] / sigm2_[2] 
           - z[chN][1] * h_[1] / sigm2_[1] 
           + z[chN][5] * h_[5] / sigm2_[5] 
           - z[chN][4] * h_[4] / sigm2_[4]) * 0.5;

  AA[1][3] += sq3 * (h_[2] / sigm2_[2] 
           -         h_[1] / sigm2_[1] 
           +         h_[5] / sigm2_[5] 
           -         h_[4] / sigm2_[4]) * 0.5;

  AA[2][0] = AA[0][2];

  AA[2][1] = AA[1][2];

  AA[2][2] += 3.0 * (z2_[2] * h_[2] / sigm2_[2] 
           +         z2_[1] * h_[1] / sigm2_[1] 
           +         z2_[5] * h_[5] / sigm2_[5] 
           +         z2_[4] * h_[4] / sigm2_[4]) * 0.5;

  AA[2][3] += 3.0 * (z[chN][2] * h_[2] / sigm2_[2] 
           +         z[chN][1] * h_[1] / sigm2_[1] 
           +         z[chN][5] * h_[5] / sigm2_[5] 
           +         z[chN][4] * h_[4] / sigm2_[4])   * 0.5;

  AA[3][0] = AA[0][3];
  AA[3][1] = AA[1][3];
  AA[3][2] = AA[2][3];
  AA[3][3] += 3.0 * (0.5 * h_[2] / sigm2_[2] 
           +  0.5 *        h_[1] / sigm2_[1] 
           +  0.5 *        h_[5] / sigm2_[5] 
           +  0.5 *        h_[4] / sigm2_[4]);
}



void BmnMwpcHitFinderSRC::FillFreeCoefVectorXUV(Int_t ichNum ,  Double_t* F, Float_t** XVU_, Float_t** z, Float_t* sigm2_, Int_t* h_) {
  // F - vector to be filled
  // XVU_ - coordinates of segment in chamber (Is it correct definition?)
  // segIdx - index of current segment
  // z - local z-positions of planes(layers)
  // sigm2_ - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)
  F[0] += 
    2 * XVU_[ichNum][0]  * z[ichNum][0] * h_[0] / sigm2_[0] 
    + 
    XVU_[ichNum][1]   * z[ichNum][1] * h_[1] / sigm2_[1] 
    + 
    XVU_[ichNum][2]  * z[ichNum][2] * h_[2] / sigm2_[2] 
    + 
    2 * XVU_[ichNum][3]  * z[ichNum][3] * h_[3] / sigm2_[3] 
    + 
    XVU_[ichNum][4]  * z[ichNum][4] * h_[4] / sigm2_[4] 
    + 
    XVU_[ichNum][5]  * z[ichNum][5] * h_[5] / sigm2_[5];

  F[1] += 2 * XVU_[ichNum][0]  * h_[0] / sigm2_[0] + XVU_[ichNum][1]  * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * h_[2] / sigm2_[2] + 2 * XVU_[ichNum][3]  * h_[3] / sigm2_[3] + XVU_[ichNum][4]  * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * h_[5] / sigm2_[5];
  F[2] += (-XVU_[ichNum][1]  * z[ichNum][1] * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * z[ichNum][2] * h_[2] / sigm2_[2] - XVU_[ichNum][4]  * z[ichNum][4] * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * z[ichNum][5] * h_[5] / sigm2_[5]);
  F[3] +=  (-XVU_[ichNum][1]  * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * h_[2] / sigm2_[2] - XVU_[ichNum][4]  * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * h_[5] / sigm2_[5]);

  F[2]=F[2]*sq3;
  F[3]=F[3]*sq3;
}


void BmnMwpcHitFinderSRC::InverseMatrix(Double_t** AA, Double_t** bb) {
  // Gaussian algorithm for 4x4 matrix inversion 
  Double_t factor;
  Double_t temp[4];
  // Set b to I
  for (Int_t i1 = 0; i1 < 4; i1++){
    for (Int_t j1 = 0; j1 < 4; j1++){
      if (i1 == j1) bb[i1][j1] = 1.0;
      else bb[i1][j1] = 0.0;
    }
  }
  for (Int_t i1 = 0; i1 < 4; i1++) {
    for (Int_t j1 = i1 + 1; j1 < 4; j1++) {
      if (fabs(AA[i1][i1]) < fabs(AA[j1][i1])) {
        for (Int_t l1 = 0; l1 < 4; l1++) temp[l1] = AA[i1][l1];
        for (Int_t l1 = 0; l1 < 4; l1++) AA[i1][l1] = AA[j1][l1];
        for (Int_t l1 = 0; l1 < 4; l1++) AA[j1][l1] = temp[l1];
        for (Int_t l1 = 0; l1 < 4; l1++) temp[l1] = bb[i1][l1];
        for (Int_t l1 = 0; l1 < 4; l1++) bb[i1][l1] = bb[j1][l1];
        for (Int_t l1 = 0; l1 < 4; l1++) bb[j1][l1] = temp[l1];
      }
    }
    factor = AA[i1][i1];
    for (Int_t j1 = 4 - 1; j1>-1; j1--) {
      bb[i1][j1] /= factor;
      AA[i1][j1] /= factor;
    }
    for (Int_t j1 = i1 + 1; j1 < 4; j1++) {
      factor = -AA[j1][i1];
      for (Int_t k1 = 0; k1 < 4; k1++) {
          AA[j1][k1] += AA[i1][k1] * factor;
          bb[j1][k1] += bb[i1][k1] * factor;
      }
    }
  } // i1
  for (Int_t i1 = 3; i1 > 0; i1--) {
    for (Int_t j1 = i1 - 1; j1>-1; j1--) {
      factor = -AA[j1][i1];
      for (Int_t k1 = 0; k1 < 4; k1++) {
          AA[j1][k1] += AA[i1][k1] * factor;
          bb[j1][k1] += bb[i1][k1] * factor;
      }
    }
  } // i1
}    //end inverse

void BmnMwpcHitFinderSRC::Finish() {
  printf("MWPC hit finder: write hists to file... ");
  fOutputFileName = TString("hMWPChits.root");
  TFile file(fOutputFileName, "RECREATE");
  if(fDoTest) fList.Write();
  file.Close();
  printf("done\n");
  delete fMwpcGeometrySRC;	

  cout << "Work time of the MWPC hit finder: " << workTime << " s" << endl;
}

ClassImp(BmnMwpcHitFinderSRC)


