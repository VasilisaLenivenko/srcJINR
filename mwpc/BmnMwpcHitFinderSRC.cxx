// @(#)bmnroot/mwpc:$Id$
// Author: Pavel Batyuk <pavel.batyuk@jinr.ru> and V.Lenivenko 2017-02-10

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcHitFinderSRC                                                        //
//                                                                            //
// Implementation of an algorithm developed by                                //
// Vasilisa Lenivenko  and Vladimir Palichik                                  //
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

BmnMwpcHitFinderSRC::BmnMwpcHitFinderSRC(Bool_t isExp) :
  fEventNo(0),
  fUseDigitsInTimeBin(kTRUE),
  expData(isExp) {
  fInputBranchName = "MWPC";
  fOutputBranchName = "BmnMwpcHit";
  fOutputBranchName1 = "BmnMwpcSegment";
  fOutputBranchName2 = "BmnMwpcTrack";
  thDist = 1.;
  nInputDigits = 3;
  nTimeSamples = 3;
  kBig = 100;
  }

BmnMwpcHitFinderSRC::~BmnMwpcHitFinderSRC() {

}

InitStatus BmnMwpcHitFinderSRC::Init() {
  if (!expData) return kERROR;
  if (fVerbose) cout << " BmnMwpcHitFinderSRC::Init() " << endl;

  FairRootManager* ioman = FairRootManager::Instance();
  fBmnMwpcDigitArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
  if (!fBmnMwpcDigitArray)
  {
    cout<<"BmnMwpcHitFinderSRC::Init(): branch "<<fInputBranchName<<" not found! Task will be deactivated"<<endl;
    SetActive(kFALSE);
    return kERROR;
  }

  fBmnMwpcHitArray = new TClonesArray(fOutputBranchName);
  ioman->Register(fOutputBranchName.Data(), "MWPC", fBmnMwpcHitArray, kTRUE);

  fBmnMwpcSegmentsArray = new TClonesArray(fOutputBranchName1);
  ioman->Register(fOutputBranchName1.Data(), "MWPC", fBmnMwpcSegmentsArray, kTRUE);

  fBmnMwpcTracksArray = new TClonesArray(fOutputBranchName2);
  ioman->Register(fOutputBranchName2.Data(), "MWPC", fBmnMwpcTracksArray, kTRUE);

  fMwpcGeometry = new BmnMwpcGeometry();//fMwpcGeometry = new BmnMwpcGeometry(6);
  kNChambers = 4;//fMwpcGeometry->GetNChambers(); 
  kNPlanes = fMwpcGeometry->GetNPlanes(); // 6
  kNWires = fMwpcGeometry->GetNWires();

  cout<<" kNChambers "<<kNChambers<<" kNPlanes "<<kNPlanes<<" kNWires "<<kNWires<<endl;

  TVector3 Ch1Cent = fMwpcGeometry->GetChamberCenter(0);
  TVector3 Ch2Cent = fMwpcGeometry->GetChamberCenter(1);
  TVector3 Ch3Cent = fMwpcGeometry->GetChamberCenter(2);
  TVector3 Ch4Cent = fMwpcGeometry->GetChamberCenter(3);

  ZCh =  new Float_t[kNChambers];

  ZCh[0] = Ch1Cent.Z(); 
  ZCh[1]=  Ch2Cent.Z();
  ZCh[2] = Ch1Cent.Z(); 
  ZCh[3] = Ch2Cent.Z();

  for(Int_t i = 0; i < kNChambers; i++){ 
    cout<<" ZCh["<<i<<"]= "<< ZCh[i]<<endl;
  }

  kZmid =  new Float_t[kNChambers];

  kZmid[0]= -75.75; // ( Ch1Cent.Z()-Ch2Cent.Z() )*0.5;////cm
  kZmid[1]=  75.75; // -( Ch1Cent.Z()-Ch2Cent.Z() )*0.5;//
  kZmid[2]= -50.; //  ( Ch3Cent.Z()-Ch4Cent.Z() )*0.5;
  kZmid[3]=  50.; //  -( Ch3Cent.Z()-Ch4Cent.Z() )*0.5;

  for(Int_t i = 0; i < kNChambers; i++){ 
    cout<<" kZmid["<<i<<"]= "<< kZmid[i]<<endl;
  }

  cout<<"  dZ(ch1-ch2) = "<< -( Ch1Cent.Z()-Ch2Cent.Z() )<<endl;

  // kZ_to_pole = - (Ch2Cent.Z() + kZmid1);
  // cout<<" kZmid1 = "<<kZmid1<<" kZmid2 = "<<kZmid2 << " kZ_to_pole "<< - (Ch2Cent.Z() + kZmid1)<<endl;   

  hNp_best_ch1 =  new TH1D("hNp_best_ch1", " Np in ch1; point; ", 6, 1., 7.); fList.Add(hNp_best_ch1);
  hNp_best_ch2 =  new TH1D("hNp_best_ch2", " Np in ch2; point; ", 6, 1., 7.); fList.Add(hNp_best_ch2);
  hNbest_Ch1 =    new TH1D("hNbest_Ch1", " Nbest_Ch1 ", 6, 0.,6.);  fList.Add(hNbest_Ch1);
  hNbest_Ch2 =    new TH1D("hNbest_Ch2", " Nbest_Ch2 ", 6, 0.,6.);  fList.Add(hNbest_Ch2);
  hChi2_ch1_2 =  new TH1D("hChi2_ch1_2", " Chi2_ch1_2 ", 500, 0., 500.);fList.Add(hChi2_ch1_2);

  kMinHits = 4;
  kChi2_Max = 20.;

  dw = fMwpcGeometry->GetWireStep();//0.25; // [cm] // wires step
  dw_half = 0.5*dw;
  sq3 = sqrt(3.);
  sq12 = sqrt(12.);
  sigma = dw/sq12;
  kMiddlePl = 47.25;

  shift = new Float_t*[kNChambers];
  kZ_loc = new Float_t*[kNChambers];// kZ1_loc = new Float_t[kNPlanes];
  z_gl = new Float_t*[kNChambers];// z_gl1 = new Float_t[kNPlanes];

  sigm2 = new Float_t[kNPlanes]; 
  ipl = new Int_t[kNPlanes];
  
  

 
  shift1_2 = new Float_t[4];

 
  z2 = new Float_t[kNPlanes];

  cout<<" 1 "<<endl;
  
  
  kPln = new Int_t*[kNChambers];
  iw = new Int_t*[kNChambers];
  iw_Ch = new Int_t*[kNChambers];
  wire_Ch = new Int_t**[kNChambers];
  xuv_Ch = new Float_t**[kNChambers];
  Wires_Ch = new Int_t**[kNChambers]; //Wires_Ch1 = new Int_t*[kNPlanes];
  clust_Ch = new Int_t**[kNChambers];// clust_Ch1 = new Int_t*[kNPlanes];
  XVU_Ch = new Float_t**[kNChambers];// XVU_Ch1 = new Float_t*[kNPlanes];
  Nhits_Ch = new Int_t*[kNChambers];// Nhits_Ch1 = new Int_t[kBig];
  Nseg_Ch = new Int_t[kNChambers];// Nseg_Ch1
  Nbest_Ch = new Int_t[kNChambers];//Nbest_Ch1 = 0;
  ind_best_Ch = new Int_t*[kNChambers];
  best_Ch_gl  = new Int_t*[kNChambers];
  Chi2_ndf_Ch = new Double_t*[kNChambers];
  Chi2_ndf_best_Ch = new Double_t*[kNChambers];
  par_ab_Ch = new Double_t**[kNChambers]; // par_ab_Ch1 = new Double_t*[4];
  XVU = new Float_t*[kNChambers];// XVU1 = new Float_t[kNPlanes];
  XVU_cl = new Float_t*[kNChambers];// XVU_cl1 = new Float_t[kNPlanes];
 

  cout<<" 2 "<<endl;

  for(Int_t ii = 0; ii < kNChambers; ii++){ 

    kPln[ii] = new Int_t[kNPlanes];
    iw[ii] = new Int_t[kNPlanes];
    iw_Ch[ii] = new Int_t[kNPlanes]; 
    shift[ii] = new Float_t[4];// shift1 = new Float_t[4];
    kZ_loc[ii] = new Float_t[kNPlanes];
    z_gl[ii] = new Float_t[kNPlanes];
    Nhits_Ch[ii] = new Int_t[kBig];

    wire_Ch[ii]= new Int_t*[kNWires];
    xuv_Ch[ii] = new Float_t*[kNWires];
    Wires_Ch[ii] = new Int_t*[kNPlanes];
    clust_Ch[ii] = new Int_t*[kNPlanes];
    XVU_Ch[ii] = new Float_t*[kNPlanes];
    par_ab_Ch[ii] = new Double_t*[4];
    XVU[ii] = new Float_t[kNPlanes];
    XVU_cl[ii] = new Float_t[kNPlanes];
    ind_best_Ch[ii] = new Int_t[5];//ind_best_Ch1 = new Int_t[5];
    best_Ch_gl[ii]  = new Int_t[5];
    Chi2_ndf_Ch[ii] = new Double_t[kBig];
    Chi2_ndf_best_Ch[ii] = new Double_t[5];

  }
  cout<<" 3 "<<endl;
 
  for(Int_t ii = 0; ii < kNChambers; ii++){ 
    for(Int_t iWire=0; iWire < kNWires; iWire++){
      wire_Ch[ii][iWire] = new Int_t[kNPlanes];
      xuv_Ch[ii][iWire] = new Float_t[kNPlanes];    
    }
  }

  for(Int_t ii = 0; ii < kNChambers; ii++){ 
    for(Int_t iPla=0; iPla < kNPlanes; iPla++){
      Wires_Ch[ii][iPla] = new Int_t[kBig];
      clust_Ch[ii][iPla] = new Int_t[kBig];
      XVU_Ch[ii][iPla] = new Float_t[kBig];
    }
    for(Int_t i=0; i<4; i++){
      par_ab_Ch[ii][i] = new Double_t[kBig];
    }
  }


  cout<<" 4 "<<endl;
  
    
  // par_ab_Ch1 = new Double_t*[4];
  // par_ab_Ch2 = new Double_t*[4];
  par_ab_Ch1_2 = new Double_t*[4];
  matrA = new Double_t*[4];
  matrb = new Double_t*[4];
  
  
  for(Int_t ii=0; ii<4; ii++){
    // par_ab_Ch1[ii] = new Double_t[100];
    //  par_ab_Ch2[ii] = new Double_t[100];
    par_ab_Ch1_2[ii] = new Double_t[5];
    matrA[ii] = new Double_t[4];
    matrb[ii] = new Double_t[4];         
  }

  
  // ind_best_Ch1 = new Int_t[5];
  //  ind_best_Ch2 = new Int_t[5];
  //   best_Ch1_gl  = new Int_t[5];
  //   best_Ch2_gl  = new Int_t[5];
    ind_best_Ch1_2 = new Int_t[5];
    //   Chi2_ndf_Ch1 = new Double_t[kBig];
    //  Chi2_ndf_Ch2 = new Double_t[kBig];
    // Chi2_ndf_best_Ch1 = new Double_t[5];
    //  Chi2_ndf_best_Ch2 = new Double_t[5];
    Chi2_match = new Double_t[5];
    Chi2_ndf_Ch1_2 = new Double_t[5];
    
    //  Ch :                     1                                    2
    //                      v+  u-  x-  v-  u+  x+         v+  u-  x-   v-  u+  x+
    //    kPln[12] =      {  4,  5,  0,  1,  2,  3,         7, 11,  6,  10,  8,  9 };  //run6-II   r.1397-last
    cout<<" 5 "<<endl;
       
    kPln[0][0] = 4;
    kPln[0][1] = 5;
    kPln[0][2] = 0;
    kPln[0][3] = 1;
    kPln[0][4] = 2;
    kPln[0][5] = 3;//{4,5,0,1,2,3,  7,11,6,10,9,8,  0,0,0,0,0,0,  0,0,0,0,0,0};

    kPln[1][0] = 4;//1
    kPln[1][1] = 5;
    kPln[1][2] = 0;
    kPln[1][3] = 1;//4
    kPln[1][4] = 2;
    kPln[1][5] = 3;
    
    kPln[2][0] = 5;
    kPln[2][1] = 0;
    kPln[2][2] = 1;
    kPln[2][3] = 2;
    kPln[2][4] = 3;
    kPln[2][5] = 4;
    
    kPln[3][0] = 1;
    kPln[3][1] = 0;
    kPln[3][2] = 5;
    kPln[3][3] = 4;
    kPln[3][4] = 3;
    kPln[3][5] = 2;
  
    cout<<" 6 "<<endl;
    //                   x-    v-    u+    x+    v+    u-   // canonical order
    //    kZ1_loc[6] = {-0.5,  0.5,  1.5,  2.5, -2.5, -1.5}; //cm   run5  
    //    kZ2_loc[6] = {-0.5,  0.5,  1.5,  2.5, -2.5, -1.5}; //cm   run5, run6

    //for mpd detectors
    
    //        3         -1.5, -0.5,  0.5,  1.5,  2.5,  -2.5
    //        4         -1.5, -2.5,  2.5,  1.5,  0.5,  -0.5   
    
    for(Int_t ichh = 0; ichh < kNChambers; ichh++){
      for(int ii = 0; ii < 6; ii++){

	if ( ichh == 0 || ichh == 1){
	  kZ_loc[ichh][ii] = -0.5 + ii;
	  if(ii == 4) { kZ_loc[ichh][ii] = -2.5;}
	  if(ii == 5) { kZ_loc[ichh][ii] = -1.5;}
	}//if ( ich == 0 || ich == 1) {

	kZ_loc[2][0] = -1.5;
	kZ_loc[2][1] = -0.5;
	kZ_loc[2][2] =  0.5;
	kZ_loc[2][3] =  1.5;
	kZ_loc[2][4] =  2.5;
	kZ_loc[2][5] = -2.5;

	kZ_loc[3][0] = -1.5;
	kZ_loc[3][1] = -2.5;
	kZ_loc[3][2] =  2.5;
	kZ_loc[3][3] =  1.5;
	kZ_loc[3][4] =  0.5;
	kZ_loc[3][5] = -0.5;

       	z_gl[ichh][ii] =  kZmid[ichh] + kZ_loc[ichh][ii];
 
	//	cout<<" ich "<<ichh<<" ii "<<ii<<" kZ_loc "<<kZ_loc[ichh][ii]<<" z_gl "<<z_gl[ichh][ii]<<endl;	
      }
    }//ich   
    cout<< endl;

    //slope and shift  for parametrs    
    shift[0][0]=   0.;//x1_slope_sh = 0;
    shift[0][2]= 0.01;//y1_slope_sh = 0.01;
    shift[0][1]= -.40;//x1_sh = -.40;
    shift[0][3]= 7.83;//y1_sh = 7.83;
    
    shift[1][0]=    0.;//x2_slope_sh = 0;
    shift[1][2]= -.008;//y2_slope_sh = -.008;   
    shift[1][1]=   .24;//x2_sh = .24;/
    shift[1][3]=  6.67;//y2_sh = 6.67;

    shift[2][0]= 0.;//x1_slope_sh = 0;
    shift[2][2]= 0.;//y1_slope_sh = 0.01;
    shift[2][1]= 0.;//x1_sh = -.40;
    shift[2][3]= 0.;//y1_sh = 7.83;
    
    shift[3][0]= 0.;//x2_slope_sh = 0;
    shift[3][2]= 0.;//y2_slope_sh = -.008;   
    shift[3][1]= 0.;//x2_sh = .24;/
    shift[3][3]= 0.;//y2_sh = 6.67;
  

    shift1_2[0]= (shift[2][1] - shift[1][1])/(-( Ch1Cent.Z()-Ch2Cent.Z() ) );
    shift1_2[2]= (shift[2][3] - shift[1][3])/(-( Ch1Cent.Z()-Ch2Cent.Z() ) );
    shift1_2[1]= 0.5*(shift[2][1] - shift[1][1]);
    shift1_2[3]= 0.5*(shift[2][3] - shift[1][3]);

    cout<<" slope and shift  for parametrs   "<<endl;
    for(int ii=0; ii<4; ii++){
      cout<<ii<<" shift Ch0= "<< shift[0][ii]<<" Ch1= "<<shift[1][ii]<<" Ch2= "<<shift[2][ii]<<" Ch3= "<<shift[3][ii]<<" shift1_2 "<<shift1_2[ii]<<endl;
    }
    cout<<endl;

    return kSUCCESS;
}

void BmnMwpcHitFinderSRC::PrepareArraysToProcessEvent(){

  fBmnMwpcHitArray->Clear();
  fBmnMwpcSegmentsArray->Clear();
  fBmnMwpcTracksArray->Clear();
   
// Clean and initialize arrays:

  for(Int_t icha=0; icha < kNChambers; icha++){
    for(Int_t iPl=0; iPl < kNPlanes; iPl++){
      iw[icha][iPl] = 0;//  iw_Ch1[iPl] = 0;
      iw_Ch[icha][iPl] = 0;
      XVU[icha][iPl] = 0;//XVU1[iPl] = 0;
      XVU_cl[icha][iPl] = 0;// XVU_cl1[iPl] = 0;
   
    }
  }

  for(Int_t iCh = 0; iCh < kNChambers;  iCh++){       
    for(Int_t iWire=0; iWire<kNWires; iWire++){
      for(Int_t iPlane=0; iPlane<kNPlanes; iPlane++){
	wire_Ch[iCh][iWire][iPlane] = 0;//  wire_Ch1[iWire][iPlane] = 0;
	xuv_Ch[iCh][iWire][iPlane] = 0.;// xuv_Ch1[iWire][iPlane] = 0;	
      }
    }
  }//iCh

  for(Int_t iCh = 0; iCh < kNChambers;  iCh++){ 
    Nseg_Ch[iCh] = 0;
    Nbest_Ch[iCh] = 0;

    for(Int_t i=0; i< 5; i++){
      ind_best_Ch[iCh][i] = 0;//	ind_best_Ch1[i] = 0;
      best_Ch_gl[iCh][i] = -1;	//	best_Ch1_gl[i] = -1;
      Chi2_ndf_best_Ch[iCh][i] = -999.;	//	Chi2_ndf_best_Ch1[i] = -999.;
    }

    for(Int_t ii=0; ii<4; ii++){
      for(Int_t jj=0; jj<100; jj++){
	par_ab_Ch[iCh][ii][jj] = 999.;
      }
    }

    for(Int_t iBig=0; iBig<kBig; iBig++){     
      Nhits_Ch[iCh][iBig] = 0;
      Chi2_ndf_Ch[iCh][iBig] = 0;
    }
  }

  for(Int_t iPl=0; iPl<kNPlanes; iPl++){     
    sigm2[iPl] = sigma*sigma;  
    ipl[iPl] = 6;
   
   
    z2[iPl] = 0;
  }

  
      Nbest_Ch12_gl = 0;

      for(Int_t i=0; i< 5; i++){
	//	ind_best_Ch1[i] = 0;
	//	ind_best_Ch2[i] = 0;
	//	best_Ch1_gl[i] = -1;
	//	best_Ch2_gl[i] = -1;
	//	Chi2_ndf_best_Ch1[i] = -999.;
	//	Chi2_ndf_best_Ch2[i] = -999.;
	Chi2_match[i] = 999.;
	Chi2_ndf_Ch1_2[i] = 999.;
	ind_best_Ch1_2[i]= -1;
	
      }

      for(Int_t ii=0; ii<4; ii++){
	for(Int_t jj=0; jj<100; jj++){
	  //  par_ab_Ch1[ii][jj] = 999.;
	  //	  par_ab_Ch2[ii][jj] = 999.;
	}
      }

      for(Int_t ii=0; ii<4; ii++){
	for(Int_t jj=0; jj<5; jj++){
	  par_ab_Ch1_2[ii][jj] = 999.;	 
	}
      }
    
    for(Int_t ii=0; ii<4; ii++){
      for(Int_t jj=0; jj<4; jj++){
        matrA[ii][jj] = 0.;
        matrb[ii][jj] = 0.;
      }
    }
    

    for(Int_t iCh = 0; iCh < kNChambers;  iCh++){    
      for(Int_t iPlane=0; iPlane<kNPlanes; iPlane++){
	for(Int_t iBig=0; iBig<kBig; iBig++){
	  Wires_Ch[iCh][iPlane][iBig] = -1;// Wires_Ch1[iPlane][iBig] = -1;
	  clust_Ch[iCh][iPlane][iBig] = -1;// clust_Ch1[iPlane][iBig] = -1;
	  XVU_Ch[iCh][iPlane][iBig] = -999.;//XVU_Ch1[iPlane][iBig] = -999.;    
	}
      }
    }
      
      for(Int_t iBig=0; iBig < kBig; iBig++){
	//	Chi2_ndf_Ch1[iBig] = 0;
	//	Chi2_ndf_Ch2[iBig] = 0;
      }
}//PrepareArraysToProcessEvent

void BmnMwpcHitFinderSRC::Exec(Option_t* opt) {
    if (!IsActive())
        return;
    clock_t tStart = clock();
    PrepareArraysToProcessEvent();
    if (fVerbose) cout << "\n======================== MWPC hit finder exec started =====================\n" << endl;
    //if (fVerbose) 
    cout << "Event number: " << fEventNo++ << endl; 
    //  cout<<"NWires = "<<kNWires<<", NPlanes = "<<kNPlanes<<endl;
    
    Short_t st, wn, pn, ts, pl;
    
    for (Int_t iDigit = 0; iDigit < fBmnMwpcDigitArray->GetEntriesFast(); iDigit++) {

      //  BmnMwpcDigit* digit = (BmnMwpcDigit*) fBmnMwpcDigitArray->UncheckedAt(iDigit);
      //digit->SetUsing(kFALSE);
      BmnMwpcDigit* digit = (BmnMwpcDigit*) fBmnMwpcDigitArray->At(iDigit);
      st = digit->GetStation();
      wn = digit->GetWireNumber();
      pl = digit->GetPlane();
      ts = digit->GetTime();
      //  cout<<"++++++++"<<endl;
      // cout<<" st "<<st<<" wn = "<<wn<<", pl = "<<pl<<", ts = "<<ts<<endl;
	// digits[digit->GetPlane() / kNPlanes][digit->GetPlane() % kNPlanes].push_back(digit);
         	  
        pn =  kPln[st][pl];//
	// 	cout<<" st "<< st <<" pl = "<<pl<<", pn = "<<pn<<", wn = "<<wn<<endl;

	// Bool_t repeat = 0; 
	// if (iw[st][pn] > 0) {
	//   for (Int_t ix = 0; ix < iw[st][pn]; ix++) {
	//     if (wn == wire_Ch[st][ ix ][pn]  ) {
	//      repeat = 1;
	//      break;
	//     }
	//    }//ix
	//    if (repeat) continue;
       
	 wire_Ch[ st ][ iw[st][pn] ][ pn ] = wn;
         xuv_Ch[st][ iw[st][pn]  ][pn] = (wn - kMiddlePl) * dw;

	 //	 cout<<"iw["<<st<<"]["<<pn<<"]= "<<iw[st][pn]<<" wire_Ch "<<wire_Ch[ st ][ iw[st][pn] ][ pn ]<<" xuv_Ch "<<xuv_Ch[st][ iw[st][pn]  ][pn]<<endl;
	 if (pn == 0 || pn == 1 || pn == 5 ) xuv_Ch[st][iw[st][pn]][pn]= -xuv_Ch[st][iw[st][pn]][pn]; 
  
	 iw[st][pn]++;

	 iw_Ch[st][pn] = iw[st][pn];

	 //	 cout<<" st "<<st<<" pn "<<pn<<" iw[st][pn] "<<iw[st][pn]<<" iw_Ch[st][pn] "<<iw_Ch[st][pn]<<endl;
	 // cout<<" ---------------"<<endl;
	 
    }// iDigit


    cout<<" after digit "<<endl;       
   
    
    for (Int_t iChamber = 0; iChamber < kNChambers; iChamber++) {
      for(Int_t iCase= 1; iCase < 9; iCase ++){
	SegmentFinder(iChamber, Wires_Ch, clust_Ch, XVU_Ch, Nhits_Ch, iw_Ch, Nseg_Ch, wire_Ch, xuv_Ch, kMinHits, iCase, kBig);	
      }
      
      //   for (Int_t ise = 0; ise < Nseg_Ch[iChamber]; ise++) {
      //	hNp_best_ch1->Fill(Nhits_Ch[iChamber][ise]);
      //  }

    }
    cout<<"SegmentFinder: Nseg_Ch Ch0= "<<Nseg_Ch[0]<<" Ch1= "<<Nseg_Ch[1]<<" Ch2="<<Nseg_Ch[2]<<" Ch3="<<Nseg_Ch[3]<<endl;
    /*
   
	  
    if(Nseg_Ch2 > 0) {
    ProcessSegments(2,  sigma,   dw_half,  
    kZ_loc,  kMinHits,   Nseg_Ch, Nhits_Ch,   
    Wires_Ch,  clust_Ch,   XVU_Ch,  Nbest_Ch,   
    ind_best_Ch, Chi2_ndf_Ch,  Chi2_ndf_best_Ch, 
    par_ab_Ch, kNPlanes, ipl,  XVU, XVU_cl,  kChi2_Max);
     
      for (Int_t iBest = 0; iBest < Nbest_Ch2; iBest++) {	  
	if (Nhits_Ch2[ind_best_Ch2[iBest]] > 3) { 	 
	//  cout<<" iBest "<< iBest<<" Chi2_ndf_best_Ch2 "<<Chi2_ndf_best_Ch2[iBest]<<endl;
	  for (Int_t i = 0; i < 6; i++){
	//    cout<<" Ch= "<<2<<" XVU_Ch "<<XVU_Ch2[i][ind_best_Ch2[iBest]]<<endl; 	
	  }
	//  cout<<endl;	  
	}
      }
    }// if(Nseg_Ch2 > 0)

	  	  
    if(Nseg_Ch1 > 0) {ProcessSegments(1,  sigma,   dw_half,  kZ1_loc,  kMinHits,   Nseg_Ch1,  Nhits_Ch1,   Wires_Ch1,  clust_Ch1,   XVU_Ch1,  Nbest_Ch1,   ind_best_Ch1, Chi2_ndf_Ch1,  Chi2_ndf_best_Ch1, par_ab_Ch1,  
				     //   matrA, matrb, 
				     kNPlanes, ipl, XVU1, XVU_cl1,  kChi2_Max, dX_i1 );
      
     
      for (Int_t iBest = 0; iBest < Nbest_Ch1; iBest++) {	  
	if (Nhits_Ch1[ind_best_Ch1[iBest]] > 3) { 
	//  cout<<" iBest "<< iBest<<" Chi2_ndf_best_Ch1 "<<Chi2_ndf_best_Ch1[iBest]<<endl;
	  for (Int_t i = 0; i < 6; i++){
	//    cout<<" Ch= "<<1<<" XVU_Ch "<<XVU_Ch1[i][ind_best_Ch1[iBest]]<<endl; 	
	  }
//	  cout<<endl;
	}
      }
    

    }// if(Nseg_Ch1 > 0)
   
    //  cout<<endl;
    //  cout<<"ProcessSegments: Nbest_Ch1 "<<Nbest_Ch1<<" Nbest_Ch2 "<<Nbest_Ch2<<endl;
    //   cout<<endl;
     
	for (Int_t ise = 0; ise < Nbest_Ch1; ise++) {
	  //  cout<<" Ch1 ise "<<ise<<" ind "<<ind_best_Ch1[ise]<<" Chi2 "<<Chi2_ndf_best_Ch1[ise]<<" Ax "<<par_ab_Ch1[0][ise]<<" bx "<<par_ab_Ch1[1][ise]<<" Ay "<<par_ab_Ch1[2][ise]<<" by "<<par_ab_Ch1[3][ise]<<endl;//" kZ1 "<<Ch1Cent.Z()<<endl;

	
	  BmnTrack *pSeg = new ((*fBmnMwpcSegmentsArray)[fBmnMwpcSegmentsArray->GetEntriesFast()]) BmnTrack();
	  pSeg->SetChi2(Chi2_ndf_best_Ch1[ise]);
	  FairTrackParam pSegParams;
	  pSegParams.SetPosition(TVector3(par_ab_Ch1[1][ise], par_ab_Ch1[3][ise],ZCh1));
	  pSegParams.SetTx(par_ab_Ch1[0][ise]);
	  pSegParams.SetTy(par_ab_Ch1[2][ise]);
	  pSeg->SetParamFirst(pSegParams);
	}

	for (Int_t ise = 0; ise < Nbest_Ch2; ise++) {
	  //	  cout<<" Ch2 ise "<<ise<<" ind "<<ind_best_Ch2[ise]<<" Chi2 "<<Chi2_ndf_best_Ch2[ise]<<" Ax "<<par_ab_Ch2[0][ise]<<" bx "<<par_ab_Ch2[1][ise]<<" Ay "<<par_ab_Ch2[2][ise]<<" by "<<par_ab_Ch2[3][ise]<<endl;

	 
	  BmnTrack *pSeg1 = new ((*fBmnMwpcSegmentsArray)[fBmnMwpcSegmentsArray->GetEntriesFast()]) BmnTrack();
	  pSeg1->SetChi2(Chi2_ndf_best_Ch2[ise]);
	  FairTrackParam pSegParams1;
	  pSegParams1.SetPosition(TVector3(par_ab_Ch2[1][ise], par_ab_Ch2[3][ise],ZCh2));
	  pSegParams1.SetTx(par_ab_Ch2[0][ise]);
	  pSegParams1.SetTy(par_ab_Ch2[2][ise]);
	  pSeg1->SetParamFirst(pSegParams1);
	}

	

	//  cout<<endl;
	  
      hNbest_Ch1->Fill(Nbest_Ch1);	     
      hNbest_Ch2->Fill(Nbest_Ch2);

      if (Nbest_Ch1 > 0 && Nbest_Ch2 > 0){

      SegmentParamAlignment(Nbest_Ch1, ind_best_Ch1, par_ab_Ch1, shift1);
      SegmentParamAlignment(Nbest_Ch2, ind_best_Ch2, par_ab_Ch2, shift2);

     
    //  cout<<" After alignment "<<endl;
      for (Int_t ise = 0; ise < Nbest_Ch1; ise++) {
//	  cout<<" Ch1 ise "<<ise<<" ind "<<ind_best_Ch1[ise]<<" Chi2 "<<Chi2_ndf_best_Ch1[ise]<<" Ax "<<par_ab_Ch1[0][ind_best_Ch1[ise]]<<" bx "<<par_ab_Ch1[1][ind_best_Ch1[ise]]<<" Ay "<<par_ab_Ch1[2][ind_best_Ch1[ise]]<<" by "<<par_ab_Ch1[3][ind_best_Ch1[ise]]<<endl;
      }

      for (Int_t ise = 0; ise < Nbest_Ch2; ise++) {
//	  cout<<" Ch2 ise "<<ise<<" ind "<<ind_best_Ch2[ise]<<" Chi2 "<<Chi2_ndf_best_Ch2[ise]<<" Ax "<<par_ab_Ch2[0][ind_best_Ch2[ise]]<<" bx "<<par_ab_Ch2[1][ind_best_Ch2[ise]]<<" Ay "<<par_ab_Ch2[2][ind_best_Ch2[ise]]<<" by "<<par_ab_Ch2[3][ind_best_Ch2[ise]]<<endl;
      }
    //  cout<<endl;
    

      SegmentMatching( Nbest_Ch1, Nbest_Ch2, par_ab_Ch1, par_ab_Ch2, kZmid1, kZmid2, ind_best_Ch1, ind_best_Ch2, best_Ch1_gl, best_Ch2_gl, Nbest_Ch12_gl, Chi2_match);

      //  cout<<" SegmentMatching: Nbest_Ch12_gl "<<Nbest_Ch12_gl<<endl;

      }// if (Nbest_Ch1 > 0 && Nbest_Ch2 > 0){

      
     

      // ----spatial track---

      if ( Nbest_Ch12_gl > 0){ 

	//	cout<<" SegmentFit: matrA "<<matrA[0][0]<<endl;

	SegmentFit(z_gl1, z_gl2,  sigm2,   Nbest_Ch12_gl, ind_best_Ch1, ind_best_Ch2, best_Ch1_gl, best_Ch2_gl,  par_ab_Ch1_2, Chi2_ndf_Ch1_2, 	  Wires_Ch1, Wires_Ch2,  XVU_Ch1,  XVU_Ch2,   ind_best_Ch1_2, Nhits_Ch1, Nhits_Ch2);	

	for (Int_t ise = 0; ise < Nbest_Ch12_gl; ise++) {
	  //	  cout<<" before Alignment: 1-2 ise "<<ise<<" ind "<<ind_best_Ch1_2[ise]<<" Chi2 "<<Chi2_ndf_Ch1_2[ise]
	  //	      <<" Ax "<<par_ab_Ch1_2[0][ise]<<" bx "<<par_ab_Ch1_2[1][ise]<<" Ay "<<par_ab_Ch1_2[2][ise]<<" by "<<par_ab_Ch1_2[3][ise]<<endl;		  
	}  
  
	SegmentParamAlignment(Nbest_Ch12_gl, ind_best_Ch1_2, par_ab_Ch1_2, shift1_2);

	for (Int_t bst = 0; bst < Nbest_Ch12_gl; bst++) {

	  Float_t X_par_to_pole=par_ab_Ch1_2[0][bst]*kZ_to_pole+ par_ab_Ch1_2[1][bst];
	  Float_t Y_par_to_pole=par_ab_Ch1_2[2][bst]*kZ_to_pole+ par_ab_Ch1_2[3][bst];

	  //  cout<<" after Alignment: 1-2 ise "<<ise<<" ind "<<ind_best_Ch1_2[ise]<<" Chi2 "<<Chi2_ndf_Ch1_2[ise]
	  //    <<" Ax "<<par_ab_Ch1_2[0][ise]<<" bx "<<par_ab_Ch1_2[1][ise]<<" Ay "<<par_ab_Ch1_2[2][ise]<<" by "<<par_ab_Ch1_2[3][ise]<<endl;

	  //  cout<<" X_par_to_pole "<<X_par_to_pole<<" Y_par_to_pole "<< Y_par_to_pole<<" kZ_to_pole "<<kZ_to_pole<<endl;

	  hChi2_ch1_2->Fill(Chi2_ndf_Ch1_2[bst]);

	  BmnTrack *Tr = new ((*fBmnMwpcTracksArray)[fBmnMwpcTracksArray->GetEntriesFast()]) BmnTrack();
	  Tr->SetChi2(Chi2_ndf_Ch1_2[bst]);
	  FairTrackParam TrParams;
	  TrParams.SetPosition(TVector3(par_ab_Ch1_2[1][bst], par_ab_Ch1_2[3][bst],-kZ_to_pole));
	  TrParams.SetTx(par_ab_Ch1_2[0][bst]);
	  TrParams.SetTy(par_ab_Ch1_2[2][bst]);
	  Tr->SetParamFirst(TrParams);

	}//Nbest_Ch12_gl

      }//Nbest_Ch12_gl > 0)

*/
      // create a track and put in into TClonesArray:
      
      
          //   for (Int_t iPlane = 0; iPlane < kNPlanes; iPlane++)
          //  digits_filtered[iChamber][iPlane] = CheckDigits(digits[iChamber][iPlane]);

        // Z-coordinate of created hit is considered as known Z-position of planes 1 and 4 respecting to the considering chamber (1 or 2) 
      //  CreateMwpcHits(CreateHitsBy3Planes(digits_filtered[iChamber][0], digits_filtered[iChamber][1], digits_filtered[iChamber][2],
      //          fMwpcGeometry->GetZPlanePos(iChamber, 1)), fBmnMwpcHitArray, iChamber);
     //   CreateMwpcHits(CreateHitsBy3Planes(digits_filtered[iChamber][3], digits_filtered[iChamber][4], digits_filtered[iChamber][5],
     //	fMwpcGeometry->GetZPlanePos(iChamber, 4)), fBmnMwpcHitArray, iChamber);


      
      //    }
    if (fVerbose) cout << "\n======================== MWPC hit finder exec finished ====================" << endl;
    clock_t tFinish = clock();
    workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
   }//?
 
void BmnMwpcHitFinderSRC::SegmentFinder(
 Int_t chNum, Int_t*** wires_Ch, Int_t ***clust_Ch_, Float_t ***XVU_Ch_, 
 Int_t **Nhits_Ch_, Int_t **iw_Ch_, 
 Int_t *Nseg, Int_t ***wires_glob, Float_t ***xuv_glob, 
 Int_t minHits, Short_t code, Int_t kBig_ ) {
  //
 
  //xuv_glob     - coordinates of all hits
  //wires_glob   - wires of all hits
  Float_t delta = (chNum == 3) ? 1.125 : 0.75; //0.5; //condition for the nearest "lighted" wires
  Int_t cntr_sumWW_x =  95 ;//center wires Wx1+Wx2 due to x-slope
  Int_t min_sumWW_x =  (chNum == 3) ? 3 : 2; //min sum of wires to 95 

  Int_t min_sumWW =  (chNum == 3) ? 3 : 2; //min sum of wires to 95  //wide UV

  Int_t cntr_sumWW =  95 ;//center wires W1+W2 due to u,v-slope
  //Int_t min_sumWW =  2 ; //min sum of wires to 95 //narrow UV

  Int_t minHits4_5= minHits;

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

  if (Nseg[chNum] > kBig_ - 2)return;// MP

    if (iw_Ch_[chNum][x] > 0) {
    for (Int_t ix = 0; ix < iw_Ch_[chNum][x]; ix++) {
      if (iw_Ch_[chNum][v] > 0) {
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
		 
		  //		  ::out1<<" ch "<<chNum<<" case "<<code<<" iseg "<<iseg<<" ix "<<x<<" wr "<<wires_Ch[chNum][x][iseg]<<" iu "<<u<<" wr "<<wires_Ch[chNum][u][iseg]<<" iv "<<v<<" wr "<<wires_Ch[chNum][v][iseg]<<endl;

		  if (wires_Ch[chNum][x][iseg] == wires_glob[chNum][ix][x]) it_was_x = 1;
		  if ((clust_Ch_[chNum][x][iseg] == -1 || clust_Ch_[chNum][x][iseg] == 3) && (wires_Ch[chNum][x][iseg] -1 == wires_glob[chNum][ix][x])) it_was_x = 1; 
		  if ((clust_Ch_[chNum][x][iseg] ==  1 || clust_Ch_[chNum][x][iseg] == 3) && (wires_Ch[chNum][x][iseg] +1 == wires_glob[chNum][ix][x])) it_was_x = 1;
		  if (wires_Ch[chNum][v][iseg] == wires_glob[chNum][iv][v]) it_was_v = 1;
		  if ((clust_Ch_[chNum][v][iseg] == -1 || clust_Ch_[chNum][v][iseg] == 3) && wires_Ch[chNum][v][iseg] -1 == wires_glob[chNum][iv][v]) it_was_v = 1; 
		  if ((clust_Ch_[chNum][v][iseg] ==  1 || clust_Ch_[chNum][v][iseg] == 3) && wires_Ch[chNum][v][iseg] +1 == wires_glob[chNum][iv][v]) it_was_v = 1;
		  if (wires_Ch[chNum][u][iseg] == wires_glob[chNum][iu][u]) it_was_u = 1;
		  if ((clust_Ch_[chNum][u][iseg] == -1 || clust_Ch_[chNum][u][iseg] == 3) && wires_Ch[chNum][u][iseg] -1 == wires_glob[chNum][iu][u]) it_was_u = 1; 
		  if ((clust_Ch_[chNum][u][iseg] ==  1 || clust_Ch_[chNum][u][iseg] == 3) && wires_Ch[chNum][u][iseg] +1 == wires_glob[chNum][iu][u]) it_was_u = 1;				 
		  //   bef.clustering    if (wires_Ch[chNum][x][iseg] == wires_glob[chNum][ix][x] && wires_Ch[chNum][v][iseg] == wires_glob[chNum][iv][v] && wires_Ch[chNum][u][iseg] == wires_glob[chNum][iu][u]) {
		  if (it_was_x * it_was_v * it_was_u)  { it_was = 1;//{ it_was = 1;
		     break;
		  }
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

		    if ( double_u2 > 0 && Nseg[chNum] <= 99 ) {
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

		    if ( double_v2 > 0 && Nseg[chNum] <= 99 ) { 
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

		    if ( double_x2 > 0 && Nseg[chNum] <= 99 ) { 
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

		  //
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

/*


void BmnMwpcHitFinderSRC::ProcessSegments(
                             Int_t chNum, Double_t sigma_, Float_t dw_half_,
			     Float_t **z_loc, Int_t Min_hits, Int_t *Nseg, 
			     Int_t **Nhits, Int_t ***wires_Ch, Int_t ***clust_Ch_,
			     Float_t ***XVU_Ch_, Int_t *Nbest, Int_t **ind_best,			 
			     Double_t **Chi2_ndf, Double_t **Chi2_ndf_best,
			     Double_t ***par_ab,  Int_t nPlanes, Int_t* ipl_,
			     Float_t** XVU_, Float_t** XVU_cl_,  Double_t kChi2_Max_ ) {

  // cout<<"start "<< chNum <<" "<< sigma_ <<" "<< dw_half_<<" "<<z_loc[0]<<" "
  //    << Min_hits <<" "<<Nseg[chNum]<< Nhits[chNum][1]<<" "<<wires_Ch[0][0]<<" "<<Nbest[chNum]<<" "<<ind_best[chNum][0]<<" "
  //    <<Chi2_ndf[chNum][0]<<" "<<Chi2_ndf_best[chNum][0]<<" "<<par_ab[chNum][0][0]<<" "<<A_[0][0]<<" "<<b_[0][0]<<" "<<nPlanes<<ipl_[0]<<" "<<XVU_[chNum][0]<<" "<<XVU_cl_[chNum][0]<<" "<<kChi2_Max_<<" "
  //    <<endl;

  
  Float_t delta = (chNum == 3) ? 1.125 : 0.75; //0.5;

  Float_t sigma2 = sigma_ * sigma_;
  //  Float_t sigm2[nPlanes] = {sigma2, sigma2, sigma2, sigma2, sigma2, sigma2};

  //  Int_t h6[nPlanes] = {1, 1, 1, 1, 1, 1};   
  Int_t h6[6] = {1, 1, 1, 1, 1, 1};   
  Int_t Min_hits6 = Min_hits;
  for (Int_t Nhitm = nPlanes; Nhitm > Min_hits - 1; Nhitm--) {
    Int_t ifNhitm = 0;

    if (Nhitm < Min_hits6) break;

    for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {

      //  cout<<" chNum "<<chNum<<"-- iseg "<<iseg<<" Nhits "<<Nhits[chNum][iseg]<<" Chi2_ndf "<<Chi2_ndf[chNum][iseg]<<endl;

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
	  if (wires_Ch[i][iseg] > -1) { h[i] = 1;
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
	  if (wires_Ch[i][iseg] < 0 ) continue;
	  if ( abs(clust_Ch_[chNum][i][iseg]) == 1) { ifirst++; ipl_[ifirst] = i;}
	}
      }// Nhits[chNum][iseg] > 4 
      else { for (Int_t i = 0; i < nPlanes; i++){ 
	  if (wires_Ch[i][iseg] > -1)  h[i] = 1;
	  ipl_[i]=i;
	}
      }
     
    //  cout<<endl;
    //  cout<<" ::ProcessSegments "<<" iseg "<<iseg<<endl;
      for (Int_t i = 0; i < nPlanes; i++){
    //  	cout<<" Ch= "<<chNum<<" h "<<h[i]<<" XVU_Ch_ "<<XVU_Ch_[chNum][i][iseg]<<endl;   	
      }
    //  cout<<endl;
     
      //  cout<<" iseg "<<iseg<<" clust_Ch_[chNum] "<<clust_Ch_[chNum][0][iseg] <<" "<<clust_Ch_[chNum][1][iseg] <<" "<<clust_Ch_[chNum][2][iseg] <<" "<<clust_Ch_[chNum][3][iseg] <<" "<<clust_Ch_[chNum][4][iseg]<<" "<<clust_Ch_[chNum][5][iseg]
      //	<<" wires_Ch "<<wires_Ch[0][iseg]<<" "<<wires_Ch[1][iseg]<<" "<<wires_Ch[2][iseg]<<" "<<wires_Ch[3][iseg]<<" "<<wires_Ch[4][iseg]<<" "<<wires_Ch[5][iseg]
      //	  <<" XVU_Ch_ "<<XVU_Ch_[chNum][0][iseg]<<" "<<XVU_Ch_[chNum][1][iseg]<<" "<<XVU_Ch_[chNum][2][iseg]<<" "<<XVU_Ch_[chNum][3][iseg]<<" "<<XVU_Ch_[chNum][4][iseg]<<" "<<XVU_Ch_[chNum][5][iseg]
      //	   <<" ipl "<< ipl_[0] <<" "<< ipl_[1] <<" "<< ipl_[2] <<" "<< ipl_[3] <<" "<< ipl_[4] <<" "<< ipl_[5] <<endl;
    

      //linear fit
      //     Double_t A_[4][4]; //coef matrix
      // Double_t F[4]; //free coef 
      
      Float_t dX[nPlanes];
      //      Float_t XVU[nPlanes];
      //      Float_t XVU_cl[nPlanes];
      Double_t Chi2_min_iseg = 999;
      Double_t par_ab_curr[4];
      Double_t par_ab_cl[4];
      
      for (Int_t i10 = min_i[ipl_[0]]; i10 <= max_i[ipl_[0]] ; i10++){ 
	//-?	if ( ipl_[0] == 6 ) continue;
	if ( h[ipl_[0]] == 1){
	  XVU_[chNum][ipl_[0]] = XVU_Ch_[chNum][ipl_[0]][iseg] + dw_half_*i10 ;
	  sigm2[ipl_[0]] = sigma2 ;
	  if ( i10 == -1 || i10 == 1 )  sigm2[ipl_[0]] = 0.5*sigma2 ;
	}
	for (Int_t i2 = min_i[ipl_[1]]; i2 <= max_i[ipl_[1]] ; i2++){ 
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

		  Float_t dX_i[nPlanes] = {0,0,0, 0,0,0};
		  Double_t Chi2_curr_iseg = 0;

		  //Att
		 
		  for (Int_t i1 = 0; i1 < nPlanes; i1++) {
		    if (wires_Ch[i1][iseg]>-1) {
		      if (i1 == 0 || i1 == 3) dX_i[i1] = XVU_[chNum][i1] - par_ab_curr[0] * z_loc[i1] - par_ab_curr[1];
		      if (i1 == 2 || i1 == 5) dX_i[i1] = XVU_[chNum][i1] - 0.5 * (par_ab_curr[0] + sq3 * par_ab_curr[2]) * z_loc[i1] - 0.5 * (par_ab_curr[1] + sq3 * par_ab_curr[3]);
		      if (i1 == 1 || i1 == 4) dX_i[i1] = XVU_[chNum][i1] - 0.5 * (par_ab_curr[0] - sq3 * par_ab_curr[2]) * z_loc[i1] - 0.5 * (par_ab_curr[1] - sq3 * par_ab_curr[3]);
		      Chi2_curr_iseg =  Chi2_curr_iseg + dX_i[i1] * dX_i[i1] / sigm2[i1];//??
		    } 
		  }//i1
		  		 
		  //	  cout<<" iseg "<<iseg<<"  --- Chi2_curr_iseg "<<Chi2_curr_iseg<<" Chi2_min_iseg "<<Chi2_min_iseg<<endl;

		  //??	  if ( Chi2_curr_iseg > Chi2_i6) break;              // optimization
		  if ( Chi2_curr_iseg > Chi2_min_iseg * 10. ) break; // optimization
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
	      }
	    }
	  }
	}//i2
      }//i10

      Chi2_ndf[chNum][iseg] = Chi2_min_iseg; 

      // cout<<"1 iseg "<<iseg<<" Chi2_ndf "<<Chi2_ndf[chNum][iseg]<<endl;

      if (Nhits[chNum][iseg] > 4)  Chi2_ndf[chNum][iseg] = Chi2_ndf[chNum][iseg] / (Nhits[chNum][iseg] - 4);

      if (Chi2_ndf[chNum][iseg] > kChi2_Max_) {
	if (Nhits[chNum][iseg] <= Min_hits) { Nhits[chNum][iseg] = 0; continue;}
	else  { //reject most distant point
	  Float_t Max_dev = 0;
	  Int_t irej = -1;
	  for (Int_t i1 = 0; i1 < nPlanes; i1++)
	    if (wires_Ch[i1][iseg]>-1)
	      if (fabs(dX[i1]) > Max_dev) {
		irej = i1;
		Max_dev = fabs(dX[i1]);
	      }
	  wires_Ch[irej][iseg] = -1;
	  Nhits[chNum][iseg]--;// reject most distant point
	  continue;}
      }
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

      Double_t Chi2_best = 9999;
      Int_t iseg_best = -1;
      for (Int_t iseg = 0; iseg < Nseg[chNum]; iseg++) {
	if (Nhits[chNum][iseg] != Nhitm) continue;
	if (Chi2_ndf[chNum][iseg] >= Chi2_best) continue;
	Chi2_best = Chi2_ndf[chNum][iseg];
	iseg_best = iseg;
      } // iseg

      if (iseg_best == -1) continue;
     
      //  cout<<"3 Ch= "<< chNum <<" iseg_best1 "<< iseg_best <<" Chi2_best "<< Chi2_best<<endl;
	//  <<" wires_Ch "<< wires_Ch[0][iseg_best]<<" "<< wires_Ch[1][iseg_best]<<" "<< wires_Ch[2][iseg_best]<<" "<< wires_Ch[3][iseg_best]<<" "<< wires_Ch[4][iseg_best]<<" "<< wires_Ch[5][iseg_best]
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
	  if (wires_Ch[i1][iseg]>-1) {
	   
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

	//	cout<<" Ch= "<< chNum <<" iseg_best2 "<< iseg_best2 <<" Chi2_best "<< Chi2_best<<" wires_Ch "<< wires_Ch[0][iseg_best2]<<" "
	//   <<wires_Ch[1][iseg_best2]<<" "<<wires_Ch[2][iseg_best2]<<" "<<wires_Ch[3][iseg_best2]<<" "<<wires_Ch[4][iseg_best2]<<" "<<wires_Ch[5][iseg_best2]<<endl;
		
	
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
	    if (wires_Ch[i1][iseg]>-1) {
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
	    if (wires_Ch[i1][iseg]>-1) {
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
	    if (wires_Ch[i1][iseg]>-1) {
	      if( fabs(XVU_Ch_[chNum][i1][iseg] - XVU_Ch_[chNum][i1][iseg_best3]) < 3*dw_half ) Nhits[chNum][iseg] = 0;
	    }
	  }
	}
      }

    }

    //   cout<<"Ch= "<<chNum<<" Nhits[chNum](Two-Five bests) ";  
    //   	for (int iseg=0; iseg< Nseg[chNum]; iseg++){ 
    //   	  cout<<Nhits[chNum][iseg]<<" "; 
    //    	} cout<<endl;
    //	cout<<"Ch= "<<chNum<<" Nbest[chNum] "<<Nbest[chNum]<<endl;
    //	for (int iseg=0; iseg< Nseg[chNum]; iseg++){ 	
	  // cout<<" Ch= "<<chNum<< " ind "<<ind_best[chNum][iseg]<<" Chi2 "<<Chi2_ndf_best[chNum][ind_best[chNum][iseg]]<<endl;
    //	}

    
	//      << Min_hits <<" "<<Nseg[chNum]<< Nhits[chNum][1]<<" "<<wires_Ch[0][0]<<" "<<Nbest[chNum]<<" "<<ind_best[chNum][0]<<" "
	//   <<Chi2_ndf[chNum][0]<<" "<<Chi2_ndf_best[chNum][0]<<" "<<par_ab[chNum][0][0]<<" "<<A_[0][0]<<" "<<b_[0][0]<<" "<<nPlanes<<ipl_[0]<<" "<<XVU_[chNum][0]<<" "<<XVU_cl_[chNum][0]<<" "<<kChi2_Max_<<" "
	//   <<endl;
	
    
  } // Nhitm
  
}// ProcessSegments


void BmnMwpcHitFinderSRC::SegmentParamAlignment(Int_t & Nbest, Int_t *ind_best, Double_t **par_ab, Float_t *shift ){
 
  //local parameters to Global parameters

  for (Int_t iBest = 0; iBest < Nbest; iBest++) {

    // cout<<"before Alignment: iBest "<<iBest<<" Ax "<< par_ab[chNum][0][ind_best[chNum][iBest]]<<" bx "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" Ay "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" by "<< par_ab[chNum][3][ind_best[chNum][iBest]]<<endl;

    //                                     ax          alpha                                      ax^2    
    par_ab[chNum][0][ind_best[chNum][iBest]] += shift[0] +  shift[0]* par_ab[chNum][0][ind_best[chNum][iBest]]* par_ab[chNum][0][ind_best[chNum][iBest]];
    par_ab[chNum][2][ind_best[chNum][iBest]] += shift[2] +  shift[2]* par_ab[chNum][2][ind_best[chNum][iBest]]* par_ab[chNum][2][ind_best[chNum][iBest]];
    par_ab[chNum][1][ind_best[chNum][iBest]] += shift[1];
    par_ab[chNum][3][ind_best[chNum][iBest]] += shift[3];

    //cout<<"after Alignment: iBest "<<iBest<<" Ax "<< par_ab[chNum][0][ind_best[chNum][iBest]]<<" bx "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" Ay "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" by "<< par_ab[chNum][3][ind_best[chNum][iBest]]<<endl;

  }//iBest 
}//SegmentParamAlignment

void BmnMwpcHitFinderSRC::SegmentMatching(  Int_t & Nbest_Ch1_, Int_t & Nbest_Ch2_, Double_t **par_ab_Ch1_,  Double_t **par_ab_Ch2_, Float_t Zmid1, Float_t Zmid2, Int_t *ind_best_Ch1_, Int_t *ind_best_Ch2_, Int_t *best_Ch1_gl_, Int_t *best_Ch2_gl_, Int_t & Nbest_Ch12_gl_, Double_t *Chi2_match_){

  //  cout<<" Nbest_Ch1 "<<Nbest_Ch1_<<" Ch2 "<<Nbest_Ch2_<<" par_ab_Ch1 "<<par_ab_Ch1_[0][0]<<" par_ab_Ch2 "<<par_ab_Ch2_[0][0]<<" Zmid1 "<<Zmid1<<" Zmid2 "<<Zmid2<<" ind_best_Ch1 "<<ind_best_Ch1_[0]<<" ind_best_Ch2 "<<ind_best_Ch2_[0]<<" ind best for match ch1 "<<best_Ch1_gl_[0]<<" ch2 "<<best_Ch2_gl_[0]<<" Nbest_Ch12_gl "<<Nbest_Ch12_gl_<< " chi2_match "<<Chi2_match_[0]<<endl;

    Int_t best = -1;
    Int_t best3 = -1;
 
    Float_t sig_dx= 3.6; //0.8;//22;- z=0 //0.85;- Zmid
    Float_t sig_dy= 3.6; //0.7;//18;- z=0 //0.76;
    Float_t sig_dax= 0.055; //0.04; //0.063;
    Float_t sig_day= 0.055; //0.04; //0.045;
	    
    Float_t min_Chi2m = 100; // 40; //100; //400
    Float_t min_distX = 99;
    Float_t min_distY = 99;
    Float_t dAx12 = 0;
    Float_t dAy12 = 0;

    Float_t DAx12 = 0;
    Float_t DAy12 = 0;
    Float_t Min_distX = 0;
    Float_t Min_distY = 0;
	    
    for (Int_t bst1 = 0; bst1 < Nbest_Ch1_; bst1++) {
      
      //ch1                                                 zloc0 -z_i
      Float_t x1mid = par_ab_Ch1_[0][ind_best_Ch1_[bst1]] *( 0 - Zmid1) + par_ab_Ch1_[1][ind_best_Ch1_[bst1]] ;
      Float_t y1mid = par_ab_Ch1_[2][ind_best_Ch1_[bst1]] *( 0 - Zmid1) + par_ab_Ch1_[3][ind_best_Ch1_[bst1]] ;

      //  cout<<"par Ch1 bst1 "<<bst1<<" Ax "<<par_ab_Ch1_[0][ind_best_Ch1_[bst1]]<<" bx "<<par_ab_Ch1_[1][ind_best_Ch1_[bst1]]<<" Ay "<<par_ab_Ch1_[2][ind_best_Ch1_[bst1]]<<" by "<<par_ab_Ch1_[3][ind_best_Ch1_[bst1]]<<endl;
                 

      for (Int_t bst2 = 0; bst2 < Nbest_Ch2_; bst2++) {
	
	//ch2       
	Float_t x2mid =  par_ab_Ch2_[0][ind_best_Ch2_[bst2]] *( 0 - Zmid2)  + par_ab_Ch2_[1][ind_best_Ch2_[bst2]] ;
	Float_t y2mid =  par_ab_Ch2_[2][ind_best_Ch2_[bst2]] *( 0 - Zmid2)  + par_ab_Ch2_[3][ind_best_Ch2_[bst2]] ;
	
	//	cout<<"par Ch2 bst2 "<<bst2<<" Ax "<<par_ab_Ch2_[0][ind_best_Ch2_[bst2]]<<" bx "<<par_ab_Ch2_[1][ind_best_Ch2_[bst2]]<<" Ay "<<par_ab_Ch2_[2][ind_best_Ch2_[bst2]]<<" by "<<par_ab_Ch2_[3][ind_best_Ch2_[bst2]]<<endl;


	// cout<<" x2mid "<<x2mid<<" y2mid "<<y2mid<<"  ind_bst1 "<<ind_best_Ch1_[bst1]<<" ind_bst2 "<<ind_best_Ch2_[bst2]<<endl;
		  		  
	dAx12 = par_ab_Ch1_[0][ind_best_Ch1_[bst1]] - par_ab_Ch2_[0][ind_best_Ch2_[bst2]];
	dAy12 = par_ab_Ch1_[2][ind_best_Ch1_[bst1]] - par_ab_Ch2_[2][ind_best_Ch2_[bst2]];
	min_distX = x1mid - x2mid; //min
	min_distY = y1mid - y2mid; //min

	//  hdX_Zmid_Ch12->Fill(min_distX);
	//   hdY_Zmid_Ch12->Fill(min_distY);
	//   hdAx12->Fill(dAx12);
	//   hdAy12->Fill(dAy12);

	Float_t Chi2_m = ( min_distX*min_distX/(sig_dx*sig_dx) 
			   +  min_distY*min_distY/(sig_dy*sig_dy)
			   +  dAx12*dAx12 /(sig_dax*sig_dax)
			   +  dAy12*dAy12 /(sig_day*sig_day) );

	//	 cout<<" bst1 "<<bst1<<" bst2 "<<bst2<<" min_distX "<<min_distX<<" min_distY "<<min_distY <<" dAx12 "<<dAx12<<" dAy12 "<<dAy12<<" Chi2_m "<<Chi2_m<<endl;
	
	  if (Chi2_m < min_Chi2m && Chi2_m < Chi2_match_[0]){
	    //  cout<<"Matching: " << bst1 << " : " <<  bst2 <<endl;

	    DAx12=dAx12;
	    DAy12=dAy12;
	    Min_distX= min_distX;
	    Min_distY= min_distY;

	    Nbest_Ch12_gl_++;

	    for (int j = 4; j > 0; j--) {
	      best_Ch1_gl_[j] = best_Ch1_gl_[j-1];
	      best_Ch2_gl_[j] = best_Ch2_gl_[j-1];
	      Chi2_match_[j] =  Chi2_match_[j-1];
	    }
	    best_Ch1_gl_[0] = bst1;
	    best_Ch2_gl_[0] = bst2;
	    Chi2_match_[0] =  Chi2_m;
	    //  cout<<"Matching: " << Chi2_m << " : " << min_Chi2m <<endl;
	    //   cout<<"Matching: " << best_Ch1_gl_[0]  << " : " << best_Ch2_gl_[0]  <<endl;
	  }//if (Chi2_m < min_Chi2m

		

      }//bst2++     

      //   cout<<" Nbest_Ch12_gl "<< Nbest_Ch12_gl_<<" best_Ch1 "<<best_Ch1_gl[bst1]<<" best_Ch2 "<< best_Ch2_gl[bst1]<<" Chi2_match "<<Chi2_match_[bst1]<<endl;
 
    }//bst1++

     // cout<<" Nbest_Ch1 "<<Nbest_Ch1_<<" Ch2 "<<Nbest_Ch2_<<" par_ab_Ch1 "<<par_ab_Ch1_[0][0]<<" par_ab_Ch2 "<<par_ab_Ch2_[0][0]<<" Zmid1 "<<Zmid1<<" Zmid2 "<<Zmid2<<" ind_best_Ch1 "<<ind_best_Ch1_[0]<<" ind_best_Ch2 "<<ind_best_Ch2_[0]<<" ind best for match ch1 "<<best_Ch1_gl_[0]<<" ch2 "<<best_Ch2_gl_[0]<<" Nbest_Ch12_gl "<<Nbest_Ch12_gl_<< " chi2_match "<<Chi2_match_[0]<<endl;
}// SegmentMatching



void BmnMwpcHitFinderSRC::SegmentFit(Float_t *z_gl1_, Float_t *z_gl2_, 
				     Float_t *sigm2_,
				     Int_t & Nbest_Ch12_gl_,
				     Int_t *ind_best_Ch1_, Int_t *ind_best_Ch2_, 
				     Int_t *best_Ch1_gl_, Int_t *best_Ch2_gl_, 
				     Double_t **par_ab_Ch1_2_,  Double_t * Chi2_ndf_Ch1_2_, 
				     //    Int_t *h1, Int_t *h2, 
				     Int_t **Wires_Ch1_, Int_t **Wires_Ch2_,
				     // Double_t **Amatr,  Double_t **bmatr,
				     Float_t **XVU_Ch1_,Float_t **XVU_Ch2_, 
				     //   Float_t *XVU1_, Float_t *XVU2_,
				     Int_t *ind_best_Ch1_2_, Int_t *Nhits_Ch1_, Int_t *Nhits_Ch2_ ){
  for (Int_t bst = 0; bst < Nbest_Ch12_gl_; bst++) {

    Int_t bst1 = best_Ch1_gl_[bst];
    Int_t bst2 = best_Ch2_gl_[bst];
	
    Int_t best1 = ind_best_Ch1_[bst1]; // it's segment index !
    Int_t best2 = ind_best_Ch2_[bst2]; // it's segment index!

    // cout<<" SegmentFit: best1 "<<best1<<" best2 "<<best2<<endl;
    
    int h1[6] = {0,0,0,0,0,0};
    int h2[6] = {0,0,0,0,0,0};
    for(Int_t i = 0; i<6;i++){
      if(Wires_Ch1_[i][best1]> -1) h1[i] = 1;
      if(Wires_Ch2_[i][best2]> -1) h2[i] = 1;
    }

   int vkh1[6] = {0,0,0,0,0,0};
   int vkh2[6] = {0,0,0,0,0,0};
    for(Int_t i = 0; i<6;i++){
       if(Wires_Ch1_[i][best1]> -1) vkh1[i] = 1;
       if(Wires_Ch2_[i][best2]> -1) vkh2[i] = 1;
    }

    for(Int_t i = 0; i < 6;i++){
      //  sigm2_gl2[i]=(sigma*sigma);
      //  if ( clust_Ch_2[i]!=0 ) 
      sigm2_[i]= 0.00260417;
	//sigm2_[i]*0.5; // for clust 		  
    }
  
    for(Int_t i = 0; i < 6;i++){
//	cout<<" sigm2_gl1 "<<sigm2_[i]<<" sigm2_gl2 "<<sigm2_[i]<<endl;
    }
  
    //  cout<<endl;
    Amatr = new Double_t*[4];
    bmatr = new Double_t*[4];

    for(Int_t ii=0; ii<4; ii++){
      Amatr[ii] = new Double_t[4];
      bmatr[ii] = new Double_t[4];
    }

    for(Int_t im=0; im<4; im++){
      for(Int_t ii=0; ii<4; ii++){
	Amatr[im][ii] = 0.;
	bmatr[im][ii] = 0.;
      }
    }


   

    FillFitMatrix(Amatr, z_gl1_, sigm2_, h1);   //Ch_1
    FillFitMatrix(Amatr, z_gl2_, sigm2_, h2);   //Ch_2
   
   vector<vector<double>> Am1(4, vector<double>(4));
   vector<vector<double>> Am2(4, vector<double>(4));   
   
   Am1 = vkFillFitMatrix(z_gl1_, sigm2_, vkh1);   //Ch_1
   Am2 = vkFillFitMatrix(z_gl2_, sigm2_, vkh2);   //Ch_2
    
    for(int iterR = 0; iterR < 4; iterR++){
        for(int iterC = 0; iterC < 4; iterC++){

	  //  Amatr[iterR][iterC] = Am1[iterR][iterC] + Am2[iterR][iterC];
	  // cout<<"vk A["<< iterR <<"]" <<"["<< iterC <<"] " << Am1[iterR][iterC] << " : " << Am2[iterR][iterC] << " : " << Am1[iterR][iterC] + Am2[iterR][iterC]<<endl;
        }
    }
   
    for(int ff = 0; ff < 4; ff++){
        for(int tt = 0; tt < 4; tt++){
//	  cout<<"1 Amatr ["<<ff <<"] ["<<tt<<"] "<<Amatr[ff][tt]<<endl;
	}
    }
   
   
   Double_t matrF[4] = {0,0,0,0};//free coef 
   
   cout<<endl;
   for (Int_t i1 = 0; i1 < 4; i1++){
   //  cout<<"0 i1 "<<i1<<" "<< matrF[i1]<<endl;
   }
   
 //  cout<<" before FillFreeCoefVector "<<endl;
   for(Int_t i = 0; i < 6;i++){
    //   cout<<" sigm2 "<< sigm2_[i]<<" h1 "<<h1[i]<<" XVU_Ch1_ "<<XVU_Ch1_[i][best1]<<" h2 "<<h2[i]<<" XVU_Ch2_ "<<XVU_Ch2_[i][best1]<<endl;
   }  
  
   FillFreeCoefVector(matrF, XVU_Ch1_, best1, z_gl1_ , sigm2_, h1, 6);//   FillFreeCoefVector(matrF, XVU1_ , z_gl1_ , sigm2_, h1, 6);
   FillFreeCoefVector(matrF, XVU_Ch2_, best2, z_gl2_ , sigm2_, h2, 6);//  FillFreeCoefVector(matrF, XVU2_ , z_gl2_ , sigm2_, h2, 6);

  
   // cout<<endl;
   for (Int_t i1 = 0; i1 < 4; i1++){
   //   cout<<"2 i1 "<<i1<<" F "<< matrF[i1]<<endl;
   }
  // cout<<endl;
  

  

  //Gaussian algorithm for 4x4 matrix inversion 
	Double_t A0matr[4][4];	 			
	for (Int_t i1 = 0; i1 < 4; i1++){
	  for (Int_t j1 = 0; j1 < 4; j1++){
	    A0matr[i1][j1] = Amatr[i1][j1];
	  }
	}

	InverseMatrix(Amatr,bmatr);			 			  
	
	  Double_t sum;		  
	  Double_t A1[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};  
	  //	  cout<<" A1 "<<endl;

	  for (Int_t i1 = 0; i1 < 4; ++i1) 
	    for (Int_t j1 = 0; j1 < 4; ++j1) {
	      sum = 0; 
	      for (Int_t k1 = 0; k1 < 4; ++k1) {
		Double_t a0 = A0matr[i1][k1];
		Double_t b0 = bmatr[k1][j1];        
		sum += a0 * b0;       
		A1[i1][j1] = sum;
	      }
	      //  cout<<A1[i1][j1]<<" ";
	    }
	 
	  for (Int_t i1 = 0; i1 < 4; i1++){
	    for (Int_t j1 = 0; j1 < 4; j1++){
	      //   cout<<"2 Amatr[i1][j1] "<<Amatr[i1][j1]<<endl;
	    }
	  }

	  for (Int_t i1 = 0; i1 < 4; i1++){
	    for (Int_t j1 = 0; j1 < 4; j1++){
	      //  cout<<"2 bmatr[i1][j1] "<<bmatr[i1][j1]<<endl;
	    }
	  }
	  

	  for(Int_t i1 = 0 ; i1 < 4; i1++){
	    par_ab_Ch1_2_[i1][bst] = 0;
	    for(Int_t j1 = 0; j1 < 4; j1++){
	      par_ab_Ch1_2_[i1][bst] += bmatr[i1][j1]*matrF[j1];// par_ab_Ch23[i1][bst] += b[i1][j1]*F[j1];
	      //  cout<<" i1 "<<i1<<" bmatr "<<bmatr[i1][j1]<<" F "<<matrF[j1] <<endl;
	    }    
	  } // i1
	
	//  cout<<endl;
//	  cout<<" par [0] "<<par_ab_Ch1_2_[0][bst]<<" par [1] "<<par_ab_Ch1_2_[1][bst]<<" par [2] "<<par_ab_Ch1_2_[2][bst]<<" par [3] "<<par_ab_Ch1_2_[3][bst]<<endl;
//	  cout<<endl;
	  

	  Float_t sigm_gl1[12];
	  Float_t sigm_gl2[12];
	  for(Int_t i1 = 0 ; i1 < 6; i1++){
	    sigm_gl1[i1]=0.0722;
	    sigm_gl2[i1]=0.0722;
	  }

		  
	  Float_t dx_[kNPlanes];

	  Chi2_ndf_Ch1_2_[bst]=0;
	  for(Int_t i1 = 0 ; i1 < 6; i1++){
	    //XVU_Ch1
	    if(Wires_Ch1_[i1][best1]>-1){
	      if(i1==0 || i1==3) dx_[i1]=XVU_Ch1_[i1][best1]-par_ab_Ch1_2_[0][bst]*z_gl1_[i1]-par_ab_Ch1_2_[1][bst];
	      if(i1==2 || i1==5) dx_[i1]=XVU_Ch1_[i1][best1]-0.5*(par_ab_Ch1_2_[0][bst]+sq3*par_ab_Ch1_2_[2][bst])*z_gl1_[i1]-0.5*(par_ab_Ch1_2_[1][bst]+sq3*par_ab_Ch1_2_[3][bst]);
	      if(i1==1 || i1==4) dx_[i1]=XVU_Ch1_[i1][best1]-0.5*(par_ab_Ch1_2_[0][bst]-sq3*par_ab_Ch1_2_[2][bst])*z_gl1_[i1]-0.5*(par_ab_Ch1_2_[1][bst]-sq3*par_ab_Ch1_2_[3][bst]);
	      Chi2_ndf_Ch1_2_[bst]= Chi2_ndf_Ch1_2_[bst]+dx_[i1]*dx_[i1]/(sigm_gl1[i1]*sigm_gl1[i1]);

	      //  cout<<"best1 "<<best1 <<" i1 "<<i1<<" dx_ "<<dx_[i1]<<" XVU_Ch1 "<<XVU_Ch1_[i1][best1]<<" Chi2_ndf_Ch1_2 "<<Chi2_ndf_Ch1_2_[bst]<<" z_gl1 "<<z_gl1_[i1]<<endl;

	    }// if( Wires_Ch1[i1][best2]>-1){
	  }
	  //  cout<<endl;

	  for(Int_t i2 = 0 ; i2 < 6; i2++){
	    //XVU_Ch2
	    if(Wires_Ch2_[i2][best2]>-1){
	      if(i2==0 || i2==3) dx_[i2]=XVU_Ch2_[i2][best2]-par_ab_Ch1_2_[0][bst]*z_gl2_[i2]-par_ab_Ch1_2_[1][bst];
	      if(i2==2 || i2==5) dx_[i2]=XVU_Ch2_[i2][best2]-0.5*(par_ab_Ch1_2_[0][bst]+sq3*par_ab_Ch1_2_[2][bst])*z_gl2_[i2]-0.5*(par_ab_Ch1_2_[1][bst]+sq3*par_ab_Ch1_2_[3][bst]);
	      if(i2==1 || i2==4) dx_[i2]=XVU_Ch2_[i2][best2]-0.5*(par_ab_Ch1_2_[0][bst]-sq3*par_ab_Ch1_2_[2][bst])*z_gl2_[i2]-0.5*(par_ab_Ch1_2_[1][bst]-sq3*par_ab_Ch1_2_[3][bst]);
	      Chi2_ndf_Ch1_2_[bst]= Chi2_ndf_Ch1_2_[bst]+dx_[i2]*dx_[i2]/(sigm_gl2[i2]*sigm_gl2[i2]);

	      //    cout<<"best2 "<<best2 <<" i2 "<<i2<<" dx_ "<<dx_[i2]<<" XVU_Ch2 "<<XVU_Ch2_[i2][best2]<<" Chi2_ndf_Ch1_2 "<<Chi2_ndf_Ch1_2_[bst]<<" z_gl2 "<<z_gl2_[i2]<<endl;

	    }// if( Wires_Ch2[i2][best2]>-1){
	  }
	  //  cout<<endl;
	  
	  // Hist_Nhits_Ch2->Fill(Nhits_Ch2_[best2]);
	  // Hist_Nhits_Ch1->Fill(Nhits_Ch1_[best1]);
	  // Hist_Nhits_Ch12->Fill(Nhits_Ch1_[best1]+Nhits_Ch2_[best2]);

	  //  cout<<" bst: Nhits_Ch1 "<<Nhits_Ch1_[best1]<<" ch2 "<<Nhits_Ch2_[best2]<<endl;

	   if (Nhits_Ch1_[best1]+Nhits_Ch2_[best2]> 4)
	     Chi2_ndf_Ch1_2_[bst]= Chi2_ndf_Ch1_2_[bst]/(Nhits_Ch1_[best1]+Nhits_Ch2_[best2]-4);

	   //  cout<<" Chi2_ndf_Ch1_2 "<<Chi2_ndf_Ch1_2_[bst]<<endl;

	  // par_ab_Ch1_2_[1][bst] += (x1_sh + x2_sh)/2;
	  // par_ab_Ch1_2_[3][bst] += (y1_sh + y2_sh)/2;		  
	  // par_ab_Ch1_2_[0][bst] += ax12_sh + ax12_sh* par_ab_Ch1_2_[0][bst]* par_ab_Ch1_2_[0][bst];
	  // par_ab_Ch1_2_[2][bst] += ay12_sh + ay12_sh* par_ab_Ch1_2_[2][bst]* par_ab_Ch1_2_[2][bst];


	  ind_best_Ch1_2_[bst]= bst;		    
	  // Hist_Chi2_ndf_Ch12->Fill(Chi2_ndf_Ch1_2_[bst]);
	  // Hist_parAx_Ch12->Fill(par_ab_Ch1_2_[0][bst]);
	  // Hist_parBx_Ch12->Fill(par_ab_Ch1_2_[1][bst]);
	  // Hist_parAy_Ch12->Fill(par_ab_Ch1_2_[2][bst]);
	  // Hist_parBy_Ch12->Fill(par_ab_Ch1_2_[3][bst]);
	 

  }//< Nbest_Ch12_gl_
}//SegmentFit



void BmnMwpcHitFinderSRC::FillFitMatrix(Double_t** AA, Float_t* z, Float_t* sigm2_, Int_t* h_) {

  //out1<<" in FillFitMatrix "<<endl;

  // AA - matrix to be filledlayers)
  // sigm2 - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)
  // Float_t z2_[nPlanes];
  Float_t z2_[6] = {z[0] * z[0], z[1] * z[1], z[2] * z[2], z[3] * z[3], z[4] * z[4], z[5] * z[5]}; //cm


  AA[0][0] += 2 * z2_[0] * h_[0] / sigm2_[0] 
           +      z2_[2] * h_[2] / (2 * sigm2_[2]) 
           +      z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z2_[3] * h_[3] / sigm2_[3] 
           +      z2_[5] * h_[5] / (2 * sigm2_[5]) 
           +      z2_[4] * h_[4] / (2 * sigm2_[4]); //Ax

  AA[0][1] += 2 * z[0] * h_[0] / sigm2_[0] 
           +      z[2] * h_[2] / (2 * sigm2_[2]) 
           +      z[1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z[3] * h_[3] / sigm2_[3] 
           +      z[5] * h_[5] / (2 * sigm2_[5]) 
           +      z[4] * h_[4] / (2 * sigm2_[4]); //Bx

  AA[0][2] += sq3 * (z2_[2] * h_[2] / (2 * sigm2_[2]) 
           -         z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +         z2_[5] * h_[5] / (2 * sigm2_[5]) 
           -         z2_[4] * h_[4] / (2 * sigm2_[4])); //Ay

  AA[0][3] += sq3 * (z[2] * h_[2] / (2 * sigm2_[2]) 
           -         z[1] * h_[1] / (2 * sigm2_[1]) 
           +         z[5] * h_[5] / (2 * sigm2_[5]) 
		       -         z[4] * h_[4] / (2 * sigm2_[4])); //By

  AA[1][0] = AA[0][1];

  AA[1][1] +=   2 * h_[0] / sigm2_[0] 
           +  0.5 * h_[2] / sigm2_[2] + 0.5 * h_[1] / sigm2_[1] 
           +    2 * h_[3] / sigm2_[3] + 0.5 * h_[5] / sigm2_[5] 
           +  0.5 * h_[4] / sigm2_[4];

  AA[1][2] += sq3 * (z[2] * h_[2] / sigm2_[2] 
           - z[1] * h_[1] / sigm2_[1] 
           + z[5] * h_[5] / sigm2_[5] 
           - z[4] * h_[4] / sigm2_[4]) * 0.5;

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

  AA[2][3] += 3.0 * (z[2] * h_[2] / sigm2_[2] 
           +         z[1] * h_[1] / sigm2_[1] 
           +         z[5] * h_[5] / sigm2_[5] 
           +         z[4] * h_[4] / sigm2_[4])   * 0.5;

  AA[3][0] = AA[0][3];
  AA[3][1] = AA[1][3];
  AA[3][2] = AA[2][3];
  AA[3][3] += 3.0 * (0.5 * h_[2] / sigm2_[2] 
           +  0.5 *        h_[1] / sigm2_[1] 
           +  0.5 *        h_[5] / sigm2_[5] 
           +  0.5 *        h_[4] / sigm2_[4]);
}

vector<vector<double>> BmnMwpcHitFinderSRC::vkFillFitMatrix(Float_t* z, Float_t* sigm2_, Int_t* h_) {

  vector<vector<double>> AA(4, vector<double>(4));

  double z2_[6] = {z[0] * z[0], z[1] * z[1], z[2] * z[2], z[3] * z[3], z[4] * z[4], z[5] * z[5]}; //cm
  // cout<<"vkffm: h "<< h_[0] << " " << h_[1] << " " << h_[2] << " " << h_[3] << " " << h_[4] << " " << h_[5] <<endl;
  // cout<<"vkffm: z "<< z[0] << " " << z[1] << " " << z[2] << " " << z[3] << " " << z[4] << " " << z[5] <<endl;

  AA[0][0] += 2 * z2_[0] * h_[0] / sigm2_[0] 
           +      z2_[2] * h_[2] / (2 * sigm2_[2]) 
           +      z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z2_[3] * h_[3] / sigm2_[3] 
           +      z2_[5] * h_[5] / (2 * sigm2_[5]) 
           +      z2_[4] * h_[4] / (2 * sigm2_[4]); //Ax

  AA[0][1] += 2 * z[0] * h_[0] / sigm2_[0] 
           +      z[2] * h_[2] / (2 * sigm2_[2]) 
           +      z[1] * h_[1] / (2 * sigm2_[1]) 
           +  2 * z[3] * h_[3] / sigm2_[3] 
           +      z[5] * h_[5] / (2 * sigm2_[5]) 
           +      z[4] * h_[4] / (2 * sigm2_[4]); //Bx

  AA[0][2] += sq3 * (z2_[2] * h_[2] / (2 * sigm2_[2]) 
           -         z2_[1] * h_[1] / (2 * sigm2_[1]) 
           +         z2_[5] * h_[5] / (2 * sigm2_[5]) 
           -         z2_[4] * h_[4] / (2 * sigm2_[4])); //Ay

  AA[0][3] += sq3 * (z[2] * h_[2] / (2 * sigm2_[2]) 
           -         z[1] * h_[1] / (2 * sigm2_[1]) 
           +         z[5] * h_[5] / (2 * sigm2_[5]) 
		       -         z[4] * h_[4] / (2 * sigm2_[4])); //By

  AA[1][0] = AA[0][1];

  AA[1][1] +=   2 * h_[0] / sigm2_[0] 
           +  0.5 * h_[2] / sigm2_[2] + 0.5 * h_[1] / sigm2_[1] 
           +    2 * h_[3] / sigm2_[3] + 0.5 * h_[5] / sigm2_[5] 
           +  0.5 * h_[4] / sigm2_[4];

  AA[1][2] += sq3 * (z[2] * h_[2] / sigm2_[2] 
           - z[1] * h_[1] / sigm2_[1] 
           + z[5] * h_[5] / sigm2_[5] 
           - z[4] * h_[4] / sigm2_[4]) * 0.5;

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

  AA[2][3] += 3.0 * (z[2] * h_[2] / sigm2_[2] 
           +         z[1] * h_[1] / sigm2_[1] 
           +         z[5] * h_[5] / sigm2_[5] 
           +         z[4] * h_[4] / sigm2_[4])   * 0.5;

  AA[3][0] = AA[0][3];
  AA[3][1] = AA[1][3];
  AA[3][2] = AA[2][3];
  AA[3][3] += 3.0 * (0.5 * h_[2] / sigm2_[2] 
           +  0.5 *        h_[1] / sigm2_[1] 
           +  0.5 *        h_[5] / sigm2_[5] 
           +  0.5 *        h_[4] / sigm2_[4]);
  return AA;
}


void BmnMwpcHitFinderSRC::FillFreeCoefVector(Double_t* F, Float_t** XVU_Ch, Int_t ise, Float_t* z, Float_t* sigm2_, Int_t* h_, Int_t nPlanes) {
  // F - vector to be filled
  // XVU_Ch - coordinates of segment in chamber (Is it correct definition?)
  // segIdx - index of current segment
  // z - local z-positions of planes(layers)
  // sigm2_ - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)

  F[0] += 
    2 * XVU_Ch[0][ise] * z[0] * h_[0] / sigm2_[0] 
    + 
    XVU_Ch[1][ise]  * z[1] * h_[1] / sigm2_[1] 
    + 
    XVU_Ch[2][ise] * z[2] * h_[2] / sigm2_[2] 
    + 
    2 * XVU_Ch[3][ise] * z[3] * h_[3] / sigm2_[3] 
    + 
    XVU_Ch[4][ise] * z[4] * h_[4] / sigm2_[4] 
    + 
    XVU_Ch[5][ise] * z[5] * h_[5] / sigm2_[5];

  F[1] += 2 * XVU_Ch[0][ise] * h_[0] / sigm2_[0] + XVU_Ch[1][ise] * h_[1] / sigm2_[1] + XVU_Ch[2][ise] * h_[2] / sigm2_[2] + 2 * XVU_Ch[3][ise] * h_[3] / sigm2_[3] + XVU_Ch[4][ise] * h_[4] / sigm2_[4] + XVU_Ch[5][ise] * h_[5] / sigm2_[5];

  F[2] += sq3*(-XVU_Ch[1][ise] * z[1] * h_[1] / sigm2_[1] + XVU_Ch[2][ise] * z[2] * h_[2] / sigm2_[2] - XVU_Ch[4][ise] * z[4] * h_[4] / sigm2_[4] + XVU_Ch[5][ise] * z[5] * h_[5] / sigm2_[5]);

  F[3] +=  sq3*(-XVU_Ch[1][ise] * h_[1] / sigm2_[1] + XVU_Ch[2][ise] * h_[2] / sigm2_[2] - XVU_Ch[4][ise] * h_[4] / sigm2_[4] + XVU_Ch[5][ise] * h_[5] / sigm2_[5]);



}


void BmnMwpcHitFinderSRC::FillFreeCoefVectorXUV(Double_t* F, Float_t* XVU_Ch,  Float_t* z, Float_t* sigm2_, Int_t* h_) {
  // F - vector to be filled
  // XVU_Ch - coordinates of segment in chamber (Is it correct definition?)
  // segIdx - index of current segment
  // z - local z-positions of planes(layers)
  // sigm2_ - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)

  F[0] += 
    2 * XVU_Ch[0]  * z[0] * h_[0] / sigm2_[0] 
    + 
    XVU_Ch[1]   * z[1] * h_[1] / sigm2_[1] 
    + 
    XVU_Ch[2]  * z[2] * h_[2] / sigm2_[2] 
    + 
    2 * XVU_Ch[3]  * z[3] * h_[3] / sigm2_[3] 
    + 
    XVU_Ch[4]  * z[4] * h_[4] / sigm2_[4] 
    + 
    XVU_Ch[5]  * z[5] * h_[5] / sigm2_[5];

  F[1] += 2 * XVU_Ch[0]  * h_[0] / sigm2_[0] + XVU_Ch[1]  * h_[1] / sigm2_[1] + XVU_Ch[2]  * h_[2] / sigm2_[2] + 2 * XVU_Ch[3]  * h_[3] / sigm2_[3] + XVU_Ch[4]  * h_[4] / sigm2_[4] + XVU_Ch[5]  * h_[5] / sigm2_[5];
  F[2] += (-XVU_Ch[1]  * z[1] * h_[1] / sigm2_[1] + XVU_Ch[2]  * z[2] * h_[2] / sigm2_[2] - XVU_Ch[4]  * z[4] * h_[4] / sigm2_[4] + XVU_Ch[5]  * z[5] * h_[5] / sigm2_[5]);
  F[3] +=  (-XVU_Ch[1]  * h_[1] / sigm2_[1] + XVU_Ch[2]  * h_[2] / sigm2_[2] - XVU_Ch[4]  * h_[4] / sigm2_[4] + XVU_Ch[5]  * h_[5] / sigm2_[5]);

  F[2]=F[2]*sq3;
  F[3]=F[3]*sq3;


}


void BmnMwpcHitFinderSRC::InverseMatrix(Double_t** AA, Double_t** bb) {
  // Gaussian algorithm for 4x4 matrix inversion 


    Double_t factor;
    Double_t temp[4];

    // Set b to I
    for (Int_t i1 = 0; i1 < 4; i1++)
        for (Int_t j1 = 0; j1 < 4; j1++)
            if (i1 == j1) bb[i1][j1] = 1.0;
            else bb[i1][j1] = 0.0;

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
    //end inverse

   
  
}
 */

void BmnMwpcHitFinderSRC::Finish() {
  cout<<" delete "<<endl;
 
//  TCanvas* c1 = new TCanvas("c1","c1",600,600);
//   hNp_best_ch1->Draw();

//  TCanvas* c4 = new TCanvas("c4","c4",600,600);
//    hNp_best_ch2->Draw();

//  TCanvas* c2 = new TCanvas("c2","c2",600,600);
//   hNbest_Ch2->Draw();

//  TCanvas* c3 = new TCanvas("c3","c3",600,600);
//   hNbest_Ch1->Draw();

//  TCanvas* c5 = new TCanvas("c5","c5",600,600);
//    hChi2_ch1_2->Draw();
 
  // fList.Draw();

 
    delete fMwpcGeometry;

    // delete 1d arrays:

    
   
    delete [] shift1_2;
    delete [] Chi2_match;
    delete [] Chi2_ndf_Ch1_2;
    delete [] ind_best_Ch1_2;
    delete [] sigm2;
    delete [] ipl;
   
   
   
    delete [] z2;

    //  3d-arrays:
  
    for(Int_t iCh=0; iCh< kNChambers; iCh++){
      for(Int_t iWire=0; iWire< kNWires; iWire++){
	delete [] wire_Ch[iCh][iWire];	// delete [] wire_Ch1[iWire];
	delete [] xuv_Ch[iCh][iWire];//  delete [] wire_Ch2[iWire];	
      }
      for(Int_t iPlane=0; iPlane<kNPlanes; iPlane++){
	delete [] Wires_Ch[iCh][iPlane];
	delete [] clust_Ch[iCh][iPlane];
	delete [] XVU_Ch[iCh][iPlane];
      }

      for(Int_t ii=0; ii<4; ii++){
	delete [] par_ab_Ch[iCh][ii];
      }
 
      // delete 2d arrays

      delete [] kPln[iCh];
      delete [] iw[iCh];
      delete [] iw_Ch[iCh];
      delete [] wire_Ch[iCh];
      delete [] xuv_Ch[iCh];
      delete [] Wires_Ch[iCh];
      delete [] clust_Ch[iCh];
      delete [] XVU_Ch[iCh];
      delete [] Nhits_Ch[iCh];
      delete [] shift[iCh]; //  delete [] shift1;
      delete [] z_gl[iCh];
      delete [] kZ_loc[iCh];// delete [] kZ1_loc;
      delete [] ind_best_Ch[iCh];
      delete [] best_Ch_gl[iCh];
      delete [] Chi2_ndf_Ch[iCh];
      delete [] Chi2_ndf_best_Ch[iCh];
      delete [] par_ab_Ch[iCh];
      delete [] XVU[iCh];// delete [] XVU1;
      delete [] XVU_cl[iCh];// delete [] XVU_cl1;
   
         
    }
    // delete 1d arrays

    for(Int_t ii=0; ii<4; ii++){
      // delete [] par_ab_Ch1[ii];
      delete [] par_ab_Ch1_2[ii];
      delete [] matrA[ii];
      delete [] matrb[ii];
    }

    delete [] wire_Ch;
    delete [] xuv_Ch;
    delete [] kPln;
    delete [] Wires_Ch;
    delete [] clust_Ch;
    delete [] XVU_Ch;
    delete [] Nhits_Ch;
    delete [] Nseg_Ch;
    delete [] Nbest_Ch;
    delete [] shift;
    delete [] kZ_loc;
    delete [] z_gl;
    delete [] ZCh;
    delete [] kZmid;
    delete [] ind_best_Ch;
    delete [] best_Ch_gl;
    delete [] Chi2_ndf_Ch;
    delete [] Chi2_ndf_best_Ch;
    delete [] par_ab_Ch;
    delete [] XVU;
    delete [] XVU_cl;

    
    
    cout << "Work time of the MWPC hit finder: " << workTime << " s" << endl;
}

ClassImp(BmnMwpcHitFinderSRC)


