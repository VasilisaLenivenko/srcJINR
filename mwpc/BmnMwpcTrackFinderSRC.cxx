// @(#)bmnroot/mwpc:$Id$


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcTrackFinder                                                         //
//                                                                            // 
//                                                                            //
//Implementation of an algorithm developed by                                 // 
// Vasilisa Lenivenko  and Vladimir Palchik                                   //
// to the BmnRoot software                                                    //
//                                                                            //
// The algorithm serves for searching for track segments                      //
// in the MWPC of the BM@N experiment                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <Rtypes.h>

#include "BmnMwpcTrackFinderSRC.h"

static Float_t workTime = 0.0;

BmnMwpcTrackFinderSRC::BmnMwpcTrackFinderSRC(Bool_t isExp, Int_t runP) :
  fEventNo(0),
  expData(isExp) {
  fInputBranchName = "BmnMwpcSegment";
  fOutputBranchName = "BmnMwpcTrack";
  fRunPeriod = runP;
  fMwpcGeo = new BmnMwpcGeometrySRC(fRunPeriod);
  kBig = 100;
  kNumPairs = 1;
  kCh_min = 2; 
  kCh_max = 4;
  if ( fRunPeriod == 7 ) {
    kNumPairs = 2; 
    kCh_min = 0;
  }   
}

BmnMwpcTrackFinderSRC::~BmnMwpcTrackFinderSRC() {
}

void BmnMwpcTrackFinderSRC::Exec(Option_t* opt) {
  if (!IsActive()) return;
  clock_t tStart = clock();
  PrepareArraysToProcessEvent();
  if (fDebug) cout << "\n======================== MWPC track finder exec started ===================\n" << endl;
  if (fDebug) cout << "Event number: " << fEventNo++ << endl;

    // fBmnMwpcTracksArray->Clear();
   
  for (Int_t iSegment = 0; iSegment < fBmnMwpcSegmentsArray->GetEntries(); iSegment++) {
    BmnTrack* segment = (BmnTrack*) fBmnMwpcSegmentsArray->At(iSegment);

    Int_t iCh;
    Double_t Z = segment->GetParamFirst()->GetZ();
    Int_t ise = segment->GetFlag();//iSegmentID

    if ( ZCh[0] == Z) iCh = 0;
    if ( ZCh[1] == Z) iCh = 1;
    if ( ZCh[2] == Z) iCh = 2;
    if ( ZCh[3] == Z) iCh = 3;

    Nhits_Ch[iCh][ise] = segment->GetNHits();
    Chi2_ndf_Ch[iCh][ise] = segment->GetChi2();
    par_ab_Ch[iCh][1][ise] = segment->GetParamFirst()->GetX();
    par_ab_Ch[iCh][0][ise] = segment->GetParamFirst()->GetTx();
    par_ab_Ch[iCh][3][ise] = segment->GetParamFirst()->GetY();
    par_ab_Ch[iCh][2][ise] = segment->GetParamFirst()->GetTy();   
   
    XVU_Ch[iCh][0][ise] =  segment->GetParamLast()->GetX();
    XVU_Ch[iCh][1][ise] =  segment->GetParamLast()->GetY();
    XVU_Ch[iCh][2][ise] =  segment->GetParamLast()->GetZ();
    XVU_Ch[iCh][3][ise] =  segment->GetParamLast()->GetTx();
    XVU_Ch[iCh][4][ise] =  segment->GetParamLast()->GetTy();
    XVU_Ch[iCh][5][ise] =  segment->GetParamLast()->GetQp();
    Nbest_Ch[iCh]++;
  }//iSegment

    //  cout<<" iSegment ---------------------------------------------"<<endl;
  for (Int_t iChamber = kCh_min; iChamber < kCh_max; iChamber++) {

   //  cout<<" Nbest_Ch["<<iChamber<<"] = "<<Nbest_Ch[iChamber]<<endl;
   for (Int_t ise = 0; ise < Nbest_Ch[iChamber]; ise++) {
      // cout<<" iCh "<<iChamber<<" Nhits_Ch[iCh]["<<ise<<"]= "<<Nhits_Ch[iChamber][ise]<<endl;
      hChi2_ndf_Ch.at(iChamber) -> Fill(Chi2_ndf_Ch[iChamber][ise]);
      //	 cout<<" Ch= "<<iChamber<<" ise "<<ise<<" Chi2 "<<Chi2_ndf_Ch[iChamber][ise]<<" Ax "<<par_ab_Ch[iChamber][0][ise]<<" bx "<<par_ab_Ch[iChamber][1][ise]<<" Ay "<<par_ab_Ch[iChamber][2][ise]<<" by "<<par_ab_Ch[iChamber][3][ise]<<endl;
     for (Int_t i = 0; i < 6; i++){
       //   cout<<" Ch= "<<iChamber<<" i "<<i<<" XVU_Ch "<<XVU_Ch[iChamber][i][ise]<<endl;
     }
      // cout<<endl;
    }
  }//iChamber

  for (Int_t iChamber = kCh_min; iChamber < kCh_max; iChamber++) {
    if ( Nbest_Ch[iChamber] > 0){
      SegmentParamAlignment(iChamber, Nbest_Ch, par_ab_Ch, shift);
      for (Int_t ise = 0; ise < Nbest_Ch[iChamber]; ise++) {
        hpar_Ax_Ch.at(iChamber) -> Fill( par_ab_Ch[iChamber][0][ise]);
        hpar_Bx_Ch.at(iChamber) -> Fill( par_ab_Ch[iChamber][1][ise]);
        hpar_Ay_Ch.at(iChamber) -> Fill( par_ab_Ch[iChamber][2][ise]);
        hpar_By_Ch.at(iChamber) -> Fill( par_ab_Ch[iChamber][3][ise]);
      }
    }
  }

  //  cout<<" N0= "<<Nbest_Ch[0]<<" N1= "<<Nbest_Ch[1]<<" N2= "<<Nbest_Ch[2]<<" N3= "<<Nbest_Ch[3]<<endl;
  if ( Nbest_Ch[0] > 0 && Nbest_Ch[1] > 0)  SegmentMatching( 0, Nbest_Ch, par_ab_Ch, kZmid,   ind_best_Ch, Nbest_pair,  Chi2_match_pair, XVU_Ch, Nhits_Ch);
  if ( Nbest_Ch[2] > 0 && Nbest_Ch[3] > 0)  SegmentMatching( 2, Nbest_Ch, par_ab_Ch, kZmid,   ind_best_Ch, Nbest_pair,  Chi2_match_pair, XVU_Ch, Nhits_Ch);


  for (Int_t p = 0; p < kNumPairs; p++) {
    for (Int_t se = 0; se < 5; se++) {
      if (Chi2_match_pair[p][se] != 999.) {
        hChi2_match_pair.at(p) -> Fill( Chi2_match_pair[p][se]);
        //   cout<<" pair "<<p<<" se "<<se<<" Chi2_match "<<Chi2_match_pair[p][se]<<endl;
      }
    }
  }
     
     // ----spatial track---
   
  Int_t First_Chamber = 2; 
  if ( fRunPeriod == 7 ) First_Chamber = 0; 
  //  if (fDebug) cout<<" Nbest_pair[0] "<<Nbest_pair[0]<<endl;
  if ( Nbest_pair[0] > 0) SegmentFit(First_Chamber, z_gl, sigm2,  Nbest_pair, ind_best_Ch,  par_ab_pair, Chi2_ndf_pair,  XVU_Ch,  ind_best_pair,  Nhits_Ch);
  
  // if (fDebug) cout<<" Nbest_pair[1] "<<Nbest_pair[1]<<endl;
  if ( Nbest_pair[1] > 0) SegmentFit(2, z_gl, sigm2,  Nbest_pair, ind_best_Ch,  par_ab_pair, Chi2_ndf_pair,  XVU_Ch,  ind_best_pair,  Nhits_Ch);

  Double_t theta, phi;
 
  for (Int_t iPair = 0; iPair < kNumPairs; iPair++) {
    if ( Nbest_pair[iPair] > 0){
      SegmentParamAlignment(iPair, Nbest_pair,  par_ab_pair, shift_pair);
      for (Int_t itr = 0; itr < Nbest_pair[iPair]; itr++) {
	phi = TMath::ATan2(par_ab_pair[iPair][2][itr],par_ab_pair[iPair][0][itr]); // phi = arctan(tgy/tgx)
        theta = TMath::ATan2(par_ab_pair[iPair][0][itr], TMath::Cos(phi));// theta = arctan(tgx/cos(phi))

        hpar_Ax_pair.at(iPair)   -> Fill(TMath::RadToDeg()* par_ab_pair[iPair][0][itr]);
        hpar_Bx_pair.at(iPair)   -> Fill( par_ab_pair[iPair][1][itr]);
        hpar_Ay_pair.at(iPair)   -> Fill(TMath::RadToDeg()* par_ab_pair[iPair][2][itr]);
        hpar_By_pair.at(iPair)   -> Fill( par_ab_pair[iPair][3][itr]);
        hpar_theta_pair.at(iPair) -> Fill(TMath::RadToDeg()*theta);
	hpar_phi_pair.at(iPair) -> Fill(TMath::RadToDeg()*phi);
	//	cout<<" theta "<<theta<<endl;
        //cout<<" iPair "<<iPair<<" itr "<<itr<<" Chi2_ndf_pair "<<Chi2_ndf_pair[iPair][itr]<<endl;
        if ( fRunPeriod == 7 && iPair == 1 ) {
        Float_t X_par_to_pole = par_ab_pair[iPair][0][itr]*( 0 - kZ_midle_pair[iPair]) + par_ab_pair[iPair][1][itr];
        Float_t Y_par_to_pole = par_ab_pair[iPair][2][itr]*( 0 - kZ_midle_pair[iPair]) + par_ab_pair[iPair][3][itr];
        //cout<<" X_par_to_pole "<<X_par_to_pole<<" Y_par_to_pole "<< Y_par_to_pole<<endl;
	}
	 

	  Float_t X_par_to_target = par_ab_pair[iPair][0][itr]*( kZ_target - kZ_midle_pair[iPair]) + par_ab_pair[iPair][1][itr];
	  Float_t Y_par_to_target = par_ab_pair[iPair][2][itr]*( kZ_target - kZ_midle_pair[iPair]) + par_ab_pair[iPair][3][itr];
	  if ( fRunPeriod == 7 && iPair == 0 ) {
	     hX_in_target ->   Fill(X_par_to_target);
	     hY_in_target ->   Fill(Y_par_to_target);
	     hY_X_in_target -> Fill(X_par_to_target,Y_par_to_target);
	  }

	  hAx_bx_in_target -> Fill(X_par_to_target, TMath::RadToDeg()*par_ab_pair[iPair][0][itr]);
	  hAy_by_in_target -> Fill(Y_par_to_target, TMath::RadToDeg()*par_ab_pair[iPair][2][itr]);
	  
  
        BmnTrack *Tr = new ((*fBmnMwpcTracksArray)[fBmnMwpcTracksArray->GetEntriesFast()]) BmnTrack();
        Tr -> SetChi2(Chi2_ndf_pair[iPair][itr]);
        FairTrackParam TrParams;
        TrParams.SetPosition(TVector3(par_ab_pair[iPair][1][itr], par_ab_pair[iPair][3][itr],kZ_midle_pair[iPair]));
        //TrParams.SetPosition(TVector3(X_par_to_pole, Y_par_to_pole, 0);
        TrParams.SetTx(par_ab_pair[iPair][0][itr]);
        TrParams.SetTy(par_ab_pair[iPair][2][itr]);
        Tr -> SetParamFirst(TrParams);
      }//Nbest_pair[iPair]
    }//> 0
  }//iPair

  if ( Nbest_pair[0] > 0 && Nbest_pair[1] > 0 ) PairMatching(Nbest_pair, par_ab_pair, kZmid);

  if (fDebug) cout << "\n======================== MWPC track finder exec finished ==================" << endl;
  clock_t tFinish = clock();
  workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
}//Exec

//----------------------------------------------------------------------------------


void BmnMwpcTrackFinderSRC::SegmentParamAlignment(Int_t chNum, Int_t *Nbest,  Double_t ***par_ab, Float_t **shiftt ){
  //local parameters to Global parameters

  for (Int_t iBest = 0; iBest < Nbest[chNum]; iBest++) {

    // cout<<"before Alignment: iBest "<<iBest<<" Ax "<< par_ab[chNum][0][ind_best[chNum][iBest]]<<" bx "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" Ay "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" by "<< par_ab[chNum][3][ind_best[chNum][iBest]]<<endl;

    //                                     ax          alpha                                      ax^2    
    par_ab[chNum][0][iBest] += shiftt[chNum][0] +  shiftt[chNum][0]* par_ab[chNum][0][iBest]* par_ab[chNum][0][iBest];
    par_ab[chNum][2][iBest] += shiftt[chNum][2] +  shiftt[chNum][2]* par_ab[chNum][2][iBest]* par_ab[chNum][2][iBest];
    par_ab[chNum][1][iBest] += shiftt[chNum][1];
    par_ab[chNum][3][iBest] += shiftt[chNum][3];

    //cout<<"after Alignment: iBest "<<iBest<<" Ax "<< par_ab[chNum][0][ind_best[chNum][iBest]]<<" bx "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" Ay "<< par_ab[chNum][1][ind_best[chNum][iBest]]<<" by "<< par_ab[chNum][3][ind_best[chNum][iBest]]<<endl;

  }//iBest 
}//SegmentParamAlignment



void BmnMwpcTrackFinderSRC::SegmentMatching( Int_t first_Ch, Int_t *Nbest, Double_t ***par_ab, Float_t *Zmid, Int_t **best_Ch, Int_t *Nbest_pair_, Double_t **Chi2_match_, Float_t ***XVU_Ch_, Int_t **Nhits_){

  //  cout<<" SegmentMatching( "<<" first_Ch "<<first_Ch<<" Nbest "<<Nbest[first_Ch]<<" par_ab_ "<<par_ab[first_Ch][0][0]<<" "<<par_ab[first_Ch][1][0]<<" "<<par_ab[first_Ch][2][0]<<" "<<par_ab[first_Ch][3][0]<<" Zmid2 "<< Zmid[2]<<" Zmid3 "<<Zmid[3]<<endl;

  Float_t sig_dx= 3.6; //0.8;//22;- z=0 //0.85;- Zmid
  Float_t sig_dy= 3.6; //0.7;//18;- z=0 //0.76;
  Float_t sig_dax= 0.055; //0.04; //0.063;
  Float_t sig_day= 0.055; //0.04; //0.045;
	    
  Float_t min_Chi2m = 100; // 40; //100; //400
  Float_t min_distX = 99;
  Float_t min_distY = 99;
  Float_t dAx12 = 0;
  Float_t dAy12 = 0;
  Int_t  best1 = -1;
  Int_t  best2 = -1;
  Float_t Min_distX[kmaxPairs];
  Float_t Min_distY[kmaxPairs];

  for (Int_t i = 0; i < kmaxPairs ; i++) {
    //   best1[i] = -1;
    // best2[i] = -1;
    Min_distX[i] = -1;
    Min_distY[i] = -1;
  }
  
  Float_t DAx12 = 0;
  Float_t DAy12 = 0;

  Int_t Pairr = 0;
  if ( fRunPeriod == 7 && first_Ch == 2) Pairr = 1;
  Int_t Secon_Ch = first_Ch+1;

  if (Nbest_Ch[first_Ch] > 0 && Nbest_Ch[Secon_Ch] > 0){ 
   	   	    	 	  
    //   cout<<" Nbest[ "<<first_Ch<<"] = "<<Nbest_Ch[first_Ch]<<endl;
    //  cout<<" Nbest[ "<<Secon_Ch<<"] = "<<Nbest_Ch[Secon_Ch]<<endl;

    for (Int_t bst1 = 0; bst1 < Nbest[first_Ch]; bst1++) {	  
      if ( bst1 ==  best_Ch[first_Ch][Nbest_pair_[Pairr]] -1) continue;
      //  cout<<" bst1 "<<bst1<<"  best1 "<< best1<<endl;
      //ch1                                     zloc0 -z_i
      Float_t x1mid = par_ab[first_Ch][0][bst1] *( 0 - kZmid[first_Ch]) + par_ab[first_Ch][1][bst1] ;
      Float_t y1mid = par_ab[first_Ch][2][bst1] *( 0 - kZmid[first_Ch]) + par_ab[first_Ch][3][bst1] ;
      
      //   cout<<" bst1 " <<bst1<<" x1mid "<<x1mid<<" y1mid "<<y1mid<<endl;
                 
      for (Int_t bst2 = 0; bst2 < Nbest[Secon_Ch]; bst2++) {


	if ( bst2 ==  best2) continue;
	//	cout<<" bst2 "<<bst2<<"  best2 "<<best2<<endl;

	//ch2       
	Float_t x2mid =  par_ab[Secon_Ch][0][bst2] *( 0 - kZmid[Secon_Ch])  + par_ab_Ch[Secon_Ch][1][bst2] ;
	Float_t y2mid =  par_ab[Secon_Ch][2][bst2] *( 0 - kZmid[Secon_Ch])  + par_ab_Ch[Secon_Ch][3][bst2] ;
	
	//	cout<<" bst2 " <<bst2<<" x2mid "<<x2mid<<" y2mid "<<y2mid<<endl;
		  		  
	dAx12 = par_ab[first_Ch][0][bst1] - par_ab[Secon_Ch][0][bst2];
	dAy12 = par_ab[first_Ch][2][bst1] - par_ab[Secon_Ch][2][bst2];
	min_distX = x1mid - x2mid; //min
	min_distY = y1mid - y2mid; //min

	Float_t Chi2_m = ( min_distX*min_distX/(sig_dx*sig_dx)  +  min_distY*min_distY/(sig_dy*sig_dy)
			   +  dAx12*dAx12 /(sig_dax*sig_dax) +  dAy12*dAy12 /(sig_day*sig_day) );

	//	cout<<" Chi2_m "<<Chi2_m<<" min_Chi2m "<<min_Chi2m<<" Chi2_match_ "<<Chi2_match_[first_Chr][0]<<endl;

	if (Chi2_m < min_Chi2m && Nbest_pair_[Pairr] < kmaxPairs ){
	  min_Chi2m = Chi2_m; //min
	  best1 = bst1;
	  best2 = bst2;
	  best_Ch[first_Ch][Nbest_pair_[Pairr]] = bst1;
	  best_Ch[Secon_Ch][Nbest_pair_[Pairr]] = bst2;
	  Chi2_match_[Pairr][Nbest_pair_[Pairr]] =  Chi2_m;
	  Min_distX[Nbest_pair_[Pairr]]= min_distX;
	  Min_distY[Nbest_pair_[Pairr]]= min_distY;
	  Nbest_pair_[Pairr]++;

	  DAx12=dAx12;
	  DAy12=dAy12;	  
	  // cout<<"  best1 "<< best1<<"  best2 "<<best2<<endl;
	  /*
	  for (int j = kmaxPairs -1; j > 0; j--) {
	    best_Ch[first_Ch][j] =  best_Ch[first_Ch][j-1];
	    best_Ch[Secon_Ch][j] =  best_Ch[Secon_Ch][j-1];
	    Chi2_match_[Pairr][j] =  Chi2_match_[Pairr][j-1];
	  }	  
	  best_Ch[first_Ch][0] = bst1;
	  best_Ch[Secon_Ch][0] = bst2;
	  Chi2_match_[Pairr][0] =  Chi2_m;
	  */

	  //  cout<<" Nbest_pair_["<<Pairr<<"] "<< Nbest_pair_[Pairr]
	  //     <<" best_Ch[Ch"<<first_Ch<<"][Nbest_pair"<<Nbest_pair_[Pairr]-1<<"] "<< best_Ch[first_Ch][Nbest_pair_[Pairr]-1]<<" bst1 "<<bst1<<" best1 "<<best1
	  //     <<" best_Ch[Ch"<<Secon_Ch<<"][Nbest_pair"<<Nbest_pair_[Pairr]-1<<"] "<< best_Ch[Secon_Ch][Nbest_pair_[Pairr]-1]<<" bst2 "<<bst2<<" best2 "<<best2
	  //     <<" Chi2_match_["<<Pairr<<"]["<<Nbest_pair_[Pairr]-1<<"] "<< Chi2_match_[Pairr][Nbest_pair_[Pairr]-1]<<endl;
	  
	}//if (Chi2_m < min_Chi2m
      }//bst2++     
    }//bst1++
   
      for (Int_t ii = 0; ii < Nbest_pair_[Pairr]; ii++) {
	/*
	if (  Chi2_match_[Pairr][ii+1] <  Chi2_match_[Pairr][ii]){
	  best_Ch[first_Ch][ii]  = bst1;
	  best_Ch[Secon_Ch][ii]  = bst2;
	  Chi2_match_[Pairr][ii] = Chi2_m;
	}
	*/
	FillEfficiency( first_Ch, XVU_Ch, Nhits_Ch,  kMinHits, best_Ch[first_Ch][ii], Min_distX[ii], Min_distY[ii]);
	FillEfficiency( Secon_Ch, XVU_Ch, Nhits_Ch,  kMinHits, best_Ch[Secon_Ch][ii], Min_distX[ii], Min_distY[ii]);
      } 
  }// N >0

}// SegmentMatching


void BmnMwpcTrackFinderSRC::PairMatching( Int_t *Nbest_p, Double_t ***par_ab,  Float_t *Zmid){

  Double_t dist,  dist2, cut2v;
  cut2v = 100.;
  Int_t Npair_dist = 0, Npair1;
  	   	    	 	  
  for (Int_t pair1 = 0; pair1 < Nbest_p[0]; pair1++) {	  
    Npair1 = 0;
    //pair 0                                 zloc0 -z_i
    Float_t x1 = par_ab[0][0][pair1] *( kZ_target - kZmid[0] ) + par_ab[0][1][pair1] ;
    Float_t y1 = par_ab[0][2][pair1] *( kZ_target - kZmid[0] ) + par_ab[0][3][pair1] ;
                      
    for (Int_t pair2 = 0; pair2 < Nbest_p[1]; pair2++) {

      //pair 1      
      Float_t x2 =  par_ab[1][0][pair2] *( kZ_target - kZmid[1] )  + par_ab_Ch[1][1][pair2] ;
      Float_t y2 =  par_ab[1][2][pair2] *( kZ_target - kZmid[1] )  + par_ab_Ch[1][3][pair2] ;
		  		  
      Float_t dAx12 = par_ab[0][0][pair1] - par_ab[1][0][pair2];
      Float_t dAy12 = par_ab[0][2][pair1] - par_ab[1][2][pair2];
      Float_t dX = x1 - x2; 
      Float_t dY = y1 - y2; 

      dist2 = dX*dX + dY*dY;
      if ( dist2 < cut2v) { 
	dist = sqrt(dist2);

      	hdist_target -> Fill(dist);
       	hdX_target ->Fill(dX);
       	hdY_target ->Fill(dY);
       	Npair_dist++;
       	Npair1++;
     
      }
      //  cout<<"PairMatching: dX "<<dX<<" dY "<<dY<<endl;
      // cout<<endl;

    }//pair2
    if ( Npair1 > 1 ) break;
  }//pair1

   hNsecondaries -> Fill(Npair_dist);

}//PairMatching


void BmnMwpcTrackFinderSRC::SegmentFit(Int_t First_Ch, Float_t **z_gl_, Float_t *sigm2_,Int_t *Nbest_pair_, Int_t ** ind_best_Ch_, Double_t ***par_ab_pair_,  Double_t **Chi2_ndf_pair_, Float_t ***XVU_Ch_,Int_t **ind_best_pair_,Int_t **Nhits_Ch_){
  int chiEl = 0;
  if (First_Ch == 2) chiEl = 1;
  if (Nbest_pair_[chiEl] >= 10) {printf("!!! ERROR: Nbest_pair_[%d] > 10\n", chiEl); return;}
   /*
  printf("PS input: %d %8.4f %8.4f %3d %3d %8.4f %8.4f %8.4f %3d %3d\n", 
                First_Ch,                 z_gl[0][0],                sigm2_[0],                Nbest_pair_[0],
                ind_best_Ch_[0][0],       par_ab_pair_[chiEl][0][0], Chi2_ndf_pair_[chiEl][0], XVU_Ch_[0][0][0], 
                ind_best_pair_[chiEl][0], Nhits_Ch_[0][0]);
   */
  Int_t Pair1 = 0;
  if (fRunPeriod == 7 && First_Ch == 2) Pair1 = 1;

   for (Int_t bst = 0; bst < Nbest_pair_[Pair1]; bst++) { // for (Int_t bst = 0; bst < Nbest_Ch12_gl_; bst++) {

     //  cout<<"------- bst "<<bst<<endl;

   Int_t fir = First_Ch;
   Int_t sec = First_Ch+1;

   Int_t  best1 = ind_best_Ch_[fir][bst];  //  Int_t bst1 = best_Ch1_gl_[bst];
   Int_t  best2 = ind_best_Ch_[sec][bst]; //  Int_t bst2 = best_Ch2_gl_[bst];
	
   // Int_t best1 = ind_best_Ch1_[bst1]; // it's segment index !
   // Int_t best2 = ind_best_Ch2_[bst2]; // it's segment index!

   // cout<<" SegmentFit: best1 "<<best1<<" best2 "<<best2<<" bst "<<bst<<" fir "<<fir<<" sec "<<sec<<endl;
    
   int h1[6] = {0,0,0,0,0,0};
   int h2[6] = {0,0,0,0,0,0};
    for(Int_t i = 0; i<6;i++){
      if ( XVU_Ch_[fir][i][best1]  > -999.) h1[i] = 1; //  if(Wires_Ch1_[i][best1]> -1) h1[i] = 1;
      if ( XVU_Ch_[sec][i][best2]  > -999.) h2[i] = 1; //   if(Wires_Ch2_[i][best2]> -1) h2[i] = 1;
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
	//	cout<<" Amatr "<<Amatr[im][ii]<<" bmatr "<< bmatr[im][ii]<<endl;
      }
    }
    for (Int_t i1 = 0; i1 < 6; i1++){
      //  cout<<" fir "<<fir<<" z_gl_ "<<z_gl_[fir][i1]<<" h1 "<<h1[i1]<<endl;
      // cout<<" sec "<<fir<<" z_gl_ "<<z_gl_[sec][i1]<<" h1 "<<h2[i1]<<endl;
    }

    FillFitMatrix(fir, Amatr, z_gl_, sigm2_, h1); //FillFitMatrix(Amatr, z_gl1_, sigm2_, h1);//Ch_1
    FillFitMatrix(sec, Amatr, z_gl_, sigm2_, h2);//FillFitMatrix(Amatr, z_gl2_, sigm2_, h2);//Ch_2
   
    //    for(int ff = 0; ff < 4; ff++){
    //        for(int tt = 0; tt < 4; tt++){
    //	  cout<<"1 Amatr ["<<ff <<"] ["<<tt<<"] "<<Amatr[ff][tt]<<endl;
    //	}
    //    }
   
   Double_t matrF[4] = {0,0,0,0};//free coef 
   
 //  cout<<" before FillFreeCoefVector "<<endl;
   for(Int_t i = 0; i < 6;i++){
     //  cout<<" sigm2 "<< sigm2_[i]<<" h1 "<<h1[i]<<" XVU_Ch1_ "<<XVU_Ch_[fir][i][best1]<<" h2 "<<h2[i]<<" XVU_Ch2_ "<<XVU_Ch_[sec][i][best2]<<endl;
   }  
  
   FillFreeCoefVector(fir, matrF, XVU_Ch, best1, z_gl_ , sigm2_, h1);
   // cout<<" between "<<endl;
   for (Int_t i1 = 0; i1 < 4; i1++){
     //  cout<<"b matrF["<<i1<<"]= "<< matrF[i1]<<endl;
   }

   FillFreeCoefVector(sec, matrF, XVU_Ch, best2, z_gl_ , sigm2_, h2);

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

	  
	  for(Int_t i1 = 0 ; i1 < 4; i1++){
	    par_ab_pair_[Pair1][i1][bst] = 0;//  par_ab_Ch1_2_[i1][bst] = 0;
	    for(Int_t j1 = 0; j1 < 4; j1++){
	      par_ab_pair_[Pair1][i1][bst] += bmatr[i1][j1]*matrF[j1];//  par_ab_Ch1_2_[i1][bst] += bmatr[i1][j1]*matrF[j1];// par_ab_Ch23[i1][bst] += b[i1][j1]*F[j1];
	    //  cout<<" i1 "<<i1<<" bmatr "<<bmatr[i1][j1]<<" F "<<matrF[j1] <<endl;
	    }    
	  } // i1	
	  //	  cout<<endl;
	  // cout<<" Pair "<<Pair1<<" par [0] "<<par_ab_pair_[Pair1][0][bst]<<" par [1] "<<par_ab_pair_[Pair1][1][bst]<<" par [2] "<<par_ab_pair_[Pair1][2][bst]<<" par [3] "<<par_ab_pair_[Pair1][3][bst]<<endl;
	  //  cout<<endl;
	 
	  Float_t sigm_gl1[12];
	  Float_t sigm_gl2[12];
	  for(Int_t i1 = 0 ; i1 < 6; i1++){
	    sigm_gl1[i1]=0.0722;
	    sigm_gl2[i1]=0.0722;
	  }
	  
	  Float_t dx_[kNPlanes];
	  Chi2_ndf_pair_[Pair1][bst]=0;

	  for(Int_t i1 = 0 ; i1 < 6; i1++){
	    dx_[i1] = 0.;
	    
	    if ( XVU_Ch_[fir][i1][best1]  > -999.){// if(Wires_Ch1_[i1][best1]>-1){
	      if(i1==0 || i1==3) dx_[i1]=XVU_Ch_[fir][i1][best1]-par_ab_pair_[Pair1][0][bst]*z_gl_[fir][i1]-par_ab_pair_[Pair1][1][bst];
	      if(i1==2 || i1==5) dx_[i1]=XVU_Ch_[fir][i1][best1]-0.5*(par_ab_pair_[Pair1][0][bst]+sq3*par_ab_pair_[Pair1][2][bst])*z_gl_[fir][i1]-0.5*(par_ab_pair_[Pair1][1][bst]+sq3*par_ab_pair_[Pair1][3][bst]);
	      if(i1==1 || i1==4) dx_[i1]=XVU_Ch_[fir][i1][best1]-0.5*(par_ab_pair_[Pair1][0][bst]-sq3*par_ab_pair_[Pair1][2][bst])*z_gl_[fir][i1]-0.5*(par_ab_pair_[Pair1][1][bst]-sq3*par_ab_pair_[Pair1][3][bst]);
	      Chi2_ndf_pair_[Pair1][bst]= Chi2_ndf_pair_[Pair1][bst]+dx_[i1]*dx_[i1]/(sigm_gl1[i1]*sigm_gl1[i1]);

	      //   cout<<"best1 "<<best1 <<" i1 "<<i1<<" dx_ "<<dx_[i1]<<" XVU_Ch1 "<<XVU_Ch_[fir][i1][best1]<<" Chi2_ndf_Ch1_2 "<<Chi2_ndf_pair_[Pair1][bst]<<" z_gl1 "<<z_gl_[fir][i1]<<endl;

	    }// if( Wires_Ch1[i1][best2]>-1){
	  }
	  //  cout<<endl;

	  for(Int_t i2 = 0 ; i2 < 6; i2++){
	   
	    if ( XVU_Ch_[sec][i2][best1]  > -999.){// if(Wires_Ch2_[i2][best2]>-1){
	      if(i2==0 || i2==3) dx_[i2]=XVU_Ch_[sec][i2][best2]-par_ab_pair_[Pair1][0][bst]*z_gl_[sec][i2]-par_ab_pair_[Pair1][1][bst];
	      if(i2==2 || i2==5) dx_[i2]=XVU_Ch_[sec][i2][best2]-0.5*(par_ab_pair_[Pair1][0][bst]+sq3*par_ab_pair_[Pair1][2][bst])*z_gl_[sec][i2]-0.5*(par_ab_pair_[Pair1][1][bst]+sq3*par_ab_pair_[Pair1][3][bst]);
	      if(i2==1 || i2==4) dx_[i2]=XVU_Ch_[sec][i2][best2]-0.5*(par_ab_pair_[Pair1][0][bst]-sq3*par_ab_pair_[Pair1][2][bst])*z_gl_[sec][i2]-0.5*(par_ab_pair_[Pair1][1][bst]-sq3*par_ab_pair_[Pair1][3][bst]);
	      Chi2_ndf_pair_[Pair1][bst]= Chi2_ndf_pair_[Pair1][bst]+dx_[i2]*dx_[i2]/(sigm_gl2[i2]*sigm_gl2[i2]);

	      //   cout<<"best2 "<<best2 <<" i2 "<<i2<<" dx_ "<<dx_[i2]<<" XVU_Ch2 "<<XVU_Ch_[sec][i2][best2]<<" Chi2_ndf_Ch1_2 "<<Chi2_ndf_pair_[Pair1][bst]<<" z_gl2 "<<z_gl_[sec][i2]<<endl;

	    }// if( Wires_Ch2[i2][best2]>-1){
	  }
	  //  cout<<endl;

	  if (Nhits_Ch_[fir][best1]+Nhits_Ch_[sec][best2]> 4)
	    Chi2_ndf_pair_[Pair1][bst]= Chi2_ndf_pair_[Pair1][bst]/(Nhits_Ch_[fir][best1]+Nhits_Ch_[sec][best2]-4);
	  //  cout<<" Chi2_ndf_pair "<<Chi2_ndf_pair_[Pair1][bst]<<endl;

	  ind_best_pair_[Pair1][bst]= bst;	    
	 
  }//< Nbest_Ch12_g_l

}//SegmentFit



void BmnMwpcTrackFinderSRC::FillFitMatrix(Int_t chN, Double_t** AA, Float_t** z, Float_t* sigmm2, Int_t* h_) {

  // AA - matrix to be filledlayers)
  // sigm2 - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)
  // Float_t z2_[nPlanes];
  Float_t z2_[6] = {z[chN][0] * z[chN][0], z[chN][1] * z[chN][1], z[chN][2] * z[chN][2], z[chN][3] * z[chN][3], z[chN][4] * z[chN][4], z[chN][5] * z[chN][5]}; //cm

  //  cout<<" chN "<<chN<<endl;
  for (Int_t i1 = 0; i1 < 6; i1++){
    // cout<<" z "<<z[chN][i1]<<" z2_ "<<z2_[i1]<<endl;
  }

  AA[0][0] += 2 * z2_[0] * h_[0] / sigmm2[0] 
           +      z2_[2] * h_[2] / (2 * sigmm2[2]) 
           +      z2_[1] * h_[1] / (2 * sigmm2[1]) 
           +  2 * z2_[3] * h_[3] / sigmm2[3] 
           +      z2_[5] * h_[5] / (2 * sigmm2[5]) 
           +      z2_[4] * h_[4] / (2 * sigmm2[4]); //Ax

  AA[0][1] += 2 * z[chN][0] * h_[0] / sigmm2[0] 
           +      z[chN][2] * h_[2] / (2 * sigmm2[2]) 
           +      z[chN][1] * h_[1] / (2 * sigmm2[1]) 
           +  2 * z[chN][3] * h_[3] / sigmm2[3] 
           +      z[chN][5] * h_[5] / (2 * sigmm2[5]) 
           +      z[chN][4] * h_[4] / (2 * sigmm2[4]); //Bx

  AA[0][2] += sq3 * (z2_[2] * h_[2] / (2 * sigmm2[2]) 
           -         z2_[1] * h_[1] / (2 * sigmm2[1]) 
           +         z2_[5] * h_[5] / (2 * sigmm2[5]) 
           -         z2_[4] * h_[4] / (2 * sigmm2[4])); //Ay

  AA[0][3] += sq3 * (z[chN][2] * h_[2] / (2 * sigmm2[2]) 
           -         z[chN][1] * h_[1] / (2 * sigmm2[1]) 
           +         z[chN][5] * h_[5] / (2 * sigmm2[5]) 
		       -         z[chN][4] * h_[4] / (2 * sigmm2[4])); //By

  AA[1][0] = AA[0][1];

  AA[1][1] +=   2 * h_[0] / sigmm2[0] 
           +  0.5 * h_[2] / sigmm2[2] + 0.5 * h_[1] / sigmm2[1] 
           +    2 * h_[3] / sigmm2[3] + 0.5 * h_[5] / sigmm2[5] 
           +  0.5 * h_[4] / sigmm2[4];

  AA[1][2] += sq3 * (z[chN][2] * h_[2] / sigmm2[2] 
           - z[chN][1] * h_[1] / sigmm2[1] 
           + z[chN][5] * h_[5] / sigmm2[5] 
           - z[chN][4] * h_[4] / sigmm2[4]) * 0.5;

  AA[1][3] += sq3 * (h_[2] / sigmm2[2] 
           -         h_[1] / sigmm2[1] 
           +         h_[5] / sigmm2[5] 
           -         h_[4] / sigmm2[4]) * 0.5;

  AA[2][0] = AA[0][2];

  AA[2][1] = AA[1][2];

  AA[2][2] += 3.0 * (z2_[2] * h_[2] / sigmm2[2] 
           +         z2_[1] * h_[1] / sigmm2[1] 
           +         z2_[5] * h_[5] / sigmm2[5] 
           +         z2_[4] * h_[4] / sigmm2[4]) * 0.5;

  AA[2][3] += 3.0 * (z[chN][2] * h_[2] / sigmm2[2] 
           +         z[chN][1] * h_[1] / sigmm2[1] 
           +         z[chN][5] * h_[5] / sigmm2[5] 
           +         z[chN][4] * h_[4] / sigmm2[4])   * 0.5;

  AA[3][0] = AA[0][3];
  AA[3][1] = AA[1][3];
  AA[3][2] = AA[2][3];
  AA[3][3] += 3.0 * (0.5 * h_[2] / sigmm2[2] 
           +  0.5 *        h_[1] / sigmm2[1] 
           +  0.5 *        h_[5] / sigmm2[5] 
           +  0.5 *        h_[4] / sigmm2[4]);

}


void BmnMwpcTrackFinderSRC::FillFreeCoefVector(Int_t ichNum, Double_t* F, Float_t*** XVU_, Int_t ise, Float_t** z, Float_t* sigmm2, Int_t* h_) {
  // F - vector to be filled
  // XVU_ - coordinates of segment in chamber (Is it correct definition?)
  // segIdx - index of current segment
  // z - local z-positions of planes(layers)
  // sigmm2 - square of sigma
  // h_ - array to include/exclude planes (h_[i] = 0 or 1)

  for(Int_t i = 0; i < 6;i++){
    // cout<<" sigm2 "<< sigmm2[i]<<" h "<<h_[i]<<" XVU_ "<<XVU_[ichNum][i][ise]<<endl;
  }  

  F[0] += 
    2 * XVU_[ichNum][0][ise] * z[ichNum][0] * h_[0] / sigmm2[0] 
    + 
    XVU_[ichNum][1][ise]  * z[ichNum][1] * h_[1] / sigmm2[1] 
    + 
    XVU_[ichNum][2][ise] * z[ichNum][2] * h_[2] / sigmm2[2] 
    + 
    2 * XVU_[ichNum][3][ise] * z[ichNum][3] * h_[3] / sigmm2[3] 
    + 
    XVU_[ichNum][4][ise] * z[ichNum][4] * h_[4] / sigmm2[4] 
    + 
    XVU_[ichNum][5][ise] * z[ichNum][5] * h_[5] / sigmm2[5];

  F[1] += 2 * XVU_[ichNum][0][ise] * h_[0] / sigmm2[0] + XVU_[ichNum][1][ise] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * h_[2] / sigmm2[2] + 2 * XVU_[ichNum][3][ise] * h_[3] / sigmm2[3] + XVU_[ichNum][4][ise] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * h_[5] / sigmm2[5];

  F[2] += sq3*(-XVU_[ichNum][1][ise] * z[ichNum][1] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * z[ichNum][2] * h_[2] / sigmm2[2] - XVU_[ichNum][4][ise] * z[ichNum][4] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * z[ichNum][5] * h_[5] / sigmm2[5]);

  F[3] +=  sq3*(-XVU_[ichNum][1][ise] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * h_[2] / sigmm2[2] - XVU_[ichNum][4][ise] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * h_[5] / sigmm2[5]);

}

void BmnMwpcTrackFinderSRC::InverseMatrix(Double_t** AA, Double_t** bb) {
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

void BmnMwpcTrackFinderSRC::FillEfficiency(Int_t ChN, Float_t ***XVU_Ch_, Int_t **Nhits_C, Int_t MinHits, Int_t ind_best, Float_t min_distX, Float_t min_distY ){ 
  //cout<<" FillEfficiency : ChN "<<ChN<<" min_distX "<<min_distX<<" min_distY "<<min_distY<<" Nhits_Ch["<<ChN<<"]["<<ind_best<<"] "<<Nhits_Ch[ChN][ ind_best ]<<endl;
        // segIdx - index of current segment ch2 or ch3 // Int_t ind_best_Ch[5]
	//4p&4p -> all matched / Efficiency per layer

    if (fabs(min_distX)< 5. && fabs(min_distY)< 5.5) { 
      for(int i1 = 0 ; i1 < 6; i1++){
	//	cout<<" XVU_Ch_["<<ChN<<"]["<<i1<<"]["<< ind_best<<"] = "<<XVU_Ch_[ChN][i1][ ind_best]<<endl;
	if( XVU_Ch_[ChN][i1][ ind_best ] > -999.  && Nhits_Ch[ChN][ ind_best ] == MinHits)  continue;//segIdx[ChN][j]

	Denom_Ch.at(ChN)->Fill(i1);
			 
	if(XVU_Ch_[ChN][i1][ ind_best ] > -999.){
	  Nomin_Ch.at(ChN)->Fill(i1);		    
	}
      }// i1 
     }//min_distX  
}// FillEfficiency



InitStatus BmnMwpcTrackFinderSRC::Init() {
  if (!expData)
      return kERROR;
  if (fDebug) cout << "BmnMwpcTrackFinderSRC::Init()" << endl;
  FairRootManager* ioman = FairRootManager::Instance();

  fBmnMwpcSegmentsArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
  if (!fBmnMwpcSegmentsArray)
  {
    cout<<"BmnMwpcTrackFinderSRC::Init(): branch "<<fInputBranchName<<" not found! Task will be deactivated"<<endl;
    SetActive(kFALSE);
    return kERROR;
  }

  fBmnMwpcTracksArray = new TClonesArray(fOutputBranchName.Data());
  ioman ->Register(fOutputBranchName.Data(), "MWPC", fBmnMwpcTracksArray, kTRUE);

  fMwpcGeo   = new BmnMwpcGeometrySRC(fRunPeriod); // period number 7
  kNChambers = fMwpcGeo -> GetNChambers(); 
  kNPlanes   = fMwpcGeo -> GetNPlanes(); // 6
  if (fDebug) printf("fRunPeriod: %d %d %d\n", fRunPeriod, kNChambers, kNPlanes);

  kMinHits = 4;
  kChi2_Max = 20.;
  kmaxPairs = 10;

  dw = fMwpcGeo->GetWireStep();//0.25; // [cm] // wires step
  dw_half = 0.5*dw;
  sq3 = sqrt(3.);
  sq12 = sqrt(12.);
  sigma = dw/sq12;
  kMiddlePl = 47.25;
  kZ_target = -648.4;//cm
  cout<<" kZ_target "<<kZ_target<<" kNumPairs "<<kNumPairs<<endl;

  ChCent          = new TVector3[kNChambers];
  ZCh             = new Float_t[kNChambers];
  kZmid           = new Float_t[kNChambers];
  shift           = new Float_t*[kNChambers];
  shift_pair      = new Float_t*[kNumPairs];
  kZ_midle_pair   = new Float_t[kNumPairs];
  XVU_Ch          = new Float_t**[kNChambers];
  par_ab_Ch       = new Double_t**[kNChambers];
  par_ab_pair     = new Double_t**[kNumPairs];// par_ab_Ch1_2 = new Double_t*[4];
  kPln            = new Int_t*[kNChambers];
  kZ_loc          = new Float_t*[kNChambers];
  z_gl            = new Float_t*[kNChambers];
  Chi2_match_pair = new Double_t*[kNumPairs];
  Chi2_ndf_pair   = new Double_t*[kNumPairs];
  ind_best_pair   = new Int_t*[kNumPairs]; //ind_best_Ch1_2 = new Int_t[5]; 
  Chi2_ndf_Ch     = new Double_t*[kNChambers];
  Nhits_Ch        = new Int_t*[kNChambers];
  sigm2           = new Float_t[kNPlanes]; 
  ipl             = new Int_t[kNPlanes];
  Nbest_pair      = new Int_t[kNumPairs];
  Nbest_Ch        = new Int_t[kNChambers];
  ind_best_Ch     = new Int_t*[kNChambers];

  for(Int_t i = 0; i < kNChambers; i++){// for(Int_t i = kCh_min; i < kCh_max; i++){
    TH1D *h;
    h = new TH1D(Form("par_Ax_Ch%d",i), Form("par_Ax_Ch%d",i), 100, -.4, .4);      fList.Add(h); hpar_Ax_Ch.push_back(h);
    h = new TH1D(Form("par_Bx_Ch%d",i), Form("par_Bx_Ch%d",i), 100, -10., 10.0);   fList.Add(h); hpar_Bx_Ch.push_back(h);
    h = new TH1D(Form("par_Ay_Ch%d",i), Form("par_Ay_Ch%d",i), 100, -.4, .4);      fList.Add(h); hpar_Ay_Ch.push_back(h);
    h = new TH1D(Form("par_By_Ch%d",i), Form("par_By_Ch%d",i), 100, -10., 10.0);   fList.Add(h); hpar_By_Ch.push_back(h);
    h = new TH1D(Form("Chi2_ndf_Ch%d",i), Form("Chi2_ndf_Ch%d",i), 100, 0., 20.0); fList.Add(h); hChi2_ndf_Ch.push_back(h);
    h = new TH1D(Form("Nomin_Ch%d",i), Form("Nomin_Ch%d",i), 6, 0., 6.);           fList.Add(h); Nomin_Ch.push_back(h);
    h = new TH1D(Form("Denom_Ch%d",i), Form("Denom_Ch%d",i), 6, 0., 6.);           fList.Add(h); Denom_Ch.push_back(h);
    h = new TH1D(Form("Efficiency_Ch%d",i), Form("Efficiency_Ch%d",i), 6, 0., 6.); fList.Add(h); Eff_Ch.push_back(h);
    shift[i]        = new Float_t[4];
    ChCent[i]       = fMwpcGeo -> GetChamberCenter(i);
    ZCh[i]          = ChCent[i].Z();
    shift[i][0]     = fMwpcGeo -> GetTx(i);
    shift[i][2]     = fMwpcGeo -> GetTy(i);
    shift[i][1]     = ChCent[i].X();
    shift[i][3]     = ChCent[i].Y();
    kPln[i]         = new Int_t[kNPlanes];
    kZ_loc[i]       = new Float_t[kNPlanes];
    z_gl[i]         = new Float_t[kNPlanes];
    XVU_Ch[i]       = new Float_t*[kNPlanes];
    par_ab_Ch[i]    = new Double_t*[4];
    Nhits_Ch[i]     = new Int_t[kBig];
    Chi2_ndf_Ch[i]  = new Double_t[kBig];
    ind_best_Ch[i]  = new Int_t[kmaxPairs];
  }

   //----- hists booking -----
  for(Int_t i = 0; i < kNChambers; i++){// for(Int_t i = kCh_min; i < kCh_max; i++){
    if (i== 0 || i== 2) { kZmid[i] = (ZCh[i]     - ZCh[i + 1]) *  0.5;}
    if (i== 1 || i== 3) { kZmid[i] = (ZCh[i - 1] - ZCh[i])     * -0.5;}
    printf("Chamber %d Z: %f Zmid: %f\n", i, ZCh[i], kZmid[i]);
  }

  for (int i=0; i < kNumPairs; ++i){ 
    TH1D *h;
    h = new TH1D(Form("Chi2_match_pair%d",i), Form("Chi2_match_pair%d",i), 100, 0., 100.0);  fList.Add(h); hChi2_match_pair.push_back(h);
    h = new TH1D(Form("par_Ax_pair%d",i), Form("slopeX pair%d; ; Events",i), 100, -2.3, 2.3);    fList.Add(h); hpar_Ax_pair.push_back(h);
    h = new TH1D(Form("par_Bx_pair%d",i), Form("posX pair%d; [cm]; Events",i), 100, -10., 10.0);   fList.Add(h); hpar_Bx_pair.push_back(h);
    h = new TH1D(Form("par_Ay_pair%d",i), Form("slopeY pair%d; ; Events",i), 100, -2.3, 2.3);    fList.Add(h); hpar_Ay_pair.push_back(h);
    h = new TH1D(Form("par_By_pair%d",i), Form("posY pair%d; [cm]; Events",i), 100, -10., 10.0);   fList.Add(h); hpar_By_pair.push_back(h);
    h = new TH1D(Form("theta_pair%d",i),Form("theta_pair%d; degrees; Events",i), 160,  0., 8.);     fList.Add(h); hpar_theta_pair.push_back(h);
    h = new TH1D(Form("phi_pair%d",i), Form("phi_pair%d;   degrees; Events",i), 380, -190., 190.);  fList.Add(h); hpar_phi_pair.push_back(h);
    Chi2_match_pair[i] = new Double_t[kmaxPairs];
    Chi2_ndf_pair[i]   = new Double_t[kmaxPairs];
    ind_best_pair[i]   = new Int_t[kmaxPairs];
    par_ab_pair[i]     = new Double_t*[4];
    shift_pair[i]      = new Float_t[4];
  }

  hdX_target =       new TH1D("dX_target", " dX(pair0-pair1) in target;[cm]; Events ", 100, -20.,20.);                   fList.Add(hdX_target);
  hdY_target =       new TH1D("dY_target", " dY(pair0-pair1) in target;[cm]; Events ", 100, -20.,20.);                   fList.Add(hdY_target);
  hX_in_target =     new TH1D("X_in_target", " posX pair0 in target;[cm]; Events ", 100, -10.,10.);                   fList.Add(hX_in_target);
  hY_in_target =     new TH1D("Y_in_target", " posY pair0 in target;[cm]; Events ", 100, -10.,10.);                   fList.Add(hY_in_target);
  hAx_bx_in_target=  new TH2D("Ax_bx_in_target", "slopeX vs posX in target; posX[cm]; slopeX", 100, -10.,10., 100, -2.3, 2.3);  fList.Add(hAx_bx_in_target);
  hAy_by_in_target=  new TH2D("Ay_by_in_target", "slopeY vs posY in target; posY[cm]; slopeY", 100, -10.,10., 100, -2.3, 2.3);  fList.Add(hAy_by_in_target);
  hY_X_in_target  =  new TH2D("Y_X_in_target", "posY vs posX (pair0) in target; X[cm]; Y[cm]", 100, -10.,10., 100,  -10., 10.); fList.Add(hY_X_in_target);
  hdist_target    =  new TH1D("dist_target","dist_target;[cm];Events", 100, 0., 10.0);                                       fList.Add(hdist_target);
  hNsecondaries   =  new TH1D("Nsecondaries", "Nsecondaries", 11, 0., 11.);                                      fList.Add(hNsecondaries);

  Int_t i1 = 0; 
  for(Int_t i = 0; i < kNumPairs; i++){ 
     i1=i;
    if ( fRunPeriod == 6 && i == 0 ) i1=2; 
    if ( i == 1) i1=2;
    kZ_midle_pair[i] =  ZCh[i1] + kZmid[i1+1];

    shift_pair[i][0]= (shift[i1+1][1] - shift[i1][1])/( ZCh[i1+1] - ZCh[i1] );
    shift_pair[i][2]= (shift[i1+1][3] - shift[i1][3])/( ZCh[i1+1] - ZCh[i1] );
    shift_pair[i][1]=  0.5*(shift[i1+1][1] + shift[i1][1]);
    shift_pair[i][3] = 0.5*(shift[i1+1][3] + shift[i1][3]);

    cout<<" i "<<i<<" kZ_midle_pair[i] "<<kZ_midle_pair[i]<<" i1 "<<i1<<" i1+1 "<<i1+1<<" -( ZCh[i1]- ZCh[i1+1] )= "<<-( ZCh[i1]- ZCh[i1+1] )<<endl;
  }

    for(Int_t ichh = 0; ichh < kNChambers; ichh++){// for(Int_t ichh = kCh_min; ichh < kCh_max; ichh++){
      for(int ii = 0; ii < 6; ii++){

	if ( fRunPeriod == 6 ){

	  if ( ichh == 0 || ichh == 1){
	    kZ_loc[ichh][ii] = -0.5 + ii;
	    if(ii == 4) { kZ_loc[ichh][ii] = -2.5;}
	    if(ii == 5) { kZ_loc[ichh][ii] = -1.5;}
	  }
	}

	if ( fRunPeriod == 7 ){
	if (ichh == 0 || ichh == 1){
	  kZ_loc[0][0] = -1.5;
	  kZ_loc[0][1] = -0.5;
	  kZ_loc[0][2] =  0.5;
	  kZ_loc[0][3] =  1.5;
	  kZ_loc[0][4] =  2.5;
	  kZ_loc[0][5] = -2.5;

	  kZ_loc[1][0] = -1.5;
	  kZ_loc[1][1] = -2.5;
	  kZ_loc[1][2] =  2.5;
	  kZ_loc[1][3] =  1.5;
	  kZ_loc[1][4] =  0.5;
	  kZ_loc[1][5] = -0.5;
	}
	}
	if ( ichh == 2 || ichh == 3){
	  kZ_loc[ichh][ii] = -0.5 + ii;
	  if(ii == 4) { kZ_loc[ichh][ii] = -2.5;}
	  if(ii == 5) { kZ_loc[ichh][ii] = -1.5;}
	}//if ( ich == 0 || ich == 1)

       	z_gl[ichh][ii] =  kZmid[ichh] + kZ_loc[ichh][ii];
 
	//	cout<<" ich "<<ichh<<" ii "<<ii<<" kZ_loc "<<kZ_loc[ichh][ii]<<" z_gl "<<z_gl[ichh][ii]<<endl;	
      }
    }//ich   
    cout<< endl;


   for(Int_t ii = 0; ii < kNChambers; ii++){// for(Int_t ii = kCh_min; ii < kCh_max; ii++){ 
      for(Int_t iPla=0; iPla < kNPlanes; iPla++){
	XVU_Ch[ii][iPla] = new Float_t[kBig];
      }
      for(Int_t i=0; i<4; i++){
	par_ab_Ch[ii][i] = new Double_t[kBig];
      }
    }

    for(Int_t ip=0; ip < kNumPairs; ip++){
      for(Int_t i4=0; i4<4; i4++){
	par_ab_pair[ip][i4] = new Double_t[kmaxPairs];	
      }
    }

    return kSUCCESS;
}//Init


void BmnMwpcTrackFinderSRC::PrepareArraysToProcessEvent(){

  fBmnMwpcTracksArray->Clear();

  // Clean and initialize arrays:

 for(Int_t iCh = 0; iCh < kNChambers; iCh++){// for(Int_t iCh = kCh_min; iCh < kCh_max;  iCh++){
    Nbest_Ch[iCh] = 0;

     for(Int_t iPlane=0; iPlane<kNPlanes; iPlane++){
	for(Int_t iBig=0; iBig<kBig; iBig++){
	  XVU_Ch[iCh][iPlane][iBig] = -999.;
	}//iBig
     }//iPlane

     for(Int_t ii=0; ii<4; ii++){
      for(Int_t jj=0; jj<100; jj++){
	par_ab_Ch[iCh][ii][jj] = 999.;
      }
    }

     for(Int_t iBig=0; iBig<kBig; iBig++){     
       Nhits_Ch[iCh][iBig] = 0;
       Chi2_ndf_Ch[iCh][iBig] = 0;
     }

     for(Int_t i=0; i< kmaxPairs; i++){
       ind_best_Ch[iCh][i] = -1;    
     }
  }//iCh 

  for(Int_t iPl=0; iPl<kNPlanes; iPl++){     
    sigm2[iPl] = sigma*sigma;  
    ipl[iPl] = 6;
  }

  for(Int_t ip=0; ip < kNumPairs; ip++){
    Nbest_pair[ip] = 0; //Nbest_Ch12_gl = 0;

      for(Int_t i4=0; i4 < 4; i4++){
	for(Int_t i5=0; i5 < kmaxPairs; i5++){	
	  par_ab_pair[ip][i4][i5] =999.;//par_ab_Ch1_2[ii][jj] = 999.;	
	}
      }
    for(Int_t i=0; i <kmaxPairs; i++){	
      Chi2_match_pair[ip][i] = 999.;
      Chi2_ndf_pair[ip][i] = 999.;
      ind_best_pair[ip][i]= -1;
    }
  }

}//PrepareArraysToProcessEvent

void BmnMwpcTrackFinderSRC::Finish() {
  delete fMwpcGeo;
  printf("MWPC track finder: write hists to file... ");
  fOutputFileName = TString("hMWPCtracks.root");
  TFile file(fOutputFileName, "recreate");
  
  for(Int_t iCh = 0; iCh < kNChambers; iCh++){
    Eff_Ch.at(iCh)->Sumw2();
    Eff_Ch.at(iCh)->Divide(Nomin_Ch.at(iCh),Denom_Ch.at(iCh),100,1);
  }

  fList.Write();
  printf("done\n");
    
  cout << "Work time of the MWPC track finder: " << workTime << " s" << endl;
}//Finish

ClassImp(BmnMwpcTrackFinderSRC)

