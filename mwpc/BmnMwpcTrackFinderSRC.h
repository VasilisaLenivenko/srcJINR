// @(#)bmnroot/mwpc:$Id$


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcTrackFinder                                                         //
//                                                                            //
//                                                                            //
// The algorithm serves for searching for track segments                      //
// in the MWPC of the BM@N experiment                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BMNMWPCTRACKFINDERSRC_H
#define BMNMWPCTRACKFINDERSRC_H 

#include <TMath.h>
#include <TNamed.h>
#include <TClonesArray.h>
#include <TString.h>
#include "FairTask.h"
#include "BmnMwpcTrack.h"
#include "BmnMwpcHit.h"
#include "BmnMwpcGeometrySRC.h"
#include "BmnEnums.h"
#include "BmnMath.h"
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <sstream>

class TH1D;
class TH2D;

using namespace std;
using namespace TMath;

class BmnMwpcTrackFinderSRC : public FairTask {
public:

    BmnMwpcTrackFinderSRC() {};
    BmnMwpcTrackFinderSRC(Bool_t, Int_t);
    virtual ~BmnMwpcTrackFinderSRC();
    
    virtual InitStatus Init();

    virtual void Exec(Option_t* opt);

    virtual void Finish();
   
private:
    Bool_t expData;
    Bool_t fDebug = 1;
    UInt_t fEventNo; // event counter
    Int_t fRunPeriod;
    
    TString fInputBranchName;
    TString fOutputBranchName;
    TString fOutputFileName;
   
    /** Input array of MWPC hits **/
    TClonesArray* fBmnMwpcSegmentsArray;
    
    /** Output array of MWPC tracks **/
    TClonesArray* fBmnMwpcTracksArray; 
    
    BmnMwpcGeometrySRC* fMwpcGeo;

    vector<TH1D*> hpar_Ax_Ch, hpar_Bx_Ch, hpar_Ay_Ch, hpar_By_Ch, hChi2_match_pair, hChi2_ndf_Ch, hpar_Ax_pair, hpar_Bx_pair, hpar_Ay_pair, hpar_By_pair, hpar_theta_pair, hpar_phi_pair, Nomin_Ch, Denom_Ch, Eff_Ch;
    TH1D *hdX_target, *hdY_target, *hX_in_target, *hY_in_target, *hdist_target, *hNsecondaries;
    TH2D *hAx_bx_in_target, *hAy_by_in_target, *hY_X_in_target;

    Short_t kNChambers;
    Short_t kNPlanes;
    Int_t kBig;
    Int_t kCh_min;
    Int_t kCh_max;
    Int_t kNumPairs;
    TVector3 *ChCent;

    Float_t *kZmid;
    Float_t *ZCh;
    Float_t **kZ_loc;
    Float_t *kZ_midle_pair;
 
    Float_t **z_gl;
    Float_t **shift; 
    Float_t **shift_pair;
    
    Float_t kZ_to_pole;
    Float_t kZ_target;
    Int_t kMinHits;
    Int_t kmaxPairs;
    Double_t kChi2_Max;

    Float_t dw;
    Float_t dw_half;
    Double_t sq3;
    Double_t sq12;
    Double_t sigma;
    Short_t kMiddlePl;
   
  
    Int_t **kPln;
    Float_t *sigm2;
    Int_t *ipl;
    Double_t **matrA;
    Double_t **matrb;
    Double_t **Amatr;
    Double_t **bmatr;

    Double_t ***par_ab_Ch;
    Double_t ***par_ab_pair;
    Float_t ***XVU_Ch;
    Double_t **Chi2_match_pair;
    Double_t **Chi2_ndf_pair;
    Double_t **Chi2_ndf_Ch;
    Int_t **ind_best_pair;
    Int_t **ind_best_Ch;
    Int_t *Nbest_pair;
    Int_t *Nbest_Ch;
    Int_t **Nhits_Ch;

    TList fList;

    void PrepareArraysToProcessEvent();
    void SegmentParamAlignment(Int_t, Int_t *,  Double_t ***, Float_t **);
    void SegmentMatching(Int_t,Int_t *, Double_t ***, Float_t *, Int_t **,  Int_t *,Double_t **, Float_t ***, Int_t **);
    void SegmentFit(Int_t, Float_t **, Float_t *, Int_t *,  Int_t **, Double_t ***, Double_t **, Float_t ***, Int_t **, Int_t **);
    void FillFitMatrix(Int_t, Double_t** , Float_t** , Float_t* , Int_t*);
    void FillFreeCoefVector(Int_t , Double_t* , Float_t*** , Int_t , Float_t**, Float_t*, Int_t*);
    void InverseMatrix(Double_t**, Double_t**);
    void PairMatching( Int_t *, Double_t ***,  Float_t *);
    void FillEfficiency(Int_t, Float_t ***, Int_t **, Int_t, Int_t, Float_t, Float_t);
         
    ClassDef(BmnMwpcTrackFinderSRC, 1)
};

#endif
