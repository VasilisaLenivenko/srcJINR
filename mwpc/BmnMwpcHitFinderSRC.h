// @(#)bmnroot/mwpc:$Id$


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcHitFinderSRC                                                        //
//                                                                            //
// Implementation of an algorithm developed by                                // 
// Vasilisa Lenivenko                                                         //
// to the BmnRoot software                                                    //
//                                                                            //
// The algorithm serves for searching for hits                                //
// in the Multi Wire Prop. Chambers of the BM@N experiment                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BMNMWPCHITFINDERSRC_H
#define	BMNMWPCHITFINDERSRC_H

#include <map>
#include <algorithm>
#include <Rtypes.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <TString.h>
#include  "FairTask.h"
#include  "BmnMwpcHit.h"
#include  "BmnMwpcDigit.h"
#include  "BmnMwpcGeometrySRC.h"
#include  "FairTask.h"
#include "TList.h"
#include "BmnTrack.h"
#include <vector>

#include "TH1D.h"
#include "TH2D.h"

class TH1D;
class TH2D;

using namespace std;

class BmnMwpcHitFinderSRC : public FairTask {
public:
  /** Default constructor **/
  BmnMwpcHitFinderSRC() {};

  /** Constructor to be used **/
  BmnMwpcHitFinderSRC(Bool_t, Int_t);

  /** Destructor **/
  virtual ~BmnMwpcHitFinderSRC();

  /** Virtual method Init **/
  virtual InitStatus Init();

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);

  /** Virtual method Finish **/
  virtual void Finish();

private:
  Bool_t expData;
  Bool_t fDebug   = 1;
  Bool_t fDoTest  = 1;
  UInt_t fEventNo; // event counter

  TString fInputBranchName;
  TString fOutputBranchName;

  /** Input array of MWPC digits **/
  TClonesArray* fBmnMwpcDigitArray;
  /** Output array of MWPC hits **/
  TClonesArray* fBmnMwpcSegmentsArray; 

  TString fOutputFileName;

  Int_t nInputDigits;  // max. number of found digits per plane
  Int_t nTimeSamples;  // 

  BmnMwpcGeometrySRC* fMwpcGeometrySRC;

  TList fList;

  Int_t fRunPeriod;
  Short_t kNChambers, kNPlanes, kNWires;
  Int_t kNumPairs;
  Int_t kCh_min, kCh_max;
  Int_t kBig;
  TVector3 *ChCent;
  Float_t *Zmid;
  Float_t *ChZ;
//  vector<TVector3> ChCent;
//  vector<Double_t> ChZ, Zmid;
  vector<TH1D*> hNp_best_Ch, hNbest_Ch, hOccupancy, hTime;
  Int_t kMinHits;
  Int_t kmaxSeg;
  Int_t kChMaxAllWires;
  Double_t kChi2_Max;
  Float_t dw, dw_half;
  Double_t sq3, sq12, sigma;
  Short_t kMiddlePl;

// Arrays 
  Int_t    *Nseg_Ch; // Int_t Nseg_Ch1;
  Int_t    *Nbest_Ch;//Int_t Nbest_Ch1;
  Float_t  *sigm2;
  Float_t  *z2;
  Int_t    *ipl;
  Float_t  *dX_i;
  Float_t  **kZ_loc;
  Float_t  **z_gl;
  Int_t    **kPln;
  Int_t    **iw;
  Int_t    **iw_Ch;
  Int_t    **Nhits_Ch;//  Int_t *Nhits_Ch1; 
  Float_t  **XVU;// Float_t *XVU1;
  Float_t  **XVU_cl;// Float_t *XVU_cl1;
  Int_t    **ind_best_Ch;// Int_t *ind_best_Ch1;
  Int_t    **best_Ch_gl;// Int_t *best_Ch1_gl;
  Double_t **Chi2_ndf_Ch;// Double_t *Chi2_ndf_Ch1;
  Double_t **Chi2_ndf_best_Ch;//Double_t *Chi2_ndf_best_Ch1;
  Int_t    ***wire_Ch;  // Int_t **wire_Ch1;
  Float_t  ***xuv_Ch;  //  Float_t **xuv_Ch1;
  Int_t    ***Wires_Ch;// Int_t **Wires_Ch1;   
  Int_t    ***clust_Ch;// Int_t **clust_Ch1;
  Float_t  ***XVU_Ch;//  Float_t **XVU_Ch1;  
  Double_t ***par_ab_Ch;// Double_t **par_ab_Ch1;

  Double_t **matrA;
  Double_t **matrb;
  Double_t **Amatr;
  Double_t **bmatr;


  //functions for Vasilisa method:
  void PrepareArraysToProcessEvent();
  void SegmentFinder(Int_t, Int_t***, Int_t***,  Float_t***,  Int_t**, Int_t** , Int_t *, Int_t ***, Float_t ***, Int_t, Short_t , Int_t);
  void ProcessSegments(Int_t, Double_t ,Float_t , Float_t **, Int_t , Int_t * , Int_t **,Int_t ***,Int_t ***, Float_t ***,Int_t * , Int_t **, Double_t **, Double_t **, Double_t ***, Int_t ,Int_t* , Float_t**,  Float_t**, Double_t);
  void FillFreeCoefVectorXUV(Int_t, Double_t*, Float_t**,  Float_t**, Float_t*, Int_t*);
  void FillFreeCoefVector(Int_t, Double_t*, Float_t**, Int_t, Float_t*, Float_t*, Int_t*, Int_t);
  void FillFitMatrix(Int_t, Double_t **, Float_t **, Float_t *, Int_t *);
  void InverseMatrix(Double_t**, Double_t**);

  ClassDef(BmnMwpcHitFinderSRC, 1);
};

#endif	

