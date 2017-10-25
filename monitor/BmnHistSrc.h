#ifndef BMNHISTSRC_H
#define BMNHISTSRC_H 1
#include <TNamed.h>
#include <exception>
#include <stdlib.h>
#include <vector>
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFolder.h"
#include "TString.h"
#include "THttpServer.h"

#include "BmnHist.h"
#include "BmnTrigWaveDigit.h"

#define SRC_ROWS 13
#define SRC_COLS  3

using namespace std;

class BmnHistSrc : public BmnHist {
public:

    BmnHistSrc(TString title = "SRC", TString path = "");
    virtual ~BmnHistSrc();
    void Reset();
    void Register(THttpServer *serv);
    void SetDir(TFile *outFile = NULL, TTree *recoTree = NULL);
    void DrawBoth();
    void FillFromDigi(DigiArrays *fDigiArrays);
    BmnStatus  SetRefRun(Int_t id);
    
private:
    vector<TString> Names;
    vector<vector<TH1F*> > hists;
    TCanvas *canvas;
    vector<PadInfo*> canPads;
    const TString trigNames[13] = 
    {
        "BC1", "BC2", "BC3", "BC4", "VC",
        "X1_Left", "X1_Right", "X2_Left", "X2_Right",
        "Y1_Left", "Y1_Right", "Y2_Left", "Y2_Right"
    };

    ClassDef(BmnHistSrc, 1)
};

#endif /* BMNHISTSRC_H */

