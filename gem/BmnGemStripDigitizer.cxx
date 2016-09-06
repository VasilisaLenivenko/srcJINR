#include "BmnGemStripDigitizer.h"
#include "CbmMCTrack.h"

#include "BmnGemStripStationSet_RunSummer2016.h"

static Float_t workTime = 0.0;
static int entrys = 0;

BmnGemStripDigitizer::BmnGemStripDigitizer()
: fOnlyPrimary(kFALSE), fStripMatching(kTRUE) {

    fInputBranchName = "StsPoint";
    fOutputDigitsBranchName = "BmnGemStripDigit";
    fOutputDigitMatchesBranchName = "BmnGemStripDigitMatch";

    fVerbose = 1;

    fCurrentConfig = BmnGemStripConfiguration::None;
}

BmnGemStripDigitizer::~BmnGemStripDigitizer() {

}

InitStatus BmnGemStripDigitizer::Init() {

    if(fVerbose) cout << "\nBmnGemStripDigitizer::Init()\n ";

    //if GEM configuration is not set -> return a fatal error
    if(!fCurrentConfig) Fatal("BmnGemStripDigitizer::Init()", " !!! Current GEM config is not set !!! ");

    if(fVerbose && fOnlyPrimary) cout << "  Only primary particles are processed!!! " << endl;

    if(fVerbose && fStripMatching) cout << "  Strip Matching is activated!!! " << endl;
    else cout << "  Strip Matching is deactivated!!! " << endl;

    FairRootManager* ioman = FairRootManager::Instance();

    fBmnGemStripPointsArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
    fMCTracksArray = (TClonesArray*) ioman->GetObject("MCTrack");

    fBmnGemStripDigitsArray = new TClonesArray(fOutputDigitsBranchName);
    ioman->Register(fOutputDigitsBranchName, "GEM", fBmnGemStripDigitsArray, kTRUE);

    if(fStripMatching) {
        fBmnGemStripDigitMatchesArray = new TClonesArray("BmnMatch");
        ioman->Register(fOutputDigitMatchesBranchName, "GEM", fBmnGemStripDigitMatchesArray, kTRUE);
    }

    if(fVerbose) cout << "BmnGemStripDigitizer::Init() finished\n\n";
    return kSUCCESS;
}

void BmnGemStripDigitizer::Exec(Option_t* opt) {
    clock_t tStart = clock();
    fBmnGemStripDigitsArray->Clear();

    if(fStripMatching) {
        fBmnGemStripDigitMatchesArray->Clear();
    }

    if (!fBmnGemStripPointsArray) {
        Error("BmnGemStripDigitizer::Exec()", " !!! Unknown branch name !!! ");
        return;
    }

    if(fVerbose) {
        cout << "\n BmnGemStripDigitizer::Exec(), event = " << entrys << "\n";
        cout << " BmnGemStripDigitizer::Exec(), Number of BmnGemStripPoints = " << fBmnGemStripPointsArray->GetEntriesFast() << "\n";
    }

    ProcessMCPoints();

    if(fVerbose) cout << " BmnGemStripDigitizer::Exec() finished\n\n";
    entrys++;
    clock_t tFinish = clock();
    workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
}

void BmnGemStripDigitizer::ProcessMCPoints() {

    FairMCPoint* GemStripPoint;
    Int_t NNotPrimaries = 0;

    BmnGemStripStationSet *StationSet = 0;
    switch (fCurrentConfig) {
        case BmnGemStripConfiguration::RunSummer2016:
            StationSet = new BmnGemStripStationSet_RunSummer2016();
            if(fVerbose) cout << "   Current Configuration : RunSummer2016" << "\n";
            break;
        default:
            StationSet = 0;
    }

    for(UInt_t ipoint = 0; ipoint < fBmnGemStripPointsArray->GetEntriesFast(); ipoint++) {
        GemStripPoint = (FairMCPoint*) fBmnGemStripPointsArray->At(ipoint);

        if(fOnlyPrimary) {
            CbmMCTrack *track = (CbmMCTrack*)fMCTracksArray->UncheckedAt(GemStripPoint->GetTrackID());
            if(!track) continue;
            if(track->GetMotherId() != -1)  {
                NNotPrimaries++;
                continue;
            }
        }

        Double_t x = -GemStripPoint->GetX(); // invert because in current geometry +x -left, -x - right
        Double_t y = GemStripPoint->GetY();
        Double_t z = GemStripPoint->GetZ();

        Double_t px = -GemStripPoint->GetPx(); // invert because in current geometry +x -left, -x - right
        Double_t py = GemStripPoint->GetPy();
        Double_t pz = GemStripPoint->GetPz();

        Double_t dEloss = GemStripPoint->GetEnergyLoss()*1e6; // in keV
        Int_t refId = ipoint;

        StationSet->AddPointToDetector(x, y, z, px, py, pz, dEloss, refId);
    }

    Int_t NAddedPoints = StationSet->CountNAddedToDetectorPoints();
    if(fVerbose && fOnlyPrimary) cout << "   Number of not primaries points : " << NNotPrimaries << "\n";
    if(fVerbose) cout << "   Processed MC points  : " << NAddedPoints << "\n";

    for(Int_t iStation = 0; iStation < StationSet->GetNStations(); ++iStation) {
        BmnGemStripStation *station = StationSet->GetGemStation(iStation);

        for(Int_t iModule = 0; iModule < station->GetNModules(); ++iModule) {
            BmnGemStripModule *module = station->GetModule(iModule);

            for(Int_t iLayer = 0; iLayer < module->GetNStripLayers(); ++iLayer) {
                BmnGemStripLayer layer = module->GetStripLayer(iLayer);

                for(Int_t iStrip = 0; iStrip < layer.GetNStrips(); ++iStrip) {
                    new ((*fBmnGemStripDigitsArray)[fBmnGemStripDigitsArray->GetEntriesFast()])
                    BmnGemStripDigit(iStation, iModule, iLayer, iStrip, layer.GetStripSignal(iStrip));

                    if(fStripMatching) {
                        new ((*fBmnGemStripDigitMatchesArray)[fBmnGemStripDigitMatchesArray->GetEntriesFast()])
                            BmnMatch(layer.GetStripMatch(iStrip));
                    }
                }
            }
        }
    }
    if(StationSet) delete StationSet;
}

void BmnGemStripDigitizer::Finish() {
    cout << "Work time of the GEM digitizer: " << workTime << endl;
}

ClassImp(BmnGemStripDigitizer)
