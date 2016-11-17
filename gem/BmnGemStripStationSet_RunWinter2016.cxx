#include "BmnGemStripStationSet_RunWinter2016.h"
#include "BmnGemStripStation_RunWinter2016.h"

BmnGemStripStationSet_RunWinter2016::BmnGemStripStationSet_RunWinter2016(BmnGemStripConfiguration::GEM_CONFIG config) {

    NStations = 6;

    // !!! Classical coordinate system is used !!!
    XStationPositions = new Double_t[NStations];
    YStationPositions = new Double_t[NStations];
    ZStationPositions = new Double_t[NStations];

    fCurrentConfig = config;

    switch(fCurrentConfig) {
        case BmnGemStripConfiguration::RunWinter2016:
            for (Int_t iStat = 0; iStat < NStations; iStat++) {
                XStationPositions[iStat] = -BmnGemStationPositions_RunWinter2016::XStationPositions[iStat]; //inverted : (bm@n x-coord -> classical x-coord)
                YStationPositions[iStat] = BmnGemStationPositions_RunWinter2016::YStationPositions[iStat];
                ZStationPositions[iStat] = BmnGemStationPositions_RunWinter2016::ZStationPositions[iStat];
            }
            break;
            
        default:
            for (Int_t iStat = 0; iStat < NStations; iStat++) {
                XStationPositions[iStat] = -BmnGemStationPositions_RunWinter2016::XStationPositions[iStat]; //inverted : (bm@n x-coord -> classical x-coord)
                YStationPositions[iStat] = BmnGemStationPositions_RunWinter2016::YStationPositions[iStat];
                ZStationPositions[iStat] = BmnGemStationPositions_RunWinter2016::ZStationPositions[iStat];
            }
            break;
    }

    BeamHoleRadiuses = new Double_t[NStations];

    DefineBeamHoleRadiuses();

    BuildStations();
}

BmnGemStripStationSet_RunWinter2016::~BmnGemStripStationSet_RunWinter2016() {

    if (XStationPositions) {
        delete [] XStationPositions;
        XStationPositions = NULL;
    }
    if (YStationPositions) {
        delete [] YStationPositions;
        YStationPositions = NULL;
    }
    if (ZStationPositions) {
        delete [] ZStationPositions;
        ZStationPositions = NULL;
    }
    if (BeamHoleRadiuses) {
        delete [] BeamHoleRadiuses;
        BeamHoleRadiuses = NULL;
    }

    for (Int_t i = 0; i < NStations; i++) {
        if (GemStations[i]) {
            delete GemStations[i];
            GemStations[i] = NULL;
        }
    }
    if (GemStations) {
        delete [] GemStations;
        GemStations = NULL;
    }

}

void BmnGemStripStationSet_RunWinter2016::DefineBeamHoleRadiuses() {
    for (UInt_t iStation = 0; iStation < NStations; iStation++) {
        BeamHoleRadiuses[iStation] = 0.0;
    }
    BeamHoleRadiuses[6] = 4.0; // real hole (without a frame) in the plane (163x45)
}

void BmnGemStripStationSet_RunWinter2016::BuildStations() {
    GemStations = new BmnGemStripStation* [NStations];

    for (Int_t iStation = 0; iStation < NStations; iStation++) {
        GemStations[iStation] =
                new BmnGemStripStation_RunWinter2016(iStation,
                XStationPositions[iStation], YStationPositions[iStation], ZStationPositions[iStation],
                BeamHoleRadiuses[iStation], fCurrentConfig);
    }
}

ClassImp(BmnGemStripStationSet_RunWinter2016)