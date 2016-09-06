#ifndef BMNGEMSTRIPCONFIGURATION_H
#define	BMNGEMSTRIPCONFIGURATION_H

namespace BmnGemStripConfiguration {
    enum GEM_CONFIG {
        None = 0, // empty config
        RunSummer2016 // 7 stations (summer 2016) - (geometry file: GEMS_RunSummer2016.root)
    };
}

namespace BmnGemStripPositions_RunSummer2016 {
    //BM@N coordinate system is used
    const Int_t NStations = 7;
    const Double_t XStationPositions[NStations] = {+0.0, -1.5, +1.5, -1.5, +1.5, -1.5, +0.0};
    const Double_t YStationPositions[NStations] = {+0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0};
    const Double_t ZStationPositions[NStations] = {5.0+0.1, 41.8+0.7, 62.8+0.7, 83.8+0.7, 104.8+0.7, 125.8+0.7, 146.8+1.5}; //All z-positions are start pos. of the sensitive vol.
}

#endif