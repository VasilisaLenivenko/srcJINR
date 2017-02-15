/**
 * @brief 
 * Macro for creating geometry root-file for MWPC
 * @author Sergey Merts
 * @date 15.02.2017
 * 
 */

#include "TGeoManager.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"
//#include "BmnMwpcGeometry.h"

using namespace TMath;

TGeoManager* gGeoMan = NULL;

void create_rootgeom_MWPC_RunWinter2016() {
    
    // Load the necessary FairRoot libraries 
    gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
    bmnloadlibs(); // load libraries

    BmnMwpcGeometry* mwpcGeo = new BmnMwpcGeometry();
    
    //Number of active planes with wires
    const Int_t NofPlanes = mwpcGeo->GetNPlanes();

    //Detector's position
    const Double_t MWPC0_Xpos = mwpcGeo->GetChamberCenter(0).X();
    const Double_t MWPC0_Ypos = mwpcGeo->GetChamberCenter(0).Y();
    const Double_t MWPC0_Zpos = mwpcGeo->GetChamberCenter(0).Z();
    
    const Double_t MWPC1_Xpos = mwpcGeo->GetChamberCenter(1).X();
    const Double_t MWPC1_Ypos = mwpcGeo->GetChamberCenter(1).Y();
    const Double_t MWPC1_Zpos = mwpcGeo->GetChamberCenter(1).Z();

    //           2   
    //        A______B
    //       / \     |\
    //   3  /   \    | \  1
    //     /     \   |  \   
    //     \      \  |  /   
    //   4  \      \ | /  0
    //       \______\|/
    //                C
    //           5   


    //                        X     
    //     _____ <----1----> _____ <-2-> _____
    //    |     |     |     |     |     |     |
    //    |     |     |     |     |     |     |
    //    |     |     |     |     |     |     |
    //    |     |     |     |     |     |     |
    //    |_____|_____|_____|_____|_____|_____|
    //                        X'
    //
    //   1 - GapWireOnePlane
    //   2 - GapWireTwoPlanes


    //Detector's construct parameters   
    const Double_t ABsize = mwpcGeo->GetPlaneWidth(); //cm
    const Double_t BCsize = mwpcGeo->GetPlaneHeight(); //cm
    const Double_t GapZsize = mwpcGeo->GetPlaneStep(); //cm

    const Double_t XSizeOfActiveVolume = ABsize;
    const Double_t YSizeOfActiveVolume = BCsize;
    const Double_t ZSizeOfActiveVolume = GapZsize * NofPlanes; // cm

    //FIXME!!! make sure! sizes of frame  
    const Double_t XWidthChamber = XSizeOfActiveVolume + 1.0;
    const Double_t YHeightChamber = YSizeOfActiveVolume + 1.0;
    const Double_t ZLengthChamber = ZSizeOfActiveVolume + 1.0;

    // -------   Load media from media file   -----------------------------------
    FairGeoLoader* geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    TString geoPath = TString(gSystem->Getenv("VMCWORKDIR")) + "/geometry/";
    TString medFile = geoPath + "media.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    gGeoMan = gGeoManager;
    // --------------------------------------------------------------------------

    // -------   Geometry file name (output)   ----------------------------------
    const TString geoDetectorName = "MWPC";
    const TString geoDetectorVersion = "RunWinter2016";
    TString geoFileName = geoPath + geoDetectorName + "_" + geoDetectorVersion + ".root";
    // --------------------------------------------------------------------------  

    // -----------------   Get and create the required media    -----------------
    FairGeoMedia* geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    FairGeoMedium* mAir = geoMedia->getMedium("air");
    if (!mAir) Fatal("Main", "FairMedium air not found");
    geoBuild->createMedium(mAir);
    TGeoMedium* pMedAir = gGeoMan->GetMedium("air");
    if (!pMedAir) Fatal("Main", "Medium air not found");

    FairGeoMedium* mArCO27030 = geoMedia->getMedium("arco27030");
    if (!mArCO27030) Fatal("Main", "FairMedium arco27030 not found");
    geoBuild->createMedium(mArCO27030);
    TGeoMedium* pMedArCO27030 = gGeoMan->GetMedium("arco27030");
    if (!pMedArCO27030) Fatal("Main", "Medium arco27030 not found");

    // --------------------------------------------------------------------------

    // --------------   Create geometry and top volume  -------------------------
    gGeoMan = (TGeoManager*) gROOT->FindObject("FAIRGeom");
    gGeoMan->SetName(geoDetectorName + "_geom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    top->SetMedium(pMedAir);
    gGeoMan->SetTopVolume(top);
    //gGeoMan->SetTopVisible(1);
    // --------------------------------------------------------------------------

    // Define TOP Geometry
    TGeoVolume* MWPC0Top = new TGeoVolumeAssembly(geoDetectorName);
    MWPC0Top->SetMedium(pMedAir);
    MWPC0Top->SetTransparency(50);
    
    TGeoVolume* MWPC1Top = new TGeoVolumeAssembly(geoDetectorName);
    MWPC1Top->SetMedium(pMedAir);
    MWPC1Top->SetTransparency(50);

    //Transformations (translations, rotations and scales)
    TGeoTranslation *DetPos0_trans = new TGeoTranslation("DetPos0_trans", MWPC0_Xpos, MWPC0_Ypos, MWPC0_Zpos);
    TGeoTranslation *DetPos1_trans = new TGeoTranslation("DetPos1_trans", MWPC1_Xpos, MWPC1_Ypos, MWPC1_Zpos);

    //Solids (shapes)  
    //hexagon which contains active wire planes 
    TGeoPgon *MWPCContainerS = new TGeoPgon("MWPCContainerS", 0, 360, NofPlanes, 2);
    MWPCContainerS->DefineSection(0, -ZLengthChamber / 2.0, 0.0, YSizeOfActiveVolume / 2.0 + 1.0);
    MWPCContainerS->DefineSection(1, ZLengthChamber / 2.0, 0.0, YSizeOfActiveVolume / 2.0 + 1.0);
    //active wire plane
    TGeoBBox *MWPCActivePlaneS = new TGeoBBox("MWPCActivePlaneS", XSizeOfActiveVolume / 2.0, YSizeOfActiveVolume / 2.0, 0.1/*GapZsize / 2.0*/);

    //Volumes
    TGeoVolume *MWPC0ContainerV = new TGeoVolume("MWPC0ContainerV", MWPCContainerS);
    MWPC0ContainerV->SetMedium(pMedAir);
    MWPC0ContainerV->SetVisibility(kTRUE);
    MWPC0ContainerV->SetTransparency(60);

    TGeoVolume *MWPC0ActivePlaneV = new TGeoVolume("MWPC0ActivePlaneV", MWPCActivePlaneS);
    MWPC0ActivePlaneV->SetMedium(pMedArCO27030);
    MWPC0ActivePlaneV->SetLineColor(kBlue);
    MWPC0ActivePlaneV->SetTransparency(40);
    
    TGeoVolume *MWPC1ContainerV = new TGeoVolume("MWPC1ContainerV", MWPCContainerS);
    MWPC1ContainerV->SetMedium(pMedAir);
    MWPC1ContainerV->SetVisibility(kTRUE);
    MWPC1ContainerV->SetTransparency(60);

    TGeoVolume *MWPC1ActivePlaneV = new TGeoVolume("MWPC1ActivePlaneV", MWPCActivePlaneS);
    MWPC1ActivePlaneV->SetMedium(pMedArCO27030);
    MWPC1ActivePlaneV->SetLineColor(kBlue);
    MWPC1ActivePlaneV->SetTransparency(40);

    //Adding volumes to the TOP Volume
    top->AddNode(MWPC0Top, 1, DetPos0_trans);
    MWPC0Top->AddNode(MWPC0ContainerV, 1);
    for (Int_t iPlane = 1; iPlane <= NofPlanes; ++iPlane) {
        TGeoTranslation t1(0.0, 0.0, GapZsize * (iPlane - (NofPlanes + 1) / 2.0));
        TGeoRotation r1("r1", 0.0, 0.0, (iPlane - 1) * 60.0);
        MWPC0ContainerV->AddNode(MWPC0ActivePlaneV, iPlane, new TGeoCombiTrans(t1, r1));
    }
    
    top->AddNode(MWPC1Top, 1, DetPos1_trans);
    MWPC1Top->AddNode(MWPC1ContainerV, 1);
    for (Int_t iPlane = 1; iPlane <= NofPlanes; ++iPlane) {
        TGeoTranslation t1(0.0, 0.0, GapZsize * (iPlane - (NofPlanes + 1) / 2.0));
        TGeoRotation r1("r1", 0.0, 0.0, (iPlane - 1) * 60.0);
        MWPC1ContainerV->AddNode(MWPC1ActivePlaneV, iPlane, new TGeoCombiTrans(t1, r1));
    }

    top->SetVisContainers(kTRUE);

    // ---------------   Finish   -----------------------------------------------
    gGeoMan->CloseGeometry();
    gGeoMan->CheckOverlaps(0.001);
    gGeoMan->PrintOverlaps();
    gGeoMan->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();
    top->Draw("ogl");
}