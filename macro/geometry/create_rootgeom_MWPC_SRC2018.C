/**
 * @brief 
 * Macro for creating geometry root-file for MWPC
 * @author Maria Patsyuk
 * @date 15.02.2018
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

void create_rootgeom_MWPC_SRC2018() {
    
    // Load the necessary FairRoot libraries 
    gROOT->LoadMacro("$VMCWORKDIR/macro/run/bmnloadlibs.C");
    bmnloadlibs(); // load libraries

    BmnMwpcGeometrySRC* mwpcGeo = new BmnMwpcGeometrySRC();
    
    //Number of active planes with wires
    const Int_t NofPlanes = mwpcGeo->GetNPlanes();

    //Number of chambers used:
    Short_t Ncham = mwpcGeo->GetNChambers();
    
    //Detector's position and angles
    Double_t MWPC_Xpos[4];
    Double_t MWPC_Ypos[4];
    Double_t MWPC_Zpos[4];
    TVector3 OXprime[4];
    TVector3 OYprime[4];
    TVector3 OZprime[4];
    
    for(Int_t nCham =0; nCham<Ncham; nCham++){
      MWPC_Xpos[nCham] = mwpcGeo->GetChamberCenter(nCham).X();
      MWPC_Ypos[nCham] = mwpcGeo->GetChamberCenter(nCham).Y();
      MWPC_Zpos[nCham] = mwpcGeo->GetChamberCenter(nCham).Z();
      OXprime[nCham] = mwpcGeo->GetAxisPrime(nCham,0);
      OYprime[nCham] = mwpcGeo->GetAxisPrime(nCham,1);
      OZprime[nCham] = mwpcGeo->GetAxisPrime(nCham,2);
    }
    /*    
    cout<<"0 chamber center positions: x = "<<MWPC_Xpos[0]<<", y = "<<MWPC_Ypos[0]<<", z = "<<MWPC_Zpos[0]<<endl;
    cout<<"1 chamber center positions: x = "<<MWPC_Xpos[1]<<", y = "<<MWPC_Ypos[1]<<", z = "<<MWPC_Zpos[1]<<endl;
    cout<<"2 chamber center positions: x = "<<MWPC_Xpos[2]<<", y = "<<MWPC_Ypos[2]<<", z = "<<MWPC_Zpos[2]<<endl;
    cout<<"3 chamber center positions: x = "<<MWPC_Xpos[3]<<", y = "<<MWPC_Ypos[3]<<", z = "<<MWPC_Zpos[3]<<endl;
    cout<<" no rotation: 90 0 90 90 0 0"<<endl;
    cout<<"0 chamber ox prime: theta = "<<OXprime[0].Theta()/TMath::Pi()*180<<", phi = "<<OXprime[0].Phi()/TMath::Pi()*180<<endl;
    cout<<"0 chamber oy prime: theta = "<<OYprime[0].Theta()/TMath::Pi()*180<<", phi = "<<OYprime[0].Phi()/TMath::Pi()*180<<endl;
    cout<<"0 chamber oz prime: theta = "<<OZprime[0].Theta()/TMath::Pi()*180<<", phi = "<<OZprime[0].Phi()/TMath::Pi()*180<<endl;
    cout<<"1 chamber ox prime: theta = "<<OXprime[1].Theta()/TMath::Pi()*180<<", phi = "<<OXprime[1].Phi()/TMath::Pi()*180<<endl;
    cout<<"1 chamber oy prime: theta = "<<OYprime[1].Theta()/TMath::Pi()*180<<", phi = "<<OYprime[1].Phi()/TMath::Pi()*180<<endl;
    cout<<"1 chamber oz prime: theta = "<<OZprime[1].Theta()/TMath::Pi()*180<<", phi = "<<OZprime[1].Phi()/TMath::Pi()*180<<endl;
    cout<<"2 chamber ox prime: theta = "<<OXprime[2].Theta()/TMath::Pi()*180<<", phi = "<<OXprime[2].Phi()/TMath::Pi()*180<<endl;
    cout<<"2 chamber oy prime: theta = "<<OYprime[2].Theta()/TMath::Pi()*180<<", phi = "<<OYprime[2].Phi()/TMath::Pi()*180<<endl;
    cout<<"2 chamber oz prime: theta = "<<OZprime[2].Theta()/TMath::Pi()*180<<", phi = "<<OZprime[2].Phi()/TMath::Pi()*180<<endl; 
    cout<<"3 chamber ox prime: theta = "<<OXprime[3].Theta()/TMath::Pi()*180<<", phi = "<<OXprime[3].Phi()/TMath::Pi()*180<<endl;
    cout<<"3 chamber ox prime: theta = "<<OYprime[3].Theta()/TMath::Pi()*180<<", phi = "<<OYprime[3].Phi()/TMath::Pi()*180<<endl;
    cout<<"3 chamber ox prime: theta = "<<OZprime[3].Theta()/TMath::Pi()*180<<", phi = "<<OZprime[3].Phi()/TMath::Pi()*180<<endl;
    */
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

    const Double_t ZLengthChamber = ZSizeOfActiveVolume + 1.;
    
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
    const TString geoDetectorVersion = "SRC2018";
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

    FairGeoMedium* mMWPCgas = geoMedia->getMedium("DCH_MWPC_gas");
    if (!mMWPCgas) Fatal("Main", "FairMedium DCH_MWPC_gas not found");
    geoBuild->createMedium(mMWPCgas);
    TGeoMedium* pMWPCgas = gGeoMan->GetMedium("DCH_MWPC_gas");
    if (!pMWPCgas) Fatal("Main", "Medium DCH_MWPC_gas not found");
    
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
    TGeoVolume* MWPC = new TGeoVolumeAssembly(geoDetectorName);
    MWPC->SetMedium(pMedAir);
    MWPC->SetTransparency(50);
    
    // Transformations (translations, rotations and scales)
    TGeoTranslation DetPos0_trans("DetPos0_trans", MWPC_Xpos[0], MWPC_Ypos[0], MWPC_Zpos[0]);
    TGeoTranslation DetPos1_trans("DetPos1_trans", MWPC_Xpos[1], MWPC_Ypos[1], MWPC_Zpos[1]);
    TGeoTranslation DetPos2_trans("DetPos2_trans", MWPC_Xpos[2], MWPC_Ypos[2], MWPC_Zpos[2]);
    TGeoTranslation DetPos3_trans("DetPos3_trans", MWPC_Xpos[3], MWPC_Ypos[3], MWPC_Zpos[3]);

    // check that phi angles are in (0, 360):
    Double_t phiX[4], phiY[4], phiZ[4];
    for(Int_t i=0; i<4; i++){
      phiX[i] = 45.;//OXprime[i].Phi()/TMath::Pi()*180.;
      phiY[i] = 120.;//OYprime[i].Phi()/TMath::Pi()*180.;
      phiZ[i] = OZprime[i].Phi()/TMath::Pi()*180.;
      if(phiX[i] < 0.)  { phiX[i] = 360. - phiX[i];}
      if(phiX[i] > 360.){ phiX[i] = phiX[i] - 360.;}
      if(phiY[i] < 0.)  { phiY[i] = 360. - phiY[i];}
      if(phiY[i] > 360.){ phiY[i] = phiY[i] - 360.;}
      if(phiZ[i] < 0.)  { phiZ[i] = 360. - phiZ[i];}
      if(phiZ[i] > 360.){ phiZ[i] = phiZ[i] - 360.;}
    }
    
    // Rotations of each chamber:
    TGeoRotation Det0_rot("Det0_rot", OXprime[0].Theta()/TMath::Pi()*180., phiX[0],
			  OYprime[0].Theta()/TMath::Pi()*180., phiY[0],
			  OZprime[0].Theta()/TMath::Pi()*180., phiZ[0]);
    TGeoRotation Det1_rot("Det1_rot", OXprime[1].Theta()/TMath::Pi()*180., phiX[1],
			  OYprime[1].Theta()/TMath::Pi()*180., phiY[1],
			  OZprime[1].Theta()/TMath::Pi()*180., phiZ[1]);
    TGeoRotation Det2_rot("Det2_rot", OXprime[2].Theta()/TMath::Pi()*180., phiX[2],
			  OYprime[2].Theta()/TMath::Pi()*180., phiY[2],
			  OZprime[2].Theta()/TMath::Pi()*180., phiZ[2]);
    TGeoRotation Det3_rot("Det3_rot", OXprime[3].Theta()/TMath::Pi()*180., phiX[3],
			  OYprime[3].Theta()/TMath::Pi()*180., phiY[3],
			  OZprime[3].Theta()/TMath::Pi()*180., phiZ[3]);
    
    //Solids (shapes)  
    //hexagon which contains active wire planes 
    TGeoPgon *MWPCContainerS = new TGeoPgon("MWPCContainerS", 0, 360, NofPlanes, 2);
    MWPCContainerS->DefineSection(0, -ZLengthChamber / 2.0, 0.0, YSizeOfActiveVolume / 2.0 + 1.0);
    MWPCContainerS->DefineSection(1, ZLengthChamber / 2.0, 0.0, YSizeOfActiveVolume / 2.0 + 1.0);
    
    //active wire plane
    TGeoBBox *MWPCActivePlaneS = new TGeoBBox("MWPCActivePlaneS", XSizeOfActiveVolume / 2.0, YSizeOfActiveVolume / 2.0, 0.1);

    //Volumes
    TGeoVolume *MWPCContainerV = new TGeoVolume("MWPCContainerV", MWPCContainerS);
    MWPCContainerV->SetMedium(pMedAir);
    MWPCContainerV->SetVisibility(kTRUE);
    MWPCContainerV->SetTransparency(75);

    TGeoVolume *MWPCActivePlaneV = new TGeoVolume("MWPCActivePlaneV", MWPCActivePlaneS);
    MWPCActivePlaneV->SetMedium(pMWPCgas);//(pMedArCO27030);
    MWPCActivePlaneV->SetLineColor(kBlue);
    MWPCActivePlaneV->SetTransparency(40);
    
    for (Int_t iPlane = 1; iPlane <= NofPlanes; ++iPlane) {
        TGeoTranslation t0(0.0, 0.0, GapZsize * (iPlane - (NofPlanes + 1) / 2.0));
        TGeoRotation r0("r0", 0.0, 0.0, (iPlane - 1) * 60.0);
        MWPCContainerV->AddNode(MWPCActivePlaneV, iPlane - 1, new TGeoCombiTrans(t0, r0));
    }
    
    MWPC->AddNode(MWPCContainerV, 0, new TGeoCombiTrans(DetPos0_trans, 0));//Det0_rot));
    MWPC->AddNode(MWPCContainerV, 1, new TGeoCombiTrans(DetPos1_trans, 0));//Det1_rot));
    MWPC->AddNode(MWPCContainerV, 2, new TGeoCombiTrans(DetPos2_trans, 0));//Det2_rot));
    MWPC->AddNode(MWPCContainerV, 3, new TGeoCombiTrans(DetPos3_trans, 0));//Det3_rot));
    
    //Adding volumes to the TOP Volume
    top->AddNode(MWPC, 0);
    
    top->SetVisContainers(kTRUE);

    // ---------------   Finish   -----------------------------------------------
    gGeoMan->CloseGeometry();
    gGeoMan->CheckOverlaps(0.001);
    gGeoMan->PrintOverlaps();
    gGeoMan->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();
    //top->Draw("ogl");
}
