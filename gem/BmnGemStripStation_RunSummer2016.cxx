#include "BmnGemStripStation_RunSummer2016.h"

BmnGemStripStation_RunSummer2016::BmnGemStripStation_RunSummer2016(Int_t iStation,
                       Double_t xpos_station, Double_t ypos_station, Double_t zpos_station,
                       Double_t beamradius) {

    //All linear values must be entered in units of cm

    StationNumber = iStation;

    //Station shift position
    XPosition = xpos_station;
    YPosition = ypos_station;
    ZPosition = zpos_station;

    BeamHoleRadius = beamradius;


    //Parameters for each station
        //ZSizeModuleWithFrames : GEM thickness (frames+sensitive)
        //ZSizeModule : sensitive volume thickness
        //XSize : x-size of the station
        //YSize : y-size of the station
        //ZSize : z-size of the station

    switch (iStation) {
        case 0: //one plane 10x10
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = +90.0;
            LowerLayerPitch = 0.04;
            UpperLayerPitch = 0.04;
            ZSizeModuleWithFrames = 1.3;
            ZSizeModule = 0.9;
            EDriftDirection = ForwardZAxisEDrift;
            break;

        case 1: //one plane 66x41
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = +15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 2.3;
            ZSizeModule = 0.9;
            EDriftDirection = ForwardZAxisEDrift;
            break;

        case 2: //one plane 66x41
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = -15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 2.3;
            ZSizeModule = 0.9;
            EDriftDirection = BackwardZAxisEDrift;
            break;

        case 3: //one plane 66x41
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = +15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 2.3;
            ZSizeModule = 0.9;
            EDriftDirection = ForwardZAxisEDrift;
            break;

        case 4: //one plane 66x41
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = -15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 2.3;
            ZSizeModule = 0.9;
            EDriftDirection = BackwardZAxisEDrift;
            break;

        case 5: //one plane 66x41
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = +15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 2.3;
            ZSizeModule = 0.9;
            EDriftDirection = ForwardZAxisEDrift;
            break;

        case 6: //one plane 163x45
            LowerLayerStripAngle = 0.0;
            UpperLayerStripAngle = -15.0;
            LowerLayerPitch = 0.08;
            UpperLayerPitch = 0.08;
            ZSizeModuleWithFrames = 3.9;
            ZSizeModule = 0.9;
            EDriftDirection = BackwardZAxisEDrift;
            break;
    }

    //Planes SIZES

    //station 0 (plane 10x10)
    XModuleSize_Plane10x10 = 10.0;
    YModuleSize_Plane10x10 = 10.0;

    //station 1-5 (plane 66x41)
    XModuleSize_Plane66x41 = 66.0;
    YModuleSize_Plane66x41 = 41.0;

    XHotZoneSize_Plane66x41 = 15.0;
    YHotZoneSize_Plane66x41 = 10.0;

    //station 6 (plane 163x45)
    XModuleSize_Plane163x45 = 163.2*0.5;
    YModuleSize_Plane163x45 = 45.0;

    XRectHotZoneSize_Plane163x45 = 40.0;
    YRectHotZoneSize_Plane163x45 = 15.0;

    XSlopeHotZoneSize_Plane163x45[0] = 37.99;
    XSlopeHotZoneSize_Plane163x45[1] = 42.01;

    YSlopeHotZoneSize_Plane163x45 = 15.0;


    //Assembling a station -----------------------------------------------------
    if(StationNumber == 0) {
        NModules = 1; //1 modules in the station
        XSize = XModuleSize_Plane10x10;
        YSize = YModuleSize_Plane10x10;
        ZSize = ZSizeModuleWithFrames;
        BuildModules_One10x10Plane();
    }

    if(StationNumber >= 1 && StationNumber <= 5) {
        NModules = 1; //1 modules in the station
        XSize = XModuleSize_Plane66x41;
        YSize = YModuleSize_Plane66x41;
        ZSize = ZSizeModuleWithFrames;
        BuildModules_One66x41Plane();
    }

    if(StationNumber == 6) {
        NModules = 2; //2 modules in the station (left+right)
        XSize = 2*XModuleSize_Plane163x45;
        YSize = YModuleSize_Plane163x45;
        ZSize = ZSizeModuleWithFrames;
        BuildModules_One163x45Plane();
    }
    //end assembling the station -----------------------------------------------
}

BmnGemStripStation_RunSummer2016::~BmnGemStripStation_RunSummer2016() {
    for(Int_t i = 0; i < NModules; ++i) {
        if(Modules[i]) {
            delete Modules[i];
            Modules[i] = NULL;
        }
    }
    if(Modules) {
        delete [] Modules;
        Modules = NULL;
    }
}

Int_t BmnGemStripStation_RunSummer2016::GetPointModuleOwnership(Double_t xcoord, Double_t ycoord, Double_t zcoord) {

    for(Int_t imodule = 0; imodule < NModules; ++imodule) {
        if( Modules[imodule]->IsPointInsideModule(xcoord, ycoord) ) return imodule;
    }

    return -1;
}
//------------------------------------------------------------------------------

void BmnGemStripStation_RunSummer2016::BuildModules_One10x10Plane() {

    Modules = new BmnGemStripModule* [NModules];

    //module
    Modules[0] = new BmnGemStripModule(ZPosition, EDriftDirection);

        //lower strip layer --------------------------------------------------------
        BmnGemStripLayer lower_layer(0, LowerStripLayer,
                                     XModuleSize_Plane10x10, YModuleSize_Plane10x10,
                                     XPosition-XModuleSize_Plane10x10*0.5, YPosition-YModuleSize_Plane10x10*0.5,
                                     LowerLayerPitch, LowerLayerStripAngle);

        lower_layer.SetStripNumberingOrder(RightToLeft);
        lower_layer.SetStripNumberingBorders(LeftTop, RightBottom);
        //--------------------------------------------------------------------------

        //upper strip layer --------------------------------------------------------
        BmnGemStripLayer upper_layer(0, UpperStripLayer,
                                     XModuleSize_Plane10x10, YModuleSize_Plane10x10,
                                     XPosition-XModuleSize_Plane10x10*0.5, YPosition-YModuleSize_Plane10x10*0.5,
                                     UpperLayerPitch, UpperLayerStripAngle);

        upper_layer.SetStripNumberingOrder(LeftToRight);
        upper_layer.SetStripNumberingBorders(LeftTop, RightBottom);
        //--------------------------------------------------------------------------

    Modules[0]->AddStripLayer(lower_layer);
    Modules[0]->AddStripLayer(upper_layer);

    return;
}

void BmnGemStripStation_RunSummer2016::BuildModules_One66x41Plane() {

    Modules = new BmnGemStripModule* [NModules];

    //module
    Modules[0] = new BmnGemStripModule(ZPosition, EDriftDirection);

    //zone 0 (big) -------------------------------------------------------------

        //lower strip layer (zone 0) -------------------------------------------
        BmnGemStripLayer lower_layer_zone0(0, LowerStripLayer,
                                           XModuleSize_Plane66x41, YModuleSize_Plane66x41,
                                           XPosition-XModuleSize_Plane66x41*0.5, YPosition-YModuleSize_Plane66x41*0.5,
                                           LowerLayerPitch, LowerLayerStripAngle);

        if(UpperLayerStripAngle >= 0.0) {
            lower_layer_zone0.SetStripNumberingOrder(LeftToRight);
            lower_layer_zone0.SetStripNumberingBorders(LeftTop, RightBottom);
        }
        else {
            lower_layer_zone0.SetStripNumberingOrder(LeftToRight);
            lower_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);
        }

        //dead zone for layers in zone 0
        const Int_t NPoints_DeadZone_zone0 = 4;
        Double_t XPoints_DeadZone_zone0[NPoints_DeadZone_zone0];
        Double_t YPoints_DeadZone_zone0[NPoints_DeadZone_zone0];

        if(UpperLayerStripAngle >= 0.0) {
            XPoints_DeadZone_zone0[0] = XPosition-XModuleSize_Plane66x41*0.5;
            XPoints_DeadZone_zone0[1] = XPosition-XModuleSize_Plane66x41*0.5;
            XPoints_DeadZone_zone0[2] = XPosition-XModuleSize_Plane66x41*0.5 + XHotZoneSize_Plane66x41;
            XPoints_DeadZone_zone0[3] = XPosition-XModuleSize_Plane66x41*0.5 + XHotZoneSize_Plane66x41;

            YPoints_DeadZone_zone0[0] = YPosition-YModuleSize_Plane66x41*0.5;
            YPoints_DeadZone_zone0[1] = YPosition-YModuleSize_Plane66x41*0.5 + YHotZoneSize_Plane66x41;
            YPoints_DeadZone_zone0[2] = YPosition-YModuleSize_Plane66x41*0.5 + YHotZoneSize_Plane66x41;
            YPoints_DeadZone_zone0[3] = YPosition-YModuleSize_Plane66x41*0.5;
        }
        else {
            XPoints_DeadZone_zone0[0] = XPosition-XModuleSize_Plane66x41*0.5;
            XPoints_DeadZone_zone0[1] = XPosition-XModuleSize_Plane66x41*0.5;
            XPoints_DeadZone_zone0[2] = XPosition-XModuleSize_Plane66x41*0.5 + XHotZoneSize_Plane66x41;
            XPoints_DeadZone_zone0[3] = XPosition-XModuleSize_Plane66x41*0.5 + XHotZoneSize_Plane66x41;

            YPoints_DeadZone_zone0[0] = YPosition+YModuleSize_Plane66x41*0.5 - YHotZoneSize_Plane66x41;
            YPoints_DeadZone_zone0[1] = YPosition+YModuleSize_Plane66x41*0.5;
            YPoints_DeadZone_zone0[2] = YPosition+YModuleSize_Plane66x41*0.5;
            YPoints_DeadZone_zone0[3] = YPosition+YModuleSize_Plane66x41*0.5 - YHotZoneSize_Plane66x41;
        }

        lower_layer_zone0.AddDeadZone(NPoints_DeadZone_zone0, XPoints_DeadZone_zone0, YPoints_DeadZone_zone0);
        //----------------------------------------------------------------------

        //upper strip layer (zone 0) -------------------------------------------
        BmnGemStripLayer upper_layer_zone0(0, UpperStripLayer,
                                           XModuleSize_Plane66x41, YModuleSize_Plane66x41,
                                           XPosition-XModuleSize_Plane66x41*0.5, YPosition-YModuleSize_Plane66x41*0.5,
                                           UpperLayerPitch, UpperLayerStripAngle);

        if(UpperLayerStripAngle >= 0.0) {
            upper_layer_zone0.SetStripNumberingOrder(LeftToRight);
            upper_layer_zone0.SetStripNumberingBorders(LeftTop, RightBottom);
        }
        else {
            upper_layer_zone0.SetStripNumberingOrder(LeftToRight);
            upper_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);
        }

        upper_layer_zone0.AddDeadZone(NPoints_DeadZone_zone0, XPoints_DeadZone_zone0, YPoints_DeadZone_zone0);
        //----------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //zone 1 (hot) -------------------------------------------------------------

        BmnGemStripLayer lower_layer_zone1;
        BmnGemStripLayer upper_layer_zone1;

        if(UpperLayerStripAngle >= 0.0) {
            //lower strip layer (zone 1) -------------------------------------------
            lower_layer_zone1 = BmnGemStripLayer(1, LowerStripLayer,
                                               XHotZoneSize_Plane66x41, YHotZoneSize_Plane66x41,
                                               XPosition-XModuleSize_Plane66x41*0.5, YPosition-YModuleSize_Plane66x41*0.5,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone1.SetStripNumberingOrder(LeftToRight);
            lower_layer_zone1.SetStripNumberingBorders(LeftTop, RightBottom);
            //----------------------------------------------------------------------

            //upper strip layer (zone 1) -------------------------------------------
            upper_layer_zone1 = BmnGemStripLayer(1, UpperStripLayer,
                                               XHotZoneSize_Plane66x41, YHotZoneSize_Plane66x41,
                                               XPosition-XModuleSize_Plane66x41*0.5, YPosition-YModuleSize_Plane66x41*0.5,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone1.SetStripNumberingOrder(LeftToRight);
            upper_layer_zone1.SetStripNumberingBorders(LeftTop, RightBottom);
            //----------------------------------------------------------------------
        }
        else {
            //lower strip layer (zone 1) -------------------------------------------
            lower_layer_zone1 = BmnGemStripLayer(1, LowerStripLayer,
                                               XHotZoneSize_Plane66x41, YHotZoneSize_Plane66x41,
                                               XPosition-XModuleSize_Plane66x41*0.5, YPosition+YModuleSize_Plane66x41*0.5-YHotZoneSize_Plane66x41,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone1.SetStripNumberingOrder(LeftToRight);
            lower_layer_zone1.SetStripNumberingBorders(LeftTop, RightBottom);
            //----------------------------------------------------------------------

            //upper strip layer (zone 1) -------------------------------------------
            upper_layer_zone1 = BmnGemStripLayer(1, UpperStripLayer,
                                               XHotZoneSize_Plane66x41, YHotZoneSize_Plane66x41,
                                               XPosition-XModuleSize_Plane66x41*0.5, YPosition+YModuleSize_Plane66x41*0.5-YHotZoneSize_Plane66x41,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone1.SetStripNumberingOrder(LeftToRight);
            upper_layer_zone1.SetStripNumberingBorders(LeftBottom, RightTop);
            //----------------------------------------------------------------------
        }
    //--------------------------------------------------------------------------

    Modules[0]->AddStripLayer(lower_layer_zone0);
    Modules[0]->AddStripLayer(upper_layer_zone0);
    Modules[0]->AddStripLayer(lower_layer_zone1);
    Modules[0]->AddStripLayer(upper_layer_zone1);

    return;
}

void BmnGemStripStation_RunSummer2016::BuildModules_One163x45Plane() {
    Modules = new BmnGemStripModule* [NModules];


//module 0 (left) --------------------------------------------------------------
    {
        Modules[0] = new BmnGemStripModule(ZPosition, EDriftDirection);

            //zone 0 (big) ---------------------------------------------------------

            //lower strip layer (zone 0) -------------------------------------------
            BmnGemStripLayer lower_layer_zone0(0, LowerStripLayer,
                                               XModuleSize_Plane163x45, YModuleSize_Plane163x45,
                                               XPosition-XModuleSize_Plane163x45, YPosition-YModuleSize_Plane163x45*0.5,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone0.SetStripNumberingOrder(RightToLeft);
            lower_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for lower layer in zone 0
            const Int_t NPoints_DeadZone_lower_zone0 = 4;
            Double_t XPoints_DeadZone_lower_zone0[NPoints_DeadZone_lower_zone0];
            Double_t YPoints_DeadZone_lower_zone0[NPoints_DeadZone_lower_zone0];

            XPoints_DeadZone_lower_zone0[0] = XPosition-XRectHotZoneSize_Plane163x45;
            XPoints_DeadZone_lower_zone0[1] = XPosition-XRectHotZoneSize_Plane163x45;
            XPoints_DeadZone_lower_zone0[2] = XPosition+0.0;
            XPoints_DeadZone_lower_zone0[3] = XPosition+0.0;

            YPoints_DeadZone_lower_zone0[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_lower_zone0[1] = YPosition-YModuleSize_Plane163x45*0.5 + YRectHotZoneSize_Plane163x45;
            YPoints_DeadZone_lower_zone0[2] = YPosition-YModuleSize_Plane163x45*0.5 + YRectHotZoneSize_Plane163x45;
            YPoints_DeadZone_lower_zone0[3] = YPosition-YModuleSize_Plane163x45*0.5;

            lower_layer_zone0.AddDeadZone(NPoints_DeadZone_lower_zone0, XPoints_DeadZone_lower_zone0, YPoints_DeadZone_lower_zone0);
            //----------------------------------------------------------------------

            //upper strip layer (zone 0) -------------------------------------------
            BmnGemStripLayer upper_layer_zone0(0, UpperStripLayer,
                                               XModuleSize_Plane163x45, YModuleSize_Plane163x45,
                                               XPosition-XModuleSize_Plane163x45, YPosition-YModuleSize_Plane163x45*0.5,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone0.SetStripNumberingOrder(RightToLeft);
            upper_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for upper layer in zone 0
            const Int_t NPoints_DeadZone_upper_zone0 = 4;
            Double_t XPoints_DeadZone_upper_zone0[NPoints_DeadZone_upper_zone0];
            Double_t YPoints_DeadZone_upper_zone0[NPoints_DeadZone_upper_zone0];

            XPoints_DeadZone_upper_zone0[0] = XPosition-XSlopeHotZoneSize_Plane163x45[0];
            XPoints_DeadZone_upper_zone0[1] = XPosition-XSlopeHotZoneSize_Plane163x45[1];
            XPoints_DeadZone_upper_zone0[2] = XPosition+0.0;
            XPoints_DeadZone_upper_zone0[3] = XPosition+0.0;

            YPoints_DeadZone_upper_zone0[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone0[1] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone0[2] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone0[3] = YPosition-YModuleSize_Plane163x45*0.5;

            upper_layer_zone0.AddDeadZone(NPoints_DeadZone_upper_zone0, XPoints_DeadZone_upper_zone0, YPoints_DeadZone_upper_zone0);
            //----------------------------------------------------------------------

            //lower strip layer (zone 1) -------------------------------------------
            BmnGemStripLayer lower_layer_zone1(1, LowerStripLayer,
                                               XRectHotZoneSize_Plane163x45, YRectHotZoneSize_Plane163x45,
                                               XPosition-XRectHotZoneSize_Plane163x45, YPosition-YModuleSize_Plane163x45*0.5,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone1.SetStripNumberingOrder(RightToLeft);
            lower_layer_zone1.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for lower layer in zone 1 (hole)
            const Int_t NPoints_DeadZone_lower_zone1 = 4;
            Double_t XPoints_DeadZone_lower_zone1[NPoints_DeadZone_lower_zone1];
            Double_t YPoints_DeadZone_lower_zone1[NPoints_DeadZone_lower_zone1];

            XPoints_DeadZone_lower_zone1[0] = XPosition-BeamHoleRadius;
            XPoints_DeadZone_lower_zone1[1] = XPosition-BeamHoleRadius;
            XPoints_DeadZone_lower_zone1[2] = XPosition+0.0;
            XPoints_DeadZone_lower_zone1[3] = XPosition+0.0;

            YPoints_DeadZone_lower_zone1[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_lower_zone1[1] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_lower_zone1[2] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_lower_zone1[3] = YPosition-YModuleSize_Plane163x45*0.5;

            lower_layer_zone1.AddDeadZone(NPoints_DeadZone_lower_zone1, XPoints_DeadZone_lower_zone1, YPoints_DeadZone_lower_zone1);
            //----------------------------------------------------------------------

            //upper strip layer (zone 1) -------------------------------------------
            BmnGemStripLayer upper_layer_zone1(1, UpperStripLayer,
                                               XSlopeHotZoneSize_Plane163x45[1], YSlopeHotZoneSize_Plane163x45,
                                               XPosition-XSlopeHotZoneSize_Plane163x45[1], YPosition-YModuleSize_Plane163x45*0.5,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone1.SetStripNumberingOrder(RightToLeft);
            upper_layer_zone1.SetStripNumberingBorders(LeftTop, RightTop);

            //dead zone for upper layer in zone 1 (hole)
            const Int_t NPoints_DeadZone_upper_zone1_hole = 4;
            Double_t XPoints_DeadZone_upper_zone1_hole[NPoints_DeadZone_upper_zone1_hole];
            Double_t YPoints_DeadZone_upper_zone1_hole[NPoints_DeadZone_upper_zone1_hole];

            XPoints_DeadZone_upper_zone1_hole[0] = XPosition-BeamHoleRadius;
            XPoints_DeadZone_upper_zone1_hole[1] = XPosition-BeamHoleRadius;
            XPoints_DeadZone_upper_zone1_hole[2] = XPosition+0.0;
            XPoints_DeadZone_upper_zone1_hole[3] = XPosition+0.0;

            YPoints_DeadZone_upper_zone1_hole[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone1_hole[1] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_upper_zone1_hole[2] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_upper_zone1_hole[3] = YPosition-YModuleSize_Plane163x45*0.5;

            upper_layer_zone1.AddDeadZone(NPoints_DeadZone_upper_zone1_hole, XPoints_DeadZone_upper_zone1_hole, YPoints_DeadZone_upper_zone1_hole);

            //dead zone for upper layer in zone 1 (triangle)
            const Int_t NPoints_DeadZone_upper_zone1_triang = 3;
            Double_t XPoints_DeadZone_upper_zone1_triang[NPoints_DeadZone_upper_zone1_triang];
            Double_t YPoints_DeadZone_upper_zone1_triang[NPoints_DeadZone_upper_zone1_triang];

            XPoints_DeadZone_upper_zone1_triang[0] = XPosition-XSlopeHotZoneSize_Plane163x45[1];
            XPoints_DeadZone_upper_zone1_triang[1] = XPosition-XSlopeHotZoneSize_Plane163x45[1];
            XPoints_DeadZone_upper_zone1_triang[2] = XPosition-XSlopeHotZoneSize_Plane163x45[0];

            YPoints_DeadZone_upper_zone1_triang[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone1_triang[1] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone1_triang[2] = YPosition-YModuleSize_Plane163x45*0.5;

            upper_layer_zone1.AddDeadZone(NPoints_DeadZone_upper_zone1_triang, XPoints_DeadZone_upper_zone1_triang, YPoints_DeadZone_upper_zone1_triang);

            //----------------------------------------------------------------------

            Modules[0]->AddStripLayer(lower_layer_zone0);
            Modules[0]->AddStripLayer(upper_layer_zone0);
            Modules[0]->AddStripLayer(lower_layer_zone1);
            Modules[0]->AddStripLayer(upper_layer_zone1);

    }
//--------------------------------------------------------------------------

//module 1 (right) --------------------------------------------------------------
    {
        Modules[1] = new BmnGemStripModule(ZPosition, EDriftDirection);

            //zone 0 (big) ---------------------------------------------------------

            //lower strip layer (zone 0) -------------------------------------------
            BmnGemStripLayer lower_layer_zone0(0, LowerStripLayer,
                                               XModuleSize_Plane163x45, YModuleSize_Plane163x45,
                                               XPosition+0.0, YPosition-YModuleSize_Plane163x45*0.5,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone0.SetStripNumberingOrder(RightToLeft);
            lower_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for lower layer in zone 0
            const Int_t NPoints_DeadZone_lower_zone0 = 4;
            Double_t XPoints_DeadZone_lower_zone0[NPoints_DeadZone_lower_zone0];
            Double_t YPoints_DeadZone_lower_zone0[NPoints_DeadZone_lower_zone0];

            XPoints_DeadZone_lower_zone0[0] = XPosition+0.0;
            XPoints_DeadZone_lower_zone0[1] = XPosition+0.0;
            XPoints_DeadZone_lower_zone0[2] = XPosition+XRectHotZoneSize_Plane163x45;
            XPoints_DeadZone_lower_zone0[3] = XPosition+XRectHotZoneSize_Plane163x45;

            YPoints_DeadZone_lower_zone0[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_lower_zone0[1] = YPosition-YModuleSize_Plane163x45*0.5 + YRectHotZoneSize_Plane163x45;
            YPoints_DeadZone_lower_zone0[2] = YPosition-YModuleSize_Plane163x45*0.5 + YRectHotZoneSize_Plane163x45;
            YPoints_DeadZone_lower_zone0[3] = YPosition-YModuleSize_Plane163x45*0.5;

            lower_layer_zone0.AddDeadZone(NPoints_DeadZone_lower_zone0, XPoints_DeadZone_lower_zone0, YPoints_DeadZone_lower_zone0);
            //----------------------------------------------------------------------

            //upper strip layer (zone 0) -------------------------------------------
            BmnGemStripLayer upper_layer_zone0(0, UpperStripLayer,
                                               XModuleSize_Plane163x45, YModuleSize_Plane163x45,
                                               XPosition+0.0, YPosition-YModuleSize_Plane163x45*0.5,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone0.SetStripNumberingOrder(RightToLeft);
            upper_layer_zone0.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for upper layer in zone 0
            const Int_t NPoints_DeadZone_upper_zone0 = 4;
            Double_t XPoints_DeadZone_upper_zone0[NPoints_DeadZone_upper_zone0];
            Double_t YPoints_DeadZone_upper_zone0[NPoints_DeadZone_upper_zone0];

            XPoints_DeadZone_upper_zone0[0] = XPosition+0.0;
            XPoints_DeadZone_upper_zone0[1] = XPosition+0.0;
            XPoints_DeadZone_upper_zone0[2] = XPosition+XSlopeHotZoneSize_Plane163x45[0];
            XPoints_DeadZone_upper_zone0[3] = XPosition+XSlopeHotZoneSize_Plane163x45[1];

            YPoints_DeadZone_upper_zone0[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone0[1] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone0[2] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone0[3] = YPosition-YModuleSize_Plane163x45*0.5;

            upper_layer_zone0.AddDeadZone(NPoints_DeadZone_upper_zone0, XPoints_DeadZone_upper_zone0, YPoints_DeadZone_upper_zone0);
            //----------------------------------------------------------------------

            //lower strip layer (zone 1) -------------------------------------------
            BmnGemStripLayer lower_layer_zone1(1, LowerStripLayer,
                                               XRectHotZoneSize_Plane163x45, YRectHotZoneSize_Plane163x45,
                                               XPosition+0.0, YPosition-YModuleSize_Plane163x45*0.5,
                                               LowerLayerPitch, LowerLayerStripAngle);

            lower_layer_zone1.SetStripNumberingOrder(RightToLeft);
            lower_layer_zone1.SetStripNumberingBorders(LeftBottom, RightTop);

            //dead zone for lower layer in zone 1 (hole)
            const Int_t NPoints_DeadZone_lower_zone1 = 4;
            Double_t XPoints_DeadZone_lower_zone1[NPoints_DeadZone_lower_zone1];
            Double_t YPoints_DeadZone_lower_zone1[NPoints_DeadZone_lower_zone1];

            XPoints_DeadZone_lower_zone1[0] = XPosition+0.0;
            XPoints_DeadZone_lower_zone1[1] = XPosition+0.0;
            XPoints_DeadZone_lower_zone1[2] = XPosition+BeamHoleRadius;
            XPoints_DeadZone_lower_zone1[3] = XPosition+BeamHoleRadius;

            YPoints_DeadZone_lower_zone1[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_lower_zone1[1] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_lower_zone1[2] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_lower_zone1[3] = YPosition-YModuleSize_Plane163x45*0.5;

            lower_layer_zone1.AddDeadZone(NPoints_DeadZone_lower_zone1, XPoints_DeadZone_lower_zone1, YPoints_DeadZone_lower_zone1);
            //----------------------------------------------------------------------

            //upper strip layer (zone 1) -------------------------------------------
            BmnGemStripLayer upper_layer_zone1(1, UpperStripLayer,
                                               XSlopeHotZoneSize_Plane163x45[1], YSlopeHotZoneSize_Plane163x45,
                                               XPosition+0.0, YPosition-YModuleSize_Plane163x45*0.5,
                                               UpperLayerPitch, UpperLayerStripAngle);

            upper_layer_zone1.SetStripNumberingOrder(RightToLeft);
            upper_layer_zone1.SetStripNumberingBorders(LeftBottom, RightBottom);

            //dead zone for upper layer in zone 1 (hole)
            const Int_t NPoints_DeadZone_upper_zone1_hole = 4;
            Double_t XPoints_DeadZone_upper_zone1_hole[NPoints_DeadZone_upper_zone1_hole];
            Double_t YPoints_DeadZone_upper_zone1_hole[NPoints_DeadZone_upper_zone1_hole];

            XPoints_DeadZone_upper_zone1_hole[0] = XPosition+0.0;
            XPoints_DeadZone_upper_zone1_hole[1] = XPosition+0.0;
            XPoints_DeadZone_upper_zone1_hole[2] = XPosition+BeamHoleRadius;
            XPoints_DeadZone_upper_zone1_hole[3] = XPosition+BeamHoleRadius;

            YPoints_DeadZone_upper_zone1_hole[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone1_hole[1] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_upper_zone1_hole[2] = YPosition-YModuleSize_Plane163x45*0.5 + BeamHoleRadius;
            YPoints_DeadZone_upper_zone1_hole[3] = YPosition-YModuleSize_Plane163x45*0.5;

            upper_layer_zone1.AddDeadZone(NPoints_DeadZone_upper_zone1_hole, XPoints_DeadZone_upper_zone1_hole, YPoints_DeadZone_upper_zone1_hole);

            //dead zone for upper layer in zone 1 (triangle)
            const Int_t NPoints_DeadZone_upper_zone1_triang = 3;
            Double_t XPoints_DeadZone_upper_zone1_triang[NPoints_DeadZone_upper_zone1_triang];
            Double_t YPoints_DeadZone_upper_zone1_triang[NPoints_DeadZone_upper_zone1_triang];

            XPoints_DeadZone_upper_zone1_triang[0] = XPosition+XSlopeHotZoneSize_Plane163x45[1];
            XPoints_DeadZone_upper_zone1_triang[1] = XPosition+XSlopeHotZoneSize_Plane163x45[0];
            XPoints_DeadZone_upper_zone1_triang[2] = XPosition+XSlopeHotZoneSize_Plane163x45[1];

            YPoints_DeadZone_upper_zone1_triang[0] = YPosition-YModuleSize_Plane163x45*0.5;
            YPoints_DeadZone_upper_zone1_triang[1] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;
            YPoints_DeadZone_upper_zone1_triang[2] = YPosition-YModuleSize_Plane163x45*0.5 + YSlopeHotZoneSize_Plane163x45;

            upper_layer_zone1.AddDeadZone(NPoints_DeadZone_upper_zone1_triang, XPoints_DeadZone_upper_zone1_triang, YPoints_DeadZone_upper_zone1_triang);

            //----------------------------------------------------------------------

            Modules[1]->AddStripLayer(lower_layer_zone0);
            Modules[1]->AddStripLayer(upper_layer_zone0);
            Modules[1]->AddStripLayer(lower_layer_zone1);
            Modules[1]->AddStripLayer(upper_layer_zone1);

    }
//--------------------------------------------------------------------------

    return;
}

ClassImp(BmnGemStripStation_RunSummer2016)