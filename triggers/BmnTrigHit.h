/* 
 * File:   BmnTrigHit.h
 * Author: segarra
 *
 * Created on Feb 14, 2018, 10:51 AM
 */

/*
TODO: add slewing, make sure that spacing of layers is only thickness of bars 10cm, add T0

BmnLandHit->GetDetectorID()  			returns plane of LAND (5 = veto, 0-4 = LAND)
BmnLandHit->GetRefIndex()			returns bar index in plane (0-20)
BmnLandHit->GetX() / GetY() / GetZ()		returns global position of LAND in lab frame
BmnLandHit->GetTimeStamp()			returns time rel to T0 with slewing corr
BmnLandHit->GetEnergy()				returns calibrated energy deposit for hit
BmnLandHit->GetDx() / GetDy() / GetDz()		returns error in position calculation

*/


#ifndef BMNTRIGGERHIT_H
#define BMNTRIGGERHIT_H

#include <math.h>
#include <iostream>
#include "BmnHit.h"
#include <TObject.h>
using namespace std;

// class TClonesArray;

class BmnTrigHit : public BmnHit {
public:

    /** Default constructor **/
    BmnTrigHit();

    /** Constructor to use **/
    BmnTrigHit(Float_t BC3, Float_t BC4 );

    /** Destructor **/
    virtual ~BmnTrigHit();
    
    Float_t BC3_Peak() const{
	return fBC3;
    }

    Float_t BC4_Peak() const{
	return fBC4;
    }
    void SetBC3(Float_t BC3) {
	fBC3 = BC3;
    }

    void SetBC4(Float_t BC4) {
	fBC4 = BC4;
    }

private:

	Float_t fBC3;
	Float_t fBC4;

    ClassDef(BmnTrigHit, 1);
};

#endif /* BMNLANDHIT_H */

