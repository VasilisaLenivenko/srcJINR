#ifndef BMNLANDCAND_H
#define BMNLANDCAND_H

#include <math.h>
#include <iostream>
#include "BmnHit.h"
#include <TObject.h>

class BmnLANDHit;

class BmnLANDCand : public BmnHit {
public:

    /** Default constructor **/
    BmnLANDCand();

    /** Constructor to use **/
    BmnLANDCand(BmnLANDHit const &, Float_t, Int_t);

    /** Destructor **/
    virtual ~BmnLANDCand();
    
    Float_t GetTime() const{
	return fTime;
    }
    Float_t GetEnergy() const{
	return fEnergy;
    }

private:
    Float_t fTime;
    Float_t fEnergy;
    Int_t fVeto;
    ClassDef(BmnLANDCand, 1);
};

#endif /* BMNLANDCAND_H */
