
#include "BmnLANDHit.h"

BmnLANDHit::BmnLANDHit() : BmnHit() {
}

BmnLANDHit::BmnLANDHit(Int_t plane, Int_t bar, TVector3 pos, TVector3 dpos, Float_t time, Float_t energy)
: BmnHit(plane,pos,dpos,bar),
  fTime(time),
  fEnergy(energy)
{
}

BmnLANDHit::~BmnLANDHit() {
}

Float_t BmnLANDHit::GetEnergy() const
{
	return fEnergy;
}

Float_t BmnLANDHit::GetTime() const
{
	return fTime;
}

ClassImp(BmnLANDHit)
