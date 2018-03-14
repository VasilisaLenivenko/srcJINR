#include "BmnLANDCand.h"
#include "BmnLANDHit.h"

BmnLANDCand::BmnLANDCand(): BmnHit() {
}

BmnLANDCand::BmnLANDCand(BmnLANDHit const &a_hit, Float_t a_energy, Int_t a_veto)
: BmnHit(a_hit.GetDetId(),
         TVector3(a_hit.GetX(), a_hit.GetY(), a_hit.GetZ()),
         TVector3(a_hit.GetDx(), a_hit.GetDy(), a_hit.GetDz()),
         a_hit.GetRefIndex()),
  fTime(a_hit.GetTime()),
  fEnergy(a_energy),
  fVeto(a_veto)
{
}

BmnLANDCand::~BmnLANDCand()
{
}

ClassImp(BmnLANDCand)
