//--------------------------------------------------------------------------------------------------------------------------------------
#ifndef __BmnTrigHitProducer_H
#define __BmnTrigHitProducer_H 1

#include <TClonesArray.h>
#include <FairTask.h>
#include <TVector3.h>
#include "BmnTrigHit.h"
#include "BmnTrigDigit.h"
#include "BmnTrigWaveDigit.h"

class TRandom2;
class TEfficiency;
class TH1D;
class TH2D;
//class BmnTof1GeoUtils;
//--------------------------------------------------------------------------------------------------------------------------------------
class BmnTrigHitProducer : public FairTask 
{
protected:
	TClonesArray		*BC3_digits;   // Exp input
	TClonesArray		*BC4_digits;   // Exp input
	TClonesArray		*T0_digits;
	TClonesArray		*aTrigHits;     //! output
	
	Bool_t			fOnlyPrimary;
	Bool_t			fUseMCData;
public:
	BmnTrigHitProducer(const char *name = "Trigger HitProducer", Bool_t useMCdata = true, Int_t verbose = 1, Bool_t DoTest = false);
	virtual ~BmnTrigHitProducer();

	virtual InitStatus 	Init();
	virtual void 		Exec(Option_t* opt);
	virtual void 		Finish();


private:

ClassDef(BmnTrigHitProducer, 2);
};

#endif
//--------------------------------------------------------------------------------------------------------------------------------------

