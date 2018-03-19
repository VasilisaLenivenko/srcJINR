//--------------------------------------------------------------------------------------------------------------------------------------
#include<assert.h>
#include<map>

#include <TRandom2.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TVector3.h>
#include "BmnTrigWaveDigit.h"
#include "FairLogger.h"

#include "CbmMCTrack.h"

#include "BmnTrigHitProducer.h"

using namespace std;


ClassImp(BmnTrigHitProducer)
	//--------------------------------------------------------------------------------------------------------------------------------------
BmnTrigHitProducer::BmnTrigHitProducer(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test)
	:  FairTask(name,verbose), fOnlyPrimary(false), fUseMCData(false), aTrigHits(nullptr), BC3_digits(nullptr), BC4_digits(nullptr){

	}
//--------------------------------------------------------------------------------------------------------------------------------------
BmnTrigHitProducer::~BmnTrigHitProducer() 
{
	//delete pGeoUtils;
}
//--------------------------------------------------------------------------------------------------------------------------------------
InitStatus 		BmnTrigHitProducer::Init() 
{
	FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "Begin [BmnTrigHitProducer::Init].");

	if(fOnlyPrimary) cout<<" Only primary particles are processed!!! \n"; // FIXME NOT used now ADDD

	BC3_digits = (TClonesArray*) FairRootManager::Instance()->GetObject("TQDC_BC3");
	if (!BC3_digits)
	{
		cout<<"BmnTrigHitProducer::Init(): branch BC3 not found! Task will be deactivated"<<endl;
		SetActive(kFALSE);
		return kERROR;
	}
	BC4_digits = (TClonesArray*) FairRootManager::Instance()->GetObject("TQDC_BC4");
	if (!BC3_digits)
	{
		cout<<"BmnTrigHitProducer::Init(): branch BC4 not found! Task will be deactivated"<<endl;
		SetActive(kFALSE);
		return kERROR;
	}
	T0_digits = (TClonesArray*) FairRootManager::Instance()->GetObject("T0");
	if (!T0_digits)
	{
		cout<<"BmnTrigHitProducer::Init(): branch T0 not found! Task will be deactivated"<<endl;
		SetActive(kFALSE);
		return kERROR;
	}
	

	aTrigHits = new TClonesArray("BmnTrigHit");
	FairRootManager::Instance()->Register("BmnTrigHit", "Triggers", aTrigHits, kTRUE);
	

	FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "Initialization [BmnTrigHitProducer::Init] finished succesfully.");

	return kSUCCESS;
}

void 		BmnTrigHitProducer::Exec(Option_t* opt) {
	

	
	if (!IsActive())
		return;
	clock_t tStart = clock();

	if (fVerbose) cout << endl << "======================== Triggers exec started ====================" << endl;
	

	aTrigHits->Clear();
	Int_t nT0Digits = T0_digits->GetEntriesFast();
	Int_t nBC3 = BC3_digits->GetEntriesFast();
	Int_t nBC4 = BC4_digits->GetEntriesFast();
	
	double fEn3 = 0.;
	double fEn4 = 0.;
	if( nT0Digits!=0 ){
		if( (nBC3!=0) && (nBC4!=0)  ){
				BmnTrigWaveDigit* signalBC3 = (BmnTrigWaveDigit*) BC3_digits->At(0);
				fEn3 += signalBC3->GetPeak();	
		
				BmnTrigWaveDigit* signalBC4 = (BmnTrigWaveDigit*) BC4_digits->At(0);
				fEn4 += signalBC4->GetPeak();	


			BmnTrigHit *pHit = new ((*aTrigHits)[aTrigHits->GetEntriesFast()]) BmnTrigHit(fEn3,fEn4 );
			pHit->SetBC3(fEn3);
			pHit->SetBC4(fEn4);
		}
	}	




    	//BmnTriggersHit *pHit = new ((*aLandHits)[aLandHits->GetEntriesFast()]) BmnTriggersHit(digLand->GetPlane(), digLand->GetBar(), poslab, dpos,digLand->GetTime(), digLand->GetEnergy());

	clock_t tFinish = clock();

	if (fVerbose) cout << "======================== Triggers exec finished ====================" << endl;
}
//--------------------------------------------------------------------------------------------------------------------------------------

void BmnTrigHitProducer::Finish() {
}


