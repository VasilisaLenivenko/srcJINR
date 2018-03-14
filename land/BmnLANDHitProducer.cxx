//--------------------------------------------------------------------------------------------------------------------------------------
#include <cassert>
#include <TRandom2.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TVector3.h>

#include "FairLogger.h"

#include "CbmMCTrack.h"

#include "BmnLANDHitProducer.h"
#include "BmnLANDCand.h"

using namespace std;

namespace{
	bool istokenish(int a_c)
	{
		return isalnum(a_c) || '_' == a_c || '.' == a_c;
	}
	
	std::vector<std::string> tokenize(std::string const &a_str)
	{
		std::vector<std::string> result;
		unsigned i = 0;
		for (;;){
			for (;;){
				if (a_str.size() == i) {
					return result;
				}
				if (istokenish(a_str.at(i))){
					break;
				}
					++i;
				}
			unsigned start = i;
			for (; i < a_str.size() && istokenish(a_str.at(i)); ++i);
			result.push_back(a_str.substr(start, i-start));
		}
	}
}

static Float_t workTime = 0.0;

ClassImp(BmnLANDHitProducer)
	//--------------------------------------------------------------------------------------------------------------------------------------
BmnLANDHitProducer::BmnLANDHitProducer(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test)
	:  FairTask(name,verbose), fLandDigits(nullptr), fUseMCData(false), fLandHits(nullptr){
	}
//--------------------------------------------------------------------------------------------------------------------------------------
BmnLANDHitProducer::~BmnLANDHitProducer() 
{
	//delete pGeoUtils;
}
//--------------------------------------------------------------------------------------------------------------------------------------
InitStatus BmnLANDHitProducer::Init() 
{
	FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "Begin [BmnLANDHitProducer::Init].");

	fLandDigits = (TClonesArray*) FairRootManager::Instance()->GetObject("LAND");
	if (!fLandDigits) {
		cout<<"BmnLANDHitProducer::Init(): branch LAND not found! Task will be deactivated"<<endl;
		SetActive(kFALSE);
		return kERROR;
	}

	fLandHits = new TClonesArray("BmnLANDHit");
	FairRootManager::Instance()->Register("BmnLandHit", "LAND", fLandHits, kTRUE);
	fLandCands = new TClonesArray("BmnLANDCand");
	FairRootManager::Instance()->Register("BmnLandCand", "LAND", fLandCands, kTRUE);
	
	if (!fUseMCData) {
		TString vScintMap = "neuland_sync_0.txt";
		SetVelMap(vScintMap);
	}

	FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "Initialization [BmnLANDHitProducer::Init] finished succesfully.");

	return kSUCCESS;
}

void BmnLANDHitProducer::Exec(Option_t* opt)
{
	if (!IsActive())
		return;

//std::cout << "LAND exec\n";
	clock_t tStart = clock();
	if (fVerbose) cout << endl << "======================== LAND exec started ====================" << endl;

	fLandHits->Clear();
	BmnTrigDigit* digit_T0 = (BmnTrigDigit*) fLandDigits->At(0);

	// Where to put LAND from the target.
	float offset_x = -130.9;
	float offset_z = 1425;
	float angle = atan(offset_x / offset_z);

	for (Int_t iDig = 1; iDig < fLandDigits->GetEntriesFast(); ++iDig) {
		BmnLANDDigit* digit = (BmnLANDDigit*) fLandDigits->At(iDig);

		int p = digit->GetPlane();
		int b = digit->GetBar();		

		TVector3 pos;
		pos.SetXYZ(digit->GetX(), digit->GetY(), (p * 10. + 5.));

		// Assume t res is 500ps for now.
		float err_along = m_vscint[p][b].vscint * sqrt(2.0) * 0.5;
		float err_z = 10. / sqrt(12.);

		// Use correct error depending on plane orientation.
		float err_x;
		float err_y;
		switch (p) {
			case 0:
			case 2:
			case 4:
			case 5:
				err_x = err_z;
				err_y = err_along;
				break;
			case 1:
			case 3:
				err_x = err_along;
				err_y = err_z;
				break;
			default:
				std::cerr << "Invalid LAND plane " << p << ".\n";
				return;
		}

		// Transform error to lab frame.
		float lab_err_x = sqrt(pow(err_x,2) * pow(TMath::Cos(5.2*TMath::DegToRad()),2) + pow(err_y,2) * pow(TMath::Sin(5.2*TMath::DegToRad()),2));
		float lab_err_y = sqrt(pow(err_x,2) * pow(TMath::Sin(5.2*TMath::DegToRad()),2) + pow(err_y,2) * pow(TMath::Cos(5.2*TMath::DegToRad()),2));

		TVector3 dpos;
		dpos.SetXYZ(lab_err_x, lab_err_y, err_z);			

		TVector3 poslab;
		poslab.SetXYZ(pos.X() + offset_x, pos.Y(), pos.Z() + offset_z);
		poslab.RotateY(angle * TMath::DegToRad());

		BmnLANDHit *pHit = new ((*fLandHits)[fLandHits->GetEntriesFast()])
			BmnLANDHit(digit->GetPlane(), digit->GetBar(), poslab, dpos, digit->GetTime(), digit->GetEnergy());

		// TODO: apply slewing correction for T0 time
		pHit->SetTimeStamp(digit->GetTime() - digit_T0->GetTime());
	}

	// TODO: Split LAND and VETO hits.

	/*
	 * Clusterize, clusterification, whatever.
	 */
#define CLUSTERIZE_DIST 27.5
#define VETO_DIST 20.
	// Cluster is an array of digits.
	struct Cluster {
		BmnLANDHit const *hit_array[100];
		size_t hit_num;
	};
	// Array of clusters.
	Cluster *cluster_array[100];
	size_t cluster_num = 0;
	// Leak check.
	int alloc_num = 0;
	for (Int_t hit_i = 0; hit_i < fLandHits->GetEntriesFast(); ++hit_i) {
		BmnLANDHit *hit = (BmnLANDHit *)fLandHits->At(hit_i);

		if (5 == hit->GetDetectorID()) {
			// Skip VETO hits.
			continue;
		}

//std::cout << "Hit=" << hit_i << '\n';
		auto x = hit->GetX();
		auto y = hit->GetY();
		auto z = hit->GetZ();

		// Look for all clusters to which this hit can be added.
		Cluster *cand_array[100];
		size_t cand_num = 0;
		for (size_t cluster_i = 0; cluster_i < cluster_num; ++cluster_i) {
			auto const cluster = cluster_array[cluster_i];
			for (size_t hit_j = 0; hit_j < cluster->hit_num; ++hit_j) {
				auto h = cluster->hit_array[hit_j];
				auto dx = x - h->GetX();
				auto dy = y - h->GetY();
				auto dz = z - h->GetZ();
				if (dx * dx + dy * dy + dz * dz < CLUSTERIZE_DIST * CLUSTERIZE_DIST) {
					cand_array[cand_num++] = cluster;
					break;
//std::cout << " Cand=" << cluster_i << '\n';
				}
			}
		}
//std::cout << "Cands=" << cand_num << '\n';

		Cluster *cluster;
		if (0 == cand_num) {
			// Create cluster.
			cluster = cluster_array[cluster_num++] = new Cluster;
			cluster->hit_num = 0;
			++alloc_num;
		} else {
			// Merge all candidate clusters into the first candidate.
			cluster = cand_array[0];
			for (size_t cand_i = 1; cand_i < cand_num; ++cand_i) {
				auto cluster_src = cand_array[cand_i];
				assert(cluster != cluster_src);
				// Move hits from source cluster into destination cluster.
				for (size_t hit_j = 0; hit_j < cluster_src->hit_num; ++hit_j) {
					cluster->hit_array[cluster->hit_num++] = cluster_src->hit_array[hit_j];
				}
				cluster_src->hit_num = 0;
			}
			// Remove cleared out clusters.
			size_t cluster_dst_i = 1;
			for (size_t cluster_src_i = 1; cluster_src_i < cluster_num; ++cluster_src_i) {
				auto cluster_src = cluster_array[cluster_src_i];
				cluster_array[cluster_dst_i] = cluster_src;
				if (0 != cluster_src->hit_num) {
					++cluster_dst_i;
				} else {
					delete cluster_src;
					--alloc_num;
				}
			}
			cluster_num = cluster_dst_i;
		}

		// Put hit in cluster.
		cluster->hit_array[cluster->hit_num++] = hit;
	}

//std::cout << "Clusters=" << cluster_num << " Allocs=" << alloc_num << '\n';
	for (size_t cluster_i = 0; cluster_i < cluster_num; ++cluster_i) {
		auto cluster = cluster_array[cluster_i];
//std::cout << " Digits=" << cluster->hit_num << '\n';
		BmnLANDHit const *hit_first = nullptr;
		float energy = 0.;
		for (size_t hit_j = 0; hit_j < cluster->hit_num; ++hit_j) {
			auto const hit = cluster->hit_array[hit_j];
			if (!hit_first || hit->GetTime() < hit_first->GetTime()) {
				hit_first = hit;
			}
			energy += hit->GetEnergy();
//std::cout
//	<< "  (x=" << hit->GetX()
//	<< ", y=" << hit->GetY()
//	<< ", z=" << hit->GetZ() << ") t=" << hit->GetTime() << " q=" << hit->GetEnergy() << "\n";
		}
//std::cout << " Time=" << hit_first->GetTime() << " Energy=" << energy << '\n';
		delete cluster;
		--alloc_num;
		// Look if there's a VETO in front.
		int veto = 0;
		for (Int_t hit_j = 0; hit_j < fLandHits->GetEntriesFast(); ++hit_j) {
			auto const hit = (BmnLANDHit *)fLandHits->At(hit_j);
			if (5 != hit->GetDetectorID()) {
				// Skip LAND hits.
				continue;
			}
			double dx = hit->GetX() - hit_first->GetX();
			double dy = hit->GetY() - hit_first->GetY();
			if (dx * dx + dy * dy < VETO_DIST * VETO_DIST) {
				veto = 1;
				break;
			}
		}
//std::cout << " VETO=" << veto << '\n';
		BmnLANDCand *cand = new((*fLandCands)[fLandCands->GetEntriesFast()])
		    BmnLANDCand(*hit_first, energy, veto);
	}
	assert(0 == alloc_num);

	clock_t tFinish = clock();
	workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
	if (fVerbose) cout << "======================== LAND exec finished ====================" << endl;
}
//--------------------------------------------------------------------------------------------------------------------------------------

void BmnLANDHitProducer::Finish() {
	/*if (fDoTest) {
		FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "[BmnLANDHitProducer::Finish] Update  %s file. ", fTestFlnm.Data());
		TFile *ptr = gFile;
		TFile file(fTestFlnm.Data(), "RECREATE");
		fList.Write();
		file.Close();
		gFile = ptr;
		if (!fUseMCData) 
			for (Int_t i = 0; i < fNDetectors; i++)
				pDetector[i] -> SaveHistToFile(fTestFlnm.Data());
	}*/

	cout << "Work time of the LAND hit finder: " << workTime << endl;
}


bool BmnLANDHitProducer::SetVelMap(TString a_vscint_filename)
{
	auto path = std::string(getenv("VMCWORKDIR")) + "/input/" + a_vscint_filename.Data();
	std::ifstream in(path.c_str());
	if (!in.is_open()){
		std::cerr << " Could not open vscint map for LAND " << path << "\n";
		return false;
	}
	for (unsigned line_no = 1;; ++line_no){
		std::string line;
		std::getline(in,line);
		if( !in.good() ) break;
		int globbar;
		auto const token = tokenize(line);
		if (token.size() < 4) continue;
		char const *p;
		char *end;
		p = token.at(0).c_str();
		globbar = strtol(p, &end, 10);
		if (end == p || globbar < 1 || globbar > 120) continue;
		--globbar;
		int land_plane = globbar / LAND_BAR_N;
		int land_bar = globbar % LAND_BAR_N;
		m_vscint[land_plane][land_bar].vscint = strtod(token.at(3).c_str(), NULL);
	}
	in.close();	
	return true;
}
