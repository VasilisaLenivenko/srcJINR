// Two proton events (March 11, 2018)
// If there are any questions, sent me an email (reynier@mit.edu)
#include "TApplication.h"
#include "TRint.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "BmnGemStripHit.h"
#include "BmnTofHit.h"
#include "BmnMwpcHit.h"
#include "BmnTrigHit.h"
#include "SRCEvent.h"
#include "TStopwatch.h"
using namespace std;

const double mp   = 0.938272;
const double Ac12 = 12.;
const double mC12 = 11.1770; // 12 * 0.931494028 GeV
double Pbeam      = 4.; //GeV/c
double Ebeam_tmp=TMath::Sqrt( pow(Ac12*Pbeam,2) + pow(mC12,2) );

TLorentzVector vc12_ForReconstructions( TVector3(0,0,Ac12*Pbeam), Ebeam_tmp  ); // un-smeared Nucleus beam 4-momentum (DUBNA lab frame)
TLorentzVector v_proton_beam_Standing ( TVector3(0,0,0         ), mp         ); // Standing proton (DUBNA lab frame)
TLorentzVector v_proton_beam_UnSMEARED( TVector3(0,0,0         ), mp         ); // Standing proton (DUBNA lab frame)

SRCEvent * currentEvent;
// =============================================================================================================================
int main(int argc, char ** argv){
	if (argc !=3){
		cerr << "Wrong number of arguments. Instead use\n"<< "\ttrack_arms /path/to/run/reco/file /path/to/output/file\n";
		return -1;
	}

	v_proton_beam_UnSMEARED.Boost(-vc12_ForReconstructions.BoostVector()); // Boost from DUBNA to GSI frame

	// Graphical Cuts from George
	TCutG *cut;
	TFile *f = new TFile("/nica/user/s/segarrae/software/build-bmnroot/bin/theta_P2_boost_vs_theta_P1_boost_4GeVpercperu.root");
	cut = (TCutG*)f->Get("mycut");

	TCutG *cut1;
	TFile *f1 = new TFile("/nica/user/s/segarrae/software/build-bmnroot/bin/P1_boost_vs_theta_P1_boost_4GeVpercperu.root");
	cut1= (TCutG*)f1->Get("mycut");

	// Defining histograms and other objects

	// Hit coordinates in the TOF400 detectors
	TH2F * XY_hits_tof400 = new TH2F("XY_hits_tof400","XY hits tof400;X [cm];Y [cm]",200,-500,500,200,-80, 80);
	TH2F * XZ_hits_tof400 = new TH2F("XZ_hits_tof400","XZ hits tof400;X [cm];Z [cm]",200,-500,500,200,200,600);
	TH2F * YZ_hits_tof400 = new TH2F("YZ_hits_tof400","YZ hits tof400;Y [cm];Z [cm]",200,- 80, 80,200,200,600);

	// Angular distributions
	TH1F * h1_thetaL  = new TH1F("h1_thetaL" ,"#theta_{L};#theta_{L};Counts",200, 20, 40);
	TH1F * h1_thetaR  = new TH1F("h1_thetaR" ,"#theta_{R};#theta_{R};Counts",200, 20, 40);
	TH1F * h1_phiL    = new TH1F("h1_phiL"   ,"#phi_{L};#phi_{L};Counts"    ,200,-30, 30);
	TH1F * h1_phiR    = new TH1F("h1_phiR"   ,"#phi_{R};#phi_{R};Counts"    ,200,150,210);

	// Correlation plots
	TH2F * h2_thetaLR = new TH2F("h2_thetaLR","#theta_{L} vs. #theta_{R};#theta_{L};#theta_{R}"  , 200 , 20, 40 , 200, 20, 40);
	TH2F * h2_pLR     = new TH2F("h2_pLR"    ,"P_{L} vs. P_{R};P_{L} [GeV];P_{R} [GeV]"          , 200 ,  0,  5 , 200,  0 , 5);

	// Left and right arms momenta
	TH1F * h1_PL      = new TH1F("h1_PL"     ,"momentum, left arm;P_{L} [GeV];Counts"                  , 200, 0,5);
	TH1F * h1_PL_x    = new TH1F("h1_PL_x"   ,"momentum, left arm, x component;P_{L}^{X} [GeV];Counts" , 200,-3,3);
	TH1F * h1_PL_y    = new TH1F("h1_PL_y"   ,"momentum, left arm, y component;P_{L}^{Y} [GeV];Counts" , 200,-2,2);
	TH1F * h1_PL_z    = new TH1F("h1_PL_z"   ,"momentum, left arm, z component;P_{L}^{Z} [GeV];Counts" , 200, 0,5);
	TH1F * h1_PR      = new TH1F("h1_PR"     ,"momentum, right arm;P_{R} [GeV];Counts"                 , 200, 0,5);
	TH1F * h1_PR_x    = new TH1F("h1_PR_x"   ,"momentum, right arm, x component;P_{R}^{X} [GeV];Counts", 200,-3,3);
	TH1F * h1_PR_y    = new TH1F("h1_PR_y"   ,"momentum, right arm, y component;P_{R}^{Y} [GeV];Counts", 200,-2,2);
	TH1F * h1_PR_z    = new TH1F("h1_PR_z"   ,"momentum, right arm, z component;P_{R}^{Z} [GeV];Counts", 200, 0,5);

	int nhitL, nhitR;
	double tofhitL_x, tofhitL_y, tofhitL_z, tofhitL_t;
	double tofhitR_x, tofhitR_y, tofhitR_z, tofhitR_t;
	double thetaL_labf, thetaR_labf, phiL_labf, phiR_labf, dL, dR, ToFL, ToFR;
	double thetaL_Crf , thetaR_Crf , phiL_Crf , phiR_Crf;
	double Emiss, Pmiss, Pmiss_x, Pmiss_y, Pmiss_z;
	double PL_labf, PL_labf_x, PL_labf_y, PL_labf_z, PR_labf, PR_labf_x, PR_labf_y, PR_labf_z, EL, ER;
	double PL_Crf , PL_Crf_x , PL_Crf_y , PL_Crf_z , PR_Crf , PR_Crf_x , PR_Crf_y , PR_Crf_z;
	double mand_s_Crf ,mand_u_Crf ,mand_t_Crf ;
	double mand_s_labf,mand_u_labf,mand_t_labf;
	double BC3_en, BC4_en;
	TVector3 L_tof400_xyz, R_tof400_xyz;
	TLorentzVector P4L        , P4R        ;
	TLorentzVector P4L_boosted, P4R_boosted;

	// Set up the input file
	TFile * infile = NULL;
	infile = new TFile(argv[1]);
	if (!infile){cerr << "Could not open file " << argv[1] <<"\n"<< "\tBailing out\n";	return -2;}
	else{cerr << "Successfully opened file " << argv[1] << " and saved it to address " << infile << "\n";}

	// -----------------------------------------------------------------------------
	// Set up the output file

	TFile * outfile = new TFile(argv[2],"RECREATE");
	TTree * outtree = new TTree("t","Output tree for the track_arms program");

	// TOF400 hits
	outtree->Branch("tofhitL_x"  ,&tofhitL_x  ,"tofhitL_x/D"  );
	outtree->Branch("tofhitL_y"  ,&tofhitL_y  ,"tofhitL_y/D"  );
	outtree->Branch("tofhitL_z"  ,&tofhitL_z  ,"tofhitL_z/D"  );
	outtree->Branch("tofhitL_t"  ,&tofhitL_t  ,"tofhitL_t/D"  );
	outtree->Branch("tofhitR_x"  ,&tofhitR_x  ,"tofhitR_x/D"  );
	outtree->Branch("tofhitR_y"  ,&tofhitR_y  ,"tofhitR_y/D"  );
	outtree->Branch("tofhitR_z"  ,&tofhitR_z  ,"tofhitR_z/D"  );
	outtree->Branch("tofhitR_t"  ,&tofhitR_t  ,"tofhitR_t/D"  );
	outtree->Branch("dL"         ,&dL         ,"dL/D"         );
	outtree->Branch("dR"         ,&dR         ,"dR/D"         );
	outtree->Branch("nhitL"	     ,&nhitL	  ,"nhitL/I"	  );
	outtree->Branch("nhitR"	     ,&nhitR	  ,"nhitR/I"	  );


	// Momenta in the lab frame
	outtree->Branch("PL_labf"    ,&PL_labf    ,"PL_labf/D"    );
	outtree->Branch("PL_labf_x"  ,&PL_labf_x  ,"PL_labf_x/D"  );
	outtree->Branch("PL_labf_y"  ,&PL_labf_y  ,"PL_labf_y/D"  );
	outtree->Branch("PL_labf_z"  ,&PL_labf_z  ,"PL_labf_z/D"  );
	outtree->Branch("PR_labf"    ,&PR_labf    ,"PR_labf/D"    );
	outtree->Branch("PR_labf_x"  ,&PR_labf_x  ,"PR_labf_x/D"  );
	outtree->Branch("PR_labf_y"  ,&PR_labf_y  ,"PR_labf_y/D"  );
	outtree->Branch("PR_labf_z"  ,&PR_labf_z  ,"PR_labf_z/D"  );

	// Momenta in the Carbon rest frame
	outtree->Branch("PL_Crf"     ,&PL_Crf     ,"PL_Crf/D"     );
	outtree->Branch("PL_Crf_x"   ,&PL_Crf_x   ,"PL_Crf_x/D"   );
	outtree->Branch("PL_Crf_y"   ,&PL_Crf_y   ,"PL_Crf_y/D"   );
	outtree->Branch("PL_Crf_z"   ,&PL_Crf_z   ,"PL_Crf_z/D"   );
	outtree->Branch("PR_Crf"     ,&PR_Crf     ,"PR_Crf/D"     );
	outtree->Branch("PR_Crf_x"   ,&PR_Crf_x   ,"PR_Crf_x/D"   );
	outtree->Branch("PR_Crf_y"   ,&PR_Crf_y   ,"PR_Crf_y/D"   );
	outtree->Branch("PR_Crf_z"   ,&PR_Crf_z   ,"PR_Crf_z/D"   );

	// Angles in the lab frame
	outtree->Branch("thetaL_labf",&thetaL_labf,"thetaL_labf/D");
	outtree->Branch("thetaR_labf",&thetaR_labf,"thetaR_labf/D");
	outtree->Branch("phiL_labf"  ,&phiL_labf  ,"phiL_labf/D"  );
	outtree->Branch("phiR_labf"  ,&phiR_labf  ,"phiR_labf/D"  );

	// Angles in the Carbon rest frame
	outtree->Branch("thetaL_Crf" ,&thetaL_Crf ,"thetaL_Crf/D" );
	outtree->Branch("thetaR_Crf" ,&thetaR_Crf ,"thetaR_Crf/D" );
	outtree->Branch("phiL_Crf"   ,&phiL_Crf   ,"phiL_Crf/D"   );
	outtree->Branch("phiR_Crf"   ,&phiR_Crf   ,"phiR_Crf/D"   );

	// Mandelstam variables in the lab frame
	outtree->Branch("mand_s_labf",&mand_s_labf,"mand_s_labf/D");
	outtree->Branch("mand_u_labf",&mand_u_labf,"mand_u_labf/D");
	outtree->Branch("mand_t_labf",&mand_t_labf,"mand_t_labf/D");

	// Mandelstam variables in the Carbon rest frame
	outtree->Branch("mand_s_Crf" ,&mand_s_Crf ,"mand_s_Crf/D" );
	outtree->Branch("mand_u_Crf" ,&mand_u_Crf ,"mand_u_Crf/D" );  
	outtree->Branch("mand_t_Crf" ,&mand_t_Crf ,"mand_t_Crf/D" );

	// Missing momentum
	outtree->Branch("Pmiss"      ,&Pmiss      ,"Pmiss/D"      );
	outtree->Branch("Pmiss_x"    ,&Pmiss_x    ,"Pmiss_x/D"    );
	outtree->Branch("Pmiss_y"    ,&Pmiss_y    ,"Pmiss_y/D"    );
	outtree->Branch("Pmiss_z"    ,&Pmiss_z    ,"Pmiss_z/D"    );

	outtree->Branch("BC3_en"    ,&BC3_en    ,"BC3_en/D"    );
	outtree->Branch("BC4_en"    ,&BC4_en    ,"BC4_en/D"    );
	// -----------------------------------------------------------------------------
	// Set up the tree
	TClonesArray * tofData  = new TClonesArray("BmnTofHit"     );
	TClonesArray * trigData  = new TClonesArray("BmnTrigHit"   );
	TClonesArray * mwpcData = new TClonesArray("BmnMwpcHit"    );
	TClonesArray * gemData  = new TClonesArray("BmnGemStripHit");

	TTree * intree = NULL;
	intree = (TTree*) infile->Get("cbmsim");

	if (!intree){cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";return -3;}
	else{cerr << "Successfully loaded tree at address " << intree << "\n";}

	const int nEvents = intree->GetEntries();
	intree->SetBranchAddress("BmnTof1Hit"    ,&tofData );
	intree->SetBranchAddress("BmnTrigHit"       ,&trigData );
	//intree->SetBranchAddress("BmnMwpcHit"    ,&mwpcData);
	//intree->SetBranchAddress("BmnGemStripHit",&gemData );

	// Shift from downstream coordinate system
	// to target-centered coordinate system:
	double shift_coord = 650.; // cm

	// ========================================================================================
	// IMPORTANT:
	// Since we only have tof400 information at the moment,
	// we are assuming that the interaction vertex is centerred
	// in the very center of the target.
	// ========================================================================================

	cout << "===================================" << endl;
	cout << "Beginning loop over events:" << endl;
	// Loop over events
	for (int event=0 ; event<nEvents ; event++){
		if (event % 10000 == 0) cerr << "Working on event " << event << " out of " << nEvents << endl;
		intree->GetEvent(event);	

		L_tof400_xyz.SetXYZ(0.,0.,0.);
		R_tof400_xyz.SetXYZ(0.,0.,0.);

		nhitL = 0;
		nhitR = 0;

		for(int j=0; j<trigData->GetEntriesFast() ;j++){
			BmnTrigHit * thisHit = (BmnTrigHit*)trigData->At(j);
			BC3_en = thisHit->BC3_Peak();
			BC4_en = thisHit->BC4_Peak();
		}


		for (int t=0 ; t<tofData->GetEntriesFast() ; t++){
			BmnTofHit * thisHit = (BmnTofHit*)tofData->At(t);
			if(thisHit->GetTimeStamp() > 25) continue;
			
			// Hits on the left arm
			if (thisHit->GetX() > 0.){
				L_tof400_xyz.SetXYZ(thisHit->GetX() ,thisHit->GetY() ,(thisHit->GetZ() + shift_coord));
				XY_hits_tof400 -> Fill(L_tof400_xyz.X(),L_tof400_xyz.Y());
				XZ_hits_tof400 -> Fill(L_tof400_xyz.X(),L_tof400_xyz.Z());
				YZ_hits_tof400 -> Fill(L_tof400_xyz.Y(),L_tof400_xyz.Z());
				ToFL = thisHit -> GetTimeStamp();
				nhitL++;
				
				tofhitL_x = L_tof400_xyz.X();
				tofhitL_y = L_tof400_xyz.Y();
				tofhitL_z = L_tof400_xyz.Z();
				tofhitL_t = ToFL;


			}

			// Hits on the right arm
			else{	
				R_tof400_xyz.SetXYZ(thisHit->GetX() ,thisHit->GetY() ,(thisHit->GetZ() + shift_coord));
				XY_hits_tof400 -> Fill(R_tof400_xyz.X(),R_tof400_xyz.Y());
				XZ_hits_tof400 -> Fill(R_tof400_xyz.X(),R_tof400_xyz.Z());
				YZ_hits_tof400 -> Fill(R_tof400_xyz.Y(),R_tof400_xyz.Z());
				ToFR = thisHit -> GetTimeStamp();
				nhitR++;
			
				tofhitR_x = R_tof400_xyz.X();
				tofhitR_y = R_tof400_xyz.Y();
				tofhitR_z = R_tof400_xyz.Z();
				tofhitR_t = ToFR;
			
			}
			outtree->Fill();
		
		}
		/*
		// Ignoring events that have more than two proton candidates
		if ((nhitL!=1) || (nhitR!=1)) continue;

		// -------------------------------
		// ToF400 hits
		tofhitL_x = L_tof400_xyz.X();
		tofhitL_y = L_tof400_xyz.Y();
		tofhitL_z = L_tof400_xyz.Z();
		tofhitL_t = ToFL;
		tofhitR_x = R_tof400_xyz.X();
		tofhitR_y = R_tof400_xyz.Y();
		tofhitR_z = R_tof400_xyz.Z();
		tofhitR_t = ToFR;

		// ----------------------------------------------------------------------------------
		// Calculating several quantities

		thetaL_labf = L_tof400_xyz.Theta();
		thetaR_labf = R_tof400_xyz.Theta();
		phiL_labf   = L_tof400_xyz.Phi  ();
		phiR_labf   = R_tof400_xyz.Phi  ();
		if(phiR_labf<0) phiR_labf += 2*3.14159;

		// The momentum calculation below is given by:
		//        ________________      __________________________________________
		//       |     m^2             |                  m^2
		// P = /\| --------------  = /\| ----------------------------------------- , where 3.3 = 1/(speed of light [meters/nanosec])
		//       |  1/beta^2 - 1       |  (tof[nanosec]/(3.3*path[meters]))^2 - 1

		// -------------------------------
		// Left arm variables
		//
		dL = sqrt(pow(L_tof400_xyz.X(),2) + pow(L_tof400_xyz.Y(),2) + pow(L_tof400_xyz.Z(),2))/100.; // [meter]
		PL_labf   = sqrt( pow(mp,2) / (pow(ToFL/(3.3*dL),2) - 1) );
		PL_labf_x = PL_labf*cos(phiL_labf)*sin(thetaL_labf);
		PL_labf_y = PL_labf*sin(phiL_labf)*sin(thetaL_labf);
		PL_labf_z = PL_labf*cos(thetaL_labf);
		EL   = sqrt(pow(PL_labf,2)+pow(mp,2));
		P4L        .SetXYZM(PL_labf_x,PL_labf_y,PL_labf_z,mp);
		P4L_boosted.SetXYZM(PL_labf_x,PL_labf_y,PL_labf_z,mp);
		P4L_boosted.Boost(-vc12_ForReconstructions.BoostVector()); // Boosting to Carbon rest frame

		PL_Crf   = P4L_boosted.Rho();
		PL_Crf_x = P4L_boosted.X();
		PL_Crf_y = P4L_boosted.Y();
		PL_Crf_z = P4L_boosted.Z();	

		// -------------------------------
		// Right arm variables
		dR = sqrt(pow(R_tof400_xyz.X(),2) + pow(R_tof400_xyz.Y(),2) + pow(R_tof400_xyz.Z(),2))/100.; // [meter]
		PR_labf   = sqrt( pow(mp,2) / (pow(ToFR/(3.3*dR),2) - 1) );
		PR_labf_x = PR_labf*cos(phiR_labf)*sin(thetaR_labf);
		PR_labf_y = PR_labf*sin(phiR_labf)*sin(thetaR_labf);
		PR_labf_z = PR_labf*cos(thetaR_labf);
		ER   = sqrt(pow(PR_labf,2)+pow(mp,2));
		P4R        .SetXYZM(PR_labf_x,PR_labf_y,PR_labf_z,mp);
		P4R_boosted.SetXYZM(PR_labf_x,PR_labf_y,PR_labf_z,mp);
		P4R_boosted.Boost(-vc12_ForReconstructions.BoostVector()); // Boosting to Carbon rest frame

		PR_Crf   = P4R_boosted.Rho();
		PR_Crf_x = P4R_boosted.X();
		PR_Crf_y = P4R_boosted.Y();
		PR_Crf_z = P4R_boosted.Z();

		// Angles in the boosted frame
		thetaL_Crf = P4L_boosted.Theta();
		thetaR_Crf = P4R_boosted.Theta();
		phiL_Crf   = P4L_boosted.Phi  ();
		phiR_Crf   = P4R_boosted.Phi  ();

		// -------------------------------
		// Missing momentum
		Pmiss_x = P4L_boosted.X() + P4R_boosted.X();
		Pmiss_y = P4L_boosted.Y() + P4R_boosted.Y();
		Pmiss_z = P4L_boosted.Z() + P4R_boosted.Z() + v_proton_beam_UnSMEARED.P();
		Pmiss   = sqrt(pow(Pmiss_x,2) + pow(Pmiss_y,2) + pow(Pmiss_z,2));
		Emiss   = sqrt(pow(Pmiss,2) + pow(mp,2));
		// ----------------------------------------------------------------------------------
		// Filling histograms

		h1_thetaL -> Fill(thetaL_labf*180./3.14159);
		h1_thetaR -> Fill(thetaR_labf*180./3.14159);
		h1_phiL   -> Fill(phiL_labf  *180./3.14159);
		h1_phiR   -> Fill(phiR_labf  *180./3.14159);
		h2_thetaLR-> Fill(thetaL_labf*180./3.14159,thetaR_labf*180./3.14159);

		h1_PL     -> Fill(PL_labf  );
		h1_PL_x   -> Fill(PL_labf_x);
		h1_PL_y   -> Fill(PL_labf_y);
		h1_PL_z   -> Fill(PL_labf_z);

		h1_PR     -> Fill(PR_labf  );
		h1_PR_x   -> Fill(PR_labf_x);
		h1_PR_y   -> Fill(PR_labf_y);
		h1_PR_z   -> Fill(PR_labf_z);

		h2_pLR    -> Fill(PL_labf,PR_labf);

		// ----------------------------------------------------------------------------------
		// Calculating Mandelstam variables

		// Frame in which proton is moving, Carbon is at rest (Carbon rest frame)
		mand_s_Crf     = (P4L_boosted+P4R_boosted)*(P4L_boosted+P4R_boosted);
		mand_u_Crf     = (v_proton_beam_UnSMEARED-P4R_boosted)*(v_proton_beam_UnSMEARED-P4R_boosted);
		mand_t_Crf     = (v_proton_beam_UnSMEARED-P4L_boosted)*(v_proton_beam_UnSMEARED-P4L_boosted);

		// Frame in which proton is standing, Carbon is moving (lab frame)
		mand_s_labf = (P4L+P4R)*(P4L+P4R);
		mand_u_labf = (v_proton_beam_Standing-P4R)*(v_proton_beam_Standing-P4R);
		mand_t_labf = (v_proton_beam_Standing-P4L)*(v_proton_beam_Standing-P4L);

		// Mandelstam variables are Lorentz-invariant, so the values above should be
		// the same before and after the boost. Let's make sure that's the case:
		if(             (mand_s_Crf-mand_s_labf>0.001)||
				(mand_u_Crf-mand_u_labf>0.001)||
				(mand_t_Crf-mand_t_labf>0.001)
		  ) cout << "There seems to be a problem with the Mandelstam variables. Check it!!!" << endl;

		// Applying graphical cuts
		if (!(cut->IsInside(thetaR_labf*180./TMath::Pi(),thetaL_labf*180./TMath::Pi())))   continue;
		//if (!(cut1->IsInside(PL_labf,thetaL_labf*180./TMath::Pi())))  continue;
		//if (!(cut1->IsInside(PR_labf,thetaR_labf*180./TMath::Pi())))      continue;
		if ((PL_Crf+PR_Crf)<4.4 || (PL_Crf+PR_Crf)<4.8)  continue; 
		
		// ----------------------------------------------------------------------------------
		// Saving informormation onto output tree
		outtree->Fill();
		*/
	}

	/*
	TCanvas * c1 = new TCanvas("c1");	XY_hits_tof400 -> Draw("COLZ");
	TCanvas * c2 = new TCanvas("c2");       XZ_hits_tof400 -> Draw("COLZ");
	TCanvas * c3 = new TCanvas("c3");       YZ_hits_tof400 -> Draw("COLZ");

	TCanvas * c4 = new TCanvas("c4","Angles");
	c4 -> Divide(2,2);
	c4 -> cd(1);	h1_thetaL -> Draw();
	c4 -> cd(2);	h1_thetaR -> Draw();
	c4 -> cd(3);	h1_phiL   -> Draw();
	c4 -> cd(4);	h1_phiR   -> Draw();

	TCanvas * c5 = new TCanvas("c5","Momenta");
	c5 -> Divide(4,2);
	c5 -> cd(1);	h1_PL   -> Draw();
	c5 -> cd(2);    h1_PL_x -> Draw();
	c5 -> cd(3);    h1_PL_y -> Draw();
	c5 -> cd(4);    h1_PL_z -> Draw();
	c5 -> cd(5);    h1_PR   -> Draw();
	c5 -> cd(6);    h1_PR_x -> Draw();
	c5 -> cd(7);    h1_PR_y -> Draw();
	c5 -> cd(8);    h1_PR_z -> Draw();

	TCanvas * c6 = new TCanvas("c6","Correlations");
	c6 -> Divide(2,1);
	c6 -> cd(1);	h2_thetaLR -> Draw("COLZ");
	c6 -> cd(2);	h2_pLR     -> Draw("COLZ");

	// ----------------------------------------------------------------------------------
	// Writing information onto output file
	*/
	infile ->Close();	// Closing input file
	outtree->Write();	// Writing tree onto output file
	//c1     ->Write();
	//c2     ->Write();
	//c3     ->Write();
	//c4     ->Write();
	//c5     ->Write();
	//c6     ->Write();
	outfile->Close();	// Closing output file

	cout << "===================================" << endl;
	cout << "Notation in the output tree:\n----------------------------" << endl;
	cout << " labf = Laboratory frame" << endl;
	cout << " Crf  = Carbon rest frame" << endl;
	cout << "===================================" << endl;

	return 0;
}
// =============================================================================================================================







