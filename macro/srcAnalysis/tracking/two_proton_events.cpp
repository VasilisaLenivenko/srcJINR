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

	v_proton_beam_UnSMEARED.Boost(-vc12_ForReconstructions.BoostVector()); // Boost from lab frame to Carbon rest frame

	// -----------------------------------------------------------------------------------------
	// Event selection cuts
	const float min_tof = 0.5; // ns
	const float max_tof = 25.; // ns
	
	// Graphical Cuts from George
	//TCutG *cut;
	//TFile *f = new TFile("theta_P2_boost_vs_theta_P1_boost_3p5GeVpercperu.root");
	//cut = (TCutG*)f->Get("mycut");

	//TCutG *cut1;
	//TFile *f1 = new TFile("P1_boost_vs_theta_P1_boost_3p5GeVpercperu.root");
	//cut1= (TCutG*)f1->Get("mycut");

	// -----------------------------------------------------------------------------------------
	// Defining histograms and other objects
	int min_t_L, min_t_R;

	// Hit coordinates in the TOF400 detectors
	TH2F * XY_hits_tof400 = new TH2F("XY_hits_tof400","XY hits tof400;X [cm];Y [cm]",200,-500,500,200,-80, 80);
	TH2F * XZ_hits_tof400 = new TH2F("XZ_hits_tof400","XZ hits tof400;X [cm];Z [cm]",200,-500,500,200,200,600);
	TH2F * YZ_hits_tof400 = new TH2F("YZ_hits_tof400","YZ hits tof400;Y [cm];Z [cm]",200,- 80, 80,200,200,600);

	int nhitL, nhitR;
	double tofhitL_x  [20], tofhitL_y  [20], tofhitL_z[20], tofhitL_t[20];
	double tofhitR_x  [20], tofhitR_y  [20], tofhitR_z[20], tofhitR_t[20];
	double thetaL_labf[20], thetaR_labf[20], phiL_labf[20], phiR_labf[20], dL[20], dR[20], ToFL, ToFR;
	double thetaL_Crf [20], thetaR_Crf [20], phiL_Crf [20], phiR_Crf [20];
	double Pmiss, Pmiss_x, Pmiss_y, Pmiss_z;
	double PL_labf[20], PL_labf_x[20], PL_labf_y[20], PL_labf_z[20], PR_labf[20], PR_labf_x[20], PR_labf_y[20], PR_labf_z[20];
	double PL_Crf [20], PL_Crf_x [20], PL_Crf_y [20], PL_Crf_z [20], PR_Crf [20], PR_Crf_x [20], PR_Crf_y [20], PR_Crf_z [20];
	double mand_s_Crf ,mand_u_Crf ,mand_t_Crf ;
	double mand_s_labf,mand_u_labf,mand_t_labf;

	TVector3 L_tof400_xyz, R_tof400_xyz;
	TLorentzVector P4L        [20], P4R        [20];
	TLorentzVector P4L_boosted[20], P4R_boosted[20];

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
	outtree->Branch("nhitL"      ,&nhitL      ,"nhitL/I"             );
        outtree->Branch("nhitR"      ,&nhitR      ,"nhitR/I"             );
	outtree->Branch("tofhitL_x"  , tofhitL_x  ,"tofhitL_x[nhitL]/D"  );
	outtree->Branch("tofhitL_y"  , tofhitL_y  ,"tofhitL_y[nhitL]/D"  );
	outtree->Branch("tofhitL_z"  , tofhitL_z  ,"tofhitL_z[nhitL]/D"  );
	outtree->Branch("tofhitL_t"  , tofhitL_t  ,"tofhitL_t[nhitL]/D"  );
	outtree->Branch("tofhitR_x"  , tofhitR_x  ,"tofhitR_x[nhitR]/D"  );
	outtree->Branch("tofhitR_y"  , tofhitR_y  ,"tofhitR_y[nhitR]/D"  );
	outtree->Branch("tofhitR_z"  , tofhitR_z  ,"tofhitR_z[nhitR]/D"  );
	outtree->Branch("tofhitR_t"  , tofhitR_t  ,"tofhitR_t[nhitR]/D"  );
	outtree->Branch("dL"         , dL         ,"dL[nhitL]/D"         );
	outtree->Branch("dR"         , dR         ,"dR[nhitR]/D"         );

	// Momenta in the lab frame
	outtree->Branch("PL_labf"    , PL_labf    ,"PL_labf[nhitL]/D"    );
	outtree->Branch("PL_labf_x"  , PL_labf_x  ,"PL_labf_x[nhitL]/D"  );
	outtree->Branch("PL_labf_y"  , PL_labf_y  ,"PL_labf_y[nhitL]/D"  );
	outtree->Branch("PL_labf_z"  , PL_labf_z  ,"PL_labf_z[nhitL]/D"  );
	outtree->Branch("PR_labf"    , PR_labf    ,"PR_labf[nhitR]/D"    );
	outtree->Branch("PR_labf_x"  , PR_labf_x  ,"PR_labf_x[nhitR]/D"  );
	outtree->Branch("PR_labf_y"  , PR_labf_y  ,"PR_labf_y[nhitR]/D"  );
	outtree->Branch("PR_labf_z"  , PR_labf_z  ,"PR_labf_z[nhitR]/D"  );

	// Momenta in the Carbon rest frame
	outtree->Branch("PL_Crf"     , PL_Crf     ,"PL_Crf[nhitL]/D"     );
	outtree->Branch("PL_Crf_x"   , PL_Crf_x   ,"PL_Crf_x[nhitL]/D"   );
	outtree->Branch("PL_Crf_y"   , PL_Crf_y   ,"PL_Crf_y[nhitL]/D"   );
	outtree->Branch("PL_Crf_z"   , PL_Crf_z   ,"PL_Crf_z[nhitL]/D"   );
	outtree->Branch("PR_Crf"     , PR_Crf     ,"PR_Crf[nhitR]/D"     );
	outtree->Branch("PR_Crf_x"   , PR_Crf_x   ,"PR_Crf_x[nhitR]/D"   );
	outtree->Branch("PR_Crf_y"   , PR_Crf_y   ,"PR_Crf_y[nhitR]/D"   );
	outtree->Branch("PR_Crf_z"   , PR_Crf_z   ,"PR_Crf_z[nhitR]/D"   );

	// Angles in the lab frame
	outtree->Branch("thetaL_labf", thetaL_labf,"thetaL_labf[nhitL]/D");
	outtree->Branch("thetaR_labf", thetaR_labf,"thetaR_labf[nhitR]/D");
	outtree->Branch("phiL_labf"  , phiL_labf  ,"phiL_labf[nhitL]/D"  );
	outtree->Branch("phiR_labf"  , phiR_labf  ,"phiR_labf[nhitR]/D"  );

	// Angles in the Carbon rest frame
	outtree->Branch("thetaL_Crf" , thetaL_Crf ,"thetaL_Crf[nhitL]/D" );
	outtree->Branch("thetaR_Crf" , thetaR_Crf ,"thetaR_Crf[nhitR]/D" );
	outtree->Branch("phiL_Crf"   , phiL_Crf   ,"phiL_Crf[nhitL]/D"   );
	outtree->Branch("phiR_Crf"   , phiR_Crf   ,"phiR_Crf[nhitR]/D"   );

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

	// -----------------------------------------------------------------------------
	// Set up the tree
	TClonesArray * tofData  = new TClonesArray("BmnTofHit"     );
	//TClonesArray * trigData = new TClonesArray("BmnTrigHit"    );
	//TClonesArray * mwpcData = new TClonesArray("BmnMwpcHit"    );
	//TClonesArray * gemData  = new TClonesArray("BmnGemStripHit");

	TTree * intree = NULL;
	intree = (TTree*) infile->Get("cbmsim");

	if (!intree){cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";return -3;}
	else{cerr << "Successfully loaded tree at address " << intree << "\n";}

	const int nEvents = intree->GetEntries();
	intree->SetBranchAddress("BmnTof1Hit"    ,&tofData );
	//intree->SetBranchAddress("BmnTrig"       ,&trigData);
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
	int ctr = 0;

	cout << "===================================" << endl;
	cout << "Beginning loop over events:" << endl;
	// Loop over events
	for (int event=0 ; event<nEvents ; event++){
		if (event % 10000 == 0) cerr << "Working on event " << event << " out of " << nEvents << endl;
		intree->GetEvent(event);	

		L_tof400_xyz.SetXYZ(0.,0.,0.);
		R_tof400_xyz.SetXYZ(0.,0.,0.);

		nhitL   = 0;
		nhitR   = 0;
		min_t_L = 0;
		min_t_R = 0;

		if(tofData->GetEntriesFast()==0) ctr++;

		for (int t=0 ; t<tofData->GetEntriesFast() ; t++){
			BmnTofHit * thisHit = (BmnTofHit*)tofData->At(t);	

			// --------------------------------------------------------------------------------------------------------
			// Hits on the left arm
			if ( (thisHit->GetX()>0.) && (thisHit->GetTimeStamp()>min_tof) && (thisHit->GetTimeStamp()<max_tof) ){

				if(nhitL==0) min_t_L = 0;
				else{ if(tofhitL_t[min_t_L]>thisHit->GetTimeStamp()) min_t_L = nhitL;}

				L_tof400_xyz.SetXYZ(thisHit->GetX() ,thisHit->GetY() ,(thisHit->GetZ() + shift_coord));
				XY_hits_tof400 -> Fill(L_tof400_xyz.X(),L_tof400_xyz.Y());
				XZ_hits_tof400 -> Fill(L_tof400_xyz.X(),L_tof400_xyz.Z());
				YZ_hits_tof400 -> Fill(L_tof400_xyz.Y(),L_tof400_xyz.Z());
				ToFL = thisHit -> GetTimeStamp();

				tofhitL_x  [nhitL] = L_tof400_xyz.X();
                		tofhitL_y  [nhitL] = L_tof400_xyz.Y();
                		tofhitL_z  [nhitL] = L_tof400_xyz.Z();
                		tofhitL_t  [nhitL] = ToFL;

				// Angles in the lab frame
				thetaL_labf[nhitL] = L_tof400_xyz.Theta();
               			phiL_labf  [nhitL] = L_tof400_xyz.Phi  ();

				// Path-length
				dL[nhitL] = sqrt(pow(L_tof400_xyz.X(),2) + pow(L_tof400_xyz.Y(),2) + pow(L_tof400_xyz.Z(),2))/100.; // [meter]
                		
				// Momentum components in the lab frame assuming proton
				PL_labf  [nhitL] = sqrt( pow(mp,2) / (pow(ToFL/(3.3*dL[nhitL]),2) - 1) );
                		PL_labf_x[nhitL] = PL_labf[nhitL]*cos(phiL_labf[nhitL])*sin(thetaL_labf[nhitL]);
                		PL_labf_y[nhitL] = PL_labf[nhitL]*sin(phiL_labf[nhitL])*sin(thetaL_labf[nhitL]);
                		PL_labf_z[nhitL] = PL_labf[nhitL]*cos(thetaL_labf[nhitL]);
                		
				// 4-vector momenta, in lab frame and in Carbon rest frame (boosted)
                		P4L        [nhitL].SetXYZM(PL_labf_x[nhitL],PL_labf_y[nhitL],PL_labf_z[nhitL],mp);
                		P4L_boosted[nhitL].SetXYZM(PL_labf_x[nhitL],PL_labf_y[nhitL],PL_labf_z[nhitL],mp);
                		P4L_boosted[nhitL].Boost(-vc12_ForReconstructions.BoostVector()); // Boosting to Carbon rest frame

				// Momentum components in the Carbon rest frame
                		PL_Crf  [nhitL] = P4L_boosted[nhitL].Rho();
                		PL_Crf_x[nhitL] = P4L_boosted[nhitL].X();
                		PL_Crf_y[nhitL] = P4L_boosted[nhitL].Y();
                		PL_Crf_z[nhitL] = P4L_boosted[nhitL].Z();
				
				// Angles in the Carbon rest frame
				thetaL_Crf[nhitL] = P4L_boosted[nhitL].Theta();	
			        phiL_Crf  [nhitL] = P4L_boosted[nhitL].Phi  ();

				nhitL++;
			}

			// --------------------------------------------------------------------------------------------------------
			// Hits on the right arm
			else if ( (thisHit->GetX()<0.) && (thisHit->GetTimeStamp()>min_tof) && (thisHit->GetTimeStamp()<max_tof) ){	

				if(nhitR==0) min_t_R = 0;
                                else{ if(tofhitR_t[min_t_R]>thisHit->GetTimeStamp()) min_t_R = nhitR;}

				R_tof400_xyz.SetXYZ(thisHit->GetX() ,thisHit->GetY() ,(thisHit->GetZ() + shift_coord));
				XY_hits_tof400 -> Fill(R_tof400_xyz.X(),R_tof400_xyz.Y());
				XZ_hits_tof400 -> Fill(R_tof400_xyz.X(),R_tof400_xyz.Z());
				YZ_hits_tof400 -> Fill(R_tof400_xyz.Y(),R_tof400_xyz.Z());
				ToFR = thisHit -> GetTimeStamp();
				
				tofhitR_x  [nhitR] = R_tof400_xyz.X();
                		tofhitR_y  [nhitR] = R_tof400_xyz.Y();
                		tofhitR_z  [nhitR] = R_tof400_xyz.Z();
                		tofhitR_t  [nhitR] = ToFR;

				// Angles in the lab frame
				thetaR_labf[nhitR] = R_tof400_xyz.Theta();
                		phiR_labf  [nhitR] = R_tof400_xyz.Phi  ();
                		if(phiR_labf[nhitR]<0) phiR_labf[nhitR] += 2*3.14159;
	
				// Path-length				
				dR[nhitR] = sqrt(pow(R_tof400_xyz.X(),2) + pow(R_tof400_xyz.Y(),2) + pow(R_tof400_xyz.Z(),2))/100.; // [meter]
                		
				// Momentum components in the lab frame assuming proton
				PR_labf  [nhitR] = sqrt( pow(mp,2) / (pow(ToFR/(3.3*dR[nhitR]),2) - 1) );
                		PR_labf_x[nhitR] = PR_labf[nhitR]*cos(phiR_labf[nhitR])*sin(thetaR_labf[nhitR]);
                		PR_labf_y[nhitR] = PR_labf[nhitR]*sin(phiR_labf[nhitR])*sin(thetaR_labf[nhitR]);
                		PR_labf_z[nhitR] = PR_labf[nhitR]*cos(thetaR_labf[nhitR]);
                		
				// 4-vector momenta, in lab frame and in Carbon rest frame (boosted)
                		P4R        [nhitR].SetXYZM(PR_labf_x[nhitR],PR_labf_y[nhitR],PR_labf_z[nhitR],mp);
                		P4R_boosted[nhitR].SetXYZM(PR_labf_x[nhitR],PR_labf_y[nhitR],PR_labf_z[nhitR],mp);
                		P4R_boosted[nhitR].Boost(-vc12_ForReconstructions.BoostVector()); // Boosting to Carbon rest frame

				// Momentum components in the Carbon rest frame
                		PR_Crf  [nhitR] = P4R_boosted[nhitR].Rho();
                		PR_Crf_x[nhitR] = P4R_boosted[nhitR].X();
                		PR_Crf_y[nhitR] = P4R_boosted[nhitR].Y();
                		PR_Crf_z[nhitR] = P4R_boosted[nhitR].Z();

				// Angles in the Carbon rest frame
                                thetaR_Crf[nhitR] = P4R_boosted[nhitR].Theta();
			        phiR_Crf  [nhitR] = P4R_boosted[nhitR].Phi  ();
				if(phiR_Crf[nhitR]<0) phiR_Crf [nhitR] += 2*3.14159;

				nhitR++;
			}
		}		
		
		// =========================================================================================================================
		// The momentum calculation below is given by:
		//        ________________      __________________________________________
		//       |     m^2             |                  m^2
		// P = /\| --------------  = /\| ----------------------------------------- , where 3.3 = 1/(speed of light [meters/nanosec])
		//       |  1/beta^2 - 1       |  (tof[nanosec]/(3.3*path[meters]))^2 - 1
		// =========================================================================================================================

		// Ignoring events that have more than two proton candidates
		// if ((nhitL!=1) || (nhitR!=1)) continue;

		// -------------------------------
		// Missing momentum
		Pmiss_x = P4L_boosted[min_t_L].X() + P4R_boosted[min_t_R].X();
		Pmiss_y = P4L_boosted[min_t_L].Y() + P4R_boosted[min_t_R].Y();
		Pmiss_z = P4L_boosted[min_t_L].Z() + P4R_boosted[min_t_R].Z() + v_proton_beam_UnSMEARED.P();
		Pmiss   = sqrt(pow(Pmiss_x,2) + pow(Pmiss_y,2) + pow(Pmiss_z,2));	

		// ----------------------------------------------------------------------------------
		// Calculating Mandelstam variables

		// Frame in which proton is moving, Carbon is at rest (Carbon rest frame)
		mand_s_Crf     = (P4L_boosted[min_t_L]+P4R_boosted[min_t_R])*(P4L_boosted[min_t_L]+P4R_boosted[min_t_R]);
		mand_u_Crf     = (v_proton_beam_UnSMEARED-P4R_boosted[min_t_R])*(v_proton_beam_UnSMEARED-P4R_boosted[min_t_R]);
		mand_t_Crf     = (v_proton_beam_UnSMEARED-P4L_boosted[min_t_L])*(v_proton_beam_UnSMEARED-P4L_boosted[min_t_L]);

		// Frame in which proton is standing, Carbon is moving (lab frame)
		mand_s_labf = (P4L[min_t_L]+P4R[min_t_R])*(P4L[min_t_L]+P4R[min_t_R]);
		mand_u_labf = (v_proton_beam_Standing-P4R[min_t_R])*(v_proton_beam_Standing-P4R[min_t_R]);
		mand_t_labf = (v_proton_beam_Standing-P4L[min_t_L])*(v_proton_beam_Standing-P4L[min_t_L]);

		// Mandelstam variables are Lorentz-invariant, so the values above should be
		// the same before and after the boost. Let's make sure that's the case:
		if(             (mand_s_Crf-mand_s_labf>0.001)||
				(mand_u_Crf-mand_u_labf>0.001)||
				(mand_t_Crf-mand_t_labf>0.001)
		  ) cout << "There seems to be a problem with the Mandelstam variables. Check it!!!" << endl;

		// Applying graphical cuts
		//if (!(cut->IsInside(thetaR_labf,thetaL_labf)))   continue;
		//if (!(cut1->IsInside(PL_labf,thetaL_labf)))      continue;
		//if (!(cut1->IsInside(PR_labf,thetaR_labf)))      continue;
		//if ((PL_Crf+PR_Crf)<4.4 || (PL_Crf+PR_Crf)<4.8)  continue; 
		
		// ----------------------------------------------------------------------------------
		// Saving informormation onto output tree
		if((nhitL!=0)&&(nhitR!=0)){
			outtree->Fill();
		}
	}

	TCanvas * c1 = new TCanvas("c1");	XY_hits_tof400 -> Draw("COLZ");
	TCanvas * c2 = new TCanvas("c2");       XZ_hits_tof400 -> Draw("COLZ");
	TCanvas * c3 = new TCanvas("c3");       YZ_hits_tof400 -> Draw("COLZ");

	// ----------------------------------------------------------------------------------
	// Writing information onto output file

	infile ->Close();	// Closing input file
	outtree->Write();	// Writing tree onto output file
	c1     ->Write();
	c2     ->Write();
	c3     ->Write();
	outfile->Close();	// Closing output file

	cout << "===================================" << endl;
	cout << "Notation in the output tree:\n----------------------------" << endl;
	cout << " labf = Laboratory frame" << endl;
	cout << " Crf  = Carbon rest frame" << endl;
	cout << "===================================" << endl;

	cout << "Number of events in which there were no tof400 hits: " << ctr << endl << endl;

	return 0;
}
// =============================================================================================================================







