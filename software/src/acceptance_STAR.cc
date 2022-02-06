#include "acceptance.h"
#include "part.h"
#include "parametermap.h"
#include "resonances.h"

using namespace std;

CAcceptance_STAR::CAcceptance_STAR(CparameterMap *parmapin) : CAcceptance(){
	ETAMIN=-0.9; // Don't bother calling Acceptance Routine if outside these boundaries.
	ETAMAX=0.9;
	PTMIN=200.0; 
	PTMAX=1600.0;
	parmap=parmapin;
	CENTRALITY=parmap->getI("STARCENTRALITY",0);   // CENTRALITY=0 is most central
	//  cen=0 0-5%	
	//  cen=1 5-10%	
	//  cen=2 10-20%	
	//  cen=3 20-30%	
	//  cen=4 30-40%	
	//  cen=5 40-50%	
	//  cen=6 50-60%	
	//  cen=7 60-70%	
	//  cen=8 70-80%
	
	/*
	// initialize arrays for ID of non-identified (from Manuel Calderon)
	for (size_t i = 0; i<20; ++i)
	for (size_t j = 0; j<30; ++j) {
		ManuelData[i][j] = 0;
	}

	// read actual data from file
	char filename[100],identstring[100];
// Hi : 0 - 20%   (i.e. central events)
// Me : 20 - 50%
// Lo : 50 - 80%  (i.e. peripheral events)
	//Unfortunately, I only have the Hi CENTRALITY arrays from Manuel
	if(CENTRALITY==0 || CENTRALITY==1 || CENTRALITY==2)
		sprintf(identstring,"PiMinusHi");
	else if(CENTRALITY==3 || CENTRALITY==4 || CENTRALITY==5)
		sprintf(identstring,"PiMinusHi");
	else
		sprintf(identstring,"PiMinusHi");
  sprintf(filename,"../acceptancedata/manuel/full_field/efficiency%s.txt",identstring);
	ifstream ifs(filename);
	//assert(ifs);
	int eta_index, pt_index;
	double value, error;
	while (!ifs.eof()) {
		ifs >> eta_index >> pt_index >> value >> error;
		ManuelData[eta_index][pt_index] = value;

		// I'm not using an array for the error, but if needed, it can be
		// treated just like the ManuelData array.
	}
	ifs.close();
	*/
    
}

void CAcceptance_STAR::CalcAcceptance(bool &accept,double &efficiency,CPart *part){
	double eta,pt,pmag;
	double dca[4];
	int pid=part->resinfo->code,starpid;
	if(abs(pid)==211) starpid=1;
	else if(abs(pid)==321) starpid=2;
	else if(pid==-2212) starpid=3;
	else if(pid==2212) starpid=4;
	else{
		if(abs(pid)!=2112 && abs(pid)!=311 && abs(pid)!=111 && abs(pid)!=22){
			printf("CAcceptance_STAR::CalcAcceptance, pid=%d isn't in STAR list\n",pid);
			exit(1);
		}
		accept=false;
		efficiency=0.0;
		return;
	}
	part->CalcDCA(dca);
	accept=false;
	efficiency=0.0;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	if(eta>ETAMIN && eta<ETAMAX && dca[0]<2.0){
		if(pt>PTMIN && pt<PTMAX){
			accept=true;
			star_acc_eff(starpid,0.001*pt,eta,accept,efficiency);
		}
	}
}

void CAcceptance_STAR::CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part){
	double eta,pt,pmag;
	double dca[4];
	int pid=part->resinfo->code;
	if(abs(pid)!=2112 && abs(pid)!=311 && abs(pid)!=111 && abs(pid)!=22){
		printf("CAcceptance_STAR::CalcAcceptance, pid=%d isn't in STAR list\n",pid);
		exit(1);
	}
	part->CalcDCA(dca);
	accept=false;
	efficiency=0.0;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	if(eta>ETAMIN && eta<ETAMAX && pt>PTMIN && pt<PTMAX){
		accept=true;
	}
	efficiency=ScottEff(eta,pt);
}

// slightly modified for 2013
void CAcceptance_STAR::star_acc_eff(int pid,double pt,double eta,bool &accept,double &eff){
	int cen=CENTRALITY;
	// Routine to return acceptance and efficiency
	// Gary D. Westfall
	// January, 2014
	// Uses results published by STAR PRC 79 034909 2009, page 034909-13 for functional forms
	// Uses embedding results from Hui Wang for Run 10, 200 GeV
	// Determined TOF matching efficiency looking at pion spectra
	//  and double checked with kaons and protons
	// Input
	//  pid=1 pion (either sign)
	//  pid=2 kaon (with sign)
	//  pid=3 anti-proton
	//  pid=4 proton
	//  cen=0 0-5%	
	//  cen=1 5-10%	
	//  cen=2 10-20%	
	//  cen=3 20-30%	
	//  cen=4 30-40%	
	//  cen=5 40-50%	
	//  cen=6 50-60%	
	//  cen=7 60-70%	
	//  cen=8 70-80%	
	//  basic acceptance cuts are
	//   |eta| < 1.0 for TPC only
	//   |eta| < 0.9 for TPC+TOF
	//   full phi acceptance
	//  pt dependence is calculated in detail

	float f1;
	float P0,P1,P2,P3;
	float p,theta;
	// Pion parameters
	float P0_pion[9]={0.714,0.745,0.779,0.810,0.830,0.842,0.851,0.862,0.874};
	float P1_pion[9]={0.083,0.082,0.079,0.087,0.105,0.085,0.121,0.109,0.121};
	float P2_pion[9]={1.906,2.021,2.031,2.331,2.915,2.498,3.813,3.052,3.287};
	// Kaon parameters
	float P0_kaon[9]={0.634,0.646,0.691,0.711,0.748,0.742,0.783,0.758,0.748};
	float P1_kaon[9]={0.226,0.220,0.215,0.205,0.204,0.194,0.202,0.198,0.218};
	float P2_kaon[9]={1.387,1.490,1.447,1.533,1.458,1.485,1.372,1.559,1.665};
	float P3_kaon[9]={0.021,0.026,0.021,0.025,0.019,0.027,0.016,0.029,0.042};
	// pbar parameters
	float P0_pbar[9]={0.707,0.735,0.770,0.803,0.821,0.836,0.849,0.891,0.889};
	float P1_pbar[9]={0.191,0.185,0.196,0.194,0.191,0.200,0.199,0.156,0.197};
	float P2_pbar[9]={2.440,2.695,3.336,3.684,4.021,5.037,4.850,2.701,4.538};
	// p parameters
	float P0_p[9]={0.720,0.750,0.785,0.818,0.837,0.852,0.865,0.885,0.921};
	float P1_p[9]={0.201,0.197,0.201,0.193,0.202,0.202,0.192,0.192,0.178};
	float P2_p[9]={3.247,3.467,4.124,4.267,5.414,5.621,5.080,4.950,3.957};

	eff=0.0;

	// Check the input data
	
	if((pid < 1) || (pid > 4)){
		return;
	}
	if((pt < 0.2) || (pt > 3.0)){
		return;
	}
	if(fabs(eta) > 1.0){
		return;
	}
	if((cen < 0) || (cen > 8)){
		return;
	}
	
	// The basic parameters are OK, check the details
	
	// Pions
	if(pid == 1){
		theta=2.0*atan(exp(-eta));
		p=pt/sin(theta);
		if((pt > 0.2) && (p < 1.6)){
			accept=true;
			P0=P0_pion[cen];
			P1=P1_pion[cen];
			P2=P2_pion[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2));
			if(pt > 0.60){
				if(fabs(eta) < 0.9){
					eff=0.59*eff;
				} else {
					eff=0.0;
				}
			}
		}
	}

	// Kaons
	if(pid == 2){
		theta=2.0*atan(exp(-eta));
		p=pt/sin(theta);
		if((pt > 0.20) && (p < 1.6)){
			accept=true;
			P0=P0_kaon[cen];
			P1=P1_kaon[cen];
			P2=P2_kaon[cen];
			P3=P3_kaon[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2))+P3*pt;
			if(pt > 0.6){
				if(fabs(eta) < 0.9){
					eff=0.59*eff;
				} else {
					eff=0.0;
				}
			}
		}
	}

	// pbar
	if(pid == 3){
		theta=2.0*atan(exp(-eta));
		p=pt/sin(theta);
		if((pt > 0.4) && (p < 3.0)){
			accept=true;
			P0=P0_pbar[cen];
			P1=P1_pbar[cen];
			P2=P2_pbar[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2));
			if(pt > 1.0){
				if(fabs(eta) < 0.9){
					eff=0.59*eff;
				} else {
					eff=0.0;
				}
			}
		}
	}

	// p
	if(pid == 4){
		theta=2.0*atan(exp(-eta));
		p=pt/sin(theta);
		if((pt > 0.4) && (p < 3.0)){
			accept=true;
			P0=P0_p[cen];
			P1=P1_p[cen];
			P2=P2_p[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2));
			if(pt > 1.0){
				if(fabs(eta) < 0.9){
					eff=0.59*eff;
				} else {
					eff=0.0;
				}
			}
		}
	}
}

/*
double CAcceptance_STAR::ManuelEff(double eta, double pt) {
	if (fabs(eta)>1.0) {
		//cout << "Eta out of range " << eta << endl;
		return 0;
	}
	if (pt>3000 || pt<0) {
		//cout << "pT out of range " << pt << endl;
		return 0;
	}
	eta+=1.0;
	size_t eta_index = static_cast<int>(eta/0.1);
	size_t pt_index  = static_cast<int>(pt/100.0);
	if (eta_index>=20) {
		// this should not happen if the eta_index
		// is properly done! But let's be safe...
		cout << "Eta index out of range! " << eta_index << " " << eta << endl;
		return 0;
	}
	if (pt_index>=30) {
		// ditto
		cout << "pT index out of range! " << pt_index << " " << pt << endl;
		return 0;
	}
	return ManuelData[eta_index][pt_index];
}
*/

double CAcceptance_STAR::ScottEff(double eta, double pt) {
	// Very Cheap efficiency for non-indentified particles
	if (fabs(eta)>1.0) {
		//cout << "Eta out of range " << eta << endl;
		return 0;
	}
	if(pt>400.0){
		return 0.84;
	}
	else{
		if(pt<200.0) {
			return 0.0;
		}
		else{
			return 0.35+(0.84-0.35)*(pt-200.0)/200.0;
		}
	}
	 // approximate efficiency
}

