#ifndef __BALANCE_ARRAYS_H__
#define __BALANCE_ARRAYS_H__

#include "commondefs.h"

using namespace std;

class CBFNumer{
public:
	string name;
	long long int npairs;
	vector<double> Bqinv,Bqout,Bqlong,Bqside,By,By1,Beta,Bphi,Betas,Beta1;  //q = 1/2 rel. momentum in CM frame
	vector<double> Cqinv,Cqout,Cqlong,Cqside,Cy,Cy1,Ceta,Cphi,Cetas,Ceta1;
	vector<vector<double>> Byphi,Cyphi; // Array for opp. signs, binned by both y and phi
	double Dy,Deta,Dqinv,Dphi;
	int Netabins,Nybins,Nqbins,Nphibins;
	CparameterMap *parmap;
	CBFNumer(CparameterMap *parmapset);
	void Reset();
	void Increment(CPart *parta,CPart *partb,double effa,double effb);
	void WriteNumer(string dirname,string numertype,bool NoQ);
	void Print();
	static CAcceptance *acceptance;
};

class CBFDenom{
public:
	double Nplus,Nminus;
	double dNdy;
	CparameterMap *parmap;
	CBFDenom(CparameterMap *parmapset);
	void Reset();
	void Increment(CPart *parta,double eff);
};

class CBalanceArrays{
public:
	CB3D *b3d;
	double NSAMPLE_HYDRO2UDS;
	int NSAMPLE_UDS2BAL,Nchi,NPHI,NEVENTS;
	bool FROM_UDS,NoKsNoPhi;  // FROM_UDS=true means BF from uds charges, if false, then brute forcd
	double BF_YMAX,BF_YMIN,DELY,BF_PHICUT;
	CPartMap bfpartmap;
	CResList *reslist;
	CparameterMap *parmap;
	string qualifier;
	string bf_results_dirname;
	CAcceptance *acceptance;
	string acceptance_description;
	CBalanceArrays(CB3D *b3dset);
	void Reset();
	void InitArrays();
	double gammap,gammapnorm,v2perfect,v2perfectnorm,v2,normtest,v2norm,v2prime,v2primenorm,dNdeta;
	double GAMMAP_SS,GAMMAP_OS;
	long long int NORM_GAMMAP_SS,NORM_GAMMAP_OS;
	CBFNumer *numer_pipi,*numer_piK,*numer_pip,*numer_KK;
	CBFNumer *numer_Kp,*numer_pp,*numer_allcharges;
	CBFDenom *denom_pi,*denom_K,*denom_p,*denom_allcharges;
	CBFNumer *bf_pipi,*bf_piK,*bf_pip,*bf_KK,*bf_Kp,*bf_pp,*bf_allcharges;
	CBFNumer *numer_allcharges_phi0,*numer_allcharges_phi45,*numer_allcharges_phi90;
	CBFNumer *bf_allcharges_phi0,*bf_allcharges_phi45,*bf_allcharges_phi90;
	CBFDenom *denom_allcharges_phi0,*denom_allcharges_phi45,*denom_allcharges_phi90;
	double rapdist[10];
	
	void ProcessPartMap();
	void ProcessBFPartMap();
	void ProcessV2Perfect();
	
	void IncrementNumer(CPart *parta,CPart *partb);
	void IncrementDenom(CPart *part);
	void IncrementGammaP(CPart *parta,CPart *partb,double effa,double effb);
	
	void SetQualifier(string qualifier_set);
	void PrintBFNumer();
	void PrintBFDenom();
	void CreateBFArrays();
	void ConstructBFs();
	void ConstructBF(CBFNumer *numer,CBFDenom *denom,CBFNumer *bf,double doublecount,bool NoQ);
	void WriteNumers();
	void WriteBFs();
	void WriteGammaP();
	void WriteDenoms();
	double GetMinv(CPart *parta,CPart *partb);
};

#endif