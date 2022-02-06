#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include "commondefs.h"
class CResList;

using namespace std;

typedef multimap<double,CResInfo *> CMassMap;

class CBranchInfo{
public:
	vector<CResInfo *> resinfoptr; //pointers for resinfo
	double branching;
	void Copy(CBranchInfo *branch);
	CBranchInfo();
};

class CMerge{
public:
	CMerge(CResInfo *resinfo,double branching, int L);
	CResInfo *resinfo;
	int L;
	double branching;
	CMerge *next;
};

class CResInfo{
public:
	int ires;
	double mass;
	double spin;
	double width;
	double netchi,netchi0;
	double minmass;
	string name;
	int code;
	int charge;
	int strange;
	int baryon;
	int q[3];
	int up,down;
	int G_Parity;
	bool decay; //false if stable, true if can decay. check if true
	CBranchList branchlist;
	CBranchList finalproductslist;
	CBranchInfo	*bptr_minmass;
	void Print();
	void DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	void DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	bool CheckForDaughters(int code);
	void FindFinalProducts(double taumax);
	void PrintFinalProducts();
	bool FindContent(int codecheck,double weight0,double taumax,double &weight);
	bool FindContentPairs(int codecheck1,int codecheck2,double weight0,double taumax,double &weight);
	bool CheckForNeutral();
	double GenerateMass();
	double GenerateThermalMass(double maxweight, double T);
	double ChiInt(double T,double vmax); // Integral used by ChiOmega
	double ChiTilde(double T,double vmax); // Integral used by ChiOmega
	CResInfo();
	static CRandy *randy;
	static CResList *reslist;
	static double **ChiA; // stored array used by ChiOmegaInt
};

class CResList{
public:
	CResList();
	~CResList();
	CResList(CparameterMap* parmap_in);
	CResInfoMap resmap;
	CMassMap massmap;
	CResInfo *GetResInfoPtr(int ID);
	void ReadResInfo();
	void CalcEoSandChi(double T,double &P,double &epsilon,double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi);
	void CalcEoSandChiandQdens(double T,double &P,double &epsilon,
double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi,double &udcontent,double &strangecontent);
	double CalcBalanceNorm(int pid,int pidprime,double taumax);
	void CalcConductivity(double T,double &P,double &epsilon,double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma);
	void freegascalc_onespecies(double m,double T,double &e,double &p,double &dens,double &sigma2,double &dedt);
	void freegascalc_onespecies_finitewidth(double m, double m1, double m2, double T,double width,double &epsilon,double &P,double &dens,double &sigma2,double &dedt, double &maxweight);
	double GetLambda(double T,double P,double epsilon);
	void freegascalc_onespecies(double m,double T,double &e,double &p,double &dens,double &sigma2,double &dedt,double &Jtot);
	void FindFinalProducts(double taumax);
	bool finalproductsfound;
	CparameterMap *parmap;
	CMerge ***MergeArray;
	double **SigmaMaxArray;
	double RESWIDTH_ALPHA;
	bool RESONANCE_DECAYS;
	bool USEPOLEMASS;
	//void freegascalc_onespecies_offshell(CResInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt);
	// freezeout props
	double Tf,epsilonf,Pf,lambdaf,nf;
	vector<double> densityf,maxweightf;
	Eigen::Matrix3d chif,chiinvf;
	double triangle(double m0,double m1,double m2);
	static CB3D *b3d;
	static CBalance *cb;
	static CSampler *sampler;
};

#endif
