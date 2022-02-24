#ifndef __MUTINFO_H__
#define __MUTINFO_H__
#include "commondefs.h"
using namespace std;


class CMuTInfo{
public:
	CMuTInfo(double tau_set);
	double tau,volume;
	double Pxpi,Pypi,Epi,PxK,PyK,EK,PxB,PyB,EB;
	double epsilonpi,epsilonK,epsilonB,rhopi,rhoK,rhoB,rhoBS;
	double Tpi,TK,TB,mupi,muK,muB,muBS;
	double Txxpi,Tyypi,Txypi;
	double TxxK,TyyK,TxyK;
	double TxxB,TyyB,TxyB;
	double Uxpi,Uypi,UxK,UyK,UxB,UyB;
	int Npi,NK,NB,NBS;
	bool sufficientN;
	
	void UpdateNPE(CB3DCell *cell);
	void CalcAllMuTU();
	static void GetMuT(double mass,double degen,double rho_target,double epsilon_target,
		double &T,double &mu);
	static bool GetMuT_Baryon(double rhoB_target,double rhoBS_target,double epsilon_target,
		double &T,double &muB,double &muBS);
	static void GetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
		double &Ux,double &Uy,double &epsilon);
	void Print();
	//void MuTCalc_PionsWithBose();
	static CB3D *b3d;
	static double DELTAU;
	static int NTAU;
	static int NETEVENTS;
	static int NMINCALC;
	static vector<CResInfo*> Bresinfo;
};

#endif