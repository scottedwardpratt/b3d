#ifndef __MUCALC_H__
#define __MUCALC_H__
#include "b3d.h"
#include "resonances.h"
#include "part.h"
using namespace std;

class CLocalInfo;

class CLocalSpeciesInfo{
public:
	CResInfo *resinfo;
	double degen;
	vector<int> pid;
	vector<double> T,ur,mu,Pr,E,M;
	vector<int> N;
	void Initialize(CResInfo *resinfo_set,string name_set);
	void Zero();
	void MuCalc();
	void MuCalc_PionsWithBose();
	void MuPrint();
	void Initialize();
	static CB3D *b3d;
};

class CLocalInfo{
public:
	static int NRBINS,NETEVENTS,IMUCALC;
	static double DELR2;
	static CB3D *b3d;
	static bool printing;
	multimap<int,CLocalSpeciesInfo *> mucalc_species;
	CLocalSpeciesInfo pion;
	CLocalSpeciesInfo kaon;
	CLocalSpeciesInfo nucleon;
	CLocalSpeciesInfo delta;
	CLocalSpeciesInfo lambda;
	CLocalSpeciesInfo sigma;
	CLocalSpeciesInfo sigmastar;
	CLocalSpeciesInfo xsi;
	CLocalSpeciesInfo xsistar;
	CLocalSpeciesInfo omega;
	void Initialize();
	void MuCalc();
};

#endif