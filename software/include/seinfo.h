#ifndef __SEINFO_H__
#define __SEINFO_H__

#include "commondefs.h"

using namespace std;

class CB3D;

class CSEInfo{
public:
	CSEInfo(CB3D *b3dset);
	// vectors hold information for different times
	vector<double> Pbar,Tzz,epsilon,nhadrons,K0,F0,uperpbar;
	// this is information for latest time
	void Zero(); // sets arrays to zero
	void SECalc();
	void Print();
	double R,TAU0,TAUMAX,DELTAU,ETAOVERS,RMAX;
	int NTAU;
	int NETEVENTS;
	CB3D *b3d;
};

#endif