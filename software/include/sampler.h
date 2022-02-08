#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#include "commondefs.h"
#include "resonances.h"
#include "b3d.h"

//#define __SAMPLER_WRITE_XY__

using namespace std;
class CSampler;

class CSampler{
public:
	CparameterMap *parmap;
	CRandy *randy;
	CResList *reslist;
	bool VISCOUSCORRECTIONS,TRIANGLE_FORMAT,JAKI_FORMAT,OSU_FORMAT;
	FILE *xyfptr;
	int NSAMPLE;

	CSampler(CB3D *b3d); // Constructor
	~CSampler();

	void ReadHyperElements2D_OSU();
	vector<CHyperElement> hyper;
	int GenHadronsFromHyperSurface();
	double cummulative_N,cummulative_random;
	double ETAMAX;
	int NBOSE;
	void SetPiByHand(double pixx,double pixy,double pixz,double piyy,double piyz,double pizz);
	double GetLambda(); // Calculates lambda in terms of T and densities
	void CalcPiFromParts();
	double Tf,epsilonf,Pf,sf,lambdaf,nhadronsf;
	Eigen::Matrix3d chif,chiinvf;
	vector<double> densityf,maxweightf;
	double GetLambda(double T,double P,double epsilon);
	CB3D *b3d;
	double MEANPT,MEANPT_NORM,NPI,NP;
	int nevents,ievent,nparts,nelements,nvertices;
	CLog *samplerlog;
	char message[500];
};

#endif
