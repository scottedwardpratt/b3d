#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#include "commondefs.h"
#include "resonances.h"
#include "b3d.h"

//#define __SAMPLER_WRITE_XY__

using namespace std;
class CSampler;
/*
class Cvertex2D;
class CvolumeElement2D{
public:
	double T,muB,muE,muS,XUD,XS,Xscale;
	double nhadrons,nbaryonsS0,nbaryonsS1,nbaryonsS2,nbaryonsS3;
	double nmesonsS0,nmesonsS1,nmesonsS2;
	vector<double> *density; // densities of various species
	FourVector Omega; // volume elements (four-vectors)
	double OmegaMax,udotdOmega;
	Cvertex2D *vertex[3]; //vertices of triangle
	double x,y,tau,eta,ux,uy,uz;
	double **pitilde; // shear tensor
	CvolumeElement2D();
	void Initialize();
	void CalcEquilibriumQuantities();
	void CopyBulkQuantities(CvolumeElement2D *element);
	void CopyEquilibriumQuantities(CvolumeElement2D *element);
	double P,epsilon,h,lambda; // prefactor used for viscous corrections
	void FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz);
	void Print();
	void CalcOmegaMax();
	void CalcOmegaMax3D();
	int MakeParts();
	static CSampler *sampler;
};

class Cvertex2D{
public:
	double r[3],ux,uy;
	void Print(){
		printf("vertex: r=(%g,%g,%g), u=(%g,%g,%g)\n",r[0],r[1],r[2],sqrt(1.0+ux*ux+uy*uy),ux,uy);
	};
};
*/

class CSampler{
public:
	CparameterMap *parmap;
	CRandy *randy;
	CResList *reslist;
	//vector<CvolumeElement2D> volume_element;
	//vector<Cvertex2D> vertex;
	bool VISCOUSCORRECTIONS,TRIANGLE_FORMAT,JAKI_FORMAT,OSU_FORMAT;
	FILE *xyfptr;
	int NSAMPLE;

	CSampler(CB3D *b3d); // Constructor
	~CSampler();

	/*
	void ReadVolumeElements2D();
	void ReadVolumeElements2D_triangles();
	void ReadVolumeElements2D_Jaki();
	void ReadVolumeElements2D_center();
	*/
	void ReadHyperElements2D_OSU();
	//void ReadVolumeElements3D();
	//vector<CvolumeElement2D> element;
	vector<CHyperElement> hyper;
	int GenHadronsFromHyperSurface();
	double cummulative_N,cummulative_random;
	double ETAMAX;
	int NBOSE;
	//void SetVolumeElementsByHand(double T,int nelements_set,double elementvg);
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
};

#endif
