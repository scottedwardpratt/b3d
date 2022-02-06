#ifndef __HYPER_H__
#define __HYPER_H__

#include "commondefs.h"

// Info for HyperSurface
class CHyperElement{
public:
	CHyperElement();
	~CHyperElement(){};
	double tau,x,y,dOmegaX,dOmegaY,dOmega0,ux,uy,T;
	vector<double> *density,*maxweight;
	double epsilon,P,h,nhadrons,lambda;
	double pitildexx,pitildexy,pitildeyy;
	double udotdOmega;
	void GetP(CResInfo *resinfo,FourVector &p,double &mass,double maxweight);
	void Copy(CHyperElement *oldhyper);
	int MakeParts();
	void Print();
	// Note dOmegaX and dOmegaY = d\Omega_\mu
	// (subscript mu, so p.Omega = p0*dOmega0+px*dOmegaX+py*dOmegaY)
	static CSampler *sampler;
	static CHydroBalance *hb;
};

#endif