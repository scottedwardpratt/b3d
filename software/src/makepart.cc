#include "constants.h"
#include "hyper.h"
#include "resonances.h"
#include "randy.h"
#include "misc.h"
#include "b3d.h"
#include "sampler.h"

using namespace std;

void CHyperElement::GetP(CResInfo *resinfo,FourVector &p,double &mass,double mw){
	bool VISCOUSCORRECTIONS=true;
	CResList *reslist=resinfo->reslist;
	bool reflect;
	CRandy *randy=(reslist->b3d)->randy;
	double pdotdOmega,nhatnorm,nhatdotp,wreflect;
	double pitilde[4][4];
	FourVector dOmegaprime,dOmega;
	int alpha,beta;
	FourVector pnoviscous,u;
	double m;
	FourVector nhat={0.0,0.0,0.0,0.0};
	pitilde[1][1]=pitildexx;
	pitilde[1][2]=pitilde[2][1]=pitildexy;
	pitilde[1][3]=pitilde[3][1]=0.0;
	pitilde[2][2]=pitildeyy;
	pitilde[2][3]=pitilde[3][2]=0.0;
	pitilde[3][3]=-pitildexx-pitildeyy;
	
	if(mw<0.0 || resinfo->width<1.0){
		m=resinfo->mass;
		randy->generate_boltzmann(m,T,pnoviscous);
	}
	else{
		m=resinfo->GenerateThermalMass(mw,T);
		randy->generate_boltzmann(m,T,pnoviscous);
	}
	mass=m;
	if(VISCOUSCORRECTIONS){
		for(alpha=1;alpha<4;alpha++){
			p[alpha]=pnoviscous[alpha];
			for(beta=1;beta<4;beta++){
				p[alpha]+=pitilde[alpha][beta]*pnoviscous[beta]/(h*lambda);
			}
		}
		p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+m*m);
	}
	else{
		for(alpha=0;alpha<4;alpha++)
			p[alpha]=pnoviscous[alpha];
	}
	u[1]=ux; u[2]=uy;	u[3]=0.0;
	u[0]=sqrt(1.0+ux*ux+uy*uy);
	dOmega[0]=dOmega0; dOmega[1]=dOmegaX; dOmega[2]=dOmegaY; dOmega[3]=0.0;
	Misc::BoostToCM(u,dOmega,dOmegaprime);  //dOmegaprime is dOmega in fluid (u=0) frame
	pdotdOmega=p[0]*dOmegaprime[0]-p[1]*dOmegaprime[1]-p[2]*dOmegaprime[2];
	//printf("u=(%g,%g), dOmega=(%g,%g,%g), dOmegaprime=(%g, %g,%g)\n",
	//ux,uy,dOmega[0],dOmega[1],dOmega[2],dOmegaprime[0],dOmegaprime[1],dOmegaprime[2]);
	
	wreflect=pdotdOmega/(p[0]*dOmegaprime[0]);
	reflect=false;
	if(wreflect<0.0)
		reflect=true;
	if(wreflect<1.0){
		if(wreflect<randy->ran())
			reflect=true;
	}
	if(reflect){
		nhatnorm=sqrt(dOmegaprime[1]*dOmegaprime[1]+dOmegaprime[2]*dOmegaprime[2]);
		nhat[1]=dOmegaprime[1]/nhatnorm;
		nhat[2]=dOmegaprime[2]/nhatnorm;
		nhatdotp=nhat[1]*p[1]+nhat[2]*p[2];
		p[1]-=2.0*nhat[1]*nhatdotp;
		p[2]-=2.0*nhat[2]*nhatdotp;
	}
	
	Misc::Boost(u,p);
	if(abs(resinfo->code)==2212 || abs(resinfo->code)==2112){
		sampler->MEANPT+=sqrt(p[1]*p[1]+p[2]*p[2]);
		sampler->MEANPT_NORM+=1.0;
		sampler->NP+=1.0;
	}
	if(abs(resinfo->code)==211)
		sampler->NPI+=1.0;
}
