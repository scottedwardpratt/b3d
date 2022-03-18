#include "hyper.h"
#include "sampler.h"
#include "part.h"
#include "randy.h"
#include "misc.h"

//#define __XYREFLECT__

CSampler* CHyperElement::sampler=NULL;
char *CHyperElement::message=new char[500];

CHyperElement::CHyperElement(){
	//
}

void CHyperElement::Copy(CHyperElement *oldhyper){
	tau=oldhyper->tau;
	x=oldhyper->x; y=oldhyper->y;
	dOmegaX=oldhyper->dOmegaX;
	dOmegaY=oldhyper->dOmegaY;
	dOmega0=oldhyper->dOmega0;
	ux=oldhyper->ux;
	uy=oldhyper->uy;
	udotdOmega=oldhyper->udotdOmega;
	pitildexx=oldhyper->pitildexx;
	pitildexy=oldhyper->pitildexy;
	pitildeyy=oldhyper->pitildeyy;
	T=oldhyper->T;
}

void CHyperElement::Print(){
	sprintf(message,"HyperElement Info:\n");
	sprintf(message,"%stau=%g, T=%g, x=%g, y=%g, ux=%g, uy=%g, dOmega0=%g, dOmegaX=%g, dOmegaY=%g, udotdOmega=%g\n",
	message,tau,T,x,y,ux,uy,dOmega0,dOmegaX,dOmegaY,udotdOmega);
	sprintf(message,"%spitildexx=%g, pitildeyy=%g, pitildexy=%g\n",
	message,pitildexx,pitildeyy,pitildexy);
	sprintf(message,"%s---------------------------------------------------------------\n",message);
	CLog::Info(message);
}

int CHyperElement::MakeParts(){
	int nparts=0;
	CPart *part;
	double bweight,mass,eta,rapidity;
	int ires,nsample=sampler->NSAMPLE;
	FourVector plab;
	double delN,r[3];
	double delNtot=nhadrons*udotdOmega*nsample;
	CResInfoMap *resmap=&(sampler->reslist->resmap);
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	CRandy *randy=sampler->randy;
	double ETAMAX=sampler->ETAMAX;
	if(sampler->cummulative_N+delNtot > sampler->cummulative_random){
		for(rpos=resmap->begin();rpos!=resmap->end();rpos++){
			resinfo=rpos->second;
			if(resinfo->code==22){
				rpos++;
				resinfo=rpos->second;
			}
			ires=resinfo->ires;
			delN=(*density)[ires]*udotdOmega*nsample;
			sampler->cummulative_N+=delN;
			while(sampler->cummulative_N>sampler->cummulative_random){
				//printf("howdy, code=%d\n",resinfo->code);
				if(abs(resinfo->code)==3112)
					sampler->b3d->ncheck+=1;
				if(abs(resinfo->code)==3212){
					sampler->b3d->ncheck1+=1;
					resinfo->Print();
				}
				if(abs(resinfo->code)==3222)
					sampler->b3d->ncheck2+=1;
				GetP(resinfo,plab,mass,(*maxweight)[ires]);
				part=sampler->b3d->GetDeadPart();	
#ifdef __XY_REFLECT__	
				if(randy->ran()<0.5){
					r[1]=-r[1];
					plab[1]=-plab[1];
				}
				if(randy->ran()<0.5){
					r[2]=-r[2];
					plab[2]=-plab[2];
				}
#endif
				r[0]=tau; r[1]=x; r[2]=y;
				bweight=1.0;
				rapidity=atanh(plab[3]/plab[0]);
				eta=(1.0-2.0*randy->ran())*ETAMAX;
				rapidity+=eta;
				if(fabs(eta)>ETAMAX){
					sprintf(message,"eta=%g ??, ETAMAX=%g\n",eta,ETAMAX);
					CLog::Fatal(message);
				}
				part->InitBalance(resinfo->code,r[1],r[2],r[0],eta,plab[1],plab[2],mass,rapidity,bweight,-1);
				/*
				if(part->balanceID<0){
					for(int a=0;a<3;a++){
						for(int b=0;b<3;b++){
							sampler->chitot(a,b)+=resinfo->q[a]*resinfo->q[b];
						}
					}
				}
				*/
				sampler->cummulative_random-=log(randy->ran());
				nparts+=1;
			}
		}	
	}
	else
		sampler->cummulative_N+=delNtot;
	return nparts;
}
