#include "sampler.h"
#include "hyper.h"
#include "constants.h"
#include "misc.h"
#include "randy.h"
#include "part.h"

using namespace std;

CSampler::CSampler(CB3D *b3dset){
	b3d=b3dset;
	samplerlog=b3d->b3dlog;
	parmap=&(b3d->parmap);
	xyfptr=NULL;
#ifdef __SAMPLER_WRITE_XY__
	xyfptr=fopen("xy.txt","w");
#endif
	VISCOUSCORRECTIONS=true;
	TRIANGLE_FORMAT=false;
	//CvolumeElement2D::sampler=this;
	CHyperElement::sampler=this;
	randy=b3d->randy;
	reslist=b3d->reslist;
	cummulative_N=0.0;
	nevents=0;
	cummulative_random=-log(randy->ran());
	Tf=b3d->parmap.getD("FREEZEOUT_TEMP",155);
	ETAMAX=b3d->parmap.getD("B3D_ETAMAX",1.0);
	NBOSE=b3d->parmap.getI("B3D_NBOSE",1);
	NSAMPLE=b3d->parmap.getI("B3D_NSAMPLE",1);
	TRIANGLE_FORMAT=b3d->parmap.getB("SAMPLER_TRIANGLE_FORMAT",false);
	OSU_FORMAT=b3d->parmap.getB("SAMPLER_OSU_FORMAT",true);
	densityf.clear();
	maxweightf.clear();
	densityf.resize(reslist->resmap.size());
	maxweightf.resize(reslist->resmap.size());
	reslist->CalcEoSandChi(Tf,Pf,epsilonf,nhadronsf,densityf,maxweightf,chif);
	chiinvf=chif.inverse();
	sf=(epsilonf+Pf)/Tf;
	lambdaf=GetLambda(Tf,Pf,epsilonf);
	MEANPT=MEANPT_NORM=0.0;
}

CSampler::~CSampler(){
	if(xyfptr!=NULL){
		fclose(xyfptr);
		xyfptr=NULL;
	}
}

int CSampler::GenHadronsFromHyperSurface(){
	double Omega0Sum=0.0,udotdOmegaSum=0.0;
	nevents+=1;
	int ielement;
	nparts=0;
	for(ielement=0;ielement<nelements;ielement++){
		Omega0Sum+=hyper[ielement].dOmega0;
		udotdOmegaSum+=hyper[ielement].udotdOmega;
		nparts+=hyper[ielement].MakeParts();
	}
	sprintf(message,"Event %4d sampling: %d parts created from Hyper Surface\n",nevents-1,int(b3d->PartMap.size()));
	samplerlog->Info(message);
	return nparts;
}

double CSampler::GetLambda(double T,double P,double epsilon){
	int n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp,J,nfact,sign,alpha;
	//double dIpptest=0.0,dp=4.0,p,e,Ipptest=0.0,Ptest=0.0;
	double lambdafact;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			dIpp=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			
			Ipp+=dIpp;
		}
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);
	return lambdafact;
}

void CSampler::ReadHyperElements2D_OSU(){
	string filename;
	CHyperElement *elem;
	//double PIbulk,dOmegaMax;
	double u0,ux,uy,x,y,udotdOmega,dOmega0,dOmegaX,dOmegaY,pitildexx,pitildeyy,pitildexy,tau,dummy1,dummy2;
	int ielement,initarraysize=1000;
	char dummy[300];
	hyper.clear();
	nelements=0;
	b3d->TotalVolume=0.0;
	bool oldhyperformat=parmap->getB("B3D_OLDHYPERFORMAT",false);
	string hyperdata_home=parmap->getS("B3D_HYPERDATA_HOME","hyperdata");
	filename=hyperdata_home+"/"+b3d->qualifier+"/hyper.txt";
	sprintf(message,"opening %s\n",filename.c_str());
	samplerlog->Info(message);
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,200,fptr);	fgets(dummy,200,fptr);
	ielement=0;
	while(!feof(fptr)){
		if(int(hyper.size())==ielement)
			hyper.resize(hyper.size()+initarraysize);
		elem=&hyper[ielement];
		elem->T=b3d->parmap.getD("FREEZEOUT_TEMP",155);
		if(oldhyperformat){
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&tau,&x,&y,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,&dummy1,&dummy2,&pitildexx,&pitildeyy,&pitildexy);
			dOmegaX=-dOmegaX; dOmegaY=-dOmegaY;
		}
		else{
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&tau,&x,&y,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,&pitildexx,&pitildeyy,&pitildexy);
		}
		if(!feof(fptr)){
			u0=sqrt(1.0+ux*ux+uy*uy);
			udotdOmega=u0*dOmega0-ux*dOmegaX-uy*dOmegaY;
			if(udotdOmega<0.0){
				sprintf(message,"udotdOmega<0 =%g!!!\n",udotdOmega);
				udotdOmega=u0*dOmega0+ux*dOmegaX+uy*dOmegaY;
				sprintf(message,"%sudotdOmega<0 =%g!!!\n",message,udotdOmega);
				sprintf(message,"%stau=%g, x=%g y=%g\n",message,tau,x,y);
				sprintf(message,"%sux=%g, uy=%g\n",message,ux,uy);
				sprintf(message,"%sdOmega=(%g,%g,%g)\n",message,dOmega0,dOmegaX,dOmegaY);
				samplerlog->Fatal(message);
			}
			elem->tau=tau;
			elem->dOmega0=dOmega0*2.0*b3d->ETAMAX;
			elem->dOmegaX=dOmegaX*2.0*b3d->ETAMAX;
			elem->dOmegaY=dOmegaY*2.0*b3d->ETAMAX;
			elem->udotdOmega=udotdOmega*2.0*b3d->ETAMAX;
		
			elem->x=x;
			elem->y=y;
			elem->ux=ux;
			elem->uy=uy;
			elem->pitildexx=1000.0*pitildexx;
			elem->pitildeyy=1000.0*pitildeyy;
			elem->pitildexy=1000.0*pitildexy;

			elem->epsilon=epsilonf;
			elem->density=&densityf;
			elem->T=Tf;
			elem->P=Pf;
			elem->h=elem->P+elem->epsilon;
			elem->lambda=lambdaf;
			elem->nhadrons=nhadronsf;
			elem->maxweight=&maxweightf;
			if(udotdOmega<0.0){
				sprintf(message,"udotdOmega<0\n");
				elem->Print(samplerlog);
				samplerlog->Fatal(message);
			}
			ielement+=1;
			b3d->TotalVolume+=udotdOmega;
		}
	}
	nelements=ielement;
	sprintf(message,"Exiting ReadHyperElements2D_OSU() happily, TotalVolume=%g\n",b3d->TotalVolume);
	samplerlog->Info(message);
}
