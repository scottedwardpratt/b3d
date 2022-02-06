#include "sampler.h"
#include "hyper.h"
#include "constants.h"
#include "misc.h"
#include "randy.h"
#include "part.h"

using namespace std;
//CSampler* CvolumeElement2D::sampler=NULL;

CSampler::CSampler(CB3D *b3dset){
	printf("sampler check a\n");
	b3d=b3dset;
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
	printf("sampler check b\n");
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

/*
int CSampler::GenHadronsFromHyperSurface(){
double Omega0Sum=0.0;
nevents+=1;
int ielement,nparts=0;
for(ielement=0;ielement<nelements;ielement++){
Omega0Sum+=element[ielement].Omega[0];
nparts+=element[ielement].MakeParts();
}
printf("Event %4d sampling: %d parts created\n",nevents,int(b3d->PartMap.size()));
return nparts;
}
*/
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
	printf("Event %4d sampling: %d parts created from Hyper Surface\n",nevents-1,int(b3d->PartMap.size()));
	return nparts;
}

/*
void CSampler::CalcPiFromParts(){
	double pi[4][4]={{0.0}};
	CPartMap::iterator ppos;
	CPart *part;
	double pressure,Ptest=0.0,etest=0.0;
	int alpha,beta;
	double volume=nelements*element[0].Omega[0];
	for(ppos=b3d->DeadPartMap.begin();ppos!=b3d->PartMap.end();ppos++){
		part=ppos->second;
		pressure=(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3])/(3.0*part->p[0]);
		Ptest+=pressure;
		etest+=part->p[0];
		for(alpha=1;alpha<4;alpha++){
			pi[alpha][alpha]-=pressure;
			for(beta=1;beta<4;beta++){
				pi[alpha][beta]+=part->p[alpha]*part->p[beta]/part->p[0];
			}
		}
	}
	printf("----From particle momenta  -----\n");
	printf("epsilon=%10.4f, P=%10.4f\n",etest/volume,Ptest/volume);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pi[alpha][beta]=pi[alpha][beta]/volume;
			printf("%9.6f ",pi[alpha][beta]/element[0].P);
		}
		printf("\n");
	}
	printf("---- From pi_ij  input -------\n");
	printf("epsilon=%10.4f, P=%10.4f\n",element[0].epsilon,element[0].P);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%9.6f ",element[0].pitilde[alpha][beta]/element[0].P);
		}
		printf("\n");
	}
}*/


/*
void CSampler::ReadVolumeElements2D(){
	if(TRIANGLE_FORMAT){
		ReadVolumeElements2D_triangles();
		printf("triangle format\n");
	}
	else
	if(OSU_FORMAT){
		ReadHyperElements2D_OSU();
		printf("OSU format\n");
	}
	else{
		printf("Sampler format screwed up\n");
		exit(1);
	}
}
*/

/*	
void CSampler::ReadVolumeElements2D_triangles(){
	string filename;
	CvolumeElement2D *elem;
	double H;
	double dTxxoverH,dTxyoverH,dTxzoverH,dTyyoverH,dTyzoverH,dTzzoverH;
	int ivertex,iv1,iv2,iv3,ielement,initarraysize=1000;
	
	vertex.clear();
	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/vertices2D.txt";
	FILE *fptr=fopen(filename.c_str(),"r");
	
	fscanf(fptr,"%d",&ivertex);
	while(!feof(fptr)){
		if(ivertex==int(vertex.size()))
			vertex.resize(vertex.size()+initarraysize);
		fscanf(fptr,"%lf %lf %lf %lf %lf",&vertex[ivertex].r[0],&vertex[ivertex].r[1],&vertex[ivertex].r[2],
		&vertex[ivertex].ux,&vertex[ivertex].uy);
		fscanf(fptr,"%d",&ivertex);
	}
	fclose(fptr);
	nvertices=vertex.size();

	element.clear();
	nelements=0;
	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/triangles2D.txt";
	fptr=fopen(filename.c_str(),"r");
	
	ielement=0;
	fscanf(fptr,"%d %d %d",&iv1,&iv2,&iv3);
	double pi33overPbar=0.0,pi33overPbarnorm=0.0;
	while(!feof(fptr)){
		if(int(element.size())==ielement)
			element.resize(element.size()+initarraysize);
		elem=&element[ielement];
		fscanf(fptr,"%lf %lf %lf",
		&(elem->Omega[0]),&(elem->Omega[1]),&(elem->Omega[2]));

		fscanf(fptr,"%lf",&elem->T);
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf",
		&dTxxoverH,&dTxyoverH,&dTxzoverH,&dTyyoverH,&dTyzoverH,&dTzzoverH);
		elem->vertex[0]=&vertex[iv1];
		elem->vertex[1]=&vertex[iv2];
		elem->vertex[2]=&vertex[iv3];
		elem->epsilon=epsilonf;
		elem->density=&densityf;
		elem->P=Pf;
		H=Pf+epsilonf;
		elem->h=H;
		elem->lambda=lambdaf;
		elem->nhadrons=nhadronsf;
		//elem->FillOutShearTensor(H*dTxxoverH,H*dTxyoverH,H*dTxzoverH,H*dTyyoverH,H*dTyzoverH,H*dTzzoverH);
		elem->pitilde[1][1]=H*dTxxoverH;
		elem->pitilde[1][2]=elem->pitilde[2][1]=H*dTxyoverH;
		elem->pitilde[1][3]=elem->pitilde[3][1]=H*dTxzoverH;
		elem->pitilde[2][2]=H*dTyyoverH;
		elem->pitilde[2][3]=elem->pitilde[3][2]=H*dTyzoverH;
		elem->pitilde[3][3]=H*dTzzoverH;
		elem->CalcOmegaMax();
		pi33overPbar+=H*dTzzoverH*nhadronsf*elem->OmegaMax/Pf;
		pi33overPbarnorm+=nhadronsf*elem->OmegaMax;
		ielement+=1;
		fscanf(fptr,"%d %d %d",&iv1,&iv2,&iv3);
	}
	pi33overPbar=pi33overPbar/pi33overPbarnorm;
	nelements=ielement;
	printf("data read in, nelements=%d, pi33/P=%g\n",nelements,pi33overPbar);
}
*/

/*
void CSampler::ReadVolumeElements3D(){
	string filename;
	CvolumeElement2D *elem;
	double pixx, pixy, pixz, piyy, piyz, pizz;
	
	int ielement, initarraysize=1000;
	element.clear();
	ielement=0;
	nelements=0;

	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/hydro3D.txt";
	printf("\n***opening %s***\n\n", filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"r");

	while(!feof(fptr)){
		if(int(element.size())==ielement)
			element.resize(element.size()+initarraysize);
		elem=&element[ielement];

		fscanf(fptr,"%lf %lf %lf %lf",&elem->tau,&elem->eta,&elem->x,&elem->y);
		fscanf(fptr,"%lf %lf %lf",&(elem->ux),&(elem->uy),&(elem->uz));
		fscanf(fptr,"%lf %lf %lf %lf",&(elem->Omega[0]),&(elem->Omega[1]),&(elem->Omega[2]),&(elem->Omega[3]));
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf",&pixx,&pixy,&pixz,&piyy,&piyz,&pizz);
		
		printf("%lf %lf %lf %lf",elem->tau,elem->eta,elem->x,elem->y);
		printf("%lf %lf %lf",elem->ux,elem->uy,elem->uz);
		printf("%lf %lf %lf %lf",elem->Omega[0],elem->Omega[1],elem->Omega[2],elem->Omega[3]);
		printf("%lf %lf %lf %lf %lf %lf\n",pixx,pixy,pixz,piyy,piyz,pizz);
		
		elem->epsilon=epsilonf;
		elem->density=&densityf;
		elem->P=Pf;
		elem->lambda=lambdaf;
		elem->nhadrons=nhadronsf;

		elem->pitilde[1][1]=pixx;
		elem->pitilde[1][2]=elem->pitilde[2][1]=pixy;
		elem->pitilde[1][3]=elem->pitilde[3][1]=pixz;
		elem->pitilde[2][2]=piyy;
		elem->pitilde[2][3]=elem->pitilde[3][2]=piyz;
		elem->pitilde[3][3]=pizz;

		printf("[%lf %lf %lf]\n",elem->pitilde[1][1],elem->pitilde[1][2],elem->pitilde[1][3]);
		printf("[%lf %lf %lf]\n",elem->pitilde[2][1],elem->pitilde[2][2],elem->pitilde[2][3]);
		printf("[%lf %lf %lf]\n",elem->pitilde[3][1],elem->pitilde[3][2],elem->pitilde[3][3]);

		printf("ielement=%d\n",ielement);
		ielement+=1;
	}
	nelements=ielement;
}
*/

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
			/*
			dIpptest=0.0;
			for(p=0.5*dp;p<3000;p+=dp){
			e=sqrt(m*m+p*p);
			dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*( (2.0/3.0)*(p*p/e) - (2.0/15.0)*pow(p,4)/pow(e,3) );
			//dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*((2.0/15.0)*pow(p,4)/pow(e,3) );
			}
			Ipptest+=dIpptest;
			*/
			
			Ipp+=dIpp;
			//Ptest+=Ipptest;
			//printf("dIpptest=%g =? %g, ratio=%g\n",dIpptest,dIpp,dIpptest/dIpp);
		}
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);
	//printf("P=%g, epsilon=%g\n",P,epsilon);
	//printf("Ipp=%g =? %g, lambdafact=%g\n",Ipptest,2.0*P-4.0*Ipp,lambdafact);
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
	printf("opening %s\n",filename.c_str());
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
				printf("udotdOmega<0 =%g!!!\n",udotdOmega);
				udotdOmega=u0*dOmega0+ux*dOmegaX+uy*dOmegaY;
				printf("udotdOmega<0 =%g!!!\n",udotdOmega);
				printf("tau=%g, x=%g y=%g\n",tau,x,y);
				printf("ux=%g, uy=%g\n",ux,uy);
				printf("dOmega=(%g,%g,%g)\n",dOmega0,dOmegaX,dOmegaY);
				exit(1);
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
				printf("udotdOmega<0\n");
				elem->Print();
				Misc::Pause();
			}
			ielement+=1;
			b3d->TotalVolume+=udotdOmega;
		}
	}
	nelements=ielement;
	printf("Exiting ReadHyperElements2D_OSU() happily, TotalVolume=%g\n",b3d->TotalVolume);
}
