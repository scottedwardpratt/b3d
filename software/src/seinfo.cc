#include "seinfo.h"
#include "b3d.h"
#include "part.h"
#include "resonances.h"
#include "constants.h"
#include "sampler.h"

CSEInfo::CSEInfo(CB3D *b3dset){
	b3d=b3dset;
	DELTAU=1.0;
	RMAX=b3d->parmap.getD("B3D_XYMAX",12);
	TAUMAX=b3d->parmap.getD("B3D_TAUCOLLMAX",30);
	TAU0=b3d->parmap.getD("SEINFO_TAU0",5);
	NTAU=lrint(TAUMAX-TAU0);
	ETAOVERS=0.0;
	NETEVENTS=0;
	R=RMAX-(TAUMAX-TAU0);
	sprintf(message,"For SEINFO: TAU0=%g, TAUMAX=%g, RMAX=%g, R=%g, DELTAU=%g\n",TAU0,TAUMAX,RMAX,R,TAUMAX-TAU0);
	CLog::Info(message);
	epsilon.resize(NTAU+1);
	Tzz.resize(NTAU+1);
	Pbar.resize(NTAU+1);
	nhadrons.resize(NTAU+1);
	uperpbar.resize(NTAU+1);
	K0.resize(NTAU+1);
	F0.resize(NTAU+1);
	Zero();
}

void CSEInfo::Zero(){
	int itau;
	NETEVENTS=0;
	for(itau=0;itau<=NTAU;itau++){
		Tzz[itau]=Pbar[itau]=epsilon[itau]=nhadrons[itau]=K0[itau]=F0[itau]=0.0;
	}
}

void CSEInfo::SECalc(){
	double e,pmag2,pz,et,r,volume;
	CPartMap::iterator ppos;
	CPart *part;
	int itau;
	itau=lrint((b3d->tau-TAU0)/DELTAU);
	for(ppos=b3d->PartMap.begin();ppos!=b3d->PartMap.end();ppos++){
		part=ppos->second;
		part->Propagate(b3d->tau);
		r=sqrt(part->r[1]*part->r[1]+part->r[2]*part->r[2]);
		if(r<R){
			pmag2=part->p[1]*part->p[1]+part->p[2]*part->p[2];
			et=sqrt(part->msquared+pmag2);
			pz=et*sinh(part->y-part->eta);
			e=sqrt(et*et+pz*pz);
			pmag2+=pz*pz;
			Tzz[itau]+=pz*pz/e;
			Pbar[itau]+=pmag2/(3.0*e);
			epsilon[itau]+=e;
			K0[itau]+=pmag2*(1.0-0.2*pmag2/(e*e))/(3.0*e);
			F0[itau]+=pmag2*pmag2/(15.0*e*e);
			nhadrons[itau]+=1;
			uperpbar[itau]+=(part->p[1]*part->r[1]+part->p[2]*part->r[2])/R;
		}
	}
	volume=2.0*b3d->ETAMAX*b3d->tau*PI*R*R*NETEVENTS;
	sprintf(message,"tau=%5.2f, pizz=%9.6f, <uperp>=%g\n",b3d->tau,-(Pbar[itau]-Tzz[itau])/volume,uperpbar[itau]/nhadrons[itau]);
	CLog::Info(message);
}

void CSEInfo::Print(){
	int itau;
	double eta;
	double alpha,dpizzoveralpha_dt,p0=0.0,p1=0.0,p2=0.0;
	double tau,volume;
	char filename[150];
	sprintf(filename,"model_output/%s/%s/tij.txt",b3d->run_name.c_str(),b3d->qualifier.c_str());
	FILE *fptr=fopen(filename,"w");
	sprintf(message,"#tau  epsilon   Pbar    Tzz   nhadrons  eta      K0      F0     -pizz   alpha\n");
	CLog::Info(message);
	fprintf(fptr,"#tau  epsilon   Pbar    Tzz   nhadrons  eta      K0      F0     -pizz    alpha\n");
	for(itau=0;itau<=NTAU;itau++){
		tau=TAU0+itau*DELTAU;
		volume=2.0*b3d->ETAMAX*tau*PI*R*R*NETEVENTS;
		eta=0.75*tau*(Pbar[itau]-Tzz[itau])/volume;
		alpha=sqrt(TAU0/tau)*sqrt(F0[itau]/F0[0]);
		sprintf(message,"%5.2f %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %7.2f %7.4f %7.4f\n",
		tau,epsilon[itau]/volume,Pbar[itau]/volume,Tzz[itau]/volume,nhadrons[itau]/volume,eta,K0[itau]/volume,F0[itau]/volume,(Pbar[itau]-Tzz[itau])/volume,alpha);
		CLog::Info(message);
		fprintf(fptr,"%5.2f %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %7.2f %7.4f %7.4f\n",
		tau,epsilon[itau]/volume,Pbar[itau]/volume,Tzz[itau]/volume,nhadrons[itau]/volume,eta,K0[itau]/volume,F0[itau]/volume,(Pbar[itau]-Tzz[itau])/volume,alpha);
		if(itau==0)
			p0=(Pbar[itau]-Tzz[itau])/alpha;
		if(itau==1)
			p1=(Pbar[itau]-Tzz[itau])/alpha;
		if(itau==2)
			p2=(Pbar[itau]-Tzz[itau])/alpha;
	}
	dpizzoveralpha_dt=(2.0*p1-1.5*p0-0.5*p2)/DELTAU;
	sprintf(message,"&&&&&&&&&eta/s=%g D(pizz/alpha)/Dt=%g &&&&&&&&&&&&\n",ETAOVERS,dpizzoveralpha_dt);
	CLog::Info(message);
	fprintf(fptr,"#&&&&&&&&&eta/s=%g D(pizz/alpha)/Dt=%g &&&&&&&&&&&&\n",ETAOVERS,dpizzoveralpha_dt);
	fclose(fptr);
}
