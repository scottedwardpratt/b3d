#include "balancearrays.h"
#include "part.h"
#include "resonances.h"
#include "parametermap.h"
#include "misc.h"
#include "acceptance.h"
#include "constants.h"

using namespace std;

CLog *CBFNumer::barrayslog=NULL;
CLog *CBFDenom::barrayslog=NULL;
char *CBFNumer::message=NULL;
char *CBFDenom::message=NULL;
CAcceptance *CBFNumer::acceptance=NULL;

CBFNumer::CBFNumer(CparameterMap *parmapset){
	int iy;
	npairs=0;
	parmap=parmapset;
	Netabins=parmap->getD("BF_NETABINS",50);
	Nybins=parmap->getD("BF_NYBINS",50);
	Nqbins=parmap->getD("BF_NQINVBINS",100);
	Nphibins=parmap->getD("BF_NPHIBINS",180);
	Deta=parmap->getD("BF_DETA",0.1);
	Dy=parmap->getD("BF_DELY",0.1);
	Dqinv=parmap->getD("BF_DQINV",10);
	Dphi=360.0/double(Nphibins);
	
	Bqinv.resize(Nqbins,0);
	Bqout.resize(Nqbins,0);
	Bqside.resize(Nqbins,0);
	Bqlong.resize(Nqbins,0);
	Beta.resize(Netabins,0);
	Beta1.resize(Netabins,0);
	By.resize(Nybins,0);
	By1.resize(Nybins,0);
	Betas.resize(Netabins,0);
	Bphi.resize(Nphibins,0);
	
	Cqinv.resize(Nqbins,0);
	Cqout.resize(Nqbins,0);
	Cqside.resize(Nqbins,0);
	Cqlong.resize(Nqbins,0);
	Ceta.resize(Netabins,0);
	Ceta1.resize(Netabins,0);
	Cy.resize(Nybins,0);
	Cy1.resize(Nybins,0);
	Cetas.resize(Netabins,0);
	Cphi.resize(Nphibins,0);
	
	Byphi.resize(Nybins);
	Cyphi.resize(Nybins);
	for(iy=0;iy<Nybins;iy++){
		Byphi[iy].resize(Nphibins,0);
		Cyphi[iy].resize(Nphibins,0);
	}
}

void CBFNumer::Reset(){
	Bqinv.assign(Nqbins,0);
	Bqout.assign(Nqbins,0);
	Bqside.assign(Nqbins,0);
	Bqlong.assign(Nqbins,0);
	Beta.assign(Netabins,0);
	Beta1.assign(Netabins,0);
	By.assign(Nybins,0);
	By1.assign(Nybins,0);
	Bphi.assign(Nphibins,0);
	Betas.assign(Netabins,0);
	Cqinv.assign(Nqbins,0);
	Cqout.assign(Nqbins,0);
	Cqside.assign(Nqbins,0);
	Cqlong.assign(Nqbins,0);
	Ceta.assign(Netabins,0);
	Ceta1.assign(Netabins,0);
	Cy.assign(Nybins,0);
	Cy1.assign(Nybins,0);
	Cphi.assign(Nphibins,0);
	Cetas.assign(Netabins,0);
	for(int iy=0;iy<Nybins;iy++){
		Byphi[iy].assign(Nphibins,0);
		Cyphi[iy].assign(Nphibins,0);
	}
}

void CBFNumer::Increment(CPart *parta,CPart *partb,double effa,double effb){
	int ibin,iphi,iy,pida=parta->resinfo->code,pidb=partb->resinfo->code;
	double qinv,qout,qside,qlong,deleta,dely,delphi,deletas;
	double QaQb,CaCb,cphi;
	Misc::outsidelong(parta->p,partb->p,qinv,qout,qside,qlong,deleta,dely,delphi);
	cphi=cos(delphi*PI/180.0);
	qout=fabs(qout); qside=fabs(qside); qlong=fabs(qlong);
	QaQb=(parta->resinfo->charge*parta->bweight)*(partb->resinfo->charge*partb->bweight);
	QaQb*=effa*effb;
	CaCb=parta->bweight*partb->bweight*effa*effb;
	deletas=fabs(parta->eta-partb->eta);
	npairs+=1;
	
	if(dely<0.0 || deleta<0.0 || qout<0.0 || qside<0.0 || qlong <0.0 || qinv<0.0 || deletas<0.0){
		sprintf(message,"bad sign: dely=%g, deleta=%g, deletas=%g, q=(%g,%g,%g,%g)\n",dely,deleta,deletas,qout,qside,qlong,qinv);
		barrayslog->Fatal(message);
	}
	
	ibin=floorl(qinv/Dqinv);	
	if(ibin>=0 && ibin<Nqbins){
		Bqinv[ibin]-=QaQb;
		Cqinv[ibin]+=CaCb;
	}
	
	ibin=floorl(qout/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqout[ibin]-=QaQb;
		Cqout[ibin]+=CaCb;
	}
	
	ibin=floorl(qside/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqside[ibin]-=QaQb;
		Cqside[ibin]+=CaCb;
	}
	
	ibin=floorl(qlong/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqlong[ibin]-=QaQb;
		Cqlong[ibin]+=CaCb;
	}
	
	ibin=floorl(deleta/Deta);
	if(ibin>=0 && ibin<Netabins){
		Beta[ibin]-=QaQb;
		Beta1[ibin]-=QaQb*cphi;
		Ceta[ibin]+=CaCb;
		Ceta1[ibin]+=CaCb*cphi;
	}
	
	ibin=floorl(dely/Dy);
	if(ibin>=0 && ibin<Nybins){
		By[ibin]-=QaQb;
		By1[ibin]-=QaQb*cphi;
		Cy[ibin]+=CaCb;
		Cy1[ibin]+=CaCb*cphi;
	}
	iy=ibin;
	
	ibin=floorl(deletas/Dy);
	if(ibin>=0 && ibin<Nybins){
		Betas[ibin]-=QaQb;
		Cetas[ibin]+=CaCb;
	}
	
	//only increment if dely is in measurable range
	if(fabs(dely)<acceptance->GetDelYMax(pida,pidb)){
		double phia=atan2(parta->p[2],parta->p[1]);
		if(sin(2.0*phia)<0.0){
			delphi=-delphi;
		}
		ibin=floorl((180.0+delphi)/Dphi);
		if(ibin>=0 && ibin<Nphibins){
			Bphi[ibin]-=QaQb;
			Cphi[ibin]+=CaCb;
		}
		iphi=ibin;
		if(iphi<Nphibins && iy<Nybins){
			Byphi[iy][iphi]-=QaQb;
			Cyphi[iy][iphi]+=CaCb;
		}
	}
}

CBFDenom::CBFDenom(CparameterMap *parmapset){
	parmap=parmapset;
	Nplus=Nminus=dNdy=0.0;
}

void CBFDenom::Reset(){
	Nplus=Nminus=dNdy=0.0;
}

void CBFDenom::Increment(CPart *part,double eff){
	int charge=part->resinfo->charge;
	if(charge==1)
		Nplus+=eff;
	else if(charge==-1)
		Nminus+=eff;
	else if (fabs(charge)>1){
		sprintf(message,"charge in CBFDenom::Increment > 1 !! = %d", charge);
		barrayslog->Info(message);
	}
}

void CBFNumer::WriteNumer(string dirname,string numertype,bool NoQ){
	string filename;
	FILE *fptr;
	int ibin,jbin;

	string command="mkdir -p "+dirname+"/"+name;
	system(command.c_str());
	
	if(!NoQ){
		filename=dirname+"/"+name+"/"+numertype+"_qinv.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ibin=0;ibin<Nqbins;ibin++){
			fprintf(fptr,"%7.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
			(0.5+ibin)*Dqinv,Bqinv[ibin],Bqout[ibin],Bqside[ibin],Bqlong[ibin],
			Cqinv[ibin],Cqout[ibin],Cqside[ibin],Cqlong[ibin]);
		}
		fclose(fptr);
	}
	
	filename=dirname+"/"+name+"/"+numertype+"_y.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nybins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Dy,By[ibin],Cy[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_y1.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nybins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Dy,By1[ibin],Cy1[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_eta.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Deta,Beta[ibin],Ceta[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_eta1.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Deta,Beta1[ibin],Ceta1[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_etas.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Deta,Betas[ibin],Betas[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_phi.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nphibins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",-180.0+(0.5+ibin)*Dphi,Bphi[ibin],Cphi[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_Byphi.txt";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Nybins= %d Nphibins= %d\n",Nybins,Nphibins);
	for(ibin=0;ibin<Nybins;ibin++){
		for(jbin=0;jbin<Nphibins;jbin++)
			fprintf(fptr,"%12.4e ",Byphi[ibin][jbin]);
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_Cyphi.txt";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Nybins= %d Nphibins= %d\n",Nybins,Nphibins);
	for(ibin=0;ibin<Nybins;ibin++){
		for(jbin=0;jbin<Nphibins;jbin++)
			fprintf(fptr,"%12.4e ",Cyphi[ibin][jbin]);
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
}

void CBFNumer::Print(){
	int ibin;
	sprintf(message,"----- By ------\n");
	for(ibin=0;ibin<Nybins;ibin++){
		sprintf(message,"%s%6.1f %10.3e\n",message,(0.5+ibin)*Dy,By[ibin]);
	}
	barrayslog->Info(message);
	sprintf(message,"----- Beta -----\n");
	for(ibin=0;ibin<Netabins;ibin++){
		sprintf(message,"%s%6.1f %9.5f %9.5f\n",message,(0.5+ibin)*Deta,Beta[ibin],Beta1[ibin]);
	}
	barrayslog->Info(message);
	sprintf(message,"----- Bphi -----\n");
	for(ibin=0;ibin<Nphibins;ibin++){
		sprintf(message,"%s%3d: %6.1f %10.3e\n",message,ibin,-180.0+(0.5+ibin)*Dphi,Bphi[ibin]);
	}
	barrayslog->Info(message);
	
}
