#include "mucalc.h"

CB3D *CLocalSpeciesInfo::b3d=NULL;
CB3D *CLocalInfo::b3d=NULL;
int CLocalInfo::NRBINS=0;
int CLocalInfo::NETEVENTS=0;
int CLocalInfo::IMUCALC=0;
double CLocalInfo::DELR2=0.0;
bool CLocalInfo::printing=false;

void CLocalSpeciesInfo::Initialize(){
	degen=0;
	int nrbins=CLocalInfo::NRBINS;
	T.resize(nrbins);
	mu.resize(nrbins);
	ur.resize(nrbins);
	N.resize(nrbins);
	E.resize(nrbins);
	Pr.resize(nrbins);
	M.resize(nrbins);
	Zero();
}

void CLocalSpeciesInfo::Zero(){
	int ir2;
	for(ir2=0;ir2<CLocalInfo::NRBINS;ir2++){
		T[ir2]=mu[ir2]=ur[ir2]=-1.0;
		Pr[ir2]=E[ir2]=N[ir2]=0;
	}
}

void CLocalInfo::Initialize(){
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(211,&pion));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-211,&pion));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(111,&pion));
	pion.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(321,&kaon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-321,&kaon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(311,&kaon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-311,&kaon));
	kaon.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(2212,&nucleon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-2212,&nucleon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(2112,&nucleon));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-2112,&nucleon));
	nucleon.Initialize();
	

	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(1114,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-1114,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(2114,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-2114,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(2214,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-2214,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(2224,&delta));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-2224,&delta));
	delta.Initialize();
		
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3122,&lambda));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3122,&lambda));
	lambda.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3222,&sigma));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3222,&sigma));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3212,&sigma));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3212,&sigma));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3112,&sigma));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3112,&sigma));
	sigma.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3224,&sigmastar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3224,&sigmastar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3214,&sigmastar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3214,&sigmastar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3114,&sigmastar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3114,&sigmastar));
	sigmastar.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3312,&xsi));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3312,&xsi));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3322,&xsi));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3322,&xsi));
	xsi.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3314,&xsistar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3314,&xsistar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3324,&xsistar));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3324,&xsistar));
	xsistar.Initialize();
	
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(3334,&omega));
	mucalc_species.insert(pair<int,CLocalSpeciesInfo *>(-3334,&omega));
	omega.Initialize();
	
	

}

void CLocalSpeciesInfo::MuCalc(){
	double delr2=CLocalInfo::DELR2;
	int nrbins=CLocalInfo::NRBINS;
	int ir2,netevents=CLocalInfo::NETEVENTS,ntry;
	double Tguess,Etot,u0,e,rho0,rhotarget,etarget,dedT,drhodT,accuracy;
	double dT,de,ddedT,dens,p,sigma2,volume;
	double tau=b3d->tau;
	double mass=resinfo->mass;
	
	for(ir2=0;ir2<nrbins;ir2++){
		if(T[ir2]>0)
			Tguess=T[ir2];
		else
			Tguess=150.0;
		if(N[ir2]>10){
			ur[ir2]=Pr[ir2]/E[ir2];
			ur[ir2]=ur[ir2]/sqrt(1.0-ur[ir2]*ur[ir2]);
		}
		if(N[ir2]>1){
			Etot=sqrt(E[ir2]*E[ir2]-Pr[ir2]*Pr[ir2]);
			u0=sqrt(1.0+ur[ir2]*ur[ir2]);
			volume=u0*PI*delr2*2.0*b3d->ETAMAX*tau;
			
			Etot=Etot/(volume*netevents);
			rhotarget=N[ir2]/(volume*netevents);
			etarget=Etot/rhotarget;
			if(etarget<mass){
				printf("etarget=%g, mass=%g, N[%d]=%d\n",etarget,resinfo->mass,ir2,N[ir2]);
				exit(1);
			}
			if(((etarget-mass)/mass)<0.05){
				T[ir2]=2.0*(etarget-mass)/3.0;
			}
			else{
				ntry=0;
				do{
					ntry+=1;
					e=dedT=drhodT=rho0=0.0;
					accuracy=0.0;
					b3d->reslist->freegascalc_onespecies(mass,Tguess,de,p,dens,sigma2,ddedT);
					rho0+=degen*dens;
					e+=degen*de;
					dedT+=degen*ddedT;
					drhodT+=degen*de/(Tguess*Tguess);
					e=e/rho0;
					accuracy=fabs(e-etarget);
					dT=(etarget-e)/(dedT/rho0-e*drhodT/rho0);
					if(fabs(dT)>0.5*Tguess)
						dT=0.5*Tguess*dT/fabs(dT);
					Tguess+=dT;
				}while(accuracy>0.1);
				T[ir2]=Tguess;
			}
			//printf("T[%d]=%g\n",ir2,T[ir2]);
			mu[ir2]=T[ir2]*log(rhotarget/rho0);
			/*
			if(name=="pion")
				printf("tau=%g, netevents=%d, N[%d]=%d, rhotarget=%g, mu=%g\n",
			b3d->tau,netevents,ir2,N[ir2],rhotarget,mu[ir2]);
			*/
		}
	}
}

void CLocalSpeciesInfo::MuCalc_PionsWithBose(){
	int nbose=b3d->NBOSE;
	double dalpha,alpha=0.0;
	double delr2=CLocalInfo::DELR2;
	int nrbins=CLocalInfo::NRBINS;
	int ir2,netevents=CLocalInfo::NETEVENTS,ntry;
	double Tguess,Etot,u0,rhotarget,etarget,dedT,ddedT,drhodT,accuracy,dT;
	double dedalpha,drhodalpha,denom,dele,delrho,e,rho,x,de,dens,sigma2,p,volume;
	double tau=b3d->tau;
	double mass=resinfo->mass;
	int n;

	for(ir2=0;ir2<nrbins;ir2++){
		if(N[ir2]>0){
			if(T[ir2]>0)
				Tguess=T[ir2];
			else
				Tguess=150.0;
			ur[ir2]=Pr[ir2]/E[ir2];
			ur[ir2]=ur[ir2]/sqrt(1.0-ur[ir2]*ur[ir2]);
			if(N[ir2]>1){
				Etot=sqrt(E[ir2]*E[ir2]-Pr[ir2]*Pr[ir2]);
				u0=sqrt(1.0+ur[ir2]*ur[ir2]);
				volume=u0*PI*delr2*2.0*b3d->ETAMAX*tau;
				etarget=Etot/(netevents*volume);
				rhotarget=N[ir2]/(netevents*volume);
				etarget=etarget/rhotarget;
				ntry=0;
				do{
					ntry+=1;
					accuracy=0.0;
					e=rho=dedT=drhodT=dedalpha=drhodalpha=0.0;
					for(n=1;n<=nbose;n++){
						b3d->reslist->freegascalc_onespecies(mass,Tguess/double(n),de,p,dens,sigma2,ddedT);
						x=exp(n*alpha);
						e+=degen*de*x;
						rho+=degen*dens*x;
						dedT+=degen*x*ddedT/double(n);
						drhodT+=degen*x*n*de/(Tguess*Tguess);
						dedalpha+=degen*de*x*n;
						drhodalpha+=degen*dens*x*n;
					}
					dedT=(dedT/rho)-e*drhodT/(rho*rho);
					dedalpha=(dedalpha/rho)-e*drhodalpha/(rho*rho);
					e=e/rho;
					denom=dedT*drhodalpha-dedalpha*drhodT;
					dele=etarget-e;
					delrho=rhotarget-rho;
					if(ntry>30){
						printf("MuCalc:: ntry=%d, ir1=%d, etarget=%g, dele=%g, rhotarget=%g, delrho=%g, T=%g, mu=%g\n",
						ntry, ir2,etarget,dele,rhotarget,delrho,Tguess,alpha*Tguess);
						//Misc::Pause();
					}
					accuracy=fabs(dele/etarget)+fabs(delrho/rhotarget);
					dT=(drhodalpha*dele-dedalpha*delrho)/denom;
					if(fabs(dT)>0.5*Tguess){
						dT=0.5*dT*Tguess/fabs(dT);
					}
					Tguess+=dT;
					dalpha=(dedT*delrho-drhodT*dele)/denom;
					if(fabs(dalpha)>0.5){
						dalpha=dalpha*0.5/fabs(dalpha);
					}
					alpha+=dalpha;
				}while(accuracy>0.01*rhotarget);
				T[ir2]=Tguess;
				mu[ir2]=alpha*T[ir2];
			}
		}
	}
}

void CLocalInfo::MuCalc(){
	double r2,r,pr,e,eperp;
	double tau=b3d->tau;
	double r2max=CLocalInfo::DELR2*CLocalInfo::NRBINS;
	CLocalSpeciesInfo *lsi;
	int ir2;
	CResInfo *resinfo;
	CPart *part;
	CPartMap::iterator ppos;
		
	CPartMap *pmap=&(b3d->PartMap);
	ppos=pmap->begin();
	while(ppos!=pmap->end()){
		part=ppos->second;
		if(part->active){
			part->Propagate(tau);
			r2=part->r[1]*part->r[1]+part->r[2]*part->r[2];
			ir2=lrint(floor(NRBINS*r2/r2max));
			if(ir2<NRBINS){
				r=sqrt(r2);
				pr=(part->r[1]*part->p[1]+part->r[2]*part->p[2])/r;
				eperp=sqrt(part->msquared+part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				e=eperp*cosh(part->y-part->eta);
				
				resinfo=part->resinfo;
				lsi=NULL;
				
				if(abs(resinfo->code)==211 || resinfo->code==111){
					lsi=&pion;
				}
				else if(abs(resinfo->code)==321 || abs(resinfo->code)==311){
					lsi=&kaon;
				}
				if(resinfo->baryon!=0){
					if(abs(resinfo->code)==2112 || abs(resinfo->code)==2212){
						lsi=&nucleon;
					}
					else if(abs(resinfo->code)==1114 || abs(resinfo->code)==2114 || (resinfo->code)==2214 || (resinfo->code)==2224){
						lsi=&delta;
					}
					else if(abs(resinfo->code)==3122){
						lsi=&lambda;
					}
					else if(abs(resinfo->code)==3222 || abs(resinfo->code)==3212 || (resinfo->code)==3112){
						lsi=&sigma;
					}
					else if(abs(resinfo->code)==3224 || abs(resinfo->code)==3214 || (resinfo->code)==3114){
					lsi=&sigmastar;
					}
					else if(abs(resinfo->code)==3312 || abs(resinfo->code)==3322){
						lsi=&xsi;
					}
					else if(abs(resinfo->code)==3314 || abs(resinfo->code)==3324){
					lsi=&xsistar;
					}
					else if(abs(resinfo->code)==3334){
						lsi=&omega;
					}
				}
				if(lsi!=NULL){
					lsi->N[ir2]+=1;
					lsi->E[ir2]+=e;
					lsi->Pr[ir2]+=pr;
					lsi->M[ir2]+=part->GetMass();
				}
			}
		}
		ppos++;
	}
	
	pion.MuCalc();
	kaon.MuCalc();
	nucleon.MuCalc();
	delta.MuCalc();
	lambda.MuCalc();
	sigma.MuCalc();
	sigmastar.MuCalc();
	xsi.MuCalc();
	xsistar.MuCalc();
	omega.MuCalc();

	if(printing){
		pion.MuPrint();
		kaon.MuPrint();
		nucleon.MuPrint();
		delta.MuPrint();
		lambda.MuPrint();
		sigma.MuPrint();
		sigmastar.MuPrint();
		xsi.MuPrint();
		xsistar.MuPrint();
		omega.MuPrint();
	}	
}

void CLocalSpeciesInfo::MuPrint(){
	double rho;
	FILE *fptr;
	char command[120];
	char filename[120];
	sprintf(command,"mkdir -p mucalc_results/%s",b3d->run_name.c_str());
	system(command);
	sprintf(filename,"mucalc_results/%s/%s_tau%g.txt",b3d->run_name.c_str(),resinfo->name.c_str(),b3d->tau);
	fptr=fopen(filename,"w");
	fprintf(fptr,"#------- tau=%6.3f --------\n",b3d->tau);
	fprintf(fptr,"#r(fm)     T      mu    u_r      rho       N\n");
	for(int ir2=0;ir2<CLocalInfo::NRBINS;ir2++){
		rho=N[ir2]/(PI*CLocalInfo::DELR2*2.0*b3d->ETAMAX*b3d->tau*CLocalInfo::NETEVENTS);
		fprintf(fptr,"%6.3f %7.2f %7.2f %7.2f %9.7f %7d\n",
		sqrt(CLocalInfo::DELR2*(0.5+ir2)),T[ir2],mu[ir2],ur[ir2],rho,N[ir2]);
	}
	fclose(fptr);
}
