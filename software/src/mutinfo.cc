#include "mutinfo.h"
#include "part.h"
#include "b3d.h"
#include "cell.h"
#include "resonances.h"

CB3D *CMuTInfo::b3d=NULL;
int CMuTInfo::NETEVENTS=0;
int CMuTInfo::NTAU=0;
double CMuTInfo::DELTAU=0.0;

CMuTInfo::CMuTInfo(double tau_set){
	tau=tau_set;
	volume=tau*2.0*b3d->ETAMAX*b3d->DXY*b3d->DXY;
	Pxpi=Pypi=PxK=PyK=PxB=PyB=0.0;
	Txxpi=Tyypi=Txypi=0.0;
	TxxK=TyyK=TxyK=0.0;
	TxxB=TyyB=TxyB=0.0;
	Npi=NK=NBS=0;
	Tpi=TK=TB=145.0;
	mupi=muK=muB=muBS=0.0;
}

void CMuTInfo::UpdateNPE(CB3DCell *cell){
	CResInfo *resinfo;
	int S;
	CPartMap::iterator ppos;
	CPart *part;
	for(ppos=cell->partmap.begin();ppos!=cell->partmap.end();ppos++){
		part=ppos->second;
		resinfo=part->resinfo;
		if(resinfo->code==111 || abs(resinfo->code)==211){
			Npi+=1;
			Pxpi+=part->p[1];
			Pypi+=part->p[2];
			Epi+=cosh(part->eta)*part->p[0]-sinh(part->eta)*part->p[3];
			Txxpi+=part->p[1]*part->p[1]/part->p[0];
			Tyypi+=part->p[2]*part->p[2]/part->p[0];
			Txypi+=part->p[1]*part->p[2]/part->p[0];
		}
		else if(abs(part->resinfo->code)==321 || abs(part->resinfo->code)==311){
			NK+=1;
			PxK+=part->p[1];
			PyK+=part->p[2];
			EK+=cosh(part->eta)*part->p[0]-sinh(part->eta)*part->p[3];
			TxxK+=part->p[1]*part->p[1]/part->p[0];
			TyyK+=part->p[2]*part->p[2]/part->p[0];
			TxyK+=part->p[1]*part->p[2]/part->p[0];
		}
		else if(resinfo->baryon!=0){
			NB+=1;
			PxB+=part->p[1];
			PyB+=part->p[2];
			EB+=cosh(part->eta)*part->p[0]-sinh(part->eta)*part->p[3];
			S=abs(resinfo->strange);
			NBS+=S;
			TxxB+=part->p[1]*part->p[1]/part->p[0];
			TyyB+=part->p[2]*part->p[2]/part->p[0];
			TxyB+=part->p[1]*part->p[2]/part->p[0];
		}
	}
}

void CMuTInfo::Print(){
	char message[500];
	sprintf(message,"-------- MuT Info, tau=%g ----------\n",tau);
	sprintf(message,"%sNpi=%d, Epi/N=%g, Pxpi/Npi=%g, Pypi/Npi=%g\n",
	message,Npi,Epi/Npi,Pxpi/Npi,Pypi/Npi);
	sprintf(message,"%sNK=%d, EK/N=%g, PxK/NK=%g, PyK/NK=%g\n",
	message,NK,EK/NK,PxK/NK,PyK/NK);
	CLog::Info(message);
}

void CMuTInfo::CalcMuTU(){
	double Txx,Tyy,Txy,T00,T0x,T0y,gamma,degen;
	Txx=Txxpi/(volume*double(NETEVENTS));
	Tyy=Tyypi/(volume*double(NETEVENTS));
	Txx=Txxpi/(volume*double(NETEVENTS));
	T00=Epi/(volume*double(NETEVENTS));
	T0x=Pxpi/(volume*double(NETEVENTS));
	T0y=Pypi/(volume*double(NETEVENTS));
	GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Uxpi,Uypi,epsilonpi);
	gamma=sqrt(1.0+Uxpi*Uxpi+Uypi*Uypi);
	rhopi=double(Npi)/(gamma*volume*double(NETEVENTS));
	degen=3.0;
	GetMuT(PionMass,degen,rhopi,epsilonpi,Tpi,mupi);
}

void CMuTInfo::GetMuT(double mass,double degen,double rho_target,double epsilon_target,double &T,double &mu){
	double &x,rho,epsilon,P,sigma2,x;
	double dT,dx,dedx,dedT,drhodx,drhodT;
	do{
		EqofState::fregascalc_onespecies(mass,T,P,epsilon,rho,sigma2,dedT);

		epsilon*=degen;
		dedx=epsilon;
		epsilon*=x;

		dedT*=degen;

		rho*=degen;
		drhodx=rho;
		rho*=x;

		drhodT=degen*epsilon/(T*T);
		drhodT*=x;
	}
}

/*void CMuTInfo::CalcMuTU(){
	double E,volume,epsilon,U0,Ux,Uy,Pmag;
	double Mxx,Mxy,Myy;
	volume=b3d->DXY*b3d->DXY*tau*b3d->DETA;



	printf("Inside CalcMuTU\n");
}*/

void CMuTInfo::GetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
double &Ux,double &Uy,double &epsilon){
	double Qx,Qy,dQxdUx,dQxdUy,dQydUx,dQydUy;
	double gamma,Det,dUx,dUy;
	double A,B,C,dAdUx,dAdUy,dBdUx,dBdUy,dCdUx,dCdUy;
	bool success=false;

	Ux=T0x/T00;
	Uy=T0y/T00;
	do{
		gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
		A=gamma*T00;
		B=(1+gamma/(1.0+gamma))*(T0x*Ux+T0y*Uy);
		C=-(1.0/(1.0+gamma))*(Ux*Txx*Ux+2.0*Ux*Txy*Uy+Uy*Tyy*Uy);
		Qx=gamma*T0x+(A+B+C)*Ux;
		Qy=gamma*T0y+(A+B+C)*Uy;
		if(Qx*Qx+Qy*Qy>1.0E-8){
			dAdUx=A*Ux/(gamma*gamma);
			dAdUy=A*Uy/(gamma*gamma);
			dBdUx=(1.0+(gamma/(1.0+gamma)))*T0x+(Ux/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dBdUy=(1.0+(gamma/(1.0+gamma)))*T0y+(Uy/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dCdUx=-C*Ux/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Txx*Ux+Txy*Uy);
			dCdUy=-C*Uy/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Tyy*Uy+Txy*Ux);
			dQxdUx=(Ux/gamma)*T0x+Ux*(dAdUx+dBdUx+dCdUx);
			dQxdUy=(Uy/gamma)*T0x+Ux*(dAdUy+dBdUy+dCdUy);
			dQydUx=(Ux/gamma)*T0y+Uy*(dAdUx+dBdUx+dCdUx);
			dQydUy=(Uy/gamma)*T0y+Uy*(dAdUy+dBdUy+dCdUy);
			dQxdUx+=(A+B+C);
			dQydUy+=(A+B+C);
			Det=dQxdUx*dQydUy-dQxdUy*dQydUx;
			dUx=-(dQydUy*Qx-dQxdUy*Qy)/Det;
			dUy=-(-dQydUx*Qx+dQxdUx*Qy)/Det;
			Ux+=dUx;
			Uy+=dUy;
		}
		else
			success=true;
	}while(!success);
	gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
	epsilon=gamma*T00*gamma-2.0*gamma*T0x*Ux-2.0*gamma*T0y*Uy
	+Ux*Txx*Ux+Uy*Tyy*Uy+2.0*Ux*Txy*Uy;
}

void CMuTinfo::

/*
	
void CMuTInfo::FindMuTUxUy(double tau,int N,double E,double M,double Px,double Py,double degen,double &T,double &mu,double &ux,double &uy){
	double u0,e,dedT,Tguess,p,dens,sigma2,ddedT,rho0,Etot,rhotarget,etarget;
	double mass,volume,drhodT,accuracy,de,dT;
	int ntry;
	Tguess=T;
	if(Tguess!=Tguess || Tguess<20.0 || Tguess>180.0)
		Tguess=150.0;
	if(N==0){
		mu=T=ux=uy=0.0;
	}
	else if(N==1){
		ux=Px/E;
		uy=Py/E;
		u0=1.0/sqrt(1.0-ux*ux-uy*uy);
		ux*=u0;
		uy*=u0;
		mu=T=0.0;
	}
	else{
		ux=Px/E;
		uy=Py/E;
		u0=1.0/sqrt(1.0-ux*ux-uy*uy);
		ux*=u0;
		uy*=u0;
		Etot=sqrt(E*E-Px*Px-Py*Py);
		volume=u0*b3d->DXY*b3d->DXY*2.0*b3d->ETAMAX*tau;
		mass=M/N;
		etarget=Etot/N;
		rhotarget=N/(volume*NETEVENTS);
			
		if(((etarget-mass)/mass)<0.1){
			T=2.0*(etarget-mass)/3.0;
			if(T<5.0)
				T=5.0;
			b3d->reslist->freegascalc_onespecies(mass,T,de,p,dens,sigma2,ddedT);
			rho0=degen*dens;
		}
		else{
			ntry=0;
			do{
				ntry+=1;
				e=dedT=drhodT=rho0=0.0;
				accuracy=0.0;
				b3d->reslist->freegascalc_onespecies(mass,Tguess,de,p,dens,sigma2,ddedT);
				rho0=degen*dens;
				e=degen*de;
				dedT=degen*ddedT;
				drhodT=degen*de/(Tguess*Tguess);
				e=e/rho0;
				accuracy=fabs((e-etarget)/mass);
				dT=(etarget-e)/(dedT/rho0-e*drhodT/rho0);
				if(fabs(dT)>0.5*Tguess)
					dT=0.5*Tguess*dT/fabs(dT);
				Tguess+=dT;
			}while(accuracy>1.0E-6);
			T=Tguess;
		}
		mu=T*log(rhotarget/rho0);
	}
}

void CMuTInfo::FindMuTInfo_pi(int itau){
	double tau=(1+itau)*DELTAU;
	FindMuTUxUy(tau,Npi[itau],Epi[itau],Mpi[itau],Pxpi[itau],Pypi[itau],3,Tpi,mupi,uxpi,uypi);
}
void CMuTInfo::FindMuTInfo_K(int itau){
	double tau=(1+itau)*DELTAU;
	FindMuTUxUy(tau,NK[itau],EK[itau],MK[itau],PxK[itau],PyK[itau],4,TK,muK,uxK,uyK);
}


}
*/
