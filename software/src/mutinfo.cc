#include "mutinfo.h"
#include "part.h"
#include "b3d.h"
#include "cell.h"
#include "log.h"
#include "resonances.h"

CB3D *CMuTInfo::b3d=NULL;
int CMuTInfo::NETEVENTS=0;
int CMuTInfo::NMINCALC=2;
int CMuTInfo::NXY=0;
double CMuTInfo::DXY=0.0;
vector<vector<double>> CMuTInfo::taumin{};
vector<CResInfo *> CMuTInfo::Bresinfo{};
vector<double> CMuTInfo::massB{938.27,1189.37,1314.83,1115.68,1232.0,1385,1530,1672.43};
vector<double> CMuTInfo::degenB{8,12,8,4,16,12,8,4};

CMuTInfo::CMuTInfo(double tau_set){
	tau=tau_set;
	Epi=EK=0.0;
	Pxpi=Pypi=PxK=PyK=0.0;
	Txxpi=Tyypi=Txypi=0.0;
	TxxK=TyyK=TxyK=0.0;
	Npi=NK;
	Tpi=TK=145.0;

	NB.resize(8);
	TB.resize(8);
	muB.resize(8);
	TxxB.resize(8);
	TyyB.resize(8);
	TxyB.resize(8);
	PxB.resize(8);
	PyB.resize(8);
	EB.resize(8);
	epsilonB.resize(8);
	rhoB.resize(8);
	UxB.resize(8);
	UyB.resize(8);

	for(int btype=0;btype<8;btype++){
		NB[btype]=0;
		EB[btype]=PxB[btype]=PyB[btype]=0.0;
		TxxB[btype]=TyyB[btype]=TxyB[btype]=0.0;
		muB[btype]=0.0;
		TB[btype]=100.0;
		epsilonB[btype]=0.0;
	}
	sufficientN=false;
}

void CMuTInfo::Print(){
	char message[500];
	sprintf(message,"-------- MuT Info, tau=%g ----------\n",tau);

	sprintf(message,"%sNpi=%d, Epi/N=%g, Pxpi/Npi=%g, Pypi/Npi=%g\n",
		message,Npi,Epi/Npi,Pxpi/Npi,Pypi/Npi);
	sprintf(message,"%sTpi=%g, mupi=%g\n",message,Tpi,mupi);

	sprintf(message,"%sNK=%d, EK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
		message,NK,EK/NK,PxK/NK,PyK/NK);
	sprintf(message,"%sTK=%g, muK=%g\n",message,TK,muK);

	CLog::Info(message);

	for(int btype=0;btype<8;btype++){
		sprintf(message,"btype=%d\n",btype);
		sprintf(message,"%sNB=%d, EK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
			message,NB[btype],EB[btype]/NB[btype],PxB[btype]/NB[btype],PyB[btype]/NB[btype]);
		sprintf(message,"%sTK=%g, muK=%g\n",message,TB[btype],muB[btype]);
		CLog::Info(message);
	}

}

void CMuTInfo::CalcAllMuTU(){
	double Txx,Tyy,Txy,T00,T0x,T0y,gamma,degen;
	int btype;

	double volume=4.0*tau*2.0*b3d->ETAMAX*DXY*DXY*double(NETEVENTS);   // factor or 4 due to combining quadrants
	if(Npi>NMINCALC){
		Txx=Txxpi/volume;
		Tyy=Tyypi/volume;
		Txy=Txypi/volume;
		T00=Epi/volume;
		T0x=Pxpi/volume;
		T0y=Pypi/volume;
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Uxpi,Uypi,epsilonpi);
		gamma=sqrt(1.0+Uxpi*Uxpi+Uypi*Uypi);
		rhopi=double(Npi)/(gamma*volume);
		degen=3.0;
		GetMuT(PionMass,degen,rhopi,epsilonpi,Tpi,mupi);
	}
	else{
		Tpi=-1.0;
		mupi=0.0;
	}

	if(NK>NMINCALC){
		Txx=TxxK/volume;
		Tyy=TyyK/volume;
		Txy=TxyK/volume;
		T00=EK/volume;
		T0x=PxK/volume;
		T0y=PyK/volume;
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxK,UyK,epsilonK);
		gamma=sqrt(1.0+UxK*UxK+UyK*UyK);
		rhoK=double(NK)/(gamma*volume);
		degen=4.0;
		GetMuT(KaonMass,degen,rhoK,epsilonK,TK,muK);
	}
	else{
		TK=-1.0;
		muK=0.0;
	}
	for(btype=0;btype<8;btype++){
		if(NB[btype]>NMINCALC){
			Txx=TxxB[btype]/volume;
			Tyy=TyyB[btype]/volume;
			Txy=TxyB[btype]/volume;
			T00=EB[btype]/volume;
			T0x=PxB[btype]/volume;
			T0y=PyB[btype]/volume;
			GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxB[btype],UyB[btype],epsilonB[btype]);
			gamma=sqrt(1.0+UxB[btype]*UxB[btype]+UyB[btype]*UyB[btype]);
			rhoB[btype]=double(NB[btype])/(gamma*volume);
			GetMuT(massB[btype],degenB[btype],rhoB[btype],epsilonB[btype],TB[btype],muB[btype]);
		}
		else{
			TB[btype]=-1.0;
			muB[btype]=0.0;
		}
	}
}

void CMuTInfo::GetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
double &Ux,double &Uy,double &epsilon){
	double Qx,Qy,dQxdUx,dQxdUy,dQydUx,dQydUy;
	double gamma,Det,dUx,dUy;
	double A,B,C,dAdUx,dAdUy,dBdUx,dBdUy,dCdUx,dCdUy,dUmag;
	bool success=false;
	int ntry=0;

	Ux=T0x/T00;
	Uy=T0y/T00;
	do{
		gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
		A=-gamma*T00;
		B=(1.0+gamma/(1.0+gamma))*(T0x*Ux+T0y*Uy);
		C=-(1.0/(1.0+gamma))*(Ux*Txx*Ux+2.0*Ux*Txy*Uy+Uy*Tyy*Uy);
		Qx=gamma*T0x-Ux*Txx-Uy*Txy+(A+B+C)*Ux;
		Qy=gamma*T0y-Uy*Tyy-Ux*Txy+(A+B+C)*Uy;
		if(Qx*Qx+Qy*Qy>1.0E-10){
			dAdUx=A*Ux/(gamma*gamma);
			dAdUy=A*Uy/(gamma*gamma);
			dBdUx=(1.0+(gamma/(1.0+gamma)))*T0x+(Ux/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dBdUy=(1.0+(gamma/(1.0+gamma)))*T0y+(Uy/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dCdUx=-C*Ux/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Txx*Ux+Txy*Uy);
			dCdUy=-C*Uy/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Tyy*Uy+Txy*Ux);
			dQxdUx=-Txx+(Ux/gamma)*T0x+Ux*(dAdUx+dBdUx+dCdUx);
			dQxdUy=-Txy+(Uy/gamma)*T0x+Ux*(dAdUy+dBdUy+dCdUy);
			dQydUx=-Txy+(Ux/gamma)*T0y+Uy*(dAdUx+dBdUx+dCdUx);
			dQydUy=-Tyy+(Uy/gamma)*T0y+Uy*(dAdUy+dBdUy+dCdUy);
			dQxdUx+=(A+B+C);
			dQydUy+=(A+B+C);
			Det=dQxdUx*dQydUy-dQxdUy*dQydUx;
			dUx=-(dQydUy*Qx-dQxdUy*Qy)/Det;
			dUy=-(-dQydUx*Qx+dQxdUx*Qy)/Det;
			if(fabs(dUx*dUx+dUy*dUy)>0.5){
				dUmag=sqrt(dUx*dUx+dUy*dUy);
				dUx*=0.5/dUmag;
				dUy*=0.5/dUmag;
			}
			Ux+=dUx;
			Uy+=dUy;
			ntry+=1;
		}
		else
			success=true;
	}while(!success && ntry<200);
	if(!success){
		char message[300];
		sprintf(message,"Yikes, No Convergence in CMuTInfo::GetEpsilonU\n");
		sprintf(message,"%sQx=%g, Qy=%g, Ux=%g, Uy=%g\n",message,Qx,Qy,Ux,Uy);
		sprintf(message,"%sT00=%g, T0x=%g, T0y=%g\n",message,T00,T0x,T0y);
		sprintf(message,"%sTxx=%g, Tyy=%g, Txy=%g\n",message,Txx,Tyy,Txy);
		CLog::Fatal(message);
	}
	gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
	epsilon=gamma*T00*gamma-2.0*gamma*T0x*Ux-2.0*gamma*T0y*Uy
	+Ux*Txx*Ux+Uy*Tyy*Uy+2.0*Ux*Txy*Uy;
}

void CMuTInfo::GetMuT(double mass,double degen,double rho_target,double epsilon_target,double &T,double &mu){
	double E,dEdT,ETarget=epsilon_target/rho_target,epsilon0,dedT,sigma2,P,rho0,dT;
	int ntry=0;
	char message[100];
	do{
		ntry+=1;
		CResList::freegascalc_onespecies(mass,T,epsilon0,P,rho0,sigma2,dedT);
		E=epsilon0/rho0;
		dEdT=dedT/rho0-epsilon0*epsilon0/(rho0*rho0*T*T);
		dT=(ETarget-E)/dEdT;
		if(fabs(dT)>0.6*T)
			dT=0.6*T*dT/fabs(dT);
		T+=dT;
	}while(fabs(dT)>1.0E-5 && ntry<20);
	if(ntry==20){
		sprintf(message,"CMuTInfo::GetMuT did not converge!!!, T=%g, dT=%g\n",T,dT);
		CLog::Info(message);
	}
	CResList::freegascalc_onespecies(mass,T,epsilon0,P,rho0,sigma2,dedT);
	mu=log(rho_target/(rho0*degen));
}

/*
bool CMuTInfo::GetMuT_Baryon(double rhoB_target,double rhoBS_target,double epsilon_target,double &T,double &muB,double &muBS){
	char message[500];
	CResInfo *resinfo;
	Eigen::Vector3d rho,rho_target,x,dx,drho;
	Eigen::Matrix3d drhodx;
	vector<string> name{"nucleon","Lambda","Cascade"};
	double mass,degen,nstrange;
	rho_target(0)=epsilon_target/rhoB_target;
	rho_target(1)=rhoB_target;
	rho_target(2)=rhoBS_target;
	double factor,dfactordxB,dfactordxBS;
	double e,dedt,dedx1,dedx2;
	double P,sigma2,epsilon0,rho0,dedt0,accuracy,D;
	int ntry=0,ispecies,nmax=1000;
	bool success=false;


	T=100.0;
	x(0)=T;
	x(1)=x(2)=0.0;

	do{
		ntry+=1;
		T=x(0);
		rho.setZero();
		drhodx.setZero();
		e=dedt=dedx1=dedx2=0.0;
		for(ispecies=0;ispecies<int(Bresinfo.size());ispecies++){
			resinfo=Bresinfo[ispecies];
			degen=2.0*(resinfo->spin+1.0);
			mass=resinfo->mass;
			nstrange=fabs(resinfo->strange);
			CResList::freegascalc_onespecies(mass,T,epsilon0,P,rho0,sigma2,dedt0);

			factor=degen*exp(x(1)+nstrange*x(2));
			dfactordxB=factor;
			dfactordxBS=factor*nstrange;

			e+=epsilon0*factor;
			rho(1)+=rho0*factor;
			rho(2)+=rho0*factor*nstrange;

			dedt+=dedt0*factor;
			dedx1+=epsilon0*dfactordxB;
			dedx2+=epsilon0*dfactordxBS;

			drhodx(1,0)+=(epsilon0/(T*T))*factor;
			drhodx(1,1)+=rho0*dfactordxB;
			drhodx(1,2)+=rho0*dfactordxBS;

			drhodx(2,0)+=(epsilon0/(T*T))*factor*nstrange;
			drhodx(2,1)+=rho0*dfactordxB*nstrange;
			drhodx(2,2)+=rho0*dfactordxBS*nstrange;
		}
		rho(0)=e/rho(1);
		drhodx(0,0)=(dedt/rho(1))-e*drhodx(1,0)/(rho(1)*rho(1));
		drhodx(0,1)=(dedx1/rho(1))-e*drhodx(1,1)/(rho(1)*rho(1));
		drhodx(0,2)=(dedx2/rho(1))-e*drhodx(1,2)/(rho(1)*rho(1));

		drho=rho_target-rho;
		dx=drhodx.colPivHouseholderQr().solve(drho);

		if(dx(0)!=dx(0) || dx(1)!=dx(1) || dx(2)!=dx(2)){
			sprintf(message,"FAILURE in GetMuT_Baryon, ntry=%d, dx!=dx\n",ntry);
			sprintf(message,"%srho_target=(%g,%g,%g)\n",message,rho_target(0),rho_target(1),rho_target(2));
			sprintf(message,"%srho=(%g,%g,%g)\n",message,rho(0),rho(1),rho(2));
			sprintf(message,"%sdx=(%g,%g,%g)\n",message,dx(0),dx(1),dx(2));
			sprintf(message,"%sTB=%g, muB=%g, muBS=%g\n",message,x[0],x[1],x[2]);
			CLog::Fatal(message);
		}

		double dxmax=0.5*x(0);
		if(fabs(dx(0))>dxmax){
			D=fabs(dx(0));
			dx=dx*dxmax/D;
		}
		dxmax=1.0;
		if(fabs(dx(1))>dxmax){
			D=fabs(dx(1));
			dx=dx*dxmax/D;
		}
		dxmax=1.0;
		if(fabs(dx(2))>dxmax){
			D=fabs(dx(2));
			dx=dx*dxmax/D;
		}

		x=x+dx;

		accuracy=0.001*drho(0)*drho(0)+drho(1)*drho(1)+drho(2)*drho(2);

	}while(accuracy>1.0E-10 && ntry<nmax);
	if(ntry<nmax){
		success=true;
	}
	else{
		sprintf(message,"In CMuTInfo::GetMuT_Baryon, ntry=%d\n",ntry);
		sprintf(message,"%sT=%g, muB=%g, muBS=%g, E/rhoB=%g=?%g, rhoB=%g=?%g, rhoBS=%g=?%g\n",
			message,T,muB,muBS,rho[0],rho_target[0],rho[1],rho_target[1],rho[2],rho_target[2]);
		CLog::Info(message);
	}
	
	T=x(0);
	muB=x(1);
	muBS=x(2);
	return success;
}
*/

void CMuTInfo::GetIxIy(double x,double y,int &ix,int &iy){
	ix=lrint(floor(fabs(x)/DXY));
	iy=lrint(floor(fabs(y)/DXY));
}

int CMuTInfo::GetBtype(int pid){
	int btype=-1;
	pid=abs(pid);
	if(pid==2122 || pid==2112)
		btype=0;
	else if(pid==3222 || pid==3212 || pid==3112)
		btype=1;
	else if(pid==3322 || pid==3312)
		btype=2;
	else if(pid==3112)
		btype=3;
	else if(pid==2224 || pid==2214 || pid==2114 || pid==1114)
		btype=4;
	else if(pid==3224 || pid==3214 || pid==3114)
		btype=5;
	else if(pid==3324 || pid==3314)
		btype=6;
	else if(pid==3334)
		btype=7;
	return btype;
}