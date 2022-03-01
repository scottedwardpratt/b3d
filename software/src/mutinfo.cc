#include "mutinfo.h"
#include "part.h"
#include "b3d.h"
#include "cell.h"
#include "log.h"
#include "resonances.h"

CB3D *CMuTInfo::b3d=NULL;
int CMuTInfo::NETEVENTS=0;
int CMuTInfo::NTAU=0;
int CMuTInfo::NMINCALC=10;
double CMuTInfo::DELTAU=0.0;
int CMuTInfo::NXY=0;
int CMuTInfo::DXY=0.0;
vector<vector<double>> CMuTInfo::taumin{};
vector<CResInfo *> CMuTInfo::Bresinfo{};

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
	sufficientN=false;
}

 void CMuTInfo::Init(CB3D *b3dset){
 	b3d=b3dset;
 	DELTAU=b3d->MUTCALC_DELTAU;
 	NT
		CMuTInfo::DELTAU=MUTCALC_DELTAU;
		CMuTInfo::NTAU=TAUCOLLMAX/MUTCALC_DELTAU;
		CMuTInfo::NETEVENTS=0;
		muTinfo.resize(CMuTInfo::NTAU);
		for(itau=0;itau<CMuTInfo::NTAU;itau++){
			muTinfo[itau].resize(2*NXY);
			for(ix=0;ix<2*NXY;ix++)
				muTinfo[itau][ix].resize(2*NXY);
			for(ix=0;ix<2*NXY;ix++){
				for(iy=0;iy<2*NXY;iy++){
					muTinfo[itau][ix][iy]=new CMuTInfo((itau+0.5)*CMuTInfo::DELTAU);
				}
			}
		}
		CResInfoMap::iterator rpos;
		CResInfo *resinfo;
		for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();++rpos){
			resinfo=rpos->second;
			if(resinfo->baryon>0){
				CMuTInfo::Bresinfo.push_back(resinfo);
			}
		}
 }

void CMuTInfo::UpdateNPE(CB3DCell *cell){
	CResInfo *resinfo;
	int S;
	double gamma,gammav,E;
	CPartMap::iterator ppos;
	CPart *part;
	for(ppos=cell->partmap.begin();ppos!=cell->partmap.end();ppos++){
		part=ppos->second;
		resinfo=part->resinfo;
		if(resinfo->code==111 || abs(resinfo->code)==211){
			gamma=cosh(part->eta);
			gammav=sinh(part->eta);
			Npi+=1;
			Pxpi+=part->p[1];
			Pypi+=part->p[2];
			E=gamma*part->p[0]-gammav*part->p[3];
			Epi+=E;
			Txxpi+=part->p[1]*part->p[1]/E;
			Tyypi+=part->p[2]*part->p[2]/E;
			Txypi+=part->p[1]*part->p[2]/E;
		}
		else if(abs(part->resinfo->code)==321 || abs(part->resinfo->code)==311){
			gamma=cosh(part->eta);
			gammav=sinh(part->eta);
			NK+=1;
			PxK+=part->p[1];
			PyK+=part->p[2];
			E=gamma*part->p[0]-gammav*part->p[3];
			EK+=E;
			TxxK+=part->p[1]*part->p[1]/E;
			TyyK+=part->p[2]*part->p[2]/E;
			TxyK+=part->p[1]*part->p[2]/E;
		}
		else if(resinfo->baryon!=0){
			gamma=cosh(part->eta);
			gammav=sinh(part->eta);
			NB+=1;
			PxB+=part->p[1];
			PyB+=part->p[2];
			E=gamma*part->p[0]-gammav*part->p[3];
			EB+=E;
			S=abs(resinfo->strange);
			NBS+=S;
			TxxB+=part->p[1]*part->p[1]/E;
			TyyB+=part->p[2]*part->p[2]/E;
			TxyB+=part->p[1]*part->p[2]/E;
		}
	}
}

void CMuTInfo::Print(){
	char message[500];
	sprintf(message,"-------- MuT Info, tau=%g ----------\n",tau);
	sprintf(message,"%sNpi=%d, Epi/N=%g, Pxpi/Npi=%g, Pypi/Npi=%g\n",
	message,Npi,Epi/Npi,Pxpi/Npi,Pypi/Npi);
	sprintf(message,"%sNK=%d, EK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
	message,NK,EK/NK,PxK/NK,PyK/NK);
	sprintf(message,"%sNB=%d, NBS=%d, EB/NB=%g, PxB/NB=%g, PyB/NB=%g\n",
	message,NB,NBS,EB/NB,PxB/NB,PyB/NB);
	sprintf(message,"%sTpi=%g, TK=%g, TB=%g\nmupi=%g, muK=%g, muB=%g, muBS=%g\n",
		message,Tpi,TK,TB,mupi,muK,muB,muBS);
	CLog::Info(message);
}

void CMuTInfo::CalcAllMuTU(){
	char message[100];
	bool success;
	double Txx,Tyy,Txy,T00,T0x,T0y,gamma,degen;
	if(Npi>NMINCALC){
		Txx=Txxpi/(volume*double(NETEVENTS));
		Tyy=Tyypi/(volume*double(NETEVENTS));
		Txy=Txypi/(volume*double(NETEVENTS));
		T00=Epi/(volume*double(NETEVENTS));
		T0x=Pxpi/(volume*double(NETEVENTS));
		T0y=Pypi/(volume*double(NETEVENTS));
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Uxpi,Uypi,epsilonpi);
		gamma=sqrt(1.0+Uxpi*Uxpi+Uypi*Uypi);
		rhopi=double(Npi)/(gamma*volume*double(NETEVENTS));
		degen=3.0;
		GetMuT(PionMass,degen,rhopi,epsilonpi,Tpi,mupi);
	}
	else{
		Tpi=-1.0;
		mupi=0.0;
	}
	if(NK>NMINCALC){
		Txx=TxxK/(volume*double(NETEVENTS));
		Tyy=TyyK/(volume*double(NETEVENTS));
		Txy=TxyK/(volume*double(NETEVENTS));
		T00=EK/(volume*double(NETEVENTS));
		T0x=PxK/(volume*double(NETEVENTS));
		T0y=PyK/(volume*double(NETEVENTS));
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxK,UyK,epsilonK);
		gamma=sqrt(1.0+UxK*UxK+UyK*UyK);
		rhoK=double(NK)/(gamma*volume*double(NETEVENTS));
		degen=4.0;
		GetMuT(KaonMass,degen,rhoK,epsilonK,TK,muK);
	}
	else{
		TxxK=-1.0;
		muK=0.0;
	}
	if(NB>NMINCALC && NBS>0){ 
		Txx=TxxB/(volume*double(NETEVENTS));
		Tyy=TyyB/(volume*double(NETEVENTS));
		Txy=TxyB/(volume*double(NETEVENTS));
		T00=EB/(volume*double(NETEVENTS));
		T0x=PxB/(volume*double(NETEVENTS));
		T0y=PyB/(volume*double(NETEVENTS));
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxB,UyB,epsilonB);
		gamma=sqrt(1.0+UxB*UxB+UyB*UyB);
		rhoB=double(NB)/(gamma*volume*double(NETEVENTS));
		rhoBS=double(NBS)/(gamma*volume*double(NETEVENTS));
		success=GetMuT_Baryon(rhoB,rhoBS,epsilonB,TB,muB,muBS);
		if(!success){
			sprintf(message,"GetMuT_Baryon Failed\n");
			CLog::Info(message);
			Print();
		}

	}
	else{
		TB=-1.0;
		muB=0.0;
		muBS=0.0;
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
void CMuTInfo::GetMuT(double mass,double degen,double rho_target,double epsilon_target,double &T,double &mu){
	double x,rho,epsilon,P,sigma2,depsilon,drho,epsilon0,rho0;
	double dT,dmu,dedmu,dedT,drhodmu,drhodT,accuracy;
	bool success=false;
	int ntry=0;
	double Det;
	T=0.5*((epsilon_target/rho_target)-135.0);
	mu=0.0;
	x=exp(mu);
	do{
		ntry+=1;
		CResList::freegascalc_onespecies(mass,T,epsilon0,P,rho0,sigma2,dedT);
		epsilon=epsilon0*degen*x;
		rho=rho0*degen*x;
		accuracy=pow((rho-rho_target)/rho_target,2)+pow((epsilon-epsilon_target)/epsilon_target,2);
		if(accuracy>1.0E-8){
			dedmu=epsilon;
			dedT*=degen*x;
			drhodT=epsilon/(T*T);
			drhodmu=rho;
			Det=dedT*drhodmu-dedmu*drhodT;
			depsilon=epsilon_target-epsilon;
			drho=rho_target-rho;

			dT=(drhodmu*depsilon-dedmu*drho)/Det;
			dmu=-(drhodT*depsilon-dedT*drho)/Det;

			if(fabs(dT)>0.5*T){
				dmu=(0.5*T/fabs(dT))*dmu;
				dT=(0.5*T/fabs(dT))*dT;
			}
			if(fabs(dmu)>1.0){
				dT=(1.0/fabs(dmu))*dT;
				dmu=(1.0/fabs(dmu))*dmu;
			}

			T+=dT;
			mu+=dmu;
			
		}
		else
			success=true;
	}while(!success && ntry<40);
	if(ntry==40){
		printf("failure, accuracy=%g\n",accuracy);
	}
}
*/

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
