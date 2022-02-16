#include "b3d.h"
#include "randy.h"
#include "misc.h"
#include "resonances.h"

using namespace std;

void TestGetMuT(double mass,double degen,double rho_target,double epsilon_target,double &T,double &mu){
	double x,rho,epsilon,P,sigma2,depsilon,drho,epsilon0,rho0;
	double dT,dx,dedx,dedT,drhodx,drhodT,accuracy;
	bool success=false;
	int ntry=0;
	double Det;
	printf("epsilon_target=%g, rho_target=%g\n",epsilon_target,rho_target);
	x=exp(mu);
	do{
		ntry+=1;
		CResList::freegascalc_onespecies(mass,T,epsilon0,P,rho0,sigma2,dedT);
		epsilon=epsilon0*degen*x;
		rho=rho0*degen*x;
		accuracy=pow((rho-rho_target)/rho_target,2)+pow((epsilon-epsilon_target)/epsilon_target,2);
		printf("T=%g, epsilon=%g, rho=%g, epsilon/rho=%g, accuracy=%g\n",T,epsilon,rho,epsilon/rho,accuracy);
		if(accuracy>1.0E-8){
			dedx=epsilon0*degen;
			dedT*=degen*x;
			drhodT=epsilon/(T*T);
			drhodx=rho0*degen;
			Det=dedT*drhodx-dedx*drhodT;
			depsilon=epsilon_target-epsilon;
			drho=rho_target-rho;

			dT=(drhodx*depsilon-dedx*drho)/Det;
			dx=-(drhodT*depsilon-dedT*drho)/Det;

			printf("dT=%g, dx=%g\n",dT,dx);

			if(fabs(dT)>0.5*T){
				dx=(0.5*T/dT)*dx;
				dT=0.5*T;
			}

			T+=dT;
			x+=dx;
		}
		else
			success=true;
	}while(!success && ntry<20);

	mu=log(x);
}

void TestGetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
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
		//printf("--------------------\n");
		//printf("Qx=%g, Qy=%g\n",Qx,Qy);
		//printf("Ux=%g, Uy=%g\n",Ux,Uy);
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
			if(fabs(dUx*dUx+dUy*dUy)>1.0){
				dUmag=sqrt(dUx*dUx+dUy*dUy);
				dUx*=1.0/dUmag;
				dUy*=1.0/dUmag;
			}
			Ux+=dUx;
			Uy+=dUy;
			ntry+=1;
		}
		else
			success=true;
	}while(!success && ntry<50);
	/*
	if(success){
		printf("SUCCESS!!!!! ntry=%d\n",ntry);
		int alpha,beta;
		double **Ttest,**Ttestboosted;
	Ttest=new double *[4];
	Ttestboosted=new double *[4];
	for(alpha=0;alpha<4;alpha++){
		Ttest[alpha]=new double[4];
		Ttestboosted[alpha]=new double[4];
	}
	FourVector u;	
		u[1]=Ux;
		u[2]=Uy;
		u[3]=0.0;
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		for(alpha=0;alpha<3;alpha++){
			for(beta=0;beta<3;beta++)
				Ttest[alpha][beta]=Ttestboosted[alpha][beta]=0.0;
		}
		for(alpha=0;alpha<3;alpha++){
			for(beta=0;beta<3;beta++)
				Ttest[alpha][beta]=Ttestboosted[alpha][beta]=0.0;
		}
		Ttest[0][0]=T00;
		Ttest[0][1]=Ttest[1][0]=T0x;
		Ttest[0][2]=Ttest[2][0]=T0y;
		Ttest[1][1]=Txx;
		Ttest[2][2]=Tyy;
		Ttest[1][2]=Ttest[2][1]=Txy;
		BoostToCM(u,Ttest,Ttestboosted);
		printf("SE tensor originally:\n");
		for(alpha=0;alpha<3;alpha++){
			for(beta=0;beta<3;beta++)
				printf("%8.4f ",Ttest[alpha][beta]);
			printf("\n");
		}
		printf("SE tensor boosted:\n");
		for(alpha=0;alpha<3;alpha++){
			for(beta=0;beta<3;beta++)
				printf("%8.4f ",Ttestboosted[alpha][beta]);
			printf("\n");
		}

	}
	*/
	gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
	epsilon=gamma*T00*gamma-2.0*gamma*T0x*Ux-2.0*gamma*T0y*Uy
	+Ux*Txx*Ux+Uy*Tyy*Uy+2.0*Ux*Txy*Uy;
}

int main(){
	int nparts=2000000;
	double Ux,Uy,epsilon,rho=1.0,weight,pdotu,umag,degen=3.0;
	CRandy *randy=new CRandy(time(NULL));
	double mass=139.57,mu,T=100.0;
	double T0x=0.0,T0y=0.0,T00=0.0;
	double Txx=0.0,Txy=0.0,Tyy=0.0;
	int ipart;
	FourVector p,u;
	u[1]=0.0;
	u[2]=0.0;
	u[3]=0.0;
	printf("Enter ux, uy: ");
	scanf("%lf %lf",&u[1],&u[2]);
	umag=sqrt(u[1]*u[1]+u[2]*u[2]);
	u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
	for(ipart=0;ipart<nparts;ipart++){
		randy->generate_boltzmann(mass,T,p);
		
		weight=(1.0+(u[1]*p[1]+u[2]*p[2])/(u[0]*p[0]));
		if(weight<randy->ran()){
			pdotu=p[1]*u[1]+p[2]*u[2];
			double oldpmag=sqrt(p[1]*p[1]+p[2]*p[2]);
			p[1]=p[1]-2.0*pdotu*u[1]/(umag*umag);
			p[2]=p[2]-2.0*pdotu*u[2]/(umag*umag);
			double newpmag=sqrt(p[1]*p[1]+p[2]*p[2]);
			if(fabs(newpmag-oldpmag)>1.0E-4)
				printf("newpmag=%g, oldpmag=%g\n",newpmag,oldpmag);
		}
		
		Misc::Boost(u,p);
		//printf("E=%g\n",p[0]);
		T00+=p[0];
		T0x+=p[1];
		T0y+=p[2];
		Txx+=p[1]*p[1]/p[0];
		Tyy+=p[2]*p[2]/p[0];
		Txy+=p[1]*p[2]/p[0];
	}
	T00*=u[0]*rho/double(nparts);
	T0x*=u[0]*rho/double(nparts);
	T0y*=u[0]*rho/double(nparts);
	Txx*=u[0]*rho/double(nparts);
	Tyy*=u[0]*rho/double(nparts);
	Txy*=u[0]*rho/double(nparts);
	TestGetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Ux,Uy,epsilon);
	printf("epsilon=%g, Ux=%g, Uy=%g\n",epsilon,Ux,Uy);

	double epsilon0,epsilon1,epsilon2,P,rho0,sigma2,dedT;
	CResList::freegascalc_onespecies(mass,T-0.1,epsilon0,P,rho0,sigma2,dedT);
	CResList::freegascalc_onespecies(mass,T+0.1,epsilon2,P,rho0,sigma2,dedT);
	CResList::freegascalc_onespecies(mass,T+0.1,epsilon1,P,rho0,sigma2,dedT);
	printf("dedT=%g =? %g\n",dedT,(epsilon2-epsilon0)/0.2);

	mu=3.6;
	TestGetMuT(mass,degen,rho,epsilon,T,mu);
	printf("T=%g, mu=%g\n",T,mu);
	return 0;
}
