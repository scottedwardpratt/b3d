#include "b3d.h"
#include "randy.h"
#include "misc.h"
#include "resonances.h"
#include "commondefs.h"

using namespace std;

void TestGetMuT_Baryon(double rhoB_target,double rhoBS_target,double epsilon_target,double &T,double &muB,double &muBS){
	Eigen::Vector3d rho,rho_target,x,dx,drho;
	Eigen::Matrix3d drhodx;
	const int NSPECIES=2;
	vector<string> name{"nucleon","Lambda","Cascade"};
	vector<double> mass{939.0,1116.0,1314.86};
	vector<int> degen{8,4,8};
	vector<int> nstrange{0,1,2};
	rho_target(0)=epsilon_target/rhoB_target;
	rho_target(1)=rhoB_target;
	rho_target(2)=rhoBS_target;
	double factor,dfactordxB,dfactordxBS;
	double xB,xBS,e,dedt,dedx1,dedx2;
	double P,sigma2,epsilon0,rho0,dedt0,accuracy,D;
	int ntry=0,ispecies,nmax=200;

	x[0]=T;
	x[1]=muB;
	x[2]=muBS;

	//printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	//printf("rho_target=(%g,%g,%g)\n",rho_target(0),rho_target(1),rho_target(2));

	do{
		ntry+=1;
		T=x(0);
		xB=x(1);
		xBS=x(2);
		rho.setZero();
		drhodx.setZero();
		e=dedt=dedx1=dedx2=0.0;
		for(ispecies=0;ispecies<NSPECIES;ispecies++){
			CResList::freegascalc_onespecies(mass[ispecies],T,epsilon0,P,rho0,sigma2,dedt0);

			factor=degen[ispecies]*exp(x(1))*exp(x(2)*nstrange[ispecies]);
			dfactordxB=factor;
			dfactordxBS=factor*nstrange[ispecies];

			e+=epsilon0*factor;
			rho(1)+=rho0*factor;
			rho(2)+=rho0*factor*nstrange[ispecies];

			dedt+=dedt0*factor;
			dedx1+=epsilon0*dfactordxB;
			dedx2+=epsilon0*dfactordxBS;

			drhodx(1,0)+=(epsilon0/(T*T))*factor;
			drhodx(1,1)+=rho0*dfactordxB;
			drhodx(1,2)+=rho0*dfactordxBS;

			drhodx(2,0)+=(epsilon0/(T*T))*factor*nstrange[ispecies];
			drhodx(2,1)+=rho0*dfactordxB*nstrange[ispecies];
			drhodx(2,2)+=rho0*dfactordxBS*nstrange[ispecies];
		}
		rho(0)=e/rho(1);
		drhodx(0,0)=(dedt/rho(1))-e*drhodx(1,0)/(rho(1)*rho(1));
		drhodx(0,1)=(dedx1/rho(1))-e*drhodx(1,1)/(rho(1)*rho(1));
		drhodx(0,2)=(dedx2/rho(1))-e*drhodx(1,2)/(rho(1)*rho(1));

		drho=rho_target-rho;
		dx=drhodx.colPivHouseholderQr().solve(drho);

		/*
		printf("rho_target=(%g,%g,%g)\n",rho_target(0),rho_target(1),rho_target(2));
		printf("rho=(%g,%g,%g)\n",rho(0),rho(1),rho(2));
		printf("dx=(%g,%g,%g)\n",dx(0),dx(1),dx(2));
		printf("x=(%g,%g,%g)\n",x(0),x(1),x(2));
		printf("TB=%g, muB=%g, muBS=%g\n",x[0],x[1],x[2]);
		printf("-----------------------------------------\n");
		*/
		if(dx(0)!=dx(0) || dx(1)!=dx(1) || dx(2)!=dx(2)){
			printf("FAILURE, ntry=%d, dx!=dx\n",ntry);
			printf("rho_target=(%g,%g,%g)\n",rho_target(0),rho_target(1),rho_target(2));
			printf("rho=(%g,%g,%g)\n",rho(0),rho(1),rho(2));
			printf("dx=(%g,%g,%g)\n",dx(0),dx(1),dx(2));
			printf("x=(%g,%g,%g)\n",x(0),x(1),x(2));
			printf("TB=%g, muB=%g, muBS=%g\n",x[0],x[1],x[2]);
			printf("drhodx=\n");
			cout << drhodx << endl;
			exit(1);
		}
		
		if(x(0)+dx(0)<0.3*T){
			D=fabs(dx(0));
			dx=dx*(0.3*T/D);
		}
		if(dx(0)>20.0){
			D=fabs(dx(0));
			dx=dx*20.0/D;
		}
		
		double dxmax=1.0;
		if(fabs(dx(1))>dxmax){
			D=fabs(dx(1));
			dx=dx*dxmax/D;
		}
		if(fabs(dx(2))>dxmax){
			D=fabs(dx(2));
			dx=dx*dxmax/D;
		}

		

		//printf("----------------------------------------\n");
		//printf("rho=(%g,%g,%g)\n",rho(0),rho(1),rho(2));
		//printf("dx=(%g,%g,%g)\n",dx(0),dx(1),dx(2));
		x=x+dx;
		//printf("x=(%g,%g,%g)\n",x(0),x(1),x(2));

		accuracy=0.001*drho[0]*drho[0]+drho[1]*drho[1]+drho[2]*drho[2];

	}while(accuracy>1.0E-10 && ntry<nmax);
	if(ntry==nmax || accuracy!=accuracy){
		printf("FAILURE, ntry=nmax\n");
		printf("rho_target=(%g,%g,%g)\n",rho_target(0),rho_target(1),rho_target(2));
		printf("rho=(%g,%g,%g)\n",rho(0),rho(1),rho(2));
		printf("dx=(%g,%g,%g)\n",dx(0),dx(1),dx(2));
		printf("x=(%g,%g,%g)\n",x(0),x(1),x(2));
		printf("TB=%g, muB=%g, muBS=%g\n",x[0],x[1],x[2]);
		exit(1);
	}
	T=x[0];
	muB=x[1];
	muBS=x[2];
	printf("ntry=%d, accuracy=%g\n",ntry,accuracy);
	//printf("TB=%g, muB=%g, muBS=%g\n",T,muB,muBS);
}


int main(){
	double epsilon,T=100.0,muB,muBS,rhoB,rhoBS;
	CRandy randy(1234);
	for(int itest=0;itest<10000000;itest++){
		rhoB=0.0+0.5*randy.ran();
		rhoBS=0.5*rhoB*randy.ran();
		epsilon=rhoB*(15.0+939.0+300*randy.ran())+rhoBS*177.0;
		muB=muBS=0.0;
		T=100.0;
		TestGetMuT_Baryon(rhoB,rhoBS,epsilon,T,muB,muBS);
		printf("%6d: T=%7.3f, muB=%5.3f, muBS=%5.3f\n",itest,T,muB,muBS);
	}
	return 0;
}
