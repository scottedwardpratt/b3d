#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int main(){
	int ix,iy,NR=30,ir,N;
	double delR=1.0,DXY=1.0,x,y,tau,r;
	double Ux,Uy,mu,T,rho,epsilon;
	FILE *fptr,*output;
	vector<double> muB,muK,mupi,Tpi,TK,TB,Upi,UK,UB,rhopi,rhoK,rhoB,epi,eK,eB;
	muB.resize(NR,0.0); muK.resize(NR,0.0); mupi.resize(NR,0.0);
	TB.resize(NR,0.0); TK.resize(NR,0.0); Tpi.resize(NR,0.0);
	Upi.resize(NR,0.0); UK.resize(NR,0.0); UB.resize(NR,0.0);
	rhopi.resize(NR,0.0); rhoK.resize(NR,0.0); rhoB.resize(NR,0.0);
	epi.resize(NR,0.0); eK.resize(NR,0.0); eB.resize(NR,0.0);
	vector<int> npts;
	npts.resize(NR,0);
	char filename[200],dummy[200];

	for(tau=1;tau<47.88888;tau+=1){
		

		//Pions
		for(ir=0;ir<NR;ir++){
			npts[ir]=0;
			mupi[ir]=Tpi[ir]=Upi[ir]=rhopi[ir]=epi[ir]=0.0;
		}

		sprintf(filename,"mucalc_results/mutinfo_pi_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon);
			if(!feof(fptr)){
				x=(ix+0.5)*DXY;
				y=(iy+0.5)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(floor(r/delR));
				if(ir<NR && N>0){
					npts[ir]+=1;
					mupi[ir]+=mu;
					Tpi[ir]+=T;
					Upi[ir]+=sqrt(Ux*Ux+Uy*Uy);
					rhopi[ir]=rho;
					epi[ir]=epsilon;
				}
			}
		}while(!feof(fptr));
		fclose(fptr);

		sprintf(filename,"mucalc_results/muTvsR_pi_tau%g.txt",tau);
		output=fopen(filename,"w");
		for(ir=0;ir<NR;ir++){
			if(npts[ir]>0){
				mupi[ir]=mupi[ir]/double(npts[ir]);
				Upi[ir]=Upi[ir]/double(npts[ir]);
				Tpi[ir]=Tpi[ir]/double(npts[ir]);
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,npts[ir],Tpi[ir],Upi[ir],mupi[ir],rhopi[ir],epi[ir]);
			}
			else{
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,0,0.0,0.0,0.0,0.0,0.0);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Kaons
		for(ir=0;ir<NR;ir++){
			npts[ir]=0;
			muK[ir]=TK[ir]=UK[ir]=rhoK[ir]=eK[ir]=0.0;
		}
		
		sprintf(filename,"mucalc_results/mutinfo_K_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon);
			//printf("ix=%d, iy=%d\n",ix,iy);
			if(!feof(fptr)){
				x=(ix+0.5)*DXY;
				y=(iy+0.5)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(floor(r/delR));
				if(ir<NR && N>0){
					//printf("ir=%d\n",ir);
					npts[ir]+=1;
					muK[ir]+=mu;
					TK[ir]+=T;
					UK[ir]+=sqrt(Ux*Ux+Uy*Uy);
					rhoK[ir]=rho;
					eK[ir]=epsilon;
				}
			}
		}while(!feof(fptr));
		fclose(fptr);

		sprintf(filename,"mucalc_results/muTvsR_K_tau%g.txt",tau);
		output=fopen(filename,"w");
		for(ir=0;ir<NR;ir++){
			if(npts[ir]>0){
				muK[ir]=muK[ir]/double(npts[ir]);
				UK[ir]=UK[ir]/double(npts[ir]);
				TK[ir]=TK[ir]/double(npts[ir]);
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,npts[ir],TK[ir],UK[ir],muK[ir],rhoK[ir],eK[ir]);
			}
			else{
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,0,0.0,0.0,0.0,0.0,0.0);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Baryons and Hyperons

		for(int btype=0;btype<8;btype++){
			for(ir=0;ir<NR;ir++){
				npts[ir]=0;
				muB[ir]=TB[ir]=UB[ir]=rhoB[ir]=eB[ir]=0.0;
			}
			sprintf(filename,"mucalc_results/mutinfo_B%d_tau%g.txt",btype,tau);
			fptr=fopen(filename,"r");
			fgets(dummy,200,fptr);
			do{
				fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon);
				if(!feof(fptr)){
					x=(ix+0.5)*DXY;
					y=(iy+0.5)*DXY;
					r=sqrt(x*x+y*y);
					ir=lrint(floor(r/delR));
					if(ir<NR && N>0){
						npts[ir]+=1;
						muB[ir]+=mu;
						TB[ir]+=T;
						UB[ir]+=sqrt(Ux*Ux+Uy*Uy);
						rhoB[ir]=rho;
						eB[ir]=epsilon;
					}
				}
			}while(!feof(fptr));
			fclose(fptr);


			sprintf(filename,"mucalc_results/muTvsR_B%d_tau%g.txt",btype,tau);
			output=fopen(filename,"w");
			for(ir=0;ir<NR;ir++){
				if(npts[ir]>0){
					muB[ir]=muB[ir]/double(npts[ir]);
					TB[ir]=TB[ir]/double(npts[ir]);
					UB[ir]=UB[ir]/double(npts[ir]);
					r=(ir+0.5)*delR;
					fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,npts[ir],TB[ir],UB[ir],muB[ir],rhoB[ir],eB[ir]);
				}
				else{
					r=(ir+0.5)*delR;
					fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f %7.4f\n",r,0,0.0,0.0,0.0,0.0,0.0);
				}
				npts[ir]=0;
			}
			fclose(output);
		}

	}
	return 0;
}
