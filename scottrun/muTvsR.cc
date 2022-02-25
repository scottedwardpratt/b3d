#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int main(){
	int ix,iy,itau,NR=20,ir;
	double delR=1.5,DXY=1.5,x,y,tau,r;
	double Ux,Uy,mu,T,muS;
	int N,NS;
	FILE *fptr,*output;
	vector<double> muB,muBS,muK,mupi,Tpi,TK,TB,Upi,UK,UB;
	muB.resize(NR,0.0); muBS.resize(NR,0.0); muK.resize(NR,0.0); mupi.resize(NR,0.0);
	TB.resize(NR,0.0); TK.resize(NR,0.0); Tpi.resize(NR,0.0);
	Upi.resize(NR,0.0); UK.resize(NR,0.0); UB.resize(NR,0.0);
	vector<int> npts;
	npts.resize(NR,0);
	char filename[200],dummy[200];

	for(itau=10;itau<61;itau+=10){
		tau=0.5*itau;

		//Pions

		sprintf(filename,"mucalc_results/mutinfo_pi_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu);
			if(!feof(fptr)){
				x=-36.0+(0.5+ix)*DXY;
				y=-36.0+(0.5+iy)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(r/delR);
				if(ir<NR){
					npts[ir]+=1;
					mupi[ir]+=mu;
					Tpi[ir]+=T;
					Upi[ir]+=sqrt(Ux*Ux+Uy*Uy);
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
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f\n",r,npts[ir],Tpi[ir],Upi[ir],mupi[ir]);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Kaons
		
		sprintf(filename,"mucalc_results/mutinfo_K_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu);
			//printf("ix=%d, iy=%d\n",ix,iy);
			if(!feof(fptr)){
				x=-36.0+(0.5+ix)*DXY;
				y=-36.0+(0.5+iy)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(r/delR);
				if(ir<NR){
					//printf("ir=%d\n",ir);
					npts[ir]+=1;
					muK[ir]+=mu;
					TK[ir]+=T;
					UK[ir]+=sqrt(Ux*Ux+Uy*Uy);
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
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f\n",r,npts[ir],TK[ir],UK[ir],muK[ir]);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Baryons and Hyperons

		sprintf(filename,"mucalc_results/mutinfo_B_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf",&ix,&iy,&N,&NS,&T,&Ux,&Uy,&mu,&muS);
			if(!feof(fptr)){
				x=-36.0+(0.5+ix)*DXY;
				y=-36.0+(0.5+iy)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(r/delR);
				if(ir<NR){
					npts[ir]+=1;
					muB[ir]+=mu;
					muBS[ir]+=muS;
					TB[ir]+=T;
					UB[ir]+=sqrt(Ux*Ux+Uy*Uy);
				}
			}
		}while(!feof(fptr));
		fclose(fptr);

		sprintf(filename,"mucalc_results/muTvsR_B_tau%g.txt",tau);
		output=fopen(filename,"w");
		for(ir=0;ir<NR;ir++){
			if(npts[ir]>0){
				muB[ir]=muB[ir]/double(npts[ir]);
				muBS[ir]=muBS[ir]/double(npts[ir]);
				TB[ir]=TB[ir]/double(npts[ir]);
				UB[ir]=UB[ir]/double(npts[ir]);
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.2f %7.4f %7.4f %7.4f\n",r,npts[ir],TB[ir],UB[ir],muB[ir],muBS[ir]);
			}
			npts[ir]=0;
		}
		fclose(output);

	}
	return 0;
}
