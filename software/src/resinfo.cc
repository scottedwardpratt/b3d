#include <cmath>
#include <cstdlib>
#include "resonances.h"
#include "constants.h"
#include "parametermap.h"
#include "randy.h"
#include "resonances.h"

using namespace std;
//using namespace boost::math;

CResList *CResInfo::reslist=NULL;
double **CResInfo::ChiA=NULL;
CRandy *CResInfo::randy=new CRandy(-1234);
char *CResInfo::message=new char[500];

CResInfo::CResInfo(){
	minmass=0.0;
	branchlist.clear();
	netchi=0.0;
}

void CResInfo::DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	double r,bsum;
	int ibody,ibranch;
	CBranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	r=randy->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			sprintf(message,"In DecayGetResInfo: bsum too large, = %g\n",bsum);
			CLog::Fatal(message);
		}
	}while(bsum<r);
	nbodies=bptr->resinfoptr.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfoptr[ibody];
	}
}

bool CResInfo::CheckForDaughters(int codecheck){
	//checks to see if any decay daughters match code, for code=0, checking for charged parts
	int ibody,nbodies,ibranch;
	bool exists=false;
	CResInfo *daughter;
	CBranchInfo *bptr;
	CBranchList::iterator bpos;
	if(codecheck!=0){
		if(code==codecheck){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->code==codecheck){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				bpos++;
			}while(bpos!=branchlist.end());
		}
	}
	else{
		if(charge!=0){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->charge!=0){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				++bpos;
			}while(bpos!=branchlist.end());
		}
	}
	return exists;
}

void CResInfo::DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	nbodies=bptr_minmass->resinfoptr.size();
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr_minmass->resinfoptr[ibody];
	}
}

bool CResInfo::CheckForNeutral(){
	bool neutral=true;
	if(charge!=0 || strange!=0 || baryon!=0)
		neutral=false;
	return neutral;
}

double CResInfo::GenerateMass(){
	double m,r,weight=0.0,rho,k,lor,gamma;
	double alpha=reslist->RESWIDTH_ALPHA;;
	if(decay){
		double m1=branchlist[0]->resinfoptr[0]->mass;
		double m2=0.0;
		for(int n=1;n<int(branchlist[0]->resinfoptr.size());n++){
			m2+=branchlist[0]->resinfoptr[n]->mass;
		}
		double kr=sqrt(pow((mass*mass-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*mass);
		do{
			r = randy->ran();
			m = ((width/2)*tan(PI*(r - .5))) + mass;
			if ((m < (m1+m2))||(m>2.0*mass)) continue;
			k=sqrt(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*m);
			gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
			rho=(2.0/(width*PI))*(0.25*gamma*gamma)/((0.25*gamma*gamma)+(mass-m)*(mass-m));
			lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-m,2.0));
			weight = rho/(lor*8.0);
			r=randy->ran();
		}while(r>weight || m<minmass);
	}
	else
		m=mass;
	return m;
}

double CResInfo::GenerateThermalMass(double maxweight, double T){
	double m,kr,weight,r1,r2,k,gamma,rho,lor,k2,m1,m2,k2mr;
	double alpha=reslist->RESWIDTH_ALPHA;
	if(decay || maxweight<0.0){
		m1=branchlist[0]->resinfoptr[0]->mass;
		m2=0.0;
		for(int n=1;n<int(branchlist[0]->resinfoptr.size());n++){
			m2+=branchlist[0]->resinfoptr[n]->mass;
		}
		k2mr = gsl_sf_bessel_Kn(2,mass/T); // K2 for resmass
		kr=pow(mass*mass-m1*m1-m2*m2,2) - (4*m1*m1*m2*m2);
		kr = (1/(2*mass))*sqrt(kr); // k at resonant mass
		int i = 0; // for use in while loop
		while(i==0){
			r1 = randy->ran(); // get random numbers
			r2 = randy->ran(); // between [0, 1]
			m = ((width/2)*tan(PI*(r1 - .5))) + mass;// generate random mass value proportional to the lorentz distribution
			if ((m < minmass) ) continue;
			// throw out values out of range
			k=sqrt(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*m);
			gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
			rho=(2.0/(width*PI))*(0.25*gamma*gamma)/((0.25*gamma*gamma)+(mass-m)*(mass-m));
			lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-m,2.0));
			k2 = gsl_sf_bessel_Kn(2,(m/T)); // K2 value
			weight = rho*k2*m*m/(lor*k2mr*mass*mass*maxweight);
			if (r2 < weight) i = 1; // success
		}
	}
	else m=mass;
	if(m<minmass){
		sprintf(message,"In GenerateThermalMass, m=%g, but minmass=%g\n",m,minmass);
		CLog::Fatal(message);
	}
	return m; 
}

void CResInfo::Print(){
	sprintf(message,"+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",code,mass,minmass,name.c_str());
	CLog::Info(message);
	sprintf(message,"Gamma=%g, Spin=%g, Decay=%d\n",width,spin,int(decay));
	CLog::Info(message);
	sprintf(message,"Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
	CLog::Info(message);
}

void CResInfo::FindFinalProducts(double taumax){
	// all decay products -- decaying stops if tau_decay > taumax
	finalproductslist.clear();
	CBranchList blist;
	CBranchInfo *bptr,*fbptr;
	CResInfo *resinfo;
	bool foundsplit=false;
	long unsigned int ibranch,iibranch,iires,iiires;
	double netbranching=0.0;
	finalproductslist.clear();
	bptr=new CBranchInfo();
	if(decay && (HBARC/width)<taumax){
		finalproductslist=branchlist;
		ibranch=0;
		foundsplit=true;
		ibranch=0;
		while(ibranch<finalproductslist.size()){
			foundsplit=false;
			bptr=finalproductslist[ibranch];
			iires=0;
			do{
				if(bptr->resinfoptr[iires]->decay && (HBARC/bptr->resinfoptr[iires]->width)<taumax){
					foundsplit=true;
					resinfo=bptr->resinfoptr[iires];
					for(iibranch=1;iibranch<resinfo->finalproductslist.size();iibranch++){
						fbptr=new CBranchInfo();
						fbptr->Copy(bptr);
						finalproductslist.push_back(fbptr);
						fbptr->resinfoptr[iires]=resinfo->finalproductslist[iibranch]->resinfoptr[0];
						for(iiires=1;iiires<resinfo->finalproductslist[iibranch]->resinfoptr.size();iiires++){
							fbptr->resinfoptr.push_back(resinfo->finalproductslist[iibranch]->resinfoptr[iiires]);
						}
						fbptr->branching=bptr->branching*resinfo->finalproductslist[iibranch]->branching;
						//bcheck+=resinfo->finalproductslist[iibranch]->branching;
						//netbranching+=fbptr->branching;
					}
					iibranch=0;
					bptr->resinfoptr[iires]=resinfo->finalproductslist[iibranch]->resinfoptr[0];
					for(iiires=1;iiires<resinfo->finalproductslist[iibranch]->resinfoptr.size();iiires++){
						bptr->resinfoptr.push_back(resinfo->finalproductslist[iibranch]->resinfoptr[iiires]);
					}
					bptr->branching=bptr->branching*resinfo->finalproductslist[iibranch]->branching;
				}
				iires+=1;
			}while(iires<bptr->resinfoptr.size());
			if(!foundsplit)
				ibranch+=1;
		}
	}
	
	netbranching=0.0;
	for(ibranch=0;ibranch<finalproductslist.size();ibranch++){
		netbranching+=finalproductslist[ibranch]->branching;
		int netq[3]={0};
		for(iires=0;iires<finalproductslist[ibranch]->resinfoptr.size();iires++){
			for(int iq=0;iq<3;iq++)
				netq[iq]+=finalproductslist[ibranch]->resinfoptr[iires]->q[iq];
		}
		double netB=0.0,netQ=0.0;
		netB=double(q[0]+q[1]+q[2])/3.0-double(netq[0]+netq[1]+netq[2])/3.0;
		netQ=(double(2.0*q[0]-q[1]-q[2])/3.0)-double(2*netq[0]-netq[1]-netq[2])/3.0;
		if(fabs(netB)>1.0E-7 || fabs(netQ)>1.0E-7){
			sprintf(message,"In FindFinalProducs, charge not conserved for %5d: netQ=%g, netB=%g\n",code,netQ,netB);
			CLog::Fatal(message);
		}
	}
	if(decay && (HBARC/width)<taumax && fabs(netbranching-1.0)>1.0E-5){
		sprintf(message,"oops, netbranching for final states=%g, pid=%d\n",netbranching,code);
		CLog::Fatal(message);
	}
}

bool CResInfo::FindContent(int codecheck,double weight0,double taumax,double &weight){
	// finds how many hadrons of type codecheck result from decays
	bool foundpart=false;
	CBranchInfo *bptr;
	CResInfo *resinfo1;
	unsigned long int ibranch,ibody1;
	if(!reslist->finalproductsfound)
		reslist->FindFinalProducts(taumax);
	weight=0.0;
	if(decay && (HBARC/width)<taumax){
		for(ibranch=0;ibranch<finalproductslist.size();ibranch++){
			bptr=finalproductslist[ibranch];
			for(ibody1=0;ibody1<bptr->resinfoptr.size();ibody1++){
				resinfo1=bptr->resinfoptr[ibody1];
				if(resinfo1->code==codecheck){
					weight+=weight0*bptr->branching;
					foundpart=true;
				}
			}
		}
	}
	else{
		if(code==codecheck){
			weight=weight0;
			foundpart=true;
		}
	}
	return foundpart;	
}

bool CResInfo::FindContentPairs(int codecheck1,int codecheck2,double weight0,double taumax,double &weight){
	// finds how many hadrons of type codecheck result from decays
	bool foundpair=false;
	CBranchInfo *bptr;
	CResInfo *resinfo1,*resinfo2;
	unsigned long int ibranch,ibody1,ibody2;
	if(!reslist->finalproductsfound)
		reslist->FindFinalProducts(taumax);
	weight=0.0;
	for(ibranch=0;ibranch<finalproductslist.size();ibranch++){
		bptr=finalproductslist[ibranch];
		for(ibody1=0;ibody1<bptr->resinfoptr.size();ibody1++){
			resinfo1=bptr->resinfoptr[ibody1];
			if(abs(resinfo1->code)==abs(codecheck1)){
				for(ibody2=0;ibody2<bptr->resinfoptr.size();ibody2++){
					if(ibody1!=ibody2){
						resinfo2=bptr->resinfoptr[ibody2];
						if(abs(resinfo2->code)==abs(codecheck2)){
							int sign=1;
							if(resinfo1->code*resinfo2->code<0)
								sign=-1;
							weight+=weight0*bptr->branching*sign;
							foundpair=true;
						}
					}
				}
			}
		}
	}
	return foundpair;
}

void CResInfo::PrintFinalProducts(){
	CBranchInfo *bptr;
	CResInfo *resinfo1;
	unsigned long int ibranch,ibody1;
	sprintf(message,"Final Decay Products of %s:\n",name.c_str());
	double netbranching=0.0;
	for(ibranch=0;ibranch<finalproductslist.size();ibranch++){
		bptr=finalproductslist[ibranch];
		sprintf(message,"%s____ BRANCH %ld _____ branching=%g\n",message,ibranch,bptr->branching);
		netbranching+=bptr->branching;
		for(ibody1=0;ibody1<bptr->resinfoptr.size();ibody1++){
			resinfo1=bptr->resinfoptr[ibody1];
			sprintf(message,"%s%6d ",message,resinfo1->code);
		}
		sprintf(message,"%s\n",message);	
	}
	sprintf(message,"%s---- NET BRANCHING=%g =? 1.0 ----\n",message,netbranching);
	CLog::Info(message);
}

void CResInfo::SetBtype(){
	Btype=-1;
	if(abs(baryon)!=0){
		Btype=abs(strange);
		if(spin>1.0)
			Btype+=4;
		if(Btype==1){
			if(name[0]=='L')
				Btype=3;
		}
		//if(Btype==4){
		//	printf("%2d:%s\n",Btype,name.c_str());
		//}
	}
}
