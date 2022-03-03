#include <cmath>
#include <cstdlib>
#include "resonances.h"
#include "constants.h"
#include "parametermap.h"
#include "randy.h"

using namespace std;
//using namespace boost::math;

CSampler* CResList::sampler=NULL;
CB3D* CResList::b3d=NULL;

CResList::CResList(CparameterMap* parmap_in){
	parmap=parmap_in;
	RESWIDTH_ALPHA=parmap->getD("RESONANCE_ALPHA",0.5);
	RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	USEPOLEMASS=parmap->getB("USEPOLEMASS",false);
	CResInfo::reslist=this;
	finalproductsfound=false;
	ReadResInfo();
}

CResList::~CResList(){
	CResInfo *resinfo;
	CResInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		resmap.erase(rpos);
		delete resinfo;
		rpos=resmap.begin();
	}
}

CMerge::CMerge(CResInfo *resinfo_in,double branching_in, int L_in){
	resinfo=resinfo_in;
	branching=branching_in;
	L = L_in;
	next=NULL;
}

CBranchInfo::CBranchInfo(){
}

void CResList::freegascalc_onespecies(double m,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	char message[100];
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=T*T;
	t3=t2*T;
	z=m/T;
	if(z>1000.0){
		P=epsilon=dens=dedt=0.0;
		sprintf(message,"freegascalc_onespecies: z is huge=%g, m=%g, t=%g\n",z,m,T);
		CLog::Info(message);
	}
	else{
		if(z<0.0){
			sprintf(message,"freegascalc_onespecies: negative z=%g,m=%g,T=%g ___\n",z,m,T);
			CLog::Info(message);
		}
		k0=gsl_sf_bessel_Kn(0,z);
		k1=gsl_sf_bessel_Kn(1,z);
		P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
		dens=P/T;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
		Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(T,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

void CResList::freegascalc_onespecies_finitewidth(double resmass, double m1, double m2, double T, double width,
double &epsilon,double &P,double &dens,double &sigma2,
double &dedt,double &maxweight){

	double kr,k,E,E0,dE,gamma,rho,percent,dp,closest;
	double sum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
	double n0,resn0,lor,weight;
	double alpha=RESWIDTH_ALPHA;
	dE=1.0;
	percent=0.001;
	dp=0.002;
	if (width<1.0){
		dE=0.1;
		percent=0.0005;
		dp=0.001;
	}
	closest=1000.0;
	E0=m1+m2;
	maxweight=-1.0;
	kr=sqrt(pow((resmass*resmass-m1*m1-m2*m2),2.0)-4.0*m1*m1*m2*m2)/(2.0*resmass);

	for(E=(m1+m2+0.5*dE);E<2.0*resmass;E+=dE){

		k=sqrt(pow((E*E-m1*m1-m2*m2),2.0)-(4.0*m1*m1*m2*m2))/(2.0*E);
		gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
		rho=(2.0/(width*PI))*0.25*gamma*gamma/((0.25*gamma*gamma)+(resmass-E)*(resmass-E));
		sum+=rho*dE;

		n0=gsl_sf_bessel_Kn(2,E/T)*E*E*T/(2*PI*PI*pow(HBARC,3.0));
		resn0=gsl_sf_bessel_Kn(2,resmass/T)*resmass*resmass*T/(2*PI*PI*pow(HBARC,3.0));
		lor=(width/(2.0*PI))/(0.25*width*width+(resmass-E)*(resmass-E));
		weight=n0*rho/(resn0*lor);

		if(weight>maxweight) maxweight=weight;
		if (abs(sum-percent)<closest) closest=abs(sum-percent);
		else{
			freegascalc_onespecies(E,T,epsilon,P,dens,sigma2,dedt);
			esum+=epsilon*rho*(E-E0);
			psum+=P*rho*(E-E0);
			dsum+=dens*rho*(E-E0);
			sigsum+=sigma2*rho*(E-E0);
			dedtsum+=dedt*rho*(E-E0);
			closest=1000.0;
			percent+=dp;
			E0=E;
		}
	}

	epsilon=esum/sum;
	P=psum/sum;
	dens=dsum/sum;
	sigma2=sigsum/sum;
	dedt=dedtsum/sum;
}

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,NResonances;
	double mothermass,bsum,netm,qR2,bmax;
	int ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,LDecay;
	int netq,netb,nets;
	string name, filename;
	CResInfo *resinfoptr=NULL;
	CBranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	FILE * decayinfofile;
	char dummy[200],cname[200];
	filename=parmap->getS("B3D_RESONANCES_INFO_FILE",string("../resinfo/resonances_standardhadrons.txt"));
	sprintf(message,"will read res info from %s\n",filename.c_str());
	CLog::Info(message);
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	MergeArray=new CMerge **[NResonances];
	SigmaMaxArray=new double *[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new CMerge *[NResonances];
		SigmaMaxArray[ires]=new double[NResonances];
		for(jres=0;jres<NResonances;jres++){
			MergeArray[ires][jres]=NULL;
			SigmaMaxArray[ires][jres]=0.0;
		}
	}
	for(ires=0;ires<NResonances;ires++){
		resinfoptr=new CResInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfoptr->code,&resinfoptr->mass,&resinfoptr->charge,&resinfoptr->baryon, &resinfoptr->strange,&resinfoptr->spin,&resinfoptr->G_Parity,&decay,&resinfoptr->width);

		resinfoptr->q[0]=resinfoptr->baryon+resinfoptr->charge;
		resinfoptr->q[1]=2.0*resinfoptr->baryon+resinfoptr->strange-resinfoptr->charge;
		resinfoptr->q[2]=-resinfoptr->strange;
		resinfoptr->netchi0+=resinfoptr->charge*resinfoptr->charge;
		resinfoptr->Btype=-1;
		if(!resinfoptr->decay)
			resinfoptr->netchi+=abs(resinfoptr->charge);
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfoptr->name=cname;
		resinfoptr->name.erase(0,resinfoptr->name.find_first_not_of(' '));
		resinfoptr->name.erase(resinfoptr->name.find_last_not_of(' ')+1);
		resinfoptr->SetBtype();
		resinfoptr->decay=bool(decay);
		if(!RESONANCE_DECAYS)
			resinfoptr->decay=false;
		resinfoptr->ires=ires;
		resinfoptr->branchlist.clear();
		resmap.insert(CResInfoPair(resinfoptr->code,resinfoptr));
	} 
	fclose(resinfofile);

	filename=parmap->getS("B3D_RESONANCES_DECAYS_FILE",string("../resinfo/decays_pdg_weak.txt"));
	sprintf(message,"will read decay info from %s\n",filename.c_str());
	CLog::Info(message);
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfoptr=GetResInfoPtr(mothercode);
		resinfoptr->minmass=1.0E10;
		bsum=0.0;
		bmax=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			bptr->resinfoptr.clear();
			resinfoptr->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfoptr.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfoptr[ibody]->charge;
				netb+=bptr->resinfoptr[ibody]->baryon;
				nets+=bptr->resinfoptr[ibody]->strange;
				netm+=bptr->resinfoptr[ibody]->mass;
			} 
			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				sprintf(message,"Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				sprintf(message,"MOTHER=%d: (ichannel=%d, nbodies=%d):\n",resinfoptr->code,ichannel,nbodies);
				sprintf(message,"DAUGHTERS: ");
				for(ibody=0;ibody<nbodies;ibody++){
					sprintf(message,"%s %d ",message,bptr->resinfoptr[ibody]->code);
				}
				sprintf(message,"%s\n",message);
				CLog::Fatal(message);
			}

			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			for(ibody=0;ibody<nbodies;ibody++){
				if(!bptr->resinfoptr[ibody]->decay)
					resinfoptr->netchi+=abs(bptr->resinfoptr[ibody]->charge)*bptr->branching;
			}
			//store two body decays only
			if(nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				if(resinfoptr->mass>bptr->resinfoptr[0]->mass+bptr->resinfoptr[1]->mass){
					qR2=triangle(resinfoptr->mass,
					bptr->resinfoptr[0]->mass,bptr->resinfoptr[1]->mass);
					SigmaMaxArray[ires1][ires2]+=
						(bptr->branching*4.0*PI*HBARC*HBARC)*(2.0*resinfoptr->spin+1.0)/
							((2.0*bptr->resinfoptr[0]->spin+1.0)*(2.0*bptr->resinfoptr[1]->spin+1.0)*qR2);
					if(ires2!=ires1)
						SigmaMaxArray[ires2][ires1]=SigmaMaxArray[ires1][ires2];
				}
			}
			bsum+=bptr->branching;
			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfoptr->minmass){
				resinfoptr->minmass=netm;
				resinfoptr->bptr_minmass=bptr;
			}
			// switch places to make sure first branch has largest 
			if(bptr->branching>bmax){
				bmax=bptr->branching>bmax;
				if(ichannel>0){
					firstbptr=resinfoptr->branchlist[0];
					resinfoptr->branchlist[0]=bptr;
					resinfoptr->branchlist[ichannel]=firstbptr;
				}
			}
		}
	}
	fclose(decayinfofile);
}

CResInfo* CResList::GetResInfoPtr(int code){
	CResInfoMap::iterator rpos;
	rpos=resmap.find(code);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		sprintf(message,"Warning GetResInfoPtr() can't find match for PID=%d\n",code);
		CLog::Fatal(message);
		return NULL;
	}
}

void CResList::CalcConductivity(double T,double &P,double &epsilon,double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi,Eigen::Matrix3d &sigma){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,m1=0.0,m2=0.0,degen,s;
	double width,minmass,maxweighti;
	double pi,epsiloni,densi,sigma2i,dedti,Ji;
	Eigen::Matrix3d sigmai(3,3);
	int a,b,n,ires;
	chi.setZero();
	sigma.setZero();
	P=epsilon=s=nh=0.0;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		ires=resinfoptr->ires;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			if(USEPOLEMASS){
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti,Ji);
				maxweighti=-1.0;
			}
			else{
				width=resinfoptr->width;
				minmass=resinfoptr->minmass;
				if(resinfoptr->decay){
					m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
					m2=0.0;
					for(n=1;n<int(resinfoptr->branchlist[0]->resinfoptr.size());n++){
						m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
					}
				}
				if((minmass>0.0)&&(width>1.0E-3)){
					freegascalc_onespecies_finitewidth(m,m1,m2,T,width,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
				}
				else{
					freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti,Ji);
					maxweighti=-1.0;
				}
			}
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			nh+=densi*degen;
			density[ires]=densi*degen;
			maxweight[ires]=maxweighti;
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi(a,b)+=densi*degen*resinfoptr->q[a]*resinfoptr->q[b];
					sigma(a,b)+=Ji*degen*resinfoptr->q[a]*resinfoptr->q[b];
				}
			}
		}
		else{
			density[ires]=0.0;
		}
	}
}

void CResList::CalcEoSandChi(double T,double &P,double &epsilon,
double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,m1=0.0,m2=0.0,degen,s;
	double width,minmass,maxweighti;
	double pi,epsiloni,densi,sigma2i,dedti;
	double netchi=0.0,netchi0=0.0;
	int a,b,n,ires,nres=resmap.size();
	chi.setZero();
	P=epsilon=s=nh=0.0;
	density.resize(nres);
	maxweight.resize(nres);
	for(ires=0;ires<nres;ires++){
		density[ires]=0.0;
	}
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		ires=resinfoptr->ires;
		
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			if(USEPOLEMASS){
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
				maxweighti=-1.0;
			}
			else{
				width=resinfoptr->width;
				minmass=resinfoptr->minmass;
				if(resinfoptr->decay){
					m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
					m2=0.0;
					for(n=1;n<int(resinfoptr->branchlist[0]->resinfoptr.size());n++){
						m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
					}
				}
				if((minmass>0.0)&&(width>1.0E-3)){
					freegascalc_onespecies_finitewidth(m,m1,m2,T,width,
					epsiloni,pi,densi,sigma2i,dedti,maxweighti);
					if(densi!=densi){
						resinfoptr->Print();
						exit(1);
					}
				}
				else{
					freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
					maxweighti=-1.0;
				}
			}
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			nh+=densi*degen;
			density[ires]=densi*degen;
			maxweight[ires]=maxweighti;
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi(a,b)+=densi*degen*resinfoptr->q[a]*resinfoptr->q[b];
				}
			}
		}
		else{
			density[ires]=0.0;
		}
		netchi+=resinfoptr->netchi*density[ires];
		netchi0+=resinfoptr->netchi0*density[ires];
	}
}

void CResList::CalcEoSandChiandQdens(double T,double &P,double &epsilon,
double &nh,vector<double> &density,vector<double> &maxweight,Eigen::Matrix3d &chi,double &strangecontent,double &udcontent){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,m1=0.0,m2=0.0,degen,s;
	double width,minmass,maxweighti;
	double pi,epsiloni,densi,sigma2i,dedti;
	double Nud,Nstrange;
	bool special;
	strangecontent=udcontent=0.0;
	double netchi=0.0,netchi0=0.0;
	int a,b,n,ires,nres=resmap.size();
	chi.setZero();
	P=epsilon=s=nh=0.0;
	density.resize(nres);
	maxweight.resize(nres);
	
	for(ires=0;ires<nres;ires++){
		density[ires]=0.0;
	}
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		Nstrange=Nud=0.0;
		resinfoptr=rpos->second;
		ires=resinfoptr->ires;
		special=false;
		size_t found = resinfoptr->name.find("phi");
		if(found!=std::string::npos){ 
			Nstrange=2.0;
			special=true;
		}
		found = resinfoptr->name.find("eta");
		if(found!=std::string::npos){ 
			Nstrange=2.0/3.0;
			Nud=4.0/3.0;
			special=true;
		}
		if(!special){
			Nstrange=abs(resinfoptr->strange);
			if(resinfoptr->baryon!=0)
				Nud=3-Nstrange;
			else{
				if(resinfoptr->code!=22)
					Nud=2-Nstrange;
			}
		}		
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			if(USEPOLEMASS){
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
				maxweighti=-1.0;
			}
			else{
				width=resinfoptr->width;
				minmass=resinfoptr->minmass;
				if(resinfoptr->decay){
					m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
					m2=0.0;
					for(n=1;n<int(resinfoptr->branchlist[0]->resinfoptr.size());n++){
						m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
					}
				}
				if((minmass>0.0)&&(width>1.0E-3)){
					freegascalc_onespecies_finitewidth(m,m1,m2,T,width,
					epsiloni,pi,densi,sigma2i,dedti,maxweighti);
					if(densi!=densi){
						resinfoptr->Print();
						exit(1);
					}
				}
				else{
					freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
					maxweighti=-1.0;
				}
			}
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			nh+=densi*degen;
			strangecontent+=Nstrange*densi*degen;
			udcontent+=Nud*densi*degen;
			density[ires]=densi*degen;
			maxweight[ires]=maxweighti;
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi(a,b)+=densi*degen*resinfoptr->q[a]*resinfoptr->q[b];
				}
			}
		}
		else{
			density[ires]=0.0;
		}
		netchi+=resinfoptr->netchi*density[ires];
		netchi0+=resinfoptr->netchi0*density[ires];
	}
	strangecontent=strangecontent/s;
	udcontent=udcontent/s;
}

double CResList::triangle(double m0,double m1,double m2){
	double answer,m0sq,m1sq,m2sq;
	if(m0<m1+m2) {
		cout << "Disaster with triangle"<<endl;
		exit(1);
	}
	m0sq=m0*m0;m1sq=m1*m1;m2sq=m2*m2;
	answer=(m0sq*m0sq+m1sq*m1sq+m2sq*m2sq-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq))/(4.0*m0sq);
	return answer;
}

double CResList::GetLambda(double T,double P,double epsilon){
	int i,n;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp=0.0;
	//double Ipptest=0.0,Ptest=0.0,dIpptest=0.0,dp=4.0,p,e;
	double J,nfact,sign,alpha;
	double lambdafact;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			dIpp=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			/*
			dIpptest=0.0;
			for(p=0.5*dp;p<3000;p+=dp){
			e=sqrt(m*m+p*p);
			dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*( (2.0/3.0)*(p*p/e) - (2.0/15.0)*pow(p,4)/pow(e,3) );
			//dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*((2.0/15.0)*pow(p,4)/pow(e,3) );
			}
			Ipptest+=dIpptest;
			*/
			
			Ipp+=dIpp;
			//Ptest+=Ipptest;
		}
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);

	return lambdafact;
}

void CResList::freegascalc_onespecies(double mass,double T,double &epsiloni,double &Pi,double &densi,double &sigma2i,double &dedti,double &Ji){
	/*
	Ji=\frac{1}{2\pi^2} \int d^3p~e^{-E/T} \frac{p^2}{3*T*E^2}
	*/
	int n;
	double eta=mass/T;
	double J[100]={0.0};
	double expeta=exp(-eta),kappa;
	J[1]=-gsl_sf_expint_Ei(-eta)/expeta;
	for(n=2;n<100;n++){
		J[n]=(-eta*J[n-1]+1.0)/double(n-1);
	}
	Ji=1.0/eta;
	kappa=1.0;
	for(n=1;n<50;n++){
		kappa*=-(1.5-n)/double(n);
		Ji+=J[2*n]*kappa;
	}
	Ji*=expeta;
	freegascalc_onespecies(mass,T,epsiloni,Pi,densi,sigma2i,dedti);
	double Jcon=1.0/(2.0*PI*PI*pow(HBARC,3));
	Ji=densi-mass*mass*mass*Ji*Jcon;
	Ji=Ji/(3.0*T);
	
	double p,E,Jtest=0.0,dp=1.0;
	for(p=0.5*dp;p<2000;p+=dp){
		E=sqrt(p*p+mass*mass);
		Jtest+=Jcon*exp(-E/T)*dp*p*p*p*p/(3.0*E*E*T);
	}
	
}


void CResList::FindFinalProducts(double taumax){
	// all decay products -- decaying stops if tau_decay > taumax
	CResInfoMap::iterator rpos;
	CMassMap::iterator mpos;
	CResInfo *resinfo;
	massmap.clear();
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfo=rpos->second;
		massmap.insert(pair<double,CResInfo*>(resinfo->mass,resinfo));
	}
	for(mpos=massmap.begin();mpos!=massmap.end();mpos++){
		resinfo=mpos->second;
		resinfo->FindFinalProducts(taumax);
	}
	finalproductsfound=true;
}

double CResList::CalcBalanceNorm(int pid,int pidprime,double taumax){
	// ideal norm of B_hh'
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	CMassMap::iterator mpos;
	Eigen::Vector3d rho,rhoprime;
	Eigen::Matrix3d chitest,unity;
	double dens,densprime,weight,norm;
	int a,ires;

	chiinvf=chif.inverse();
	for(a=0;a<3;a++){
		rho(a)=rhoprime(a)=0.0;
	}
	
	norm=dens=densprime=0.0;
	double netq[3]={0.0};
	int iq;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfo=rpos->second;
		ires=resinfo->ires;
		for(iq=0;iq<3;iq++)
			netq[iq]+=densityf[ires]*resinfo->q[iq];
		
		if(resinfo->FindContent(pid,1.0,taumax,weight)){
			dens+=weight*densityf[ires];
			for(a=0;a<3;a++)
				rho(a)+=weight*densityf[ires]*resinfo->q[a];
		}
		if(resinfo->FindContent(pidprime,1.0,taumax,weight)){
			densprime+=weight*densityf[ires];
			for(a=0;a<3;a++)
				rhoprime(a)+=weight*densityf[ires]*resinfo->q[a];
		}
		if(resinfo->FindContent(-pid,1.0,taumax,weight)){
			dens+=weight*densityf[ires];
			for(a=0;a<3;a++)
				rho(a)-=weight*densityf[ires]*resinfo->q[a];
		}
		if(resinfo->FindContent(-pidprime,1.0,taumax,weight)){
			densprime+=weight*densityf[ires];
			for(a=0;a<3;a++)
				rhoprime(a)-=weight*densityf[ires]*resinfo->q[a];
		}
		
		if(resinfo->FindContentPairs(pid,pidprime,1.0,taumax,weight)){
			norm-=weight*densityf[ires];
		}
		
	}
	norm+=double((rho.transpose())*(chiinvf*rhoprime));
	norm=norm/densprime;
	
	return norm;
}

void CBranchInfo::Copy(CBranchInfo *oldbranch){
	resinfoptr=oldbranch->resinfoptr; //pointers for resinfo
	branching=oldbranch->branching;
}
