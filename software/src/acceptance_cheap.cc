//
//  acceptance_CHEAP.cc
//
#include "acceptance.h"
#include "part.h"
#include "resonances.h"
#include "parametermap.h"

using namespace std;

CAcceptance_CHEAP::CAcceptance_CHEAP(CparameterMap *parmapin) : CAcceptance(){
	ETAMAX=parmapin->getD("B3D_ETAMAX",6.0);
	ETAMIN=-ETAMAX;
	ptmin=200.0;
	ptmax=20000000.0;
}

void CAcceptance_CHEAP::CalcAcceptance(bool &accept,double &efficiency,CPart *part){
	double pt;
	double eta,pmag;
	
	accept=false;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	//printf("pt=%g\n",pt);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	/*
	int pid=part->resinfo->code;
	double y,gammav,m,A0,ctau_kaon=3.7,ctau_pion=7.8,lmin=1.0;
	if(dca[0]<1.5){
	efficiency=0.0;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	//printf("pt=%g\n",pt);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	//y=atanh(part->p[3]/part->p[0]);
	m=part->resinfo->mass;
	 
	if(eta>ETAMIN && eta<ETAMAX && pt>ptmin && pt<ptmax){
	if(pid==211 || pid==-211){
	A0=0.759;
	accept=true;
	m=sqrt(part->p[0]*part->p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_pion));
	}
	if(pid==2212 || pid==-2212){
	if(pid==2212) A0=0.221;
	else A0=0.246;
	accept=true;
	efficiency=A0;
	}
	if(pid==321 || pid==-321){
	A0=0.45;
	accept=true;
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	m=sqrt(part->p[0]*part->p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_kaon));
	}
	if(pt>600.0) efficiency*=0.65;
	}
	}
	else{
	printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	}
	if(m!=m){
	part->Print();
	exit(1);
	}
	 
	*/
	
	/*
	accept=false;
	efficiency=0.0;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	//y=atanh(part->p[3]/part->p[0]);
	//if(dca[0]>0.000001) printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	if(pt>ptmin && pt<ptmax && eta>ETAMIN && eta<ETAMAX){
	accept=true;
	efficiency=0.8;
	}
	**/
	if(pt<10000000 && fabs(eta)<ETAMAX){
		accept=true;
		efficiency=1.0;
	}
	else{
		accept=false;
		efficiency=0.0;
	}
}

void CAcceptance_CHEAP::CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part){
	double pt,eta,pmag;
	accept=false;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	/* accept=false;
	int pid=part->resinfo->code;
	double gammav,m;
	double ctau_kaon=3.7,ctau_pion=7.8,lmin=1.0;
	double A0;
	if(dca[0]<1.5){
	efficiency=0.0;
	//printf("pt=%g\n",pt);
	//y=atanh(part->p[3]/part->p[0]);
	m=part->resinfo->mass;
	 
	if(eta>ETAMIN && eta<ETAMAX && pt>ptmin && pt<ptmax){
	if(pid==211 || pid==-211){
	A0=0.759;
	accept=true;
	m=sqrt(part->p[0]*part->p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_pion));
	}
	if(pid==2212 || pid==-2212){
	if(pid==2212) A0=0.221;
	else A0=0.246;
	accept=true;
	efficiency=A0;
	}
	if(pid==321 || pid==-321){
	A0=0.45;
	accept=true;
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	m=sqrt(part->p[0]*part->p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_kaon));
	}
	if(pt>600.0) efficiency*=0.65;
	}
	}
	else{
	printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	}
	if(m!=m){
	part->Print();
	exit(1);
	}
	 
	*/
	
	/*
	accept=false;
	efficiency=0.0;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	//y=atanh(part->p[3]/part->p[0]);
	//if(dca[0]>0.000001) printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	if(pt>ptmin && pt<ptmax && eta>ETAMIN && eta<ETAMAX){
	accept=true;
	efficiency=0.8;
	}
	**/
	
	if(pt<10000000 && fabs(eta)<ETAMAX){
		accept=true;
		efficiency=1.0;
	}
	else{
		accept=false;
		efficiency=0.0;
	}
	
	
}
