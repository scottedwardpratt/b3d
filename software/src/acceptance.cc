//
//  acceptance.cc
//  

#include "acceptance.h"

using namespace std;

CAcceptance::CAcceptance(){
	ETAMIN=-1.0;
	ETAMAX=1.0;
	PTMIN=0.0;
	PTMAX=100000000.0;
	CENTRALITY=0.0;
}

CAcceptance::CAcceptance(CparameterMap *parmapin){
	parmap=parmapin;
	ETAMIN=-1.0;
	ETAMAX=1.0;
	PTMIN=0.0;
	PTMAX=100000000.0;
	CENTRALITY=0.0;
}

void CAcceptance::CalcAcceptance(bool &accept,double &efficiency,CPart *part){
	sprintf(message,"hmmmmmm, should not be here in dummy routing for CalcAcceptance\n");
	CLog::Fatal(message);
	if(part==NULL){
		accept=true; efficiency=1.0; 
		accept=false;
	}
}

void CAcceptance::CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part){
	sprintf(message,"hmmmmmm, should not be here in dummy routing for CalcAcceptanceNoID\n");
	CLog::Fatal(message);
	if(part==NULL){
		accept=true; efficiency=1.0; 
		accept=false;
	}
}

double CAcceptance::GetDelYMax(int pida,int pidb){
	if(abs(pida)==0 && abs(pidb)==0){
		return 0.0;
	}
	else
		return 2.0;
}
