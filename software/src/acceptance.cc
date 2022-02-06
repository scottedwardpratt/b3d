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
	printf("hmmmmmm, should not be here in dummy routine for CalcAcceptance\n");
	if(part!=NULL){
		accept=true;
		efficiency=1.0;
	}
	else
		printf("Calculating acceptance for null part\n");
}

void CAcceptance::CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part){
	printf("hmmmmmm, should not be here in dummy routing for CalcAcceptance\n");
  accept=true;
	efficiency=1.0;
	if(part!=NULL){
		accept=true;
		efficiency=1.0;
	}
	else
		printf("Calculating acceptance for null part\n");
}

double CAcceptance::GetDelYMax(int pida,int pidb){
	if(abs(pida)==0 && abs(pidb)==0){
		return 0.0;
	}
	else
		return 2.0;
}
