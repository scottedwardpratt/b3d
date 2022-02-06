#ifndef __ACCEPTANCE_H__
#define __ACCEPTANCE_H__

#include "commondefs.h"

class CAcceptance{
public:
	CAcceptance();
	CAcceptance(CparameterMap *parmapin);
	double ETAMIN,ETAMAX,PTMIN,PTMAX;
	int CENTRALITY;
	CparameterMap *parmap;
	virtual void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	virtual void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	virtual double GetDelYMax(int pida,int pidb);
};

class CAcceptance_CHEAP : public CAcceptance{
public:
	CAcceptance_CHEAP(CparameterMap *parmapin);
	double ptmin,ptmax;
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	double GetDelYMax(int pida,int pidb){
		if(abs(pida)==0 && abs(pidb)==0){
			return 0.0;
		}
		else
			return 2.0;
	}
};

class CAcceptance_ALICE : public CAcceptance{
public:
	CAcceptance_ALICE(CparameterMap *parmapin);
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);  // ignores cuts in rapidity
	void CalcAcceptance_Realistic(bool &accept,double &efficiency,CPart *part); // cuts in rapidity
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	double GetDelYMax(int pida,int pidb);
};

class CAcceptance_ALICE_Perfect : public CAcceptance{
public:
	CAcceptance_ALICE_Perfect(CparameterMap *parmapin);
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);  // ignores cuts in rapidity
	void CalcAcceptance_Realistic(bool &accept,double &efficiency,CPart *part); // cuts in rapidity
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	double GetDelYMax(int pida,int pidb);
};

class CAcceptance_STAR : public CAcceptance{
public:
	CAcceptance_STAR(CparameterMap *parmapin);
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	void star_acc_eff(int pid,double pt,double eta,bool &accept,double &eff);
	
	// For efficiencies of non-Identified parts
  //double ManuelEff(double eta, double pt);
	double ScottEff(double eta, double pt);
	double GetDelYMax(int pida,int pidb){
		if(abs(pida)==0 && abs(pidb)==0){
			return 0.0;
		}
		else
			return 2.0;
	}
	
};


#endif
