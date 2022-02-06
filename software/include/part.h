#ifndef __PART_H__
#define __PART_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <map>
#include <unordered_map>
#include <sys/stat.h>
#include "commondefs.h"

using namespace std;
//using namespace boost;

class CPart{
public:
	CPart();
	CPart(int keyset);
	~CPart();
	CB3DCell *cell,*nextcell;
	double tau0,tau_lastint,tauexit,taudecay;
	double y,eta,msquared,weight;
	int badmother;
	FourVector p,r;
	double eta0,phi0;
	int listid,nscatt,balanceID,key;
	double bweight;
	int actionmother; //refers to action from which particle was created
	CResInfo *resinfo;
	bool active;
	CPartMap *currentmap; // PartList for a Cell, or b3d->DeadPartList
	
	void InitBalance(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,double bweight,int balanceid);
	void Propagate(double tau);
	void FindDecay();
	void FindCellExit();
	void Init(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,double bweight);
	void Init_NoColls(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,double bweight);
	void Setp0();
	void SetMass();
	void Copy(CPart *part);
	void CopyPositionInfo(CPart *part);
	void CopyMomentumInfo(CPart *part);
	double GetMass();
	void CyclicReset();
	void SetInitialKey();
	void CheckRapidity();

	void Kill();
	void FindActions();
	void SubtractAction(CAction *actionptr);
	void AddAction(CAction *actionptr);
	void KillActions();
	void Print();
	void CheckMap(CPartMap *expectedpartmap);
	void CheckCell();
	void ChangeMap(CPartMap *newmap);
	void BjorkenTranslate();
	void BjorkenUnTranslate();
	void Boost(FourVector &u);
	void BoostP(FourVector &u);
	void BoostR(FourVector &u);
	void CalcDCA(double *dca);
	void BoostRap(double dely);
	//~CPart();

	// These are the actions involving these particles
	CActionMap actionmap;

	CB3DCell *FindCell();

	static CB3D *b3d;
	static CBalance *cb;
	double GetEta(double tau);
	double GetPseudoRapidity();
	double GetRapidity();
	double GetMT();
	void SetY();
	void SetEta(double neweta);
	void GetHBTPars(double &t,double &rout,double &rside,double &rlong);
	void FindCollisions();

	CPartMap::iterator GetPos(CPartMap *pmap);
	CPartMap::iterator DeleteFromMap(CPartMap *partmap);
	void ChangeCell(CB3DCell *newcell);
	void RemoveFromCell();

	void ChangePartMap(CPartMap *newmap);
	CPartMap::iterator DeleteFromCurrentMap();
	void AddToMap(CPartMap *newmap);
	void AddToMap(CPartMap::iterator guess,CPartMap *newmap);
};

class CB3DBinaryPartInfo{
public:
	int ID;
	double tau,x,y,eta;
	double px,py,rapidity,bweight,weight;
	void Print(){
		printf("ID=%d, tau=%g, x=%g, y=%g, eta=%g\n",ID,tau,x,y,eta);
		printf("px=%g,py=%g, rapidity=%g, bweight=%g\n",px,py,rapidity,weight);
	}
};

class CB3DBinaryBalancePartInfo{
public:
	int ID,balanceID;
	double px,py,rapidity,bweight;
	void Print(){
		printf("ID=%d, balanceID=%d\n",ID,balanceID);
		printf("px=%g,py=%g, rapidity=%g,bweight=%g\n",px,py,rapidity,bweight);
	}
};

#endif
