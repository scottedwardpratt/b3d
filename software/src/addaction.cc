#include "b3d.h"
#include "part.h"
#include "cell.h"

void CB3D::AddAction_Activate(CPart *part){
	CActionMap::iterator epos;
	part->active=false;
	CAction *action;
	if(BJORKEN && fabs(part->eta)>ETAMAX){
		sprintf(message,"CB3D::AddAction_Activate, eta out of bounds, =%g\n",fabs(part->eta));
		CLog::Fatal(message);
	}
	action=GetDeadAction();
	if(action->currentmap==&ActionMap){
		sprintf(message,"don't even try, action wasn't dead, key=%d\n",int(action->key));
		CLog::Fatal(message);
	}
	action->type=0;
	action->tau=part->tau0;
	action->MoveToActionMap();
	action->partmap.insert(CPartPair(part->key,part));
	if(action->tau<tau){
		sprintf(message,"trying to AddAction_Activate at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
		CLog::Fatal(message);
	}
	part->actionmap.insert(CActionPair(action->key,action));
}

void CB3D::AddAction_Decay(CPart *part,double taudecay){
	CAction *action=GetDeadAction();
	action->tau=taudecay;
	action->partmap.clear();
	action->type=1;
	action->MoveToActionMap();
	action->partmap.insert(CPartPair(part->key,part));
	part->actionmap.insert(CActionPair(action->key,action));
	if(action->tau<tau){
		part->Print();
		part->cell->Print();
		sprintf(message,"CB3D::AddAction_Decay, trying to AddAction_Decay at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
		CLog::Fatal(message);
	}
}

void CB3D::AddAction_ExitCell(CPart *part){
	CAction *action;
	action=GetDeadAction();
	if(part->tauexit<TAUCOLLMAX){
		CActionMap::iterator epos=DeadActionMap.begin();
		action=epos->second;
		action->type=6;
		action->tau=part->tauexit;
		action->MoveToActionMap();
		action->partmap.insert(CPartPair(part->key,part));
		part->actionmap.insert(CActionPair(action->key,action));
		if(action->tau<tau-1.0E-10){
			part->Print();
			part->cell->Print();
			sprintf(message,"CB3D::AddAction_ExitCell, trying to AddAction_ExitCell at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
			CLog::Fatal(message);
		}
	}
}

void CB3D::AddAction_Collision(CPart *part1,CPart *part2,double taucoll,double pibsquared){
	CAction *action=GetDeadAction();
	action->type=2;
	action->tau=taucoll;
	action->pibsquared=pibsquared;
	action->MoveToActionMap();
	if(action->tau<tau){
		action->Print();
		sprintf(message,"trying to AddAction_Collision at earler time!!!  tau=%g\n",tau);
		CLog::Fatal(message);
	}
	action->partmap.insert(CPartPair(part1->key,part1));
	action->partmap.insert(CPartPair(part2->key,part2));

	part1->actionmap.insert(CActionPair(action->key,action));
	part2->actionmap.insert(CActionPair(action->key,action));
}

void CB3D::AddAction_DensCalc(double taucalc){
	CAction *action;
	action=GetDeadAction();
	action->type=4;
	action->tau=taucalc;
	action->MoveToActionMap();
	action->partmap.clear(); 
	if(action->tau<tau){
		action->Print();
		sprintf(message,"trying to AddAction_DensCalc at earler time!!!  tau=%g\n",tau);
		CLog::Fatal(message);
	}
}

void CB3D::AddAction_MuTCalc_UpdateNPE(double taucalc){
	CAction *action;
	action=GetDeadAction();
	action->type=5;
	action->tau=taucalc;
	action->MoveToActionMap();
	action->partmap.clear(); 
	if(action->tau<tau){
		action->Print();
		sprintf(message,"trying to AddAction_MuTCalc at earler time!!!  tau=%g\n",tau);
		CLog::Fatal(message);
	}
}

void CB3D::AddAction_SECalc(double taucalc){
	CAction *action;
	action=GetDeadAction();
	action->type=7;
	action->tau=taucalc;
	action->MoveToActionMap();
	action->partmap.clear(); 
	if(action->tau<tau){
		action->Print();
		sprintf(message,"trying to AddAction_MuTCalc at earler time!!!  tau=%g\n",tau);
		CLog::Fatal(message);
	}
}
