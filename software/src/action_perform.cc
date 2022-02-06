#ifndef __ACTION_PERFORM_CC__
#define __ACTION_PERFORM_CC__

#include "b3d.h"
#include "part.h"

void CAction::Perform(){
	CPartMap::iterator ppos;
	CPart *part;
	//printf("performing action of type %d, tau=%g\n",type,tau);

	b3d->nactions+=1;
	Kill();
	b3d->tau=tau;

	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Propagate(tau);
		if(part->currentmap != &(b3d->PartMap)){
			printf("wrong map for part in Action List\n");
			part->Print();
			part->active=true;
			part->ChangeMap(&(b3d->PartMap));
		}
	}
	if(tau+1.0E-4<b3d->tau){
		printf("FATAL:: action earlier than tau!!!!, b3d->tau=%15.10e, action tau=%15.10e\n",b3d->tau,tau);
		exit(1);
	}
	if(type==6){
		PerformExitCell();
	}
	else if(type==0){
		PerformActivate();
	}
	else if(type==1){
		PerformDecay();
	}
	else if(type==2){
		PerformCollide();
	}
	else if(type==4){
		PerformDensCalc();
	}
	else if(type==5){
		PerformMuTCalc();
	}
	else if(type==7){
		PerformSECalc();
	}
	else{
		printf("FATAL: action type = %d is unknown, exiting\n",type);
		exit(1);
	}
	//b3d->CheckPartMap();
	//printf("---------\n");
	
}

#endif
