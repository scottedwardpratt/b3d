#include "action.h"
#include "part.h"

void CAction::Perform(){
	CPartMap::iterator ppos;
	CPart *part;
	CB3DCell *cell;

	//sprintf(message,"Performing Action of type %d\n",type);
	//CLog::Info(message);

	b3d->nactions+=1;
	Kill();
	b3d->tau=tau;

	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Propagate(tau);
		if(part->currentmap != &(b3d->PartMap)){
			sprintf(message,"wrong map for part in Action List\n");
			CLog::Info(message);
			part->Print();
			part->active=true;
			part->ChangeMap(&(b3d->PartMap));
			cell=part->FindCell();
			if(part->cell!=cell){
				part->ChangeCell(cell);
			}
		}
	}
	if(tau+1.0E-4<b3d->tau){
		sprintf(message,"FATAL:: action earlier than tau!!!!, b3d->tau=%15.10e, action tau=%15.10e\n",b3d->tau,tau);
		CLog::Fatal(message);
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
		PerformMuTCalcUpdateNPE();
	}
	else if(type==7){
		PerformSECalc();
	}
	else{
		sprintf(message,"FATAL: action type = %d is unknown, exiting\n",type);
		CLog::Fatal(message);
	}

	
}
