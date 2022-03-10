#include "action.h"
#include "part.h"
#include "cell.h"
#include "resonances.h"

void CAction::PerformDecay(){
	CPart *mother,*dptr;
	CPartMap::iterator ppos;
	int ibody,nbodies;
	double mtot,mt,etamax=b3d->ETAMAX,mothermass;
	double deleta;
	ppos=partmap.begin();
	mother=ppos->second;
	b3d->GetDeadParts(product);

	mothermass=mother->GetMass();
	if(mother->cell!=NULL && mother->cell!=mother->FindCell() && tau<b3d->TAUCOLLMAX){
		mother->CheckRapidity();
		mother->cell->Print();
		mother->FindCell()->Print();
		mother->Print();
		sprintf(message,"Cells don't match for decaying mother\n");
		CLog::Fatal(message);
	}
	if(tau>b3d->TAUCOLLMAX || mother->cell==NULL){
		while(b3d->BJORKEN && (mother->eta<-etamax || mother->eta>etamax)){
			if(mother->eta<-etamax) deleta=2.0*etamax*ceil((-etamax-mother->eta)/(2.0*etamax));
			if(mother->eta>etamax) deleta=-2.0*etamax*ceil((mother->eta-etamax)/(2.0*etamax));
			mother->eta+=deleta;
			mother->y+=deleta;
			mt=mother->GetMT();
			mother->p[0]=mt*cosh(mother->y);
			mother->p[3]=mt*sinh(mother->y);
			mother->r[0]=tau*cosh(mother->eta);
			mother->r[3]=tau*sinh(mother->eta);
		}
	}
	int ntry=0;
	do{
		mtot=0.0;
		if(ntry<25)
			mother->resinfo->DecayGetResInfoPtr(nbodies,daughterresinfo);
		else{
			mother->resinfo->DecayGetResInfoPtr_minmass(nbodies,daughterresinfo);
		}
		for(ibody=0;ibody<nbodies;ibody++){
			mtot+=daughterresinfo[ibody]->mass;
		}
		if(ntry>25){
			mother->Print();
			sprintf(message,"FATAL: action_perform_decay, ntry too big, mothermass=%g\n",mother->GetMass());
			CLog::Fatal(message);
		}
		ntry++;
	}while(mtot>mothermass);
	for(ibody=0;ibody<nbodies;ibody++){
		product[ibody]->resinfo=daughterresinfo[ibody];
	}
	
	b3d->Decay(mother,nbodies,product);
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for(ibody=0;ibody<nbodies;ibody++){
		dptr=product[ibody];
		dptr->active=true;
		dptr->bweight=mother->bweight;
		dptr->balanceID=mother->balanceID;
		dptr->nscatt=0;
		dptr->tau_lastint=tau;
		dptr->actionmother=b3d->nactions;
		dptr->ChangeCell(dptr->FindCell());
		if(dptr->currentmap!=&(b3d->PartMap))
			dptr->ChangeMap(&(b3d->PartMap));
		dptr->FindActions();
	}
	mother->Kill();
	b3d->ndecay+=1;
}
