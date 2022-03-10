#include "action.h"
#include "part.h"

void CAction::PerformCollide(){
	int colltype,iproduct,nproducts;
	double sigmatot;
	CPart *part1,*part2,*part;
	CPartMap::iterator ppos;
	CB3DCell *cell;
	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	b3d->GetDeadParts(product);

	sigmatot=sigma_scatter+sigma_merge+sigma_annihilation+sigma_inel;
	double r=b3d->randy->ran();
	if(r<sigma_scatter/sigmatot){
		colltype=b3d->Collide_Scatter(part1,part2,nproducts,product);
	}
	else if(r<(sigma_scatter+sigma_merge)/sigmatot){
		colltype=b3d->Collide_Merge(part1,part2,sigma_merge,dsigma_merge,nproducts,product);
	}
	else if(r<(sigma_scatter+sigma_merge+sigma_annihilation)/sigmatot){
		colltype=b3d->Collide_Annihilate(part1,part2,nproducts,product);
	}
	else{
		if(!b3d->INELASTIC){
			sprintf(message,"In PeformCollide: Trying to Perform Inelastic Collision???");
			CLog::Fatal(message);
		}
		colltype=b3d->Collide_Inelastic(part1,part2,nproducts,product);
	}

	//colltype=b3d->Collide(part1,part2,nproducts,product,pibsquared);

	if((part1->balanceID<0 && part2->balanceID>=0) || (part1->balanceID>=0 && part2->balanceID<0)){
		colltype=-2;
	}

	if(colltype==0 || nproducts==0){
		b3d->npass+=1;
		part1->actionmother=b3d->nactions;
		part2->actionmother=b3d->nactions;
	}
	else if(colltype==-2){
		if(part1->balanceID>=0 && part2->balanceID<0){
			part1->CopyMomentumInfo(product[0]);
			part1->SetMass();
		}
		if(part2->balanceID>=0 && part1->balanceID<0){
			part2->CopyMomentumInfo(product[1]);
			part2->SetMass();
		}
		part1->tau_lastint=tau;
		if(part1->currentmap!=&(b3d->PartMap))
			part1->ChangeMap(&(b3d->PartMap));
		part1->active=true;
		part1->actionmother=b3d->nactions;
		part1->FindActions();
		part2->tau_lastint=tau;
		if(part2->currentmap!=&(b3d->PartMap))
			part2->ChangeMap(&(b3d->PartMap));
		part2->active=true;
		part2->actionmother=b3d->nactions;
		part2->FindActions();
	}
	else{
		part1->Kill();
		part2->Kill();
		for(iproduct=0;iproduct<nproducts;iproduct++){
			part=product[iproduct];
			part->SetMass();
			part->active=true;
			part->weight=1.0;
			part->tau_lastint=tau;
			part->actionmother=b3d->nactions;
			cell=part->FindCell();
			part->ChangeCell(cell);
			part->ChangeMap(&b3d->PartMap);
			part->FindActions();
		}
		for(iproduct=nproducts;iproduct<5;iproduct++){
			part=product[iproduct];
			part->Kill();
		}
	}
}
