#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__

#include "b3d.h"
#include "part.h"

void CAction::PerformCollide(){
	int colltype,iproduct,nproducts;
	CPart *part1,*part2,*part;
	CPartMap::iterator ppos;
	CB3DCell *cell;
	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	part1->actionmother=b3d->nactions;
	part2->actionmother=b3d->nactions;
	b3d->GetDeadParts(product);
	colltype=b3d->Collide(part1,part2,nproducts,product,pibsquared);
	if((part1->balanceID<0 && part2->balanceID>=0) || (part1->balanceID>=0 && part2->balanceID<0)){
		colltype=-2;
	}
	if(colltype==0 || nproducts==0){
		b3d->npass+=1;
		return;
	}
	if(colltype==1)
		b3d->nmerge+=1;
	if(colltype==2)
		b3d->nscatter+=1;
	if(colltype==-2)
		b3d->nbscatter+=1;
	if(colltype==3)
		b3d->ninelastic+=1;
	if(colltype==4)
		b3d->nannihilate+=1;
	
	if(colltype==-2){
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
	}

	if(colltype!=-2){
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
	
	//printf("check colltype=%d, partmap size=%d\n",colltype,int(b3d->PartMap.size()));
}

#endif