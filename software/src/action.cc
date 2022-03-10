#include  "b3d.h"
#include "action.h"

CB3D *CAction::b3d=NULL;
char *CAction::message=new char[500];

CAction::CAction(){
}

CAction::CAction(int keyset){
	key=keyset;
	listid=key;
	tau=0.0;
	type=-1;
	currentmap=&b3d->DeadActionMap;
	b3d->DeadActionMap.insert(CActionPair(key,this));
	partmap.clear();
	b3d->nactionstot+=1;
}

// type=0(creation) 1(decay) 2(collision) 3(VizWrite) 4(DensCalc)

CActionMap::iterator CAction::GetPos(CActionMap *emap){
	pair<CActionMap::iterator,CActionMap::iterator> epospair;
	CActionMap::iterator epos;
	epospair=emap->equal_range(key);
	epos=epospair.first;
	while(epos->second!=this && epos!=epospair.second){
		++epos;
	}
	if(epos->second!=this){
		sprintf(message,"CAction::GetPos cannot find this action\n");
		CLog::Fatal(message);
		return emap->end();
	}
	else
		return epos;
}

void CAction::MoveToActionMap(){
	CActionMap::iterator epos,eepos;
	if(currentmap==&b3d->ActionMap){
		sprintf(message,"trying to move action to ActionMap even though action is already in ActionMap\n");
		sprintf(message,"%sWrong current map\n",message);
		CLog::Fatal(message);
	}
	if(currentmap==&b3d->DeadActionMap){
		epos=GetPos(currentmap);
		if(epos!=currentmap->end()){
			b3d->DeadActionMap.erase(epos);
			partmap.clear();
		}
		else{
			sprintf(message,"cannot find epos for action in DeadActionMap!!!\n");
			sprintf(message,"%sDeadActionMap.size=%d\n",message,int(b3d->DeadActionMap.size()));
			CLog::Fatal(message);
		}
		key=tau;
		AddToMap(&b3d->ActionMap);
	}
}

void CAction::Kill(){
	CPart *part;
	CPartMap::iterator ppos;
	if(currentmap==&b3d->ActionMap){
		CActionMap::iterator eepos,epos=GetPos(currentmap);
		if(epos==currentmap->end()){
			sprintf(message,"in CAction::Kill(), not in map\n");
			sprintf(message,"%sb3d->ActionMap.size()=%d\n",message,int(b3d->ActionMap.size()));
			CLog::Info(message);
			ppos=partmap.begin();
			part=ppos->second;
			part->Print();
			epos=GetPos(&b3d->DeadActionMap);
			if(epos!=b3d->DeadActionMap.end()){
				sprintf(message,"found action in DeadActionMap\n");
			}
			sprintf(message,"%snot in DeadActionMap either\n",message);
			CLog::Fatal(message);
		}
		currentmap->erase(epos);
		b3d->nactionkills+=1;
		key=0;
		AddToMap(b3d->DeadActionMap.end(),&b3d->DeadActionMap);
		//AddToMap(&b3d->DeadActionMap);
		ppos=partmap.begin();
		while(ppos!=partmap.end()){
			part=ppos->second;
			epos=part->actionmap.begin();
			while(epos!=part->actionmap.end()){
				eepos=epos;
				++eepos;
				if(epos->second==this)
					part->actionmap.erase(epos);
				epos=eepos;
			}
			++ppos;
		}
	}
}

void CAction::AddToMap(CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(CActionPair(key,this));
}

void CAction::AddToMap(CActionMap::iterator guess,CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(guess,CActionPair(key,this));
}

void CAction::AddPart(CPart *part){
	partmap.insert(CPartPair(part->key,part));
}

void CAction::Print(){
	sprintf(message,"___________ type=%d, tau=%g, nparts=%d ___________\n",type,tau,int(partmap.size()));
	CLog::Info(message);
	CPartMap::iterator ppos;
	CPart *part;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Print();
	}
	sprintf(message,"_________________________________________________\n");
	CLog::Info(message);
}

void CAction::CheckPartList(){
	CPart *part;
	CPartMap::iterator ppos,ppos2;
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		part=ppos->second;
		ppos2=part->GetPos(&(b3d->PartMap));
		if(ppos2==b3d->PartMap.end()){
			part->Print();
			sprintf(message,"%s____________ CAction::CheckPartList FATAL, action type=%d ________________\n",message,type);
			sprintf(message,"%siterator not in expected pmap\n",message);
			CLog::Fatal(message);
		}
		++ppos;
	}
}

void CAction::PerformDensCalc(){
	int itau;
	itau=lrint(floor((tau-0.001)/b3d->DENSWRITE_DELTAU));
	if(itau>=b3d->DENSWRITE_NTAU){
		sprintf(message,"trying to perform CAction::DensCalc() for itau>=DENSWRITE_NTAU, =%d",itau);
		CLog::Fatal(message);
	}
	CB3DCell *cell;
	int ix,iy,ieta;
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			for(ieta=0;ieta<2*b3d->NETA;ieta++){
				cell=b3d->cell[ix][iy][ieta];
				cell->dens[itau]+=cell->partmap.size();
			}
		}
	}
}

