#include "b3d.h"
#include "part.h"
#include "misc.h"
#include "resonances.h"
#include "cell.h"
#include "constants.h"
#include "randy.h"
#include "action.h"

void CB3D::WriteAnnihilationData(){
	if(BARYON_ANNIHILATION){
		double total=0.0;
		int iitau,imax=lrint(TAUCOLLMAX);
		for(iitau=0;iitau<imax;iitau++){
			sprintf(message,"%6.2f %g\n",(iitau+0.5),annihilation_array[iitau]);
			CLog::Info(message);
			total+=annihilation_array[iitau];
		}
		sprintf(message,"%g total annihilations, nbaryons=%d, annihilation fraction=%g\n",total,nbaryons,2.0*total/double(nbaryons));
		CLog::Info(message);
	}
}

void CB3D::PerformAllActions(){
	int iitau;
	if(DENSWRITE){
		for(iitau=0;iitau<DENSWRITE_NTAU;iitau++){
			AddAction_DensCalc((iitau+1.0)*DENSWRITE_DELTAU);
		}
	}
	CAction *action;
	nscatter=nbscatter=ndecay=npass=nmerge=nswallow=nexit=nactivate=0;
	ninelastic=ncheck=nactionkills=nbaryons=ncheck1=ncheck2=0;
	ncollisions=oldncollisions=nannihilate=nregenerate=0;
	tau=0.0;
	nactions=0;	
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Perform();
		epos=ActionMap.begin();
	}
}

void CB3D::KillAllActions(){
	CAction *action;
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Kill();
		epos=ActionMap.begin();
	}
	nactions=0;
}

void CB3D::KillAllParts(){
	CPart *part;
	CPartMap::iterator ppos;
	
	ppos=PartMap.begin();
	while(ppos!=PartMap.end()){
		part=ppos->second;
		if(part->currentmap!=&PartMap){
			sprintf(message,"Fatal: KillAllParts:  currentpart not listed as PartMap\n");
			CLog::Info(message);
			part->Print();
			sprintf(message,"PartMap.size=%d, DeadPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()));
			part->currentmap=&PartMap;
			CLog::Fatal(message);
			//Misc::Pause();
		}
		part->Kill();
		ppos=PartMap.begin();
	}
	
	// Recheck 
	for(ppos=DeadPartMap.begin();ppos!=DeadPartMap.end();++ppos){
		part=ppos->second;
		if(part->currentmap!=&DeadPartMap){
			sprintf(message,"particle in dead part map has wrong current map\n");
			CLog::Info(message);
			part->Print();
			sprintf(message,"PartMap.size=%d, DeadPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()));
			CLog::Fatal(message);
		}
	}

	int ix,iy,ieta;
	if(COLLISIONS){
		CPartMap *partmap;
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					partmap=&(cell[ix][iy][ieta]->partmap);
					ppos=partmap->begin();
					while(ppos!=partmap->end()){
						part=ppos->second;
						part->currentmap=&DeadPartMap;
						partmap->erase(ppos);
						//part->Kill();
						ppos=partmap->begin();
					}
				}
			}
		}
	}
}

void CB3D::PrintActionMap(CActionMap *actionmap){
	CActionMap::iterator epos;
	CAction *action;
	int iaction=0;
	sprintf(message,"_________________ ACTIONMAP %d actions _________________________\n",int(actionmap->size()));
	CLog::Info(message);
	for(epos=actionmap->begin();epos!=actionmap->end();++epos){
		iaction+=1;
		action=epos->second;
		action->Print();
	}
}

/*
void CB3D::FindAllCollisions(){
	double taucoll,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	CActionMap::iterator epos;
	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
	}
	ppos1=PartMap.begin();
	part1=ppos1->second;
	ppos2=ppos1; ++ppos2;
	while(ppos2!=PartMap.end()){
		part1=ppos1->second; part2=ppos2->second;
		if(part1->balanceID<0 || part2->balanceID<0)
			FindCollision(part1,part2,taucoll);
		ppos1=ppos2;
		++ppos2;
	}
}
*/

void CB3D::PrintPartList(){
	CPartMap::iterator ppos2,ppos1=PartMap.begin();
	while(ppos1!=PartMap.end()){
		sprintf(message,"%d ",ppos1->second->listid);
		ppos2=ppos1; ++ppos2;
		if(ppos2!=PartMap.end()){
			if(ppos1->second->actionmother!=ppos2->second->actionmother)
				sprintf(message,"%s| ",message);
		}
		++ppos1;
	}
	sprintf(message,"%s\n",message);
	CLog::Info(message);
}

void CB3D::ListFutureCollisions(){
	CActionMap::iterator epos=ActionMap.begin();
	CAction *action;
	CPartMap::iterator p1,p2;
	sprintf(message,"------------------- LIST OF FUTURE COLLISIONS ---------------------\n");
	while(epos!=ActionMap.end()){
		action=epos->second;
		if(action->type==2){
			p1=action->partmap.begin();
			p2=p1; ++p2;
			sprintf(message,"%s%d  %d  will collide at %g\n",message,p1->second->listid,p2->second->listid,double(action->tau));
		}
		epos++;
	}
	CLog::Info(message);
}

// Note part2 is fake and part1 is real
void CB3D::SplitPart(CPart *part1,CPart *part2){
	double oldeta,mt,g1,g2;
	CB3DCell *ccell;
	part2->Copy(part1); // does not change reality or weights
	if(BJORKEN){
		oldeta=part1->eta;
		part1->eta=-ETAMAX+2.0*ETAMAX*randy->ran();
		part1->y+=(part1->eta-oldeta);
		part1->r[3]=part1->tau0*sinh(part1->eta);
		part1->r[0]=sqrt(part1->tau0*part1->tau0+part1->r[3]*part1->r[3]);
		mt=sqrt(part1->msquared+part1->p[1]*part1->p[1]+part1->p[2]*part1->p[2]);
		part1->p[3]=mt*sinh(part1->y);
		part1->Setp0();
	}
	else{
		randy->ran_gauss2(g1,g2);
		part1->r[1]+=0.5*g1;
		part1->r[2]+=0.5*g2;
		oldeta=part1->eta;
		part1->eta+=0.5*randy->ran_gauss()/tau;
		part1->y+=(part1->eta-oldeta);
		part1->r[3]=part1->tau0*sinh(part1->eta);
		part1->r[0]=sqrt(part1->tau0*part1->tau0+part1->r[3]*part1->r[3]);
		mt=sqrt(part1->msquared+part1->p[1]*part1->p[1]+part1->p[2]*part1->p[2]);
		part1->p[3]=mt*sinh(part1->y);
		part1->Setp0();
	}
	
	ccell=part1->FindCell();
	part1->ChangeCell(ccell);
	if(part1->currentmap!=&PartMap)
		part1->ChangeMap(&PartMap);
	
	ccell=part2->FindCell();
	part2->ChangeCell(ccell);
	if(part2->currentmap!=&PartMap)
		part2->ChangeMap(&PartMap);
}

CPart* CB3D::GetDeadPart(){
	if(DeadPartMap.size()==0){
		for(int ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++){
			new CPart(npartstot);
		}
		sprintf(message,"made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
		CLog::Info(message);
	}
	return DeadPartMap.begin()->second;
}

void CB3D::GetDeadParts(CPart *&part1,CPart *&part2){
	int ipart;
	while(DeadPartMap.size()<2){
		for(ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++)
			new CPart(npartstot);
		sprintf(message,"made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
		CLog::Info(message);
	}
	CPartMap::iterator ppos=DeadPartMap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
}

void CB3D::GetDeadParts(array<CPart*,5> &product){
	int ipart;
	while(DeadPartMap.size()<5){
		for(ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++)
			new CPart(npartstot);
		sprintf(message,"made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
		CLog::Info(message);
	}
	CPartMap::iterator ppos=DeadPartMap.begin();
	for(ipart=0;ipart<5;ipart++){
		product[ipart]=ppos->second;
		++ppos;
	}
}

CAction* CB3D::GetDeadAction(){
	if(DeadActionMap.size()==0){
		for(int iaction=0;iaction<DELNACTIONSTOT*NSAMPLE;iaction++)
			new CAction(nactionstot);		
		//sprintf(message,"created %d new actions, nactionstot=%d\n",DELNACTIONSTOT*NSAMPLE,nactionstot);
		//CLog::Info(message);
	}
	return DeadActionMap.begin()->second;
}

void CB3D::CheckPartMap(){
	CPartMap::iterator iter;
	CPart *part;
	for(iter=PartMap.begin();iter!=PartMap.end();iter++){
		part=iter->second;
		if(part->currentmap!=&PartMap){
			part->Print();
			sprintf(message,"----- FAILED CheckPartMap-----\n");
			CLog::Fatal(message);
		}
	}
}

int CB3D::CountBaryons(){
	CPartMap::iterator iter;
	CPart *part;
	nbaryons=0;
	for(iter=PartMap.begin();iter!=PartMap.end();iter++){
		part=iter->second;
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
	}
	return nbaryons;
}

void CB3D::InitMuTCalc(){
	int ix,iy,ntau;
	CMuTInfo::b3d=this;
	CMuTInfo::NXY=parmap.getI("B3D_MUTCALC_NXY",24);
	CMuTInfo::DXY=parmap.getD("B3D_MUTCalc_DXY",1.0);
	CMuTInfo::NMINCALC=parmap.getD("B3D_MUTCALC_NMINCALC",5);
	CMuTInfo::taumin.resize(NXY);
	//CMuTInfo::massB.resize(8);
	//CMuTInfo::degenB.resize(8);
	for(ix=0;ix<NXY;ix++){
		CMuTInfo::taumin[ix].resize(NXY);
		for(iy=0;iy<NXY;iy++)
			CMuTInfo::taumin[ix][iy]=0.0;
	}
	CMuTInfo::NETEVENTS=0;
	ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
	muTinfo.resize(ntau);
	for(itau=0;itau<ntau;itau++){
		muTinfo[itau].resize(NXY);
		for(ix=0;ix<NXY;ix++){
			muTinfo[itau][ix].resize(NXY);
			for(iy=0;iy<NXY;iy++)
				muTinfo[itau][ix][iy]=new CMuTInfo((itau+0.5)*MUTCALC_DELTAU);
		}
	}
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->baryon>0){
			CMuTInfo::Bresinfo.push_back(resinfo);
		}
	}
	CMuTInfo::massB[0]=939.0;
}

void CB3D::CalcMuTU(){
	int ix,iy,iitau,ntau;
	ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
	for(iitau=0;iitau<ntau;iitau++){
		for(ix=0;ix<CMuTInfo::NXY;ix++){
			for(iy=0;iy<CMuTInfo::NXY;iy++){
				muTinfo[iitau][ix][iy]->CalcAllMuTU();
			}
		}
	}
}
