#include "b3d.h"
#include "part.h"
#include "cell.h"
#include "resonances.h"
#include "constants.h"

bool CB3D::CheckKinematics(CPart *part1,CPart *part2,
	double &Minv2,double &pibsquared,double &taucoll){
	double denom;
	double p1dotp2=0.0,p1dotr=0.0,p2dotr=0.0,m1squared,m2squared,rsquared=0.0;
	double t1,t2,z1,z2,tau1,tau2;
	FourVector r,r1,r2,p1,p2;
	const int g[4]={1,-1,-1,-1};
	int alpha;
	bool bjtranslate=false;
	bool possible=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	m1squared=part1->msquared;
	m2squared=part2->msquared;
	
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=part1->r[alpha]-part2->r[alpha];
		rsquared+=g[alpha]*r[alpha]*r[alpha];
		p1dotp2+=g[alpha]*part1->p[alpha]*part2->p[alpha];
		p1dotr+=g[alpha]*part1->p[alpha]*r[alpha];
		p2dotr+=g[alpha]*part2->p[alpha]*r[alpha];
	}
	Minv2=m1squared+m2squared+2.0*p1dotp2;
		
	denom=p1dotp2*p1dotp2-m1squared*m2squared;
	pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*m2squared-p2dotr*p2dotr*m1squared)/denom);

	for(alpha=0;alpha<4;alpha++){
		p1[alpha]=part1->p[alpha];
		p2[alpha]=part2->p[alpha];
		r1[alpha]=part1->r[alpha];
		r2[alpha]=part2->r[alpha];
	}

	t1=p1[0]*(p1dotr*m2squared-p2dotr*p1dotp2)/denom;
	t2=-p2[0]*(p2dotr*m1squared-p1dotr*p1dotp2)/denom;

	if(t1+r1[0]>0 && t2+r2[0]>0 && (t1+t2)>0.0){
		z1=r1[3]+(p1[3]/p1[0])*t1;
		z2=r2[3]+(p2[3]/p2[0])*t2;
		t1+=r1[0];
		t2+=r2[0];
			/*
			if(pibsquared<sigmamax && tau<20){
				printf("b=%g, tau1=%g, tau2=%g actionmother1=%d actiomother2=%d\n",
					sqrt(pibsquared/PI),part1->tau0,part2->tau0,part1->actionmother, part2->actionmother);
				printf("bmax=%g\n",sqrt(sigmamax/PI));
			}
			*/
		if(fabs(z1)<t1 && fabs(z2)<t2){
			tau1=sqrt(t1*t1-z1*z1);
			tau2=sqrt(t2*t2-z2*z2);
			taucoll=0.5*(tau1+tau2);
			if(taucoll>part1->tau0 && taucoll>part2->tau0 && taucoll<part1->tauexit && taucoll<part2->tauexit)
				possible=true;
		}
	}
	if(bjtranslate)
		part1->BjorkenUnTranslate();
	return possible;
}

double CB3D::GetSigma(CPart *part1,CPart *part2,double Minv2,
		double &sigma_scatter,double &sigma_merge,double &sigma_annihilation,double &sigma_inel,
		vector<double> &dsigma_merge){
	double sigmatot=0,MR,M,Gamma,b,jR,j1,j2,qR2,q2,q3,q4,tan2delta,G,dsigma;
	int L_merge,ir1,ir2;
	CMerge *merge;
	
	if((part1->balanceID>=0 && part2->balanceID<0) || (part2->balanceID>=0 && part1->balanceID<0)){
		sigma_scatter=SIGMABF;
	}
	else{
		sigma_scatter=SIGMADEFAULT;
	}

	if(BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		sigma_annihilation=GetAnnihilationSigma(part1,part2);
	}
	else
		sigma_annihilation=0.0;

	if(SIGMAINELASTIC)
		sigma_inel=SIGMAINELASTIC;
	else
		sigma_inel=0.0;

	ir1=part1->resinfo->ires;
	ir2=part2->resinfo->ires;
	merge=reslist->MergeArray[ir1][ir2];
	dsigma_merge.clear();
	while(merge!=NULL){
		j1=part1->resinfo->spin;
		j2=part2->resinfo->spin;
		M=sqrt(Minv2);
		if(M>merge->resinfo->minmass){
			q2=Misc::triangle2(Minv2,part1->msquared,part2->msquared);
			Gamma=merge->resinfo->width;
			b=merge->branching;
			jR=merge->resinfo->spin;
			MR=merge->resinfo->mass;
			L_merge = merge->L;
			qR2=Misc::triangle2(MR*MR,part1->msquared,part2->msquared);
			q3=pow(q2/qR2,(2*L_merge + 1)/2);
			q4=pow(q2/qR2,(2*L_merge)/2);
		//q3=q4=1.0;
			G=Gamma*(MR/M)*q3*1.2/(1.0+0.2*q4);
			tan2delta=pow(0.5*G/(M-MR),2);
			dsigma=b*((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
				*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))));
			dsigma_merge.push_back(dsigma);
			sigma_merge+=dsigma;
		}
		merge=merge->next;
	}
	sigmatot=sigma_scatter+sigma_inel+sigma_annihilation+sigma_merge;
	return sigmatot;
}

bool CB3D::FindCollision(CPart *part1,CPart *part2,double &taucoll){
	if((part1->balanceID>=0) && (part2->balanceID>=0)){
		return false;
	}
	bool checkpossible;
	double sigma_scatter=0.0,sigma_merge=0.0,sigma_annihilation=0.0,sigma_inel=0.0;
	vector<double> dsigma_merge;
	bool collide=false;
	double pibsquared,Minv2,sigmatot;
	dsigma_merge.clear();

	checkpossible=CheckKinematics(part1,part2,Minv2,pibsquared,taucoll);
	if(checkpossible){
		sigmatot=GetSigma(part1,part2,Minv2,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
			dsigma_merge);
		if(pibsquared<sigmatot/double(NSAMPLE)){
			AddAction_Collision(part1,part2,taucoll,pibsquared,
						sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
						dsigma_merge);
			collide=true;
		}
	}
	return collide;
}


	

	
	

	/*
	resinfo1=part1->resinfo;
	resinfo2=part2->resinfo;
	if(part1->balanceID<0 && part2->balanceID<0){
		sigmamax=SIGMADEFAULT+reslist->SigmaMaxArray[resinfo1->ires][resinfo2->ires]+SIGMAINELASTIC;
	}
	else{
		sigmamax=SIGMABF;
	}
	if(BARYON_ANNIHILATION && resinfo1->baryon*resinfo2->baryon<0)
		sigmamax+=30.0; // cuts off annihilation Xsection at 300 mb
		*/
		

void CB3D::FindAllCollisions(){
	double taucoll;
	int nfound=0;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	CActionMap::iterator epos;
	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
		part1->active=true;
		part1->ChangeCell(part1->FindCell());
		if(part1->cell!=NULL){
			part1->FindCellExit();
		}
		if(part1->resinfo->decay)
			part1->FindDecay();
	}
	ppos2=PartMap.begin();
	++ppos2;
	while(ppos2!=PartMap.end()){
		part2=ppos2->second;
		ppos1=PartMap.begin();
		while(ppos1!=ppos2){
			part1=ppos1->second;
			if(part1->balanceID<0 || part2->balanceID<0){
				if(part1->cell!=NULL & part2->cell!=NULL){
					if(FindCollision(part1,part2,taucoll))
						nfound+=1;
				}
			}
			++ppos1;
		}
		++ppos2;
	}
}
