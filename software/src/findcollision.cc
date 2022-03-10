#include "b3d.h"
#include "part.h"
#include "cell.h"
#include "resonances.h"
#include "constants.h"

double CB3D::GetPiBSquared(CPart *part1,CPart *part2){
	double denom;
	CB3DCell *cell1=part1->cell;
	CB3DCell *cell2=part2->cell;
	double p1dotp2=0.0,p1dotr=0.0,p2dotr=0.0,m1squared,m2squared,rsquared=0.0;
	double tau1,eta1,y1,mt;
	FourVector p1,p2,r1,r2,r;
	double pibsquared;
	const int g[4]={1,-1,-1,-1};
	int alpha;
	bool bjtranslate=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	for(alpha=0;alpha<4;alpha++){
		r1[alpha]=part1->r[alpha];
		r2[alpha]=part2->r[alpha];
		p1[alpha]=part1->p[alpha];
		p2[alpha]=part2->p[alpha];
	}
	m1squared=part1->msquared;
	m2squared=part2->msquared;
	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		tau1=part1->tau0;
		if(cell1->ieta==0){
			eta1=part1->eta+2.0*ETAMAX; y1=part1->y+2.0*ETAMAX;
		}
		else{
			eta1=part1->eta-2.0*ETAMAX; y1=part1->y-2.0*ETAMAX;
		}
		r1[0]=tau1*cosh(eta1);
		r1[3]=tau1*sinh(eta1);
		mt=part1->resinfo->mass;
		mt=sqrt(mt*mt+p1[1]*p1[1]+p1[2]*p1[2]);
		p1[0]=mt*cosh(y1);
		p1[3]=mt*sinh(y1);
	}
	
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=r1[alpha]-r2[alpha];
		rsquared+=g[alpha]*r[alpha]*r[alpha];
		p1dotp2+=g[alpha]*p1[alpha]*p2[alpha];
		p1dotr+=g[alpha]*p1[alpha]*r[alpha];
		p2dotr+=g[alpha]*p2[alpha]*r[alpha];
	}
		
	denom=p1dotp2*p1dotp2-m1squared*m2squared;
	pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*m2squared-p2dotr*p2dotr*m1squared)/denom);
	if(bjtranslate)
		part1->BjorkenUnTranslate();
	return pibsquared;
}

double CB3D::GetSigma(CPart *part1,CPart *part2,
		double &sigma_scatter,double &sigma_merge,double &sigma_annihilation,double &sigma_inel,
		vector<double> &dsigma_merge){
	double sigma=0,MR,M,Gamma,b,jR,j1,j2,qR2,q2,q3,q4,p1dotp2,P2,tan2delta,G,vrel,dsigma;
	int L_merge,alpha,ir1,ir2;
	const double g[4]={1,-1,-1,-1};
	bool bjtranslate=false;
	CMerge *merge;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	
	if((part1->balanceID>=0 && part2->balanceID<0) || (part2->balanceID>=0 && part1->balanceID<0)){
		sigma_scatter=SIGMABF;
	}
	else{
		sigma_scatter=SIGMADEFAULT/double(NSAMPLE);
	}

	if(BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		sigma_annihilation=GetAnnihilationSigma(part1,part2,vrel);
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
		p1dotp2=0.0;
		for(alpha=0;alpha<4;alpha++){
			p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
		}
		P2=part1->msquared+part2->msquared+2.0*p1dotp2;
		M=sqrt(P2);
		if(M>merge->resinfo->minmass){
			q2=Misc::triangle2(P2,part1->msquared,part2->msquared);
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
				*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))))/double(NSAMPLE);
			dsigma_merge.push_back(dsigma);
			sigma_merge+=dsigma;
		}
		merge=merge->next;
	}

	if(bjtranslate)
		part1->BjorkenUnTranslate();
	sigma=sigma_scatter+sigma_inel+sigma_annihilation+sigma_merge;
	return sigma;
}

bool CB3D::FindCollision(CPart *part1,CPart *part2,double &taucoll){
	if((part1->balanceID>=0) && (part2->balanceID>=0)){
		return false;
	}
	double denom,sigmamax;
	double sigma_scatter=0.0,sigma_merge=0.0,sigma_annihilation=0.0,sigma_inel=0.0;
	vector<double> dsigma_merge;
	CResInfo *resinfo1,*resinfo2;
	bool collide=false;
	CB3DCell *cell1=part1->cell;
	CB3DCell *cell2=part2->cell;
	double p1dotp2=0.0,p1dotr=0.0,p2dotr=0.0,m1squared,m2squared,rsquared=0.0;
	double tau1,tau2,eta1,y1,mt;
	FourVector p1,p2,r1,r2,r;
	double t1,t2,z1,z2;
	double pibsquared;
	const int g[4]={1,-1,-1,-1};
	int alpha;
	dsigma_merge.clear();
	for(alpha=0;alpha<4;alpha++){
		r1[alpha]=part1->r[alpha];
		r2[alpha]=part2->r[alpha];
		p1[alpha]=part1->p[alpha];
		p2[alpha]=part2->p[alpha];
	}
	m1squared=part1->msquared;
	m2squared=part2->msquared;
	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		tau1=part1->tau0;
		if(cell1->ieta==0){
			eta1=part1->eta+2.0*ETAMAX; y1=part1->y+2.0*ETAMAX;
		}
		else{
			eta1=part1->eta-2.0*ETAMAX; y1=part1->y-2.0*ETAMAX;
		}
		r1[0]=tau1*cosh(eta1);
		r1[3]=tau1*sinh(eta1);
		mt=part1->resinfo->mass;
		mt=sqrt(mt*mt+p1[1]*p1[1]+p1[2]*p1[2]);
		p1[0]=mt*cosh(y1);
		p1[3]=mt*sinh(y1);
	}
	
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=r1[alpha]-r2[alpha];
		rsquared+=g[alpha]*r[alpha]*r[alpha];
		p1dotp2+=g[alpha]*p1[alpha]*p2[alpha];
		p1dotr+=g[alpha]*p1[alpha]*r[alpha];
		p2dotr+=g[alpha]*p2[alpha]*r[alpha];
	}
		
	denom=p1dotp2*p1dotp2-m1squared*m2squared;
	pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*m2squared-p2dotr*p2dotr*m1squared)/denom);

	resinfo1=part1->resinfo;
	resinfo2=part2->resinfo;
	
	sigmamax=GetSigma(part1,part2,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
		dsigma_merge);

	/*
	if(part1->balanceID<0 && part2->balanceID<0){
		sigmamax=SIGMADEFAULT+reslist->SigmaMaxArray[resinfo1->ires][resinfo2->ires]+SIGMAINELASTIC;
	}
	else{
		sigmamax=SIGMABF;
	}
	if(BARYON_ANNIHILATION && resinfo1->baryon*resinfo2->baryon<0)
		sigmamax+=30.0; // cuts off annihilation Xsection at 300 mb
		*/
		
	if(pibsquared<sigmamax/double(NSAMPLE)){
		
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
				/*
				if(tau<20)
				printf("-----------------\ntau=%g, taucoll=%g\ncell exit times are %g, %g\n",
					tau,taucoll,part1->tauexit,part2->tauexit);
					*/
				if(taucoll>tau && taucoll<part1->tauexit && taucoll<part2->tauexit && taucoll<TAUCOLLMAX){
					AddAction_Collision(part1,part2,taucoll,pibsquared,
						sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
						dsigma_merge);
						//if(tau<20)
						//printf("XXXXXX found one\n");
					collide=true;
				}
			}
		}
	}
	return collide;
}

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
