#include "b3d.h"
#include "part.h"
#include "cell.h"
#include "resonances.h"
#include "randy.h"
#include "constants.h"
#include "misc.h"

int CB3D::Annihilate(CPart *part1,CPart *part2,int &ndaughters,array<CPart*,5> &daughter){
	CPart *dptr;
	int netq,nets,nK0bar,nK0,nKplus,nKminus,npi0,npiplus,npiminus,npions,nkaons,npaircheck,qpions;
	FourVector *pa,*pb,pc,u;
	double ma,mb,q,cthet,sthet,phi;
	bool bjtranslate=false;
	int idaughter,iK,ipair,alpha;
	double Minv;
	double MM,P[4]={0.0},PP[4],T;
	const double g[4]={1.0,-1.0,-1.0,-1.0};
	CPartMap::iterator ppos;
	CB3DCell *newcell;
	ndaughters=0;
	if(BJORKEN && fabs(part1->eta)>ETAMAX)
		bjtranslate=true;

	netq = part1->resinfo->charge+part2->resinfo->charge;
	nets = part1->resinfo->strange+part2->resinfo->strange;
	nkaons=abs(nets);
	RECHARGE:
	nK0bar=nK0=nKplus=nKminus=npiplus=npiminus=npi0=0;
	for(iK=0;iK<nkaons;iK++){
		if(randy->ran()>0.5){
			if(nets>0)
				nK0+=1;
			else
				nK0bar+=1;
		}
		else{
			if(nets>0)
				nKplus+=1;
			else
				nKminus+=1;
		}
	}
	qpions=netq-nKplus+nKminus;
	if(qpions>0)
		npiplus+=qpions;
	else
		npiminus-=qpions;
	npions=5-nkaons;
	npi0=npions-npiplus-npiminus;
	if(netq!= nKplus+npiplus-nKminus-npiminus){
		sprintf(message,"charges don't add up\n");
		CLog::Fatal(message);
	}
	if(nets!= nKplus+nK0-nKminus-nK0bar){
		sprintf(message,"charges don't add up\n");
		CLog::Fatal(message);
	}

	if(npi0<0){
		goto RECHARGE;
	}
	npaircheck=0;
	if(npi0>=2)
		npaircheck=1;
	if(npi0>=4)
		npaircheck=2;
	for(ipair=0;ipair<npaircheck;ipair++){
		if(randy->ran()<0.66666666667){
			npiplus+=1;
			npiminus+=1;
			npi0-=2;
		}
	}
	ndaughters=npi0+npiplus+npiminus+nKplus+nKminus+nK0+nK0bar;
	if(ndaughters != 5){
		sprintf(message,"annihilation doesn't go to 5 particles\n");
		CLog::Fatal(message);
	}
	Minv=0.0;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part1->p[alpha]+part2->p[alpha];
		Minv+=g[alpha]*P[alpha]*P[alpha];
	}
	Minv=sqrt(Minv);
	T=Minv; // T will be KE of emitted particles
	idaughter=0;
	while(nK0bar>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(311);
		nK0bar-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nK0>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-311);
		nK0-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nKplus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(321);
		nKplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}while(nKminus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-321);
		nKminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npi0>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(111);
		T-=daughter[idaughter]->resinfo->mass;
		npi0-=1;
		idaughter+=1;
	}
	while(npiplus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(211);
		npiplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npiminus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-211);
		npiminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	T=0.5*T/double(npions+nkaons); // Pick temperature of 0.5*KE/particles
	
	// Do a 2-body decay for the last 2 particles
	DO_OVER:
	for(alpha=0;alpha<4;alpha++)
		PP[alpha]=0;
	for(idaughter=0;idaughter<ndaughters-2;idaughter++){
		dptr=daughter[idaughter];
		randy->generate_boltzmann(dptr->resinfo->mass,T,dptr->p);
		for(alpha=0;alpha<4;alpha++){
			PP[alpha]-=dptr->p[alpha];
		}
	}
	// PP is momentum of remaining two particles
	ma=daughter[ndaughters-2]->resinfo->mass;
	mb=daughter[ndaughters-1]->resinfo->mass;
	PP[0]=Minv+PP[0];
	if(PP[0]<ma+mb)
		goto DO_OVER;
	MM=0;
	for(alpha=0;alpha<4;alpha++)
		MM+=g[alpha]*PP[alpha]*PP[alpha];
	if(MM<(ma+mb)*(ma+mb))
		goto DO_OVER;
	MM=sqrt(MM);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=PP[alpha]/MM;
	pa=&daughter[ndaughters-2]->p;
	pb=&daughter[ndaughters-1]->p;
	cthet=1.0-2.0*randy->ran();
	sthet=sqrt(1.0-cthet*cthet);
	phi=2.0*PI*randy->ran();
	q=sqrt(Misc::triangle(MM,ma,mb));
	(*pa)[3]=q*cthet;
	(*pa)[1]=q*sthet*cos(phi);
	(*pa)[2]=q*sthet*sin(phi);
	(*pb)[3]=-(*pa)[3];
	(*pb)[2]=-(*pa)[2];
	(*pb)[1]=-(*pa)[1];
	(*pa)[0]=sqrt(ma*ma+q*q);
	(*pb)[0]=sqrt(mb*mb+q*q);
	Misc::Boost(u,*pa,pc);
	for(alpha=0;alpha<4;alpha++)
		(*pa)[alpha]=pc[alpha];
	Misc::Boost(u,*pb,pc);
	for(alpha=0;alpha<4;alpha++)
		(*pb)[alpha]=pc[alpha];

	for(alpha=0;alpha<4;alpha++){
		u[alpha]=P[alpha]/Minv;
		P[alpha]=0.0;
	}
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		pa=&daughter[idaughter]->p;
		Misc::Boost(u,*pa,pc);
		for(alpha=0;alpha<4;alpha++){
			(*pa)[alpha]=pc[alpha];
			P[alpha]+=(*pa)[alpha];
		}
	}
	Minv=0.0;
	for(alpha=0;alpha<4;alpha++)
		Minv+=g[alpha]*P[alpha]*P[alpha];
	Minv=sqrt(Minv);
	double rbar[4],etabar;

	rbar[1]=0.5*(part1->r[1]+part2->r[1]);
	rbar[2]=0.5*(part1->r[2]+part2->r[2]);
	etabar=0.5*(part1->eta+part2->eta);
	rbar[0]=tau*cosh(etabar);
	rbar[3]=tau*sinh(etabar);
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		dptr=daughter[idaughter];
		dptr->tau0=tau;
		dptr->eta=etabar;
		for(alpha=0;alpha<4;alpha++)
			dptr->r[alpha]=rbar[alpha];
		dptr->active=true;
		dptr->SetMass();
		dptr->SetY();
		if(bjtranslate)
			dptr->CyclicReset();
		newcell=dptr->FindCell();
		dptr->ChangeCell(newcell);

		dptr->ChangeMap(&PartMap);
		if(fabs(dptr->eta)>ETAMAX){
			sprintf(message,"eta out of range\n");
			CLog::Fatal(message);
		}
		if(dptr->p[0]<0.0){
			sprintf(message,"dptr->p[0]=%g\n",dptr->p[0]);
			CLog::Fatal(message);
		}
	}
	return ndaughters;
}

double CB3D::GetAnnihilationSigma(CPart *part1,CPart *part2,double &vrel){
	const double g[4]={1,-1,-1,-1};
	double Plab,p1dotp2,triangle,sigma_annihilation,rstrange,m1squared,m2squared;
	int alpha,ix,iy;
	double taumin1,taumin2;
	CMuTInfo::GetIxIy(part1->r[1],part1->r[2],ix,iy);
	taumin1=CMuTInfo::taumin[ix][iy];
	CMuTInfo::GetIxIy(part2->r[1],part2->r[2],ix,iy);
	taumin2=CMuTInfo::taumin[ix][iy];
	if(tau>taumin1 && tau>taumin2){
		return 0.0;
	}
	part1->SetMass(); part2->SetMass();
	m1squared=part1->msquared;
	m2squared=part2->msquared;
	p1dotp2=0.0;
	for(alpha=0;alpha<4;alpha++){
		p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
	}
	//Plab=sqrt((p1dotp2*p1dotp2/(part2->msquared))-part1->msquared);
	triangle=p1dotp2*p1dotp2-m1squared*m2squared;
	Plab=0.5*(m1squared+m2squared)*triangle/(m1squared*m2squared);
	Plab=sqrt(Plab);
	sigma_annihilation=6.7*pow(Plab/1000.0,-0.7)/double(NSAMPLE);
	rstrange=0.5*sqrt(sigma_annihilation);
	rstrange*=pow(ANNIHILATION_SREDUCTION,abs(part1->resinfo->strange))+pow(ANNIHILATION_SREDUCTION,abs(part2->resinfo->strange));
	sigma_annihilation=rstrange*rstrange;
	vrel=sqrt(triangle)/(part1->p[0]*part2->p[0]);
	return sigma_annihilation;
}

bool CB3D::CancelAnnihilation(CPart *part1,CPart *part2){
	double mupi,muK,muB,muBS,mutot,netS;
	double reduction_factor=1.0;
	CB3DCell *cell1,*cell2;
	cell1=part1->cell;
	cell2=part2->cell;
	if(cell1!=NULL && cell2!=NULL){
		int ix1,iy1,ix2,iy2,iitau;
		CMuTInfo::GetIxIy(part1->r[1],part1->r[2],ix1,iy1);
		CMuTInfo::GetIxIy(part1->r[1],part1->r[2],ix2,iy2);
		iitau=lrint(tau/MUTCALC_DELTAU);
		if(iitau<muTinfo.size()){
			if(muTinfo[iitau][ix1][iy1]->sufficientN && muTinfo[iitau][ix2][iy2]->sufficientN){
				mupi=0.5*(muTinfo[iitau][ix1][iy1]->mupi+muTinfo[iitau][ix2][iy2]->mupi);
				muK=0.5*(muTinfo[iitau][ix1][iy1]->muK+muTinfo[iitau][ix2][iy2]->muK);
				muB=0.5*(muTinfo[iitau][ix1][iy1]->muB+muTinfo[iitau][ix2][iy2]->muB);
				muBS=0.5*(muTinfo[iitau][ix1][iy1]->muK+muTinfo[iitau][ix2][iy2]->muBS);
				netS=fabs(part1->resinfo->strange+part2->resinfo->strange);
				mutot=(5.0-netS)*mupi+netS*muK-2.0*muB-netS*muBS;
				reduction_factor=1.0-exp(mutot);
				if(randy->ran()>reduction_factor){
					return false;
				}
				else{
					return true;
				}
			}
		}
	}
	return false;
}
