#ifndef __DECAY_CC__
#define __DECAY_CC__

#include "b3d.h"
#include "randy.h"
#include "resonances.h"
#include "part.h"
#include "misc.h"
#include "constants.h"

void CB3D::Decay(CPart *mother,int &nbodies,array<CPart *,5> &daughter){
	int ibody,alpha;
	double mtot;
	CPart *dptr;
	vector<double> mass(6);
	vector<FourVector> p(5);
	FourVector u,pprime;

	mass[0]=mother->GetMass();
	
	/* Make masses for daughters */
	do{
		mtot=0.0;
		for(ibody=0;ibody<nbodies;ibody++){
			if(daughter[ibody]->resinfo->decay)
				mass[ibody+1]=daughter[ibody]->resinfo->GenerateMass();
			else
				mass[ibody+1]=daughter[ibody]->resinfo->mass;
			mtot+=mass[ibody+1];
		}
	}while(mtot>mass[0]);
	
	if(nbodies==2){
		decay_nbody->SetMasses2(mass[0],mass[1],mass[2]);
		decay_nbody->GenerateMomenta2(p[0],p[1]);
	}
	else if(nbodies==3){
		decay_nbody->SetMasses3(mass[0],mass[1],mass[2],mass[3]);
		decay_nbody->GenerateMomenta3(p[0],p[1],p[2]);
	}
	else if(nbodies==4){
		decay_nbody->SetMasses4(mass[0],mass[1],mass[2],mass[3],mass[4]);
		decay_nbody->GenerateMomenta4(p[0],p[1],p[2],p[3]);
	}
	else{
		decay_nbody->SetMasses(nbodies,mass);
		decay_nbody->GenerateMomenta(p);
	}
	
	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=mother->p[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::lorentz(u,p[ibody],pprime);
		dptr->active=true;
		for(alpha=0;alpha<4;alpha++)
			dptr->p[alpha]=pprime[alpha];
		dptr->CopyPositionInfo(mother);
		//dptr->msquared=pow(dptr->resinfo->mass,2);
		dptr->SetMass();
		//dptr->Setp0();
		dptr->SetY();
		dptr->eta0=mother->eta0;
		dptr->phi0=mother->phi0;
	}
}

#endif
