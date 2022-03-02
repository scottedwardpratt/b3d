#include "b3d.h"
#include "cell.h"
#include "mutinfo.h"
#include "part.h"
#include "mutinfo.h"
#include "resonances.h"

void CAction::PerformMuTCalcUpdateNPE(){
	int ix,iy,itau;
	itau=lrint(floor(tau/b3d->MUTCALC_DELTAU));
	CPartMap::iterator ppos;
	CPart *part;
	CResInfo *resinfo;
	CMuTInfo *mti;
	int S;
	double gamma,gammav,E,px,py;

	for(ppos=b3d->PartMap.begin();ppos!=b3d->PartMap.end();++ppos){
		part=ppos->second;
		part->Propagate(tau);
		resinfo=part->resinfo;
		CMuTInfo::GetIxIy(part->r[1],part->r[2],ix,iy);
		if(ix<CMuTInfo::NXY && iy<CMuTInfo::NXY){
			if(b3d->tau>CMuTInfo::taumin[ix][iy]){
				px=part->p[1];
				py=part->p[2];
				if(part->r[1]<0.0)
					px=-px;
				if(part->r[2]<0.0)
					py=-py;
				if(ix<CMuTInfo::NXY && iy<CMuTInfo::NXY){
					mti=b3d->muTinfo[itau][ix][iy];

					if(resinfo->code==111 || abs(resinfo->code)==211){
						gamma=cosh(part->eta);
						gammav=sinh(part->eta);
						mti->Npi+=1;
						mti->Pxpi+=px;
						mti->Pypi+=py;
						E=gamma*part->p[0]-gammav*part->p[3];
						mti->Epi+=E;
						mti->Txxpi+=px*px/E;
						mti->Tyypi+=py*py/E;
						mti->Txypi+=px*py/E;
					}
					else if(abs(part->resinfo->code)==321 || abs(part->resinfo->code)==311){
						gamma=cosh(part->eta);
						gammav=sinh(part->eta);
						mti->NK+=1;
						mti->PxK+=px;
						mti->PyK+=py;
						E=gamma*part->p[0]-gammav*part->p[3];
						mti->EK+=E;
						mti->TxxK+=px*px/E;
						mti->TyyK+=py*py/E;
						mti->TxyK+=px*py/E;
					}
					else if(resinfo->baryon!=0){
						gamma=cosh(part->eta);
						gammav=sinh(part->eta);
						mti->NB+=1;
						mti->PxB+=px;
						mti->PyB+=py;
						E=gamma*part->p[0]-gammav*part->p[3];
						mti->EB+=E;
						S=abs(resinfo->strange);
						mti->NBS+=S;
						mti->TxxB+=px*px/E;
						mti->TyyB+=py*py/E;
						mti->TxyB+=px*py/E;
					}
				}
			}
		}
	}
	b3d->FindAllCollisions();
}