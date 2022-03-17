#include "action.h"
#include "cell.h"
#include "mutinfo.h"
#include "part.h"
#include "resonances.h"

void CAction::PerformMuTCalcUpdateNPE(){
	int ix,iy,itau,btype,pid;
	itau=lrint(floor(tau/b3d->MUTCALC_DELTAU));
	CPartMap::iterator ppos;
	CPart *part;
	CResInfo *resinfo;
	CMuTInfo *mti;
	double gamma,gammav,E,px,py,eta,t,x,y;

	for(ppos=b3d->PartMap.begin();ppos!=b3d->PartMap.end();++ppos){
		part=ppos->second;
		
		resinfo=part->resinfo;
		pid=abs(resinfo->code);
		btype=CMuTInfo::GetBtype(pid);
		if(pid==111 || pid==211 || pid==311 || pid==321 || btype>=0){
			eta=part->GetEta(tau);
			t=tau*cosh(eta);
			x=part->r[1]+(t-part->r[0])*part->p[1]/part->p[0];
			y=part->r[2]+(t-part->r[0])*part->p[2]/part->p[0];
			CMuTInfo::GetIxIy(x,y,ix,iy);
			if(ix<CMuTInfo::NXY && iy<CMuTInfo::NXY){
				if(b3d->tau>CMuTInfo::taumin[ix][iy]){
					px=part->p[1];
					py=part->p[2];
					if(x<0.0)
						px=-px;
					if(y<0.0)
						py=-py;
					if(ix<CMuTInfo::NXY && iy<CMuTInfo::NXY){
						mti=b3d->muTinfo[itau][ix][iy];

						if(resinfo->code==111 || abs(resinfo->code)==211){
							gamma=cosh(eta);
							gammav=sinh(eta);
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
							gamma=cosh(eta);
							gammav=sinh(eta);
							mti->NK+=1;
							mti->PxK+=px;
							mti->PyK+=py;
							E=gamma*part->p[0]-gammav*part->p[3];
							mti->EK+=E;
							mti->TxxK+=px*px/E;
							mti->TyyK+=py*py/E;
							mti->TxyK+=px*py/E;
						}
						else if(btype>=0){
							gamma=cosh(eta);
							gammav=sinh(eta);
							mti->NB[btype]+=1;
							mti->PxB[btype]+=px;
							mti->PyB[btype]+=py;
							E=gamma*part->p[0]-gammav*part->p[3];
							mti->EB[btype]+=E;
							mti->TxxB[btype]+=px*px/E;
							mti->TyyB[btype]+=py*py/E;
							mti->TxyB[btype]+=px*py/E;
						}
					}
				}
			}
		}
	}
}