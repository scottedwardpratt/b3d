#ifndef __ACTION_PERFORM_REGEN_CC__
#define __ACTION_PERFORM_REGEN_CC__

#include "b3d.h"
#include "cell.h"
#include "mutinfo.h"

void CAction::PerformMuTCalcUpdateNPE(){
	int ix,iy,ieta,itau;
	itau=lrint(floor(tau/CMuTInfo::DELTAU));
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			for(ieta=0;ieta<2*b3d->NETA;ieta++){
				b3d->muTinfo[itau][ix][iy]->UpdateNPE(b3d->cell[ix][iy][ieta]);
			}
		}
	}
}

#endif
