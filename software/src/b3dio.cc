#ifndef __B3DIO_CC__
#define __B3DIO_CC__

#include "b3d.h"
#include "part.h"
#include "resonances.h"
#include "randy.h"
#include "cell.h"
#include "misc.h"

using namespace std;

double CB3D::WriteOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	double dnchdy=0;
	int ipart;
	CPart *part;
	CPartMap::iterator ppos;
	
	int nparts=PartMap.size();
	sprintf(message,"writing %d particles to %s\n",nparts,oscarfilename.c_str());
	b3dlog->Info(message);
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"b3d output\n");
		}
	}
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parmap.getD("GLAUBER_B",0.0),
	parmap.getD("GLAUBER_B",0.0));
	ppos=PartMap.begin();
	for(ipart=0;ipart<nparts;ipart++){
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->code;
			bpart.tau=part->tau0;
			bpart.x=part->r[1];
			bpart.y=part->r[2];
			bpart.eta=part->eta;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			bpart.weight=part->weight;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else
			fprintf(oscarfile,"%5d %5d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %g\n",
		ipart,part->resinfo->code,part->p[1],part->p[2],part->p[3],part->p[0],sqrt(part->msquared),part->r[1],part->r[2],part->r[3],part->r[0],part->weight);
		if(ppos==PartMap.end()){
			sprintf(message,"ppos shouldn't be here\n");
			b3dlog->Fatal(message);
		}
		++ppos;
	}
	return dnchdy/(2.0*ETAMAX);
}

void CB3D::ReadOSCARHeader(){
	int ndead=3,idead;
	char dummy[200];
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.txt";
	if(BINARY_RW)
		oscarfile=fopen(oscarfilename.c_str(),"rb");
	else{
		oscarfile=fopen(oscarfilename.c_str(),"r");
		for(idead=0;idead<ndead;idead++)
			fgets(dummy,200,oscarfile);
	}
}

int CB3D::ReadOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	CResInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,ID;
	int nparts_read,ipart=0;
	int ievent_read;
	double bmin,bmax; // impact parameter
	CPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			sprintf(message,"trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			b3dlog->Fatal(message);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			tau0=bpart.tau;
			r[1]=bpart.x;
			r[2]=bpart.y;
			eta=bpart.eta;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			resinfo=reslist->GetResInfoPtr(ID);
			mass=resinfo->mass;
			weight=bpart.weight;
		}
		else{
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
			&ipart,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0],&weight);
			tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
			eta=asinh(r[3]/tau0);
			rapidity=asinh(p[3]/p[0]);
		}
		mother->Init(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight);
	}
	return nparts_read;
}

double CB3D::WriteBalanceParts(int ievent){
	CB3DBinaryBalancePartInfo bpart;
	double sigma=0;
	int nsigma=0;
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"b3d output\n");
		}
	}
	int nparts=PartMap.size();
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parmap.getD("GLAUBER_B",0.0),
	parmap.getD("GLAUBER_B",0.0));
	double dnchdy=0,rapidity;
	int ipart;
	CPart *part;
	CPartMap::iterator ppos;
	ipart=0;
	ppos=PartMap.begin();
	do{
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->code;
			bpart.balanceID=part->balanceID;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			sigma+=part->y*part->y;
			nsigma+=1;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else{
			rapidity=0.5*log((part->p[0]+part->p[3])/(part->p[0]-part->p[3]));
			if(bpart.balanceID!=-1){
				fprintf(oscarfile,"%5d %5d %7d %12.6e %12.6e %12.6e\n",
				ipart,part->resinfo->code,part->balanceID,part->p[1],part->p[2],rapidity);
			}
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	sprintf(message,"WriteBalance -- sigma=%g\n",sqrt(sigma/double(nsigma)));
	b3dlog->Info(message);
	return dnchdy/(2.0*ETAMAX);
}

int CB3D::ReadBalanceParts(int ievent){
	CB3DBinaryBalancePartInfo bpart;
	CResInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,ID,balanceID;
	int nparts_read,ipart=0;
	int ievent_read;
	double bmin,bmax; // impact parameter
	CPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			sprintf(message,"trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			b3dlog->Fatal(message);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			balanceID=bpart.balanceID;
		}
		else{
			sprintf(message,"should only work for binary rw\n");
			b3dlog->Fatal(message);
			//fscanf(oscarfile,"%d %d %lf %lf %lf",&ipart,&ID,&p[1],&p[2],&rapidity);
		}
		r[1]=r[2]=0.0;
		eta=rapidity;
		while(fabs(eta)>ETAMAX)
			eta-=2.0*ETAMAX*fabs(eta)/eta;
		tau0=10.0;
		weight=1.0;
		resinfo=reslist->GetResInfoPtr(ID);
		mass=resinfo->mass;
		if(resinfo->charge!=0 || resinfo->decay){
			mother->InitBalance(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight,balanceID);
		}
	}
	return nparts_read;
}

void CB3D::WriteDens(){
	string densfilename="model_output/"+run_name+"/"+qualifier+"/dens.txt";
	FILE *densfile = fopen(densfilename.c_str(),"w");
	fprintf(densfile,"#ix iy  dens[itau=0] dens[itau=1]...\n");
	double dxy;
	int ix,iy,ieta,iitau;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			fprintf(densfile,"%3d %3d",ix,iy);
			for(iitau=0;iitau<DENSWRITE_NTAU;iitau++){
				dxy=0.0;
				for(ieta=0;ieta<2*NETA;ieta++){
					dxy+=cell[ix][iy][ieta]->dens[iitau];
				}
				fprintf(densfile," %6.0f",dxy);
			}
			fprintf(densfile,"\n");
		}
	}
	fclose(densfile);
}

#endif
