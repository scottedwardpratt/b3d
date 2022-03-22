#include <cmath>
#include <cstdlib>
#include "resonances.h"
#include "constants.h"
#include "parametermap.h"
#include "randy.h"
#include "misc.h"

using namespace std;

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,NResonances;
	double mothermass,netm,qR2,bmax;
	int ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,LDecay;
	int netq,netb,nets;
	string name, filename;
	CResInfo *resinfoptr=NULL;
	CBranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	FILE * decayinfofile;
	char dummy[200],cname[200];
	filename=parmap->getS("B3D_RESONANCES_INFO_FILE",string("../resinfo/resonances_standardhadrons.txt"));
	sprintf(message,"will read res info from %s\n",filename.c_str());
	CLog::Info(message);
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	MergeArray=new CMerge **[NResonances];
	SigmaMaxArray=new double *[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new CMerge *[NResonances];
		SigmaMaxArray[ires]=new double[NResonances];
		for(jres=0;jres<NResonances;jres++){
			MergeArray[ires][jres]=NULL;
			SigmaMaxArray[ires][jres]=0.0;
		}
	}
	for(ires=0;ires<NResonances;ires++){
		resinfoptr=new CResInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfoptr->code,&resinfoptr->mass,&resinfoptr->charge,&resinfoptr->baryon, &resinfoptr->strange,&resinfoptr->spin,&resinfoptr->G_Parity,&decay,&resinfoptr->width);

		resinfoptr->q[0]=resinfoptr->baryon+resinfoptr->charge;
		resinfoptr->q[1]=2.0*resinfoptr->baryon+resinfoptr->strange-resinfoptr->charge;
		resinfoptr->q[2]=-resinfoptr->strange;
		resinfoptr->netchi0+=resinfoptr->charge*resinfoptr->charge;
		resinfoptr->Btype=-1;
		if(!resinfoptr->decay)
			resinfoptr->netchi+=abs(resinfoptr->charge);
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfoptr->name=cname;
		resinfoptr->name.erase(0,resinfoptr->name.find_first_not_of(' '));
		resinfoptr->name.erase(resinfoptr->name.find_last_not_of(' ')+1);
		resinfoptr->SetBtype();
		resinfoptr->decay=bool(decay);
		if(!RESONANCE_DECAYS)
			resinfoptr->decay=false;
		resinfoptr->ires=ires;
		resinfoptr->branchlist.clear();
		resmap.insert(CResInfoPair(resinfoptr->code,resinfoptr));
	} 
	fclose(resinfofile);

	filename=parmap->getS("B3D_RESONANCES_DECAYS_FILE",string("../resinfo/decays_pdg_weak.txt"));
	sprintf(message,"will read decay info from %s\n",filename.c_str());
	CLog::Info(message);
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfoptr=GetResInfoPtr(mothercode);
		resinfoptr->minmass=0.0;
		bmax=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			bptr->resinfoptr.clear();
			resinfoptr->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfoptr.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfoptr[ibody]->charge;
				netb+=bptr->resinfoptr[ibody]->baryon;
				nets+=bptr->resinfoptr[ibody]->strange;
				netm+=bptr->resinfoptr[ibody]->mass;
			} 
			if(ichannel==0){
				resinfoptr->minmass=netm;
				if(resinfoptr->minmass < resinfoptr->mass-300.0)
					resinfoptr->minmass=resinfoptr->mass-300.0;
				if(resinfoptr->minmass < resinfoptr->mass-5.0*resinfoptr->width)
					resinfoptr->minmass=resinfoptr->mass-5.0*resinfoptr->width;
				if(resinfoptr->minmass>resinfoptr->mass-0.1){
					resinfoptr->minmass=resinfoptr->mass-0.1;
					resinfoptr->Print();
				}
			}
			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				sprintf(message,"Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				sprintf(message,"MOTHER=%d: (ichannel=%d, nbodies=%d):\n",resinfoptr->code,ichannel,nbodies);
				sprintf(message,"DAUGHTERS: ");
				for(ibody=0;ibody<nbodies;ibody++){
					sprintf(message,"%s %d ",message,bptr->resinfoptr[ibody]->code);
				}
				sprintf(message,"%s\n",message);
				CLog::Fatal(message);
			}

			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			for(ibody=0;ibody<nbodies;ibody++){
				if(!bptr->resinfoptr[ibody]->decay)
					resinfoptr->netchi+=abs(bptr->resinfoptr[ibody]->charge)*bptr->branching;
			}
			//store two body decays only
			if(nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				if(resinfoptr->mass>bptr->resinfoptr[0]->mass+bptr->resinfoptr[1]->mass){
					qR2=Misc::triangle(resinfoptr->mass,
					bptr->resinfoptr[0]->mass,bptr->resinfoptr[1]->mass);
					SigmaMaxArray[ires1][ires2]+=
						(bptr->branching*4.0*PI*HBARC*HBARC)*(2.0*resinfoptr->spin+1.0)/
							((2.0*bptr->resinfoptr[0]->spin+1.0)*(2.0*bptr->resinfoptr[1]->spin+1.0)*qR2);
					if(ires2!=ires1)
						SigmaMaxArray[ires2][ires1]=SigmaMaxArray[ires1][ires2];
				}
			}
			//if the total mass is smaller than the minimum required mass, replace it
			//if(netm<resinfoptr->minmass){
			//	resinfoptr->minmass=netm;
				//resinfoptr->bptr_minmass=bptr;
			//}
			// switch places to make sure first branch has largest 
			if(bptr->branching>bmax){
				bmax=bptr->branching>bmax;
				if(ichannel>0){
					firstbptr=resinfoptr->branchlist[0];
					resinfoptr->branchlist[0]=bptr;
					resinfoptr->branchlist[ichannel]=firstbptr;
				}
			}
		}
	}
	fclose(decayinfofile);
}
