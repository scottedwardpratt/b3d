#include "b3d.h"
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include "b3d.h"
//#include "hist.h"
#include "pow.h"
#include "part.h"
#include "resonances.h"
#include "cell.h"
#include "sampler.h"
#include "balancearrays.h"
#include "randy.h"
#include "hyper.h"
#include "action.h"

CB3D::CB3D(){
};

CB3D::CB3D(string run_name_set){
	run_name=run_name_set;
	string parsfilename,dirname;
	dirname="model_output/"+run_name;
	parsfilename="model_output/fixed_parameters.txt";
	sprintf(message,"reading %s\n",parsfilename.c_str());
	CLog::Info(message);
	parmap.ReadParsFromFile(parsfilename);
	parsfilename=dirname+"/parameters.txt";
	sprintf(message,"reading %s\n",parsfilename.c_str());
	CLog::Info(message);
	parmap.ReadParsFromFile(parsfilename);
	CopyParMapPars();
	string command="mkdir -p model_output/"+run_name;
	system(command.c_str());
	if(BFCALC){
		parmap.ReadParsFromFile("udsdata/udsparameters.txt");
		balancearrays=new CBalanceArrays(this);
	}
	ibalmax=0;
	npartstot=nactionstot=0;
	CResList::b3d=this;
	CPart::b3d=this;
	tau=0.0;
	chitotH.setZero();
	chitotQ.setZero();
	npartstot=nactionstot=0;
	PartMap.clear();
	DeadPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	randy=new CRandy(-1234);
	CAction::b3d=this;
	CPart::b3d=this;
	CB3DCell::b3d=this;
	oscarfile=NULL;
	reslist=new CResList(&parmap);
	sampler=new CSampler(this);
	decay_nbody=new CDecay_NBody(randy);
}

void CB3D::CopyParMapPars(){
	NACTIONSMAX=parmap.getI("B3D_NACTIONSMAX",100000);
	NPARTSMAX=parmap.getI("B3D_NPARTSMAX",200000);
	TAUCOLLMAX=parmap.getD("B3D_TAUCOLLMAX",50.0);
	DENSWRITE=parmap.getB("B3D_DENSWRITE",false);
	DENSWRITE_NTAU=parmap.getI("B3D_DENSWRITE_NTAU",20);
	DENSWRITE_DELTAU=parmap.getD("B3D_DENSWRITE_DELTAU",1.0);
	input_dataroot=parmap.getS("B3D_INPUT_DATAROOT","udsdata/b3d");
	output_dataroot=parmap.getS("B3D_OUTPUT_DATAROOT","udsdata/b3d");
	NSAMPLE=parmap.getI("B3D_NSAMPLE",1);
	NSAMPLE_UDS2BAL=parmap.getI("NSAMPLE_UDS2BAL",1);
	ERROR_PRINT=parmap.getB("B3D_ERROR_PRINT",true);
	XYMAX=parmap.getD("B3D_XYMAX",15);
	ETAMAX=parmap.getD("B3D_ETAMAX",1.0);
	NETA=parmap.getI("B3D_NETA",10);
	NXY=parmap.getI("B3D_NXY",10);
	SIGMAMAX=parmap.getD("B3D_SIGMAMAX",30);
	NRINGSMAX=parmap.getI("B3D_NRINGSMAX",200);
	NPRCELLSMAX=parmap.getI("B3D_PR_NPRCELLSMAX",100000);
	COLLISIONS=parmap.getB("B3D_COLLISIONS",true);
	BJORKEN=parmap.getB("B3D_BJORKEN",true);
	BFCALC=parmap.getB("B3D_BFCALC",true);
	SIGMADEFAULT=parmap.getD("B3D_SIGMADEFAULT",1.0);
	SIGMAINELASTIC=parmap.getD( "B3D_SIGMAINELASTIC",1.0);
	INELASTIC=parmap.getB( "B3D_INELASTIC",false);
	NBOSE=parmap.getI("B3D_NBOSE",1);
	Q0=parmap.getD( "B3D_INELASTIC_Q0", 0);
	COOPERFRYE_CREATENEGPARTS=parmap.getB("B3D_COOPERFRYE_CREATENEGPARTS","false");
	COOPERFRYE_WMAX=parmap.getI("B3D_COOPERFRYE_WMAX",1);
	HYDRO_OCTANT_SYMMETRY=parmap.getI("HYDRO_OCTANT",2);
	HYDRO_PURE_BJORKEN=parmap.getB("HYDRO_PURE_BJORKEN",false);
	BARYON_ANNIHILATION=parmap.getB("B3D_BARYON_ANNIHILATION",false);
	DELNPARTSTOT=parmap.getD("B3D_DELNPARTSTOT",1000);
	DELNACTIONSTOT=parmap.getD("B3D_DELNACTIONSTOT",2000);
	USE_OLD_SAMPLER=parmap.getB("B3D_USE_OLD_SAMPLER",false);
	BINARY_RW=parmap.getB("B3D_BINARY_RW",false);
	MUTCALC=parmap.getB("B3D_MUTCALC",false);
	MUTCALC_DELTAU=parmap.getD("B3D_MUTCALC_DELTAU",0.5);
	ANNIHILATION_SREDUCTION=parmap.getD("B3D_ANNIHILATION_SREDUCTION",1.0);
	RESONANCE_DECAYS=parmap.getB("B3D_RESONANCE_DECAYS",true);
	SIGMAMAX=SIGMAMAX/double(NSAMPLE);
	SIGMABF=parmap.getD("B3D_SIGMABF",2.3);
	NPARTSMAX*=NSAMPLE;
	NACTIONSMAX*=NSAMPLE;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);
	BALANCE_CALC=parmap.getB("B3D_BALANCE_CALC",false);
}

void CB3D::InitCascade(){
	// First initialize cells
	int ix,iy,ieta,jx,jy,jeta,iitau,imax;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CB3DCell *c;
	CInelasticList::b3d=this;
	CInelasticList::UseFile = false;
	CInelasticList::UseInelasticArray = false;
	if(INELASTIC)
		inelasticlist = new CInelasticList();
	
	cell.resize(2*NXY);
	for(ix=0;ix<2*NXY;ix++){
		xmin=-XYMAX+ix*DXY;
		xmax=xmin+DXY;
		cell[ix].resize(2*NXY);
		for(iy=0;iy<2*NXY;iy++){
			ymin=-XYMAX+iy*DXY;
			ymax=ymin+DXY;
			cell[ix][iy].resize(2*NETA);
			for(ieta=0;ieta<2*NETA;ieta++){
				etamin=-ETAMAX+DETA*ieta;
				etamax=etamin+DETA;
				cell[ix][iy][ieta]=new CB3DCell(xmin,xmax,ymin,ymax,etamin,etamax);
			}
		}
	}
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				c->ix=ix; c->iy=iy; c->ieta=ieta;
				c->creflection=NULL;
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						for(jeta=ieta-1;jeta<=ieta+1;jeta++){
							if(jx>=0 && jy>=0 && jeta>=0 && jx<2*NXY && jy<2*NXY && jeta<2*NETA){
								c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=cell[jx][jy][jeta];
							}
							else c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=NULL;
						}
					}
				}
			}
		}
	}
	if(BJORKEN){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				c=cell[ix][iy][0];
				c->ireflection=-1;
				c->creflection=cell[ix][iy][2*NETA-1];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][0]=cell[jx][jy][2*NETA-1];
						}
					}
				}				
				c=cell[ix][iy][2*NETA-1];
				c->ireflection=1;
				c->creflection=cell[ix][iy][0];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][2]=cell[jx][jy][0];
						}
					}
				}
			}
		}
	}
	if(DENSWRITE){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					c=cell[ix][iy][ieta];
					c->dens.resize(DENSWRITE_NTAU);
					for(iitau=0;iitau<DENSWRITE_NTAU;iitau++){
						c->dens[iitau]=0.0;
					}
				}
			}
		}
	}
	if(BARYON_ANNIHILATION){
		imax=lrint(TAUCOLLMAX);
		annihilation_array.resize(imax);
		for(int i=0;i<imax;i++){
			annihilation_array[i]=0.0;
		}
	}
	if(MUTCALC || BARYON_ANNIHILATION){
		InitMuTCalc();
	}
}

void CB3D::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	string command="mkdir -p model_output/"+run_name+"/"+qualifier;
	system(command.c_str());
	if(sampler!=NULL)
		sampler->nevents=0;
	if(oscarfile!=NULL){
		fclose(oscarfile);
		oscarfile=NULL;
	}
	if(BFCALC){
		balancearrays->SetQualifier(qualifier);
	}
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.txt";
}

void CB3D::Reset(){
	int iitau,ntau;
	double taucalc;
	KillAllParts();
	KillAllActions();
	tau=0.0;
	nactions=0;
	//npartstot=0;
	if(MUTCALC){
		ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
		for(iitau=1;iitau<ntau;iitau++){
			taucalc=iitau*MUTCALC_DELTAU;
			AddAction_MuTCalc_UpdateNPE(taucalc);
		}
		CMuTInfo::NETEVENTS+=NSAMPLE;
	}
}

CB3D::~CB3D(){
	if(oscarfile!=NULL)
		fclose(oscarfile);
}
