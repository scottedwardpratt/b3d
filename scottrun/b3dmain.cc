#include "b3d.h"
#include "sampler.h"
#include "balancearrays.h"
#include "qualifier.h"
#include "randy.h"
#include "misc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 4) {
		printf("Usage: b3d run_name ievent0 ieventf\n");
		exit(-1);
  }
  char message[500];
	long long int nparts,ninit;
	long long int nscatter,nmerge,nannihilate,nregen,nbaryons,norm,npass,nexit;
	int ievent,iqual,nevents;
	string run_name=argv[1];
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		nscatter=nmerge=nparts=ninit=nannihilate=nregen=nbaryons=npass=nexit=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		sprintf(message,"_________________ iqual=%d, nevents=%d ________________\n",iqual,nevents);
		CLog::Info(message);
		b3d->sampler->ReadHyperElements2D_OSU();
		b3d->ReadMuTInfo();
		for(ievent=ievent0;ievent<=ieventf;ievent++){
			sprintf(message,"------ beginning, ievent=%d --------\n",ievent);
			CLog::Info(message);
			b3d->Reset();
			ninit+=b3d->sampler->GenHadronsFromHyperSurface(); // Generates particles from hypersurface
			b3d->PerformAllActions();
			nscatter+=b3d->nscatter;
			nmerge+=b3d->nmerge;
			nbaryons+=b3d->nbaryons;
			nannihilate+=b3d->nannihilate;
			nregen+=b3d->nregenerate;
			npass+=b3d->npass;
			nexit+=b3d->nexit;
			nparts+=b3d->PartMap.size();
			nbaryons+=b3d->CountBaryons();
			//printf("ncheck=%lld, ncheck1=%lld, ncheck2=%lld\n",b3d->ncheck,b3d->ncheck1,b3d->ncheck2);
		}
		norm=nevents*b3d->NSAMPLE;
		sprintf(message,"<Nparts>=%8.2f, initial<Nparts>=%8.2f\n",double(nparts)/norm,double(ninit)/norm);
		sprintf(message,"%s<Nscatter>=%9.2f, <Nmerge>=%9.2f,  <NB>=%7.3f, <Nannihilate>=%7.4f, <Nregen>=%7.4f\n",
			message,double(nscatter)/norm,double(nmerge)/norm,double(nbaryons)/norm,double(nannihilate)/norm,double(nregen)/norm);
		sprintf(message,"%s<Nexit/part>=%g\n",message,double(nexit)/double(nparts));
		CLog::Info(message);
		sprintf(message,"Npass=%g\n",double(npass)/norm);
		CLog::Info(message);
		b3d->WriteMuTInfo();
	}
	return 0;
}
