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
	long long int ncolls,nannihilate,nregen,nbaryons,norm,npass;
	int ievent,iqual,nevents;
	string run_name=argv[1];
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		ncolls=nparts=ninit=nannihilate=nregen=nbaryons=npass=0;
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
			ncolls+=b3d->nscatter+b3d->nmerge;
			nbaryons+=b3d->nbaryons;
			nannihilate+=b3d->nannihilate;
			nregen+=b3d->nregenerate;
			npass+=b3d->npass;
			nparts+=b3d->PartMap.size();
			nbaryons+=b3d->CountBaryons();
		}
		norm=nevents*b3d->NSAMPLE;
		sprintf(message,"<Nparts>=%8.2f, initial<Nparts>=%8.2f\n",double(nparts)/norm,double(ninit)/norm);
		sprintf(message,"%s<Ncolls>=%9.2f, <NB>=%7.3f, <Nannihilate>=%7.4f, <Nregen>=%7.4f\n",
			message,double(ncolls)/norm,double(nbaryons)/norm,double(nannihilate)/norm,double(nregen)/norm);
		CLog::Info(message);
		sprintf(message,"Npass=%g\n",double(npass)/norm);
		CLog::Info(message);
		b3d->WriteMuTInfo();
	}
	return 0;
}
