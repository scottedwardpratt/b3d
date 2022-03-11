#ifndef __B3D_H__
#define __B3D_H__

#include "commondefs.h"
#include "parametermap.h"
#include "inelastic.h"
#include "mutinfo.h"
#include "decay_nbody.h"
#include "acceptance.h"
#include "resonances.h"
#include "log.h"

using namespace std;

class CB3D{
public:
	CparameterMap parmap;
	CChargeMap chargemap;
	CPartMap DeadPartMap;
	CPartMap PartMap;		//!< A C++ map for active CPart objects in the model.

	/*!
	This map is used to schedule and organize the various actions that the model must perform in time order. It contains all actions (as CAction objects) that have yet to occur, and the map's key is the boost-invariant time \f$\tau\f$ at which the action is scheduled to occur.
	*/
	CActionMap ActionMap;
	CActionMap DeadActionMap; // action objects not in line to be processed
	CResList *reslist;	//!< The CResList instance for the mucalodel (dynamically allocated).
	CInelasticList *inelasticlist;	//!< The CInelasicList instance for the model (dynamically allocated).
	
	int NXY;	//!< Determines size of mesh. The mesh size is \f$(2NXY,2NXY, 2NETA)\f$.
	int NETA;
	int NRINGSMAX;
	int NPRCELLSMAX;
	int COOPERFRYE_WMAX;
	bool COOPERFRYE_CREATENEGPARTS;
	bool USE_OLD_SAMPLER;
	bool RESONANCE_DECAYS,BFCALC;
	int HYDRO_OCTANT_SYMMETRY;
	double XYMAX,ETAMAX,DXY,DETA;
	CSampler *sampler;
	CHYDROtoB3D *hydrotob3d;
	bool BJORKEN,COLLISIONS,INELASTIC,HYDRO_PURE_BJORKEN,DENSWRITE,BARYON_ANNIHILATION,MUTCALC,SECALC;
	double ANNIHILATION_SREDUCTION;  // reduces annihilation cross section based on amount of strangeness
	int DENSWRITE_NTAU;
	int NBOSE;
	char message[200];
	double DENSWRITE_DELTAU,MUTCALC_DELTAU;
	vector<vector<vector<CB3DCell *> > > cell;
	vector<vector<vector<CMuTInfo *>>> muTinfo;
	vector<double> annihilation_array;
	CDecay_NBody *decay_nbody;
	
	void ReadCharges(int ichargefile);
	void GenHadronsFromCharges();
	void GenHadronsFromCharge(int balanceID,CCharge *charge);
	void TestChargeConservation(int pid);
	
	
	void ReadHydroInput();
	CPart *GetDeadPart();
	void GetDeadParts(CPart *&part1,CPart *&part2);
	void GetDeadParts(array<CPart*,5> &product);
	CAction *GetDeadAction();
	//vector<CLocalInfo *> localinfo;
	//
	vector<long long int> phicount; // for testing purposes
	void AnalyzeCharges(); // for testing purposes
	
	//!Constructor.
	/*!
	The constructor, which initializes all other elements of the model run.
	
	\param[in] run_name_set This is the "run name" read in from the command line.
	*/
	CB3D(); // this is a constructor which does nothing but create an object
	CB3D(string run_name_set); // this gets all arrays ready
	CBalanceArrays *balancearrays,*bfbalancearrays;
	void CopyParMapPars(); // copies parmap to B3D variables (usually similar name)
	void InitCascade();
	void InitAnalysis();
	~CB3D();	//!< Destructor.
	double tau,TAUCOLLMAX;
	int itau;
	//
	// READ IN FROM PARAMETER FILE
	int NACTIONSMAX;
	int NPARTSMAX,nbaryons,npartstot,nactionstot;
	int NSAMPLE,NSAMPLE_UDS2BAL;
	int DELNPARTSTOT,DELNACTIONSTOT;
	bool BINARY_RW;
	double SIGMAMAX,SIGMADEFAULT, SIGMABF,SIGMAINELASTIC, Q0; // cross sections in sq. fm
	double RESWIDTH_ALPHA; // sets spectral function for res widths
	string input_dataroot;
	string output_dataroot;
	string run_name,qualifier;

	string oscarfilename;
	FILE *oscarfile;
	//
	void SetQualifier(string qualifier_set);
	void MovePartsToFinalMap();
	double WriteOSCAR(int ievent);  // returns dnch/deta
	void ReadOSCARHeader();
	int ReadOSCAR(int ievent);
	double WriteBalanceParts(int ievent); // returns dnch/deta
	//void ReadBalanceParts();
	int ReadBalanceParts(int ievent);
	void WriteDens();
	void WriteAnnihilationData();
	
	void FindAllCollisions();
	void FindAllCellExits();
	void PerformAllActions();
	void Reset();
	void KillAllActions();
	void KillAllParts();
	void SplitPart(CPart *part1,CPart *part2);
	void CheckPartMap();

	void AddAction_Activate(CPart *part);
	void AddAction_Decay(CPart *part,double taudecay);
	void AddAction_Collision(CPart *part1,CPart *part2,double tau,double pibsquared,
						double sigma_scatter,double sigma_merge,double sigma_annihilation,
						double sigma_inel,vector<double> dsigma_merge);
	void AddAction_ResetCollisions(double taureset);
	//void AddAction_SwallowParticles(double tau_breakup);
	void AddAction_ExitCell(CPart *part);
	void AddAction_DensCalc(double tauwrite);
	void AddAction_MuTCalc_UpdateNPE(double taucalc);
	void AddAction_MuTCalc();
	void AddAction_SECalc(double taucalc);

	void ListFutureCollisions();
	void PrintPartList();
	void WriteMuTInfo();
	void ReadMuTInfo();
	void WriteWeights();
	void IncrementWeightArrays();
	int CountBaryons();
	void InitMuTCalc();

	bool FindCollision(CPart *part1,CPart *part2,double &taucoll);
	void Decay(CPart *mother,int &nbodies,array<CPart *,5> &daughter);
	double GetSigma(CPart *part1,CPart *part2,double Minv2,
		double &sigma_scatter,double &sigma_merge,double &sigma_annihilation,double &sigma_inel,
		vector<double> &dsigma_merge);

	CRandy *randy;

	void PrintActionMap(CActionMap *actionmap);

	bool CheckKinematics(CPart *part1,CPart *part2,double &Minv2,
		double &pibsquared,double &taucoll);
	int Collide(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product,double pibsquared); // will collide if sigma>scompare

	int Collide_Scatter(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product);
	int Collide_Merge(CPart *part1,CPart *part2,double sigma_merge,vector<double> &dsigma_merge,
		int &nproducts,array<CPart*,5> &product);
	int Collide_Annihilate(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product);
	int Collide_Inelastic(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product);


	void Scatter(CPart *part1,CPart *part2,CPart *part3,CPart *part4);
	bool Merge(CPart *part1,CPart *part2,CPart *part3,CResInfo *resinfo);
	void InelasticScatter(CPart *part1, CPart *part2,CPart *part3,CPart *part4,CInelasticInfo inelinfo);
	double GetAnnihilationSigma(CPart *part1,CPart *part2);
	bool CancelAnnihilation(CPart *part1,CPart *part2);
	int Annihilate(CPart *part1,CPart *part2,int &nproducts,array<CPart *,5> &product);

	void CheckActions();
	bool ERROR_PRINT;

	long long int nscatter,nbscatter,n,nmerge,nswallow,npass,nexit;
	long long int nactivate,nannihilate,nregenerate,nactionkills;
	long long int nactions,ninelastic, ncollisions,oldncollisions,ndecay,ncheck;

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	
	bool BALANCE_DECAY,BALANCE_CALC;
	int ibalmax;
	
	// These are used for Analysis
	CPart **partarray;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	CRegenerate *regen;
	int DecayParts(int nparts);
	double CalcSpectra_PHENIX();
	void CalcSpectra_PHENIXppbar();
	double CalcSpectra_ALICE();
	void Calc3DSpectra();
	void CalcBalance();
	void CalcBalanceQGP();
	void CalcV2_STAR();
	void CalcV2_ALICE();
	void CalcHBT_STAR();
	void CalcHBT_ALICE();
	void CalcRealityDiff();
	void CalcGamma();
	void CalcMuTU();
	double legendre(int ell,double x);
	void Consolidate(string run_name);
	void CBalanceArray(CB3D *b3d);
	void IncrementChiTotFromCharges();
	void IncrementChiTotFromHadrons();
	void DeleteCharges();
	double TotalVolume;
	Eigen::Matrix3d chitotQ,chitotH;
};

#endif
