#ifndef __B3D_ACTION_H__
#define __B3D_ACTION_H__

#include "b3d.h"
#include "part.h"
#include "cell.h"
#include "resonances.h"

//!An action in the CB3D model.
/*!
\version 1.0
\author Scott Pratt
\date March 2011

This class handles any actions that the model takes during execution. Examples of "actions" that the model takes are a resonance decaying, a particle crossing a cell boundary, a collision, new particles being generated, etc. In this way, a complex system of interacting particles is reduced to a scheduled list of actions. Scheduling is handled using a C++ map container of CAction objects, keyed by the boost-invariant time tau (\f$\tau\f$) at which they are scheduled to occur. Note that this map is revised consistently, as future actions often change dramatically as a result of the current action.

Actions are allocated and are moved from the map of future actions (CB3D::ActionMap) to the list of completed actions (CB3D::DeadActionMap) once they have been performed.
*/
class CAction{
public:
	double tau;
	double pibsquared,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel;
	vector<double> dsigma_merge;
	int listid;
	double key;
	int type; // =0 for activation, 1 for decay, 2 for collision, ....  6 for ExitCell
	// These are the particles in the action
	CPartMap partmap;

	void Kill();
	void AddPart(CPart *partptr);
	void Print();

	void Perform();
	void PerformDensCalc();
	void PerformMuTCalcUpdateNPE();
	void PerformSECalc();
	void PerformActivate();
	void PerformExitCell();
	void PerformDecay();
	void PerformCollide();
	void PerformCollide_BALANCE();
	void PerformResetCollisions();
	array <CResInfo *,5> daughterresinfo;
	//void PerformSwallowParticles();
	CAction();
	CAction(int keyset);
	~CAction();

	static CB3D *b3d;
	static char *message;

	CActionMap::iterator GetPos(CActionMap *actionmap);
	void MoveToActionMap();
	void RemoveFromActionMap();
	void AddToMap(CActionMap *newmap);
	void AddToMap(CActionMap::iterator guess,CActionMap *newmap);
	void CheckPartList();
	CActionMap *currentmap;
	array<CPart *,5> product;
};

#endif