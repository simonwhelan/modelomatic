/*
 * modelomatic.h
 *
 *  Created on: May 3, 2012
 *      Author: Simon Whelan
 */

#ifndef MODELOMATIC_H_
#define MODELOMATIC_H_
#ifdef ISUNIX
#include "optimise.h"
#include "TreeList.h"
#include "model.h"
#include "interface.h"
#else
#include "../PhyloLib/optimise.h"
#include "../PhyloLib/TreeList.h"
#include "../PhyloLib/model.h"
#include "../PhyloLib/interface.h"
#endif

enum Lcorrection { L_EQU, L_EMP, L_NA };	// Uses likelihood correction from equal unosbserved character frequencies; empirical unobserved character frequencies; not applicable (e.g. nucleotides and codons).

struct SModelDetails {
	string Name;				// Model name
	double OrilnL;				// Unadjusted model likelihood
	double lnL;					// Adjusted likelihood
	double TreeLength;			// Tree length
	double NoPar;				// Number of parameters
	double AIC;					// Normalised AIC
	Lcorrection Correction;				// Likelihood correction used
	EDataType DataType;			// Data type
};

SModelDetails DoModelRun(CBaseModel *M, int NoPar, Lcorrection, double Adj = 0.0);

int GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
int GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
int GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
int GetCODModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
bool GetModels(string file = "modelomatic.ini");

double GetAIC(double lnL, int NoPar) { return 2*(NoPar-lnL); }
#endif /* MODELOMATIC_H_ */
