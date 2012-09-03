/*
 * ModelAssess.h
 *
 *  Created on: May 3, 2012
 *      Author: Simon Whelan
 */

#ifndef MODELASSESS_H_
#define MODELASSESS_H_
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


struct SModelDetails {
	string Name;				// Model name
	double OrilnL;				// Unadjusted model likelihood
	double lnL;					// Adjusted likelihood
	double TreeLength;			// Tree length
	double NoPar;				// Number of parameters
	double AIC;					// Normalised AIC
	EDataType DataType;			// Data type
};

SModelDetails DoModelRun(CBaseModel *M, int NoPar, double Adj = 0.0);

void GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
void GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
void GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);
void GetCODModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out = cout);

double GetAIC(double lnL, int NoPar) { return 2*(NoPar-lnL); }
#endif /* MODELASSESS_H_ */
