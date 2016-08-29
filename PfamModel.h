/*
 * PfamModel.h
 *
 *  Created on: May 12, 2016
 *      Author: simon
 *
 *      ---
 *
 *      This file defines the spatial heterogeneity models based on Pfam HMMs
 *
 *      The expected input into the program is an MSA and tree (details of their names below)
 */

#ifndef PFAMMODEL_H_
#define PFAMMODEL_H_

#include "optimise.h"
#include "TreeList.h"
#include "model.h"
#include "interface.h"
#include "Leaphy.h"

/////////////// DING FUNCTIONS //////////////////
vector <int> ReadNameFile(string NameFile, CData *D);	// Read a set of names and return a vector of indexes specifying those names
template <class TVecCon> vector <TVecCon> VecCon(vector <TVecCon> A, vector <TVecCon> B) {
	vector <TVecCon> Ret;
	Ret.reserve(A.size() + B.size());
	Ret.insert(Ret.end(),A.begin(),A.end());
	Ret.insert(Ret.end(),B.begin(),B.end());
	return Ret;
}
// Finds the branch in the tree that specifies the novel branch of that hypothesis
int FindSpecialSplit(CTree *T,vector <int> S1,vector <int> S2,vector <int> S3,vector <int> S4,vector <int> S5,vector <int> S6);
// Adds branch labels onto the tree such that 0 is the novel branch, and (1,4) are the subclades provides by (S1,S4);
void PaintTree(CTree *T, int KeyBra, vector <int> S1, vector <int> S2, vector <int> S3, vector <int> S4);
void PaintTree(CTree *T, vector <vector <int> > Subtrees);


void PfamModelAnalysis(string RunFile);

class CPfamModel : public CBaseModel	{
public:
	CPfamModel(CData *Data, CTree *Tree, string ProbFileName, bool DoExt = false);		// Data/Tree/File name of frequencies
	~CPfamModel();
	// Variables
	vector <vector <double> > m_vvdSiteFreq;						// The storage of the sitewise probability vectors from the ProbFileName (derived from pfam)
	CQPar *m_pFreqFactor;			// The frequency factor	(used in experimental Pfam model
	// Implementation
	double lnL(bool ForceReal = false);	// perform a likelihood calculation (if over-ridden, be careful other likelihood functions are too e.g. DoBralnL(...) )
	vector <double> GetDerivatives(double CurlnL = -BIG_NUMBER, bool *OK = NULL);		// Calculate the processes derivatives
	void PrepareFastCalc() { return; };
	// Fast branch calculations
	double FastBranchOpt(double CurlnL, double Tol = 1.0E-7, bool *Conv = NULL, int NoIter = 5, bool CheckPars = true);	// Controller function
	void BranchOpt(int First,int NTo, int NFr, double *BestlnL,double tol);
	// Space update functions for model (used in stepwise addition routines
	void Leaf_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoFullUpdate = false);
	void Bran_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoFullUpdate = false);
	void PreparePT(int Br) { } ;												// This needs to be blank cos everything is done within each loop
};

class CPfamProcess : public CBaseProcess {
public:
	CPfamProcess(CData *Data, CTree *Tree, vector <vector <double> > *Freqs);		// Data/Tree/Pointer to a permanent store of the frequencies (to be held in the model above)
	~CPfamProcess();
	CBaseProcess *RateProcessCopy();								// Returns a pseudoprocess with a seperate rate category
	// Variables
	vector <CQMat *> m_vpSiteQ;	// Array of Q matrices for the process (better because then it doesn't involve rebuilding the matrix each time)

	// Implementation
	bool Likelihood(bool ForceReal = false);		// Main likelihood function
	virtual void MakeZeroRateLikelihood();				// Makes likelihood for Rate = 0 models
	virtual void MakeMaxRateLikelihood();				// Makes likelihood for Rate = inf model (garbage collector models)
	void PartialL(CTree *pTree, int iNoTo, int iNoFr, int Branch, bool LeafFirst = true, bool BlockWriteback = false);
	void BranchNodePartialL(int SpFrom, int SpTo, double *PT, bool First);
	void LeafNodePartialL(CTree *pTree, int LeafNode, int Branch, int SpFlag, double *PT,bool First);		// Usual leaf node calculation
	void TransScale(int NodeTo,int NodeFr,bool First, bool Partial = true, bool ForceReal = false);	// Transfer (sum) values from NodeFr to NodeTo
	void DoScale(int Node, bool ForceRealScale = false);					// Performs scaling on the fly; ForceRealScale = true will treat all nodes as real
	void CleanScale(int NodeNum, bool ForceAll = true);						// Removes all scales in a node
	void PrepareBraDer();						// Prepare the process for get branch derivatives function
	// Branch derivatives
	bool GetBraDer(CProb *ModelL);                          // Gets the branch derivatives. ModelL = array of partial likelihoods for whole model
	void Branch_dT(int NTo, int NFr, int Branch, CTree *pTree, int First, int *BrError);
	void LeafNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	void BranNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	double PartialGrad(int site,double Total,int SiteScale);
	void FullLeafNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoCompleteUpdate = false);			// Loop function going through sites
	void FullBranNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoCompleteUpdate = false);	// Loop function going through sites
	void LeafNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoCompleteUpdate = false);		// Sitewise
	void BranNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoCompleteUpdate = false); // Sitewise
	// Fast branch calculation stuff
	void GetBranchPartL(CProb **P, int NTo, int NFr, int Br);	// Get the sitewise likelihoods

	// Variables
	int m_iCurSite;		// Used in partial likelihood computations
};


////////////////////////////////////////// Models with an additional parameter that fits the frequencies to the data //////////////////////////////////////////////////////
class CExpPfamModel : public CPfamModel {
public:
	CExpPfamModel(CData *Data, CTree *Tree, string ProbFileName, double StartE = 1.0);
	~CExpPfamModel() { /* BLANK */ };
	vector <vector <double> > m_vvdSiteFreq;	// Store of original pfam frequencies
	// Functions
	vector <double *> GetOptPar(bool ExtBranch = true, bool IntBranch = true, bool Parameters = true, bool Eqm = false);
};

class CExpPfamProcess : public CPfamProcess {
public:
	CExpPfamProcess(CData *Data, CTree *Tree, vector <vector <double> >*Freqs, double StartE = 1.0);
	~CExpPfamProcess();
	CBaseProcess *RateProcessCopy();								// Returns a pseudoprocess with a seperate rate category
	// Likelihood function that updates the frequencies as well
	bool Likelihood(bool ForceReal = false);
	// Variables
	CQPar *m_pFreqFactor;
	vector <vector <double> > m_vvdOriSiteFreq;
	vector <double> m_vdBackgroundFreq;
	vector <vector <double> > m_vvdFreqFact;
	// The frequencies used in the model are m_vdBackgroundFreq[i] * exp(m_vvdFreqFact[site][i] * m_pFreqFactor->Val());
	//  When m_pFreqFactor == 0 -> Background frequencies
	//  When m_pFreqFactor == 1 -> Pfam frequencies

};



#endif /* PFAMMODEL_H_ */
