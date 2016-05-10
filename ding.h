/*
 * ding.h
 *  //////// CLASS DEFINITIONS USED FOR DING'S WORK ///////////
 *  Created on: Apr 29, 2016
 *      Author: simon
 */

#ifndef DING_H_
#define DING_H_

int GetRoot(CTree * T, vector <int> RootSeqs);		// Get root branch

// Heterogeneous empirical amino acid model
class CHeteroEMP : public CBaseModel	{
public:
	CHeteroEMP(CData *Data, CTree *Tree, string Name, double *S_ij, vector <vector <int> > SubTrees, vector <int> Root);	// Constructor
	~CHeteroEMP();																						// Destructor
};

class CHeteroEmpProc : public CBaseProcess {
public:
	CHeteroEmpProc(CData *D, CTree *T, string Name, double *S_ij, vector< vector <double> > Freq, int RootFreq);		// Constructor
	~CHeteroEmpProc();																					// Destructor
	// So gamma distribution works (These functions are changed so they returns a CHeteroEmpProc model
	CBaseProcess *RateProcessCopy();								// Returns a pseudoprocess with a seperate rate category
	bool Likelihood(bool ForceReal);								// Likelihood function so matrices are correctly remade

private:
	// Functions
	void CreateEMPmodel(double *S_ij,vector <double> Freq);													// Very similar copy to CAAProcess version
	void MakeHetQMat(int ProcNum);																			// Builds the empirical models ready for use
	bool Make_PT(int Branch, bool RedoRate = false);														// Makes the P(t) matrix for the likelihood function based on the m_vpHetero* stuff below and the labelled tree
	vector <double> RootEqm();																				// Returns the root equilibrium for likelihood calculations (specified later)
	void AssignRootEqm(int root);																			// Assigns the root equilibrium
	// Variables
	vector <CBaseEqm *> m_vpHeteroEqm;																	// Vector of the heterogeneous equilibriums
	vector <CQMat *> m_vpHeteroQMat;																	// Vector of the heterogeneous Q matrices
	int m_iHeteroRootEqm;																				// Identified that states what the root equilibrium is
	int m_iRootEquilibrium;																				// The number of the process that the root belongs to
};

#endif /* DING_H_ */
