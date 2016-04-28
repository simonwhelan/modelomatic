/* **************************************************************
 *  Code for new process stuff
 * ************************************************************** */


#ifndef NEW_PROCESS_STUFF
#define NEW_PROCESS_STUFF

#include "process.h"

///////////////////////////////////////////////////////
// Amino acid model stuff
vector <double> vWAGVal();		// Returns the 190 parameters of the WAG model
vector <double> vWAGFreq();		// Returns the 20 frequencies of the WAG model
vector <double> vJTTVal();		// Returns the 190 parameters of the JTT model
vector <double> vJTTFreq();		// Returns the 20 frequencies of the JTT model
vector <double> vMTREVVal();		// Returns the 190 parameters of the mtREV model
vector <double> vMTREVFreq();	// Returns the 20 frequencies of the mtREV model
vector <double> vDAYVal();		// Returns the 190 parameters of the Dayhoff model
vector <double> vDAYFreq();		// Returns the 20 frequencies of the Dayhoff model
vector <double> vcpVal();		// Returns the 190 parameters of the cpREV model
vector <double> vcpFreq();		// Returns the 20 frequencies of the cpREV model

/* ********************** Class defining s parameters in covarion style models ************* */
class CCSCovPar : public CQPar	{
public:
	CCSCovPar(string Name, int DataChar, int HiddenChar, double Value, bool Optimise = true, double Low = 0.0, double Up = MAX_PAR_VALUE, ParOp Oper=MULTIPLY);
	void AddCSstates(int From, int To, string sABET, int iABET = 1);
	vector <int> From()	{ return m_viCSstatesFrom; }
	vector <int> To()	{ return m_viCSstatesTo; }
private:
	int m_iHiddenChar;
	vector <int> m_viCSstatesFrom;
	vector <int> m_viCSstatesTo;
	ostream &Output(ostream &os);
};

/* ************************* Class defining substitution processes in THMMs ******************** */
class CTHMMSubProc	{
public:
	// Constructor function
	CTHMMSubProc(int ProcNum,int DataChar, double *QMat, int QChar,CPar *Rate, CBaseEqm *Eqm);
	// Destructor function
	~CTHMMSubProc();
	// Interaction functions
	void Scale(double Rate = -1.0);										// Scale the process so the mean rate is that of m_pRatePar
	void Output(ostream &os = cout, char delim = '\t', string lead = "");		// Output details of the process
	double CalcRate();													// Calculate the current rate of the process
	double Rate() { return m_pRatePar->Val(); }							// Returns the processes rate
	double Rate(double R) { if(!m_pRatePar->Special()) { m_pRatePar->SetVal(R); } return m_pRatePar->Val(); }			// Sets the processes rate
	void SetOptRate(bool Val)	{ m_pRatePar->SetOptimise(Val); }		// Sets whether the rate is optimised
	bool Opt() { return m_pRatePar->Opt(); }							// Whether the rate is optimised
	bool SpecialRate() { return m_pRatePar->Special(); }				// Returns whether the rate parameter is special
	vector <double> Eqm() { return m_pEqm->SubEqm(m_iHiddenState); }	// Returns the equilibrium distribution of the subprocess
private:
	// Private variables
	int m_iChar;					// Number of observable characters in the process
	int m_iHiddenState;				// The number of the state the process corresponds with
	vector <double *> m_vpQ;		// The QMatrix of size m_iChar*m_iChar for the substitution process
	CPar *m_pRatePar;				// The rate paramter
	CBaseEqm *m_pEqm;				// The equilibrium associated with the process
};

/* ************************* Class defining temporal Hidden Markov Models (THMM) **************************** */
// Definition of different THMM equilibrium classes
//enum THMM_EQM_TYPE {equ,obs,complex,complex_GC};
// equ			= equiprobable (e.g. Jukes and Cantor)
// obs			= single frequencies produced from data
// complex		= seperate frequency vectors for each hidden process
// complex_GC	= single frequency vector with different GC content for each process

/* ********* general THMM process that data specific ones are based on ********************** */
class CTHMMProcess : public CBaseProcess {
public:
	// Constructor
	CTHMMProcess(CData *D, CTree *T, string Name, int NoCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE Eqm = obs, bool OptFreq = false, EHIDDEN_TYPE DoHidden = H_same, ERateTypes RateTypes = same);
	// Destructor
	~CTHMMProcess();
	/* *********************** Interaction functions ***************** */
	// Control optimisation of covarion probabilities
	void OptimiseHiddenProbs(int ProcNum = -1)		{ SetHiddenProbsOpt(true,ProcNum); };	// Optimise the Hidden state probabilities
	void NoOptimiseHiddenProbs(int ProcNum = -1)	{ SetHiddenProbsOpt(false,ProcNum); };	// Do not optimise the hidden state probabilities
	void SetHiddenProbsOpt(bool Optimise,int ProcNum = -1);									// Set the optimise probabilities
	// Function for  that makes a state fixed (i.e the fixed state in the covarion process
	void EnforceFixedState(int State2Fix = 1);
	void MakeSimpleCovarionModelEqm() { int i; FOR(i,(int) m_vpEqm.size()) { m_vpEqm[i]->SetDoBasicNonReversibleCovarion(true); } }
	double SetSubProcRate(int ProcNum, double Rate);				// Set the rate of a subprocess
	double SubProcRate(int Num) { assert(InRange(Num,0,m_iHiddenChar)); return m_vpSubProcs[Num]->Rate(); }
	CTHMMSubProc *SubProc(int Num) { assert(InRange(Num,0,m_iHiddenChar)); return m_vpSubProcs[Num]; }	// Returns a pointer to a subprocess

	// Functions for rate optimisation and finding out about rates
	void SetOptRates(int ProcNum = -1,bool Val = true);				// Set whether to optimise the rates
	bool OptRate(int N) { assert(InRange(N,0,m_iHiddenChar)); return SubProc(N)->Opt(); }
	void OptimiseRates(int ProcNum = -1)	{ SetOptRates(ProcNum,true); };		// Optimise the rates of the individual substitution processes (-1 = all)
	void NoOptimiseRates(int ProcNum = -1)	{ SetOptRates(ProcNum,false); };	// Do not optimise the rates of the individual substitution processes (-1 = all)
	bool DoRates() { return m_bDoRateScaling; }						// Whether to do rate scaling
	void SetRateScaling(bool V) { m_bDoRateScaling = V; if(V == false) { NoOptimiseRates(); } }			// Sets whether to do rate scaling
	// Functions returning interesting quantities associated with the process (overriding virtual functions)
	double OverallSubRate(bool FWP = false);						// Calculate the overall substitution rate
	double OverallTransRate();										// Calculate the overall rate of switching between hidden states

protected:
	// Implementation functions
	///////////////////////////////////////////////
	// Stuff  for parameters
	CQPar *CreateZeros();										// Create the off-diagonal zeros in the covarion matrix
	CQPar *CreateChangePar(double Value,int From=-1,int To=-1,bool DoOpt = true);	// Create the parameters describing changes between states in covarion matrix
	void CreateCovProbs(vector <double> Probs);
	CQPar *AddSubProc(int ProcNum, double InitRate, bool Optimise);	// Add the details for a substitution process
	// Functions for adding Q matrices to the model
	CQMat *Add_QMat(string Name,EDataType Type);					// Function for adding simple Q matrices based on standard data types
	CQMat *Add_QMat(string Name, int Char);							// Function for adding Q matrices of size Char
	void DoPostHocQMatUpdate();										// Use to perform rate scaling of individual processes
	bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);		// New routine that does gamma rate calculations if needed
	bool PrepareQMats(bool DoScale = true)								// Overloaded function for doing all the Q matrices
				{ vector <int> v; return PrepareQMats(v,DoScale); }


	///////////////////////////////////////////////////////
	// Variables
	// ---------
	// 1. Variables holding model information
	bool m_bUseGammaRates;						// Whether rate parameters are taken from a gamma distribution
	bool m_bDoRateScaling;						// Whether to do rate scaling
	vector <CPar *> m_vpRPars;					// Parameters describing rates of change between observable characters
	vector <CCSCovPar *> m_vpSPars;					// Parameters describing rates of change between hidden processes
	vector <CTHMMSubProc *> m_vpSubProcs;		// The substitution processes in the model
	CQPar *m_Alfa;								// The alpha parameters of the gamma distribution
	int m_iNoGamCat;							// Number of categories in the gamma distribution
	// 2. Variables holding flags

};

/* ************************ DNA specific THMMs -- Only rate varying ***************************** */
// 1. General class
class CDNATHMMProcess : public CTHMMProcess {
public:
	// Constructor
	CDNATHMMProcess(CData *D, CTree *T,string Name,int NoCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE Eqm = obs, bool OptFreq = false, EHIDDEN_TYPE DoHidden = H_same, ERateTypes SeperateRates = same);
	// Destructor
	~CDNATHMMProcess();
	// Implementation
	void CreateHKYModel();
};
// 2. Class derived for the WARS model
class CDNAWARSProcess : public CDNATHMMProcess {
public:
	// Constructor
	CDNAWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv);
	// Destructor
	~CDNAWARSProcess();
	// Implementation
	void DoParameterisation();		// Applies the WARS parameters to the standardised model
	bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);	// New routine that performs DoParameterisation() in addition to usual routine
	bool PrepareQMats(bool DoScale = true)								// Overloaded function for doing all the Q matrices
				{ vector <int> v; return PrepareQMats(v,DoScale); }
protected:
	// Variables
	// 1. Parameters
	CQPar *m_pInv;					// The proportion of invariant sites
	CQPar *m_SigAlfa;				// The rate of switching within gamma distributed rates classes
	CQPar *m_SigInv;				// The rate of switching between gamma rate classes and invariant rate classes
};

// This is essentially the WARS model, but allows parameters to vary on a branch by branch basis
const double CBranchWARSProcess_PenaltyIntensity = 1.0E-1;
class CDNABranchWARSProcess : public CDNAWARSProcess {
public:
	// Constructor
	CDNABranchWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv, bool VarySig, bool VarypInv);
	// Destructor
	~CDNABranchWARSProcess();
	// Function that makes the branch specific P(t) matrices
	bool Make_PT(int B, bool RedoRate = false);
	void PrepareBraDer();			// Makes sure QP is calculated correctly
	// Output function
	ostream &Output(ostream &os);									// Overriding output function
	// The root eqm
	vector <double> RootEqm();										// Get's the root eqm from m_RootInv
	// The penalty function
	double Penalty();												// Returns the penalty function
	double PenaltyTraverse(int NodeTo, int NodeFrom);				// Does the tree traversal to calculate the penalty
	void ShuffleParameters() { Error("Haven't shuffled THMM parameters..."); }
protected:
	vector <CQPar *> m_BraSigAlfa;		// The branch wise values of SigAlfa
	vector <CQPar *> m_BraSigInv;		// The branch wise values of SigInv
	vector <CQPar *> m_BrapInv;			// The branch wise values of pInv
	CQPar *m_RootInv;					// The root value of pInv

};

/* ************************ AA specific THMMs *********************************************** */
// 1. General class
class CAATHMMProcess : public CTHMMProcess {
public:
	// Constructor
	CAATHMMProcess(CData *D, CTree *T,string Name,int NoCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE Eqm = obs, bool OptFreq = false, EHIDDEN_TYPE DoHidden = H_same, ERateTypes SeperateRates = same);
	// Destructor
	~CAATHMMProcess();
	// Implementation
	void CreateEMPModel(vector <double> S_ij, vector <double> Freq, bool AddF);
};
// 2. Class derived for the WARS model
class CWARSProcess : public CAATHMMProcess {
public:
	// Constructor
	CWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv);
	// Destructor
	~CWARSProcess();
	// Implementation
	void DoParameterisation();		// Applies the WARS parameters to the standardised model
	bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);	// New routine that performs DoParameterisation() in addition to usual routine
	bool PrepareQMats(bool DoScale = true)								// Overloaded function for doing all the Q matrices
				{ vector <int> v; return PrepareQMats(v,DoScale); }
	// Output function
	ostream &Output(ostream &os);									// Overriding output function
protected:
	// Variables
	// 1. Parameters
	CQPar *m_pInv;					// The proportion of invariant sites
	CQPar *m_SigAlfa;				// The rate of switching within gamma distributed rates classes
	CQPar *m_SigInv;				// The rate of switching between gamma rate classes and invariant rate classes
};

// This is essentially the WARS model, but allows parameters to vary on a branch by branch basis
class CBranchWARSProcess : public CWARSProcess {
public:
	// Constructor
	CBranchWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv, bool VarySig, bool VarypInv);
	// Destructor
	~CBranchWARSProcess();
	// Function that makes the branch specific P(t) matrices
	bool Make_PT(int B, bool RedoRate = false);						// Do the PT matrices
	void PrepareBraDer();											// Make sure the QP matrices are correctly calculated
	// Output function
	ostream &Output(ostream &os);									// Overriding output function
	// The root eqm
	vector <double> RootEqm();										// Get's the root eqm from m_RootInv
	// The penalty function
	double Penalty();												// Returns the penalty function
	double PenaltyTraverse(int NodeTo, int NodeFrom);				// Does the tree traversal to calculate the penalty
	void ShuffleParameters();				// Useful for mixture models where parameters combinations cause problems...
protected:
	vector <CQPar *> m_BraSigAlfa;		// The branch wise values of SigAlfa
	vector <CQPar *> m_BraSigInv;		// The branch wise values of SigInv
	vector <CQPar *> m_BrapInv;			// The branch wise values of pInv
	CQPar *m_RootInv;					// The root value of pInv

};


/* ************************ Full DNA specific THMMs ********************************************** */
class CFullDNATHMMProcess : public CTHMMProcess	{
public:
	// Constructor
	CFullDNATHMMProcess(CData *D, CTree *T,string Name,int NoCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE Eqm = obs, bool OptFreq = false, EHIDDEN_TYPE DoHidden = H_same, ERateTypes SeperateRates = same);
	// Destructor
	~CFullDNATHMMProcess();
	/* *********************** Interaction functions ***************** */
	// Function for adding parameters describing substitutions between characters to the model
	// if HiddenProc == -1 then apply to all processes
	CQPar *AddKappa(int HiddenProc = -1);													// Add a kappa parameter
	CQPar *AddDNASubPar(int S1, int S2, int HiddenProc = -1);									// Add a parameter describing changes between observable state S1 and S2
	void AddSeperateKappas();																// Add a different kappa parameter for each process
	void AddREV(int HiddenProc = -1);														// Add a REV model

	// Output function
	ostream &Output(ostream &os);									// Overriding output function
	// Functions for dealing with labelling issues
	void ShuffleParameters();				// Useful for mixture models where parameters combinations cause problems...
};

class C2StateCovProc : public CFullDNATHMMProcess {
public:
	// Constructor function
	C2StateCovProc(CData *D, CTree *T, string Name, vector <double> Probs, ETHMM_EQM_TYPE eqm_type = obs);
	// Destructor function
	~C2StateCovProc()  { /* BLANK */ };
};

/* ************************************************************************************
 * 				Coevolution models -- started 7ii11
 * ************************************************************************************
 * Model will create a model where the n*n state-space represents pairs of characters
 * The rate of a pair of characters i,k -> j,l = Q_{ik,jl}
 * This model is then produced from
 * (i) an n-state empirical rate matrix (e.g. WAG/JTT) named S
 * (ii) an n-state equilibrium distribution where the pairs of characters are independent
 * 			pi_ik = pi_i * pi_k
 * (iii) an n*n-state equilibrium where the pairs of characters are dependent
 * 		pi_ik
 *
 * ***
 *
 * The Q matrix can then be computed as
 * 		q_{ik->jl} =	{	0 								 	for i!=j && k!=l
 * 						{	S_{i->j} * ( pi_{jk} / pi_{k} ) 	for i!=j && k==l
 * 						{	S_{k->l} * { pi_{jk} / pi_{j} ) 	for i==j && k!=l
 * 						{ -row_sum(.)							for i==j && k==l
 *
 * The strength of the interaction may be measured by the relative disequilibrium between
 * 	pi_ik and pi_i:
 * 		RD = Sum_{i,k} [ pi_{ik} - (pi_{i} * pi_{k}) ]^{2}
 *  Or maybe sqrt(RD)...
 *
 *  *** */
/* ---
 * Initial Pairwise Coevolution Process
 *
 *	The joint probability pi_{ik} will be created from:
 *		a. The independent probabilities pi_{i} and pi_{k}
 *		b. A 'propensity matrix', R, describing the relative propensity for different amino acids pair to be interacting
 *		c. A variable psi, which describes the intensity of the interaction
 *				psi = 0		-> no interaction
 *				psi < 0		-> negative interaction (opposite of the propensity matrix)
 *				psi > 0		-> positive interaction (inferred from propensity matrix)
 *
 *	The formula for creating the the joint probability is pi_{ik} = pi_{i} * pi_{k} * e^(psi * r_{ik}) * c
 *	Note: the c is a scaling factor for a particular amount of coevolution and will be used to scale the pairwise process
 * ---
 */

class CPairwiseCoevoProcess: public CBaseProcess {
public:
	// Constructor / destructor
	CPairwiseCoevoProcess(CData *D, vector <double> Left, vector <double> Right, CTree *T, vector <double> * R = NULL, double init_psi = 0); // R = Propensity matrix; init_psi = Psi starting value;
	~CPairwiseCoevoProcess();
	// Equilibrium creation
	vector <double> MakeCoevoEqm(double Psi = -BIG_NUMBER);

protected:
	// Functions controlling the QMatrices
	vector <double> RootEqm();				// Created the equilibrium distribution for likelihood calculations based on Left/Right/Psi/R; Also created the correct rate
	bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);			// Prepares the Q matrices in a process
	int m_iOneChar;						// Size of state-space for a single character
	CSimpleEqm *OneCharEqmLeft, *OneCharEqmRight;	// The eqm distribution of the single character states (Left = first character; Right = second character).
	CQPar *m_pPsi;						// Intensity of coevolution parameter
	vector <double> m_vRMat;			// The 'propensity matrix'
	double m_dCoevoScale;				// The coevolution scaling factor
	vector <CQPar *> m_vpOneEqm;		// The equilibrium distribution for the one state model

};

class CAACoevoProcess : public CPairwiseCoevoProcess {
public:
	// Constructor
	CAACoevoProcess(CData *Data,vector <double> Left, vector <double> Right, CTree *Tree, vector <double> *R = NULL, double initial_psi = 0, bool DoF = true); // DoF = whether to do +F a.a. model (true).
	// Destructor
	~CAACoevoProcess();
};

#endif
