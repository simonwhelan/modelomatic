////////////////////////////////////////////////////////////////
// Header file for models
// ----------------------
// Models contain a number of processes
//
// The model is responsible for


#ifndef MODEL_HEADER
#define MODEL_HEADER


#include "process.h"
#include "new_process.h"
#include "parsimony.h"

struct SBestModel	{
	double m_dlnL;
	vector <double> m_dParVal;
	vector <CTree> m_vTree;
	int m_iOptObs;
};

enum EParType {branch_par,model_par};				// Defines the type of optimisation a parameter undergoes
enum ECalcType { cML,cMP,cRMSD };					// Defines the type of calculations to be performed
class CBaseModel	{
public:
	CBaseModel(CData *D,CTree *T,string Name = "Unnamed model");
	~CBaseModel();
	// Variables
	vector <CBaseProcess *> m_vpProc;		// The list of proceses
	vector <double>	m_vdProbProc;			// The probability of the individual processes
	vector <CPar *> m_vpPar;				// The parameters in the processes (not to be deleted. That is done by the processes themselves
	vector <CPar *> m_vpAllOptPar;			// All the optimised parameters including branches
	CProb *m_arL;							// The vector of sitewise likleihoods for the model (used for calculating GetLogLikelihood
	CData *m_pData;							// Pointer to the data
	CTree *m_pTree;							// Pointer to the current tree
	CTree *m_pSubTree;						// Pointer to current subtree
	string m_sName;							// Name of the model

	////////////////////////////////////////////////////////////
	// Assign data to a model
	void MakeNewData(CData *Data, bool Overwrite = false);

	////////////////////////////////////////////////////////////
	// Some access functions
	int NoSeq() { return m_pData->m_iNoSeq; }	// The number of sequences
	int NoPar()	{ return (int)m_vpPar.size(); }		// The Number of parameters
	void GetParVec(vector <CPar*> *Par);		// Get a vector of parameters
	void RemovePar(string Name);				// Removes a parameter
	void RedoScale(bool Force = false);			// Rescale all the parameters in the model
	void RandomiseParameters(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm);	// Randomise the optimised parameter values
	bool IsViable();							// Checks whether all the memory and
	////////////////////////////////////////////////////////////
	// Some calculation based functions
	void CreateProcessSpace(bool force = false);	// Function that creates space for all the sub functions
	void FinalInitialisation();			// Function that does any extra initialisations required before calculations
	void PrepareFastCalc();				// Prepare for fast computation
	void CleanFastCalc(bool Force = false);				// Clean the fast calc stuff
	void PreparelnL(bool ForceRemake = false);			// Prepares the Q matrices in the model and sets the rate
	virtual double lnL(bool ForceReal = false);	// perform a likelihood calculation (if over-ridden, be careful other likelihood functions are too e.g. DoBralnL(...) )
	double CalculateL(CBaseProcess *Process = NULL, bool DoFulllnL = true);			// Calculate the overall log-likelihood from a single process
	vector <double> SitewiseL(CBaseProcess *Process = NULL, bool DoFulllnL = true);	// Calculate the per site log-likelihood from a single process
	vector <double> BranchDerivatives();
	SBestModel ModelScores(double BestlnL);				// Outputs the best model
	double RestoreBestModel(SBestModel *M);	// Restores the best model
	vector <double> PartialPW() { assert(IsSubTree()); return m_vdPartialTreeDist; }
	bool ForceSeperateParOpt() { return m_bDoSepParOpt; };	// Difficult models insist that they to be optimised more carefully...
	int OptNum(int ON = -1) { if(ON > 0) { m_iOptNum = ON; } return m_iOptNum; }
	// Clever function for fast optimisation ofbranch lengths
	double FastBranchOpt(double CurlnL, double Tol = 1.0E-7, bool *Conv = NULL, int NoIter = 5, bool CheckPars = true);	// Controller function
	void SingleBranchOpt(int Br, double *BestlnL, double tol);				// Do a branch optimise for a single branch, calculates partial likelihoods so inefficient when traversing a tree
	// Space update functions for model (used in stepwise addition routines
	void Leaf_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoFullUpdate = false);
	void Bran_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoFullUpdate = false);

	////////////////////////////////////////////////////////////
	// Functions required for tree rearrangement
	void MakeNewTree(CTree *T, bool Overwrite);				// Replaces this tree
	vector <int> PrepareBranchCP(int BranchCP, int Depth);	// Calculates partial likelihoods from a BRANCH centre point and puts them in forward and backward space. Returns the nodes for which the partial likelihoods prepared
	vector <int> PrepareNodeCP(int NodeCP, int Depth);		// Calculates partial likelihoods from a NODE centre pointand puts them in forward and backward space. Returns the nodes for which the partial likelihoods prepared
	int NumLeafCP();										// Returns the number of the CP has produced
	void ApplySubTree(CTree *Tree, bool UseExtBra = true, bool Overwrite = false);	// Once correct leaf mapping is done, this will apply a new sub-tree to all processes FOR CALCULATIONS ONLY. Will overwrite the tree if Overwrite == true
	void FixSmallBranches();				// Resets zero branches in trees to a small number
	bool CheckSameTree();					// Checks whether all processes are using the same tree
	void CleanCPMapping();					// Goes through all the processes and removes the centre-point mapping
	void BuildOriSubTree(CTree *);			// Construct the original subtree and return it;
	bool IsSubTreeSame(CTree *SubTree);		// Compares subtree with original tree to see if they're the same (memory check only)
	bool IsSubTree();						// Whether any processes are working from subtrees
	CTree *Tree() { if(IsSubTree()) { return m_pSubTree; } return m_pTree; }
	vector <int> LeafMap()		{ return m_viLeafMap; }
	void WriteLeafMap(vector <int> NewLeafMap);		// Write LeafMap to this and associated models
	vector <double> ExtBra()	{ return m_vdExtBra; }
	vector <int> NodeFr()		{ return m_viNodeFrom; }
	vector <int> NodesCovered()	{ return m_viCPNodesCovered; }
	vector <int> NodesBelow(int N) { return m_vviNodesBelow[N]; }
	vector <vector <int> > NodesBelow() { return m_vviNodesBelow; }
	void PreparePT(int Br);												// Prepare PT matrix for a specific branch
	// Pairwise preperation
	void PreparePairwiseCalc(int Seq1, int Seq2, CTree *T);				// Prepare for a pairwise calculation
	// Triplet calc preperation
	void PrepareTripletCalc(vector <int> LeafMap, CTree *T,int SPR);	// Prepare for a triplet calculations; assumes space in interior nodes is right
	// SPR preperation
	int PrepareSPR(int Branch, int LinkNum, int *OriBr, LazyType DoLazy = lazy);		// Prepares model for SPR calculations; returns spare node for computations

	// Process handling functions
	void MakeGammaModel(int ProcessNum, int NoCat, double InitAlpha = INITIAL_GAMMA);	// Takes m_vpProc[ProcessNum] and turns it into a gamma distributed rate model
	void MakeInvariantSitesModel(int ProcessNum);				// Adds an invariant sites class to the model
	void MakeGarbageCollectorModel(int ProcessNum);				// Adds a site class to the model with infinite rate (serves to collect poorly characterised columns).
	void LockModel() { m_bLockModel = true; }					// Lock all parameters in a model so they're not optimised
	void UnlockModel() { m_bLockModel = false; }				// Allows parameters in a model to be optimised
	bool Locked() { return m_bLockModel; }						// Returns whether the model is locked or not
	virtual void CompressModel() { };							// If available, compresses the model for tree estimation
	virtual void UncompressModel() { };							// Uncompresses the model
	// Functions dealing with probabilities and other bits and pieces of different processes
	void CleanProcProbs();										// Remove process probabilities from the parameters held in m_vpPar
	void PrepareProcessProbs(bool OptProbs = true);				// Add process probabilities to the optimised parameters
	vector <CPar *> GetProcessProbs();							// Get the distribution of probabilities for the different processes
	double ModelRate(double R = -1.0);							// Returns the models rate (or sets it to R if R >= 0.0)
	vector <double> Rates();									// Returns the rate distribution across processes

	// Functions associated with TreeHMMs
	void MakeTreeHMM();											// Makes a TreeHMM from the current data

	// Output detail and functions
	string Name() { return m_sName; }
	bool OutputDetail() { return m_bOutputDetail; }
	void ComplexOutput() { m_bOutputDetail = true; }
	void SimpleOutput() { m_bOutputDetail = false; }
	vector <string> ParNameOut();								// Output the name of parameters
	vector <double> ParValOut();								// Output the values of parameters
	bool MainModel() { return m_bMainModel; }					// Flag as to whether this is a main model or not
	void FastBranchDetails(ostream &os = cout) { os << "\nFastBranchOptimiser details: TotalCalls = " << m_iFastBralnL_Calls << "; TotallnL = " << m_iFastBralnL << "; BracketlnL = " << m_iFastBralnL_Bracket; if(m_iFastBralnL_Bracket!=0) { cout << "; ratio = " << m_iFastBralnL/m_iFastBralnL_Bracket; } }

	////////////////////////////////////////////////////////////
	// Some space allocation functions
	void NewDataTree(CData *D, CTree *T);
	void ZeroSpace();						// Sets all the space in the model to 0.0
	CDNAProcess *AddDNAProcess(CData *Data,CTree *Tree, DNAProc Model, string name = "");
	CAAProcess  *AddAAProcess(CData *Data, CTree *Tree, AAProc Model, bool AddF);
	CCodonProcess *AddCodonProcess(CData *D, CTree *T, CodonProc Model, ECodonEqm CE, int GenCode);
	string AddFullDNATHMMProcess(int NoProc, EHIDDEN_TYPE DoHidden, ETHMM_EQM_TYPE DoFreqs, bool DoKappa, ERateTypes DoRates,vector <double> Rates);	// Returns a suitable name for the model
	string AddFullDNATHMMProcess(int NoProc, EHIDDEN_TYPE DoHidden, ETHMM_EQM_TYPE DoFreqs, bool DoKappa, ERateTypes DoRates);							// Returns a suitable name for the model
	string AddAATHMMProcess(double Alfa, int CatGam, double SigAlpha, double pInv, double SigInv, bool VarySig, bool VarypInv); // Returns the name of a THMM_AA
	string AddDNATHMMProcess(double Alfa, int CatGam, double SigAlfa, double pInv, double SigInv, bool VarySig,bool VarypInv);
	string AddAACoevoProcess(CData *JointData, vector <double> Left, vector <double> Right, vector <double> *PropMat, double init_psi);		// Returns the name of the coevolution process
	CBaseProcess *AddCoevoProcess(CData *Data, CTree *Tree, vector <double> *R =NULL, int init_psi = 0);
	vector <CPar *> CreateOptPar();			// Transfers the optimised parameters to m_vpPar;
	// Optimisation interaction functions
	vector <double *> GetOptPar(bool ExtBranch = true, bool IntBranch = true, bool Parameters = true, bool Eqm = false);
	int CountOptPar(bool ExtBranch = true, bool InBranch = true, bool Parameters = true, bool Eqm = false);
	vector <double> GetDerivatives(double CurlnL = -BIG_NUMBER, bool *OK = NULL);		// Calculate the processes derivatives
	double GetNumDerivative(double *Par, double lnL);
	virtual CBaseModel *PreOptModel()	{ return NULL; }								// Returns the 'preoptimisation' model for complex models (e.g. THMMs)
	virtual void ApplyPreOptModel(CBaseModel *PreOpt)		{  }						// Applies the pre-opt model to the data
	bool ReplaceParValue(CPar *Par);				// Searches for *Par in current model and replaces it (if it exists)
	// Space access functions
	void OutSpace(ostream &os,int Proc, int Node, int PBeg = -1, int PEnd = -1);
	void OutPT(ostream &os, int Proc = -1, int Branch = -1);
	void SetRandomSeed(long v)	{ m_lOriRandomSeed = v; PlantSeeds(v);}
	long GetOriRandomSeed()		{ return m_lOriRandomSeed; }
	long GetCurrentRandomSeed() { long v; GetSeed(&v); return v; }
	// Controllers for doing RMSD and MP calculations instead of likelihoods
	ECalcType DoRMSDCalc() { if(IsRMSDCalc()) { return m_CalcType; } m_CalcType = cRMSD; return m_CalcType; }		// Do RMSD calcs instead of proper likelihoods
	ECalcType DoParsimony(double Limit = PROP_SITES_PARS);					// Do parsimony calculations
	ECalcType DoLikelihoodCalc() { if(IsLikelihoodCalc()) { return m_CalcType; } m_CalcType = cML; return m_CalcType; }	// Do full likelihood computation for likelihoods
	bool IsRMSDCalc()		{ if(m_CalcType == cRMSD) { return true; } return false; }
	bool IsLikelihoodCalc()	{ if(m_CalcType == cML) { return true; } return false; }
	bool IsParsimony()	{ if(m_CalcType == cMP) { return true; } return false; }
	void ResetParsimony(double Limit = PROP_SITES_PARS) { m_Pars->CreateLimitMask((int)((double)m_pData->m_iNoSeq * Limit)); }
	double GetFullParsimony();
	// Compression values
	bool Compressed() { return m_bCompressedModel; }
	void SetCompressed()	{ m_bCompressedModel = true; }
	void SetUncompressed()	{ m_bCompressedModel = false; }
	/////////////////////////////////////////////////////////////////
	// Functions for doing simulation
	void DoSimulation(int NoSeq, int Size,vector <vector <int> > *Seqs,bool KeepAnc);
	void EvolveSequences(int NodeFrom, int NodeTo, vector <vector <int> > *Seqs, vector <int> *Procs);

	////////////////////////////////////////////////////////////////
	// Functions for obtaining partial likelihood information
	//  for robust counting procedures
	void PartialLOutput(string File, int DoBranch = -1);	// File = output file; DoBranch = output branch (-1 = all)

	////////////////////////////////////////////////////////////////
	// General function for altering likelihood calculations
	double (*pLikelihood)(CBaseModel *M);			// Function that can be used to adjust the likelihood
protected:
	bool m_bDoSepParOpt;							// Whether the model insists that it needs extra optimisations
	// Functions relating to checking and/or optimising tree branches
	void BranchOpt(int First,int NTo, int NFr, double *BestlnL,double tol);	// Recursive function
	virtual void DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol,bool AllowUpdate = true);
			// Do the actual optimisation. (NB: virtualised so that other model/branch specific parameters can be optimsised too)
	virtual double DoBralnL(int B, int NL,int NR);						// Do calculations for a branch
	// Functions relating to partial likelihood output
	void RecPartLOut(int First,int NTo, int NFr, ostream &out, int Branch);
	void DoPartLOut(int NTo, int NFr, int Br, ostream &out,int Branch);
	// Models associated with this model
	vector <CBaseModel *> m_vpAssociatedModels;		// These models are associated with the current model. Only one can be main model
private:
	//////////////////////////// Private variables ////////////////////////////////////
	// Models associated with this model
	bool m_bMainModel;								// Flag as to whether this is a main model or not
	bool m_bTreeHMM;								// Flag as to whether model is a tree HMM
	// Model associated parameters (Only basic parameters stored here)
	CPar *m_ModelRate;								// Rate of the model
	// The random seed
	long m_lOriRandomSeed;
	// Variables for optimisation algorithms
	bool m_bCompressedModel;						// Whether the model has been compressed
	ECalcType m_CalcType;							// The type of calculation to be performed
	bool m_bOptReady;								// Whether the model is ready for optimisation
	bool m_bLockModel;								// Whether parameter values will ever be optimised
	vector <bool>		m_vbDoBranchDer;			// Whether processes require branch derivative calculations
	// Variables relating to centre point mapping
	vector <int> m_viCPNodesCovered;				// Internal nodes created by the centre point
	vector <int> m_viLeafMap;						// Leaf nodes created by the centre point
	vector <vector <int> > m_vviNodesBelow;			// Lists of nodes below each of the nodes in m_viLeafMap
	vector <int> m_viNodeFrom;						// The nodes the LeafNodes came from
	vector <double> m_vdExtBra;						// External branch lengths (leading to leaf nodes) created by CP
	// Variables relating to pairwise distances
	vector <double> m_vdDist;						// Likelihood distances
	vector <double> m_vdPartialTreeDist;			// Partial tree distances for subtrees
	// Variables relating to parsimony calculations...
	CParsimony *m_Pars;								// The parsimony object
	// Other bits and pieces
	bool m_bOutputDetail;							// Amount of output detail
	// Preoptimiser stuff
	EModel m_PreOptModel;							// Model
	// Optimiser information
	int m_iOptNum;									// Number of times the optimiser should be run
	int m_iFastBralnL;								// Number of times FastBranchlnL is run;
	int m_iFastBralnL_Bracket;						// Number of those runs that are associated with bracketing
	int m_iFastBralnL_Calls;						// Number of DoBraOpt calls
	//////////////////////////// Private functions ////////////////////////////////////
	void CleanPar();						// Clean the parameter vector
	void CleanMemory();						// Clean the memory
	// Calculation functions
	bool FormMixtureSitewiseL();													// Calculate likelihoods as though they are from a mixture model
	// Functions relating to centrepointing and partial likelihood computations
	void ApplyCPMapping(vector <int> LeafMapping, vector <int> NodeFrom);	// Applies the leaf mapping to all processes
	// Output functions
	friend ostream &operator<<(ostream &os, CBaseModel &Model);
};

ostream &operator<<(ostream &os, CBaseModel &Model);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General function that returns a pointer to class CBaseModel with the model you want
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CBaseModel * GetMyModel(EModel ModelChoice, CData *Data, CTree *Tree);

////////////////////////////////////////////////////////////////
//			DNA models.
////////////////////////////////////////////////////////////////

// The model definitions
class CJC : public CBaseModel	{
public:
	CJC(CData *Data, CTree *Tree);
};

class CFEL : public CBaseModel	{
public:
	CFEL(CData *Data, CTree *Tree);
};

class CK2P : public CBaseModel	{
public:
	CK2P(CData *Data, CTree *Tree);
};
class CHKY : public CBaseModel	{
public:
	CHKY(CData *Data,CTree *Tree);
};

class C2RateHKY : public CBaseModel {
public:
	C2RateHKY(CData *Data, CTree *Tree);
};

class C2RateJC : public CBaseModel {
public:
	C2RateJC(CData *Data,CTree *Tree);
};
class C3RateFEL : public CBaseModel {
public:
	C3RateFEL(CData *Data, CTree *Tree);
};

class CREV : public CBaseModel	{
public:
	CREV(CData *Data,CTree *Tree);
	void CompressModel();
	void UncompressModel();
};


////////////////////////////////////////////////////////////////////////////
// The RY model
class CRY : public CBaseModel {
public:
	CRY(CData *D, CTree *T);
};

//////////////////////////////////////////////////////////////
// THMM DNA models with everything varying across evolution
// ---
//	See Whelan. 2008. MBE.
// ---
// There are four factors
// H = HiddenChanges:
//		H_none = no transition between states (mixture model);
//		H_same  = same rate of transition for all processes;
//		H_diff = different rates between all processes)
// F = Nucleotide frequencies:
//		equ = equiprobable model used in all frequencies
//		simple = one set of frequencies for all processes
//		complex = different sets of frequencies in all processes;
// K = Kappa:
//		false = one Kappa for all processes; true = different Kappa for all processes;
// R = Rates:
//		false = same rate for all processes; true = different rate for each process;
class CTHMM_FULL : public CBaseModel {
public:
	// Constructor defaults to the most complex case
	CTHMM_FULL(CData *D, CTree *T, int NoProc, ERateTypes DoRates = varyall, vector <double> *InitRates = NULL, ETHMM_EQM_TYPE DoFreq = complex, bool DoKappa = true, EHIDDEN_TYPE DoHidden = H_diff);
	// Implementation
	void CompressModel();
	void UncompressModel();
	CBaseModel *PreOptModel();								// Returns the 'preoptimisation' model for complex models. This will be HKY
	void ApplyPreOptModel(CBaseModel *PreOpt);								// Applies the pre-opt model to the data to get starting parameters
private:
	// Internal variable
	bool m_bRates;
};

// A general THMM for DNA data based around the HKY and +dG models
// ---
// The model has 4 rate categories, linked together by a THMM
// Essentially the model of Galtier (if I remember correctly)
class CHKYdG_THMM : public CTHMM_FULL {
public:
	// Constructor function
	CHKYdG_THMM(CData *D, CTree *T, int NoCat);
	// Destructor function
	~CHKYdG_THMM() { };

};

//////////////////////////////////////////////////////////////
// THMM DNA models (only rate varying)
// ---
// Uses the HKY model as a base
// OriAlpha		= Starting values of alpha in gamma distribution (-1 = no gamma)
// NoCatAlpha	= Number of categories in gamma distribution (1 = no gamma)
// ProbInv		= Probability of invariant site class ( < 0 -> no invariant sites)
// SigAlpha		= Rate of change between
//
class CTHMMDNA: public CBaseModel {
public:
	// Constructor function
	CTHMMDNA(CData *D, CTree *T, double OriAlpha = 0.5, int NoCatAlpha = 4, double ProbInv = 0.15, double SigAlpha = 0.1, double SigInv = 0.15);
	// Destructor function
	~CTHMMDNA() { /* empty */ };
	// Implementation
};
// As above but allows things to vary over branches
#define ALLOW_BDNA_PEN 0 // Whether to incorporate a penalty function into the BAA code
class CTHMMBDNA : public CBaseModel {
public:
	// Constructor function
	CTHMMBDNA(CData *D, CTree *T, double OriAlpha = 0.5, int NoCatAlpha = 4, double ProbInv = 0.15, double SigAlpha = 0.1, double SigInv = 0.15, bool VarySig = true, bool VarypInv = true);
	// Destructor function
	~CTHMMBDNA() { /* empty */ };
	// Implementation
#if ALLOW_BDNA_PEN == 1
	double DoBralnL(int B, int NL,int NR);						// Do calculations for a branch
	double lnL(bool ForceReal = false);	// perform a likelihood calculation
#endif
};


//////////////////////////////////////////////////////////////
// THMM AA models
// ----------
// MatType		= WAG/JTT/mtREV		Decides what replacement matrix is being used for a particular model
// OriAlpha		= Starting values of alpha in gamma distribution (-1 = no gamma)
// NoCatAlpha	= Number of categories in gamma distribution (1 = no gamma)
// ProbInv		= Probability of invariant site class ( < 0 -> no invariant sites)
// SigAlpha		= Rate of change between
//
class CTHMMAA : public CBaseModel {
public:
	// Constructor function
	CTHMMAA(CData *D, CTree *T, int MatType, double OriAlpha = 0.5, int NoCatAlpha = 4, double ProbInv = 0.15, double SigAlpha = 0.1, double SigInv = 0.15);
	// Destructor function
	~CTHMMAA() { /* empty */ };
	// Implementation
	int m_MatType, m_NoCatAlpha;
	double m_OriAlpha, m_ProbInv, m_SigAlpha, m_SigInv;
	CBaseModel* PreOptModel();					// Get initial parameter values
	void ApplyPreOptModel(CBaseModel *PreOpt);	// Apply the preoptimised model values to the current model
};

// As above but allows things to vary over branches
#define ALLOW_BAA_PEN 0 // Whether to incorporate a penalty function into the BAA code
class CTHMMBAA : public CBaseModel {
public:
	// Constructor function
	CTHMMBAA(CData *D, CTree *T, int MatType, double OriAlpha = 0.5, int NoCatAlpha = 4, double ProbInv = 0.15, double SigAlpha = 0.1, double SigInv = 0.15, bool VarySig = true, bool VarypInv = true);
	// Destructor function
	~CTHMMBAA() { /* empty */ };
	// Implementation
	int m_MatType, m_NoCatAlpha;
	double m_OriAlpha, m_ProbInv, m_SigAlpha, m_SigInv;
	CBaseModel* PreOptModel();					// Get initial parameter values
	void ApplyPreOptModel(CBaseModel *PreOpt);	// Apply the preoptimised model values to the current model
	double DoBralnL(int B, int NL,int NR);						// Do calculations for a branch
#if ALLOW_BAA_PEN == 1
	double lnL(bool ForceReal = false);	// perform a likelihood calculation
#endif
protected:
	bool m_bAllowParBraOpt;			// Whether or not to allow a sigma to be estimated on a branch-by-branch basis
	void DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol,bool AllowUpdate = true); // Branch specific optimiser for sigmas
};

// A standard THMM using 4 discrete gamma distributed rates that change through time
class CWAGdG_THMM : public CTHMMAA {
public:
	// Constructor function
	CWAGdG_THMM(CData *D, CTree *T, int MatType, double OriAlpha = 0.5, int NoCatAlpha = 4);
	// Destructor function
	~CWAGdG_THMM();
	// Implementation
};

// Exon detector
class CExonDetector : public CBaseModel {
public:
	CExonDetector(CData *D, CTree *T, int NoProc,double LowRange = 0.0, double UpRange = 1.0,vector <double> *Probs = NULL, bool OptProbs = true);
};

// DNA covarion models
class CClassicCovJC : public CBaseModel {
public:
	CClassicCovJC(CData *D, CTree *T,double OriProb = 0.75);
};
class CCovJC : public CBaseModel {
public:
	CCovJC(CData *D, CTree *T,double OriProb = 0.75);
};
class CClassicCovHKY : public CBaseModel {
public:
	CClassicCovHKY(CData *D, CTree *T, double OriProb = 0.75);
};
class CCovHKY : public CBaseModel {
public:
	CCovHKY(CData *D, CTree *T, double OriProb = 0.75);
};
class CCovREV : public CBaseModel {
public:
	CCovREV(CData *D, CTree *T, double OriProb = 0.75);
};

/////////////////////////////////////////////////////////////////
//		Amino acid models
/////////////////////////////////////////////////////////////////
class CEQU : public CBaseModel	{
public:
	CEQU(CData *Data, CTree *Tree, bool AddF = true);
};

class CWAG : public CBaseModel	{
public:
	CWAG(CData *Data, CTree *Tree, bool AddF = true);
};

class CJTT : public CBaseModel	{
public:
	CJTT(CData *Data, CTree *Tree, bool AddF = true);
};

class CDAY : public CBaseModel	{
public:
	CDAY(CData *Data, CTree *Tree, bool AddF = true);
};

class CMTREV: public CBaseModel	{
public:
	CMTREV(CData *Data, CTree *Tree, bool AddF = true);
};

// Creates a generic Q matrix
class CEMP : public CBaseModel {
public:
	CEMP(CData *Data, CTree *Tree, string Name, bool AddF, double *S_ij, double *pi_j);
};

//////////////////////////////////////////////////////////////////
//		Codon models
//////////////////////////////////////////////////////////////////

class CCodonM0 : public CBaseModel {
public:
	CCodonM0(CData *Data, CTree *Tree, ECodonEqm CE = F3X4, int GenCode = 0);
};

//////////////////////////////////////////////////////////////////
//		Pseudo-codon models with site specific models or branches
// ---
// ModelPar contains three numbers that specify how model parameters are shared between codon positions. { 0 , 1, 0} indicated that pos1 and pos3 have the same model {0} that is seperate to pos2 {1}. There's a strong requirement that the numbers must be in order and maximum value 2;
// BranchPar contains three numbers that specify how branch parameters are shared between codon positions.
class CSiteCodon : public CBaseModel {
public:
	// Constructor
	CSiteCodon(CData *Data, CTree *Tree, vector <int> ModelPar, vector <int> BranchPar, EModel CoreModel, bool WithGamma);
	// Destructor
	~CSiteCodon();
	// Interaction functions
	double lnL(bool ForceReal = false);				// Over-ride of the virtual function. Basically just calls NormaliseParameters, then calculates the likelihood as normal
	void DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol,bool AllowUpdate = true);
				// Override to call DoBraOpt. Just calls NormaliseParameters and then lets the function do its thing


private:
	// Variables
	vector <CData *> m_vpDataSites;					// Stores the data for site 0, 1, and 2; Will either have all three sites or be empty
	vector <CTree *> m_vpTreeSites;					// Stores the trees for site 0,1, and 2; These may be pointers to a previous sites tree. For example, m_vpTreeSites[0] == m_vpTreeSites[2]. So be careful with memory changes!
	vector <int> m_viModelMap;						// Map of the models between sites. Always of size 3

	//Functions
	bool NormaliseParameters();		// Enforces the same parameters between the models sharing the same site number in m_viModelMap. If SetOptToo == true, then it will set the optimiser off for those parameters as well.
};

/////////////////////////////////////////////////////////////////
// 		Coevolution models
/////////////////////////////////////////////////////////////////

// Generic base class
class CBaseCoevo : public CBaseModel {
public:
	// Construction and destruction
	CBaseCoevo(CData *JointData, CData* Data1, CData* Data2, CTree *Tree);
	~CBaseCoevo();
	// Likelihood functions for working with single and double process
	double lnL(bool ForceReal = false);	// perform a likelihood calculation
	double DoBralnL(int B, int NL,int NR);	// Do calculations for a single branch

protected:
	// Variables
	CData * m_pSingleData1;			// Holds data set 1
	CData * m_pSingleData2;			// Holds data set 2
	CBaseProcess *m_pSingleProc1;	// Holds the process for data set 1
	CBaseProcess *m_pSingleProc2;	// Holds the process for data set 2

};

class CWAGCoevo : public CBaseCoevo {
public:
	// Construction and destructions
	CWAGCoevo(CData *JointData, CData *Data1, CData* Data2, CTree *Tree);
	~CWAGCoevo();

private:
};


#endif
