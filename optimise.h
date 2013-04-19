// Optimiser header file

#ifndef __OPT_HEADER
#define __OPT_HEADER

#include "TreeList.h"
#include "model.h"


// Structure holding all the information required to make rearrangements on a phylogeny using CopyTree
struct TreeArrange	{
	CTree Tree;					// The tree
	int NodeNum;				// The node that the rearrangement was made on
	vector <int> LeafList;		// The list of leaves
	vector <int> NodCovList;	// The nodes covered
	int TabuDist;				// Tabu distance
	double lnL;					// Its likelihood
};

enum SNAPchange { none, tabu, normal, star };
// Details about the node under consideration
// NodeNum = Node number in tree
// Length = Length statistic used for sorting (Currently average branch length of collapseable nodes)
// ColSize = Size the node collapse down to
// LeafList = List of nodes describing in the collapsed tree (Leaf < NoSp = real_sequence);
// LeafFrom = The node from which each LeafList comes from
struct NodeLength { int NodeNum; float Length; int ColSize; vector <int> LeafList; vector <int> LeafFrom; };

void CopyTree(CTree *OriTree, CTree *SubTree,int NodeNum, vector <int> LeafList, vector <int> FromList);

double GoldenSection(double OrilnL, double *x, CPar *Par,CBaseModel *M);
double lnsrch(vector <double *> x,double fold,vector <double> g, double p[], double pold[], double *f, bool Do_GS,CBaseModel *Model);
double MulD_Optimise(double OrilnL,double gtol ,double ltol,vector <double *> x,CBaseModel *Model,int NoIterations, bool DoBasicOutput = false,bool OnlyBranches = false, int OptTol = 2, bool NewOne = true, double lnL2Best = BIG_NUMBER,int NoBranchOpt = NUM_FAST_BRA_OPT, bool AllowOnlyParOpt = true,bool TryReallyHard = true);
double DoOnlyParOpt(double OrilnL,double gtol ,double ltol,vector <double *> x,CBaseModel *Model,int NI,int OptTol, int MaxIter = 3);
double PredictlnL(double *Last, int n);	// Function that predicts what the optimal likelihood will be after examining n previous likelihoods
bool CheckAllPar(CBaseModel *M, double lnL, vector <double *> x, double Tol, ostream &os = cout);
bool HardCheckOpt(CBaseModel *M, double lnL, double *x, double Tol, int ParNum,ostream &os = cout);	//	Hard checks that parameter *x is at an optima within +- Tol
double SubSetlnsrch(double Prob, vector <double *> x,double *step_xi, vector <double> g,double lnL,CBaseModel *M);		// Optimises a subset of x, each with probability Prob (i.e. binomial)
double Get2ndDer(CPar *Par,CBaseModel *M, double lnL = -BIG_NUMBER); // Get the numerical second derivative of a parameter

// Tree search routines routines
///////////////////////////////////////////
// General
double DoMPHillClimb(CBaseModel *M,CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, double LimitMPScore = 0,bool DoOutput = true);
// SNAP specific
double TreeSNAP(CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, CBaseModel *Model, bool DoOutput = true, LazyType DoLazy = LAZY_SUBTREE);
void GetParsimonySNAPList(vector <double> *Score,CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, CBaseModel *Model);		// Get list of proposed changes
SNAPchange DoScoreSubTree(int NodeNum, double *BestlnL, CTreeList *TL, CBaseModel *Model, vector <TreeArrange> *TA, LazyType DoLazy);
SNAPchange DoSimpleSubTree(int NodeNum, double *BestlnL, CTreeList *TL, CBaseModel *Model, vector <TreeArrange> *TA, LazyType DoLazy);
vector <STreeScore> DoSubTreeCalcs(double OrilnL, vector <STreeScore> *OriList, CTreeList *TL, CBaseModel *M,int SetSize, bool DoFullOpt, LazyType DoLazy,double PropDiff = BIG_NUMBER, bool ForceSetSize = true);
// SPR specific routines
double TreeSPR(CBaseModel *M, bool DoOutput, double OrilnL = -BIG_NUMBER,bool TryHard=false, bool DoAll=true, int Only1SPR=10, LazyType DoLazy = LAZY_SUBTREE);	// SPR driving routine; if(!Only1SPR) will continue until no improvement found
double DoSimpleSPR(CBaseModel *M, bool DoOutput = true, double OrilnL = -BIG_NUMBER,bool TryHard=false, bool DoAll=true, LazyType DoLazy = LAZY_SUBTREE);					// The routine that actually applies SPR for parsimony to the tree
double DoScoreSPR(CBaseModel *M, bool DoOutput = true, double OrilnL = -BIG_NUMBER,bool TryHard=false, bool DoAll=true, LazyType DoLazy = LAZY_SUBTREE);					// The routine that actually applies SPR for ML and RMSD to the tree
////////////////////////////////////////////////////
// Likelihood optimiser routines
double PreOpt(CBaseModel *Model, bool DoOutput);		// Pre-optimiser that gets branch lengths and initial parameters for models
double FullOpt(CBaseModel *Model, bool DoPar = true, bool DoBra = true, bool DoFreq = false, double CurlnL = -BIG_NUMBER, bool FullLikAcc = true, int MaxNoIterations = DEFAULT_OPTNUM, double lnL2Beat = -BIG_NUMBER, double lnL_tol = FULL_LIK_ACC, bool DoOutput = false, bool TightFullOpt = false);
double LazyOpt(CBaseModel *Model, bool DoPar = true, bool DoBra = true, bool DoFreq = false, double CurlnL = -BIG_NUMBER, bool FullLikAcc = true, int MaxNoIterations = 2, double lnL2Beat = -BIG_NUMBER, bool DoOutput=false);
double LazyBraOpt(CBaseModel *Model, double CurlnL = -BIG_NUMBER,int NoIterations = 2, double Tol = 0.001);
double IntBraOpt(CBaseModel *Model, double CurlnL = -BIG_NUMBER, bool FullLikAcc = true, int MaxNoIterations = DEFAULT_OPTNUM);

// RMSD optimiser
double FMOpt(CBaseModel *Model, vector <double> PWDist);	// Optimise tree according to FitchMargoliash least-squares

// BioNJ
string DoBioNJ(vector <double> PWdists, vector <string> Names, bool DoNumbers = false);
string DoRandomBioNJ(CBaseModel *M, CData *D, bool DoJCDists = false);

// Step wise addition pseudo-importance sampler
void GetSAStartTree(CTree *T, CBaseModel *M, bool DoOutput = true, LazyType DoLazy = LAZY_SA, int DoMP = 0);
double GetTree(CBaseModel *M, CTreeList *Four, CTreeList *Five, CTreeList *Six,LazyType DoLazy = LAZY_SA, bool DoOutput = true);
void AddSequences(bool DoOutput, vector <int> *SeqRem, CBaseModel *M,CTree *Tree, double DespProb = 0.0, bool AllowRandomSnap = true,LazyType DoLazy = LAZY_SA);
void GetNewTree(int NumRandomSteps,CTree *Tree,CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, LazyType DoLazy = LAZY_SA);
void GetRMSDSampTree(CBaseModel *M, CTreeList *Four, CTreeList *Five, CTreeList *Six,LazyType DoLazy = LAZY_SA,double ProbParsimony = PROB_DO_MP_SA,bool DoOutput = true);
void GetSATree(vector <STreelnL *> *Trees, int SeqAdd, CBaseModel *M,double ProbSubOptBra = 0.0,LazyType DoLazy = LAZY_SA);
void BraSACalc(int Seq, vector <STreelnL*>*, CBaseModel *M,bool DoTripleOpt, double *lnL2Beat = NULL, int SPR=-1, int OriBr = -1,LazyType DoLazy = LAZY_SA, CTree *OriT = NULL);
void Branch_SA(int Seq, vector <STreelnL*>*, int First, int NTo, int NFr, CBaseModel *M,bool DoTripleOpt,double *lnL2Beat = NULL, int SPR=-1, int OriBr = -1,LazyType DoLazy = LAZY_SA, CTree *OriT = NULL);
void DoSA(int Seq, vector <STreelnL*>*, int First, int NTo, int NFr, int Br, CBaseModel *M, bool IsExtBra,bool DoTripleOpt, double *lnL2Beat = NULL,int SPR=-1, int OriBr = -1,LazyType DoLazy = LAZY_SA, CTree *OriT = NULL);

// Some STreelnL functions
void SortList(vector <STreelnL *> *L);
STreelnL GetBestList(vector <STreelnL *> *L);
void CleanList(vector <STreelnL *> *L);


// Pairwise distance routines
vector <double> GetPW(CBaseModel *Model, vector <double> *CurrentPW = NULL, bool GetJCDist = false, bool GetStdErr = false);					// Get pairwise distances
double GetDist(CBaseModel *M, CData *D,int i, int j,double Current = -1.0);			// Get a specific pairwise distance
double GetDistVar(CBaseModel *M, CData *D,int Seq1, int Seq2, double CurrentDist);	// Get the variance of a pairwise distance estimate

// Graphing functions
void Parameter1DGraph(CBaseModel *M, CPar *Par,double low=0.0,double high = 10.0, double step = 1.0, bool OptAll = true,ostream &os = cout);
void Parameter2DGraph(CBaseModel *M, CPar *Par1,CPar *Par2, double l1=0.0,double h1=10.0,double s1 = 1.0,double l2=0.0,double h2=10.0,double s2 = 1.0, bool OptAll = true, ostream &os = cout);

#endif
