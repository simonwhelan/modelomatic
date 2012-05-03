
// Predefined lists of tree topologies

#include "tree.h"
#include "data.h"
#include "model.h"

#ifndef __TREELIST_H
#define __TREELIST_H

class CDataSummary;

//////////////////////////////////////////////////////////////////////
// Class for phylogenetic data
class CPhyloDat {
public:
	static int SimNum;					// Number of the simulation
	// Main functions
	CPhyloDat();						// Constructor
	~CPhyloDat();						// Destructor
	// General Interaction
	bool NeedMenu();					// Whether to get info from menu
	// Data file
	void SetIn(string S);				// Set input file
	string In() { return m_sInFile;	}	// Returns input file name
	void CleanData();					// Cleans the datafile
	void ClearDataPointer() { m_pData = NULL; };			// Clears the data pointer without deleting it...
	EDataType DataType();				// Returns the datatype
	int Size();							// Length of sequences
	int TrueSize();						// Length of true sequences
	int NoSeq();						// Number of sequences
	void GetData();						// General get data function that decides whether to get single data or phylogenomic data
	void GetSingleData();				// Gets sequence data for a single gene (or concatenation)

	CData *pData() { return m_pData; }	// Returns a pointer to the data

	// Simulation functions
	void DoSimulation(string ModelInputFile, string DataOutputFile, int NumberOfSims);
	void SimulateData(bool Force=true,bool MakeData=false, CDataSummary *DS = NULL,vector <vector <bool> > *GapMask = NULL);	// Function to simulate data; if(Force) then it will overwrite what is in m_pData;

	// Partial likelihood generation functions
	void DoPartL(CData *Dat, string model_file, string out_file);	// Highest level function that calculates and output the partial likelihoods per branch

	// Multiple data set management
	void GetMultiData();				// Gets sequence data for multiple genes
	void CleanMultiData();				// Cleans the multiple data set
	int NoDataSets();					// Returns the number of data sets
	void SetMultiSet(int i);			// Sets the data for data analysis to
	int CurrentMulti() { return m_iCurrentMulti; }	// Returns the current m_iCurrentMulti
	int GetMultiRFDist(int Gene1,int Gene2);	// Calculate the RF distance between the trees for Gene1 and Gene2

	// Tree file
	void SetTreeFile(string S);			// Set the tree file
	bool IsTree();						// Whether a tree exists
	string TreeFile() { return m_sTreeFile; }	// Returns tree file
	void CleanTree();					// Cleans tree memory
	void GetTree();						// Gets the tree
	CTree *pTree() { return m_pTree; }	// Returns pointer to tree
	CTree *CreateTree(string T);			// Create a tree
	CTree *CreateTree();					// Creates a blank tree

	// Output file
	void SetOut(string S);				// Set output file
	string Out(){ return m_sOutFile; }	// Returns output file name
	void GetOutputFile();				// Gets a valid output file name

	// Substition model stuff
	void ReadModel(string InputFile = "\0",CData *D = NULL); 	// Reads a model and parameters from file
	void CleanModel();					// Cleans the model
	void CreateModel();					// Creates the model
	CBaseModel *pModel() { return m_pModel; }							// Returns a pointer to model
	void SetpModel(CBaseModel *M, bool AllowForce = false) { assert(AllowForce || m_pModel == NULL); m_pModel = M; };	// Sets the model (if not already done!)
	EModel Model() { return m_Model; }									// Returns the type of model it is
	double ModelAcc() { return m_dModelOptTol; }						// Returns the accuracy to which the model is originally optimised
	void SetModelAcc(double acc) { assert(acc >= FULL_LIK_ACC); m_dModelOptTol = acc; }
	void SetModel(EModel M) { m_Model = M; }							// Set the model to use
	bool SetModel(string s);											// Set model from a string
	void SetGamma(bool V) { m_bDoGamma = V; }							// Set m_bDoGamma
	void SetInv(bool V) { m_bDoInv = V; }								// Set m_bDoInv
	void SetGamCat(int V) { if(!InRange(V,0,10)) { Error("Gamma cat out of range\n"); } m_iGamCat = V; }	// Set the number of gamma categories
	void SetDoF(bool F) { m_bF = F; }									// Set the value of Do F
	bool DoGam()	{ return m_bDoGamma; }
	bool DoInv()	{ return m_bDoInv; }
	bool DoF();															// Whether the model has +F option (for amino acid models)
	int NoGamCat()	{ return m_iGamCat; }
	void SetRandomSeed(long l)  { if(m_pModel == NULL) { Error("\nTrying to PhyloDat::SetRandomSeed before model is specified\n"); } m_pModel->SetRandomSeed(l); };
	long GetOriRandomSeed()		{ if(m_pModel == NULL) { Error("\nTrying to PhyloDat::GetOriRandomSeed before model is specified\n"); } return m_pModel->GetOriRandomSeed(); };
	long GetCurrentRandomSeed()	{ if(m_pModel == NULL) { Error("\nTrying to PhyloDat::GetCurrentRandomSeed before model is specified\n"); } return m_pModel->GetCurrentRandomSeed(); };
	void LikelihoodSurfaces(ostream & os = cout, vector <CPar *> *Pars = NULL);	// Output the likelihood surface for each parameter in turn in the model (or in vector <CPar *> Pars)
	bool Compress()   { return m_bAllowCompress; }
	void DoCompress() { m_bAllowCompress = true; }
	void NoCompress() { m_bAllowCompress = false; }

	// Tree estimation stuff
	void SetSA(bool V) { m_bDoSA = V; }			// Sets the DoSA variable
	void SetBioNJ(bool V) { m_bDoBioNJ = V; }	// Sets the DoBioNJ variable
	bool SA()		{ return m_bDoSA; }			// Whether to do stepwise addition
	bool bionj()	{ return m_bDoBioNJ; }		// Whether to do bionj
	void SetStartTree(string s);				// Set the tree estimation stuff
	// THMM stuff
	void GetTHMMStdIn();						// Get THMM stuff from stdin
	string GetFullDNATHMMString(string Str);			// Get THMM stuff from string for full DNA models
	string GetDNATHMMString(string Str);		// Get THMM stuff from string for DNA models
	string GetAATHMMString(string Str);			// Get THMM stuff from string for AA models
	bool VerifyTHMM();							// Checks that the THMM is valid
	string THMM_Details();						// Function that provides a suffix name for a THMM
	// Codon stuff
	void GetCodonString(string s);
	// What to optimise
	bool OptBra() { return m_bOptBra; }
	bool OptPar() { return m_bOptPar; }
	bool OptFre() { return m_bOptFre; }
	// Output function
	friend ostream &operator<<(ostream &os, CPhyloDat Dat);
protected:
	// Variables
	CBaseModel *m_pModel;				// The model
	EModel m_Model;						// The model choice
	bool m_bDoGamma;					// Do gamma distribution
	int m_iGamCat;						// Number of categories in the gamma distribution
	bool m_bDoInv;						// Do invariant sites
	bool m_bDoGarbage;					// Do Garbage collector
	bool m_bF;							// Whether the model has the +F option (a.a. models only)
	bool m_bDoTreeHMM;					// Whether to use a tree HMM model rather than a simple
	// Optimise parameters
	bool m_bOptBra, m_bOptPar, m_bOptFre;	// The parameters to optimise in initial optimisations
	bool m_bAllowCompress;					// Whether the optimiser will compress the model for tree-search
	// Some THMM specific stuff
	int THMM_NoCat;						// Number of hidden states
	ERateTypes THMM_R;					// What type of variable Rates to include
	bool THMM_F;						// Whether to include variable frequencies
	bool THMM_K;						// Whether to include variable kappa
	EHIDDEN_TYPE THMM_H;				// The type of transitions between hidden states allowed by the model
	// Some codon specific stuff
	int m_iGenCode;						// The genetic code
	ECodonEqm m_CodonEqm;				// The frequency parameters
	// Other things
	int m_iPseudoDataSize;				// Data length; recovered from file when reading a model
	int m_iPseudoDataSeq;				// Num seq; recovered from file when reading a model
	EDataType m_PseudoDataType;			// Data type
	CData *m_pData;						// The data
	CTree *m_pTree;						// The tree
	string m_sInFile;					// The data file
	string m_sOutFile;					// The output file
	string m_sTreeFile;					// The tree file
	bool m_bDoSA;							// Whether stepwise addition is performed
	bool m_bDoBioNJ;						// Whether bionj is performed
	double m_dModelOptTol;				// The degree to which the likelihood of the model is optimised.
	// Information for multigene data sets
	int m_iCurrentMulti;				// Current MultiDataSet
	vector <CData *> m_vpMultiDat;		// Vector containing the sequence data multiple data sets
	vector <CTree *> m_vpMultiTree;	// Vector containing the trees for multigene data sets
	vector <double> m_vdMultiScore;	// Vector containing the scores for multigene data sets
	vector <string > m_vsMultiName;	// Vector containing the names for multigene data sets
};

/////////////////////////////////////////////////////////////////////
// Class defining some summary statistics of a sequence alignment
class CDataSummary {
public:
	// Constructor
	CDataSummary();		// Blank constructor
	CDataSummary(string Name,EDataType Type, vector <vector <int> > *Seqs, vector <string> *Names);	// Construct from sequences
	CDataSummary(string Name,CData *D);																// Construct from a data object
	// Underlying simple constructor
	void MakeDataSummary(string Name,EDataType Type, vector <vector <int> > *Seqs, vector <string> *Names);
	friend ostream &operator <<(ostream &os, CDataSummary DS);
	int NoStats() { if(m_dlnL > 0) { return (3*NumStates(m_DataType)) + 1; } return (3*NumStates(m_DataType)) + 4; }
	// Attributes
	string m_sName;									// Name of the data
	bool m_bDoFullOutput;							// Whether to output the header too
	EDataType m_DataType;							// The type of sequence data
	int m_iSize, m_iNoSeq;							// The size of the data matrix
	int m_iNoDataPatterns;							// Number of unique data patterns
	double m_dlnL,m_dlnL_multi,m_ddelta;						// Test statistics for Goldman test
	vector <string> m_vsNames;						// Names of the sequences
	vector <double> m_vdFrq;						// The frequencies of the different states (this really need to be between species...)
	vector <double> m_vdFrqStdErr;					// Gets the standard error for the frequencies
	vector <int> m_vdStatesObs;					// Distribution of the number of observed states at sites;
};

//////////////////////////////////////////////////////////////////////
// General tabuing functions
bool IsTabu(CTree *Tree, bool OnlyRadTabuOld,bool AddToTabu = true, CBaseModel *M = NULL, double lnL = -BIG_NUMBER, bool IsOptima = false);	// Checks whether calculations need to be done for this tree.
int MinTabuRFDist(CTree * Tree, bool CompareOnly2Old = false);
void AgeTabuTrees();									// Age all the trees in the tabu list
bool OldTabuTrees();
void VivifyTabuTrees();								// Make all the tabu trees young again
void RemoveTabuTree(CTree *Tree);
bool IsSameTree(CTree *T1, CTree *T2);

struct STabuTree	{
	CTree Tree;									// The tree
	double lnL;									// The trees optimal likelihood
	vector <double> vPar;									// The parameter values that gave the optimal likelihood
	vector <int> List;		// The identifying functions of the tree
	bool Optima;								// Whether the tree is an optima
};

struct STreeScore { int TreeNum; double Lik; int TabuDist; };
void SortSTreeScore(vector <STreeScore> *x);


class CTreeList	{
public:
	CTreeList(int NoSp);
	// The trees store
	vector <CTree> m_Trees;				// Vector of trees to be examined
	list <STreeScore> m_TreeEst;			// List of tree estimates made
	list <STreeScore>::iterator m_ListPos;	// Iterator for the STreeList structure
	int m_iNoTree;
	int m_ariLeafSp[6];						// Whether the external nodes in a tree are actual sequences (-1 = Partial likelihood, 0 > n = species number);

};

// Structure storing trees and their likelihoods
// Used in the Stepwise addition algorithm and will eventually be used in the TABU search.
struct STreelnL {
	CTree Tree;			// Tree examined
	double lnL;			// likelihood score
	double score;		// Score associated with SA algorithm
	int BraNum;			// Branch sequence was added to
	int SPRID;		// Identifier for SPR routine
	vector <int> BraLinks;	// The Branches that the sequence added branch associates with (up to 4)
	bool IsOpt;			// Whether this branch represents an optima
	bool DoOpt;			// Whether to optimise this branch again
};
void OutSTreelnL(vector <STreelnL*> *TL, bool FullDetail = false, ostream &os= cout);;

const char Sp4Tree[3][13]	= {	"((1,4),2,3);",
								"(1,(2,4),3);",
								"(1,2,(3,4));"};
const char Sp5Tree[15][19] = {	"(((1,4),5),2,3);",
								"((1,4),(2,5),3);",
								"((1,4),2,(3,5));",
								"(((1,5),4),2,3);",
								"((1,(4,5)),2,3);",
								"((1,5),(2,4),3);",
								"(1,((2,4),5),3);",
								"(1,(2,4),(3,5));",
								"(1,((2,5),4),3);",
								"(1,(2,(4,5)),3);",
								"((1,5),2,(3,4));",
								"(1,(2,5),(3,4));",
								"(1,2,((3,4),5));",
								"(1,2,((3,5),4));",
								"(1,2,(3,(4,5)));"};
const char Sp6Tree[105][24] = {	"((((1,4),5),6),2,3);",
								"(((1,4),5),(2,6),3);",
								"(((1,4),5),2,(3,6));",
								"((((1,6),4),5),2,3);",
								"(((1,(4,6)),5),2,3);",
								"((((1,4),6),5),2,3);",
								"(((1,4),(5,6)),2,3);",
								"(((1,4),6),(2,5),3);",
								"((1,4),((2,5),6),3);",
								"((1,4),(2,5),(3,6));",
								"(((1,6),4),(2,5),3);",
								"((1,(4,6)),(2,5),3);",
								"((1,4),((2,6),5),3);",
								"((1,4),(2,(5,6)),3);",
								"(((1,4),6),2,(3,5));",
								"((1,4),(2,6),(3,5));",
								"((1,4),2,((3,5),6));",
								"(((1,6),4),2,(3,5));",
								"((1,(4,6)),2,(3,5));",
								"((1,4),2,((3,6),5));",
								"((1,4),2,(3,(5,6)));",
								"((((1,5),4),6),2,3);",
								"(((1,5),4),(2,6),3);",
								"(((1,5),4),2,(3,6));",
								"((((1,5),6),4),2,3);",
								"(((1,5),(4,6)),2,3);",
								"((((1,6),5),4),2,3);",
								"(((1,(5,6)),4),2,3);",
								"(((1,(4,5)),6),2,3);",
								"((1,(4,5)),(2,6),3);",
								"((1,(4,5)),2,(3,6));",
								"(((1,6),(4,5)),2,3);",
								"((1,((4,5),6)),2,3);",
								"((1,((4,6),5)),2,3);",
								"((1,(4,(5,6))),2,3);",
								"(((1,5),6),(2,4),3);",
								"((1,5),((2,4),6),3);",
								"((1,5),(2,4),(3,6));",
								"((1,5),((2,6),4),3);",
								"((1,5),(2,(4,6)),3);",
								"(((1,6),5),(2,4),3);",
								"((1,(5,6)),(2,4),3);",
								"((1,6),((2,4),5),3);",
								"(1,(((2,4),5),6),3);",
								"(1,((2,4),5),(3,6));",
								"(1,(((2,6),4),5),3);",
								"(1,((2,(4,6)),5),3);",
								"(1,(((2,4),6),5),3);",
								"(1,((2,4),(5,6)),3);",
								"((1,6),(2,4),(3,5));",
								"(1,((2,4),6),(3,5));",
								"(1,(2,4),((3,5),6));",
								"(1,((2,6),4),(3,5));",
								"(1,(2,(4,6)),(3,5));",
								"(1,(2,4),((3,6),5));",
								"(1,(2,4),(3,(5,6)));",
								"((1,6),((2,5),4),3);",
								"(1,(((2,5),4),6),3);",
								"(1,((2,5),4),(3,6));",
								"(1,(((2,5),6),4),3);",
								"(1,((2,5),(4,6)),3);",
								"(1,(((2,6),5),4),3);",
								"(1,((2,(5,6)),4),3);",
								"((1,6),(2,(4,5)),3);",
								"(1,((2,(4,5)),6),3);",
								"(1,(2,(4,5)),(3,6));",
								"(1,((2,6),(4,5)),3);",
								"(1,(2,((4,5),6)),3);",
								"(1,(2,((4,6),5)),3);",
								"(1,(2,(4,(5,6))),3);",
								"(((1,5),6),2,(3,4));",
								"((1,5),(2,6),(3,4));",
								"((1,5),2,((3,4),6));",
								"((1,5),2,((3,6),4));",
								"((1,5),2,(3,(4,6)));",
								"(((1,6),5),2,(3,4));",
								"((1,(5,6)),2,(3,4));",
								"((1,6),(2,5),(3,4));",
								"(1,((2,5),6),(3,4));",
								"(1,(2,5),((3,4),6));",
								"(1,(2,5),((3,6),4));",
								"(1,(2,5),(3,(4,6)));",
								"(1,((2,6),5),(3,4));",
								"(1,(2,(5,6)),(3,4));",
								"((1,6),2,((3,4),5));",
								"(1,(2,6),((3,4),5));",
								"(1,2,(((3,4),5),6));",
								"(1,2,(((3,6),4),5));",
								"(1,2,((3,(4,6)),5));",
								"(1,2,(((3,4),6),5));",
								"(1,2,((3,4),(5,6)));",
								"((1,6),2,((3,5),4));",
								"(1,(2,6),((3,5),4));",
								"(1,2,(((3,5),4),6));",
								"(1,2,(((3,5),6),4));",
								"(1,2,((3,5),(4,6)));",
								"(1,2,(((3,6),5),4));",
								"(1,2,((3,(5,6)),4));",
								"((1,6),2,(3,(4,5)));",
								"(1,(2,6),(3,(4,5)));",
								"(1,2,((3,(4,5)),6));",
								"(1,2,((3,6),(4,5)));",
								"(1,2,(3,((4,5),6)));",
								"(1,2,(3,((4,6),5)));",
								"(1,2,(3,(4,(5,6))));"};

#endif
