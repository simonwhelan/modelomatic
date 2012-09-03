 /* **********************************************************
 * Header file defining what Markov processes are
 * Each process needs with associated with it the following:
 * 	- Data
 * 	- Tree
 *	- A calculation function
 * Each process can be used to calculate partial likelihoods up to
 * given nodes in the tree or get full likelihoods (usual calcs + mixture models)
 * ********************************************************* */

#ifndef PROCESS_HEADER
#define PROCESS_HEADER

#include "Leaphy.h"
#include "data.h"
#include "tree.h"

static int ProcID = 0;
static int QMatID = 0;

class CBaseProcess;

// Some definitions for THMM models
#define DEBUG_HMP_MODEL 1
#define HMP_MODEL_2_SIMPLE 0	// Forces the THMM model to mimic a standard model, useful for debugging.

// Some other definitions
enum ECodonEqm { cEQU,F1X4,F3X4,F64 };				// Defines the different types of codon frequencies

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino acid model information
const double dWAGVal[190] = {0.551571,0.509848,0.635346,0.738998,0.147304,5.429420,1.027040,0.528191,0.265256,0.0302949,0.908598,3.035500,1.543640,0.616783,0.0988179,1.582850,0.439157,0.947198,6.174160,0.021352,5.469470,1.416720,0.584665,1.125560,0.865584,0.306674,0.330052,0.567717,0.316954,2.137150,3.956290,0.930676,0.248972,4.294110,0.570025,0.249410,0.193335,0.186979,0.554236,0.039437,0.170135,0.113917,0.127395,0.0304501,0.138190,0.397915,0.497671,0.131528,0.0848047,0.384287,0.869489,0.154263,0.0613037,0.499462,3.170970,0.906265,5.351420,3.012010,0.479855,0.0740339,3.894900,2.584430,0.373558,0.890432,0.323832,0.257555,0.893496,0.683162,0.198221,0.103754,0.390482,1.545260,0.315124,0.174100,0.404141,4.257460,4.854020,0.934276,0.210494,0.102711,0.0961621,0.0467304,0.398020,0.0999208,0.0811339,0.049931,0.679371,1.059470,2.115170,0.088836,1.190630,1.438550,0.679489,0.195081,0.423984,0.109404,0.933372,0.682355,0.243570,0.696198,0.0999288,0.415844,0.556896,0.171329,0.161444,3.370790,1.224190,3.974230,1.071760,1.407660,1.028870,0.704939,1.341820,0.740169,0.319440,0.344739,0.967130,0.493905,0.545931,1.613280,2.121110,0.554413,2.030060,0.374866,0.512984,0.857928,0.822765,0.225833,0.473307,1.458160,0.326622,1.386980,1.516120,0.171903,0.795384,4.378020,0.113133,1.163920,0.0719167,0.129767,0.717070,0.215737,0.156557,0.336983,0.262569,0.212483,0.665309,0.137505,0.515706,1.529640,0.139405,0.523742,0.110864,0.240735,0.381533,1.086000,0.325711,0.543833,0.227710,0.196303,0.103604,3.873440,0.420170,0.398618,0.133264,0.428437,6.454280,0.216046,0.786993,0.291148,2.485390,2.006010,0.251849,0.196246,0.152335,1.002140,0.301281,0.588731,0.187247,0.118358,7.821300,1.800340,0.305434,2.058450,0.649892,0.314887,0.232739,1.388230,0.365369,0.314730};
const double dWAGFreq[20] = {0.0866279,0.043972,0.0390894,0.0570451,0.0193078,0.0367281,0.0580589,0.0832518,0.0244313,0.048466,0.086209,0.0620286,0.0195027,0.0384319,0.0457631,0.0695179,0.0610127,0.0143859,0.0352742,0.0708956};
const double dJTTVal[190] = {58,54,45,81,16,528,56,113,34,10,57,310,86,49,9,105,29,58,767,5,323,179,137,81,130,59,26,119,27,328,391,112,69,597,26,23,36,22,47,11,17,9,12,6,16,30,38,12,7,23,72,9,6,56,229,35,646,263,26,7,292,181,27,45,21,14,54,44,30,15,31,43,18,14,33,479,388,65,15,5,10,4,78,4,5,5,40,89,248,4,43,194,74,15,15,14,164,18,24,115,10,102,21,16,17,378,101,503,59,223,53,30,201,73,40,59,47,29,92,285,475,64,232,38,42,51,32,33,46,245,25,103,226,12,118,477,9,126,8,4,115,18,10,55,8,9,52,10,24,53,6,35,12,11,20,70,46,209,24,7,8,573,32,24,8,18,536,10,63,21,71,298,17,16,31,62,20,45,47,11,961,180,14,323,62,23,38,112,25,16};
const double dJTTFreq[20] = {7.674791534958306e-2,5.169091033818206e-2,4.26450085290017e-2,5.1543910308782054e-2,1.980300396060079e-2,4.075200815040162e-2,6.182991236598247e-2,7.315191463038292e-2,2.2944004588800915e-2,5.3760910752182145e-2,9.190391838078367e-2,5.867591173518234e-2,2.382600476520095e-2,4.01260080252016e-2,5.090091018018203e-2,6.876491375298274e-2,5.856491171298234e-2,1.4261002852200569e-2,3.210200642040128e-2,6.600491320098263e-2}; //renormalised
const double dDAYVal[190] = {27,98,32,120,0,905,36,23,0,0,89,246,103,134,0,198,1,148,1153,0,716,240,9,139,125,11,28,81,23,240,535,86,28,606,43,10,65,64,77,24,44,18,61,0,7,41,15,34,0,0,73,11,7,44,257,26,464,318,71,0,153,83,27,26,46,18,72,90,1,0,0,114,30,17,0,336,527,243,18,14,14,0,0,0,0,15,48,196,157,0,92,250,103,42,13,19,153,51,34,94,12,32,33,17,11,409,154,495,95,161,56,79,234,35,24,17,96,62,46,245,371,26,229,66,16,53,34,30,22,192,33,136,104,13,78,550,0,201,23,0,0,0,0,0,27,0,46,0,0,76,0,75,0,24,8,95,0,96,0,22,0,127,37,28,13,0,698,0,34,42,61,208,24,15,18,49,35,37,54,44,889,175,10,258,12,48,30,157,0,28};
const double dDAYFreq[20] = {0.087127,0.040904,0.040432,0.046872,0.033474,0.038255,0.049530,0.088612,0.033618,0.036886,0.085357,0.080482,0.014753,0.039772,0.050680,0.069577,0.058542,0.010494,0.029916,0.064718};
const double dmtREVVal[190] = {23.18,26.95,13.24,17.67,1.90,794.38,59.93,103.33,58.94,1.90,1.90,220.99,173.56,55.28,75.24,9.77,1.90,63.05,583.55,1.90,313.56,120.71,23.03,53.30,56.77,30.71,6.75,28.28,13.90,165.23,496.13,113.99,141.49,582.40,49.12,1.90,96.49,1.90,27.10,4.34,62.73,8.34,3.31,5.98,12.26,25.46,15.58,15.16,1.90,25.65,39.70,1.90,2.41,11.49,329.09,8.36,141.40,608.70,2.31,1.90,465.58,313.86,22.73,127.67,19.57,14.88,141.88,1.90,65.41,1.90,6.18,47.37,1.90,1.90,11.97,517.98,537.53,91.37,6.37,4.69,15.20,4.98,70.80,19.11,2.67,1.90,48.16,84.67,216.06,6.44,90.82,54.31,23.64,73.31,13.43,31.26,137.29,12.83,1.90,60.97,20.63,40.10,50.10,18.84,17.31,387.86,6.04,494.39,69.02,277.05,54.11,54.71,125.93,77.46,47.70,73.61,105.79,111.16,64.29,169.90,480.72,2.08,238.46,28.01,179.97,94.93,14.82,11.17,44.78,368.43,126.40,136.33,528.17,33.85,128.22,597.21,1.90,21.95,10.68,19.86,33.60,1.90,1.90,10.92,7.08,1.90,32.44,24.00,21.71,7.84,4.21,38.58,9.99,6.48,1.90,191.36,21.21,254.77,38.82,13.12,3.21,670.14,25.01,44.15,51.17,39.96,465.58,16.21,64.92,38.73,26.25,195.06,7.64,1.90,1.90,1.90,19.00,21.14,2.53,1.90,1222.94,91.67,1.90,387.54,6.35,8.23,1.90,204.54,5.37,1.90};
const double dmtREVFreq[20] = {0.072,0.019,0.039,0.019,0.006,0.025,0.024,0.056,0.028,0.088,0.169,0.023,0.054,0.061,0.054,0.072,0.086,0.029,0.033,0.043};
const double dLGVal[190] = {0.425093,0.276818,0.751878,0.395144,0.123954,5.076149,2.489084,0.534551,0.528768,0.062556,0.969894,2.807908,1.695752,0.523386,0.084808,1.038545,0.363970,0.541712,5.243870,0.003499,4.128591,2.066040,0.390192,1.437645,0.844926,0.569265,0.267959,0.348847,0.358858,2.426601,4.509238,0.927114,0.640543,4.813505,0.423881,0.311484,0.149830,0.126991,0.191503,0.010690,0.320627,0.072854,0.044265,0.008705,0.108882,0.395337,0.301848,0.068427,0.015076,0.594007,0.582457,0.069673,0.044261,0.366317,4.145067,0.536518,6.326067,2.145078,0.282959,0.013266,3.234294,1.807177,0.296636,0.697264,0.159069,0.137500,1.124035,0.484133,0.371004,0.025548,0.893680,1.672569,0.173735,0.139538,0.442472,4.273607,6.312358,0.656604,0.253701,0.052722,0.089525,0.017416,1.105251,0.035855,0.018811,0.089586,0.682139,1.112727,2.592692,0.023918,1.798853,1.177651,0.332533,0.161787,0.394456,0.075382,0.624294,0.419409,0.196961,0.508851,0.078281,0.249060,0.390322,0.099849,0.094464,4.727182,0.858151,4.008358,1.240275,2.784478,1.223828,0.611973,1.739990,0.990012,0.064105,0.182287,0.748683,0.346960,0.361819,1.338132,2.139501,0.578987,2.000679,0.425860,1.143480,1.080136,0.604545,0.129836,0.584262,1.033739,0.302936,1.136863,2.020366,0.165001,0.571468,6.472279,0.180717,0.593607,0.045376,0.029890,0.670128,0.236199,0.077852,0.268491,0.597054,0.111660,0.619632,0.049906,0.696175,2.457121,0.095131,0.248862,0.140825,0.218959,0.314440,0.612025,0.135107,1.165532,0.257336,0.120037,0.054679,5.306834,0.232523,0.299648,0.131932,0.481306,7.803902,0.089613,0.400547,0.245841,3.151815,2.547870,0.170887,0.083688,0.037967,1.959291,0.210332,0.245034,0.076701,0.119013,10.649107,1.702745,0.185202,1.898718,0.654683,0.296501,0.098369,2.188158,0.189510,0.249313};
const double dLGFreq[20] = {7.906592093407908e-2,5.594094405905594e-2,4.197695802304198e-2,5.3051946948053055e-2,1.2936987063012939e-2,4.0766959233040766e-2,7.158592841407159e-2,5.733694266305734e-2,2.2354977645022357e-2,6.215693784306216e-2,9.908090091909909e-2,6.45999354000646e-2,2.2950977049022953e-2,4.2301957698042306e-2,4.403995596004404e-2,6.11969388030612e-2,5.328694671305329e-2,1.2065987934012068e-2,3.4154965845034156e-2,6.914693085306915e-2}; // renormalised

////////////////////////////////////////////////////////
// Data compression routine
vector <int> GetCompressedData(CTree *T,CData *D);	// Returns vector of size (m_iNoIntNode * m_iSize) describing the compression (see implementation for more details)
void ColSorter(int NT, int NF, CTree *T, CData *D, vector < vector <int> > &ColInfo,bool First = true);
vector <int> DoCompress(vector <vector <int> > &CV, CTree *T, CData *D);

///////////////////////////////////////////////////////////////////////////
// Classes expanding on CPar for Q matrix specific parameters
///////////////////////////////////////////////////////////////////////////

// Standard Q matrix parameters
class CQPar : public CPar	{
public:
	CQPar(string Name, int Char, double Value, bool Optimise=true, double Lowbound = 0.0, double Upbound=MAX_PAR_VALUE,ParOp Oper=MULTIPLY);	// Core parameter construction routine
	// Variables mapping to the Q matrix
	int m_iChar;
	vector <int> m_viQMap;									// Elements of the Q matrix it applies to
	vector <int> m_viApply2Q;								// The list of Q matrix IDs that the parameter applies to (empty = all)
	// Functions
	void AddQij(int i, int j, bool Rev = true);				// Adds coordinates of Q to
	void Par2Q(double *QMat, int QMatID);					// Applies the parameter to the Q matrix of number QNum
	// Updated output function
	ostream &Output(ostream &os);
};

// The alpha parameter of the gamma distribution
#define INITIAL_GAMMA MODEL_INITIAL_GAMMA
#define MAX_ALFA 100
class CGammaPar : public CPar	{
public:
	CGammaPar(string Name, int Char, CPar* Rate,double Value = INITIAL_GAMMA, bool Optimise=true, double Lowbound = 0.01, double Upbound=MAX_ALFA,ParOp Oper=REPLACE);	// Constructor
	~CGammaPar();							// Destructor
	// public functions
	void GlobalApply();						// How the parameter is applied
	void AddRateToGamma(CBaseProcess *Proc);// Control this rate using the gamma distribution
	int NoCat() { return m_iNoCat; }		// Returns the number of rate categories
	ostream &Output(ostream &os);			// Output function
private:
	// Private variables
	int m_iNoCat;						// Number of categories in the gamma distribution
	vector <CQPar *> m_arpRates;		// Rate parameters
	CPar *m_pProcRate;					// Pointer to the rate of the process
	CQPar *m_pProb;						// Probability of the base process
};

// This is the Q matrix class
/////////////////////////////////////////////////////////////////////////
// This contains the space for a single Q matrix and the functions
// for turning it into a P(t) matrix
class CQMat {
public:
	CQMat(EDataType Type,string Name = "Unnamed Q");
	CQMat(int Char, EDataType Type = OTHER,string Name = "Unnamed Q");
	// Access functions
	double *Q(int i=0,int j=0)	{ return &m_ardQMat[(i*m_iChar)+j]; };		// Access to the Q matrix
	int Char()					{ return m_iChar; };						// Access to the number of characters
	vector <double> Eqm();													// Access to the processes equilibrium
	EDataType Type()			{ return m_DataType; };						// Access to the datatype
	int ID()	{ return m_iQMatID; }										// Returns the identifier of the matrix
	void OutQRoot(ostream &os=cout) { int i; FOR(i,m_iChar) { os << m_ardRoot[i] << " "; } };	// Outputs the eigen values
	void OutQ(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardQMat[i] << " "; } } // Outputs the Q matrix
	void OutU(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardU[i] << " "; } } // Outputs the U matrix
	void OutV(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardV[i] << " "; } } // Outputs the V matrix
	void OutRoot(ostream &os=cout)	{ int i; FOR(i,m_iChar)  { os << m_ardRoot[i] << " "; } }
	bool IsLocked() { return FlipBool(m_bAllowModelUpdate); }							// Returns whether the process is locked
	void Lock()		{ m_bAllowModelUpdate = false; }						// Locks the model so the Q matrix and the eigen vectors can't be changed
	void Unlock()	{ m_bAllowModelUpdate = true; }							// Unlock the model so Q and eigens can change
	// Linear algebra functions
	void Decompose(double *Eqm,bool Scale = true, bool REV = true,double Rate = 1.0);			// Decompose the matrix to its eigen vectors and values (Eqm contains equilibrium distribution if known)
	void Decompose(vector <double> Eqm,bool Scale = true,bool REV = true,double Rate = 1.0);	// Decompose the matrix to its eigen vectors and values (Eqm contains equilibrium distribution if known)
	bool MakePT(double T,double *PT);										// Returns the P(t) matrix
	virtual double OverallSubRate(bool ForceWholeProcess = false);			// Returns the overall substitution rate of the process
	virtual double OverallTransRate() { return 0.0; }						// Returns the overall transition rate between hidden states
	virtual double TransRateFrom(int From) { return 0.0; };					// Returns the overall transition rate from a hidden state
	virtual double TransRateTo(int To) { return 0.0; };						// Returns the overall transition rate to a hidden state
	virtual void ScaleQ(double Rate);												// Scale the Q matrix to the overall rate Rate

	// Parameter application functions
	void InitQ(double Val);								// Initialise Q matrix to a particular value
	void ApplyPar2Q(CQPar *Par);						// Apply the parameter to the Q matrix
	void DoQDiag();															// Make all rows sum to 0
	friend ostream &operator<<(ostream &os,CQMat &Mat) { return Mat.Output(os); };

protected:
	// Variables
	bool m_bScaleReady;					// Whether the process is ready to be scaled (requires the eqm distribution to be verified)
	double *m_ardEqm;					// The equilibrium of the process
	double m_dScale;					// Current overall rate of process
	bool m_bAlwaysI;					// When the rate is set to zero, P(t) will always produce an identity matrix
	// Variables
	int m_iQMatID;						// Unique identifier for the Q matrix
	string m_sName;						// Name for the process
	EDataType m_DataType;				// The datatype stored
	int m_iChar;						// The number of characters in the data
	int m_iChar2;						// The number of characters in the data squared
	bool m_bAllowModelUpdate;			// Whether the model can be updated
	double *m_ardQMat;					// The Q matrix
	double *m_ardU,*m_ardV,*m_ardRoot;	// The eigen information
	double *m_ardRootEqm;				// The square root of the eqm distribution (for fast eigen decomposition);
	// Memory management functions
	int GetQMatID()	{ return QMatID++; }									// Get the identifier for a Q matrix
	void MakeSpace(int Char,EDataType Type,string Name);
	void MakeSpace(EDataType Type,string Name)	{ MakeSpace(NumStates(Type),Type,Name); }
	void CleanSpace();
	virtual ostream &Output(ostream &os);									// Output function
};

///////////////////////////////////////////////////////////////////////////
// This is the Q matrix class for covarion models
// --
// The only difference between this and the class CQMat is how scaling
//  is performed
///////////////////////////////////////////////////////////////////////////
class CTHMMQMat : public CQMat	{
public:
	// Constructor
	CTHMMQMat(int DataChar, int m_iChar, EDataType Type,string Name = "Unnamed Q");
	// Destructor
	~CTHMMQMat() { /* empty */ };
	// Interaction functions
	double OverallSubRate(bool FWP);		// Returns the overall rate of the *SUBSTITUTION* process -- ignores transitions between hidden states
	double OverallTransRate();				// Returns the overall rate of transitions in the *HIDDEN* process -- ignores substitions between observable states
	double TransRateFrom(int From);			// Returns the overall transition rate from a state
	double TransRateTo(int To);				// Returns the overall transition rate to a state
	void ScaleQ(double Rate);				// Scales only the substitution process in Q
private:
	int m_iDataChar;							// Size of the data (used for calculating
};

///////////////////////////////////////////////////////////////////////////
// Eqm defining classes - Used for applying a defined Eqm to a CQMat
///////////////////////////////////////////////////////////////////////////
// equ			= equiprobable (e.g. Jukes and Cantor)
// obs			= single frequencies produced from data
// complex		= seperate frequency vectors for each hidden process
// complex_GC	= single frequency vector with different GC content for each process

class CBaseEqm	{
public:
	// Core functions
	CBaseEqm(int Char, vector <CQPar *> *ProcPar);			// The constructor function
	~CBaseEqm();											// The destructor function
	// Interaction functions
	void ApplyEqm2QMat(double *Q, int MatID);				// Applies the equilibrium distribution to these matrices (based on m_viQMatID)
	bool IsID(int ID);										// Whether the equilibrium distribution applies to a particular matrix
	void AddMatID(int ID);									// Add matrix ID that equilibrium applies to
	void SetOpt(bool Optimise);								// Whether to optimise the Eqm distribution
	int Char() { return m_iChar; }							// Returns the number of characters in teh eqm distribution
	virtual vector <double> Eqm() {							// Returns the equilibrium distribution for the whole process
		vector <double> v; Error("Empty eqm distribution"); return v; }
	// HMM based eqm functions
	virtual vector <double> SubEqm(int HiddenState) {		// Returns the equilibrium distribution for a substitution state
		vector <double> v; Error("Empty CBaseEqm::SubEqm() function"); return v; }
	virtual vector <double> TransEqm() {					// Returns the equilibrium for the transition matrix
		vector <double> v; Error("Empty CBaseEqm::TransEqm() function"); return v; }
	vector <double *> OptimiserValues();					// Returns a vector containing the optimiser values for the eqm distribution
	void SetDoBasicNonReversibleCovarion(bool Val) { DoBasicNonReversibleCovarion = Val; }
	virtual void Shuffle()	{};										// Blank shuffle routine for when there are multiple eqms
	virtual void ResetEqm(vector <double> Eqm, bool RandomFactor) { Error("Empty CBaseEqm::ResetEqm function...\n"); };// Resets the eqm parameters to the values in the vector
protected:
	// Variables
	int m_iChar;											// Number of character states in the equilibrium
	vector <int> m_viQMatID;								// The CQMat objects that the eqm process is applied to
	vector <CPar *> m_vpEqmPar;								// Pointers to the parameters that define the equilibrium distribution
	vector <CQPar *> *m_pProcPar;							// Pointer to process parameters
	double *m_ardQModifier;									// The modifier matrix (element-wise)
	// Elements redundant for some kinds of models
	bool DoBasicNonReversibleCovarion;					// Produce a model that relies on \pi^{k}_i == \pi^{l}_i i.e. standard covarion model
	// Checks as to whether a recalculation of m_ardQModifier is required
	vector <double> m_vdOldParVals;							// Old parameters values (used for whether updating is required)
	bool CheckQ();											// Checks whether the equilibrium distribution has changed
	// The function that actually applies the distribution
	virtual void AdapterFunc() { Error("Empty adaptor function"); }	// Creates m_ardQModifier matrix for the equilibrium distribution

};

/* ********************** Class describing a simple eqm distribution *********************** */
class CSimpleEqm : public CBaseEqm {
public:
	// Core functions
	CSimpleEqm(int Char, vector <CQPar *> *ProcPar, CData *Data);			// Constructor
	CSimpleEqm(int Char, vector <CQPar *> *ProcPar, vector <double> Eqm);	// Constructor
	~CSimpleEqm();															// Destructor
	// Interaction
	void ResetEqm(vector <double> Eqm, bool RandomFactor);					// Resets the equlibrium distribution
	vector <double> Eqm();									// Returns equilibrium distribution
private:
	void ConstructSimpleEqm(vector <double> eqm);			// Construct the equilibrium
	void AdapterFunc();										// Applies the simple eqm distribution to Q
};

/* ************************ Class describing codon equ distribution ************************* */
class CCodonEqm : public CBaseEqm {
public:
	// Core functions
	CCodonEqm(ECodonEqm CE, int GenCode, vector <CQPar *> *ProcPar,CData *D);			// Constructor
	~CCodonEqm();																// Destructor
	// Interaction
	ECodonEqm m_EqmType;														// Stores the eqm type
	int m_iGenCode;																// Store the genetic code
	vector <double> Eqm();														// Returns eqm distribution
private:
	ECodonEqm m_CE;
	void ConstructCodonEqm(vector <double> eqm);								// Construct the eqm
	void AdapterFunc();															// Applies the eqm distribution to Q
	CData *m_pData;																// Pointer to data
};

/* ************************ Class describing covarion eqm process *************************** */
class CCovEqm : public CBaseEqm {
public:
	// Core functions
	// Two constructors to create simple covarion equilibriums (1 observable distribution for all hidden states)
	CCovEqm(int Char, vector <CQPar *> *ProcPar, CData *Data, vector <CPar *> *StateProbs);			// Constructor
	CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <double> Eqm, vector <CPar *> *StateProbs);	// Constructor
	// Constructor functions to create complex covarion equilibriums (>1 observable distribution for all hidden states)
	CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <vector <double> > Eqm, vector <vector <int> > EqmMap, vector <CPar *> *StateProbs, bool AllowOpt = false);
	// Destructor function
	~CCovEqm();															// Destructor
	// Interaction
	vector <double> Eqm();								// Returns equilibrium distribution for the whole process
	vector <double> SubEqm(int HiddenState);			// Returns the equilibrium distribution for a specific hidden state (size == m_iDataChar);
	vector <double> TransEqm();							// Returns the equilibrium distribution for the C matrix (describes transitions between hidden states; size == m_iNoCovState)
	void Shuffle();										// Shuffle the states that the eqm applies to
	void ResetEqm(vector <double> Eqm, bool RandomBit);	// Resets the eqm
private:
	// General variables
	int m_iNoCovState;								// Number of covarion states
	int m_iDataChar;									// Nunber of characters in the data (m_iChar == m_iNoCovStates * m_iDataChar)
	vector <CPar *> m_vpStateProbs;						// Probability of each of the m_iDataChar states
	// Variables relating to observed equilibrium distributions (i.e. those of characters)
	int m_iNoCovObsEqm;									// Number of different equilibrium distributions for the hidden-states
	vector <vector <int> > m_vviEqm2State;				// Mapping of these eqm distributions onto hidden states
	// Functions
	void ConstructSimpleCovEqm(vector <double> eqm,vector <CPar *> *StateProbs);	// Construct the equilibrium
	void AdapterFunc();									// Applies the simple eqm distribution to Q
};

///////////////////////////////////////////////////////////////////////////
// Class defining the space for an individual site at a node
///////////////////////////////////////////////////////////////////////////
class CSite {
public:
	CSite(int *Char);						// Basic constructor
	CSite(const CSite &Site);						// Copy constructor
	~CSite();								// Basic destructor
	// Interaction functions
	inline void Overwrite(CSite &Site,int SiteNum);	// Function that makes a site a copy of another
	inline void CopyVals(CSite *Site, bool ForceReal = false);	// Copies values in Site->m_pSpacePointer and Site->m_pScalePointer to this
	inline bool IsReal()	{ return m_bReal; }		// Whether computations should be performed for a site
	inline double *Sp() { return m_pSpacePointer;}	// Returns the space for the site
	inline double *RealSp() { return m_ardSpace; }	// Returns the real space regardless (for derivative calculations)
	inline int *Sc()	 { return m_pScalePointer;}	// Returns the space for the scaling factor
	inline int *RealSc() { return &m_iScale; }		// Returns the real scaling factor regardless
	inline int OriSite()	{ return m_iOriSite; }	// Returns the site a sequence was based on
	inline int Copies()	{ return m_iCopyNum; }	// Returns the number of copies of a site
	inline void ResetSite();						// Resets the site to a real site
	inline void ZeroSite(bool ForceAll = true);	// Sets all values in a site to zero
	friend ostream &operator<<(ostream &os, const CSite &Site);	// Output function
private:
	// Variables
	int *m_Char;							// The size
	double *m_pSpacePointer;				// Pointer to the space used
	double *m_ardSpace;						// The space for the node calculations
	int m_iScale;							// Value holding the scaling factor
	int *m_pScalePointer;					// The Scaling factor
	int m_iOriSite;							// if(m_bReal == false) The site it has been copied from
	int m_iCopyNum;							// if(m_bReal == true) the number of copies of the site
	bool m_bReal;							// Whether the node is actually real
	CSite *m_pOrSite;						// Pointer to original site
};

// This is the base process class
/////////////////////////////////////////////////////////////////////////
// Provides basic functionality for all likelihood calculations
// It will include information held by the individual processes, including:
//
//	- Q matrix information
//  - Partial likelihood information
//
// The process will
//	- Scale the process to mean rate 1 or leave it unscaled (for special cases)
//
// The model will expect to be able to:
//  - call the process calculation functions and receive back a final vector of likelihoods

class CBaseProcess {
public:
	// Constructor
	CBaseProcess(CData *Data, CTree *Tree,string Name = "Unnamed process");
	// Destructor
	~CBaseProcess();

	////////////////////////////////////////////////////////////////////////////////////////
	// Main interaction functions for calculations and so on
	////////////////////////////////////////////////////////////////////////////////////////
	void PrepareFastCalc(vector <int> *C = NULL);					// Prepares the columns for fast calculations from vector (if blank it makes its own
	void DecompressSpace();																// Reverts to normal memory allocation

	bool PrepareLikelihood(bool DoQ, bool ForceRemake);		// Prepare for the main likelihood function
	bool CreatePTMats(int Bra = -1);				// Creates the PT matrices ready for for branch [Bra]; if(Bra == -1) Do all branches...
	bool Likelihood(bool ForceReal = false);		// Main likelihood function
	virtual double Penalty() { return 0.0; }		// Returns a penalty value for penalized likelihood (if appropriate)
	double LogL();													// Calculate log likelihood of process (used after Likelihood is called);
	CProb &Lsum(int Site);						// Get the sitewise likelihood

	// Data and tree memory manipulation routines
	void ApplyNewData(CData *D, bool RedoSpace = true);	// Apply new data to the process (RedoSpace will delete existing space and replace it)
	void ApplyNewTree(CTree *T);						// Apply new tree to the process (RedoSpace will delete existing space and replace it)
	void CleanCalcSpace();											// Clean space used by calculatios
	void MakeCalcSpace(bool AllowRemake = false);					// Make the calculation space (or remake it if necessary) from character size

	// Functions for getting branch derivatives
	virtual void PrepareBraDer();						// Prepare the process for get branch derivatives function
	bool GetBraDer(CProb *ModelL);				// Gets the branch derivatives. ModelL = array of partial likelihoods for whole model


	// Parameter access functions
	////////////////////////////////////////////////////////////////////////////////////////
	CQPar * Kappa(EDataType Type);					// Adds kappa to a DNA (4x4) or Codon (64x64) process
	int Size() { return m_iSize; }						// Returns the size
	int PatOcc(int i) { assert(m_pData != NULL); assert(InRange(i,0,m_pData->m_iSize)); return m_pData->m_ariPatOcc[i]; }
	inline string Name()	{ return m_sName; }
	inline int NoPar()		{ return (int)m_vpPar.size(); }			// Returns the number of parameters
	inline int Char() { return m_iChar; };				// Returns the number of states in the model
	inline int HiddenChar() { return m_iHiddenChar; }	// Returns the number of hidden states in the process
	inline int NoSeq() { return m_pData->m_iNoSeq; }	// Returns the number of sequences in the process
	inline int DataChar() { return NumStates(m_DataType); }	// Returns the data character
	bool PseudoProcess() { return m_bPseudoProcess; }	// Returns whether process is a pseudoprocess
	CQPar *RatePar() { return m_pRate; }			// Returns a pointer to the rate parameter
	CQPar *ProbPar() { return m_pProcProb; }		// Returns a pointer to the probability parameters
	CPar *pPar(int i){ assert(i<NoPar()); return m_vpPar[i]; }	// Return a parameter
	CPar *AddQPar(CQPar *P) { assert(P!=NULL); m_vpPar.push_back(P); return P; }	// Adds a parameter
	CTree *Tree() { if(IsSubTree()) { return m_pSubTree; } return m_pTree; }	// Returns a pointer to the tree
	vector <CBaseEqm *> EqmPar() { return m_vpEqm; }			// Returns the Eqm parameters
	vector<double> Eqm(int QMatNum) { return m_vpQMat[QMatNum]->Eqm(); };	// Returns the numerical value of equilibrium
	virtual vector <double> RootEqm();				// Returns the root equilibrium for likelihood calculations
	void OutEqm(int QMat = -1, ostream &os = cout);	// Returns the eqm distribution
	double Rate(double NewRate = -1.0,bool MakeRateOpt = false);// Returns the rate and if NewRate >= 0 sets m_pRate = NewRate
	bool MaxRate() { return m_bMaxRate;}						// Whether the process is meant to be max rate
	bool MakeMaxRate() { m_bMaxRate = true; }					// Enforces max rate
	bool MakeNormalRate() { m_bMaxRate = false;	}  				// Goes back to a 'normal' rate
	CQPar *AddRatePar2Opt();									// Adds rate parameter to the optimised parameters
	void RemovePar(string Name);						// Removes a parameter from the process

	CTree *MainTree(CTree *T = NULL);							// Pointer to underlying tree
	CData *MainData() { return m_pData; }						// Pointer to underlying data
	inline bool IsSubTree() { if(m_pSubTree == NULL) { return false; } return true; }
	inline bool IsCompressed() { return m_bCompressedSpace; }

	// Functions dealing with the processes probability
	double Prob() { return m_pProcProb->Val() / *m_piProbFactor; }					// Returns the probability of the process
	void SetProbFactor(int i) { *m_piProbFactor = i; }

	// Functions for shuffling parameters
	virtual void ShuffleParameters() { cout << "\nDoing an empty shuffle..."; exit(-1); };				// Useful for mixture models where parameters combinations cause problems...

	// Some useful functions for debugging
	vector <double> GetQMatSums();						// Returns the sum of the off diagonal entries in the Q matrix (useful for detecting silly changes)

	// Functions for getting copies of the process (used in model.cxx)
	// These *MUST* be overloaded when other virtual functions are used
	virtual CBaseProcess *RateProcessCopy();								// Returns a pseudoprocess with a seperate rate category
	virtual CBaseProcess *GammaRateProcessCopy();							// Returns a pseudoprocess with a seperate rate category and tied probability (used for gamma distribution)

	// interaction functions for subtree assignment and calculations
	////////////////////////////////////////////////////////////////////////////////////////
	void ApplyCPMapping(vector <int> LeafMap, vector <int> NodeFr, bool DoPartial = true);	// Functions that applies CP mapping
	void CleanCPMapping();												// Function that cleans everything up and removes the partial likelihood mappings
	// Functions that allocate for a specific tree
	void ApplySubTree(CTree *Tree);							// Applies a subtree for likelihood computation
	void CleanSubTree();									// Clears subtree from likelihood computation
	bool OldGetBraDer(CTree *Tree, bool Init = true);			// Calculate branch derivatives (returns whether it was successful)
	// Function that delivers the likelihood of a site
	CProb &L(int Site)	{ assert(IsProb(Prob())); return m_ardL[Site].Multiply(Prob(),false); }; 			// The final likelihoods
	CProb &ModelL(int Site) { return m_arModelL[Site]; }
	bool LOkay() { if(MATCH_PAML == 1) { return true; } return FlipBool(m_bFailedL); };			// Returns whether the likelihood has computed
	bool IsGamma() { return m_bIsGamma;	}					// returns whether the process is gamma distributed
	void MakeGamma() { m_bIsGamma = true; }					// Flags that the process has a gamma distributed rate associated with it
	// Miscellaneous space access and output functions
	void OutQ(int QNum = -1,ostream &os = cout);								// Output the Q matrices for the process
	ostream &SiteOut(ostream &os, int Node, int PosBeg = -1, int PosEnd = -1);	// Output function to output details of a site
	void OutPT(ostream &os, int Branch);										// Output the PT for a branch
	void ZeroSpace();															// Function that removes all values currently held in the space
	vector <double> AggregateQ(int Position);				// Returns the aggregated Q matrix for the Position-th character
//	vector <double> AggregatePT(int Position);				// Returns the aggregated P(t) matrix for the Position-th character

	ostream &SiteL(int Site, ostream &os = cout);								// Output likelihood of a site
	void SetOutputDetail(bool V) { m_bDetailedOutput = V; }						// Set level of detail in general output
	void OutPartL(int Site, int Node, ostream &out = cout);						// Formatted output of partial likelihood for a node and site
	vector <double> GetPartL(int Site, int Node);								// Returns the partial likelihood for a node and a site
	int GetPartLScale(int Site, int Node);										// Returns the partial likelihood scale factor for a node and a site
	vector <double> GetPT(int Branch);											// Returns the PT matrix
	// Space updating routines
	virtual bool Make_PT(int Branch, bool RedoRate = false);				// Makes the P(t) matrix for the likelihood function
	void ScaleQ() { int i; FOR(i,(int) m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); } };	// Scale Q matrices ready for calcs
	void LeafNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoCompleteUpdate = false);
	void BranNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoCompleteUpdate = false);
	// Space functions that calculate partial likelihoods from a centre point from SNAP paper
	void PreparePartialL(int Node, int NodeFrom,int Node2Trans=-1);	// Prepare a partial likelihood: if(Node2Trans != -1) { transfer space info to Node2Trans space; }
	// Function for FastBranchCalc
	void GetBranchPartL(CProb **P, int NTo, int NFr, int Br);	// Get the sitewise likelihoods
	// Reset the calculation space...
	void ResetCalcSpace();
protected:
	// Critical variables
	int m_iChar;						// Number of characters in the process
	int m_iHiddenChar;					// Number of hidden states in the process (only matters for THMMs)
	int m_iDataChar;					// Number of characters in the observable data (m in the model)
	int m_iSize;						// Length of the sequences
	int m_iSiCh;						// Is m_iChar * m_iLength
	int m_iChar2;						// Is m_iChar * m_iChar
	string m_sName;						// Name of process
	vector <CPar *> m_vpCovProbs;				// The probability parameters of each hidden state (size == m)
	EDataType m_DataType;				// The type of data the process is meant to apply to
	int m_iSpaceSize;					// Size allocated in space
	int m_iSpaceNoSeq;					// NoSeq allocated in space
	// Some flags of various things
	bool m_bPseudoProcess;				// Whether the process is real or not (i.e. needs no recomputing of QMat and m_vdEqm and only scales)
	bool m_bFailedL;					// Likelihood computations failed or is 0.0
	bool m_bIsGamma;					// Whether the process has a gamma distributed rate associated with it.
	bool m_bBraDerReady;				// Flag to check whether branch derivative calculations were ready
	bool m_bDetailedOutput;				// Flag as to whether to do detailed output
	bool m_bModelPerBranch;				// Flag highlighting whether parameters need estimating per branch
	bool m_bDoStandardDecompose;		// Whether eigen vectors/values for the model are produced in the usual manner

	////////////////////////////////////////////////////////////////////////////////////////
	//			Parameter related definitions - Including some standard model stuff
	//			Mostly used for preparing and defining models
	//
	// ********* NOTE: Parameters are held in two arrays: m_vpPar and m_vpEqmPar ***********
	////////////////////////////////////////////////////////////////////////////////////////

	// Some Q related values including parameters
	vector <CQMat *> m_vpQMat;			// Vector of the Q matrices used by the process (1 for most, but many when process varies across branches)
	vector <double> m_vdEqm;			// Pointer to the equilibrium distribution at the root (as used for likelihood computations)
	int m_iRootQ;						// The Q matrix of the process at the root
	string m_sABET;						// The string defining the processes alphabet
	vector <int> m_viQ2Bra;				// Vector describing which Q matrices describe which branches (if empty then use m_vpQMat[0] for all branches
	int m_iABET_length;					// The length of each character in the alphabet
	vector <CQPar *> m_vpPar;			// The parameters in the process
	CQPar *m_pRate;						// The rate parameter -- Usually set to one, but varies when rate varies between processes
	bool m_bMaxRate;					// Whether the process is meant to be max rate;
	CQPar *m_pProcProb;					// The probability of the process occurring
	int *m_piProbFactor;				// The number of processes a whole probability is divided between (important for gamma rates)
	double m_dBaseVal;					// The value that a Q is initialised to (doesn't matter is its scaled...)
	int m_iProcID;

	// Functions for creating Q matrices

	virtual bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);			// Prepares the Q matrices in a process
	virtual bool PrepareQMats(bool DoScale = true)								// Overloaded function for doing all the Q matrices
				{ vector <int> v; return PrepareQMats(v,DoScale); }
	virtual void DoPostHocQMatUpdate() {};								// Function available if necessary

	// Variables and functions relating to the eqm distribution of Q matrices
	vector <CBaseEqm *> m_vpEqm;										// The equilibrium functions
	vector <double> SimpleEqm(int QMatID,bool AllowSumNotOne = false);	// Function returning a vector of the eqm distribution for a Q matrix

	// Functions relating to parameters and Q matrix
	void CleanPar();					// Removes all the parameters in m_arpPar and their mappings
	void CleanQ();						// Removes all Q matrices and their mappings
	void MakeABET(EDataType Type);		// Define it as a standard type
	int GetProcID() { return ProcID++; }// Gets the unique process ID;

	// Functions for adding Q matrices to the model
	virtual CQMat *Add_QMat(string Name,EDataType Type);					// Function for adding simple Q matrices based on standard data types
	virtual CQMat *Add_QMat(string Name, int Char);							// Function for adding Q matrices of size Char
	virtual CQMat *Add_CodRedQMat(string Name, int Char);			// Function for adding Q matrices for codon data with the genetic code enforced
	// Functions for mapping Q matrices back to the likelihood function
	int QMat4Bra(int Branch);										// Returns the number of the QMat object for a specific branch

	// The variables carrying the information about the process
	// Uses structure CSite for vector space
	bool m_bCheckSpace;					// Checks that the space has been initialised
	double *m_ardPT;					// The PT matrices
	double *m_ardQP;					// The Q * P(t) matrices

	bool m_bCompressedSpace;			// Whether the space has been compressed or not
	bool m_bDoingPartial;				// Set to true when making a partial likelihood (stops anything being transferred to m_iStartCalc);
	vector <CSite> m_vSpace;			// Space used for calculations
	vector <CSite> m_vBackSp;			// Space used for backwards calculations (for branch derivatives)

	CProb *m_ardL;						// The final likelihoods
	CProb *m_arModelL;					// Pointer to the final model likelihoods (set in PrepareBraDer and calculated in CBaseModel::Sitewise
	// Space access helpers
	inline int InitNodePos(int Node)	{ return (Node * m_iSize); }
	inline int PartLNode()				{ return m_pTree->NoNode(); }

	// Functions that define access to calculation space
	void CleanSpace();												// Cleans all the currently used space
	void MakeBasicSpace(int Char);									// Make the basic properties of the process
	void MakeBasicSpace(EDataType Type);							// Make basic properties under a certain standard data type

	// Space access functions
	inline double *PT(int Br) { return &m_ardPT[m_iChar2 * Br]; };					// The PT matrix for a branch
	inline double *QP(int Br) { return &m_ardQP[m_iChar2 * Br]; };					// The Q * PT matrix for a branch
	inline double *Fd(int Node,int Pos) { return QkFd(InitNodePos(Node) + Pos); };	// Get partial likelihood vector for specific node
	inline double *QkFd(int NodePos) { return m_vSpace[NodePos].Sp(); };			// Quick access function to Forward space

	inline bool FdReal(int Node, int Pos) {return QkFdReal(InitNodePos(Node)+Pos);};// Returns whether forward space needs calculating
	inline bool QkFdReal(int NodePos) { return m_vSpace[NodePos].IsReal(); };		// Quick checker to see if space needs calculating
	inline bool BkReal(int Node, int Pos) {return QkBkReal(InitNodePos(Node)+Pos);};// Returns whether back space needs calculating
	inline bool QkBkReal(int NodePos) { return m_vBackSp[NodePos].IsReal(); };		// Quick checker to see if back space needs calculating
	inline int OrFd(int NodePos) { return m_vSpace[NodePos].OriSite(); }			// Return the original label for a node in Fdspace
	inline int OrBk(int NodePos) { return m_vBackSp[NodePos].OriSite(); }			// Return the original lable for a node in the BkSpace
	inline int FdCopies(int Node, int Pos) { return QkFdCopies(InitNodePos(Node) + Pos); }	// Returns the number of copies of a site
	inline int QkFdCopies(int NodePos) { return m_vSpace[NodePos].Copies(); }				// Returns the number of copies of a site
	inline double *Bk(int Node,int Pos){ return QkBk(InitNodePos(Node) + Pos); };	// Get the Backwards likelihood vector for specific node
	inline double *QkBk(int NodePos) { return m_vBackSp[NodePos].Sp(); };			// Quick access function to backward space
	inline int *FdScale(int NodeNum, int Pos) { return QkFdSc(InitNodePos(NodeNum)+Pos); };	// Access to the forwards space held in m_ardSpace
	inline int *BkScale(int NodeNum, int Pos) { return QkBkSc(InitNodePos(NodeNum)+Pos); };	// Access to the backwards space held in m_ardBackSpace
	inline int *QkFdSc(int NodePos) { return m_vSpace[NodePos].Sc(); }; 			// Quick forward scale access
	inline int *QkBkSc(int NodePos)	{ return m_vBackSp[NodePos].Sc(); };			// Quick backward scale access
	// Force space functions
	inline double *ForceRealFd(int Node,int Pos) { return QkForceRealFd(InitNodePos(Node) + Pos); }	// Force use of real backward space (required for derivative calculations)
	inline double *QkForceRealFd(int NodePos) { return m_vSpace[NodePos].RealSp(); }
	inline int *ForceRealFdSc(int Node, int Pos) { return QkForceRealFdSc(InitNodePos(Node) + Pos); }
	inline int *QkForceRealFdSc(int NodePos) { return m_vSpace[NodePos].RealSc(); }
	inline double *ForceRealBk(int Node,int Pos) { return QkForceRealBk(InitNodePos(Node) + Pos); }	// Force use of real backward space (required for derivative calculations)
	inline double *QkForceRealBk(int NodePos) { return m_vBackSp[NodePos].RealSp(); }
	inline int *ForceRealBkSc(int Node, int Pos) { return QkForceRealBkSc(InitNodePos(Node) + Pos); }
	inline int *QkForceRealBkSc(int NodePos) { return m_vBackSp[NodePos].RealSc(); }

	void CopyNode(int NodeFr, int NTo);										// Copy all space in NodeFr to NTo
	void CleanScale(int NodeNum, bool ForceAll = true);						// Removes all scales in a node
	void TransScale(int NodeTo,int NodeFr,bool First, bool Partial = true, bool ForceReal = false);	// Transfer (sum) values from NodeFr to NodeTo
	void DoScale(int Node, bool ForceRealScale = false);					// Performs scaling on the fly; ForceRealScale = true will treat all nodes as real

	inline int *LScale(int Site) { return FdScale(PartLNode(),Site); }		// Final node and m_ardPerSitelnL scaling factor
	inline double *PartL(int Site) { return Fd(PartLNode(),Site); };		// The Partial likelihoods


	/////////////////////////////////////////////////////////////////////////////////////////
	// Variables and functions for calculating partial likelihoods upto subnodes, &c
	bool m_bSubTreeActive;									// Flag stating that the likelihood function is expecting to work with a subtree
	vector <int> m_viLeafMap;								// The mapping of the full data onto the subtrees nodes (0<i<m_iNoSeq) == Real sequences; -1 = Partial likelihood from Space i;
	CTree *m_pSubTree;										// The tree being used with the partial likelihood

	/////////////////////////////////////////////////////////////////////////////////////////
	// Likelihood computation functions
	// --------------------------------
	// These are more complicated than neccessary, but at least they work!

	// Normal calculations
	void PartialL(CTree *pTree, int iNoTo, int iNoFr, int Branch, bool LeafFirst = true, bool BlockWriteback = false);
	void BranchNodePartialL(int SpFrom, int SpTo, double *PT, bool First);
	void LeafNodePartialL(CTree *pTree, int LeafNode, int Branch, int SpFlag, double *PT,bool First);		// Usual leaf node calculation
	void MakeZeroRateLikelihood();				// Makes likelihood for Rate = 0 models
	void MakeMaxRateLikelihood();				// Makes likelihood for Rate = inf model (garbage collector models)

	// Branch derivative calculations
	void Branch_dT(int NTo, int NFr, int Branch, CTree *pTree, int First, int *BrError);
	void LeafNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	void BranNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	double PartialGrad(int site,double Total,int SiteScale);

	// Some pointers to other data based other stuff
	CTree *m_pTree;				// Pointer to the tree the model should be using
	CData *m_pData;				// Pointer to the data

	// Some virtual functions that will vary considerably within the calculations
	virtual void Data2PartL(int Char, double *PT, double *RetSpace, vector <double> *eqm);
	virtual void Vec_by_Data(int Char, double *Vec);		// Multiplies a vector by the characters (i.e. set all values not of type Char to 0); Returns the sum of these values
	virtual double Sum_Vec(int Char, double *Vec, vector <double> eqm);	// Returns the sum of elements * eqm matching Char
	void *Make_PT(double T,double *PT,double *U,double *V,double *Root);	// Makes the P(t) matrix													// Makes the Q matrix
	// Output functions
	friend ostream &operator<<(ostream &os, CBaseProcess &Proc) { return Proc.Output(os); }
	virtual ostream &Output(ostream &os);

	///////////////////////////////////////////////////////////////////////////////////
	// Some useful functions for building the equilibrium distribution of models
	// AddSimpleEqm adds a simple eqm distribution with only a single eqm distribution.
	void AddSimpleEqm(bool Opt = true);										// Add a simple eqm distribution from sequences
	void AddSimpleEqm(int Char, vector <double> Freq, bool Opt = true);		// Add a simple eqm distribution from Freq
	void AddSimpleCovEqm(vector <CPar *> *StateProbs, bool ForceEqu = false,bool Opt = false);						// Add a simple eqm for Covarion models
	void AddSimpleCovEqm(int Char, vector <double> Freq, vector <CPar *> *StateProbs, bool ForceEqu = false, bool Opt = false);	// Add a simple eqm for Covarion models
	// AddComplexEqm adds multiple eqm distributions
	// 1. Add seperate eqm for each hidden state. The initial eqm distributions are guessed from looking at the data
	void AddComplexCovEqm(int Char, vector <CPar *> *StateProbs, bool AllowOpt = true);
	// 2. Add seperate eqm for each hidden state. The initial eqm distributions are taken from Freq
	void AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <CPar *> *StateProbs, bool AllowOpt);
	// 3. Equilibrium distribution specified in Freq (the actual frequencies) and FreqMap (the list of states that each vector in Freq maps to)
	void AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <vector <int> > FreqMap,vector <CPar *> *StateProbs, bool AllowOpt); // Adds an equilibrium distribution with different probabilities for each CS state
	// AddCodonEqm adds a codon equilibrium distribution
	void AddCodonEqm(int GenCode, int Char, ECodonEqm CE, bool Opt);
};


///////////////////////////////////////////////////////////////////////////////////////
// Some basic process types that the processes can be assigned

enum DNAProc	{pRY,pJC,pFEL,pK2P,pHKY,pREV,pCOV,pDNA_THMM,pBRACOV,DNA_OTHER};
enum AAProc		{pEQU,pJTT,pWAG,pDAY,pMTREV,AA_OTHER};
enum CodonProc  {pM0};

class CDNAProcess : public CBaseProcess	{
public:
	CDNAProcess(CData *Data, CTree *Tree, DNAProc Model, string name);		// Constructor to produce a defaulted model
	~CDNAProcess();
	// Specific routines for building different types of model
	vector <CQPar *> MakeDNAREV(EDataType Type);		// Make a full reversible model
};

class CAAProcess : public CBaseProcess {
public:
	CAAProcess(CData *Data, CTree *Tree, AAProc Model, bool AddF = true);								// Constructor to produce a defaulted model
	CAAProcess(CData *D, CTree *T, string Name, bool AddF, double *S_ij, double *pi_j);			// Constructor to produce a general empirical model
	~CAAProcess();
	// General empirical model creation routine
	void CreateEMPmodel(double *S_ij,double *Freq,bool DoF);
	// Model specific routines
	void MakeEQU(bool AddF);
	void MakeDAY(bool AddF);
	void MakeWAG(bool AddF);
	void MakeJTT(bool AddF);
	void MakeMTREV(bool AddF);
};



class CCodonProcess : public CBaseProcess {
public:
	// Constructor function
	CCodonProcess(CData *Data, CTree *Tree, CodonProc Model,ECodonEqm CE,int GenCode = 0);
	// Destructor function
	~CCodonProcess();
	// Parameter specific bits and pieces
	CQPar * AddOmega(int GenCode);
	void AddMultiChangeZeros(int GenCode);
	// Functions for getting overall rates of certain types of events
	double ObsSynRate();					// Gets the observed rate of synonymous substitutions
	double ObsNonsynRate();					// Gets the observed rate of non-synonymous substitutions
	// Other interaction functions
	int GenCode() { return m_iGenCode; }
	// Output
	ostream &Output(ostream &os);
private:
	int	m_iGenCode;							// Stores the genetic code
};


#endif
