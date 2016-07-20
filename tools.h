///////////////////// Tools.h -- Tools header //////////////////
////////////////////////////////////////////////////////////////
//  Includes implementation of some useful tools, including
//  1) Eigen routines
//  2) Conversion and runtime routines
//  3) Basic parameters with internal bounding


#ifndef TOOLS_VAL_FILE
#define TOOLS_VAL_FILE
#define _CRT_SECURE_NO_DEPRECATE

#define DEVELOPER_BUILD 0

#include <cassert>
#include <iostream>
#include <map>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <ctype.h>
#include <string.h>
#include <string>
#include <sstream>
#include <cfloat>
#include <cctype>
#include <ctime>
#include <vector>
#include <istream>
#include <cmath>

using namespace std;

// Math functions that have trouble porting
int my_isnan(double x);
int my_isinf(double x);

#define DO_MEMORY_CHECK 0
#if DO_MEMORY_CHECK
class CMemChecker {
public:
	CMemChecker() {
		CountCBaseModel = 0;
		CountCBaseProcess = 0;
		CountCProb = 0;
		CountCPar = 0;
		CountCBaseEqm = 0;
		CountCQMat = 0;
		CountCSite = 0;
		CountCTree = 0;
		CountCNode = 0;
		CountCData = 0;
	}
	int CountCBaseModel;
	int CountCBaseProcess;
	int CountCProb;
	int CountCPar;
	int CountCBaseEqm;
	int CountCQMat;
	int CountCSite;
	int CountCTree;
	int CountCNode;
	int CountCData;
};
#endif

// Define which models exist
#define NoDNAModels 11
#define NoAAModels 8
#define NoCodonModel 1
enum EModel {	JC,FEL,K2P,HKY,REV,								// Standard nucleotide
				RY_model, COV_HKY,COV_REV,THMM_FULLDNA, HKYdG_THMM,	THMM_DNA, // Unusual nucleotide
				EQU,WAG,JTT,DAY,mtREV,							// Standard amino acid
				THMM_AA,WAGdG_THMM,								// Unusual amino acid
				Coevo_WAG,										// Coevolution models
				CodonM0, CodonM0_DrDc, CodonEMPRest, CodonEMPUnrest,			// Codon
				UNKNOWN};
const string sModelNames[] = {	"JC","FEL","K2P","HKY","REV",
						"RY","CovHKY",",CovREV","THMM_FULLDNA","HKYdG_THMM","THMM_DNA",
						"EQU","WAG","JTT","DAY","mtREV",
						"THMM_AA","WAGdG_THMM",
						"Coevo_WAG",
						"CodonM0", "CodonM0_DrDc","CodonEMPRest", "CodonEMPUnrest",
						"Unknown"};
enum ERateTypes {same,varyall,gammarates};			// Defines the type of rate variation occurring in THMMs
ostream &operator <<(ostream &os,EModel M);


#define MAX_SPACE 1000				// Maximum number of characters in model
#define BIG_NUMBER 100000000		// A large number
#define MIN_PROB 1.0E-3			// Minimum value a probability can take
#define MAX_CHAR 1000				// Maximum number of states in a process
#define FUNC_COUNTERS 1				// If 1 then number of matrix calculations is logged
#define MATCH_PAML 1			// Whether to match PAML likelihood calcultions or not
								// This affects:	i.  Columns of an alignment with a single character in
								//					ii. Cases where the likelihood fails (rather than -BIG_NUMBER, each site gets assigned 10^-80)
#define GRAD_MIN 1.0E-4			// Minimum value allowed for a gradient


const string RY_ABET = "RY";			// Purine/Pyrimidine or MatchMatch/MatchDelete or DeleteDelete/DeleteMatch
const string RY2_ABET= "RR RY YY YR";
const string DNA_ABET = "ACGT";
const string DNA2_ABET = "AAACAGATCACCCGCTGAGCGGGTTATCTGTT";
const string AA_ABET  = "ARNDCQEGHILKMFPSTWYV";
const string AA2_ABET = "AAARANADACAQAEAGAHAIALAKAMAFAPASATAWAYAVRARRRNRDRCRQRERGRHRIRLRKRMRFRPRSRTRWRYRVNANRNNNDNCNQNENGNHNINLNKNMNFNPNSNTNWNYNVDADRDNDDDCDQDEDGDHDIDLDKDMDFDPDSDTDWDYDVCACRCNCDCCCQCECGCHCICLCKCMCFCPCSCTCWCYCVQAQRQNQDQCQQQEQGQHQIQLQKQMQFQPQSQTQWQYQVEAERENEDECEQEEEGEHEIELEKEMEFEPESETEWEYEVGAGRGNGDGCGQGEGGGHGIGLGKGMGFGPGSGTGWGYGVHAHRHNHDHCHQHEHGHHHIHLHKHMHFHPHSHTHWHYHVIAIRINIDICIQIEIGIHIIILIKIMIFIPISITIWIYIVLALRLNLDLCLQLELGLHLILLLKLMLFLPLSLTLWLYLVKAKRKNKDKCKQKEKGKHKIKLKKKMKFKPKSKTKWKYKVMAMRMNMDMCMQMEMGMHMIMLMKMMMFMPMSMTMWMYMVFAFRFNFDFCFQFEFGFHFIFLFKFMFFFPFSFTFWFYFVPAPRPNPDPCPQPEPGPHPIPLPKPMPFPPPSPTPWPYPVSASRSNSDSCSQSESGSHSISLSKSMSFSPSSSTSWSYSVTATRTNTDTCTQTETGTHTITLTKTMTFTPTSTTTWTYTVWAWRWNWDWCWQWEWGWHWIWLWKWMWFWPWSWTWWWYWVYAYRYNYDYCYQYEYGYHYIYLYKYMYFYPYSYTYWYYYVVAVRVNVDVCVQVEVGVHVIVLVKVMVFVPVSVTVWVYVV";
const string COD_ABET = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT";
const string GAP_ABET = "-X.?";
const string RNA_ABET = "ACGU";
const string RNASTEM_ABET = "AAACAGAUCACCCGCUGAGCGGGUUAUCUGUU";


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Some data types and their related functions
// Definitions: RY=Purine/Pyrimidine; DNA; AA; COD=Codon(64 states); COD_RED=ReducedCodon(After genetic code);
// 				RNA = standard 4 state RNA; RNASTEM = Basepairs in RNA; RNAGAP
enum EDataType { RY, RY2, DNA, DNA2, AA, AA2, COD, COD_RED, RNA, RNASTEM, GAP, OTHER,NONE };
const double MIN_DATA_PERCENT=0.15;								// Minimum amount of sequence data matching datatype for program to run
ostream &operator<<(ostream &os, EDataType Type);				// Outstream operator
EDataType GetDataType(string s);								// Translate a string to EDataType
// Data type functions
int LenStates(EDataType Type);											// Returns the length of the character in 'string' units
int NumStates(EDataType Type);									// Get the number of states in a Datatype
string DataStates(EDataType Type);								// Get the string of a data type
string State(EDataType Type,int State);							// Get state State of a data type
int FindState(EDataType Type,string ToFind,int Start = 0);		// Find what state the ToFind is in a datatype
string GetPos(string Seq, int Pos, EDataType Type);				// Get position Pos of a sequence and return it as a string
string GetPos(string Seq, int Pos, int PosLength);				// Get position Pos of a sequence and return it as a string
int CodonDiff(string Codon1,string Codon2);						// Number of differences between two codons
// Other various useful functions
EDataType GuessDataType(string Data);									// Guesses the data type
double RateDataType(EDataType Type, string Seq);		// Scores the fit of data to a particular type (0 = no characters in common; 1 = all characters in common)
bool IsDataType(EDataType Type, string Seq, bool AllowGaps = false);	// Checks whether is of a certain dataype
bool IsDataType(EDataType Type, char Seq);								// Checks whether a character is of a certain datatype
bool IsTsByChar(char a, char b);										// Compares two chars and determines whether a transition
bool IsGap(char c);														// Is character c a gap? yes: true; no: false
bool IsTs(string a, string b, EDataType Type);							// Compares two states (inc. codons) and determines whether transition (only one change allowed)
vector <string> Tree2Names(string Tree,int NoSeq);						// Breaks up a tree and spits out names as a vector
vector <double> GetFreqFromSeq(string Seq, string ABET);				// Function to get state frequencies from a sequence (string vs string)
vector <double> GetFreqFromSeq(string Seq, EDataType DT);				// Function to get state frequencies from a sequence (string vs datatype)
vector <double> GetFreqFromSeq(vector <int> Seq, int iChar);			// Function to get state frequencies from a sequence (vector vs size)
// Various macro and function definitions
ostream& Error(string str = " ", ostream  &os = cout);
#define diff(a,b) fabs(a - b) > FLT_EPSILON?(1):(0)		// Checks difference to floating accuracy
#define tdiff(a,b,tol) fabs(a - b) > tol?(1):(0)		// Checks difference to tol accuracy
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define SHFT(a,b,c,d) (a) = (b); (b) = (c); (c) = (d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define pos(i,j,n)      ((i)*(n)+(j))
static double sqrarg;
static double maxarg1, maxarg2;
#define FMAX(a,b) ( maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
#define FOR(i,n) for(i=0; i<n; i++)				// Forwards for-loop
#define rFOR(i,n) for(i=n-1; i >= 0; i--)		// Backwards for-loop
#define IFOR(iter, list) for(iter = list.begin(); iter!= list.end(); iter++)

#define bIFOR(iter, list) for(iter = list.end(); iter != list.begin(); iter--);
#define pIFOR(iter, list) for(iter = list->begin(); iter != list->end(); iter++)
#define GET_MEM(var,type,size) var = new type[size]; assert(var != NULL);
#define DEL_MEM(var) if(var != NULL) { delete [] var; } var = NULL;
#define GET_LINE(file,var) var[0] = '\0'; while(var[0] == '\0' || var[0] == '\n') { file.getline(var,sizeof(var)); assert(!file.eof()); }
// Get no hash line, i.e. skips lines that start with '#'
#define GET_NHLINE(file,var) var = "#"; while(var[0] == '\r' || var[0] == '#' || var[0] =='\n' || var[0] == '\0') { getline(file,var); if(file.eof() == 1) { cout << "\nUnexpected end of file\n\n"; ex(-1); }}
#define FINOPEN(name,file) ifstream name(file); if(name.fail()) { cout << "\nCouldn't open file <"<<file<<">\n\n"; exit(-1); }
vector <string> Tokenise(string line);	// Tokenise a string
vector <string> Tokenise(string line, string Delim);		// Tokenise a string according to delimiter Delim
string EatWhiteSpace(string line);		// Remove all w/s from a string
string GetDataLine(ifstream *in);

// Functions taken from PAML relating to gamma distributed rates
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define MODEL_INITIAL_GAMMA 0.25
double PointChi2(double prob, double v);
double LnGamma(double alpha);
double IncompleteGamma(double x, double alpha, double ln_gamma_alpha);
double LBinormal(double h, double k, double r);
double PointNormal(double prob);
double CDFNormal(double x);
int DiscreteGamma(double freqK[],double rK[], double alpha, double beta,
	int K, int median);
// Some probability functions taken from PAML
double probBinomial (int n, int k, double p);
double Binomial(double n, int k, double *scale);
double probBetaBinomial (int n, int k, double p, double q);

////////////////////////////////////////////////////////////////////////////////////
// Other probability distributon stuff
double SampleExp(double lambda);						// Obtain a sample from the exponential distributions
double GetRiemannZeta(double a); 						// Calculate the mean of a RiemannZeta function numerically (lots of rubbish to screen)
double MeanZipf(double a, int MaxLength);		 		// What's the mean length from a zipf with parameter a, conditional on the max length being MaxLength
vector <double> ProbZipf(double a, int Maxlength); 		// Returns the probability distribution of a Zipf distribution truncated at Maxlength
int Zipf(double a, int MaxLength);						// Sample from a Zipf distribution (copied from DAWG)
	/*
// Other required functions
int cmatinv( complex x[], int n, int m, double space[]);
complex cexp (complex a);
complex cby (complex a, complex b);
*/

///////////////////////////////////////////////////////////////////////////////////////
// Linear algebra tools
///////////////////////////////////////////////////////////////////////////////////////
// Norms
double FrobeniusNorm(double * a,double * b, int n);
double FrobeniusNorm(vector <double> &a,vector <double> &b);
// Function that checks detailed balance equations for reversibility
bool CheckReversibility(vector <double> Mat, vector <double> Eqm, bool Output = true);
bool CheckReversibility(int n, double *Mat, double *Eqm,bool Output = true);
bool CheckReversibility(int n, double *Mat, vector <double> Eqm, bool Output = true);
// Checks whether a matrix is a probability matrix
bool IsProbMat(int n, double *M);
bool IsProbMat(vector <double> M);
// Adjusts a vector of probabilities so that none are true zeros
void FixProb(vector <double> *Prob);

// This implementation of eigen comes from PAML
int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0,
               double Root[], double U[], double V[]);
int eigenRealSym(double A[], int n, double Root[], double Offdiag[]);
// Matrix outputter
ostream &VecOut(int n, double *Vec, char delim = '\t',ostream &os = cout);
ostream &MatOut(int n, double *Mat, char delim = '\t',ostream &os = cout);
ostream &MatOut(int n,double **Mat,char delim = '\t', ostream &os = cout);
template <class MatType> ostream &MatOut(int n, vector <MatType> Mat, char delim = '\t', ostream &os = cout) {
	int i,j;
	assert(n*n == Mat.size());
	FOR(i,n) {
		os << endl;
		FOR(j,n) { os << Mat[(i*n)+j] << delim; }
	}
	return os;
}

// Standard matrix multipliers
void VMat(double *a,double *b, double *c, int n);					// a[1*SIZE], b[SIZE*SIZE], c[1*SIZE]  ......  c = a*b
void MulMat(double *a, double *b, double *c, int n, int m, int k);
void IMat(double *I, int n);		// Makes a n * n matrix an identity matrix
void IMat(double **I,int n);		// Makes a n * n matrix an identity matrix
void ZeroMat(double *Mat,int n);	// Sets a n * n matrix to zero

// Tools
bool FlipBool(bool V);								// Swap the value of a bool
bool FlipBin(int i);								// Swap the value of a binary integer

int GetOption(string message, int low, int high);	// Get a numerical option
bool FileExist(string File);						// Tests whether a file exists
bool GetYesNoOption(string message);				// Returns yes[1] or no[0] from a menu
string GetInFileName();								// Finds an existing input file
string GetOutFileName(string SuggestedFile = "\0");	// Finds or clears an output file

// Stream to string variables...
template<class T>
 string toString(const T& t, bool *ok = NULL)
 {
         ostringstream stream;     // line A
         stream << t;              // line B
         if(ok != NULL)
               *ok = (stream.fail() == false);  // line C
         return stream.str();      // Line D
}

// Reverse sort
template <class TSort> void RSort(vector <TSort> *List)	{
	int i,j;
	TSort temp;
	FOR(i,(int)List->size()-1)	{
		for(j=(int)List->size()-1;i<j;j--)	{
			if(List->at(j-1)<List->at(j)) {
				temp = List->at(j-1); List->at(j-1) = List->at(j); List->at(j) = temp;
}	}	}	}

// Normal sort
template <class TSort> void Sort(vector <TSort> *List)	{
	int i,j;
	TSort temp;
	FOR(i,(int)List->size()-1)	{
		for(j=(int)List->size()-1;i<j;j--)	{
			if(List->at(j-1)>List->at(j)) {
				temp = List->at(j-1); List->at(j-1) = List->at(j); List->at(j) = temp;
}	}	}	}

template <class TComp> bool Compare(vector <TComp> *List1, vector <TComp> *List2) {
	int i;
	bool same = true;
	if(List1->size() != List2->size()) { return false; }
	FOR(i,(int)List1->size()) { if(fabs((double) (List1->at(i) - List2->at(i))) > 1.0E-6) { same = false; break; } }
	return same;
}

template <class TSumType> TSumType Sum(vector <TSumType> *Vec) {
	TSumType Sum = 0;
	int i;
	FOR(i,(int)Vec->size()) { Sum += Vec->at(i); }
	return Sum;
}

template <class TQuantType> double Quantile(vector <TQuantType> Vec, double quantile) {
	int Pos;
	double Weight1;
	assert(quantile >= 0.0 && quantile <= 1.0);
	Sort(&Vec);
	Pos = (int) (Vec.size() * quantile); Weight1 = 1.0 - ((double)(Vec.size() * quantile) - (double) Pos);
//	cout << "\nGetting position: " << Pos << " weight " << Weight1 << " == " << Vec[Pos] * Weight1;
//	cout << "\nGetting position: " << Pos+1 << " weight " << 1-Weight1 << " == " << Vec[Pos] * (1.0 - Weight1);
	return (Vec[Pos]* Weight1) + (Vec[Pos+1] * (1 - Weight1));
}

template <class TOutVec> ostream& operator<<(ostream & os, vector <TOutVec> Vec) { int i; FOR(i,(int)Vec.size()) { os << Vec[i] << " "; } return os; }

template <class TRange> bool InRange(TRange Val, TRange LowerBound, TRange UpperBound) { return ( !(Val < LowerBound) && ( Val < UpperBound) ); }

template <class TIsIn> bool IsIn(TIsIn Val, vector <TIsIn> Vec) {
	if(Vec.empty()) { return false; }
	if(find(Vec.begin(),Vec.end(),Val) == Vec.end()) { return false; }
	return true;
}

// Displays a prompt for input, and assigns a default value if nothing is entered.
template <class T> void read(T &val, const std::string &prompt, const T &def_val)	{
	std::string inputLine;
	std::cerr << prompt << " [" << def_val << "] ";
	getline(std::cin,inputLine);
	if(inputLine.length() > 0)	{
		std::istringstream input(inputLine.c_str());
		if(! (input >> val)) { val = def_val; }
	}
	else val = def_val;
}
template <class T> void read(T &val, const std::string &prompt) {
	bool Flag = true;
	std::string inputLine;
	std::cerr << prompt;
	getline(std::cin,inputLine);
	while(Flag)     {
		if(inputLine.length() > 0)      {
			std::istringstream input(inputLine.c_str());
			if(input >> val) { Flag = false; }
	}	}
}


// String manipulations
string int_to_string(int num);
string double_to_string(double num);
void number(int a, char *ret);
string find_and_replace(string source, string find, string replace);
int wildcmp(const char *wild, const char *string); // Wildcard string comparison

// Vector manipulations
vector <double> NormaliseVector(vector <double> Vec, double VectorSumTo = 1.0);
int DiscreteRand(vector <double> *Probs);					// Get a random number (0, Probs.size()) with proportion Probs
double GetRMSD(vector <double> Obs, vector <double> Exp);	// Obs: pairwise distance; Exp: distance from tree

// Some useful definitions
#define MAX_PAR_VALUE 100		// fabs(MAX_PAR_VALUE) are the bounds for all parameters
#define SMALL_PROB 1.0E-3		// When a probability becomes small enough to be ignored

//////////////////////////////////////////////////////////////////////////////////////////
// Scaled probability class
//////////////////////////////////////////////////////////////////////////////////////////
// Simple class containing scaled probabilities - used for likelihood computations
// Contains probabilities in the form RealValue = m_dValue * pow(10,-m_iScale));
#define L_SCALE_VAL (1.0E-5)	// Value at which scaling takes place
#define L_SCALE_NUM 5
#define LOG10 2.3025850929940459		// Value of Ln(10) used in likelihood computations
inline bool Double_Zero(double value) { if(value > DBL_MIN) { return false; } else { return true; } }

class CProb {
public:
	CProb(double InitVal = 0.0);						// Constructor from double
	CProb(CProb &Prob);									// Constructor from CProb object
	CProb(double Val, int Scal);						// Constructor from Value and Scaling factor
	virtual ~CProb();											// Explicit blank constructor
	// Basic operators
	friend ostream &operator<<(ostream &os,CProb &Prob);// Output operator
	CProb &operator=(CProb &Prob);						// Copy operator
	CProb &Zero() { m_dValue = 0.0; m_iScale = 0; return *this; }	// Zeros the probability
	CProb &Assign(double Value);						// Assigns Value to the probability
	CProb &Assign(CProb &Prob);							// Assign Prob to the probability
	CProb &Assign(double Val, int Scal);				// Assign from Value and Scaling factor
	// Return value functions
	double Prob();										// Returns the raw probability (subject to underflow)
	double LogP();										// Returns the log of the probability
	bool IsZero();										// Returns whether the probability is zero or not
	// Partial return value functions
	int Scale() { return m_iScale; }					// Scale Value
	double ScalVal() { return m_dValue; }				// Returns the Scaled Value
	// Functions for numerical operations (if Overwrite is true it overwrites current CProb with new value
	CProb &Multiply(double Value,bool Overwrite);		// returns this probability multiplied by value
	CProb &Multiply(CProb &Prob,bool Overwrite);		// returns this probability multiplied by Prob
	CProb &Add(double Value,bool Overwrite);			// returns this probability added to Value
	CProb &Add(CProb &Prob,bool Overwrite);				// returns Adds Prob
	CProb &Divide(double Val, bool Overwrite);			// return this probability divided by value
	CProb &Divide(CProb &Prob, bool Overwrite);			// return this probability divided by Prob
	friend void VMat(CProb *a, double *b,CProb *c, int n);	// a[1*SIZE], b[SIZE*SIZE], c[1*SIZE]  ......  c = a*b
#if DEVELOPER_BUILD == 1
private:
#endif
	void DoScale();										// Performs the scaling
	void MatchScales(CProb *Prob, bool DoHigh);			// Matches the scales between two probabilities; this is dangerous as information can be lost)
	double m_dValue;									// The normal part of the probability
	int m_iScale;										// The 10 to the minus bit of the probability
};
bool IsProb(double V); // Check whether a double is a probability

//////////////////////////////////////////////////////////////////////////////////////////
// Scaled Value class
//////////////////////////////////////////////////////////////////////////////////////////
// Very simple class used in parameter class, that contains a value scaled to 1.0
#define PAR_TO_SCALE 1.0		// The value that parameters stop being scaled
class CScaledParVal {
public:
	bool Scaled() { if(fabs(m_dScaler-1.0) < FLT_EPSILON) { return false; } return true; }
	double Scaler()		   { return m_dScaler; }								// Returns the scaler
	double *pScaledValue() { return &m_dValue; }								// Returns a pointer to m_dValue
	double Value();						// Returns the real value
	// Set new value (returned as Value) by adjusting m_dScaler and fixing m_dValue
	double SoftNewVal(double Val); 	// Set new value (returned as Value) by adjusting m_dValue and fixing m_dScaler
	double NewVal(double Val); 		// Changes so that Val = m_dValue * m_dScaler
	double RescaleVal(double Val);	// Set value <- Val and  rescale so m_dValue = 1.
	CScaledParVal(double Val = 1.0,bool AllowScale = true) {						// Will automatically produce value scaled to 1.0
		m_bAllowScale = AllowScale;
		RescaleVal(Val);
	}
private:
	bool m_bAllowScale;
	double m_dValue;
	double m_dScaler;
};

//////////////////////////////////////////////////////////////////////////////////////////
// The parameter class
//////////////////////////////////////////////////////////////////////////////////////////

enum ParOp { REPLACE, MULTIPLY, DIVIDE, ADD, SUBTRACT };
ostream & operator<<(ostream &os, ParOp Par);
class CPar; // Predeclare class for functions

// Function used by CPar for probabilities
void ProbabilityScale(vector <CPar*> *Parameters, bool Value2Scale,bool first = false,bool ReorderProbs = false);

// Currently scaling is only required for probabilities
class CPar {
public:
	// Constructor routines
	CPar(string Name, double Value, bool Optimise=true, double Upbound = 0.0, double Lowbound=MAX_PAR_VALUE,ParOp Oper=MULTIPLY);	// Core parameter construction routine
	// Destructor routine
	virtual ~CPar();
	////////////////////////////////////////////////////
	// Interaction functions
	// 1. Value and scale related
	double Val()	{ UpdatePar(); return m_dRealValue; }				// Returns the real value of the parameter
	double DoOper(double OldValue);										// Apply the paramter's operator to OldValue
	double *OptimiserValue() { return m_dScaledValue.pScaledValue(); }	// Returns the value passed to the optimiser function
	double Scal()	{ UpdatePar(); return m_dScaledValue.Value(); }		// Returns the real scaled value
	double Scaler() { UpdatePar(); return m_dScaledValue.Scaler(); }	// Returns the scaling value;
	void Rescale(bool Force = false)	{ if(Force || (pDoUpdate == NULL && m_ardBounds[1] > 1.1)) { m_dScaledValue.RescaleVal(Val()); } }				// Rescale the value (but not probabilities)
	double SetVal(double Value, bool DoUpdate=true, bool Force=false, bool Normalise=true);	// Set the real value (Normalise only matters for probabilities or other tied parameters
	void SetBounds(double Low, double High) { assert(High+DBL_EPSILON > Low); m_ardBounds[0] = Low; m_ardBounds[1] = High; CheckBound(); }
	bool CheckBound(bool ForceBounds = true);						// Check the bounds of a parameter (ForceBounds will reset value if out of bounds)
	bool CheckLowBound(bool ForceBounds = true);					// Checks the lower bound of a parameter
	bool CheckUpBound(bool ForceBounds = true);						// Checkes the upper bound of a parameter
	double LowBound() { return m_ardBounds[0]; }					// Returns the lower bound
	double UpBound() { return m_ardBounds[1]; }						// Returns the upper bound
	double BoundDist();													// Get distance from bounds
	bool StoreOptBounds(double Low, double Up) { assert(Up + FLT_EPSILON > Low); m_ardOptBounds[0] = max(Low,m_ardBounds[0]); m_ardOptBounds[1] = min(Up,m_ardBounds[1]);  }
	double OptLow() { return m_ardOptBounds[0]; }
	double OptUp() { return m_ardOptBounds[1]; }

	// 2. The update function -- Probably the most important parameter function
	bool UpdatePar(bool force = false, bool RedoScale = false);		// Function to decide whether to update parameters
	bool NeedUpdatePar() { if(pDoUpdate == NULL) { return false; } return true; }
	bool IsProb() { if(pDoUpdate == &ProbabilityScale) { return true; } else { return false; } }
	// 3. Other general parameter related
	string Name(string Name = "\0") { if(Name != "\0") { m_sName = Name; } return m_sName; }
	// 4. Locks
	bool CheckLocked()	{ return m_bLockedScale; }
	bool Unlock();													// Unlock the parameter
	// 5. Optimisation related
	bool Opt()		{ return m_bOpt; }
	bool SetOptimise(bool Opt) { if(Special()) { m_bOpt = false; } m_bOpt = Opt; return Opt; }		// Set whether to optimise the parameter or not
	// 6. Other interactions
	bool Special()				{ return m_bSpecial; }
	bool SetSpecial(bool Sp)	{ m_bSpecial = Sp; return Sp; }		// Sets the special values
	void PreOptCheck();												// Function to be run before parameters are to be optimised
	void SetUpdate(void (*func) (vector <CPar*>*P,bool V2C,bool First,bool Other)) { pDoUpdate = func;  };
	bool SetAllowScale(bool Sc)	{ m_bAllowScale = Sc; return m_bAllowScale; }
	virtual void GlobalApply()	{ };								// Applies the parameter globally to all processes in a model (empty except for CTiePar);
	// 7. Access to the derivatives (does gradient checking and bounding dynamically)
	double grad() { return m_dGrad; };								// The 1st derivative of the parameter
	double grad(double g);											// Sets the gradient (does gradient checking and bounding dynamically)
	bool DoHardDer() { return m_bDoHardDer;} 						// Whether the derivative needs extra attention
	void SetHardDer(bool V) { m_bDoHardDer = V; }					// Sets whether do do a complex derivative calculation
	// Accessible parameters
	vector <CPar*> m_arpPar;	// Parameters whose values are linked to this one
	bool m_bOutDetail;			// Whether to output lots of detail about parameter
	virtual ostream &Output(ostream &os);			// Output function
	bool IsBranch() { return m_bIsBranch; }						// Returns whether a parameter is a branch
	bool SetAsBranch() { m_bIsBranch = true; return true; }		// Set parameter as a branch
	// Stuff for numerical derivatives
	bool DoNumDer()					{ return m_bNumDer; }
	bool NumericalDerivative()		{ m_bNumDer = true; }
	bool AnalyticDerivative()		{ m_bNumDer = false; }
	bool InitialiseDerivativeType()	{ m_bNumDer = FlipBool(IsBranch()); return m_bNumDer; }
	// overloaded Operators
	CPar &operator=(CPar &Par);
private:
	///////////////// Core values ////////////////////
	string m_sName;				// The parameters
	double m_dRealValue;		// The real value of the parameter used in modelling
	double m_dGrad;				// 1st derivative associated with parameter
	bool m_bNumDer;				// Whether the parameter is calculated using a numerical derivative
	bool m_bDoHardDer;			// Whether the numerical derivative should be calculated using three points rather than two (more accurate)
	ParOp m_Operator;			// What the parameter does
	CScaledParVal m_dScaledValue;	// The value of the parameter used in optimisation (see class definition above)
	double m_ardBounds[2];		// The [0] Upbound and the [1] Lowbound
	double m_ardOptBounds[2];	// The [0] lower and [1] upper bounds used in the last round of optimisation. Good starting points for next round
	bool m_bOpt;				// Whether to optimise the parameter
	bool m_bAllowScale;			// Whether Scaler is allowed
	bool m_bSpecial;			// Whether the parameter is categorised as special (often used for probability scaling)
	bool m_bLockedScale;		// The scale values become locked to prevent rescaling during optimisation
	bool m_bIsBranch;			// Whether the parameter is a branch
	///////// Other variables and functions //////////
	double m_dOldReal;			// Old RealValue of the parameter   } Used to check whether needs updating
	double m_dOldScale;			// Old ScaledValue of the parameter }
	// Update function(s) for the parameter (e.g. fixes probabilities and so on)
	void (*pDoUpdate) (vector <CPar*> *Parameters, bool Value2Scale, bool First, bool Other);
	friend void ProbabilityScale(vector <CPar*> *Parameters, bool Value2Scale,bool first,bool ReorderProbs);
};

ostream &operator<<(ostream &os, CPar &Par);

/////////////////////////////////////////////////
// Profile HMM stuff
string GetProfilePath(string Seq, string ABET);

#if FUNC_COUNTERS == 1
	void OutputCounters(ostream &os = cout);		// Output the counters
#endif

/* -----------------------------------------------------------------------
 * Name            : rngs.h  (header file for the library file rngs.c)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 09-22-98
 * -----------------------------------------------------------------------
 */

#if !defined( _RNGS_ )
#define _RNGS_

double Random(void);
void   PlantSeeds(long x);
void   GetSeed(long *x);
void   PutSeed(long x);
void   SelectStream(int index);
void   TestRandom(void);

#endif

int RandInt(int begin,int end);
inline double RandDouble(double Low, double High)	{ return Low + (fabs(High - Low) * Random()); }


#endif
