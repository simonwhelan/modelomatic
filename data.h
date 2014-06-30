////////////////////////////////////////////
// Data.h
////////////////////////////////////////////
// Contains all the Data class structures

#ifndef __DATA_HEADER
#define __DATA_HEADER

#include "Leaphy.h"

#define ALLOW_CODON_REDUCTION 1		// Whether the code allows reducing the states of codon models (used for debugging);

// Structure describing the data
// Code is located at the beginning of Data.cxx
// None is included for initialisations
class CData	{
public:
    // Constructor functions
	// Inputs data from a specific file starting at point Pos
	CData(string file, EDataType SpecType=NONE, bool AllowFail = false, streampos FilePos = 0);
	// Inputs data from a set of arrays (for pairwise distances)
	CData(int NoSeq,int Size, vector <string> InSeq, vector <string> InName,EDataType Type = NONE);
	// Blank data array
	CData(int NoSeq, int Size, EDataType Type, vector <string> *Names = NULL);
	// Takes strings into object
    void InputData(EDataType Type, vector <string> cInputSeq, vector <string> cInputName, vector <int> SiteLabels, bool AllowFail = false);
    // Takes a NEXUS file and it's tags into the object (Work Still in Progress...)
    void InputNexus(string File);
    // Destructor
	~CData();						// Clears memory
	void Clean();					// Does the memory clearing
    // Attributes
	bool Valid() {return m_bValid;}// Whether it is valid data
    int m_iNoSeq;					// Number of sequences
    int m_iSize;					// Length of sequences
	EDataType m_DataType;			// The type of data
	string m_sABET;					// The data's alphabet
	int m_iChar;					// Number of characters in the data
    vector <vector <int> > m_ariSeq;// Vector array of sequences[m_iNoSeq][m_iSize]
	vector <int> m_viNoChange;		// Vector of when no changes occur (0,m_iChar-1) for no changes of character; -1: change
    vector <int> m_ariPatOcc;		// Array of how often the patterns occurred (this is the one analysed
    int m_iName;					// Whether we've got names for the sequences
	vector <string> m_vsName;		// Vector of names
	vector <double> m_vFreq;		// Array of state frequencies[SIZE]
	vector <int> m_viSiteLabels;	// Used when there are multiple site types (MultiGene option)
	vector <vector <int> > m_vviCoevoMapping; 	// Used when the data is a Coevo model to allow mapping between data sets
	// Useful interaction functions
	vector <double> GetFreq(int Seq = -1);	// Get the distribution of frequencies for sequence Seq [-1: all sequences];
	string Seq(int SeqNum)	{ return m_vsTrueSeq[SeqNum]; }	// Returns the sequence
	// Functions to do data pattern stuff
	void OutDataDist(string FileName, EDataType Type);					// Output the distribution of character patterns to file <FileName> of type Type.
	void CreatePatterns(int Pos, vector <int> *iters,EDataType Type, vector <vector <int> > *pats);	// Recursive function to get the patterns
	// Bootstrapping functions
	vector <int> m_viRealPatOcc;	// Array of the original data's PatOcc;
	void Bootstrap();				// Bootstrap the data
	bool IsBootstrap() { return FlipBool(m_viRealPatOcc.empty()); }	// Whether the data is a bootstrap
	void CleanBootstrap();			// Cleans up the bootstrap data
	void GetGapMask(vector <vector <bool> > *Mask);	// Gets the mask of which sites in an alignment are gaps
	double MultilnL();				// Get the multinomial probability of the data
	// Other true data stuff
	vector <int> m_ariPatMap;				// m_ariPatMap[m_iTrueSize] = Mapping between position i in true sequence to pattern in analysed sequence.
									// May be used to obtain real data
	int m_iTrueSize;
	vector <string> m_vsTrueSeq;	// vector of real sequences
    // Implementation
	////////////////////////////////////////////////////////////
	void RemoveInvariantSites();	// Removes invariant sites from the alignment
	void RemoveSparseSeqs(bool sparse=true,CTree *Tree = NULL, bool do_stdout=true);		// Removes all sequences with no data (sparse=tree: without at least MIN_DATA_PERCENT% data); if Tree != NULL then will also remove the sequences from the tree
	// Distances
	double PropDiff(int S1, int S2, bool IgnoreGaps);	// Count differences between sequences S1 and S2;
	double PoissonDist(int S1, int S2);					// Calculate the poisson distance between S1 and S2
	double PoissonVar(int S1, int S2);					// Calculate the poisson variance between S1 and S2
    // file data input functions
    char *GetName(char *string, char *name);
	string GetName(string name);
	int RemoveSeq(int RemSeq);
	int RemoveSeq(int RemSeq,CTree *TREE);
	void CondenseGaps();			// Remove sites where all sequences are gaps
	void DNA2RY();					// Changes normal DNA into a 2 state RY form
	void CleanToDNACodon();			// Takes a DNA data set and cleans it up so it is codon compatible
	void AddColumn(vector <int> Column, int NoOccurs = 1);			// Add a new column to a data set, along with the number of times it occurs
	// Some useful functions for guessing parameter values
	double GetPar(vector <string> s1, vector <string> s2);
	double GuessKappa();					// Guesses Kappa for nucleotide models
	// Some functions for switching between different data types
	void Translate(int GenCode = 0);		// Translate DNA sequences into amino acid sequences using the specified genetic code
	void MakeCodonData();					// Changes DNA data into Codon data. No imposition of the genetic code
	void ReduceCodonData(int GenCode);		// Reduces the state space of codon data to a more manageable size
	void ExpandCodonData(int GenCode);		// Expands the state space of codon data back to a 64 state model
	bool GetCodonPositions(bool First, bool Second, bool Third);            // Rewrite the data structure with a subset of the codon positions. Requires the DNA formatted data
	int GenCode() { return m_iGenCode; }	// Returns the genetic code assigned to the data
	// Functions for getting likelihood scaling between amino acid and codon models, among other adjustments
	double GetAminoToCodonlnLScale(int GenCode, int *df);		// Decides what adjustment to do and returns the value
	double GetAdjustmentScore(CData *AA_Data, CData *COD_Data, vector <double> AAFreq, vector <double> CodFreq, int GenCode = 0); // Returns eqn (6) from Seo and Kishino (Syst Biol. 2008) for a set of frequencies
	double GetRYToCodonlnLScale(int GenCode, int *df);			// Returns the equivalent to eqn (6) from Seo and Kishino, Syst Biol. 2008, but for RY to codons
	double GetNT2MissinglnLScale();							// Returns the probability of observing some nucleotide data based purely on the frequencies (NB *NOT* the multinomial, but individual nucleotides)
	double OldGetAminoToCodonlnLScale(int GeneticCode); // Old amino acid version for error checking.
	// Some stuff relating to tree HMMs
	CData *MakeMatchData();				// Returns the data for a match state TreeHMM
	CData *MakeDeleteData();			// Returns the data for a delete state TreeHMM
	// Code for deailing with pairwise coevolution models
	CData * CreateAllPairs(CData *SecondData = NULL, bool Output=true);			// Function that creates all pairs of characters from existing data
	bool CreateAllPairMapping(CData *Data1, CData *Data2);
	void OutRealData(ostream &os = cout);										// Output the sequences in sequential
	int CountMSAChars();					// Function to count number of characters actually used in an analysis
private:
	bool m_bValid;			// Whether the data is valid
	int m_iGenCode;			// Stores the genetic code
};

// Overloaded ostream operator
ostream& operator<<(ostream& os, const CData& DATA);

// File input functions
bool ReadData(string File, vector <string > &Names, vector <string > &Seqs, bool AllowFail = true);		// Read data from FileName
// file data input functions
char *GetName(char *string, char *name);
string GetName(string name);

// Function that produces corrects a set of codon frequencies so that they match a set of amino acid frequencies
// For a codon c_{i,j} and the corresponding amino acid a_c_{i,j} we know
//   a_c_{i,j} = \sum_{c_{i,j} \in a_c_{i,j}} c_{i,j}
// For a given set of c_{i,j} they are normalised to produce the correct a_c_{i,j}
vector <double> EnforceAAFreqOnCodon(vector<double> CodFreq, vector <double> AAFreq, int GCode);

// Definitions of the genetic code
const int NumGenCode = 12;;
// The genetic codes are:
//	0:  Universal code
//	1:  Vertebrate mt
//	2:  Yeast mt
//	3:  Mould mt
//	4:  Invertebrate mt
//	5:  Ciliate nuclear
//	6:  Echinoderm mt
//	7:  Euplotid mt
//	8:  Alternative yeast nuclear
//	9:  Ascidian mt
//  10: Blepharisma nuclear
//	11: Fake code where everything codes
const string GenCodeName[] = {
		"Universal",					// [0]
		"Vertebrate mt",				// [1]
		"Yeast mt",						// [2]
		"Mould mt",						// [3]
		"Invertebrate mt",				// [4]
		"Ciliate nuclear",				// [5]
		"Echinoderm mt",				// [6]
		"Euplotid mt",					// [7]
		"Alternative yeast nuclear",	// [8]
		"Ascidian mt",					// [9]
		"Blepharisma",					// [10]
		"Fake with everything coding"	// [11]
};

const int GenCodes[][64] = {
	{	11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16,-1,15,-1,15,12, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16,-1,15,-1,15,12, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15,12, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,16,16,16,16,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16,15,15,15,15,12, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		 5,18, 5,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
	},{
		 2, 2,11, 2,16,16,16,16,15,15,15,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15, 4, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,15,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 7,15, 7,15,12, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
	},{
		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
		-1,18, 5,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
	},{
		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
}	};

#endif
