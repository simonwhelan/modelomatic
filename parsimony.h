////////////////////////////////////////////////////////////////
// Header file detailing the parsimony based algorithms
////////////////////////////////////////////////////////////////

#include "Leaphy.h"
#include "data.h"
#include "tree.h"

class CParsimony	{
public:
	// Constructors
	CParsimony(CData *D, CTree *T);										// Usual constructor
	CParsimony(EDataType Type, vector <vector <int> > *Data, CTree *T);	// Constructor for other data
	// Destructors
	~CParsimony();
	// Simple interaction
	int Size()  { return m_iSize; }
	int NoSeq() { return m_iNoSeq; }
	CTree *Tree() { return m_pTree; }
	// Scoring functions
	int Score(bool FullScore = false);					// Get parsimony score; FullScore forces calculation without a mask
	vector <int> SitewiseScore(bool PenaliseGaps);		// Get sitewise parsimony score; if(PenaliseGaps) { all gap characters add 1 to score)
	// Masking functions
	int CreateZeroMask();					// Masks out all sites with no changes
	int CreatePropMask(double Prop);		// Masks out all sites except the (Prop) slowest evolving 
	int CreateCountMask(int Count);		// Masks out all sites except the Count slowest evolving
	int CreateLimitMask(int MaxChange);	// Mask out all sites except those with fewer than MaxChange changes
	void GetMask(vector <bool> &Mask) { Mask = m_vbMask; };		// Gets the mask for parsimony
	void SetMask(vector <bool> &Mask) { m_vbMask = Mask; };		// Sets the mask for parsimony
private:
	// Variables
	int m_iNoSeq;								// Number of sequences
	int m_iNoNode;								// Number of internal nodes (min 1);
	int m_iSize;								// Length of sequences
	int m_iChar;								// Characters
	CTree *m_pTree;
	CData *m_pData;
	bool ***m_BrNodes;						// Data held at the branch nodes [m_iNoSeq][m_iSize][m_iChar]
	vector <vector <int> > m_vviLfNodes;	// Data held at the leaf nodes
	vector <bool> m_vbMask;						// Whether to do calculation at a site
	// Implementation
	void CreateMemory(int NoSeq, int Size, EDataType Type, vector <vector <int> > *Data);
	void CleanMemory();						// Function to clean the memory	
	// The Step1 functions
	int Step1(int NT, int NF);							// Do the first round of parsimony stuff
	int Step1(int NT, int NF, vector <int> *SW);		// Do the first round of parsimony stuff and remember the score for each site
	inline int Form2Leaf(int C1, int C2, bool *bPar);
	inline int FormLeafBranch(int C1, bool *bC2, bool *bPar);
	inline int Form2Branch(bool *bC1, bool *bC2, bool *bPar);
};

