
#ifndef __TREE_HEADER
#define __TREE_HEADER

#include "Leaphy.h"

#define TREE_START 1			// Number that input trees first species start with
#define SET_BRANCH 0.1			// Default branch length

enum ENodeType	{ branch, leaf, root};

struct SSplit{ int BrLabel; vector <int> Left, Right; };

class CNode		{
public:
    // Member variables
    ////////////////////////
	int m_iNoLinks;
	int m_iInternalNodeNum;
    vector <int> m_viLink;
    vector <int> m_viBranch;
    ENodeType m_NodeType;
    // Member functions
    //////////////////////
	// Note parent node ALWAYS should be specified in the last branchx and linkx where applicable
	// General constructor
	CNode(int NoLinks = -1,int *LinkList = NULL);
    // Constructor for internal nodes
    CNode(int linka, int linkb, int linkc, int brancha, int branchb, int branchc, int IntVal = -1);
    // Constructor for external nodes
    CNode(int linka, int brancha, int IntVal = -1);
	// Other constructor
	CNode(vector <int> Links, vector <int> Branches, int IntVal);
	// Copy constructor
	CNode(const CNode &Node);
    // Destructor function
    virtual ~CNode();
	// Member functions specific for bifurcating trees where linkc is defined as a parent
	void SetNode(int linka, int linkb, int linkc, int brancha, int branchb, int branchc, int IntVal);
	void SetNode(int la, int lb, int ba, int bb, int IntVal);
    void SetNode(int linka, int brancha, int IntVal);
	// General functions for assigning links (No parent assumed)
	void SetNode(int NoLinks, int *LinkList);
	// Cleaning operation
	void CleanNode();
	// Operator functions
    CNode &operator=(const CNode &);
};

ostream& operator<<(ostream& os, const CNode &Node);

class CTree
{
public:
    // Member functions
    /////////////////////////////////

    // Constructor functions
    CTree(string TREE, int NoSeq, bool AllowFail = false, CData *Data = NULL, bool AllowSubTree = false);			// Basic constructor
	CTree(string TREE, bool GetTreeFromFile,CData *Data =NULL, bool AllowFail = false, bool AllowSubTree = false);	// Risky constructor (may be wrong tree for data)
	void CreateTree(string TREE,int NoSeq,bool CheckVar = true, bool AllowFail = false,bool AllowSubTree = false,CData *D = NULL);// Underlying construction function
    // Blank constructor
    CTree(CData *D = NULL);
	// Copy Constructor
	CTree(const CTree &Tree);
    // Destructor function
    virtual ~CTree();
	// Memory functions
	void CleanTree();
    void GetMemory(int NoSeq);

	// Access functions
	////////////////////////////////////////
	// Main tree
    inline bool IsRooted()  { return m_bRooted; }		// Returns whether the tree is rooted
    inline int Root()		{ return m_iRootNode; }			// Returns the root node index
	inline int NoNode()		{ return m_iNoNode; }		// Number of nodes in tree
	inline int NoBra()		{ return m_iNoBra; }		// Returns number of branches in tree
	inline int NoOptBra()	{ return m_iNoOptBra; }		// Returns the number of optimised branches in tree
	inline int NoSeq()		{ return m_iNoSeq; }		// Returns number of sequences in tree
	inline int StartCalc()	{ return m_iStartCalc; }	// Returns where the calculation starts
	bool OldTree()			{ return m_bOldTree; }		// Returns whether the tree is old or not
	bool SetOldTree(bool V)	{ m_bOldTree = V; return V; }	// Sets the m_bOldTree value
	int BestStartCalc();								// Sets the m_iStartCalc to the best value
	int SetStartCalc(int S) { m_iStartCalc=S;return S;}	// Forces m_iStartCalc to a value -- power user function because can cause problems
	bool FastCalcOK()		{ return m_bFastCalcOK; }	// Returns whether the fast calc is okay or needs resetting
	void SetFastCalcOK(bool V)	{ m_bFastCalcOK = V; }	// Sets the m_bFastCalcOK flag
	int SetOutBra(int V)	{ m_iOutBra = V; }			// Sets m_iOutBra
	void OutBra()			{ m_iOutBra = 1; }			// Sets m_iOutBra to output branches
	void NoOutBra()			{ m_iOutBra = 0; }			// Sets m_iOutBra to not output branches
	void OutBraNum()		{ m_iOutBra = 2; }			// Sets m_iOutBra to output branch numbers
	void OutName()			{ m_bOutName = true; }		// Set for outputting names
	void NoOutName()		{ m_bOutName = false; }		// Set for not outputting names
	void OutLabel()			{ m_bOutLabel = true; }		// Set for outputting labels on tree
	void NoOutLabel()		{ m_bOutLabel = false; }	// Set for not outputting labels on tree
	void CreateBranchLabels(bool Force = true);			// Create labels on branches that correspond with their numbers
	bool Valid() { return m_bValid; }					// Returns whether a valid constructor has been run
	void ForceReady() { m_bReady = true; }				// Makes the tree ready.
	vector <string> Names() { return m_vsName; }		// Returns the vector of names
	void SetNames(vector <string > NewNames, bool Overwrite = false);	// Set the names in the tree (e.g. for output)
	bool BranchLabels() { if(m_viBranchLabels.empty()) { return false; } return true; }
	vector <int> Labels() { return m_viBranchLabels; }
	int NoLabels() { return m_iNoBraLabels; }
	// Branch
	vector <double> Branches();					// Returns a vector of the branch lengths
	double B(int Branch);						// Returns the value of branch length
	double *OptimiserB(int Branch);				// Returns the pointer to the value to be used in optimisation
	double TreeLength();						// Returns the tree length;
	CPar *pBra(int B)		{ return m_vpBra[B]; }		// Returns the parameter for a branch
	bool OptB(int Branch)	{ return m_vpBra[Branch]->Opt(); }								// Returns whether the branch is optimised or not (also u
	double BScale(int B)	{ return m_vpBra[B]->Scaler(); }							// Returns the scaling factor for a branch=
	void BRescale(int B)	{ m_vpBra[B]->Rescale(); }									// Rescales the value in a branch
	bool SetOptB(int B, bool V) { return m_vpBra[B]->SetOptimise(V); }					// Sets whether a branch is optimised or not
	double SetB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Set branch length Branch to Value
	double MulB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Multiply branch length Branch by Value
	double AddB(int Branch, double Value, bool Update = false, bool Rescale = false);		// Add Value to branch length Branch
	double QuadB(int Branch);																// Get the average value of Branch and the (upto) 4 surrounding branches
	void UpdateB(int Branch) { m_vpBra[Branch]->UpdatePar(); }								// Update the branch parameter
	inline int BraLink(int B, int L) { return m_ariBraLinks[B][L]; }						// Get the L link for branch B
	inline void ReplaceBraLink(int B,int L,int Val) { m_ariBraLinks[B][L]=Val;};			// Set the L link for branch B to Val; returns Val;

	// Node
	inline int NoLinks(int N)	{ return (int) m_Node[N]->m_viLink.size(); }	// Returns the # links in a node
	int NodeLink(int N,int L)	{ return m_Node[N]->m_viLink[L]; }		// Returns link L of node N
	bool IsNodeLink(int N, int Val) { return IsIn(Val,m_Node[N]->m_viLink); }	// Returns whether a Val is in Node N m_viLink
	void ReplaceNodeLink(int N, vector<int>L);							// Replace m_viLink of node N with L;
	void ReplaceNodeLinkElement(int N, int Element, int Val);			// Replace a single element of m_viLink in Node N with Val
	int NodeBra(int N, int B)	{ return m_Node[N]->m_viBranch[B]; }	// Returns Branchlink B of Node N
	void ReplaceNodeBra(int N, vector <int>B);							// Replace m_viBranch of Node N with B
	void ReplaceNodeBraElement(int N, int Element, int Val);			// Replace a single element of m_viBranch in Node N with Val
	ENodeType NodeType(int N)	{ return m_Node[N]->m_NodeType; }		// Returns the NodeType of Node N
	bool NodeNull(int N) { if(m_Node[N] == NULL) { return true; } return false; }	// Returns whether a node is NULL or not
	int NoLeafLink(int N);												// Returns the number of leaf links for Node N
	void AssignNodeType(int N, ENodeType Type);							// Set the NodeType

	// Output functions
	friend ostream& operator<<(ostream& os, CTree &Tree);		// Standard output routine
	bool OutDetail(ostream &os = cout, bool ForceExit = false);	// Detailed output routine (nodes, branches, and tree)

	// Tree modification routines
	void Unroot();											// If rooted, this function permanently unroots it
	void OrderNode(int NodeNum = -1,bool DoBraToo = true);	// Orders tree nodes to ascending numerical value (-1 does all)
	int CutBranch(int Branch,int Link);					// Cuts a branch: returns the branch on priority subtree where they used to be attached -- Subtrees are maintained; Link specifies which one is given priority
	int RemoveBraNode(int Node);						// Removes a node from a tree: returns the branch where link used to be attached -- NB: One Link must be set to -1
	int RemoveLeafNode(int RemNode);
	int AddSeq(int SeqNum, int Branch, double BranchProp =0.5);	// Adds a sequence to the tree
							// SeqNum = sequence added; Branch = branch added to; BranchProp = relative position in branch (from low to high node nums)
	int GetStart(bool replace = true);			// Find branch in tree from which to recurse
	double GetTreeLength(bool first = true, int NTo = -1, int NFr = -1);	// Get remaining stuff in a tree
	void BuildOriSubTree(CTree *T, vector <bool> NodesBool);										// Returns subtree from an array of bools describing which nodes it covers
	void BuildOriSubTree(CTree *T, vector <int> LeafMap, vector <int> NCover, vector <int> NFrom);	// Returns subtree from LeafMap and NodesCovered
	void BuildOriSubTree(CTree *T, vector <int> NodesCovered);										// Returns subtree from the nodes it covers
	void PerformSPR(int OriB, int OriL, int IBr,CTree Triplet,vector <int> LeafList);
	void BuildBraLinks(bool Verify = true);					// Builds the branch links from nodes

	// Function to create a consistent output for comparison
	vector <int> ConstOut();
	int NodeDist(int Node1, int Node2, int NodeFrom = -1);		// Count number of branches separating nodes i and j;
	int BranchDist(int Br1, int Br2, bool AllowZero = false);	// Count the number of branches separating branches Br1 and Br2; !AllowZero means that very short branches will not be considered real branches
	// Functions for comparing trees using Robinson-Foulds distance
	int GetRFDist(CTree &Tree);			// Standard comparison between two trees of same number of taxa
	bool IsCompatible(CTree &SubTree);	// Compares tree <SubTree> with #Seq <= *this->#Seq to check whether they are compatible (High level function accepting tree objects)

	// Functions for adding sequences to a subtree tree based on an existing full tree (Greedy algorithm for maximising tree length)


	// Functions for testing properties of trees
	bool IsNode(int Node);  // Used for assert statements
	bool IsBra(int Branch); // Used for assert statements
	bool IsCutTree();		// Whether tree has had sequences removed...
	bool GoodBra(int Branch);	// Is an active branch in the tree
	bool GoodNode(int Node);	// Is an active node in the tree

	// Tree split-based functions
	void BuildSplits();									// Build the set of splits associated with a tree. Current implementation always forces the rebuild
	SSplit GetSplit(int Bra);						// Return the split set for branch Bra
	int BranchSets(int BranchNum, vector <int> *Left, vector <int> *Right);	// Find the sets of leaf sequences that the branch splits
	int FindBra(int Node_i, int Node_j);	// Find branch linking Nodes i and j, returns -1 if none
		// Returns the value of total number in the Left set (i.e. Left = m_ariBraLinks[0] )
	void GetBraSets(int NTo, int NFr, vector <int> *List, bool First = true);
	void OutSplits(ostream &os = cout);
	// Functions to get pairwise distances from a tree
	vector <double> GetTreePW();
	vector <double> GetAllTreePW();
	void PWDistSub(int NodeTo, int NodeFrom,vector <double> *d,bool DoInternalNodes = false);
	vector <double> GetPartialTreeDist(vector <int> LeafMap, vector <vector <int> > NBelow);
	vector <double> GetSubTreePW(vector <int> LeafMap, vector <vector <int> > NBelow, vector <double> BaseDist);

	// Function to get RMSD between tree pairwise distances and observed pairwise distances
	vector <double> MeasureSeqStability(vector <double> PWDist);	// Gets RMSD by sequence. If SubSet != NULL then only the comparisons between the sequences specified are considered
	vector <double> MeasureBraStability(vector <double> PWDist);	// Gets RMSD by branch.
	vector <double> MeasureNodeStability(vector <double> PWDist);	// Gets RMSD by three branches around a node.

	double RMSD(vector <double> PWDist);		// Get the overall RMSD score for a tree
	double SubTreeRMSD(vector<int> LeafMap, vector <vector <int> > NBelow, vector <double> PWDist, vector <double> SubPWDist);	// Get the RMSD score for a subtree.

	// Function to get all of the nodes of a certain depth from a particular node; bool GetLess == true, will also return leaf nodes when they fall within this range
	vector <int> GetNodesOfDepth(int InitNode, int NodeDepth, bool GetLess, vector <int> *NodesFrom = NULL, vector <int> *NodeCov = NULL, vector <double> *ExtBra = NULL ,int NodeFr = -1, bool First = true);
	// Centering point functions used for snap
	vector <int> BranchCP(int CP, int depth, vector <int> *NodeFr, vector <int> *NodesCovered = NULL, vector <double> *ExtBra = NULL);
	vector <int> NodeCP(int Node, int depth, vector <int> *NodeFr, vector <int> *NodesCovered = NULL, vector <double> *ExtBra = NULL);
	void ReplaceTreeCP(CTree *NT,vector <int> LeafMap,vector <int> NCover,bool VerifyBranchLinks = true);		// Overwrites the current tree (from CPs) and puts NewTree in its place
    // Operator= function
    void operator=(const CTree &);

	vector <vector <int> > GetKnotClusters(vector <bool> IntNodeChanges, int ChangeRad = 2);		// Function that takes an array[NoNodes()] of which branches change (or zero lengths, or whatever)
					// and returns a set of clusters of overlapping changes of SNAP radius = ChangeRad
	///////////////////////////////////////////////////////
	// Rearrangement functions
	void GetSPRList(vector <CTree *> *TreeList, int MinDist = 1, int MaxDist = BIG_NUMBER);
private:
    // Private member variable
    ///////////////////////
	int m_iRootNode;			// The root node (-1 if not rooted)
	bool m_bRooted;				// Whether the tree is rooted
    int m_iNoNode;				// The # of nodes in tree
    int m_iMemNode;				// Number of nodes stored in memory
    int m_iNoBra;				// The # of branches in tree
	int m_iNoOptBra;			// The # of branches optimised in the current tree
    int m_iNoSeq;				// The # of sequences in tree
	int m_iStartCalc;			// Leaf node from which to start calculations
	vector <CPar *> m_vpBra;	// The Branch length parameters
    CNode ** m_Node;			// The nodes in the tree
    bool m_bReady;				// Whether tree is ready
    int m_iOutBra;				// Out branch [0=no branches,1=branches,2=branch numbers]
	int m_bOutName;				// Out Names
	bool m_bOutLabel;			// Whether to output tree labels or not
	vector <string> m_vsName;	// The names in the sequence
    int **m_ariBraLinks;		// the nodes linked to each branch [branch_number][links]
	bool m_bOldTree;			// Whether the tree has already been through a round of SNAP
	bool m_bFastCalcOK;			// Whether the tree is okay for FastCalc computations
	bool m_bValid;				// Flag to identify whether a valid constructor has been run
	vector <int> m_viBranchLabels;	// Labels for different branch types (e.g. parameter per branch models)
	int m_iNoBraLabels;				// Number of unique branch labels
	vector <SSplit> m_vSplits;	// Vector containing the splits on a tree

	// Private member functions
    ////////////////////////////////////////////////////////////

	// Find closest: gets closest nodes in a string
    void find_closest(string *tree, int *c1, int *c2, int *p,int n_p, double *bra,int *label, int *IntVal);
	void GetBraDetails(string s, int *node, double *bra, int *label);
	ostream& OutNode(int FromNodeNum, int ToNode, ostream &os);
	ostream& OutBranch(int ToNode, int FromNode, ostream &os);
	double DoBranch(string *tree,int *pos, int *IntVal);
	void BuildBranches(int NT = -1, int NF = -1);	// Builds the branchs of a tree from nodes recursively
	bool ValidateTree(bool AllowExit = true);

	////////////////////////////////////////////
	// Rearrangement functions
	void BuildBraSA(CTree *T, int Seq, vector <CTree *> *TList, int OriBr, int SPRID, int MinDist, int MaxDist);
	void Branch_SA(CTree *T, int Seq, vector <CTree*> *TList, int First, int NTo, int NFr, int OriBr, int SPRID, int MinDist, int MaxDist);
	void DoSA(CTree *T, int Seq, vector <CTree*> *TList, int First, int NTo, int NFr, int Br, bool IsExtBra, int OriBr, int SPRID, int MinDist,int MaxDist);
};

void ExpandNode(int Node, string *String, int Stringpos, CTree *TREE);
int IsTreeGap(char Check);
// Functions for finding the greedy subtree
CTree FindGreedySubTree(CTree *FullTree, int NoSeq); 		// Driver function: uses the greedy algorithm to identify the optimal subtree
double TravAddGreedy(CTree *CurT, int To, int Fr, int Seq2Add, vector <double> *PWDist, int *BestBra); // Function to traverse a tree and test which position maximises distance
double GetDist(int Add, int a, int b,vector <double> *PWDist);	// Check how much improvement a given sequence added to the tree would provide
void GreedySeq2Tree(int Bra,int Seq2Add, CTree *CurTree, vector <double> *PWDists);  // Adds the sequence to the tree

bool SplitsCompatible(vector <SSplit> Splits1, int S1_seq, vector <SSplit> Splits2, int S2_seq);	// Low level function just comparing a set of splits
bool CompareSplit(SSplit S1, SSplit S2);															// Simple function for comparing splits

#endif
