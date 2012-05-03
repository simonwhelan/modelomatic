////////////////////////////////////////////////////////////////
// The implementation of the Parsimony class
////////////////////////////////////////////////////////////////

#include "parsimony.h"

////////////////////////////////////////////////////////////////
// 1. Constructor for known data
CParsimony::CParsimony(CData *D, CTree *T)	{
	// Some initialisation
	m_pData = D; m_pTree = T;
	// Some verification
	if(D->m_iNoSeq != T->NoSeq()) { Error("\nNumber of sequences in tree and data don't match in CParsimony::CParsimony constructor...\n\n"); }
	// Get memory
	CreateMemory(D->m_iNoSeq,D->m_iSize,D->m_DataType,NULL);
	// Create simple mask
	CreateZeroMask();
}
// 2. Constructor for uncharacterised data
CParsimony::CParsimony(EDataType Type,vector <vector <int> > *Data,CTree *T)	{
	// Some initialisation
	m_pData = NULL; m_pTree = T; m_iNoSeq = T->NoSeq();
	// Some verification
	if((int)Data->size() != T->NoSeq()) { Error("\nNumber of sequences in tree and data don't match in CParsimony::CParsimony constructor...\n\n"); }
	// Get memory
	CreateMemory(T->NoSeq(),(int)Data[0].size(),Type,Data);
	// Create simple mask
	CreateZeroMask();
}

//////////////////////////////////////////////////////////////
// Destructor function
CParsimony::~CParsimony()	{
	CleanMemory();
}

void CParsimony::CleanMemory()	{
	int i,j;
	m_pTree = NULL;
	m_pData = NULL;
	m_vviLfNodes.clear();
	if(m_BrNodes != NULL) { FOR(i,m_iNoNode) { FOR(j,m_iSize) { DEL_MEM(m_BrNodes[i][j]); } DEL_MEM(m_BrNodes[i]); } DEL_MEM(m_BrNodes); }
}	

//////////////////////////////////////////////////////////////
// Function to create memory
// ---
//  1. The number of sequences specifies the size of m_vviBrNodes;
//  2. The length of the sequences
//  3. The number of characters in the data
//  4. The EDataType; used to define m_iChar
//  5. If Data != NULL then m_vviLfNodes is assigned from the vector;
void CParsimony::CreateMemory(int NoSeq, int Size, EDataType Type, vector <vector <int> > *Data)	{
	int i,j;
	// Some basic numbers
	m_iNoSeq = NoSeq; 
	m_iSize = Size;
	m_iChar = NumStates(Type);
	m_iNoNode = max(1,m_iNoSeq - 2); 
	if(m_iChar < 2) { Error("\nInitialising parsimony with too few states...\n"); }
	// Get the leaf nodes
	if(Data != NULL) { m_vviLfNodes = *Data; }
	// Get the internal nodes
	// 1. Old version using bools
	GET_MEM(m_BrNodes,bool**,m_iNoNode);
	FOR(i,m_iNoNode) { GET_MEM(m_BrNodes[i],bool *,m_iSize); FOR(j,m_iSize) { GET_MEM(m_BrNodes[i][j],bool,m_iChar); } }
}

////////////////////////////////////////////////////////
// Simple function to create a mask for calculations
int CParsimony::CreateZeroMask()	{
	int i,j,Seq,Sites = 0;
	m_vbMask.clear(); m_vbMask.assign(m_iSize,false); 
	if(m_pData != NULL) {
		FOR(i,m_iSize) { 
			Seq = -1;
			FOR(j,m_iNoSeq) { 
				if(m_pData->m_ariSeq[j][i] == m_iChar) { continue; }
				if(Seq == -1) { Seq = m_pData->m_ariSeq[j][i]; continue; }
				if(m_pData->m_ariSeq[j][i] != Seq) { Sites += m_pData->m_ariPatMap[i]; m_vbMask[i] = true; break; }
			}
		}
	} else {
		Error("\nHaven't done simple mask for fake data yet...\n"); 
	}
	return Sites;
}

/////////////////////////////////////////////////////////////////
// Function to mask a certain proportion of fast evolving sites
int CParsimony::CreatePropMask(double Prop)	{ 
	if(m_pData == NULL) { Error("\nHaven't done CreatePropMask for fake data yet...\n");  }
	return CreateCountMask((int) (Prop * m_pData->m_iTrueSize)); 
}
/////////////////////////////////////////////////////////////////
// Function to mask all but the slowest "Count" sites
int CParsimony::CreateCountMask(int Count)	{
	int i,j,Sites = 0;
	vector <int> Scores = SitewiseScore(true),Sorted;
	if(m_pData == NULL) { Error("\nHaven't done CreatePropMask for fake data yet...\n");  }
	// Initialise
	m_vbMask.clear(); m_vbMask.assign(m_iSize,false); 
	// Get the sorted list
	Sorted = Scores; Sort(&Sorted);
	if((int)Sorted.size() != m_iSize) { Error("\nMismatch between number of sites and stuff in CParsimony::CreateCountMask()...\n"); }
	FOR(i,m_iSize) { 
		if(Sorted[i] == 0)  { continue; }
		if(i > 0) { if(Sorted[i] == Sorted[i-1]) { continue; } }
		FOR(j,m_iSize) { if(Scores[j] == Sorted[i]) { Sites += m_pData->m_ariPatOcc[j]; m_vbMask[j] = true; } }
		if(Sites > Count) { break; }
	}
	return Sites;
}
//////////////////////////////////////////////////////////////////
// Function to mask to limit of MaxChange changes
int CParsimony::CreateLimitMask(int MaxChange)	{
	int i,Sites = 0;
	vector <int> Scores = SitewiseScore(true);
	if(m_pData == NULL) { Error("\nHaven't done CreatePropMask for fake data yet...\n");  }
	// Initialise
	m_vbMask.clear(); m_vbMask.assign(m_iSize,false); 
	// Get the limited sites
	FOR(i,(int)Scores.size()) { if(Scores[i] <= MaxChange && Scores[i] > 0) { Sites += m_pData->m_ariPatOcc[i]; m_vbMask[i] = true; } }
	return  Sites;
}
/*	**************************************************
		Functions to do the calculations
	************************************************** */

// Main function for getting the score
int CParsimony::Score(bool FullScore)	{
	int i,j,k,Start = m_pTree->GetStart();
	vector <bool> OldMask;
	// Check everything is okay...
	assert(m_BrNodes != NULL);
	// Check a mask exists
	if(m_vbMask.empty()) { m_vbMask.assign(m_iSize,true); }
	// Clean the memory
	FOR(i,m_iNoNode) { FOR(j,m_iSize) { FOR(k,m_iChar) { m_BrNodes[i][j][k] = false; } } }
	// Deal with the mask is FullScore
	if(FullScore) { OldMask = m_vbMask; CreateZeroMask(); }
	// Do the calculations
	if(Start >= m_pTree->NoSeq()) { Error("\nTrying to do parsimony starting from a branch node...\n"); }
	// Deal with the mask
	if(FullScore) { m_vbMask = OldMask; }
	return Step1(Start,-1);
}

// Main function for getting the score
vector <int> CParsimony::SitewiseScore(bool PenaliseGaps)	{
	int i,j,k,Start = m_pTree->GetStart();
	vector <int> SWscore;
	// Check everything is okay...
	assert(m_BrNodes != NULL);
	// Check a mask exists
	if(m_vbMask.empty()) { m_vbMask.assign(m_iSize,true); }
	// Clean the memory
	FOR(i,m_iNoNode) { FOR(j,m_iSize) { FOR(k,m_iChar) { m_BrNodes[i][j][k] = false; } } }
	// Do the calculations
	if(Start >= m_pTree->NoSeq()) { Error("\nTrying to do parsimony starting from a branch node...\n"); }
	Step1(Start,-1,&SWscore);
	// Deal with gaps if required
	if(PenaliseGaps) {
		if(m_pData == NULL) { Error("\nDunno how to deal with fake data in CParsimony::SitewiseScore(true)\n"); }
		FOR(i,m_iSize) { 
			FOR(j,m_iNoSeq) { 
				if(m_pData->m_ariSeq[j][i] == m_iChar) { SWscore[i] += m_pData->m_ariPatOcc[i]; }
	}	}	}
	return SWscore;
}

/* *************************************************
		Tree traversals to get the scores
   ************************************************* */


// Post order tree traversal to get the score
int CParsimony::Step1(int NT, int NF)	{
	int i,j=0,Score = 0;
	int Links[2];
	// Do the tree traversal...
	FOR(i,m_pTree->NoLinks(NT)) { 
		if(m_pTree->NodeLink(NT,i) == NF)		{ continue; }
		Links[j++] = m_pTree->NodeLink(NT,i); 
		if(m_pTree->NodeLink(NT,i) < m_iNoSeq)	{ continue; }
		Score += Step1(m_pTree->NodeLink(NT,i),NT);	
	}
	// Do calculations for final node 
	if(NF == -1) { 
		if(m_pTree->NodeType(NT) == branch) { Error("\nError: started parsimony from an internal node..."); }
		NF = m_pTree->NodeLink(NT,0);
		if(m_pData != NULL) {	// Do normal data
			if(NF < m_iNoSeq)	{	// Do for a pair of sequences
				assert(m_iNoSeq == 2);
				FOR(i,m_iSize) { if(m_vbMask[i]) { Score += m_pData->m_ariPatOcc[i] * Form2Leaf(0,1,m_BrNodes[0][i]); } }
			}
			else { // Do for a normal set of sequences
				NF -= m_iNoSeq;
				FOR(i,m_iSize)	{
					if(!m_vbMask[i]) { continue; }
					// Do intersection
					if(m_pData->m_ariSeq[NT][i] != m_iChar) { 
						if(m_BrNodes[NF][i][m_pData->m_ariSeq[NT][i]]) { FOR(j,m_iChar) { m_BrNodes[NF][i][j] = false; } } 
						else { Score += m_pData->m_ariPatOcc[i]; }
						m_BrNodes[NF][i][m_pData->m_ariSeq[NT][i]] = true;
					}
		}	}	}	else {				// Do faked data
			Error("\nParsimony for faked data not ready yet...\n"); 
		}
		return Score; 
	}
	if(j != 2) { Error("\nUnexpected Multifurcation...\n"); }
	// Do the parsimony step for real data 
	NT -= m_iNoSeq;				// Adjust node to for the correct space
	if(Links[0] > Links[1]) { i = Links[1]; Links[1] = Links[0]; Links[0] = i; }
	if(m_pData != NULL)	{
		// Do step 1 parsimony
		if(Links[0] < m_iNoSeq)	{
			if(Links[1] < m_iNoSeq)	{	// Case both children are leaf nodes
				FOR(i,m_iSize) { if(m_vbMask[i]) { Score += m_pData->m_ariPatOcc[i] * Form2Leaf(m_pData->m_ariSeq[Links[0]][i],m_pData->m_ariSeq[Links[1]][i],m_BrNodes[NT][i]); } }
			} else {					// Case Link[0] is a leaf node; Link[1] is a branch node
				FOR(i,m_iSize) { if(m_vbMask[i]) { Score += m_pData->m_ariPatOcc[i] * FormLeafBranch(m_pData->m_ariSeq[Links[0]][i],m_BrNodes[Links[1] - m_iNoSeq][i],m_BrNodes[NT][i]); } }
		}	} else {	
			FOR(i,m_iSize) { if(m_vbMask[i]) { Score += m_pData->m_ariPatOcc[i] * Form2Branch(m_BrNodes[Links[0] - m_iNoSeq][i],m_BrNodes[Links[1] - m_iNoSeq][i],m_BrNodes[NT][i]); } }
		} 
	} else {
	// Do the parsimony step for fake data...
		Error("\nParsimony for faked data not ready yet...\n"); 
	}
	return Score;
}

// Post order tree traversal to get the score per site
int CParsimony::Step1(int NT, int NF,vector <int> *SWscore)	{
	int i,j=0,Score = 0;
	int Links[2];
	if(SWscore == NULL) { Error("\nNeed vector to do CParsimony::Step1 sitewise scores..."); }
	else if(SWscore->empty()){ if(m_pData != NULL) { SWscore->assign(m_iSize,0); } }
	// Do the tree traversal...
	FOR(i,m_pTree->NoLinks(NT)) { 
		if(m_pTree->NodeLink(NT,i) == NF)		{ continue; }
		Links[j++] = m_pTree->NodeLink(NT,i); 
		if(m_pTree->NodeLink(NT,i) < m_iNoSeq)	{ continue; }
		Score += Step1(m_pTree->NodeLink(NT,i),NT,SWscore);	
	}
	// Do calculations for final node 
	if(NF == -1) { 
		if(m_pTree->NodeType(NT) == branch) { Error("\nError: started parsimony from an internal node..."); }
		NF = m_pTree->NodeLink(NT,0);
		if(m_pData != NULL) {	// Do normal data
			if(NF < m_iNoSeq)	{	// Do for a pair of sequences
				assert(m_iNoSeq == 2);
				FOR(i,m_iSize) { SWscore->at(i) += m_pData->m_ariPatOcc[i] * Form2Leaf(0,1,m_BrNodes[0][i]); }
			}
			else { // Do for a normal set of sequences
				NF -= m_iNoSeq;
				FOR(i,m_iSize)	{
					// Do intersection
					if(m_pData->m_ariSeq[NT][i] != m_iChar) { 
						if(m_BrNodes[NF][i][m_pData->m_ariSeq[NT][i]]) { FOR(j,m_iChar) { m_BrNodes[NF][i][j] = false; } } 
						else { SWscore->at(i) += m_pData->m_ariPatOcc[i]; }
						m_BrNodes[NF][i][m_pData->m_ariSeq[NT][i]] = true;
					}
		}	}	}	else {				// Do faked data
			Error("\nParsimony for faked data not ready yet...\n"); 
		}
		return Score; 
	}
	if(j != 2) { Error("\nUnexpected Multifurcation...\n"); }
	// Do the parsimony step for real data 
	NT -= m_iNoSeq;				// Adjust node to for the correct space
	if(Links[0] > Links[1]) { i = Links[1]; Links[1] = Links[0]; Links[0] = i; }
	if(m_pData != NULL)	{
		// Do step 1 parsimony
		if(Links[0] < m_iNoSeq)	{
			if(Links[1] < m_iNoSeq)	{	// Case both children are leaf nodes
				FOR(i,m_iSize) { SWscore->at(i) += m_pData->m_ariPatOcc[i] * Form2Leaf(m_pData->m_ariSeq[Links[0]][i],m_pData->m_ariSeq[Links[1]][i],m_BrNodes[NT][i]); }
			} else {					// Case Link[0] is a leaf node; Link[1] is a branch node
				FOR(i,m_iSize) { SWscore->at(i) += m_pData->m_ariPatOcc[i] * FormLeafBranch(m_pData->m_ariSeq[Links[0]][i],m_BrNodes[Links[1] - m_iNoSeq][i],m_BrNodes[NT][i]); }
		}	} else {	
			FOR(i,m_iSize) { SWscore->at(i) += m_pData->m_ariPatOcc[i] * Form2Branch(m_BrNodes[Links[0] - m_iNoSeq][i],m_BrNodes[Links[1] - m_iNoSeq][i],m_BrNodes[NT][i]); }
		} 
	} else {
	// Do the parsimony step for fake data...
		Error("\nParsimony for faked data not ready yet...\n"); 
	}
	return Score;
}
/* ************************************************************
		The union and intersection functions for parsimony 
   ************************************************************ */

// Both children are leaf nodes...
// Form bPar from characters C1 and C2
int CParsimony::Form2Leaf(int C1, int C2, bool *bPar)	{
	int i;
	// Deal with gaps
	if(C1 == m_iChar) { 
		if(C2 == m_iChar) { FOR(i,m_iChar) { bPar[i] = true; } return 0; }
		bPar[C2] = true; return 0;
	} else if(C2 == m_iChar) { bPar[C1] = true; return 0; }
	// Do the normal stuff
	if(C1 == C2) { bPar[C1] = true; return 0; }
	bPar[C1] = true;
	bPar[C2] = true;
	return 1;
}

// One child is a leaf node; the other is a branch node
// Form bPar from C1 leaf node and bC2 a branch node
int CParsimony::FormLeafBranch(int C1, bool *bC2, bool *bPar) {
	int i;
	if(C1 == m_iChar) { FOR(i,m_iChar) { bPar[i] = bC2[i]; } return 0; }
	if(bC2[C1]) { bPar[C1] = true; return 0; }
	FOR(i,m_iChar) { bPar[i] = bC2[i]; }
	bPar[C1] = true; 
	return 1;
}

// Does step 1 of parsimony on children bC1 and bC2 to form bPar
// Returns false if no intersection found
int CParsimony::Form2Branch(bool *bC1, bool *bC2, bool *bPar)	{
	int i;
	bool DoUnion = true;
	FOR(i,m_iChar) { if(bC1[i] && bC2[i]) { bPar[i] = true; DoUnion = false; } }
	if(!DoUnion) { return 0; }
	FOR(i,m_iChar) { if(bC1[i] || bC2[i]) { bPar[i] = true; } }
	return 1;
}


