//
// Implementation of importance sampling on tree

#include "optimise.h"
#include "TreeList.h"
#define DEBUG_ADDSEQ 0

#if IS_LNL_CALC
extern vector <double> PWDists;
extern int TABU_RADIUS;
extern int EXIT_OBS,OptObs;
extern double PROB_RAN_SEQ_REM;
#else
vector <double> PWDists;
int TABU_RADIUS;
int EXIT_OBS,OptObs;
double PROB_RAN_SEQ_REM;
#define DO_NOT_RUN 1
#endif
const int NumResamp = 1;

///////////////////////////////////////////////////////////////////////////////////////////////
// This is the main resampling function
double GetTree(CBaseModel *M, CTreeList *Four, CTreeList *Five, CTreeList *Six, LazyType DoLazy,bool DoOutput)	{
	int i;
	double BestLnL = -BIG_NUMBER, lnL;
	CTree BestTree;
	BestTree = *M->m_pTree;
#if DO_NOT_RUN == 1
	Error("\nSet IS_LNL_CALC to 1 in Leaphy.h file when compiling a program that does likelihood computation");
#endif
	FOR(i,NumResamp)	{
		// Reset to best current tree
		if(i != 0) { M->MakeNewTree(&BestTree,true); }
		// Resample tree space
		// ------------------------
		GetRMSDSampTree(M,Four,Five,Six,LAZY_SA,PROB_DO_MP_SA,DoOutput);
		// Store the new tree
		if(!DoLazy) { lnL = FullOpt(M); } else { lnL = LazyOpt(M); }
		if(lnL > BestLnL)  {
			BestTree = *M->m_pTree;
			BestLnL = lnL;
	}	}
	M->MakeNewTree(&BestTree,true);
	return BestLnL;
}

/////////////////////////////////////////////////////////////////////////////////////
// This is the routine that gets the starting stepwise addition tree
// --------------------------------------------------------------------------------
// 1. Get the stepwise addition tree using Maximum likelihood or Maximum parsimony
void GetSAStartTree(CTree *RetTree, CBaseModel *M, bool DoOutput, LazyType DoLazy, int DoParsimony)	{
	int i,BraNum = M->NoSeq(),InNode = M->NoSeq();
	double CurlnL;
	bool OriMLCalc = M->IsLikelihoodCalc();
	stringstream sT;
	vector <int> Seq2Add, Lnk;
	CTree T(M->m_pData);
	assert(RetTree != NULL);
	// Do some output
	if(DoOutput)	{ cout << "\n---------------- Building stepwise addition starting tree ----------------\n"; }
	// Initialise Seq2Add vector
	FOR(i,M->NoSeq()) { Seq2Add.push_back(i); } random_shuffle(Seq2Add.begin(),Seq2Add.end());
	// For initial tree from first three sequences
	sT << "(" << Seq2Add[0]+1 << ":0.1," << Seq2Add[1]+1 << ":0.1," << Seq2Add[2]+1 << ":0.1);";
	T.CreateTree(sT.str(),M->NoSeq(),false);
	FOR(i,3) { Seq2Add.erase(Seq2Add.begin()); }
	// Now build the remaining nodes and branches
	FOR(i,(int)Seq2Add.size()) {
		// Get space nodes and branches; Do some checking
		while(T.BraLink(BraNum,0) != -1 || T.BraLink(BraNum,1) != -1)	{ BraNum++; }
		while(T.NoLinks(InNode) != 0)									{ InNode++; }
		assert(T.NoLinks(Seq2Add[i]) == 0 && T.BraLink(Seq2Add[i],0) == -1 && T.BraLink(Seq2Add[i],1) == -1);
		// Do the nodes
		Lnk.clear(); Lnk.push_back(InNode); T.ReplaceNodeLink(Seq2Add[i],Lnk);
		Lnk.clear(); Lnk.push_back(Seq2Add[i]); Lnk.push_back(-1); Lnk.push_back(-1); T.ReplaceNodeLink(InNode,Lnk);
		Lnk.clear(); Lnk.push_back(Seq2Add[i]); T.ReplaceNodeBra(Seq2Add[i],Lnk);
		Lnk.clear(); Lnk.push_back(Seq2Add[i]); Lnk.push_back(BraNum); Lnk.push_back(-1); T.ReplaceNodeBra(InNode,Lnk);
		// Assign node types
		T.AssignNodeType(Seq2Add[i],leaf);
		T.AssignNodeType(InNode,branch);
		// Move pointers on
		BraNum++; InNode++;
	}
	T.ForceReady(); T.GetStart(); T.BuildBraLinks(false);
	// Optimise the tree
	M->MakeNewTree(&T,true);
	// Organise the calculation type
	if(M->IsRMSDCalc()) { Error("\nCannot do GetSAStartTree for RMSD calculations...\n"); }
	if(DoParsimony == 1) { M->DoParsimony(); }
	else { M->DoLikelihoodCalc(); }
	// New faster resampling...
	switch(DoLazy) {
		case fullopt:	CurlnL = FullOpt(M,true,true); break;
		case lazy:		CurlnL = LazyOpt(M,true,true); break;
		case randomtree:	CurlnL = -BIG_NUMBER; break;
	};
	// Add the sequences back
	AddSequences(DoOutput,&Seq2Add,M,M->m_pTree,0,false,DoLazy);
	// Do output
	if(DoOutput) { cout << "\n\n--------------------------------------------------------------------------\n"; }
	// Reset the calculation type
	if(OriMLCalc) { M->DoLikelihoodCalc(); }
	else { M->DoParsimony(); }
}

/////////////////////////////////////////////////////////////////////////////////////
// This is the driving routine for the bionj random sampler
// ---------------------------------------------------------------
string DoRandomBioNJ(CBaseModel *M, CData *D, bool DoJC)	{
	vector <double > NewDists;
	string RetVal;
	D->Bootstrap();
	NewDists = GetPW(M,NULL,false);
	RetVal = DoBioNJ(NewDists,D->m_vsName,true);
	D->CleanBootstrap();
	return RetVal;
}

/////////////////////////////////////////////////////////////////////////////////////
// This is the driving routine for the RMSD based sampler
// ---------------------------------------------------------------
// This randomly picks leaves in proportion to the RMSD fit between tree and observed distances
//  and removes them.
// SNAP is then used to perform local rearrangements
// Quasi-stepwise addition is used to propose a new tree
// //////////////////////////////////////////////////////////
// The function that controls how sequences are removed from the tree is in two parts
// Pr(SeqRemoved) = Normalised(f(RMSD) * g(NumberOfTimesSeqSampled))
//  f(RMSD) is the normalised RMSD values, so larger RMSD are more frequently removed from the tree
//  g(x) is 1 - normalised frequency of times a sequence has been pick; more frequently picked sequences are less likely to be chosen
// Arguments:
// DoLazy: whether to do lazy optimisation
// ProbParsimony: probability of doing a parsimony stepwise addition
// DoOutput: whether to do output
void GetRMSDSampTree(CBaseModel *M, CTreeList *Four, CTreeList *Five, CTreeList *Six, LazyType DoLazy, double ProbParsimony,bool DoOutput)	{
	int i,j,NodStore = 0,Seq, DespCount = 0, Seq2Remove, NearestT, RandomCount = 0;
	int DespThreshold = DESP_THRESHOLD;
	double RanNum, BestLnL, ProbSubOptBra = 0.0;
	vector <int> SeqRem,TreeCheck;
	vector <double> ProbSeq, ProbSeqRem;
	static bool SetInit = true;
	static vector <int> NumSeqRem;			// Number of times each sequence has been removed
	static vector <int> BestTreeVal;		// The unique values associated with the current best tree (used to reset NumSeqRem)
	CTree OriTree, *T = M->m_pTree;
	OriTree = *T;
	bool flag = true, Treefirst = true, OriMLCalc = M->IsLikelihoodCalc();

	// RMSD isn't allowed for resampling
	if(M->IsRMSDCalc()) { Error("\nCannot do RMSD resampling of trees..."); }

	// Set the Probability of a suboptimal solution
	ProbSubOptBra = min(max (0.0, (double) OptObs * DESPERATE_SUBOPT_PROB),0.33);

	// Do first Opt
	if(!DoLazy) { BestLnL = FullOpt(M); } else { BestLnL = LazyOpt(M); }
	// On first call to function initialise SeqRem
	if(SetInit == true) {
		BestTreeVal = T->ConstOut();
		FOR(i,T->NoSeq()) { NumSeqRem.push_back(0); } FOR(i,T->NoSeq()) { ProbSeqRem.push_back(1.0); }
		SetInit = false;
	} else {
	// Otherwise do a check to see whether NumSeqRem needs to be reset
		TreeCheck = T->ConstOut();
		if(Compare(&TreeCheck,&BestTreeVal) == false) {
			BestTreeVal = TreeCheck;
			FOR(i,T->NoSeq()) { NumSeqRem[i] = 0; }
	}	}

	ProbSeqRem.resize(T->NoSeq());
	while(flag == true)	{
		if(DoOutput) { cout << "\n\n---------- Resampling treespace using RMSD based random plucking ----------\n"; }
		if(Treefirst == false) { *T = OriTree; } else { Treefirst = false; }
		////////////////////////////////////////////////////////////////////////
		// Code for random plucking based on RMSD
		// Define probability of removing particular sequences
		ProbSeq = OriTree.MeasureSeqStability(PWDists);
		ProbSeq = NormaliseVector(ProbSeq);
		// Adjust the probabilities of sequence removal for how often a sequence has been sampled
		FOR(i,M->m_pTree->NoSeq()) { ProbSeqRem[i] = (double) NumSeqRem[i]; }
		ProbSeqRem = NormaliseVector(ProbSeqRem);
		FOR(i,M->m_pTree->NoSeq()) { ProbSeq[i] *= 1- ProbSeqRem[i]; }
		ProbSeq = NormaliseVector(ProbSeq);
#if DEBUG_ADDSEQ == 1
		cout << "\nProcessed probs = " << ProbSeq << " Sum = " << Sum(ProbSeq);
#endif
		SeqRem.clear();
		if(DoOutput) { cout << "Removed sequences:"; }
		// Randomly remove some sequences from the tree
		Seq2Remove = 0;
		FOR(i,T->NoSeq()) {
			RanNum = Random();
			if(RanNum < PROB_RAN_SEQ_REM) { Seq2Remove++; }
		}
		if(Seq2Remove < 3) { Seq2Remove = 3; }
		FOR(i,Seq2Remove)	{
			if(T->NoSeq() - SeqRem.size() == 3) { break; }
			// If a sequence has been chosen to be removed randomly choose which one w.r.t. RMSD
			Seq = DiscreteRand(&ProbSeq);
			// Debug check to ensure it hasn't already been removed...
			FOR(j,(int)SeqRem.size()) { if(SeqRem[j] == Seq) { Error("Trying to multiply remove a sequence..."); } }
			NumSeqRem[Seq]++;
			if(DoOutput) { cout << " " << Seq+1; }
			SeqRem.push_back(Seq);
			T->RemoveLeafNode(Seq);
			// Reset the probabilities
			ProbSeq[Seq] = 0;
			ProbSeq = NormaliseVector(ProbSeq);
		}
		random_shuffle(SeqRem.begin(),SeqRem.end());
		// Set whether to do parsimony resampling
		if(Random() < ProbParsimony) {
			if(DoOutput) { cout << "\nDoing parsimony stepwise addition"; }
			M->DoParsimony(); DespThreshold = PARS_DESP_THRESHOLD;
		} else {
			if(DoOutput) { cout << "\nDoing likelihood stepwise addition"; }
			M->DoLikelihoodCalc();
		}
		// Find node to start calculation from and get rid of any compression
		T->GetStart(); M->CleanFastCalc(true);
		RanNum = Random();
		if(RanNum < PROB_DO_IMPSNAP && M->IsParsimony() && DoLazy != randomtree)	{
			// ReSNAP the reduced tree to try and do big rearrangements
			// If only three sequences then just do optimisation
			if(T->NoSeq() - SeqRem.size() == 3)	{
				if(DoLazy = fullopt) { BestLnL = FullOpt(M); } else if(DoLazy = lazy) { BestLnL = LazyOpt(M); }
			// Otherwise do switch routine
			} else {
				if(DoOutput) { cout << "\nResampled space: now doing SNAP on subtree" << flush; }
				TreeSNAP(Four,Five,Six,M,false);
			}	} else { if(!DoLazy) { BestLnL = FullOpt(M); } else { BestLnL = LazyOpt(M); } }
		// Add them back again in a random order
		AddSequences(DoOutput,&SeqRem,M,T,ProbSubOptBra,true,DoLazy);
		// Do MP hill-climb if required
		if(DO_MP_HILLCLIMB > Random() && M->IsParsimony()) { DoMPHillClimb(M,Four,Five,Six); }
		NearestT = MinTabuRFDist(T);
		if(NearestT < TABU_RADIUS && NearestT >= 0) { flag = true; } else { break; }
		if(flag == true) {
			if(DoOutput) { cout << " -- min RF dist " << NearestT * 2<< " .. Tree tabu"; }
			PROB_RAN_SEQ_REM *= 1.1; if(PROB_RAN_SEQ_REM > 1.0) { PROB_RAN_SEQ_REM = 1.0; }
			// If the resampling isn't desperate then just resample as before
			if(DespCount < DespThreshold) { DespCount++; ProbSubOptBra += DESPERATE_SUBOPT_PROB;
			} else { // Otherwise random SNAP
				if(DoOutput) { cout << "\n\n---------- Become desperate: using random SNAPs ----------";  }
				DespCount = 0;	// Reset DespCount for progression to random trees
				RandomCount = 1;
				if(DoOutput) { cout << "\n\tProcessing: "; }
				while(flag ==true) {
					if(DoOutput) { cout << "." << flush; }
					*T  = OriTree;
					if(DespCount < 5) { // Do some random SNAPs
						int RandStepSize = (RandomCount++) + (int) ( (double) (T->NoNode() - T->NoSeq()) * (double) PROB_NODE_SNAP);
						GetNewTree(RandStepSize,T,Four,Five,Six,DoLazy);
						if(DO_MP_HILLCLIMB > Random()) { DoMPHillClimb(M,Four,Five,Six); }
						NearestT = MinTabuRFDist(T);
						if(NearestT < TABU_RADIUS) { flag = true; } else { flag = false; break; }
						if(RandomCount == 25) { if(DoOutput) { cout << " forced break" << flush; } flag = false; break; }
					} else { // Randomly sample tree space
						cout << "\nVery desperate! Random tree: ";
						GetSAStartTree(T, M, false, randomtree, false);
						NearestT = MinTabuRFDist(T);
						cout << "RF = " << NearestT;
						if(NearestT < TABU_RADIUS && NearestT >= 0) { flag = true; } else { break; }
					}
					DespCount++;
	}	}	}	}
	T->GetStart();
	T = NULL;
	if(DoOutput) { cout << " --- min RF dist " << NearestT * 2<< " ... Tree accepted!" << flush; }
	if(OriMLCalc) { M->DoLikelihoodCalc(); } else { M->DoParsimony(); }
}

void AddSequences(bool DoOutput, vector <int> *SeqRem, CBaseModel *M,CTree *Tree, double DespProb, bool AllowRandomSnap, LazyType DoLazy)	{
	int i;
	double BestlnL;
	/////////////////////////////////////////////////////////////
	// The advanced sequence adding approach
	vector <STreelnL*> CandTreeList;	// List of candidate trees
	STreelnL *StartTree;

	if(DoOutput) { cout << "\nAdding sequences: "; }
	// Initialise to starting tree
	StartTree = new STreelnL;
	StartTree->Tree = *M->m_pTree;
	if(DoLazy != randomtree) { StartTree->lnL = M->lnL(); } else { StartTree->lnL = -BIG_NUMBER; }
	CandTreeList.push_back(StartTree); StartTree = NULL;
	////////////////////////////////////////////////////////////////////
	// Loop through sequences to add
	FOR(i,(int)SeqRem->size())	{
#if DEBUG_ADDSEQ == 1
		cout << "\n\nAdding sequence " << i+1 << "/" << SeqRem->size() << " (ProbSubOptTree = " << DespProb << ");";
#endif
		if(DoOutput && i%20 == 0) { cout << "\n\t"; }
		if(DoOutput) { cout << " " << SeqRem->at(i) +1 << flush; }
		GetSATree(&CandTreeList,SeqRem->at(i),M,DespProb, DoLazy);
	}
	// Find minimum likelihood
	BestlnL = 0;
	FOR(i,(int)CandTreeList.size())	{ if(CandTreeList[i]->lnL < BestlnL) { BestlnL = CandTreeList[i]->lnL; } }
	// Check whether all trees in list are tabu
	FOR(i,(int)CandTreeList.size())	{ if(IsTabu(&CandTreeList[i]->Tree,false,false,M) == false) { *Tree = CandTreeList[i]->Tree; break; } }
	*Tree = CandTreeList[0]->Tree;
	// Clean up after this
	FOR(i,(int)CandTreeList.size())	{ delete CandTreeList[i]; CandTreeList[i] = NULL; }
	CandTreeList.clear();
	StartTree = NULL;
	return;
}

// Highest level Stepwise Addition routine
void GetSATree(vector <STreelnL *> *Trees, int SeqAdd, CBaseModel *M,double ProbSubOptBra, LazyType DoLazy)	{
	vector <STreelnL*> List, BestTrees;
	vector <int> NC;
	int TreeNum,i,j,Node;
	int StorePerRound = SA_STORE_PER_ROUND_ML;
	CTree T(M->m_pData);
	SBestModel BestM;
	double lnL2Beat;
	if(M->IsParsimony()) { StorePerRound = SA_STORE_PER_ROUND_MP; }
#if DEBUG_ADDSEQ == 1
	cout << "\nProducing " << Trees->size() << " lists from tree" << flush;
	FOR(i,Trees->size()) { Trees->at(i)->Tree.OutBra(); cout << "\n\tTree["<<i<<"]: " << Trees->at(i)->Tree << " = " << Trees->at(i)->lnL; }
#endif
	FOR(TreeNum,(int)Trees->size())	{
		// Get the start tree
		M->MakeNewTree(&Trees->at(TreeNum)->Tree,true);
		if(M->IsLikelihoodCalc()) {
			switch(DoLazy) {
				case fullopt:	Trees->at(TreeNum)->lnL = FullOpt(M); break;
				case lazy:		Trees->at(TreeNum)->lnL = LazyBraOpt(M); break;
				case randomtree:	Trees->at(TreeNum)->lnL = Random(); break;
			};
		} else { Trees->at(TreeNum)->lnL = M->lnL(); }

#if DEBUG_ADDSEQ == 1
		M->m_pTree->OutBra(); Trees->at(TreeNum)->Tree.OutBra();
		cout << "\n\nDoing tree: " << TreeNum << ": " << Trees->at(TreeNum)->lnL << " == " << M->lnL(true) << "\nOldTree: " << Trees->at(TreeNum)->Tree << "\nNewTree: " << *M->m_pTree  << flush;
#endif
		// Do calculation from i
		BraSACalc(SeqAdd,&List,M,false,NULL,-1,-1,DoLazy);
#if DEBUG_ADDSEQ == 1
		cout << "\n\tHave starting list["<<List.size()<<"]";   FOR(i,List.size()) { List[i]->Tree.OutBra(); cout << "\n\t\t" << i << "\t" << List[i]->Tree << " = " << List[i]->lnL; }
#endif
		// Sort The list
		SortList(&List);
		// Now fully optimise the required branches
		FOR(i,(int)List.size()) { if(i==SA_OPT_NUM_ACC) { break; } List[i]->DoOpt = true; }
		// Redo the calculation
		lnL2Beat = List[0]->lnL;
		BraSACalc(SeqAdd,&List,M,true,&lnL2Beat,-1,-1,DoLazy);
#if DEBUG_ADDSEQ == 1
		cout << "\n\tHave list["<<List.size()<<"]";   FOR(i,List.size()) { List[i]->Tree.OutBra(); cout << "\n\t\t" << i << "\t" << List[i]->Tree << " = " << List[i]->lnL; }
#endif
		// Create the full optimised trees and store the them in the BestTrees list
		Node = M->m_pTree->NodeLink(SeqAdd,0);		// Find a spare internal node to use for tree reconstruction
		FOR(i,min((int)List.size(),SA_OPT_NUM_ACC) ) {
			T = *M->m_pTree;
			T.AddSeq(SeqAdd,List[i]->BraNum);
			Node = T.NodeLink(SeqAdd,0);
			FOR(j,T.NoLinks(Node))	{ T.SetB(T.FindBra(Node,List[i]->BraLinks[j]),List[i]->Tree.B(j)); }
			List[i]->Tree = T;
			BestTrees.push_back(List[i]); List[i] = NULL;
		}
		// Clean the List
		CleanList(&List);
	}
	// The BestTrees vector is now sorted and the best StorePerRound trees are used as the basis of the next round of SA
	SortList(&BestTrees);
#if DEBUG_ADDSEQ == 1
		cout << "\n\tHave list of BestTrees[size="<<BestTrees.size()<<"]";   FOR(i,BestTrees.size()) {
			cout << "\n\t\t" << i << "\t" << BestTrees[i]->Tree << " = " << BestTrees[i]->lnL; }
#endif
	CleanList(Trees);
	// Store the best trees
	if(ProbSubOptBra < DBL_EPSILON)	{ // Do normal storage
		FOR(i,min(StorePerRound,(int)BestTrees.size())) { Trees->push_back(BestTrees[i]); BestTrees[i] = NULL; }
	} else {	// Otherwise accept poorer trees
		vector <double> Probs;
		double Cur = 1.0, RanNum;
		FOR(i,(int)BestTrees.size()) { Probs.push_back(Cur); if(i == 1) { Cur = ProbSubOptBra; } }
		FOR(i,min(StorePerRound,(int)BestTrees.size()))	{
			RanNum = Random();
			Probs = NormaliseVector(Probs,1.0);
			Cur = 0.0; FOR(j,(int)Probs.size()) { Cur += Probs[j]; if(Cur > RanNum) { Probs[j] = 0.0; break; } }
			Trees->push_back(BestTrees[j]);
			BestTrees[j] = NULL;
		}
	}
#if DEBUG_ADDSEQ == 1
	Trees->at(0)->Tree.OutBra();
	cout << "\n\tHave list of trees for next round [size="<<Trees->size()<<"]";   FOR(i,Trees->size()) { cout << "\n\t\t" << i << "\t" << Trees->at(i)->Tree << " = " << Trees->at(i)->lnL; }
#endif
	FOR(i,(int)BestTrees.size() ) { if(BestTrees[i] != NULL) { delete BestTrees[i]; } BestTrees[i] = NULL; } BestTrees.clear();
}

//////////////////////////////////////////////////////////////////////////
// Fast Step-wise addition
// This steps through a tree and calculates the SA of all the branches
void BraSACalc(int SAdd, vector <STreelnL*> *TreeList, CBaseModel *M,bool DoTripleOpt,double *lnL2Beat, int SPR, int OriBr, LazyType DoLazy, CTree *OriT)	{
	int Nod;
	string tree;
//	cout << "\n--- Into BraSACalc, SAdd: " << SAdd << " ---";
	// Get starting likelihood
	if(SPR == -1) { M->m_pTree->GetStart(); }
	// Do full likelihood C to initialise Bk and Fd space
	if(M->IsLikelihoodCalc() && DoLazy != randomtree) { M->PrepareFastCalc(); M->lnL(true); }
	// Check entry conditions (slightly weaker than the AddSeq check...)
	assert(!M->m_pTree->NodeNull(SAdd) && SAdd < M->m_pTree->NoSeq());
	assert(M->m_pTree->NodeType(SAdd) == leaf && M->m_pTree->NoLinks(SAdd) != 0);
	Nod = M->m_pTree->NodeLink(SAdd,0);
	assert(M->m_pTree->NodeType(Nod) == branch && M->m_pTree->NoLinks(Nod) == 3);
	assert(M->m_pTree->NodeLink(Nod,1) == -1 && M->m_pTree->NodeLink(Nod,2) == -1);
	// Run the stepwise addition
	Branch_SA(SAdd,TreeList,-1,M->m_pTree->StartCalc(),-1,M,DoTripleOpt,lnL2Beat,SPR,OriBr,DoLazy,OriT);
}

// In-order tree traversal for SA
void Branch_SA(int SAdd, vector <STreelnL*> *TreeList, int First,int NTo, int NFr, CBaseModel *M,bool DoTripleOpt, double *lnL2Beat, int SPR,int OriBr, LazyType DoLazy, CTree *OriT)	{
	int i;
//	cout << "\n>> Branch_SA: NTo: " << NTo << ", NFr: " << NFr;
	// Always perform the calculations in the first place
	if(M->m_pTree->NodeType(NTo) == leaf)	{ // Do the leaf calculations (Should only be the first calc)
		if(NFr == -1) { // Its the starting node
			DoSA(SAdd,TreeList,First,NTo,M->m_pTree->NodeLink(NTo,0),M->m_pTree->NodeBra(NTo,0),M,true,DoTripleOpt,lnL2Beat,SPR,OriBr,DoLazy,OriT);
		} else { // For other nodes
			DoSA(SAdd,TreeList,First,NFr,NTo,M->m_pTree->NodeBra(NTo,0),M,true,DoTripleOpt,lnL2Beat,SPR,OriBr,DoLazy,OriT);
		}
	} else { // Do the internal calculations
		FOR(i,M->m_pTree->NoLinks(NTo))	{ if(M->m_pTree->NodeLink(NTo,i) == NFr || M->m_pTree->NodeLink(NTo,i) == -1) { break; } }
		assert(i != M->m_pTree->NoLinks(NTo));
		// If the node from isn't a leaf node do the internal calculation (i.e. avoids first node)
		if(M->m_pTree->NodeType(M->m_pTree->NodeLink(NTo,i)) != leaf)	{
			DoSA(SAdd,TreeList,First,NTo,NFr, M->m_pTree->NodeBra(NTo,i), M,false,DoTripleOpt,lnL2Beat,SPR,OriBr,DoLazy,OriT);
	}	}
	// Do the looping
	First = 0;
	FOR(i,M->m_pTree->NoLinks(NTo))	{
		if(M->m_pTree->NodeLink(NTo,i) == NFr || M->m_pTree->NodeLink(NTo,i) == -1) { continue; }
		Branch_SA(SAdd,TreeList,First,M->m_pTree->NodeLink(NTo,i),NTo,M,DoTripleOpt,lnL2Beat,SPR,OriBr,DoLazy,OriT);
		First = 1;
}	}

//  0) First = The number of the visit to this node (0 = first, 1 = second, -1 = First branch on tree)
//	1) NTo = Node the likelihood is directed to
//	2) NFr = Node the likelihood is directed from
//	3) Br  = Branch being examined

void DoSA(int SAdd, vector <STreelnL*> *TreeList, int First,int NTo, int NFr, int Br, CBaseModel *M, bool IsExtBra,bool DoTripleOpt, double *lnL2Beat, int SPR, int OriBr, LazyType DoLazy, CTree *OriT)	{
	int ListNum;
	bool Flag = false;
	double T = M->m_pTree->B(M->m_pTree->FindBra(NTo,NFr)), lnL;
	static CTree Triplet("(1:0.1,2:0.1,3:0.1);",3);
	vector <int> L;
	STreelnL *TL = NULL;

	if(OriT != NULL)    {
		int BrDist = OriT->BranchDist(Br,OriBr);
		if(BrDist < SPR_LOWBOUND || BrDist > SPR_UPBOUND) {
			if(M->IsLikelihoodCalc() && DoLazy != randomtree)	{
				if(Br < 3) { M->PreparePT(Br); }
				if(IsExtBra) { M->Leaf_update(NTo,NFr,Br,M->m_pTree,First); } else { M->Bran_update(NTo,NFr,Br,M->m_pTree,First); }
			}
			return;
	}	}

//	cout << "\nDoSA; NTo: " << NTo << ", NFr: " << NFr << ", SPR: " << SPR;

	// Firstly check whether the sequence is in the list, and if it is whether it needs optimising
	FOR(ListNum,(int)TreeList->size()) { if(TreeList->at(ListNum)->BraNum == Br && TreeList->at(ListNum)->SPRID == SPR) { if(TreeList->at(ListNum)->DoOpt == true) { L = TreeList->at(ListNum)->BraLinks; Flag = true; } break; } }
	// If its new then add it to the list
	if(ListNum == (int)TreeList->size()) {
		TL = new STreelnL; TL->BraNum = Br; TL->BraLinks.push_back(NTo); TL->BraLinks.push_back(NFr); TL->BraLinks.push_back(SAdd);	L = TL->BraLinks; TL->SPRID = SPR;
		Flag = true;
	}
	// If calculations are required then do them
	if(Flag && OriBr != Br)	{
		if(M->IsParsimony()) {
			if(SPR != -1) { Error("\nShouldn't do SPR tree DoSA for parsimony...");  }
			CTree TempTree; TempTree = *M->Tree();
			M->Tree()->AddSeq(SAdd,Br);
			lnL = M->lnL();
//			cout << "\nScore: " << lnL;
			if(ListNum == (int)TreeList->size()) {
				TL->lnL = lnL;
				TL->Tree = Triplet; TL->DoOpt = false;
				TreeList->push_back(TL);
				TL = NULL;
			} else { TreeList->at(ListNum)->lnL = lnL; TreeList->at(ListNum)->Tree = Triplet; }
			// Restore the old tree
			*M->Tree() = TempTree;
		} else if(DoLazy != randomtree) {
			//////////////////////////////////////////////////////////////////
			// Triplet calculations
			// 1. Set up the partial likelihoods and mapping in Model
//			cout << "\nDoing Triplet calculation: " << L;
			M->PrepareTripletCalc(L,&Triplet,min(SPR,OriBr));
			M->ApplySubTree(&Triplet,false,false);
			// 2. Initialise the branches
			Triplet.SetB(0,T/2); Triplet.SetB(1,T/2); Triplet.SetB(2,0.1);
			// 3. Do the calculations and store the result
			// Note for optimisation only branches are ever optimised...
			lnL = M->lnL();
			// Don't optimise bloody awful trees (defined as within 1.5% or 5% of likelihood depending on optimisation strength)
//			cout << "\nTripleCalc: lnL: " << lnL;
			if(lnL2Beat != NULL) {
//				cout << " cf. lnL2Beat: " << *lnL2Beat << "; ratio: " << lnL / *lnL2Beat << flush;
				if(
						((lnL / *lnL2Beat) > FULL_SPR_RATIO && DoLazy == fullopt) ||
						((lnL / *lnL2Beat) > LAZY_SPR_RATIO && DoLazy == lazy)
					) { DoTripleOpt = false; } }
//			DoTripleOpt = false;
			if(DoTripleOpt) {
//				cout << "\nCompetitive optimise: " << DoLazy << " from " << M->lnL(); if(lnL2Beat != NULL) { cout << " (vs " << *lnL2Beat << ")... "; } cout << flush;
				switch(DoLazy) {
					case fullopt:		lnL = M->FastBranchOpt(lnL,0.001,NULL,10,false); break;
					case lazy:			lnL = M->FastBranchOpt(lnL,1.0,NULL,3,false); break;
					case randomtree:	lnL = Random(); break;
				};
			}else { if(DoLazy == randomtree) { lnL = Random(); } else { lnL = M->lnL(); } }
//			cout << " -> " << M->lnL(true); // << " cf. lnL: " << lnL << flush;
//			cout << " Diff: " << *lnL2Beat - M->lnL(true);
			if(ListNum == (int)TreeList->size()) {
				TL->lnL = lnL;
				TL->Tree = Triplet; TL->DoOpt = false;
				TreeList->push_back(TL);
				TL = NULL;
			} else { TreeList->at(ListNum)->lnL = lnL; TreeList->at(ListNum)->Tree = Triplet; }
			// 4. Clean up
			M->CleanCPMapping();
//			cout << " ... lnL: " << lnL;
		} else { // Do a random sample likelihood
			CTree TempTree; TempTree = *M->Tree();
			M->Tree()->AddSeq(SAdd,Br);
			if(ListNum == TreeList->size()) {
				TL->lnL = Random();
				TL->Tree = Triplet; TL->DoOpt = false;
				TreeList->push_back(TL);
				TL = NULL;
			} else { TreeList->at(ListNum)->lnL = Random(); TreeList->at(ListNum)->Tree = Triplet; }
			// Restore the old tree
			*M->Tree() = TempTree;
		}
	}
	// Finally update the space
	if(M->IsLikelihoodCalc() && DoLazy != randomtree)	{
		if(Br < 3) { M->PreparePT(Br); }
		if(IsExtBra) { M->Leaf_update(NTo,NFr,Br,M->m_pTree,First); } else { M->Bran_update(NTo,NFr,Br,M->m_pTree,First); }
	}
}


//////////////////////////////////////////////////////////////////////////
// Crude tree resampling routine using multiple SNAPS
//////////////////////////////////////////////////////////////////////////

void GetNewTree(int Steps,CTree *InTree,CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, LazyType DoLazy)	{
	int ChangeNum;
	CTreeList *TL = NULL;
	vector <int> NodeFrom, NodeCov, LeafList;
	if(Steps < MIN_NUM_RANSNAPS) { Steps = MIN_NUM_RANSNAPS; }
	// Do the random snaps
	FOR(ChangeNum,Steps)	{
		LeafList = InTree->NodeCP(RandInt(InTree->NoSeq(),InTree->NoNode()),2,&NodeFrom,&NodeCov);
		switch(LeafList.size()) {
		case 4: TL = FourSp; break;
		case 5: TL = FiveSp; break;
		case 6: TL = SixSp; break;
		default: Error("\nUnknown number of leaves in GetNewTree");
		}
		InTree->ReplaceTreeCP(&TL->m_Trees[RandInt(0,(int)TL->m_Trees.size())],LeafList,NodeCov);
	}
	TL = NULL;
}

// Some STreelnL useful structure functions
void SortList(vector <STreelnL *> *L)	{
	int i,j;
	if(L->size() == 1) { return; }
	STreelnL *TempTL;
	vector <int> Good;
	// Sort The list
	TempTL = NULL;
	FOR(i,(int)L->size())	{ for(j=i;j<(int)L->size();j++)	{ if(L->at(j)->lnL > L->at(i)->lnL) { TempTL = L->at(i); L->at(i) = L->at(j); L->at(j) = TempTL; TempTL = NULL; } } }
	// Sort out the case when multiple members of the list have the same score
	if(L->at(0)->lnL == L->at(1)->lnL)	{
		vector <STreelnL *> Temp;
		FOR(i,(int)L->size()) { if(L->at(0)->lnL < L->at(i)->lnL) { break; } Good.push_back(i); }
		random_shuffle(Good.begin(),Good.end());
		FOR(i,(int)Good.size()) { Temp.push_back(L->at(Good[i])); }
		FOR(i,(int)Good.size()) { L->at(i) = Temp[i]; Temp[i] = NULL; }

	}
	assert(TempTL == NULL);
}
STreelnL GetBestList(vector <STreelnL *> *L) {
	int i,best = -1;
	double bestlnL = -BIG_NUMBER;
	FOR(i,(int)L->size()) {
		if(L->at(i)->lnL > bestlnL) { bestlnL = L->at(i)->lnL; best = i; }
	}
	return *L->at(best);
}

void CleanList(vector <STreelnL *> *L)	{
	int i;
	// Clean the List
	FOR(i,(int)L->size()) { if(L->at(i) != NULL) { delete L->at(i); } L->at(i) = NULL; }
	L->clear();
}
void OutSTreelnL(vector <STreelnL*> *TL, bool FullDetail, ostream &os)	{
	int i;
	FOR(i,(int)TL->size())	{
		os << "\nTree["<<i<<"] ID(Br:"<< TL->at(i)->BraNum << ",SPRID:"<<TL->at(i)->SPRID<<"): " << TL->at(i)->Tree << " = " << TL->at(i)->lnL;
		if(TL->at(i)->DoOpt == true) { os << " * "; }
		if(FullDetail) { os << "; score: " << TL->at(i)->score << "; IsOpt: " << (bool) TL->at(i)->IsOpt << "; DoOpt: " << (bool) TL->at(i)->DoOpt; }
}	}

