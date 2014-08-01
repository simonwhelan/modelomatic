// Optimise.cxx
//
// Contains:
//
// 1) The SNAP optimisation routine for each node (inc. variants)
// 2) The dfpmin numerical optimisation routine
//

#include "optimise.h"
#include "TreeList.h"
#if IS_LNL_CALC
extern CPhyloDat PhyDat;
extern bool WarningMulD;
extern vector <double> PWDists;
extern vector <STabuTree> TabuTrees;
extern int TABU_RADIUS, OptObs, EXIT_OBS;
extern bool ALLOW_PREDICTLNL;
#else
CPhyloDat PhyDat;
bool WarningMulD;
extern vector <double> PWDists;
vector <STabuTree> TabuTrees;
extern int TABU_RADIUS, EXIT_OBS, OptObs;
#define DO_NOT_RUN 1
#endif
const int NumResamp = 1;

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter;
	extern int LFunc_Log_Counter, SubLFunc_Log_Counter,SPR_Log_Counter;
#endif

#define SPR_DEBUG 0				// Debugging of SPR routine
#define DEBUG_MULD_OPT 0		// [0: None; 1: Gentle; 2: Hard] Debugging of MulD_Optimiser routine (outputs stuff to opt.output)
#define DEBUG_MULD_OPT_SEP_FILES 0	// Whether to output stuff to opt.output.Run#
int DebugOptNum=0;
const string OptFile = "opt.output";
#define ALLOW_ONLYPAROPT 1		// Whether the DoOnlyParOpt(...) function is used

/////////////////////////////////////////////////////////////////////////////////
// Standard likelihood optimising routine
double FullOpt(CBaseModel *Model, bool DoPar, bool DoBra, bool DoFreq, double CurlnL,bool FullLikAcc, int NoIterations, double lnL2Beat, double lnL_tol, bool DoOutput,bool TightFullOpt)	{
	int i;
//	cout << "\nInto FullOpt for Model = <" << Model->m_sName << ">  DoPar: " << DoPar << ", DoBra: " << DoBra << ", DoFreq: " << DoFreq << endl << flush;
//	if(DoOutput) { cout << "\nOutput"; } else { cout << "\nNo output"; }
	double ACC = lnL_tol, gtol = FULL_GTOL;
	bool OnlyBra = false;
	// Deal with parsimony
	if(Model->IsParsimony()) { return Model->lnL(true); }
	// Do everything else
	if(FullLikAcc == false) { ACC = PART_LIK_ACC; }
	if(Model->IsRMSDCalc()) { ACC = RMSD_ACC; }
#if DO_NOT_RUN == 1
	Error("\nSet IS_LNL_CALC to 1 in Leaphy.h file when compiling programs that do likelihood computation\n");
#endif
	// Get the parameter pointers.
	vector <double *> vPar = Model->GetOptPar(DoBra,DoBra,DoPar,DoFreq);
	if( (!DoPar && DoBra) || (DoBra && vPar.size() == Model->Tree()->NoBra()) ) { OnlyBra = true; }
	Model->PrepareFastCalc();
	// If rqd, get the initial likelihood
	if(fabs(CurlnL + BIG_NUMBER) < DBL_EPSILON)	{
		CurlnL = Model->lnL();
		// Ensure this is a good starting likelihood
		if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) {
			Model->FixSmallBranches();
			CurlnL = Model->lnL();
			if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) { Model->m_pTree->OutBra(); cout << "\nTree:\n" << *Model->m_pTree << "\nLikelihood: " << Model->lnL(); Error("\nFullOpt(...) failed due to weird likelihood..."); }
	}	}
	// Get the optimised likelihood
	if(ACC > gtol) { gtol = ACC; }
//	cout << "\nInto optimiser for " << Model->Name() << " = " << CurlnL << " cf. " << Model->lnL() << " cff. " << Model->lnL(true) << flush;
//	if(fabs(CurlnL - Model->lnL(true)) > 0.001) { cout << "\nAnd it's broken already...\n\n"; exit(-1); }
//	cout << "\nModel: " << Model->Name() << ", Options: [" << DoPar << "," << DoBra << "," << DoFreq << "," << CurlnL << "," << FullLikAcc << "," << NoIterations << "," << lnL2Beat << "," << lnL_tol << "]";
//	DoOutput = true; cout << "\nOptimising: " << flush; FOR(i,(int)vPar.size()) { cout << "\t" << *vPar[i] << flush; }
	if(!vPar.empty()) {
//		cout << "\nAttempting to optimise..." << flush;
		if(NoIterations == DEFAULT_OPTNUM) { NoIterations = Model->OptNum(); }
		CurlnL = MulD_Optimise(CurlnL,gtol,ACC,vPar,Model,NoIterations,DoOutput,OnlyBra,2,true,-lnL2Beat,7,true,TightFullOpt); }
	// Clean up and return
	FOR(i,(int)vPar.size()) { vPar[i] = NULL; }
	FOR(i,(int)Model->m_vpAllOptPar.size()) { Model->m_vpAllOptPar[i]->GlobalApply(); }
	Model->Tree()->OutBra();
//	cout << "\nlnL: " << Model->lnL();
//	cout << "\nTree: " << *Model->Tree();
	assert(CurlnL < 0);
	return CurlnL;
}

/////////////////////////////////////////////////////////////////////////////////
// Standard likelihood optimising routine but only does a bit...
double LazyOpt(CBaseModel *Model, bool DoPar, bool DoBra, bool DoFreq, double CurlnL,bool FullLikAcc, int NoIterations, double lnL2Beat, bool DoOutput)	{
	int i;
	double ACC = FULL_LIK_ACC; double GTOL = FULL_GTOL;
	bool OnlyBra = false;
	// Deal with parsimony
	if(Model->IsParsimony()) { return Model->lnL(true); }
	// Do everything else
	if(FullLikAcc == false) { ACC = PART_LIK_ACC; GTOL = PART_GTOL; }
	if(Model->IsRMSDCalc()) { ACC = RMSD_ACC; }
#if DO_NOT_RUN == 1
	Error("\nSet IS_LNL_CALC to 1 in Leaphy.h file when compiling programs that do likelihood computations");
#endif
	// Get the parameter pointers.
	vector <double *> vPar = Model->GetOptPar(DoBra,DoBra,DoPar,DoFreq);
	if( (!DoPar && DoBra) || (DoBra && vPar.size() == Model->Tree()->NoBra()) ) { OnlyBra = true; }
	Model->PrepareFastCalc();
	// If rqd, get the initial likelihood
	if(fabs(CurlnL + BIG_NUMBER) < DBL_EPSILON)	{
		CurlnL = Model->lnL();
		// Ensure this is a good starting likelihood
		if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) {
			Model->FixSmallBranches();
			CurlnL = Model->lnL();
			cout << "\nFinished new likelihood..." << CurlnL << flush;
			if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) { Model->m_pTree->OutBra(); cout << "\nTree:\n" << *Model->m_pTree << "\nLikelihood: " << Model->lnL(); Error("\nLazyOpt(...) failed due to weird likelihood..."); }
	}	}
	// Get the optimised likelihood
	CurlnL = MulD_Optimise(CurlnL,GTOL,ACC,vPar,Model,NoIterations,DoOutput,OnlyBra,2,true,-lnL2Beat,2,false,false);
//	CurlnL = MulD_Optimise(CurlnL,GTOL,ACC,vPar,Model,NoIterations,DoOutput,OnlyBra,2,true,-lnL2Beat,2);
	// Clean up and return
	FOR(i,(int)vPar.size()) { vPar[i] = NULL; }
//	if(CurlnL > 0) { cout << "\n --- RETURNING CurlnL > 0 ???"; }
	assert(CurlnL < 0);
	return CurlnL;
}

////////////////////////////////////////////////////////////////////////////////
// Preoptimiser optimising routine. Runs a simpler model to get good
//  starting parameters for a more complex model
double PreOpt(CBaseModel *Model,bool DoOutput)	{
	double likelihood;
	CBaseModel *PreOptModel = Model->PreOptModel();
	if(PreOptModel == NULL)	{ return -BIG_NUMBER; }			// If no preoptimisation required, don't do any.
	PreOpt(PreOptModel,DoOutput);
//	cout << "\nDoing pre-opt with model: " << PreOptModel->Name();
	likelihood = PreOptModel->lnL();
/*	cout << "\nBefore";
	int i; FOR(i,PreOptModel->NoPar()) {
		if(PreOptModel->m_vpPar[i]->Name().find("Alpha") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
		if(PreOptModel->m_vpPar[i]->Name().find("P(Inv)") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
		if(PreOptModel->m_vpPar[i]->Name().find("Sigma") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
	}
	cout << "\nLikelihood " << likelihood;
*/	likelihood = FullOpt(PreOptModel,true,true,false,likelihood,true,100,-BIG_NUMBER,FULL_LIK_ACC,DoOutput);
/*	cout << "\nAfter";
	FOR(i,PreOptModel->NoPar()) {
		if(PreOptModel->m_vpPar[i]->Name().find("Alpha") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
		if(PreOptModel->m_vpPar[i]->Name().find("P(Inv)") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
		if(PreOptModel->m_vpPar[i]->Name().find("Sigma") != string::npos) { cout << "\nPar["<<i<<"]: " << *PreOptModel->m_vpPar[i]; }
	}
	cout << "\nLikelihood " << likelihood;
*/	Model->ApplyPreOptModel(PreOptModel);
	delete PreOptModel;
	return likelihood;
}

/////////////////////////////////////////////////////////////////////////////////
// Likelihood optimising routine for branches; doesn't work rigorously
double LazyBraOpt(CBaseModel *Model, double CurlnL,int NoIterations, double Tol)	{
	int i;
	double ACC = FULL_LIK_ACC;
	bool OnlyBra = false;
	// Deal with parsimony
	if(Model->IsParsimony()) { return Model->lnL(true); }
#if DO_NOT_RUN == 1
	Error("\nSet IS_LNL_CALC to 1 in Leaphy.h file when compiling programs that do likelihood computations");
#endif
	// Get the parameter pointers.
	vector <double *> vPar = Model->GetOptPar(true,true,false,false);
	Model->PrepareFastCalc();
	// If rqd, get the initial likelihood
	if(fabs(CurlnL + BIG_NUMBER) < DBL_EPSILON)	{
		CurlnL = Model->lnL();
		// Ensure this is a good starting likelihood
		if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) {
			Model->FixSmallBranches();
			CurlnL = Model->lnL();
			cout << "\nFinished new likelihood..." << CurlnL << flush;
			if(fabs(CurlnL+BIG_NUMBER) < DBL_EPSILON) { Model->m_pTree->OutBra(); cout << "\nTree:\n" << *Model->m_pTree << "\nLikelihood: " << Model->lnL(); Error("\nLazyOpt(...) failed due to weird likelihood..."); }
	}	}
	// Get the optimised likelihood
	CurlnL = Model->FastBranchOpt(CurlnL,Tol,NULL,NoIterations);
	// Clean up and return
	FOR(i,(int)vPar.size()) { vPar[i] = NULL; }
	return CurlnL;
}

///////////////////////////////////////////////////////////////////////////////
// Standard Fitch-Margoliash branch length optimising routine
double FMOpt(CBaseModel *Model, vector <double> D)	{
	int i;
	bool DoingRMSDCalc = Model->IsRMSDCalc();
	double RMSD_score;
	Model->DoRMSDCalc();		// Tell model to do RMSD calcs
	RMSD_score = Model->lnL();	// Initialise score
	vector <double *> RMSD_Par = Model->GetOptPar(true,true,false,false);	// Get branches parameters only
	RMSD_score = MulD_Optimise(RMSD_score,FULL_GTOL,RMSD_ACC,RMSD_Par,Model,50,false,false,2);				// Do optimisation
	if(!DoingRMSDCalc) { Model->DoLikelihoodCalc();	} // Tell model to return to calculations
	FOR(i,(int)RMSD_Par.size()) { RMSD_Par[i] = NULL; }
	return RMSD_score;
}


/* **************************************************************************** */
//			GENERAL ROUTINES FOR HILL CLIMBING
/* **************************************************************************** */
double DoMPHillClimb(CBaseModel *M, CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, double LimitMPScore, bool DoOutput)	{
	int count = 0;
	bool OriML = M->IsLikelihoodCalc();
	double BestPars = -BIG_NUMBER, CurPars;
	// Check some entry conditions
	if(M->IsRMSDCalc()) { Error("\nCannot run DoMPHillClimb(...) for RMSD tree\n\n"); }
	// Do some initialisation
	if(OriML) { M->DoParsimony(); } else { M->ResetParsimony(); }

	CurPars = M->lnL();
	if(DoOutput) { cout << "\nDoing parsimony hillclimb\n-------------------------\nInitial: " << -CurPars; }
	while(CurPars - BestPars > LimitMPScore)	{
		BestPars = CurPars;
		CurPars = TreeSNAP(FourSp,FiveSp,SixSp,M,false,lazy);
		CurPars = DoSimpleSPR(M,false,CurPars);
		if(DoOutput) { if(count++ %5 == 0) { cout << "\n\t"; } else { cout << " -> "; } cout << CurPars; }
	}
	if(OriML) { M->DoLikelihoodCalc(); }
	return BestPars;
}

/* **************************************************************************** */

//			ROUTINES TO DO SNAP REARRANGEMENTS

/* **************************************************************************** */
///////////////////////////////////////////////////////////////////////////////
// Routine to do the SNAP algorithm
// --------------------------------
// Returns optimal likelihood (e.g double value less than 0) when new trees
// If tree already examined, returns 1;
int RMSD_SUBSET = LOOSE_RMSD_SUBSET, PARS_SUBSET = LOOSE_PARS_SUBSET, OPT_NUM_ACC = LOOSE_OPT_NUM_ACC;
double PERCENT_CHANGE = LOOSE_PERCENT_CHANGE;

double TreeSNAP(CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, CBaseModel *Model, bool DoOutput, LazyType OriDoLazy)	{
	bool ChangeFlag = true, *DoNodes, LastRun = false, SubChange, HadKnot = false, DoFullOpt = false, HadTabu = false, DoTheTabu, DoSNAP = true, AllowFast = true;
	LazyType DoLazy = OriDoLazy;
	int i, OutCount=0,NumSp = Model->m_pTree->NoSeq(), NumSubTree = 0, RunCount = 0, OptCount = 0, NoChange, ActualChange, NodeNum, VisibleRunCount = 1, NumSkipped, KnotCount;
	RMSD_SUBSET = LOOSE_RMSD_SUBSET; PARS_SUBSET = LOOSE_PARS_SUBSET; PERCENT_CHANGE = LOOSE_PERCENT_CHANGE; OPT_NUM_ACC = LOOSE_OPT_NUM_ACC;
	double Cur_lnL,Best_lnL, last_lnL = -BIG_NUMBER, LastParsScore = -BIG_NUMBER,CurParsScore;
	SNAPchange ChangeType;
	vector <TreeArrange> TA;								// Vector of potential tree arrangements
	vector <TreeArrange>::iterator TA_i;					// Iterator for topologies
	vector <vector <int> > Knot_Nodes;						// Each vector<int> contains mutually exclusive "knots" in the tree for considering seperately.
	vector <bool> boolTA(Model->Tree()->NoNode(),false);

#if DO_QUICK_SNAP == 1
	vector <double> Stability, NodeScores;
	vector <int> NodeOrder(Model->Tree()->NoSeq()-2,0);
	int j,NodeCount;
#endif

	CTree InTree, BestTree,TestTree,FirstTree,nodeBestTree;	// Too many tree objects?
	CTreeList *pTLst;									// Tree list passed
	double nodeBestlnL;
	InTree = *Model->m_pTree;
	GET_MEM(DoNodes,bool,NumSp-2); FOR(i,NumSp-2) { DoNodes[i] = true; }
	// Check entry conditions
	if(!Model->CheckSameTree()) { Error("Model contains mixed trees... SNAP can't cope with this (yet...)\n\n"); }
	// Initialise
	if(DoLazy == lazy && Model->IsParsimony()) { DoLazy = fullopt; }
	BestTree = InTree; // Optimise the model and initialise some variables
	switch(DoLazy) {
		case fullopt:	Best_lnL = Cur_lnL = FullOpt(Model); break;
		case lazy:		Best_lnL = Cur_lnL = LazyOpt(Model); break;
		case randomtree:	Error("\nHaven't sorted random in TreeSnap yet"); break;
	};
	// Do some output
	if(DoOutput) {
		cout << "\n\n~ ~ ~\n\n\tPerforming SNAP"; if(ALLOW_SNAP_SKIP != 0) { cout << " (allow skipping)"; } cout << ": " << Best_lnL << endl;
		cout << "\t\tKey: . = no improvement; + = better tree found; * = resolving star\n\t\t     ^ = fast step; _ = fast step failed; ~ = no fast step\n\t\t     t = tabu tree; | = Skipping node; L = Two leaf nodes so skipping\n";
	}
	///////////////////////////////////////////////////////////////////////////////
	// Loop to apply tree estimation algorithm
	while(ChangeFlag== true)	{
		// Never allow too many snaps
		if(RunCount >= max(100,NumSp) && HadKnot == true) { if(DoOutput) { cout << "\n\t\tBreaking SNAP because of high numbers of rearrangements"; } break; }
		DoFullOpt = false; ChangeFlag = false;	// Set flags
		if(RunCount % 5 == 0) { AllowFast = true; }
//		if(NumSp <= 6) { DoFullOpt = true; AllowFast = false; }	// Do full search for small numbers of sequences
		// 1.) Loop through node list and examine sub-trees
		if(DoOutput) { cout << "\n\tRun " << VisibleRunCount << " [" << NumSp - 2 << "]:\t" << flush; }
		// Set some counters and stuff
		nodeBestTree = InTree; nodeBestlnL = Best_lnL;
		OutCount=0; NoChange = 0; NumSkipped = 0;
		FOR(i,Model->m_pTree->NoNode()) { boolTA[i] = false; }
		// Set the DoNodes if examining Knots
		if(Knot_Nodes.size() > 0)	{
			if(DoOutput) { cout << "\n\tDoing knotfix:\t"; }
			FOR(i,NumSp-2) { DoNodes[i] = false; }
			FOR(i,(int)Knot_Nodes[0].size())	{ if(Knot_Nodes[0][i] < NumSp) { continue; } DoNodes[Knot_Nodes[0][i] - NumSp] = true; }
		}
		// This is to identify when parsimony is no longer a good predictor of the likelihood function
		if(Model->IsLikelihoodCalc() && PARS_SUBSET < 105) {
			CurParsScore = Model->GetFullParsimony();
			if(LastParsScore - CurParsScore >= 0) { PARS_SUBSET = 106; OPT_NUM_ACC = FULL_OPT_NUM_ACC; }
			LastParsScore = CurParsScore;
		}

#if DO_QUICK_SNAP == 1
		// Node order stuff for fast calculations...
		if(Model->IsParsimony()) { FOR(i,NumSp-2) { NodeOrder[i] = i+NumSp; } }
		else {
			Stability = Model->Tree()->MeasureNodeStability(PWDists);
			GetParsimonySNAPList(&Stability,FourSp,FiveSp,SixSp,Model);
			NodeScores = Stability;
			RSort(&NodeScores);
			FOR(i,NumSp - 2) {
				FOR(j,NumSp - 2) {
					if(fabs(Stability[j] - NodeScores[i]) < DBL_EPSILON) {
						NodeOrder[i] = j + NumSp; NodeScores[i] = -1.0; break;
		}	}	}	}
		// See if a big step can be made using parsimony
		DoSNAP = true;
		if(Model->IsLikelihoodCalc() && !IsSameTree(Model->Tree(),&InTree) && AllowFast) {
			Cur_lnL = LazyBraOpt(Model,-BIG_NUMBER,5);
			if(Cur_lnL > Best_lnL) {
				if(DoOutput) { cout << "^ " << flush; for(i=1;i<(int)NodeOrder.size();i++) { if(i%80==0 && i > 1) { break; } if(i%20 == 0) { cout << " "; } cout << " "; } }
				DoSNAP = false;
				ChangeFlag = true;
			}
			else {
				if(DoOutput) { cout << "_" << flush; }
				AllowFast =false;
				*Model->Tree() = InTree;
				Cur_lnL = Best_lnL;
		}	} else { *Model->Tree() = InTree; if(DoOutput) { cout << "~" << flush; } }
#endif
		/////////////////////////////////////////////////
		// Enter the main loop
#if DO_QUICK_SNAP == 1
		FOR(NodeCount,(int)NodeOrder.size())	{
			NodeNum = NodeOrder[NodeCount];
			LastRun = true;
#else
		for(NodeNum = NumSp;NodeNum < Model->m_pTree->NoNode(); NodeNum++)	{
#endif
			// If the big step has been made using parsimony then don't bother with SNAP
			if(!DoSNAP) { break; }
			// Some output and counters
			if(DoOutput && OutCount % 20 == 0 && OutCount > 0) { cout << " "; }
			if(DoOutput && OutCount % 80 == 0 && OutCount > 0) { cout << "\n\t\t\t"; }
			OutCount++;
			// For subtrees: skip nodes with -1 in their links
			if(Model->m_pTree->IsNodeLink(NodeNum,-1)) { continue; }
			if(DoNodes[NodeNum - InTree.NoSeq()] == true || ALLOW_SNAP_SKIP == 0)	{
				// i) Skip nodes that only run to two leaf nodes. They'll be covered by other rearrangements
#if DO_LEAF_SNAP == 0
				if(Model->m_pTree->NoLeafLink(NodeNum) > 1) { if(DoOutput) { cout << "L"; } continue; }
#endif
				// ii) Prepare tree for calculation
				Model->PrepareNodeCP(NodeNum,2);
				switch(Model->NumLeafCP())	{
				case 4: pTLst = FourSp; break;
				case 5: pTLst = FiveSp; break;
				case 6: pTLst = SixSp; break;
				default: Error("Unlisted type of tree..."); break;
				};
				// i) Estimate sub-tree likelihoods
				if(Model->IsLikelihoodCalc() || Model->IsRMSDCalc()) { 	ChangeType = DoScoreSubTree(NodeNum, &Best_lnL, pTLst, Model, &TA,DoLazy); }
				else { ChangeType = DoSimpleSubTree(NodeNum, &Best_lnL, pTLst, Model, &TA,DoLazy); }
				switch(ChangeType)	{
				case normal:	if(DoOutput) { cout << "+" << flush; } SubChange = true; boolTA[NodeNum] = true; break;
				case star:		if(DoOutput) { cout << "*" << flush; } SubChange = true; boolTA[NodeNum] = true; break;
				case tabu:		if(DoOutput) { cout << "t" << flush; } SubChange = false; HadTabu = true; break;
				case none:		if(DoOutput) { cout << "." << flush; } SubChange = false; break;
				};

				if(SubChange)	{
					NoChange++;
					ChangeFlag = true;
					if(Best_lnL > nodeBestlnL)	{	// If the new tree is better than any others store it
						nodeBestTree = InTree; nodeBestlnL = Best_lnL;
				}	}
#if DO_QUICK_SNAP == 1
			if(!Model->IsParsimony() && (ChangeType == normal || ChangeType == star)) {
				if(DoOutput) { for(i=NodeCount+1;i<(int)NodeOrder.size();i++) { if(i%80==0) { break; } if(i%20 == 0) { cout << " "; } cout << " "; } }
				break;
			}
#endif
			} else { if(DoOutput) { cout << "|"; NumSkipped++; } }
		}
		// Clean the leaf mapping
		Model->CleanCPMapping();
		// Make changes to the tree and specify which nodes to look at next
		// Only nodes that are within 4 of changed nodes are examined again (i.e. <= 2 SNAPS)
		ActualChange = 0;
		if(!TA.empty() || !DoSNAP)	{	// TA was sorted in DoSubTree
			if(DoSNAP) {
#if DO_QUICK_SNAP == 1
				FOR(i,NodeCount) { DoNodes[NodeOrder[i] - NumSp] = false; }
#else
				// Set all nodes to false
				FOR(i,NumSp-2) { DoNodes[i] = false; }
#endif
				// All nodes that can change are upgraded to true
				FOR(i,(int)TA.size()) { 	DoNodes[TA[i].NodeNum - Model->m_pData->m_iNoSeq] = true; }
				// Do some rearrangements
				while(!TA.empty())	{
					TestTree = *Model->m_pTree;	// Get a copy of the current tree
					// Make the change
					TestTree.ReplaceTreeCP(&TA[0].Tree,TA[0].LeafList,TA[0].NodCovList,false);
					// If the tree hasn't been seen before accept it and
					// Remove all other NodeNum (centre of SNAP) within 5 steps of the best
					if(IsTabu(&TestTree,true,false,Model) == false || DO_TABU == 0)	{
						// Set to true all DoNodes within 4 of changed node
						FOR(i,NumSp-2)	{ if(Model->m_pTree->NodeDist(TA[0].NodeNum,Model->m_pTree->NoSeq()+i) < 5) { DoNodes[i] = true; } }
						// Apply the change to the tree
						Model->m_pTree->ReplaceTreeCP(&TA[0].Tree,TA[0].LeafList,TA[0].NodCovList,false);
						ActualChange++;
						// If its the first change then store it in FirstTree (in case combinations of moves fail)
						if(ActualChange == 1) {  FirstTree = *Model->m_pTree; }
						// Remove the TA that overlap it.
						for(TA_i = TA.begin(), TA_i++; TA_i < TA.end(); TA_i++)	{ // If within the distance then cut it
							if(InTree.NodeDist(TA[0].NodeNum,TA_i->NodeNum) < 6)	{ TA.erase(TA_i); TA_i--;}
					}	}
					// Remove the original one now.
					TA.erase(TA.begin());
				}
			} else { TA.clear(); NoChange = 1; } // If(!DoSNAP) then note a change has been made
			Model->CleanCPMapping();
			// Assess the new topology and store the tree
			if(DoLazy != fullopt) {	// Do full optimisation
					if(RunCount % STEPS_BETWEEN_PAROPT ==0 && RunCount > 0) { Cur_lnL = LazyOpt(Model); }
					else													{ Cur_lnL = LazyBraOpt(Model); }
			} else {								// Do quick optimisation
				Cur_lnL = FullOpt(Model,false,true);
			}
			// If this likelihood is worse then try resorting to only a single change in topology
			if(Best_lnL > Cur_lnL + FULL_LIK_ACC) { *Model->m_pTree = FirstTree; Cur_lnL = FullOpt(Model); DoFullOpt = true; }
			if(DoOutput) { cout << " : " << Cur_lnL << flush; }
		} else { if(DoOutput) { cout << " : " << Cur_lnL << flush; } }
		InTree = *Model->m_pTree;	// Update the current tree
		// 4.) Output current tree
		if(DoOutput)	{
			ofstream out(PhyDat.Out().c_str(), ios::app);
			out.setf(ios::fixed); out.precision(4);
			InTree.OutBra();
			out << "\n" << VisibleRunCount << ":=\t" <<  Cur_lnL << "\t" << BestTree << flush;
			out.close();
		}
		////////////////////////////////////////////////////////////////////
		// Manage exit conditions
		////////////////////////////////////////////////////////////////////
		// Before anything, lets knock on some counters
		RunCount++; VisibleRunCount++;
		// 0. If doing a full opt and tabu dist is < 2 then exit
#if DO_TABU == 1
		if(MinTabuRFDist(Model->m_pTree,true) < 2) {
			if(DoOutput) { cout << " skipping opt..."; }
			ChangeFlag = false; DoFullOpt = false; break;
		}
#endif
		// 1. If no changes and some nodes have been skipped then check them (except for knots)
		if(NoChange == 0 && LastRun == false && ALLOW_SNAP_SKIP != 0 && NumSkipped > 0 && Knot_Nodes.empty()) {
			// Check the marked nodes before exiting
			FOR(i,NumSp - 2) { if(DoNodes[i] == true) { DoNodes[i] = false; } else { NoChange++; DoNodes[i] = true; } }
			LastRun = ChangeFlag = true;
			if(KnotCount == -1) { if(DoOutput) { cout << " too knotted"; } break; }
			if(DoOutput) { cout << " R"; }
			continue;
		}
		// 2. Otherwise, if no changes then exit; not if examining knots
		else if(NoChange == 0) {
			if(Knot_Nodes.empty()) {	// If not working on knots
				if(DoLazy == lazy || DoLazy == randomtree) {
					// If already seen tabu trees and very close to them return before proper optimisation
#if DO_TABU == 1
					if(MinTabuRFDist(Model->m_pTree,true) < 2 || HadTabu) {
						if(DoOutput) { cout << " skipping opt..."; }
						ChangeFlag = false; DoFullOpt = false; break;
					}
#endif
					// Otherwise there might be another good tree nearby
					DoLazy = fullopt; FOR(i,NumSp - 2) { DoNodes[i] = true; }
					Best_lnL = Cur_lnL = FullOpt(Model); DoFullOpt = true;
					if(DoOutput) { cout << " tightening opt = " << Best_lnL << flush; }
					// Tighten the set of trees examined...
					if(PARS_SUBSET < FULL_PARS_SUBSET) { PARS_SUBSET = FULL_PARS_SUBSET; }
					RMSD_SUBSET = FULL_RMSD_SUBSET;
					OPT_NUM_ACC = FULL_OPT_NUM_ACC;
					ChangeFlag = true;
					continue;
				}
				if(DoOutput) { cout << "  New optima!"; } break;
#if DO_TABU == 1
				IsTabu(Model->Tree(),true,true,Model,Best_lnL,true);
#endif
			}			// Break if no knots
			else { Knot_Nodes.erase(Knot_Nodes.begin()); if(DoOutput) { cout << "  Knot resolved!"; } KnotCount = -1; FOR(i,NumSp-2) { DoNodes[i] = true; } ChangeFlag = true; continue; }	// Remove Knot if no possible changes have been found
		}
		// 3. If an change has been made and an improvement in likelihood is visible
		else if(Cur_lnL - Best_lnL > FULL_LIK_ACC) {
			if(DoOutput) { cout << " +" << ActualChange << " = " << Cur_lnL - Best_lnL; }
#if DO_TABU == 1
			i = MinTabuRFDist(Model->m_pTree,true);
			if(i>=0 && DoOutput && !Model->m_pTree->IsCutTree() && OldTabuTrees()) { cout << " ("<<2*MinTabuRFDist(Model->m_pTree,true) << ")"; }
#endif
			if(Knot_Nodes.empty()) { KnotCount = 0; }
			if(Cur_lnL - Best_lnL > FULL_LIK_ACC * 10) { OptCount =0; }
			LastRun = false;
		}
		// 4. Otherwise if there is a change, but no improvement in likelihood we have knots
		else if(fabs(Cur_lnL - Best_lnL) <= 1.0E-4 && Model->IsLikelihoodCalc())	{
			HadKnot = true;
			if(DoOutput) { cout << " +" << ActualChange << " Knot"; }
			if(KnotCount == -1) { cout << " and stuck"; break; }		// Leave if the tree is hopelessly knotted
			if(Knot_Nodes.empty())	{ // New knot untangling code
				Knot_Nodes = Model->Tree()->GetKnotClusters(boolTA);
				KnotCount = 0;
			} else {
				if(++KnotCount > 10) {
					if(DoOutput) { cout << " limit reached"; }
					Knot_Nodes.erase(Knot_Nodes.begin());
					if(Knot_Nodes.empty()) { FOR(i,NumSp-2) { DoNodes[i] = true; } KnotCount = -1; continue; }
					if(DoOutput) { cout << " (trying next)"; }
					KnotCount = 0; continue;
			}	}
		// 5. Otherwise the likelihood has decreased:
		} else {
			// Otherwise give up... It looks like an error
			if(DoOutput) { cout << " -" << ActualChange << " bad tree"; }
			Cur_lnL = Best_lnL; *Model->m_pTree = BestTree;
			// If already tried the reverse the exit
			if(LastRun == true) { if(DoOutput) { cout << " giving up"; } break; }
			// Check the marked nodes before exiting
			NoChange = 0;
			FOR(i,NumSp - 2) { if(DoNodes[i] == true) { DoNodes[i] = false; } else { NoChange++; DoNodes[i] = true; } }
			if(NoChange == 0) { if(DoOutput) { cout << " done" << flush; } break; }
			LastRun = ChangeFlag = true;
			if(DoOutput) { cout << " R"; }
			continue;
		}
		/////////////////////////////////////////////////////
		// Test whether tree is tabu; if so and there is skipping, then check other directions
#if ALLOW_TABU_PART_OPT == 1
		DoTheTabu = true; if(Model->IsParsimony()) { DoTheTabu = false; }
		if(IsTabu(Model->m_pTree,true,DoTheTabu,Model,Cur_lnL)) {
#else
		DoTheTabu = true; if(Model->IsParsimony() || !DoFullOpt) { DoTheTabu = false; }
		if(IsTabu(Model->m_pTree,true,DoFullOpt,Model,Cur_lnL)) {
#endif
			if(DoOutput) { cout << " Tabu"; } break;
		}	// Exit if complete
		Best_lnL = Cur_lnL; BestTree = *Model->m_pTree;
	}	// Clear up memory
	DEL_MEM(DoNodes);
	if(DoOutput) {
		if(DoFullOpt)	{ cout << "\n\n\tReturning best lnL: " << Best_lnL << endl; }
		else			{ cout << "\n\n\tReturning optimised lnL: " << Best_lnL << endl; }
	}
	return Best_lnL;
}

////////////////////////////////////////////////////////////////////
// Function that gets the proposed SNAP changes
void GetParsimonySNAPList(vector <double> *Score,CTreeList *FourSp,CTreeList *FiveSp,CTreeList *SixSp, CBaseModel *Model)	{
	bool ChangeFlag = true, LastRun = false;
	int ActualChange, NumSp = Model->m_pTree->NoSeq(), NumSubTree = 0, RunCount = 0, OptCount = 0, NodeNum;
	if(Score->size() != NumSp - 2) { Error("\nNeed to pass Score.size() == NumSp - 2 to GetParsimonySNAPList(...)\n\n"); }
	if(Model->IsParsimony() || Model->IsRMSDCalc()) { return; }
	// Some more initialisation
	Model->DoParsimony();
	vector <TreeArrange> DoTA,TA;								// Vector of potential tree arrangements
	vector <TreeArrange>::iterator DoTA_i;					// Iterator for topologies
	CTree TestTree,InTree;
	InTree = *Model->Tree();
	CTreeList *pTLst;									// Tree list passed
	double Best_lnL = Model->lnL(), OriMP = Best_lnL;
	for(NodeNum = NumSp;NodeNum < Model->m_pTree->NoNode(); NodeNum++)	{
		// Some output and counters
		// For subtrees: skip nodes with -1 in their links
		if(Model->m_pTree->IsNodeLink(NodeNum,-1)) { continue; }
		// i) Skip nodes that only run to two leaf nodes. They'll be covered by other rearrangements
#if DO_LEAF_SNAP == 0
		if(Model->m_pTree->NoLeafLink(NodeNum) > 1) { Score->at(NodeNum - NumSp) = 0.05 * Random(); continue; }
#endif
		// ii) Prepare tree for calculation
		Model->PrepareNodeCP(NodeNum,2);
		switch(Model->NumLeafCP())	{
		case 4: pTLst = FourSp; break;
		case 5: pTLst = FiveSp; break;
		case 6: pTLst = SixSp; break;
		default: Error("Unlisted type of tree..."); break;
		};
		// i) Estimate sub-tree likelihoods
		DoSimpleSubTree(NodeNum, &Best_lnL, pTLst, Model, &TA,fullopt);
		if(!TA.empty()) { Score->at(NodeNum - NumSp) = (100 * (TA[0].lnL - OriMP)) + Random(); DoTA.push_back(TA[0]); }
		TA.clear();
	}

	// Do some rearrangements to get a single tree to propose...
	Model->CleanCPMapping(); ActualChange = 0;
	while(!DoTA.empty())	{
		TestTree = *Model->m_pTree;	// Get a copy of the current tree
		// Make the change
		TestTree.ReplaceTreeCP(&DoTA[0].Tree,DoTA[0].LeafList,DoTA[0].NodCovList,false);
		// If the tree hasn't been seen before accept it and
		// Remove all other NodeNum (centre of SNAP) within 5 steps of the best
		if(IsTabu(&TestTree,true,false,Model) == false)	{
			// Apply the change to the tree
			Model->m_pTree->ReplaceTreeCP(&DoTA[0].Tree,DoTA[0].LeafList,DoTA[0].NodCovList,false);
			ActualChange++;
			// Remove the DoTA that overlap it.
			for(DoTA_i = DoTA.begin(), DoTA_i++; DoTA_i < DoTA.end(); DoTA_i++)	{ // If within the distance then cut it
				if(InTree.NodeDist(DoTA[0].NodeNum,DoTA_i->NodeNum) < 6)	{ DoTA.erase(DoTA_i); DoTA_i--;}
		}	}
		// Remove the original one now.
		DoTA.erase(DoTA.begin());
	}
	Model->DoLikelihoodCalc();
}



/////////////////////////////////////////////////////////////////////
// Do subtree search based on something really simple
// --------------------------------------------------
// This does searching when the branches don't need optimising
// Currently this is only done for parsimony
SNAPchange DoSimpleSubTree(int NodeNum, double *OrilnL, CTreeList *TL, CBaseModel *Model, vector <TreeArrange> *TA, LazyType DoLazy)	{
	int i,j,BestTree = -1, OriTreeNum,MaxTabuDist, SNAP_AWAY_TABU_RADIUS = 2  * TABU_RADIUS;
	bool flag = true, tflag = false, okay = true, IsStar = false, HaveTabu = false;
	double lnL, lnL2Beat = *OrilnL, baseRMSD = 0,BestLnL = -BIG_NUMBER;
	vector <STreeScore> lnL_score;
	vector <STreeScore>::iterator vI, vOri;
	STreeScore Score;
	vector <int> StarFix;
	CTree TestTree,OldTree;
	vector <TreeArrange>::iterator TA_i;
#if ENCOURAGE_SNAP_AWAY == 0
	SNAP_AWAY_TABU_RADIUS = 0;
#endif
	// Initialisation
	Model->BuildOriSubTree(&OldTree);
	TestTree = *Model->m_pTree;
	// Find which of the trees is the original
	FOR(i,(int)TL->m_Trees.size()) {
		if(IsSameTree(&OldTree,&TL->m_Trees[i])) { OriTreeNum = i; }
		Score.TreeNum = i;
		// Form the tree
		Model->Tree()->ReplaceTreeCP(&TL->m_Trees[i],Model->LeafMap(),Model->NodesCovered(),false);
		// Get the score
		Score.Lik = Model->lnL(true);
		lnL_score.push_back(Score);
		// Restore the tree
		*Model->Tree() = TestTree;
	}
	SortSTreeScore(&lnL_score);
	////////////////////////////////////////////////////////////////////////////
	// Make sure the original tree is put last of trees of the same likelihood
	vOri = lnL_score.end() - 1;
	IFOR(vI,lnL_score)	{
		// Get the placement of the original
		if(vI->TreeNum == OriTreeNum) { vOri = vI; continue; }
		// Do the swap if required
		if(vOri->Lik > vI->Lik && vOri != lnL_score.end() - 1) {
			vI--; j = vI->TreeNum; lnL = vI->Lik; vI->TreeNum = vOri->TreeNum; vI->Lik = vOri->Lik; vOri->TreeNum = j; vOri->Lik = lnL;
			break;
	}	}
	// If the best tree is the original, then return no change
	if(OriTreeNum == lnL_score[0].TreeNum) { return none; }
	////////////////////////////////////////////////////////////////////////////
	// Only do clever SNAP stuff (e.g. tabu, &c) for the full tree
	if(!Model->m_pTree->IsCutTree() && OldTabuTrees())	{
		MaxTabuDist = 0;
		// If required encourage SNAP away from already visited trees
		// Will only accept trees that improve the likelihood, but if faced with multiple choices
		//  it will choose the one that maximises the TabuDist *NOT* the likelihood!
		FOR(i,(int)lnL_score.size()) {
			// Once reached the original tree then stop
			if(lnL_score[i].TreeNum == OriTreeNum) { break; }
			// Form the tree
			TestTree = *Model->m_pTree;
			TestTree.ReplaceTreeCP(&TL->m_Trees[lnL_score[i].TreeNum],Model->LeafMap(),Model->NodesCovered(),false);
			// Get the TabuDist
			lnL_score[i].TabuDist = MinTabuRFDist(&TestTree,true);
			if(lnL_score[i].TabuDist == 0) { HaveTabu = true; }
			if(IsTabu(&TestTree,true,false,Model)) { lnL_score[i].TabuDist = -1; }
			// If the tree isn't tabu then add
			if(lnL_score[i].TabuDist > MaxTabuDist) {		// Catch new best TabuDist
				// Always take the highest likelihood when ENCOURAGE_SNAP_AWAY == 0 or after OptObs > EXIT_OBS - 3
				if(MaxTabuDist != 0 && (ENCOURAGE_SNAP_AWAY == 0 || OptObs < EXIT_OBS - COMP_SNAP_AWAY_2_OPT_OBS)) { break;  }
				BestLnL = lnL_score[i].Lik; StarFix.clear(); StarFix.push_back(i); MaxTabuDist = lnL_score[i].TabuDist;
				if(lnL_score[i].TabuDist > SNAP_AWAY_TABU_RADIUS) { break; }
			} else if(lnL_score[i].TabuDist > 0 && lnL_score[i].TabuDist == MaxTabuDist && lnL_score[i].Lik +  FULL_LIK_ACC > BestLnL) {
				// If the best likelihood then store only it
				if(lnL_score[i].Lik > BestLnL) { BestLnL = lnL_score[i].Lik; StarFix.clear(); StarFix.push_back(i); }
				// Else add the tree to StarFix
				if(fabs(BestLnL-lnL_score[i].Lik) < FULL_LIK_ACC) { StarFix.push_back(i); }
		}	}
		if(MaxTabuDist == 0) { return tabu; } // Skip tabu trees
		if(StarFix.size() > 1)	{ random_shuffle(StarFix.begin(),StarFix.end()); }
		BestLnL = lnL_score[StarFix[0]].Lik; BestTree = lnL_score[StarFix[0]].TreeNum;
//		cout << "\n-----------------------------------------------------------------\nTrees:"; FOR(i,lnL_score.size()) { cout << "\n\tTree["<<i<<"]: " << TL->m_Trees[lnL_score[i].TreeNum] << " = " << lnL_score[i].Lik << "("<< lnL_score[i].TabuDist << ")"; if(IsIn(i,StarFix)) { cout << " *In StarFix*"; } if(i==StarFix[0]) { cout << " !BEST!"; } } cout << "\nBest: " << TL->m_Trees[BestTree] << " = " << BestLnL << " cf. Original: " << OldTree << " = " << *OrilnL;
	} else { // Still have to avoid tabu trees
		MaxTabuDist = -1;
		FOR(i,(int)lnL_score.size()) {
			if(lnL_score[i].TreeNum == OriTreeNum) { if(HaveTabu) { return tabu; } else { return none; } }
			// Form the tree
			TestTree = *Model->m_pTree;
			TestTree.ReplaceTreeCP(&TL->m_Trees[lnL_score[i].TreeNum],Model->LeafMap(),Model->NodesCovered(),false);
			// Check whether tabu
			if(IsTabu(&TestTree,true,false,Model)) { continue; }
			BestLnL = lnL_score[i].Lik; BestTree = lnL_score[i].TreeNum; break;
		}
//		cout << "\n-----------------------------------------------------------------\nTrees:"; FOR(i,lnL_score.size()) { cout << "\n\tTree["<<i<<"]: " << TL->m_Trees[lnL_score[i].TreeNum] << " = " << lnL_score[i].Lik ; } cout << "\nBest: " << TL->m_Trees[BestTree] << " = " << BestLnL << " cf. Original: " << OldTree << " = " << *OrilnL;
	}
	// Store the best tree
	/////////////////////////////////////////////////////
	TreeArrange NewTA;
	if(!InRange(BestTree,0,(int)TL->m_Trees.size())) { Error("Error: trying to assign best tree too large"); }
	NewTA.lnL = BestLnL; NewTA.NodeNum = NodeNum;
	NewTA.Tree = TL->m_Trees[BestTree]; NewTA.TabuDist = MaxTabuDist;
	NewTA.LeafList = Model->LeafMap(); NewTA.NodCovList = Model->NodesCovered();
	if(TA->empty()) { TA->push_back(NewTA); }
	else	{
		pIFOR(TA_i,TA)	{
			// Do sorting by likelihood when not worrying about TABU dist
			if(ENCOURAGE_SNAP_AWAY == 0 || OptObs < EXIT_OBS - COMP_SNAP_AWAY_2_OPT_OBS || MaxTabuDist> SNAP_AWAY_TABU_RADIUS) { if(BestLnL > TA_i->lnL) { break; } }
			else { if(MaxTabuDist > TA_i->TabuDist || (MaxTabuDist == TA_i->TabuDist && BestLnL > TA_i->lnL) ) { break; } }
		}
		TA->insert(TA_i,NewTA);
	}
	if(StarFix.size() > 1) { return star; } return normal;
}

////////////////////////////////////////////////////////////////////////////////////////
// Section that defines how subtree searching is done
// ---------------------------------------------------
// DoScoreSubTree does calculations where the branches need optimising
// Currently this includes ML and RMSD style calculations
/////////////////////////////////////////////////////////////////////////////////////////

SNAPchange DoScoreSubTree(int NodeNum, double *OrilnL, CTreeList *TL, CBaseModel *Model, vector <TreeArrange> *TA, LazyType DoLazy)	{
	int i,j,BestTree = -1, OriTreeNum,MaxTabuDist, SNAP_AWAY_TABU_RADIUS = 2  * TABU_RADIUS, pars_set = PARS_SUBSET;
	bool flag = true, tflag = false, okay = true, IsStar = false, HaveTabu = false;
	double lnL, lnL2Beat = *OrilnL, baseRMSD = 0,BestLnL = -BIG_NUMBER;
	vector <STreeScore> lnL_score;
	vector <STreeScore>::iterator vI, vOri;
	STreeScore Score;
	vector <int> StarFix;
	CTree TestTree,OldTree;
	vector <TreeArrange>::iterator TA_i;
#if ENCOURAGE_SNAP_AWAY == 0
	SNAP_AWAY_TABU_RADIUS = 0;
#endif
	// Initialisation
	Model->BuildOriSubTree(&OldTree);
	// Find which of the trees is the original
	FOR(i,(int)TL->m_Trees.size()) { if(IsSameTree(&OldTree,&TL->m_Trees[i])) { OriTreeNum = i; break; } }
	// For larger trees check RMSD stuff
	if((int)TL->m_Trees.size() > RMSD_SUBSET && !Model->m_pTree->IsCutTree()) {
		// Loop through trees and get the RMSD stuff
		FOR(i,(int)TL->m_Trees.size()) {
			Score.TreeNum = i;
			Model->ApplySubTree(&TL->m_Trees[i],true);
			// Get the RMSD score
			Score.Lik = -TL->m_Trees[i].SubTreeRMSD(Model->LeafMap(),Model->NodesBelow(),PWDists,Model->PartialPW());
			lnL_score.push_back(Score);
		}
		SortSTreeScore(&lnL_score);
	} else if((int)TL->m_Trees.size() > PARS_SUBSET && !Model->m_pTree->IsCutTree()) {
		TestTree = *Model->m_pTree;
		Model->DoParsimony();
		// Find which of the trees is the original
		FOR(i,(int)TL->m_Trees.size()) {
			Score.TreeNum = i;
			// Form the tree
			Model->Tree()->ReplaceTreeCP(&TL->m_Trees[i],Model->LeafMap(),Model->NodesCovered(),false);
			// Get the score
			Score.Lik = Model->lnL(true);
			lnL_score.push_back(Score);
			// Restore the tree
			*Model->Tree() = TestTree;
		}
		lnL_score[OriTreeNum].Lik = -BIG_NUMBER;
		// Sort scores and ensure that all scores with the same best MP score are included
		SortSTreeScore(&lnL_score);
		FOR(i,(int)lnL_score.size()) { if(diff(lnL_score[0].Lik,lnL_score[i].Lik) != 0) { break; } }
		pars_set = max(i,PARS_SUBSET);
		// Reset the calculation type
		Model->DoLikelihoodCalc();

	} else { Score.Lik = 0.0; FOR(i,(int)TL->m_Trees.size()) { Score.TreeNum = i; lnL_score.push_back(Score); } }
	//////////////////////////////////////////////////////
	// Get the likelihood estimate guess
	//////////////////////////////////////////////////////
	lnL_score = DoSubTreeCalcs(*OrilnL,&lnL_score,TL,Model,min(RMSD_SUBSET,pars_set),false,DoLazy,PERCENT_CHANGE,false);
//	cout << "\nExamined " << lnL_score.size() << "/" << TL->m_Trees.size() << " approx_trees" << flush;

	//////////////////////////////////////////////////////
	// Do accurate likelihood calculations.
	lnL_score = DoSubTreeCalcs(*OrilnL,&lnL_score,TL,Model,OPT_NUM_ACC,true,DoLazy,BIG_NUMBER,true);
//	cout << " and " << lnL_score.size() << " accurate trees" << flush;

	// Add the original tree if required (original tree should alway be available for comparison)
	///////////////////////////////////////////////////////////////////////////////////////////////
	FOR(i,(int)lnL_score.size()) {
		if(lnL_score[i].TreeNum == OriTreeNum) {
			if(fabs(*OrilnL - lnL_score[i].Lik) > 1.0E-4) { lnL_score[i].Lik = *OrilnL; SortSTreeScore(&lnL_score); }
			break;
	}	}
	if(i == lnL_score.size()) {
		Score.Lik = *OrilnL; Score.TreeNum = OriTreeNum;
		lnL_score.push_back(Score); SortSTreeScore(&lnL_score);
	}

	////////////////////////////////////////////////////////////////////////////
	// Make sure the original tree is put last of trees of the same likelihood
	vOri = lnL_score.end() - 1;
	IFOR(vI,lnL_score)	{
		// Get the placement of the original
		if(vI->TreeNum == OriTreeNum) { vOri = vI; continue; }
		// Do the swap if required
		if(vOri->Lik > vI->Lik && vOri != lnL_score.end() - 1) {
			vI--; j = vI->TreeNum; lnL = vI->Lik; vI->TreeNum = vOri->TreeNum; vI->Lik = vOri->Lik; vOri->TreeNum = j; vOri->Lik = lnL;
			break;
	}	}
	// If the best tree is the original, then return no change
	if(OriTreeNum == lnL_score[0].TreeNum) { return none; }
	////////////////////////////////////////////////////////////////////////////
	// Only do clever SNAP stuff (e.g. tabu, &c) for the full tree
	if(!Model->m_pTree->IsCutTree() && OldTabuTrees())	{
		MaxTabuDist = 0;
		// If required encourage SNAP away from already visited trees
		// Will only accept trees that improve the likelihood, but if faced with multiple choices
		//  it will choose the one that maximises the TabuDist *NOT* the likelihood!
		FOR(i,(int)lnL_score.size()) {
			// Once reached the original tree then stop
			if(lnL_score[i].TreeNum == OriTreeNum) { break; }
			// Form the tree
			TestTree = *Model->m_pTree;
			TestTree.ReplaceTreeCP(&TL->m_Trees[lnL_score[i].TreeNum],Model->LeafMap(),Model->NodesCovered(),false);
			// Get the TabuDist
			lnL_score[i].TabuDist = MinTabuRFDist(&TestTree,true);
			if(lnL_score[i].TabuDist == 0) { HaveTabu = true; }
			if(IsTabu(&TestTree,true,false,Model)) { lnL_score[i].TabuDist = -1; }
			// If the tree isn't tabu then add
			if(lnL_score[i].TabuDist > MaxTabuDist) {		// Catch new best TabuDist
				// Always take the highest likelihood when ENCOURAGE_SNAP_AWAY == 0 or after OptObs > EXIT_OBS - 3
				if(MaxTabuDist != 0 && (ENCOURAGE_SNAP_AWAY == 0 || OptObs < EXIT_OBS - COMP_SNAP_AWAY_2_OPT_OBS)) { break;  }
				BestLnL = lnL_score[i].Lik; StarFix.clear(); StarFix.push_back(i); MaxTabuDist = lnL_score[i].TabuDist;
				if(lnL_score[i].TabuDist > SNAP_AWAY_TABU_RADIUS) { break; }
			} else if(lnL_score[i].TabuDist > 0 && lnL_score[i].TabuDist == MaxTabuDist && lnL_score[i].Lik +  FULL_LIK_ACC > BestLnL) {
				// If the best likelihood then store only it
				if(lnL_score[i].Lik > BestLnL) { BestLnL = lnL_score[i].Lik; StarFix.clear(); StarFix.push_back(i); }
				// Else add the tree to StarFix
				if(fabs(BestLnL-lnL_score[i].Lik) < FULL_LIK_ACC) { StarFix.push_back(i); }
		}	}
		if(MaxTabuDist == 0) { return tabu; } // Skip tabu trees
		if(StarFix.size() > 1)	{ random_shuffle(StarFix.begin(),StarFix.end()); }
		BestLnL = lnL_score[StarFix[0]].Lik; BestTree = lnL_score[StarFix[0]].TreeNum;
//		cout << "\n-----------------------------------------------------------------\nTrees:"; FOR(i,lnL_score.size()) { cout << "\n\tTree["<<i<<"]: " << TL->m_Trees[lnL_score[i].TreeNum] << " = " << lnL_score[i].Lik << "("<< lnL_score[i].TabuDist << ")"; if(IsIn(i,StarFix)) { cout << " *In StarFix*"; } if(i==StarFix[0]) { cout << " !BEST!"; } } cout << "\nBest: " << TL->m_Trees[BestTree] << " = " << BestLnL << " cf. Original: " << OldTree << " = " << *OrilnL;
	} else { // Still have to avoid tabu trees
		MaxTabuDist = -1;
		FOR(i,(int)lnL_score.size()) {
			if(lnL_score[i].TreeNum == OriTreeNum) { if(HaveTabu) { return tabu; } else { return none; } }
			// Form the tree
			TestTree = *Model->m_pTree;
			TestTree.ReplaceTreeCP(&TL->m_Trees[lnL_score[i].TreeNum],Model->LeafMap(),Model->NodesCovered(),false);
			// Check whether tabu
			if(IsTabu(&TestTree,true,false,Model)) { continue; }
			BestLnL = lnL_score[i].Lik; BestTree = lnL_score[i].TreeNum; break;
		}
//		cout << "\n-----------------------------------------------------------------\nTrees:"; FOR(i,lnL_score.size()) { cout << "\n\tTree["<<i<<"]: " << TL->m_Trees[lnL_score[i].TreeNum] << " = " << lnL_score[i].Lik ; } cout << "\nBest: " << TL->m_Trees[BestTree] << " = " << BestLnL << " cf. Original: " << OldTree << " = " << *OrilnL;
	}
	// Store the best tree
	/////////////////////////////////////////////////////
	TreeArrange NewTA;
	if(!InRange(BestTree,0,(int)TL->m_Trees.size())) { Error("Error: trying to assign best tree too large"); }
	NewTA.lnL = BestLnL; NewTA.NodeNum = NodeNum;
	NewTA.Tree = TL->m_Trees[BestTree]; NewTA.TabuDist = MaxTabuDist;
	NewTA.LeafList = Model->LeafMap(); NewTA.NodCovList = Model->NodesCovered();
	if(TA->empty()) { TA->push_back(NewTA); }
	else	{
		pIFOR(TA_i,TA)	{
			// Do sorting by likelihood when not worrying about TABU dist
			if(ENCOURAGE_SNAP_AWAY == 0 || OptObs < EXIT_OBS - COMP_SNAP_AWAY_2_OPT_OBS || MaxTabuDist> SNAP_AWAY_TABU_RADIUS) { if(BestLnL > TA_i->lnL) { break; } }
			else { if(MaxTabuDist > TA_i->TabuDist || (MaxTabuDist == TA_i->TabuDist && BestLnL > TA_i->lnL) ) { break; } }
		}
		TA->insert(TA_i,NewTA);
	}
	if(StarFix.size() > 1) { return star; } return normal;
}

vector <STreeScore> DoSubTreeCalcs(double CurlnL, vector <STreeScore> *OriList, CTreeList *TL, CBaseModel *M,int SetSize, bool DoFullOpt, LazyType DoLazy, double PropDiff, bool ForceSetSize)	{
	bool FoundBetter = false;
	int i, BestTreeCount  = 0;
	double lnL2Beat = CurlnL, lnL, ThresholdScore;
	STreeScore Score; Score.TabuDist = -1;
	vector <STreeScore> NewList;
	if(PropDiff > 1.0) { ThresholdScore = OriList->at(0).Lik * PropDiff; } else { ThresholdScore = -BIG_NUMBER; }
	FOR(i,(int)OriList->size())	{
		if(ForceSetSize && i == SetSize) { break; }
		assert(OriList->at(i).TreeNum < (int)TL->m_Trees.size());
		M->ApplySubTree(&TL->m_Trees[OriList->at(i).TreeNum]);		// Apply the new subtree
		// Get some likelihoods
		DoLazy = fullopt; DoFullOpt = true;
		if(DoFullOpt)	{
			if(DoLazy == fullopt) { lnL = FullOpt(M,false,true,false,-BIG_NUMBER,true,DEFAULT_OPTNUM,lnL2Beat); }
			else		{ lnL = LazyOpt(M,false,true,false,-BIG_NUMBER,true,DEFAULT_OPTNUM,lnL2Beat); }
		} else { lnL = M->lnL(); }
		// If its is a good value put into the list of trees to further examine
		if(lnL + FULL_LIK_ACC > lnL2Beat) { FoundBetter = true; lnL2Beat = lnL; BestTreeCount = 0;} BestTreeCount++;
		Score.Lik = lnL; Score.TreeNum = OriList->at(i).TreeNum;
		// Store the likelihood
		NewList.push_back(Score);
		// Some new stuff for fast rearrangements
		if(M->NoSeq() > 6) {
			if(DO_QUICK_SNAP == 1 && DoFullOpt && lnL - CurlnL > 0.1)  { // Only need one good tree in full opt
				SortSTreeScore(&NewList);
				return NewList;
			}
			if(BestTreeCount >= SetSize && OriList->at(i).Lik > ThresholdScore) { break; }
	}	}
	SortSTreeScore(&NewList);
	return NewList;
}

/* **************************************************************************** */

//		ROUTINES TO DO SPR REARRANGEMENTS

/* **************************************************************************** */

//////////////////////////////////////////////////////////////////////
// SPR routine -- fast and uses routines for adding sequences
// -----------
// Variables:
//  - M			= model used for calculations
//  - DoOutput	= whether any output is thrown to screen during calculations
//  - OrilnL	= Original likelihood (if -BIG_NUMBER then recalculated)
//  - TryHard	= if(false) = Only top SA_OPT_NUM_ACC potential trees are triplet optimised
//				= if(true)  = SA_OPT_NUM_ACC and the top tree for every subtree are triplet optimised
//  - DoAll		= if(false) = only subtrees rearrangements are calculated
//				= if(true)  = subtrees and leaf rearrangements are calculated
//  - MaxRoundsSPR = Maximum number of SPR rounds on allowed on the tree
// ---
// TODO: There are instances where the same rearrangement is calculated multiple times. This should be avoided
//////////////////////////////////////////////////////////////////////
// This is the driver routine
double TreeSPR(CBaseModel *M, bool DoOutput, double OrilnL,bool TryHard, bool DoAll, int MaxRoundsSPR, LazyType DoLazy)	{
	int RoundNum, NoImp = 0;
	double BestlnL = OrilnL,OldlnL;
	CTree T; T = *M->m_pTree;
	if(fabs(OrilnL+BIG_NUMBER) < DBL_EPSILON) { OrilnL = FullOpt(M); } OldlnL = BestlnL = OrilnL;
	// Do some output
	if(DoOutput) { cout << "\n\n~ ~ ~\n\n\tPerforming SPR: " << BestlnL << endl; }
	FOR(RoundNum,MaxRoundsSPR)	{
		if(DoOutput) { cout << "\nRound "<< RoundNum << "/" << MaxRoundsSPR << "; lnL: " << BestlnL; }
		if(M->IsLikelihoodCalc() || M->IsRMSDCalc()) { BestlnL = DoScoreSPR(M,DoOutput,BestlnL,TryHard,DoAll,DoLazy); }
		else { BestlnL = DoSimpleSPR(M,DoOutput,BestlnL,TryHard,DoAll,DoLazy); }
		// Manage return statements
		assert(BestlnL < 0);
		if(IsSameTree(&T,M->m_pTree)) {
			if(DoOutput) { cout << " best"; }
			break;  // No change -- optimum found
		}
		if(DoOutput) { cout << " + " << BestlnL - OldlnL; }
		if(OldlnL > BestlnL) { if(DoOutput) { cout << " worse"; } BestlnL = OldlnL; *M->m_pTree = T; break; }
		OldlnL = BestlnL;
		T = *M->m_pTree;
		NoImp = 0;
	}
	return BestlnL;
}

// This is the SPR subroutine that does the work when trees have no branch by branch score
// Currently this only includes parsimony
double DoSimpleSPR(CBaseModel *M, bool DoOutput, double OrilnL,bool TryHard, bool DoAll, LazyType Lazy)	{
	int i,Best;
	double BestlnL, lnL;
	bool OriMLCalc = M->IsLikelihoodCalc();
	CTree OT; OT = *M->Tree();
	vector <CTree *> TL;
	if(M->IsRMSDCalc()) { Error("\nCannot do DoSimpleSPR when base method is RMSD..."); }
	M->DoParsimony();
	BestlnL = M->lnL();
	M->Tree()->GetSPRList(&TL,1,5);
	Best = -1; FOR(i,(int)TL.size()) {
		*M->Tree() = *TL[i];
		lnL = M->lnL();
		if(lnL > BestlnL) { Best = i; BestlnL = lnL; }
	}
	if(Best != -1) { *M->Tree() = *TL[Best]; } else { *M->Tree() = OT; }
	FOR(i,(int)TL.size()) { delete TL[i]; }
	return BestlnL;
}
// This is the SPR subroutine that does the work for score based routines
// Currently this includes ML and RMSD
#define DO_BRA_STABILITY 0

double DoScoreSPR(CBaseModel *M, bool DoOutput, double OrilnL,bool TryHard, bool DoAll, LazyType DoLazy)	{
	int i,j,Br,BrRef, Link,Seq2Add,OriBr,iOptlnL,count = 0;
	bool flag = false;
	double BestlnL = OrilnL, lnL, temp;
	vector <int> BranchList;
	vector <double> Stability;
	vector <STreelnL*> List;
	STreelnL *TempTreelnL = NULL;
	vector <int> RepeatList;
	double dOptlnL;
	CTree T;
	assert(!M->IsSubTree());
	T = *M->m_pTree;

	if(fabs(OrilnL+BIG_NUMBER) < DBL_EPSILON) { FullOpt(M); }
	M->CleanFastCalc(true);
	if(M->IsRMSDCalc()) { Error("\nThe SPR routine doesn't yet work for RMSD calculations..."); }

	// Set up order branches will be visited by SPR using a branch stability measure
	FOR(i,T.NoBra()) { BranchList.push_back(i); }
#if DO_BRA_STABILITY == 0
	random_shuffle(BranchList.begin(),BranchList.end());
#else
	Stability = T.MeasureBraStability(PWDists);
	assert(Stability.size() == BranchList.size());
	FOR(i,(int)Stability.size()-1)	{
		for(j=(int)Stability.size()-1;i<j;j--)	{
			if(Stability[j-1]<Stability[j]) {
				temp = Stability[j-1]; Stability[j-1] = Stability[j]; Stability[j] = temp;
				Br = BranchList[j-1]; BranchList[j-1] = BranchList[j]; BranchList[j] = Br;
	}	}	}
#endif

	// The SPR routines
	///////////////////////////////////////////////////
	// 1. Get the rough likelihoods
	FOR(BrRef,T.NoBra())	{			// Loop through all branches in the tree
		Br = BranchList[BrRef];
		FOR(Link,2)	{				// Loop though links in the branch
			if(T.BraLink(Br,Link) >= T.NoSeq()) {	// If its a subtree then perform SPR
//				cout << "\nPrepareSPR(" << Br << ")";
				Seq2Add = M->PrepareSPR(Br,Link,&OriBr,DoLazy);
				if(Seq2Add == -1)	{ *M->m_pTree = T; continue; }
				BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,OriBr,max(DoLazy,lazy),&T);

			} else if(DoAll == true)	{			// If its a single sequence then do addition if DoAll == true
				Seq2Add = M->m_pTree->BraLink(Br,Link);
				M->m_pTree->RemoveLeafNode(Seq2Add);
				M->m_pTree->GetStart();
				BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,-1,max(DoLazy,lazy),&T);
			}
			// Restore the old tree
			*M->m_pTree = T;
			// If a tree likelihood is better and doing lazy optimisation, then accept it
			FOR(i,List.size()  && DoLazy == lazy) {
				if(List[i]->lnL > BestlnL) {
					// Check the better tree isn't tabu
					if(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2) >= T.NoSeq()) {	// Do Subtree rearrangement
						M->m_pTree->PerformSPR(List[i]->SPRID/2,List[i]->SPRID%2,List[i]->BraNum,List[i]->Tree,List[i]->BraLinks);
					} else {	// Do sequence addition
						M->m_pTree->RemoveLeafNode(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2));
						M->m_pTree->AddSeq(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2),List[i]->BraNum);
					}
					if(IsTabu(M->m_pTree,true,false,M)) { *M->m_pTree = T; flag = false; }
					else { // Wipe the list clean and store only the best tree
						// Note the tree is now changed here, so no rewrite of tree can be allowed!
						TempTreelnL = new STreelnL; *TempTreelnL = GetBestList(&List);
						CleanList(&List); List.push_back(TempTreelnL);
						TempTreelnL = NULL;
						flag = true; break;
					}
				}
				if(flag == true) { break; }
			}
			// Do Output
#if FUNC_COUNTERS == 1
			SPR_Log_Counter++;
#endif
			if(DoOutput) { if( (count) % 80 == 0) { cout << "\n\t"; } if((count++)%10 == 0) { cout << " "; } cout << "." << flush; }
			if(flag == true) { break; }

		}
		if(flag == true) { break; }
	}
	// 2. Sort the list
	SortList(&List);;
//	cout << "\nSortedTrees[size:"<<List.size()<<"]: CurBest: " << OrilnL;
//	M->FastBranchDetails(); exit(-1);
//	FOR(i,min(50,List.size())) { cout << "\nTree["<<i<<"]: " << List[i]->Tree << " = " << List[i]->lnL << " SPRID= " << List[i]->SPRID << flush; }
	// Only do the remaining steps if it's actually needed...
	if(flag == false) {
		// 3. Organise to fully optimise the top SA_OPT_NUM_ACC proposed trees
		FOR(i,SA_OPT_NUM_ACC)	{ List[i]->DoOpt = true; RepeatList.push_back(List[i]->SPRID); }
		// 4. If TryHard==true then also optimise the best choice for each possible subtree arrangement
		if(TryHard)	{
			for(Br=M->m_pTree->NoSeq();Br<M->m_pTree->NoBra();Br++)	{	// Loop through branches in tree
				if(M->m_pTree->BraLink(Br,1) < M->m_pTree->NoSeq()) { continue; }
				FOR(Link,2)	{				// Loop though links in the branch
					dOptlnL = -BIG_NUMBER; iOptlnL = -1;
					FOR(i,(int)List.size()) { if(List[i]->SPRID == (Br*2)+Link && List[i]->lnL > dOptlnL) { dOptlnL = List[i]->lnL; iOptlnL = i; } }
					if(iOptlnL == -1) { continue; }
					List[iOptlnL]->DoOpt = true; RepeatList.push_back(List[iOptlnL]->SPRID);
				}	}	}
		if(DoOutput) { count = 0; }
		// 5. Do the full optimisations. Note only those in RepeatList are done
		FOR(Br,T.NoBra())	{			// Loop through all branches in the tree
			FOR(Link,2)	{				// Loop though links in the branch
				if(!IsIn((Br*2)+Link,RepeatList)) { continue; }
				if(T.BraLink(Br,Link) >= T.NoSeq()) {	// If its a subtree then perform SPR
					Seq2Add = M->PrepareSPR(Br,Link,&OriBr,DoLazy);
					if(Seq2Add == -1)	{ *M->m_pTree = T; continue; }
					BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,OriBr,DoLazy,&T);
				} else if(DoAll == true)	{			// If its a single sequence then do addition if DoAll == true
					Seq2Add = M->m_pTree->BraLink(Br,Link);
					M->m_pTree->RemoveLeafNode(Seq2Add);
					M->m_pTree->GetStart();
					BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,-1,DoLazy,&T);
				}
				// Do Output
				if(DoOutput) { if( (count) % 80 == 0) { cout << "\n\t"; } if((count++)%10 == 0) { cout << " "; } cout << "*" << flush; }
				// Restore the old tree
				*M->m_pTree = T;
			}	}
		// 6. Sort again
		SortList(&List);

#if SPR_DEBUG == 1
		cout << "\nBestTrees[size:"<<List.size()<<"]: CurBest: " << OrilnL;
		FOR(i,min(50,List.size())) { cout << "\nTree["<<i<<"]: " << List[i]->Tree << " = " << List[i]->lnL << " SPRID= " << List[i]->SPRID << flush; }
#endif
		// 5. Try the best tree that isn't tabu and see whether any improvement is found
		FOR(i,min(SA_OPT_NUM_ACC,(int)List.size())) {
			if(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2) >= T.NoSeq()) {	// Do Subtree rearrangement
				M->m_pTree->PerformSPR(List[i]->SPRID/2,List[i]->SPRID%2,List[i]->BraNum,List[i]->Tree,List[i]->BraLinks);
			} else {	// Do sequence addition
				M->m_pTree->RemoveLeafNode(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2));
				M->m_pTree->AddSeq(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2),List[i]->BraNum);
			}
			lnL = M->FastBranchOpt(M->lnL(true),0.01,NULL,5,false);
#if SPR_DEBUG == 1
			cout << "\nTrying tree["<<i<<"]: " << *M->m_pTree << " = Before: " << List[i]->lnL << ", branch_opt: " << lnL << " cf. full: " << M->lnL();
#endif
			if(IsTabu(M->m_pTree,true,false,M)) { *M->m_pTree = T; continue; }
			else if(lnL > BestlnL) { break; }
			else { *M->m_pTree = T; }
		}
	}
	// 6. Optimise the best SPR tree
	if(DoLazy == fullopt)	{ BestlnL = FullOpt(M,true,true,false,List[0]->lnL,true,DEFAULT_OPTNUM,OrilnL); }
	else					{ BestlnL = LazyOpt(M,true,true,false,List[0]->lnL,true,DEFAULT_OPTNUM,OrilnL); }
	if(DoOutput) { for(i = count %80; i< 80;i++) { cout << " "; } cout << " : " << BestlnL << flush; }
	// 7. Manage the return
	if(IsTabu(M->m_pTree,true,false,M,BestlnL)) { cout << " tabu" << flush; *M->m_pTree = T; }
	else if(BestlnL < OrilnL) { *M->m_pTree = T; }
	assert(BestlnL < 0);
	return BestlnL;
}

// Ori ScoreSPR function
/*
double DoScoreSPR(CBaseModel *M, bool DoOutput, double OrilnL,bool TryHard, bool DoAll, LazyType DoLazy)	{
	int i,j,Br,BrRef, Link,Seq2Add,OriBr,iOptlnL,count = 0;
	bool flag = false;
	double BestlnL = OrilnL, lnL, temp;
	vector <int> BranchList;
	vector <double> Stability;
	vector <STreelnL*> List;
	STreelnL *TempTreelnL = NULL;
	vector <int> RepeatList;
	double dOptlnL;
	CTree T;
	assert(!M->IsSubTree());
	T = *M->m_pTree;

	if(fabs(OrilnL+BIG_NUMBER) < DBL_EPSILON) { FullOpt(M); }
	M->CleanFastCalc(true);
	if(M->IsRMSDCalc()) { Error("\nThe SPR routine doesn't yet work for RMSD calculations..."); }

	// Set up order branches will be visited by SPR using a branch stability measure
	FOR(i,T.NoBra()) { BranchList.push_back(i); }
#if DO_BRA_STABILITY == 0
	random_shuffle(BranchList.begin(),BranchList.end());
#else
	Stability = T.MeasureBraStability(PWDists);
	assert(Stability.size() == BranchList.size());
	FOR(i,(int)Stability.size()-1)	{
		for(j=(int)Stability.size()-1;i<j;j--)	{
			if(Stability[j-1]<Stability[j]) {
				temp = Stability[j-1]; Stability[j-1] = Stability[j]; Stability[j] = temp;
				Br = BranchList[j-1]; BranchList[j-1] = BranchList[j]; BranchList[j] = Br;
	}	}	}
#endif

	// The SPR routines
	///////////////////////////////////////////////////
	// 1. Get the rough likelihoods
	FOR(BrRef,T.NoBra())	{			// Loop through all branches in the tree
		Br = BranchList[BrRef];
		FOR(Link,2)	{				// Loop though links in the branch
			if(T.BraLink(Br,Link) >= T.NoSeq()) {	// If its a subtree then perform SPR
//				cout << "\nPrepareSPR(" << Br << ")";
				Seq2Add = M->PrepareSPR(Br,Link,&OriBr,DoLazy);
				if(Seq2Add == -1)	{ *M->m_pTree = T; continue; }
				BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,OriBr,max(DoLazy,lazy),&T);

			} else if(DoAll == true)	{			// If its a single sequence then do addition if DoAll == true
				Seq2Add = M->m_pTree->BraLink(Br,Link);
				M->m_pTree->RemoveLeafNode(Seq2Add);
				M->m_pTree->GetStart();
				BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,-1,max(DoLazy,lazy),&T);
			}
			// Restore the old tree
			*M->m_pTree = T;
			// If a tree likelihood is better and doing lazy optimisation, then accept it
			FOR(i,List.size()  && DoLazy == lazy) {
				if(List[i]->lnL > BestlnL) {
					// Check the better tree isn't tabu
					if(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2) >= T.NoSeq()) {	// Do Subtree rearrangement
						M->m_pTree->PerformSPR(List[i]->SPRID/2,List[i]->SPRID%2,List[i]->BraNum,List[i]->Tree,List[i]->BraLinks);
					} else {	// Do sequence addition
						M->m_pTree->RemoveLeafNode(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2));
						M->m_pTree->AddSeq(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2),List[i]->BraNum);
					}
					if(IsTabu(M->m_pTree,true,false,M)) { *M->m_pTree = T; flag = false; }
					else { // Wipe the list clean and store only the best tree
						TempTreelnL = new STreelnL; *TempTreelnL = GetBestList(&List);
						CleanList(&List); List.push_back(TempTreelnL);
						TempTreelnL = NULL;
						flag = true; break;
			} 	}	}
			// Do Output
#if FUNC_COUNTERS == 1
			SPR_Log_Counter++;
#endif
			if(DoOutput) { if( (count) % 80 == 0) { cout << "\n\t"; } if((count++)%10 == 0) { cout << " "; } cout << "." << flush; }
		}
		if(flag == true) { break; }
	}
	// 2. Sort the list
	SortList(&List);;
	cout << "\nSortedTrees[size:"<<List.size()<<"]: CurBest: " << OrilnL;
//	M->FastBranchDetails(); exit(-1);
	FOR(i,min(50,List.size())) { cout << "\nTree["<<i<<"]: " << List[i]->Tree << " = " << List[i]->lnL << " SPRID= " << List[i]->SPRID << flush; }
	// Only do the remaining steps if it's actually needed...
	if(flag == false) {
		// 3. Organise to fully optimise the top SA_OPT_NUM_ACC proposed trees
		FOR(i,SA_OPT_NUM_ACC)	{ List[i]->DoOpt = true; RepeatList.push_back(List[i]->SPRID); }
		// 4. If TryHard==true then also optimise the best choice for each possible subtree arrangement
		if(TryHard)	{
			for(Br=M->m_pTree->NoSeq();Br<M->m_pTree->NoBra();Br++)	{	// Loop through branches in tree
				if(M->m_pTree->BraLink(Br,1) < M->m_pTree->NoSeq()) { continue; }
				FOR(Link,2)	{				// Loop though links in the branch
					dOptlnL = -BIG_NUMBER; iOptlnL = -1;
					FOR(i,(int)List.size()) { if(List[i]->SPRID == (Br*2)+Link && List[i]->lnL > dOptlnL) { dOptlnL = List[i]->lnL; iOptlnL = i; } }
					if(iOptlnL == -1) { continue; }
					List[iOptlnL]->DoOpt = true; RepeatList.push_back(List[iOptlnL]->SPRID);
				}	}	}
		if(DoOutput) { count = 0; }
		// 5. Do the full optimisations. Note only those in RepeatList are done
		FOR(Br,T.NoBra())	{			// Loop through all branches in the tree
			FOR(Link,2)	{				// Loop though links in the branch
				if(!IsIn((Br*2)+Link,RepeatList)) { continue; }
				if(T.BraLink(Br,Link) >= T.NoSeq()) {	// If its a subtree then perform SPR
					Seq2Add = M->PrepareSPR(Br,Link,&OriBr,DoLazy);
					if(Seq2Add == -1)	{ *M->m_pTree = T; continue; }
					BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,OriBr,DoLazy,&T);
				} else if(DoAll == true)	{			// If its a single sequence then do addition if DoAll == true
					Seq2Add = M->m_pTree->BraLink(Br,Link);
					M->m_pTree->RemoveLeafNode(Seq2Add);
					M->m_pTree->GetStart();
					BraSACalc(Seq2Add,&List,M,true,&BestlnL,(Br*2)+Link,-1,DoLazy,&T);
				}
				// Do Output
				if(DoOutput) { if( (count) % 80 == 0) { cout << "\n\t"; } if((count++)%10 == 0) { cout << " "; } cout << "*" << flush; }
				// Restore the old tree
				*M->m_pTree = T;
			}	}
		// 6. Sort again
		SortList(&List);

#if SPR_DEBUG == 1
		cout << "\nBestTrees[size:"<<List.size()<<"]: CurBest: " << OrilnL;
		FOR(i,min(50,List.size())) { cout << "\nTree["<<i<<"]: " << List[i]->Tree << " = " << List[i]->lnL << " SPRID= " << List[i]->SPRID << flush; }
#endif
		// 5. Try the best tree that isn't tabu and see whether any improvement is found
		FOR(i,min(SA_OPT_NUM_ACC,(int)List.size())) {
			*M->m_pTree = T;
			if(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2) >= T.NoSeq()) {	// Do Subtree rearrangement
				M->m_pTree->PerformSPR(List[i]->SPRID/2,List[i]->SPRID%2,List[i]->BraNum,List[i]->Tree,List[i]->BraLinks);
			} else {	// Do sequence addition
				M->m_pTree->RemoveLeafNode(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2));
				M->m_pTree->AddSeq(T.BraLink(List[i]->SPRID/2,List[i]->SPRID%2),List[i]->BraNum);
			}
			lnL = M->FastBranchOpt(M->lnL(true),0.01,NULL,5,false);
#if SPR_DEBUG == 1
			cout << "\nTrying tree["<<i<<"]: " << *M->m_pTree << " = Before: " << List[i]->lnL << ", branch_opt: " << lnL << " cf. full: " << M->lnL();
#endif
			if(IsTabu(M->m_pTree,true,false,M)) { *M->m_pTree = T; continue; }
			else if(lnL > BestlnL) { break; }
			else { *M->m_pTree = T; }
	}	} else {	// Otherwise apply the good tree found and optimise its branch lengths a bit


	}
	// 6. Optimise the best SPR tree
	if(DoLazy == fullopt)	{ BestlnL = FullOpt(M,true,true,false,List[0]->lnL,true,DEFAULT_OPTNUM,OrilnL); }
	else					{ BestlnL = LazyOpt(M,true,true,false,List[0]->lnL,true,DEFAULT_OPTNUM,OrilnL); }
	if(DoOutput) { for(i = count %80; i< 80;i++) { cout << " "; } cout << " : " << BestlnL << flush; }
	// 7. Manage the return
	if(IsTabu(M->m_pTree,true,false,M,BestlnL)) { cout << " tabu" << flush; *M->m_pTree = T; }
	else if(BestlnL < OrilnL) { *M->m_pTree = T; }
	assert(BestlnL < 0);
	return BestlnL;
} */

#define GS_DELTA 5
double GoldenSection(double OrilnL, double *x, CPar *Par,CBaseModel *M)	{
	int i;
	double dx = 1.0E-3;
	double p_old = *x, lnL_old = -fabs(OrilnL), x1,x2=*x,x3,x1_lnL = 1.0,x2_lnL = -fabs(OrilnL),x3_lnL = 1.0, xi,temp;
//	cout << "\n\n---\nDoing GoldenSection"; Par->m_bOutDetail = true;
//	cout << "\n--> " << *Par;
//	cout << "\nPar = " << Par->Name() << "; x2: " << x2 << " == " << x2_lnL << " cf. " << M->lnL() << " Bounds("<<Par->LowBound()<<","<<Par->UpBound()<<")";
//	cout << "\nPar: " << Par->Val() << " cf " << *x;
	// Get left bracketing
//	cout << "\nGoing left...";
	cout << "\nGolden section bounding could/should be improved";
	while(x2_lnL < x1_lnL)	{
//		cout << "-" << flush;
		*x = x1 = x2 - max((fabs(x2) * dx),dx); Par->Val();
		if(!Par->CheckLowBound()) { x1 = Par->Val(); x1_lnL = M->lnL();  break; }
		x1_lnL = M->lnL();
		dx *= GS_DELTA;
//		cout << "\nleft: " << x1 << " == " << x1_lnL;
		if(x1_lnL > x2_lnL)	{
			x3 = x2; x3_lnL = x2_lnL;
			x2 = x1; x2_lnL = x1_lnL;
			x1_lnL = 1.0;
		}
	}
	if(x3_lnL > 0.0)	{
		dx = 1.0E-4;
		// Get right bracketing
		while(x2_lnL < x3_lnL)	{
//			cout <<"+" << flush;
			*x = x3 = x2 + max((fabs(x2) * dx),dx); Par->Val();
			if(!Par->CheckUpBound()) { x3 = Par->Val(); x3 = M->lnL(); break; }
			x3_lnL = M->lnL();
//			cout << "\nRight: " << x3 << " == " << x3_lnL;
			if(x3_lnL > x2_lnL)	{
				x1 = x2; x1_lnL = x2_lnL;
				x2 = x3; x2_lnL = x3_lnL;
				x3_lnL = 1.0;
			}
			dx *= GS_DELTA;
	}	}
/*
	cout << "\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
	*x = x2; cout << "\nx2: " << x2 << " == " << x2_lnL << " cf. " << M->lnL();
	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
	cout << "\nPar: " << Par->Val();
	cout << "\n--> " << *Par;
 	cout << "\n==> " << M->lnL();
	if(fabs(x2_lnL - M->lnL()) > 0.01)  { exit(-1); } */
//	exit(-1);
	// Do various checks to make sure it looks like a hill
	if(!Par->CheckBound()) { *x = p_old; return lnL_old; }
	if(x1_lnL > max(x2_lnL,x3_lnL)) {
		// Deal with the U shape case
		if(x3_lnL > x2_lnL) {
			if(x1_lnL > x3_lnL) { *x = x1; GoldenSection(x1_lnL,x,Par,M); }
			else				{ *x = x3; GoldenSection(x3_lnL,x,Par,M); }
		} else { *x = x1; return x1_lnL; }
	}
	else if (x3_lnL > x2_lnL)		{ *x = x3; return x3_lnL; }
	// Now do golden section search to find optimal value
	assert(x2_lnL + FLT_EPSILON > x1_lnL && x2_lnL + FLT_EPSILON> x3_lnL);
	FOR(i,20)	{
//		cout << "." << flush;
		if(fabs(x3 - x2) > fabs(x2 - x1))	{ *x = xi = x2 + ( ( x3 - x2) * GOLDEN_NUMBER); }
		else								{ *x = xi = x2 - ( ( x2 - x1) * GOLDEN_NUMBER); }
		temp = M->lnL();
		if(temp > x2_lnL) {
			if(fabs(x3 - x2) > fabs(x2 - x1))	{ x1 = x2; x1_lnL = x2_lnL; x2 = xi; x2_lnL = temp; }
			else								{ x3 = x2; x3_lnL = x2_lnL; x2 = xi; x2_lnL = temp; }
		} else		  {
			if(fabs(x3 - x2) > fabs(x2 - x1))	{ x3 = xi; x3_lnL = temp; }
			else								{ x1 = xi; x1_lnL = temp; }
		}
//		cout << "\n\t\tLeft: " << x2_lnL - x1_lnL <<"; Right: " << x2_lnL - x3_lnL;
//		cout << "\n\t["<<i<<"] ---\n\t\tx1: " << x1 << " = " << x1_lnL;
//		cout << "\n\t\tx2: " << x2 << " = " << x2_lnL;
//		*x = x2; cout << " cf " << M->lnL();
//		cout << "\n\t\tx3: " << x3 << " = " << x3_lnL;
		if(fabs(x2_lnL - x3_lnL) < FULL_LIK_ACC && fabs(x2_lnL - x1_lnL) < FULL_LIK_ACC && i > 2)	{ break; }
	}
//	cout << "\n---\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
//	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
//	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
	if(lnL_old > x2_lnL) { *x = p_old; return lnL_old; }
//	cout << " done: i=" << i << ", left= " << x2_lnL - x1_lnL << "; right= " << x2_lnL - x3_lnL << flush; // exit(-1);
//	if(Par->Name().find("Inv") != string::npos) { exit(-1); }
	*x = x2;
	return x2_lnL;
}

// Line searching routine
/// From NR in C Section 9.7 pp 385
/////////////////////////////////////////
// Arguments passed to function:
// -----------------------------
// n = number of dimensions
// xold[n]	= n-dimensional starting point
// fold 	= function value at xold
// g[n] 	= gradients at xold
// p[n] 	= the resulting direction of movement in space
// x[n] 	= the new point along direction p from xold
// f 		= the new function value
// stpmax	= limits the length of steps
// check	= Rather than return a value (0 = OK:1 = BAD)
//////////////////////////////////////////

#define INIT_ALAM 1.0

double lnsrch(vector <double *> x,double fold,vector <double> g, double p[], double pold[], double *f, bool Do_GS,CBaseModel *Model)	{
   	int i, j, n = (int)x.size(), num_run = 10;
    double a, alam = INIT_ALAM, alam2 = INIT_ALAM, alamin, b, disc, f2 = -BIG_NUMBER,fold2, rhs1, rhs2, slope, temp,
    	test, tmplam;
	double a0 = 0,a1 = INIT_ALAM ,a2 = -1,	v0 = *f, v1 = *f, v2 = -1, best_a, best_f, ori_f = *f;
	double BestAlam = -BIG_NUMBER;
	double BestlnL = BIG_NUMBER;
	// Get old parameters and impose a maximum step
	FOR(i,n) { pold[i] = *x[i]; }
    for(slope = 0.0,i = 0;i < n;i++)	{ slope += g[i] * p[i]; }
    test = 0.0;
    for(i = 0; i < n;i++)	{
	    temp = fabs(p[i]) / FMAX(fabs(pold[i]),1.0);
		if(temp > test) { test = temp; }
    }
#if DEBUG_MULD_OPT > 1
	cout << "\nEntering lnsrch (n="<<n<<"): lnL = " << fold << " cf. " << Model->lnL();
#endif
	if(test < 0.99) { alamin = FLT_EPSILON/test; } else { alamin = 100 * FLT_EPSILON; }



/*
	double checker;
	cout << "\nEntering lnsrch (n="<<n<<"): lnL = " << fold << " cf. " << Model->lnL();
	cout << "\nOriginal: ";
	cout << "\nOrips:    "; FOR(i,n) { cout << "\t" << pold[i]; }
	checker = 0; cout << "\nSteps:    "; FOR(i,n) { cout << "\t" << alam*p[i]; checker += fabs(alam*p[i]); }
	cout << "\nTotal step: " << checker ;
	cout << "\nThere's clearly a problem here with the chosen step sizes... Even when close to the optima the step sizes are at least an order of magnitude larger than one should expect for maximum steps...";
	cout << "\nConsider implementing the logged parameters bit";
	/*
	cout << "\nDirection:"; FOR(i,n) {
		if(fabs(p[i]) < 3*DX) { cout << "\tdone"; continue; }
		if(p[i]>0) 	{ *x[i] = pold[i] + (3*DX); checker = Model->lnL(); }
		else 		{ *x[i] = pold[i] - (3*DX); checker = Model->lnL(); }
		cout << "\n\t" << Model->m_vpAllOptPar[i]->Name() << " ; grad = " << Model->m_vpAllOptPar[i]->grad() << " ; step = " << alam*p[i] << " ; fold = " << fold << " cf. new: " << checker;
		cout << "\t" << checker + fold;
		*x[i] = pold[i];
	}
*/
    // Start main loop
    for(;;)	{
/*    	// Hard debug code...
    	cout << "\nChecking each step (alam = " << alam << ")";
    	FOR(i,n) { *x[i] = pold[i]; }
    	cout << "\nReset: ori=" << fold << " cf. " << -Model->lnL(true) << " -- diff: " << fold + Model->lnL(true);
    	FOR(i,n) {
    		cout << "\n\tPar[" << i << "] " << *x[i] << " + " << (alam * p[i]) << " -> " << pold[i] + (alam * p[i]);
    		*x[i] = pold[i] + (alam * p[i]); cout << " == " << -Model->lnL(true) << " ; imp: " << Model->lnL(true) + fold;
    		*x[i] = pold[i];
    	}

    	cout << "\nOrips:"; FOR(i,n) { cout << "\t" << pold[i]; }
    	cout << "\nSteps:"; FOR(i,n) { cout << "\t" << alam*p[i]; }
    	FOR(i,n) {
    		double left, middle, right;
    		*x[i] = pold[i] - 0.01;
    		left = -Model->lnL();
    		*x[i] = pold[i];
    		middle = -Model->lnL();
    		*x[i] = pold[i] + 0.01;
  			right = -Model->lnL();
    		cout << "\nPar["<<i<<"] = " << pold[i] << " += " << p[i] <<": " << left << " (" << middle << ") " << right;
    		if(middle < left && middle < right) { cout << " +++"; }
    		*x[i] = pold[i];
    	}
*/		// Add the value
		FOR(i,n) { *x[i] = pold[i] + (alam * p[i]); }
		// Perform calculation
		*f = -Model->lnL(); v2 = v1; v1 = *f; a2 = a1; a1 = alam;
//		cout << "\n\t\talam: " << alam << ":  " << *f;
#if DEBUG_MULD_OPT > 1
		cout << "\n\t\tlnsrch: alam="<<alam << "; func=" << *f; if(x.size() == 1) { cout << " x= " << *x[0]; }
#endif
 	    // If reached convergence where no improvement in likelihood can be found. Either alam too small or likelihood improvement too small.
		if(alam < alamin || fabs(fold - *f) < FULL_LIK_ACC)	{
			FOR(i,n) { *x[i] = pold[i]; } *f = fold;
//			cout << "\nConverged alam = " << alam;
//				cout << "\nNew:"; FOR(i,n) { cout << "\t" << pold[i]; }
//			cout << " Converged -- Returning: " << fold << " == " << Model->lnL() << " *f: " << *f;
//			cout << "\nAnd:"; FOR(i,n) { cout << "\t" << pold[i]; }
//			cout << "\nDiff: " << fabs(Model->lnL() + fold);
//			if(fabs(Model->lnL() + fold) > 0.00001) { cout << "\n\nYUCK!";  exit(-1); }
			return -Model->lnL();
		// Only exit if there is an increase in likelihood...
		} else if(fold - *f > FULL_LIK_ACC) {
//			cout << "\nlnL increase -- fold: " << fold << " - " << *f << " = " << fold - *f << " cf. Best (lnL: " << BestlnL << "; alam: " << alam << ")";
//			cout << "\nImprovement over best: " << BestlnL - *f;
			if(*f < BestlnL) { BestlnL =  *f; BestAlam = alam; } // If going well continue with next step
			else {	// Revert to the best alam
				alam = BestAlam;
				FOR(i,n) { *x[i] = pold[i] + (alam * p[i]); }
				*f = BestlnL;
				// DEBUG CHECKING
/*				cout << "\nDebug check in alam return";
				if(fabs(*f - -Model->lnL()) > 1.0E-6) {
					double number = -Model->lnL();
					cout << "\nError in alam restoral..."; cout << "\n\t\tReturning from lnsrch: fp: " << *f << "; lnL: " << number << "; Imp: " << ori_f - *f;
					cout << "\nDiff = " << *f - number << " == fabs() " << fabs(*f-number);

					exit(-1);
				}
*/				return alam;
			}
//			double BestAlam = -BIG_NUMBER;
//			double BestlnL = -BIG_NUMBER;
//			cout << "\n\t\tReturning from lnsrch: fp: " << *f << "; lnL: " << Model->lnL() << "; Imp: " << ori_f - *f;
//			cout << " return 1: " << *f;
//			return alam;
		}
		// Otherwise adjust the alam
//		else	{
	        // Need catch to find when the search for a new value is getting silly
			if(fabs(alam - INIT_ALAM) < DBL_EPSILON) { tmplam = -slope / (2.0 * (*f - fold - slope)); if(tmplam > 0.95 * alam) { tmplam = 0.95 * alam; } }
			else if(fabs(*f - BIG_NUMBER) < DBL_EPSILON || fabs(f2 - BIG_NUMBER) < DBL_EPSILON) { tmplam = 0.5 * alam; }
			else	{
				rhs1 = *f - fold - alam * slope;
				rhs2 = f2 - fold2 - alam2 * slope;
				// Check the divides for zeros
				a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
				b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
				if(a == 0.0) { tmplam = -slope / (2.0 * b); }
				else	{
					disc = b * b - 3.0 * a * slope;
					if(disc < 0.0) { disc = 0; }
					tmplam = (-b + sqrt(disc)) / (3.0 * a);
				}
				if(tmplam > 0.5 * alam) { tmplam = 0.5 * alam; }
			}
//		}
	    alam2 = alam;
	    f2 = *f;
	    fold2 = fold;
	    alam = max(tmplam, 0.1 * alam);
	}
	cout << "\n FAILURE EXIT";
	return 0.0;
}

///////////////////////////////////////////////////////////////////////
// Break the optimisation into little bits to try and do it more quickly
// Works when there a couple of odd gradients.
double SubSetlnsrch(double Prob, vector <double *> x,double *step_xi, vector <double> g,double lnL,CBaseModel *M)	{
	int i, n = (int)x.size();
	vector <int> DoGrad; FOR(i,n) { DoGrad.push_back(i); } random_shuffle(DoGrad.begin(),DoGrad.end());
	vector <double *> sub_x; vector <double> sub_g(n,0);
	double temp, temp2, *sub_xi,func,*pold;
	GET_MEM(sub_xi,double,n); GET_MEM(pold,double,n);
	cout << "\nDoGrad: " << DoGrad;
	while(!DoGrad.empty())	{	// Run linesearch for various subsets of parameters
		temp = 0.0;temp2 = Random();
		FOR(i,n)	{
			sub_x.push_back(x[DoGrad[0]]); pold[i] = *x[DoGrad[i]];
			sub_xi[i] = step_xi[DoGrad[i]]; sub_g[i] = g[DoGrad[i]];
			DoGrad.erase(DoGrad.begin());
			temp += probBinomial(n,i,Prob); if((i > 1 && temp > temp2) || DoGrad.empty()) { break; }
		}
		cout << "\nDoing new linesearch: " << sub_x.size();
		cout << "\nx=  "; FOR(i,(int)sub_x.size()) { cout << " " << *sub_x[i]; }
		cout << "\nxi= ";  FOR(i,(int)sub_x.size()) { cout << " " << sub_xi[i]; }
		lnsrch(sub_x,lnL,sub_g,sub_xi,pold,&func,true,M); lnL = func;
		cout << "\nNew lnL: " << lnL;
		FOR(i,(int)sub_x.size()) { sub_x[i] = NULL; } sub_x.clear(); sub_g.clear();
	}
	cout << "\n>>> New improved lnL: " << lnL;
	DEL_MEM(sub_xi); DEL_MEM(pold);
	return lnL;
}

// dfpmin specific definitions
// From NR in C, pp 428-9
/////////////////////////////////////////////////////
// Arguments passed to function:
// -----------------------------
// p[n] 	= Starting points of the parameters
// n		= Number of parameters
// iter		= Number of iterations performed
// fret		= The minimum value of the function
// (*dfunc)	= The function to calculate gradients
// (*func)	= The function to minimise
// gtol		= The tolerence of the gradients
//////////////////////////////////////////////////////
// 19x06: Now added a variable lnL2Beat used by the function to decide to bail out if it's blatantly not going to beat it.
//		lnL2Beat == -BIG_NUMBER when nothing is to be beated (remember MulD_Optimise minimises the -lnL)

#if DEBUG_MULD_OPT > 0
#define FREEALL delete [] dg; delete[]  hdg; for(i=0;i<n;i++) { delete [] hessin[i];}; DEL_MEM(sub_xi); \
	delete [] hessin; delete [] pold; delete [] xi; DEL_MEM(grad_delta); DEL_MEM(oldxi);  \
	optout <<"\nDoing hard check: "; CheckAllPar(Model,-fp,x,1.0E-4,optout); optout.close();
#else
#define FREEALL DEL_MEM(dg); DEL_MEM(hdg); for(i=0;i<n;i++) { DEL_MEM(hessin[i]);}; DEL_MEM(sub_xi); \
	DEL_MEM(hessin); DEL_MEM(pold); DEL_MEM(xi); DEL_MEM(grad_delta); DEL_MEM(oldxi); \
	//cout << "\nDoing hard check: "; CheckAllPar(Model,-fp,x,1.0E-4);
	//cout << "\n\tThere were " << its << " iterations: " << -fp << flush; /*  << "\n\t  grads ="; FOR(i,n) { cout << g[i] << " "; } cout << "\n\t  Paras ="; FOR(i,n) { cout << *x[i] << " "; }*/ \
	//exit(-1);
#endif

double MulD_Optimise(double OrilnL,double gtol ,double ltol,vector <double *> x,CBaseModel *Model,int NI, bool DoBasicOutput,bool OnlyBranches, int OptTol, bool NewOne, double lnL2Beat, int NoBranchOpt,bool AllowOnlyParOpt,bool TryReallyHard)	{
	int i, its, j, NumberIter = NI, LikTol = 0, HessWarning = 0, n = (int)x.size();
    double den, fold, fac, fad, fae, fp, sum = 0.0, sumdg, sumxi, temp, test, max_g, last_improvement = BIG_NUMBER;
	double inc, temp_lnL;	// Some values describing the increases in likelihood
    double *dg, *hdg, **hessin, *pold, *xi, *oldxi, fret, alpha, *sub_xi;
    double step_max = BIG_STEP_MAX, *grad_delta, Last5[5] = {BIG_NUMBER,BIG_NUMBER,BIG_NUMBER,BIG_NUMBER,BIG_NUMBER};
	double PredictedlnL;
	vector <double> g, temp_g; g.assign(x.size(),0.0);
	vector <double *> temp_x;
	bool flag, ResetHess, Do_GS, GradOK, DoingSubset = false;
	if(n == 0)
	if(n <= 0) { return -BIG_NUMBER; }

	if(DoBasicOutput) { cout << "\n\tOptimising likelihood for model " << Model->Name() << "; this may take some time..." << flush; }

    // Allocate memory
	GET_MEM(dg,double,n);  GET_MEM(hdg,double,n); GET_MEM(sub_xi,double,n);
	GET_MEM(hessin,double*,n); FOR(i,n) { hessin[i] = NULL; GET_MEM(hessin[i],double,n); }
	GET_MEM(pold,double,n); GET_MEM(xi,double,n); GET_MEM(grad_delta,double,n);
	GET_MEM(oldxi,double,n);

//	cout << "\nDifference fp: " << fp << " and " << Model->lnL() << " cf. " << Model->lnL(true);

	// Calculate starting values
    fold = BIG_NUMBER;
	fret = fp = -OrilnL;
	// Do some initial branch optimising to get something parabola like
	if(Model->Locked() || (!Model->IsRMSDCalc() && n > 1))	{
		flag = false;
		fret = fp = -Model->FastBranchOpt(-fret,ltol,&flag,NoBranchOpt); // Do a round of fast branch optimisation
	}
	// If only doing branches and has converged
	if(Model->Locked() || (OnlyBranches == true && flag ==true)) {
		FREEALL;
		return -fp;
	}
	// Do only parameter optimisation if allowed. TODO: Should be removed?
	if(AllowOnlyParOpt || Model->ForceSeperateParOpt())	{ fret = fp = -DoOnlyParOpt(fp,gtol,ltol,x,Model,NI,OptTol,2); }
	// Set up hessian and initialise steps
	g = Model->GetDerivatives(-fp,&GradOK);
    FOR(i,n) {
		IMat(hessin,n);
		if(fabs(*x[i]) > 1 && fabs(g[i]) < 0.2) { xi[i] = -g[i] * fabs(*x[i]); } else { xi[i] = -g[i]; }
		sum += *x[i] * *x[i];
    }
	// Iterate
	FOR(its,NumberIter)	{
//		cout << "\n\t--- Iter: " << its << ": " << fp; //  << " cf. " << Model->lnL(true) << " (" << fabs(fp+Model->lnL(true)) << ")" << flush ;
		if(DoBasicOutput) {
//			cout << "\n\t["<<its<<"] " <<fold << " -> " << fp << flush; //  << " cf. " << Model->lnL() << flush;
			if(its == 0) { cout << "\n\t\t" << its << ": " << flush; }
			else if (its % 60 == 0) { cout << ": " << fp << "\n\t\t"<<its<<": "; }
			//if(fold < fp - FULL_LIK_ACC) {
			if(fp - fold > FULL_LIK_ACC * 100) {
				cout.precision(12);
				cout << "\ndiff ("<<fold << "-" << fp << "): " << fabs(fold-fp) << ";";
				Error(" ... Error... likelihood decreased in value???");
			}
			cout << "." << flush;
		}

		// Some initialisation
		ResetHess = false; Do_GS = false; step_max = BIG_STEP_MAX;
		// Decide whether to do Golden Section search in line search
		if(its < THOROUGH_LINE_SEARCH)  { Do_GS = true; }
		if(fabs(fold - fp) > 0.25) { step_max = SMALL_STEP_MAX; }
		fold = fp;

		// --------------------------------- Perform the line search ----------------------------------
		// Debug information
/*		cout << "\n<<<<<< ITER " << its << ": " << fp << " >>>>>>>>";
		cout << "\n\tPnames:"; FOR(i,n) {
			if(Model->m_vpAllOptPar[i]->Special()) { cout << "*SPECIAL*"; exit(-1); }
			cout << "\t" << Model->m_vpAllOptPar[i]->Name(); }
		cout << "\n\tPreals:"; FOR(i,n) { cout<< "\t" << Model->m_vpAllOptPar[i]->Val(); }
		cout << "\n\tP:     "; FOR(i,n) { cout << "\t" << *x[i]; }
		cout << "\n\tG:     "; FOR(i,n) { cout << "\t" << g[i]; }*/
		// Before performing line-search decide whether only a subset of parameters are worth examining
		// This is decided when there are very large gradients relative to all other gradients.
		double GOOD_IMPROVE = 0.5;
		double BIG_GRADIENT_MULTIPLIER = 0.33;
//		cout << "\nLast improvement: " << last_improvement;
		if((last_improvement > GOOD_IMPROVE || NumberIter - its < 3) && its < 10) { // If the last improvement was good then try looking at subsets of parameters
			max_g = 0.0; FOR(i,(int)g.size()) { max_g = max(fabs(g[i]),max_g); } max_g = max(max_g,10);
			FOR(i,(int)g.size()) { if(fabs(max_g * BIG_GRADIENT_MULTIPLIER) > fabs(g[i])) { g[i] = xi[i] = 0.0; DoingSubset = true; } } // Set small gradients to step size zero; Will need to reset the Hessian
		} else { // If we're close to an optima then gradient should agree with direction. Sign of a screwed up Hessian when it doesn't
			if(DoingSubset) { ResetHess = true; DoingSubset = false; }	// Finished with subsets, so reset the Hessian and get going again
			flag = false;
			FOR(i,(int)g.size()) {
				if(g[i] > 0 && xi[i] > 0) { xi[i] *= -1; flag = true; } // if(fabs(xi[i]) > 0.01) { xi[i] = -0.01;  } }
				if(g[i] < 0 && xi[i] < 0) { xi[i] *= -1; flag = true; } // if(fabs(xi[i]) > 0.01) { xi[i] = 0.01; } }
			}
			if(flag == true) { HessWarning ++; }
		}
		// Now adjust the step sizes to conform to step_max
		max_g = -BIG_NUMBER; FOR(i,(int)g.size()) { if(fabs(xi[i]) > max_g) { max_g = fabs(xi[i]); } }
		if(max_g > step_max) {
			max_g = step_max / max_g;
			FOR(i,(int)g.size()) { xi[i] *= max_g; }
		}



//		if(DoingSubset) { cout << "\nWorking with subset: "; FOR(i,(int)g.size()) { if(fabs(g[i]) > FLT_EPSILON) { cout << Model->m_vpAllOptPar[i]->Name() << " "; } } }
//		double temp_fp_s = fp; cout << "\n\tLinesearch --  fp: " << fp;

		last_improvement = fp;
		alpha = lnsrch(x,fp,g,xi,pold,&fret,Do_GS,Model); fp = fret;
		last_improvement = last_improvement - fp;

//		cout << " -> fp: " << fp << " imp(" << temp_fp_s - fp << ")";

		if(alpha < DBL_EPSILON || GradOK == false) { ResetHess = true; }
		if(its == NumberIter - 1) { break; }	// No point doing all this if about to step out.
		// Do the last five shuffle
		FOR(i,4) { Last5[i] = Last5[i+1]; } Last5[4] = fp;
		// Save the old gradients and step sizes
		FOR(i,n) { dg[i] = g[i]; oldxi[i] = xi[i]; }
		FOR(i,n)	{ xi[i] = *x[i] - pold[i]; }
		// Test for convergence
		for(test = 0.0, i = 0; i < n; i++)	{
			temp = fabs(xi[i]) / max(fabs(*x[i]),1);
			if(temp - test > DBL_EPSILON) 	{ test = temp; }
		}
		// Get the new derivatives
#if DEBUG_MULD_OPT > 1
		cout << "\n <<<<<<<<<<<<<<<<<<<<<<<<<<< UPDATE DERIVATIVES: exp: " << -fp << "; lnL: " << Model->lnL() << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
#endif
//		cout << "\nGetting derivatives: " << Model->lnL();
		g = Model->GetDerivatives(-fp,&GradOK);
//		cout << " ... Done derivatives: " << Model->lnL();
#if DEBUG_MULD_OPT > 1
		cout << "\n <<<<<<<<<<<<<<<<<<<<<<<<<<< DONE DERIVATIVES: exp: " << -fp << "; lnL: " << Model->lnL() << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
#endif
		// max_g = 0.0; FOR(i,(int)g.size()) { max_g = max(fabs(g[i]),max_g); }
		// Test for convergence on zero gradient
		test = 0.0; den = max(fret,1.0);
	    FOR(i,n)	{
	        temp = fabs(g[i]) * max(fabs(*x[i]), 1.0) / den;
	        if(temp > test)	{ test = temp; }
		}


		///////////////////////////// Checking standard exit conditions /////////////////////////////////
		// When the likelihood looks good but the gradients haven't converged try this
		// Break out when gradient passes appropriate test and when likelihood no longer decreasing]
		inc = fold - fp;
#if DEBUG_MULD_OPT > 0
		optout << "\nChecking return conditions; cur_lnL: " << fp << " cf. old: "<< fold <<"\n\tfabs(fp-fold)="<<fabs(fp-fold)<<" cf. " << ltol; if(fabs(fp - fold) < ltol) { optout << " *OKAY*"; }
		optout << "\n\tmax_g= " << test << " cf. gtol=" << gtol; if(test <= gtol) { optout << " *okay*"; }
		optout << "\n\tPredictedlnL: " << PredictlnL(Last5,5) << "; Last5: "; FOR(i,5) { optout << Last5[i] << " "; }
		if(inc < 0) { optout << "\nError: likelihood increasing in MulD_Optimise(...): fold: " << fold << ", fp: " << fp << "; diff= " << inc; exit(-1); }
#endif
		if(inc < ltol) {	// If there's no discernable increase in likelihood
			temp_lnL = -Model->lnL(true);
			if(fabs(temp_lnL - fp) > 1.0E-6) {	// If there's an error detected reset the likelihood. Hopefully this won't happen too often...
				 if(fabs(temp_lnL - fp) > 1.0E-1) { ResetHess = true; } // If it's a large error reset the hessian matrix
				 fp = temp_lnL;
				 LikTol = 0;
			} else {
				if(n==1) { LikTol = OptTol; } else { LikTol++; }
				if(fabs(inc) < FLT_EPSILON) { ResetHess = true; }
		} 	} else { LikTol = 0; }
		// Check return conditions
		if((test <= gtol && inc < ltol) || LikTol > OptTol)	{
			// Do a round of fast branch optimisation
			PredictedlnL = PredictlnL(Last5,5);
			// Check the predicted likelihood
			if((LikTol > OptTol) || fabs(PredictedlnL) < DBL_EPSILON) {
				if(TryReallyHard)	{
					FOR(i,20)	{
						fold = fp;
						fp = -Model->FastBranchOpt(-fp,FULL_LIK_ACC,NULL,10);			// Do a round of fast branch optimisation
#if DEBUG_MULD_OPT > 0
						optout << "\nHard["<<i<<"]: Branches-> " << fp;
#endif
					if(AllowOnlyParOpt || Model->ForceSeperateParOpt()) { fp = -DoOnlyParOpt(fp,gtol,ltol,x,Model,NI,OptTol); }	// Do a round of parameter optimisation
#if DEBUG_MULD_OPT > 0
						optout << " Pars-> " << fp;
#endif
						if(fabs(fp - fold) < FULL_LIK_ACC) { break; }
					}
					if(fabs(fold - fp) < ltol)	{
#if DEBUG_MULD_OPT > 0
						optout << "\n>>> Returning naturally. lnL= " << fp << " cf. old= " << fold << "; diff= " << fp - fold;
#endif
						// Final step is to check whether the
//		cout << "\n---Return B: fp: " << -fp << " cf. " << Model->lnL(true);
						FREEALL return -fp;
				}	} else {
#if DEBUG_MULD_OPT > 0
						optout << "\n>>> Returning naturally. lnL= " << fp << " cf. old= " << fold << "; diff= " << fp - fold;
#endif
						// Final step is to check whether the
#if DEBUG_MULD_OPT_SEP_FILES == 1 && DEBUG_MULD_OPT > 1
		cout << " " << -fp << "("<<its<<")";
#endif
//		cout << "\n---Return C: fp: " << -fp << " cf. " << Model->lnL(true);
						FREEALL return -fp;
				}
			}
			// Reset the hessian...
			ResetHess = true;
		}
		///////////////////////////// Checking lnL2Best exit conditition ///////////////////////////////
		// Return condition: use previous lnL to predict what the likelihood is going to be.
		if(lnL2Beat < fp && its > 10) {
			PredictedlnL = PredictlnL(Last5,5);
			if(fabs(PredictedlnL) > DBL_EPSILON && PredictedlnL - BIG_LNL_DIFF > lnL2Beat) {
#if DEBUG_MULD_OPT > 0
			optout << "\n>>> Returning due to lnL prediction";
#endif
			FREEALL;
#if DEBUG_MULD_OPT_SEP_FILES == 1 && DEBUG_MULD_OPT > 1
			cout << " " << -fp << "("<<its<<")";
#endif
			return -fp;
		}	}

//		cout << "\nThe last 5 likelihoods: "; FOR(i,5) { cout << "\t" << Last5[i]; }
		// Set warning if gradients and steps are not decreasing.
		// This is an assumption for using the hessian
		// Reset the hessian when this occurs
/*		flag = true;
		if(inc > 0.5) {
			FOR(i,n)	{
				if((fabs(dg[i] - g[i])  < 1.0E-3 && g[i] > 1.0E-3))	{
//					cout << "\nHess warning " << i<< " dg[i]= " << dg[i] << " cf. g[i]= " << g[i] << flush;
					FOR(j,n) { if(i!=j) { hessin[i][j] = hessin[j][i] = 0.0; } }
					flag = false;
		}	}	}
		if(flag == false) { HessWarning++; if(HessWarning >=5) { ResetHess = true; } } else { HessWarning = 0; }
*/
		// Reset the hessian if Expected convergence on condition that
		// i) haven't got it, or ii) expecting to exit next iteration
		// iii) The warning regarding the dot product has passed a threshold
		if(ResetHess || (its % 100 == 0 && its > 0))	{ // Otherwise reset the hessian
#if DEBUG_MULD_OPT > 0
			optout << "\n\tResetting hessian... lnL = " << Model->lnL();
#endif
			// Rescale the parameters
			Model->RedoScale();
			if(!Model->IsRMSDCalc())	{
#if DEBUG_MULD_OPT > 1
				cout << "\n <<<<<<<<<<<<<<<<<<<<<<<<<<< DOING RESCALE UPDATES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
#endif
				j = 1;
				if(TryReallyHard)	{ j = 20; }
				FOR(i,j)	{
					fold = fp;
//					cout << "\nDoing FBA: " << Model->lnL();
					fp = -Model->FastBranchOpt(-fp,FULL_LIK_ACC,NULL,10);			// Do a round of fast branch optimisation
//					cout << "\nDone FBA: " << Model->lnL();
#if DEBUG_MULD_OPT > 0
					optout << "\nHard["<<its<<"]: Branches-> " << fp;
#endif
			if(AllowOnlyParOpt || Model->ForceSeperateParOpt()) { fp = -DoOnlyParOpt(fp,gtol,ltol,x,Model,NI,OptTol); }	// Do a round of parameter optimisation
#if DEBUG_MULD_OPT > 0
					optout << " Pars-> " << fp;
#endif
					if(fabs(fp - fold) < FULL_LIK_ACC) { break; }
				}
#if DEBUG_MULD_OPT > 1
	cout << "\n <<<<<<<<<<<<<<<<<<<<<<<<<<< DONE RESCALE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
#endif
				max_g = 0.0; FOR(i,(int)g.size()) { max_g = max(fabs(g[i]),max_g); }
				g = Model->GetDerivatives(-fp,&GradOK);			// Get the new derivatives
			}
			// Reset the hessian
			FOR(i,n) {
				IMat(hessin,n);
				if(*x[i] > 1) { xi[i] = -g[i] * fabs(*x[i]); } else { xi[i] = -g[i]; }
				HessWarning = 0;
		}	} else {
			// Compute difference of gradients
			FOR(i,n) { dg[i] = g[i] - dg[i]; }
	        // And difference times current matrix
			FOR(i,n) {
				hdg[i] = 0.0;
				FOR(j,n) { hdg[i] += hessin[i][j] * dg[j]; }
			}
			// Calculate dot products for the denominators
			fac = fae = sumdg = sumxi = 0.0;
			FOR(i,n)	{
				fac += dg[i] * xi[i];
				fae += dg[i] * hdg[i];
				sumdg += SQR(dg[i]);
				sumxi += SQR(xi[i]);
			}
		    // Skip update if fac not sufficiently positive
		    if(fac * fac > FLT_EPSILON * sumdg * sumxi)	{
		        fac = 1.0 / fac;
		        fad = 1.0 / fae;
		        // The vector tht makes BFGS different from DFP
		        FOR(i,n)	{ dg[i] = fac * xi[i] - fad * hdg[i]; }
		        // The BFGS updating formula
				FOR(i,n)	{
					FOR(j,n)	{
			            hessin[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] +fae * dg[i] * dg[j];
			}	}	}
		    // Now calculate the next direction to go
			FOR(i,n)	{
			    xi[i] = 0.0;
				FOR(j,n)	{ xi[i] -= hessin[i][j] * g[j]; }
			}
			// Ensure its goes in the direction of the gradient
			//	This basically happens when the hessian is screwing up because the
			//	function isn't looking quadratic, defined here as being > 2.5
			FOR(i,n) {
				if(*x[i] > 0.5 && ((xi[i] < 0 && g[i] < 0) || (xi[i] > 0 && g[i] > 0)) ) {
					FOR(j,n) { hessin[i][j] = hessin[j][i] = 0.0; } hessin[i][i] = 1.0;	// Partially reset hessian
					if(*x[i] > 1 && fabs(g[i]) < 0.2) { xi[i] = -g[i] * fabs(*x[i]); } else { xi[i] = -g[i]; }
	}	}	}	}
    FREEALL
	WarningMulD = true;
#if DEBUG_MULD_OPT > 0
			optout << "\n>>> Returning due to excess of optimisation iterations";
#endif
//	cout << "\nDone " << its << " (all) iterations...";
#if DEBUG_MULD_OPT_SEP_FILES == 1 && DEBUG_MULD_OPT > 1
			cout << " " << -fp << "("<<its<<")";
#endif
	return -fp;
}

/////////////////////////////////// Likelihood predictor ///////////////////////////////////////
double PredictlnL(double *Last, int n)	{
	int i,count = 0;
	double OrilnL = Last[n-1], lnL = Last[n-1], dx1, dx2 = 0,temp = BIG_NUMBER;
	assert(n >3);
	if(!ALLOW_PREDICTLNL) { return 0.0; }
	cout << "\nDOUPEE!!!"; exit(-01);

	if(fabs(lnL) < DBL_EPSILON || Last[0] + DBL_EPSILON > BIG_NUMBER) { return 0.0; }
//	cout << "\n\t\tlast 5: "; FOR(i,n) { cout << Last[i] << " "; }
	// Get 1st derivative
	dx1 = Last[1] / Last[0];
	// Get 2nd derivative (average of n-2 estimates)
	FOR(i,n-2)	{
		if(fabs(Last[i+1]) < DBL_EPSILON || Last[i+1] > Last[i]){ return 0; }
		dx2 += (Last[i+2] * Last[i]) / (Last[i+1] * Last[i+1]);
	}
	dx2 /= (double) n-2;	// Make average
//	cout << " dx1= " << dx1 << "; dx2= " << dx2;
//	cout << "\nOri_lnL: " << lnL;
	// Do the prediction
	if(dx1 >= 1.0) { return 0; } // Return a big likelihood
	while(fabs(lnL - temp) > FULL_LIK_ACC && dx1 < 1.0)	{
		temp = lnL; dx1 *= dx2;
		if(dx1 > 1.0) { break; }
		lnL *= dx1;
//		cout << "\n\tstep " << count << " dx1: " << dx1 << "; dx2: "<< dx2 << " lnL: " << lnL;
		if(count ++ > 20) { if(OrilnL - lnL > FULL_LIK_ACC) { return lnL; } else { cout << "%"; return 0; } } // Ensure the loop doesn't get stuck
	}
//	cout << " ... return: " << lnL;
	return lnL;
}

/////////////////////////////////////////////////////////////////////////
// Optimise only parameters -- not branches
//	equivalent of Model->DoFastBraOpt(..)
double DoOnlyParOpt(double OrilnL,double gtol ,double ltol,vector <double *> x,CBaseModel *Model,int NI,int OptTol,int MaxIter)	{
	int i,Iter =0;
	double Ret = -fabs(OrilnL), Old_lnL;
	vector <int> Order;
	if(ALLOW_ONLYPAROPT == 0 && !Model->ForceSeperateParOpt()) { return -fabs(OrilnL); }
	// Loop through the optimised values and optimise them one at a time
	FOR(i,(int)Model->m_vpAllOptPar.size()) {
//		cout << "\nParameter["<<i<<"] " << Model->m_vpAllOptPar[i]->Name() << " = " <<  Model->m_vpAllOptPar[i]->Val() << " == " << *x[i];

		if(Model->IsLikelihoodCalc() && Model->m_vpAllOptPar[i]->Name().find("Branch") != string::npos)	{ continue; }

		if(Model->m_vpAllOptPar[i]->Name().find("Inv") == string::npos) { continue; }
		Order.push_back(i);
	}
	if(Order.empty()) { return -fabs(Ret); }
//	cout << "\n================================== Doing Only Par =======================";
	FOR(Iter,MaxIter)	{
		Old_lnL = Ret;
		random_shuffle(Order.begin(),Order.end());
//		cout << "\n>>>>>>>>>>>>>>>>>>>>>>> Start " << Iter << ": " << Ret << flush;
		FOR(i,(int)Order.size())	{
			// Skip branches
			Ret = GoldenSection(fabs(Ret),x[Order[i]],Model->m_vpAllOptPar[Order[i]],Model);
//			cout << "\nOrderDone: " << Ret << " cf. " << Model->lnL();
		}
//		cout << "\n\tDone Par" << flush;
		if(Iter % 2 == 0) { Ret = Model->FastBranchOpt(Ret); }
//		cout << "\n\tIter " << Iter << ": " << Old_lnL << " -> " << Ret << " cf. " << Model->lnL() << flush;
		if(fabs(Old_lnL - Ret) < FULL_LIK_ACC) { break; }
	}
//	cout << "\n\tDone -> " << Ret << " cf. " << Model->lnL() << flush;
//	exit(-1);
	return -fabs(Ret);
}

/////////////////////////////////////////////////////////////////////////
// Check that a parameter is at a true optima

bool CheckAllPar(CBaseModel *M, double lnL, vector <double *> x, double Tol, ostream &os, bool ForceShow)	{
	int i;
	bool RetVal = true;
	if(os != cout) { os << "\nHard checking parameter estimates:"; }
	FOR(i,(int)x.size()) { if(!HardCheckOpt(M,lnL,x[i],Tol,i,os,ForceShow)) { RetVal = false; }	}
	return RetVal;
}

bool HardCheckOpt(CBaseModel *M, double lnL, double *x, double Tol, int ParNum,ostream &os, bool ForceShow)	{
	double l_lnL,r_lnL, x_ori = *x;
	double i;
	// Get left likelihood
	*x = x_ori - Tol; l_lnL = M->lnL();
	// Get right likelihood
	*x = x_ori + Tol; r_lnL = M->lnL();
	*x = x_ori;	// Reset value
	// Do some debug output if required
	if(os!= cout || ForceShow) {
		int prec = os.precision(); os.precision(8);
		os << "\n\tPar ["<<ParNum <<"] " << M->m_vpAllOptPar[ParNum]->Name() << ": " << *x << " == " << lnL;
		for(i=1;i<1000;i*=10) {
			*x = x_ori - (Tol* i); l_lnL = M->lnL();	// Get left likelihood
			*x = x_ori + (Tol* i); r_lnL = M->lnL();	// Get right likelihood
			*x = x_ori;									// Reset value
			os << "\n\t\t["<<i<<"] left_lnL(par="<< x_ori - (Tol* i) <<"): "<< l_lnL << "(" << lnL - l_lnL << "); ";
			os << "right_lnL (par="<<x_ori+ (Tol* i) << "): " << r_lnL << "(" << lnL - r_lnL << ")";
		}
		os.precision(prec);
		if(lnL + FULL_LIK_ACC >= l_lnL && lnL + FULL_LIK_ACC >= r_lnL) { os << "\tOkay"; } else {
			os << "\tFailed";
//			Parameter1DGraph(M,M->m_vpAllOptPar[ParNum],max(0,x_ori - (0.01*25)),x_ori + (0.01*25),0.01,false,os);
		}
	}
	// Return appropriate value
	if(lnL + FULL_LIK_ACC >= l_lnL && lnL + FULL_LIK_ACC >= r_lnL) { return true; }
 	return false;
}

///////////////////////////////////////////////////////////
// Function to get pairwise distances from data
// These pairwise distances are used in the sequence plucking algorithms
// Parameters are *NOT* optimised
// ---
// Args:
//	Model ==		Current model
//  CurrentPW ==	Current pairwise distance matrix
//  GetJCDist ==	Whether to only get JC distances (and StdErrs)
//  GetStdErr ==	Whether to get the standard errors for the distances
vector <double> GetPW(CBaseModel *Model, vector <double> *CurrentPW, bool GetJCDist, bool GetStdErr)	{
	int i,j;
	vector <double> Dists(Model->m_pData->m_iNoSeq*Model->m_pData->m_iNoSeq, -1.0);
	if(GetJCDist == false)	{	/////////// Get full likelihood distances //////////
		bool DoingRMSD = Model->IsRMSDCalc();
		Model->DoLikelihoodCalc();
		CData *OldDat = Model->m_pData;
		CTree *OldTree = Model->m_pTree;
		if(CurrentPW!=NULL) { if(CurrentPW->empty()) { CurrentPW = NULL; } }
		if(CurrentPW == NULL)	{
			FOR(i,OldDat->m_iNoSeq)	{
				Dists[(i*OldDat->m_iNoSeq)+i] = 0.0;
				FOR(j,i)	{
					Dists[(i*OldDat->m_iNoSeq)+j] = Dists[(j*OldDat->m_iNoSeq)+i] = GetDist(Model,OldDat,i,j); }
		}	} else {
			assert(CurrentPW->size() == OldDat->m_iNoSeq * OldDat->m_iNoSeq);
			FOR(i,OldDat->m_iNoSeq)	{
				Dists[(i*OldDat->m_iNoSeq)+i] = 0.0;
				FOR(j,i)	{ Dists[(i*OldDat->m_iNoSeq)+j] = Dists[(j*OldDat->m_iNoSeq)+i] = GetDist(Model,OldDat,i,j,CurrentPW->at((i*OldDat->m_iNoSeq)+j)); }
		}	}
		// The way this is incorporated is not the most efficient as several things are set up multiple times
		if(GetStdErr)	{
			FOR(i,OldDat->m_iNoSeq)	{
				FOR(j,i)	{ Dists[(i*OldDat->m_iNoSeq)+j] = sqrt(GetDistVar(Model,OldDat,i,j,Dists[(j*OldDat->m_iNoSeq)+i])); }
		}	}
		// Restore models original values are return
		Model->m_pTree = OldTree;
		Model->m_pData = OldDat;
		if(DoingRMSD) { Model->DoRMSDCalc(); }
	} else {					//////////// Get only JC distances (using formula) ///////
		FOR(i,Model->m_pData->m_iNoSeq)	{
			Dists[(i*Model->m_pData->m_iNoSeq)+i] = 0.0;
			FOR(j,i)	{ Dists[(i*Model->m_pData->m_iNoSeq)+j] = Dists[(j*Model->m_pData->m_iNoSeq)+i] = Model->m_pData->PoissonDist(i,j); }
		}
		if(GetStdErr)	{
			FOR(i,Model->m_pData->m_iNoSeq)	{
				FOR(j,i)	{ Dists[(i*Model->m_pData->m_iNoSeq)+j] = sqrt(Model->m_pData->PoissonVar(i,j)); }
	}	}	}
	FOR(i,(int)Dists.size()) {
		if(Dists[i] < 0) { Dists[i] = 0.0; } if(isnan(Dists[i])) { Dists[i] = 20; }
	}
	return Dists;
}

void SortSTreeScore(vector <STreeScore> *x)	{
	int i,j;
	STreeScore Temp;
	FOR(i,(int)x->size())	{
		for(j=i;j<(int)x->size();j++)	{
			if(x->at(j).Lik > x->at(i).Lik) { Temp = x->at(i); x->at(i) = x->at(j); x->at(j) = Temp; }
	}	}


}

double GetDist(CBaseModel *M, CData *D,int Seq1, int Seq2, double CurrentDist)	{
	double dist;
	if(Seq1==Seq2) { return 0.0; }
	static CTree T("(1,2);",2);
	M->PreparePairwiseCalc(Seq1,Seq2,&T);
	if(CurrentDist >= 0) { T.SetB(0,CurrentDist); }
	GoldenSection(M->lnL(),T.pBra(0)->OptimiserValue(),T.pBra(0),M);
	dist = T.B(0);
	M->CleanCPMapping();
	return dist;
}
double GetDistVar(CBaseModel *M, CData *D,int Seq1, int Seq2, double CurrentDist)	{
	double dist;
	if(Seq1==Seq2) { return 0.0; }
	static CTree T("(1,2);",2);
	M->PreparePairwiseCalc(Seq1,Seq2,&T);
	if(CurrentDist >= 0) { T.SetB(0,CurrentDist); } else { Error("\nNeed to provide valid distance to GetDistVar...\n\n"); }
	// Get the stderr
	dist = Get2ndDer(T.pBra(0),M,M->lnL());
	M->CleanCPMapping();
	return dist;
}


//////////////////////////////////////////////
// Derivative functions
double Get2ndDer(CPar *Par, CBaseModel *M, double lnL)	{
	assert(Par != NULL && M != NULL);
	double x1,x2,x3,dx, p_old = Par->Val();
	if(lnL - FLT_EPSILON < -BIG_NUMBER) { x2 = M->lnL(); }
	else								{ x2 = lnL; }
	dx = max(1.0E-4,Par->Val() * 1.0E-4);
	if(!Par->CheckBound()) { Error("\nTrying to get derivative on boundary..."); }
	// Do left
	Par->SetVal(p_old - dx,true,true);
	x1 = M->lnL();
	// Do right
	Par->SetVal(p_old + dx,true,true);
	x3 = M->lnL();
	// Reset parameter
	Par->SetVal(p_old); if(Par->Val() != p_old) { Error("\nParameter doesn't behave normally in Get2ndDer()...\n\n"); }
	// Get the derivative
	cout.flags(ios::fixed); cout .precision(8);
//	cout << "\nx1: " << p_old - dx <<" == " << x1 << "\nx2: " << p_old <<" == " << x2 << "\nx3: " << p_old + dx <<" == " << x3;
//	cout << "\nGrad: " << (x3 + x1 - (x2*2)) / (dx * dx) << " cf. " << (dx * dx) / ((x2*2) - x1 - x3);
	if(x1 > x2) { Error("\nIn Get2ndDer() have x1 > x2...\n"); }
	if(x3 > x2) { Error("\nIn Get2ndDer() have x3 > x2...\n"); }
	return (dx * dx) / ((x2*2) - x1 - x3);
}

//////////////////////////////////////////////
// Graphing functions

void Parameter1DGraph(CBaseModel *M, CPar *Par,double low,double high, double step, bool OptAll,ostream &os)	{
	double i;
	Par->SetOptimise(false);
	os << "\nParameter 1D graph for " << Par->Name() << ": ";
	for(i=low;i<high;i+=step)	{
		os << "\n" << i << flush;
		Par->SetVal(i,OptAll,OptAll);
		if(OptAll == true) { os << "\t" << FullOpt(M,true,true); } else { os << "\t" << M->lnL(); }
}	}
void Parameter2DGraph(CBaseModel *M, CPar *Par1,CPar *Par2, double l1,double h1,double s1,double l2,double h2,double s2, bool OptAll,ostream &os)	{
	double i,j;
	Par1->SetOptimise(false); Par2->SetOptimise(false);
	os << "\nParameter 2D graph for " << Par1->Name() << " and " << Par2->Name() << ": ";
	os << "\n\nX";
	for(j=l2;j<h2;j+=s2)	{ os << "\t" << j; }
	for(i=l1;i<h1;i+=s1)	{
		os << "\n" << i;
		Par1->SetVal(i,true,true);
		for(j=l2;j<h2;j+=s2)	{
			Par2->SetVal(j,true,true);
			if(OptAll == true) { os << "\t" << FullOpt(M,true,true); } else { os << "\t" << M->lnL(); }
}	}	}

