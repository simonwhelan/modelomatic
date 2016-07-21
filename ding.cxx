/*	///////////////////////////////////////////////////////////////////////
		ModelOMatic -- Program for assessing the fit of lots of models
			Simon Whelan, Uppsala University

		ModelOMatic <data_file> <tree_file> <output_file> <Genetic_code; default=universal> <fast>

	Option definitions
	<data_file>		Any sequence data file in phylip or sequential format
	<tree_file>		The tree
	<output_file>	Output file
	<Genetic_code> 	a number between 0 and something
	<fast>			if fast is present it will only optimise the branches and alpha parameter for the first model of a data type and then fix for all subsequent
	/////////////////////////////////////////////////////////////////////// */

//#include "./modelomatic.h"
//#include <time.h>
//#include <set>


#if DO_MEMORY_CHECK == 1
extern CMemChecker memory_check;
#endif

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter, MakeQ_Log_Counter, MakePT_Log_Counter, LFunc_Log_Counter, SubLFunc_Log_Counter, SPR_Log_Counter;
#endif

#include "optimise.h"
#include "TreeList.h"
#include "model.h"
#include "interface.h"
#include "Leaphy.h"
#include "ding.h"
#include "PfamModel.h"

////////////////////////////////////////////////////////////////
// External variables
extern vector <STabuTree> TabuTrees;
extern vector <double> PWDists;
extern int TABU_RADIUS;

#define DO_PFAM 1		// Whether to do the Pfam analysis in PfamModel.cxx

int ding()	{
		int Ret = 0;
		int GeneticCode = 0;
		int count = 0;
		// Stuff from Leaphy
		int i,j,k,l,NumModelReruns = 1;
		long RandomSeed = 0;
		string temp_string, Name, outfilestring;
		vector <string> Toks;

#if DO_PFAM == 1
		PfamModelAnalysis();
		exit(-1);
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Stuff for Ding's work on tree support
	// -----------------
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Files expected
	// Tree files: mito55_Simon_H1.tre, mito55_Simon_H2.tre, mito55_Simon_H3.tre, mito55_Simon_H0.tre
	// Data file: mito55_exRog.phy
	// Eukaryote name file: mito55_Euk.names
	// Rickettsia name file: mito55_Rik.names
	// Others name file: mito55_Oth.names
	// Proteobacteria name file: mito55_Pro.names

	// Work in progress files expected
	// Tree files: mito55_exRog_H1.tre mito55_exRog_H2.tre mito55_exRog_H3.tre

	cout << "\nInitialising data and hypotheses" << flush;

	// Input
	//   Data
	CData DingDat("MitoData.phy",AA);
	cout << "\nData read successfully ("<<DingDat.m_iNoSeq << " x " << DingDat.m_iSize << ")";
	// Name files
	vector <int> EukNames = ReadNameFile("mito55_Euk.names",&DingDat);
	vector <int> RikNames = ReadNameFile("mito55_Rik.names",&DingDat);
	vector <int> OthNames = ReadNameFile("mito55_Oth.names",&DingDat);
	vector <int> ProNames = ReadNameFile("mito55_Pro.names",&DingDat);
	vector <int> RootNames = ReadNameFile("mito55_root.names",&DingDat);
	vector <vector <int> > SubNames;
	assert(EukNames.size() + RikNames.size() + OthNames.size() + ProNames.size() == DingDat.m_vsName.size());
	SubNames.push_back(EukNames); SubNames.push_back(RikNames); SubNames.push_back(OthNames); SubNames.push_back(ProNames);
	vector <int> Temp;
	// Build pairwise sets of these four => 4C2 = 6
	Temp = VecCon(EukNames,RikNames); sort(Temp.begin(),Temp.end());
	vector <int> EukRik = Temp;
	Temp = VecCon(EukNames,OthNames); sort(Temp.begin(),Temp.end());
	vector <int> EukOth = Temp;
	Temp = VecCon(EukNames,ProNames); sort(Temp.begin(),Temp.end());
	vector <int> EukPro = Temp;
	Temp = VecCon(RikNames,OthNames); sort(Temp.begin(),Temp.end());
	vector <int> RikOth = Temp;
	Temp = VecCon(RikNames,ProNames); sort(Temp.begin(),Temp.end());
	vector <int> RikPro = Temp;
	Temp = VecCon(OthNames,ProNames); sort(Temp.begin(),Temp.end());
	vector <int> OthPro = Temp;



	cout << "\nNames read successfully" << flush;
	bool DoGamma = true;
	// 	 Tree H1
	cout << "\nReading trees: H1" << flush;
	CTree OriT_H1("mito55_exRog_H1.tre", true, &DingDat);
	OriT_H1.OutBra(); OriT_H1.OutName();
	CTree T_H1 = OriT_H1;
	int H1_Split = FindSpecialSplit(&T_H1,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
//	PaintTree(&T_H1, H1_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H1(&DingDat,&T_H1,"LG",true,(double*)dLGVal,(double*)dLGFreq);
//	if(DoGamma) { LG_H1.MakeGammaModel(0,4,0.63); }
	cout.precision(10);
	CHeteroEMP LGHet_H1(&DingDat,&T_H1,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);


//	CWAG LG_H1(&DingDat,&T_H1);
	//   Tree H2
	cout << " ... H2" << flush;
	CTree OriT_H2("mito55_exRog_H2.tre", true, &DingDat);
	OriT_H2.OutBra(); OriT_H2.OutName();
	CTree T_H2 = OriT_H2;
	int H2_Split = FindSpecialSplit(&T_H2,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
//	PaintTree(&T_H2, H2_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H2(&DingDat,&T_H2,"LG",true,(double*)dLGVal,(double*)dLGFreq);
//	if(DoGamma) { LG_H2.MakeGammaModel(0,4,0.63); }
//	CWAG LG_H2(&DingDat,&T_H2);
	CHeteroEMP LGHet_H2(&DingDat,&T_H2,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);

	//   Tree H3
	cout << " ... H3" << flush;
	CTree OriT_H3("mito55_exRog_H3.tre", true, &DingDat);
	OriT_H3.OutBra(); OriT_H3.OutName();
	CTree T_H3 = OriT_H3;
	int H3_Split = FindSpecialSplit(&T_H3,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
//	PaintTree(&T_H3, H3_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H3(&DingDat,&T_H3,"LG",true,(double*)dLGVal,(double*)dLGFreq);
//	if(DoGamma) { LG_H3.MakeGammaModel(0,4,0.63); }
//	CWAG LG_H3(&DingDat,&T_H3);
	CHeteroEMP LGHet_H3(&DingDat,&T_H3,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);

	cout << " ... done" << flush;
	cout << "\nTrees and data input" << flush;

#define EXAMINE_AA_CONTENT 1
#if EXAMINE_AA_CONTENT == 1
	cout << "\nExamining AA content across the tree";
	// Branches and splits
	vector <SSplit> CompSet1;
	CTree CompTree = OriT_H1;
	CompSet1 = CompTree.BuildSplits();
	vector <double> X_BranchComp(CompSet1.size(),0), LeftComp,RightComp, LeftFreq, RightPred;
	vector <int> LeftSet,RightSet;
	double Count, CompCount = -1;
	double X;
	// Loop through the branches. The Larger set is always assigned to left and used to predict the composition of the smaller set
	FOR(i,(int)CompSet1.size()) {
		// Set up
		LeftComp.clear(); LeftComp.assign(20,0);
		RightComp.clear(); RightComp.assign(20,0);
		if(CompSet1[i].Left.size() > CompSet1[i].Right.size()) { LeftSet = CompSet1[i].Left; RightSet = CompSet1[i].Right; }
		else { RightSet = CompSet1[i].Left; LeftSet = CompSet1[i].Right; }
		// Get the count vectors
		// Left (biggest)
		FOR(j,LeftSet.size()) {
			FOR(k,DingDat.m_iSize) {
				if(DingDat.m_ariSeq[LeftSet[j]][k] != 20) { LeftComp[DingDat.m_ariSeq[LeftSet[j]][k]] += DingDat.m_ariPatOcc[k]; }
		}	}
		// Right
		FOR(j,RightSet.size()) {
			FOR(k,DingDat.m_iSize) {
				if(DingDat.m_ariSeq[RightSet[j]][k] != 20) { RightComp[DingDat.m_ariSeq[RightSet[j]][k]] += DingDat.m_ariPatOcc[k]; }
		}	}
		// Error checks
		Count = 0; FOR(j,20) { Count += LeftComp[j] + RightComp[j]; }
		if(CompCount < 0) { CompCount = Count; }
		else { assert(fabs(CompCount-Count) < 1.0E-3); }
		// Do the statistics test (Obs = RightComp, Exp = RightPred) => X = (Obs - Exp)^2 / Exp
		RightPred = NormaliseVector(LeftComp,Sum(&RightComp));
		Count = 0; FOR(j,20) { Count += pow(RightComp[j] - RightPred[j],2) / RightPred[j]; }
		X_BranchComp[i] = Count;
//		cout << "\nPartition["<<i<<"]: " << LeftSet << " | " << RightSet;
//		cout << "\nLeft      " << LeftComp;
//		cout << "\nRight     " << RightComp << " == " << Sum(&RightComp);
//		cout << "\nRightPred " << RightPred << " == " << Sum(&RightPred);
//		cout << "\nX^2       " << Count;
	}
	// Put the values on a tree and print the tree out
	CompTree.m_vdBootStrap = X_BranchComp;

	CompTree.OutBra(); CompTree.OutName();
	cout << "\n\nThe tree: " << endl << CompTree;
	exit(-1);
#endif


#define BUILD_NULL_TREES 0
#if BUILD_NULL_TREES == 1
	// Get the bifurcating tree
	CTree OriT_H0("mito55_exRog_H0.tre", true, &DingDat);
	// Find the 'special branch'
	int SpecB = FindSpecialSplit(&OriT_H0,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch


#endif


/*
	// Match branches across the hypotheses
	vector <SSplit> BSet1, BSet2, BSet3;	// Storage for the splits
	vector <vector <int> > BMap; vector <int> NewBMap;		// BMap is a vector [#Branches][2] Where the first index is the branch in H1 and the [2] is the corresponding branch in other trees; NewBMap[3] is the branch that varies across all trees (the one that determines the tree
	// Build sets
	BSet1 = OriT_H1.BuildSplits();
	BSet2 = OriT_H2.BuildSplits();
	BSet3 = OriT_H3.BuildSplits();
	// Now build the mapping
	vector <int> AddRow(2,-1);
	vector <bool> Checker;
	NewBMap.assign(2,-1);
	// Set 2
	cout << "\nNames";
	FOR(i,DingDat.m_vsName.size()) { cout << "\n["<<i<<"]:\t" << DingDat.m_vsName[i]; }
	Checker.clear(); Checker.assign((int)BSet2.size(),true);
	FOR(i,(int)BSet1.size())	{
		AddRow[0] = AddRow[1] = -1;
		cout << "\nBranch["<<i<<"]:  " << BSet1[i].Left << " | " << BSet1[i].Right;
		// Find Split i in BSet2
		FOR(j,(int)BSet2.size()) {
			if(BSet1[i].Left == BSet2[j].Left && BSet1[i].Right == BSet2[j].Right) { cout << " --> " << j; assert(Checker[j]); AddRow[0] = j; Checker[j] = false; break; }
		}
		if(j == BSet2.size()) {
			cout << "\n<><><><><><> CHECK <><><><><><><><><><";
			FOR(j,(int)BSet2.size()) {
					cout << "\nBranch2["<<j<<"]: " << BSet2[j].Left << " | " << BSet2[j].Right;
			}
			cout << "\n<><><><><><> CHECK <><><><><><><><><><";
		}
	}
	cout << "\nChecker: " << Checker;
	FOR(j,(int)BSet2.size()) {
		cout << "\nChecker["<<j<<"]: " << Checker[j] << flush;
		if(Checker[j] == true) { assert(NewBMap[0] == -1); NewBMap[0] = j; }
	}
	assert(NewBMap[0] != -1);
	// Set 3
	Checker.clear(); Checker.assign((int)BSet2.size(),true);
	FOR(i,(int)BSet1.size())	{
		AddRow[0] = AddRow[1] = -1;
		// Find Split i in BSet3
		FOR(j,(int)BSet3.size()) {
			if(BSet1[i].Left == BSet3[j].Left && BSet1[i].Right == BSet3[j].Right) { assert(Checker[j]); AddRow[0] = j; Checker[j] = false; break; }
		}
	}
	FOR(j,(int)BSet3.size()) {
		if(Checker[j] == true) { assert(NewBMap[1] == -1); NewBMap[1] = j; }
	}
	assert(NewBMap[1] != -1);
*/

	cout << "\nEuk: " << EukNames;
	cout << "\nPro: " << ProNames;
	cout << "\nRik: " << RikNames;
	cout << "\nOth: " << OthNames;


	//cout << "\n"

	bool DoNormal = true;
	bool DoNormalGamma = true;
	bool DoHet = true;
	bool DoHetGamma = true;
	bool DoHetFreq = false;

	cout << "\nInitial Likelihoods"; cout.precision(10);

	CBaseModel *M;
	CBaseModel *MH;
	cout << "\n------------------------- H1------------------------------";
	M = &LG_H1;
	MH = &LGHet_H1;
	if(DoNormal)	{
		cout << "\nStart " << M->Name() << " " << M->lnL() << flush;
		cout << "\nOptimising...";
//		FullOpt(M);
		cout << "\nLikelihood: "<< M->lnL();
//		FOR(i,5) { cout << "\nSite["<<i<<"]: " << M->m_arL[i].LogP(); FOR(j,M->m_vpProc.size()) { cout << "  " << M->m_vpProc[j]->m_ardL[i]; } }
		cout << "\nTree\n" << *M->Tree();

	}
	if(DoNormalGamma)	{
		M->MakeGammaModel(0,4,0.63);
		cout << "\n\nAdded gamma";
		cout << "\nStart " << M->Name() << ": " << M->lnL() << flush;
		cout << "\nOptimising...";
//		FullOpt(M);
		cout << "\nLikelihood+dG: "<< M->lnL();
//		FOR(i,5) { cout << "\nSite["<<i<<"]: " << M->m_arL[i].LogP(); FOR(j,M->m_vpProc.size()) { cout << "  " << M->m_vpProc[j]->m_ardL[i]; } }
		cout << "\nTree\n" << *M->Tree();
//		cout << "\nModel+dG\n" << *M;
	}
	if(DoHet)	{
		cout << "\nStart " << MH->Name() << ": " << MH->lnL() << flush;
		cout << "\nOptimising...";
//		FullOpt(MH);
		cout << "\nLikelihood: "<< MH->lnL();
//		FOR(i,5) { cout << "\nSite["<<i<<"]: " << M->m_arL[i].LogP(); FOR(j,MH->m_vpProc.size()) { cout << "  " << MH->m_vpProc[j]->m_ardL[i]; } }
		cout << "\nTree\n" << *MH->Tree();
	}
	if(DoHetGamma)	{
		MH->MakeGammaModel(0,4,0.63);
		cout << "\n\nAdded gamma" << flush;
		cout << "\nStart " << MH->Name() << ": " << MH->lnL() << flush;
		cout << "\nOptimising...";
//		FullOpt(MH,true,true,DoHetFreq);
		cout << "\nLikelihood+dG: "<< MH->lnL() << flush;
//		FOR(i,5) { cout << "\nSite["<<i<<"]: " << M->m_arL[i].LogP(); FOR(j,MH->m_vpProc.size()) { cout << "  " << MH->m_vpProc[j]->m_ardL[i]; } }
		cout << "\nTree\n" << *M->Tree();
//		cout << "\nModel+dG\n" << *MH << flush;
	}
	M =NULL;
	MH = NULL;

	exit(-1);

	cout << "\n------------------------- H2------------------------------";
	M = &LG_H2;
	MH = &LGHet_H2;
	if(DoNormal)	{
		cout << "\nStart: " << M->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(M);
		cout << "\nLikelihood: "<< M->lnL();
		cout << "\nTree\n" << *M->Tree();

	}
	if(DoNormalGamma)	{
		M->MakeGammaModel(0,4,0.63390);
		cout << "\n\nAdded gamma";
		cout << "\nStart+dG: " << M->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(M);
		cout << "\nLikelihood+dG: "<< M->lnL();
		cout << "\nModel+dG\n" << *M;
	}
	if(DoHet)	{
		cout << "\nStart: " << MH->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(MH);
		cout << "\nLikelihood: "<< MH->lnL();
		cout << "\nTree\n" << *MH->Tree();

	}
	if(DoHetGamma)	{
		MH->MakeGammaModel(0,4,0.63390);
		cout << "\n\nAdded gamma";
		cout << "\nStart+dG: " << MH->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(MH);
		cout << "\nLikelihood+dG: "<< MH->lnL();
		cout << "\nModel+dG\n" << *MH;
	}
	M =NULL;
	MH = NULL;

	cout << "\n------------------------- H3------------------------------";
	M = &LG_H3;
	MH = &LGHet_H3;
	if(DoNormal)	{
		cout << "\nStart: " << M->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(M);
		cout << "\nLikelihood: "<< M->lnL();
		cout << "\nTree\n" << *M->Tree();

	}
	if(DoNormalGamma)	{
		M->MakeGammaModel(0,4,0.63390);
		cout << "\n\nAdded gamma";
		cout << "\nStart+dG: " << M->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(M);
		cout << "\nLikelihood+dG: "<< M->lnL();
		cout << "\nModel+dG\n" << *M;
	}
	if(DoHet)	{
		cout << "\nStart: " << MH->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(MH);
		cout << "\nLikelihood: "<< MH->lnL();
		cout << "\nTree\n" << *MH->Tree();

	}
	if(DoHetGamma)	{
		MH->MakeGammaModel(0,4,0.63390);
		cout << "\n\nAdded gamma";
		cout << "\nStart+dG: " << MH->lnL() << flush;
		cout << "\nOptimising...";
		FullOpt(MH);
		cout << "\nLikelihood+dG: "<< MH->lnL();
		cout << "\nModel+dG\n" << *MH;
	}
	M =NULL;
	MH = NULL;


	exit(-1);
	////////////////// END DING WORK //////////////////

	return Ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DING SPECIFIC FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read a set of names and return a vector of indexes specifying those names
vector <int> ReadNameFile(string NameFile, CData *D)	{
	int i;
	vector <int> RetSet;
	vector <bool> Checkers((int)D->m_vsName.size(),true);
	string store;
	vector <string> Toks;
//	cout << "\nReading " << NameFile;
	//Open file
	FINOPEN(FileIn,NameFile.c_str());
	getline(FileIn,store);
	while(!FileIn.eof())	{
		Toks = Tokenise(store);
		assert((int)Toks.size() == 1);
//		cout << "\n\t" << Toks[0];
		FOR(i, (int)D->m_vsName.size()) {
			if(Toks[0].find(D->m_vsName[i]) != string::npos) {
//				cout << " found ["<<i<<"] == " << D->m_vsName[i];
				assert(Checkers[i]); Checkers[i] = false; RetSet.push_back(i); break; }
		}
		if(i == (int)D->m_vsName.size()) { cout << "\nCannot find: " << Toks[0] << " in name list"; exit(-1); }
		// Move on
		getline(FileIn,store);
		if(FileIn.eof()) { break; }
	}

	FileIn.close();
	sort(RetSet.begin(),RetSet.end());
	return RetSet;
}

/////////////////////// Find the branch that's unique to a hypothesis //////////////////
// NB: This is poorly written code and not generalised. Will need rewriting for more general questions
int FindSpecialSplit(CTree *T,vector <int> S1,vector <int> S2,vector <int> S3,vector <int> S4,vector <int> S5,vector <int> S6)	{
	int i,j,RetSplit = -1;
	vector <vector <int> > SplitList;
	SplitList.push_back(S1); SplitList.push_back(S2); SplitList.push_back(S3); SplitList.push_back(S4); SplitList.push_back(S5); SplitList.push_back(S6);
	vector <SSplit> Splits = T->BuildSplits();
	FOR(i,Splits.size())	{
		FOR(j,6) {
			if(Splits[i].Left == SplitList[j]) { break; }
		}
		if(j == 6) { continue; } 	// Jump if not matches
		FOR(j,6) {
			if(Splits[i].Right == SplitList[j]) { break; }
		}
		if(j == 6) { Error("\nSplit list only matches one side...\n"); } 	// Jump if not matches
//		cout << "\nHave match on branch[" << i << "]: " << T->B(i);
		RetSplit = i; break;
	}
	assert(RetSplit != -1);
	return RetSplit;
}

//////////////////////// Label branches for heterogeneous model ///////////////////
void PaintTree(CTree *T, int KeyBra, vector <int> S1, vector <int> S2, vector <int> S3, vector <int> S4)	{
	int i,j;
	// Set initial labels to 0 (will be only KeyBra at end)
	T->SetLabels(0);

	T->PaintSubtree(S1,1);
	T->PaintSubtree(S2,2);
	T->PaintSubtree(S3,3);
	T->PaintSubtree(S4,4);

	// Error checking outputting
//	cout << "\nInto PaintTree"; T->OutLabel(); T->OutBra(); T->OutName(); cout << "Tree\n" << *T; cout << "\nBranchLabels: " << T->BranchLabels();

}

void PaintTree(CTree *T, vector <vector <int> > Subtrees)	{
	int i;
	// Check the Subtree labels are unique, so only one species occurs in each group
	vector <int> vTemp; FOR(i,(int)Subtrees.size()) { vTemp.insert(vTemp.end(),Subtrees[i].begin(),Subtrees[i].end()); }
	sort(vTemp.begin(),vTemp.end()); vector <int>::iterator last = unique(vTemp.begin(),vTemp.end()); vTemp.erase(last,vTemp.end());
	// Paint the tree with the last label first
	T->SetLabels((int)Subtrees.size());
	// Now paint the rest of the tree
	FOR(i,(int)Subtrees.size()) { T->PaintSubtree(Subtrees[i],i); }
}

 // Returns the root branch from a set of Root sequences that specify the branch
int GetRoot(CTree * T, vector <int> RootSeqs)	{
	int i;
	vector <SSplit> Splits = T->BuildSplits();
//	cout << "\nRoot" << RootSeqs;
	FOR(i,Splits.size()) {
		if(Splits[i].Left == RootSeqs || Splits[i].Right == RootSeqs) {
//			cout << "\nFound match: ";
//			cout << "\nLeft:  " << Splits[i].Left;
//			cout << "\nRight: " << Splits[i].Right;
			return i;
		}
	}
	cout << "\nError: Root split failure..."; exit(-1);
	return -1;
}

//////////////////////////////////////////////////////////////////////////////////
// Heterogeneous model implementation
//////////////////////////////////////////////////////////////////////////////////
CHeteroEMP::CHeteroEMP(CData *Data, CTree *Tree, string Name, double *S_ij, vector <vector <int> > SubTrees, vector <int> Root) : CBaseModel(Data,Tree)	{
	int i,j,k;
	vector <vector <double> > Freqs;
	vector <double> vdTemp;
	CBaseProcess *HetProc;
	// Init stuff
	m_sName = Name;

	// Start by painting the tree (Also checks validity of SubTrees)
	PaintTree(Tree,SubTrees);
	// Create the vector of frequencies
	FOR(i,(int)SubTrees.size()) {
		vdTemp.assign(20,0);
		FOR(j,(int)SubTrees[i].size())	{
			FOR(k,Data->m_iSize) {
				if(Data->m_ariSeq[SubTrees[i][j]][k] == 20) { continue; }	// Gaps don't count
				vdTemp[Data->m_ariSeq[SubTrees[i][j]][k]] += Data->m_ariPatOcc[k];
		}	}
		vdTemp = NormaliseVector(vdTemp);
		Freqs.push_back(vdTemp);
		vdTemp.clear();
	}
	Freqs.push_back(Data->m_vFreq);	// Final category has same as global frequencies

//	Freqs.clear(); FOR(i,5) { Freqs.push_back(Data->m_vFreq); }
//	cout << "\nFrequencies: "; FOR(i,5) { cout << "\nFreq["<<i<<"]: " << Freqs[i]; }


	// Build the process with the new equilibriums
	int RootFreq = 2;
	cout << "\nRoot Frequency linked to grouping["<<2<<"] = " << SubTrees[RootFreq];
	HetProc = new CHeteroEmpProc(Data,Tree,"HeteroProcess",S_ij,Freqs, RootFreq);
	m_vpProc.push_back(HetProc);
	HetProc = NULL;
	Tree->SetStartCalc(GetRoot(Tree,Root));
	SetFastBranchOpt(true);
	FinalInitialisation();
	cout << "\HetRoot set to " << Tree->StartCalc();
}

// Destructor function
CHeteroEMP::~CHeteroEMP()	{
}

//////////////////////////////////////////////////////////////////////////////////
// Model parameters for optimisation
vector <double *> CHeteroEMP::GetOptPar(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm)	{
	int i,j;
	vector <double *> RetVal = CBaseModel::GetOptPar(ExtBranch,IntBranch,Parameters,false);	// Always falsify the eqm for the het models so the correct thing is initialised

	if(Eqm) {
		cout << "\nAdding frequencies for optimisation -- this is unstable and untested...";
		FOR(i,m_vpProc[0]->NoPar())	{
			if(m_vpProc[0]->pPar(i)->Name().find("Freq_") != string::npos && !m_vpProc[0]->pPar(i)->Special()) {
				RetVal.push_back(m_vpProc[0]->pPar(i)->OptimiserValue());
				m_vpAllOptPar.push_back(m_vpProc[0]->pPar(i));
		}	}
	}
//	cout << "\nBuilt the original parameters -- size " << RetVal.size();
//	FOR(i,m_vpAllOptPar.size()) { cout << "\nPar["<<i<<"]: = " << m_vpAllOptPar[i]->Name() << " = " << *RetVal[i]; }
/*
	int ParToTest = 118;
	double OriVal = *RetVal[ParToTest];
	cout << "\nOriginal lnL for par = " << *RetVal[ParToTest] << " = " << lnL(true);
	cout << "\n------------------- Varying parameter " << ParToTest <<": " << *m_vpAllOptPar[ParToTest] << " -----------";
	for(double dodo =0.1;dodo < 1; dodo +=0.05)	{
		*RetVal[ParToTest] = dodo;
		cout << "\n" << *RetVal[ParToTest] << " = " << lnL(true);
	}
	*RetVal[ParToTest] = OriVal;
	cout << "\nOriginal lnL for par = " << *RetVal[ParToTest] << " = " << lnL(true);
	exit(-1);*/
	return RetVal;
}

/*
 * ORIGINAL STUFF
 * 	// Get values fo optimising parameters
	if(Parameters == true && !Locked())	{
		FOR(i,(int)m_vpPar.size())	{
			if(m_vpPar[i]->Opt() == false || m_vpPar[i]->Special()) { continue; }
			string::size_type loc = m_vpPar[i]->Name().find("Freq",0);
			if(Eqm == false && loc != string::npos) { continue; }
			OptVal.push_back(m_vpPar[i]->OptimiserValue());
			m_vpAllOptPar.push_back(m_vpPar[i]);
	}	}
 *
 */

//////////////////////////////////////////////////////////////////////////////////
// Heterogeneous process implementation
//////////////////////////////////////////////////////////////////////////////////
CHeteroEmpProc::CHeteroEmpProc(CData *D, CTree *T, string Name, double *S_ij, vector< vector <double> > Freq, int RootFreq) : CBaseProcess(D,T)	{
	int i,j;
	CSimpleEqm *Eqm;
	CQMat *QMat;
	string fName;
	// Initialise some basic stuff
	m_iRootEquilibrium = -1;
	MakeBasicSpace(AA);
	m_DataType = AA;
	m_sName = Name;

	// Some minor validation. That's mostly done at the model level
	assert(D != NULL); assert(D->m_DataType == AA);
	if(D->m_DataType != AA) { Error("Trying to initialise AA model with data that doesn't look like amino acids...\n\n"); }
	FOR(i,Freq.size())	{ assert(Freq[i].size() == 20 && fabs(Sum(&Freq[i]) - 1) < 1.0E-6); }	// Frequencies valid
	// Create the equilibriums and prepare the QMatrices
	FOR(i,Freq.size()) {
		// Add Eqm
		Eqm = new CSimpleEqm(20,&m_vpPar,Freq[i]);
		Eqm->SetOpt(false);
		FOR(j,20) {
			fName = Eqm->EqmPar(j)->Name();
			fName.insert(fName.find("Freq") +4 ,"_" + int_to_string(i));
			Eqm->EqmPar(j)->Name(fName);
		}
		m_vpHeteroEqm.push_back(Eqm);
		Eqm = NULL;
		// Add QMat
		QMat = new CQMat(AA,"Mat["+int_to_string(i)+"]"); QMat->InitQ(m_dBaseVal);
		m_vpHeteroQMat.push_back(QMat);
		QMat = NULL;

	}
	AssignRootEqm(RootFreq);
	// Create the rest as though it were a standard EMP model
	Add_QMat(Name,AA);
	CreateEMPmodel(S_ij,Freq[Freq.size()-1]);
	// Go through the matrices and prepare them
	FOR(i,m_vpHeteroQMat.size()) {
//		cout << "\n------------------------------------------------------------------------------------------------";
//		cout << "\nMatrix["<<i<<"]: =name= " << *m_vpHeteroQMat[i];
		MakeHetQMat(i);
	}

}

bool DoEMPoutput =false;
void CHeteroEmpProc::CreateEMPmodel(double *S_ij,vector <double> Freq)	{
	int i,j,count = 0;
	double Max = -1;
	string Name;
	CQPar *Par;
	// Correct S_ij to ensure that none are over the maximum value allowed by a parameter
	if(DoEMPoutput) { cout << "\nThe matrix used is: "; }
	FOR(i,m_iChar) {
		if(DoEMPoutput) { cout << "\n"; }
		FOR(j,i)	{
			if(DoEMPoutput) { cout << S_ij[count] << "\t"; }
			Name = GetPos(m_sABET,i,m_iABET_length/m_iChar) + " <-> " + GetPos(m_sABET,j,m_iABET_length/m_iChar);
			Par = new CQPar(Name,m_iChar,S_ij[count++],false,0.0,BIG_NUMBER);
			Par->AddQij(i,j);
			Par->SetOptimise(false);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}
	if(DoEMPoutput) { PrepareQMats(false); m_vpQMat[0]->OutQ(); }
	// Do the equilibrium distribution
	AddSimpleEqm();
	PrepareQMats();
	m_vpQMat[0]->Lock();
}

CHeteroEmpProc::~CHeteroEmpProc()	{
	int i;
	if(!m_bPseudoProcess) {
		FOR(i,(int)m_vpHeteroEqm.size()) { if(m_vpHeteroEqm[i] != NULL) { delete m_vpHeteroEqm[i]; m_vpHeteroEqm[i] = NULL; } }
		FOR(i,(int)m_vpHeteroQMat.size()) { if(m_vpHeteroQMat[i] != NULL) { delete m_vpHeteroQMat[i]; m_vpHeteroQMat[i] = NULL; } }
	}
}

// Makes the heterogeneous Q matrix
void CHeteroEmpProc::MakeHetQMat(int ProcNum)	{
	int i;
	// Some basic error checking
	assert(m_vpHeteroQMat.size() == m_vpHeteroEqm.size() && InRange(ProcNum,0,(int)m_vpHeteroQMat.size()));
	assert(!m_vpHeteroQMat[ProcNum]->IsLocked());
	// Make the process
	m_vpHeteroQMat[ProcNum]->InitQ(m_dBaseVal);
//	cout << "\nInitialised" << flush;
	FOR(i,m_vpPar.size()) {
		m_vpPar[i]->UpdatePar();
		m_vpHeteroQMat[ProcNum]->ApplyPar2Q(m_vpPar[i]);
	}
//	cout << "\nApplied "<< m_vpPar.size() << " parameters" << flush;
	m_vpHeteroEqm[ProcNum]->ApplyEqm2QMat(m_vpHeteroQMat[ProcNum]->Q(),m_vpHeteroQMat[ProcNum]->ID());
//	cout << "\nApplied Frequencies\n" << *m_vpHeteroQMat[ProcNum] << endl << flush;
	m_vpHeteroQMat[ProcNum]->DoQDiag();
	m_vpHeteroQMat[ProcNum]->Decompose(m_vpHeteroEqm[ProcNum]->Eqm());

	m_vpHeteroQMat[ProcNum]->Lock();
//	cout << "\nBuilt" << flush;
}

bool CHeteroEmpProc::Make_PT(int Branch, bool RedoRate)	{
	bool RetVal;
	Tree()->UpdateB(Branch);
	if(RedoRate) {
		m_vpHeteroQMat[Tree()->BranchLabel(Branch)]->ScaleQ(m_pRate->Val()); }

	RetVal = m_vpHeteroQMat[Tree()->BranchLabel(Branch)]->MakePT(Tree()->B(Branch),PT(Branch));
//	cout << "\nBranch["<<Branch<<"] uses Q["<<Tree()->BranchLabel(Branch)<<"]\n"; OutPT(cout,Branch);
/*
	if(Branch < 2) {
		cout << "\nMaking PT for branch[" << Branch<<"] = " << Tree()->BranchLabel(Branch) << " rate= " << m_pRate->Val() << " ; REDO = " << RedoRate <<" ; PT: ";
		int i; FOR(i,5) { cout << "  " << GetPT(Branch)[i]; }
	}
*/

	return RetVal;
}

bool CHeteroEmpProc::Likelihood(bool ForceReal)	{
	int i;
	// Rebuild the Q matrices (TODO: This needs flagging for efficiency)
	if(!m_bIsProcessCopy) {
		FOR(i,m_vpHeteroQMat.size()) {
			m_vpHeteroQMat[i]->Unlock();
			MakeHetQMat(i);
	}	}
//	FOR(i,m_vpHeteroEqm.size()) { cout << "\nEqm["<<i<<"]: " << m_vpHeteroEqm[i]->Eqm(); }
	FOR(i,m_vpHeteroQMat.size()) { m_vpHeteroQMat[i]->ScaleQ(m_pRate->Val()); }
	return CBaseProcess::Likelihood(ForceReal);
}

void CHeteroEmpProc::AssignRootEqm(int root) {
	assert(m_iRootEquilibrium == -1);
	assert(InRange(root,0,(int)m_vpHeteroEqm.size()));
	m_iRootEquilibrium = root;
}

vector <double> CHeteroEmpProc::RootEqm()	{
	assert(m_iRootEquilibrium != -1);
	return m_vpHeteroEqm[m_iRootEquilibrium]->Eqm();
}

// Functions dealing with rate process copies
CBaseProcess * CHeteroEmpProc::RateProcessCopy()	{
	CHeteroEmpProc *NewProc;
	static int CopyNum = 0;
	string Name;
	int i,count;
	vector <vector <double> > Freq;
	vector <double> storeFreq;
	double S_ij[190];
	// Do naming and initialise
	Name = "Pseudo" + m_sName + "(CopyID=" + int_to_string(CopyNum++) + ")";

	// Get the right bits and pieces for the process copy
	FOR(i,m_vpHeteroEqm.size()) { storeFreq = m_vpHeteroEqm[i]->Eqm(); Freq.push_back(storeFreq); }

	count = 0 ; FOR(i,(int)m_vpPar.size()) { if(m_vpPar[i]->Name().find("<->") != string::npos) { assert(count < 190); S_ij[count++] = m_vpPar[i]->Val(); } }
	assert(count == 190);

	NewProc = new CHeteroEmpProc(m_pData, m_pTree,Name,S_ij,Freq,m_iRootEquilibrium);
	// Fill in some stuff
	NewProc->CleanPar();
	NewProc->CleanQ();
	// Clean up the m_vpQMat and m_viQ2Bra since they need to match the original processes
	FOR(i,NewProc->m_vpQMat.size()) { delete NewProc->m_vpQMat[i]; } NewProc->m_vpQMat.clear();
	FOR(i,(int)m_vpQMat.size())	{ NewProc->m_vpQMat.push_back(m_vpQMat[i]); }
	// New stuff
	FOR(i,(int)m_vpHeteroEqm.size()) { delete NewProc->m_vpHeteroEqm[i]; NewProc->m_vpHeteroEqm[i] = m_vpHeteroEqm[i]; }
	FOR(i,(int)m_vpHeteroQMat.size()) { delete NewProc->m_vpHeteroQMat[i]; NewProc->m_vpHeteroQMat[i] = m_vpHeteroQMat[i]; }
	NewProc->m_bPseudoProcess = true;
	NewProc->m_iHiddenChar = m_iHiddenChar;
	NewProc->m_iDataChar = m_iDataChar;
	FOR(i,(int)m_vpCovProbs.size()) { NewProc->m_vpCovProbs.push_back(m_vpCovProbs[i]); }
	FOR(i,(int)m_vpEqm.size()) { NewProc->m_vpEqm.push_back(m_vpEqm[i]); }
	NewProc->m_bMaxRate = m_bMaxRate;
	NewProc->m_bIsProcessCopy = true;



	// Return it
	return NewProc;
}





