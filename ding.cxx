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

////////////////////////////////////////////////////////////////
// External variables
extern vector <STabuTree> TabuTrees;
extern vector <double> PWDists;
extern int TABU_RADIUS;
/////////////// DING FUNCTIONS //////////////////
vector <int> ReadNameFile(string NameFile, CData *D);	// Read a set of names and return a vector of indexes specifying those names
template <class TVecCon> vector <TVecCon> VecCon(vector <TVecCon> A, vector <TVecCon> B) {
	vector <TVecCon> Ret;
	Ret.reserve(A.size() + B.size());
	Ret.insert(Ret.end(),A.begin(),A.end());
	Ret.insert(Ret.end(),B.begin(),B.end());
	return Ret;
}
// Finds the branch in the tree that specifies the novel branch of that hypothesis
int FindSpecialSplit(CTree *T,vector <int> S1,vector <int> S2,vector <int> S3,vector <int> S4,vector <int> S5,vector <int> S6);
// Adds branch labels onto the tree such that 0 is the novel branch, and (1,4) are the subclades provides by (S1,S4);
void PaintTree(CTree *T, int KeyBra, vector <int> S1, vector <int> S2, vector <int> S3, vector <int> S4);


int ding()	{
		int Ret = 0;
		int GeneticCode = 0;
		int count = 0;
		// Stuff from Leaphy
		int i,j,k,l,NumModelReruns = 1;
		long RandomSeed = 0;
		string temp_string, Name, outfilestring;
		vector <string> Toks;

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
	CData DingDat("mito55_exRog.phy",AA);
	cout << "\nData read successfully ("<<DingDat.m_iNoSeq << " x " << DingDat.m_iSize << ")";
	// Name files
	vector <int> EukNames = ReadNameFile("mito55_Euk.names",&DingDat);
	vector <int> RikNames = ReadNameFile("mito55_Rik.names",&DingDat);
	vector <int> OthNames = ReadNameFile("mito55_Oth.names",&DingDat);
	vector <int> ProNames = ReadNameFile("mito55_Pro.names",&DingDat);
	assert(EukNames.size() + RikNames.size() + OthNames.size() + ProNames.size() == DingDat.m_vsName.size());
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
	// 	 Tree H1
	cout << "\nReading trees: H1" << flush;
	CTree OriT_H1("mito55_exRog_H1.tre", true, &DingDat);
	CTree T_H1 = OriT_H1;
	int H1_Split = FindSpecialSplit(&T_H1,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
	PaintTree(&T_H1, H1_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H1(&DingDat,&T_H1,"LG",false,(double*)dLGVal,(double*)dLGFreq);
//	CWAG LG_H1(&DingDat,&T_H1);
	//   Tree H2
	cout << " ... H2" << flush;
	CTree OriT_H2("mito55_exRog_H2.tre", true, &DingDat);
	CTree T_H2 = OriT_H2;
	int H2_Split = FindSpecialSplit(&T_H2,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
	PaintTree(&T_H2, H2_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H2(&DingDat,&T_H2,"LG",false,(double*)dLGVal,(double*)dLGFreq);
//	CWAG LG_H2(&DingDat,&T_H2);
	//   Tree H3
	cout << " ... H3" << flush;
	CTree OriT_H3("mito55_exRog_H3.tre", true, &DingDat);
	CTree T_H3 = OriT_H3;
	int H3_Split = FindSpecialSplit(&T_H3,EukRik,EukOth,EukPro,RikOth,RikPro,OthPro);	// The number of the unique branch
	PaintTree(&T_H3, H3_Split, EukNames, RikNames, OthNames, ProNames);
	CEMP LG_H3(&DingDat,&T_H3,"LG",false,(double*)dLGVal,(double*)dLGFreq);
//	CWAG LG_H3(&DingDat,&T_H3);
	cout << " ... done" << flush;
	cout << "\nTrees and data input" << flush;

#define EXAMINE_AA_CONTENT 0
#if EXAMINE_AA_CONTENT == 1
	cout << "\nExamining AA content";
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

	cout << "\n------------------------- H1------------------------------";
	//cout << "\n"


	cout << "\nInitial Likelihoods"; cout.precision(10);
	cout << "\nH1: " << LG_H1.lnL() << flush;
//	FullOpt(&LG_H1);
//	cout << " -> " << LG_H1.lnL();
	cout << "\nH2: " << LG_H2.lnL() << flush;
	cout << "\nH3: " << LG_H3.lnL() << flush;

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


}
