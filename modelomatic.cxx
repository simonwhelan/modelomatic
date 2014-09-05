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

#include "./modelomatic.h"
#include <time.h>
#include "ini/cpp/INIReader.h"
#include <set>

#define CHECK_LNL_OUT 1

#define DEVELOPER_VERSION_MAIN 1
#if DEVELOPER_VERSION_MAIN == 1
#define VERSION_NUMBER "1.1 (development)"
#else
#define VERSION_NUMBER "1.0 (release)"
#endif

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter, MakeQ_Log_Counter, MakePT_Log_Counter, LFunc_Log_Counter, SubLFunc_Log_Counter, SPR_Log_Counter;
#endif
vector <double> GetPWVar(CBaseModel *Model, vector <double> *CurrentPW);
double GetDistVar(CBaseModel *M, CData *D,int Seq1, int Seq2, double CurrentDist);
extern bool ALLOW_PREDICTLNL;

/////////////////////////////////////////////////////////////////////////////////
// Main tree estimation routine
SBestModel DoTreeEstimation(CBaseModel *M, CData *D, CTree *T);

// Pointer to likelihood function that changes depending on the type of calculation performed
vector <STabuTree> TabuTrees;
vector <double> PWDists;
CPhyloDat PhyDat;
int TABU_RADIUS = DEFAULT_TABU_RADIUS, EXIT_OBS, OptObs;
double PROB_RAN_SEQ_REM = DEFAULT_PROB_RAN_SEQ_REM;
bool AllowPreOpt = true, DoItFast = true, DoItTrim = false, ModelOut = false;
int TrimTree = 10;		// Number of sequences defaulted to by the DoItTrim option
void DoInstructions();
bool RunHadError = false;

bool WarningMulD;

// File specific global variables
int DebugOutput = 0;

int n_models=13;
string model_names[] = {"JTT","WAG","LG","DAY","mtREV","mtMam","mtArt","rtREV","cpREV","BLOSUM62","VT","HIVb","HIVw"};
double* smat[] = {(double*)dJTTVal,(double*)dWAGVal,(double*)dLGVal,(double*)dDAYVal,(double*)dmtREVVal,(double*)dmtMAMVal,(double*)dmtArtVal,(double*)drtREVVal,(double*)dcpREVVal,(double*)dBLOSUM62Val,(double*)dVTVal,(double*)dHIVbVal,(double*)dHIVwVal};
double* freq[] = {(double*)dJTTFreq,(double*)dWAGFreq,(double*)dLGFreq,(double*)dDAYFreq,(double*)dmtREVFreq,(double*)dmtMAMFreq,(double*)dmtArtFreq,(double*)drtREVFreq,(double*)dcpREVFreq,(double*)dBLOSUM62Freq,(double*)dVTFreq,(double*)dHIVbFreq,(double*)dHIVwFreq};

std::map <string,double**> aa_model_map;
std::set <string> codon_model_set;

int main(int argc, char *argv[])	{
		int GeneticCode = 0;
		int count = 0;
		int RY_count, DNA_count, AA_count, COD_count;
		vector <SModelDetails> Models;

		// Stuff from Leaphy
		ALLOW_PREDICTLNL = false;
		int i,j,k,NumModelReruns = 1;
		long RandomSeed = 0;
		string temp_string, Name, outfilestring;
		vector <string> Toks;
		TabuTrees.clear();
		WarningMulD = false;
		SBestModel Result;
		bool DoSurface = false, DoHardOptCheck = false, bDoBioNJ = true;
		string TreeName, DataName;	// Some files strings for Trim option
		ofstream *out = NULL;		// Model output pointer

	string InTree = "file=";

	// Some initial verification
	if(min(LOOSE_RMSD_SUBSET,FULL_RMSD_SUBSET) < 105 && min(LOOSE_PARS_SUBSET,FULL_PARS_SUBSET) < 105) { Error("\nTrying to choose subsets of trees to examine based both on RMSD and parsimony"); }

	// Get information
	if(!InRange(argc,2,8)) {
		cout << "ModelOMatic <data_file> <tree> <output_file> <genetic_code> <normal/fast/trim> <models_out>\n";
		cout << "\n---";
		cout << "\n\t<data_file>:  \tInput data in sequential or interleaved format";
		cout << "\n\t<tree_file>:  \tEither input tree file in Newick or use 'bionj' to build a distance tree";
		cout << "\n\t<output_file> \tWhere the output from the program will go";
		cout << "\n\t<genetic_code>\tThe genetic code used for codon models and amino acid models [default = Universal]";
		cout << "\n\t<normal/fast/trim>        \tOption controlling how analyses will be done [default = normal]";
		cout << "\n\t\t\t\t\tnormal = full ML estimation for each tree";
		cout << "\n\t\t\t\t\tfast = Only do full MLE for first tree of each data type; after that just model parameters";
		cout << "\n\t\t\t\t\ttrim = Only use ten species with greatest tree coverage to perform analysis";
		cout << "\n\t\t\t\t\ttrim=5 = Only use 5 species with the greatest tree coverage to perform the analysis. NB: 5 can vary";
		cout << "\n\t<models_out>  \tyes/no Option to output full model details to file model.out [default = no]";
		cout << "\nExiting...\n";
		exit(-1);
	}
	cout << "\n---------------------------------------------------------------------------------------------";
	cout << "\n  ModelOMatic (v" << VERSION_NUMBER << ").\n\tA program for choosing substitutions models for phylogenetic inference.\n\tWritten by Simon Whelan.\n\tContributions from James Allen, Ben Blackburne and David Talavera.";
	cout << "\n---------------------------------------------------------------------------------------------";

	///////////////////////////////////////////////////////////////////////////////////////////
	// Create the data structures
	// 0. Input the raw data
	cout << "\nData: <" << argv[1] << ">" << flush;
	PhyDat.SetIn(argv[1]); PhyDat.GetData();
	assert(PhyDat.pData()->m_DataType == DNA);
	cout << ": " << PhyDat.pData()->m_iNoSeq << " sequences of length " << PhyDat.pData()->m_iTrueSize << " (DataMatrix: " << PhyDat.pData()->m_iNoSeq << " x " << PhyDat.pData()->m_iSize << ")" << flush;
	PhyDat.pData()->CleanToDNACodon();	// Make sure gaps are codon compatible early, otherwise it causes problems with Trim and other functions
	PhyDat.pData()->RemoveSparseSeqs(true,NULL);
	// 3. Set genetic code is done first so translation will work for bionj tree
	if(argc>4) {
		assert(InRange(atoi(argv[4]),0,NumGenCode));
		GeneticCode = atoi(argv[4]);
	}
	// 1. Create a tree (default: bionj)
	cout << "\nCreating start tree ... " << flush;
	bDoBioNJ = true;
	if(argc > 2) { if(strcmp(argv[2],"bionj")) { bDoBioNJ = false; } }
	CTree Tree;
	clock_t start,end;
    start = clock();

	if(bDoBioNJ) {
		// Clean data as required
		// PhyDat.pData()->RemoveSparseSeqs(true,NULL); //PhyDat.pData()->CondenseGaps();
		// Create a bionj starting tree
        CData tmp_AA_Data = *PhyDat.pData();
        tmp_AA_Data.Translate(GeneticCode);
		CEQU EQU_PW(&tmp_AA_Data,NULL);
		PWDists = GetPW(&EQU_PW,NULL,true);  // Get pairwise distances
		if(PhyDat.pData()->m_iNoSeq > 2) {
			CTree T_bionj(DoBioNJ(PWDists,PhyDat.pData()->m_vsName,true),tmp_AA_Data.m_iNoSeq);
			Tree = T_bionj;
			cout << " estimated using bionj" << flush;
		} else {
			temp_string ="(" + PhyDat.pData()->m_vsName[0] + ":" + double_to_string(min(PWDists[1],0.5)) + "," + PhyDat.pData()->m_vsName[1] + ":0.0);";
			CTree Pair_tree(temp_string,2,false,PhyDat.pData());
			Tree = Pair_tree;
			cout << " no tree for a pair" << flush;
		}


	} else {
		// Take tree from file
		InTree += argv[2];
		PhyDat.SetStartTree(InTree);
		if(PhyDat.pTree() == NULL) { Error("Failed to read tree from file <" + (string) argv[2] + ">\n\n"); }
		PhyDat.pTree()->Unroot();
		Tree = *PhyDat.pTree();
		cout << " taken from <" << InTree <<"> successfully" << flush;
		//PhyDat.pData()->CondenseGaps();
	}
    end = clock();
    cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)"<<flush;
    start=clock();
	// 2. Set output
	if(argc > 3) {
		outfilestring = argv[3];
		PhyDat.SetOut(argv[3]);
	} else {
		outfilestring = argv[1]; outfilestring += ".output";
		cout << "\nTrying to work with: " << outfilestring;
		PhyDat.SetOut(outfilestring.c_str());
	}
	// 4. If needed do the DoItTrim option
   	// Check whether if is meant to be running fast
    if(argc>5) {
    	TreeName = argv[5];
    	// Organize the optimisation settings
    	if(TreeName.find("normal") != string::npos) { DoItFast = false; TreeName = find_and_replace(TreeName,"normal",""); }
    	else if(TreeName.find("fast") != string::npos) { DoItFast = true; TreeName = find_and_replace(TreeName,"fast",""); }
    	// Do the trim settings
    	if(TreeName.find("trim") != string::npos) {
    		DoItTrim = true; TreeName = find_and_replace(TreeName,"trim","");
    		if(!TreeName.empty()) {
    			if(TreeName[0] == '=') { // Get the trim length
    				TreeName.erase(0,1);
    				TrimTree = atoi(TreeName.c_str());
    				if(TreeName.empty()) { cout << "\nError with number of trimmed sequences. Needs to be set through (e.g.) trim=5. You have set: trim= without a number"; exit(-1); }
    				if(!InRange(TrimTree,0,PhyDat.pData()->m_iNoSeq)) { cout << "\nTrying to trim tree down to " << TrimTree << " sequences (from string: " << TreeName << "), but there are only " << PhyDat.pData()->m_iNoSeq << " sequences...\n"; exit(-1); }
    				FOR(i,TreeName.size()) {
    					if(!isdigit(TreeName[i])) { break; }
    				}
    				TreeName.erase(0,i);
    			}
    		}
    	}
    	// There should be nothing left in the string
   		if(!TreeName.empty()) { cout << "\nParsing of trim/fast option has resulted in left over characters: " << TreeName << "\n\n"; exit(-1); }
   	}
    cout << "\nOptimisation settings: "; if(DoItFast) { cout << "fast "; } else { cout << "normal "; } if(DoItTrim) { cout << "trim=" << TrimTree; }
    // 5. Check model output options
    if(argc>6) {
    	TreeName = argv[6];
    	if(TreeName.find("yes") != string::npos) { ModelOut = true; }
    	else if(TreeName.find("no") != string::npos) { ModelOut = false; }
    	else { cout << "\nUnknown option for models_out: " << TreeName << "\n\n"; exit(-1); }
    }
    // 6. Some file validation
    if(FileExist(outfilestring)) { cout << "\nError: main output file <" << outfilestring << "> already exists. Please delete this and related files before continuing.\n"; exit(-1); }
    if(ModelOut) { if(FileExist((string)outfilestring + (string)".model.out")) { cout << "\nError: models_out file <" << outfilestring + (string)".model.out" << "> already exists. Please delete this and related files before continuing.\n"; exit(-1); } }
    // 7. Do trimming if required
	if(DoItTrim && Tree.NoSeq() > TrimTree)	{
		bool Check;
		// Build file names
		TreeName = outfilestring + ".trim.tree";
		DataName = outfilestring + ".trim.data";
		if(FileExist(TreeName)) { cout << "\n\tError: <" << TreeName << "> exists already. Please delete before continuing.\n"; exit(-1); }
		if(FileExist(DataName)) { cout << "\n\tError: <" << DataName << "> exists already. Please delete before continuing.\n"; exit(-1); }
		cout << "\nTRIMMING: \tTree contains " << Tree.NoSeq() << " sequences (>TrimTree=" << TrimTree << ") ..." << flush;
		// Very briefly optimise the tree under JC
		cout << "\n          \tLoose optimisation of start tree under JC ..." << flush;
		CJC *TrimJC = NULL; TrimJC = new CJC(PhyDat.pData(),&Tree);
		TrimJC->lnL(); LazyOpt(TrimJC,false,true,false);
		delete TrimJC;
		cout << " done" << flush;

		// Now build the greedy tree
		cout << "\n          \tObtaining greedy start tree with " << TrimTree << " sequences ..." << flush;
		Tree  = FindGreedySubTree(&Tree,TrimTree);
		cout << " done" << flush;
		// Output tree to appropriate place
		Tree.SetNames(PhyDat.pData()->m_vsName,true);
		Tree.OutName(); Tree.OutBra();
		ofstream trout(TreeName.c_str());
		trout << Tree << flush;
		trout.close();
		// Output the new data
		// Needs to be slightly clever so that columns of '-' are not included
		vector <string> NewSeq(TrimTree,"");
		vector <int> Seqs2Out(TrimTree,-1), temp;
		Tree.BuildSplits();
		FOR(i,Tree.NoSeq()) { if(Tree.BraLink(i,0) != -1) { break; } }
		assert(i!=Tree.NoSeq());
		temp = Tree.GetSplit(i).Right;
		Seqs2Out = Tree.GetSplit(i).Left; Seqs2Out.insert(Seqs2Out.end(),temp.begin(),temp.end());
		assert((int)Seqs2Out.size() == TrimTree);
		FOR(j,(int)PhyDat.pData()->m_iTrueSize) {
			FOR(i,TrimTree) {
				if(GAP_ABET.find(PhyDat.pData()->m_vsTrueSeq[Seqs2Out[i]][j]) == string::npos) { break; } // Real sequence exists so output
			}
			if(i != TrimTree) { // Sequence exists at that site so add to data for output
				FOR(i,TrimTree) {
					NewSeq[i].push_back(PhyDat.pData()->m_vsTrueSeq[Seqs2Out[i]][j]);
		}	}	}
		ofstream datout(DataName.c_str());
		datout << TrimTree << "  " << NewSeq[0].size() << endl << endl;
		FOR(i,TrimTree) { datout << PhyDat.pData()->m_vsName[Seqs2Out[i]] << " \t" << NewSeq[i] << "\n\n" << flush; }
		datout.close();
		// Read back the data
		cout << "\n          \tReinitialising objects for trimmed data ..." << flush;
		PhyDat.CleanData();	PhyDat.SetIn(DataName); PhyDat.GetData();
		//CData NewData(DataName,DNA); PhyDat.
		// Finally read back the tree
		CTree SmallTree(TreeName,true,PhyDat.pData()); Tree = SmallTree;
		cout << " done" << flush;
		cout << "\n          \tNew files available: data = <" << DataName << ">; tree = <" << TreeName << ">";
		// Create
	}
	// 8. Check modelomatic.ini to decide what models to run
	cout << "\nScanning for modelomatic.ini file ...";
	GetModels();
	cout << " done";
	// 9. Create the other data sets
	CData NT_Data = *PhyDat.pData();
	CData AA_Data = *PhyDat.pData();
	CData AA_Temp = *PhyDat.pData();  AA_Temp.Translate(GeneticCode);	// Error check the translation for stop codons and so on
	CData COD_Data = *PhyDat.pData(); //COD_Data.MakeCodonData();
	CData RY_Data = *PhyDat.pData(); //RY_Data.DNA2RY();

//	cout << "\nRY:  \t" << RY_Data.CountMSAChars();
//	cout << "\nDNA: \t" << PhyDat.pData()->CountMSAChars();
//	cout << "\nAA:  \t" << AA_Temp.CountMSAChars();
//	cout << "\nCoD:  \t" << COD_Data.CountMSAChars();
	// Set output
	PhyDat.SetOut(outfilestring);

	cout << "\nWorking with genetic code: " << GenCodeName[GeneticCode];
	if(PhyDat.pData()->m_iNoSeq == 2) { cout << "\nWorking with 2 sequences so cannot apply gamma distributed rates-across-sites"; }
	cout << "\n>>> Doing model analysis <<< \n" << flush;

	///////////////////////////////////////////////////////////////////////////////////////////
	// Development stuff goes here
#if DEVELOPER_VERSION_MAIN
	CBaseModel *ModelPointer = NULL;
	double cVal,rVal;
	cout << "\nDoing codon based analysis"; cout.precision(10); // exit(-1);

	// Some general settings
	bool OptM0 = true;
	bool OptDrDc = true;
	bool Opt1Dr2Dc = true;
	bool Opt2Dr1Dc = true;
	bool Simulation = true;
	// Some optimiser stuff
	bool DoBra = true;
	bool DoPar = true;

	vector <int> RadMat(20*20,-1);
	FINOPEN(Radin, sRadicalFileName.c_str());
	FOR(i,20)	{
		FOR(j,i)	{
			Radin >> RadMat[(i*20)+j];
			if(!InRange(RadMat[(i*20)+j],0,2)) { cout << "\nError reading Radical Matrix from: " << sRadicalFileName << " at ["<<i << "," << j << "] = " << RadMat[(i*20)+j] << "\nMatrix so far: " << MatOut(20,RadMat); exit(-1); }
			RadMat[(j*20)+i] = RadMat[(i*20)+j];
		}
	}
	Radin.close();
//		cout << "\nRadical Matrix" << endl <<  MatOut(20, RadMat);
//		cout << "\n\nDone";

	CData NT1 = *PhyDat.pData(); NT1.GetCodonPositions(true,false,false);
	CData NT2 = *PhyDat.pData(); NT2.GetCodonPositions(false,true,false);
	CData NT3 = *PhyDat.pData(); NT3.GetCodonPositions(false,false,true);
/*	CData NT12 = *PhyDat.pData(); NT12.GetCodonPositions(true,true,false);
	CData NT13 = *PhyDat.pData(); NT13.GetCodonPositions(true,false,true);
	CData NT23 = *PhyDat.pData(); NT23.GetCodonPositions(false,true,true);
	CData NT123 = *PhyDat.pData(); NT123.GetCodonPositions(true,true,true);*/

//	int ShowSeq = RandInt(0,NT1.m_iNoSeq-1);
//	cout << "\nOriginal data:		   " << PhyDat.pData()->m_iNoSeq << " " << PhyDat.pData()->m_iTrueSize << "\t" << PhyDat.pData()->m_vsTrueSeq[ShowSeq].substr(0,15);;


	/*
	// Checking rate idea...
	cout << "\n\nParameters:";
	FOR(i,M0Test->m_vpProc[0]->NoPar()) {
		cout << "\nPar[" << i<<"]: " << *M0Test->m_vpProc[0]->pPar(i);
	}
	M0Test->lnL(true);
	cout << "\nALPHABET\n";
	FOR(i,64) {
		cout << i << ":" << State(COD,i) << "[" << GenCodes[0][i] << "] ";
	}
	cout << "\nEqm: " << M0Test->m_vpProc[0]->Eqm(0);
	double SetOmegaVal;

	M0Test->m_vpProc[0]->pPar(14)->SetVal(2.5);

	// Omega = 0.1
	SetOmegaVal = 0.1;
	M0Test->m_vpProc[0]->pPar(13)->SetVal(SetOmegaVal);
	cout << "\nHave Omega = " << SetOmegaVal << " == " << *M0Test->m_vpProc[0]->pPar(13) << " ; lnL = " << M0Test->lnL(true);
	// Observed rates
	cout << "\nObsSynRate:    " << M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nObsNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(true);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(true) / M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nTotal rate:    "<< M0Test->m_vpProc[0]->NonsynRate(true) + M0Test->m_vpProc[0]->SynRate(true);
	// Actual rates
/*	cout << "\nSynRate:    " << M0Test->m_vpProc[0]->SynRate(false);
	cout << "\nNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(false);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(false) / M0Test->m_vpProc[0]->SynRate(false); */
//	cout << "\n\nQMat: "; M0Test->m_vpProc[0]->OutQ();
//	cout << "\n\nEqm: " << M0Test->m_vpProc[0]->Eqm(0);
//	cout << "\nModel: " << *M0Test;
/*
	// Omega = 0.5
	SetOmegaVal = 0.5;
	M0Test->m_vpProc[0]->pPar(13)->SetVal(SetOmegaVal);
	cout << "\nHave Omega = " << SetOmegaVal << " == " << *M0Test->m_vpProc[0]->pPar(13) << " ; lnL = " << M0Test->lnL(true);
	// Observed rates
	cout << "\nObsSynRate:    " << M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nObsNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(true);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(true) / M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nTotal rate:    "<< M0Test->m_vpProc[0]->NonsynRate(true) + M0Test->m_vpProc[0]->SynRate(true);
	// Actual rates
	cout << "\nSynRate:    " << M0Test->m_vpProc[0]->SynRate(false);
	cout << "\nNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(false);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(false) / M0Test->m_vpProc[0]->SynRate(false);

	// Omega = 1.0
	SetOmegaVal = 1.0;
	M0Test->m_vpProc[0]->pPar(13)->SetVal(SetOmegaVal);
	cout << "\nHave Omega = " << SetOmegaVal << " == " << *M0Test->m_vpProc[0]->pPar(13) << " ; lnL = " << M0Test->lnL(true);
	// Observed rates
	cout << "\nObsSynRate:    " << M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nObsNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(true);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(true) / M0Test->m_vpProc[0]->SynRate(true);
	cout << "\nTotal rate:    "<< M0Test->m_vpProc[0]->NonsynRate(true) + M0Test->m_vpProc[0]->SynRate(true);
	// Actual rates
	cout << "\nSynRate:    " << M0Test->m_vpProc[0]->SynRate(false);
	cout << "\nNonsynRate: " << M0Test->m_vpProc[0]->NonsynRate(false);
	cout << "\nRatio:         " << M0Test->m_vpProc[0]->NonsynRate(false) / M0Test->m_vpProc[0]->SynRate(false);
*/
	double testval = 0.0;

	// ************* Do Basic M0 Model *******************
	cout << "\nCreating M0";
	CCodonM0 *M0calc = NULL;
	CData Cod0 = *PhyDat.pData();
	M0calc = new CCodonM0(&Cod0,&Tree);
	ModelPointer = M0calc;
	testval = ModelPointer->lnL(true);
	if(OptM0) {  FullOpt(ModelPointer,DoPar,DoBra,false,-BIG_NUMBER,true,100,-BIG_NUMBER,FULL_LIK_ACC,true); } else { ModelPointer->lnL(true); }
	if(OptM0)	{
	cout << "\n>>>>>>>>>>>>>> FINAL DETAILS " << flush;
	cout << "\n" << *ModelPointer;
	FOR(i,ModelPointer->m_vpProc.size())	{
		cVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 0);		// Conservative
		rVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 1);		// Radical
		cout << "\n\tProc["<<i<<"] expectedObservations:\tConservative: " << cVal << "\tRadical: " << rVal << "\tDr/Dc: " << rVal/cVal;
	} 	}
	ModelPointer = NULL;



	// ************* Do Basic DrDc Model *******************
	cout << "\n>>>>>>>>>>>>>> Created DrDc model" << flush;
	CCodonDrDc *M0New = NULL;
	CData Cod2 = *PhyDat.pData();
	M0New = new CCodonDrDc(&Cod2,&Tree);
	ModelPointer = M0New;
	if(OptDrDc) {  FullOpt(ModelPointer,DoPar,DoBra,false,-BIG_NUMBER,true,100,-BIG_NUMBER,FULL_LIK_ACC,true); } else { ModelPointer->lnL(true); }
	if(OptDrDc) { cout << "\n>>>>>>>>>>>>>> FINAL DETAILS " << flush;
	cout << "\n" << *ModelPointer;
	FOR(i,ModelPointer->m_vpProc.size())	{
		cVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 0);		// Conservative
		rVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 1);		// Radical
		cout << "\n\tProc["<<i<<"] expectedObservations:\tConservative: " << cVal << "\tRadical: " << rVal << "\tDr/Dc: " << rVal/cVal;
	} 	}
	ModelPointer = NULL;

	if(Simulation)	{
	cout << "\n>>>>>>>>>>>>> TRYING SIMULATION" << flush;
int CompCount;
const int MIXTURECAT = 4;
// Type 1
//const double Dr[] = {0.2, 0.6, 0.2, 0.6};
//const double Dc[] = {0.4, 0.4, 0.8, 0.8};
// Type 2
const double Dr[] = {0.6, 0.6, 0.2, 0.6};
const double Dc[] = {0.4, 0.8, 0.8, 0.8};

const int SimSize[] = {200, 300, 0, 0};

//	CTree SimTree("(1:0.1,3:0.1,(((14:0.1,(11:0.1,15:0.1):0.1):0.1,(2:0.1,10:0.1):0.1):0.1,(12:0.1,((5:0.1,((8:0.1,16:0.1):0.1,(4:0.1,13:0.1):0.1):0.1):0.1,(6:0.1,(7:0.1,9:0.1):0.1):0.1):0.1):0.1):0.1);",16);
	CTree SimTree("(1:0.3,3:0.3,(((14:0.3,(11:0.3,15:0.3):0.3):0.3,(2:0.3,10:0.3):0.3):0.3,(12:0.3,((5:0.3,((8:0.3,16:0.3):0.3,(4:0.3,13:0.3):0.3):0.3):0.3,(6:0.3,(7:0.3,9:0.3):0.3):0.3):0.3):0.3):0.3);",16);
	CTree TreeStore; TreeStore = Tree;
	*M0New->m_pTree = SimTree;
	int NoSeq = 16, Size = 0;
	vector <int> Blob;
	vector <vector <int> > NewSeq, FullSeq;
	vector <double> Probs, Rates, FullRates; // Rates = RawProcessRates; FullRates = Those used for scaling
	double SumRate;

	// Get memory sorted
	FOR(i,MIXTURECAT) { Size += SimSize[i]; }
	Blob.reserve(Size);
	FOR(i,NoSeq) {
		FullSeq.push_back(Blob);
	}
	Blob.clear();

	// Create the relative probabilities & scales
	FOR(CompCount,MIXTURECAT)	{
		Probs.push_back((double)SimSize[CompCount]/(double)Size);
	}
	cout << "\nDoing simulation of length " << Size;
	FOR(CompCount,MIXTURECAT)	{
		cout << "\n\tCat[" << CompCount << "]: P(" << CompCount<<") =  " << Probs[CompCount] <<  "; Dc = " << Dc[CompCount] << " & Dr = " << Dr[CompCount];
		assert(NewSeq.empty());
		M0New->m_vpProc[0]->GetPar("Omega_Conservative")->SetVal(Dc[CompCount]);
		M0New->m_vpProc[0]->GetPar("Omega_Radical")->SetVal(Dr[CompCount]);
		M0New->m_vpProc[0]->GetPar("Kappa")->SetVal(2.5);

		// Copied from DrDc mixer
		// ---
		// Go through and apply all global parameters
		FOR(i,(int)M0New->NoPar()) { M0New->m_vpPar[i]->GlobalApply(); }
		// Prepare for likelihood computations
		M0New->m_vpProc[0]->PrepareLikelihood(true,true,false);
//		cout << "\n==================== SIMULATION MATRIX " << CompCount << "=================";
//		M0New->m_vpProc[0]->OutQ(0);
//		cout << "\n//";
//		cout << "\nEqm: " << M0New->m_vpProc[0]->Eqm(0);
		// Get the current rate of the individual processes (they currently scaled correctly relative to one another, just not the right mean rate)
		Rates.push_back(M0New->m_vpProc[0]->CalcRate());
		FullRates.push_back(Rates[CompCount] * Probs[CompCount]);

	}

	SumRate = Sum(&FullRates);
	cout << "\nProbs:     " << Probs;
//	cout << "\nRates:     " << Rates;
	FOR(i,MIXTURECAT) { FullRates[i] = Rates[i] / SumRate; }
	cout << "\nRates: " << FullRates;

	SumRate = 0.0; FOR(i,MIXTURECAT) { SumRate += Probs[i] * FullRates[i]; }
	cout << " == TotalRate: " << SumRate;

	// Do the actual simulation

	FOR(CompCount,MIXTURECAT)	{

		assert(NewSeq.empty());
		*M0New->m_pTree = SimTree;
//		cout << "\nBranchesProc[CompCount]: ";
		FOR(i,M0New->m_pTree->NoBra()) {
			M0New->m_pTree->SetB(i,M0New->m_pTree->B(i) * FullRates[CompCount]);
//			if(i <4) { cout << "\t" << SimTree.B(i) << "->" << M0New->m_pTree->B(i); }
		}
		M0New->m_vpProc[0]->GetPar("Omega_Conservative")->SetVal(Dc[CompCount]);
		M0New->m_vpProc[0]->GetPar("Omega_Radical")->SetVal(Dr[CompCount]);
		M0New->m_vpProc[0]->GetPar("Kappa")->SetVal(2.5);
//		cout << "\nSimulating from model for size: " << SimSize[CompCount] <<   flush; // *M0New;
		if(SimSize[CompCount]> 0) {
			M0New->DoSimulation(NoSeq,SimSize[CompCount],&NewSeq,false);
//			cout << "\nSimulation done... " << flush;
			FOR(i,NoSeq) {
				FullSeq[i].insert(FullSeq[i].end(),NewSeq[i].begin(),NewSeq[i].end());
			}
			// Output the individual sim for a category
			string OutName = "sim.data." + int_to_string(CompCount+1);
			ofstream chouder(OutName.c_str());
			chouder << NoSeq << "  " << SimSize[CompCount];
			FOR(i,NoSeq) {
				chouder << "\n\nSeq_" << i << "\t ";
				FOR(j,SimSize[CompCount]) {
					assert(InRange(NewSeq[i][j],0,62));
					if(NewSeq[i][j] == 61) { chouder << "---"; }
					else {
						FOR(k,3) { chouder << Cod2.m_sABET[(NewSeq[i][j]*3)+k]; }
					}
				}
			}
			chouder << "\n";
			chouder.close();

			// Clear the sequences
			NewSeq.clear();
		}
	}
	cout << "\n>>>>>>>>>>>>> Simulation done";
//	cout << "\nHave data " << Cod2.m_iChar;
//	cout << "\nABET[" << Cod2.m_sABET.size() <<" ]: " << Cod2.m_sABET;
	ofstream simout("sim.data");
//	cout << "\nHave " << FullSeq.size() << " of sizes"; FOR(i,NoSeq) { cout << "\n[" << i << "] " << FullSeq[i].size() << " " << FullSeq[i]; }
	assert(FullSeq[0].size() == Size);
	simout << NoSeq << "  " << Size;
	FOR(i,NoSeq) {
		simout << "\n\nSeq_" << i << "\t ";
		FOR(j,Size) {
			assert(InRange(FullSeq[i][j],0,62));
			if(FullSeq[i][j] == 61) { simout << "---"; }
			else {
				FOR(k,3) { simout << Cod2.m_sABET[(FullSeq[i][j]*3)+k]; }
			}
		}
	}
	simout << "\n";
	simout.close();
	Tree = TreeStore;
	}
//	exit(-1);	// Set here for testing...

	//////////////////////////////////////////
	// Mixture models

	// 2Dr1Dc
	CCodon2Dr1Dc *M0NewTester = NULL;
	CData Cod10 = *PhyDat.pData();
	M0NewTester = new CCodon2Dr1Dc(&Cod10,&Tree);
	ModelPointer = M0NewTester;
	cout << "\n>>>>>>>>>>>>>> Created new model" << flush;

//	cout << "\n>>>>>>>>>>>>>> Likelihood computation" << flush;

	Tree.OutBra();
	cout << "\nTree: " << Tree;

	testval = M0NewTester->lnL(true);
	if(Opt2Dr1Dc) {  FullOpt(ModelPointer,DoPar,DoBra,false,-BIG_NUMBER,true,100,-BIG_NUMBER,FULL_LIK_ACC,true); } else { ModelPointer->lnL(true); }

	if(Opt2Dr1Dc) { cout << "\n>>>>>>>>>>>>>> FINAL DETAILS " << flush;

	cout << "\n" << *ModelPointer;
	FOR(i,ModelPointer->m_vpProc.size())	{
		cVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 0);		// Conservative
		rVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 1);		// Radical
		cout << "\n\tProc["<<i<<"] expectedObservations:\tConservative: " << cVal << "\tRadical: " << rVal << "\tDr/Dc: " << rVal/cVal;
	}	}
	ModelPointer = NULL;

	// 1Dr2Dc
	CCodon1Dr2Dc *M0NewTester2 = NULL;
	CData Cod11 = *PhyDat.pData();
	M0NewTester2 = new CCodon1Dr2Dc(&Cod11,&Tree);
	ModelPointer = M0NewTester2;
	cout << "\n>>>>>>>>>>>>>> Created new model" << flush;

//	cout << "\n>>>>>>>>>>>>>> Likelihood computation" << flush;

	testval = ModelPointer->lnL(true);

	if(Opt1Dr2Dc) {  FullOpt(ModelPointer,DoPar,DoBra,false,-BIG_NUMBER,true,200,-BIG_NUMBER,FULL_LIK_ACC,true); } else { ModelPointer->lnL(true); }

	if(Opt1Dr2Dc) { cout << "\n>>>>>>>>>>>>>> FINAL DETAILS " << flush;
	cout << "\n" << *ModelPointer;
	FOR(i,ModelPointer->m_vpProc.size())	{
		cVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 0);		// Conservative
		rVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 1);		// Radical
		cout << "\n\tProc["<<i<<"] expectedObservations:\tConservative: " << cVal << "\tRadical: " << rVal << "\tDr/Dc: " << rVal/cVal;
	}	}
	ModelPointer = NULL;


	exit(-1);	// Set here for testing...

	// 2Dr2Dc
	CCodon2Dr2Dc *M0NewTester3 = NULL;
	CData Cod12 = *PhyDat.pData();
	M0NewTester3 = new CCodon2Dr2Dc(&Cod12,&Tree);
	ModelPointer = M0NewTester3;
	cout << "\n>>>>>>>>>>>>>> Created new model" << flush;

//	cout << "\n>>>>>>>>>>>>>> Likelihood computation" << flush;

	testval = ModelPointer->lnL(true);
	FullOpt(ModelPointer,true,false,false,-BIG_NUMBER,true,50,-BIG_NUMBER,FULL_LIK_ACC,true);

	cout << "\n>>>>>>>>>>>>>> FINAL DETAILS " << flush;
	cout << "\n" << *ModelPointer;
	FOR(i,ModelPointer->m_vpProc.size())	{
		cVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 0);		// Conservative
		rVal = GetAminoAcidCountFromCodon( ModelPointer->m_vpProc[i]->GetQMat(0), 0, RadMat, 1);		// Radical
		cout << "\n\tProc["<<i<<"] expectedObservations:\tConservative: " << cVal << "\tRadical: " << rVal << "\tDr/Dc: " << rVal/cVal;
	}
	ModelPointer = NULL;


	exit(-1);


	exit(-1);

	CREV *RevTest = NULL;
	double Lsum = 0.0;
	RevTest = new CREV(&NT1,&Tree);
	RevTest->lnL(true);
	cout << "\n-------------------------------------------- Individual codon positions -----------------------------";
	FullOpt(RevTest); cout <<	"\nPos1 " << FullOpt(RevTest);  cout << " == " << RevTest->lnL(true); Lsum += RevTest->lnL(true);
//	cout << "\n" << *RevTest;
	delete RevTest;

	RevTest = new CREV(&NT2,&Tree);
	RevTest->lnL(true);
	FullOpt(RevTest); cout <<	"\nPos2 " << FullOpt(RevTest);  cout << " == " << RevTest->lnL(true); Lsum += RevTest->lnL(true);
//	cout << "\n" << *RevTest;
	delete RevTest;

	RevTest = new CREV(&NT3,&Tree);
	RevTest->lnL(true);
	FullOpt(RevTest); cout <<	"\nPos3 " << FullOpt(RevTest);  cout << " == " << RevTest->lnL(true); Lsum += RevTest->lnL(true);
//	cout << "\n" << *RevTest;
	delete RevTest;

	cout << "\nFull " << Lsum;


//	exit(-1);
/*
	cout << "\nNT1    				" << NT1.m_iNoSeq << " " << NT1.m_iTrueSize << "\t" << NT1.m_vsTrueSeq[ShowSeq].substr(0,5);
	cout << "\nNT2    				" << NT2.m_iNoSeq << " " << NT2.m_iTrueSize << "\t" << NT2.m_vsTrueSeq[ShowSeq].substr(0,5);;
	cout << "\nNT3    				" << NT3.m_iNoSeq << " " << NT3.m_iTrueSize << "\t" << NT3.m_vsTrueSeq[ShowSeq].substr(0,5);;
	cout << "\nNT12    			" << NT12.m_iNoSeq << " " << NT12.m_iTrueSize << "\t" << NT12.m_vsTrueSeq[ShowSeq].substr(0,10);;
	cout << "\nNT13    			" << NT13.m_iNoSeq << " " << NT13.m_iTrueSize << "\t" << NT13.m_vsTrueSeq[ShowSeq].substr(0,10);;
	cout << "\nNT23    			" << NT23.m_iNoSeq << " " << NT23.m_iTrueSize << "\t" << NT23.m_vsTrueSeq[ShowSeq].substr(0,10);;
	cout << "\nNT123   			" << NT123.m_iNoSeq << " " << NT123.m_iTrueSize << "\t" << NT123.m_vsTrueSeq[ShowSeq].substr(0,15);;
*/


	cout << "\nTrying to initialise new model object" << flush;
	int iModPos[3] = {0,1,-1}, iBraPos[3] = {0,1,-1};
	vector <int> vModPos(3,0), vBraPos(3,0);
	FOR(i,3) { vModPos[i] = iModPos[i]; vBraPos[i] = iBraPos[i]; }
	CSiteCodon *CodonModel;

	CodonModel = new CSiteCodon(PhyDat.pData(),&Tree,vModPos,vBraPos,REV,false);

	cout << "\nModel initialised" << flush;
	cout << "\nCodonModel initialised, with lnL: " << CodonModel->lnL(true);
	cout << "\n\nDone!" << flush;

//	CodonModel->FastBranchOpt(CodonModel->lnL(true));	cout << "\nFinished branch opt: " << CodonModel->lnL(true);

	cout << "\nTrying optimiser..." << flush;
	double IThink = FullOpt(CodonModel);
	cout << "\nFull Opt gives lnl: " << CodonModel->lnL(true) << " cf. " << IThink << flush;
	Models.push_back(DoModelRun(CodonModel,15,L_EQU,0));
	cout << "\nAnd after model run " << CodonModel->lnL(true);

//	cout << "\n------------ Models -------------\n" << *CodonModel;
	// DoModelRun(CBaseModel *M, int NoPar, Lcorrection Lcor, double Adj
exit(-1);

	GetFullCodonModels(PhyDat.pData(),&Tree,&Models,GeneticCode, *out);

	exit(-1);

#endif

	///////////////////////////////////////////////////////////////////////////////////////////
	// Do the models
	// ---
	// Open up the model output file if needed
	if(ModelOut) {
		TreeName = (string)outfilestring + (string)".model.out";
		out = new ofstream(TreeName.c_str());
		//ofstream out(TreeName.c_str());
	}
	RY_count = GetRYModels(&RY_Data,&Tree,&Models,GeneticCode, *out);
        cout<<"\rRY Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
	DNA_count = GetNTModels(&NT_Data,&Tree,&Models,GeneticCode, *out);
        cout<<"\rNT Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
    if(RY_count != DNA_count) { Error("DNA[" + toString(DNA_count) + "] and RY[" + toString(RY_count) + "] character counts do not match... Fatal internal error\n"); }
	AA_count = GetAAModels(&AA_Data,&Tree,&Models,GeneticCode, *out);
        cout<<"\rAA Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
    if(AA_count * 3 != DNA_count) { Error("DNA[" + toString(DNA_count) + "] and AA[" + toString(AA_count) + "] character counts do not match... Fatal internal error\n"); }
	COD_count = GetCODModels(&COD_Data,&Tree,&Models,GeneticCode, *out);
        cout<<"\rCodon Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
    if(COD_count * 3 != DNA_count) { Error("DNA[" + toString(DNA_count) + "] and Codon[" + toString(COD_count) + "] character counts do not match... Fatal internal error\n"); }
    if(ModelOut) { assert(out != NULL); out->close(); cout << "\nOutputting models to <" << (string)outfilestring << ".model.out>"; delete out; }
    cout << "\nOutputting results to <" << outfilestring << ">";
    ofstream output(outfilestring.c_str());
	double minAIC = Models[0].AIC;
	for (i=1; i<(int)Models.size(); i++) {
		if (Models[i].AIC < minAIC) { minAIC = Models[i].AIC; };
	}
	output << "\nModel Information\n---\nModel#\tName\tDataType\tTreeLength\tOrilnL\tCorrectionType\tAdjlnL\tNoPar\tAIC\tDeltaAIC";
	FOR(i,(int)Models.size()) {
		output << "\nModel["<<i<<"]\t" << Models[i].Name << "\t" << Models[i].DataType << "\t" << Models[i].TreeLength << "\t" << Models[i].OrilnL;
		switch(Models[i].Correction) {
		case L_EQU:
			output << "\tEQUAL";
			break;
		case L_EMP:
			output << "\tEMPIRICAL";
			break;
		case L_NA:
			output << "\tNA";
			break;
		default: Error("Error in switch at end of main(...) ... ");
		}
		output << "\t" << Models[i].lnL << "\t" << Models[i].NoPar << "\t" << Models[i].AIC << "\t" << Models[i].AIC - minAIC;

	}
	output.close();
	if(RunHadError) { cout << "\nWARNING: Run had a potential error, but appeared to recover"; }
	cout << "\nSuccessful exit\n";
	return 0;
}

// Functions implementing the models
///////////////////////////////////////////////////
// ToDo:
// 1. Can probably be smarter in the order models are estimated and the way parameters are shared between them
int GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	SModelDetails RY_Model; RY_Model.DataType = RY;
	int df = -1, Ret;
	Lcorrection Correction = L_EQU;
	// Make data
	CData RY_Data = *Data; RY_Data.DNA2RY();
	Ret = RY_Data.CountMSAChars();
	// 1. Get RYmodel
	CRY RY(&RY_Data,Tree);
	// Get the correction
	double RY2Cod_Adj = Data->GetRYToCodonlnLScale(GeneticCode,&df);
	if(df > 0) { Correction = L_EMP; }
	// Optimise
	RY.lnL();
	RY_Model.OrilnL = FullOpt(&RY);
	RY_Model.lnL = RY_Model.OrilnL + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 1 + df;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	RY_Model.Correction = Correction;
	Models->push_back(RY_Model);
    cout<<"."<<flush;
    if(ModelOut) { os << RY << endl << flush; }	// Output model details
	// 2. RY+dG
	RY.MakeGammaModel(0,4);
	RY.lnL();
	RY_Model.OrilnL = FullOpt(&RY,true,FlipBool(DoItFast),FlipBool(DoItFast));
	RY_Model.lnL = RY_Model.OrilnL + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 2 + df;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	RY_Model.Correction = Correction;
	Models->push_back(RY_Model);
	cout<<"."<<flush;
	if(ModelOut) { os << RY << endl << flush; }	// Output model details
	return Ret;
}

int GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	CBaseModel *Model;
	double Alpha;
	int count = 2,Ret = Data->CountMSAChars();
	bool Fast = DoItFast;
	CTree PlainTree, GammaTree;
	// 1. JC
	CJC *JC; JC = new CJC(Data,Tree); Model = JC;
	if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
	Models->push_back(DoModelRun(Model,0,L_NA));
	if(Fast) { PlainTree = *Model->Tree(); }
	if(ModelOut) { os << *JC << endl << flush; }	// Output model details
    if(Model->NoSeq() > 2) {
    	Model->MakeGammaModel(0,4);
    	if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
    	Models->push_back(DoModelRun(Model,1,L_NA));
    	Alpha = JC->m_vpPar[0]->Val();
    	if(Fast) { GammaTree = *Model->Tree(); }
    	if(ModelOut) { os << *JC << endl << flush; }	// Output model details
    }
        cout<<"."<<flush;
	Model = NULL;
//	cout << "\nTree JCdG: \t" << Tree->TreeLength() << "\t" << JC->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	delete JC;
	DoItFast = Fast;
	// 2. FEL
	CFEL *FEL; FEL = new CFEL(Data,Tree); Model = FEL;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,3,L_NA));
	if(ModelOut) { os << *FEL << endl << flush; }	// Output model details
//	delete FEL; FEL = new CFEL(Data,Tree); Model = FEL;
//	cout << "\nTree FEL:  \t" << Tree->TreeLength() << "\t" << FEL->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	if(Model->NoSeq() > 2) {
		Model->MakeGammaModel(0,4,Alpha);
		if(DoItFast) { *Model->Tree() = GammaTree; }
		Models->push_back(DoModelRun(Model,4,L_NA));
		if(ModelOut) { os << *FEL << endl << flush; }	// Output model details
	}
        cout<<"."<<flush;
	Model = NULL;
//	cout << "\nTree FELdG:\t" << Tree->TreeLength() << "\t" << FEL->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	delete FEL;
	// 3. K2P
	CK2P *K2P; K2P= new CK2P(Data,Tree); Model = K2P;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,1,L_NA));
	if(ModelOut) { os << *K2P << endl << flush; }	// Output model details
	if(Model->NoSeq() > 2) {
		Model->MakeGammaModel(0,4,Alpha);
		if(DoItFast) { *Model->Tree() = GammaTree; }
		Models->push_back(DoModelRun(Model,2,L_NA));
		if(ModelOut) { os << *K2P << endl << flush; }	// Output model details
	}
        cout<<"."<<flush;
	Model = NULL;
	delete K2P;
	// 4. HKY
	CHKY *HKY; HKY = new CHKY(Data,Tree); Model = HKY;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,4,L_NA));
	if(ModelOut) { os << *HKY<< endl << flush; }	// Output model details
	if(Model->NoSeq() > 2) {
		Model->MakeGammaModel(0,4,Alpha);
		if(DoItFast) { *Model->Tree() = GammaTree; }
		Models->push_back(DoModelRun(Model,5,L_NA));
		if(ModelOut) { os << *HKY << endl << flush; }	// Output model details
	}
        cout<<"."<<flush;
	Model = NULL;
	delete HKY;
	// 5. GTR
	CREV *REV; REV = new CREV(Data,Tree); Model = REV;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,8,L_NA));
	if(ModelOut) { os << *REV<< endl << flush; }	// Output model details
	if(Model->NoSeq() > 2) {
		Model->MakeGammaModel(0,4,Alpha);
		if(DoItFast) { *Model->Tree() = GammaTree; }
		Models->push_back(DoModelRun(Model,9,L_NA));
		if(ModelOut) { os << *REV<< endl << flush; }	// Output model details
	}
        cout<<"."<<flush;
	Model = NULL;
	delete REV;
	return Ret;
}

// Need something that chooses an appropriate set of models based on the genetic code
int GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	double Alpha = 1.0;
	vector <double> ModelF(20,0.05);
	CData AmA = *Data; AmA.Translate(GeneticCode);
	CBaseModel *Model;
	CEMP *EMP;
	bool First = true, Fast = DoItFast;
	CTree PlainTree, GammaTree;
	int Ret = AmA.CountMSAChars();
     // The correction information
	int df = -1;
	Lcorrection Correction = L_EQU;
	double AA2Cod_Adj = Data->GetAminoToCodonlnLScale(GeneticCode,&df);
	if(df > 0) { Correction = L_EMP; }
	// Do the analyses and apply the correction ad hoc

    if (aa_model_map.find("EQU")!=aa_model_map.end()){
                // 1. EQU
//    			if(First) { DoItFast = false; }
                CEQU *EQU; EQU = new CEQU(&AmA,Tree,false); Model = EQU;
                if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                Models->push_back(DoModelRun(Model,0+df,Correction,AA2Cod_Adj));
                if(First) { PlainTree = *Model->Tree(); }
                if(ModelOut) { os << *EQU << endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4);
                	if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                	Models->push_back(DoModelRun(Model,1+df,Correction,AA2Cod_Adj));
                	if(First) { GammaTree = *Model->Tree(); First = false; DoItFast = Fast;}
                	Alpha = Model->m_vpPar[0]->Val();
                	if(ModelOut) { os << *EQU << endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete EQU;
                // 2. EQU+F
                EQU = new CEQU(&AmA,Tree,true); Model = EQU;
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,19+df,Correction,AA2Cod_Adj));
                if(ModelOut) { os << *EQU << endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4,Alpha);
                	if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                	Models->push_back(DoModelRun(Model,20+df,Correction,AA2Cod_Adj));
                	if(ModelOut) { os << *EQU << endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete EQU;
        }
        for(std::map<string,double**>::iterator iter = aa_model_map.begin(); iter != aa_model_map.end(); ++iter){
                if (iter->first == "EQU") continue;
                string name=iter->first;
                double* mySMat=iter->second[0];
                double* myFreq=iter->second[1];
                //model frequencies
//                if(First) { DoItFast = false; }
                EMP = new CEMP(&AmA,Tree,name,false,mySMat,myFreq); Model = EMP;
                if(!First && DoItFast) { *Model->Tree() = PlainTree; } else if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                Models->push_back(DoModelRun(Model,0+df,Correction,AA2Cod_Adj));
                if(First) { PlainTree = *Model->Tree(); }
                if(ModelOut) { os << *EMP << endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4,Alpha);
                	if(!First && DoItFast) { *Model->Tree() = GammaTree; } else if(DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                	Models->push_back(DoModelRun(Model,1+df,Correction,AA2Cod_Adj));
                	if(First) { GammaTree = *Model->Tree(); First = false; DoItFast = Fast; }
                	if(ModelOut) { os << *EMP << endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
                //+F
                EMP = new CEMP(&AmA,Tree,name,true,mySMat,myFreq); Model = EMP;
                if(ModelOut) { os << *EMP << endl << flush; }	// Output model details
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,19+df,Correction,AA2Cod_Adj));
                if(ModelOut) { os << *EMP << endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4,Alpha);
                	if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                	Models->push_back(DoModelRun(Model,20+df,Correction,AA2Cod_Adj));
                	if(ModelOut) { os << *EMP << endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
        }
        return Ret;
}

// Note, running with codon data regularly into Codon models doesn't work quite right...
int GetCODModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models,int GeneticCode, ostream &os) {
	// Initialise
	CData CoD = *Data;
	CBaseModel *Model;
	CCodonM0 *M0;
	int i,NoF64 = 0,Ret = -1;
	bool First = true, Fast = DoItFast;
	CTree PlainTree, GammaTree;
	// Do the calculations
	FOR(i,64) { if(GenCodes[GeneticCode][i] >= 0) { NoF64++; } }
        if (codon_model_set.count("F0")){
                // 1. F0
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,cEQU,GeneticCode); Model = M0;
                Ret = CoD.CountMSAChars();
                // if(First) { DoItFast = false; }
                if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                Models->push_back(DoModelRun(Model,2,L_NA));
                if(First) { PlainTree = *Model->Tree(); }
                if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4);
                	if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                	Models->push_back(DoModelRun(Model,3,L_NA));
                	if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                	if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }
        if (codon_model_set.count("F1X4")){
                        // 1. F1X4
                        CoD = *Data;
                        M0 = new CCodonM0(&CoD,Tree,F1X4,GeneticCode); Model = M0;
                        if(Ret == -1) { Ret = CoD.CountMSAChars(); }
//                        if(First) { DoItFast = false; }
                        if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                        if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                        Models->push_back(DoModelRun(Model,6,L_NA));
                        if(First) { PlainTree = *Model->Tree(); }
                        if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                        if(Model->NoSeq() > 2) {
                        	Model->MakeGammaModel(0,4);
                        	if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                        	if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                        	Models->push_back(DoModelRun(Model,7,L_NA));
                        	if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                        	if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                        }
                        cout<<"."<<flush;
                        Model = NULL;
                        delete M0;
        }
        if (codon_model_set.count("F3X4")){
                // 1. F3X4
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,F3X4,GeneticCode); Model = M0;
                if(Ret == -1) { Ret = CoD.CountMSAChars(); }
//                if(First) { DoItFast = false; }
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                Models->push_back(DoModelRun(Model,11,L_NA));
                if(First) { PlainTree = *Model->Tree(); }
                if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4);
                	if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                	if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                	Models->push_back(DoModelRun(Model,12,L_NA));
                	if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                	if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }
        if (codon_model_set.count("F64")){
                // 1. F64
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,F64,GeneticCode); Model = M0;
                if(Ret == -1) { Ret = CoD.CountMSAChars(); }
//                if(First) { DoItFast = false; }
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                Models->push_back(DoModelRun(Model,2 + NoF64,L_NA));
                if(First) { PlainTree = *Model->Tree(); }
                if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                if(Model->NoSeq() > 2) {
                	Model->MakeGammaModel(0,4);
                	if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                	if(First && DoItFast) { LazyBraOpt(Model,Model->lnL(),5,MATIC_BRANCH_ACC); }
                	Models->push_back(DoModelRun(Model,3 +NoF64,L_NA));
                	if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                	if(ModelOut) { os << *Model<< endl << flush; }	// Output model details
                }
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }
        if(Ret == -1) {
        	CoD = *Data;
        	CoD.MakeCodonData();
        	//CoD.ReduceCodonData(GeneticCode);
        	Ret = CoD.CountMSAChars();
        }
        return Ret;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Full codon model run, including site specific nucleotide models and empirical codon models
/* Old code for getting codon models
 *

	int i,j,k, i1,j1;
	vector <double> Mat, Freq;
	vector <vector <double> > Mat2;
	vector <string> ABET, Toks;
	vector <int> CharMap(61,-1);
	string str;
	// Assign the bigger matrix
	Mat.assign(61,-1);
	FOR(i,61) { Mat2.push_back(Mat); }
	Mat.clear();

	// Temporary stuff to read in empirical codon models
	ifstream readmodel("ECMrest.dat");
	for(i=1;i<61;i++) {
		getline(readmodel,str);
		Toks = Tokenise(str);
		assert(Toks.size() == i);
		FOR(j,i) { Mat.push_back(atof(Toks[j].c_str())); }
	}
	getline(readmodel,str);
	getline(readmodel,str);
	getline(readmodel,str);
	Toks = Tokenise(str);
	assert(Toks.size() == 61);
	FOR(i,61) { Freq.push_back(atof(Toks[i].c_str())); }
	getline(readmodel,str);
	getline(readmodel,str);
	// L1 ABET
	getline(readmodel,str);
	Toks = Tokenise(str);
	FOR(i,Toks.size()) { ABET.push_back(Toks[i]); }
	// L2 ABET
	getline(readmodel,str);
	Toks = Tokenise(str);
	FOR(i,Toks.size()) { ABET.push_back(Toks[i]); }
	// L3 ABET
	getline(readmodel,str);
	Toks = Tokenise(str);
	FOR(i,Toks.size()) { ABET.push_back(Toks[i]); }
	// L4 ABET
	getline(readmodel,str);
	Toks = Tokenise(str);
	FOR(i,Toks.size()) { ABET.push_back(Toks[i]); }
	assert(ABET.size() == 61);
	readmodel.close();
	// Test output
	k=0;
	cout << "\n\n--- Original ---";
	FOR(i,61) {
		cout << "\n" << ABET[i] << " : " << Freq[i];
		FOR(j,i) {
			Mat2[i][j] = Mat2[j][i] = Mat[k];
			cout << "\t[" << ((i-1) * 61) + j << "]" << Mat[k];
			k++;
		}
	}

//	cout << "\n\n--- Processed ---";
//	FOR(i,61) { cout << "\n" << Mat2[i]; }

	// Do the new matrix
	FOR(i,64) {
		if(GenCodes[0][i] == -1) { continue; }
		str = State(COD,i);
//		cout << "\n" << str;
		FOR(j,61) {
			if(strcmp(ABET[j].c_str(),str.c_str()) == 0) {
//				cout << " == [" << j << "]: " << ABET[j];
				CharMap[i] = j; break; }
		}

	}


	// Output
	cout << "\nFrequencies:\n";
	FOR(i,64) {
		if(GenCodes[0][i] == -1) { continue; }
		assert(CharMap[i] != -1);
		cout << Freq[CharMap[i]] << ",";
	}
	cout << "\n\nMatrix:\n";
	FOR(i,64) {
		if(GenCodes[0][i] == -1) { continue; }
		i1 = CharMap[i];
		FOR(j,i) {
			if(GenCodes[0][j] == -1) { continue; }
//			cout << "\n" << ABET[CharMap[i]] << " -> " << ABET[CharMap[j]] << " = ";
			cout << Mat2[i1][CharMap[j]] << " ,";
		}
	}


 *
 */

int GetFullCodonModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out)	{
	int i, j, k;
	CEMPCodonREST *NewModel1; CEMPCodonUNREST *NewModel2;
	CData CodData = *Data;
	assert(GeneticCode == 0);

	CCodonM0 *CodModel;
	CodData = *Data;
	CodModel = new CCodonM0(&CodData,Tree,F64,0);
	cout << "\nOld lnL: " << CodModel->lnL(true);
	DoModelRun(CodModel,63,L_NA);
	cout << "\nOptimised lnL: " << CodModel->lnL(true);

	CodData = *Data;
	NewModel1 = new CEMPCodonREST(&CodData, Tree, false, 0);
	cout << "\nRestrained model made..." << flush;
	cout << "\nNew lnL: " << NewModel1->lnL(true);
	DoModelRun(NewModel1,63,L_NA);
	cout << "\nOptimised lnL: " << NewModel1->lnL(true);
	delete NewModel1;

	CodData = *Data;
	NewModel2 = new CEMPCodonUNREST(&CodData, Tree, false, 0);
	cout << "\nUnrestrained model made..." << flush;
	cout << "\nNew lnL: " << NewModel2->lnL(true);
	DoModelRun(NewModel2,63,L_NA);
	cout << "\nOptimised lnL: " << NewModel2->lnL(true);
	delete NewModel2;


	cout << "\nSuccessful exit...";
	exit(-1);

}

SModelDetails DoModelRun(CBaseModel *M, int NoPar, Lcorrection Lcor, double Adj) {
	SModelDetails ModDet;
	ModDet.DataType = M->m_pData->m_DataType;
	double CurlnL;
	int NoIter = 5 + (NoPar * 1.5);
	/* --- OLD VERSION OF THE DoItFast OPTION ---
	M->lnL();
	ModDet.OrilnL = FullOpt(M,true,FlipBool(DoItFast)) ;
	*/
	// --- NEW VERSION OF THE DoItFast OPTION ---
	CurlnL = M->lnL();
//	cout << "\n-----------------------------------------------\nModel: " << M->Name();
	if(DoItFast) {
		if(M->m_pData->m_DataType == DNA || M->m_pData->m_DataType == COD || M->m_pData->m_DataType == COD_RED) { NoIter = max(10,NoIter); }
//		cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>> LazyBraOpt1 <<<<<<<<<<<<<<<<<<<<<<< ";
		CurlnL = LazyBraOpt(M,CurlnL,1,MATIC_BRANCH_ACC);
//		cout << " CurlnL: " << CurlnL;
//		cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>> LazyParOpt <<<<<<<<<<<<<<<<<<<<<<< ";
		CurlnL = LazyOpt(M,true,false,false,CurlnL,false,NoIter);
//		cout << " CurlnL: " << CurlnL;
//		cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>> LazyBraOpt2 <<<<<<<<<<<<<<<<<<<<<<< ";
		ModDet.OrilnL = LazyBraOpt(M,CurlnL,1,MATIC_BRANCH_ACC);
//		cout << " CurlnL: " << CurlnL;
	} else {
		ModDet.OrilnL = FullOpt(M,true,true);
	}
//	cout << "\n\t ... final: ret= " << ModDet.OrilnL << " cf. real= " << M->lnL(true);

//	cout << "\nFinished opt: " << ModDet.OrilnL << " cf. actual lnL calc: " << M->lnL(true);
//	double FullOpt(CBaseModel *Model, bool DoPar, bool DoBra, bool DoFreq, double CurlnL,bool FullLikAcc, int NoIterations, double lnL2Beat, double lnL_tol, bool DoOutput,bool TightFullOpt)	{
//	ModDet.OrilnL = FullOpt(M,true,FlipBool(DoItFast),false,-BIG_NUMBER,true, DEFAULT_OPTNUM,-BIG_NUMBER, FULL_LIK_ACC,true,false);
	ModDet.lnL = ModDet.OrilnL + Adj;
	ModDet.NoPar = NoPar; ModDet.AIC = GetAIC(ModDet.lnL,ModDet.NoPar);
	ModDet.TreeLength = M->Tree()->GetTreeLength();
	ModDet.Name = M->Name();
	ModDet.Correction = Lcor;
#if CHECK_LNL_OUT == 1
	// Error checking. This no longer throws a terminal error, but instead makes the program rerun
	if(fabs(ModDet.OrilnL - M->lnL(true)) > 0.001)	{
			cout << "\nModel: " << M->Name() << " obtained " << ModDet.lnL << " cf. " << M->lnL() << " cff. " << M->lnL(true) << "\n\n";
			cout << "\nModel details: " << *M; exit(-1);
	}
#endif
	return ModDet;
}

bool GetModels(string file) {
	//set up aa_model_map
    //
	int count;
    for (int i=0; i < n_models; i++){
            double** vals = (double**) malloc(2*sizeof(double*));
            vals[0]=smat[i];
            vals[1]=freq[i];
            aa_model_map[model_names[i]]=vals;
    }
    codon_model_set.insert("F0");
    codon_model_set.insert("F1X4");
    codon_model_set.insert("F3X4");
    codon_model_set.insert("F64");

    //read conf
    string confFile="modelomatic.ini";
    char const* home_c = getenv("HOME");
    string homedir = (home_c == NULL) ? std::string() : std::string(home_c);
    string dotConfFile=homedir.append("/.modelomatic.ini");
    INIReader reader(confFile);
    INIReader dotReader(dotConfFile);
    if (reader.ParseError() > 0) {
    	std::cout << "Error parsing " << confFile << "\n" ;
        return 1;
    }
    if (dotReader.ParseError() > 0) {
    	std::cout << "Error parsing " << dotConfFile << "\n" ;
    	return 1;
    }
    count = 0;
    for (int i=0; i < n_models; i++){
    	//default true, ini in cwd takes priority
    	if (!reader.GetBoolean("Amino Acid",model_names[i],dotReader.GetBoolean("Amino Acid",model_names[i],true))){
    		aa_model_map.erase(model_names[i]);
    		if(count == 0) { cout << "\n\tAmino acid models skipping " << model_names[i]; count++; }
    		else { cout << ", " << model_names[i]; }
    	}
    }
    //handle EQU separately
    if (reader.GetBoolean("Amino Acid","EQU",dotReader.GetBoolean("Amino Acid","EQU",true))){
    	aa_model_map["EQU"]=NULL;
    }
    std::set <string> new_codon_model_set;
    std::set<string>::iterator it;
    count = 0;
    for( it = codon_model_set.begin(); it != codon_model_set.end(); it++ ) {
//        	cout << "\nGoing through codon models";
    	if (reader.GetBoolean("Codon",*it,dotReader.GetBoolean("Codon",*it,true))){
//        		cout << "\nInserting " << *it;
    		new_codon_model_set.insert(*it);
    	}else {
    		if(count == 0) { cout << "\n\tCodon models skipping " << *it; count++; }
    		else { cout << ", " << *it; }
    	}
    }
    codon_model_set=new_codon_model_set;

	return true;
}

