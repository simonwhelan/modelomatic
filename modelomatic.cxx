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
#define VERSION_NUMBER "1.0a"
#define DEVELOPER_VERSION_MAIN 1

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
		int i,j,NumModelReruns = 1;
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
	PhyDat.pData()->RemoveSparseSeqs(true,NULL);
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
        tmp_AA_Data.Translate();
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
	// 3. Set genetic code if required
	if(argc>4) {
		assert(InRange(atoi(argv[4]),0,NumGenCode));
		GeneticCode = atoi(argv[4]);
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
	PhyDat.pData()->CleanToDNACodon();
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
	cout << "\nDoing codon based analysis"; cout.precision(10);

	CData NT1 = *PhyDat.pData(); NT1.GetCodonPositions(true,false,false);
	CData NT2 = *PhyDat.pData(); NT2.GetCodonPositions(false,true,false);
	CData NT3 = *PhyDat.pData(); NT3.GetCodonPositions(false,false,true);
/*	CData NT12 = *PhyDat.pData(); NT12.GetCodonPositions(true,true,false);
	CData NT13 = *PhyDat.pData(); NT13.GetCodonPositions(true,false,true);
	CData NT23 = *PhyDat.pData(); NT23.GetCodonPositions(false,true,true);
	CData NT123 = *PhyDat.pData(); NT123.GetCodonPositions(true,true,true);*/

//	int ShowSeq = RandInt(0,NT1.m_iNoSeq-1);
//	cout << "\nOriginal data:		   " << PhyDat.pData()->m_iNoSeq << " " << PhyDat.pData()->m_iTrueSize << "\t" << PhyDat.pData()->m_vsTrueSeq[ShowSeq].substr(0,15);;

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

	/*
	cout << "\nTrying to initialise new model object";
	int iModPos[3] = {0,1,2}, iBraPos[3] = {0,1,2};
	vector <int> vModPos(3,0), vBraPos(3,0);
	FOR(i,3) { vModPos[i] = iModPos[i]; vBraPos[i] = iBraPos[i]; }
	CSiteCodon *CodonModel;
	CodonModel = new CSiteCodon(PhyDat.pData(),&Tree,vModPos,vBraPos,REV,false);


	cout << "\nCodonModel initialised, with lnL: " << CodonModel->lnL(true);
	cout << "\n\nDone!" << flush;
//			CSiteCodon::CSiteCodon(CData *D, CTree *T, vector <int> ModelPar, vector <int> BranchPar) : CBaseModel(D,T)	{

//	CodonModel->FastBranchOpt(CodonModel->lnL(true));	cout << "\nFinished branch opt: " << CodonModel->lnL(true);

	cout << "\nTrying optimiser...";
//	double IThink = FullOpt(CodonModel);
//	cout << "\nFull Opt gives lnl: " << CodonModel->lnL(true) << " cf. " << IThink;
	Models.push_back(DoModelRun(CodonModel,15,L_EQU,0));
	cout << "\nAnd after model run " << CodonModel->lnL(true);

//	cout << "\n------------ Models -------------\n" << *CodonModel;
	// DoModelRun(CBaseModel *M, int NoPar, Lcorrection Lcor, double Adj
	*/
	cout <
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
const double dECMrest[1830] = {0.366267 ,10.762877 ,0.982893 ,4.929483 ,13.863648 ,1.010212 ,3.91936 ,0 ,0 ,0 ,0 ,2.13474 ,0 ,0 ,6.967655 ,0 ,0 ,5.869857 ,0 ,12.678802 ,38.414969 ,0 ,0 ,0 ,11.715476 ,28.307876 ,6.03574 ,19.578982 ,8.423762 ,0 ,0 ,0 ,0.754055 ,0 ,0 ,0 ,0 ,7.784768 ,0 ,0 ,0 ,21.337943 ,0 ,0 ,0.1209 ,0 ,0 ,21.000358 ,0 ,0 ,0 ,22.577271 ,0 ,36.292705 ,6.01197 ,0 ,0 ,0 ,16.623529 ,0 ,0 ,0 ,39.399322 ,1.792245 ,26.637668 ,3.324581 ,2.64762 ,0 ,0 ,0 ,4.96235 ,0 ,0 ,0 ,1.903175 ,0 ,0 ,0 ,0 ,0.292691 ,0 ,0 ,0 ,2.231468 ,0 ,0 ,0 ,0.437779 ,0 ,0 ,4.623936 ,0 ,0 ,0.975074 ,0 ,0 ,0 ,8.155285 ,0 ,0 ,0 ,2.231662 ,0 ,2.082128 ,1.303951 ,0 ,0 ,0 ,0.501054 ,0 ,0 ,0 ,2.801564 ,0 ,0 ,0 ,0.578118 ,18.2989 ,10.140728 ,1.281784 ,6.404436 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.105602 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.918904 ,0 ,0 ,2.715692 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,14.959151 ,5.369463 ,0 ,0 ,0 ,5.004762 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.44038 ,28.08616 ,1.086327 ,0 ,0 ,0 ,0 ,0.925575 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.381868 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.43455 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.774544 ,0 ,0 ,12.053827 ,0 ,0 ,0 ,0 ,0 ,0 ,1.23845 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.794388 ,0 ,11.295213 ,28.557939 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,6.578976 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.021682 ,35.480963 ,17.461627 ,8.407091 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,92.372238 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,6.415757 ,0 ,0 ,0 ,1.693662 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.03363 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.016899 ,0 ,0 ,0 ,1.699302 ,0 ,0 ,14.930266 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,133.296291 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,19.920977 ,0 ,0 ,0 ,3.144296 ,0 ,13.919752 ,79.48373 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.834666 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.019726 ,0 ,0 ,0 ,0.334165 ,99.459951 ,14.603857 ,30.80475 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.134459 ,0 ,0 ,0 ,7.419058 ,0 ,0 ,0 ,1.39368 ,0 ,0 ,0 ,3.023939 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,6.381841 ,0 ,0 ,0 ,3.230751 ,0 ,0 ,0 ,1.389735 ,0 ,0 ,0 ,2.25077 ,0 ,0 ,7.811205 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.570123 ,0 ,0 ,0 ,1.81254 ,0 ,0 ,0 ,0.477616 ,0 ,0 ,0 ,1.462945 ,0 ,5.650116 ,22.078564 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.301225 ,0 ,0 ,0 ,1.583116 ,0 ,0 ,0 ,1.263838 ,0 ,0 ,0 ,0.779565 ,56.803876 ,8.905484 ,8.432339 ,2.811398 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.485369 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.686074 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.006827 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.650088 ,0 ,0 ,2.090641 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.686782 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,10.15357 ,3.779053 ,0 ,0 ,0 ,6.091654 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.383706 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.71061 ,10.605899 ,1.296294 ,0 ,0 ,0 ,0 ,3.064069 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.923972 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.774445 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.312255 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.120189 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.866565 ,0 ,0 ,5.090423 ,0 ,0 ,0 ,0 ,0 ,0 ,25.461549 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,6.034856 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,9.441919 ,0 ,9.529141 ,20.715347 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.107629 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.159185 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.693026 ,21.910225 ,6.516319 ,5.512586 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.43241 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.091041 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.874 ,0 ,0 ,0 ,2.985093 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.285543 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.25247 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.561629 ,0 ,0 ,0 ,2.303487 ,0 ,0 ,4.120953 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.388514 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.803738 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.39381 ,0 ,0 ,0 ,6.644971 ,0 ,19.084271 ,18.064826 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,10.59078 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.04115 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.042624 ,0 ,0 ,0 ,1.541379 ,20.5181 ,9.48852 ,13.246936 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,25.294298 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.737489 ,0 ,0 ,0 ,1.277861 ,0 ,0 ,0 ,6.291148 ,0 ,0 ,0 ,0.702411 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,10.137337 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.914972 ,0 ,0 ,0 ,0.208522 ,0 ,0 ,0 ,4.308415 ,0 ,0 ,0 ,0.542717 ,0 ,0 ,3.531461 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.078444 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.164805 ,0 ,0 ,0 ,0.476105 ,0 ,0 ,0 ,6.166554 ,0 ,0 ,0 ,0.302501 ,0 ,11.898141 ,21.657664 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.096281 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.473623 ,0 ,0 ,0 ,0.352915 ,0 ,0 ,0 ,3.682092 ,0 ,0 ,0 ,0.503385 ,26.045078 ,6.669955 ,8.901167 ,0 ,0.228361 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.087801 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.294702 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.245405 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.444964 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.367553 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,19.942818 ,0 ,0 ,0 ,0 ,6.033932 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.312811 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,9.880382 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,15.793595 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.704601 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.675838 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.615472 ,0 ,13.681751 ,0 ,0 ,0 ,0 ,0 ,0 ,17.103904 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.30348 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,21.863158 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,12.991861 ,65.097021 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,38.2291 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.769355 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,5.540052 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.206054 ,44.407955 ,13.889102 ,17.057443 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.047654 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.632945 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.552851 ,0 ,0 ,0 ,0 ,0 ,0 ,0.144278 ,0 ,0 ,9.339206 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.825082 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,3.026086 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.810856 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.662001 ,0 ,0.073268 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4.719489 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.617091 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.104727 ,0 ,0 ,0 ,0 ,0 ,0.582084 ,0 ,0 ,0 ,6.191481 ,44.777964 ,0.677177 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,10.116588 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,75.752638 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.36272 ,0 ,0 ,0 ,0 ,0 ,2.952332 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.786043 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.266373 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.923392 ,0 ,0 ,5.17793 ,0 ,0 ,1.913571 ,0 ,0 ,1.558523 ,0 ,0 ,0.010896 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,7.911096 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,20.877218 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.022101 ,0 ,0 ,0 ,0 ,0 ,8.126914 ,0 ,0 ,5.369644 ,0 ,24.748755 ,4.756288 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1.682029 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.090751 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,2.261813 ,0 ,6.610894 ,0 ,0 ,0 ,1.658051 ,0 ,0 ,3.347364 ,1.31561 ,11.192024 ,5.427076};
const double dECMunrest[1830] = {0.413957 ,12.931524 ,1.038682 ,2.075154 ,17.634677 ,0.524647 ,1.523251 ,0.259909 ,0.256174 ,1.243753 ,0.089476 ,1.346385 ,0.623366 ,0.067554 ,8.373639 ,0.199589 ,0.822133 ,1.144639 ,0.224397 ,19.459275 ,28.470047 ,0.878163 ,0.296066 ,0.291009 ,2.057449 ,31.915858 ,13.991583 ,12.116657 ,5.815294 ,0.273855 ,1.197614 ,0.93684 ,1.316696 ,0.021698 ,0.087905 ,0.432918 ,0.334224 ,4.127466 ,0.928876 ,0.643789 ,0.441084 ,4.699375 ,2.871848 ,0.622868 ,0.403913 ,1.868194 ,0.600415 ,7.316623 ,0.454347 ,0.511968 ,0.140701 ,0.765495 ,0.309656 ,43.916187 ,1.28699 ,1.29386 ,0.586685 ,0.389004 ,4.971507 ,2.397415 ,0.275865 ,0.563566 ,4.943439 ,1.644621 ,28.579806 ,1.477696 ,0.667397 ,0.020837 ,0.056267 ,0.256283 ,2.260424 ,0.016701 ,0.154441 ,0.351882 ,0.780008 ,0.021539 ,0.317466 ,0.269862 ,0.02302 ,0.175491 ,0.116405 ,0.014351 ,0.048902 ,1.066285 ,0.311759 ,0.075878 ,0.012331 ,0.238533 ,0.062369 ,0.018727 ,5.784672 ,0.275355 ,0.21917 ,0.497593 ,0.208924 ,0.853779 ,0.433759 ,1.376727 ,0.419831 ,0.305714 ,0.304581 ,0.90581 ,0.306375 ,2.230691 ,1.428293 ,0.221536 ,0.025405 ,0.04144 ,0.239388 ,0.296014 ,0.034573 ,0.090674 ,1.266654 ,0.116039 ,0.03114 ,0.093542 ,0.359183 ,16.415611 ,13.60931 ,1.155098 ,3.215817 ,0.377447 ,0.461383 ,1.759732 ,1.660885 ,0.066519 ,0.002854 ,0.374419 ,1.24214 ,0.232701 ,0.003786 ,1.177718 ,0.390087 ,0.010543 ,0.443617 ,0.078926 ,0.143479 ,2.203539 ,0.571469 ,0.188587 ,0.02192 ,0.969387 ,0.10308 ,0.000272 ,0.025544 ,1.282522 ,0.169022 ,0.000453 ,0.000272 ,0.127899 ,0.163406 ,0 ,0.751904 ,0.285384 ,0.985145 ,2.285052 ,0.206851 ,0.001302 ,0.617464 ,1.766961 ,0.029088 ,0.000289 ,1.418383 ,1.148387 ,0.083261 ,0.000048 ,0.105885 ,0.834976 ,0.003859 ,17.923045 ,2.503268 ,0.899262 ,0.46097 ,0.324525 ,2.691226 ,0.082373 ,0.001125 ,0.076938 ,0.898794 ,0.411115 ,0.00506 ,0.140943 ,1.293228 ,0.040109 ,0.000469 ,0.323588 ,0.193328 ,2.473439 ,38.685701 ,1.106179 ,1.042457 ,0.074161 ,0.012002 ,0.096104 ,1.018005 ,0.027138 ,0.083475 ,0.345545 ,0.493508 ,0.005195 ,0.009225 ,0.021048 ,0.480521 ,0.000448 ,0.042006 ,0.002597 ,1.05204 ,0.089835 ,0.184327 ,0.286252 ,0.012394 ,0.532092 ,0.303613 ,0.000192 ,0.02623 ,0.870967 ,0.435243 ,0.053805 ,0.000384 ,1.337629 ,0.033244 ,0.000288 ,0 ,0.270274 ,0.008455 ,0.000096 ,0.074846 ,0.730306 ,0.529594 ,0.091276 ,15.426733 ,0.00041 ,0.095418 ,0.561523 ,0.002705 ,0.015001 ,0.527094 ,0.962295 ,0.011476 ,0 ,0.104681 ,0.273876 ,0.001066 ,0 ,0.021723 ,0.193459 ,0.000246 ,0.063366 ,0.364129 ,0.715716 ,0.054021 ,15.127582 ,33.45378 ,0.247195 ,0.000959 ,0.220816 ,0.78734 ,0.368538 ,0.044125 ,0.079712 ,1.050363 ,0.011511 ,0.001247 ,0.060432 ,1.035015 ,0.002014 ,0.000096 ,0.081631 ,0.292855 ,0.435109 ,0.061583 ,0.107434 ,0.775254 ,40.922701 ,23.496083 ,10.291826 ,3.06551 ,0.016636 ,0.034519 ,0.009981 ,0.360575 ,0.002079 ,0.002911 ,0.00104 ,31.949764 ,0.354545 ,15.287419 ,0.530881 ,0.121855 ,0.001248 ,0.015596 ,0.000416 ,4.72584 ,0.526904 ,0.61732 ,1.35386 ,0.626603 ,0.154017 ,0.091073 ,0.352807 ,0.296358 ,0.863191 ,2.721043 ,0.000068 ,0 ,1.481721 ,0.014243 ,0 ,1.113671 ,1.320875 ,4.087042 ,0.018146 ,0 ,0.238839 ,0.029786 ,0 ,0.14894 ,2.478358 ,2.029914 ,0.170088 ,0.001702 ,0.563799 ,0.798346 ,0.002398 ,12.677657 ,0.011032 ,0.065212 ,3.8322 ,0.000981 ,0 ,0.017774 ,2.985421 ,0 ,2.418859 ,1.360141 ,18.531553 ,0.025006 ,0 ,0.003064 ,0.40868 ,0 ,0.299978 ,2.20843 ,4.718199 ,0.316372 ,0.0103 ,0.897101 ,1.164525 ,0.00988 ,30.574631 ,35.563093 ,1.769879 ,0.00318 ,0.560831 ,0.682254 ,0.000589 ,0.00106 ,0.004593 ,0.344601 ,17.450524 ,0.387941 ,6.798629 ,0.772349 ,0 ,0.000118 ,0.028736 ,0.076434 ,1.164262 ,0.277265 ,0.340811 ,2.460976 ,0.063413 ,0.124232 ,0.044676 ,0.500433 ,39.595443 ,27.244097 ,12.868484 ,0.43362 ,0.01454 ,0.000679 ,0.019568 ,0.94456 ,0.000679 ,0.000408 ,0.004348 ,0.409432 ,0.002718 ,0.000951 ,0.009512 ,4.40061 ,0.203017 ,1.635281 ,1.328581 ,1.205738 ,0.059927 ,0.105449 ,0.438648 ,1.07379 ,0.040087 ,0.016442 ,0.181411 ,0.918198 ,0.005436 ,0.040902 ,0.134394 ,0.002143 ,0.18788 ,0.022384 ,0 ,0.000119 ,0.89749 ,0.000595 ,0.00006 ,0.00006 ,0.300513 ,0.000952 ,0.00006 ,0.136148 ,3.765697 ,1.434579 ,0.18401 ,0.071794 ,0.617694 ,0.422671 ,0.057328 ,0.016609 ,0.648769 ,0.118229 ,0.036254 ,0.112395 ,0.544828 ,0.640495 ,0.068758 ,8.808275 ,0.000035 ,0.034916 ,0.133032 ,0.00014 ,0.00007 ,0.073657 ,0.733202 ,0.00007 ,0 ,0.0173 ,0.04032 ,0.000105 ,0.048882 ,2.469249 ,2.236557 ,0.089308 ,0.033372 ,0.330036 ,0.703795 ,0.044742 ,0.040917 ,0.149771 ,0.495176 ,0.014388 ,0.041969 ,0.563956 ,0.794366 ,0.015019 ,8.319767 ,15.484088 ,0.060922 ,0.00033 ,0.039362 ,0.25938 ,0.000857 ,0.000396 ,0.000264 ,0.831215 ,0.004747 ,0.00033 ,0.00956 ,0.322807 ,1.260239 ,0.39019 ,1.22002 ,3.330793 ,0.529045 ,0.058614 ,0.079119 ,0.664866 ,0.065669 ,0.026571 ,0.007978 ,0.822643 ,0.424607 ,0.016022 ,0.039428 ,0.663877 ,21.437257 ,15.874441 ,4.802017 ,2.258321 ,0.193544 ,0.308851 ,0.975582 ,1.164866 ,0.109664 ,0.003464 ,0.251423 ,0.668173 ,0.308372 ,0.00442 ,0.832592 ,0.361403 ,0.009174 ,0.033994 ,0.038612 ,3.081614 ,0.106069 ,0.536567 ,0.452442 ,0.890822 ,0.105112 ,0.00541 ,0.184749 ,0.440217 ,0.025678 ,0.000858 ,0.034978 ,0.313495 ,0.002679 ,0.00014 ,0.03956 ,0.053168 ,2.729647 ,0.426684 ,0.167041 ,0.003354 ,1.013424 ,0.095207 ,0.000158 ,0.000718 ,2.005334 ,0.008986 ,0.000388 ,0 ,0.032447 ,0.011379 ,0 ,0.121389 ,0.975727 ,0.928319 ,0.100928 ,0.003672 ,0.976175 ,0.644235 ,0.000096 ,0.009981 ,0.915505 ,0.150527 ,0.000236 ,0.001359 ,0.135255 ,0.015686 ,0.000066 ,0.892903 ,0.159248 ,0.372719 ,1.563846 ,0.149578 ,0.000838 ,0.360196 ,1.447392 ,0.033581 ,0.00041 ,0.822643 ,1.19579 ,0.129098 ,0.000078 ,0.031984 ,0.262465 ,0.003312 ,0.343314 ,0.566034 ,3.036767 ,0.159217 ,0.001791 ,0.294294 ,0.99033 ,0.034437 ,0.00104 ,0.528692 ,2.171984 ,0.007891 ,0.000136 ,0.010061 ,0.1763 ,0.00422 ,12.984714 ,3.195698 ,0.604859 ,0.371685 ,0.188097 ,3.529593 ,0.092988 ,0.000559 ,0.042595 ,0.928969 ,0.030474 ,0.012481 ,0.015763 ,1.702817 ,0.001171 ,0.000042 ,0.059065 ,0.066966 ,0.706782 ,0.046739 ,0.130873 ,0.746792 ,0.147873 ,0.000769 ,0.03156 ,0.661776 ,0.00707 ,0 ,0.000368 ,0.162525 ,0.012366 ,0 ,0.000491 ,0.129492 ,2.738552 ,14.214694 ,1.091224 ,0.923644 ,0.041435 ,0.004382 ,0.08462 ,2.735831 ,0.077341 ,0.318574 ,0.381796 ,0.461729 ,0.063353 ,0.00442 ,0.126382 ,0.983926 ,0.000589 ,0.517896 ,0.002262 ,1.274118 ,0.004982 ,0.002219 ,0.039265 ,1.236457 ,0.069754 ,0.087384 ,0.335924 ,0.286547 ,0.000068 ,0.000123 ,0.000589 ,0.741135 ,0.000357 ,0.000351 ,0.002374 ,0.905488 ,0.027912 ,0.157605 ,0.18909 ,0.0055 ,0.565968 ,0.373948 ,0.000045 ,0.019665 ,2.159328 ,0.498504 ,0.038559 ,0 ,2.214735 ,0.001768 ,0.000485 ,0 ,0.536362 ,0.169264 ,0.00004 ,0.018814 ,0.696016 ,0.634927 ,0.000562 ,0.013077 ,1.394892 ,0.971148 ,0.022062 ,0.000416 ,1.610248 ,0.020471 ,0.000236 ,0.000136 ,1.191514 ,0.103766 ,0.000132 ,0.140546 ,0.802516 ,0.938586 ,0.01514 ,6.74656 ,0.000482 ,0.068501 ,0.759132 ,0.003771 ,0.049848 ,0.728891 ,2.866325 ,0.023862 ,0 ,0.200205 ,0.187832 ,0.006402 ,0.000078 ,0.004629 ,0.837302 ,0.000081 ,0.001165 ,0.190037 ,1.544626 ,0.080124 ,0.017824 ,0.417565 ,1.950083 ,0.014101 ,0.000624 ,0.083744 ,3.409178 ,0.001531 ,0 ,0.000893 ,1.524024 ,0.000066 ,0.227527 ,0.702594 ,1.251589 ,0.043198 ,11.561124 ,19.373137 ,0.216244 ,0.001329 ,0.268662 ,0.64316 ,0.324162 ,0.120121 ,0.108501 ,2.359362 ,0.013852 ,0.005328 ,0.045374 ,2.323091 ,0.001795 ,0.000253 ,0.451134 ,0.536332 ,0.27295 ,0.000091 ,0.040907 ,0.42058 ,0.346978 ,0.073021 ,0.013526 ,1.272426 ,0.001871 ,0 ,0 ,0.324227 ,0.004077 ,0.000417 ,0.000772 ,0.879741 ,0.502355 ,0.018617 ,0.106257 ,0.396521 ,24.400553 ,9.612709 ,5.140068 ,0.713961 ,0.2804 ,0.00176 ,0.909861 ,0.46677 ,0.000112 ,0.000187 ,0.002449 ,0.938439 ,0.298738 ,0.165587 ,1.179053 ,0.194319 ,0.000084 ,0.037724 ,0.000283 ,0.551493 ,0.004982 ,0.000289 ,0.097554 ,0.556204 ,0.000384 ,0.000164 ,0.003549 ,0.247245 ,0.000753 ,0.010419 ,0.01531 ,0.165784 ,0.00006 ,0.00007 ,0.00033 ,0.639125 ,0.067431 ,0.053748 ,0.375097 ,1.693998 ,0.036022 ,0.046588 ,0.658541 ,0.002863 ,1.110726 ,0.176224 ,0.001122 ,0 ,0.518903 ,0.006459 ,0 ,0.000821 ,2.392606 ,0.012227 ,0.054025 ,0 ,0.055425 ,0.018603 ,0.000283 ,0.016193 ,0.566759 ,0.336277 ,0.000843 ,0.00009 ,0.636916 ,0.395771 ,0.000096 ,0.008942 ,0.676049 ,0.703728 ,0.021435 ,0 ,0.240505 ,0.047268 ,0 ,0.030585 ,1.228004 ,0.239725 ,0.037461 ,0.027973 ,1.797671 ,1.398079 ,0.022214 ,5.214689 ,0.001158 ,0.293884 ,0.744 ,0.194964 ,0.00061 ,0.000447 ,0.716441 ,0.000553 ,0.069567 ,1.149846 ,1.558784 ,0.684968 ,0.000078 ,0.000295 ,0.410199 ,0.000081 ,0.003029 ,0.053714 ,0.634203 ,0.08406 ,0.001254 ,0.001441 ,0.703092 ,0.001918 ,0.143481 ,0.097165 ,0.590833 ,0.043222 ,0 ,0.000119 ,0.185038 ,0.000066 ,0.136937 ,0.403521 ,0.968146 ,0.235747 ,0.672077 ,0.496552 ,1.526141 ,0.238907 ,25.64086 ,18.675227 ,0.04873 ,0.005905 ,0.021649 ,1.322234 ,0.000534 ,0.000112 ,0.000374 ,0.380374 ,0.037657 ,0.326035 ,0.047437 ,2.382451 ,0.000156 ,0 ,0.027479 ,0.069995 ,0.042812 ,0.000091 ,0.003377 ,0.388156 ,0.002508 ,0.000288 ,0.001312 ,0.418802 ,0.042005 ,0.004656 ,0.011768 ,0.241197 ,0.000272 ,0.00006 ,0.000175 ,0.169777 ,0.089148 ,0.03512 ,0.030223 ,0.629243 ,0.585071 ,0.086923 ,0.083552 ,1.572808 ,21.171295 ,13.93595 ,13.981617 ,0.708204 ,0.002362 ,0.000524 ,0.008844 ,3.879355 ,0.001957 ,0.002528 ,0.009008 ,0.560642 ,0.000657 ,0.002062 ,0.006305 ,11.271062 ,0.313779 ,0.569076 ,2.581654 ,0.606945 ,0.002808 ,0.000145 ,0.019586 ,0.732111 ,0.000192 ,0.000082 ,0.001151 ,0.255147 ,0 ,0 ,0.000353 ,2.047159 ,0.03173 ,0.020529 ,0.509595 ,0.469281 ,0.006776 ,0.053017 ,0.078381 ,3.415722 ,0.028569 ,0.081861 ,0.662578 ,0.324696 ,0.001557 ,0.074536 ,0.072521 ,0.00193 ,0.215492 ,0.056145 ,0 ,0.002515 ,2.862216 ,0.010298 ,0.001027 ,0.000205 ,0.381287 ,0.002062 ,0.0000097 ,0.154207 ,7.487816 ,0.350473 ,0.458622 ,0.006757 ,0.491034 ,0.06932 ,0.000375 ,0.000537 ,0.76009 ,0.031642 ,0.000192 ,0.001871 ,0.648523 ,0.005026 ,0.000589 ,0.043756 ,2.720153 ,1.08517 ,0.105361 ,0.042565 ,0.119062 ,0.245019 ,0.007112 ,0.092405 ,2.620833 ,0.968351 ,0.131929 ,0.009773 ,0.201477 ,0.064227 ,0.016776 ,5.637509 ,0.000193 ,0.013582 ,0.311885 ,0.000628 ,0.011433 ,0.089643 ,3.0925 ,0.003635 ,0.000103 ,0.00416 ,0.314379 ,0.001358 ,0.911193 ,3.487342 ,1.399673 ,0.953353 ,0.000466 ,0.044112 ,0.475496 ,0.032237 ,0.000717 ,0.001537 ,0.602428 ,0.000384 ,0.000832 ,0.003287 ,0.701399 ,0.001413 ,0.173122 ,0.477439 ,2.26249 ,0.091844 ,0.055025 ,0.04128 ,0.362328 ,0.011491 ,0.415718 ,0.682579 ,2.211488 ,0.152215 ,0.016839 ,0.051048 ,0.276276 ,0.024883 ,13.945647 ,18.744445 ,0.20148 ,0.000295 ,0.072138 ,0.318595 ,0.05625 ,0.00302 ,0.004868 ,2.898791 ,0.032116 ,0.000219 ,0.051709 ,0.501456 ,2.577578 ,0.575921 ,0.430124 ,6.546163 ,0.150804 ,0.000091 ,0.0096 ,0.338207 ,0.009225 ,0.000192 ,0.000246 ,0.708395 ,0.003743   ,0 ,0 ,0.336828 ,0.555104 ,0.072866 ,0.01158 ,1.89715 ,0.276499 ,0.004101 ,0.049842 ,0.153418 ,0.718719 ,0.046583 ,0.054049 ,2.600428 ,0.080337 ,0.00205 ,0.032177 ,0.164143 ,25.313949 ,10.956917 ,8.832621 ,0.042942 ,0.487814 ,0.140864 ,0.042762 ,0.053363 ,0.393545 ,0.177757 ,0.000599 ,0.060909 ,0.34689 ,0.26951 ,0.00018 ,0.025454 ,0.255975 ,0.300354 ,0.00018 ,0.079236 ,3.193996 ,0.332396 ,0.340541 ,0.025873 ,0.310416 ,0.037372 ,0.000479 ,0.051686 ,0.604541 ,0.081931 ,0.001377 ,0.011978 ,0.622089 ,0.308859 ,0.000299 ,0.043721 ,0.184704 ,0.136432 ,0.015452 ,0.034737 ,0.380788 ,0.106187 ,0.000419 ,0.029826 ,0.127388 ,0.014314 ,0.00012 ,0.06043 ,0.582801 ,0.116309 ,0.000359 ,0.467136 ,0.046586 ,0.072867 ,0.74787 ,0.208041 ,0.000116 ,0.042119 ,0.599932 ,0.424321 ,0.000116 ,0.200499 ,0.50647 ,0.573013 ,0.000116 ,0.294368 ,0.352731 ,0.486281 ,0.265766 ,0.09207 ,3.088481 ,0.03214 ,0.000116 ,0.004177 ,0.337879 ,0.043917 ,0.000116 ,0.003713 ,0.534375 ,0.01665 ,0.000116 ,0.004815 ,0.485469 ,0.257122 ,0.010037 ,0.042815 ,0.223763 ,0.176829 ,0.000058 ,0.027035 ,0.370309 ,0.082323 ,0.000058 ,0.007948 ,0.13094 ,0.236179 ,0.000116 ,0.00963 ,0.571505 ,24.177765 ,1.496628 ,0.187826 ,0.004438 ,0.872084 ,5.709933 ,0.150107 ,0.516594 ,2.297807 ,1.094736 ,1.049877 ,0.009261 ,5.964323 ,0.403435 ,0.000096 ,0.246094 ,0.000386 ,2.066955 ,0.01447 ,0.002122 ,0.128111 ,3.121656 ,0.072352 ,0.064924 ,0.544474 ,0.557015 ,0 ,0 ,0.000579 ,0.496046 ,0.000096 ,0 ,0.000289 ,1.265487 ,0.007428 ,0.001061 ,0.348448 ,4.380679 ,0.084218 ,0.150589 ,1.825878 ,1.720919 ,0.000289 ,0.003376 ,0.011673 ,0.81372 ,0.000096 ,0.000772 ,0.002219 ,0.070809 ,0.341791 ,0.020909 ,1.118208 ,0.25165 ,0.000929 ,0.409721 ,3.82958 ,1.562591 ,0.487317 ,0.002044 ,10.892976 ,0.028715 ,0.083636 ,0 ,0.115975 ,0.042282 ,0 ,0.080383 ,1.144692 ,0.423382 ,0.002137 ,0.138742 ,2.487136 ,0.371063 ,0.196358 ,0.004089 ,0.685812 ,0.004832 ,0.001208 ,0.000186 ,0.342813 ,0.001022 ,0.000186 ,0.08568 ,1.023142 ,0.283432 ,0.003717 ,0.5585 ,3.45257 ,1.262618 ,0.771771 ,0.001951 ,0.880032 ,0.002974 ,0.001951 ,0.000279 ,0.396897 ,0.001394 ,0.000186 ,0.524116 ,0.047115 ,15.982573 ,0.000459 ,0.725836 ,0.768607 ,0.040706 ,0.349846 ,2.395833 ,4.919174 ,0.199748 ,0.000115 ,9.818281 ,0.514392 ,0.681575 ,0 ,0.000344 ,0.527005 ,0 ,0.002179 ,0.114551 ,1.519211 ,0.07373 ,0.151589 ,0.43252 ,2.075226 ,0.078776 ,0.003211 ,0.032106 ,1.903571 ,0.00172 ,0.000115 ,0.000344 ,0.266943 ,0 ,0.002179 ,0.670108 ,1.029128 ,0.117533 ,0.505677 ,1.978448 ,4.658653 ,0.250545 ,0.00172 ,0.356038 ,2.163175 ,0.007109 ,0.000115 ,0.001261 ,0.451899 ,0.000115 ,0.213967 ,0.058136 ,17.424222 ,35.359077 ,0.340991 ,0.00339 ,0.156468 ,1.773102 ,2.018483 ,0.196983 ,0.361164 ,4.693944 ,0.04789 ,0.094762 ,0.063994 ,12.169045 ,0.000424 ,0 ,0.093491 ,0.172403 ,0.508308 ,0.000593 ,0.049246 ,1.253945 ,0.487542 ,0.101204 ,0.021444 ,2.846677 ,0.006103 ,0 ,0.000085 ,0.660622 ,0.000678 ,0.000085 ,0 ,0.277929 ,0.330481 ,0.000254 ,0.04238 ,1.120532 ,1.93135 ,0.111968 ,0.103323 ,3.994417 ,0.012121 ,0.000085 ,0.002797 ,0.976185 ,0.001356 ,0.000085 ,0.000254 ,0.563995 ,0.069588 ,0.688169 ,35.921779 ,23.65509 ,11.510965 ,0.004316 ,0.325422 ,0.034794 ,0.000135 ,0.061227 ,0.806071 ,0.173298 ,0.001483 ,0.030883 ,1.372357 ,0.070533 ,0.181794 ,0.012138 ,0.309643 ,0.324613 ,0.000405 ,0.032232 ,0.483076 ,0.105597 ,0.001483 ,0.01025 ,0.215779 ,0.004585 ,0.00027 ,0.16615 ,0.516927 ,0.183143 ,0.077815 ,0.021578 ,0.674176 ,0.265948 ,0.000674 ,0.003641 ,0.071612 ,0.013216 ,0.000135 ,0.297371 ,1.514097 ,0.586111 ,0.011733 ,0.061497 ,0.309374 ,0.183952 ,0.045044 ,0.086986 ,1.273908 ,0.530278 ,0.004316 ,0.628438 ,0.068342 ,0.121937 ,1.323765 ,0.493179 ,0.138244 ,0.051005 ,0.053908 ,0.105226 ,0.074859 ,0.060421 ,0.037822 ,0.362681 ,0.035939 ,0.145245 ,0.167138 ,0.496706 ,0.139046 ,0.170041 ,0.136849 ,0.321642 ,0.078312 ,0.055163 ,0.369273 ,0.24749 ,0.195073 ,0.070308 ,0.032564 ,0.215788 ,0.061833 ,0.146501 ,0.176476 ,1.13591 ,0.019696 ,0.017106 ,0.113701 ,0.504866 ,0.010122 ,0.009966 ,0.066384 ,0.137634 ,0.028249 ,0.031466 ,0.044178 ,0.446015 ,0.03586 ,0.073682 ,0.105305 ,0.381671 ,0.039313 ,0.058616 ,0.092044 ,0.277072 ,0.038763 ,1.421996 ,1.140051 ,0.195833 ,0.081312 ,0.410046 ,0.069758 ,0.210115 ,0.099519 ,0.000932 ,0.011489 ,0.360038 ,0.415775 ,0.000776 ,0.053874 ,1.089585 ,0.144388 ,0.248099 ,0.122186 ,2.007768 ,0.069555 ,0.000311 ,0.298091 ,0.431456 ,0.163174 ,0.000155 ,0.05791 ,0.634065 ,0.002795 ,0 ,0.000155 ,0.239715 ,0.232728 ,0.012886 ,0.023909 ,0.726908 ,0.020649 ,0.000466 ,0.002639 ,0.54324 ,0.028878 ,0 ,0.014439 ,0.099985 ,0.806554 ,0.000466 ,0.028567 ,1.899865 ,0.197641 ,0.007608 ,0.105884 ,0.420899 ,0.754233 ,0.001087 ,0.126844 ,1.598824 ,0.05486 ,0.822999 ,0.75198 ,0.343463 ,0.143447 ,1.811753 ,56.838378 ,0.264556 ,0.61266 ,0.00284 ,0.000062 ,0.271628 ,0.871568 ,0.000062 ,0.000062 ,0.012473 ,0.423774 ,0.000432 ,0.000247 ,0.010929 ,5.159075 ,0.042112 ,1.418466 ,2.130663 ,0.576476 ,0.000123 ,0 ,0.026551 ,0.261069 ,0 ,0 ,0.001358 ,0.136647 ,0 ,0 ,0.000741 ,27.219895 ,1.373823 ,1.045943 ,12.035723 ,0.287373 ,0.000062 ,0 ,0.042112 ,1.041312 ,0 ,0 ,0.008892 ,0.201112 ,0 ,0 ,0.000741 ,2.299543 ,0.00599 ,0.024699 ,0.748196 ,0.033529 ,0.798582 ,1.184813 ,0.011609 ,0.010188 ,0.249707 ,0.040012 ,0.271072 ,0.645571 ,0.005996 ,0.143614 ,0.043609 ,0.000094 ,0.010539 ,0.371541 ,0.036489 ,0.000141 ,0.0126 ,0.166706 ,0.056771 ,0.000047 ,0.064968 ,1.006045 ,0.800602 ,0.059956 ,0.016628 ,0.761771 ,0.107218 ,0.001312 ,0.006558 ,0.285635 ,0.079161 ,0.000141 ,0.01602 ,0.200806 ,0.027355 ,0.000094 ,0.198558 ,2.30611 ,1.341144 ,0.207692 ,0.006464 ,0.055366 ,0.034428 ,0 ,0.016535 ,0.512203 ,0.074851 ,0.000234 ,0.018315 ,0.172889 ,0.032039 ,0.000141 ,0.047403 ,0.946136 ,0.371072 ,0.065905 ,4.483501 ,0.503397 ,0.056677 ,0.656004 ,0.198277 ,0.052602 ,0.56162 ,1.183337 ,0.076912 ,0.151858 ,0.004963 ,0.060699 ,0.33378 ,0.077878 ,0.011134 ,0.00299 ,0.586818 ,0.015652 ,0.001718 ,0.00579 ,0.506841 ,0.009289 ,1.065537 ,0.478019 ,3.062807 ,1.292935 ,0.007508 ,0.002163 ,0.250748 ,0.040275 ,0.005535 ,0.000382 ,0.112999 ,0.003499 ,0.001527 ,0.000064 ,0.127696 ,0.003054 ,16.560966 ,5.651603 ,12.455337 ,11.161511 ,0.005472 ,0.006808 ,0.159955 ,0.038557 ,0.019406 ,0.001654 ,0.220908 ,0.03563 ,0.001082 ,0.000191 ,0.175861 ,0.003881 ,0.762425 ,0.140867 ,1.108802 ,0.529619 ,0.177833 ,0.337279 ,0.611887 ,0.158873 ,0.694091 ,0.27499 ,0.240632 ,0.632947 ,0.395942 ,18.541946 ,0.675537 ,0.148268 ,0.000093 ,0.013122 ,0.263567 ,0.079948 ,0.000093 ,0.008032 ,0.497293 ,0.140516 ,0.00009 ,0.08331 ,0.225554 ,1.075187 ,0.056038 ,0.67937 ,1.119411 ,0.12632 ,0.000467 ,0.00976 ,0.638696 ,0.025217 ,0.00014 ,0.004063 ,0.324368 ,0.030121 ,0 ,0.000934 ,0.143832 ,1.036474 ,0.083684 ,0.07355 ,2.016257 ,0.084945 ,0 ,0.010974 ,0.064397 ,0.109975 ,0.000093 ,0.005884 ,0.483563 ,0.113991 ,0 ,0.018773 ,0.164659 ,0.565566 ,0.026338 ,0.068927 ,0.86397 ,0.481042 ,4.317981 ,0.27809 ,0.030074 ,0.034137 ,0.773935 ,0.045951 ,0.786871 ,0.733587 ,2.395822 ,16.011531 ,1.204356 };

const double dECMrestFreq[61] = {0.030001,0.020171,0.026344,0.023006,0.015528,0.020201,0.012142,0.013424,0.010372,0.011679,0.008195,0.010142,0.013551,0.023441,0.020102,0.025576,0.016076,0.011703,0.020211,0.011097,0.010642,0.0101,0.011843,0.010007,0.0048,0.014148,0.007837,0.008311,0.0076,0.017386,0.028839,0.014467,0.033223,0.024532,0.031878,0.028198,0.015908,0.028307,0.018853,0.019005,0.015796,0.022982,0.010191,0.016852,0.010901,0.018938,0.022747,0.019047,0.015782,0.015965,0.00975,0.011131,0.008956,0.01188,0.007029,0.01188,0.006025,0.016387,0.021383,0.015425,0.022103};
const double dECMunrestFreq[61] = {0.03109,0.020321,0.026699,0.022276,0.01312,0.017882,0.010682,0.012656,0.009746,0.013701,0.006788,0.01031,0.012814,0.023762,0.02118,0.024759,0.017168,0.01104,0.02073,0.010671,0.011165,0.010408,0.012199,0.010425,0.004809,0.014604,0.008158,0.008491,0.007359,0.016798,0.028497,0.015167,0.034527,0.025285,0.030606,0.02859,0.016516,0.026817,0.018288,0.018907,0.016386,0.023659,0.010223,0.016883,0.010921,0.018419,0.022626,0.01902,0.016697,0.017237,0.010366,0.010761,0.008721,0.011798,0.007415,0.012744,0.006441,0.016195,0.021349,0.015717,0.021414};

int GetFullCodonModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &out)	{
	int i, j, k;
	k = 0;
	FOR(i,61) {
		FOR(j,i) {
			cout << "\n[" << k << "]: " << dECMrest[k] << " : " << dECMunrest[k];
			k++;
		}
	}
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

