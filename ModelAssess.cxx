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

// Note: this version was sent to David T @ Thu 14 March 2013

#include "./ModelAssess.h"
#include <time.h>
#include "ini/cpp/INIReader.h"
#include <set>

#define CHECK_LNL_OUT 1
#define VERSION_NUMBER "1.0a"

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter, MakeQ_Log_Counter, MakePT_Log_Counter, LFunc_Log_Counter, SubLFunc_Log_Counter, SPR_Log_Counter;
#endif
vector <double> GetPWVar(CBaseModel *Model, vector <double> *CurrentPW);
double GetDistVar(CBaseModel *M, CData *D,int Seq1, int Seq2, double CurrentDist);

/////////////////////////////////////////////////////////////////////////////////
// Main tree estimation routine
SBestModel DoTreeEstimation(CBaseModel *M, CData *D, CTree *T);

// Pointer to likelihood function that changes depending on the type of calculation performed
vector <STabuTree> TabuTrees;
vector <double> PWDists;
CPhyloDat PhyDat;
int TABU_RADIUS = DEFAULT_TABU_RADIUS, EXIT_OBS, OptObs;
double PROB_RAN_SEQ_REM = DEFAULT_PROB_RAN_SEQ_REM;
bool AllowPreOpt = true, DoItFast = false, DoItTrim = false, ModelOut = false;
int TrimTree = 10;		// Number of sequences defaulted to by the DoItTrim option
void DoInstructions();

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

		// Stuff from Leaphy
		int i,j,NumModelReruns = 1;
		long RandomSeed = 0;
		string temp_string, Name;
		vector <string> Toks;
		TabuTrees.clear();
		WarningMulD = false;
		SBestModel Result;
		bool DoSurface = false, DoHardOptCheck = false;
		string TreeName, DataName;	// Some files strings for Trim option

        //set up aa_model_map
        //

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
        string confFile="modelassess.ini";
        char const* home_c = getenv("HOME");
        string homedir = (home_c == NULL) ? std::string() : std::string(home_c);
        string dotConfFile=homedir.append("/.modelassess.ini");
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
                                if(count == 0) { cout << "\nSkipping " << model_names[i]; count++; }
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
                        if (reader.GetBoolean("Codon",*it,dotReader.GetBoolean("Codon",*it,true))){
                                new_codon_model_set.insert(*it);
                        }else {
                        	if(count == 0) { cout << "\nSkipping " << *it; count++; }
                        	else { cout << ", " << *it; }
                        }
                }
                codon_model_set=new_codon_model_set;



	string InTree = "file=";

	// Some initial verification
	if(min(LOOSE_RMSD_SUBSET,FULL_RMSD_SUBSET) < 105 && min(LOOSE_PARS_SUBSET,FULL_PARS_SUBSET) < 105) { Error("\nTrying to choose subsets of trees to examine based both on RMSD and parsimony"); }

	// Get information
	if(!InRange(argc,4,8)) {
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
		cout << "\n\t<models_out>  \tyes/no Option to output full model details to file model.out [default = no]";
		cout << "\nExiting...\n";
		exit(-1);
	}
	cout << "\n---------------------------------------------------------------------------------------------";
	cout << "\n  ModelOMatic (v" << VERSION_NUMBER << ").\n\tA program for choosing substitutions models for phylogenetic inference.\n\tWritten by Simon Whelan.\n\tContributions from James Allen, Ben Blackburne and David Talavera.";
	cout << "\n---------------------------------------------------------------------------------------------";

	///////////////////////////////////////////////////////////////////////////////////////////
	// Create the data structures
	// 0. Set output
	PhyDat.SetOut(argv[3]);
	// 1. Input the raw data
	cout << "\nData: <" << argv[1] << ">" << flush;
	PhyDat.SetIn(argv[1]); PhyDat.GetData();
	assert(PhyDat.pData()->m_DataType == DNA);
	cout << ": " << PhyDat.pData()->m_iNoSeq << " sequences of length " << PhyDat.pData()->m_iTrueSize << " (DataMatrix: " << PhyDat.pData()->m_iNoSeq << " x " << PhyDat.pData()->m_iSize << ")" << flush;
	// 2. Create a tree
	CTree Tree;
	clock_t start,end;
	cout << "\nCreating start tree ... " << flush;
    start = clock();
	if(!strcmp(argv[2],"bionj")) {
		// Create a bionj starting tree
                CData tmp_AA_Data = *PhyDat.pData();
                tmp_AA_Data.Translate();
		CEQU EQU_PW(&tmp_AA_Data,NULL);
		PWDists = GetPW(&EQU_PW,NULL,true);  // Get pairwise distances
		CTree T_bionj(DoBioNJ(PWDists,PhyDat.pData()->m_vsName,true),tmp_AA_Data.m_iNoSeq);
		cout << " estimated using bionj" << flush;
		Tree = T_bionj;
	} else {
		// Take tree from file
		InTree += argv[2];
		PhyDat.SetStartTree(InTree);
		if(PhyDat.pTree() == NULL) { Error("Failed to read tree from file <" + (string) argv[2] + ">\n\n"); }
		PhyDat.pTree()->Unroot();
		Tree = *PhyDat.pTree();
		cout << " taken from <" << InTree <<"> successfully" << flush;
	}
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)"<<flush;
        start=clock();

	// 3. If needed do the DoItTrim option
   	// Check whether if is meant to be running fast
    if(argc>5) {
    	if(strcmp(argv[5],"normal") == 0) { DoItFast = false; }
    	else {
    		TreeName = argv[5];
    		if(TreeName.find("fast") != string::npos) { DoItFast = true; TreeName = find_and_replace(TreeName,"fast",""); }
    		if(TreeName.find("trim") != string::npos) { DoItTrim = true; TreeName = find_and_replace(TreeName,"trim",""); }
    		if(!TreeName.empty()) { cout << "\nParsing of trim/fast option has resulted in left over characters: " << TreeName << "\n\n"; exit(-1); }
    	}
   	}
    cout << "\nOptimisation settings: "; if(!DoItFast && !DoItTrim) { cout << "normal"; } else { if(DoItFast) { cout << "fast "; } if(DoItTrim) { cout << "trim"; } }
    // 4. Check model output options
    if(argc>6) {
    	TreeName = argv[6];
    	if(TreeName.find("yes") != string::npos) { ModelOut = true; }
    	else if(TreeName.find("no") != string::npos) { ModelOut = false; }
    	else { cout << "\nUnknown option for models_out: " << TreeName << "\n\n"; exit(-1); }
    }
    // 5. Some file validation
    if(FileExist(argv[3])) { cout << "\nError: main output file <" << argv[3] << "> already exists. Please delete this and related files before continuing.\n"; exit(-1); }
    if(FileExist((string)argv[3] + (string)".model.out")) { cout << "\nError: models_out file <" << (string)argv[3] + (string)".model.out" << "> already exists. Please delete this and related files before continuing.\n"; exit(-1); }

	if(DoItTrim && Tree.NoSeq() > TrimTree)	{
		bool Check;
		// Build file names
		TreeName = (string)argv[3] + ".trim.tree";
		DataName = (string)argv[3] + ".trim.data";
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
	// 3. Create the other data sets
	CData NT_Data = *PhyDat.pData();
	CData AA_Data = *PhyDat.pData();
	CData AA_Temp = *PhyDat.pData();  AA_Temp.Translate(GeneticCode);	// Error check the translation for stop codons and so on
	CData COD_Data = *PhyDat.pData(); COD_Data.MakeCodonData();
	CData RY_Data = *PhyDat.pData();

	// Set output
	PhyDat.SetOut(argv[3]);

	// Set genetic code if required
	if(argc>4) {
		assert(InRange(atoi(argv[4]),0,NumGenCode));
		GeneticCode = atoi(argv[4]);
	}
	cout << "\nWorking with genetic code: " << GenCodeName[GeneticCode];

	cout << "\n>>> Doing model analysis <<< \n";
	///////////////////////////////////////////////////////////////////////////////////////////
	// Do the models

//	cout << "\n>>>>>>>>>>>>>>>>>> DOING DEBUG STUFFING <<<<<<<<<<<<<<<<<<<<<<<<";

	TreeName = (string)argv[3] + (string)".model.out";
	ofstream out(TreeName.c_str());
	vector <SModelDetails> Models;
	GetRYModels(&RY_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rRY Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
	GetNTModels(&NT_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rNT Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
	GetAAModels(&AA_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rAA Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
	GetCODModels(&COD_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rCodon Done ";
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();
    out.close();
    cout << "\nOutputting results to <" << argv[3] << ">";
    ofstream output(argv[3]);
	double minAIC = Models[0].AIC;
	for (i=1; i<(int)Models.size(); i++) {
		if (Models[i].AIC < minAIC) { minAIC = Models[i].AIC; };
	}
	output << "\nModel Information\n---\nModel#\tName\tDataType\tTreeLength\tOrilnL\tAdjlnL\tNoPar\tAIC\tDeltaAIC";
	FOR(i,(int)Models.size()) {
		output << "\nModel["<<i<<"]\t" << Models[i].Name << "\t" << Models[i].DataType << "\t" << Models[i].TreeLength << "\t" << Models[i].OrilnL << "\t" << Models[i].lnL << "\t" << Models[i].NoPar << "\t" << Models[i].AIC << "\t" << Models[i].AIC - minAIC;

	}
	output.close();
	cout << "\nSuccessful exit\n";
	return 0;
}

// Functions implementing the models
///////////////////////////////////////////////////
// ToDo:
// 1. Can probably be smarter in the order models are estimated and the way parameters are shared between them
void GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	SModelDetails RY_Model; RY_Model.DataType = RY;
	int df = -1;
	string NameAdd = "";
	// Make data
	CData RY_Data = *Data; RY_Data.DNA2RY();
	// 1. Get RYmodel
	CRY RY(&RY_Data,Tree);
	// Get the correction
	double RY2Cod_Adj = Data->GetRYToCodonlnLScale(GeneticCode,&df);

	if(df > 0) { NameAdd = "(+EMP)"; }
	// Optimise
	RY.m_sName += NameAdd;
	RY.lnL();
	RY_Model.OrilnL = FullOpt(&RY);
	RY_Model.lnL = RY_Model.OrilnL + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 1 + df;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	Models->push_back(RY_Model);
    cout<<"."<<flush;
    if(os != cout) { os << RY << endl << flush; }	// Output model details
	// 2. RY+dG
	RY.MakeGammaModel(0,4);
	RY.lnL();
	RY_Model.OrilnL = FullOpt(&RY,true,FlipBool(DoItFast),FlipBool(DoItFast));
	RY_Model.lnL = RY_Model.OrilnL + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 2 + df;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	Models->push_back(RY_Model);
	cout<<"."<<flush;
	if(os != cout) { os << RY << endl << flush; }	// Output model details
}

void GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	CBaseModel *Model;
	double Alpha;
	int count = 2;
	bool Fast = DoItFast;
	CTree PlainTree, GammaTree;
	// 1. JC
	DoItFast = false;
	CJC *JC; JC = new CJC(Data,Tree); Model = JC;
	Models->push_back(DoModelRun(Model,0));
	if(Fast) { PlainTree = *Model->Tree(); }
	if(os != cout) { os << *JC << endl << flush; }	// Output model details
//	cout << "\nTree JC:   \t" << Tree->TreeLength() << "\t" << JC->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,1));
	Alpha = JC->m_vpPar[0]->Val();
	if(Fast) { GammaTree = *Model->Tree(); }
	if(os != cout) { os << *JC << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
//	cout << "\nTree JCdG: \t" << Tree->TreeLength() << "\t" << JC->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	delete JC;
	DoItFast = Fast;
	// 2. FEL
	CFEL *FEL; FEL = new CFEL(Data,Tree); Model = FEL;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,3));
	if(os != cout) { os << *FEL << endl << flush; }	// Output model details
//	delete FEL; FEL = new CFEL(Data,Tree); Model = FEL;
//	cout << "\nTree FEL:  \t" << Tree->TreeLength() << "\t" << FEL->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	Model->MakeGammaModel(0,4,Alpha);
	if(DoItFast) { *Model->Tree() = GammaTree; }
	Models->push_back(DoModelRun(Model,4));
	if(os != cout) { os << *FEL << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
//	cout << "\nTree FELdG:\t" << Tree->TreeLength() << "\t" << FEL->lnL(true) << " cf. " << Models->at(count++).OrilnL;
	delete FEL;
	// 3. K2P
	CK2P *K2P; K2P= new CK2P(Data,Tree); Model = K2P;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,1));
	if(os != cout) { os << *K2P << endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	if(DoItFast) { *Model->Tree() = GammaTree; }
	Models->push_back(DoModelRun(Model,2));
	if(os != cout) { os << *K2P << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete K2P;
	// 4. HKY
	CHKY *HKY; HKY = new CHKY(Data,Tree); Model = HKY;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,4));
	if(os != cout) { os << *HKY<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	if(DoItFast) { *Model->Tree() = GammaTree; }
	Models->push_back(DoModelRun(Model,5));
	if(os != cout) { os << *HKY << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete HKY;
	// 5. GTR
	CREV *REV; REV = new CREV(Data,Tree); Model = REV;
	if(DoItFast) { *Model->Tree() = PlainTree; }
	Models->push_back(DoModelRun(Model,8));
	if(os != cout) { os << *REV<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	if(DoItFast) { *Model->Tree() = GammaTree; }
	Models->push_back(DoModelRun(Model,9));
	if(os != cout) { os << *REV<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete REV;
}

// Need something that chooses an appropriate set of models based on the genetic code
void GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	string NameAdd;
	double Alpha = 1.0;
	vector <double> ModelF(20,0.05);
	CData AmA = *Data; AmA.Translate(GeneticCode);
	CBaseModel *Model;
	CEMP *EMP;
	bool First = true, Fast = DoItFast;
	CTree PlainTree, GammaTree;
     // The correction information
	int df = -1;
	double AA2Cod_Adj = Data->GetAminoToCodonlnLScale(GeneticCode,&df);
	if(df > 0) { NameAdd = "(+EMP)"; }
	// Do the analyses and apply the correction ad hoc

    if (aa_model_map.find("EQU")!=aa_model_map.end()){
                // 1. EQU
    			if(First) { DoItFast = false; }
                CEQU *EQU; EQU = new CEQU(&AmA,Tree,false); Model = EQU; Model->m_sName += NameAdd;
                Models->push_back(DoModelRun(Model,0+df,AA2Cod_Adj));
                if(First) { PlainTree = *Model->Tree(); }
                if(os != cout) { os << *EQU << endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4);
                Models->push_back(DoModelRun(Model,1+df,AA2Cod_Adj));
                if(First) { GammaTree = *Model->Tree(); First = false; DoItFast = Fast;}
                Alpha = Model->m_vpPar[0]->Val();
                if(os != cout) { os << *EQU << endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete EQU;
                // 2. EQU+F
                EQU = new CEQU(&AmA,Tree,true); Model = EQU; Model->m_sName += NameAdd;
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,19+df,AA2Cod_Adj));
                if(os != cout) { os << *EQU << endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                Models->push_back(DoModelRun(Model,20+df,AA2Cod_Adj));
                if(os != cout) { os << *EQU << endl << flush; }	// Output model details
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
                if(First) { DoItFast = false; }
                EMP = new CEMP(&AmA,Tree,name,false,mySMat,myFreq); Model = EMP; Model->m_sName += NameAdd;
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,0+df,AA2Cod_Adj));
                if(First) { PlainTree = *Model->Tree(); }
                if(os != cout) { os << *EMP << endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                Models->push_back(DoModelRun(Model,1+df,AA2Cod_Adj));
                if(First) { GammaTree = *Model->Tree(); First = false; DoItFast = Fast; }
                if(os != cout) { os << *EMP << endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
                //+F
                EMP = new CEMP(&AmA,Tree,name,true,mySMat,myFreq); Model = EMP; Model->m_sName += NameAdd;
                if(os != cout) { os << *EMP << endl << flush; }	// Output model details
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,19+df,AA2Cod_Adj));
                if(os != cout) { os << *EMP << endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                Models->push_back(DoModelRun(Model,20+df,AA2Cod_Adj));
                if(os != cout) { os << *EMP << endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
        }

}

// Note, running with codon data regularly into Codon models doesn't work quite right...
void GetCODModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models,int GeneticCode, ostream &os) {
	// Initialise
	CData CoD = *Data;
	CBaseModel *Model;
	CCodonM0 *M0;
	int i,NoF64 = 0;
	bool First = true, Fast = DoItFast;
	CTree PlainTree, GammaTree;
	// Do the calculations
	FOR(i,64) { if(GenCodes[GeneticCode][i] >= 0) { NoF64++; } }
        if (codon_model_set.count("F0")){
                // 1. F0
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,cEQU,GeneticCode); Model = M0;
                if(First) { DoItFast = false; }
                Models->push_back(DoModelRun(Model,2));
                if(First) { PlainTree = *Model->Tree(); }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4);
                Models->push_back(DoModelRun(Model,3));
                if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }
        if (codon_model_set.count("F1X4")){
                        // 1. F1X4
                        CoD = *Data;
                        M0 = new CCodonM0(&CoD,Tree,F1X4,GeneticCode); Model = M0;
                        if(First) { DoItFast = false; }
                        if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                        Models->push_back(DoModelRun(Model,6));
                        if(First) { PlainTree = *Model->Tree(); }
                        if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                        Model->MakeGammaModel(0,4);
                        if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                        Models->push_back(DoModelRun(Model,7));
                        if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                        if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                        cout<<"."<<flush;
                        Model = NULL;
                        delete M0;
        }
        if (codon_model_set.count("F3X4")){
                // 1. F3X4
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,F3X4,GeneticCode); Model = M0;
                if(First) { DoItFast = false; }
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,11));
                if(First) { PlainTree = *Model->Tree(); }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4);
                if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                Models->push_back(DoModelRun(Model,12));
                if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }
        if (codon_model_set.count("F64")){
                // 1. F64
                CoD = *Data;
                M0 = new CCodonM0(&CoD,Tree,F64,GeneticCode); Model = M0;
                if(First) { DoItFast = false; }
                if(!First && DoItFast) { *Model->Tree() = PlainTree; }
                Models->push_back(DoModelRun(Model,2 + NoF64));
                if(First) { PlainTree = *Model->Tree(); }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4);
                if(!First && DoItFast) { *Model->Tree() = GammaTree; }
                Models->push_back(DoModelRun(Model,3 +NoF64));
                if(First) { GammaTree = *Model->Tree(); DoItFast = Fast; First = false; }
                if(os != cout) { os << *Model<< endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete M0;
        }

}

SModelDetails DoModelRun(CBaseModel *M, int NoPar,double Adj) {
	SModelDetails ModDet;
	ModDet.DataType = M->m_pData->m_DataType;
	M->lnL();
	ModDet.OrilnL = FullOpt(M,true,FlipBool(DoItFast)) ;
//	cout << "\nFinished opt: " << ModDet.OrilnL << " cf. actual lnL calc: " << M->lnL(true);
//	double FullOpt(CBaseModel *Model, bool DoPar, bool DoBra, bool DoFreq, double CurlnL,bool FullLikAcc, int NoIterations, double lnL2Beat, double lnL_tol, bool DoOutput,bool TightFullOpt)	{
//	ModDet.OrilnL = FullOpt(M,true,FlipBool(DoItFast),false,-BIG_NUMBER,true, DEFAULT_OPTNUM,-BIG_NUMBER, FULL_LIK_ACC,true,false);
	ModDet.lnL = ModDet.OrilnL + Adj;
	ModDet.NoPar = NoPar; ModDet.AIC = GetAIC(ModDet.lnL,ModDet.NoPar);
	ModDet.TreeLength = M->Tree()->GetTreeLength();
	ModDet.Name = M->Name();
#if CHECK_LNL_OUT == 1
	if(fabs(ModDet.OrilnL - M->lnL(true)) > 0.001)	{
			cout << "\nModel: " << M->Name() << " obtained " << ModDet.lnL << " cf. " << M->lnL() << " cff. " << M->lnL(true) << "\n\n"; exit(-1);
	}
#endif
	return ModDet;
}
