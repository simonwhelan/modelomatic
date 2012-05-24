/*	///////////////////////////////////////////////////////////////////////
		ModelAssess -- Program for assessing the fit of lots of models
			Simon Whelan, University of Manchester

		ModelAssess <data_file> <tree_file> <output_file> <Genetic_code; default=universal>

	Option definitions
	<data_file>		Any sequence data file in phylip or sequential format
	<tree_file>		The tree
	<output_file>	Output file
	/////////////////////////////////////////////////////////////////////// */

#include "./ModelAssess.h"

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
bool AllowPreOpt = true;

void DoInstructions();

bool WarningMulD;

// File specific global variables
int DebugOutput = 0;

int main(int argc, char *argv[])	{
	string InTree = "file=";
	int GeneticCode = 0;

	// Stuff from Leaphy
	int i,NumModelReruns = 1;
	long RandomSeed = 0;
	string temp_string, Name;
	vector <string> Toks;
	TabuTrees.clear();
	WarningMulD = false;
	SBestModel Result;
	bool DoSurface = false, DoHardOptCheck = false;

	// Some initial verification
	if(min(LOOSE_RMSD_SUBSET,FULL_RMSD_SUBSET) < 105 && min(LOOSE_PARS_SUBSET,FULL_PARS_SUBSET) < 105) { Error("\nTrying to choose subsets of trees to examine based both on RMSD and parsimony"); }

	// Get information
	if(argc != 4) {
		Error("ModelAssess <data_file> <tree_file> <output_file>\n");
	}



	///////////////////////////////////////////////////////////////////////////////////////////
	// Create the data structures
	PhyDat.SetIn(argv[1]); PhyDat.GetData();
	assert(PhyDat.pData()->m_DataType == DNA);
	CData NT_Data = *PhyDat.pData();
	CData AA_Data = *PhyDat.pData();
	CData COD_Data = *PhyDat.pData(); COD_Data.MakeCodonData();
	CData RY_Data = *PhyDat.pData();

	///////////////////////////////////////////////////////////////////////////////////////////
	// Get tree
	CTree Tree;
	cout << "\nCreating start tree ... " << flush;
	if(!strcmp(argv[2],"bionj")) {
		// Create a bionj starting tree
		CEQU EQU_PW(&AA_Data,NULL);
		PWDists = GetPW(&EQU_PW,NULL,true);  // Get pairwise distances
		CTree T_bionj(DoBioNJ(PWDists,PhyDat.pData()->m_vsName,true),AA_Data.m_iNoSeq);
		cout << " estimated using bionj" << flush;
		Tree = T_bionj;

	} else {
		// Take tree from file
		InTree += argv[2];
		PhyDat.SetStartTree(InTree);
		if(PhyDat.pTree() == NULL) { Error("Failed to read tree from file <" + (string) argv[2] + ">\n\n"); }
		PhyDat.pTree()->Unroot();
		Tree = *PhyDat.pTree();
		cout << " taken from file successfully" << flush;
	}

	// Set output
	PhyDat.SetOut(argv[3]);

	// Set genetic code if required
	if(argc>4) { Error("Different genetic codes not implemented yet. Just live with the universal"); }

/*
	// Some debug code for Codon models
	cout << "\nTrying M0 Codon model...";
	CCodonM0 TestM0(&COD_Data, &Tree, F64, 0);
	TestM0.m_vpPar[1]->SetVal(1.0);
	TestM0.m_vpPar[0]->SetVal(1.0);
	cout << "\nlnL: " << flush;
	cout << TestM0.lnL() << flush;
	cout << " ... done" << flush;
//	TestM0.m_vpProc[0]->OutQ();
	cout << "\n" << TestM0;

	exit(-1);
	FullOpt(&TestM0,true,true,false,-BIG_NUMBER,true,DEFAULT_OPTNUM,-BIG_NUMBER,FULL_LIK_ACC,true);

	exit(-1);
*/

	cout << "\nDoing model analysis...\n";
	///////////////////////////////////////////////////////////////////////////////////////////
	// Do the models
	ofstream out("model.output");
	vector <SModelDetails> Models;
	GetRYModels(&RY_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rRY Done\n"<<flush;
	GetNTModels(&NT_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rNT Done\n"<<flush;
	GetAAModels(&AA_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rAA Done\n"<<flush;
	GetCODModels(&COD_Data,&Tree,&Models,GeneticCode, out);
        cout<<"\rCodon Done\n"<<flush;
    out.close();

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
	cout << "\n\nSuccessful exit\n";
	return 0;
}

// Functions implementing the models
///////////////////////////////////////////////////
// ToDo:
// 1. Can probably be smarter in the order models are estimated and the way parameters are shared between them
void GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	SModelDetails RY_Model; RY_Model.DataType = RY;
	// Make data
	CData RY_Data = *Data; RY_Data.DNA2RY();
	// 1. Get RYmodel
	CRY RY(&RY_Data,Tree);
	// Get the correction
	double RY2Cod_Adj = Data->GetRYToCodonlnLScale(GeneticCode);
	// Optimise
	RY.lnL();
	RY_Model.lnL = FullOpt(&RY) + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 1;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	Models->push_back(RY_Model);
    cout<<"."<<flush;
    if(os != cout) { os << RY << endl << flush; }	// Output model details
	// 2. RY+dG
	RY.MakeGammaModel(0,4);
	RY.lnL();
	RY_Model.lnL = FullOpt(&RY) + RY2Cod_Adj;
	RY_Model.Name = RY.Name();
	RY_Model.TreeLength = RY.Tree()->GetTreeLength();
	RY_Model.NoPar = 2;
	RY_Model.AIC = GetAIC(RY_Model.lnL,RY_Model.NoPar);
	Models->push_back(RY_Model);
	cout<<"."<<flush;
	if(os != cout) { os << RY << endl << flush; }	// Output model details
}

void GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	CBaseModel *Model;
	double Alpha;
	// 1. JC
	CJC *JC; JC = new CJC(Data,Tree); Model = JC;
	Models->push_back(DoModelRun(Model,0));
	if(os != cout) { os << *JC << endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,1));
	Alpha = JC->m_vpPar[0]->Val();
	if(os != cout) { os << *JC << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete JC;
	// 2. FEL
	CFEL *FEL; FEL = new CFEL(Data,Tree); Model = FEL;
	Models->push_back(DoModelRun(Model,3));
	if(os != cout) { os << *FEL << endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	Models->push_back(DoModelRun(Model,4));
	if(os != cout) { os << *FEL << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete FEL;
	// 3. K2P
	CK2P *K2P; K2P= new CK2P(Data,Tree); Model = K2P;
	Models->push_back(DoModelRun(Model,1));
	if(os != cout) { os << *K2P << endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	Models->push_back(DoModelRun(Model,2));
	if(os != cout) { os << *K2P << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete K2P;
	// 4. HKY
	CHKY *HKY; HKY = new CHKY(Data,Tree); Model = HKY;
	Models->push_back(DoModelRun(Model,4));
	if(os != cout) { os << *HKY<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	Models->push_back(DoModelRun(Model,5));
	if(os != cout) { os << *HKY << endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete HKY;
	// 5. GTR
	CREV *REV; REV = new CREV(Data,Tree); Model = REV;
	Models->push_back(DoModelRun(Model,8));
	if(os != cout) { os << *REV<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	Models->push_back(DoModelRun(Model,9));
	if(os != cout) { os << *REV<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete REV;
}

// Need something that chooses an appropriate set of models based on the genetic code
void GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode, ostream &os)	{
	// Initialise
	double Alpha;
	CData AmA = *Data; AmA.Translate(GeneticCode);
	CBaseModel *Model;
	CEMP *EMP;
        int n_models=3;
        string model_names[] = {"JTT","WAG","LG"};
        double* smat[] = {(double*)dJTTVal,(double*)dWAGVal,(double*)dLGVal};
        double* freq[] = {(double*)dJTTFreq,(double*)dWAGFreq,(double*)dLGFreq};

	// Get the correction
	double AA2Cod_Adj = Data->GetAminoToCodonlnLScale(GeneticCode);
	// 1. EQU
	CEQU *EQU; EQU = new CEQU(&AmA,Tree,false); Model = EQU;
	Models->push_back(DoModelRun(Model,0,AA2Cod_Adj));
	if(os != cout) { os << *EQU<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,1,AA2Cod_Adj));
	Alpha = Model->m_vpPar[0]->Val();
	if(os != cout) { os << *EQU<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete EQU;
	// 2. EQU+F
	EQU = new CEQU(&AmA,Tree,true); Model = EQU;
	Models->push_back(DoModelRun(Model,19,AA2Cod_Adj));
	if(os != cout) { os << *EQU<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4,Alpha);
	Models->push_back(DoModelRun(Model,20,AA2Cod_Adj));
	if(os != cout) { os << *EQU<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete EQU;
        for (int i=0; i < n_models; i++){
                //model frequencies
                EMP = new CEMP(&AmA,Tree,model_names[i],false,smat[i],freq[i]); Model = EMP;
                Models->push_back(DoModelRun(Model,0,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                Models->push_back(DoModelRun(Model,1,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
                //+F
                EMP = new CEMP(&AmA,Tree,model_names[i],true,smat[i],freq[i]); Model = EMP;
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                Models->push_back(DoModelRun(Model,19,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                Models->push_back(DoModelRun(Model,20,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
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
	FOR(i,64) { if(GenCodes[GeneticCode][i] >= 0) { NoF64++; } }
	// 1. F0
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,cEQU,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,2));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,3));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F1X4
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F1X4,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,6));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,7));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F3X4
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F3X4,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,11));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,12));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F64
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F64,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,2 + NoF64));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,3 +NoF64));
	if(os != cout) { os << *Model<< endl << flush; }	// Output model details
        cout<<"."<<flush;
	Model = NULL;
	delete M0;



}

SModelDetails DoModelRun(CBaseModel *M, int NoPar,double Adj) {
	SModelDetails ModDet;
	ModDet.DataType = M->m_pData->m_DataType;
	M->lnL();
	ModDet.OrilnL = FullOpt(M) ;
	ModDet.lnL = ModDet.OrilnL + Adj;
	ModDet.NoPar = NoPar; ModDet.AIC = GetAIC(ModDet.lnL,ModDet.NoPar);
	ModDet.TreeLength = M->Tree()->GetTreeLength();
	ModDet.Name = M->Name();
	return ModDet;
}
