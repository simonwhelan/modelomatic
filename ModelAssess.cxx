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
		Error("ModelAssess <data_file> <tree_file> <output_file>");
	}
	///////////////////////////////////////////////////////////////////////////////////////////
	// Create the data structures
	PhyDat.SetIn(argv[1]); PhyDat.GetData();
	assert(PhyDat.pData()->m_DataType == DNA);
	CData NT_Data = *PhyDat.pData();
	CData AA_Data = *PhyDat.pData(); AA_Data.Translate(GeneticCode);
	CData COD_Data = *PhyDat.pData(); COD_Data.MakeCodonData();
	CData RY_Data = *PhyDat.pData(); RY_Data.DNA2RY();

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
		PhyDat.SetStartTree(InTree); PhyDat.pTree()->Unroot();
		Tree = *PhyDat.pTree();
		cout << " taken from file successfully" << flush;
	}
	cout << "\nStart topology: " << Tree << "\n";

	// Set output
	PhyDat.SetOut(argv[3]);

	// Set genetic code if required
	if(argc>4) { Error("Different genetic codes not implemented yet. Just live with the universal"); }



	///////////////////////////////////////////////////////////////////////////////////////////
	// Do the models
	vector <SModelDetails> Models;
	GetRYModels(PhyDat.pData(),&Tree,&Models,GeneticCode);
        cout<<"\rRY Done\n"<<flush;
	GetNTModels(PhyDat.pData(),&Tree,&Models,GeneticCode);
        cout<<"\rNT Done\n"<<flush;
	GetAAModels(PhyDat.pData(),&Tree,&Models,GeneticCode);
        cout<<"\rAA Done\n"<<flush;
	GetCODModels(PhyDat.pData(),&Tree,&Models,GeneticCode);
        cout<<"\rCodon Done\n"<<flush;

	cout << "\nModel Information\n---\nModel#\tName\tDataType\tTreeLength\tlnL\tNoPar\tlnL";
	FOR(i,(int)Models.size()) {
		cout << "\nModel["<<i<<"]\t" << Models[i].Name << "\t" << Models[i].DataType << "\t" << Models[i].TreeLength << "\t" << Models[i].lnL << "\t" << Models[i].NoPar << "\t" << Models[i].AIC;

	}

	cout << "\n\nSuccessful exit\n";
	return 0;
}

// Functions implementing the models
///////////////////////////////////////////////////
// ToDo:
// 1. Can probably be smarter in the order models are estimated and the way parameters are shared between them
void GetRYModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode)	{
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
}

void GetNTModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode)	{
	// Initialise
	CBaseModel *Model;
	// 1. JC
	CJC *JC; JC = new CJC(Data,Tree); Model = JC;
	Models->push_back(DoModelRun(Model,0));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,0));
        cout<<"."<<flush;
	Model = NULL;
	delete JC;
	// 2. FEL
	CFEL *FEL; FEL = new CFEL(Data,Tree); Model = FEL;
	Models->push_back(DoModelRun(Model,3));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,4));
        cout<<"."<<flush;
	Model = NULL;
	delete FEL;
	// 3. K2P
	CK2P *K2P; K2P= new CK2P(Data,Tree); Model = K2P;
	Models->push_back(DoModelRun(Model,1));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,2));
        cout<<"."<<flush;
	Model = NULL;
	delete K2P;
	// 4. HKY
	CHKY *HKY; HKY = new CHKY(Data,Tree); Model = HKY;
	Models->push_back(DoModelRun(Model,4));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,5));
        cout<<"."<<flush;
	Model = NULL;
	delete HKY;
	// 5. GTR
	CREV *REV; REV = new CREV(Data,Tree); Model = REV;
	Models->push_back(DoModelRun(Model,8));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,9));
        cout<<"."<<flush;
	Model = NULL;
	delete REV;
}

// Need something that chooses an appropriate set of models based on the genetic code
void GetAAModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models, int GeneticCode)	{
	// Initialise
	CData AmA = *Data; AmA.Translate(GeneticCode);
	CBaseModel *Model;
	CEMP *EMP;
	// Get the correction
	double AA2Cod_Adj = Data->GetAminoToCodonlnLScale(GeneticCode);
	// 1. EQU
	CEQU *EQU; EQU = new CEQU(&AmA,Tree,false); Model = EQU;
	Models->push_back(DoModelRun(Model,0,AA2Cod_Adj));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,1,AA2Cod_Adj));
        cout<<"."<<flush;
	Model = NULL;
	delete EQU;
	// 2. EQU+F
	EQU = new CEQU(&AmA,Tree,true); Model = EQU;
	Models->push_back(DoModelRun(Model,19,AA2Cod_Adj));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,20,AA2Cod_Adj));
        cout<<"."<<flush;
	Model = NULL;
	delete EQU;
	// 3. WAG
	EMP = new CEMP(&AmA,Tree,"WAG",false,(double*) dWAGVal,(double*) dWAGFreq); Model = EMP;
	Models->push_back(DoModelRun(Model,0,AA2Cod_Adj));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,1,AA2Cod_Adj));
        cout<<"."<<flush;
	Model = NULL;
	delete EMP;
	// 3. WAG+F
	EMP = new CEMP(&AmA,Tree,"WAG",true,(double*) dWAGVal,(double*) dWAGFreq); Model = EMP;
	Models->push_back(DoModelRun(Model,19,AA2Cod_Adj));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,20,AA2Cod_Adj));
        cout<<"."<<flush;
	Model = NULL;
	delete EMP;

}

// Note, running with codon data regularly into Codon models doesn't work quite right...
void GetCODModels(CData *Data, CTree *Tree, vector <SModelDetails> *Models,int GeneticCode) {
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
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,3));
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F1X4
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F1X4,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,6));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,7));
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F3X4
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F3X4,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,11));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,12));
        cout<<"."<<flush;
	Model = NULL;
	delete M0;
	// 1. F64
	CoD = *Data;
	M0 = new CCodonM0(&CoD,Tree,F64,GeneticCode); Model = M0;
	Models->push_back(DoModelRun(Model,2 + NoF64));
	Model->MakeGammaModel(0,4);
	Models->push_back(DoModelRun(Model,3 +NoF64));
        cout<<"."<<flush;
	Model = NULL;
	delete M0;



}

SModelDetails DoModelRun(CBaseModel *M, int NoPar,double Adj) {
	SModelDetails ModDet;
	ModDet.DataType = M->m_pData->m_DataType;
	M->lnL();
	ModDet.lnL = FullOpt(M) + Adj;
	ModDet.NoPar = NoPar; ModDet.AIC = GetAIC(ModDet.lnL,ModDet.NoPar);
	ModDet.TreeLength = M->Tree()->GetTreeLength();
	ModDet.Name = M->Name();
	return ModDet;
}
