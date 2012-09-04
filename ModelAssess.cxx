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
#include <time.h>
#include "ini/cpp/INIReader.h"
#include <set>
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

int n_models=13;
string model_names[] = {"JTT","WAG","LG","DAY","mtREV","mtMam","mtArt","rtREV","cpREV","BLOSUM62","VT","HIVb","HIVw"};
double* smat[] = {(double*)dJTTVal,(double*)dWAGVal,(double*)dLGVal,(double*)dDAYVal,(double*)dmtREVVal,(double*)dmtMAMVal,(double*)dmtArtVal,(double*)drtREVVal,(double*)dcpREVVal,(double*)dBLOSUM62Val,(double*)dVTVal,(double*)dHIVbVal,(double*)dHIVwVal};
double* freq[] = {(double*)dJTTFreq,(double*)dWAGFreq,(double*)dLGFreq,(double*)dDAYFreq,(double*)dmtREVFreq,(double*)dmtMAMFreq,(double*)dmtArtFreq,(double*)drtREVFreq,(double*)dcpREVFreq,(double*)dBLOSUM62Freq,(double*)dVTFreq,(double*)dHIVbFreq,(double*)dHIVwFreq};

std::map <string,double**> aa_model_map;
std::set <string> codon_model_set;



int main(int argc, char *argv[])	{

        //set up aa_model_map
        //

        for (int i=0; i < n_models; i++){
                double* vals[] = {smat[i],freq[i]};
                aa_model_map[model_names[i]]=vals;
        }
        codon_model_set.insert("F0");
        codon_model_set.insert("F1X4");
        codon_model_set.insert("F3X4");
        codon_model_set.insert("F64");

        //read conf
        string confFile="modelassess.ini";
        ifstream ifile("modelassess.ini");
        if (ifile){
                INIReader reader(confFile);
                if (reader.ParseError() < 0) {
                        std::cout << "Can't read " << confFile << "\n" ;
                        return 1;
                }
                for (int i=0; i < n_models; i++){
                        if (!reader.GetBoolean("Amino Acid",model_names[i],true)){
                                aa_model_map.erase(model_names[i]);
                                cout << "Skipping " << model_names[i] << "\n";
                        }
                }
                std::set <string> new_codon_model_set;
                std::set<string>::iterator it;
                for( it = codon_model_set.begin(); it != codon_model_set.end(); it++ ) {
                        if (reader.GetBoolean("Codon",*it,true)){
                                new_codon_model_set.insert(*it);
                        }else {
                                cout << "Skipping " << *it << "\n";
                        }
                }
                codon_model_set=new_codon_model_set;
        }



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
		cout << " taken from file successfully" << flush;
	}
        end = clock();
        cout << " (" << (double)(end-start)/CLOCKS_PER_SEC << "s)\n"<<flush;
        start=clock();

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
        for(std::map<string,double**>::iterator iter = aa_model_map.begin(); iter != aa_model_map.end(); ++iter){
                string name=iter->first;
                double* mySMat=iter->second[0];
                double* myFreq=iter->second[1];
                //model frequencies
                EMP = new CEMP(&AmA,Tree,name,false,mySMat,myFreq); Model = EMP;
                Models->push_back(DoModelRun(Model,0,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                Model->MakeGammaModel(0,4,Alpha);
                Models->push_back(DoModelRun(Model,1,AA2Cod_Adj));
                if(os != cout) { os << *EMP<< endl << flush; }	// Output model details
                cout<<"."<<flush;
                Model = NULL;
                delete EMP;
                //+F
                EMP = new CEMP(&AmA,Tree,name,true,mySMat,myFreq); Model = EMP;
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
        if (codon_model_set.count("F0")){
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
        }
        if (codon_model_set.count("F1X4")){
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
        }
        if (codon_model_set.count("F3X4")){
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
        }
        if (codon_model_set.count("F64")){
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
