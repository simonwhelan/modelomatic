////////////////////////////////////////////////////////////////
// Implementation of various tree carrying and data managing
// routine

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

////////////////////////////////////////////////////////////////
// General data structure
int CPhyloDat::SimNum = 0;
// Constructor
CPhyloDat::CPhyloDat()	{
	m_pData = NULL; m_pModel = NULL; m_pTree = NULL;
	m_sInFile = "\0"; m_sOutFile = "\0"; m_sTreeFile = "\0";
	m_bDoGamma = DEFAULT_GAMMA; m_bDoInv = DEFAULT_INV; m_iGamCat = DEFAULT_GAMMA_CAT; m_bDoGarbage = DEFAULT_GARBAGE;
	m_bDoTreeHMM = false;
	m_dModelOptTol = 0.01; m_bF = true; m_bDoSA = m_bDoBioNJ = true;
	// Some THMM specific stuff
	THMM_NoCat = 1; THMM_R = same; THMM_F = false; THMM_K = false; THMM_H = H_none;
	// Some codon stuff
	m_iGenCode = 0; m_CodonEqm = cEQU;
	// The optimise stuff
	m_bOptPar = m_bOptBra = true; m_bOptFre = false; m_bAllowCompress = true;
	// The model
	m_Model = UNKNOWN; m_iPseudoDataSize = m_iPseudoDataSeq = -1; m_PseudoDataType = NONE;
}
// Destructor
CPhyloDat::~CPhyloDat()	{
	if(m_vpMultiDat.empty()) { CleanData(); } else { m_pModel->MakeNewData(NULL); m_pData = NULL; CleanMultiData(); }
	CleanModel();
	CleanTree();
}
// Whether I need to run the menu
bool CPhyloDat::NeedMenu()	{
	if(m_pData != NULL) { if(!m_pData->Valid()) { return true; CleanData(); } } else if(m_vpMultiDat.empty()) { return true; }
	if(m_pTree != NULL) { if(!m_pTree->Valid()) { return true; CleanTree(); } } else { if(!m_bDoBioNJ && ! m_bDoSA) { return true; } }
	if(m_Model == UNKNOWN) { return true; }
	return false;
}

// Interaction functions
bool CPhyloDat::DoF() { if(DataType() != AA) { return false; } return m_bF; }
void CPhyloDat::SetIn(string S) { if(S.size() > 0) { m_sInFile = S; } else { m_sInFile = "\0"; } }
void CPhyloDat::SetOut(string S) { if(S.size() > 0) { m_sOutFile = S; } else { m_sOutFile = "\0"; } }
void CPhyloDat::SetTreeFile(string S){ if(S.size() > 0) { m_sTreeFile = S; } else { m_sOutFile = "\0"; } }
// Data interaction functions
EDataType CPhyloDat::DataType() { if(m_pData != NULL) { return m_pData->m_DataType; } return NONE; }
int CPhyloDat::Size() { if(m_pData != NULL) { return m_pData->m_iSize; } return -1; }
int CPhyloDat::TrueSize() { if(m_pData != NULL) { return m_pData->m_iTrueSize; } return -1; }
int CPhyloDat::NoSeq() { if(m_pData != NULL) { return m_pData->m_iNoSeq; } return -1; }
void CPhyloDat::CleanTree() { if(m_pTree != NULL) { delete m_pTree; m_pTree = NULL; } m_sTreeFile = ""; }
void CPhyloDat::CleanData() { if(m_pData != NULL) { delete m_pData; m_pData = NULL; } m_sInFile = "";  }
void CPhyloDat::CleanMultiData() {
	int i;
	if(m_vpMultiDat.size() != m_vpMultiTree.size()) { Error("\nUnexpected condition where number of data sets not equal to number of trees...\n"); }
	FOR(i,(int)m_vpMultiDat.size()) { if(m_vpMultiDat[i] != NULL) { delete m_vpMultiDat[i]; } m_vpMultiDat[i] = NULL; }
	FOR(i,(int)m_vpMultiTree.size()) { if(m_vpMultiTree[i] != NULL) { delete m_vpMultiTree[i]; } m_vpMultiTree[i] = NULL; }
	m_vpMultiDat.clear();
	m_vpMultiTree.clear();
	m_vsMultiName.clear();
	m_vdMultiScore.clear();
}
void CPhyloDat::CleanModel() { if(m_pModel != NULL) { delete m_pModel; m_pModel = NULL; } }
void CPhyloDat::GetOutputFile() {
	bool GetFileName = true;
	string SugName;
	SugName = In() + ".output";
	if(!Out().empty()) {
		if(FileExist(Out())) {
			cout << "\nWarning: file <"<<Out()<<"> exists already...";
			if(GetYesNoOption("Do you want to delete it [y/n]: ")) { ofstream out(Out().c_str()); out.close(); GetFileName = false; }
		} else { GetFileName = false; }
	}
	if(GetFileName) { SetOut(GetOutFileName(SugName)); }
}
bool CPhyloDat::IsTree()	{ if(m_pTree == NULL) { return false; } else { return true; } }
// Get the data
void CPhyloDat::GetData() {
	int i;
	bool SingleData = true;
	string store;
	// If required get the input file
	if(In().empty() || !FileExist(In())) { SetIn(GetInFileName()); }
	// Check what kind of data it is
	FINOPEN(incheck,In().c_str());
	getline(incheck,store); incheck.close();
	FOR(i,(int)store.size()) { store[i] = toupper(store[i]); }
	if(store.find("MULTIGENE") == string::npos) {
//		cout << "\nGetting single data";
		GetSingleData();
	} else {
//		cout << "\nGetting multidata";
		GetMultiData();
	}
}
// Get data for a single data set
void CPhyloDat::GetSingleData()	{
	// If required get the input file
	if(In().empty() || !FileExist(In())) { SetIn(GetInFileName()); }
	while(m_pData == NULL)	{
//		cout << "\nGetting data" << flush;
		m_pData = new CData(In(),NONE,true);
//		cout << " ... done" << flush;
		if(!m_pData->Valid()) {
			cout << "\nFile does not appear to contain sequence data...";
			delete m_pData; m_pData = NULL; SetIn(GetInFileName());
	}	}
	m_iPseudoDataSeq = m_pData->m_iNoSeq;
	m_iPseudoDataSize = m_pData->m_iTrueSize;
	m_PseudoDataType = m_pData->m_DataType;
}
// Get data for a multiple data set
void CPhyloDat::GetMultiData() {
	int i, GeneCount = 0;
	vector <string > Toks;
	CData *Input;
	// Initialise file pointer
	string store; FINOPEN(in,In().c_str()); getline(in, store);
	FOR(i,(int)store.size()) { store[i] = toupper(store[i]); }
	assert(store.find("MULTIGENE") != string::npos);
	for(;;) {
		getline(in,store);
		if(in.eof()) { break; };
		if(store[0] == '#') { continue; }
		Toks = Tokenise(store);
		if((int)Toks.size() != 1) { cout << "\nLine in gene file <"<<In()<<"> has multiple tokens: "<< store << "\n\n"; Error(); }
		cout << "\n\tGetting gene[" << GeneCount+1 << "] <"<<Toks[0] << ">";
		Input = new CData(Toks[0]);
		if(!m_vpMultiDat.empty()) { if(Input->m_DataType != m_vpMultiDat[0]->m_DataType) { Error("\nNot all data sets of the same type?\n\n"); } }
		cout << "\t" << Input->m_iNoSeq << ":" << Input->m_iSize;
		m_vpMultiDat.push_back(Input); Input = NULL;
		m_vpMultiTree.push_back(NULL);
		m_vdMultiScore.push_back(-BIG_NUMBER);
		m_vsMultiName.push_back(Toks[0]);
		GeneCount ++;
		cout << " ... done";
	}
	if(GeneCount == 0) { Error("\nTried to input multiple data sets, but ended up with no data\n\n"); }
}

// Get the tree
void CPhyloDat::GetTree()	{
	string temp;
	if(m_pData == NULL) { Error("Cannot get tree before datafile is read...\n"); }
	// If required get the tree input file
	if(TreeFile().empty() || !FileExist(TreeFile())) { SetTreeFile(GetInFileName()); }
	m_pTree = new CTree(TreeFile(),true,m_pData,true);
	if(!m_pTree->Valid()) { cout << "\nInvalid tree... please try again."; CleanTree(); }
	else {
		m_pTree->OutBra();
		m_pTree->NoOutBra(); m_pTree->NoOutName();
	}
}
// Create a tree
CTree *CPhyloDat::CreateTree(string T)	{ assert(m_pData != NULL); CleanTree(); m_pTree = new CTree(T,NoSeq(),false); if(m_pModel != NULL) { m_pModel->MakeNewTree(m_pTree,false); } return m_pTree; }
CTree *CPhyloDat::CreateTree()			{ CleanTree(); m_pTree = new CTree; if(m_pModel != NULL) { m_pModel->MakeNewTree(m_pTree,false); } return m_pTree; }
// Get the model from a string
// ---
// Note: when adding models that have similar names make sure the more complex goes first.
bool CPhyloDat::SetModel(string s)	{
	vector <string> Toks = Tokenise(s,".");
	m_Model = UNKNOWN;
	m_bDoGamma = false; m_bDoInv = false;
	// Do rate variation
	if(s.find("+dG",0) != string::npos || wildcmp("*+?dG*",s.c_str())) {
		if(s.find("+dG",0) != string::npos) { m_iGamCat = DEFAULT_GAMMA_CAT; s = find_and_replace(s,"+dG",""); }
		else {
			m_iGamCat = atoi(s.c_str() + s.find("dG",0) - 1);
			if(s.c_str()[s.find("dG",0) - 2] != '+') { cout << "\nFormat for gamma distribution is incorrect in model name: " << s<< ". Should be of form  +?dG\n"; exit(-1); }
			s.erase(s.find("dG",0) - 2,4);
		}
		if(!InRange(m_iGamCat,2,10)) { Error("\nLeaphy only allows 2-10 discrete categories for the gamma distribution\n\n"); }
		m_bDoGamma = true;
	}
	// Do invariant sites
	if(Toks[0].find("+I",0) != string::npos) { m_bDoInv = true; s = find_and_replace(s,"+I",""); }
	// Do garbage collector
	if(Toks[0].find("+gC",0) != string::npos)	{ m_bDoGarbage = true; s = find_and_replace(s,"+gC",""); }
	// Do +F option
	if(Toks[0].find("+F",0) != string::npos) { m_bF = true; s = find_and_replace(s,"+F",""); } else { m_bF = false; }
	// Do Tree HMM if required
	if(s.find("+TreeHMM") != string::npos) {
		s = find_and_replace(s,"+TreeHMM",""); // cut name
		m_bDoTreeHMM = true;
	}

	// Do DNA models
	if(Toks[0].find(sModelNames[(int)COV_HKY],0) != string::npos)			{ m_Model = COV_HKY; s = find_and_replace(s,sModelNames[(int)COV_HKY],""); }
	else if(Toks[0].find(sModelNames[(int)COV_REV],0) != string::npos)	{ m_Model = COV_REV; s = find_and_replace(s,sModelNames[(int)COV_REV],""); }
	else if(Toks[0].find(sModelNames[(int)HKYdG_THMM],0) != string::npos) { m_Model = HKYdG_THMM; s = find_and_replace(s,sModelNames[(int)HKYdG_THMM],"");}
	else if(Toks[0].find(sModelNames[(int)RY_model],0) != string::npos)			{ m_Model = RY_model; s = find_and_replace(s,sModelNames[(int)RY_model],""); }
	else if(Toks[0].find(sModelNames[(int)JC],0) != string::npos)			{ m_Model = JC; s = find_and_replace(s,sModelNames[(int)JC],""); }
	else if(Toks[0].find(sModelNames[(int)FEL],0) != string::npos)		{ m_Model = FEL; s = find_and_replace(s,sModelNames[(int)FEL],""); }
	else if(Toks[0].find(sModelNames[(int)K2P],0) != string::npos)		{ m_Model = K2P; s = find_and_replace(s,sModelNames[(int)K2P],""); }
	else if(Toks[0].find(sModelNames[(int)HKY],0) != string::npos)		{ m_Model = HKY; s = find_and_replace(s,sModelNames[(int)HKY],""); }
	else if(Toks[0].find(sModelNames[(int)REV],0) != string::npos)		{ m_Model = REV; s = find_and_replace(s,sModelNames[(int)REV],""); }
	else if(Toks[0].find(sModelNames[(int)THMM_FULLDNA],0) != string::npos &&
		Toks[0].find(sModelNames[(int)THMM_AA],0) == string::npos)
	{	m_Model = THMM_FULLDNA; s = GetFullDNATHMMString(s); }
	else if(Toks[0].find(sModelNames[(int)THMM_DNA],0) != string::npos) { m_Model = THMM_DNA; s = GetDNATHMMString(s); }
	// Do protein models
	else if(Toks[0].find(sModelNames[(int)WAGdG_THMM],0) != string::npos)   { m_Model = WAGdG_THMM; s = find_and_replace(s,sModelNames[(int)WAGdG_THMM],""); }
	else if(Toks[0].find(sModelNames[(int)Coevo_WAG],0) != string::npos) { m_Model = Coevo_WAG; s = find_and_replace(s,sModelNames[(int)Coevo_WAG],""); }
	else if(Toks[0].find(sModelNames[(int)EQU],0) != string::npos)			{ m_Model = EQU; s = find_and_replace(s,sModelNames[(int)EQU],""); }
	else if(Toks[0].find(sModelNames[(int)WAG],0) != string::npos)		{ m_Model = WAG; s = find_and_replace(s,sModelNames[(int)WAG],""); }
	else if(Toks[0].find(sModelNames[(int)JTT],0) != string::npos)		{ m_Model = JTT; s = find_and_replace(s,sModelNames[(int)JTT],""); }
	else if(Toks[0].find(sModelNames[(int)DAY],0) != string::npos)		{ m_Model = DAY; s = find_and_replace(s,sModelNames[(int)COV_HKY],"DAY"); }
	else if(Toks[0].find(sModelNames[(int)mtREV],0) != string::npos)		{ m_Model = mtREV; s = find_and_replace(s,sModelNames[(int)mtREV],""); } // These REV models *MUST* be positioned after the DNA REV model
	else if(Toks[0].find(sModelNames[(int)THMM_AA],0) != string::npos)		{ m_Model = THMM_AA; s = GetAATHMMString(s); }
	// Do Codon models
	else if(Toks[0].find(sModelNames[(int)CodonM0],0) != string::npos)		{ m_Model = CodonM0; GetCodonString(s); s = find_and_replace(s,sModelNames[(int)CodonM0],""); }
	if(!s.empty()) { cout << "\nParsing of model name has resulted in left over characters: " << s << "\n\n"; exit(-1); }
	return true;
}
void CPhyloDat::ReadModel(string InputFile, CData *Data)	{
	int i, j,k,ParNum;
	string store,name;
	EDataType Type = NONE;
	vector <string> Toks;
	string Tree;
	vector <string> SeqNames;
	vector <string> ParNames;
	vector <double> ParVals;
	vector <CPar *> Pars;
	long StartingSeed = -1;
	// Initialise
	if(InputFile == "\0" || !FileExist(InputFile)) { InputFile = GetInFileName(); }
	CleanModel();
	if(m_pData != NULL && Data != NULL) { Error("\nTrying to pass data to a CPhyloData::ReadModel function when data already exists...\n\n"); }
	assert(m_pModel == NULL);
	// The file input commands
	FINOPEN(in,InputFile.c_str());
	getline(in,store,';');
	if(store.find("# Model file for Leaphy",0) == string::npos) { Error("\nTrying to input a model without correct header \"# Model file for Leaphy\"\n\n"); }
	while(!in.eof()) {
		getline(in,store,';');
		if(in.eof()) { break; }
		Toks = Tokenise(store);
		if(Toks.empty()) { continue; }
		if(Toks[0].find("#") != string::npos) { continue; }
		// 1. Get the model
		if(Toks[0].find("Model:") != string::npos) {
			if(Toks.size() < 2) { Error("\nLine with model must be of form \"Model:\t<ModelName>\"\n" + store + "\n"); }
			SetModel(Toks[1]);
			if(m_Model == UNKNOWN) { Error("\nInput model unknown: " + Toks[1] + "\n\n"); }
		}
		// 2. Get the EDataType
		else if(Toks[0].find("DataType:",0) != string::npos) {
			if(Toks.size() < 2) { Error("Line with DataType details must be of form \"DataType:\t<Type>\n"+store+"\n\n"); }
			m_PseudoDataType = Type = GetDataType(Toks[1]);
		}
		// 3. Get the data specifications
		else if(Toks[0].find("Data:") != string::npos) {
			if(Toks.size() < 3) { Error("Line with data details must be of form \"Data:\t<NoSeq>\t<Length>\"\n"+ store +"\n"); }
			m_iPseudoDataSeq = atoi(Toks[1].c_str()); m_iPseudoDataSize = atoi(Toks[2].c_str());
//			cout << "\nGot datasize: "<< m_iPseudoDataSize << "\nfrom: "<< Toks[2];
			if(m_iPseudoDataSeq < 0 || m_iPseudoDataSize < 0) { Error("Line with data details must be of form \"Data:\t<NoSeq>\t<Length>\"\n"+ store +"\n"); }
			if(m_pData != NULL) {
				if(NoSeq() != m_iPseudoDataSeq || m_pData->m_iTrueSize != m_iPseudoDataSize) {
					cout << "\nData in model file doesn't match that of current data file";
					cout << "\n\tDatafile:  " << NoSeq() << " " << m_pData->m_iTrueSize;
					cout << "\n\tModelfile: " << m_iPseudoDataSeq << " " << m_iPseudoDataSize;
					exit(-1);
		}	}	}
		// 4. Get parameter values
		else if(Toks[0].find("Parameter",0) != string::npos) {
			Toks = Tokenise(store,"\n"); store = Toks[0];
			Toks = Tokenise(store);
//			cout << "\nHave '" << store << "' -> ";
			if(Toks.size() < 4) { Error("Line with parameter details must be of form \"Paramteter[#]:\tName\t=\tValue\"\n"+ store +"\n"); }
			name = "";
			for(i=1;i<(int)store.size();i++) { if(store[i-1] == ']') { break; } }
			for(i;i<(int)store.size();i++) { if(isspace(store[i]) == 0) { break; } }
			rFOR(j,(int)store.size()-1){ if(store[j+1] == '=') { break; } }
			for(j;j>0;j--) { if(isspace(store[j-1]) == 0) { break; } }
			for(k=i;k<j;k++) { name += store[k]; }
//			cout << "'" << name << "'";
			ParNames.push_back(name);
			ParVals.push_back(atof(Toks[(int)Toks.size()-1].c_str()));
			}
		// 5. Get the tree
		else if(Toks[0].find("Tree:",0) != string::npos) {
			if(Toks.size() != 2) { Error("Line with tree details must be of form \"Tree:\t<tree_in_newick>\"; Note: There must be no spaces in the tree\n"+store+"\n\n"); }
			Tree = Toks[1];
		}
		// 6. Get the starting seed
		else if(Toks[0].find("Current_random_seed:",0) != string::npos) {
			if(Toks.size() < 2) { Error("Line with random seed must be of form \"Current_random_seed:\t<long>\""+store+"\n\n"); }
			StartingSeed = atol(Toks[1].c_str());
	}	}
	// Validate and create data
	cout << "\nCreating data of type: " << Type;
	if(m_pData == NULL && Data == NULL)  {
		if(!Tree.empty()) { SeqNames = Tree2Names(Tree,m_iPseudoDataSeq); }
		if(SeqNames.empty())	{ m_pData = new CData(m_iPseudoDataSeq,m_iPseudoDataSize,Type); }
		else					{ m_pData = new CData(m_iPseudoDataSeq,m_iPseudoDataSize,Type,&SeqNames); }
	} else if(m_pData == NULL) { m_pData = Data; }

//	cout << "\nData: Seq=" << m_pData->m_iNoSeq << "; size= " << m_pData->m_iSize;
//	cout << "SeqNames " << SeqNames;
	// Validate and create tree
	if(!Tree.empty()) { m_pTree = new CTree(Tree,m_iPseudoDataSeq,false,m_pData); }
	// Validate and create model
	CreateModel(); m_pModel->SetRandomSeed(StartingSeed);
	m_pModel->GetParVec(&Pars);
	if((int)ParNames.size() != (int)Pars.size()) { Error("\nNumber of parameters in specified model is " + int_to_string((int) Pars.size()) + "; Number of parameters in file is " + int_to_string((int)ParNames.size()) + "\n\n"); }
	FOR(i,(int)ParNames.size()) {
		ParNum = -1;
		// First guess should be good most of the time
		if(Pars[i]->Name() == ParNames[i]) { ParNum = i; }
		else { // Otherwise look for parameter
			FOR(j,(int)Pars.size()) { if(ParNames[i] == Pars[j]->Name()) { ParNum = i; break; } }
		}
		if(ParNum == -1) { cout << "\nCan't find parameter '" << ParNames[i] << " in list of Parameters: "; FOR(i,(int)Pars.size()) { cout << "\n" << i <<". " << *Pars[i]; } exit(-1); }
		Pars[ParNum]->SetVal(ParVals[i],false,true,false);
	}
	FOR(i,(int)Pars.size()) { Pars[i] = NULL; }
	Pars.clear();
	m_pModel->RedoScale(true);

//	cout << "\nNew Model: " << *m_pModel;
}

//////////////////////////////////////////////////////////////////////////////
// Data simulation routines
void CPhyloDat::DoSimulation(string ModelInputFile, string DataOutputFile, int NumberOfSims)	{
	int i,j,k,max_name =0,count=0;
	cout << "\nDoing " << NumberOfSims <<" simulations using <"<<ModelInputFile<<"> -> <"<<DataOutputFile<<">";
	ReadModel(ModelInputFile);
	SetRandomSeed((long)time(NULL));
	cout << "\nRandom seed: " << GetOriRandomSeed();

	FOR(i,m_pData->m_iNoSeq) { if((int) m_pData->m_vsName[i].size() > max_name) { max_name = (int)m_pData->m_vsName[i].size(); } }
	ofstream out(DataOutputFile.c_str());
	cout << "\n\nDoing Sims\n";
	FOR(i,NumberOfSims) {
		cout << "." << flush; if(count % 50 == 0 && count > 0) { cout << " " << count << endl; }
		SimulateData(true,true);
		out << m_pData->m_iNoSeq << " " << m_pData->m_iTrueSize << "\n";
		FOR(j,m_pData->m_iNoSeq) {
			out << m_pData->m_vsName[j];
			for(k=(int)m_pData->m_vsName[j].size();k<max_name;k++) { out << " "; }
			out << " " << m_pData->m_vsTrueSeq[j] << endl;
		}
		count ++;
	}
	cout << " " << count;
	out.close();
//	cout << "\nModel used for simulation output to sim.model\n";
//	ofstream modout("sim.model");
//	modout << *m_pModel << flush;
//	modout.close();
}

void CPhyloDat::SimulateData(bool Force,bool MakeData, CDataSummary *DS, vector <vector <bool> > *GapMask)	{
	int i,j;
	string ABET;
	vector <string> Names,Seqs;
	vector <vector <int> > vviSeqs;
	Names.clear();
	int Char = NumStates(m_PseudoDataType);
	// Check entry conditions
	if(m_pData != NULL && !Force) { Error("\nTrying to CPhyloDat::SimulateData from an object that already holds data...\n"); }
	if(m_pModel == NULL) { Error("\nTrying to CPhyloDat::SimulateData without a model specified\n"); }
	if(!m_pModel->IsViable()) { Error("\nTrying to CPhyloDat::SimulateData without a viable model\n"); }
	if(m_iPseudoDataSeq < 0 || m_iPseudoDataSize < 0) { Error("\nTrying to CPhyloDat::SimulateData when m_iPseudoDataSeq < 0 || m_iPseudoDataSize < 0\n\n"); }
	if(m_PseudoDataType == NONE) { Error("\nUnknown m_PseudoDataType in CPhyloDat::SimulateData\n\n"); }
	if(!Names.empty() && (int)Names.size() != m_iPseudoDataSeq) { Error("\nPhyloDat::SimulateData; Number of names in tree doesn't match what is expected...\n\n"); }
	if(DS == NULL && MakeData == false) { Error("\nDoing CPhyloDat::SimulateData when no valid form of output\n\n"); }
	if(GapMask != NULL) {
		if((int)GapMask->size() != m_iPseudoDataSeq) { Error("\nCPhyloDat::SimulateData GapMask doesn't have the right number of sequences\n\n"); }
		FOR(i,m_iPseudoDataSeq) { if((int)GapMask->at(i).size() != m_iPseudoDataSize) { Error("\nCPhyloDat::SimulateData GapMask doesn't have the right length of sequences\n\n"); } }
	}


	// Do data prep
	m_pModel->CleanFastCalc(true);
	CleanData();
	// Initialise
	Names = m_pTree->Names();
	if(Names.empty())	{ FOR(i,m_iPseudoDataSeq) { Names.push_back("Seq" + int_to_string(i+1)); } }
	// Create the sequences
	Seqs.clear(); vviSeqs.clear();
	m_pModel->DoSimulation(m_iPseudoDataSeq,m_iPseudoDataSize,&vviSeqs,false);
	// Do the gaps
	if(GapMask != NULL)	{ FOR(i,m_iPseudoDataSeq) { FOR(j,m_iPseudoDataSize) { if(GapMask->at(i)[j]) { vviSeqs[i][j] = Char; }  }	}	}
	// Make the sequence data
	if(MakeData)	{
		// Translate the integers to alphabet sequences
		if((int)vviSeqs.size() != (2*m_iPseudoDataSeq)-2) { Error("\nCreated the wrong number of sequences in CPhyloDat::SimulateData"); }
		FOR(i,m_iPseudoDataSeq)	{
			if((int)vviSeqs[i].size() != m_iPseudoDataSize) { Error("\nCreate sequence["+int_to_string(i)+"] of wrong length in CPhyloData::SimulateData\n\n"); }
			Seqs.push_back("");
			FOR(j,m_iPseudoDataSize) {
				if(vviSeqs[i][j] == Char) { Seqs[i] += "-"; }
				Seqs[i] += State(m_PseudoDataType,vviSeqs[i][j]); }
		}
		// Input the data
		if(m_pData != NULL) { Error("\nm_pData should be null...\n"); }
		m_pData = new CData(m_iPseudoDataSeq,m_iPseudoDataSize,Seqs,Names,m_PseudoDataType);
		m_pModel->m_pData = m_pData;
	}
	// Create the data summary
	if(DS != NULL) { DS->MakeDataSummary("Sim["+int_to_string(++SimNum)+"]",m_PseudoDataType,&vviSeqs,&Names); }
}


/*  	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
* 			Code for calculating partial likelihoods and outputting them from a set model file
*/ 	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CPhyloDat::DoPartL(CData *Dat, string model_file, string out_file) {
	cout << "\nCreating partial likelihoods from model file <" << model_file << "> -> <"<<out_file<<">";
	ReadModel(model_file, Dat);

/*	cout << "\nThe lnL is: " << m_pModel->lnL(true);
	cout << "\nPerforming BranchOpt: ";
	m_pModel->GetOptPar();
	m_pModel->FastBranchOpt(m_pModel->lnL(true));
	cout << "\nDoing rest..."; */

	// Do partial likelihood stuff
	pModel()->PartialLOutput(out_file);
}



////////////////////////////////////////////////////////////////
// Do the starting tree
void CPhyloDat::SetStartTree(string s)	{
	int i;
	string file;
	string::size_type pos;
	CTree *T;
//	cout << "\nWorking with string: " << s;
	m_bDoSA = false; m_bDoBioNJ = false;
	// Do automatic stuff
	if(s.find("SA") != string::npos && s.size() == 2)	{ m_bDoSA = true; }
	else if(s.find("bionj") != string::npos && s.size() == 5)	{ m_bDoBioNJ = true; }
	else if(s.find("all") != string::npos && s.size() == 3)	{ m_bDoSA = m_bDoBioNJ = true; }
	// Do the file
	else if(s.find("file=") != string::npos && (m_pData != NULL || !m_vpMultiDat.empty()))	{
		pos = s.find("file=") + 5;
		file = s.substr(pos);
		T = new CTree(file,true,m_pData,true);
		if(T->Valid()) { m_sTreeFile = file; m_pTree = T; } else { m_bDoSA = m_bDoBioNJ = true; delete T; }
		// Do the multiple data set stuff if required
		if(!m_vpMultiDat.empty()) {
			if(!T->Valid()) { Error("\nTrying to initialise from an invalid tree...\n\n"); }
			cout << "\nInitialising multitree to\n" << *m_pTree << flush;
			FOR(i,(int)m_vpMultiTree.size()) {
				m_vpMultiTree[i] = new CTree();
				*m_vpMultiTree[i] = *m_pTree;
				if(i%10 == 0) { cout << endl; } cout << m_vsMultiName[i] << "\t";
			}
		}
		T = NULL;
	} else { Error("\nUnknown choice of tree..."); }
}

////////////////////////////////////////////////////////////////
// Create the model
// ----------------
// Bloody messy function, but what the hell
void CPhyloDat::CreateModel()	{
	double pInv = 0.0, sigInv = 0.0; bool VSig = false, VInv = false;
	// DNA
	CJC *mJC; CFEL *mFEL; CK2P *mK2P; CHKY *mHKY; CREV *mREV;
	// Exotic nucleotide based models
	CRY *mRY; CCovHKY *mCovHKY; CCovREV *mCovREV; CTHMM_FULL *mTHMM_FULLDNA; CHKYdG_THMM *mHKYdGTHMM;
	// Protein
	CEQU *mEQU; CWAG *mWAG; CJTT *mJTT; CDAY *mDAY; CMTREV *mMTREV;
	// Exotic amino acid models
	CTHMMAA *mTHMMAA; CTHMMBAA *mTHMMBAA; CTHMMDNA *mTHMMDNA; CTHMMBDNA *mTHMMBDNA; CWAGdG_THMM *mCWAGdGTHMM;
	CWAGCoevo *mCoevoWAG;
	// Codon
	CCodonM0 *mCodonM0;
	// For multiple data sets get the largest number of sequences and length of alignment
	if(!m_vpMultiDat.empty()) {
		int i,Sp = 0, Size = 0;
		assert(m_pData == NULL);
		FOR(i,(int)m_vpMultiDat.size()) {
			if(m_vpMultiDat[i]->m_iNoSeq > Sp) { Sp = m_vpMultiDat[i]->m_iNoSeq; }
			if(m_vpMultiDat[i]->m_iSize > Size) { Size = m_vpMultiDat[i]->m_iSize; }
		}
		cout << "\nSp: " << Sp << "; Size: " << Size;
		m_pData = new CData(Sp,Size,m_vpMultiDat[0]->m_DataType);

	}
	switch(m_Model)	{
	////////////////////////// DNA models ////////////////////////////////
	case JC:
		mJC = new CJC(m_pData,m_pTree); m_pModel = mJC; mJC = NULL;
		break;
	case FEL:
		mFEL = new CFEL(m_pData,m_pTree); m_pModel = mFEL; mFEL = NULL;
		break;
	case K2P:
		mK2P = new CK2P(m_pData,m_pTree); m_pModel = mK2P; mK2P = NULL;
		break;
	case HKY:
		mHKY = new CHKY(m_pData,m_pTree); m_pModel = mHKY; mHKY = NULL;
		break;
	case REV:
		mREV = new CREV(m_pData,m_pTree); m_pModel = mREV; mREV = NULL;
		break;
	case RY_model:
		mRY = new CRY(m_pData,m_pTree); m_pModel = mRY; mRY = NULL;
		break;
	case COV_HKY:
		mCovHKY = new CCovHKY(m_pData,m_pTree); m_pModel = mCovHKY; mCovHKY = NULL;
		break;
	case COV_REV:
		mCovREV = new CCovREV(m_pData,m_pTree); m_pModel = mCovREV; mCovREV = NULL;
		break;
	case THMM_FULLDNA:
		if(THMM_F == true)	{ m_bOptFre = true; mTHMM_FULLDNA = new CTHMM_FULL(m_pData,m_pTree,THMM_NoCat,THMM_R,NULL,complex,THMM_K,THMM_H); }
		else				{ mTHMM_FULLDNA = new CTHMM_FULL(m_pData,m_pTree,THMM_NoCat,THMM_R,NULL,obs,THMM_K,THMM_H); }
		m_pModel = mTHMM_FULLDNA; mTHMM_FULLDNA = NULL;
		break;
	case THMM_DNA:
		if(THMM_R == varyall) { sigInv = 0.15; }
		if(THMM_K == true) { pInv = 0.1; }
		if(!THMM_F) {
			mTHMMDNA = new CTHMMDNA(m_pData, m_pTree,0.5,4,pInv,0.1,sigInv); m_pModel = mTHMMDNA; mTHMMDNA = NULL;
		} else {
			switch(THMM_H) {	// Set up the correct things to vary
			case H_diff: VInv = VSig = true; break;
			case H_same: VSig = true; break;
			case H_none: VInv = true; break;
			}
			mTHMMBDNA = new CTHMMBDNA(m_pData, m_pTree,0.5,4,pInv,0.1,sigInv,VSig,VInv); m_pModel = mTHMMBDNA; mTHMMBDNA = NULL;
		}
		break;
	case HKYdG_THMM:
		cout << "\nDoing CHKYdG_THMM...";
		mHKYdGTHMM = new CHKYdG_THMM(m_pData,m_pTree,4); m_pModel = mHKYdGTHMM; mHKYdGTHMM = NULL;
		break;
	////////////////////////// Protein models /////////////////////////////
	case EQU:
		mEQU = new CEQU(m_pData,m_pTree,DoF()); m_pModel = mEQU; mEQU = NULL;
		break;
	case WAG:
		mWAG = new CWAG(m_pData,m_pTree,DoF()); m_pModel = mWAG; mWAG= NULL;
		break;
	case JTT:
		mJTT = new CJTT(m_pData,m_pTree,DoF()); m_pModel = mJTT; mJTT = NULL;
		break;
	case DAY:
		mDAY = new CDAY(m_pData,m_pTree,DoF()); m_pModel = mDAY; mDAY = NULL;
		break;
	case mtREV:
		mMTREV = new CMTREV(m_pData,m_pTree,DoF()); m_pModel = mMTREV; mMTREV = NULL;
		break;
	case THMM_AA:
		if(THMM_R == varyall) { sigInv = 0.15; }
		if(THMM_K == true) { pInv = 0.1; }
		if(!THMM_F) {
			mTHMMAA = new CTHMMAA(m_pData, m_pTree, 0,0.5,4,pInv,0.1,sigInv); m_pModel = mTHMMAA; mTHMMAA = NULL;
		} else {
			switch(THMM_H) {	// Set up the correct things to vary
			case H_diff: VInv = VSig = true; break;
			case H_same: VSig = true; break;
			case H_none: VInv = true; break;
			}
			mTHMMBAA = new CTHMMBAA(m_pData, m_pTree, 0,0.5,4,pInv,0.1,sigInv,VSig,VInv); m_pModel = mTHMMBAA; mTHMMBAA = NULL;
		}
		break;
	case WAGdG_THMM:
		mCWAGdGTHMM = new CWAGdG_THMM(m_pData,m_pTree,0);
		m_pModel = mCWAGdGTHMM; mCWAGdGTHMM = NULL;
		break;
//	case rtREV:
//		cout << "\nModel not currently ready..."; exit(-1);
//		break;
	///////////////////////// Codon models /////////////////////////////////
	case CodonM0:
		mCodonM0 = new CCodonM0(m_pData,m_pTree,m_CodonEqm,m_iGenCode); m_pModel = mCodonM0; mCodonM0 = NULL;
		break;
	///////////////////////// Coevolution models ///////////////////////////
	case Coevo_WAG:
		cout << "\nCoevo_WAG creation:";

//		tempvec.assign(m_pData->m_iNoSeq,0); int i;
//		FOR(i,NoSeq()) { tempvec[i] = m_pData->m_ariSeq[i][3]; } tempvec[3] = 20;
//		m_pData->AddColumn(tempvec,2);

		cout << "\nDone...";

		mCoevoWAG = new CWAGCoevo(m_pData->CreateAllPairs(m_pData,true),m_pData, m_pData, m_pTree);
		m_pModel = mCoevoWAG; mCoevoWAG = NULL;
		break;
	default:
		Error("Unknown model in CPhyloDat::CreateModel()...\n");
	};
	// Do rates-across-sites
	if(m_bDoGamma)	{ m_pModel->MakeGammaModel(0,m_iGamCat); }
	if(m_bDoInv)	{ m_pModel->MakeInvariantSitesModel(0); }
	if(m_bDoGarbage){ m_pModel->MakeGarbageCollectorModel(0); }
	// Do Tree HMM stuff
	if(m_bDoTreeHMM) { m_pModel->MakeTreeHMM(); }

	if(!m_vpMultiDat.empty()) { m_pModel->CreateProcessSpace(); CleanData(); }
}

/////////////////////////////////////////////////////////
// Functions to get the THMM stuff
void CPhyloDat::GetTHMMStdIn()	{
	int Option;
	THMM_NoCat = GetOption("\nChoose number of hidden states",2,5);
	if(GetYesNoOption("\nAllow variation in rates across states")) { THMM_R = varyall; } else { THMM_R = same; };
	THMM_F = GetYesNoOption("\nAllow variation in nucleotide composition across states");
	THMM_K = GetYesNoOption("\nAllow variation in Ts/Tv rate ratio across states");
	Option = GetOption("\nChoose how transitions between hidden states are defined; 1: None (mixture model); 2: One set; 3: All transitions independent",1,3);
	switch(Option) {
		case 1: THMM_H = H_none; break;
		case 2: THMM_H = H_same; break;
		case 3: THMM_H = H_diff; break;
		default: Error("Unknown THMM_H type in CPhyloDat::GetTHMMStdIn");
	};
	if(!VerifyTHMM()) { Error("THMM specified is not valid..."); }
}

// The THMM models have a series of flags that are used to define the model
//	The form of the flags varies with model
//	For the THMM_FULLDNA: THMM.<#Cat>.<R><F><K><-H/H>
// The flag values represent:
//	THMM_H:	H_none = mixture model;	H_diff = different rates between hidden classes; H_same = 1 rate between classes
//	THMM_R: true = different rates between classes; false = one rate
//  THMM_K: true = different kappas between classes; false = one kappa
//	THMM_F: true = different freqs between classes; false = one set of freq
string CPhyloDat::GetFullDNATHMMString(string Str)	{
	string ErrorMessage = "\nTemporal HMM should be specified in form THMM_FULLDNA.<NoCat>.<R><F><K><-H/H>\n";
	vector <string > Toks;
	// Clean up the string
	if(Str.find("THMM") == string::npos) { Error("Passed string " + Str + " to CPhyloDat::GetDNATHMMString, but doesn't contain \"THMM\"...\n"); }
	transform(Str.begin(),Str.end(),Str.begin(),(int(*)(int)) toupper);
	Str = Str.replace(0,Str.find("THMM_FULLDNA")+12,"");
	// Working with the dots that specify how many categories

	if(Str[0] != '.' || Str[2] != '.') { Error(ErrorMessage + "The dots seems to be in the wrong place for #cats...\n\n"); }
	THMM_NoCat = atoi(Str.c_str()+ 1);
	if(!InRange(THMM_NoCat,2,5)) { Error("\nTemporal HMM can only have between 2 and 5 hidden states; you specified: " + int_to_string(THMM_NoCat) + "\n\n"); }
	Toks = Tokenise(Str,"."); if(Toks.size() == 0) { Str = ""; } else if(Toks.size() > 2) { Error("\nToo much information passed to CPhyloDat::GetFullDNATHMMString " + Str +"\n\n"); } else { Str = Toks[1]; }

	// The next bit specifies what type of parameters are in the model
	if(Str.find("-H") != string::npos)		{ THMM_H = H_none; }
	else if(Str.find("H") != string::npos)	{ THMM_H = H_diff; }
	else									{ THMM_H = H_same; }
	if(Str.find("R") != string::npos)		{ THMM_R = varyall; }
	else									{ THMM_R = same; }
	if(Str.find("K") != string::npos)		{ THMM_K = true; }
	else									{ THMM_K = false; }
	if(Str.find("F") != string::npos)		{ THMM_F = true; }
	else									{ THMM_F = false; }
	if(Str.find("ALL") != string::npos)		{ THMM_H = H_diff; THMM_R = varyall; THMM_F = true; THMM_K = true; }
	if(!VerifyTHMM()) { Error("\nTHMM specified by cmd line is not valid\n" + ErrorMessage); }
	return "";
}
// This is the controller class for the amino acid THMMs (both THMM_AA and THMM_BAA, which varies things across branches)
// The flag values represent:
//	THMM_H:	H_none = Vary pInv by branch; H_diff = vary pInv and sigma by branch; H_same = differ only sigma by branch
//	THMM_R: varyall = SigmaInv and SigmaAlfa are different; same = single Sigma
//  THMM_K: true = Allows a pInv parameter; false = no invariant site category
//	THMM_F: true = things vary over branches; false = things don't vary over branches
// --
// The format for THMM_(B)AA is: THMM_<B>AA.<S><I>.<S><I>
//	The first S represents whether there is one sigma or two
//	The first I represents whether there is an invariant sites category
// The second set only applies for THMM_BAA and describes what varies across branches
// --
// Function returns the string minus the THMM_<B>AA.<S><I>.<S><I> bit
string CPhyloDat::GetAATHMMString(string Str)	{
	vector <string> Toks;
	bool F1 = false,F2 = false;
	transform(Str.begin(),Str.end(),Str.begin(),(int(*)(int)) toupper);
	Toks = Tokenise(Str,".");
	string ErrorMessage = "\nAmino acid temporal HMM should be specified in form THMM_AA.<I><S>.<I><S>\nFirst:  <I> = invariant sites present;  <S> = seperate sigma for gamma and invariant sites; otherwise filler (e.g. none)\nSecond: <I> = P(inv) varies by branch;  <S> Sigma(s) vary by branch\n";
	// Error check
	THMM_F = false;
	if(Toks[0].find("THMM") == string::npos) { Error("Passed string " + Str + " to CPhyloDat::GetAATHMMString, but doesn't contain \"THMM\"...\n"); }
	if((int)Toks.size() > 3) { Error("Too much information passed to CPhyloDat::GetAATHMMString...\n" + ErrorMessage); }
	// Get the info in the secondtoken
	if(Toks.size() > 1) {
		if(Toks[1].find("S") != string::npos)		{ THMM_R = varyall; } else { THMM_R = same; }
		if(Toks[1].find("I") != string::npos)		{ THMM_K = true; } else { THMM_K = false; }
	} else {
		THMM_R = same; THMM_K = false;
	}
	// Get the info in the third token
	if(Toks.size() > 2) {
		if(Toks[2].find("S") != string::npos)		{ F1 = true; }
		if(Toks[2].find("I") != string::npos)		{ F2 = true; }
		THMM_F = true;
		if(F1 && F2)	{ THMM_H = H_diff; }
		else if(F1)		{ THMM_H = H_same; }
		else if(F2)		{ THMM_H = H_none; }
		else { THMM_H = H_none; THMM_F = false; }
	} else { THMM_H = H_none; THMM_F = false; }
	string str_temp = Str;
	str_temp.erase(0,Str.find("THMM_"));
	if(str_temp.find(" ") != string::npos) { str_temp = str_temp.erase(str_temp.find(" "),str_temp.size() - 1); }
	return find_and_replace(Str,str_temp,"");
}

// This is the controller class for the nucleotide THMMs (both THMM_DNA and THMM_BDNA, which varies things across branches)
// The flag values represent:
//	THMM_H:	H_none = Vary pInv by branch; H_diff = vary pInv and sigma by branch; H_same = differ only sigma by branch
//	THMM_R: varyall = SigmaInv and SigmaAlfa are different; same = single Sigma
//  THMM_K: true = Allows a pInv parameter; false = no invariant site category
//	THMM_F: true = things vary over branches; false = things don't vary over branches
// --
// The format for THMM_(B)DNA is: THMM_<B>DNA.<S><I>.<S><I>
//	The first S represents whether there is one sigma or two
//	The first I represents whether there is an invariant sites category
// The second set only applies for THMM_BAA and describes what varies across branches
// --
// Function returns the string minus the THMM_<B>AA.<S><I>.<S><I> bit
string CPhyloDat::GetDNATHMMString(string Str)	{
	vector <string> Toks;
	bool F1 = false,F2 = false;
	transform(Str.begin(),Str.end(),Str.begin(),(int(*)(int)) toupper);
	Toks = Tokenise(Str,".");
	string ErrorMessage = "\nNucleotide temporal HMM should be specified in form THMM_DNA.<I><S>.<I><S>\nFirst:  <I> = invariant sites present;  <S> = seperate sigma for gamma and invariant sites; otherwise filler (e.g. none)\nSecond: <I> = P(inv) varies by branch;  <S> Sigma(s) vary by branch\n";
	// Error check
	THMM_F = false;
	if(Toks[0].find("THMM") == string::npos) { Error("Passed string " + Str + " to CPhyloDat::GetDNATHMMString, but doesn't contain \"THMM\"...\n"); }
	if((int)Toks.size() > 3) { Error("Too much information passed to CPhyloDat::GetDNATHMMString...\n" + ErrorMessage); }
	// Get the info in the secondtoken
	if(Toks.size() > 1) {
		if(Toks[1].find("S") != string::npos)		{ THMM_R = varyall; } else { THMM_R = same; }
		if(Toks[1].find("I") != string::npos)		{ THMM_K = true; } else { THMM_K = false; }
	} else {
		THMM_R = same; THMM_K = false;
	}
	// Get the info in the third token
	if(Toks.size() > 2) {
		if(Toks[2].find("S") != string::npos)		{ F1 = true; }
		if(Toks[2].find("I") != string::npos)		{ F2 = true; }
		THMM_F = true;
		if(F1 && F2)	{ THMM_H = H_diff; }
		else if(F1)		{ THMM_H = H_same; }
		else if(F2)		{ THMM_H = H_none; }
		else { THMM_H = H_none; THMM_F = false; }
	} else { THMM_H = H_none; THMM_F = false; }
	string str_temp = Str;
	str_temp.erase(0,Str.find("THMM_"));
	if(str_temp.find(" ") != string::npos) { str_temp = str_temp.erase(str_temp.find(" "),str_temp.size() - 1); }
	return find_and_replace(Str,str_temp,"");
}

//////////////////////////////////////////////////////////////
// Functions to get the codon stuff
// Models of form CodonM0.<GenCode>.<Eqm>
void CPhyloDat::GetCodonString(string s) {
	bool Throw = false;
	vector <string> Toks;
	Toks = Tokenise(s,".");;
	if((int)Toks.size() != 3) { Error("\nCodon model must be of form <CodonModel>.<GenCode>.<EqmType>\n"); }
	// Get the genetic code
	m_iGenCode = atoi(Toks[1].c_str());
	if(!InRange(m_iGenCode,0,11)) { Throw = true; }
	// Get the equilibrium
	if(Toks[2].find("EQU") != string::npos)			{ m_CodonEqm = cEQU; }
	else if(Toks[2].find("F1X4") != string::npos)	{ m_CodonEqm = F1X4; }
	else if(Toks[2].find("F3X4") != string::npos)	{ m_CodonEqm = F3X4; }
	else if(Toks[2].find("F64") != string::npos)	{ m_CodonEqm = F64; }
	else { Throw = true; }
	// Do error
	if(Throw) { Error("\nCodon model must be of form <CodonModel>.<GenCode>.<EqmType>\nYou input: " + s + "\n"); }
}

// Checks whether the THMM is a sensible tree or not
bool CPhyloDat::VerifyTHMM()	{
	// For case when there shouldn't be a THMM
	if(THMM_NoCat == 1)	{
		if(THMM_H != H_none || THMM_R != false || THMM_K != false || THMM_F != false) { return false; }
		return true;
	}
	if(!InRange(THMM_NoCat,2,5)) { return false; }
	if(THMM_R == false && THMM_K == false && THMM_F == false) { return false; }
	return true;
}
string CPhyloDat::THMM_Details()	{
	string RetStr = "." + int_to_string(THMM_NoCat) + ".";
	if(THMM_R == varyall) { RetStr += "R"; }
	if(THMM_F == true) { RetStr += "F"; }
	if(THMM_K == true) { RetStr += "K"; }
	if(THMM_H == H_none) { RetStr += "-H"; } else if(THMM_H == H_diff) { RetStr += "H"; }
	return RetStr;
}

////////////////////////////////////////////////////////////////
// Some functions relating to multiple data sets
int CPhyloDat::NoDataSets() {
	if(m_vpMultiDat.empty()) { if(m_pData != NULL) { return 1; } else { return 0; }}
	return (int) m_vpMultiDat.size();
}
void CPhyloDat::SetMultiSet(int Set)	{
	if(!InRange(Set,0,(int)m_vpMultiDat.size())) { Error("\nTrying to set data to strange data set in PhyloDat::SetMultiSet..."); }
	assert(m_vpMultiTree.size() == m_vpMultiDat.size());
//	cout <<"\nDataSet["<<Set<<"]: " << m_vpMultiDat[Set]->m_iNoSeq << ":" << m_vpMultiDat[Set]->m_iSize << flush;
	// Do model stuff
	if(Set != -1) { m_pModel->CleanFastCalc(true); }
	m_pModel->MakeNewTree(m_vpMultiTree[Set],false);
	m_pModel->MakeNewData(m_vpMultiDat[Set],false);
	m_pModel->CleanFastCalc(true);
	// Do internal pointer
	m_iCurrentMulti = Set;
}

int CPhyloDat::GetMultiRFDist(int Gene1,int Gene2) {
	assert(InRange(Gene1,0,(int)m_vpMultiTree.size()));
	assert(InRange(Gene2,0,(int)m_vpMultiTree.size()));
	return m_vpMultiTree[Gene1]->GetRFDist(*m_vpMultiTree[Gene2]);
}

///////////////////// ostream operator /////////////////////////
ostream &operator<<(ostream &os, CPhyloDat Dat)	{
	cout << "\nTemporary CPhyloDat output:";
	cout << "\n\tInputfile:  " << Dat.In();
	cout << "\n\tOutputfile: " << Dat.Out();
	cout << "\n\tModel: " << Dat.m_Model;
	if(Dat.m_bDoGamma) { cout << "+"<<Dat.m_iGamCat << "dG"; }
	if(Dat.m_bDoInv) { cout << "+I"; }
	cout << "\nStart tree: ";
	if(Dat.IsTree()) { cout << *Dat.m_pTree; }
	else { if(Dat.m_bDoBioNJ) { cout << " bionj"; } if(Dat.m_bDoSA) { cout << " SA"; } }
	return os;
}

////////////////////////////////////////////////////////////////////////////
// Function that will characterise tree-space
void CPhyloDat::LikelihoodSurfaces(ostream & os, vector <CPar *> *Pars) {
	int i,j, count;
	bool CleanSpace = false;
	double dVal;
	cout << "\nEntering lnL Surface..." << flush;
	vector <vector <double> > vvlnLs;
	vector <vector <double> > vvVals;
	vector <double> vlnLs;
	vector <double> values, oriVals;
	os.setf(ios::fixed); os.precision(6);
	if(m_pModel == NULL || m_pTree == NULL) { cout << "\nCannot calculate Likelihood Surfaces for a model and tree that are not initialised\n\n"; }
	if(Pars == NULL) {
		Pars = new vector <CPar *>;
		FOR(i,m_pModel->NoPar()) {
			if(m_pModel->m_vpPar[i]->Name().find("P(") != string::npos || m_pModel->m_vpPar[i]->Special() == true || m_pModel->m_vpPar[i]->CheckLocked() || m_pModel->m_vpPar[i]->NeedUpdatePar() == true) { continue; }
			Pars->push_back(m_pModel->m_vpPar[i]);
	}	}

	FOR(i,(int)Pars->size()) { vvlnLs.push_back(values); vvVals.push_back(values); oriVals.push_back(Pars->at(i)->Val()); }
	// Do the likelihood surfaces
	cout << "\nOriginal likelihood: " << m_pModel->lnL(true);

	FOR(i,(int)Pars->size()) {
		cout << "\nVal Par["<<i<<"]: " << *Pars->at(i) << flush;
		for(dVal = 0; dVal < max(1.0,oriVals[i]) * 1.5; dVal += max(1.0,oriVals[i]) / 200) {
			FOR(j,(int)oriVals.size()) { Pars->at(j)->SetVal(oriVals[j],true); }
			Pars->at(i)->SetVal(dVal,true);
			vvVals[i].push_back(dVal);
			vvlnLs[i].push_back(m_pModel->lnL());
		}
	}
	FOR(j,(int)oriVals.size()) { Pars->at(j)->SetVal(oriVals[j],true); }
	cout << "\nAfterwards lnL: " << m_pModel->lnL() << "\n\n---\n";

	count = 0; FOR(i,(int)Pars->size()) { if(vvVals[i].size() > count) { count = vvVals[i].size(); } }

	FOR(i,(int)Pars->size()) { os << Pars->at(i)->Name() << ".val\t" << Pars->at(i)->Name() << ".lnL\t"; }
	FOR(i,count) {
		cout << endl;
		FOR(j,(int)Pars->size()) {
			if(i>(int) vvVals[j].size()) { os << "NA\tNA\t"; continue; }
			os << vvVals[j][i] << "\t" << vvlnLs[j][i] << "\t";
		}
	}

	if(CleanSpace) { delete Pars; Pars = NULL; }
	cout << "\nDone..." << flush;
}


////////////////////////////////////////////////////////////////
// Tabu management routines
bool IsTabu(CTree *Tree, bool OnlyRadTabuOld, bool AddToTabu,CBaseModel *M, double lnL, bool IsOptima)	{
#if DO_TABU
	int i,Val;
	STabuTree Test;
	vector <int> List;
	// Check some entry conditions
	if(Tree->IsCutTree() == true || M->IsRMSDCalc()) { return false; }
	List = Tree->ConstOut();
	FOR(i,(int)TabuTrees.size())	{
		// If only tabuing radius on old then do this bit
		if(OnlyRadTabuOld == false) {
			Val = MinTabuRFDist(Tree);
			if(Val < 0) { return false; }
			else if(Val < TABU_RADIUS) { return true; }
		} else	{ // Otherwise do the dumb thing
			assert(TabuTrees[i].List.size() == List.size());
			if(TabuTrees[i].List == List) { return true; }
	}	}
	if(AddToTabu == true) {
		Test.Tree = *Tree; Test.Tree.SetOldTree(false);
		Test.lnL = lnL;
		Test.List = List;
		Test.Optima = IsOptima;
		FOR(i,(int)M->m_vpPar.size()) { Test.vPar.push_back(M->m_vpPar[i]->Val()); }
		TabuTrees.push_back(Test);
	}
#endif
	return false;
}

int MinTabuRFDist(CTree * Tree,bool OnlyCompare2Old)	{
	if(Tree->IsCutTree() == true) { return -1; }
	int i,j,k,MinDist = BIG_NUMBER, ThisDist;
	static vector <int> *OL = NULL, *L = NULL, *R = NULL;
	if(OL == NULL)	{
		FOR(i,Tree->NoBra()) {
			OL = new vector<int>[Tree->NoBra()];
			L  = new vector<int>[Tree->NoBra()]; R  = new vector<int>[Tree->NoBra()];
	}	}
	if(TabuTrees.empty()) { return -1; }
	// Get Tree's branch sets
	FOR(i,Tree->NoBra())	{ Tree->BranchSets(i,&OL[i],&R[i]); }
	FOR(i,(int)TabuTrees.size())	{
		assert(TabuTrees[i].Tree.NoSeq() == Tree->NoSeq());
		if(!TabuTrees[i].Tree.OldTree() && OnlyCompare2Old) { continue; }
		ThisDist = Tree->NoBra();
		// Loop through trees and compare
		FOR(j,TabuTrees[i].Tree.NoBra())	{
			// Get Tabu trees branch sets
			L[j].clear(); R[j].clear();
			TabuTrees[i].Tree.BranchSets(j,&L[j],&R[j]);
			// compare the Left one to the new trees left one
			FOR(k,Tree->NoBra())	{ if(Compare(&OL[k],&L[j]) == true) { ThisDist--; break; } }
		}
		if(ThisDist < MinDist) { MinDist = ThisDist; }
	}
	return MinDist;
}

void AgeTabuTrees() {
	int i;
	FOR(i,(int)TabuTrees.size()) { TabuTrees[i].Tree.SetOldTree(true); }
}

void VivifyTabuTrees() {
	int i;
	FOR(i,(int)TabuTrees.size()) { TabuTrees[i].Tree.SetOldTree(false); }
}

void RemoveTabuTree(CTree *Tree) {
	vector <STabuTree>::iterator iT = TabuTrees.begin();
	vector <int> List = Tree->ConstOut();
	IFOR(iT,TabuTrees) {
		if(iT->List == List) { TabuTrees.erase(iT); iT--; }
	}
}

bool OldTabuTrees()	{
	int i;
	FOR(i,(int)TabuTrees.size()) { if(TabuTrees[i].Tree.OldTree()) { return true; } }
	return false;
}

/////////////////////////////////////////////////////////////////////
// Class defining some summary statistics of a sequence alignment
// Constructor NULL
CDataSummary::CDataSummary() { m_sName = "EMPTY"; m_iNoDataPatterns = 0; m_bDoFullOutput = false; m_iNoSeq = -1; m_iSize = -1; }
// Constructor 1.
CDataSummary::CDataSummary(string Name, EDataType Type, std::vector<vector<int> > *Seqs, std::vector<string> *Names)	{ MakeDataSummary(Name,Type,Seqs,Names); }
void CDataSummary::MakeDataSummary(string Name, EDataType Type, std::vector<vector<int> > *Seqs, std::vector<string> *Names)	{
	int i,j,TotalSites = 0,Char = NumStates(Type);
	vector <int> States, Pat;
	vector <vector <int> > PatObs;
	vector <vector <double> > SeqFrq;
	double Sum;
	m_sName = Name; m_DataType = Type;
	m_iNoDataPatterns = 0;
	m_bDoFullOutput = false;
	m_dlnL = m_dlnL_multi = m_ddelta = BIG_NUMBER;
	// Check some entry conditions
	if(Seqs->size() != Names->size() && (Names->size() * 2) -2 != Seqs->size() ) { Error("\nThe number of names and the number of sequences should be equal in CDataSummary::CDataSummary...\n\n"); }
	FOR(i,(int)Seqs->size()-1) { if(Seqs->at(i).size() != Seqs->at(i+1).size()) { Error("\nThe length of all sequences should be equal in CDataSummary::CDataSummary...\n\n"); } }
	// Get the names and other basic data details
	m_vsNames = *Names;
	m_iNoSeq = (int)Names->size(); m_iSize = (int)Seqs->at(0).size();
	// Get the distribution of number of characters in sites
	m_vdStatesObs.assign(min(m_iNoSeq,NumStates(Type)),0);
	SeqFrq.assign(Char,vector<double>()); FOR(i,Char) { SeqFrq[i].assign(m_iNoSeq,0); }
	m_vdFrq.assign(Char,0.0);
	m_vdFrqStdErr.assign(NumStates(m_DataType),0.0);
	// Get the info
	FOR(i,m_iSize)	{
		States.clear();
		Pat.clear();
		FOR(j,m_iNoSeq) {
			if(Seqs->at(j)[i] >= Char) { continue; }
			// Do frqs
			SeqFrq[Seqs->at(j)[i]][j] += 1.0; TotalSites++;
			// Do patterns
			Pat.push_back(Seqs->at(j)[i]);
			if(!IsIn(Seqs->at(j)[i],States)) { States.push_back(Seqs->at(j)[i]); }
		}
		// Get the number of different observed states
		m_vdStatesObs[(int)States.size() - 1] ++;
		// Get the site pattern
		if(!(IsIn(Pat,PatObs))) { PatObs.push_back(Pat); m_iNoDataPatterns++; }
	}


	// Get the frequency information
	FOR(i,m_iNoSeq) { Sum = 0.0; FOR(j,Char) { Sum += SeqFrq[j][i]; }	FOR(j,Char) { SeqFrq[j][i] /= Sum; }	}
//	cout.precision(12); cout << "\nThe sequence frequencies [total =" << TotalSites << "]: ";
//	FOR(i,4) {  cout << "\nstate["<<i<<"]: " << SeqFrq[i]; }
	FOR(i,Char) {
		Sum = 0;
		FOR(j,m_iNoSeq) {
			Sum += (double) SeqFrq[i][j];
			m_vdFrqStdErr[i] += (double) SeqFrq[i][j] * (double) SeqFrq[i][j];
		}
		m_vdFrq[i] = (double) Sum / (double)m_iNoSeq;
		m_vdFrqStdErr[i] = sqrt((m_vdFrqStdErr[i] - (Sum * m_vdFrq[i])) / (m_iNoSeq - 1));
	}
//	cout << "\nValues produced: ";
//	FOR(i,4) { cout << "\nState["<<i<<"]: mean = " << m_vdFrq[i] << "\tStdErr: " << m_vdFrqStdErr[i]; }
}
// Constructor 2.
CDataSummary::CDataSummary(string Name, CData *D) {
	int i,j,TotalSites = 0, Char = D->m_iChar;
	double Sum;
	vector <int> States, Pat;
	vector <vector <int> > PatObs;
	vector <vector <double> > SeqFrq;
	m_iNoDataPatterns = 0;
	m_sName = Name;
	m_bDoFullOutput = false;
	// Get the names and other basic data details
	m_vsNames = D->m_vsName; m_DataType = D->m_DataType;
	m_iNoSeq = D->m_iNoSeq; m_iSize = D->m_iTrueSize;
	m_dlnL = m_dlnL_multi = m_ddelta = BIG_NUMBER;

	// Get the distribution of number of characters in sites
	m_vdStatesObs.assign(min(m_iNoSeq,NumStates(D->m_DataType)),0);
	SeqFrq.assign(NumStates(m_DataType),vector<double>()); FOR(i,NumStates(m_DataType)) { SeqFrq[i].assign(m_iNoSeq,0); }
	m_vdFrq.assign(NumStates(m_DataType),0.0);
	m_vdFrqStdErr.assign(NumStates(m_DataType),0.0);
	// Get the info
	FOR(i,D->m_iSize)	{
		States.clear();
		Pat.clear();
		FOR(j,m_iNoSeq) {
			if(D->m_ariSeq[j][i] >= D->m_iChar) { continue; }
			// Do frqs
			SeqFrq[D->m_ariSeq[j][i]][j] += (double) D->m_ariPatOcc[i]; TotalSites += D->m_ariPatOcc[i];
			// Do patterns
			Pat.push_back(D->m_ariSeq[j][i]);
			if(!IsIn(D->m_ariSeq[j][i],States)) { States.push_back(D->m_ariSeq[j][i]); }
		}
		// Get the number of different observed states
		m_vdStatesObs[(int)States.size() - 1] += D->m_ariPatOcc[i];
		// Get the site pattern
		if(!(IsIn(Pat,PatObs))) { PatObs.push_back(Pat); m_iNoDataPatterns++; }
	}
	// Get the frequency information
	FOR(i,m_iNoSeq) { Sum = 0.0; FOR(j,Char) { Sum += SeqFrq[j][i]; }	FOR(j,Char) { SeqFrq[j][i] /= Sum; }	}
//	cout.precision(12); cout << "\nThe sequence frequencies [total =" << TotalSites << "]: ";
//	FOR(i,4) {  cout << "\nstate["<<i<<"]: " << SeqFrq[i]; }
	FOR(i,Char) {
		Sum = 0;
		FOR(j,m_iNoSeq) {
			Sum += (double) SeqFrq[i][j];
			m_vdFrqStdErr[i] += (double) SeqFrq[i][j] * (double) SeqFrq[i][j];
		}
		m_vdFrq[i] = (double) Sum / (double)m_iNoSeq;
		m_vdFrqStdErr[i] = sqrt((m_vdFrqStdErr[i] - (Sum * m_vdFrq[i])) / (m_iNoSeq - 1));
	}
//	cout << "\nValues produced: ";
//	FOR(i,4) { cout << "\nState["<<i<<"]: mean = " << m_vdFrq[i] << "\tStdErr: " << m_vdFrqStdErr[i]; }
}

// Output operator
ostream &operator<<(ostream &os, CDataSummary DS)	{
	int i;
	if(DS.m_dlnL < 0) { if(DS.m_dlnL_multi > 0) { Error("\nValues for likelihoods in CDataSummary don't look right...\n"); } }

	if(DS.m_bDoFullOutput) {
		os << "\nName\tNoSeq\tSize\tNoPatterns";
		FOR(i,NumStates(DS.m_DataType)) { os << "\tNumDiffStates="<<i+1; }
		FOR(i,NumStates(DS.m_DataType)) { os << "\tAveFreq("<<State(DS.m_DataType,i)<<")\tStdErrFreq("<<State(DS.m_DataType,i)<<")"; }
		if(DS.m_dlnL < 0) { os << "\tlnL\tMulti_lnL\tdelta"; }
	}
	os << endl << DS.m_sName << "\t" << DS.m_iNoSeq << "\t" << DS.m_iSize << "\t" << DS.m_iNoDataPatterns << "\t" << DS.m_vdStatesObs;
	FOR(i,NumStates(DS.m_DataType)) { os << "\t" << DS.m_vdFrq[i] << "\t" <<DS.m_vdFrqStdErr[i]; }
	if(DS.m_dlnL < 0) { os << "\t" << DS.m_dlnL << "\t" << DS.m_dlnL_multi << "\t" << DS.m_ddelta; }
	return os;
}
