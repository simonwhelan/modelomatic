// Implementation of the models used in programs
//

#include "model.h"

#if FUNC_COUNTERS == 1
	extern int LFunc_Log_Counter, SubLFunc_Log_Counter;
#endif
extern vector <double> PWDists;		// Pairwise distances

#if DO_MEMORY_CHECK
extern CMemChecker memory_check;
#endif


#define DERIVATIVE_DEBUG 0			// Whether to debug the derivative function
#define MODEL_DEBUG	0				// Turns on some checkers that might throw up where errors are occurring
#define LIKELIHOOD_FUNC_DEBUG 0		// Whether to debug the likelihood function
#define FASTBRANCHOPT_DEBUG 0		// Whether to debug the FastBranchOpt function
#define ALLOW_BRANCH_OPTIMISE 1		// Specifies whether fast branch optimisation is allowed... (Should always be 1!)

////////////////////////////////////////////////////////////////////////////
// CBaseModel implementation
////////////////////////////////////////////////////////////////////////////

// Constructor
CBaseModel::CBaseModel(CData *Data, CTree *Tree,string Name)	{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseModel++;
#endif
	// Set up some storage pointers
	m_pData = Data; m_pTree = Tree; m_sName = Name; m_pSubTree = NULL; m_Pars = NULL;
	m_arL = NULL;  m_bOptReady = false; m_CalcType = cML; m_bCompressedModel = false; m_bDoSepParOpt = false; m_PreOptModel = UNKNOWN;
	m_bMainModel = true; m_bTreeHMM = false; m_bLockModel = false;
	m_bOutputDetail = false; m_iOptNum = DEFAULT_OPTNUM;
	m_ModelRate = new CPar("Rate",1.0,false,1.0E-5);	// Set model rate to 1.0
	m_iFastBralnL_Calls = m_iFastBralnL = m_iFastBralnL_Bracket = 0;
	pLikelihood = NULL;
}
// Destructor
CBaseModel::~CBaseModel()	{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseModel--;
#endif
	pLikelihood = NULL;
	CleanPar();
	CleanMemory();
	DoBralnL(-1,-1,-1,true);		// Clean the static CProb vector
}
// Cleaning functions
// 1.) Clean parameters
void CBaseModel::CleanPar()	{
	int i;
	FOR(i,(int)m_vpPar.size())	{ m_vpPar[i] = NULL; }
	m_vpPar.clear();
}
// 2.)  Clean memory
void CBaseModel::CleanMemory()	{
	int i;
	// Do public
	FOR(i,m_vpProc.size()) { if(m_vpProc[i] != NULL) { delete m_vpProc[i]; } }
	m_vpProc.clear();
	m_vdProbProc.clear();
	m_pData = NULL; m_pTree = NULL; m_pSubTree = NULL;
	DEL_MEM(m_arL); m_sName = "Unassigned";
	// Do private
	m_bOptReady = false;
	FOR(i,(int)m_vpAllOptPar.size()) { m_vpAllOptPar[i] = NULL; }
	m_vbDoBranchDer.clear();
	m_viCPNodesCovered.clear(); m_viLeafMap.clear(); m_viNodeFrom.clear(); m_vdExtBra.clear();
	if(m_ModelRate != NULL) { delete m_ModelRate; m_ModelRate = NULL; }
	// Do parsimony
	if(m_Pars != NULL) { delete m_Pars; }
}

// Space creation
void CBaseModel::CreateProcessSpace(bool Force) { int i; FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->MakeCalcSpace(); } }

////////////////////////////////////////////////////////////
// Stuff for initialising parsimony
ECalcType CBaseModel::DoParsimony(double PropSlowestSites)	{
	if(IsParsimony()) { return m_CalcType; }
	m_CalcType = cMP;
	if(m_Pars == NULL) { m_Pars = new CParsimony(m_pData,m_pTree); }
	else { if(m_Pars->NoSeq() != m_pData->m_iNoSeq || m_Pars->Size() != m_pData->m_iSize) { delete m_Pars; m_Pars = new CParsimony(m_pData,m_pTree); } }
	if(PropSlowestSites > 0) { m_Pars->CreateLimitMask((int)((double)m_pData->m_iNoSeq * PropSlowestSites)); }
	return m_CalcType;
}
// Function to return full parsimony score
double CBaseModel::GetFullParsimony() {
	double Score;
	vector <bool> Mask;
	if(IsParsimony()) {
		m_Pars->GetMask(Mask);
		m_Pars->CreateZeroMask();
		Score = lnL();
		m_Pars->SetMask(Mask);
	} else if(IsLikelihoodCalc()) {
		DoParsimony(BIG_NUMBER);
		Score = lnL();
		DoLikelihoodCalc();
	}
	return Score;
}


// Get a vector of parameters
void CBaseModel::GetParVec(vector <CPar *> *ParVec)	{
	int i,j;
	FOR(i,(int)ParVec->size()) { ParVec->at(i) = NULL; } ParVec->clear();
	FOR(i,(int) m_vpProc.size())	{
		FOR(j,m_vpProc[i]->NoPar())	{
			if(!IsIn(m_vpProc[i]->pPar(j),*ParVec)) { ParVec->push_back(m_vpProc[i]->pPar(j)); }
	}	}
	FOR(i,(int)m_vpPar.size())	{ if(!IsIn(m_vpPar[i],*ParVec)) { ParVec->push_back(m_vpPar[i]); } }

}

// Remove a parameter
void CBaseModel::RemovePar(string Name)	{
	int i;
	vector <CPar *>::iterator iPar;
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->RemovePar(Name); }
	i = 0; 	IFOR(iPar,m_vpPar)	{
		if(m_vpPar[i]->Name().find(Name) != string::npos) { m_vpPar.erase(iPar);  if(iPar == m_vpPar.end()) { break; } }
		else { i++; }
}	}

// Basic function that gets the optimised parameters into the model
vector <CPar *> CBaseModel::CreateOptPar()	{
	int i,j;
	CleanPar();
	FOR(i,(int)m_vpProc.size())	{
		FOR(j,m_vpProc[i]->NoPar())	{
			if(m_vpProc[i]->pPar(j)->Opt() == true) { m_vpPar.push_back(m_vpProc[i]->pPar(j)); }
	}	}
/*
	cout << "\nThe optimised parameters are: ";
	FOR(i,(int)m_vpPar.size())	{ cout << endl << *m_vpPar[i]; }
*/
	return m_vpPar;
}

////////////////////////////////////////////////////////////
// Functions to control the best model

SBestModel CBaseModel::ModelScores(double BestlnL)	{
	int i,j;
	SBestModel Ret;
	vector <bool> DoTree(m_vpProc.size(),true);
	// Check entry conditions if wanted
//	if(fabs(BestlnL - lnL()) > FULL_LIK_ACC) { cout << "\nDifference in likelihood: " << fabs(BestlnL - lnL()) << " Exp: " << BestlnL << " cf. " << lnL(); }
//	assert(fabs(BestlnL - lnL()) < FULL_LIK_ACC);				// Debugging code
	// Do likelihood
	Ret.m_dlnL = BestlnL;
	// Do parameters
	Ret.m_dParVal.clear();
	FOR(i,(int)m_vpPar.size()) { Ret.m_dParVal.push_back(m_vpPar[i]->Val()); }
	// Do trees
	Ret.m_vTree.clear();
	FOR(i,(int)m_vpProc.size())	{
		if(DoTree[i] == false) { continue; }
		Ret.m_vTree.push_back(*m_vpProc[i]->Tree());
		for(j=i+1;j<(int)m_vpProc.size();j++) { if(m_vpProc[i]->Tree() == m_vpProc[j]->Tree()) { DoTree[j] = false; } }
	}
	return Ret;
}

double CBaseModel::RestoreBestModel(SBestModel *M)	{
	int i,j,TreeNum =0;
	vector <bool> DoTree(m_vpProc.size(),true);
	// check entry conditions
	assert(M->m_dParVal.size() == m_vpPar.size());
	// Do parameters
	FOR(i,(int)M->m_dParVal.size()) { if(!m_vpPar[i]->Special()) { m_vpPar[i]->SetVal(M->m_dParVal[i],true,true,false); } }
	FOR(i,(int)M->m_dParVal.size()) { if(m_vpPar[i]->Special()) { m_vpPar[i]->SetVal(M->m_dParVal[i],true,true,true); } }
	// Do trees
	FOR(i,(int)m_vpProc.size())	{
		if(DoTree[i] == false) { continue; }
		assert(TreeNum < (int)M->m_vTree.size());
		*m_vpProc[i]->Tree()= M->m_vTree[TreeNum++];
		for(j=i+1;j<(int)m_vpProc.size();j++) { if(m_vpProc[i]->Tree() == m_vpProc[j]->Tree()) { DoTree[j] = false; } }
	}
	CleanFastCalc(true);
	// Final check
/*	FOR(i,(int)M->m_dParVal.size())	{
		cout << "\nModelPar " << m_vpPar[i]->Name() << " == " << m_vpPar[i]->Val() << " cf. " << M->m_dParVal[i];
		if(fabs(m_vpPar[i]->Val() - M->m_dParVal[i]) > 1.0E-6) {
			 Error("\nModel update error...");
	}	}
*/
	if(fabs(lnL() - M->m_dlnL) > FULL_LIK_ACC * 100) { cout.precision(10); cout << "\nDifference " << lnL() << " cf. " << M->m_dlnL << " == " << lnL() - M->m_dlnL; }
//	assert(fabs(lnL() - M->m_dlnL) < FULL_LIK_ACC * 100); // The FULL_LIK_ACC * 100 is used because reinstating parameters can be messy
	return M->m_dlnL;
}

////////////////////////////////////////////////////////////
// Output function
ostream &operator<<(ostream &os, CBaseModel &Model)	{
	int i;
	// Collect parameters from processes
	vector <CPar *> ParVec; Model.GetParVec(&ParVec);
	vector <double> Values;
	os.precision(6); os.setf(ios::fixed);
	Model.Tree()->OutLabel();
	os << "# Model file for Leaphy;";
	os << "\nModel:\t" << Model.m_sName << ";";
	os << "\nDataType:\t" << Model.m_pData->m_DataType << ";";
	os << "\nData:\t" << Model.m_pData->m_iNoSeq << "\t" << Model.m_pData->m_iTrueSize << ";";
	Model.Tree()->OutBra(); Model.Tree()->OutName();
	os << "\nTree:\t" << *Model.Tree();
	os << "\nTreeLength:\t" << Model.Tree()->GetTreeLength();
	os << "\nOriginal_random_seed: " << Model.GetOriRandomSeed() << ";";
	os << "\nCurrent_random_seed: " << Model.GetCurrentRandomSeed() << ";";
	if(Model.m_pData->Valid()) {
		Values.clear();
		Values.push_back(Model.lnL(true));
		os << "\nLikelihood:\t" << Values[0];
		if(!Model.m_vpAssociatedModels.empty()) {
			FOR(i,(int)Model.m_vpAssociatedModels.size()) { Values.push_back(Model.m_vpAssociatedModels[i]->lnL(true)); Values[0] -= Values[i+1]; }
			os << "\t{"; FOR(i,(int)Values.size()) { os << "\t" << Values[i]; }
			os << "\t}";
		}
		os << ";";
	}
	FOR(i,(int)ParVec.size()) { os << "\nParameter["<<i<<"]\t" << *ParVec[i] << ";"; }
	if(Model.OutputDetail())	{
		os << "\n-----------------------------";
		os << "\nConsists of " << Model.m_vpProc.size() << " process"; if(Model.m_vpProc.size() > 1) { os << "es"; }
		FOR(i,(int)Model.m_vpProc.size())	{ Model.m_vpProc[i]->SetOutputDetail(Model.OutputDetail()); os << "\n" << *Model.m_vpProc[i]; }
	}
	// Do associated models if required...
	FOR(i,(int)Model.m_vpAssociatedModels.size()) {
		os << "\n----------------------------- AssociatedModels["<<i<<"]; Rate: "<< Model.m_vpAssociatedModels[i]->ModelRate() <<" -----------------------------\n";
		os << *Model.m_vpAssociatedModels[i];
	}

	return os;
}

// Output parameter names
vector <string> CBaseModel::ParNameOut()	{
	int i,j;
	vector <string> Names;
	vector <CPar *> ParVec;
	FOR(i,(int)m_vpProc.size())	{
		FOR(j,m_vpProc[i]->NoPar())	{
			if(!IsIn(m_vpProc[i]->pPar(j),ParVec))	{
				Names.push_back(m_vpProc[i]->pPar(j)->Name());
				ParVec.push_back(m_vpProc[i]->pPar(j));
	}	}	}
	FOR(j,NoPar())	{
		if(!IsIn(m_vpPar[j],ParVec))	{
			Names.push_back(m_vpPar[j]->Name());
			ParVec.push_back(m_vpPar[j]);
	}	}
	return Names;
}
// Output the values of parameters
vector <double> CBaseModel::ParValOut()	{
	int i,j;
	vector <CPar *> ParVec;
	vector <double> Vals;
	FOR(i,(int)m_vpProc.size())	{
		FOR(j,m_vpProc[i]->NoPar())	{
			if(!IsIn(m_vpProc[i]->pPar(j),ParVec))	{
				Vals.push_back(m_vpProc[i]->pPar(j)->Val());
				ParVec.push_back(m_vpProc[i]->pPar(j));
	}	}	}
	FOR(j,NoPar())	{
		if(!IsIn(m_vpPar[j],ParVec))	{
			Vals.push_back(m_vpPar[j]->Val());
			ParVec.push_back(m_vpPar[j]);
	}	}
	return Vals;
}

// Space output function
void CBaseModel::OutSpace(ostream &os,int Proc, int Node, int PBeg, int PEnd)	{
	int i;
	if(Proc == -1) {
		FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->SiteOut(os,Node,PBeg,PEnd); }
	}
	else {
		assert(InRange(Proc,0,(int)m_vpProc.size()));
		m_vpProc[Proc]->SiteOut(os,Node,PBeg,PEnd);
	}
}
void CBaseModel::OutPT(ostream &os, int Proc, int Branch)	{
	int i,j;
	assert(InRange(Proc,-1,(int)m_vpProc.size()));
	os << "\nPT output for: Process="; if(Proc == -1) { os << "ALL"; } else { os << Proc; }
	os << ", Branches="; if(Branch == -1) { os << "ALL"; } else { os << Branch; }
	if(Proc == -1)	{
		FOR(i,(int)m_vpProc.size())	{
			os << "\nProcess="<<i;
			assert(InRange(Branch,0,m_vpProc[i]->Tree()->NoBra()));
			if(Branch == -1) {
				FOR(j,m_vpProc[i]->Tree()->NoBra())	{ m_vpProc[i]->OutPT(os,j); }
			} else {
				m_vpProc[i]->OutPT(os,Branch);
			}
	}	} else {
		assert(InRange(Branch,0,m_vpProc[Proc]->Tree()->NoBra()));
		if(Branch == -1) {
			FOR(j,m_vpProc[Proc]->Tree()->NoBra())	{ m_vpProc[Proc]->OutPT(os,j); }
		} else {
			m_vpProc[Proc]->OutPT(os,Branch);
	}	}
}
// Zeroing function
void CBaseModel::ZeroSpace()	{
	int i;
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ZeroSpace(); }
}

///////////////////////////////////////////////////
// Function returning values for optimisation
vector <double *> CBaseModel::GetOptPar(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm)	{
	int i,j, grad_pointer = 0;
	vector <double *> OptVal;
	// Clean parameter space
//	cout << "\n----------- GetOptPar -----------";
	FOR(i,(int)m_vpAllOptPar.size()) { m_vpAllOptPar[i] = NULL; } m_vpAllOptPar.clear();
	m_vbDoBranchDer.clear();
	// Get some memory and space
	m_vbDoBranchDer.assign(m_vpProc.size(),false);
	// Get values for optimising branches
	if(ExtBranch == true) {
		FOR(i,(int)m_vpProc.size()) {
			m_vbDoBranchDer[i] = true;
			// Skip adding extra parameters when trees are the same in different processes
			if(i > 0 && m_vpProc[i]->Tree() == m_vpProc[0]->Tree()) { continue; }
			if(m_vpProc[i]->NoSeq() >2) {
				FOR(j,m_vpProc[i]->Tree()->NoSeq()) {
					m_vpProc[i]->Tree()->SetOptB(j,true);
					m_vpAllOptPar.push_back(m_vpProc[i]->Tree()->pBra(j));
					OptVal.push_back(m_vpProc[i]->Tree()->OptimiserB(j));
			}	} else { // If 2 sequences always branch zero
				m_vpProc[i]->Tree()->SetOptB(0,true);
				m_vpAllOptPar.push_back(m_vpProc[i]->Tree()->pBra(0));
				OptVal.push_back(m_vpProc[i]->Tree()->OptimiserB(0));
			}
	}	}
	if(IntBranch == true)	{
		FOR(i,(int)m_vpProc.size()) {
			if(m_vpProc[i]->NoSeq() == 2) { continue; }
			m_vbDoBranchDer[i] = true;
			// Skip adding extra parameters when trees are the same in different processes
			if(i > 0 && m_vpProc[i]->Tree() == m_vpProc[0]->Tree()) { continue; }
			for(j=m_vpProc[i]->Tree()->NoSeq();j<m_vpProc[i]->Tree()->NoBra();j++)	{
				m_vpProc[i]->Tree()->SetOptB(j,true);
				OptVal.push_back(m_vpProc[i]->Tree()->OptimiserB(j));
				m_vpAllOptPar.push_back(m_vpProc[i]->Tree()->pBra(j));
	}	}	}
	// Get values fo optimising parameters
	if(Parameters == true && !Locked())	{
		FOR(i,(int)m_vpPar.size())	{
			if(m_vpPar[i]->Opt() == false || m_vpPar[i]->Special()) { continue; }
			string::size_type loc = m_vpPar[i]->Name().find("Freq",0);
			if(Eqm == false && loc != string::npos) { continue; }
			OptVal.push_back(m_vpPar[i]->OptimiserValue());
			m_vpAllOptPar.push_back(m_vpPar[i]);
	}	}

/*	cout << "\nGetting Opt Par ["<< m_vpAllOptPar.size();
	cout << ":" << OptVal.size() <<"]";
	FOR(i,m_vpAllOptPar.size()) {
		cout << "\n\tPar["<<i<<"] " << m_vpAllOptPar[i]->Name() << " = " << m_vpAllOptPar[i]->Val() << " == " << *OptVal[i];
	}
	cout << "\n---------"; */
//	exit(-1);
	return OptVal;
}

// Count the number of optimised parameters
int CBaseModel::CountOptPar(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm)	{
	int i,j, grad_pointer = 0, NoPar = 0;
	// Get values for optimising branches
	if(ExtBranch == true) {
		FOR(i,(int)m_vpProc.size()) {
			// Skip when trees are the same in different processes
			if(i > 0 && m_vpProc[i]->Tree() == m_vpProc[0]->Tree()) { continue; }
			FOR(j,m_vpProc[i]->Tree()->NoSeq()) {
				NoPar++;
				if(m_vpProc[i]->Tree()->NoSeq() == 2) { break; }
	}	}	}
	if(IntBranch == true)	{
		FOR(i,(int)m_vpProc.size()) {
			// Skip when trees are the same in different processes
			if(i > 0 && m_vpProc[i]->Tree() == m_vpProc[0]->Tree()) { continue; }
			for(j=m_vpProc[i]->Tree()->NoSeq();j<m_vpProc[i]->Tree()->NoBra();j++)	{
				NoPar++;
	}	}	}
	// Get values fo optimising parameters
	if(Parameters == true && !Locked())	{
		FOR(i,(int)m_vpPar.size())	{
			string::size_type loc = m_vpPar[i]->Name().find("Freq",0);
			if(Eqm == false && loc != string::npos) { continue; }
			NoPar++;
	}	}
	return NoPar;
}
////////////////////////////////////////////////////////////////////////////
// Randomise the optimised parameters
void CBaseModel::RandomiseParameters(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm)	{
	int i;
	vector <double *> Par = GetOptPar(ExtBranch,IntBranch,Parameters,Eqm);
	assert(Par.size() == m_vpAllOptPar.size());
	FOR(i,(int)Par.size())	{
		*Par[i] = RandDouble(0.0,2* max(*Par[i],1.0)); Par[i] = NULL;
	}
}
////////////////////////////////////////////////////////////////////////////
// Function that checks that the memory allocation for the model is viable
bool CBaseModel::IsViable()	{
	int i;
	double Total = 0.0;
	// All processes should be described and have a probability associated with them
	FOR(i,(int)m_vpProc.size()) {
		if(m_vpProc[i] == NULL) { return false; }
		Total += m_vpProc[i]->Prob();
	}
	if(diff(Total,1.0)) { return false; }
	// Check optimised parameters and gradients, &c
//	if((int) m_vpAll
	// All processes
	return true;
}

// Function to fix small branch lengths to ensure reasonable computation
void CBaseModel::FixSmallBranches()	{
	int i,j;
	vector <CTree *> TreesDone;
	FOR(i,(int)m_vpProc.size())	{
		if(!IsIn(m_vpProc[i]->Tree(),TreesDone)) { TreesDone.push_back(m_vpProc[i]->Tree()); } else { continue; }
//		cout << "\nTrying fix small branches (<"<<100*DX<<")...";
//		cout << "\nTree: " << *m_vpProc[i]->Tree();
		FOR(j,m_vpProc[i]->Tree()->NoBra())	{
//			cout << "\nChecking branch " << j << "/" << m_vpProc[i]->Tree()->NoBra()<< " = " << m_vpProc[i]->Tree()->B(j) << flush;
			if(m_vpProc[i]->Tree()->B(j) < 100 * DX) { m_vpProc[i]->Tree()->SetB(j,0.01,true); }
		}
//		cout << "\nNewTree: " << *m_vpProc[i]->Tree();
	}
	FOR(i,(int)TreesDone.size()) { TreesDone[i] = NULL; } TreesDone.~vector();
}

bool CBaseModel::CheckSameTree()	{
	int i;
	FOR(i,(int)m_vpProc.size()) { if(m_vpProc[i]->Tree() != m_vpProc[0]->Tree()) { return false; } }
	return true;
}

/////////////////////////////////////////////////////////////////////////
// Rescale all parameters
void CBaseModel::RedoScale(bool Force)	{
	int i,j;
	vector <CTree *> Trees;
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->Rescale(Force); }
	FOR(i,(int)m_vpProc.size())	{
		FOR(j,(int)Trees.size()) { if(m_vpProc[i]->Tree() == Trees[j]) { break; } }
		if(j != Trees.size()) { continue; }
		Trees.push_back(m_vpProc[i]->Tree());
		FOR(j,m_vpProc[i]->Tree()->NoBra()) {
			m_vpProc[i]->Tree()->BRescale(j);
	}	}
}

/////////////////////////////////////////////////////////////////////////
// Derivative calculations
vector <double> CBaseModel::GetDerivatives(double CurlnL, bool *pOK)	{
	int i;
	bool OK = true, ForceNumBra = false;
	vector <double> Grads;
	vector <double> temp;
	double temp_lnL;
	if(m_vpAllOptPar.empty())  { return Grads; }
	assert(IsViable());
#if MODEL_DEBUG
	double ltemp = lnL();
	cout << "\nIn GetDerivatives: Comparing CurlnL= " << CurlnL << " cf. ReallnL= " << ltemp << " ; diff = " << fabs(CurlnL - ltemp);
	if(fabs(ltemp - CurlnL) > FULL_LIK_ACC) { exit(-1); }
	CurlnL = ltemp;
#else
	// If required, get the current likelihood
//	cout << "\n\n>>>>>>>>>>GetDerivatives: " << CurlnL;CurlnL = lnL(); cout << " cf. " << CurlnL;
//	if(CurlnL - DBL_EPSILON < -BIG_NUMBER) { CurlnL = lnL(); }
#endif
	if(IsRMSDCalc())	{	////////////////////////// Do RMSD derivatives ///////////////////////
		FOR(i,(int)m_vpAllOptPar.size()) {
			Grads.push_back(m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL)));
			if(Grads[i] > RMSD_GRAD_LIM)		{ OK = false; Grads[i] = RMSD_GRAD_LIM; }
			else if(Grads[i] < -RMSD_GRAD_LIM)	{ OK = false; Grads[i] = -RMSD_GRAD_LIM; }
	}	} else {			////////////////////////// Do likelihood derivatives /////////////////
#if ALLOW_ANALYTICAL_BRANCH_DERIVATIVES
//		cout << "\nTrying analytical derivatives...";
		if((int)m_vpAllOptPar.size() > 1) {
			int j;
			// Get the derivatives by process
			temp.assign((int)m_vpAllOptPar.size(),0.0);
			// Initialise
			FOR(i,(int)m_vpAllOptPar.size()) { m_vpAllOptPar[i]->InitialiseDerivativeType(); }

			///////////////////////////////////////////////////////////////////
			// For analytical derivatives
			// 1. Set up the Q matrices
			PreparelnL();
			// 2. Build all the partial likelihoods and get the processes sitewise likelihood
			FOR(i,(int)m_vpProc.size())	{
				if(m_vbDoBranchDer[i] == true) {
					if(m_vpProc[0]->Tree() != m_vpProc[i]->Tree()) { Error("\nHaven't checked that multiple trees work..."); exit(-1); }
					m_vpProc[i]->PrepareBraDer();
			}	}
			if(!FormMixtureSitewiseL()) { ForceNumBra = true; } // Form the mixture distribution
			// 3. Now get the branch derivatives
			FOR(i,(int)m_vpProc.size())	{
				if(m_vbDoBranchDer[i] == true) {
					if(!m_vpProc[i]->GetBraDer(m_arL)) { ForceNumBra = true; break; }
					FOR(j,m_vpProc[i]->Tree()->NoBra())	{ temp[j] += m_vpAllOptPar[j]->grad(); }
			}	}
			if(ForceNumBra == true)	{
				assert(CheckSameTree());
//				if(Tree()->NoSeq() > 2) { cout << "\n <<<<<<<<<<<<<<< DOING NUMERICAL DERIVATIVES: CurlnL = " << CurlnL << "; diff: "<< fabs(CurlnL-lnL()) <<" >>>>>>>>>>>>>>>>>>>>>>>>>>>>"; }
				FOR(i,m_vpProc[0]->Tree()->NoBra()) {
					assert(m_vpAllOptPar[i]->IsBranch());
//					cout << "\nGetting branch derivative for ["<<i<<"]: " << *m_vpAllOptPar[i];
					temp[i] = m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL));
//					cout << " ... have " << m_vpAllOptPar[i]->grad() << " = " << temp[i];
				}
				// Error check here
				temp_lnL = lnL(true);
				if(fabs(temp_lnL - CurlnL) > 0.00001) {
					// If there's an error try one more time
					CurlnL = temp_lnL;
					FOR(i,m_vpProc[0]->Tree()->NoBra()) {
						assert(m_vpAllOptPar[i]->IsBranch());
						temp[i] = m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL));
					}
					// I shouldn't let it continue, but I'll try...
					if(fabs(temp_lnL - CurlnL) > 0.001) { cout.precision(10); cout << "\nError in CModel::GetDerivatives(...): likelihoods don't match lnL()= " << CurlnL << " cf. "<< lnL() << " cf. " << lnL(true); exit(-1); }
				}
			}
			// Copy the derivatives to the store
			FOR(i,(int)m_vpAllOptPar.size()) { if(m_vpAllOptPar[i]->IsBranch()) { m_vpAllOptPar[i]->grad(temp[i]); }  }
			// Get the remaining numerical derivatives
//			cout << "\nThink likelihood should be " << CurlnL << " cf. " << lnL() << " == " << fabs(lnL() - CurlnL);
			FOR(i,(int)m_vpAllOptPar.size())	{
				// If it requires a numerical derivative then do it
				if(m_vpAllOptPar[i]->DoNumDer()) {
//					cout << "\n================== Getting derivative for " << *m_vpAllOptPar[i] << " ==============";
//					cout << "\nCurrent lnL: " << lnL() << " cf. " << CurlnL;
					m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL));
				}

			}
		} else {	// Always do single derivatives numerically
			m_vpAllOptPar[0]->grad(GetNumDerivative(m_vpAllOptPar[0]->OptimiserValue(),CurlnL));
		}
#else
		FOR(i,(int)m_vpAllOptPar.size())	{ m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL)); }
#endif
		// Store the gradients
		FOR(i,(int)m_vpAllOptPar.size()) { Grads.push_back(m_vpAllOptPar[i]->grad()); }
	}

//	static int Count =0;

	if(pOK != NULL) { *pOK = OK; }
#if DERIVATIVE_DEBUG == 1 && ALLOW_ANALYTICAL_BRANCH_DERIVATIVES == 1
			cout << "\nDerivative checker:";
			Tree()->OutBra(); Tree()->OutName(); cout << "\nTree: " << *Tree();
			double dt_check;
			FOR(i,(int)m_vpAllOptPar.size())	{
				cout << "\nCurlnL: " << CurlnL << " cf. real: " << lnL(true);
				dt_check = GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL);
				cout << "\n\tPar["<<i<<"] " << m_vpAllOptPar[i]->Name() << " = " << *m_vpAllOptPar[i]->OptimiserValue() << " -- Der: " << Grads[i] << " cf. numder " << dt_check << " diff= " << fabs(dt_check - Grads[i]);
				if( (dt_check > 0 && Grads[i] < 0) || (dt_check < 0 && Grads[i] > 0) ) { cout << " inv_sign!";  }
//				exit(-1);
			}
			if(fabs(CurlnL - lnL()) > FLT_EPSILON) { cout << "\n\n<<<<<<<<<<<< WARNING LIKELIHOOD ERROR: CurlnL = " << CurlnL << " cf. actual = " << lnL() << " >>>>>>>>>>>>>>>>>>>>>>\n\n"; }
#endif

//	cout << "\nGradients: " << Grads << "\n";
//	if(Count ++ > 4)
//		exit(-1);
	return Grads;
}


double CBaseModel::GetNumDerivative(double *x, double Old_lnL)	{
	double grad = 0.0, OldPar = *x, UseDX = max(DX,*x * DX), new_lnL, new_par;
	bool UpBound = false, LowBound = false;
	Old_lnL = -fabs(Old_lnL);		// Correct the likelihood (may not be necessary...).
	if(OldPar < 0) { UseDX *= -1; }
	// Get the derivative
#if DERIVATIVE_DEBUG == 1
	if(fabs(Old_lnL -lnL()) > 1.0E-5)	{ cout << "\nNumerical derivative starting; ori_p: " << OldPar << " == lnL: " << Old_lnL << " != actual lnL: " << lnL() << " diff: " << fabs(Old_lnL - lnL()); Error("Hmm");}
	cout << "\nNumerical derivative: ori_p: " << OldPar << " == " << Old_lnL;
	if(fabs(fabs(Old_lnL) - fabs(lnL(true))) > FLT_EPSILON) { cout << "(unstable: " << fabs(Old_lnL - lnL(true)) << ")"; }
#endif
	if(OldPar < 0.5)	{ *x = OldPar + UseDX; }
	else				{ *x = OldPar * (1 + UseDX); }
	new_lnL = lnL();
//#if DERIVATIVE_DEBUG == 1
//	cout << " -> new_p: " << *x << " == " << new_lnL; *x = OldPar; cout << " ==> " << lnL();
//#endif
	// If at the upper bound then sort it out
	if(fabs(new_lnL - Old_lnL) < FLT_EPSILON || fabs(OldPar - *x) < FLT_EPSILON)	{
//		cout.precision(12); cout << "\nDifficult derivative: " << OldPar;
		UseDX *= 2;
		while(fabs(UseDX) < 0.99)	{
			if(fabs(OldPar) < 0.5)	{ new_par = *x = OldPar + UseDX; }
			else { new_par = *x = OldPar * (1 + UseDX); }
//			cout << "\n\tUseDX = " << UseDX << "; OldPar: " << OldPar << "; New par: " << *x;
			new_lnL = lnL();
//			cout << " -> " << *x << ": " << new_lnL;
			// Correct is it's a boundary problem (e.g. original value was past bound)
			if(fabs(*x - new_par) > FLT_EPSILON)	{
				*x = OldPar;
#if DERIVATIVE_DEBUG == 1
					if(fabs(lnL() - Old_lnL) > 1.0E-6) {
						cout.precision(7);
						cout << "\nDiff(" << lnL() << " cf. " << Old_lnL << " = " << fabs(lnL() - Old_lnL);
						Error("\nIn CBaseModel::GetNumDerivative the likelihoods aren't matching...\n\n");
					}
#endif
				OldPar = *x;
			}
			if(fabs(new_lnL - Old_lnL) > FLT_EPSILON) {	break; }
			if(fabs(OldPar) < 0.5) { UseDX *= -2; } // Switch the sign regularly to try and ensure as much of parameter space as possible is investigated }
			else { UseDX *= 2; }
		}
		if(fabs(new_lnL - Old_lnL) < FLT_EPSILON || fabs(OldPar - *x) < FLT_EPSILON) {
//			cout << "\nWarning: flat gradient located... lnL diff=" << fabs(new_lnL - Old_lnL) << "; par diff=" << fabs(OldPar - *x);
//			int i; cout << "\nNew: " << *x << " = " << lnL(); FOR(i,m_vpPar.size()) {
//			cout << "[" << m_vpPar[i]->LowBound() << "," << m_vpPar[i]->UpBound() << "]";
	//		*x = OldPar; cout << "\nOld: " << *x << " = "<< lnL();cout << " ... " << m_vpPar[i]->Name() << " = " << m_vpPar[i]->Val();
			*x = OldPar; return 0.0;
		}
		if(fabs(OldPar - *x) < FLT_EPSILON)	{ Error("\nTrying to do derivative when parameters are the same...\n\n"); }
		grad = (new_lnL - Old_lnL) / (*x - OldPar);
	} else { // Otherwise get the normal gradient
		if(fabs(OldPar - *x) < FLT_EPSILON)	{ Error("\nTrying to do derivative when parameters are the same...\n\n"); }
		grad = (Old_lnL - new_lnL) / (*x - OldPar);
#if DERIVATIVE_DEBUG == 1
		cout << "\n\tPar: " << OldPar << " *x: " << *x <<"; UseDX: " << UseDX << " -> Grad: (" << Old_lnL << " - " << new_lnL << ") / " << *x - OldPar << " == " << grad;
#endif
	}
	*x = OldPar;
#if DERIVATIVE_DEBUG == 1
	if(diff(Old_lnL,lnL()))	{ cout << "\nNumerical derivative ending; ori_p: " << OldPar << " == lnL: " << Old_lnL << " != actual lnL: " << lnL(); }
#endif
	return grad;
}

//////////////////////////////////////////////
// Likelihood calculations

// Function that does any extra initialisations required before calculations
void CBaseModel::FinalInitialisation()		{
	CreateOptPar();			// Create the optimised parameters
	// Get the probability vector sorted
	if(m_vdProbProc.empty()) { m_vdProbProc.assign((int)m_vpProc.size(),1.0 / (double) m_vpProc.size());  }
	assert(m_vdProbProc.size() == m_vpProc.size());
	// Create the partial likelihood
	GET_MEM(m_arL,CProb,m_pData->m_iSize);
}

// Function that prepares the space in processes for quick computation
void CBaseModel::PrepareFastCalc()	{
	if(!ALLOW_FAST_CALC) { return; }
	vector <int> C;
	int i,j;
	// Return for cut trees and for when there is no tree, which occurs in special models
	FOR(i,(int)m_vpProc.size()) { if(m_vpProc[i]->Tree() == NULL) { return; } if(m_vpProc[i]->Tree()->IsCutTree()) { return; } }
	FOR(i,(int)m_vpProc.size()) {
		// Check whether tree has changed for processes
		if(m_vpProc[i]->Tree()->FastCalcOK() == false) {
			FOR(j,(int)m_vpProc.size()) { if(m_vpProc[i]->Tree() == m_vpProc[j]->Tree()) { m_vpProc[j]->DecompressSpace(); } }
		}
		// Compress space
		if(m_vpProc[i]->IsCompressed() == true) { continue; }
		// If required then do all trees in processes with this memory address
		if(!m_vpProc[i]->IsSubTree() && !m_vpProc[i]->MainTree()->IsCutTree()) {
			C = GetCompressedData(m_vpProc[i]->MainTree(),m_vpProc[i]->MainData());
			m_vpProc[i]->PrepareFastCalc(&C);
			// Do all other processes with same tree
			for(j=i+1; j < (int)m_vpProc.size(); j++)	{
				if(m_vpProc[j]->MainTree() == m_vpProc[i]->MainTree() ) { m_vpProc[j]->PrepareFastCalc(&C); }
			}
			m_vpProc[i]->MainTree()->SetFastCalcOK(true);
	}	}
	FOR(i,(int)m_vpAssociatedModels.size()) {
		assert(m_vpAssociatedModels[i]->m_vpAssociatedModels.empty());
		m_vpAssociatedModels[i]->PrepareFastCalc();
	}
}

void CBaseModel::CleanFastCalc(bool force)	{
	int i;
	if(!ALLOW_FAST_CALC) { return; }
	FOR(i,(int)m_vpProc.size()) {
		if(m_vpProc[i]->IsCompressed() == true && m_vpProc[i]->Tree()->FastCalcOK() == false
			&& force == false && IsSubTree() == false) { continue; }
		m_vpProc[i]->DecompressSpace();
	}
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->MainTree()->SetFastCalcOK(false); }
	FOR(i,(int)m_vpAssociatedModels.size()) {
		assert(m_vpAssociatedModels[i]->m_vpAssociatedModels.empty());
		m_vpAssociatedModels[i]->CleanFastCalc(force);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Functions for likelihood computations

// Function to prepare likelihood computations
void CBaseModel::PreparelnL(bool ForceRemake)	{
	int i;
	// Do the rate
//	cout << "\nRate " << flush;
	ModelRate();
	// Go through and apply all global parameters
//	cout << "\nPars " << flush;
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->GlobalApply(); }
	// Prepare for likelihood computations
//	cout << "\nProcesses" << flush;
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->PrepareLikelihood(true,ForceRemake); }
}

// Function that performs likelihood computations
double CBaseModel::lnL(bool ForceReal)	{
	int i;
	double logL = 0.0, TotalProb = 0.0;

//	cout << "\n\n========== Starting new likelihood calc ================" << flush;
	// New stuff dealing with submodels
	if(!m_vpAssociatedModels.empty()) {
//		cout << "\n>>> New likelihood computation";
		FOR(i,(int) m_vpAssociatedModels.size()) {
			// Fix cases where I can assign a tree ok
			if(m_vpAssociatedModels[i] == NULL) { Error("Associated model NULL...\nProbably a problem with TreeHMMs...\n"); }
//			cout << "\n\tlnL Ass["<<i<<"]: " << m_vpAssociatedModels[i]->lnL(ForceReal);
			logL += m_vpAssociatedModels[i]->lnL(ForceReal);
	}	}
	if(!m_pData->Valid()) { Error("\nTrying to do likelihood computation with invalid data\n\n"); }
	if(IsRMSDCalc())	{	////////////////// Do RMSD calc ////////////////////////////////
		if(IsSubTree()) { logL += -m_pSubTree->SubTreeRMSD(m_viLeafMap,m_vviNodesBelow,PWDists,m_vdPartialTreeDist);	}								// DoSubTree calc
		else			{ logL += - m_pTree->RMSD(PWDists); }	// Do normal calc
	} else if(CBaseModel::IsParsimony())	{
		assert(m_Pars != NULL);
		logL += -m_Pars->Score();
	} else {			////////////////// Do full likelihood calculations /////////////
		assert(m_pTree != NULL && m_pData != NULL);
		if(m_pData->m_iNoSeq != m_pTree->NoSeq()) { Error("\nMismatch between number of sequences in data and tree\n\n"); }
		// Prepare the Q mats
//		cout << "\nPreparing Likelihood" << flush;
		PreparelnL();
//		cout << "\nDone PreparelnL" << flush;
		// Get all of the processes ready individually and calculate likelihoods
		FOR(i,(int)m_vpProc.size()) {
#if MODEL_DEBUG == 1
			TotalProb += m_vpProc[i]->Prob();
#endif
			m_vpProc[i]->Likelihood(ForceReal);
 		}
		// Get probabilities of the processes; TODO: This is currently really basic
		if(!FormMixtureSitewiseL()) { cout << "\nReturning error..."; return -BIG_NUMBER; }
		// Get the log likelihood and return it
		FOR(i,m_pData->m_iSize) {
#if MATCH_PAML == 0
			if(m_arL[i].IsZero()) { return -BIG_NUMBER; }
#endif
#if DEVELOPER_BUILD == 1
			if(i < 10000) { cout << "\n\t\tSite["<<i<<"]: " << m_arL[i].m_dValue << " * 10^-" << m_arL[i].m_iScale << " == " << m_arL[i].LogP();}
#endif
//			if(i<30) { cout << "\n\tfullSite["<<i<<"]:\t" << m_arL[i] << " -> " << m_arL[i].LogP(); }
			logL += m_pData->m_ariPatOcc[i] * m_arL[i].LogP();
		}
#if FUNC_COUNTERS == 1
		LFunc_Log_Counter++;
		if(IsSubTree()) { SubLFunc_Log_Counter++; }
#endif
	}
#if MODEL_DEBUG == 1
	cout << "\n\tReturning likelihood: " << logL;
	if(fabs(1.0 - TotalProb) > 1.0E-4) { cout << "\nError in process probabilities, sum = " << TotalProb; FOR(i,m_vpProc.size()) { cout << "\n\tProc["<<i<<"] = " << m_vpProc[i]->Prob();  } cout << "\n\n"; exit(-1); }
#endif
//	cout << "\nRoot: " << m_vpProc[0]->RootEqm();
#if DEVELOPER_BUILD == 1
	cout << "\n\tReturning likelihood: " << logL << endl;  // exit(-1);
#endif
	// Apply extra stuff to the likelihood function if needed
	if(pLikelihood != NULL) {
		logL -= pLikelihood(NULL); // Called as blank. Other arguments intended to allow functionality
	}
//	if(m_bMainModel) { cout << "\n\tReturning likelihood: " << logL << endl;  }
//	cout << "\n\tReturning likelihood: " << logL << endl;  exit(-1);
	return logL;
}


////////////////////////////////////////////////////////////////////////
// Space updating functions
void CBaseModel::Leaf_update(int NTo, int NFr, int Br, CTree *T, int First,bool DoFullUpdate)	{
	int i; FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->LeafNode_Update(NTo,NFr,Br,T,First,DoFullUpdate);
}	}

void CBaseModel::Bran_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo, bool DoNFr, bool DoFullUpdate)	{
	int i; FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->BranNode_Update(NTo,NFr,Br,T,First,DoNTo,DoNFr,DoFullUpdate);
}	}

// Function to produce the sitewise likelihood
// -------------------------------------------
// Loops through the individual processes and adds their likelihoods to m_ardL and m_ariLScale
// This is equivalent to calculating them from a mixture model
// -- This function will be changed for HMMs and other dependent structures --
bool CBaseModel::FormMixtureSitewiseL()	{
	int i,ProcNum;
	// Check some entry conditions
	assert(m_arL != NULL);
#if LIKELIHOOD_FUNC_DEBUG == 1
	cout << "\n--- Forming mixture model from sitewise likelihoods --- ";
#endif
	// Do the looping through processes and then by data size
	FOR(ProcNum,(int)m_vpProc.size())	{
#if LIKELIHOOD_FUNC_DEBUG == 1
		cout << "\nProcess["<<ProcNum<<"]: Prob = " << m_vpProc[ProcNum]->Prob();
		//FOR(i,min(m_pData->m_iSize,5)) { cout << "\n\tSite["<<i<<"]: " << m_vpProc[ProcNum]->L(i) << " = " << m_vpProc[ProcNum]->L(i).LogP(); }
		i = 120;  cout << "\n\tSite["<<i<<"]: " << m_vpProc[ProcNum]->L(i) << " = " << m_vpProc[ProcNum]->L(i).LogP();
#endif
		if(m_vpProc[ProcNum]->Prob() < MIN_PROB) { continue; }
		if(ProcNum == 0) {	// For the first process transfer values
			FOR(i,m_pData->m_iSize)	{ m_arL[i] = m_vpProc[ProcNum]->L(i); }
		} else {			// Otherwise they need to be summed
			FOR(i,m_pData->m_iSize) { m_arL[i].Add(m_vpProc[ProcNum]->L(i),true); }
	}	}
	return true;
}

// Function to calculate the overall log-likelihood from a single process
double CBaseModel::CalculateL(CBaseProcess *Process, bool DoFulllnL)		{
	int i;
	double Sum = 0;
	vector <double> Liks;
	Liks = SitewiseL(Process,DoFulllnL); assert(Liks.size() == Process->Size());
	FOR(i,Liks.size()) {
//		cout << "\nSite[" << i << "]: " << Liks[i] << " * " << Process->PatOcc(i) << " = " << Process->PatOcc(i) * Liks[i];
		Sum += Process->PatOcc(i) * Liks[i]; }
	return Sum;
}
// Function to calculate the sitewise log-likelihoods from a single process
vector <double> CBaseModel::SitewiseL(CBaseProcess *Process, bool DoFulllnL) 	{
	int i;
	vector <double> SitelnL;
	if(DoFulllnL) { Process->Likelihood(true); }
	FOR(i,Process->Size()) { SitelnL.push_back(Process->L(i).LogP()); }
	return SitelnL;
}

/////////////////////////////////////////////////////////////////////////
// Functions that do the centre pointing for a model and tree
// and allow subtree calculations

// Calculates partial likelihoods from a BRANCH centre point and puts them in forward and backward space. Returns number of partial likelihoods prepared
vector <int> CBaseModel::PrepareBranchCP(int BranchCP, int Depth)	{
	int i,j;
	vector <int> temp;
	// Clean in preparation of new subtree
	CleanCPMapping();
	m_viLeafMap = m_pTree->BranchCP(BranchCP,Depth,&m_viNodeFrom,&m_viCPNodesCovered,&m_vdExtBra);
	assert(m_viLeafMap.size() == m_viNodeFrom.size());
	// Prepare the model for computation
	PreparelnL();
	// Get the leaf nodes below (including that node if rqd) for each of the m_viLeafMap;
	FOR(i,(int)m_viLeafMap.size()) {
		temp.clear();
		if(m_viLeafMap[i] < m_pTree->NoSeq()) { temp.push_back(m_viLeafMap[i]); m_vviNodesBelow.push_back(temp); continue; }
		// Find node LeafMap came from
		FOR(j,(int)m_viCPNodesCovered.size()) { if(m_pTree->IsNodeLink(m_viCPNodesCovered[j],m_viLeafMap[i]) ) { break; } }
		assert(j != m_viCPNodesCovered.size());
		m_pTree->GetBraSets(m_viLeafMap[i],m_viCPNodesCovered[j],&temp);
		sort(temp.begin(),temp.end());
		m_vviNodesBelow.push_back(temp);
	}
	// Apply the CP mapping
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom); }
	// Return the LeafMap in case its needed
	return m_viLeafMap;

}
/*vector <int> CBaseModel::PrepareBranchCP(int BranchCP, int Depth)	{
	int i,j;
	vector <int> temp;
	// Clean in preparation of new subtree
	CleanCPMapping();
	temp = m_pTree->BranchCP(BranchCP,Depth,&m_viNodeFrom,&m_viCPNodesCovered,&m_vdExtBra);
	WriteLeafMap(temp);
	assert(m_viLeafMap.size() == m_viNodeFrom.size());
	// Prepare the model for computation
	PreparelnL();
	// Get the leaf nodes below (including that node if rqd) for each of the m_viLeafMap;
	FOR(i,(int)m_viLeafMap.size()) {
		temp.clear();
		if(m_viLeafMap[i] < m_pTree->NoSeq()) { temp.push_back(m_viLeafMap[i]); m_vviNodesBelow.push_back(temp); continue; }
		// Find node LeafMap came from
		FOR(j,(int)m_viCPNodesCovered.size()) { if(m_pTree->IsNodeLink(m_viCPNodesCovered[j],m_viLeafMap[i]) ) { break; } }
		assert(j != m_viCPNodesCovered.size());
		m_pTree->GetBraSets(m_viLeafMap[i],m_viCPNodesCovered[j],&temp);
		sort(temp.begin(),temp.end());
		m_vviNodesBelow.push_back(temp);
	}
	// Apply the CP mapping
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom); }
	// Return the LeafMap in case its needed
	return m_viLeafMap;
}
*/
// Calculates partial likelihoods from a NODE centre pointand puts them in forward and backward space. Returns number of partial likelihoods prepared
vector <int> CBaseModel::PrepareNodeCP(int NodeCP, int Depth)		{
	int i,j;
	vector <int> temp;
	// Clean in preparation of new subtree
	CleanCPMapping();
	temp = m_pTree->NodeCP(NodeCP,Depth,&m_viNodeFrom,&m_viCPNodesCovered,&m_vdExtBra);
	WriteLeafMap(temp);
	assert(m_viLeafMap.size() == m_viNodeFrom.size() && m_viCPNodesCovered.size() == m_viNodeFrom.size() - 2);
	// Get the leaf nodes below (including that node if rqd) for each of the m_viLeafMap;
	FOR(i,(int)m_viLeafMap.size()) {
		temp.clear();
		if(m_viLeafMap[i] < m_pTree->NoSeq()) {
			temp.push_back(m_viLeafMap[i]);
			m_vviNodesBelow.push_back(temp);
			continue;
		}
		// Find node LeafMap came from
		FOR(j,(int)m_viCPNodesCovered.size()) { if(m_pTree->IsNodeLink(m_viCPNodesCovered[j],m_viLeafMap[i]) ) { break; } }
		assert(j != m_viCPNodesCovered.size());
		m_pTree->GetBraSets(m_viLeafMap[i],m_viCPNodesCovered[j],&temp);
		sort(temp.begin(),temp.end());
		m_vviNodesBelow.push_back(temp);
	}
	// Apply to the associated models
	FOR(i,(int)m_vpAssociatedModels.size()) {
		m_vpAssociatedModels[i]->m_viNodeFrom = m_viNodeFrom;
		m_vpAssociatedModels[i]->m_viCPNodesCovered = m_viCPNodesCovered;
		m_vpAssociatedModels[i]->m_vdExtBra = m_vdExtBra;
		m_vpAssociatedModels[i]->m_vviNodesBelow = m_vviNodesBelow;
	}
	if(IsLikelihoodCalc())	{
		// Prepare the model for computation
		PreparelnL();
		// Apply the CP mapping
		FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom,FlipBool(IsRMSDCalc())); }
		// And for associated models
		FOR(i,(int)m_vpAssociatedModels.size()) {
			m_vpAssociatedModels[i]->PreparelnL();
			FOR(j,(int)m_vpAssociatedModels[i]->m_vpProc.size()) { m_vpAssociatedModels[i]->m_vpProc[j]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom,FlipBool(IsRMSDCalc())); }
	}	}
	// Now calculate the m_vdPartialTreeDist
	m_vdPartialTreeDist = m_pTree->GetPartialTreeDist(m_viLeafMap,m_vviNodesBelow);
	FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->m_vdPartialTreeDist = m_vdPartialTreeDist; }
	// Return the LeafMap in case its needed
	return m_viLeafMap;
}

// Prepare pairwise calculations
void CBaseModel::PreparePairwiseCalc(int Seq1, int Seq2, CTree *Tree)	{
	vector <int> Map;
	// Do the mapping, &c
	CleanCPMapping();
	Map.push_back(Seq1); Map.push_back(Seq2);
	WriteLeafMap(Map);
	int i,j;
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom); }
	FOR(j,(int) m_vpAssociatedModels.size()) { 	FOR(i,(int)m_vpAssociatedModels[j]->m_vpProc.size()) { m_vpAssociatedModels[j]->m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom); } }
	ApplySubTree(Tree);
	// Now have a guess at the distance using Jukes cantor style distances
	Tree->SetB(0,m_pData->PoissonDist(Seq1,Seq2),true,true);
}

// Prepare the triplet computations
void CBaseModel::PrepareTripletCalc(vector <int> LeafMap, CTree *T, int SPR)	{
	int i;
	// Do the mapping, &c.
	CleanCPMapping();
	assert(LeafMap.size() == 3);
	// If it is using SPR, then there will be partial likelihoods in sequence 2;
	if(SPR != -1) { LeafMap[2] = -1; }
	WriteLeafMap(LeafMap);
	// Apply this to all the processes
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyCPMapping(m_viLeafMap,m_viNodeFrom,false); }
	// Do for all associated models
	FOR(i,(int)m_vpAssociatedModels.size()) {
		m_vpAssociatedModels[i]->PrepareTripletCalc(LeafMap,T,SPR);
	}
}

// Prepare the SPR computations; returns the spare node to be used for computation
int CBaseModel::PrepareSPR(int Br, int LNum, int *OriBr, LazyType DoLazy)	{
	int i,SpareNode, N;
	vector <int> vI;
	CTree *T = m_pTree;
	vector <int> Map;
	double BestlnL;
	assert(InRange(Br,0,T->NoBra()) && InRange(LNum,0,2) && !T->IsCutTree());
	if(T->BraLink(Br,LNum) < T->NoSeq()) { cout << "\nTree details: "; m_pTree->OutDetail(); Error("\nDoing SPR on a single sequence -- Don't be daft..."); }
	// Get a sequence to use as a leaf
	T->GetBraSets(T->BraLink(Br,FlipBin(LNum)),T->BraLink(Br,LNum),&vI);
	// Only when the subtree has at least 3 sequences is SPR considered
	if(vI.size() < 3) { return -1; }
	// Find where calculations should start
	sort(vI.begin(),vI.end());
	FOR(i,T->NoSeq()) { if(!IsIn(i,vI)) { SpareNode = i; break; } }
	// Prepare LeafMap
	CleanCPMapping();
	FOR(i,T->NoSeq()) { Map.push_back(i); } Map[SpareNode] = -1;
	WriteLeafMap(Map);
	// Prepare the partial likelihoods
	if(IsLikelihoodCalc()) { FOR(i,(int)m_vpProc.size())	{
		m_vpProc[i]->ApplyCPMapping(m_viLeafMap,vI,false);
		m_vpProc[i]->PreparePartialL(T->BraLink(Br,LNum),T->BraLink(Br,FlipBin(LNum)),2); }
	}
	T->SetStartCalc(vI[0]);
	assert(i!=T->NoSeq());
	// Now cut the branch and prepare the subtree
	*OriBr = T->CutBranch(Br,LNum);
	// Finally make sure the node to add has correct structure for adding sequences
	// 1. Find spare branches and prepare the node to be added
	vI.clear();
	FOR(i,T->NoNode())	{ if(T->NodeLink(i,0) == -1) { N = i; break; } } // Find node
	FOR(i,T->NoBra())	{ if(T->BraLink(i,0) == -1)	 { vI.push_back(i); } } // Find branches
	assert(vI.size() >= 2);
	T->ReplaceNodeLinkElement(SpareNode,0,N); T->ReplaceNodeBraElement(SpareNode,0,vI[0]);
	T->ReplaceBraLink(vI[0],0,SpareNode); T->ReplaceBraLink(vI[0],1,N);
	T->ReplaceNodeLinkElement(N,0,SpareNode); T->ReplaceNodeBraElement(N,0,vI[0]);
	T->ReplaceNodeBraElement(N,1,vI[1]); T->ReplaceBraLink(vI[1],0,N);
	T = NULL;
	// Now optimise the cut branch...
//	cout << "\nSPR branch opt... OriBr:" << *OriBr;
//	cout << " ... From " << BestlnL;
	BestlnL = lnL();
	if(DoLazy == fullopt) { SingleBranchOpt(*OriBr,&BestlnL,0.001); }
	if(DoLazy == lazy) { SingleBranchOpt(*OriBr,&BestlnL,1.0); }

//	cout << " -> " << BestlnL;
	return SpareNode;
}

// Calculate P(t) matrices for a specific branch
void CBaseModel::PreparePT(int Br)	{
	assert(InRange(Br,0,Tree()->NoBra()));
	int i; FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->Make_PT(Br,true); }
}

int CBaseModel::NumLeafCP()	{
	if(m_viLeafMap.empty()) { return -1; }	// Return error value
	return (int)m_viLeafMap.size();
}

// Goes through all the processes and removes the centre-point mapping
void CBaseModel::CleanCPMapping()	{
	int i; FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->CleanCPMapping(); }
	m_pSubTree = NULL;
	m_viCPNodesCovered.clear();
	m_vdPartialTreeDist.clear();
	m_viLeafMap.clear(); m_viNodeFrom.clear(); m_vviNodesBelow.clear();
	m_vdExtBra.clear();
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->CleanSubTree(); assert(m_vpProc[i]->Tree() == m_pTree); }
	FOR(i,(int)m_vpAssociatedModels.size()) { assert(m_vpAssociatedModels[i]->m_vpAssociatedModels.empty()); m_vpAssociatedModels[i]->CleanCPMapping(); }
}

// Applies a LeafMap to model and associated models
void CBaseModel::WriteLeafMap(vector <int> NewLeafMap) {
	int j;
	m_viLeafMap = NewLeafMap;
	FOR(j,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[j]->m_viLeafMap = m_viLeafMap; }
}

// Once correct leaf mapping is done, this will apply a new sub-tree to all processes
void CBaseModel::ApplySubTree(CTree *Tree, bool UseExtBra,bool Overwrite)	{
	int i,j;
	assert(!m_viLeafMap.empty());
	// Copy the trees
	m_pSubTree = Tree;
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->CleanSubTree(); assert(m_vpProc[i]->Tree() == m_pTree); }
	if(Tree->NoSeq() == 2) { UseExtBra = false; }
	if(UseExtBra)	{	// Use external branches to initialise tree
		assert(Tree->NoSeq() == m_vdExtBra.size());
		FOR(i,Tree->NoSeq())	{ Tree->SetB(i,max(m_vdExtBra[i],1.0E-5),true,true); }
		for(i=Tree->NoSeq();i<Tree->NoBra();i++) { Tree->SetB(i,DELTA_STEP_SUB,true,true); }
	}
	if(Overwrite)	{	// Actually overwrite the tree
		m_pTree->ReplaceTreeCP(Tree,m_viLeafMap, m_viCPNodesCovered);
	} else {					// Just apply as a subtree
		FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplySubTree(Tree); }
	}
	PrepareFastCalc();
	FOR(j,(int)m_vpAssociatedModels.size()) {
		assert(m_vpAssociatedModels[j]->m_vpAssociatedModels.empty());
		if(m_vpAssociatedModels[j]->Tree() != m_pTree) { m_vpAssociatedModels[j]->ApplySubTree(Tree,UseExtBra,Overwrite); }
		else {
			FOR(i,(int)m_vpAssociatedModels[j]->m_vpProc.size()) { m_vpAssociatedModels[j]->m_vpProc[i]->CleanSubTree(); assert(m_vpAssociatedModels[j]->m_vpProc[i]->Tree() == m_vpAssociatedModels[j]->m_pTree); }
			FOR(i,(int)m_vpAssociatedModels[j]->m_vpProc.size()) { m_vpAssociatedModels[j]->m_vpProc[i]->ApplySubTree(Tree); }
			m_vpAssociatedModels[j]->PrepareFastCalc();
		}
	}
}

bool CBaseModel::IsSubTree()	{
	if(m_pSubTree == NULL) { return false; } return true;
}

// Build the original subtree (for comparison purposes)
void CBaseModel::BuildOriSubTree(CTree *T)	{ m_pTree->BuildOriSubTree(T,m_viLeafMap,m_viCPNodesCovered,m_viNodeFrom); }

////////////////////////////////////////////////////////////////
// Replace the tree
void CBaseModel::MakeNewTree(CTree *T, bool Overwrite)	{
	int i;
	// Copy over memory
	if(!Overwrite) {
		if(m_pTree != NULL) { m_pTree = NULL; }
		m_pTree = T;
		if(m_bTreeHMM) { FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->MakeNewTree(m_pTree,Overwrite); } }
		FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyNewTree(m_pTree); }
	} else { // Physically write over the tree
		assert(m_pTree != NULL && T != NULL);
		*m_pTree = *T;
		if(m_bTreeHMM) { FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->MakeNewTree(m_pTree,Overwrite); } }
		FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->MainTree(m_pTree); }
	}
}
///////////////////////////////////////////////////////////////
// Replace the data
void CBaseModel::MakeNewData(CData *D, bool Overwrite)	{
	int i;
	// Copy over memory
	if(D == NULL) {	m_pData = NULL; }
	else if(!Overwrite) {
//		cout << "\nNot overwriting: " << D->m_iNoSeq << ":" << D->m_iSize << flush;
		if(m_pData != NULL) { m_pData = NULL; }
		m_pData = D;
		if(m_bTreeHMM) { FOR(i,(int)m_vpAssociatedModels.size()) { cout << "\nMultiple data sets with TreeHMMs not ready yet..."; exit(-1); m_vpAssociatedModels[i]->MakeNewData(m_pData,Overwrite); } }
	} else { // Physically write over the data
		assert(m_pData != NULL && D != NULL);
		*m_pData = *D;
		if(m_bTreeHMM) { FOR(i,(int)m_vpAssociatedModels.size()) { cout << "\nMultiple data sets with TreeHMMs not ready yet..."; exit(-1); m_vpAssociatedModels[i]->MakeNewTree(m_pTree,Overwrite); } }
	}
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->ApplyNewData(m_pData,Overwrite); }

}

////////////////////////////////////////////////////////////////
// Functions for handling process probabilities and other stuff
////////////////////////////////////////////////////////////////

vector <double> CBaseModel::Rates()	{
	int i;
	vector <double> Rat;
	FOR(i,(int)m_vpProc.size())	{ Rat.push_back(m_vpProc[i]->Rate()); }
	return Rat;
}

// Returns the overall rate of the model
double CBaseModel::ModelRate(double R)	{
	int i;
	double Rate = 0.0;
	vector <double> Probs;
//	cout << "\nProcessRatesBefore: " << Rates() << flush;
	FOR(i,(int) m_vpProc.size()) { if(m_vpProc[i]->MaxRate()) { Probs.push_back(0.0); } else { Probs.push_back(m_vpProc[i]->Prob()); } }
	Probs = NormaliseVector(Probs);
	FOR(i,(int)m_vpProc.size()) { if(m_vpProc[i]->MaxRate()) { continue; }  Rate += m_vpProc[i]->Rate() * Probs[i]; }
	if(R >= DBL_EPSILON)	{	// If rqd then reset the rate
		m_ModelRate->SetVal(R,true,true);
	}
	R = m_ModelRate->Val();
	FOR(i,(int)m_vpProc.size()) { if(m_vpProc[i]->MaxRate()) { continue; }
//		cout << "\nSetting process[" << i << "]: " << m_vpProc[i]->Rate() <<  " / " <<  Rate << ")  * " <<  R << flush;
		m_vpProc[i]->Rate( (m_vpProc[i]->Rate() / Rate) * R); }
	Rate = R;
//	cout << "\nProcessRatesAfter: " << Rates() << flush;
//	cout << "\nRate: " << R;
	return Rate;


}

// Removes the process probabilities from the model parameters
void CBaseModel::CleanProcProbs()	{
	int i;
	vector <CPar *>::iterator i_P;
	IFOR(i_P,m_vpPar) {
		FOR(i,(int)m_vpProc.size()) {
			if(*i_P == m_vpProc[i]->ProbPar()) { *i_P = NULL; i_P = m_vpPar.erase(i_P); }
}	}	}

// Gets a list of process probabilities
vector <CPar *> CBaseModel::GetProcessProbs()	{
	int i,j;
	vector <CPar *> ProcProbs;
	FOR(i,(int)m_vpProc.size()) {
		FOR(j,(int)ProcProbs.size()) { if(ProcProbs[j] == m_vpProc[i]->ProbPar()) { break; } }
		if(j == (int)ProcProbs.size()) { ProcProbs.push_back(m_vpProc[i]->ProbPar()); }
	}
	return ProcProbs;
}

// Takes the process probs and prepares them for optimisation (if rqd)
void CBaseModel::PrepareProcessProbs(bool OptProbs)	{
	int i;
	double Total;
	vector <CPar *> ProbProcs;
	// Do some initialisation
	CleanProcProbs();
	ProbProcs = GetProcessProbs();
	// If only a single process then do nothing but return
	if(ProbProcs.size() == 1) { ProbProcs[0] = NULL; return; }
	// Normalise them so they add to oneFastBranchOp
	Total = 0.0; FOR(i,(int)ProbProcs.size()) { Total += ProbProcs[i]->Val(); }
	FOR(i,(int)ProbProcs.size()) { ProbProcs[i]->SetVal(ProbProcs[i]->Val() / Total); }
	// Turn them into probabilities and add them to m_vpPar if they are to be optimised
	if(OptProbs)	{
		FOR(i,(int)ProbProcs.size()) { ProbProcs[i]->SetOptimise(true); m_vpPar.push_back(ProbProcs[i]); }
		ProbabilityScale(&ProbProcs,true,true,true);
	}
}

////////////////////////////////////////////////////////////////
// Fast optimisation function for branches
////////////////////////////////////////////////////////////////
// This function will alway require at least 2 passes of the tree
// The first one is to get reasonable sets of branches at a lower tolerance
double CBaseModel::FastBranchOpt(double CurlnL, double tol, bool *Conv, int NoIter, bool CheckPars)	{
	int i,j, Branches = 0;
	double newlnL, BestlnL = 0.0, working_tol;
	assert(CurlnL < 0);
#if FASTBRANCHOPT_DEBUG == 1
	cout << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW ROUND OF FAST BRANCH OPT >>>>>>>>>>>>>>>>>>>>" << flush;
	cout << "\nHave entered with lnL : " << CurlnL << " cf. " << lnL(true);
#endif
	// This is a fairly meaningless piece of code for trapping errors for multiple trees
	if(m_vbDoBranchDer.empty()) { Error("CBaseModel::FastBranchOpt(...) error. The vector m_vbDoBranchDer is empty. Try called GetOptPar(...) first\n\n"); }
	FOR(i,(int)m_vpProc.size())	{ ;
		if(m_vpProc[0]->Tree() != m_vpProc[i]->Tree()) { Error("\nCan only do CBaseModel::FastBranchOpt on a single tree"); }
		else if(i>0) { continue; }
		if(m_vbDoBranchDer[i] == true) { Branches+= m_vpProc[0]->Tree()->NoBra(); }
	}
	FOR(i,(int)m_vpAssociatedModels.size()) {
//		cout << "\nAssociated model[" << i<< "]: Memory addresses: " << m_vpProc[0]->Tree() << " : " << m_vpAssociatedModels[i]->m_vpProc[0]->Tree() << " ;; " << *m_vpProc[0]->Tree() << " : " << *m_vpAssociatedModels[i]->m_vpProc[0]->Tree();
		if(m_vpProc[0]->Tree() != m_vpAssociatedModels[i]->m_vpProc[0]->Tree()) { Error("\nCan only do CBaseModel::FastBranchOpt on a single tree when there are multiple models"); }
	}
	if(CheckPars) {
		// Check branches are actually being optimised!
		FOR(i,(int)m_vpAllOptPar.size())	{ if(m_vpAllOptPar[i]->Name().find("Branch") != string::npos) { break; } }
		// Return if conditions not met
		if(Branches != Tree()->NoBra() || i == (int)m_vpAllOptPar.size()) { return CurlnL; }
	}
	///////////////////////////////////////////////////////////////////
	// 1. Set up the Q matrices
	PreparelnL();
	FOR(j,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[j]->PreparelnL(); }
	// 2. Build all the partial likelihoods and get the processes sitewise likelihood
	FOR(i,(int)m_vpProc.size())	{
		if(m_vbDoBranchDer[i] == true) {
			if(m_vpProc[0]->Tree() != m_vpProc[i]->Tree()) { Error("\nHaven't checked that multiple trees work..."); exit(-1); }
			m_vpProc[i]->PrepareBraDer();
	}	}
	FOR(j,(int)m_vpAssociatedModels.size()) {
		FOR(i,(int)m_vpAssociatedModels[j]->m_vpProc.size())	{
			if(m_vbDoBranchDer[i] == true) {
				m_vpAssociatedModels[j]->m_vpProc[i]->PrepareBraDer();
	}	}	}
	// Do the optimisation
//	cout << "\n------------- Doing branch opt: actual_tol: " << tol <<" ------------";
//	if(tol > 1.0E-3) { working_tol = tol; } else { working_tol = 1.0E-3; }
	working_tol = max(1.0E-2,tol * 100);
	///////////////////////////////////////////////////////////////////
	// Only do cyclical optimisation with multiple branches
	if(Tree()->NoBra() == 1) {
		DoBraOpt(true,0,1,0,true,&CurlnL,tol,false);
		return CurlnL;
	}
	BestlnL = lnL(true);							// 1. Do the first calculation
	FOR(i,NoIter)	{
		newlnL = BestlnL;	// new_lnL hold current optimal likelihood
#if FASTBRANCHOPT_DEBUG == 1
//#if DEVELOPER_BUILD == 1
		cout << "\n\n--- Round " << i<< ". " << newlnL << " (tol: "<< working_tol << ") ---";;
		cout << "\nOriginal branches:  "; int j; FOR(j,Tree()->NoBra()) { cout << Tree()->B(j) << " "; }
		cout << flush;
#endif
		BranchOpt(-1,Tree()->StartCalc(),-1, &BestlnL,working_tol);	// 2. Run the fast optimisation routine
		BestlnL = lnL(true); // Run likelihood to update partial likelihoods and clean it up properly
		if(working_tol > tol) { working_tol = max(tol,working_tol/100); }
//#if DEVELOPER_BUILD == 1
#if FASTBRANCHOPT_DEBUG == 1
		cout << "; " << BestlnL << " == " << lnL() << "; diff = " << BestlnL - lnL();
		cout << "\nTree: " << *Tree();
		if(fabs(BestlnL - lnL()) > tol) { cout << "\nBig Error..."; exit(-1); }
#endif
		if(fabs(BestlnL - newlnL) < tol) { break; }				// 3. Control exit
	}
	if(Conv != NULL) { if(i==NoIter) { *Conv = false; } else { *Conv = true; } }
//	cout << "\nReturning: " << BestlnL << " cf. " << lnL() << " fabs: " << fabs(BestlnL - lnL()); // exit(-1);
	assert(BestlnL < 0);
	return BestlnL;
}

// Optimise the set of branches
void CBaseModel::BranchOpt(int First,int NTo, int NFr, double *BestlnL,double tol)	{
	int i,j,OriFirst = First;
//	cout << "\nInto CBaseModel::BranchOpt(" << First << ", " << NTo << "," << NFr << ", " << *BestlnL << ")";
	if(NFr == -1)	{
		rFOR(i,Tree()->NoLinks(NTo)) { if(Tree()->NodeLink(NTo,i) == -1) { continue; } BranchOpt(First,Tree()->NodeLink(NTo,i),NTo,BestlnL,tol);
	}	} else {
		// Always perform the calculations in the first place
		if(Tree()->NodeType(NFr) == leaf)	{ // Do the leaf calculations for first calc
#if ALLOW_BRANCH_OPTIMISE == 1
			DoBraOpt(First,NTo,NFr,Tree()->NodeBra(NFr,0),true,BestlnL,tol);
#else
			DoBralnL(Tree()->NodeBra(NFr,0),NFr,NTo);
			PreparePT(Tree()->NodeBra(NFr,0));
			FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true); }
			Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true);

#endif
		} else if(Tree()->NodeType(NTo) == leaf)	{ // Do the leaf calculations for other calcs
#if ALLOW_BRANCH_OPTIMISE == 1
			DoBraOpt(First,NFr,NTo,Tree()->NodeBra(NTo,0),true,BestlnL,tol);
#else
			DoBralnL(Tree()->NodeBra(NTo,0),NTo,NFr);
			PreparePT(Tree()->NodeBra(NTo,0));
			FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NFr,NTo,Tree()->NodeBra(NTo,0),Tree(),First,true); }
			Leaf_update(NFr,NTo,Tree()->NodeBra(NTo,0),Tree(),First,true);
#endif
		} else { // Do the internal calculations
			FOR(i,Tree()->NoLinks(NTo))	{ if(Tree()->NodeLink(NTo,i) == NFr || Tree()->NodeLink(NTo,i) == -1) { break; } }
			assert(i != Tree()->NoLinks(NTo));
			// If the node from isn't a leaf node do the internal calculation (i.e. avoids first node)
			if(Tree()->NodeType(Tree()->NodeLink(NTo,i)) != leaf)	{
#if ALLOW_BRANCH_OPTIMISE == 1
				DoBraOpt(First,NTo,NFr, Tree()->NodeBra(NTo,i),false,BestlnL,tol);
#else
				DoBralnL(Tree()->NodeBra(NTo,i),NFr,NTo);
				PreparePT(Tree()->NodeBra(NTo,i));
				FOR(j,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[j]->Bran_update(NTo,NFr,Tree()->NodeBra(NTo,j),Tree(),First,true,false); }
				Bran_update(NTo,NFr,Tree()->NodeBra(NTo,i),Tree(),First,true,false);
#endif
		}	}
		// Do the looping
		First = 0;
		FOR(i,Tree()->NoLinks(NTo))	{
			if(Tree()->NodeLink(NTo,i) == NFr || Tree()->NodeLink(NTo,i) == -1) { continue; }
			BranchOpt(First,Tree()->NodeLink(NTo,i),NTo,BestlnL,tol);
			First = 1;
	}	}
	if(NFr != -1) { if(Tree()->NodeType(NTo) != leaf && Tree()->NodeType(NFr) != leaf) {
		PreparePT(Tree()->FindBra(NTo,NFr));
		Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
		FOR(j,(int) m_vpAssociatedModels.size()) {
			m_vpAssociatedModels[j]->PreparePT(Tree()->FindBra(NTo,NFr));
			m_vpAssociatedModels[j]->Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
}	}	}	}

// ------------------ Single branch optimisation function, including calculation of partial likelihoods --------------------
void CBaseModel::SingleBranchOpt(int Br, double *BestlnL, double tol) {
	int i;
	bool Leaf = false;
	// Check input and set up partial likelihoods
	assert(Tree()->BraLink(Br,0) != -1 && Tree()->BraLink(Br,1) != -1);
//	cout << "\nBranch[" << Br<<"] = (" << Tree()->BraLink(Br,0) << "," << Tree()->BraLink(Br,1) << "): " << Tree()->B(Br) << " == " << lnL(true);
	FOR(i,(int)m_vpProc.size())	{
			if(Tree()->NodeType(Tree()->BraLink(Br,0)) != leaf) { m_vpProc[i]->PreparePartialL(Tree()->BraLink(Br,0),Tree()->BraLink(Br,1)); } else { Leaf = true; }
			if(Tree()->NodeType(Tree()->BraLink(Br,1)) != leaf) { m_vpProc[i]->PreparePartialL(Tree()->BraLink(Br,1),Tree()->BraLink(Br,0)); } else { Leaf = true; }
	}
//	cout << " from " << DoBralnL(Br,Tree()->BraLink(Br,0),Tree()->BraLink(Br,1));
	DoBraOpt(0,Tree()->BraLink(Br,0),Tree()->BraLink(Br,1),Br,Leaf,BestlnL,tol,false);
//	cout << " -> " << DoBralnL(Br,Tree()->BraLink(Br,0),Tree()->BraLink(Br,1));
	*BestlnL = lnL();
}


// Optimise a single branch providing the partial likelihoods are correctly assigned
#define RETURN_DOBRAOPT Par = NULL; PreparePT(Br); \
	/* cout << "\n\tReturning Branch: " << Tree()->B(Br) << "; best_lnL: " << *BestlnL << " (x1: " << x1_lnL << "; x2: " << x2_lnL << "; x3: " << x3_lnL << ")"; */ \
	if(AllowUpdate) { if(IsExtBra) { FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NTo,NFr,Br,Tree(),First,true); } Leaf_update(NTo,NFr,Br,Tree(),First,true); } \
	else { FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Bran_update(NTo,NFr,Br,Tree(),First,true,false); } Bran_update(NTo,NFr,Br,Tree(),First,true,false); }}  return;
#define GS_DELTA 5
#define DO_BRENT 0		// Whether to do Brent's algorithm in fast search

const double phi = (1 + sqrt(5)) /2;
const double resphi = 2 - phi;

#if FASTBRANCHOPT_DEBUG == 1
bool ErrorCheckInBralnL = false;
#endif

// Brent version of DoBraOpt
void CBaseModel::DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra, double *BestlnL,double tol, bool AllowUpdate)	{
	int i;
	bool HaveBound;
	double *p_x = Tree()->OptimiserB(Br),x1 = -BIG_NUMBER,x2,x3 = -BIG_NUMBER,x1_lnL = 1.0,x2_lnL = -fabs(*BestlnL),x3_lnL = 1.0, xi,temp, ori_x, ori_lnL;
	tol = max(tol,FULL_LIK_ACC);
	double dx = max(DX,tol);	// For looser optimisation have slightly larger bounds
	CPar *Par  = Tree()->pBra(Br);
	ori_x = x2 = *p_x;	ori_lnL = *BestlnL;
//	ori_lnL = *BestlnL = DoBralnL(Br,NFr,NTo);
#if FASTBRANCHOPT_DEBUG == 1
	cout << "\nBranch["<<Br<<"]=" << *p_x << " has DoBralnL: "<< DoBralnL(Br,NFr,NTo) << " cf. " << DoBralnL(Br,NFr,NTo) << " and ori_lnL: " << ori_lnL;
	if(fabs(*BestlnL - DoBralnL(Br,NFr,NTo)) > 1.0E-6) {
		cout.precision(16);
		cout << "\n ===================== ERROR IN BRANCH " << Br << " (" << Tree()->BraLink(Br,0) << "," << Tree()->BraLink(Br,1) << ") ================";
		cout << "\nTrying DoBralnL again = " << DoBralnL(Br,NFr,NTo);
//		cout << "\nOtherway round? = " << DoBralnL(Br,NTo,NFr);
		ErrorCheckInBralnL = true;	// Output from DoBralnL
		double Counter1 = DoBralnL(Br,NFr,NTo), Counter2 = lnL(true);
		double Pouncer1 = 0.0, Pouncer2 = 0.0;
		cout << "\nError in Update... BestlnL:  " << *BestlnL << " != " << Counter1 << " -> diff = " << fabs(*BestlnL - Counter1);
		cout << "\nReal likeklihood: " << Counter2 << endl;
		cout << "\nPouncer1: " << Pouncer1 << " cf. " << Pouncer2;
//		cout << "\nSitewise likelihoods under full: ";
//		FOR(i,m_pData->m_iSize) {  i= 215; cout << "\nSite[" << i<< "]: " << m_vpProc[0]->L(i).LogP() << " * " << m_pData->m_ariPatOcc[i] << " = " <<  m_vpProc[0]->L(i).LogP() * m_pData->m_ariPatOcc[i];  break; }
//		cout << "\nAnd the model: " << *this;
		exit(-1);
	}

#endif

	// ------------------------------------- Catch entry into bounds ------------------------------------
	if(!Par->CheckLowBound()) {	// Check lower bound
		x1 = x2; x1_lnL = *BestlnL;
		x2 = *p_x = Par->LowBound() + (DX * 2);
		x2_lnL = DoBralnL(Br,NFr,NTo); x2 = *p_x; /* catches bound corrections */ m_iFastBralnL_Bracket++;

//		cout << "\nCheck lower bound failed: x1=" << x1 << " : " << x1_lnL << " cf. x2=" << x2 << " : " << x2_lnL;
		if(x1_lnL + tol > x2_lnL) { // At genuine low bound
			*p_x = x1; *BestlnL = x1_lnL;
			Par->StoreOptBounds(Par->LowBound(),Par->LowBound() + DX);
//			cout << "\n\t\tReturning on low bound: " << *p_x << " = " << *BestlnL << " == " << DoBralnL(Br,NFr,NTo);
			RETURN_DOBRAOPT;
	}	}
	if(!Par->CheckUpBound()) {	// Check upper bound
		x3 = x2; x3_lnL = *BestlnL;
		x2 = *p_x = Par->UpBound() - (DX * 2);
		x2_lnL = DoBralnL(Br,NFr,NTo); x2 = *p_x; /* catches bound corrections */ m_iFastBralnL_Bracket++;
		if(x3_lnL + tol > x2_lnL) { // At genuine up bound
			*p_x = x3; *BestlnL = x3_lnL;
			Par->StoreOptBounds(Par->UpBound() - DX, Par->UpBound());
			RETURN_DOBRAOPT;
	}	}
	// ----------------------------------- Bracketing routine -------------------------------------
	// Initialise sensible bounds; this is necessary because after other optimisation the original value may fall out of bound
	Par->StoreOptBounds(min(x2-DX,Par->OptLow()),max(x2+DX,Par->OptUp()));
	// Get left bracketing
//	cout << "\nDoing left";
	HaveBound = true;
	if(fabs(Par->LowBound() - Par->OptLow()) < FLT_EPSILON) { HaveBound = false; dx = DX; } else { dx = 0.0; }	 // Set original dx
	while(x2_lnL < x1_lnL)	{
		if(HaveBound) {
			*p_x = x1 = max(Par->LowBound(),Par->OptLow() - (dx * 10));
		} else {
			*p_x = x1 = max(Par->LowBound(),x2 - (dx * 100));
		}

		x1_lnL = DoBralnL(Br,NFr,NTo); x1 = *p_x; /* catches bound corrections */ m_iFastBralnL_Bracket++;
//		cout << "\n\Left: (" << x1 << ": " << x1_lnL << "," << x2 << ": " << x2_lnL << "," << x3 << ": " << x3_lnL << ")";
		// If at low bound organise around it and check for exit conditions
		if(!Par->CheckLowBound()) {
			if(x2_lnL > x1_lnL) { break; } // Have curvature around x2. All okay
			// Otherwise need to try and check what's going on with the bounded parameter
			x3 = x2; x3_lnL = x2_lnL;
			x2 = *p_x = (x1+x3)/2;
			x2_lnL = DoBralnL(Br,NFr,NTo); m_iFastBralnL_Bracket++;
			if(x1_lnL + tol > x2_lnL ) { // True maxima at bound
				*p_x = x1; *BestlnL = x1_lnL;
				Par->StoreOptBounds(Par->LowBound(),Par->LowBound() + DX);
				RETURN_DOBRAOPT;
			}
			break;
		}
		if(x1_lnL > x2_lnL)	{	// When moving left improves likelihood, then keep cycling
			x3 = x2; x3_lnL = x2_lnL;
			x2 = x1; x2_lnL = x1_lnL;
			x1_lnL = 1.0;
		}
		if(fabs(dx) < FLT_EPSILON) { dx = DX; } else { dx *= GS_DELTA; }
	}
//	cout << "\n\tafter left: (" << x1 << ": " << x1_lnL << "," << x2 << ": " << x2_lnL << "," << x3 << ": " << x3_lnL << ")";
//	cout << "\nDoing right";
	// Get right bracketing
	HaveBound = true;
	if(fabs(Par->UpBound() - Par->OptUp()) < FLT_EPSILON) { dx = DX; HaveBound = false; } else { dx = 0.0; }	// Set original dx
	if(x3_lnL > 0.0)	{		// Search for right bound if it's needed
		while(x2_lnL < x3_lnL)	{
			if(HaveBound) {
				*p_x = x3 = max(Par->LowBound(),Par->OptUp() + (dx * 10));
			} else {
				*p_x = x3 = max(Par->LowBound(),x2 + (dx * 100));
			}
			x3_lnL = DoBralnL(Br,NFr,NTo); x3 = *p_x; /* catches bound corrections */ m_iFastBralnL_Bracket++;
//			cout << "\n\tRight: (" << x1 << ": " << x1_lnL << "," << x2 << ": " << x2_lnL << "," << x3 << ": " << x3_lnL << ")";
			if(!Par->CheckBound())	{
//				cout << "\nCheck Up bound break";
				if(x2_lnL > x1_lnL) { break; } // Have curvature around x2. All okay
				// Otherwise need to try and check what's going on with the bounded parameter
				x1 = x2; x1_lnL = x2_lnL;
				x2 = *p_x = (x1 + x3)/2;
				x2_lnL = DoBralnL(Br,NFr,NTo); m_iFastBralnL_Bracket++;
				if(x3_lnL + tol > x2_lnL) { // True maxima at bound
					*p_x = x3; *BestlnL = x3_lnL;
					Par->StoreOptBounds(Par->UpBound() - DX, Par->UpBound());
					RETURN_DOBRAOPT;
				}
				break;
			}
			if(x3_lnL > x2_lnL)	{
				x1 = x2; x1_lnL = x2_lnL;
				x2 = x3; x2_lnL = x3_lnL;
				x3_lnL = 1.0;
			}
			if(fabs(dx) < FLT_EPSILON) { dx = DX; } else { dx *= GS_DELTA; }
	}	}
//	cout << "\n\tafter Right: (" << x1 << ": " << x1_lnL << "," << x2 << ": " << x2_lnL << "," << x3 << ": " << x3_lnL << ")";
	if((x1_lnL > x2_lnL || x3_lnL > x2_lnL) || (x1 > x2 || x2 > x3)) {
		cout.precision(10);
		cout << "\nBounding failed in DoBraOpt(...)... ";
		cout << "\n\tafter Right: (" << x1 << ": " << x1_lnL << "," << x2 << ": " << x2_lnL << "," << x3 << ": " << x3_lnL << ")";
		cout << "\nReturning with no improvement and please report to Simon Whelan"; exit(-1);
		if(x1_lnL > x2_lnL) { *p_x = x1; *BestlnL = x1_lnL; } else { *p_x = x3; *BestlnL = x3_lnL; }
	}
	// Some output for debugging if needed
//	if(fabs(x1-x2) <  1.0E-9) { exit(-1); }
//	if(x1 < 0 || x2 < 0 || x3 < 0 || *p_x < 0) { cout << "\nHave detected an error for DoBraOpt: x1: " << x1 << " x2: " <<  x2 << " x3: " <<  x3 <<" x: " <<  *p_x << "!!!"; exit(-1); }
//	cout << "\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
//	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
//	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
#if DO_BRENT == 1
	// ------------------------------------ Brent's algorithm (adapted from Numerical Recipes in C) -----------------------
	assert(x2_lnL > x1_lnL && x2_lnL > x3_lnL);
	double d=0.0, e=0.0, etemp, p, q, r,x,w,v,fx,fw,fv, tol1=1.0E-5;
	x=w=v=x2; fx=fw=fv=-x2_lnL;
	FOR(i,20) {
/*		cout << "\nRound " << i;
		cout << "\n\tx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
		cout << "\n\tx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
		cout << "\n\tx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
		cout << "\n\tx: " << x << " == " << fx;
		cout << "\n\tw: " << w << " == " << fw;
		cout << "\n\tv: " << v << " == " << fv;
*/		// Test exit condition
//		cout << "\nfv - fx: " << fabs(fv-fx) << " cf. tol: " << tol << " == " << fabs(fv - fx) << " ; e: " << e;
		if(fabs(fv - fx) < tol && i > 1) { break; }
		// Construct trial parabolic fit
		xi = 0.5 * (x1 + x2);
		if(fabs(e) > 1.0E-5) {
			r=(x-w)*(fx-fv); q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r; q=2.0*(q-r);
			if(q>0.0) { p = -p; }
			q= fabs(q); etemp = e; e=d;
			// Determine whether parabola fits adequately
			if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(x1-x) || p >= q*(x2-x)) 	{ cout << " -- GS"; d=resphi*(e-(x>=xi ? x1-x : x2-x)); }
			else /* take the parabolic step */ 									{ cout << " -- para"; d=p/q; *p_x=x+d; if(*p_x-x1 < tol1 || x2-*p_x < tol1) { d = SIGN(tol1,xi-x); }}
		} else {
			d=resphi*(e=(x>=xi ? x1-x : x2-x));
		}
		// Do the function evaluation
		*p_x = (fabs(d) >= tol1 ? x+d : x + SIGN(tol1,d));
		temp = -DoBralnL(Br,NFr,NTo);
//		cout << "\n\tNewvalue -- x: " << *p_x << " = " << temp;
		// Now decide what to do with the function evaluation
		if(temp <= fx) {
			if(*p_x >= x) { x1 = x; } else { x2 = *p_x; }
			SHFT(v,w,x,*p_x);
			SHFT(fv,fw,fx,temp);
		} else {
			if(*p_x < x) { x1 = *p_x; } else { x2 = *p_x; }
			if(temp <= fw || w == x) { v=w; w=*p_x; fv=fw; fw=temp; }
			else if (temp <= fx || v == x || v ==w) { v=*p_x; fv = temp; }

		}
	}
	if(fx < fw && fx < fv) 	{ *p_x = x; *BestlnL=-fx; }
	else if(fw < fx && fw < fv)	{ *p_x = w; *BestlnL=-fw; }
	else 						{ *p_x = v; *BestlnL=-fv; }
#else
	// New Goldensection
	assert(x2_lnL + FLT_EPSILON > x1_lnL && x2_lnL + FLT_EPSILON> x3_lnL);
	FOR(i,20) {
		// New value
		if(x3 - x2 > x2 - x1) 	{ *p_x = xi = x2 + resphi * (x3-x2); }
		else					{ *p_x = xi = x2 - resphi * (x2-x1); }
		// break condition
		// if(fabs(x3_lnL - x1_lnL) < tol) { /* cout << "\nBreaking at tol=" << tol << " fabs(" << x3_lnL << " - " << x1_lnL << ")";  */ *p_x = x2 = (x1+x3)/2; break; }
		if(fabs(x3_lnL - x1_lnL) < tol) { /* cout << "\nBreaking at tol=" << tol << " fabs(" << x3_lnL << " - " << x1_lnL << ")";  */ *p_x = x2; break; }
		// Search
		temp = DoBralnL(Br,NFr,NTo);
//		cout.precision(10); cout << "\n[" << i << "] x1: " << x1 << ": " << x1_lnL << "; x2: " << x2 << ": " << x2_lnL << "; x3: " << x3 << ": " << x3_lnL;
//		cout << "\n[i="<<i<<"] xi:" << xi << ": " << temp << " range(" << x1 << "," << x3 << ")=" << x3-x1 << " ; lnL imp: " << fabs(x3_lnL - x1_lnL) << " cf. " << tol;
//		cout << "\n\tx2 [x1=] " << x2_lnL - x1_lnL << " [x3=] " << x3_lnL - x2_lnL;
		if(temp > x2_lnL) {
			if(x3 - x2 > x2 - x1) 	{ x1 = x2; x1_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
			else					{ x3 = x2; x3_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
		} else {
			if(x3 - x2 > x2 - x1)	{ x3 = xi; x3_lnL = temp; }
			else					{ x1 = xi; x1_lnL = temp; }
		}
//		cout << "[x1:" << x1 << ", x2:" << x2 << ", x3:" << x3 << ","<< max(x1_lnL,max(x2_lnL,x3_lnL)) << "]" << flush;
	}
	*p_x = x2;
#endif
/*
	cout.precision(16);
	cout << "\nFinished search: x: " << *p_x << " = " << temp << " == " << DoBralnL(Br,NFr,NTo);
	cout << ": tol= " << max(x2_lnL - x1_lnL,x2_lnL - x3_lnL);
	cout << "\n---\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
*/	// Finish by doing the calculation again to correctly update the partial likelihoods
	m_iFastBralnL_Calls++;
	if(x2_lnL > ori_lnL) { *p_x = x2; *BestlnL = x2_lnL; } else { if(fabs(x2_lnL - ori_lnL) > tol) { cout << "\nWeird... optimiser (tol=" << tol << ") made worse likelihood in CBaseModel::DoBraOpt(...)\n"; cout << " should have: " << ori_lnL << " and have " << DoBralnL(Br,NFr,NTo) << " diff = " << DoBralnL(Br,NFr,NTo) - ori_lnL; } *p_x = ori_x; *BestlnL = ori_lnL;  }
	Par->StoreOptBounds(x1,x3);
	RETURN_DOBRAOPT;
}

// Function that performs basic likelihood calculation
double CBaseModel::DoBralnL(int B, int NF,int NT, bool JustClean)	{
	int i;
	double logL = 0.0;
	static CProb **P = NULL;
	static int CProbSize = -1;
#if DEVELOPER_BUILD == 1
	cout << "\n>>>>>>>>>>>>>>>>>>>> Entering CBaseModel::DoBralnL("<<B<<","<<NF<<","<<NT<<")";
#endif
	// If required just clean up the static allocation
	if(JustClean)	{
		if(P != NULL) {
			FOR(i,CProbSize) { delete P[i]; } delete [] P;
			P = NULL; CProbSize = -1;

		}
		return 0.0;
	}
	// Get likelihoods from associated models
//	cout << "\n>>>>>>>>>>>>>>>>>>>> Entering CBaseModel::DoBralnL("<<B<<","<<NF<<","<<NT<<")";
//	if(m_bMainModel) { cout << "\n>>> QuicklnL"; }
	FOR(i,(int)m_vpAssociatedModels.size())	{
//		cout << "\n\tSublnL["<<i<<"]: " << m_vpAssociatedModels[i]->DoBralnL(B,NF,NT);
		logL += m_vpAssociatedModels[i]->DoBralnL(B,NF,NT);
	}
	// Do the main process
	if(P == NULL) {
		GET_MEM(P,CProb *,m_pData->m_iSize);
		FOR(i,m_pData->m_iSize) { P[i] = new CProb(0.0,0); }
		CProbSize = m_pData->m_iSize;
	} else if(CProbSize != m_pData->m_iSize) {
		FOR(i,CProbSize) { delete P[i]; } delete [] P;
		GET_MEM(P,CProb *,m_pData->m_iSize);
		FOR(i,m_pData->m_iSize) { P[i] = new CProb(0.0,0); }
		CProbSize = m_pData->m_iSize;
	} else { FOR(i,m_pData->m_iSize) { P[i]->Zero(); } }
	// Get the partial likelihoods
	FOR(i,(int)m_vpProc.size())	{
		if(m_vpProc[i]->Prob() < MIN_PROB) { continue; }
		m_vpProc[i]->GetBranchPartL(P,NT,NF,B);
	}
	FOR(i,m_pData->m_iSize) {
#if FASTBRANCHOPT_DEBUG == 1
		if(ErrorCheckInBralnL) {
			if(i == 0) { cout << "\nOutputting DoBralnL Sitewise"; }
			if(i== 215) {
				cout << "\n\t\tSite["<<i<<"]:  ";
				int j; FOR(j,m_pData->m_iNoSeq) {  if(m_pData->m_ariSeq[j][i]< m_pData->m_iChar) { cout << m_pData->m_sABET[m_pData->m_ariSeq[j][i]]; } else { cout << "-"; }}
				cout << " " << P[i]->LogP() << " * " << m_pData->m_ariPatOcc[i] << " = " << P[i]->LogP() * m_pData->m_ariPatOcc[i];
			}
		}
#endif
		logL += P[i]->LogP() * m_pData->m_ariPatOcc[i]; }
//	cout << "\nReturning lnL: " << logL; // << " cf. " << lnL(true);
#if DEVELOPER_BUILD == 1
	cout << "\nDoBralnL(Branch="<<B<<") returning: "<< logL << "\n\n\\\\";
#endif
//	if(m_bMainModel) { cout << "\n\tReturning: " << logL; }
	// Apply extra stuff to the likelihood function if needed
	if(pLikelihood != NULL) {
		logL -= pLikelihood(NULL); // Called as blank. Other arguments intended to allow functionality
	}
//	cout << "\n>>>>>>>>>>>>>>>>>>>> DoBralnL(Branch="<<B<<") returning: "<< logL;
	m_iFastBralnL++;
	return logL;
}

////////////////////////////////////////////////////////////////
// Process handling functions
////////////////////////////////////////////////////////////////

// Adds a gamma distribution to a model
void CBaseModel::MakeGammaModel(int PNum, int NoCat, double InitAlpha)	{
	int i;
	CGammaPar *GP;
	CBaseProcess *NewProc = NULL;
	// Check entry conditions
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->DecompressSpace(); }
	assert(PNum < (int)m_vpProc.size());
	assert(m_vpProc[PNum]->IsGamma() == false); // Check the process hasn't already been made into a gamma model
	// Get the gamma distribution parameter
	GP = new CGammaPar("Alpha",m_vpProc[PNum]->Char(),m_ModelRate, InitAlpha);
	GP->AddRateToGamma(m_vpProc[PNum]);
	// Add parameter to its process
	m_vpProc[PNum]->AddQPar(GP);
	// Add pseudoprocesses for rate
	FOR(i,NoCat-1)	{
		NewProc = m_vpProc[PNum]->GammaRateProcessCopy();
		GP->AddRateToGamma(NewProc);
		m_vpProc.push_back(NewProc);
		NewProc = NULL;
	}
	m_sName = m_sName + "+" + int_to_string(NoCat) + "dG";
	m_vpPar.push_back(GP);
	GP = NULL;
}

// Adds an invariant sites
void CBaseModel::MakeInvariantSitesModel(int PNum)	{
	CBaseProcess *NewProc = NULL;
	m_sName = m_sName + "+I";
	// Prepare the process
	assert(!m_vpProc[PNum]->PseudoProcess() && PNum < (int)m_vpProc.size());
	NewProc = m_vpProc[PNum]->RateProcessCopy();
	NewProc->Rate(0.0); NewProc->ProbPar()->SetVal(0.2);
	NewProc->ProbPar()->Name("P(Inv)");
	// Add the process to the model and sort out the probabilities
	m_vpProc.push_back(NewProc); NewProc = NULL;
	PrepareProcessProbs();
}

// Adds a garbage collector category to the model
// ---
// This is a class whose likelihood is the product of the equilibrium frequencies of the column
//  * This approach is the equivalent to having an infinite rate class
void CBaseModel::MakeGarbageCollectorModel(int PNum) {
	CBaseProcess *NewProc = NULL;
	m_sName = m_sName + "+gC";
	// Prepare the process
	assert(!m_vpProc[PNum]->PseudoProcess() && PNum < (int)m_vpProc.size());
	NewProc = m_vpProc[PNum]->RateProcessCopy();
	NewProc->Rate(BIG_NUMBER); NewProc->MakeMaxRate();
	NewProc->ProbPar()->SetVal(0.1);
	NewProc->ProbPar()->Name("P(Garbage)");
	// Add the process to the model and sort out the probabilities
	m_vpProc.push_back(NewProc); NewProc = NULL;
	PrepareProcessProbs();
}

//////////////////////////////////////////////////////////
// Simulation functions
void CBaseModel::DoSimulation(int NoSeq, int Size,vector <vector <int> > *vviSeqs,bool KeepAnc)	{
	int i,j,SeqNum,Char,DataChar;
	bool DoHidden = false;
	// Check some entry conditions
	if(NoSeq < 1 || Size < 1) { Error("\nCBaseModel::DoSimulation -- Need NoSeq("+int_to_string(NoSeq)+") and Size("+int_to_string(Size)+") to be greater than zero\n\n"); }
//	if(m_pData != NULL) { Error("\nDoing CBaseModel::DoSimulation when data already established...\n\n"); }
	if(m_pTree == NULL) { Error("\nNeed a tree to run CBaseModel::DoSimulation\n\n"); }
	FOR(i,(int)m_vpProc.size()) { if((int) m_vpProc[i]->EqmPar().size() > 1) { Error("\nCan only simulate when there is a single equilibrium..."); } }
	// Some initialisation
	vector <int> Procs(Size,-1);
	vector <double> Probs;
	vviSeqs->clear();
	FOR(i,(2*NoSeq)-2)	{ vviSeqs->push_back(Procs); }
	//////////////////////////////////////////////////
	// Get starting sequence
	SeqNum = Tree()->StartCalc();
	if(!InRange(SeqNum,0,(int)vviSeqs->size())) { Error("\nTrying to initialise CBaseModel::DoSimulation start sequence in odd place...\n\n"); }
	// Initialise Procs to processes that sequence evolve to
	FOR(i,(int)m_vpProc.size())	{ Probs.push_back(m_vpProc[i]->Prob()); }
	if(diff(Sum(&Probs),1.0)) { Error("\nProcess probabilities do not sum to one..."); }
	FOR(i,Size)	{ Procs[i] = DiscreteRand(&Probs); }
	// Doing a sort will make things easier in the long run
	sort(Procs.begin(),Procs.end());
	assert(vviSeqs!=NULL);
	int Counters[20]; FOR(i,20) { Counters[i] = 0; }
	FOR(i,Size)	{
		if(i==0) {
			if(m_vpProc[Procs[0]]->EqmPar().empty()) { Probs.clear(); FOR(j,m_vpProc[Procs[0]]->Char()) { Probs.push_back(1.0); } Probs = NormaliseVector(Probs); }
			else { Probs = m_vpProc[Procs[0]]->EqmPar()[0]->Eqm(); }
			if(m_vpProc[Procs[0]]->HiddenChar() > 1) { DoHidden = true; }
		} else if(Procs[i-1] != Procs[i]) {
			if(m_vpProc[Procs[i]]->EqmPar().empty()) { Probs.clear(); FOR(j,m_vpProc[Procs[i]]->Char()) { Probs.push_back(1.0); } Probs = NormaliseVector(Probs); }
			else { Probs = m_vpProc[Procs[i]]->EqmPar()[0]->Eqm(); }
			if(m_vpProc[Procs[i]]->HiddenChar() > 1) { DoHidden = true; }
		}
		vviSeqs->at(SeqNum)[i] = DiscreteRand(&Probs);
//		Counters[vviSeqs->at(SeqNum)[i]]++;
	}
//	FOR(i,20) { cout << "\nCounter["<<i<<"]\t" << Counters[i] << "\t" << (double)Counters[i] / (double)Size << "\t" << fabs(((double)Counters[i] / (double)Size)  - m_vpProc[Procs[0]]->EqmPar()[0]->Eqm()[i]); } cout << "\n\n";
	// Now evolve the rest...

	PreparelnL(true);
	FOR(i,Tree()->NoLinks(SeqNum)) { EvolveSequences(SeqNum,Tree()->NodeLink(SeqNum,i),vviSeqs,&Procs); }
	// Now fix the sequences into the observable types (for THMMs...)
	if(DoHidden) {
		FOR(i,Size)	{
			DataChar = m_vpProc[Procs[i]]->DataChar();
			Char=m_vpProc[Procs[i]]->Char();
			if(Char % DataChar != 0) { 	Error("\nTrying to compress states, but they don't divide evenly...\n"); }
			FOR(j,NoSeq) {
				if(vviSeqs->at(j)[i] == Char)	{ vviSeqs->at(j)[i] = DataChar; }
				else							{ vviSeqs->at(j)[i] = (vviSeqs->at(j)[i] % DataChar); }
				if(vviSeqs->at(j)[i] > DataChar) { cout << " -> " << vviSeqs->at(j)[i]; Error("This shouldn't happen..."); }
	}	}	}
//	FOR(i,20) { Counters[i] = 0; }
//	FOR(j,Tree()->NoSeq()) { FOR(i,Size) { Counters[vviSeqs->at(j)[i]]++; } }
//	FOR(i,20) { cout << "\n\tCounter["<<i<<"]\t" << Counters[i] << "\t" << (double)Counters[i] / ((double)Size * Tree()->NoSeq())<< "\t" << (((double)Counters[i] / ((double)Size * Tree()->NoSeq()))  - m_vpProc[Procs[0]]->EqmPar()[0]->Eqm()[i]); } cout << "\n";
//	exit(-1);
}
// Do in-order tree-traversal evolving sequences as it goes
void CBaseModel::EvolveSequences(int NodeFrom, int NodeTo, vector <vector <int> > *Seqs, vector <int> *Procs)	{
	int i,j,Branch,Char;
	vector <double> Probs;
	vector <double> PT;
	// Check entry
	if(!InRange(NodeFrom,0,Tree()->NoNode()) || !InRange(NodeTo,0,Tree()->NoNode())) { Error("\nError in CBaseModel::EvolveSequences, NodeFrom=" + int_to_string(NodeFrom) + " or NodeTo=" + int_to_string(NodeTo) + " out of range...\n\n"); }
	if(Seqs == NULL || Procs == NULL) { Error("\nPassed NULL to CBaseModel::EvolveSequences()...\n\n"); }
	if((int)Seqs->size() != Tree()->NoNode()) { Error("\nSpace passed to CBaseModel::EvolveSequences of incorrect size...\n\n"); }
	if(Seqs->at(NodeFrom).size() != Seqs->at(NodeTo).size() || Seqs->at(NodeTo).size() != Procs->size()) { Error("\nUnmatching space for sequences in CBaseModel::EvolveSequences...\n\n"); }
	if(Seqs->at(NodeTo)[0] != -1) { Error("\nIn CBaseModel::EvolveSequences, the NodeTo sequence seems to have already been assigned...\n\n"); }
	// Some more initialisation
	Branch = Tree()->FindBra(NodeFrom,NodeTo);
	// Do the evolution
	FOR(i,(int)Procs->size())	{
		// Create the appropriate PT (if necessary)
		if(i==0) { Char = m_vpProc[Procs->at(i)]->Char(); Probs.clear(); Probs.assign(Char,-1); m_vpProc[Procs->at(0)]->CreatePTMats(Branch); PT = m_vpProc[Procs->at(0)]->GetPT(Branch); }
		else if(Procs->at(i-1) != Procs->at(i)) { Char = m_vpProc[Procs->at(i)]->Char(); Probs.clear(); Probs.assign(Char,-1); m_vpProc[Procs->at(i)]->CreatePTMats(Branch); PT = m_vpProc[Procs->at(i)]->GetPT(Branch); }
		FOR(j,Char) { Probs[j] = PT[(Seqs->at(NodeFrom)[i] * Char) + j]; }
		Seqs->at(NodeTo)[i] = DiscreteRand(&Probs);
	}
	// Do the traversal
	if(NodeTo < Tree()->NoSeq()) { return; }
	FOR(i,Tree()->NoLinks(NodeTo)) {
		if(Tree()->NodeLink(NodeTo,i) == NodeFrom) { continue; }
		EvolveSequences(NodeTo,Tree()->NodeLink(NodeTo,i),Seqs,Procs);
	}
}

////////////////////////////////////////////////////////////////
// Partial likelihood functions -- used for robust counting

// Core function that sets up recursion and check input values
void CBaseModel::PartialLOutput(string File, int DoBranch) {
	int i;
	if(!m_vpAssociatedModels.empty()) { Error("Error: PartialLOutput not set up for models with multiple components (!m_vpAssociatedModels.empty())"); }
	ofstream out(File.c_str());
	out.precision(16);
	out << "OriginalLikelihood\t" << lnL(true);
	Tree()->CreateBranchLabels(); Tree()->OutLabel(); Tree()->OutBra(); Tree()->OutName();
	out << "\nTree:\t" << *m_pTree;
	out << "\nMapping to real data:";
	FOR(i,m_pData->m_iTrueSize) { out << "\t" << m_pData->m_ariPatMap[i]; }
	FOR(i,(int)m_vpProc.size()) {
		if(m_vpProc[i]->PseudoProcess()) { continue; }
		out << "\nProcess["<<i<<"] QMat:";
		m_vpProc[i]->OutQ(-1,out);
		out << "\nEqm:\t" << m_vpProc[i]->Eqm(0);
	}
	out << "\n//";
	RecPartLOut(-1,Tree()->StartCalc(),-1,out,DoBranch);
	out.close();
}

// Recursive function doing output
void CBaseModel::RecPartLOut(int First, int NTo, int NFr, ostream &output, int Branch) {
	int i,j,OriFirst = First;
	if(NFr == -1)	{
		rFOR(i,Tree()->NoLinks(NTo)) { if(Tree()->NodeLink(NTo,i) == -1) { continue; } RecPartLOut(First,Tree()->NodeLink(NTo,i),NTo,output,Branch);
	}	} else {
		// Always perform the calculations in the first place
		if(Tree()->NodeType(NFr) == leaf)	{ // Do the leaf calculations for first calc
			DoPartLOut(NTo, NFr, Tree()->NodeBra(NFr,0),output,Branch);
			PreparePT(Tree()->NodeBra(NFr,0));
			FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true); }
			Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true);
		} else if(Tree()->NodeType(NTo) == leaf)	{ // Do the leaf calculations for other calcs
			DoPartLOut(NTo, NFr, Tree()->NodeBra(NTo,0),output,Branch);
			PreparePT(Tree()->NodeBra(NTo,0));
			FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NFr,NTo,Tree()->NodeBra(NTo,0),Tree(),First,true); }
			Leaf_update(NFr,NTo,Tree()->NodeBra(NTo,0),Tree(),First,true);
		} else { // Do the internal calculations
			FOR(i,Tree()->NoLinks(NTo))	{ if(Tree()->NodeLink(NTo,i) == NFr || Tree()->NodeLink(NTo,i) == -1) { break; } }
			assert(i != Tree()->NoLinks(NTo));
			// If the node from isn't a leaf node do the internal calculation (i.e. avoids first node)
			if(Tree()->NodeType(Tree()->NodeLink(NTo,i)) != leaf)	{
				DoPartLOut(NTo, NFr, Tree()->NodeBra(NTo,i),output,Branch);
				PreparePT(Tree()->NodeBra(NTo,i));
				FOR(j,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[j]->Bran_update(NTo,NFr,Tree()->NodeBra(NTo,j),Tree(),First,true,false); }
				Bran_update(NTo,NFr,Tree()->NodeBra(NTo,i),Tree(),First,true,false);
		}	}
		// Do the looping
		First = 0;
		FOR(i,Tree()->NoLinks(NTo))	{
			if(Tree()->NodeLink(NTo,i) == NFr || Tree()->NodeLink(NTo,i) == -1) { continue; }
			RecPartLOut(First,Tree()->NodeLink(NTo,i),NTo,output,Branch);
			First = 1;
	}	}
	if(NFr != -1) {
		if(Tree()->NodeType(NTo) != leaf && Tree()->NodeType(NFr) != leaf) {
			PreparePT(Tree()->FindBra(NTo,NFr));
			Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
			FOR(j,(int) m_vpAssociatedModels.size()) {
				m_vpAssociatedModels[j]->PreparePT(Tree()->FindBra(NTo,NFr));
				m_vpAssociatedModels[j]->Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
	}	}	}
}

#define DO_PARTLCONFIRM 1

void CBaseModel::DoPartLOut(int NTo, int NFr, int Br, ostream &out,int Branch) {
	vector <int> L,R;
	int i,j;
	Tree()->BranchSets(Br,&L,&R);
	out << "\nBranch["<<Br<<"] (" << NTo << "," << NFr <<")";
	out << "\nSplitLeft\t" << L;
	out << "\nSplitRight\t" << R << flush;
	out << "\nNamedSplitLeft\t"; FOR(i,L.size()) { out << m_pData->m_vsName[L[i]] << " "; }
    out << "\nNamedSplitRight\t"; FOR(i,R.size()) { out << m_pData->m_vsName[R[i]] << " "; }
	PreparePT(Br);
	OutPT(out,-1,Br);
	out << "\n---";
	// Calculate test likelihood
#if DO_PARTLCONFIRM == 1
	double Likelihood = 0, sublnL, temp,Vec1[500],Vec2[500], MyPT[25000];
	int k,Char, Scale;
	CProb SiteP, tempP;
	vector <double> store, eqm;
//	cout << "\n--- Branch["<< Br << "] (" << NFr << "," << NTo << ")";
	FOR(i,m_pData->m_iSize) {
		sublnL = 0; Scale = 0;
		SiteP.Zero();
//		if(i<5 && Br == 27) { cout << "\n["<<i<<"]: " << m_vpProc[0]->GetPartL(i,NTo);cout << "\n["<<i<<"]: " << m_vpProc[0]->GetPartL(i,NFr); }
		FOR(j,(int)m_vpProc.size()) {
			Char = m_vpProc[j]->Char();
			store = m_vpProc[j]->GetPT(Br);
			FOR(k,(int)store.size()) { MyPT[k] = store[k]; }
			eqm = m_vpProc[j]->RootEqm();
			store = m_vpProc[j]->GetPartL(i,NTo);
			FOR(k,(int)store.size()) { Vec1[k] = store[k]; }
			VMat(Vec1,MyPT,Vec2,Char);
			store = m_vpProc[j]->GetPartL(i,NFr);
//			if(i<5) { cout << "\nSite["<<i<<"] Proc["<<j<<"]: " << store[0] << " * " << Vec2[0] << " * " << eqm[0]; }
			temp = 0; FOR(k,Char) { temp += store[k] * Vec2[k] * eqm[k]; }
//			if(i<5) { cout << " == " << temp << " * " << m_vpProc[j]->Prob() << " -> " << temp * m_vpProc[j]->Prob(); }

			if(j == 0) {
				SiteP.Assign(temp * m_vpProc[j]->Prob(),m_vpProc[j]->GetPartLScale(NFr,i) + m_vpProc[j]->GetPartLScale(NTo,i));
			} else {
				tempP.Zero();
				tempP.Assign(temp * m_vpProc[j]->Prob(),m_vpProc[j]->GetPartLScale(NFr,i) + m_vpProc[j]->GetPartLScale(NTo,i));
				SiteP.Add(tempP,true);
			}
		}
//		if(i<5) { cout << " == " << log(sublnL)<< " cf. " << m_arL[i].LogP(); }
		Likelihood += SiteP.LogP() * m_pData->m_ariPatOcc[i];

	}
	cout << " == " << Likelihood;

#endif
	// Get the partial likelihoods
	FOR(i,m_pData->m_iSize) {
//		if(i>3) { break; }
		out << "\nSite["<<i<<"]\t" << m_arL[i].LogP()  << "\t" << m_pData->m_ariPatOcc[i];
		FOR(j,(int)m_vpProc.size()) {
			out << "\n\tProcess["<<j<<"] Left";
			m_vpProc[j]->OutPartL(i,NTo,out);
			out << "\n\tProcess["<<j<<"] Right";
			m_vpProc[j]->OutPartL(i,NFr,out);
		}
	}
	out << "\n//";
}

////////////////////////////////////////////////////////////////
// Tree HMM stuff
void CBaseModel::MakeTreeHMM()	{
	cout << "\nMaking tree HMM...";
	int i;
	CData *MatchData, *DeleteData;
	CBaseModel *MatchModel, *DeleteModel;
	CPar *RatePar;
	RatePar = new CPar("IndelRate",0.75,true,1.0E-5);
	// Check object is correct
	assert(m_pData != NULL);
	assert(m_bTreeHMM == false);
	// Add TreeHMM to name
	m_sName = m_sName + "+TreeHMM";
	// Make the data
	MatchData = m_pData->MakeMatchData();
	DeleteData = m_pData->MakeDeleteData();
	// Make the Match model
	MatchModel = new CRY(MatchData,m_pTree); MatchModel->m_bMainModel = false;
	MatchModel->m_sName = "MatchModel" ;
	delete MatchModel->m_ModelRate; MatchModel->m_ModelRate = RatePar; m_vpPar.push_back(MatchModel->m_ModelRate);
//	MatchModel->MakeGammaModel(0,4);
	FOR(i,(int) MatchModel->m_vpPar.size()) {
		MatchModel->m_vpPar[i]->Name("Match(" + MatchModel->m_vpPar[i]->Name() + ")");
		m_vpPar.push_back(MatchModel->m_vpPar[i]);
	}
	// Make the Delete model
	DeleteModel = new CRY(DeleteData,m_pTree); DeleteModel->m_bMainModel = false;
	DeleteModel->m_sName = "DeleteModel";
	DeleteModel->m_ModelRate; DeleteModel->m_ModelRate = RatePar;
//	DeleteModel->MakeGammaModel(0,4);
	FOR(i,(int) DeleteModel->m_vpPar.size()) {
		DeleteModel->m_vpPar[i]->Name("Delete(" + DeleteModel->m_vpPar[i]->Name() + ")");
		m_vpPar.push_back(DeleteModel->m_vpPar[i]);
	}
	// Clean up the Rate parameter
	RatePar = NULL;
/*
	cout << "\nLikelihood Subs: " <<	lnL();
	cout << "\nLikelihood Matches: " << MatchModel->lnL();
	cout << "\nLikelihood Deletes: " << DeleteModel->lnL();
	cout << "\nTotal: " << lnL() + MatchModel->lnL() + DeleteModel->lnL();
*/
	// Add models to object
	m_bTreeHMM = true;
	m_vpAssociatedModels.push_back(MatchModel); MatchModel = NULL; MatchData = NULL;
	m_vpAssociatedModels.push_back(DeleteModel); MatchModel = NULL; DeleteData = NULL;

//	cout << "\nThe TreeHMM: " << lnL();
//	exit(-1);

//	cout << "\nThe model is: " << *this;
}

// Function for replacing parameter values
bool CBaseModel::ReplaceParValue(CPar *Par)	{
	int i;
	bool Found = false;
	FOR(i,NoPar()) {
		if(m_vpPar[i]->Name().find(Par->Name()) != string::npos) {
//			cout << "\nFound match: " << m_vpPar[i]->Name() << " cf. " << Par->Name();
			m_vpPar[i]->SetVal(Par->Val(),true,true,true);
			Found = true;
		}
	}
	return Found;
}

//////////////////////////////////////////////
// Controllers for DNA models

CDNAProcess * CBaseModel::AddDNAProcess(CData *Data,CTree *Tree, DNAProc Model, string Name)	{
	CDNAProcess *Proc;
	Proc = new CDNAProcess(Data,Tree,Model,Name);
	return Proc;
}

CJC::CJC(CData *Data, CTree *Tree) : CBaseModel(Data,Tree) {
	m_sName = sModelNames[(int)JC];
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pJC));
	FinalInitialisation();
}

CFEL::CFEL(CData *Data, CTree *Tree) : CBaseModel(Data,Tree) {
	m_sName = sModelNames[(int)FEL];
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pFEL));
	FinalInitialisation();
}

CK2P::CK2P(CData *Data, CTree *Tree) : CBaseModel(Data,Tree) {
	m_sName = sModelNames[(int)K2P];
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pK2P));
	FinalInitialisation();
}

CHKY::CHKY(CData *Data,CTree *Tree) : CBaseModel(Data,Tree)	{
	m_sName = sModelNames[(int)HKY];
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pHKY));
	FinalInitialisation();
}

CREV::CREV(CData *Data,CTree *Tree) : CBaseModel(Data,Tree)	{
	m_sName = sModelNames[(int)REV];
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pREV));
	// Get the probability vector sorted
	FinalInitialisation();
}
void CREV::CompressModel()	{ // Function to compress the model
	int i,j,proc;
	string::size_type loc;
	CQPar *Par;
	if(Compressed()) { return; }
	if(m_pData->m_DataType == DNA)	{ // Do DNA compression
		// Fix the transition probabilities in the matrix
		FOR(i,(int)m_vpPar.size())	{
			loc = m_vpPar[i]->Name().find("<->");
			if(loc != string::npos)	{ m_vpPar[i]->SetOptimise(false); }
		}
		// Now add a single ts/tv scaling parameters
		FOR(proc,(int)m_vpProc.size())	{
			if(m_vpProc[proc]->PseudoProcess() == true) { continue; }
			loc = m_vpProc[proc]->Name().find("REV");
			if(loc != string::npos && m_vpProc[proc]->Char() == 4) {	// Do the addition of kappa
				Par = new CQPar("PseudoKappa",4,1.0);
				FOR(i,4)	{ for(j=i+1;j<4;j++) { if(IsTs(State(DNA,i),State(DNA,j),DNA))	{ Par->AddQij(i,j); } } }
				m_vpProc[proc]->AddQPar(Par);
				m_vpPar.push_back(Par);
				Par = NULL;
	}	}	}
	SetCompressed();
}
void CREV::UncompressModel()	{ // Function to uncompress the model
	int i,j,proc;
	if(!Compressed()) { return; }
	if(m_pData->m_DataType == DNA)	{ // Do DNA compression
		// Fix the transition probabilities in the matrix
		FOR(i,(int)m_vpPar.size())	{
			if(m_vpPar[i]->Name().find("<->") != string::npos)	{ m_vpPar[i]->SetOptimise(true); }
		}
		// Now add remove the single ts/tv scaling parameters
		FOR(proc,(int)m_vpProc.size())	{
			if(m_vpProc[proc]->PseudoProcess() == true) { continue; }
			if(m_vpProc[proc]->Name().find("REV") != string::npos && m_vpProc[proc]->Char() == 4) {
				FOR(i,m_vpProc[proc]->NoPar()) {
					if(m_vpProc[proc]->pPar(i)->Name().find("PseudoKappa") != string::npos) { break; }
				}
				if(i == m_vpProc[proc]->NoPar()) { Error("\nCan't find PseudoKappa in CREV::UncompressModel()...\n\n"); }
				FOR(j,m_vpProc[proc]->NoPar()) {
					if(	m_vpProc[proc]->pPar(j)->Name().find("A<->G") != string::npos || m_vpProc[proc]->pPar(j)->Name().find("G<->A") != string::npos ||
						m_vpProc[proc]->pPar(j)->Name().find("T<->C") != string::npos || m_vpProc[proc]->pPar(j)->Name().find("C<->T") != string::npos) {
							m_vpProc[proc]->pPar(j)->SetVal(m_vpProc[proc]->pPar(j)->Val() * m_vpProc[proc]->pPar(i)->Val());
				}	}
		}	}
		// Remove the PseudoKappa Par
		RemovePar("PseudoKappa");
	}
	SetUncompressed();
}
///////////////////////////////////////////////////////////////
// The RY model
CRY::CRY(CData *D,CTree *T) : CBaseModel(D,T)	{
	m_sName = sModelNames[(int)RY_model];
	if(D->m_DataType == DNA) { D->DNA2RY(); }
	if(D->m_DataType != RY)	{ Error("\nCan only create RY model for DNA data"); }
	m_vpProc.push_back(AddDNAProcess(D,T,pRY));
	m_vpProc[0]->pPar(0)->Name(m_vpProc[0]->pPar(0)->Name().replace(0,4,"Prob"));
	m_vpProc[0]->pPar(1)->Name(m_vpProc[0]->pPar(1)->Name().replace(0,4,"Prob"));
	FinalInitialisation();
}

///////////////////////////////////////////////////////////////
// Temporal HMM models
//	Functions return the name of the model

// --- Controller for AA based THMM models ---
string CBaseModel::AddAATHMMProcess(double Alfa, int CatGam, double SigAlfa, double pInv, double SigInv, bool VarySig,bool VarypInv)	{
	string Name = "THMM_AA";
	// Check entry conditions
	assert(m_pTree != NULL && m_pData!= NULL);
	// Create the AA process
	CBaseProcess *Proc = NULL;
	CWARSProcess *AAProc = NULL;
	CBranchWARSProcess *BAAProc = NULL;
	if(fabs(SigInv) > FLT_EPSILON || fabs(pInv) > FLT_EPSILON) { Name += "."; }
	if(fabs(SigInv) > FLT_EPSILON) { Name += "S"; }
	if(fabs(pInv) > FLT_EPSILON) { Name += "I"; }
	if(VarySig || VarypInv) {
		if(fabs(pInv) < FLT_EPSILON && fabs(SigInv) < FLT_EPSILON) { Name += ".NONE"; }
		if(VarySig && VarypInv)	{ Name += ".SI"; }
		else if(VarySig)		{ Name += ".S";  }
		else if(VarypInv)		{ Name += ".I"; }
		BAAProc = new CBranchWARSProcess(m_pData, m_pTree,Name,Alfa,CatGam,SigAlfa,pInv,SigInv,VarySig,VarypInv);
		Proc = BAAProc;
		BAAProc = NULL;
	} else {
		AAProc = new CWARSProcess(m_pData, m_pTree,Name,Alfa,CatGam,SigAlfa,pInv,SigInv);
		Proc = AAProc;
		AAProc = NULL;
	}
	// Add the process to the list
	m_vpProc.push_back(Proc);
	Proc = NULL;
	// Return name
	return Name;
}

// --- Controller for DNA based THMM models ---
string CBaseModel::AddDNATHMMProcess(double Alfa, int CatGam, double SigAlfa, double pInv, double SigInv, bool VarySig,bool VarypInv)	{
	string Name = "THMM_DNA";
	// Check entry conditions
	assert(m_pTree != NULL && m_pData!= NULL);
	// Create the DNA process
	CBaseProcess *Proc = NULL;
	CDNAWARSProcess *DNAProc= NULL;
	CDNABranchWARSProcess *BDNAProc= NULL;
	if(fabs(SigInv) > FLT_EPSILON || fabs(pInv) > FLT_EPSILON) { Name += "."; }
	if(fabs(SigInv) > FLT_EPSILON) { Name += "S"; }
	if(fabs(pInv) > FLT_EPSILON) { Name += "I"; }
	if(VarySig || VarypInv) {
		if(fabs(pInv) < FLT_EPSILON && fabs(SigInv) < FLT_EPSILON) { Name += ".NONE"; }
		if(VarySig && VarypInv)	{ Name += ".SI"; }
		else if(VarySig)		{ Name += ".S";  }
		else if(VarypInv)		{ Name += ".I"; }
		BDNAProc= new CDNABranchWARSProcess(m_pData, m_pTree,Name,Alfa,CatGam,SigAlfa,pInv,SigInv,VarySig,VarypInv);
		Proc = BDNAProc;
		BDNAProc= NULL;
	} else {
		DNAProc= new CDNAWARSProcess(m_pData, m_pTree,Name,Alfa,CatGam,SigAlfa,pInv,SigInv);
		Proc = DNAProc;
		DNAProc= NULL;
	}
	// Add the process to the list
	m_vpProc.push_back(Proc);
	Proc = NULL;
	// Return name
	return Name;
}

// --- Controller for full DNA based THMM models ---
string CBaseModel::AddFullDNATHMMProcess(int NoProc, EHIDDEN_TYPE DoHidden, ETHMM_EQM_TYPE DoFreqs, bool DoKappa, ERateTypes DoRates,vector <double> Rates)	{
	int i;
	vector <double> Probs;
	string Name;

	// Check some entry conditions
	assert(m_pTree != NULL && m_pData != NULL);
	// Do the process (pseudo count of 0.25 to make sure no probabilities are too small for initial optimisation)
	if(DoRates) { FOR(i,NoProc) { Probs.push_back(1.0); } }
	else		{ FOR(i,NoProc) { Probs.push_back(0.25 + Random()); } }
	NormaliseVector(Probs);
	Name = "THMM.";
	Name = Name + int_to_string(NoProc) + ".";
	if(DoHidden == H_diff && DoFreqs == complex && DoKappa && DoRates) { Name = Name + "All"; } else {
		if(DoFreqs == complex)	{ Name = Name + "F"; }
		else if(DoFreqs == equ)	{ Name = Name + "-F"; }
		if(DoKappa)				{ Name = Name + "K"; }
		if(DoRates == varyall)	{ Name = Name + "R"; }
		if(DoRates == gammarates) { Name = Name + "dGR"; }
		if(DoHidden == H_diff)  { Name = Name + "H"; }
		else if(DoHidden == H_none) { Name = Name + "-H"; }
	}
	Probs = NormaliseVector(Probs);
	// Initialise the processes
	CFullDNATHMMProcess*CovP = NULL;
	CovP  = new CFullDNATHMMProcess(m_pData,m_pTree,Name,NoProc,Probs,DoFreqs,true,DoHidden,DoRates);
	CovP->SetRateScaling(true);
	if(DoKappa)	{ CovP->AddSeperateKappas(); } else { CovP->AddKappa(); }
	if(DoRates)	{
#if HMP_MODEL_2_SIMPLE == 1
		FOR(i,NoProc) { Probs[i] = 1.0; }
#else
		FOR(i,NoProc) { Probs[i] = Random() + (i*0.1); }
#endif
		Sort(&Probs);
		FOR(i,NoProc) { CovP->SetSubProcRate(i,Probs[i]); }
		CovP->OptimiseRates();
	}
	if(DoHidden != H_none) { CovP->OptimiseHiddenProbs(); }
	// Do Rates (if rqd)
	if(!Rates.empty())	{
		if((int) Rates.size() != NoProc) { Error("\nNumber of rates[" + int_to_string((int)Rates.size()) + "] should equal number of processes [" + int_to_string(NoProc) + "]"); }
		// Normalise rates so first is the largest
		RSort(&Rates);
		cout << "\nOriginal rates: " << Rates;
		rFOR(i,(int)Rates.size())	{
			CovP->SetSubProcRate(i,Rates[i] / Rates[0]); }
	}
	// Finish up
	m_vpProc.push_back(CovP);
	CovP = NULL;
	return Name;
}
string CBaseModel::AddFullDNATHMMProcess(int NoProc, EHIDDEN_TYPE DoHidden, ETHMM_EQM_TYPE DoFreqs, bool DoKappa, ERateTypes DoRates)	{
	vector <double> Rates;
	return AddFullDNATHMMProcess(NoProc, DoHidden, DoFreqs, DoKappa, DoRates,Rates);
}

string CBaseModel::AddAACoevoProcess(CData * JointData, vector <double> Left, vector <double> Right, vector <double> *PropMat, double init_psi) {
	CAACoevoProcess * Proc;
	Proc = new CAACoevoProcess(JointData,Left,Right,m_pTree,PropMat, init_psi);
	m_vpProc.push_back(Proc);
	Proc = NULL;
	return "Coevo_WAG";
}

///////////////////////////////////////////////////////////////////////////////
// Models being developed to investigate one off temporal events in AA
// 1. The constructor
CTHMMAA::CTHMMAA(CData *D, CTree *T, int MatType, double OriAlpha, int NoCatAlpha, double ProbInv, double SigAlpha, double SigInv) : CBaseModel(D,T)	{
	m_sName = AddAATHMMProcess(OriAlpha,NoCatAlpha,SigAlpha,ProbInv,SigInv,false,false);
	m_MatType = MatType; m_NoCatAlpha = NoCatAlpha;
	m_OriAlpha = OriAlpha; m_ProbInv = ProbInv; m_SigAlpha = SigAlpha; m_SigInv = SigInv;
	FinalInitialisation();
}
// Return the preoptimiser model (which will be the simple THMM of whatever is used)
CBaseModel *CTHMMAA::PreOptModel()	{
	if(m_vpProc.empty()) { Error("Unexpected empty processes...");}
	else if(m_vpProc.size()!=1) { Error("Preoptimisation expects only a single process");}
	if(m_pData == NULL || m_pTree == NULL) { Error("\nCannot run preopt without specifying data and tree...\n\n"); }
	CWAG *PreModel;
	PreModel = new CWAG(m_pData,m_pTree,true);
	PreModel->MakeGammaModel(0,4);
	if(m_ProbInv > FLT_EPSILON) { PreModel->MakeInvariantSitesModel(0); }
	return PreModel;
}
// Apply the preoptimiser parameters to the current model
void CTHMMAA::ApplyPreOptModel(CBaseModel *PreModel)	{
	int i, NumSigma = 0;
	double Sigma = -100;
	// Function
	FOR(i,PreModel->NoPar()) {
		if(PreModel->m_vpPar[i]->Name().find("P(Inv)") != string::npos) { ReplaceParValue(PreModel->m_vpPar[i]); }
		if(PreModel->m_vpPar[i]->Name().find("Alpha") != string::npos) { ReplaceParValue(PreModel->m_vpPar[i]); }
	}
}

///////////////////////////////////////////////////////////////////////////////
// As above, but things vary over branches
CTHMMBAA::CTHMMBAA(CData *D, CTree *T, int MatType, double OriAlpha, int NoCatAlpha, double ProbInv, double SigAlpha, double SigInv,bool VarySig, bool VarypInv) : CBaseModel(D,T)	{
	int i, NumSigma = 0;
	m_sName = AddAATHMMProcess(OriAlpha,NoCatAlpha,SigAlpha,ProbInv,SigInv,VarySig,VarypInv);
	m_MatType = MatType; m_NoCatAlpha = NoCatAlpha;
	m_OriAlpha = OriAlpha; m_ProbInv = ProbInv; m_SigAlpha = SigAlpha; m_SigInv = SigInv;
	m_bAllowParBraOpt = true;
	FinalInitialisation();
	// Need function for comparing parameters between models
	FOR(i,NoPar()) { if(m_vpPar[i]->Name().find("Sigma") != string::npos) { NumSigma ++; }	}
	if(NumSigma != m_pTree->NoBra()) { m_bAllowParBraOpt = false; }
}

// Return the preoptimiser model (which will be the simple THMM of whatever is used)
CBaseModel *CTHMMBAA::PreOptModel()	{
	if(m_vpProc.empty()) { Error("Unexpected empty processes...");}
	else if(m_vpProc.size()!=1) { Error("Preoptimisation expects only a single process");}
	if(m_pData == NULL || m_pTree == NULL) { Error("\nCannot run preopt without specifying data and tree...\n\n"); }
	CTHMMAA *PreModel;
	PreModel = new CTHMMAA(m_pData,m_pTree,m_MatType,m_OriAlpha, m_NoCatAlpha, m_ProbInv, m_SigAlpha, m_SigInv);
	return PreModel;
}
// Apply the preoptimiser parameters to the current model
void CTHMMBAA::ApplyPreOptModel(CBaseModel *PreModel)	{
	int i, NumSigma = 0;
	double Sigma = -100;

	// Need function for comparing parameters between models
	FOR(i,PreModel->NoPar()) {
		if(PreModel->m_vpPar[i]->Name().find("P(Inv)") != string::npos) { ReplaceParValue(PreModel->m_vpPar[i]); }
		if(PreModel->m_vpPar[i]->Name().find("Alpha") != string::npos) { ReplaceParValue(PreModel->m_vpPar[i]); }
		if(PreModel->m_vpPar[i]->Name().find("Sigma") != string::npos) {
			if(Sigma < 0) { Sigma = PreModel->m_vpPar[i]->Val(); }
			else		  { Sigma += PreModel->m_vpPar[i]->Val(); }
			NumSigma ++;
	}	}
	// Apply the single sigma to the sets of sigma
	Sigma /= (double) NumSigma;
	if(Sigma < 0.05) { Sigma = 0.05; }
	FOR(i,NoPar()) {
		if(m_vpPar[i]->Name().find("Sigma") != string::npos) {
			m_vpPar[i]->SetVal((Sigma - 0.005) + (0.01*Random()), true, true);
}	}	}

// Likelihood functions that may contain a penalty term
double CTHMMBAA::DoBralnL(int B, int NL, int NR, bool JustClean) {
	double Pen = 0;
	int i;
#if ALLOW_BAA_PEN == 1
	FOR(i,(int)m_vpProc.size()) { Pen += m_vpProc[0]->Penalty(); }

	return Pen + CBaseModel::DoBralnL(B,NL,NR,JustClean);
#else
//	m_vpProc[0]->Make_PT(B);
	return CBaseModel::DoBralnL(B,NL,NR,JustClean);
#endif
}
#if ALLOW_BAA_PEN == 1
double CTHMMBAA::lnL(bool FR) {
	double Pen = 0;
	int i;
	FOR(i,(int)m_vpProc.size()) { Pen += m_vpProc[0]->Penalty(); }
	return Pen + CBaseModel::lnL(FR);
}
#endif

void CTHMMBAA::DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol, bool AllowUpdate) {
	int i;
	// Start by optimising the branch length
//	cout << "\nNot dealing with branches properly...";
//	cout << "\n--- Doing DoBraOpt[" << Br << "] --- ";
	// Find the right sigma for the branch; if there's no natural corresponding sigma then it means there's probably groups of sigma so return.
	string Name = "Sigma[" + int_to_string(Br) + "]";
	FOR(i,m_vpPar.size()) { if(m_vpPar[i]->Name().find(Name) != string::npos) { break; } }
	if(i == (int)m_vpPar.size() || m_bAllowParBraOpt == false) { return CBaseModel::DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol,AllowUpdate); }
	CBaseModel::DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol,false);
//	cout << "\nHave found ["<< Br << "]: " << flush; cout << *m_vpPar[i];
	// Initialise variables
	CPar *Par  = m_vpPar[i];
	double *x = Par->OptimiserValue(),x1,x2=*x,x3,x1_lnL = 1.0,x2_lnL = -fabs(*BestlnL),x3_lnL = 1.0, xi,temp;
	tol = max(tol,FULL_LIK_ACC);
	double dx = tol * 10;

#if DEVELOPER_BUILD == 1
	cout << "\nBranch["<<Br<<"] has DoBralnL: "<< DoBralnL(Br,NFr,NTo) << " cf. " << DoBralnL(Br,NFr,NTo);
	cout << "\nReturning from CBaseModel::DoBraOpt (including branch updates)";
	RETURN_DOBRAOPT;
#endif
//	if(Br == 2) { RETURN_DOBRAOPT; }
//	cout << "\nDoing branch["<<Br<<"]: sent bestlnL: " << *BestlnL << " cf. DoBralnL(Br,NFr,NTo): " << DoBralnL(Br,NFr,NTo); //  << " cf. real " << lnL(); exit(-1);
//	if(IsSubTree()) { cout << "\nFinished for checking... " << lnL(true); exit(-1); }
//	RETURN_DOBRAOPT;
//	cout << "\n\n---\nDoing GoldenSection";
//	cout << "\nDoing FastBranchOpt ["<<Br<<"]; tol=" <<tol <<": ";
//	cout << "\nPar = " << Par->Name() << "; x2: " << x2 << " == " << x2_lnL << " == " << DoBralnL(Br,NFr,NTo);
	if(fabs(DoBralnL(Br,NFr,NTo) - *BestlnL) > 1.0E-5) {
//		cout << " ... broken"; RETURN_DOBRAOPT;
		cout << "\nLikelihoods don't match in CTHMMBAA::DoBraOpt";
//		cout << "\n\tBestlnL: \t" << *BestlnL << "\n\tDoBralnL: \t" << DoBralnL(Br,NFr,NTo); cout << "\n\tReallnL: \t" << lnL(true);
		cout << "\n\tBestlnL: \t" << *BestlnL << "\n\tDoBralnL: \t" << DoBralnL(Br,NFr,NTo);
//		m_vpProc[0]->OutPT(cout,Br);
		RETURN_DOBRAOPT;
		exit(-1);
	}
/*	else {
		cout << " ... fine";  RETURN_DOBRAOPT;
	}*/

	// Get left bracketing
	while(x2_lnL < x1_lnL)	{
		*x = x1 = x2 - max(fabs(x2) * dx,dx);
		x1_lnL = DoBralnL(Br,NFr,NTo);
//		cout << "\nx = " << x1 << " : " << x1_lnL << flush;
		if(x1_lnL > x2_lnL)	{
			x3 = x2; x3_lnL = x2_lnL;
			x2 = x1; x2_lnL = x1_lnL;
			x1_lnL = 1.0;
		}
		dx *= GS_DELTA;
		if(!Par->CheckLowBound()) { break; }
	}
	if(x3_lnL > 0.0)	{
		dx = tol * 10;
		// Get right bracketing
		while(x2_lnL < x3_lnL)	{
			*x = x3 = x2 + max(fabs(x2) * dx,dx);
			x3_lnL = DoBralnL(Br,NFr,NTo);
//			cout << "R" << flush;
//			cout << "\nRight: " << x3 << " == " << x3_lnL;
			if(x3_lnL > x2_lnL)	{
				x1 = x2; x1_lnL = x2_lnL;
				x2 = x3; x2_lnL = x3_lnL;
				x3_lnL = 1.0;
			}
			dx *= GS_DELTA;
			if(!Par->CheckUpBound()) { break; }
	}	}
//	cout << "\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
//	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
//	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
//	exit(-1);
	// Do various checks to make sure it looks like a hill
	if(!Par->CheckBound()) { *x = x2; *BestlnL = x2_lnL; RETURN_DOBRAOPT; }
	if(x1_lnL > max(x2_lnL,x3_lnL)) {
		// Deal with the U shape case
		if(x3_lnL > x2_lnL) {
			if(x1_lnL > x3_lnL) { *x = x1; DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol); }
			else				{ *x = x3; DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol); }
		} else { *x = x1; *BestlnL = x1_lnL; RETURN_DOBRAOPT; }
	}
	else if (x3_lnL > x2_lnL)		{ *x = x3; *BestlnL = x3_lnL; RETURN_DOBRAOPT; }
	// Now do golden section search to find optimal value
	assert(x2_lnL + FLT_EPSILON > x1_lnL && x2_lnL + FLT_EPSILON> x3_lnL);
	FOR(i,20)	{
//		cout << "." << flush;
		if(fabs(x3 - x2) > fabs(x2 - x1))	{ *x = xi = x2 + ( ( x3 - x2) * GOLDEN_NUMBER); }
		else								{ *x = xi = x2 - ( ( x2 - x1) * GOLDEN_NUMBER); }
		temp = DoBralnL(Br,NFr,NTo);
		if(temp > x2_lnL) {
			if(fabs(x3 - x2) > fabs(x2 - x1))	{ x1 = x2; x1_lnL = x2_lnL; x2 = xi; x2_lnL = temp; }
			else								{ x3 = x2; x3_lnL = x2_lnL; x2 = xi; x2_lnL = temp; }
		} else		  {
			if(fabs(x3 - x2) > fabs(x2 - x1))	{ x3 = xi; x3_lnL = temp; }
			else								{ x1 = xi; x1_lnL = temp; }
		}
		if(fabs(x2_lnL - x3_lnL) < tol && fabs(x2_lnL - x1_lnL) < tol && i > 2)	{ break; }
	}
//	cout << ": tol= " << max(x2_lnL - x1_lnL,x2_lnL - x3_lnL);
//	cout << "\n---\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
//	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
//	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
	*x = x2; *BestlnL = x2_lnL;
//	cout << "\nx = " << *x << " lnL: " << *BestlnL << " == " << DoBralnL(Br,NFr,NTo);
//	cout << " == " << lnL(true);
//	exit(-1);
	// Finish by doing the calculation again to correctly update the partial likelihoods
	DoBralnL(Br,NFr,NTo);
	RETURN_DOBRAOPT;



}

///////////////////////////////////////////////////////////////////////////////
// Models being developed to investigate one off temporal events in DNA
// 1. The constructor
CTHMMDNA::CTHMMDNA(CData *D, CTree *T, double OriAlpha, int NoCatAlpha, double ProbInv, double SigAlpha, double SigInv) : CBaseModel(D,T)	{
	m_sName = AddDNATHMMProcess(OriAlpha,NoCatAlpha,SigAlpha,ProbInv,SigInv,false,false);
	FinalInitialisation();
}

///////////////////////////////////////////////////////////////////////////////
// As above, but things vary over branches
CTHMMBDNA::CTHMMBDNA(CData *D, CTree *T, double OriAlpha, int NoCatAlpha, double ProbInv, double SigAlpha, double SigInv,bool VarySig, bool VarypInv) : CBaseModel(D,T)	{
	m_sName = AddDNATHMMProcess(OriAlpha,NoCatAlpha,SigAlpha,ProbInv,SigInv,VarySig,VarypInv);
	FinalInitialisation();
	cout << "\nThe model is: " << flush;
	cout << *this;
	cout << "\n------- DONE ------" << flush;
	exit(-1);

}
#if ALLOW_BDNA_PEN == 1
// Likelihood functions for penalty term
double CTHMMBDNA::DoBralnL(int B, int NL, int NR, bool JustClean) {
	double Pen = 0;
	int i;
	FOR(i,(int)m_vpProc.size()) { Pen += m_vpProc[0]->Penalty(); }
	return Pen + CBaseModel::DoBralnL(B,NL,NR, JustClean);
}
double CTHMMBDNA::lnL(bool FR) {
	double Pen = 0;
	int i;
	FOR(i,(int)m_vpProc.size()) { Pen += m_vpProc[0]->Penalty(); }
	return Pen + CBaseModel::lnL(FR);
}
#endif

/////////////////////////////////////////////////////////////////////////////
// Model mimicking the Galtier model for AA sequences with 4 discrete
//	gamma rates that vary through time
CWAGdG_THMM::CWAGdG_THMM(CData *D, CTree *T, int MatType, double OriAlpha, int NoCatAlpha) : CTHMMAA(D, T, MatType, OriAlpha, NoCatAlpha, 0, 0.1, 0) {
	m_sName = "WAGdG_THMM";
	FinalInitialisation();
}

CWAGdG_THMM::~CWAGdG_THMM() {
}
//////////////////////////////////////////////////////////////////////////////
// New models used to investigate temporal and spatial heterogeneity in DNA
// 1. The general constructor
CTHMM_FULL::CTHMM_FULL(CData *D, CTree *T, int NoProc, ERateTypes DoRates, vector <double> *Rates, ETHMM_EQM_TYPE DoFreq, bool DoKappa, EHIDDEN_TYPE DoHidden) : CBaseModel(D,T)	{
	m_bRates = false;
	if(Rates == NULL)	{ m_sName = AddFullDNATHMMProcess(NoProc,DoHidden,DoFreq,DoKappa,DoRates); }
	else				{ m_sName = AddFullDNATHMMProcess(NoProc,DoHidden,DoFreq,DoKappa,DoRates,*Rates); }
	if(DoRates != same) { m_bRates = true; }
	m_bDoSepParOpt = true;
	FinalInitialisation();
}
void CTHMM_FULL::CompressModel()	{ // Function to compress the model
	int i,proc;
	m_bDoSepParOpt = false;
	FOR(proc,(int)m_vpProc.size())	{
		FOR(i,(int)m_vpProc[proc]->NoPar())	{
			m_vpProc[proc]->pPar(i)->SetOptimise(false);
}	}	}
void CTHMM_FULL::UncompressModel() { // Function to uncompress the model
	int i,proc;
	m_bDoSepParOpt = true;
	FOR(proc,(int)m_vpProc.size())	{
		FOR(i,(int)m_vpProc[proc]->NoPar())	{
			m_vpProc[proc]->pPar(i)->SetOptimise(true);
}	}	}
// Get the preoptimiser model
CBaseModel *CTHMM_FULL::PreOptModel()	{
	CHKY *PreModel;
	if(m_pData == NULL || m_pTree == NULL) { Error("\nCannot run preopt without specifying data and tree...\n\n"); }
	PreModel = new CHKY(m_pData,m_pTree);
	// Add gamma distribution for model if required.
	if(m_bRates) { PreModel->MakeGammaModel(0,m_vpProc[0]->Char() / m_pData->m_iChar); }
	return PreModel;
}
// Apply the Preoptimised model to the current model
/////////////////////////////////////////////////////////////
// This code assumes certain standard parameter names
//
// 1. All Kappa parameters include "Kappa" in their name
// 2. All Rate parameters have "OverallRate" in their name
void CTHMM_FULL::ApplyPreOptModel(CBaseModel *PreOpt)		{
	int i,j,Pars;
	double Val;
	vector <double> Vals;
	// Do the kappa
	FOR(i,PreOpt->NoPar()) { if(PreOpt->m_vpPar[i]->Name().find("Kappa") != string::npos) { Val = PreOpt->m_vpPar[i]->Val(); break; } }
	if(i==PreOpt->NoPar()) { Error("\nCouldn't find Kappa parameter in HKY..."); }
	Pars = 0; FOR(i,NoPar())	{
		if(m_vpPar[i]->Name().find("Kappa") != string::npos) {
			if(Pars > 0) { m_vpPar[i]->SetVal(Val + (0.25 * Random()));} else { m_vpPar[i]->SetVal(Val); }
			Pars++;
	}	}
	if(Pars == 0) { Error("\nCouldn't find Kappa parameter in THMM..."); }
	// Do the frequencies
	Vals = PreOpt->m_vpProc[0]->Eqm(0);
	// Do the eqm distributions
	FOR(i,(int)m_vpProc[0]->EqmPar().size()) {
		if(i > 0)	{ m_vpProc[0]->EqmPar()[i]->ResetEqm(Vals,true); }
		else		{ m_vpProc[0]->EqmPar()[i]->ResetEqm(Vals,false); }
	}
	// Do the rates
	if(m_bRates)	{
		Vals.clear();
		// Get the rates
		Vals = PreOpt->Rates();
		for(i=1;i<(int)Vals.size();i++) { Vals[i] /= Vals[0]; }
		j = 1;
		FOR(i,NoPar())	{
			if(m_vpPar[i]->Name().find("OverallRate") == string::npos || m_vpPar[i]->Special()) { continue; }
			m_vpPar[i]->SetVal(Vals[j++] + 0.1 * (Random()));
		}
	}
}

////////////////////////////////////////////////////////////////
// A simple THMM for phylogenetic analysis based around HKY+dG
CHKYdG_THMM::CHKYdG_THMM(CData *D, CTree *T, int NoCat) : CTHMM_FULL(D,T,NoCat,gammarates,NULL,equ,false,H_same)	{
	cout << "\nMaking CHKYdG_THMM";
	m_sName = "HKYdG_THMM";
//	exit(-1);
}

////////////////////////////////////////////////////////////////
// DNA model for detecting changes in function
// ---
// NoProc:		The number of hidden processes the exon finder uses
// LowRange:	The rate of the slowest process
// UpRange:		The rate of the fastest process
//
// The rate of all other processes are filled in to be equally spaced apart
// These rates will *NOT* be optimised at any point
CExonDetector::CExonDetector(CData *D, CTree *T, int NoProc,double LowRange, double UpRange,vector <double> *PassedProbs, bool OptProbs) : CBaseModel(D,T)	{
	int i;
	double Rate,RateInc;
	vector <double> Probs;
	// Check entry conditions
	assert(D->m_DataType == DNA);
	assert(LowRange >= 0 && UpRange > LowRange);
	// Set initial probabilities
	if(PassedProbs == NULL) { FOR(i,NoProc) { Probs.push_back(Random()); } if(!OptProbs) { Error("\nFixing probabilities in ::ExonDetector when they are set to random...\n\n"); } }

	else { assert((int)PassedProbs->size() == NoProc); FOR(i,NoProc) { Probs.push_back(PassedProbs->at(i)); } }
	// Set up process
	m_sName = "THMM for testing functionality of exons with " + int_to_string(NoProc) + " hidden states";
	Probs = NormaliseVector(Probs);
	CFullDNATHMMProcess*CovP = NULL;
	CovP  = new CFullDNATHMMProcess(D,T,"REV[NoHiddenProcesses=" + int_to_string(NoProc) + "]",NoProc,Probs,obs,true);
	CovP->SetHiddenProbsOpt(OptProbs);
	// Organise rates
	RateInc = ((UpRange - LowRange) / (double) (NoProc - 1) ); i = 0;
	for(Rate = LowRange;Rate<UpRange - FLT_EPSILON;Rate += RateInc) { CovP->SetSubProcRate(i++,Rate); }
	CovP->SetRateScaling(true);
	CovP->NoOptimiseRates();
	CovP->AddREV();
	m_vpProc.push_back(CovP);
	CovP = NULL;
	FinalInitialisation();
}
/////////////////////////////////////////////////////
// DNA covarion models
CClassicCovJC::CClassicCovJC(CData *D, CTree *T,double Prob1of2)	: CBaseModel(D,T) {
	vector <double> Probs;
	m_sName = "Classic_covarion_JC";
	if(!InRange(Prob1of2,0.0,1.0)) { Error("\n\nProb1of2 in CCovJC out of range"); }
	Probs.push_back(Prob1of2); Probs.push_back(1.0-Prob1of2);
	Probs = NormaliseVector(Probs);
	C2StateCovProc *CovP = NULL;
	CovP  = new C2StateCovProc(D,T,"Classic covarion JC process",Probs,equ);
	CovP->MakeSimpleCovarionModelEqm();
	CovP->EnforceFixedState();
	CovP->SetRateScaling(false);
	m_vpProc.push_back(CovP); CovP->SetSubProcRate(1,0);
	CovP = NULL;
	FinalInitialisation();
}
CCovJC::CCovJC(CData *D, CTree *T,double Prob1of2)	: CBaseModel(D,T) {
	vector <double> Probs;
	m_sName = "Covarion_JC";
	if(!InRange(Prob1of2,0.0,1.0)) { Error("\n\nProb1of2 in CCovJC out of range"); }
	Probs.push_back(Prob1of2); Probs.push_back(1.0-Prob1of2);
	Probs = NormaliseVector(Probs);
	C2StateCovProc *CovP = NULL;
	CovP  = new C2StateCovProc(D,T,"Covarion JC process",Probs,equ);
	CovP->EnforceFixedState();
	CovP->SetRateScaling(false);
	m_vpProc.push_back(CovP); CovP->SetSubProcRate(1,0);
	CovP = NULL;
	FinalInitialisation();
}
CClassicCovHKY::CClassicCovHKY(CData *D, CTree *T, double Prob1of2) : CBaseModel(D,T) {
	vector <double> Probs;
	m_sName = "Classic_covarion_HKY";
	if(!InRange(Prob1of2,0.0,1.0)) { Error("\n\nProb1of2 in CCovJC out of range"); }
	Probs.push_back(Prob1of2); Probs.push_back(1.0-Prob1of2);
	Probs = NormaliseVector(Probs);
	C2StateCovProc *CovP = NULL;
	CovP  = new C2StateCovProc(D,T,"Classic covarion HKY process",Probs,obs);
	CovP->MakeSimpleCovarionModelEqm();
	CovP->EnforceFixedState();
	CovP->SetRateScaling(false); CovP->SetSubProcRate(1,0);
	CovP->AddKappa();
	m_vpProc.push_back(CovP);
	CovP = NULL;
	FinalInitialisation();
}
CCovHKY::CCovHKY(CData *D, CTree *T, double Prob1of2) : CBaseModel(D,T) {
	vector <double> Probs;
	m_sName = sModelNames[(int)COV_HKY];
	if(!InRange(Prob1of2,0.0,1.0)) { Error("\n\nProb1of2 in CCovHKY out of range"); }
	Probs.push_back(Prob1of2); Probs.push_back(1.0-Prob1of2);
	Probs = NormaliseVector(Probs);
	C2StateCovProc *CovP = NULL;
	CovP  = new C2StateCovProc(D,T,"Covarion HKY process",Probs,obs);
	CovP->SetRateScaling(true);
	CovP->AddKappa();
	CovP->OptimiseHiddenProbs();
	m_vpProc.push_back(CovP);
	CovP = NULL;
	FinalInitialisation();
}

CCovREV::CCovREV(CData *D, CTree *T, double Prob1of2) : CBaseModel(D,T)	{
	vector <double> Probs;
	m_sName = sModelNames[(int)COV_REV];
	if(!InRange(Prob1of2,0.0,1.0)) { Error("\n\nProb1of2 in CCovHKY out of range"); }
	Probs.push_back(Prob1of2); Probs.push_back(1.0-Prob1of2);
	Probs = NormaliseVector(Probs);
	C2StateCovProc *CovP = NULL;
	CovP  = new C2StateCovProc(D,T,"Covarion HKY process",Probs,obs);
	CovP->SetRateScaling(true);
	CovP->AddREV();
	CovP->OptimiseHiddenProbs();
	m_vpProc.push_back(CovP);
	CovP = NULL;
	FinalInitialisation();

}


/////////////////////////////////////////////////////
// Various rate models

C3RateFEL::C3RateFEL(CData *Data, CTree *Tree) : CBaseModel(Data,Tree)	{
	int i;

	m_sName = "2RateFEL";
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pFEL,"FEL_proc(1)"));
	m_vpProc.push_back(m_vpProc[0]->RateProcessCopy());
	m_vpProc.push_back(m_vpProc[0]->RateProcessCopy());


	m_vpProc[0]->Rate(0.072181); m_vpProc[0]->ProbPar()->SetVal(0.001650,true,true);
	m_vpProc[1]->Rate(0.518536); m_vpProc[1]->ProbPar()->SetVal(0.743378,true,true);
	m_vpProc[2]->Rate(2.390653); m_vpProc[2]->ProbPar()->SetVal(0.254971,true,true);

	FOR(i,3) { cout << "\nProc["<<i<<"]: rate= " << m_vpProc[i]->Rate() << "; prob= " << m_vpProc[i]->Prob(); }

	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->AddRatePar2Opt(); }

	FinalInitialisation();

}

C2RateHKY::C2RateHKY(CData *Data, CTree *Tree) : CBaseModel(Data,Tree) {
	int i;
	m_sName = "2RateHKY";
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pHKY,"HKY_proc(1)"));
	m_vpProc.push_back(m_vpProc[0]->RateProcessCopy());

	m_vpProc[0]->Rate(0.0);
	m_vpProc[1]->Rate(2.0);

	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->AddRatePar2Opt(); }

	FinalInitialisation();
}

C2RateJC::C2RateJC(CData *Data, CTree *Tree) : CBaseModel(Data,Tree) {
	int i;
	m_sName = "2RateJC";
	m_vpProc.push_back(AddDNAProcess(Data,Tree,pJC,"JC_proc(1)"));
	m_vpProc.push_back(m_vpProc[0]->RateProcessCopy());

	m_vpProc[0]->Rate(0.0);
	m_vpProc[1]->Rate(2.0);


	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->AddRatePar2Opt(); }

	FinalInitialisation();
}


//////////////////////////////////////////////
// Controllers for Amino acid models

CAAProcess *CBaseModel::AddAAProcess(CData *Data, CTree *Tree, AAProc Model, bool AddF)	{
	CAAProcess *Proc;
	Proc = new CAAProcess(Data,Tree,Model,AddF);
	return Proc;
}

CWAG::CWAG(CData *D, CTree *T, bool AddF) : CBaseModel(D,T)	{
	m_sName = sModelNames[(int)WAG]; if(AddF) { m_sName += "+F"; }
	m_vpProc.push_back(AddAAProcess(D,T,pWAG,AddF));
	FinalInitialisation();
}

CJTT::CJTT(CData *D, CTree *T, bool AddF) : CBaseModel(D,T) {
	m_sName = sModelNames[(int)JTT]; if(AddF) { m_sName += "+F"; }
	m_vpProc.push_back(AddAAProcess(D,T,pJTT,AddF));
	FinalInitialisation();
}

CDAY::CDAY(CData *D, CTree *T, bool AddF) : CBaseModel(D,T) {
	m_sName = sModelNames[(int)DAY]; if(AddF) { m_sName += "+F"; }
	m_vpProc.push_back(AddAAProcess(D,T,pDAY,AddF));
	FinalInitialisation();
}

CEQU::CEQU(CData *D, CTree *T, bool AddF) : CBaseModel(D,T)	{
	m_sName = sModelNames[(int)EQU]; if(AddF) { m_sName += "+F"; }
	m_vpProc.push_back(AddAAProcess(D,T,pEQU,AddF));
	FinalInitialisation();
}

CMTREV::CMTREV(CData *D, CTree *T, bool AddF) : CBaseModel(D,T) {
	m_sName = sModelNames[(int)mtREV]; if(AddF) { m_sName += "+F"; }
	m_vpProc.push_back(AddAAProcess(D,T,pMTREV,AddF));
	FinalInitialisation();
}

CEMP::CEMP(CData *D, CTree *T, string Name, bool AddF, double *S_ij, double *pi_j) : CBaseModel(D,T) {
	CAAProcess *Proc;
	m_sName = Name;
	if(AddF) { m_sName += "+F"; }
	Proc = new CAAProcess(D,T,m_sName,AddF, S_ij,pi_j);
	m_vpProc.push_back(Proc); Proc = NULL;
	FinalInitialisation();

}


////////////////////////////////////////////
// Controllers for Codon models

CCodonProcess *CBaseModel::AddCodonProcess(CData *D,CTree *T,CodonProc Model, ECodonEqm CE, int GenCode, string RadicalFile)	{
	CCodonProcess *Proc;
	Proc = new CCodonProcess(D,T,Model,CE,GenCode,RadicalFile);
	return Proc;
}

CCodonM0::CCodonM0(CData *D, CTree *T, ECodonEqm CE, int GenCode) : CBaseModel(D,T)	{
	m_sName = sModelNames[(int)CodonM0];
	// Do genetic code
	D->MakeCodonData();
	// Reduce the model and the data to the correct genetic code
	m_pData->ReduceCodonData(GenCode);

	m_sName += "." + int_to_string(GenCode) + ".";
	// Do frequencies
	switch(CE) {
	case cEQU: m_sName += "EQU";	break;
	case F1X4: m_sName += "F1X4";	break;
	case F3X4: m_sName += "F3X4";	break;
	case F64:  m_sName += "F64";	break;
	default:   Error("\nUnknown CE option in CCodonM0::CCodonM0\n\n");
	}
	// Add the process
	m_vpProc.push_back(AddCodonProcess(D,T,pM0,CE,GenCode));
	FinalInitialisation();
}

//////////////////////////////////////////////////////
// Controllers for Coevolution models

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Base coevolutionary model ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CBaseCoevo::CBaseCoevo(CData *JointData, CData *Dat1, CData *Dat2, CTree * Tree) : CBaseModel(JointData, Tree) {
	// Assign data objects
	m_pSingleData1 = Dat1;
	m_pSingleData2 = Dat2;
	m_pSingleProc1 = NULL; m_pSingleProc2 = NULL;

}



CBaseCoevo::~CBaseCoevo() {
	m_pSingleData1 = NULL; m_pSingleData2 = NULL;
	if(m_pSingleProc2 == m_pSingleProc1) { m_pSingleProc2 = NULL; } else { delete m_pSingleProc2; }
	delete m_pSingleProc1;
	cout << "\nNeed to fix destructor for CBaseCoevo"; exit(-1); }

double CBaseCoevo::lnL(bool ForceReal)	{
	int i;
	double Ind1, Ind2, Joint;
	// Do some basic error checking
	assert(m_pSingleData1 != NULL && m_pSingleData2 != NULL && m_pSingleProc1 != NULL && m_pSingleProc2 != NULL);
	cout << "\n <<<<<<<>>>>>> Into CBaseCoevo::lnL(...) <<<<<<<>>>>>>>";
	// Get the likelihoods for the independent models
	m_pSingleProc1->Likelihood(ForceReal);
	if(m_pSingleProc1 != m_pSingleProc2) { m_pSingleProc2->Likelihood(ForceReal); }

	cout << "\nSingle: " << m_pSingleProc1->Eqm(0);

	Ind1 = CalculateL(m_pSingleProc1,false);
	Ind2 = CalculateL(m_pSingleProc2,false);


	cout << "\nSingle process done: Ind1= " << Ind1 << " cf. Ind2= " << Ind2;
	cout << "\nNeed to create lnL for coevo model";
/*
	cout << "\n-----------------------------------------------------------------------";
	cout << "\nModel1: " << *m_pSingleProc1;
	cout << "\n-----------------------------------------------------------------------";
	cout << "\nTree: " << *m_pSingleProc1->Tree();
	cout << "\n-----------------------------------------------------------------------";
	cout << "\nData: " << *m_pSingleProc1->m_pData;
*/

	cout << "\nCalling joint model..." << flush;
	Joint = CBaseModel::lnL(ForceReal);

	cout << "\nThe Q matrix for the joint model:\n\n"; m_vpProc[0]->OutQ();
	cout << "\n\nThe eqm:\n" << m_vpProc[0]->Eqm(0);

	CProb pLeft,pRight,pProduct;
	cout << "\nJoint lnL: " << Joint << flush;
	cout << "\nMappingSize: " << m_pData->m_vviCoevoMapping.size() << flush;
	cout.precision(16);
	FOR(i,m_pData->m_iSize) {
		cout << "\nSite[" << i<<"]: Mapping[" << m_pData->m_vviCoevoMapping[i][0] <<"][" << m_pData->m_vviCoevoMapping[i][1] << "], PatOcc: " << m_pData->m_ariPatOcc[i];
		cout << "\tJoint: " << m_vpProc[0]->L(i).Prob();
		pLeft = m_pSingleProc1->L(m_pData->m_vviCoevoMapping[i][0]);
		pRight = m_pSingleProc2->L(m_pData->m_vviCoevoMapping[i][1]);
		pProduct = pLeft.Multiply(pRight,false);
		cout << "  Left: " << pLeft.Prob() << "  Right: " << pRight.Prob() << flush;
		cout << "  product: " << pProduct.Prob() << flush;
	}
	cout << "\n-------------------------------- Doing Aggregate -------------------------";
	m_vpProc[0]->AggregateQ(0);


	exit(-1);
}

double CBaseCoevo::DoBralnL(int B, int NL, int NR, bool JustClean) { // Do calculations for a single branch
	cout << "\nNeed to do CBaseCoevoDoBralnL function"; exit(-1);
}

/////////////////////////////////////////////////// CWAGCoevo /////////////////////////////////////////////////////////
CWAGCoevo::CWAGCoevo(CData *JointData,CData *Dat1, CData *Dat2, CTree *Tree) : CBaseCoevo(JointData,Dat1,Dat2,Tree) {
	// General definitions
	int InitPsi = 0.1;
	vector <double> PropMat(400,0);
	m_sName = sModelNames[(int)Coevo_WAG];
	// Create processes
//	m_pSingleProc1 = AddAAProcess(Dat1,Tree,pWAG,true);
	m_pSingleProc1 = AddAAProcess(Dat1,Tree,pEQU,true);
	if(Dat1 == Dat2) { m_pSingleProc2 = m_pSingleProc1; } else { m_pSingleProc2 = AddAAProcess(Dat2,Tree,pWAG,true); }
	cout << "\nOkay to here" << flush;
	AddAACoevoProcess(JointData,Dat1->m_vFreq, Dat2->m_vFreq, &PropMat,0.0);
	cout << "\ndone" << flush;
	cout << " Rates: " << Rates() << flush;
	cout << "\nDone Rates" << flush;
	FinalInitialisation();
}


CWAGCoevo::~CWAGCoevo()	{
	cout << "\nNeed to fix destructor for CWAGCoevo..."; exit(-1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function for get a model pointer
CBaseModel * GetMyModel(EModel ModelChoice, CData *Data, CTree *Tree)	{

	// DNA model stores
	CJC *JCModel = NULL;
	CFEL *FELModel = NULL;
	CK2P *K2PModel = NULL;
	CHKY *HKYModel = NULL;
	CREV *REVModel = NULL;
	// Amino acid models	(NB: Some of these aren't in the sModelNames string so it's a bit rough and ready)
	CEQU *EQUModel = NULL;
	CWAG *WAGModel = NULL;
	CJTT *JTTModel = NULL;
	CDAY *DAYModel = NULL;

	// Check data is of the appropriate form.
	if(ModelChoice >= JC && ModelChoice <= REV) { assert(Data->m_DataType == DNA); }
	else if(ModelChoice >= EQU && ModelChoice <= mtREV) { assert(Data->m_DataType == AA); }

	switch(ModelChoice) {
	case JC:
		JCModel = new CJC(Data,Tree);
		return JCModel;
	case FEL:
		FELModel = new CFEL(Data,Tree);
		return FELModel;
	case K2P:
		K2PModel = new CK2P(Data,Tree);
		return K2PModel;
	case HKY:
		HKYModel = new CHKY(Data,Tree);
		return HKYModel;
	case REV:
		REVModel = new CREV(Data,Tree);
		return REVModel;
	case EQU:
		EQUModel = new CEQU(Data,Tree);
		return EQUModel;
	case WAG:
		WAGModel = new CWAG(Data,Tree);
		return WAGModel;
	case JTT:
		JTTModel = new CJTT(Data,Tree);
		return JTTModel;
	case DAY:
		DAYModel = new CDAY(Data,Tree);
		return DAYModel;
	default:
		cout << "\nYou've asked for a model that I haven't implemented in GetMyModel(...) yet: " << ModelChoice << ". Bad luck :(\n\n"; exit(-1);
	}

	Error("Reached impossible point in GetMyModel(...)");
	return NULL;
/*
	const string sModelNames[] = {	"JC","FEL","K2P","HKY","REV",
							"RY","CovHKY",",CovREV","THMM_FULLDNA","HKYdG_THMM","THMM_DNA",
							"EQU","WAG","JTT","DAY","mtREV",
							"THMM_AA","WAGdG_THMM",
							"Coevo_WAG",
							"CodonM0",
							"Unknown"};
*/

}




