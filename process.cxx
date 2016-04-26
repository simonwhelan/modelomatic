////////////////////////////////////////////////////////////////
// Implementation of the process.h classes
// ---------------------------------------

#include "process.h"

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter, MakeQ_Log_Counter, MakePT_Log_Counter, LFunc_Log_Counter, SubLFunc_Log_Counter;
#endif

#if DO_MEMORY_CHECK
extern CMemChecker memory_check;
#endif

#define ANALYTIC_DERIVATIVE_DEBUG 0	// Checker for the analytical derivative functions
#define ADD_SITE_MAX 3			// Maximum site for doing analytic derivative debug
#define HARDCHECK_CALCS 0			// Hard check calculations for errors (not complete)

////////////////////////////////////////////////////////////////////////////////
// Data compression routines
////////////////////////////////
// Returns a vector <int> of size (Data->m_iSize * #Internal nodes), which is mapped so [Node#][Site];
// An entry of -1 indicates a real node. An entry of >=0 is a reference to the Site# it matches
vector <int> GetCompressedData(CTree *T,CData *D)	{
	vector <vector <int> > ColInfo;
	T->BestStartCalc();
	if(D == NULL) { Error("\nNeed proper data for GetCompressedData...\n"); }
	// Get the data columns
	ColSorter(T->StartCalc(),-1,T,D,ColInfo);
	// Process them to get the mappings and return them
	return DoCompress(ColInfo, T, D);
}

// Column sorter
// -------------
// Post-order tree traversal to do column sorting
void ColSorter(int NT, int NF, CTree *T, CData *D, vector <vector <int> > &ColInfo,bool First)	{
	int i,j,site,Node,From;
	// Initialise
	if(First == true) {
		vector <int> V;
		ColInfo.assign( (T->NoNode() - T->NoSeq()) * D->m_iSize, V);
	}
	// Traversal
	FOR(i,T->NoLinks(NT))	{
		if(T->NodeLink(NT,i) == NF || T->NodeLink(NT,i) == -1) { continue; }
		ColSorter(T->NodeLink(NT,i),NT,T,D,ColInfo,false);
	}
	// Do the processing
	Node = NF - T->NoSeq();
	From = NT - T->NoSeq();
	if(Node < 0) { return; } // Don't compress external nodes
	// 1. Get sequences from
	if(NT < T->NoSeq()) {
		for(i=0,site = Node * D->m_iSize; site < (Node+1) * D->m_iSize; site++,i++) { ColInfo[site].push_back(D->m_ariSeq[NT][i]); }
	} else {
		assert(InRange(From,0,T->NoNode()-T->NoSeq()));
		for(i= From * D->m_iSize,site = Node * D->m_iSize; site < (Node+1) * D->m_iSize; site++,i++) {
			FOR(j,(int)ColInfo[i].size()) { ColInfo[site].push_back(ColInfo[i][j]); }
}	}	}

///////////////////////////////////////////////////
// Do the compression
vector <int> DoCompress(vector <vector <int> > &CV, CTree *T, CData *D)	{
	int Node,i,j;
	vector <int> Compress;
	Compress.assign(D->m_iSize * (T->NoNode()-T->NoSeq()),-1);
	FOR(Node,T->NoNode()-T->NoSeq())	{	// Loop through nodes
//		cout << "\nDoing node: " << Node + T->NoSeq();
		for(i = D->m_iSize * Node; i < D->m_iSize * (Node + 1); i++) {
			if(Compress[i] >= 0) { continue; } // If its not real then skip it
			for(j = i+1; j < D->m_iSize * (Node + 1); j++) {
				if(Compress[j] >= 0) { continue; }
//				cout << "\nComparing CV["<<i<<"] = " << CV[i] << " cf. CV["<<j<<"] = " << CV[j];
				if(CV[i] == CV[j]) { Compress[j] = i; }
	}	}	}
	return Compress;
}

///////////////////////////////////////////////////////////////////////////////
// Q matrix parameter definition
///////////////////////////////////////////////////////////////////////////////

CQPar::CQPar(string Name, int Char, double Value, bool Optimise, double Lowbound, double Upbound,ParOp Oper) : CPar(Name,Value,Optimise,Lowbound,Upbound,Oper) {
	m_iChar = Char;
}

// Functions allocating the Q_ij elements
void CQPar::AddQij(int i, int j, bool Rev)	{
	m_viQMap.push_back((i*m_iChar)+j);
	if(Rev == true) { m_viQMap.push_back((j*m_iChar)+i); }
}

void CQPar::Par2Q(double *QMat,int QMatID)	{
	int i;
	if(QMatID == -1 || find(m_viApply2Q.begin(),m_viApply2Q.end(),QMatID) != m_viApply2Q.end() || m_viApply2Q.empty()) {
		FOR(i,(int)m_viQMap.size()) {
			assert(m_viQMap[i] < m_iChar * m_iChar);
			QMat[m_viQMap[i]] = DoOper(QMat[m_viQMap[i]]);
	}	}
}

// Q based Parameter outputter
ostream &CQPar::Output(ostream &os)	{
	int i;
	string delim = " ";
	CPar::Output(os);
	if(!m_viQMap.empty() && m_bOutDetail == true)	{
		os << delim << "QMat";
		FOR(i,(int)m_viQMap.size()) { os << delim << m_viQMap[i]; }
	}
	return os;
}

///////////////////////////////////////////////////////////////////////////////
// Alpha parameter of gamma distribution definition
///////////////////////////////////////////////////////////////////////////////

CGammaPar::CGammaPar(string Name, int Char, CPar * Rate, double Value, bool Optimise, double Lowbound, double Upbound,ParOp Oper) : CQPar(Name,Char,Value,Optimise,Lowbound,Upbound,Oper)	{
	m_iNoCat = 0; m_pProb = NULL; m_pProcRate = Rate;
}

CGammaPar::~CGammaPar()	{
	int i;
	FOR(i,m_iNoCat) { m_arpRates[i] = NULL; } m_arpRates.clear();
	m_pProb = NULL;
	m_pProcRate = NULL;
}

void CGammaPar::AddRateToGamma(CBaseProcess *Proc)	{
	// Check entry conditions
	assert(m_arpRates.size() == m_iNoCat);
	// If this is the first process store the probability
	if(m_iNoCat == 0) { m_pProb = Proc->ProbPar(); }
	m_arpRates.push_back(Proc->RatePar());
	Proc->SetProbFactor(++m_iNoCat);
}

ostream &CGammaPar::Output(ostream &os)	{
	int i;
	CPar::Output(os);
	// Check some entry conditions
	assert(m_arpPar.empty() && m_arpRates.size() == NoCat());
	os << "\n\t" << NoCat() << " Rate categories; TotalRate: " << m_pProcRate->Val();
	GlobalApply();
	FOR(i,NoCat())	{ os << "\n\tRate["<<i<<"]: " << m_arpRates[i]->Val() << "; Prob: " << m_pProb->Val() / (double) m_iNoCat; }
	return os;
}

void CGammaPar::GlobalApply()	{
	int i;
	double *P,*R,Rate=m_pProcRate->Val();
	if(Val() > MAX_ALFA) { CheckBound(true); FOR(i,NoCat())  { m_arpRates[i]->SetVal(1.0); } return; }
	GET_MEM(P,double,NoCat()); GET_MEM(R,double,NoCat());
	// Check some entry conditions
	assert(m_arpPar.empty() && m_arpRates.size() == NoCat());
	// Get the gamma distributed rates
	DiscreteGamma(P,R,Val(),Val(),NoCat(),0);
//	cout << "\nRate: " << Rate;

	FOR(i,NoCat())	{ R[i] *= Rate; m_arpRates[i]->SetVal(R[i]); }
//	cout << "\nGetting gamma -- alpha ["<<m_arpRates.size() << ":"<<NoCat()<<"] = " << Val() << " == { "; FOR(i,NoCat()) { cout << R[i] << " "; } cout << "}";
	// Tidy up
	DEL_MEM(P); DEL_MEM(R);
}

//////////////////////////////////////////////////////////////////////////
// This is the Q matrix implementation
/////////////////////////////////////////////////////////////////////////
// This contains the space for a single Q matrix and the functions
// for turning it into a P(t) matrix

// Constructor functions
CQMat::CQMat(EDataType Type,string Name)			{
#if DO_MEMORY_CHECK
	memory_check.CountCQMat++;
#endif
	MakeSpace(Type,Name); }
CQMat::CQMat(int Char, EDataType Type,string Name)	{
#if DO_MEMORY_CHECK
	memory_check.CountCQMat++;
#endif
	MakeSpace(Char,Type,Name); }

CQMat::~CQMat()	{
#if DO_MEMORY_CHECK
	memory_check.CountCQMat--;
#endif
	CleanSpace();
}

// Memory management
void CQMat::MakeSpace(int Char,EDataType Type,string Name)	{
	m_iQMatID = GetQMatID();
	m_sName = Name;
	m_bScaleReady = false;
	m_ardQMat = NULL; m_ardU = NULL; m_ardV = NULL; m_ardRoot = NULL;
	m_iChar = Char; m_iChar2 = Char*Char;
	m_DataType = Type;
	GET_MEM(m_ardQMat,double,m_iChar2);
	GET_MEM(m_ardU,double,m_iChar2);
	GET_MEM(m_ardV,double,m_iChar2);
	GET_MEM(m_ardRoot,double,m_iChar);
	GET_MEM(m_ardEqm,double,m_iChar);
	GET_MEM(m_ardRootEqm,double,m_iChar);
	m_bAlwaysI = false;
	Unlock();	// All everything to be changeable
}
void CQMat::CleanSpace()	{
	DEL_MEM(m_ardQMat); DEL_MEM(m_ardU); DEL_MEM(m_ardV); DEL_MEM(m_ardRoot); DEL_MEM(m_ardEqm); DEL_MEM(m_ardRootEqm); m_bAlwaysI = false;
}

// Access functions
vector <double> CQMat::Eqm()	{
	vector <double> Vec(&m_ardEqm[0],&m_ardEqm[m_iChar]);
	return Vec;
}

// Functions that modify the Q matrix
// Initialise the matrix to a particular value
void CQMat::InitQ(double Val)	{ if(m_bAllowModelUpdate == false) { return; } int i; FOR(i,m_iChar2) { m_ardQMat[i] = Val; } m_bScaleReady = false; if(Val > DBL_EPSILON) { m_bAlwaysI = false; } else { m_bAlwaysI = true; } }
// Applies a parameter to the matrix
void CQMat::ApplyPar2Q(CQPar *Par)	{ if(m_bAllowModelUpdate == false) { return; } /* cout << "\nApplying parameter: " << *Par; */ Par->Par2Q(m_ardQMat,m_iQMatID); }
// Make the correct diagonal
void CQMat::DoQDiag()				{
	if(m_bAllowModelUpdate == false) { return; }
	double Total,Rate = 0.0;
	int i,j;
	FOR(i,m_iChar)	{
		Total = 0;
		FOR(j,m_iChar)	{
			if(i==j) { continue; }
			if(*Q(i,j) < 0)	{
				if(*Q(i,j) < -1.0E-5) {
					cout << "\nOdd Q matrix, with Q["<<i<<","<<j<<"]: " << *Q(i,j);
					OutQ(cout);
					exit(-1);
				} else { *Q(i,j) = 0; }
			}
			Total += *Q(i,j);
		}
		*Q(i,i) = -Total;
}	}
// Get the overall rate of the process
double CQMat::OverallSubRate(bool FWP)	{
	int i;
	double Total = 0.0;
	assert(m_bScaleReady == true);
	Total = 0;
	FOR(i,m_iChar) { Total -= *Q(i,i) * m_ardEqm[i]; }
	return Total;
}
// Scale the Q matrix and the eigen roots to give the process the correct rate.
void CQMat::ScaleQ(double Rate)	{
	m_bAlwaysI = false;
	// Deal with zero rates
	assert(Rate >= 0.0);
//	cout << "\nScaling Q by " << Rate<<": "; OutQ();
	if(Rate < 1.0E-6 || m_dScale < DBL_EPSILON || fabs(Rate - MAX_PAR_VALUE) < DX) {
		if(OverallSubRate() > DBL_EPSILON && m_dScale < DBL_EPSILON) { Error("\nUnexpected values in CQMat::ScaleQ.\n"); }
		m_bAlwaysI = true;  // This is true even for covarion processes -- cannot have zero substitutions with changes of state!
		return;
	}
	assert(m_bScaleReady == true);
	int i;
	double NewRate = Rate / m_dScale;
	FOR(i,m_iChar)	{ m_ardRoot[i] *= NewRate; }
	FOR(i,m_iChar2)	{ m_ardQMat[i] *= NewRate; }
//	cout << "\nScaled Q: "; OutQ();
//	cout << "\nNew rate: " << OverallRate() << " cf. " << OverallRate(true);
	m_dScale = Rate;
}
// Decomposition and scaling factor
void CQMat::Decompose(double *Eqm,bool DoScale, bool REV,double Rate)	{
	if(m_bAllowModelUpdate == false) { return; }
	int i;
	double *Mat;
	GET_MEM(Mat,double,m_iChar2);	FOR(i,m_iChar2) { Mat[i] = m_ardQMat[i]; }
	// Do slow eigen decomposition if the process in non-reversible or the eqm not known
	if(Eqm == NULL || REV == false) {
		Error("Haven't done decomposition for CQMat::Decompose (Eqm == NULL || REV == false)");
	} else {
	// Otherwise do the quick eigen decomposition
/*		cout << "\nDecomposing with Eqm: ";
		int j;
		FOR(i,5) {
			double Total = 0.0;
//			cout << "\n\tEqm["<<i<<"]: ";
			FOR(j,20) {
				Total += Eqm[(i*20)+j]; // cout << Eqm[(i*20)+j] << " ";
			}
			cout << "  Sum["<<i<<"]: " << Total;
		}
*/
//		cout << "\nDoing eigen-decompose with eqm[" << m_iChar << "]:"; FOR(i,m_iChar) { cout << " " << Eqm[i]; }
		FOR(i,m_iChar) { m_ardEqm[i] = Eqm[i]; m_ardRootEqm[i] = sqrt(Eqm[i]); }
#if HARDCHECK_CALCS == 1
		if(!CheckReversibility(m_iChar,Q(),CQMat::Eqm())) {
			cout << "\nDetailed balance failed...\nQ mat:";  OutQ(); cout << "\n\neqm: " << CQMat::Eqm(); exit(-1); }
#endif

		eigenQREV(Mat,m_ardEqm,m_ardRootEqm,m_iChar,0,m_ardRoot,m_ardU,m_ardV);
	}
	m_bScaleReady = true;
	m_dScale = OverallSubRate();
	// If required scale to rate 1.0
	if(DoScale == true)	{ assert(Rate>=0.0); ScaleQ(Rate); }
	DEL_MEM(Mat);

}
void CQMat::Decompose(vector <double> Eqm,bool Scale,bool REV,double Rate)	{
	if(m_bAllowModelUpdate == false) { return; }
	int i;
	double *eqm;
	GET_MEM(eqm,double,m_iChar);	FOR(i,m_iChar) { eqm[i] = Eqm[i]; }
	Decompose(eqm,Scale,REV,Rate);
	DEL_MEM(eqm);
}
// Makes the P(t)
bool CQMat::MakePT(double T, double *PT)	{
//	cout << "\nMaking P(T="<<T<<")";
   int i,j,k,Char2 = m_iChar * m_iChar;
    double e1, e2, *P;
	const double *pdV,*pdU;
	if(T < -FLT_EPSILON) { T = 0; }
	if(T < DX || m_bAlwaysI)		{ // Small branch lengths create identity matrix regardless of reversibility
		IMat(PT,m_iChar);
		return true;
	}
	else if(T > MAX_BRANCH)	{
		FOR(i,m_iChar) { FOR(j,m_iChar) { PT[(i*m_iChar)+j] = m_ardEqm[j]; } }
		return true;
	}
    P = PT; for(i=0;i<Char2;i++) { *(P++)=0; }
    for(k=0;k<m_iChar;k++)		{
        P=PT; 	pdU = &m_ardU[k];// Set P for pointer arithmatic
        e1=exp(T*m_ardRoot[k]);
        for(i=0;i<m_iChar;i++)	{
            e2=*pdU*e1; pdV = &m_ardV[(k*m_iChar)]; pdU += m_iChar;
            for(j=0;j<m_iChar;j++)	{ *P++ += (e2 * *(pdV++)); }

	}	}
#if FUNC_COUNTERS
	MakePT_Log_Counter++;
#endif
#if HARDCHECK_CALCS == 1
	k = 0; bool okay = true;
	FOR(i,m_iChar) {
		e1 = 0.0;
		FOR(j,m_iChar) { if(PT[k] < -FLT_EPSILON || my_isnan(PT[k]) ) { okay = false; } e1 += PT[k++]; }
		if(tdiff(e1,1.0,1.0E-6)) { cout << "\nRow " << i << " failed == " << e1 << "; diff= " << fabs(1.0 - e1); okay = false; }
	}
	if(okay == false) {
		cout << "\nFailed CQMat::MakePT()";
		cout << "\nP("<<T<<"); overall rate: " << OverallSubRate(); MatOut(m_iChar,PT); cout << "\n\nThe QMat: "; OutQ();
		cout << "\nEqm: " << Eqm() << " == " << Sum(Eqm());
		if(!CheckReversibility(m_iChar,Q(),Eqm())) { cout << "\nDetailed balance failed...";  }
		cout << "\nThe U mat: \n"; OutU();
		cout << "\nThe U mat: \t"; OutV();
		cout << "\nThe root vector: \t"; OutRoot();
		exit(-1);
	}
#endif
//	cout << "\n\nnew P("<<T<<"; Rate = " << m_dScale<<"):"; FOR(i,m_iChar*m_iChar) { if(i%m_iChar==0) { cout << endl; } cout << PT[i] << "\t"; } cout << "\n";
//	exit(-1);
	return true;
}
// Output
ostream &CQMat::Output(ostream &os)	{
	int i;
	os << "Q matrix (ID=" << m_iQMatID << "): " << m_sName;
	FOR(i,m_iChar2) { if(i%m_iChar==0) { os << "\n\t"; } os << m_ardQMat[i] << "\t"; }
	return os;
}

///////////////////////////////////////////////////////////////////////////////
// Equilibrium class implementations
///////////////////////////////////////////////////////////////////////////////

///////////////////////////// CBaseEqm class //////////////////////////////////
CBaseEqm::CBaseEqm(int Char, vector <CQPar *> *ProcPar)	{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseEqm++;
#endif
	m_iChar = Char;
	GET_MEM(m_ardQModifier,double, m_iChar*m_iChar);
	m_pProcPar = ProcPar;
}
CBaseEqm::~CBaseEqm()			{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseEqm--;
#endif
	int i;
	DEL_MEM(m_ardQModifier);
	FOR(i,(int)m_vpEqmPar.size()) { m_vpEqmPar[i] = NULL; } m_vpEqmPar.clear();
	m_pProcPar = NULL;
}
// Add Q matrix ID for the eqm distribution to apply to
void CBaseEqm::AddMatID(int ID)	{ m_viQMatID.push_back(ID); }
// Whether the matrix of ID has this eqm distribution applied to it
bool CBaseEqm::IsID(int ID)	{
	if(find(m_viQMatID.begin(),m_viQMatID.end(),ID) != m_viQMatID.end() || m_viQMatID.empty() ) { return true; }
	return false;
}

// Whether to optimise the equilibrium distribution
void CBaseEqm::SetOpt(bool Opt)  { int i; FOR(i,(int)m_vpEqmPar.size()) { m_vpEqmPar[i]->SetOptimise(Opt); } }

// The basic application function that applies the eqm distribution (in the form of m_ardQModifier) to
//  to the CQMat object
void CBaseEqm::ApplyEqm2QMat(double *Q, int MatID)	{
	int i;
	// If the matrices ID is in m_viQMatID or if m_viQMatID is empty then apply the Q matrix
	if(IsID(MatID))	{
		if(!CheckQ()) { AdapterFunc(); }		// Create m_ardQModifier using AdapterFunc() if the parameters have changed
		FOR(i,m_iChar*m_iChar) { Q[i] *= m_ardQModifier[i]; }
}	}

// Check as to whether m_ardQModifier needs changing
bool CBaseEqm::CheckQ()	{
	int i;
	// If hasn't been done at all yet
	if(m_vdOldParVals.empty()) {
		FOR(i,(int)m_vpEqmPar.size()) { m_vdOldParVals.push_back(m_vpEqmPar[i]->Val()); }
		return false;
	}
	// Otherwise do the obvious check
	assert(m_vpEqmPar.size() == m_vdOldParVals.size());
	FOR(i,(int)m_vpEqmPar.size()) { if(fabs(m_vpEqmPar[i]->Val() - m_vdOldParVals[i]) > DBL_EPSILON) { break; } }
	if(i == m_vpEqmPar.size()) { return true; }
	FOR(i,(int)m_vpEqmPar.size()) { m_vdOldParVals[i] = m_vpEqmPar[i]->Val(); }
	return false;

}

vector <double *> CBaseEqm::OptimiserValues()	{
	int i;
	vector <double *> Vec;
	FOR(i,(int)m_vpEqmPar.size()) { Vec.push_back(m_vpEqmPar[i]->OptimiserValue()); }
	return Vec;
}


///////////////////////////// CSimpleEqm class //////////////////////////////////
CSimpleEqm::CSimpleEqm(int Char, vector <CQPar *> *ProcPar, vector <double> eqm)	: CBaseEqm(Char,ProcPar) {
	ConstructSimpleEqm(eqm);
}
CSimpleEqm::CSimpleEqm(int Char, vector <CQPar *> *ProcPar, CData *Data)	: CBaseEqm(Char,ProcPar)			{
	assert(Char == Data->m_iChar && Data->m_vFreq.size() == Char);
	ConstructSimpleEqm(Data->m_vFreq);
}
CSimpleEqm::~CSimpleEqm()	{ /* Currently empty */ }

void CSimpleEqm::AdapterFunc()	{
	int i,j;
	assert(m_vpEqmPar.size() == m_iChar);
	FOR(i,m_iChar) {
		FOR(j,m_iChar)	{
			m_ardQModifier[(i*m_iChar)+j] = m_vpEqmPar[j]->Val();
	}	}
//	FOR(i,m_iChar*m_iChar) { if(i%m_iChar == 0) { cout << endl; } cout << m_ardQModifier[i] << "\t"; }
}

void CSimpleEqm::ConstructSimpleEqm(vector <double> eqm)	{
	int i;
	EDataType Type;
	CQPar *Par;
	string Name;
	// Check entry conditions
//	cout << "\nConstructSimpleEqm [eqm.size()=" << eqm.size() << ", m_iChar=" << m_iChar << "]";
	assert(eqm.size() == m_iChar);
	assert(fabs(Sum(&eqm) - 1.0) < FLT_EPSILON);
	switch(eqm.size())	{
	case 2:		Type = RY; break;
	case 4:		Type = DNA; break;
	case 20:	Type = AA;	break;
	case 64:	Type = COD;	break;
	case 400:	Type = AA2; break;
	case 16: 	Type = DNA2; break;
	default: Error("Unknown type of eqm data");
	};
	// Check eqm distribution parameters are sensible (i.e. ban zeros)
	FixProb(&eqm);
	// Add the eqm parameters
	FOR(i,(int)eqm.size())	{
		Name = "Freq(" + State(Type,i) + ")";
		Par = new CQPar(Name,m_iChar,eqm[i],true,MIN_PROB,1.0,MULTIPLY);
		m_vpEqmPar.push_back(Par);
		m_pProcPar->push_back(Par);
		Par = NULL;
	}
	ProbabilityScale(&m_vpEqmPar,true,true,true);
}

vector <double> CSimpleEqm::Eqm() {
	vector <double> eq;
	int i;
	FOR(i,m_iChar) { eq.push_back(m_vpEqmPar[i]->Val()); }
	return eq;
}

void CSimpleEqm::ResetEqm(vector <double> New, bool RandomBit)	{
	int i,j = -1;
	if((int)New.size() != m_iChar) { Error("Error: CSimpleEqm::ResetEqm -- vector <double> New has wrong number of states...\n"); }
	New = NormaliseVector(New);
//	cout << "\nResetEqm: ";
	FOR(i,m_iChar) {
		if(m_vpEqmPar[i]->Special()) { assert(j == -1); j = i; continue; }
		if(!RandomBit) { m_vpEqmPar[i]->SetVal(New[i],false,true,false);  }
		else { m_vpEqmPar[i]->SetVal(New[i] + (0.1 * Random()),false,true,false); }
//		cout << "\n\t["<<i<<"]: " << New[i]; // << " --> " << m_vpEqmPar[i]->Val();

	}
	assert(j != -1);
	m_vpEqmPar[j]->SetVal(New[j],true,true,true);
}

/////////////////////////////////////////////////////////////////////////////
// CCodonEqm class

// Constructor
// How the equilibrium is built is highly data dependent.
CCodonEqm::CCodonEqm(ECodonEqm CE, int GenCode, std::vector<CQPar*> *ProcPar, CData *D) : CBaseEqm(D->m_iChar, ProcPar)	{
	int i,j,k;
	string cod_temp, temp,Name;
	CQPar *Par = NULL;
	vector <double> Frqs, tmpFrqs;
	vector <CPar *> VecPar;
	assert(D!=NULL);
	assert((D->m_iChar == 64 && D->m_DataType == COD) || D->m_DataType == COD_RED);
	m_CE = CE;
	// Function to create the different types of eqm
	switch(CE)	{
	case cEQU:
		break;
	case F1X4:
		Frqs.assign(4,0.0);
		FOR(i,D->m_iSize) {
			FOR(j,D->m_iNoSeq) {
				if(D->m_ariSeq[j][i] == D->m_iChar) { FOR(k,4) { Frqs[k] += 0.75 * (double) D->m_ariPatOcc[i]; } continue; } // Gaps treated as missing data
				cod_temp = D->m_sABET.substr(D->m_ariSeq[j][i]*3,3);
				FOR(k,3) {  temp = cod_temp[k]; Frqs[FindState(DNA,temp)]+= D->m_ariPatOcc[i]; }
		}	}
		Frqs = NormaliseVector(Frqs);
		// Add the eqm parameters
		FOR(i,4)	{
			Name = "Freq(" + State(DNA,i) + ")";
			Par = new CQPar(Name,m_iChar,Frqs[i],true,MIN_PROB,1.0,MULTIPLY);
			m_vpEqmPar.push_back(Par);
			m_pProcPar->push_back(Par);
			Par = NULL;
		}
		ProbabilityScale(&m_vpEqmPar,true,true,true);
		break;
	case F3X4:
		// Count the frequencies of nucleotides
		Frqs.assign(12,0.0);
		FOR(i,D->m_iSize) {
			FOR(j,D->m_iNoSeq) {
				if(D->m_ariSeq[j][i] == D->m_iChar) {
					FOR(k,12) { Frqs[k] += 0.25 * (double) D->m_ariPatOcc[i]; } continue; } // Gaps treated as missing data
				cod_temp = D->m_sABET.substr(D->m_ariSeq[j][i]*3,3);
				FOR(k,3) { temp = cod_temp[k]; Frqs[(k*4) + FindState(DNA,temp)]+= D->m_ariPatOcc[i]; }
		}	}

		// Add the eqm parameters
		tmpFrqs.assign(4,0);
		FOR(j,3) {
			FOR(i,4) { tmpFrqs[i] = Frqs[(j*4)+i]; } tmpFrqs = NormaliseVector(tmpFrqs);
			FOR(i,4)	{
				Name = "Freq(" + State(DNA,i) + "[" + int_to_string(j)+"])";
				Par = new CQPar(Name,m_iChar,tmpFrqs[i],true,MIN_PROB,1.0,MULTIPLY);
				m_vpEqmPar.push_back(Par);
				m_pProcPar->push_back(Par);
				VecPar.push_back(Par);
				Par = NULL;
			}
			ProbabilityScale(&VecPar,true,true,true);
			FOR(i,4) { VecPar[i] = NULL; } VecPar.clear();
		}
		break;
	case F64: // Or equivalent for other genetic codes...
		Frqs = NormaliseVector(D->m_vFreq);
		FOR(i,(int)Frqs.size()) {
			if(GenCodes[GenCode][i] == -1 && D->m_iChar == 64) { if(Frqs[i] > DBL_EPSILON) { Error("\nSequences have stop codon " + State(COD,i) + " in for that genetic code...\n"); } continue; }
			Name = "Freq(" + D->m_sABET.substr(i*3,3) +")\0";
			assert(InRange(Frqs[i],0.0,1.0));
			Par = new CQPar(Name,m_iChar,Frqs[i],true,MIN_PROB,1.0,MULTIPLY);
			m_vpEqmPar.push_back(Par);
			m_pProcPar->push_back(Par);
			Par = NULL;
		}
		ProbabilityScale(&m_vpEqmPar,true,true,true);
		break;
	default:
		Error("\nTrying to create unknown type of equilibrium in CCodonEqm::CCodonEqm(...)\n");
	};
	// Set the data pointer and genetic code
	m_pData = D;
	m_iGenCode = GenCode;
}

// Destructor
CCodonEqm::~CCodonEqm() {
	m_pData = NULL;
};

// Adapter function
void CCodonEqm::AdapterFunc()	{
	int i,j;
	vector <double> veqm = Eqm();
	assert((int)veqm.size()  == m_iChar);
	FOR(i,m_iChar) {
		FOR(j,m_iChar)	{
			m_ardQModifier[(i*m_iChar)+j] = veqm[j];
}	}	}

// Eqm function
vector <double> CCodonEqm::Eqm()	{
	int i,j;
	bool SkipStops = false;
	vector <double> retEqm,Frqs;
	vector <double>::iterator iFreq;
	string s,c;
	// Some initialisation and error checking
	if(m_pData != NULL) {
		if(m_pData->m_DataType == COD_RED) { if(m_pData->GenCode() != m_iGenCode) { Error("\nGenCode in data and CCodonEqm::Eqm() do not match\n"); } }
		if((int)retEqm.size() != m_pData->m_iChar) { retEqm.clear(); retEqm.assign(m_pData->m_iChar,1.0); }
	}
	assert(InRange(m_iGenCode,0,12));
	// Calculate the eqm distributions
	switch(m_CE) {
	case cEQU:
		retEqm.assign(m_pData->m_iChar,1.0/(double)m_pData->m_iChar);
		break;
	case F1X4:
		if(retEqm.empty()) { retEqm.assign(m_iChar,1.0); }
		// Count the frequencies of nucleotides
		if(Frqs.empty()) { Frqs.assign(4,1.0); }
		assert((int)m_vpEqmPar.size() == 4 && Frqs.size() == 4);
		FOR(i,4) { Frqs[i] = m_vpEqmPar[i]->Val(); }
		FOR(i,m_iChar) {
			if(i > (int) retEqm.size()) { Error("\nMix-up in CCodon::Eqm() with pos > retEqm.size()\n\n"); }
			if(m_iChar == 64 && GenCodes[m_iGenCode][i] == -1) { if(!SkipStops) { retEqm[i] = 0.0; } continue; }
			c = m_pData->m_sABET.substr(i*3,3);
			FOR(j,3) { s = c[j]; retEqm[i] *= Frqs[FindState(DNA,s)]; }
		}
		retEqm = NormaliseVector(retEqm);
		break;
	case F3X4:
		if(retEqm.empty()) { retEqm.assign(m_iChar,1.0); }
		// Count the frequencies of nucleotides
		if(Frqs.empty()) { Frqs.assign(12,0.0); }
		assert((int)m_vpEqmPar.size() == 12 && Frqs.size() == 12);
		FOR(i,12) { Frqs[i] = m_vpEqmPar[i]->Val(); }
		FOR(i,m_iChar) {
			if(i > (int) retEqm.size()) { Error("\nMix-up in CCodon::Eqm() with pos > retEqm.size()\n\n"); }
			if(m_iChar == 64 && GenCodes[m_iGenCode][i] == -1) { if(!SkipStops) { retEqm[i] = 0.0; } continue; }
			c = m_pData->m_sABET.substr(i*3,3);
			FOR(j,3) { s = c[j]; retEqm[i] *= Frqs[(j*4) + FindState(DNA,s)]; }
		}
		retEqm = NormaliseVector(retEqm);
		break;
	case F64:
		if(retEqm.empty()) { retEqm.assign(m_iChar,1.0); } assert((int)m_vpEqmPar.size() == m_iChar);
		FOR(i,m_iChar) {
			if(i > (int) retEqm.size()) { Error("\nMix-up in CCodon::Eqm() with pos > retEqm.size()\n\n"); }
			if(m_iChar == 64 && GenCodes[m_iGenCode][i] == -1) { if(!SkipStops) { retEqm[i] = 0.0; } continue; }
			retEqm[i] = m_vpEqmPar[i]->Val();
		}
		break;
	}
	return retEqm;
}

void CCodonEqm::ResetEqm(vector <double> New, bool RandomBit)	{
	int i,j = -1;
	if((int)New.size() != m_iChar) { Error("Error: CSimpleEqm::ResetEqm -- vector <double> New has wrong number of states...\n"); }
	New = NormaliseVector(New);
//	cout << "\nResetEqm: ";
	FOR(i,m_iChar) {
		if(m_vpEqmPar[i]->Special()) { assert(j == -1); j = i; continue; }
		if(!RandomBit) { m_vpEqmPar[i]->SetVal(New[i],false,true,false);  }
		else { m_vpEqmPar[i]->SetVal(New[i] + (0.1 * Random()),false,true,false); }
//		cout << "\n\t["<<i<<"]: " << New[i]; // << " --> " << m_vpEqmPar[i]->Val();

	}
	assert(j != -1);
	m_vpEqmPar[j]->SetVal(New[j],true,true,true);
}


/////////////////// Function to add a codon Equilibrium //////////////////////////
void CBaseProcess::AddCodonEqm(int GenCode,int Char, ECodonEqm CE, bool Opt)	{
	CCodonEqm *Eqm;
	if(m_pData == NULL) { Error("\nCannot create CodonEqm without data..."); }
	Eqm = new CCodonEqm(CE,GenCode,&m_vpPar,m_pData);
	Eqm->SetOpt(Opt);
	m_vpEqm.push_back(Eqm);
	Eqm = NULL;
}

///////////////////////////////////////////////////////////////////////////////
// CSite definition
//////// ///////////////////////////////////////////////////////////////////////

CSite::CSite(int *Char)	{
#if DO_MEMORY_CHECK
	memory_check.CountCSite++;
#endif
	m_Char = Char;
	GET_MEM(m_ardSpace,double,*m_Char);
	m_bReal = true; m_iOriSite = -1; m_iScale = 0;
	m_iCopyNum=1; m_pOrSite = NULL;
	m_pSpacePointer = m_ardSpace;
	m_pScalePointer = &m_iScale;
}
CSite::CSite(const CSite &Site)	{
#if DO_MEMORY_CHECK
	memory_check.CountCSite++;
#endif
	int i;
	m_Char = Site.m_Char;
	assert(Site.m_bReal == true && *m_Char != -1);
	GET_MEM(m_ardSpace,double,*m_Char);
	FOR(i,*m_Char) { m_ardSpace[i] = 0.0; }
	m_bReal = true; m_iOriSite = -1; m_iScale = 0;
	m_iCopyNum=1; m_pOrSite = NULL;
	m_pSpacePointer = m_ardSpace;
	m_pScalePointer = &m_iScale;
}
CSite::~CSite()			{
#if DO_MEMORY_CHECK
	memory_check.CountCSite--;
#endif
	ResetSite();
	if(m_bReal == true) { DEL_MEM(m_ardSpace);  }
	else if(m_pOrSite != NULL) { m_pOrSite->m_iCopyNum--; }
	m_ardSpace =  NULL; m_pOrSite = NULL; m_pSpacePointer = NULL; m_pScalePointer = NULL;
}
void CSite::Overwrite(CSite &Site,int OrSite)	{
	m_bReal = false;
	m_pSpacePointer = Site.m_ardSpace; m_pScalePointer = &Site.m_iScale;
	m_iOriSite = OrSite;
	m_iCopyNum=0;
	Site.m_iCopyNum++;
	m_pOrSite = &Site;;
}
void CSite::ResetSite()	{
	if(m_bReal == true) { return; }
	m_pScalePointer = &m_iScale; m_pSpacePointer = m_ardSpace;
	m_pOrSite->m_iCopyNum--; m_pOrSite = NULL;
	m_bReal = true; m_iOriSite = -1;
}
void CSite::ZeroSite(bool ForceAll)	{
	int i;
	if(m_bReal == false && ForceAll == false) { return; }
	FOR(i,*m_Char) { m_ardSpace[i] = 0.0; }
	*m_pScalePointer = 0; m_iScale = 0;
}
void CSite::CopyVals(CSite *Site,bool ForceReal)	{
	int i;
	if(ForceReal) { FOR(i,*m_Char) { m_ardSpace[i] = Site->m_ardSpace[i]; } m_iScale = Site->m_iScale; }
	else { FOR(i,*m_Char) { m_pSpacePointer[i] = Site->m_pSpacePointer[i]; } *m_pScalePointer = *Site->m_pScalePointer; }
}
ostream &operator<<(ostream &os, const CSite &Site)	{
	int i;
	os << "Space[Scale="<<*Site.m_pScalePointer<<"; MEM: "<< Site.m_pSpacePointer << "]:"; FOR(i,*Site.m_Char) { os << "\t" << Site.m_pSpacePointer[i]; }
	return os;
}

///////////////////////////////////////////////////////////////////////////////
// Process definition
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Core constructor
// ----------------
// Data = Sequence data
// Tree = Phylogenetic tree
// Name = Name of process (default = "Unknown process")

CBaseProcess::CBaseProcess(CData *Data, CTree *Tree,string Name, bool DoTreeSearch)	{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseProcess++;
#endif
	// Get process ID
	m_iProcID = GetProcID();
	// Declare space as empty
	m_pSubTree = NULL; m_arModelL = NULL; m_pRate = NULL;
	m_ardQP = NULL; m_ardL = NULL; m_pProcProb = NULL;
	// Get data and tree pointers
	m_pData = Data; m_pTree = Tree; m_sName = Name;
	// Initialise variables
	m_iChar = m_iChar2 = m_iDataChar = m_iSize = m_iSiCh = -1; m_iHiddenChar = 1;
	m_iSpaceSize = m_iSpaceNoSeq = 0;
	m_dBaseVal = 1.0;
	m_piProbFactor = NULL; GET_MEM(m_piProbFactor,int,1); *m_piProbFactor = 1;
	// Set flags
	m_bPseudoProcess = false; m_bIsGamma = false; m_bBraDerReady = false; m_bDetailedOutput = false;
	m_bCheckSpace = false; m_bSubTreeActive = false; m_bCompressedSpace = false; m_bMaxRate = false;
	m_bFailedL = true; m_bDoingPartial = false; m_bModelPerBranch = false; m_bDoStandardDecompose = true;
	m_bIsProcessCopy = false;
	m_bAllowTreeSearch = true;
	// Initialise data
	if(Data != NULL) { m_DataType = Data->m_DataType; } else { m_DataType = NONE; }
	m_iRootQ = -1;
}

CBaseProcess::~CBaseProcess()		{
#if DO_MEMORY_CHECK
	memory_check.CountCBaseProcess--;
#endif
	CleanPar();
	CleanSpace();
	if(m_pRate != NULL) { delete m_pRate; }
	if(!m_bIsProcessCopy) { // If a process copy then probfactor and procprob are defined elsewhere
		if(m_piProbFactor != NULL) { delete m_piProbFactor; }
		if(m_pProcProb != NULL) { delete m_pProcProb; }
	}
}

void CBaseProcess::MakeABET(EDataType Type) {
	m_sABET = DataStates(Type);
	m_iABET_length = NumStates(Type);
}

///////////////////////////////////////////////////////////////////////
// Functions changing memory
// 1. Apply new data
void CBaseProcess::ApplyNewData(CData *D, bool RedoSpace) {
	m_pData = D;
	if(D== NULL) { m_pData = NULL; }
	else if(RedoSpace) {
		MakeCalcSpace(true);
	} else {
		if(D->m_iNoSeq > m_iSpaceNoSeq) { Error("\nTrying to apply new data to a process with not enough space (Seq)\n\n"); }
		if(D->m_iSize > m_iSpaceSize) { Error("\nTrying to apply new data to a process with not enough space (Size)\n\n"); }
		m_iSize = D->m_iSize;
	}
}
// 2. Apply new tree
void CBaseProcess::ApplyNewTree(CTree *T) {
	m_pTree = T;
	if(T->NoSeq() > m_iSpaceNoSeq) { Error("\nTrying to apply new tree with not enough space\n\n"); }
}


//////////////////////// Output functions /////////////////////////
ostream &CBaseProcess::Output(ostream &os)	{
	int i;
	os<< "\n----- Process: " << m_sName << " : ID = " << m_iProcID << "; Rate = " << m_pRate->Val() << " -----";
	if(m_bPseudoProcess == true) { os << "\nPseudoprocess"; }
	else {
		os << "\n" << m_vpPar.size() << " Parameters";
		FOR(i,(int)m_vpPar.size()) { os << "\n\t" << *m_vpPar[i]; }
		os << "\n" << m_vpQMat.size() << " Q matrices";
		FOR(i,(int)m_vpQMat.size()) { os << "\n\t" << *m_vpQMat[i]; os << "\nEqm: " << m_vpQMat[i]->Eqm();  }
	}
	return os;
}

vector <double> CBaseProcess::GetPartL(int site, int Node) {
	int i;
	double *pL;
	vector <double> Return;
	bool NTreal = true;
	if(IsSubTree()) { // Deal with subtrees
		if(Node < Tree()->NoSeq()) { if(m_viLeafMap[Node] == -1) { NTreal = false; } else { Node = m_viLeafMap[Node]; } } else { NTreal = false; }
	} else { if(Node >= Tree()->NoSeq()) { NTreal = false; } }

	if(NTreal) {
		if(m_pData->m_ariSeq[Node][site] == m_iChar) { FOR(i,m_iChar) { Return.push_back(1.0); } }
		else { FOR(i,m_iChar) { if(i % DataChar() == m_pData->m_ariSeq[Node][site]) { Return.push_back(1.0); } else { Return.push_back(0.0); } } }
	} else {
		pL = ForceRealFd(Node,site);
		FOR(i,m_iChar) { Return.push_back(*pL++); }
	}
	if((int)Return.size() != m_iChar) { Error("The output vector length doesn't match the expected in CBaseProcess::GetPartL(...)\n"); }
	return Return;
}

int CBaseProcess::GetPartLScale(int Node,int Site) {
	int i, Return = 0;
	bool NTreal = true;
	if(IsSubTree()) { // Deal with subtrees
		if(Node < Tree()->NoSeq()) { if(m_viLeafMap[Node] == -1) { NTreal = false; } else { Node = m_viLeafMap[Node]; } } else { NTreal = false; }
	} else { if(Node >= Tree()->NoSeq()) { NTreal = false; } }
	if(NTreal) { return 0; } else { return *ForceRealFdSc(Node,Site); }
	return Return;
}


void CBaseProcess::OutPartL(int site, int Node, ostream &out)	{
	int i;
	double *pL;
	bool NTreal = true;
	out << "\n\tScale\t" << GetPartLScale(Node,site) << " : " << GetPartL(site,Node);
}

ostream &CBaseProcess::SiteOut(ostream &os, int Node, int PBeg, int PEnd)	{
	int NodePos;
	if(PBeg == -1) { PBeg = 0; }
	if(PEnd == -1) { PEnd = m_pData->m_iSize - 1; }
	if(Node == -1) { Node = PartLNode(); }
	NodePos = InitNodePos(Node);
	assert(InRange(Node,0,m_pTree->NoNode()+1));
	assert(InRange(PBeg,0,(int)m_pData->m_iSize) && InRange(PEnd,0,(int)m_pData->m_iSize));
	for(PBeg;PBeg < PEnd;PBeg++)	{
		os << "\nSite["<<PBeg<<"] = "; //if(IsProb(Lsum(PBeg).Prob())) { os << Lsum(PBeg).LogP(); }
		os << "\n\tFd" << m_vSpace[NodePos+PBeg] << "\n\tBk" << m_vBackSp[NodePos+PBeg];
	}
	return os;
}

ostream& CBaseProcess::SiteL(int Site, ostream &os)	{
	assert(InRange(Site,-1,m_iSize));
	if(Site == -1) {
		int i; FOR(i,m_iSize) { os << "\nSite["<<i<<"]: " << Lsum(Site).LogP(); }
	} else { os << Lsum(Site).LogP(); }
	return os;
}

void CBaseProcess::OutQ(int QNum, ostream &os)	{
	int i;
	if(QNum >= 0) {
		if(!InRange(QNum,0,(int)m_vpQMat.size())) { Error("\nTrying to output QMat outside range...\n\n"); }
		m_vpQMat[QNum]->OutQ(os);
	} else { FOR(i,(int) m_vpQMat.size()) { m_vpQMat[i]->OutQ(os); } }
}

void CBaseProcess::OutPT(ostream &os, int Branch) {
	int i;
	os << "\n\tBranch " << Branch << " P(" << Tree()->B(Branch) << ")";
	FOR(i,m_iChar2) {
		if(i % m_iChar == 0) { os << "\n\t"; }
		os << "\t" << PT(Branch)[i];
	}
}

void CBaseProcess::OutEqm(int ENum,ostream &os)	{
	int i;
	if(ENum == -1)	{
		FOR(i,(int)m_vpQMat.size()) { os << "\nEqm["<<i<<"]: " << m_vpQMat[i]->Eqm(); }
	} else {
		assert(InRange(ENum,0,(int)m_vpQMat.size()));
		os << "\nEqm["<<ENum<<"]: " << m_vpQMat[ENum]->Eqm();
	}
}

//////////////////////////////////////////////////////////
// Calculate the rate of a full specified Q matrix as sum of pi[i]*Q[i][i]
// QMat specifies which matrix to compute (default = 0);
double CBaseProcess::CalcRate(int QMat)	{
	// Error check
	assert(InRange(QMat,0,(int)m_vpQMat.size()));
	// Get the substitution rate
	return m_vpQMat[QMat]->OverallSubRate();
}

//////////////////////// Parameter functions //////////////////////
// Set rate
double CBaseProcess::Rate(double NewRate, bool MakeRateOpt)	{
	if(NewRate >= 0)	{ m_pRate->SetVal(NewRate,true); }	// Set rate if required
	if(MakeRateOpt == true) { assert(m_pRate != NULL); m_vpPar.push_back(m_pRate); } // Add rate to parameter set if required
	return m_pRate->Val();
}

CQPar *CBaseProcess::AddRatePar2Opt()	{
	int i;
//	cout << "\nProcess: " << m_sName << " = " << NoPar();
	FOR(i,NoPar()) { if(m_vpPar[i] == m_pRate) { return m_pRate; } }
	m_vpPar.push_back(m_pRate);

//	cout << " --> " << NoPar();
	return m_pRate;
}

// Remove parameter function
void CBaseProcess::RemovePar(string Name, bool AllowFail)	{
	int i = 0,count = 0;
	vector <CQPar *>::iterator iPar;
	IFOR(iPar,m_vpPar)	{
		if(m_vpPar[i]->Name().find(Name) != string::npos) {
			m_vpPar.erase(iPar); count ++;
			if(iPar == m_vpPar.end()) { break; }
		}
		else { i++; }
	}
	if(!AllowFail && count != 1) { cout << "\nError in CBaseProcess::RemovePar(" << Name << ") which removed " << i << " parameters...\n"; exit(-1); }
}

// Gets a parameter of a specific name
CQPar *CBaseProcess::GetPar(string Name, bool AllowFail) {
	int i = 0, count = 0;
	vector <CQPar *>::iterator iPar;
	CQPar *RetPar = NULL;
	IFOR(iPar,m_vpPar)	{
		if(m_vpPar[i]->Name().find(Name) != string::npos) {
			RetPar = m_vpPar[i]; count ++;
			if(iPar == m_vpPar.end()) { break; }
		}
		i++; // Only move counter on if not erased
	}
	if((!AllowFail && count != 1) || RetPar == NULL) { cout << "\nError in CBaseProcess::GetPar(" << Name << ") which found " << i << " parameters...\n"; exit(-1); }
	return RetPar;
}


//////////////////////// Process copying functions ////////////////
//
// Returns a pseudocopy of this process for adding to a vector of processes
// in the model
// The only parameter included in the model is its rate
CBaseProcess * CBaseProcess::RateProcessCopy()	{
	CBaseProcess *NewProc;
	static int CopyNum = 0;
	string Name;
	int i;
	// Do naming and initialise
	Name = "Pseudo" + m_sName + "(CopyID=" + int_to_string(CopyNum++) + ")";
	NewProc = new CBaseProcess(m_pData, m_pTree,Name);
	NewProc->MakeBasicSpace(m_iChar);
	NewProc->m_DataType = m_DataType;
	NewProc->m_bAllowTreeSearch = m_bAllowTreeSearch;
	// Do the copy
	NewProc->CleanPar();
	NewProc->CleanQ();
	FOR(i,(int)m_vpQMat.size())	{ NewProc->m_vpQMat.push_back(m_vpQMat[i]); }
	FOR(i,(int)m_viQ2Bra.size())	{ NewProc->m_viQ2Bra.push_back(m_viQ2Bra[i]); }
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

// Returns a pseudocopy of this process for adding to a gamma distribution model
// This is the same as RateProcessCopy, except it also ties the Prob() together
CBaseProcess *CBaseProcess::GammaRateProcessCopy()	{
	CBaseProcess *NewProc;
	NewProc = RateProcessCopy();
	// Replace the probability stuff
	delete NewProc->m_pProcProb; NewProc->m_pProcProb = NULL;
	delete NewProc->m_piProbFactor; NewProc->m_piProbFactor = NULL;
	NewProc->m_pProcProb = m_pProcProb; NewProc->m_piProbFactor = m_piProbFactor;
	return NewProc;
}
//////////////////////// Space functions ///////////////////////////////////////////////////////////
// The structure of this space depends on whether the model is intended for tree search or not
void CBaseProcess::MakeCalcSpace(bool AllowRemake)	{
	int i, NoBra = (2*m_pData->m_iNoSeq)-3, NoNode = (2*m_pData->m_iNoSeq)-2;
	if(!m_bAllowTreeSearch) { NoNode -= m_pData->m_iNoSeq; }
	if(!m_vSpace.empty() || !m_vBackSp.empty()) {
		if(AllowRemake == true) { CleanCalcSpace(); } else { Error("Trying to remake space when not allowed\n\n"); }
	}
	m_iSpaceSize = m_iSize = m_pData->m_iSize;
	m_iSpaceNoSeq = m_pData->m_iNoSeq;
	m_iSiCh = m_iChar * m_iSize;
	// Get space vectors (The last node is the final partial likelihood
	m_vSpace.assign(m_iSize * (NoNode + 1),CSite(&m_iChar));
	m_vBackSp.assign(m_iSize * (NoNode + 1),CSite(&m_iChar));
	GET_MEM(m_ardQP,double,NoBra * m_iChar2);
	GET_MEM(m_ardL,CProb,m_iSize);
	m_bCheckSpace = true; m_bCompressedSpace = false; m_bBraDerReady = false;

}
void CBaseProcess::MakeBasicSpace(int Char)	{
	int NoBra = (2*m_pData->m_iNoSeq)-3;
	string Name;
	// Do some assertions and initialisation
	assert(m_pData != NULL);
	m_iChar = Char;
	m_iDataChar  = m_pData->m_iChar;
	m_iChar2 = m_iChar * m_iChar;
	// Get the rate parameters and the process probability parameters
	Name = "Rate(" + int_to_string(m_iProcID) + ")";
	m_pRate = new CQPar(Name,m_iChar,1.0,false,0.0,MAX_PAR_VALUE,REPLACE);
	Name = "Prob(" + int_to_string(m_iProcID) + ")";
	m_pProcProb = new CQPar(Name,m_iChar,1.0,false,MIN_PROB,1.0,REPLACE);
	GET_MEM(m_ardPT,double,NoBra * m_iChar2);
	m_bCheckSpace = false; m_bCompressedSpace = false; m_bBraDerReady = false;
}
void CBaseProcess::MakeBasicSpace(EDataType Type)	{
	m_DataType = Type;
	MakeABET(Type);
	MakeBasicSpace(NumStates(Type));
}
//////////////////////// Function that resets the calculation space (used for parametric bootstrapping) //////
void CBaseProcess::ResetCalcSpace()	{
	if(m_pData->m_iNoSeq != Tree()->NoSeq()) { Error("\nData(" + int_to_string(m_pData->m_iNoSeq) + " and tree (" + int_to_string(Tree()->NoSeq()) + ") need to match in CBaseProcess::ResetCalcSpace(...)\n\n"); }
	MakeCalcSpace(true);
}

// Special memory check to see if tree search is allowed: -1 = Allow; -2 = NoAllow;
bool CBaseProcess::CheckAllowTreeSearch() { return m_bAllowTreeSearch; }

// Sets all values in space to zero. Useful for debugging
void CBaseProcess::ZeroSpace()	{
	vector <CSite>::iterator i_Sp;
	IFOR(i_Sp,m_vSpace)		{ i_Sp->ZeroSite(); }
	IFOR(i_Sp,m_vBackSp)	{ i_Sp->ZeroSite(); }
}

// Cleaning routines
void CBaseProcess::CleanSpace()	{
	DEL_MEM(m_ardPT);
	CleanCalcSpace();

	// Clear important pointers
	m_pTree = NULL; m_pData = NULL;
	// Reset flags
	m_bCheckSpace = false; m_bCompressedSpace = false; m_bBraDerReady = false;
	CleanCPMapping();
}
void CBaseProcess::CleanCalcSpace()	{
	int i;
	DEL_MEM(m_ardQP);
	DEL_MEM(m_ardL);
	if(!m_bIsProcessCopy) {
		FOR(i,m_vpQMat.size()) { if(m_vpQMat[i] != NULL) { delete m_vpQMat[i]; m_vpQMat[i] = NULL; } }
	} else { FOR(i,m_vpQMat.size()) { m_vpQMat[i] = NULL; } }
	m_vpQMat.clear();
	m_vSpace.clear(); m_vBackSp.clear();
	m_bCompressedSpace = false;
	if(MainTree() != NULL) { MainTree()->SetFastCalcOK(false); }
}
void CBaseProcess::CleanPar()	{
	int i;
	FOR(i,(int)m_vpPar.size()) { delete m_vpPar[i]; m_vpPar[i] = NULL; }
	if(!m_bIsProcessCopy) { FOR(i,(int)m_vpEqm.size()) { if(m_vpEqm[i] != NULL) { delete m_vpEqm[i]; m_vpEqm[i] = NULL; } } }

	m_vpPar.clear(); m_vpEqm.clear();
	m_bBraDerReady = false;
}
void CBaseProcess::CleanQ()	{
	int i;
	FOR(i,(int)m_vpQMat.size()) { delete m_vpQMat[i]; m_vpQMat[i] = NULL; }
	m_vpPar.clear();
	m_bBraDerReady = false;
}
/////////////////////////////////////////////////////////////////////////////////////
// Space access functions
void CBaseProcess::CleanScale(int NodeNum, bool ForceAll)	{
	int i, NodePos;
	assert(m_bCheckSpace == true);
	m_bBraDerReady = false;
	assert(InRange(NodeNum,-1,m_pTree->NoNode()) || (NodeNum == 2 && m_pData->m_iNoSeq == 2));
	if(NodeNum == -1) { // Clean the base scale
		FOR(i,m_iSize){ *LScale(i) = 0; }
	} else {			// Clean the node scale
		NodePos = InitNodePos(NodeNum);
		if(ForceAll == false)	{ FOR(i,m_iSize)	{ *QkFdSc(NodePos) = 0; *QkBkSc(NodePos++) = 0; } }
		else					{ FOR(i,m_iSize)	{ *QkForceRealFdSc(NodePos) = 0; *QkForceRealBkSc(NodePos++) = 0; } }
}	}

void CBaseProcess::TransScale(int NT, int NF,bool First, bool Partial,bool ForceReal)	{
#if ALLOW_SCALE == 1
	int i, NPosT, NPosF;
	assert(InRange(NT,-1,m_pTree->NoNode()+1) && InRange(NF,0,m_pTree->NoNode()));
	// if NT >= 0 then its a normal tranfer; otherwise transfer to PartLNode
	if(NT >= 0)		{ NPosT = InitNodePos(NT); } else { NPosT = InitNodePos(PartLNode()); assert(First == true); }
	// Set the other pointer
	NPosF = InitNodePos(NF);
	// Do the normal scale
	if(ForceReal == true || m_bCompressedSpace == false)	{
		if(Partial == true)	{
			if(First == true)	{
				FOR(i,m_iSize)	{
					*QkForceRealFdSc(NPosT++) = *QkForceRealFdSc(NPosF++);
			}	} else {
				FOR(i,m_iSize)	{
					*QkForceRealBkSc(NPosT) = *QkForceRealFdSc(NPosF);
					*QkForceRealFdSc(NPosT) += *QkForceRealBkSc(NPosT);
					NPosF++; NPosT++;
		}	}	} else {
			if(First == true)	{
				FOR(i,m_iSize)	{
					*QkForceRealFdSc(NPosT++) = *QkForceRealFdSc(NPosF++);
			}	} else {
				FOR(i,m_iSize)	{
					*QkForceRealFdSc(NPosT++) += *QkForceRealFdSc(NPosF++);
	}	}	}	} else {
		// Do the column sorted scale
		if(Partial == true)	{
			if(First == true)	{
				FOR(i,m_iSize)	{
					if(!QkFdReal(NPosT)) { NPosT++; NPosF++; continue; }
					*QkFdSc(NPosT++) = *QkFdSc(NPosF++);
			}	} else {
				FOR(i,m_iSize)	{
					if(!QkFdReal(NPosT)) { NPosT++; NPosF++; continue; }
					*QkBkSc(NPosT) = *QkFdSc(NPosF);
					*QkFdSc(NPosT) += *QkBkSc(NPosT);
					NPosF++; NPosT++;
		}	}	} else {
			if(First == true)	{
				FOR(i,m_iSize)	{
					if(!QkFdReal(NPosT)) { NPosT++; NPosF++; continue; }
					*QkFdSc(NPosT++) = *QkFdSc(NPosF++);
			}	} else {
				FOR(i,m_iSize)	{
					if(!QkFdReal(NPosT)) { NPosT++; NPosF++; continue; }
					*QkFdSc(NPosT++) += *QkFdSc(NPosF++);
	}	}	}	}
#endif
}
void CBaseProcess::DoScale(int Node, bool ForceRealScale)	{
#if ALLOW_SCALE == 1
	int i,j, NodePos = InitNodePos(Node);
	double BkMax,FdMax, *p_f = NULL, *p_b = NULL;
	if(Node == m_pTree->NoNode()) { return; }	// Don't scale if its from the last node
	assert(Node >= 0 && Node < m_pTree->NoNode());
	FOR(i,m_iSize)	{
		if(!QkFdReal(NodePos) && ForceRealScale == false && m_bCompressedSpace == true) { NodePos++; continue; }
		// Find Max at site
		p_f = QkForceRealFd(NodePos); p_b = QkForceRealBk(NodePos);
		BkMax = FdMax = 0.0; FOR(j,m_iChar)	{
			if(*(p_f) > FdMax) { FdMax = *(p_f); }
			if(*(p_b) > BkMax) { BkMax = *(p_b); }
			p_f++; p_b++;
		}
		// If max are small enough then scale
		if(BkMax < P_SCALE_VAL)	{
			// Allow true zero
			if(Double_Zero(BkMax) || *QkForceRealBkSc(NodePos) > 1000000000) { p_b = QkForceRealBk(NodePos); FOR(j,m_iChar) { *(p_b++) = 0.0; *QkForceRealBkSc(NodePos) = 0; } }
			else {
				// Do normal scale
				while(BkMax < 1.0)	{ BkMax *= 10; p_b = QkForceRealBk(NodePos); FOR(j,m_iChar)	{ *(p_b++) *= 10; } QkForceRealBkSc(NodePos)[0]++; }
				while(BkMax > 10.0)	{ BkMax /= 10; p_b = QkForceRealBk(NodePos); FOR(j,m_iChar)	{ *(p_b++) /= 10; } QkForceRealBkSc(NodePos)[0]--; }
		}	}
		if(FdMax < P_SCALE_VAL)	{
			// Allow true zero
			if(Double_Zero(FdMax)  || *QkForceRealFdSc(NodePos) > 1000000000) { FOR(j,m_iChar) { QkForceRealFd(NodePos)[j] = 0.0; } *QkForceRealFdSc(NodePos) = 0;  }
			else {
				// Do normal scale
				while(FdMax < 1.0)	{ FdMax *= 10; p_f = QkForceRealFd(NodePos); FOR(j,m_iChar)	{ *(p_f++) *= 10; } QkForceRealFdSc(NodePos)[0]++; }
				while(FdMax > 10.0)	{ FdMax /= 10; p_f = QkForceRealFd(NodePos); FOR(j,m_iChar)	{ *(p_f++) /= 10; } QkForceRealFdSc(NodePos)[0]--; }
		}	}
		NodePos++;
	}
	p_f = NULL; p_b = NULL;
#endif
}

void CBaseProcess::CopyNode(int NodeFr, int NTo)	{
	int i,j,NodePos1,NodePos2;
	double *p_SpTo, *p_SpFr, *p_BkTo, *p_BkFr;
	// Check entry conditions
	assert(InRange(NodeFr,0,m_pTree->NoNode()));
	assert(InRange(NTo,0,m_pTree->NoNode()));
	// Set Quick node access up
	NodePos1 = InitNodePos(NodeFr);
	NodePos2 = InitNodePos(NTo);
	FOR(i,m_iSize) {
		// Get pointers
		p_SpFr = QkFd(NodePos1);			p_SpTo = QkFd(NodePos2);
		p_BkFr = QkBk(NodePos1);			p_BkTo = QkBk(NodePos2);
		// Do scales
		*QkFdSc(NodePos2) = *QkFdSc(NodePos1);
		*QkBkSc(NodePos2) = *QkBkSc(NodePos1);
		FOR(j,m_iChar)	{
			// Do the partial likelihoods
			*(p_SpTo++) = *(p_SpFr++);
			*(p_BkTo++) = *(p_BkFr++);
		}
		// Move pointers along
		NodePos1++; NodePos2++;
	}
	// Clear pointers
	p_SpTo = NULL; p_SpFr = NULL; p_BkTo = NULL; p_BkFr = NULL;

}
/////////////////////////////////////////////////////////////////////////////////////
// Column sorting routines for fast likelihood computation
void CBaseProcess::PrepareFastCalc(vector <int> *C)	{
	int NodeBase,i;
	assert(MainTree()->FastCalcOK() == false && m_bCompressedSpace == false);
	assert(m_vSpace.size() == m_vBackSp.size());
	// Some entry conditions
	if(!ALLOW_FAST_CALC) { return; }
	if(m_vSpace.empty() || m_vBackSp.empty()) { return; }	// Can't prepare calculation when the space isn't ready...
	if(IsSubTree() || m_pTree->IsCutTree()) { return; }
//	cout << "\n <--------------------- Fast calc ------------------->";
	// If required get the compressed data (shouldn't be required because its passed by model)
	if(C == NULL) { *C = GetCompressedData(m_pTree,m_pData); }
	assert(C->size() == (m_pData->m_iSize * (m_pTree->NoNode() - m_pTree->NoSeq())));
	if(!m_bAllowTreeSearch) { NodeBase = 0; }
	else { NodeBase = m_pTree->NoSeq() * m_pData->m_iSize; }
	FOR(i,(int)C->size())	{
		if(C->at(i) == -1) { continue; }
		m_vSpace[NodeBase + i].Overwrite(m_vSpace[NodeBase + C->at(i)],NodeBase + C->at(i));
		m_vBackSp[NodeBase + i].Overwrite(m_vBackSp[NodeBase + C->at(i)],NodeBase + C->at(i));
	}
	m_bCompressedSpace = true;
}

// Revert to simple likelihood computations
void CBaseProcess::DecompressSpace()	{
	int i;
	if(m_bCompressedSpace == false) { return; }
//	cout << "\n <--------------------- Decompress ------------------->";
	FOR(i,(int)m_vSpace.size())	{ m_vSpace[i].ResetSite(); m_vBackSp[i].ResetSite(); }
	m_bCompressedSpace = false;
}

///////////////////////////////////////////////////////////////////////
// Eqm distribution functions
vector <double> CBaseProcess::SimpleEqm(int QMatID, bool AllowSumNotOne)	{
	int i,Mat = -1;
	vector <double> eqm;
	FOR(i,(int)m_vpEqm.size()) {
		if(m_vpEqm[i]->IsID(QMatID)) { assert(Mat == -1); Mat = i; }
	}
	// Get the equilibrium distribution
	if(Mat == -1) { eqm.assign(m_iChar,1.0/(double) m_iChar); } else { assert(m_iChar == m_vpEqm[Mat]->Char()); eqm = m_vpEqm[Mat]->Eqm();  }
	return eqm;
}


// Functions for adding Q matrices to the models
CQMat *CBaseProcess::Add_QMat(string Name,EDataType Type)	{ CQMat *Mat; Mat = new CQMat(Type,Name); Mat->InitQ(m_dBaseVal); m_vpQMat.push_back(Mat); return Mat; }
CQMat *CBaseProcess::Add_QMat(string Name, int Char)		{ CQMat *Mat; Mat = new CQMat(Char, OTHER, Name); Mat->InitQ(m_dBaseVal); m_vpQMat.push_back(Mat); return Mat; }
CQMat *CBaseProcess::Add_CodRedQMat(string Name, int Char) 	{ CQMat *Mat; Mat = new CQMat(Char, COD_RED, Name); Mat->InitQ(m_dBaseVal); m_vpQMat.push_back(Mat); return Mat; }
// Function to prepare the Q matrices for likelihood computations
bool CBaseProcess::PrepareQMats(vector <int> Qs2do, bool DoScale)	{
	bool RetVal = true;
	int i,j;
	if(m_pSubTree != NULL && m_vpQMat.size() > 1) { cout << "\nWarning: Haven't considered shorter PrepareQMats functions for subtrees"; Error(" "); }
	// If required do all matrices
	if(Qs2do.empty() || Qs2do.size() == m_vpQMat.size()) { Qs2do.clear(); FOR(i,(int)m_vpQMat.size()) { Qs2do.push_back(i); } }
	// Apply the matrices
//	cout << "\nAbout to do Qs" << flush;
	FOR(i,(int)Qs2do.size())	{
//		cout << "\nDoing Q["<<i<<"]: " << flush;
		// If a normal process then do all the Q stuff
		if(!m_bPseudoProcess)	{
//			cout << "\nScaleQ" << flush;
			if(m_vpQMat[i]->IsLocked()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); continue; }
			// Initialise the Q matrix
			m_vpQMat[i]->InitQ(m_dBaseVal);
			// Apply the parameters
			FOR(j,(int)m_vpPar.size()) {
//				cout << "\nParameter["<<j<<"] " << m_vpPar[j]->Name() << ": " << m_vpPar[j]->Val() << flush;
				m_vpPar[j]->UpdatePar(); m_vpQMat[i]->ApplyPar2Q(m_vpPar[j]); }
			// Apply the eqm distributions
			FOR(j,(int)m_vpEqm.size()) { m_vpEqm[j]->ApplyEqm2QMat(m_vpQMat[i]->Q(),m_vpQMat[i]->ID()); }
			// Do anything else that needs doing to the rate matrix
			DoPostHocQMatUpdate();
			// Do the diagonal
			m_vpQMat[i]->DoQDiag();
//			cout << "\nThe QMat is:\n"; m_vpQMat[i]->OutQ();
			// Decompose the matrix
			m_vpQMat[i]->Decompose(SimpleEqm(m_vpQMat[i]->ID()),DoScale,true,m_pRate->Val());

#if FUNC_COUNTERS
			MakeQ_Log_Counter++;
#endif
	}	}
//	cout << "\nCreated QMat: " << flush;
//	m_vpQMat[0]->OutQ();
//	cout << "\n\nEqm: " << m_vpEqm[0]->Eqm();
	return RetVal;
}
// Function to find the Q matrix for a specific branch
int CBaseProcess::QMat4Bra(int Branch)	{
	int Q2return = 0;
	// If require get the
	if(!m_viQ2Bra.empty())	{ FOR(Q2return,(int)m_viQ2Bra.size()) { if(m_viQ2Bra[Q2return] == Branch) { break; } } }
	assert(Q2return < (int)m_vpQMat.size());
	return Q2return;
}

// Function returning PT matrix from Q matrices
bool CBaseProcess::Make_PT(int Branch, bool RedoRate)	{
	Tree()->UpdateB(Branch);
	if(RedoRate) { m_vpQMat[QMat4Bra(Branch)]->ScaleQ(m_pRate->Val()); }
	if(m_bModelPerBranch) {
#if DEVELOPER_VERSION == 0
		cout << "\nPer branch models are for developers only";
#endif
		cout << "\nI need to do calculations per branch...";
	}
	return m_vpQMat[QMat4Bra(Branch)]->MakePT(Tree()->B(Branch),PT(Branch));
}
// Function returning the PT matrix
vector <double> CBaseProcess::GetPT(int Br)	{
	int i;
	vector <double> v_mat;
	double * p_mat;
	assert(InRange(Br,0,Tree()->NoBra()));
	p_mat = PT(Br);
	FOR(i,m_iChar2) { v_mat.push_back(p_mat[i]); }
	return v_mat;
}

// Function to get the equilibrium at the root
vector <double> CBaseProcess::RootEqm()	{
	vector <double> E;
	// Get the equilibrium of the root
	if(m_iRootQ == -1) { E = m_vpQMat[0]->Eqm(); }
	else { assert(m_iRootQ >= 0 && m_iRootQ < (int)m_vpQMat.size()); E = m_vpQMat[m_iRootQ]->Eqm(); }
	return E;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions that will be called externally
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Prepare the Q matrices for the calculations
bool CBaseProcess::PrepareLikelihood(bool DoQ, bool ForceRemake, bool DoScale)	{
	int i;
	vector <bool>  LockedProcs;
//	cout << "\nPrepareLikelihood for process " << flush;
//	cout << m_sName << flush;
	// Allow remakes if required
	if(ForceRemake) { FOR(i,(int)m_vpQMat.size()) { LockedProcs.push_back(m_vpQMat[i]->IsLocked()); m_vpQMat[i]->Unlock(); } }
	// Prepare Q matrices
	if(DoQ&& m_bDoStandardDecompose) { if(PrepareQMats(DoScale) == false) { return false; } }
	// Get the equilibrium of the root

	m_vdEqm = RootEqm();
	// Restore lock conditions
	if(ForceRemake) { FOR(i,(int)m_vpQMat.size()) { if(LockedProcs[i]) { m_vpQMat[i]->Lock(); } } }
	return true;
}

bool CBaseProcess::CreatePTMats(int Bra)	{
	int i;
	// Check entry conditions
	assert(InRange(Bra,-1,Tree()->NoBra()));
	// Do zero rates
	if(Rate() < DX) { if(Bra < 0) { FOR(i,Tree()->NoBra()) { IMat(PT(i),m_iChar); } } else { IMat(PT(Bra),m_iChar); } }
	// Do the rest of the PT
	FOR(i,(int)m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); }
	if(Bra < 0) { FOR(i,Tree()->NoBra()) { Tree()->UpdateB(i); Make_PT(i); } } else { Tree()->UpdateB(Bra); Make_PT(Bra); }
	return true;
}

// Does the likelihood computation
bool CBaseProcess::Likelihood(bool ForceReal)	{
	int i;
	static int Counter =0;
//	cout << "\n----------------------------------------------------------------------------------\nLikelihood comp: " << Counter++;
//	cout << "\n--------------------------------------------------\nLikelihood\nTree: " << *m_pTree;

//	cout << "\nDoing likelihood\nData: " << *m_pData << "\nTree: " << *m_pTree;

//	cout << "\nLikelihood starting at node: "<< m_pTree->StartCalc();

//	cout << "\nDoing likelihood: "; FOR(i,m_vpPar.size()) { cout << m_vpPar[i]->Val() << " "; }
	// Create space if not already done so
	if(m_vSpace.empty()) { MakeCalcSpace(false); }
	// For zero rate models this is really easy
	if(Rate() < DX) { MakeZeroRateLikelihood(); return true; }
	// For maximum rate model this is also really easy (Garbage Collector)
	if(m_bMaxRate) { MakeMaxRateLikelihood(); return true; }
	// Otherwise do the usual computations
	bool OldComp = m_bCompressedSpace;
	m_bBraDerReady = false;
	m_bFailedL = false;
	assert(m_bDoingPartial == false);
	// If the probability of the process is small then stop here
	if(Prob() < SMALL_PROB) { return true; }
	// Adjust the rate for the process (if required)
	FOR(i,(int)m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); }
//	cout << "\nScaled Q rate" << flush;
//	cout << "\nOverall eqm: " << m_vdEqm;
//	cout << "\nHidden eqm:  " <<  m_vpEqm[0]->TransEqm();
//	cout << "\nParameters: "; FOR(i,(int)m_vpPar.size()) { cout << "\n\tPar["<<i<<"]: " << *m_vpPar[i]; }
//	cout << "\nPrepared matrices. QMat: "; OutQ();
//	cout << "\nMaking PT matrices";
	// Make the PT matrices
	FOR(i,Tree()->NoBra()) { Tree()->UpdateB(i); Make_PT(i); }
	// Adjust m_bCompressedSpace if required
	if(ForceReal == true) { m_bCompressedSpace = false; }
	// Do the calculations starting at a point specified in the tree object
	CleanScale(-1,ForceReal);	// Clean the scaling array
	PartialL(Tree(),Tree()->StartCalc(),-1,-1,true);
	// Get the final likelihoods for the process
	FOR(i,m_iSize)	{ m_ardL[i].Assign(Lsum(i)); }
	// Return if okay
	m_bCompressedSpace = OldComp;
//	cout << "\nDone likelihood...\n\n";
	return true;
}

// calculate the log likelihood of the process
double CBaseProcess::LogL()	{
	int i; double lnL = 0.0;
	FOR(i,m_iSize) { lnL += L(i).LogP() * m_pData->m_ariPatOcc[i]; }
	return lnL;
}
// Lsum: returns the sitewise likelihood used for calculations in model
CProb &CBaseProcess::Lsum(int site)	{
	int i;
	double *p_a = PartL(site);
	static CProb dVal;
//	cout << "\nCBaseProcess::m_pData["<<site<<"/"<<m_pData->m_iSize<<"]: " << m_pData;
	assert(site < m_pData->m_iSize);
	CProb temp;
	dVal.Assign(0.0);
	m_vdEqm = RootEqm();
#if DEVELOPER_BUILD == 1
	if(site == 0) {
		cout << "\n---------------- Site["<<site<<"] --------------------\nVec: ";
		FOR(i,m_iChar) { cout << PartL(site)[i] << " "; }
		cout << "\nEqm: " << m_vdEqm << endl;
	}
#endif
	FOR(i,m_iChar)	{
		if(my_isnan(*(p_a))) {
			int j;
			cout << "\nIn LSum(site=" << site<< "): ";
			cout << "\nBroken\nThe process looks like this: " << *this;
			cout << "\nThe P(t) matrices are: ";
			FOR(j,Tree()->NoBra()) { cout << "\nBranch["<<j<<"]:"; OutPT(cout,j); }
			cout << "\n<PartL(site=" << site << "): "; FOR(j,m_iChar) { cout << PartL(site)[j] << ":"; }   cout << ">"; }
		if(my_isnan(m_vdEqm[i])) { cout << " <eqm>"; }

		temp.Assign(*(p_a++) * m_vdEqm[i],*LScale(site));
		dVal.Add(temp,true);
	}
#if DEVELOPER_BUILD == 1
	if(site == 0) {
		cout << "\n" << dVal.m_dValue << " x10^-"  << dVal.m_iScale << " == " << dVal.LogP();
	}
#endif
	p_a = NULL;

	return dVal;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Likelihood functions
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Partial likelihood functions
///////////////////////////////////////////////////////
// Function
// ===========
// Calculates the partial likelihoods at the the internal node (iNoFr)
// Function expects all PT matrices to have been already computed
//
// Args:
// ===========
// 1.) iNoTo:	(Node_to) The node number that the branch is going to
// 2.) iNoFr:	(Node_from) The node number that the branch is coming from
///////////////////////////////////////////////////////

///////////////////////////////// Full likelihood functions /////////////////////////////////

void CBaseProcess::MakeZeroRateLikelihood()	{
	int i;
	FOR(i,m_pData->m_iSize)	{
		if(m_pData->m_viNoChange[i] == -1) { m_ardL[i].Zero(); continue; }
		if(!InRange(m_pData->m_viNoChange[i],0,(int)m_vdEqm.size())) {
			cout << "\nGoing to throw assert error... Site["<<i<<"] which is: ";
			int j; FOR(j,m_pData->m_iNoSeq) { cout << m_pData->m_ariSeq[j][i] << " "; }
			cout << "\nAnd m_pData->m_viNoChange[" << i<<"]: " << m_pData->m_viNoChange[i];
			cout << "\nEqm[size=" << m_vdEqm.size() << "]: " << m_vdEqm;
			cout << "\nAnd the vector of m_viNoChange:\n" << m_pData->m_viNoChange;
		}
		assert( InRange(m_pData->m_viNoChange[i],0,(int)m_vdEqm.size()) );
		m_ardL[i].Assign(m_vdEqm[m_pData->m_viNoChange[i]]);
	}
}

void CBaseProcess:: MakeMaxRateLikelihood() {
	int i,j;
	FOR(i,m_pData->m_iSize) {
		m_ardL[i].Assign(1.0);
		FOR(j,NoSeq()) {
			if(m_pData->m_ariSeq[j][i] == m_iChar) { continue; }
			m_ardL[i].Multiply(m_vdEqm[m_pData->m_ariSeq[j][i]],true);
}	}	}

// Calculates partial likelihoods recursively
// --
// Note BlockWriteback is used for PreparePartialL where you don't what to write to iNoFr for the starting node only
void CBaseProcess::PartialL(CTree *pTree, int iNoTo, int iNoFr, int Branch, bool NodeFirst, bool BlockWriteback)	{
    int i,j,NodePos1, NodePos2;
	bool First = true;
	double *I = NULL, *p_a = NULL, *p_b = NULL;
//	cout << "\nPartialL - NT: " << iNoTo << ", iNoFr" << iNoFr << ", Branch: " << Branch << ", NodeFirst: " << NodeFirst;

	// Prepare node to have zero scaling factors
	if(iNoTo >= pTree->NoSeq()) { CleanScale(iNoTo,FlipBool(m_bCompressedSpace)); }
	// Do post order tree traversal
	FOR(i,pTree->NoLinks(iNoTo)) {
		if(pTree->NodeLink(iNoTo,i) == iNoFr || pTree->NodeLink(iNoTo,i) == -1) { continue; }		// Skip the node it came from
		PartialL(pTree,pTree->NodeLink(iNoTo,i),iNoTo,pTree->NodeBra(iNoTo,i),First);	// Traverse the tree
		First = false;
	}
	// If not writeback then return
	if(BlockWriteback) { return; }
	// Transfer information back to iNoFr
	if(iNoFr >= 0)	{ // Do normal nodes
		if(pTree->NodeType(iNoTo) == branch && iNoFr >= pTree->NoSeq())	{			// Do normal branch nodes
			BranchNodePartialL(iNoTo,iNoFr,PT(Branch),NodeFirst);
		} else if(pTree->NodeType(iNoTo) == leaf)	{								// Do normal leaf nodes
			LeafNodePartialL(pTree,iNoTo,pTree->NodeBra(iNoTo,0),iNoFr,PT(pTree->NodeBra(iNoTo,0)),NodeFirst);
		} else if(pTree->NodeType(iNoFr) == leaf && iNoTo >= pTree->NoSeq()) {	// Final calculations

			// Copy the first bit of likelihood to final node
			LeafNodePartialL(pTree,iNoFr,Branch,-1,PT(Branch),true);

			// Now do the second bit;
			NodePos1 = InitNodePos(iNoTo);	// Space from
			NodePos2 = InitNodePos(PartLNode());						// Space to (the storage node)
			if(m_bCompressedSpace) {
				FOR(i,m_iSize) {
					if(!QkFdReal(NodePos2)) { NodePos2++; continue; }
					p_a = QkFd(NodePos2); p_b = QkFd(NodePos1);
					FOR(j,m_iChar) { *(p_a++) *= *(p_b++); }
					// do the scale
					*LScale(i) += *QkFdSc(NodePos1);
					NodePos1++; NodePos2++;
			}	} else { // Uncompressed space
				FOR(i,m_iSize) {
					p_a = QkForceRealFd(NodePos2); p_b = QkForceRealFd(NodePos1);
					FOR(j,m_iChar) { *(p_a++) *= *(p_b++); }
					// do the scale
					*LScale(i) += *QkForceRealFdSc(NodePos1);
					NodePos1++; NodePos2++;
			}	}

//			cout << "\nAnswer:   "; FOR(i,20) { cout << " " << ForceRealFd(PartLNode(),0)[i]; }
			p_a = NULL; p_b = NULL;
	}	} else if(Tree()->NoSeq() == 2) {	// Do the final calculation for 2 species trees
		if(!m_viLeafMap.empty())	{ assert((int)m_viLeafMap.size() >= iNoTo); iNoTo = m_viLeafMap[iNoTo]; }
		NodePos1 = InitNodePos(PartLNode());
		FOR(i,m_iSize) { Vec_by_Data(m_pData->m_ariSeq[iNoTo][i],QkFd(NodePos1++)); }
	}
}
// Simple calculations for BranchNodes
void CBaseProcess::BranchNodePartialL(int SpFrom, int SpTo, double *PT, bool First)	{
	int i,k,NodePos1,NodePos2;
    double *p_a,*p_b, *p_c;
	// If the node comes from last link in StartCalc() then copy information to PartLNode()
	if(SpTo == Tree()->StartCalc() && !m_bDoingPartial) {
		FOR(k,Tree()->NoLinks(SpTo)) { if(Tree()->NodeLink(SpTo,k) == SpFrom) { break; } }
		assert(k != Tree()->NoLinks(SpTo));
		switch(k)	{	// Copy stuff to the internal node if required
		case 1:
			NodePos1 = InitNodePos(SpTo); NodePos2 = InitNodePos(PartLNode());
			FOR(i,m_iSize)	{ m_vBackSp[NodePos1++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));  }
			break;
		case 2:
			NodePos1 = InitNodePos(SpTo); NodePos2 = InitNodePos(PartLNode());
			FOR(i,m_iSize) { m_vSpace[NodePos1++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace)); }
			break;
		};
		SpTo = PartLNode();
	}
	// Initialise quick space access
	NodePos1 = InitNodePos(SpFrom);
	NodePos2 = InitNodePos(SpTo);
	if(First == true)	{ // Copies the probabilities
		// Do first round of calculations for the node
		if(m_bCompressedSpace)	{
			FOR(k,m_pData->m_iSize)	{
				if(!QkFdReal(NodePos2)) { NodePos1++; NodePos2++; continue; }
				VMat(QkFd(NodePos1++),PT,QkFd(NodePos2++),m_iChar);
		}	} else {	// Incompressed space
			FOR(k,m_pData->m_iSize)	{
				VMat(QkForceRealFd(NodePos1++),PT,QkForceRealFd(NodePos2++),m_iChar);
	}	}	} else { // Replaces the probabilities (This has been hard checked)
		if(m_bCompressedSpace)	{
			// Do the final round of calculations, including the backspace storage
			FOR(k,m_pData->m_iSize)	{
				// Deal with when the node isn't real
				if(!QkFdReal(NodePos2)) { NodePos1++; NodePos2++; continue; }
				p_a = QkFd(NodePos1++); p_b = QkBk(NodePos2);
				VMat(p_a,PT,p_b,m_iChar);
				p_c = QkFd(NodePos2++);
				FOR(i,m_iChar) { *(p_c++) *= *(p_b++); }
		}	} else {	// Uncompressed space
			FOR(k,m_pData->m_iSize)	{
				p_a = QkForceRealFd(NodePos1++); p_b = QkForceRealBk(NodePos2);
				VMat(p_a,PT,p_b,m_iChar);
				p_c = QkForceRealFd(NodePos2++);
				FOR(i,m_iChar) { *(p_c++) *= *(p_b++); }
	}	}	}
	// Capture the nodes scaling factors
	TransScale(SpTo,SpFrom,First);
	// Rescale if required
	DoScale(SpTo);
	p_a = NULL; p_b = NULL; p_c = NULL;
}

// LeafNode likelihood computation
// -------------------------------
// Uses adapter function to allow complex models with more states than data characters (e.g. covarion models)
// to transform into useful values

void CBaseProcess::LeafNodePartialL(CTree *pTree, int LeafNode, int Branch, int Sp, double *PT, bool First)	{
	int site,i,k,NodePos, NodePos2;
	static double TempSp[MAX_SPACE],*p_a, *p_b, *p_c;
	vector <double> teqm;
	// If comes from StartCalc(), then transfer info to final node
	// If the node comes from last link in StartCalc() then copy information to PartLNode()
	if(Sp == Tree()->StartCalc() && !m_bDoingPartial) {
		if(Tree()->NoSeq() > 2) {
			FOR(k,Tree()->NoLinks(Sp)) { if(Tree()->NodeLink(Sp,k) == LeafNode) { break; } }
			assert(k != Tree()->NoLinks(Sp));
			switch(k)	{	// Copy stuff to the internal node if required
			case 1:
				NodePos = InitNodePos(Sp); NodePos2 = InitNodePos(PartLNode());
				FOR(i,m_iSize) { m_vBackSp[NodePos++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));  }
				break;
			case 2:
				NodePos = InitNodePos(Sp); NodePos2 = InitNodePos(PartLNode());
				FOR(i,m_iSize) { m_vSpace[NodePos++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));  }
				break;
			};
		}
		Sp = -1;
	}
	// Check entry conditions
	assert(MAX_SPACE > m_iChar); assert(InRange(Sp,-1,m_pTree->NoNode()));
	// If required do the SubLeafNodePartialL
	if(IsSubTree())	{
		assert(!m_viLeafMap.empty() && LeafNode < (int)m_viLeafMap.size());
		// If rqd do the partial likelihood as a branch likelihood
		if(m_viLeafMap[LeafNode] == -1) {
			if(Sp == -1)	{ assert( First == true); BranchNodePartialL(LeafNode,PartLNode(),PT,First); }
			else			{ BranchNodePartialL(LeafNode,Sp,PT,First); }
			return;
		} else {
			LeafNode = m_viLeafMap[LeafNode]; // Rename the node and proceed with calcs
	}	}
	// Get the equilibrium distribution
	teqm = RootEqm(); // ZZXX: m_vpQMat[QMat4Bra(Branch)]->Eqm();
	// Organise space
	if(Sp >= 0)	{		// For normal space
		NodePos = InitNodePos(Sp);
	} else		{		// For final likelihood space
		NodePos = InitNodePos(PartLNode());
	}
	// Loop through sites
	// If TempSp copies to the space
	if(First == true)	{
		if(m_bCompressedSpace)	{
			FOR(site,m_pData->m_iSize)	{
				if(!QkFdReal(NodePos)) { NodePos++; continue; }							// Skip if required
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);			// Get data vector
				p_a = QkFd(NodePos++); p_b = TempSp;									// Set Pointers
				// *COPY* TempSp to the calculation space; BkSp not needed
				FOR(i,m_iChar) { *(p_a++) = *(p_b++); }
		}	} else {
			FOR(site,m_pData->m_iSize)	{
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);			// Get data vector
				p_a = QkForceRealFd(NodePos++); p_b = TempSp;									// Set Pointers
				// *COPY* TempSp to the calculation space; BkSp not needed
				FOR(i,m_iChar) { *(p_a++) = *(p_b++); }
	}	}	} else {
	// Check entry conditions
	// If TempSp multiplies with the space
		if(m_bCompressedSpace)	{
			FOR(site,m_pData->m_iSize)	{
				if(!QkFdReal(NodePos))	{ NodePos++; continue; }					// Skip if required
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);		// Get data vector
				p_a = QkFd(NodePos); p_b = TempSp; p_c = QkBk(NodePos++);			// Set Pointers
				// *MULTIPLY* TempSp and calculation space; BkSp has a copy
				FOR(i,m_iChar) {
					*(p_c++) = *(p_b);			// Copy to BkSp
					*(p_a++) *= *(p_b++);		// Multiply by SpaceTo
		}	}	} else {
			FOR(site,m_pData->m_iSize)	{
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);		// Get data vector
				p_a = QkForceRealFd(NodePos); p_b = TempSp; p_c = QkForceRealBk(NodePos++);			// Set Pointers
				// *MULTIPLY* TempSp and calculation space; BkSp has a copy
				FOR(i,m_iChar) {
					*(p_c++) = *(p_b);			// Copy to BkSp
					*(p_a++) *= *(p_b++);		// Multiply by SpaceTo

	}	}	}	}
	p_a = NULL; p_b = NULL; p_c = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Analytic branch derivative functions
///////////////////////////////////////////////////////////////////////////////////////////////////////
#define SITE 100

// Public function used to prepare calculations
void CBaseProcess::PrepareBraDer()	{
	int i;
	// Get the Q.P(t) matrices required for analytic derivatives
	if(m_pRate->Val() > DX)	{
		Likelihood(true);
		FOR(i,Tree()->NoBra()) {
			if(Tree()->GoodBra(i) == true) { MulMat(m_vpQMat[QMat4Bra(i)]->Q(),PT(i),QP(i),m_iChar,m_iChar,m_iChar); }
	}	}
	// NB: For small rates there is no point getting QP because it is zero
	m_bBraDerReady = true;
}

// Public function to get branch derivatives
bool CBaseProcess::GetBraDer(CProb *ModelL)     {
	int i, BrError = 1;
	// Initialise
	m_arModelL = ModelL;
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n\n--- GetBraDer: " << m_sName << " Rate: " << Rate() << " ---";
#endif
	// Get the derivatives
	FOR(i,Tree()->NoBra()) { Tree()->pBra(i)->grad(0.0); }	// Set them all to zero to start
	if(Tree()->NoBra() == 1) {	// For 2 species trees don't bother trying to be clever
		return false;
	} else if(Rate() < DX) {	// Rate 0 adds nothing to the derivative
		return true;
	} else {					// Otherwise do fast derivative calculations
		Branch_dT(Tree()->StartCalc(),-1,-1,Tree(),-1, &BrError);	// Do the calculations
	}
	// Go through gradients and zero those pointing towards an even more negative branch
	FOR(i,Tree()->NoBra())	{
		// Check the GoodBra
		if(Tree()->GoodBra(i) == false) { continue; }
	}
	// Tidy up
	m_arModelL = NULL;
	// Return if okay
	return FlipBool(m_bFailedL);
}


// In-order tree traversal to calculate the gradients
void CBaseProcess::Branch_dT(int NTo, int NFr, int Branch, CTree *pTree, int First, int *BrError)	{
	int i;
	// Do the first node
	if(NFr == -1) {
		rFOR(i,pTree->NoLinks(NTo)) { if(pTree->NodeLink(NTo,i) == -1) { continue; } Branch_dT(pTree->NodeLink(NTo,i),NTo,pTree->NodeBra(NTo,i),pTree,First,BrError); }
	} else {
	// Otherwise perform the branch_dT
		// Always perform the calculations in the first place
		if(pTree->NodeType(NTo) == leaf || pTree->NodeType(NFr) == leaf)	{ // Do the leaf calculations
			LeafNode_dT(NFr,NTo,Branch,pTree,First, BrError);
		} else if(pTree->NodeType(NTo) == branch) { // Do the internal calculations
			BranNode_dT(NTo,NFr,Branch,pTree,First, BrError);
		} else { Error("\nUnable to calculate derivative in CBaseProcess::Branch_dT for node type\n\n"); }
		// Do the looping
		First = 0;
		FOR(i,pTree->NoLinks(NTo))	{
			if(pTree->NodeLink(NTo,i) == NFr || pTree->NodeLink(NTo,i) == -1) { continue; }
			Branch_dT(pTree->NodeLink(NTo,i),NTo,pTree->NodeBra(NTo,i),pTree,First,BrError);
			First = 1;
	}	}
}

void CBaseProcess::LeafNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError)	{
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\nLeafNode_dT: NTo: " << NTo << "; NFr: " << NFr << "; Br: " << Br;
	cout << "\nEqm: " << m_vpQMat[QMat4Bra(Br)]->Eqm();
	cout << "\nP(t): "; OutPT(cout,Br);
#endif
	bool CharCheck,LeafSeq = true;
	int i,site,Seq, SiteScale = 0;	// Counters
	double Total, Value, *p_a = NULL, *p_b = NULL;
	static double Vec[MAX_CHAR];
	CProb NewL;
	vector <double> eqm = RootEqm();

	// Make sure leaf node is in NTo
	if(NTo > NFr) { i = NFr; NFr = NTo; NTo = i; }
	assert(pTree->NodeType(NFr) == branch); assert(m_arModelL != NULL);
	// Get Seq: the node number that the information goes to
	if(!m_viLeafMap.empty() && m_pSubTree != NULL) {
		assert((int)m_viLeafMap.size() > NTo);
		if(m_viLeafMap[NTo] == -1 || m_viLeafMap[NTo] > m_pData->m_iNoSeq) {
			LeafSeq = false; Seq = NTo;
		} else { LeafSeq=true; Seq = m_viLeafMap[NTo]; }
	} else {
		Seq = NTo;
		if(Seq < m_pData->m_iNoSeq) { LeafSeq = true; } else { LeafSeq = false; }
	}
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	Value = 0.0;
	FOR(site,m_iSize)	{
		// If not a gap then calculations are non-trivial
		/////////////////////////////////////////////////////////////////////
		if(m_pSubTree == NULL || LeafSeq == true)	{ // Gap check statement
			if(m_pData->m_ariSeq[Seq][site] != m_iDataChar) { CharCheck = true; } else { CharCheck = false; }
		} else { CharCheck = true; }
		if(CharCheck == true) {
			////////////////////////////////////////////////////////////////////////////////
			// Get calculation of Vec = leafnode * QP
			SiteScale = 0;
			if(LeafSeq == true) {
				Data2PartL(m_pData->m_ariSeq[Seq][site],QP(Br),Vec,&eqm);
#if ANALYTIC_DERIVATIVE_DEBUG == 1
				if(site < ADD_SITE_MAX) { cout << "\nSite["<<site<<"]: " << m_pData->m_ariSeq[Seq][site] << "\n\tData2PartL: "; FOR(i,m_iChar) { cout << Vec[i] << " "; } }
#endif
			}
			else				{ VMat(ForceRealFd(NTo,site),QP(Br),Vec,m_iChar); SiteScale += *ForceRealFdSc(NTo,site); }
			////////////////////////////////////////////////////////////////////////////////
			// Get calculation of total = sum(Vec[i] = Vec[i] * BranchNode[i] * Eqm[i]);
			Total = 0.0;
			p_a = ForceRealFd(NFr,site);
			FOR(i,m_iChar)	{
				Vec[i] = Vec[i] * *(p_a++) * eqm[i];
				Total += Vec[i];
			}
			SiteScale += *ForceRealFdSc(NFr,site);
			////////////////////////////////////////////////////////////////////////////////
			// No do the derivative calculation
			Value += PartialGrad(site,Total,SiteScale - ModelL(site).Scale());
#if ANALYTIC_DERIVATIVE_DEBUG == 1
			if(site < ADD_SITE_MAX) { cout << "\nTotal = " << Value; }
#endif

	}	}
	// update the space
	if(First != 1) { LeafNode_Update(NTo,NFr,Br,pTree,First); }
	// Clean up
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n\tBranch[" << Br <<"] = " << Tree()->B(Br) << ": Value=" << Value << " Prob()=" << Prob() << " BScale=" << Tree()->BScale(Br) << " -- ret: " << -(Value * Prob() * Tree()->BScale(Br));
#endif
	Tree()->pBra(Br)->grad(-(Value * Prob() * Tree()->BScale(Br)));
	p_a = NULL; p_b = NULL;

}
void CBaseProcess::BranNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError)	{
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\nBranNode_dT: NTo: " << NTo << "; NFr: " << NFr << "; Br: " << Br;
	cout << "\nP(t): "; OutPT(cout,Br);
#endif
	int i,site,SiteScale;	// Counters
	double Value,Total, *p_a = NULL, *p_c = NULL, *p_d = NULL;
	static double Vec[MAX_CHAR];
	vector <double> eqm = RootEqm();
	assert(InRange(NFr,pTree->NoSeq(),pTree->NoNode()) || First == -1);
	assert(pTree->NodeType(NTo) == branch); assert(m_arModelL != NULL);
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	Value = 0.0;
	FOR(site,m_iSize)	{
		SiteScale = 0;
		// Check whether everything below only consists of gaps
		FOR(i,m_iChar)	{ if(diff(ForceRealFd(NFr,site)[i],1.0) == 1) { break; } }
		if(i != m_iChar)	{	// Not a gap
			// Do the derivative calculations
			/////////////////////////////////////////////////////////////////
			switch(First)	{
			case -1:
//				cout << "\nWarning not branch node for first calculation...";
			case 0:
			case 1:
				VMat(ForceRealFd(NFr,site),QP(Br),Vec,m_iChar);
				SiteScale = *ForceRealFdSc(NFr,site);
				break;
			default:
				Error("Unknown first...");
			};
#if ANALYTIC_DERIVATIVE_DEBUG == 1
			vector <int> Left,Right;
			Tree()->BranchSets(Br,&Left,&Right);

			if(site < ADD_SITE_MAX) {
				cout << "\nSite["<<site<<"] pat_occ: " << m_pData->m_ariPatOcc[site] <<"; Left: ";
				FOR(i,(int)Left.size()) { cout << m_pData->m_ariSeq[Left[i]][site]; }
				cout << "; Right: ";
				FOR(i,(int)Right.size()) { cout << m_pData->m_ariSeq[Right[i]][site]; }
				cout << "\n\tOri: ";
				FOR(i,m_iChar) { cout << ForceRealFd(NFr,site)[i] << " "; }
				cout << "\n\tVec: ";
				FOR(i,m_iChar) { cout << Vec[i] << " "; }
			}
#endif
			// For internal branches a complete calculation is required
			SiteScale += *ForceRealFdSc(NTo,site);
			p_a = ForceRealFd(NTo,site);
			FOR(i,m_iChar) { Vec[i] *= *(p_a++); }
			Total = 0; FOR(i,m_iChar)	{ Total += Vec[i] * eqm[i]; }
			// Do the actual calculation
			Value += PartialGrad(site,Total,SiteScale - ModelL(site).Scale());
	} 	}
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n\tBranch[" << Br <<"] = " << Tree()->B(Br) << ": Value=" << Value << " Prob()=" << Prob() << " BScale=" << Tree()->BScale(Br) << " -- ret: " << -(Value * Prob() * Tree()->BScale(Br));
#endif
	Tree()->pBra(Br)->grad(-(Value * Prob() * Tree()->BScale(Br)));
	// Update the space
	BranNode_Update(NTo, NFr, Br, pTree, First);
	// Tidy up
	p_a = NULL; p_c = NULL; p_d = NULL;
}

// Update routines that allow the backward calculations
void CBaseProcess::LeafNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoCompleteUpdate)	{
	bool CharCheck,LeafSeq = true;
	int i,site,Seq, SiteScale = 0;	// Counters
	double Vec[MAX_CHAR], Value, *p_a = NULL, *p_b = NULL;
	CProb NewL;
	vector <double> eqm = m_vpQMat[QMat4Bra(Br)]->Eqm();
	// Make sure leaf node is in NTo
	if(NTo > NFr) { i = NFr; NFr = NTo; NTo = i; }
	// Deal with the case when going to leafnode from StartCalc()
	if(NFr == PartLNode()) { assert(First == -1); First = 1; }
	else { assert(pTree->NodeType(NFr) == branch); }
	// Get Seq: the node number that the information goes to
	if(IsSubTree()) {
		assert((int)m_viLeafMap.size() > NTo);
		if(m_viLeafMap[NTo] == -1 || m_viLeafMap[NTo] > m_pData->m_iNoSeq) {
			LeafSeq = false; Seq = NTo;
		} else { LeafSeq=true; Seq = m_viLeafMap[NTo]; }
	} else {
		Seq = NTo;
		if(Seq < m_pData->m_iNoSeq) { LeafSeq = true; } else { LeafSeq = false; }
	}
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	Value = 0.0;

	// Adjust First so more updating is done if DoCompleteUpdate == true
	if(DoCompleteUpdate == true)	{ if(First == 0) { First = -1; } else if(First == 1) { First = 0; } }

	// Note: You can't speed up these computations the same way as you can speed up forwards
	FOR(site,m_iSize)	{
		// If not a gap then calculations are non-trivial
		/////////////////////////////////////////////////////////////////////
		if(m_pSubTree == NULL || LeafSeq == true)	{ // Gap check statement
			if(m_pData->m_ariSeq[Seq][site] != m_iDataChar) { CharCheck = true; } else { CharCheck = false; }
		} else { CharCheck = true; }
		if(CharCheck == true) {
			// Do the updating procedure
			/////////////////////////////////////////////////////////////////
			// i) Obtain vP(t) in Vec
			SiteScale = 0;
			if(LeafSeq == true) { Data2PartL(m_pData->m_ariSeq[Seq][site],PT(Br),Vec,&eqm); }
			else				{ VMat(ForceRealFd(NTo,site),PT(Br),Vec,m_iChar); SiteScale += *ForceRealFdSc(NTo,site); }
			// ii) Update appropriate memory according to First
			// I think that if first use Fd Space
			switch(First)	{
			case -1:	// If origin branch update Fd to complete partial likelihood
				// Get new Fd Space and back space
				p_a = ForceRealFd(NFr,site); p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar)	{ *(p_a++) = *(p_b) * Vec[i]; *(p_b++) = Vec[i]; }
				*ForceRealFdSc(NFr,site) = SiteScale + *ForceRealBkSc(NFr,site);
				*ForceRealBkSc(NFr,site) = SiteScale;
				break;
			case 0:		// If the first link in an internal node update Fd in NFr to full partial likelihood
				p_a = ForceRealFd(NFr,site); p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar)	{ *(p_a++) = *(p_b++) * Vec[i]; }
				*ForceRealFdSc(NFr,site) = SiteScale + *ForceRealBkSc(NFr,site);
				break;
			case 1:		// If the second link in an internal node
				break;
			default:
				Error("Unknown first...");
			};
		// If a gap then nothing contributed to derivatives and
		//   updates easy because vP(t)_i = 1;
		/////////////////////////////////////////////////////////////////////
		} else {
			// ii) Update appropriate memory according to First
			switch(First)	{
			case -1:	// If origin branch update Fd to complete partial likelihood
				p_a = ForceRealFd(NFr,site); p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar)	{ *(p_a++) = *(p_b++); }
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar)	{ *(p_b++) = 1; }
				*ForceRealBkSc(NFr,site) = 0;
				break;
			case 0:		// If the first link in an internal node update Fd in NFr to full partial likelihood
				// Update code - Correct from 4sp tree
				p_a = ForceRealFd(NFr,site); p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar)	{ *(p_a++) = *(p_b++); }
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				break;
			case 1:		// If the second link in an internal node
				break;
			default:
				Error("Unknown first...");
			};
		}
	}
	DoScale(NFr,true);
}

///////////////////////////////////////////////////////////////////////
// Update routine for branch nodes
// -------------------------------
// This is quite complicated
// Run type 1:
// ===========
// If DoNTo == DoNFr == true, then forward and backward nodes are adjusted at the same time
// This is fine if no other branches change in the tree
// Run type 2:
// ===========
// An alternative way for performing calculations is:
// 1. adjust the To nodes in an in-order tree traversal fashion
// 2. adjust the Fr nodes in a post-order tree traversal fashion
// This means that the function needs to be called at different times, complicating its implementation.
// But it does allow branches to be changed dynamically


void CBaseProcess::BranNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoNTo, bool DoNFr, bool DoCompleteUpdate)	{
	int i,site,SiteScale,SiteScale2;	// Counters
	double *p_a = NULL, *p_b = NULL, *p_c = NULL, *p_d = NULL;
	double Vec[MAX_CHAR], Vec2[MAX_CHAR];
	assert(InRange(NFr,pTree->NoSeq(),pTree->NoNode()) || First == -1);
	assert(pTree->NodeType(NTo) == branch);
//	cout << "\nDoing BranNode_Update...";
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	if(DoNTo && DoNFr)	{
		FOR(site,m_iSize)	{
			// Do the updating procedure
			/////////////////////////////////////////////////////////////////
			// i) Update appropriate memory according to First
			switch(First)	{
			case -1:	// Do the update for the StartCalc() node -- Also requires update of BackSp in StartCalc()
				assert(NFr == pTree->StartCalc());
				// Prepare Bk space for node from
				VMat(ForceRealFd(NTo,site),PT(Br),Vec2,m_iChar);	// Get the vector of partial likelihoods NodeTo -> NodeFr
				SiteScale2 = *ForceRealFdSc(NTo,site);
				// Update the BackSpace
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) *= Vec2[i]; }
				*ForceRealBkSc(NFr,site) += SiteScale2;
				// Update ForwardSpace node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_a = ForceRealFd(NTo,site); p_b = ForceRealBk(NFr,site);
				p_c = ForceRealFd(NFr,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{
					*(p_a++) = *(p_d) * Vec[i];
					*(p_d++) = Vec[i];
					*(p_c++) = *(p_b++);
				}
				*ForceRealFdSc(NTo,site) = *ForceRealBkSc(NTo,site) + SiteScale;
				*ForceRealBkSc(NTo,site) = SiteScale;
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				// Update the BackSpace
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) = Vec2[i]; }
				*ForceRealBkSc(NFr,site) = SiteScale2;
				break;
			case 0:		// If the first link in an internal node update Fd in NFr to full partial likelihood
				// Prepare Bk space for node from
				VMat(ForceRealFd(NTo,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeTo -> NodeFr
				SiteScale = *ForceRealFdSc(NTo,site);
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) *= Vec[i]; }
				*ForceRealBkSc(NFr,site) += SiteScale;
				// Update node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_a = ForceRealFd(NTo,site); p_b = ForceRealBk(NFr,site); p_c = ForceRealFd(NFr,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{
					*(p_a++) = *(p_d) * Vec[i];
					*(p_d++) = Vec[i];
					*(p_c++) = *(p_b++);
				}
				*ForceRealFdSc(NTo,site) = *ForceRealBkSc(NTo,site) + SiteScale;
				*ForceRealBkSc(NTo,site) = SiteScale;
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				break;
			case 1:		// If the second link in an internal node
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Obtain vP(t) in Vec
				SiteScale = *ForceRealFdSc(NFr,site);
				// Do the update
				p_a = ForceRealFd(NTo,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{
					*(p_a++) = Vec[i] * *(p_d);
					*(p_d++) = Vec[i];
				}
				*ForceRealFdSc(NTo,site) = SiteScale + *ForceRealBkSc(NTo,site);
				*ForceRealBkSc(NTo,site) = SiteScale;
				break;
			default:
				Error("Unknown first...");
			};
		}
		DoScale(NFr,true); DoScale(NTo,true);
	} else if(DoNTo)	{
		///////////////////////////////////////////////////////
		// Only update NTo
		FOR(site,m_iSize)	{
			// Do the updating procedure
			/////////////////////////////////////////////////////////////////
			// i) Update appropriate memory according to First
			switch(First)	{
			case -1:	// Do the update for the StartCalc() node -- Also requires update of BackSp in StartCalc()
				assert(NFr == pTree->StartCalc());
				// Update ForwardSpace node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_a = ForceRealFd(NTo,site); p_b = ForceRealBk(NFr,site);
				p_c = ForceRealFd(NFr,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{
					*(p_a++) = *(p_d) * Vec[i];
					*(p_d++) = Vec[i];
				}
				*ForceRealFdSc(NTo,site) = *ForceRealBkSc(NTo,site) + SiteScale;
				*ForceRealBkSc(NTo,site) = SiteScale;
				// Update the BackSpace
				break;
			case 0:		// If the first link in an internal node update Fd in NFr to full partial likelihood
				// Update node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_a = ForceRealFd(NTo,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{
					*(p_a++) = *(p_d) * Vec[i];
					*(p_d++) = Vec[i];
				}
				*ForceRealFdSc(NTo,site) = *ForceRealBkSc(NTo,site) + SiteScale;
				*ForceRealBkSc(NTo,site) = SiteScale;
				break;
			case 1:		// If the second link in an internal node
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Obtain vP(t) in Vec
				SiteScale = *ForceRealFdSc(NFr,site);
				// Do the update
				p_a = ForceRealFd(NTo,site); p_d = ForceRealBk(NTo,site);
				FOR(i,m_iChar)	{ *(p_a++) = Vec[i] * *(p_d); *(p_d++) = Vec[i]; }
				*ForceRealFdSc(NTo,site) = SiteScale + *ForceRealBkSc(NTo,site);
				*ForceRealBkSc(NTo,site) = SiteScale;
				break;
			default:
				Error("Unknown first...");
			};
		}
		DoScale(NTo,true);
	} else {
		if(DoCompleteUpdate == true) { if(First == 0) { First = -1; } else if(First == 1) { First = 0; } }
		///////////////////////////////////////////////////////
		// Only update NFr
		FOR(site,m_iSize)	{
			// Do the updating procedure
			/////////////////////////////////////////////////////////////////
			// i) Update appropriate memory according to First
			switch(First)	{
			case -1:	// Do the update for the StartCalc() node -- Also requires update of BackSp in StartCalc()
				// Prepare Bk space for node from
				VMat(ForceRealFd(NTo,site),PT(Br),Vec2,m_iChar);	// Get the vector of partial likelihoods NodeTo -> NodeFr
				SiteScale2 = *ForceRealFdSc(NTo,site);
				// Update the BackSpace
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) *= Vec2[i]; }
				*ForceRealBkSc(NFr,site) += SiteScale2;
				// Update ForwardSpace node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_b = ForceRealBk(NFr,site); p_c = ForceRealFd(NFr,site);
				FOR(i,m_iChar)	{ *(p_c++) = *(p_b++); 	}
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				// Update the BackSpace
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) = Vec2[i]; }
				*ForceRealBkSc(NFr,site) = SiteScale2;
				break;
			case 0:		// If the first link in an internal node update Fd in NFr to full partial likelihood
				// Prepare Bk space for node from
				VMat(ForceRealFd(NTo,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeTo -> NodeFr
				SiteScale = *ForceRealFdSc(NTo,site);
				p_b = ForceRealBk(NFr,site);
				FOR(i,m_iChar) { *(p_b++) *= Vec[i]; }
				*ForceRealBkSc(NFr,site) += SiteScale;
				// Update node to and finish updating node from
				VMat(ForceRealFd(NFr,site),PT(Br),Vec,m_iChar);	// Get the vector of partial likelihoods NodeFr -> NodeTo
				SiteScale = *ForceRealFdSc(NFr,site);
				p_b = ForceRealBk(NFr,site); p_c = ForceRealFd(NFr,site);
				FOR(i,m_iChar)	{ *(p_c++) = *(p_b++); }
				*ForceRealFdSc(NFr,site) = *ForceRealBkSc(NFr,site);
				break;
			case 1:		// If the second link in an internal node
				break;
			default:
				Error("Unknown first...");
			};
		}
		DoScale(NFr,true);
	}
}

double CBaseProcess::PartialGrad(int site,double Total,int SiteScale)	{

	if(m_pData->IsBootstrap()) { if(m_pData->m_ariPatOcc[site] == 0) { return 0.0; } }
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	if(site < ADD_SITE_MAX) {
		cout << "\n\t\tDoing site: " << site << " [occ=" << m_pData->m_ariPatOcc[site]<<"]: Total= " << Total << "*10^" << SiteScale << "; Partial= " << ModelL(site) << "; return value: " << (Total / ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site] * pow((double)10,-SiteScale) );
	}
#endif
	if(fabs(Total) < FLT_EPSILON) { return 0; }
	if(abs(SiteScale) < 15)	{	// If scaling okay, then do normal calculations
		if(SiteScale == 0)	{
			if(fabs(ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site]) < DBL_EPSILON) { m_bFailedL = true; return -BIG_NUMBER; }
			return (Total / ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site]);
		} else {
			if(fabs(ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site]) < DBL_EPSILON) {
//				cout.precision(16);	cout << "\nFailed site["<<site<<"]: numerator: " << Total << ", denominator: " << ModelL(site).ScalVal() << " * " << m_pData->m_ariPatOcc[site];
				m_bFailedL = true; return -BIG_NUMBER; }
			return (Total / ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site] * pow((double)10,-SiteScale) );
	}	} else {					// Do extreme values
		if(SiteScale > 0)	{ return 0.0; }
		else				{
			if(Total / ModelL(site).ScalVal() > 0) { return GRAD_LIM; } else { return -GRAD_LIM; }
	}	}
	return 0;
}


// Function that calculates Node[NTo] * PT(Br), then multiplies elementwise by Node[NFr]
// --
// Note the current implementation works at the double level, whereas normal likelihood function works with CProbs.
// This difference is a potential source of numerical instability.

void CBaseProcess::GetBranchPartL(CProb **arpP, int NT, int NF, int B)	{
//	cout << "\n--> Entered CBaseProcess::GetBranchPartL(CProb **arpP, int NT, int NF, int B)";
	bool NTreal = true, NFreal = true;
	int i,SiteScale = 0,site;
	double *p_a = NULL, Total = 0.0;
	static double V[MAX_SPACE];
	CProb Pr;
	vector <double> eqm = RootEqm();
//	cout << "\nIn CBaseProcess::GetBranchPartL: eqm: " << eqm;
	// Do Garbage Collector rate
	//////////////////////////////////////////////
	if(MaxRate()) {
		FOR(site,m_pData->m_iSize) {
			Pr.Assign(1.0);
			FOR(i,NoSeq()) {
				if(m_pData->m_ariSeq[i][site] == m_iChar) { continue; }
				Pr.Multiply(eqm[m_pData->m_ariSeq[i][site]],true);

			}
			Pr.Multiply(Prob(),true);
			arpP[site]->Add(Pr,true);
	}	}
	//////////////////////////////////////////////
	// Do normal processes
	else if(Rate() > DX)	{
		m_pTree->UpdateB(B); Make_PT(B,true);

		// Organise whether NT and NF are real sequences
		if(IsSubTree()) { // Deal with subtrees
			if(NT < Tree()->NoSeq()) { if(m_viLeafMap[NT] == -1) { NTreal = false; } else { NT = m_viLeafMap[NT]; } } else { NTreal = false; }
			if(NF < Tree()->NoSeq()) { if(m_viLeafMap[NF] == -1) { NFreal = false; } else { NF = m_viLeafMap[NF]; } } else { NFreal = false; }
		} else { if(NT >= Tree()->NoSeq()) { NTreal = false; } if(NF >= Tree()->NoSeq()) { NFreal = false; } }
		//////////////////////////////////////////////////////
		// Do the calculations
		FOR(site,m_pData->m_iSize)	{
			SiteScale = 0; Total = 0.0;
			// Get first vector of calc
			if(NTreal)	{
				Data2PartL(m_pData->m_ariSeq[NT][site],PT(B),V,&eqm);
//				FOR(i,m_iChar) { if(i == m_pData->m_ariSeq[NT][site] || m_pData->m_ariSeq[NT][site] == m_iChar) { cout << "\t1"; } else { cout << "\t0"; } }
			} else {
//				FOR(i,m_iChar)	{ cout << "\t" << ForceRealFd(NT,site)[i]; }
				VMat(ForceRealFd(NT,site),PT(B),V,m_iChar); SiteScale += *ForceRealFdSc(NT,site);
			}
//			int j; 	cout << "\n\tRight:\t"; if(NFreal) { FOR(i,m_iChar) { if(i == m_pData->m_ariSeq[NF][site] || m_pData->m_ariSeq[NF][site] == m_iChar) { cout << "\t1"; } else { cout << "\t0"; } } } else { FOR(i,m_iChar)	{ cout << "\t" << ForceRealFd(NF,site)[i]; } } cout << "\n\tLeft * P(t):"; FOR(j,m_iChar) { cout << "\t" << V[j]; }

			// Get calculation of total = sum(Vec[i] = Vec[i] * BranchNode[i] * Eqm[i]);
			if(NFreal)	{
				Total = Sum_Vec(m_pData->m_ariSeq[NF][site],V,eqm);
			} else {	// Do partial likelihoods
				p_a = ForceRealFd(NF,site);
				FOR(i,m_iChar)	{ Total += V[i] * *(p_a++) * eqm[i]; }
				SiteScale += *ForceRealFdSc(NF,site);
			}
			Total *= Prob();

			// Assign the likelihood
			Pr.Assign(Total,SiteScale);
			arpP[site]->Add(Pr,true);
	}	} else {
	/////////////////////////////////////////////////
	// Do zero rate processes
		FOR(site,m_pData->m_iSize)	{
			if(m_pData->m_viNoChange[site] != -1) {
				Pr.Assign(eqm[m_pData->m_viNoChange[site]] * Prob());
				arpP[site]->Add(Pr,true);
	}	}	}
#if DEVELOPER_BUILD == 1
		cout << "\nFinished with Eqm: " << eqm;
#endif

}

///////////////////////////////// Data to partial likelihood adapter functions /////////////////
// 1. Simple routine that multiplies a leaf node by its PT matrix

/////////////////////////////////////////////////////////////////////////////////
// Functions translating data into the vectors used by model
void CBaseProcess::Data2PartL(int Char,double *PT, double *RetSpace, vector <double> *eqm)	{
	int i;
	double *p_PT = &PT[Char], *p_Ret = RetSpace;
	assert(Char >=0 && Char <= m_iDataChar);
	// If a gap then set the multiple to one
	if(Char == m_iDataChar)	{
		FOR(i,m_iChar) { *(p_Ret++) = 1.0; }
	} else {
		if(m_iHiddenChar == 1) { FOR(i,m_iChar) { *(p_Ret++) = *(p_PT); p_PT += m_iChar; } }
		else {
			int state;
			FOR(i,m_iChar) { *(p_Ret++) = 0.0; }
			// If a real character
			FOR(state,m_iHiddenChar)	{
				p_PT = &PT[(state * m_iDataChar)+Char];
				p_Ret = RetSpace;
				FOR(i,m_iChar) {
					*(p_Ret++) += *(p_PT); p_PT += m_iChar;
	}	}	}	}
	p_PT = NULL; p_Ret = NULL;
	return;
}
void CBaseProcess::Vec_by_Data(int Char,double *V)	{
	int i;
	if(Char == m_iDataChar) { return; }
	double Val = V[Char];
	FOR(i,m_iChar) { V[i] = 0.0; }
	for(i=Char;i<m_iChar;i+=m_iDataChar) { V[i] = Val; }
}

double CBaseProcess::Sum_Vec(int Char, double *Vec, vector <double> eqm)	{
	int i;
	double Total = 0.0;
	assert(InRange(Char,0,m_iChar+1));
	// Do the simple case
	if(Char != m_iDataChar)	{ for(i=Char;i<m_iChar;i+=m_iDataChar) { Total += Vec[i] * eqm[i]; } }
	// Deal with gaps
	else { FOR(i,m_iChar) { Total += Vec[i] * eqm[i]; } }
	return Total;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Functions for preparing partial likelihoods from centre point
/////////////////////////////////////////////////////////////////////////////////////////////
//
// Takes the current tree and
// Produces partial likelihoods
//
/////////////////////////// Functions controlling the partial likelihood mapping

// Function that prepares Nodes from a CP and their partial likelihoods
// Assumes that the ::PrepareLikelihood(...) function has been called to prepare the model
void CBaseProcess::ApplyCPMapping(vector <int> LeafMap, vector <int> NodeFr, bool DoPartial)	{
	int i;
	vector <int> Map;
	// Check the entry conditions
	assert(m_bSubTreeActive == false);
	// Clean up anything remaining
	CleanCPMapping();
	// Do the partial likelihoods
	// 1.) For pairwise calcs don't worry too much
	if(LeafMap.size() < 3)	{ m_viLeafMap = LeafMap; }
	// 2.) For proper trees
	else	{
		// Get the correct P(t) matrices
		FOR(i,(int)LeafMap.size())	{
			if(LeafMap[i] < m_pTree->NoSeq()) { Map.push_back(LeafMap[i]); continue; }	// If its a normal sequence then continue
			if(DoPartial) { PreparePartialL(LeafMap[i],NodeFr[i],i); }
			else {	// Copy directly from the node
				CopyNode(LeafMap[i],i);
			}
			Map.push_back(-1);
		}
		m_viLeafMap = Map;
}	}

// Partial likelihood calculator
// Note that the way PartialL is set up, it writes to the NodeFrom too...
void CBaseProcess::PreparePartialL(int Node, int NodeFrom, int Node2Trans)	{
	int i;
	// Check entry conditions
	assert(InRange(Node,0,m_pTree->NoNode()));
	assert(InRange(NodeFrom,0,m_pTree->NoNode()));
	// Decompress the space; going to need to anyway!
	DecompressSpace();
	// Adjust the rate for the process (if required)
	FOR(i,(int)m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); }
	// Make the PT matrices
	FOR(i,Tree()->NoBra()) { Tree()->UpdateB(i); Make_PT(i); }
	// Calculate the partial likelihood
	m_bDoingPartial = true;
	PartialL(m_pTree,Node,NodeFrom,m_pTree->FindBra(Node,NodeFrom),true,true);
	DoScale(Node);
	m_bDoingPartial = false;
	// If required copy the partial likelihood to new space
	if(Node2Trans != -1)	{ CopyNode(Node,Node2Trans); }
	// Set flags
	m_bSubTreeActive = true;
}

// Cleaning function
void CBaseProcess::CleanCPMapping()	{
	if(m_bSubTreeActive == true) {
		CleanSubTree();
		m_viLeafMap.clear();
		m_bSubTreeActive = false;
}	}

// Uses a specific tree for calculations
// TODO: The following two routines will be altered to allow multiple Q matrix processes
void CBaseProcess::ApplySubTree(CTree *Tree)	{
	// Entry conditions
	if(m_pSubTree != NULL) { CleanSubTree(); }
	if(m_vpQMat.size() > 1) { Error("\nDon't know how to deal with multiple Q matrices in process for Subtree calculations..."); }
	assert(Tree->NoSeq() == m_viLeafMap.size());
	// Put in the subtree
	m_pSubTree = Tree;
	// Deal with space if required
	if(Tree->NoNode() > MainTree()->NoSeq()) { DecompressSpace(); }
}


// Extracts a subtree from data
void CBaseProcess::CleanSubTree()	{
	m_pSubTree = NULL;			// Remove the subtree
	if(m_vpQMat.size() > 1) { Error("\nDon't know how to deal with multiple Q matrices in process for Subtree calculations..."); }
}

// MainTree function
CTree *CBaseProcess::MainTree(CTree *T) {
	if(T != NULL) {
		if(m_pTree == NULL) { m_pTree = T; }
		else				{ *m_pTree = *T; }
	}
	return m_pTree;
}

////////////////////////////////////////////////////////////////////////////////
// Debug function for CBaseProcess
vector <double> CBaseProcess::GetQMatSums()	{
	int i,j,k;
	double *p,total;
	vector <double> Sums;
//	cout << "\nCBaseProcess::GetQMatSums() for process " << m_sName;
	FOR(i,(int)m_vpQMat.size())	{
//		cout << "\n\tQMat["<<i<<"]: ";
		total = 0; p = m_vpQMat[i]->Q();
		FOR(j,m_iChar) {
			FOR(k,m_iChar) {
				if(j==k) { *p++; continue; }
				total += *(p++);
		}	}
		Sums.push_back(total);
//		cout << total;
	}
	return Sums;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions dealing with the particulars of models
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add kappa to a model.
CQPar * CBaseProcess::Kappa(EDataType Type)	{
	int i,j;
	CQPar *Par;
	double KappaVal;
#if ALLOW_KAPPA_GUESS == 1
	KappaVal = m_pData->GuessKappa();
#else
	KappaVal = INITIAL_KAPPA;
#endif
	// Check entry conditions
	assert(Type == DNA || Type == COD || Type == COD_RED);
	Par = new CQPar("Kappa",m_iChar,KappaVal);
	// Assign the parameters
	FOR(i,m_iChar)	{ for(j=i+1;j<m_iChar;j++)	{
		if(IsTs( m_sABET.substr(i*LenStates(Type),LenStates(Type)) , m_sABET.substr(j*LenStates(Type),LenStates(Type)),Type)) { Par->AddQij(i,j); } } }
	m_vpPar.push_back(Par);
	return Par;
}

vector <CQPar *> CDNAProcess::MakeDNAREV(EDataType Type)	{
	CQPar *Par = NULL;
	vector <CQPar *> vPars;
	int i,j;
	double KappaVal = INITIAL_KAPPA;
#if ALLOW_KAPPA_GUESS == 1
	KappaVal = m_pData->GuessKappa();
#endif
	string Name;
	FOR(i,m_iChar)	{
		for(j=i+1;j<m_iChar;j++)	{
			Name = GetPos(m_sABET,i,m_iABET_length/m_iChar) + "<->" + GetPos(m_sABET,j,m_iABET_length/m_iChar);
			if(IsTs(State(Type,i),State(Type,j),Type))	{ Par = new CQPar(Name,m_iChar,KappaVal); }
			else										{ Par = new CQPar(Name,m_iChar,1.0); }
			Par->AddQij(i,j);
			m_vpPar.push_back(Par);
			vPars.push_back(Par);
			Par = NULL;
	}	}

	vPars[vPars.size() - 1]->SetOptimise(false);

	return vPars;
}

// Add a simple equilibrium process to a model
void CBaseProcess::AddSimpleEqm(bool Opt)	{ AddSimpleEqm(m_pData->m_iChar,m_pData->m_vFreq,Opt); }

void CBaseProcess::AddSimpleEqm(int Char, vector <double> Freq, bool Opt)	{
	CSimpleEqm *Eqm;
	Eqm = new CSimpleEqm(Char,&m_vpPar,Freq);
	m_vpEqm.push_back(Eqm);
	Eqm->SetOpt(Opt);
	Eqm = NULL;
}

//////////////////////////////////////////////////////////////////////////
// Covarion style process character (observable) distribution functions

// Add a simple eqm to a Covarion process
void CBaseProcess::AddSimpleCovEqm(vector <CPar *> *StateProbs,bool ForceEqu, bool Opt)	{
	AddSimpleCovEqm(m_pData->m_iChar * (int)StateProbs->size(),m_pData->m_vFreq, StateProbs,ForceEqu,Opt); }
void CBaseProcess::AddSimpleCovEqm(int Char, vector <double> Freq, vector <CPar *> *StateProbs, bool ForceEqu, bool Opt)	{
	int i;
	assert(Freq.size() * StateProbs->size() == Char);
	if(ForceEqu == true) { FOR(i,(int)Freq.size()) { Freq[i] = 1.0 / (double) Freq.size(); }  }
	CCovEqm *Eqm;
	Eqm = new CCovEqm(Char,&m_vpPar,Freq,StateProbs);
	Eqm->SetOpt(Opt);
	m_vpEqm.push_back(Eqm);
	Eqm = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////
// Functions for adding a complex eqm to a Covarion process
// --------------------------------------------------------
// 1. Function that adds a seperate equilibrium distribution for each hidden state. These are taken from the data
void CBaseProcess::AddComplexCovEqm(int Char, std::vector<CPar*> *StateProbs, bool AllowOpt)	{
	int i, j, HiddenChar = (int) StateProbs->size();
	vector <vector <double> > RealFreq,Frq;
	vector <double> Temp(DataChar(),0.0);
	if(AllowOpt == false) { Error("\nTrying to run CBaseProcess::AddComplexCovEqm(int Char, std::vector<CPar*> *StateProbs, bool AllowOpt) with AllowOpt == false.\nNot a good idea because the routine isn't well thought out...\n"); }
	if(HiddenChar >= m_pData->m_iNoSeq) { Error("\nWarning: creating too many equilibrium distributions (" + int_to_string(HiddenChar) + ") from " + int_to_string(m_pData->m_iNoSeq) + " sequences"); }
	FOR(i,m_pData->m_iNoSeq)	{ RealFreq.push_back(m_pData->GetFreq(i)); }
	// Get some starting frequencies
/*	FOR(i,HiddenChar)	{
		FOR(j,DataChar())	{
			Temp[j] = RealFreq[RandInt(0,NoSeq())][j];
		}
		NormaliseVector(Temp);
		Frq.push_back(NormaliseVector(Temp));
	}
*/
	FOR(i,HiddenChar)	{
		FOR(j,DataChar())	{
			Temp[j] = m_pData->m_vFreq[j];
		}
		NormaliseVector(Temp);
		Frq.push_back(NormaliseVector(Temp));
	}
	AddComplexCovEqm(Char,Frq,StateProbs,AllowOpt);
}

// 2. Function that adds a seperate equilibrium distribution for each hidden state, taken from Freq.
void CBaseProcess::AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <CPar *> *StateProbs, bool AllowOpt)	{
	int i,HiddenChar = (int) StateProbs->size();
	vector <vector <int> > FreqMap;	// Set up freq map
	vector <int> FM1;
	FOR(i,HiddenChar)	{ FM1.push_back(i); FreqMap.push_back(FM1); FM1.clear(); }
	AddComplexCovEqm(Char,Freq,FreqMap,StateProbs,AllowOpt);
}
// Variables:
//	Char = Number of characters in process
//	Freq = a vector of vectors that contain the equilibrium distribution of observable states.
void CBaseProcess::AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <vector <int> > FreqMap,vector <CPar *> *StateProbs, bool AllowOpt)	{
	CCovEqm *Eqm;
	Eqm = new CCovEqm(Char,&m_vpPar,Freq,FreqMap,StateProbs,AllowOpt);
	m_vpEqm.push_back(Eqm);
	Eqm = NULL;
}

// Functions for creating aggregate matrices, useful for checking investigating property and error checking
vector <double> CBaseProcess::AggregateQ(int Position) {
	int i,j,count = 0;
	Position =1;
	string SUB_ABET, temp;
	EDataType SubType;
	int OneChar = pow((int)m_iChar,(double) 1/LenStates(m_DataType));
	vector <double> Q(OneChar*OneChar,0);
	// Do entry checking
	assert(InRange(Position,0,2));
	if(OneChar == m_iChar) { Error("\nTrying to aggregate a single character process"); }
	// Get SUB_ABET
	switch(m_DataType) {
	case RY2: SUB_ABET = RY_ABET; SubType = RY; break;
	case DNA2: SUB_ABET = DNA_ABET; SubType = DNA; break;
	case AA2: SUB_ABET = AA_ABET; SubType = AA; break;
	case COD: SUB_ABET = DNA_ABET; SubType = DNA; break;
	default:
		Error("\nTrying to do aggregated process of something unexpected...\n");

	}
	// Create Q Matrices
	PrepareQMats();
	// Aggregate
	cout << "\nAggregating matrix for position["<<Position<<"]: ";
	FOR(i,m_iChar) {
		FOR(j,m_iChar) {
			if(i==j) { continue; }
			if(State(m_DataType,i)[Position] == State(m_DataType,j)[Position]) { continue; }
			cout << "\n\t" << count++ << "["<<i<<","<<j<<"] " << State(m_DataType,i) << " <-> " << State(m_DataType,j);
			cout << "   " << State(m_DataType,i)[Position] << " cf. " << State(m_DataType,j)[Position];
			temp = State(m_DataType,i)[Position];
			cout << " ::: Pos[" << FindState(SubType,temp);
			temp = State(m_DataType,j)[Position];
			cout  << "][" << FindState(SubType,temp) << "]";
			Q[(FindState(SubType,temp)*OneChar) + FindState(SubType,temp)] += *m_vpQMat[0]->Q(i,j);
		}
	}

	cout << "\nMatrix flat: " << Q;
	cout << "\nAnd the Matrix...";
	MatOut(OneChar,Q);
	exit(-1);
}



////////////////////////////////////////////////////////////////////////////////
// Amino acid model constructors

const bool DoEMPoutput = false;

void CAAProcess::CreateEMPmodel(double *S_ij,double *Freq,bool AddF)	{
	int i,j,count = 0;
	double Max = -1;
	string Name;
	vector <double> F;
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
	if(AddF == true) {
		if(DoEMPoutput) { FOR(i,20) { F.push_back(m_pData->m_vFreq[i]); } }
		AddSimpleEqm();
	}
	else { FOR(i,20) { F.push_back(Freq[i]); } AddSimpleEqm(20,F); }
	if(DoEMPoutput) { cout << "\nThe eqm: " << F << "\nSum == " << Sum(&F) << endl;  }
	PrepareQMats();
	if(DoEMPoutput) { cout << "\nActual:  " << m_vpQMat[0]->Eqm(); exit(-1); }
	m_vpQMat[0]->Lock();
}

// Function for making the equiprobable model
void CAAProcess::MakeEQU(bool AddF)	{
	if(AddF == true) { AddSimpleEqm(); }
	PrepareQMats();
	m_vpQMat[0]->Lock();

}

// Function for making the WAG model
void CAAProcess::MakeWAG(bool AddF)	{
	CreateEMPmodel((double*) dWAGVal,(double*)dWAGFreq,AddF);
}

// Function for making the JTT model
void CAAProcess::MakeJTT(bool AddF)	{
	CreateEMPmodel((double*) dJTTVal,(double*) dJTTFreq,AddF);
}

// Function for making the Dayhoff model
void CAAProcess::MakeDAY(bool AddF)	{
	CreateEMPmodel((double*) dDAYVal,(double*) dDAYFreq,AddF);
}

void CAAProcess::MakeMTREV(bool AddF) {
	CreateEMPmodel((double*) dmtREVVal,(double*) dmtREVFreq,AddF);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Some generic process constructors
////////////////////////////////////////////////////////////////////////////////////////////////////////

/* *************************** Basic DNA processes ************************************ */
CDNAProcess::CDNAProcess(CData *Data, CTree *Tree, DNAProc Model, string name) : CBaseProcess(Data,Tree)	{
	assert(Data != NULL);

	// Firstly deal with the RY model
	if(Data->m_DataType == RY) {
		if(Model != pRY) { Error("\nTrying to create a normal DNA model with RY data...\n\n"); }
		m_DataType = RY;
		MakeBasicSpace(RY);
		m_sName = "RY model";
		AddSimpleEqm();
		if(name.size() < 1) { Add_QMat("RY proc",RY); } else { Add_QMat(name,RY); }
		return;
	}

	if(Data->m_DataType != DNA)	{ Error("Trying to initialise DNA model with data that doesn't look like nucleotides...\n\n"); }

	m_DataType = DNA;
	MakeBasicSpace(DNA);

	// Define the model
	switch(Model)	{
	case pJC:
		m_sName = "JC";
		if(name.size() < 1) { Add_QMat("JC",DNA); } else { Add_QMat(name,DNA); }
		break;
	case pK2P:
		m_sName = "K2P";
		if(name.size() < 1) { Add_QMat("K2P",DNA); } else { Add_QMat(name,DNA); }
		Kappa(DNA);
		break;
	case pFEL:
		m_sName = "FEL";
		if(name.size() < 1) { Add_QMat("FEL",DNA); } else { Add_QMat(name,DNA); }
		AddSimpleEqm();
		break;
	case pHKY:
		m_sName = "HKY";
		if(name.size() < 1) { Add_QMat("HKY",DNA); } else { Add_QMat(name,DNA); }
		Kappa(DNA);
//		Kappa(DNA)->SetVal(3.0);
//	int i;FOR(i,4) {	Data->m_vFreq[i] = (double)(i+1)/10; } cout << "\nFreq: " << Data->m_vFreq;
		AddSimpleEqm();
		break;
	case pREV:
		m_sName = "REV";
		if(name.size() < 1) { Add_QMat("REV",DNA); } else { Add_QMat(name,DNA); }
		MakeDNAREV(DNA);
		AddSimpleEqm();
		break;
	default:
		Error("Attempting to create DNA process from unknown model...");
		break;
	}
}

/* *************************** Basic amino acid processes ************************************ */

CAAProcess::CAAProcess(CData *D, CTree *T, AAProc Model, bool AddF) : CBaseProcess(D,T) {
	assert(D != NULL); assert(D->m_DataType == AA);
	if(D->m_DataType != AA) { Error("Trying to initialise AA model with data that doesn't look like amino acids...\n\n"); }
	MakeBasicSpace(AA);
	m_DataType = AA;
	// Define the model
	switch(Model)	{
	case pEQU:
		Add_QMat("EQU",AA);
		MakeEQU(AddF);
		break;
	case pJTT:
		Add_QMat("JTT",AA);
		MakeJTT(AddF);
		break;
	case pWAG:
		Add_QMat("WAG",AA);
		MakeWAG(AddF);
		break;
	case pDAY:
		Add_QMat("DAY",AA);
		MakeDAY(AddF);
		break;
	case pMTREV:
		Add_QMat("mtREV",AA);
		MakeMTREV(AddF);
		break;
	default:
		Error("Attempting to create AA process from unknown model...");
		break;
	}
}

CAAProcess::CAAProcess(CData *D, CTree *T, string Name, bool AddF, double *S_ij, double  *pi_j) : CBaseProcess(D,T) {
	assert(D != NULL); assert(D->m_DataType == AA);
	if(D->m_DataType != AA) { Error("Trying to initialise AA model with data that doesn't look like amino acids...\n\n"); }
	MakeBasicSpace(AA);
	m_DataType = AA;
	m_sName = Name;
	Add_QMat(Name,AA);
	CreateEMPmodel(S_ij,pi_j,AddF);
}

/* ******************************* Basic codon processes ************************************** */
CCodonProcess::CCodonProcess(CData *D, CTree *T, CodonProc Model, ECodonEqm CE, int GenCode, string RadicalFile) : CBaseProcess(D,T,"CodonProcess")	{
	int i, j, k, count;
	CQPar *Par = NULL;
	string sFrom, sTo,Name;
	vector <int> RadMat(20*20,-1);
	assert(D != NULL); assert(D->m_DataType == COD_RED);
	MakeBasicSpace(D->m_iChar); m_sABET = D->m_sABET;
	// Store the genetic code
	m_iGenCode = GenCode;
	// Do the model equilibrium distribution
	AddCodonEqm(GenCode,m_pData->m_iChar,CE,false);
	// Define the model parameters
	cout << "\nAdding codon process" << flush;

	if(Model == pM0DrDc) {
		cout << "\nMaking random Radical matrix -> <Random.mat>";
		ofstream outrand("Random.mat");
		FOR(i,20) {
			outrand << "\n";
			FOR(j,i) { outrand << RandInt(0,2) << "  "; }
		}
		outrand.close();
		// Input RadicalFile
		cout << "\nInputting Radical amino acids from file: <" << RadicalFile << ">";
		FINOPEN(Radin, RadicalFile.c_str());
		FOR(i,20)	{
			FOR(j,i)	{
				Radin >> RadMat[(i*20)+j];
				if(!InRange(RadMat[(i*20)+j],0,2)) { cout << "\nError reading Radical Matrix from: " << RadicalFile << " at ["<<i << "," << j << "] = " << RadMat[(i*20)+j] << "\nMatrix so far: " << MatOut(20,RadMat); exit(-1); }
				RadMat[(j*20)+i] = RadMat[(i*20)+j];
			}
		}
		Radin.close();
//		cout << "\nRadical Matrix" << endl <<  MatOut(20, RadMat);
//		cout << "\n\nDone";

	}

	switch(Model)	{
	case pM0:
		m_sName = "CodonM0";
		Add_CodRedQMat("M0",D->m_iChar);
		AddMultiChangeZeros(GenCode);
		AddOmega(GenCode);
		Kappa(COD);
		break;
	case pCodonEMPRest:
		cout << "\nMaking RESTRAINED empirical codon model";
		m_sName = "CodonEMPRest";
		Add_CodRedQMat("CodonEMP",D->m_iChar);
		count = 0;
		FOR(i,m_iChar) {
			FOR(j,i)	{
				// Get the names of the states
				sFrom = sTo = "";
				FOR(k,3) { sFrom = sFrom + m_sABET[(j*3)+k]; sTo = sTo + m_sABET[(i*3)+k]; }
				Name = sFrom + "<->" + sTo;

//					cout << "\nCreating transition " << Name << " == " << dECMrest[count] << flush;

				Par = new CQPar(Name,m_iChar,dECMrest[count++],false,0.0,BIG_NUMBER);
					Par->AddQij(i,j);
				Par->SetOptimise(false);
				m_vpPar.push_back(Par);
					Par = NULL;
		}	}
		assert(count == 1830);
		break;
	case pCodonEMPUnrest:
		cout << "\nMaking UNRESTRAINED empirical codon model" << flush;
		m_sName = "CodonEMPUnrest";
		Add_CodRedQMat("CodonEMP",D->m_iChar);
		count = 0;
		FOR(i,m_iChar) {
			FOR(j,i)	{
				// Get the names of the states
				sFrom = sTo = "";
				FOR(k,3) { sFrom = sFrom + m_sABET[(j*3)+k]; sTo = sTo + m_sABET[(i*3)+k]; }
				Name = sFrom + "<->" + sTo;
//					cout << "\nCreating transition " << Name << " == " << dECMunrest[count] << flush;
				Par = new CQPar(Name,m_iChar,dECMunrest[count++],false,0.0,BIG_NUMBER);
					Par->AddQij(i,j);
				Par->SetOptimise(false);
				m_vpPar.push_back(Par);
					Par = NULL;
		}	}
		assert(count == 1830);
		break;
	case pAAEMP:	// Amino acid empirical codon model
		Add_CodRedQMat("CodonAAEMP",D->m_iChar);
		FOR(i,m_iChar) {
			FOR(j,i)	{
				// Get the names of the states
				sFrom = sTo = "";
				FOR(k,3) { sFrom = sFrom + m_sABET[(j*3)+k]; sTo = sTo + m_sABET[(i*3)+k]; }
				Name = sFrom + "<->" + sTo;
//					cout << "\nCreating transition " << Name << " == " << dECMunrest[count] << flush;
				Par = new CQPar(Name,m_iChar,dECMunrest[count++],false,0.0,BIG_NUMBER);
					Par->AddQij(i,j);
				Par->SetOptimise(false);
				m_vpPar.push_back(Par);
					Par = NULL;
		}	}
		assert(count == 1830);
	case pM0DrDc:
		cout << "\nMaking new Dr/Dc matrix with input from <"<<RadicalFile<<">!";
		m_sName = "M0_DrDc";
		// Read in Radical Mat


		Add_CodRedQMat("M0_DrDc",D->m_iChar);
		AddMultiChangeZeros(GenCode);
		AddDrDcOmega(GenCode,RadMat,0);
		AddDrDcOmega(GenCode,RadMat,1);
		Kappa(COD);
		break;

	default:
		Error("\nAttempting to create Codon process from unknown model...");
		exit(-1);
	}
//	PrepareLikelihood(true,true);
//	CreatePTMats();
//	cout << "\nMaking P(t) matrices and yields Q Mat of:";
//	OutQ();
//	cout << "\n\/\/";

}
/////////////////////////////////////////////////////////////////////////////////////////
// Function for adding omega in M0 type manner
CQPar * CCodonProcess::AddOmega(int GenCode)	{
	int i,j, CurChar_i, CurChar_j;
	CQPar *Par;
	// Check entry conditions
	assert(InRange(GenCode,0,11));
	// Initialise
	Par = new CQPar("Omega",m_iChar,INITIAL_OMEGA,true,MIN_OMEGA);
	// Assign the parameters
	if(m_iChar == 64) { // When the Genetic code hasn't been applied
		assert(m_vpQMat[0]->Char() == 64 && m_pData->m_DataType == COD);
		FOR(i,m_iChar)	{ for(j=i+1;j<m_iChar;j++) { if(GenCodes[GenCode][i] != GenCodes[GenCode][j]) { Par->AddQij(i,j); } } }
	} else {
		assert(m_vpQMat[0]->Char() < 64 && m_pData->m_DataType == COD_RED);
		CurChar_i = 0;
		FOR(i,64)	{
			CurChar_j = CurChar_i+1;
			if(GenCodes[GenCode][i] == -1) { continue; }
			for(j=i+1;j<64;j++) {
				if(GenCodes[GenCode][j] == -1) { continue; }
				if(GenCodes[GenCode][i] != GenCodes[GenCode][j]) { Par->AddQij(CurChar_i,CurChar_j); }
				CurChar_j++;
			}
			CurChar_i++;
		}
	}
	m_vpPar.push_back(Par);
	return Par;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Function for adding omega in Dr/Dc type manner
// Radical Matrix consists of 0/1
CQPar * CCodonProcess::AddDrDcOmega(int GenCode, vector <int> RadMat, int Val2Add)	{
	int i,j, CurChar_i, CurChar_j;
	string Name = "Omega";
	CQPar *Par;
	// Check entry conditions
	assert(InRange(GenCode,0,11));
	FOR(i,20*20) { assert(InRange(RadMat[i],-1,2)); }
	assert(InRange(Val2Add,0,2));
	// Initialise
	if(Val2Add == 0) { Name += "_Conservative"; } else { Name += "_Radical"; }
	Par = new CQPar(Name,m_iChar,INITIAL_OMEGA + RandDouble(-0.05,0.05),true,MIN_OMEGA);

//	Par = new CQPar(Name,m_iChar,INITIAL_OMEGA,true,MIN_OMEGA);


	// Some DEBUG code
//	if(Val2Add == 0) { Par->SetVal(10); }
//	else { Val2Add = 25; }

	// Assign the parameters
	if(m_iChar == 64) { // When the Genetic code hasn't been applied
		assert(m_vpQMat[0]->Char() == 64 && m_pData->m_DataType == COD);
		FOR(i,m_iChar)	{ for(j=i+1;j<m_iChar;j++) { if(GenCodes[GenCode][i] != GenCodes[GenCode][j] && GenCodes[GenCode][j] != -1 && RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] == Val2Add) { Par->AddQij(i,j); } } }
	} else {
		assert(m_vpQMat[0]->Char() < 64 && m_pData->m_DataType == COD_RED);
		CurChar_i = 0;
		FOR(i,64)	{
			CurChar_j = CurChar_i+1;
			if(GenCodes[GenCode][i] == -1) { continue; }
			for(j=i+1;j<64;j++) {
				if(GenCodes[GenCode][j] == -1) { continue; }
				if(GenCodes[GenCode][i] != GenCodes[GenCode][j] && RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] == Val2Add) { Par->AddQij(CurChar_i,CurChar_j); }
				CurChar_j++;
			}
			CurChar_i++;
		}
	}
	m_vpPar.push_back(Par);

//	cout << "\nAdded Parameter: " << Par->Name();
//	cout << "\n" << Par->m_viQMap;

	return Par;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Set appropriate matrix values to zero
void CCodonProcess::AddMultiChangeZeros(int GenCode)	{
	CQPar *Par = new CQPar("CodonZero",m_iChar,0.0,false,0,0,REPLACE);
	int i,j,pos_i,pos_j,k,count;
	string a,b;
	pos_i = 0;

	FOR(i,64)		{
		if(GenCodes[GenCode][i] == -1 ) { if(m_pData->m_DataType != COD_RED) { pos_i++; } continue; }
		pos_j = pos_i + 1;
		for(j=i+1;j<64;j++)	{
			if(GenCodes[GenCode][j] == -1 ) { if(m_pData->m_DataType != COD_RED) { pos_j++; } continue; }
			a = m_sABET.substr(pos_i*3,3);
			b = m_sABET.substr(pos_j*3,3);
			assert(a.size() == b.size());
			count = 0; FOR(k,3) { if(a[k] != b[k]) { count++; } }
			if(count > 1) { Par->AddQij(pos_i,pos_j); pos_j++; continue; }
			if(m_pData->m_DataType == COD && (GenCodes[GenCode][i] == -1 || GenCodes[GenCode][j] == -1)) { Par->AddQij(pos_i,pos_j); pos_j++; continue; }
			pos_j++;
		}
		assert(pos_j == m_iChar);
		pos_i++;
	}
	assert(pos_i == m_iChar);
	m_vpPar.push_back(Par);
	Par = NULL;
}

// Old
/* void CCodonProcess::AddMultiChangeZeros(int GenCode)	{
	CQPar *Par = new CQPar("CodonZero",m_iChar,0.0,false,0,0,REPLACE);
	int i,j,pos_i,pos_j,k,count;
	string a,b;
	pos_i = 0;
	FOR(i,m_iChar)		{
		pos_j = 0;
		if(m_pData->m_DataType == COD_RED && GenCodes[GenCode][i] == -1 ) { continue; }
		for(j=i+1;j<m_iChar;j++)	{
			if(m_pData->m_DataType == COD_RED && GenCodes[GenCode][j] == -1 ) { continue; }
			a = m_sABET.substr(pos_i*3,3);
			b = m_sABET.substr(pos_j*3,3);
			assert(a.size() == b.size());
			count = 0; FOR(k,3) { if(a[k] != b[k]) { count++; } }
			if(count > 1) { Par->AddQij(i,j); continue; }
			if(m_pData->m_DataType == COD && (GenCodes[GenCode][i] == -1 || GenCodes[GenCode][j] == -1)) { Par->AddQij(i,j); continue; }
			pos_j++;
		}
		pos_i++;
	}
	m_vpPar.push_back(Par);
	Par = NULL;
}
*/
////////////////// Functions to calculate various rates ///////////////////////////////////////////////////////////
double CCodonProcess::NonsynRate(bool Observed) {
	int i,j,pos_i, pos_j, Correction =0;
	double Rate = 0.0;
	vector <double> eq;
	// Do some error checking
	if((int)m_vpQMat.size() > 1) { cout << "\nHaven't worked out how to do rates for multiple matrices in a single process"; }
	assert(InRange(m_iGenCode,0,12));
	// Get the rate
	CCodonProcess::PrepareQMats();
	if(!Observed) { eq.assign(m_iChar,1.0); } else { eq = RootEqm(); }
	pos_i = 0; FOR(i,64) {
		if(GenCodes[m_iGenCode][i] == -1) { Correction++; continue; }
		for(j=i+1,pos_j = i+1-Correction;j<64;j++) {
			if(GenCodes[m_iGenCode][i] == -1) { if(m_iChar != 64) { pos_i--; } break; }
			if(GenCodes[m_iGenCode][j] == -1) { if(m_iChar == 64) { pos_j++; } continue; }
			if(GenCodes[m_iGenCode][i] != GenCodes[m_iGenCode][j]) {
				Rate += eq[pos_i] * *m_vpQMat[0]->Q(pos_i,pos_j);
				Rate += eq[pos_j] * *m_vpQMat[0]->Q(pos_j,pos_i);
			}
			pos_j++;
		}
		pos_i++;
	}
	return Rate;
}

double CCodonProcess::SynRate(bool Observed) {
	int i,j,pos_i, pos_j, Correction = 0;
	double Rate = 0.0;
	vector <double> eq;
	// Do some error checking
	if((int)m_vpQMat.size() > 1) { cout << "\nHaven't worked out how to do rates for multiple matrices in a single process"; }
	assert(InRange(m_iGenCode,0,12));
	// Get the rate
	CCodonProcess::PrepareQMats();
	if(!Observed) { eq.assign(m_iChar,1.0); } else { eq = RootEqm(); }

	pos_i = 0; FOR(i,64) {
		if(GenCodes[m_iGenCode][i] == -1) { Correction++; continue; }
		for(j=i+1,pos_j = i+1-Correction;j<64;j++) {
			if(GenCodes[m_iGenCode][j] == -1) { if(m_iChar == 64) { pos_j++; } continue; }
//			cout << "\nComparing ["<<i<<":" << pos_i <<"][" << j << ":" << pos_j <<"] eqm = " << eq[pos_i] << " == " << *m_vpQMat[0]->Q(pos_i,pos_j) << ": " << GenCodes[m_iGenCode][i] << " cf. " << GenCodes[m_iGenCode][j];
			if(GenCodes[m_iGenCode][i] == GenCodes[m_iGenCode][j]) {
//				cout << " ... adding " << eq[i] * *m_vpQMat[0]->Q(pos_i,pos_j) << " and " << eq[j] * *m_vpQMat[0]->Q(pos_j,pos_i);
				Rate += eq[pos_i] * *m_vpQMat[0]->Q(pos_i,pos_j);
				Rate += eq[pos_j] * *m_vpQMat[0]->Q(pos_j,pos_i);
			}
			pos_j++;
		}
		pos_i++;
	}
//	cout << "\nReturning rate: " << Rate;
//	exit(-1);
	return Rate;
}

ostream &CCodonProcess::Output(ostream &os)	{
	int i;
	os<< "\n----- CodonProcess: " << m_sName << " : ID = " << m_iProcID << "; Rate = " << m_pRate->Val() << "; Prob = " << Prob() << " -----";
	if(m_bPseudoProcess == true) { os << "\nPseudoprocess"; }
	else {
		os << "\n" << m_vpPar.size() << " Parameters";
		FOR(i,(int)m_vpPar.size()) { os << "\n\t" << *m_vpPar[i]; }
		os << "\nObservedRateNonsynonymous:\t" << NonsynRate(true);
		os << "\nObservedRateSynonymous:\t" << SynRate(true);
//		os << "\n" << m_vpQMat.size() << " Q matrices";
//		FOR(i,(int)m_vpQMat.size()) { os << "\n\t" << *m_vpQMat[i]; os << "\nEqm: " << m_vpQMat[i]->Eqm();  }
/*		FOR(i,(int)m_vpQMat.size()) {
			os << "\n\t" << *m_vpQMat[i];
		}
		os << "\nEqm: "; FOR(i,64) { int j; FOR(j,64) { os << "\n" << m_vpQMat[0]->Eqm()[i];  } }
*/	}
/*
	string BET = "TCAG";
	cout << "\nSets: ";
	int j,k;
	FOR(i,4) { FOR(j,4) { FOR(k,4) {
				cout << endl << BET[i] << " " << BET[j] << " " << BET[k];
	}	}	}
*/

	return os;
}


