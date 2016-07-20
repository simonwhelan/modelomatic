//////////////// Implementation of some useful tools //////////////////
#include "tools.h"

#if DO_MEMORY_CHECK
CMemChecker memory_check;
#endif

// Implementation of troublesome math functions
int my_isnan(double x)  {
#ifdef __INTEL_COMPILER
  return isnan(x);
#else
  return std::isnan(x);
#endif
}

int my_isinf(double x) {
#ifdef __INTEL_COMPILER
  return isinf(x);
#else
  return std::isinf(x);
#endif
}

///////////////////////////////////////////////////////////////////////
// Some tools
bool ALLOW_PREDICTLNL = true;	// Whether to allow predicted lnLs in optimisation
#if FUNC_COUNTERS == 1
int LFunc_Log_Counter = 0;
int MakeQ_Log_Counter = 0;
int MakePT_Log_Counter = 0;
int Matrix_Log_Counter = 0;
int SubLFunc_Log_Counter = 0;
int SPR_Log_Counter = 0;

void OutputCounters(ostream &os)	{
	os << "\nCounter logs:\n\tLFunc:    " << LFunc_Log_Counter << "\n\tMakeQ:    "<< MakeQ_Log_Counter;
	os << "\n\tMakePT:   " << MakePT_Log_Counter << "\n\tMatrix:   " << Matrix_Log_Counter << "\n\tSubLFunc: " << SubLFunc_Log_Counter;
	os << "\n\tSPRBranches: " << SPR_Log_Counter;
}

#endif

#define PAR_DEBUG 1				// Whether debugging code for parameters is on (0: off; 1: on);
#define TOOLS_DEBUG 0			// Whether debugging code is on  (0: off;1: on);
#define HARD_DEBUG_PROBS 0		// Debugs probability functions by checking for nans. Only works under UNIX

//////////////////////////////////////////////////////////////////////////
// Data based routines
//////////////////////////////////////////////////////////////////////////

EDataType GetDataType(string s)	{
	if(s.find("RY2",0)!=string::npos) { return RY2; }
	if(s.find("DNA2",0)!=string::npos) { return DNA2; }
	if(s.find("AA2",0)!=string::npos) { return AA2; }
	if(s.find("RY",0)!=string::npos) { return RY; }
	if(s.find("DNA",0)!=string::npos) { return DNA; }
	if(s.find("AA",0)!=string::npos) { return AA; }
	if(s.find("COD",0)!=string::npos) { return COD; }
	if(s.find("RNASTEM",0)!=string::npos) { return RNASTEM; }
	if(s.find("RNA",0)!=string::npos) { return RNA; }
	return NONE;
}

string DataStates(EDataType Type)	{
	switch(Type)	{
	case RY:	return RY_ABET;
	case RY2:	return RY2_ABET;
	case DNA:	return DNA_ABET;
	case DNA2:	return DNA2_ABET;
	case AA:	return AA_ABET;
	case AA2: 	return AA2_ABET;
	case COD:	return COD_ABET;
	case GAP:	return GAP_ABET;
	case RNA:	return RNA_ABET;
	case RNASTEM: return RNASTEM_ABET;
	case COD_RED: Error("\nCan not get ABET for reduced codon set...");
	default: Error("Trying to return unspecified type of data\n"); break;
	}
	return "Unknown";
}
int LenStates(EDataType Type)	{
	switch(Type) {
	case RY:
	case DNA:
	case AA:
	case RNA:
	case GAP:	return 1;
	case RNASTEM: return 2;
	case RY2:
	case DNA2:
	case AA2: return 2;
	case COD_RED:
	case COD:	return 3;
	default: Error("Trying to get state size of unspecified data type\n"); break;
	}
	return -1;
}
int NumStates(EDataType Type)	{ return (int) DataStates(Type).size() / LenStates(Type);  }
string State(EDataType Type,int i)	{
	int pos; string Ret;
	if(i == NumStates(Type)) { FOR(pos,LenStates(Type)) { Ret += GAP_ABET[0]; } return Ret; }
	return GetPos(DataStates(Type),i,Type);
}
int FindState(EDataType Type,string ToFind,int Start)	{
	int i;
	transform(ToFind.begin(),ToFind.end(),ToFind.begin(),(int(*)(int)) toupper);
	for(i = Start; i< NumStates(Type); i++) { if(State(Type,i) == ToFind) { break; } }
	return i;
}
string GetPos(string Seq, int Pos, EDataType Type) {
	int length = LenStates(Type);
	return GetPos(Seq,Pos,length);
}
string GetPos(string Seq, int Pos, int PosLength)	{
	assert(Pos * PosLength < (int)Seq.size());
	return Seq.substr(Pos*PosLength,PosLength);
}
int CodonDiff(string c1, string c2)	{
	int i, diff = 0;
	assert(c1.size() == 3 && c2.size() == 3);
	FOR(i,3) { if(c1[i] != c2[i]) { diff++; } }
	return diff;
}
bool IsDataType(EDataType Type, string Seq, bool AllowGaps)	{
	int i;
	bool Check;
	transform(Seq.begin(),Seq.end(),Seq.begin(),(int(*)(int)) toupper);
	FOR(i,(int)Seq.size() / LenStates(Type)) {
		Check = true;
		if(FindState(DNA,GetPos(Seq,i,Type),0) == NumStates(Type)) {
			if(AllowGaps == true) {
				if(FindState(GAP,GetPos(Seq,i,Type),0) == NumStates(GAP)) { Check = false; break; } }
			else { Check = false; break; }
	}	}
	return Check;
}
bool IsDataType(EDataType Type, char Seq)	{
	int i;
	FOR(i,NumStates(Type)) { if(strcmp(State(Type,i).c_str(),&Seq) == 0) { return true; } }
	return false;
}
double RateDataType(EDataType Type, string Seq)	{
	int i;
	int Score = 0, count = 0;;
	transform(Seq.begin(),Seq.end(),Seq.begin(),(int(*)(int)) toupper);
	FOR(i,(int)Seq.size() / LenStates(Type)) {
		if(Seq[i] == 'N') { continue; } // N is a special type because it's used as non-specified for DNA
		if(FindState(Type,GetPos(Seq,i,Type),0) != NumStates(Type)) { Score += 1.0; } count++;
	}
	if(count == 0) { return -1.0; }
	return (double)Score / (double) count;
}
EDataType GuessDataType(string Data)	{
	double DNA_score = RateDataType(DNA,Data), AA_score = RateDataType(AA,Data);
	if(DNA_score < 0) { return NONE; }
	if(DNA_score * 1.1 > AA_score && DNA_score > MIN_DATA_PERCENT) { return DNA; }
	if(AA_score > MIN_DATA_PERCENT && AA_score  > DNA_score)	{ return AA; }
	return NONE;
}
bool IsTsByChar(char a, char b) {
	char c;
	if(a > b) { c = a; a = b; b = c; }
	if((a == 'A' && b == 'G') || (a == 'C' && b == 'T')) { return true; }
	return false;
}

bool IsGap(char c) {
	int i;
	FOR(i,(int)GAP_ABET.size()) { if(c == GAP_ABET[i]) { return true; } }
	return false;
}


bool IsTs(string a, string b, EDataType Type)	{
	int i;
	assert(Type == DNA || Type == COD || Type == COD_RED);
	if(Type == DNA) {	// The DNA transition function
		assert(a.size() == 1 && b.size() == 1);
		return IsTsByChar(a.c_str()[0],b.c_str()[0]);
	} else {			// The codon transition function
		if(CodonDiff(a,b) == 1)	{
			FOR(i,3) { if(a[i] != b[i]) { break; } }
			assert(i!=3);
			return IsTsByChar(a.c_str()[i],b.c_str()[i]);
	}	}
	return false;
}
ostream &operator<<(ostream &os, EDataType Type) {
	switch(Type)	{
	case RY:	os << "RY";		break;
	case RY2:	os << "RY2";	break;
	case AA:	os << "AA";		break;
	case AA2:	os << "AA2";	break;
	case DNA:	os << "DNA";	break;
	case DNA2:	os << "DNA2";	break;
	case COD:	os << "Codon";	break;
	case COD_RED:	os << "CodonReduced"; break;
	case RNA: 	os << "RNA"; 	break;
	case RNASTEM:	os << "RNASTEM"; break;
	case GAP:	os << "Gap";	break;
	case OTHER:	os << "Other";	break;
	case NONE:	os << "Unknown"; break;
	default: Error("Unknown EDataType passed to ostream\n");
	}
	return os;
}

// Get frequencies from sequence
// Note gaps are *NOT* counted
vector <double> GetFreqFromSeq(string Seq, string ABET) {
	int i,j;
	vector <double> fqs(ABET.size(),0);
	FOR(i,(int)Seq.size()) {
		FOR(j,(int)ABET.size()) {
			if(Seq[i] == ABET[j]) { fqs[j]++; break; }
	} }
	return NormaliseVector(fqs);
}
vector <double> GetFreqFromSeq(string Seq, EDataType DT) { return GetFreqFromSeq(Seq, DT); }
vector <double> GetFreqFromSeq(vector <int> Seq, int iChar) {
	int i;
	vector <double> fqs(iChar,0);
	FOR(i,Seq.size()) {
		if(Seq[i] < iChar) { fqs[Seq[i]] ++; }
	}
	return NormaliseVector(fqs);
}

/////////////////////////////////////////////////////////////////////
// Tree2Names
vector <string> Tree2Names(string Tree,int NoSeq)	{
	int i,Left = 0;
	vector <string> Names;
	string Temp;
	FOR(i,(int)Tree.size()) {
		// Get sequence
		if(Tree[i] == '(') { Left ++; }
		if((Tree[i] == '(' || Tree[i] == ',') && isalnum(Tree[i+1])) {
			Temp = "";
			for(i=i+1;i<(int)Tree.size();i++) {
				Temp += Tree[i];
				if(Tree[i+1] == ':' || Tree[i+1] == ',' || Tree[i+1] == ')') { Names.push_back(Temp); break; }
	}	}	}
	if(Left + 2 != NoSeq && NoSeq != 2) { Error("\nNumber of '(' in tree is wrong: obs= " + int_to_string(Left) + " cf. exp= " + int_to_string(NoSeq-2) + "\n"); }
	if(NoSeq != (int)Names.size()) { Error("\nNumber names in tree is wrong: obs= " + int_to_string((int)Names.size()) + " cf. exp= " + int_to_string(NoSeq) + "\n"); }
	return Names;
}

/////////////////////////////////////////////////////////////////////
// Error function
ostream& Error(string str, ostream  &os) { os << str; assert(0); exit(-1); };

// Vector manipulations
///////////////////////////
vector <double> NormaliseVector(vector <double> Vec, double VectorSumTo) {
	double sum = Sum(&Vec) / VectorSumTo;
	int i;
	if(sum > DBL_EPSILON) { FOR(i,(int)Vec.size()) { assert(Vec[i] >= 0); Vec[i] /= sum; } }
	else {					FOR(i,(int)Vec.size()) { assert(Vec[i] < DBL_EPSILON); Vec[i] = 1 / (double) Vec.size(); } }
	return Vec;
}

int DiscreteRand(vector <double> *Probs) {
	double count = 0, RandNum = Random();
	int i;
	FOR(i,(int)Probs->size()) { count += Probs->at(i); if(count > RandNum) { break; }}
	return i;
}

// Fitch and Margoliash
double GetRMSD(vector <double> Obs, vector <double> Exp) {
	int i;
	double RMSDVal = 0, O = 0 , E = 0;
	assert(Obs.size() == Exp.size());
	FOR(i,(int)Obs.size()) { if(Obs[i] < DBL_EPSILON) { continue; } RMSDVal += pow(Obs[i]-Exp[i],2) / (Obs[i]*Obs[i]); }
	return RMSDVal;
}

// Norm calculators
/////////////////////////////
double FrobeniusNorm(double * a,double * b, int n)	{
	int i;
	double total = 0.0;
	FOR(i,n) { total += pow(a[i] - b[i],2); }
	return sqrt(total);
}
double FrobeniusNorm(vector <double> &a,vector <double> &b)	{
	int i;
	double total = 0.0;
	assert((int)a.size() == (int)b.size());
	FOR(i,(int) a.size()) { total += pow(a[i] - b[i],2); }
	return sqrt(total);
}

// String manipulations
/////////////////////////////
string int_to_string(int num) { stringstream ss; ss << num; string str = ss.str(); return str; }
string double_to_string(double num) { stringstream ss; ss << num; string str = ss.str(); return str; }

void number(int a, char *ret) { stringstream ss; ss << a; string str = ss.str(); strcpy(ret,str.c_str()); }

string find_and_replace(string source, string find, string replace) {
	int pos=(int)source.find(find);
	while(pos != string::npos) { source.replace(pos,(int)find.size(),replace); pos=(int)source.find(find); }
	return source;
}

int wildcmp(const char *wild, const char *string) {
	// Written by Jack Handy - jakkhandy@hotmail.com
	// for comparisons of form "if (wildcmp("bl?h.*", "blah.jpg"))"
	//	where ? is single char and * is a run of char
	const char *cp = NULL, *mp = NULL;
	while ((*string) && (*wild != '*')) { if ((*wild != *string) && (*wild != '?')) { return 0; } wild++; string++; }
	while (*string) {
		if (*wild == '*') {
			if (!*++wild) { return 1; }
			mp = wild;
			cp = string+1;
		} else if ((*wild == *string) || (*wild == '?')) {
			wild++; string++;
		} else {
			wild = mp; string = cp++;
	}	}
	while (*wild == '*') { wild++; }
	return !*wild;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the parameter class
//////////////////////////////////////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &os, ParOp Par) {
	switch(Par)	{
	case REPLACE:	os << "Replace";	break;
	case MULTIPLY:	os << "Multiply";	break;
	case DIVIDE:	os << "Divide";		break;
	case ADD:		os << "Add";		break;
	case SUBTRACT:	os << "Subtract";	break;
	default: Error("Unknown parameter type in ostream");
	}
	return os;
}

// Constructor
CPar::CPar(string Name,double Value, bool Optimise, double LowBound, double UpBound,ParOp Type) {
#if DO_MEMORY_CHECK
	memory_check.CountCPar++;
#endif
#if PAR_DEBUG == 1
	if(my_isnan(Value)) { Error("\nTrying to assign parameter " + Name + " with nan...\n\n"); }
#endif
	m_dOldReal = m_dRealValue = Value;
	m_dOldScale = m_dScaledValue.RescaleVal(Value);
	m_sName = Name;
	// Initialisation stuff
	m_bOpt = Optimise; m_bLockedScale = false; m_bOutDetail = false;
	m_bSpecial = false; bool m_bAllowScale = true; m_bDoHardDer = false;
	m_dGrad = 0.0; m_bIsBranch = false;
	// Set default scaling routines
	pDoUpdate = NULL;
	m_Operator = Type;
	// Do the bounds (these can SetVal() so need to be at the end)
	SetBounds(LowBound,UpBound);
	StoreOptBounds(LowBound,UpBound);	// Default bounds are the optimisation bounds
}

// Destructor
CPar::~CPar() {
#if DO_MEMORY_CHECK
	memory_check.CountCPar--;
#endif
	int i;
	pDoUpdate = NULL;
	FOR(i,(int)m_arpPar.size()) { m_arpPar[i] = NULL; }
}

// Copy operator
CPar & CPar::operator=(CPar &Par)	{
#if PAR_DEBUG == 1
	if(my_isnan(Par.Val())) {
		cout << "\nHave nan in CPar::operator=(...): " << Par.m_sName << " = " << Par.Val();
		cout << "\nm_dRealValue: " << Par.m_dRealValue << "; m_dScaledValue: (Scaler=" << Par.m_dScaledValue.Scaler()<< "; Value=" << Par.m_dScaledValue.Value() << ")";
		Error("\nTrying to assign parameter " + Par.Name() + " with nan...\n\n");
	}
#endif
	assert(m_arpPar.empty());
	m_bOutDetail = Par.m_bOutDetail;
	m_sName = Par.m_sName;
	m_dOldReal = m_dRealValue = Par.m_dRealValue;
	m_dOldScale = Par.m_dOldScale;
	m_dScaledValue = Par.m_dScaledValue;
	m_ardBounds[0] = Par.m_ardBounds[0]; m_ardBounds[1] = Par.m_ardBounds[1];
	m_bOpt = Par.m_bOpt;
	m_bSpecial = Par.m_bSpecial;
	m_bAllowScale = Par.m_bAllowScale;
	m_bLockedScale = Par.m_bLockedScale;
	pDoUpdate = Par.pDoUpdate;
	m_Operator = Par.m_Operator;
	m_ardOptBounds[0] = Par.m_ardOptBounds[0];
	m_ardOptBounds[1] = Par.m_ardOptBounds[1];
	return *this;
}

double CPar::DoOper(double Value)	{
#if PAR_DEBUG == 1
	if(my_isnan(m_dRealValue) || my_isnan(m_dScaledValue.Value()) || my_isnan(Value)) {
		cout << "\nHave nan in CPar::DoOper=(...): " << Name() << " = " << Val();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		cout << "\nValue: " << Value;
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
	switch(m_Operator)	{
	case REPLACE:	return Val();
	case MULTIPLY:	return Value * Val();
	case DIVIDE:	return Value / Val();
	case ADD:		return Value + Val();
	case SUBTRACT:	return Value - Val();
	default:
		Error("Unknown type of operator...");
	}
	return 0;
}

// Parameter outputter
ostream &operator<<(ostream &os, CPar &Par)	{ return Par.Output(os); }
ostream &CPar::Output(ostream &os)	{
	int i;
	UpdatePar();
	string delim = " ";
	os << Name() << delim << "=" << delim << Val();
	if(m_bOutDetail == true) {
		os << delim << m_Operator;
		os << delim<< "("<< m_ardBounds[0] << ","<< m_ardBounds[1] << ")";
		if(Opt() == true) { os << delim << "Optimised"; } else { os << delim << "NotOptimised"; }
		os << delim << "ScaledValue="<<delim << Scal();
		if(!m_arpPar.empty()) { os << delim << "LinkedPar:"; FOR(i,(int)m_arpPar.size()) { os << delim << m_arpPar[i]->m_sName; } }
	}
	return os;
}

static int countery = 0;

double CPar::SetVal(double Value, bool Update, bool Force, bool Normalise) {
//	cout << "\nInto SetVal(" << Value << "," << Update << "," << Force << "," << Normalise << ")";
#if PAR_DEBUG == 1
	if(my_isnan(m_dRealValue) || my_isnan(Value)) {
		cout << "\nHave nan in CPar::SetVal(...): " << Name() << " = " << Val();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		cout << "\nValue: " << Value;
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
	m_dOldReal = m_dRealValue;
	m_dRealValue = Value;
	// If its a probability then the other parameters need adjusting
	// They are adjusted so that the set value is maintained.
	if(pDoUpdate == &ProbabilityScale && Normalise) {
//		cout << "\nUpdating " << Name();
		assert(m_arpPar.size() > 1);
		double total = 0.0;
		int i;
		assert(Value < 1.0 + DBL_EPSILON && Value > 0.0 - DBL_EPSILON);
		FOR(i,(int)m_arpPar.size()) { if(m_arpPar[i] == this) { continue; } total += m_arpPar[i]->Val(); }
		total /= (1.0 - Value);
		FOR(i,(int)m_arpPar.size()) { if(m_arpPar[i] == this) { continue; } m_arpPar[i]->m_dOldReal = m_arpPar[i]->m_dRealValue /= total; }
		total = 0.0; FOR(i,(int)m_arpPar.size()) { total += m_arpPar[i]->Val(); }
		assert(fabs(1.0-total) < FLT_EPSILON);
	}
//	if(Force == true && Update == true && pDoUpdate == &ProbabilityScale) { cout << "\n\nIn tools.cxx CPar::SetVal(...): doing a probability forced assignment with update... Is this really what you intend...\n\n"; }
	if(Force == true) { m_dOldReal = m_dRealValue; }
//	if(pDoUpdate == &ProbabilityScale) { cout << " 1: " << m_dRealValue; }
	if(Update == true) { UpdatePar(Force,Force); }
//	if(pDoUpdate == &ProbabilityScale) { cout << " 2: " << m_dRealValue; }
#if PAR_DEBUG == 1
	if(my_isnan(Val()) || my_isnan(Value)) {
		cout << "\nHave nan when exiting CPar::SetVal(...): " << Name() << " = " << Val();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		cout << "\nValue: " << Value;
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
//	if(pDoUpdate == &ProbabilityScale) { cout << " ... " << m_dRealValue; }
	return m_dRealValue;
}

// Bound checker (bounds defined elsewhere); if ForceBounds == true then will return parameter equal to or
//  within bounds.
bool CPar::CheckBound(bool ForceBounds) {
#if PAR_DEBUG == 1
	if(my_isnan(m_dRealValue)) {
		cout << "\nHave nan in CPar::CheckBound(...): " << Name();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
	static bool CheckOkay = true;
	if(CheckOkay == false) { return true; }
	CheckOkay = false;
	if(m_dRealValue < m_ardBounds[0] - 1.0E-6) {
		if(ForceBounds == true) { SetVal(m_ardBounds[0]); }
		CheckOkay = true;
		return false;
	} else if(m_dRealValue > m_ardBounds[1] + 1.0E-6) {
		if(ForceBounds == true) { SetVal(m_ardBounds[1]); }
		CheckOkay = true;
		return false;
	}
	CheckOkay = true;
	return true;
}

bool CPar::CheckLowBound(bool ForceBounds)	{
	static bool CheckOkay = true;
	if(CheckOkay == false) { return true; }
	CheckOkay = false;
//	if(m_dRealValue < m_ardBounds[0] - 1.0E-6) {
	if(fabs(m_dRealValue - m_ardBounds[0]) < 1.0E-6 || m_dRealValue < m_ardBounds[0] - 1.0E-6) {
//		cout << "\nLow Bound: " << m_dRealValue << " cf. " << m_ardBounds[0];
		if(ForceBounds == true) { SetVal(m_ardBounds[0]); }
		CheckOkay = true;
		return false;
	}
	CheckOkay = true;
	return true;
}
bool CPar::CheckUpBound(bool ForceBounds)	{
	static bool CheckOkay = true;
	if(CheckOkay == false) { return true; }
	CheckOkay = false;
//	if(m_dRealValue > m_ardBounds[1] + 1.0E-6) {
	if(fabs(m_dRealValue - m_ardBounds[1]) < 1.0E-6 || m_dRealValue > m_ardBounds[1] - 1.0E-6) {
		if(ForceBounds == true) { SetVal(m_ardBounds[1]); }
		CheckOkay = true;
		return false;
	}
	CheckOkay = true;
	return true;
}

double CPar::BoundDist()	{
	double dist;
	if(fabs(Val() - m_ardBounds[0]) < fabs(Val() - m_ardBounds[1])) { // Closests to lower bound
		dist = Val() - m_ardBounds[0];
	} else {
		dist = m_ardBounds[1] - Val();
	}
	return dist;
}

// Update parameters -- This will also include scaling when required
// arguments:
// bool ForceChange:	whether to allow probability values to be changed
// bool RedoScale:		whether the m_dScaledValue should be rescaled
bool CPar::UpdatePar(bool ForceChange, bool RedoScale) {
#if PAR_DEBUG == 1
	if(my_isnan(m_dRealValue) || my_isnan(m_dScaledValue.Value())) {
		cout << "\nHave nan when entering CPar::UpdatePar(...): " << Name();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
	bool OldLocked = m_bLockedScale;
	// When no change to the values has occurred then return
	if(pDoUpdate == NULL &&(fabs(m_dOldReal - m_dRealValue) < FLT_EPSILON && fabs(m_dOldScale - m_dScaledValue.Value()) < FLT_EPSILON && pDoUpdate == NULL && ForceChange == false)) {
		if(RedoScale == true) { m_dScaledValue.RescaleVal(m_dRealValue); }	// Reset the scaling if required
		CheckBound();
		return false;
	}
	// Check there is a scaling routine; if not then its a standard parameter
	// These are scaled to 1.0, so that: m_dValue = m_dScaler * m_dScaledValue
	if(pDoUpdate==NULL)	{
		assert(m_arpPar.empty());
		// Correct the scaling as appropriate
		if(m_dOldReal != m_dRealValue || ForceChange == true) {		// If the value has changed
			m_dScaledValue.NewVal(m_dRealValue);
		} else {								// If the scaledvalue has changed
			m_dRealValue = m_dScaledValue.Value();
		}
		if(RedoScale == true) { m_dScaledValue.RescaleVal(m_dRealValue); }	// Reset the scaling if required
	} else { // Otherwise there is a pDoUpdate routine that should be used
		int i;
		bool Changed = false;
		// Check whether any parameters have changed
		if(m_dOldReal != m_dRealValue && m_dOldScale != m_dScaledValue.Value()) { Changed = true; }
		if(Changed == false) {
			FOR(i,(int)m_arpPar.size()) { if(m_arpPar[i]->m_dOldScale != m_arpPar[i]->m_dScaledValue.Value()) { Changed = true; break; } }
		}
		if(ForceChange)	{ pDoUpdate(&m_arpPar,true,false,false);
		} else {
			if(Changed)	{ m_bLockedScale = true; pDoUpdate(&m_arpPar,false,false,false); }
			if(RedoScale) { pDoUpdate(&m_arpPar,false,false,true); }
		}
	}
	m_dOldReal = m_dRealValue; m_dOldScale = m_dScaledValue.Value();
	CheckBound();
#if PAR_DEBUG == 1
	if(my_isnan(m_dRealValue) || my_isnan(m_dScaledValue.Value())) {
		cout << "\nHave nan when exiting CPar::UpdatePar(...): " << Name();
		cout << "\nm_dRealValue: " << m_dRealValue << "; m_dScaledValue: (Scaler=" << m_dScaledValue.Scaler()<< "; Value=" << m_dScaledValue.Value() << ")";
		Error("\nTrying to assign parameter " + Name() + " with nan...\n\n");
	}
#endif
	return true;
}

//////////////////////////////////////////////////////////
// Function for setting parameter derivatives
double CPar::grad(double g)	{
//	cout << "\nSetting gradient: " << m_sName << " == " << Val() << " (" << m_ardBounds[0] << "," << m_ardBounds[1] << ")" << " grad = " << g;
//	if(!CheckUpBound()) { cout << " .. UP .."; }
//	if(!CheckLowBound()) { cout << " .. LOW .."; }
	if		(!CheckUpBound(true) && g < 0)	{ m_dGrad = 0.0; }
	else if (!CheckLowBound(true) && g > 0) { m_dGrad = 0.0; }
	else if (fabs(g) < GRAD_MIN)			{ m_dGrad = 0.0; }
	else { m_dGrad = g; }

//	cout << " --> " << m_dGrad;

	return m_dGrad;
}

double CScaledParVal::Value()		  {
	if(my_isnan(m_dValue)) { cout << "\nHave nan in value...\n\n"; exit(-1); }
	return m_dValue * m_dScaler; }						// Returns the real value
// Set new value (returned as Value) by adjusting m_dScaler and fixing m_dValue
double CScaledParVal::SoftNewVal(double Val) {
	assert(m_dScaler > DBL_EPSILON);
	if(m_dValue < FLT_EPSILON && diff(Val,m_dValue) == 1) { cout << "\nError: SoftNewVal"; }
	m_dScaler =  Val / m_dValue; return Value(); }
// Set new value (returned as Value) by adjusting m_dValue and fixing m_dScaler
double CScaledParVal::NewVal(double Val) { assert(m_dScaler > DBL_EPSILON); m_dValue = Val / m_dScaler; return Value(); }	// Changes so that Val = m_dValue * m_dScaler
double CScaledParVal::RescaleVal(double Val) { 											// Set value <- Val and  rescale so m_dValue = 1.
	if(fabs(Val) > PAR_TO_SCALE && m_bAllowScale == true) { m_dValue = 1.0; m_dScaler = Val; }
	else { m_dValue = Val; m_dScaler = 1.0; }
	if(my_isnan(Value()) || my_isnan(m_dValue)) { cout << "\nDone RescaleVal... have Value() == nan\n\n"; exit(-1); }
	return Val;
}


// The scaling routines
////////////////////////////////////////////////////
// These are assigned to the pDoUpdate routines
// If Val2Scale == true  -> Use value to obtain scaling factor
// If Val2Scale == false -> Use scaling factor to obtain real value
//
// Function transforming probabilities
// Unscaled[i] = exp(Scaled[i]) / ( 1.0 + Sum_j( exp(Scaled[j]) ) ) where j = 1,...,n;
void ProbabilityScale(vector <CPar*> *P, bool Value2Scale,bool first, bool ReorderProbs) {
	int i,j,non_opt = -1;
	vector <double> PVal;
	double total,Sum;
	bool ResetScales = false; 		// When a scale has run out of bounds, then correct it
	if(P->size() <2) { cout << "\nCannot make parameters of size 1 or less"; assert(0); exit(-1); }
	// First use of function is to creates a set of probability parameters from P
	if(first == true) {
		// Check the values are set up correctly
		FOR(i,(int)P->size()) { PVal.push_back(P->at(i)->m_dRealValue); }
		PVal = NormaliseVector(PVal,1.0); FixProb(&PVal);
		// Set up all the probabilities
		FOR(i,(int)P->size()) {
			P->at(i)->m_dOldReal = P->at(i)->m_dRealValue = PVal[i];
			assert(P->at(i)->m_arpPar.empty());
			FOR(j,(int)P->size()) { P->at(i)->m_arpPar.push_back(P->at(j)); }
		}
		Value2Scale = true;
	}
	// Reset the parameter that is used as the base for scaling
	if(first == true || ReorderProbs == true)	{
		// Take the largest probability as the one that will not be optimised
		total=0.0; j=0;
		FOR(i,(int)P->size()) {
			if(P->at(i)->m_dRealValue > total) { j=i; total = P->at(i)->m_dRealValue; }
			P->at(i)->SetOptimise(true); P->at(i)->SetSpecial(false);
			// Now set the correct bounds
			P->at(i)->SetUpdate(&ProbabilityScale);
				P->at(i)->SetBounds(0,1);
		}
		P->at(j)->SetOptimise(false); P->at(j)->SetSpecial(true);
	}
	FOR(i,(int)P->size()) { assert(!P->at(i)->m_arpPar.empty()); }
	// **************************************************************************************** //
	if(Value2Scale == true) {
	// Do the required scaling of Values, find the non-optimised probability, and set total
		PVal.clear();
//		cout << "\nSetting up Value2Scale";
		FOR(i,(int)P->size()) {
			PVal.push_back(P->at(i)->m_dRealValue);
			if(P->at(i)->Special() == true) {
//				cout << "\tHave found that " << i << " is special...";
				assert(non_opt == -1); non_opt = i;
				if(P->at(i)->m_dRealValue < FLT_EPSILON) { cout << "\nProblem 1" << flush; exit(-1); P->at(i)->m_dRealValue = FLT_EPSILON; }
				total = 1/ P->at(i)->m_dRealValue;
		}	}
		assert(non_opt != -1);
		FOR(i,(int)PVal.size()) {
			if(i!= non_opt) {
				if( (PVal[i] * total) < FLT_EPSILON)	{ P->at(i)->m_dScaledValue.NewVal(-100000); }
				else { P->at(i)->m_dScaledValue.NewVal(log(PVal[i] * total)); }
			} else { P->at(i)->m_dScaledValue.NewVal(0.0); }
			P->at(i)->m_dOldScale = P->at(i)->m_dScaledValue.Value();
		}
	}	else {
	// **************************************************************************************** //
	// Do the unscaling from Scale to Values
		PVal.clear();total = 1.0;
//		cout << "\nProbabilityScale->Values: ["; FOR(i,(int)P->size()) { cout << P->at(i)->m_dScaledValue.Value() << ","; } cout << "]";
//		cout << "\nBounds: "; FOR(i,(int)P->size()) { cout << "(" << P->at(i)->UpBound() << "," << P->at(i)->LowBound() << ")"; }
		FOR(i,(int)P->size()) {
			PVal.push_back(P->at(i)->m_dScaledValue.Value());
			if(P->at(i)->Special() == true) { assert(non_opt = -1); non_opt = i; }
		}
		FOR(i,(int)PVal.size()) { if(i==non_opt) { continue; } total += PVal[i] = exp(PVal[i]); }

		if(total < FLT_EPSILON) { cout << "\nWarning total reset..." << flush; exit(-1);total = FLT_EPSILON; }
		FOR(i,(int)PVal.size()) {
			if(i==non_opt)	{ PVal[i] = 1/total; }
			else			{ PVal[i] = PVal[i]/total; }
//			if(PVal[i] < MIN_PROB) { PVal[i] = MIN_PROB + FLT_EPSILON; ResetScales = true; }		// This is commented out because resetting probabilities during optimisation is really bad
		}
		PVal = NormaliseVector(PVal);
		Sum = 0;
		FOR(i,(int)PVal.size())	{
			P->at(i)->m_dRealValue = PVal[i];
			P->at(i)->m_dOldReal = P->at(i)->m_dRealValue;
			Sum += P->at(i)->m_dRealValue;
		}
//		cout << "\n\t -> ["; FOR(i,(int)P->size()) { cout << P->at(i)->m_dRealValue << ","; } cout << "]";
		if(fabs(Sum - 1.0) > 1.0e-6) {
			cout << "\nSum: " << Sum;
			Error("Odd frequencies... need to rescale\n"); }
	}
	// Assign the probabilities
	Sum = 0; FOR(i,(int)PVal.size()) { 	Sum += P->at(i)->m_dRealValue; }
	// Redo the scaling if required
	if(ResetScales) { ProbabilityScale(P,true,false,true); }

//	cout << "\nReturning probs: "; FOR(i,P->size()) { cout << P->at(i)->m_dRealValue << " "; }

}

int RandInt(int Begin, int End)	{
	return (int) (Begin + ( (End - Begin) * Random() ) );
}

/////////////////////////////////////////////////////////////
// Translates a sequence into match and delete stuff
string GetProfilePath(string Seq, string ABET)	{
	int i,j;
	string NewSeq;
	NewSeq.clear();
	FOR(i,(int)Seq.size()-1) {
		FOR(j,(int)ABET.size()) { if(Seq[i] == ABET[j]) { break; } }
		if(j != (int) ABET.size()) {
			FOR(j,(int)ABET.size()) { if(Seq[i+1] == ABET[j]) { break; } }
			if(j != (int) ABET.size()) { NewSeq = NewSeq + "MM "; }
			else { NewSeq = NewSeq + "MD "; }
		} else {
			FOR(j,(int)ABET.size()) { if(Seq[i+1] == ABET[j]) { break; } }
			if(j != (int) ABET.size()) { NewSeq = NewSeq + "DM "; }
			else { NewSeq = NewSeq + "DD "; }
		}
	}
	return NewSeq;
}
/////////////////////////////////////////////////////////////////////////
// Various bits and pieces

ostream &operator <<(ostream &os, EModel M)	{
	switch(M)	{
	case JC:	os << "Jukes and Cantor (JC)"; break;
	case FEL:	os << "Felsenstein (FEL)"; break;
	case K2P:	os << "Kimura 2-parameter (K2P)"; break;
	case HKY:	os << "Hasegawa, Kishino, and Yano (HKY)"; break;
	case REV:	os << "General time reversible (REV)"; break;
	case RY_model: os << "RY model (RY)"; break;
	case COV_HKY:	os << "Covarion model with HKY (COV_HKY)"; break;
	case COV_REV:	os << "Covarion model with REV (COV_REV)"; break;
	case HKYdG_THMM: os << "HKY+dG with rate categories linked by a THMM"; break;
	case THMM_FULLDNA:	os << "Full DNA temporal hidden Markov model (THMM_FULLDNA)"; break;
	case EQU:	os << "Equiprobable (EQU)"; break;
	case WAG:	os << "Whelan and Goldman (WAG)"; break;
	case JTT:	os << "Jones, Taylor, and Thornton (JTT)"; break;
	case DAY:	os << "Dayhoff et al. (DAY)"; break;
	case mtREV:	os << "mitochondrial REV (mtREV)"; break;
	case THMM_AA:	os << "AA temporal hidden Markov model (THMM_AA)"; break;
	case WAGdG_THMM: os << "WAG+dG with rate categories linked by a THMM (WAGdG_THMM)"; break;
//	case rtREV:	os << "Retroviral REV (rtREV)"; break;
	default:
		Error("Trying to output name of unknown model...",cerr);
	};
	return os;
}


////////////////////////////////////////////////////////////////////////////
// Some probability functions -- mostly from PAML

double Binomial(double n, int k, double *scale)
{
/* calculates n choose k. n is any real number, and k is integer.
   If(*scale!=0) the result should be c+exp(*scale).
*/
   double c=1,i,large=1e99;

   *scale=0;
   if((int)k!=k)
      Error("k is not a whole number in Binomial.");
   if(n<0 && k%2==1) c=-1;
   if(k==0) return(1);
   if(n>0 && (k<0 || k>n)) return (0);

   if(n>0 && (int)n==n) k=min(k,(int)n-k);
   for (i=1; i<=k; i++) {
      c*=(n-k+i)/i;
      if(c>large)
         { *scale+=log(c); c=1; }
   }
   return(c);
}

double probBinomial (int n, int k, double p)
{
/* calculates  {n\choose k} * p^k * (1-p)^(n-k)
*/
   double C, up, down;

   if (n<40 || (n<1000&&k<10)) {
      for (down=min(k,n-k),up=n,C=1; down>0; down--,up--) C*=up/down;
      if (fabs(p-.5)<1e-6) C *= pow(p,(double)n);
      else                 C *= pow(p,(double)k)*pow((1-p),(double)(n-k));
   }
   else  {
      C = exp((LnGamma(n+1.)-LnGamma(k+1.)-LnGamma(n-k+1.))/n);
      C = pow(p*C,(double)k) * pow((1-p)*C,(double)(n-k));
   }
   return C;
}


double probBetaBinomial (int n, int k, double p, double q)
{
/* This calculates beta-binomial probability of k succeses out of n trials,
   The binomial probability parameter has distribution beta(p, q)

   prob(x) = C1(-a,k) * C2(-b,n-k)/C3(-a-b,n)
*/
   double a=p,b=q, C1,C2,C3,scale1,scale2,scale3;

   if(a<=0 || b<=0) return(0);
   C1=Binomial(-a, k, &scale1);
   C2=Binomial(-b, n-k, &scale2);
   C3=Binomial(-a-b, n, &scale3);
   C1*=C2/C3;
   if(C1<0)
      Error("error in probBetaBinomial");
   return C1*exp(scale1+scale2-scale3);
}



////////////////////////////////////////////////////////////////////////////
// Functions for matrix multiplications
////////////////////////////////////////////////////////////////////////////

ostream &VecOut(int n, double *Vec, char delim,ostream &os)	{
	int i;
	FOR(i,n) { os << Vec[i] << delim; }
	return os;
}

ostream &MatOut(int n, double *Mat,char delim, ostream &os)	{
	int i;
	FOR(i,n*n)	{ if(i%n == 0) { os << endl; } os << Mat[i] << delim; }
	return os;
}
ostream &MatOut(int n, vector <double> Mat,char delim, ostream &os)	{
	int i;
	FOR(i,n*n)	{ if(i%n == 0) { os << endl; } os << Mat[i] << delim; }
	return os;
}
ostream &MatOut(int n,double **Mat,char delim, ostream &os)	{
	int i,j;
	FOR(i,n)	{
		os << endl;
		FOR(j,n)	{
			os << Mat[i][j] << delim;
	}	}
	return os;
}
void VMat(double *a,double *b, double *c, int n)	{
// a[1*SIZE], b[SIZE*SIZE], c[1*SIZE]  ......  c = a*b
    static int i,j;
    static double *p_a,*p_b,*p_c;
	p_a = a; p_b = b; p_c = c;
	static double t;
    for(i=0;i<n;i++, p_c++)	{
        p_a = a;
        for(j=0, t=0;j<n;j++) { t += (*(p_a++))*(*(p_b++)); }
        *p_c = t;
    }
	p_a = NULL; p_b = NULL; p_c = NULL;
#if FUNC_COUNTERS == 1
	Matrix_Log_Counter++;
#endif
}

void MulMat(double *a, double *b, double *c, int n, int m, int k)	{
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b */
   int i,j,i1;
   double t;
   FOR (i,n)  FOR(j,k) {
      for (i1=0,t=0.0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
      c[i*k+j] = t;
	}
#if FUNC_COUNTERS == 1
	Matrix_Log_Counter++;
#endif

}

void IMat(double *I, int n)		{ int i; ZeroMat(I,n); FOR(i,n)	{ I[(i*n)+i] = 1.0; } }
void IMat(double **I, int n)		{ int i,j; FOR(i,n) { FOR(j,n) { I[i][j] = 0.0; } I[i][i] = 1.0; } }
void ZeroMat(double *M, int n)	{ int i; FOR(i,n*n) { M[i] = 0.0; } }

//////////////////////////////////////////////////////////////////////////////////////////
// Scaled probability class
//////////////////////////////////////////////////////////////////////////////////////////
// Simple class containing scaled probabilities - used for likelihood computations
// Contains probabilities in the form RealValue = pow(m_dValue,(-10*m_iScale));

// Useful probability function
bool IsProb(double V)	{ if(InRange((double) V,(double) -FLT_EPSILON,(double)1.0+FLT_EPSILON)) { return true; }
#if HARD_DEBUG_PROBS
cout.precision(12);
	cout << "\nValue not a probability: " << V;
#endif
	return false;
}
bool FlipBool(bool V)	{ if(V == true) { return false; } return true; }
bool FlipBin(int i)		{ assert(i==0||i==1); if(i==0) { return 1; } return 0; }
// Constructor
CProb::CProb(double InitVal)	{
#if DO_MEMORY_CHECK
	memory_check.CountCProb++;
#endif
	m_dValue = 0.0; assert(IsProb(InitVal)); Assign(InitVal); }
CProb::CProb(CProb &Prob)		{
#if DO_MEMORY_CHECK
	memory_check.CountCProb++;
#endif
	m_dValue = 0.0; Assign(Prob); }
CProb::CProb(double Val, int Sc){
#if DO_MEMORY_CHECK
	memory_check.CountCProb++;
#endif
	m_dValue = 0.0; Assign(Val,Sc); }

CProb::~CProb()	{
#if DO_MEMORY_CHECK
	memory_check.CountCProb--;
#endif

}

///////////////////// Private functions //////////////////////////////////////////////////////

void CProb::DoScale()	{
	// Allow true zero
	if(Double_Zero(m_dValue) || m_iScale > 1000000000)  {
		m_dValue = 0.0; m_iScale = 0; return;
	}
	// Otherwise adjust so that double m_dValue hold only a single digit
	while(m_dValue < 1.0)	{ m_dValue *= 10; m_iScale++; }
	while(m_dValue > 10.0)	{ m_dValue /= 10; m_iScale--; }
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue)) { cout << "\nReturning CProb::DoScale(): m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif

}

///////////////////// Public functions ///////////////////////////////////////////////////////

// Ostream operator
ostream &operator<<(ostream &os,CProb &Prob)	{ os << Prob.m_dValue << "*10^-" << Prob.m_iScale; return os; }

// Return functions
double CProb::Prob()	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue)) { cout << "\nReturning Prob(): m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	if(m_iScale > 400) { return 0.0; }
	return m_dValue * pow(10.0,(double) -m_iScale);
}
double CProb::LogP()	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue)) { cout << "\nReturning LogP(): m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	if(IsZero()) { if(MATCH_PAML == 0) { return -BIG_NUMBER; } else { return log(pow((double)10,-80)); }} return (log(m_dValue) + (-LOG10 * m_iScale));
}
bool CProb::IsZero()	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue)) { cout << "\nReturning IsZero(): m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	if(!Double_Zero(m_dValue)) { return false; } return true;
}
// Copy operator
CProb &CProb::operator=(CProb &Prob)	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Prob.m_dValue)) { cout << "\nReturning CProb::operator=(): Prob.m_dValue= " << Prob.m_dValue << "; Prob.m_iScale= " << Prob.m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	m_dValue = Prob.m_dValue;
	m_iScale = Prob.m_iScale;
	return *this;
}

// Assign operator
CProb &CProb::Assign(double Value)		{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Value)) { cout << "\nReturning Assign(double Value): Value= " << Value << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	if(!(Value <= 1.0 + FLT_EPSILON && Value >= 0.0 - FLT_EPSILON)) {
		cout << "\nAbout to fail assert statement: Value = " << Value;
	}
	assert(Value <= 1.0 + FLT_EPSILON && Value >= 0.0 - FLT_EPSILON); m_dValue = Value; m_iScale = 0; DoScale(); return *this;
}
CProb &CProb::Assign(CProb &Prob)		{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Prob.m_dValue)) { cout << "\nReturning CProb::Assign(CProb &): Prob.m_dValue= " << Prob.m_dValue << "; Prob.m_iScale= " << Prob.m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	m_dValue = Prob.m_dValue; m_iScale = Prob.m_iScale; DoScale(); return *this;
}
CProb &CProb::Assign(double Val, int Sc){
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Val)) { cout << "\nReturning Assign(double Val, int Sc): Val= " << Val << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	m_dValue = Val; m_iScale = Sc; if(!IsProb(Prob())) {
		cout << "\nBroken Prob(): " << Prob() << " assigning from: " << Val << " with scale " << Sc << endl << flush;
		exit(-1);
	}
	assert(IsProb(Prob())); DoScale(); return *this;
}
// double functions for numerical operations
CProb &CProb::Multiply(double Value,bool Overwrite)	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Value)) { cout << "\nReturning CProb::Multiply(double Value, bool): Value= " << Value << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	CProb Prob(Value); return Multiply(Prob,Overwrite);
}
CProb &CProb::Add(double Value, bool Overwrite)		{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Value)) { cout << "\nReturning CProb::Add(double Value, bool): Value= " << Value << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	CProb Prob(Value); return Add(Prob,Overwrite);
}
CProb &CProb::Divide(double Value, bool Overwrite)	{
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Value)) { cout << "\nReturning CProb::Divide(double Value, bool): Value= " << Value << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	CProb Prob(Value); return Divide(Prob,Overwrite);
}

// CProb functions for numerical operations
CProb &CProb::Multiply(CProb &Prob, bool Overwrite) {
	static CProb New;
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Prob.m_dValue)) { cout << "\nReturning CProb::Multiply(CProb &,bool): Prob.m_dValue= " << Prob.m_dValue << "; Prob.m_iScale= " << Prob.m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	if(IsZero() || Prob.IsZero()) { if(Overwrite == true) { Zero(); } New.Zero(); }
	else {
		New = *this;
//		cout << "\nCProb::Multiply -- New: " << &New << " cf *this: " << this << " cf Prob: " << &Prob << " multiplication: " << New.m_dValue << " * " << Prob.m_dValue;
		// Do the multiplication
		New.m_dValue *= Prob.m_dValue; New.m_iScale += Prob.m_iScale; New.DoScale();
		// Replace if required
		if(Overwrite == true) { *this = New; }
	}
	return New;
}
CProb &CProb::Add(CProb &Prob, bool Overwrite) {
	static CProb New;
	CProb added(Prob);
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Prob.m_dValue)) { cout << "\nReturning CProb::Add(CProb &,bool): Prob.m_dValue= " << Prob.m_dValue << "; Prob.m_iScale= " << Prob.m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	// Return Prob if zero
	if(IsZero()) { New = Prob; if(Overwrite) { *this = New; } }
	else if(Prob.IsZero()) { New = *this; }
	else {
		New = *this;
		// Match the scales, do the addition and rescale
		New.MatchScales(&added,true); New.m_dValue += added.m_dValue; New.DoScale();
		// Replace if required
#if TOOLS_DEBUG == 1
		if(Double_Zero(New.m_dValue)) { cout .precision(16); cout << "\n Adding: " << Prob << " + " << *this << " ... zero: " << New << "__"; }
#endif
		if(Overwrite) { *this = New; }
	}
	return New;
}
CProb &CProb::Divide(CProb &Prob, bool Overwrite)	{
	static CProb New;
	New = *this;
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue) || my_isnan(Prob.m_dValue)) { cout << "\nReturning CProb::Divide(CProb &,bool): Prob.m_dValue= " << Prob.m_dValue << "; Prob.m_iScale= " << Prob.m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
	// Do the division
	Error("\nHaven't done division...");
	// Replace if required
	if(Overwrite == true) { *this = New; }
	return New;

}

// Matches the scales between two probabilities; this is dangerous as information can be lost)
// The argument DoHigh decides whether the highest or lowest scale factor takes precedence
void CProb::MatchScales(CProb *Prob, bool DoHigh)	{
	bool DoThis = false;
	int Diff = abs(Prob->m_iScale - m_iScale);
	// Find the probability with the highest scaling factor
	if(Diff == 0) { return; }
	if(m_iScale > Prob->m_iScale)	{ DoThis = true; }
	if(DoHigh == false)				{ DoThis = FlipBool(DoThis); }
	// If (DoThis == true) then this Probability is rescaled, otherwise the other probability is rescaled
	if(DoThis == true)	{ m_iScale -= Diff; m_dValue *= pow(10.0,(double) -Diff); }
	else				{ Prob->m_iScale -= Diff; Prob->m_dValue *= pow(10.0,(double) -Diff); }
#if HARD_DEBUG_PROBS == 1
	if(my_isnan(m_dValue)) { cout << "\nReturning CProb::MatchScales(CProb* Prob): Prob->m_dValue= " << Prob->m_dValue << "; Prob->m_iScale= " << Prob->m_iScale << "; m_dValue= " << m_dValue << "; m_iScale= " << m_iScale; exit(-1); }
#endif
}

// Ensures that all values in a probability vector have no true zeros
void FixProb(vector <double> *Prob)	{
	bool Fixed = false;
	int i;
	FOR(i,(int)Prob->size()) { if(Prob->at(i) < MIN_PROB) { Prob->at(i) = MIN_PROB; Fixed = true; } }
	if(Fixed) { *Prob = NormaliseVector(*Prob,1.0); }
}

// Reads data from file
string GetDataLine(ifstream *in)	{
	string line = "*";
	string::size_type loc;
	do {
		line = "#";
		while(line[0] == '\r' || line[0] == '#' || line[0] =='\n' || line[0] == '\0' || line.find_first_not_of(' ') == std::string::npos) {
			getline(*in,line);
			if(in->eof() == 1) { cout << "\nUnexpected end of file in GetDataLine(...). Try adding a new line at the end of the data file\n\n"; cout << "\nLine was: " << line; exit(-1); }
		}
		loc = line.find("*");
	} while (loc != string::npos);
	return line;
}

vector <string> Tokenise(string line) {
	string buf;
	stringstream in(line);
	vector <string> Toks;
	Toks.~vector();
	while(in >> buf) { Toks.push_back(buf); }
	return Toks;
}
vector <string> Tokenise(string line, string Delim)	{
	size_t i = 0, j,j1;
	vector <string> Toks;
	while(i != (int)line.size())	{
		j = line.find(Delim,i+1);
		if(j == string::npos) { j = j1 = (int)line.size(); } else { j1 = j+1; }
		Toks.push_back(line.substr(i,j-i));
		i = j1;
	}
	return Toks;
}

string EatWhiteSpace(string line)	{
	string::iterator it, temp;
	IFOR(it,line)	{ if(isspace(*it)) { temp = it; it = it - 1; line.erase(temp); } }
	return line;
}

/* *********** Eigen routines *************************************

	All code taken directly from PAML package (Yang 1997)

	Some small alterations have been made

  ************************************************************ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

				The simple eigen routine that assumes a symmetric matrix

   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */


int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0,
               double Root[], double U[], double V[])	{
/*
   This finds the eigen solution of the rate matrix Q for a time-reversible
   Markov process, using the algorithm for a real symmetric matrix.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
   where S is symmetrical, all elements of pi are positive, and U*V = I.
   pi_sqrt[n-npi0] has to be calculated before calling this routine.

   [U 0] [Q_0 0] [U^-1 0]    [Root  0]
   [0 I] [0   0] [0    I]  = [0     0]

   Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
*/
   int i,j, inew,jnew, nnew=n-npi0, status;
   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */
   if(pi_sqrt==NULL)  { Error("pi_sqrt should be calculated before");; }
   FOR(i,n) { if(pi_sqrt[i] < FLT_EPSILON) { pi_sqrt[i] = FLT_EPSILON; } }
   if(npi0==0) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);
      status=eigenRealSym(U, n, Root, V);
      for(i=0;i<n;i++) for(j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0;i<n;i++) for(j=0;j<n;j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i] > DBL_EPSILON) {
            for(j=0,jnew=0; j<i; j++)
               if(pi[j] > DBL_EPSILON) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew]
                                    = Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }
      status=eigenRealSym(U, nnew, Root, V);
      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i] > DBL_EPSILON) ? Root[inew--] : 0;
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i] > DBL_EPSILON) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j] > DBL_EPSILON) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[j];
                  jnew--;
               }
               else
                  V[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i] > DBL_EPSILON) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j] > DBL_EPSILON) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[i];
                  jnew--;
               }
               else
                  U[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }
   Root[0]=0;
   return(status);
}


/* eigen solution for real symmetric matrix */

void HouseholderRealSym(double a[], int n, double d[], double e[]);
int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
void EigenSort(double d[], double U[], int n);

int eigenRealSym(double A[], int n, double Root[], double Offdiag[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return,
   A has the right vectors and Root has the eigenvalues. work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
   and then using the QL algorithm with implicit shifts.

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, Offdiag);
   status=EigenTridagQLImplicit(Root, Offdiag, n, A);
   EigenSort(Root, A, n);

   return(status);
}


void EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0;i<n-1;i++) {
      p=d[k=i];
      for (j=i+1;j<n;j++)
         if (d[j] >= p) p=d[k=j];
      if (k != i) {
         d[k]=d[i];
         d[i]=p;
         for (j=0;j<n;j++) {
            p=U[j*n+i];
            U[j*n+i]=U[j*n+k];
            U[j*n+k]=p;
         }
      }
   }
}

void HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0) ? -sqrt(h) : sqrt(h);
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      }
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the
   tridiagonal matrix, or the output from HouseholderRealSym() to get the
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
               else if(bb==0)             r=0;
               else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}

#undef SIGN

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Various functions taken directly from PAML for gamma distributed rates
//////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////// discrete gamma-distribution related //////////////////////////////////////
int DiscreteGamma(double freqK[], double rK[],double alfa, double beta, int K, int median)
{
/* discretization of gamma distribution with equal proportions in each
   category
*/
   int i;
   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
         freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
         freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}

double LBinormal(double h1, double h2, double r)
{
/* L(h1,h2,r) = prob(x>h1, y>h2), where x and y are standard binormal,
   with r=corr(x,y),  error < 2e-7.
      Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
      bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107.
*/
   int i;
   double x[]={0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
   double w[]={0.018854042, 0.038088059, 0.0452707394,0.038088059,0.018854042};
   double Lh=0, r1, r2, r3, rr, aa, ab, h3, h5, h6, h7, h12;


   h12=(h1*h1+h2*h2)/2;
   if (fabs(r)>=0.7) {
      r2=1-r*r;   r3=sqrt(r2);
      if (r<0) h2*=-1;
      h3=h1*h2;   h7=exp(-h3/2);
      if (fabs(r)!=1) {
         h6=fabs(h1-h2);   h5=h6*h6/2; h6/=r3; aa=.5-h3/8;  ab=3-2*aa*h5;
         Lh = .13298076*h6*ab*(1-CDFNormal(h6))
            - exp(-h5/r2)*(ab+aa*r2)*0.053051647;
         for (i=0; i<5; i++) {
            r1=r3*x[i];  rr=r1*r1;   r2=sqrt(1-rr);
            Lh-=w[i]*exp(-h5/rr)*(exp(-h3/(1+r2))/r2/h7-1-aa*rr);
         }
      }
      if (r>0) Lh = Lh*r3*h7+(1-CDFNormal(max(h1,h2)));
      else if (r<0) Lh = (h1<h2?CDFNormal(h2)-CDFNormal(h1):0) - Lh*r3*h7;
   }
   else {
      h3=h1*h2;
      if (r!=0)
         for (i=0; i<5; i++) {
            r1=r*x[i]; r2=1-r1*r1;
           Lh+=w[i]*exp((r1*h3-h12)/r2)/sqrt(r2);
         }
      Lh=(1-CDFNormal(h1))*(1-CDFNormal(h2))+r*Lh;
   }
   return (Lh);
}

double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}


double LnGamma(double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
               +.083333333333333)/x;
}

double PointNormal(double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}

double CDFNormal(double x)
{
/* Hill ID  (1973)  The normal integral.  Applied Statistics, 22:424-427.
   Algorithm AS 66.   (error < ?)
   adapted by Z. Yang, March 1994.  Hill's routine does not look good, and I
   haven't consulted
      Adams AG  (1969)  Algorithm 39.  Areas under the normal curve.
      Computer J. 12: 197-198.
*/
    int invers=0;
    double p, limit=10, t=1.28, y=x*x/2;

    if (x<0) {  invers=1;  x*=-1; }
    if (x>limit)  return (invers?0:1);
    if (x<t)
       p = .5 - x * (    .398942280444 - .399903438504 * y
                   /(y + 5.75885480458 - 29.8213557808
                   /(y + 2.62433121679 + 48.6959930692
                   /(y + 5.92885724438))));
    else
       p = 0.398942280385 * exp(-y) /
           (x - 3.8052e-8 + 1.00000615302 /
           (x + 3.98064794e-4 + 1.98615381364 /
           (x - 0.151679116635 + 5.29330324926 /
           (x + 4.8385912808 - 15.1508972451 /
           (x + 0.742380924027 + 30.789933034 /
           (x + 3.99019417011))))));
    return (invers ? p : 1-p);
}

double PointChi2(double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      cout << "\nerr incomplete gamma";
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;


   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}

////////////////////////////////////////////////////////////////////
// Some other probability distribution stuff

// Sample from exponential distribution with mean 1/lambda
double SampleExp(double lambda)	{  return (-log(1-Random())/lambda); }

// Sample from zipf (as used in DAWG)
// (Probably won't be used and just set a discrete distribution early)
int Zipf(double a, int MaxLength)
{
	double b = pow(2.0, a-1.0);
	double x=MaxLength+1,t;
	do {
		do {
			x = floor(pow(Random(), -1.0/(a-1.0)));
			t = pow(1.0+1.0/x, a-1.0);
		} while( Random()*x*(t-1.0)*b > t*(b-1.0));
	} while (!InRange((int)x,0,MaxLength+1));
	return (int)x;
}
// Gets the probability of a zipf law distribution truncated at Maxlength
vector <double> ProbZipf(double a, int MaxLength) {
	int i;
	vector <double> dist(MaxLength,0);
	FOR(i,MaxLength) { dist[i] = pow(i+1,-a); }
	return NormaliseVector(dist);

}
// What's the mean length from a zipf with parameter a, conditional on the max length being MaxLength
double MeanZipf(double a, int MaxLength) {
	int i;
	double mean = 0.0;
	vector <double> lengths = ProbZipf(a,MaxLength);

//	cout << "\nLengths: " << lengths;

	FOR(i,MaxLength) { mean += lengths[i] * (double)(i+1); }

//	cout << "\nMean length: " << mean;
	return mean;

}

// Riemann-Zeta
double GetRiemannZeta(double a) {
	int k=2;
	double total=1,diff;
	while(k) {
		diff = pow((double)k,-a);
		total += diff;
		cout << "\nk: " << k << "\tcumulative: " << total << "\tchange" << diff;
		k++;
	}

	return total;
}

////////////////////////////////////////////////////////////////////
// Some menu and file based stuff

bool FileExist(string File)	{ ifstream in(File.c_str()); if(in) { return true; } in.close(); return false; }

bool GetYesNoOption(string message)	{
	int Count = 0; string opt = "x";
	cout << message << " [y/n]?  ";
	while(opt != "y" && opt != "n")	{
		Count++; if(Count>100) { Error("Option broken...\n"); }
		getline(cin,opt);
	}
	if(opt == "y" || opt == "yes") { return 1; }
	return 0;
}

int GetOption(string message, int low, int high)	{
	int Opt = low - 1, Count = 0;
	string Input;
	while(Opt < low || Opt > high)	{
		if(Count++ > 10) { cout << "Too indecisive... Exiting\n\n"; exit(-1); }
		cout << message << " ["<<low<<" - "<<high<<"]:  ";
		getline(cin,Input);
		Opt = atoi(Input.c_str());
	}
	return Opt;
}

string GetInFileName()	{
	int Count = 0;
	string retval = "\0", file;
	while(retval == "\0")	{
		if(Count++ > 10) { cout << "Too indecisive... Exiting\n\n"; exit(-1); }
		cout << "\nPlease give me an input file name [ls for directory listing; q for exit ]: ";
		getline(cin,file);
		if(file.empty()) { continue; }
		if(FileExist(file)) { return file; }
		else if("q" == file) { exit(-1); }
		else if("ls" == file) {
			cout << "\nDirectory contains:\n";
			if(system("ls") == 1) { system("cls"); system("dir | more"); };
		} else { cout << "\nThe file <" <<file<<"> doesn't exist"; }
	}
	return retval;
}

string GetOutFileName(string SuggestedFile)	{
	int Count = 0;
	string retval = "\0", file;
	while(retval == "\0")	{
		if(Count++ > 10) { cout << "Too indecisive... Exiting\n\n"; exit(-1); }
		while(file.empty())	{
			if(Count++ > 10) { cout << "Too indecisive... Exiting\n\n"; exit(-1); }
			read(file,"Please give me an output file name ",SuggestedFile);
			if(file.empty()) { file = SuggestedFile; }
		}
		if(FileExist(file)) {
			cout << "\nWarning: file <"<<file<<"> exists already...";
			if(GetYesNoOption("Do you want to delete it [y/n]: ")) {
				ofstream out("file"); out.close(); retval = file;
			} else { file = "\0"; }
		} else { retval = file; }
	}
	return retval;
}

/////////////////////////////////////////////////////////////////////
// Some other basic Linear Algebra stuff
bool IsProbMat(int n, double *M)	{
	int i; vector <double> Mat;
	FOR(i,n*n) { Mat.push_back(M[i]); }
	return IsProbMat(Mat);
}
bool IsProbMat(vector <double> M)	{
	int i,j,n = (int) sqrt((double)M.size());
	double Total;
	if(n*n != M.size()) { Error("Invalid matrix passed to IsProbMat"); }
	FOR(i,n)	{
		Total = 0.0;
		FOR(j,n) { Total += M[(i*n)+j]; }
		if(diff(1.0,Total)) { return false; }
	}
	return true;
}

/////////////////////////////////////////////////////////////////////
// Function that checks detailed balance equations for a matrix
bool CheckReversibility(vector <double> Mat, vector <double> Eqm, bool Output)	{
	bool RetVal = true;
	int i,j,n = (int)Eqm.size();
	assert(Mat.size() == n * n);
	FOR(i,n)	{
		for(j=i+1;j<n;j++)	{
//			if(i<4 && j < 4) { cout << "\nPi["<<i<<"] * Q["<<i<<","<<j<<"] == Pi["<<j<<"] * Q["<<j<<","<<i<<"]: " << Eqm[i] * Mat[(i*n)+j] <<" == " << Eqm[j] * Mat[(j*n)+i] << "; diff = " << fabs(Eqm[i] * Mat[(i*n)+j] - Eqm[j] * Mat[(j*n)+i]);  }
			if(tdiff(Eqm[i] * Mat[(i*n)+j],Eqm[j] * Mat[(j*n)+i],1.0E-6)) {
				cout << "\nPi["<<i<<"] * Q["<<i<<","<<j<<"] != Pi["<<j<<"] * Q["<<j<<","<<i<<"]: " << Eqm[i] * Mat[(i*n)+j] <<" == " << Eqm[j] * Mat[(j*n)+i] << "; diff = " << fabs(Eqm[i] * Mat[(i*n)+j] - Eqm[j] * Mat[(j*n)+i]); RetVal = false;
			}
	}	}
	return RetVal;
}

bool CheckReversibility(int n, double *Mat, double *Eqm,bool Output) {
	int i; vector <double> M,E; FOR(i,n*n) { M.push_back(Mat[i]); } FOR(i,n) { E.push_back(Eqm[i]); }
	return CheckReversibility(M,E,Output);
}
bool CheckReversibility(int n, double *Mat, vector <double> Eqm,bool Output) {
	int i; vector <double> M; FOR(i,n*n) { M.push_back(Mat[i]); }
	return CheckReversibility(M,Eqm,Output);
}

/* -------------------------------------------------------------------------
 * This is an ANSI C library for multi-stream random number generation.
 * The use of this library is recommended as a replacement for the ANSI C
 * rand() and srand() functions, particularly in simulation applications
 * where the statistical 'goodness' of the random number generator is
 * important.  The library supplies 256 streams of random numbers; use
 * SelectStream(s) to switch between streams indexed s = 0,1,...,255.
 *
 * The streams must be initialized.  The recommended way to do this is by
 * using the function PlantSeeds(x) with the value of x used to initialize
 * the default stream and all other streams initialized automatically with
 * values dependent on the value of x.  The following convention is used
 * to initialize the default stream:
 *    if x > 0 then x is the state
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively.
 *
 * The generator used in this library is a so-called 'Lehmer random number
 * generator' which returns a pseudo-random number uniformly distributed
 * 0.0 and 1.0.  The period is (m - 1) where m = 2,147,483,647 and the
 * smallest and largest possible values are (1 / m) and 1 - (1 / m)
 * respectively.  For more details see:
 *
 *       "Random Number Generators: Good Ones Are Hard To Find"
 *                   Steve Park and Keith Miller
 *              Communications of the ACM, October 1988
 *
 * Name            : rngs.c  (Random Number Generation - Multiple Streams)
 * Authors         : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 09-22-98
 * -------------------------------------------------------------------------
 */

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */

static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */


   double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed
 * between 0.0 and 1.0.
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0)
    seed[stream] = t;
  else
    seed[stream] = t + MODULUS;

  return ((double) seed[stream] / MODULUS);
}


   void PlantSeeds(long x)
/* ---------------------------------------------------------------------
 * Use this function to set the state of all the random number generator
 * streams by "planting" a sequence of states (seeds), one per stream,
 * with all states dictated by the state of the default stream.
 * The sequence of planted states is separated one from the next by
 * 8,367,782 calls to Random().
 * ---------------------------------------------------------------------
 */
{
  const long Q = MODULUS / A256;
  const long R = MODULUS % A256;
        int  j;
        int  s;

  initialized = 1;
  s = stream;                            /* remember the current stream */
  SelectStream(0);                       /* change to stream 0          */
  PutSeed(x);                            /* set seed[0]                 */
  stream = s;                            /* reset the current stream    */
  for (j = 1; j < STREAMS; j++) {
    x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
    if (x > 0)
      seed[j] = x;
    else
      seed[j] = x + MODULUS;
   }
}


   void PutSeed(long x)
/* ---------------------------------------------------------------
 * Use this function to set the state of the current random number
 * generator stream according to the following conventions:
 *    if x > 0 then x is the state (unless too large)
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively
 * ---------------------------------------------------------------
 */
{
  char ok = 0;

  if (x > 0)
    x = x % MODULUS;                       /* correct if x is too large  */
  if (x < 0)
    x = ((unsigned long) time((time_t *) NULL)) % MODULUS;
  if (x == 0)
    while (!ok) {
      printf("\nEnter a positive integer seed (9 digits or less) >> ");
      scanf("%ld", &x);
      ok = (0 < x) && (x < MODULUS);
      if (!ok)
        printf("\nInput out of range ... try again\n");
    }
  seed[stream] = x;
}


   void GetSeed(long *x)
/* ---------------------------------------------------------------
 * Use this function to get the state of the current random number
 * generator stream.
 * ---------------------------------------------------------------
 */
{
  *x = seed[stream];
}


   void SelectStream(int index)
/* ------------------------------------------------------------------
 * Use this function to set the current random number generator
 * stream -- that stream from which the next random number will come.
 * ------------------------------------------------------------------
 */
{
  stream = ((unsigned int) index) % STREAMS;
  if ((initialized == 0) && (stream != 0))   /* protect against        */
    PlantSeeds(DEFAULT);                     /* un-initialized streams */
}


   void TestRandom(void)
/* ------------------------------------------------------------------
 * Use this (optional) function to test for a correct implementation.
 * ------------------------------------------------------------------
 */
{
  long   i;
  long   x;
  double u;
  char   ok = 0;

  SelectStream(0);                  /* select the default stream */
  PutSeed(1);                       /* and set the state to 1    */
  for(i = 0; i < 10000; i++)
    u = Random();
  GetSeed(&x);                      /* get the new state value   */
  ok = (x == CHECK);                /* and check for correctness */

  SelectStream(1);                  /* select stream 1                 */
  PlantSeeds(1);                    /* set the state of all streams    */
  GetSeed(&x);                      /* get the state of stream 1       */
  ok = ok && (x == A256);           /* x should be the jump multiplier */
  if (ok)
    printf("\n The implementation of rngs.c is correct.\n\n");
  else
    printf("\n\a ERROR -- the implementation of rngs.c is not correct.\n\n");
}

