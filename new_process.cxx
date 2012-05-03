/* ***********************************************************
 *  Implementation of new process stuff
 * *********************************************************** */

#include "process.h"
#include "new_process.h"

// Note: DEBUG functions are located in process.h

#if FUNC_COUNTERS == 1
	extern int Matrix_Log_Counter, MakeQ_Log_Counter, MakePT_Log_Counter, LFunc_Log_Counter, SubLFunc_Log_Counter;
#endif

//////////////////////////////////////////////////////
// Functions that return values of amino acid models
vector <double> vWAGVal()	{
	int i;
	vector <double> Ret;
	double dWAGVal[190] = {0.551571,0.509848,0.635346,0.738998,0.147304,5.429420,1.027040,0.528191,0.265256,0.0302949,0.908598,3.035500,1.543640,0.616783,0.0988179,1.582850,0.439157,0.947198,6.174160,0.021352,5.469470,1.416720,0.584665,1.125560,0.865584,0.306674,0.330052,0.567717,0.316954,2.137150,3.956290,0.930676,0.248972,4.294110,0.570025,0.249410,0.193335,0.186979,0.554236,0.039437,0.170135,0.113917,0.127395,0.0304501,0.138190,0.397915,0.497671,0.131528,0.0848047,0.384287,0.869489,0.154263,0.0613037,0.499462,3.170970,0.906265,5.351420,3.012010,0.479855,0.0740339,3.894900,2.584430,0.373558,0.890432,0.323832,0.257555,0.893496,0.683162,0.198221,0.103754,0.390482,1.545260,0.315124,0.174100,0.404141,4.257460,4.854020,0.934276,0.210494,0.102711,0.0961621,0.0467304,0.398020,0.0999208,0.0811339,0.049931,0.679371,1.059470,2.115170,0.088836,1.190630,1.438550,0.679489,0.195081,0.423984,0.109404,0.933372,0.682355,0.243570,0.696198,0.0999288,0.415844,0.556896,0.171329,0.161444,3.370790,1.224190,3.974230,1.071760,1.407660,1.028870,0.704939,1.341820,0.740169,0.319440,0.344739,0.967130,0.493905,0.545931,1.613280,2.121110,0.554413,2.030060,0.374866,0.512984,0.857928,0.822765,0.225833,0.473307,1.458160,0.326622,1.386980,1.516120,0.171903,0.795384,4.378020,0.113133,1.163920,0.0719167,0.129767,0.717070,0.215737,0.156557,0.336983,0.262569,0.212483,0.665309,0.137505,0.515706,1.529640,0.139405,0.523742,0.110864,0.240735,0.381533,1.086000,0.325711,0.543833,0.227710,0.196303,0.103604,3.873440,0.420170,0.398618,0.133264,0.428437,6.454280,0.216046,0.786993,0.291148,2.485390,2.006010,0.251849,0.196246,0.152335,1.002140,0.301281,0.588731,0.187247,0.118358,7.821300,1.800340,0.305434,2.058450,0.649892,0.314887,0.232739,1.388230,0.365369,0.314730};
	FOR(i,190) { Ret.push_back(dWAGVal[i]); }
	return Ret;
}
vector <double> vWAGFreq() {
	int i;
	vector <double> Ret;
	double dWAGFreq[20] = {0.0866279,0.043972,0.0390894,0.0570451,0.0193078,0.0367281,0.0580589,0.0832518,0.0244313,0.048466,0.086209,0.0620286,0.0195027,0.0384319,0.0457631,0.0695179,0.0610127,0.0143859,0.0352742,0.0708956};
	FOR(i,20) { Ret.push_back(dWAGFreq[i]); }
	return Ret;
}
vector <double> vJTTVal()	{
	int i;
	vector <double> Ret;
	double dJTTVal[190] = {58,54,45,81,16,528,56,113,34,10,57,310,86,49,9,105,29,58,767,5,323,179,137,81,130,59,26,119,27,328,391,112,69,597,26,23,36,22,47,11,17,9,12,6,16,30,38,12,7,23,72,9,6,56,229,35,646,263,26,7,292,181,27,45,21,14,54,44,30,15,31,43,18,14,33,479,388,65,15,5,10,4,78,4,5,5,40,89,248,4,43,194,74,15,15,14,164,18,24,115,10,102,21,16,17,378,101,503,59,223,53,30,201,73,40,59,47,29,92,285,475,64,232,38,42,51,32,33,46,245,25,103,226,12,118,477,9,126,8,4,115,18,10,55,8,9,52,10,24,53,6,35,12,11,20,70,46,209,24,7,8,573,32,24,8,18,536,10,63,21,71,298,17,16,31,62,20,45,47,11,961,180,14,323,62,23,38,112,25,16};
	FOR(i,190) { Ret.push_back(dJTTVal[i]); }
	return Ret;
}
vector <double> vJTTFreq() {
	int i;
	vector <double> Ret;
	double dJTTFreq[20] = {0.0767479,0.0516909,0.042645,0.0515439,0.019803,0.040752,0.0618299,0.0731519,0.022944,0.0537609,0.0919039,0.0586759,0.023826,0.040126,0.0509009,0.0687649,0.0585649,0.014261,0.032102,0.0660049};
	FOR(i,20) { Ret.push_back(dJTTFreq[i]); }
	return Ret;
}
vector <double> vDAYVal()	{
	int i;
	vector <double> Ret;
	double dDAYVal[190] = {27,98,32,120,0,905,36,23,0,0,89,246,103,134,0,198,1,148,1153,0,716,240,9,139,125,11,28,81,23,240,535,86,28,606,43,10,65,64,77,24,44,18,61,0,7,41,15,34,0,0,73,11,7,44,257,26,464,318,71,0,153,83,27,26,46,18,72,90,1,0,0,114,30,17,0,336,527,243,18,14,14,0,0,0,0,15,48,196,157,0,92,250,103,42,13,19,153,51,34,94,12,32,33,17,11,409,154,495,95,161,56,79,234,35,24,17,96,62,46,245,371,26,229,66,16,53,34,30,22,192,33,136,104,13,78,550,0,201,23,0,0,0,0,0,27,0,46,0,0,76,0,75,0,24,8,95,0,96,0,22,0,127,37,28,13,0,698,0,34,42,61,208,24,15,18,49,35,37,54,44,889,175,10,258,12,48,30,157,0,28};
	FOR(i,190) { Ret.push_back(dDAYVal[i]); }
	return Ret;
}
vector <double> vDAYFreq() {
	int i;
	vector <double> Ret;
	double dDAYFreq[20] = {0.087127,0.040904,0.040432,0.046872,0.033474,0.038255,0.049530,0.088612,0.033618,0.036886,0.085357,0.080482,0.014753,0.039772,0.050680,0.069577,0.058542,0.010494,0.029916,0.064718};
	FOR(i,20) { Ret.push_back(dDAYFreq[i]); }
	return Ret;
}
vector <double> vMTREVVal()	{
	int i;
	vector <double> Ret;
	double dMTREVVal[190] = {23.18,26.95,13.24,17.67,1.90,794.38,59.93,103.33,58.94,1.90,1.90,220.99,173.56,55.28,75.24,9.77,1.90,63.05,583.55,1.90,313.56,120.71,23.03,53.30,56.77,30.71,6.75,28.28,13.90,165.23,496.13,113.99,141.49,582.40,49.12,1.90,96.49,1.90,27.10,4.34,62.73,8.34,3.31,5.98,12.26,25.46,15.58,15.16,1.90,25.65,39.70,1.90,2.41,11.49,329.09,8.36,141.40,608.70,2.31,1.90,465.58,313.86,22.73,127.67,19.57,14.88,141.88,1.90,65.41,1.90,6.18,47.37,1.90,1.90,11.97,517.98,537.53,91.37,6.37,4.69,15.20,4.98,70.80,19.11,2.67,1.90,48.16,84.67,216.06,6.44,90.82,54.31,23.64,73.31,13.43,31.26,137.29,12.83,1.90,60.97,20.63,40.10,50.10,18.84,17.31,387.86,6.04,494.39,69.02,277.05,54.11,54.71,125.93,77.46,47.70,73.61,105.79,111.16,64.29,169.90,480.72,2.08,238.46,28.01,179.97,94.93,14.82,11.17,44.78,368.43,126.40,136.33,528.17,33.85,128.22,597.21,1.90,21.95,10.68,19.86,33.60,1.90,1.90,10.92,7.08,1.90,32.44,24.00,21.71,7.84,4.21,38.58,9.99,6.48,1.90,191.36,21.21,254.77,38.82,13.12,3.21,670.14,25.01,44.15,51.17,39.96,465.58,16.21,64.92,38.73,26.25,195.06,7.64,1.90,1.90,1.90,19.00,21.14,2.53,1.90,1222.94,91.67,1.90,387.54,6.35,8.23,1.90,204.54,5.37,1.90};
	FOR(i,190) { Ret.push_back(dMTREVVal[i]); }
	return Ret;
}
vector <double> vMTREVFreq() {
	int i;
	vector <double> Ret;
	double dMTREVFreq[20] = {0.072,0.019,0.039,0.019,0.006,0.025,0.024,0.056,0.028,0.088,0.169,0.023,0.054,0.061,0.054,0.072,0.086,0.029,0.033,0.043};
	FOR(i,20) { Ret.push_back(dMTREVFreq[i]); }
	return Ret;
}

vector <double> vcpREVVal() {
	int i;
	vector <double> Ret;
	double dcpREVVal[190] = {105,227,357,175,43,4435,669,823,538,10,157,1745,768,400,10,499,152,1055,3691,10,3122,665,243,653,431,303,133,379,66,715,1405,331,441,1269,162,19,145,136,168,10,280,92,148,40,29,197,203,113,10,396,286,82,20,66,1745,236,4482,2430,412,48,3313,2629,263,305,345,218,185,125,61,47,159,202,113,21,10,1772,1351,193,68,53,97,22,726,10,145,25,127,454,1268,72,327,490,87,173,170,285,323,185,28,152,117,219,302,100,43,2440,385,2085,590,2331,396,568,691,303,216,516,868,93,487,1202,1340,314,1393,266,576,241,369,92,32,1040,156,918,645,148,260,2151,14,230,40,18,435,53,63,82,69,42,159,10,86,468,49,73,29,56,323,754,281,1466,391,142,10,1971,89,189,247,215,2370,97,522,71,346,968,92,83,75,592,54,200,91,25,4797,865,249,475,317,122,167,760,10,119};
	FOR(i,190) { Ret.push_back(dcpREVVal[i]); }
	return Ret;
}
vector <double> vcpREVFreq() {
	int i;
	vector <double> Ret;
	double dcpREVFreq[20] = {0.0755,0.0621,0.0410,0.0371,0.0091,0.0382,0.0495,0.0838,0.0246,0.0806,0.1011,0.0504,0.0220,0.0506,0.0431,0.0622,0.0543,0.0181,0.0307,0.0660};
	FOR(i,20) { Ret.push_back(dcpREVFreq[i]); }
	return Ret;
}

/* ********************** Class defining s parameters in covarion style models ************* */
// Constructor currently the same
CCSCovPar::CCSCovPar(string Name, int DataChar, int HiddenChar, double Value, bool Opt, double L, double U, ParOp Oper) : CQPar(Name,DataChar, Value, Opt, L, U, Oper)	{ m_iHiddenChar = HiddenChar; };

void CCSCovPar::AddCSstates(int From, int To,string sABET, int iABET)	{
	int i,j, DataChar = m_iChar / m_iHiddenChar;
	if(iABET != 1) { cout << "\nHaven't done multiple ABET yet..."; }
	if(From == -1)	{ assert(m_viCSstatesFrom.empty()); }
	if(To == -1)	{ assert(m_viCSstatesTo.empty()); }
	m_viCSstatesFrom.push_back(From);
	m_viCSstatesTo.push_back(To);
	FOR(i,m_iChar) {
		for(j=i+1;j<m_iChar;j++)	{
			if(i / DataChar == j / DataChar) { continue; }	// Ignore the on-diagonal elements of the matrix
			if(From != -1 && i / DataChar != From) { continue; }
			if(To != -1 && j / DataChar != To) { continue; }
			if(sABET[i] != sABET[j]) { continue; }
			AddQij(i,j,true);
	}	}}

ostream &CCSCovPar::Output(ostream &os)	{
	int i;
	CQPar::Output(os);
	os << "\n\t\tCS states (total=" << m_iHiddenChar << ") linked by parameters: ";
	os << " (From: ";
	if(m_viCSstatesFrom.size() == 1 && m_viCSstatesFrom[0] == -1) { os << "all"; }
	else { FOR(i,(int)m_viCSstatesFrom.size()) { os << m_viCSstatesFrom[i] << " "; } }
	os << "; To: ";
	if(m_viCSstatesTo.size() == 1 && m_viCSstatesTo[0] == -1) { os << "all"; }
	else { FOR(i,(int)m_viCSstatesTo.size()) { os << m_viCSstatesTo[i] << " "; } }
	os << ")";
	return os;
}

/////////////////////////////////////////////////////////////////////////////
//			CTHMMQMat -- class for describing Q in covarion models
/////////////////////////////////////////////////////////////////////////////
// Constructor
CTHMMQMat::CTHMMQMat(int DataChar,int Char, EDataType Type, string Name) : CQMat(Char,Type,Name)	{
	m_iDataChar = DataChar;
}
// Overall rate calculator -- if(ForceWholeProcess == false) -> ignores transitions between classes
double CTHMMQMat::OverallSubRate(bool ForceWholeProcess)	{
	if(ForceWholeProcess) { return CQMat::OverallSubRate(); }
	int i, j,k,StateNo = Char() / m_iDataChar;
	double Total = 0.0;
//	cout << "\nCalculating rate of Covarion process: "; OutQ();
//	cout << "\nEqm: " << Eqm();
	// Do some error checking with debug code
	assert(m_bScaleReady == true);
	assert(Char() % m_iDataChar == 0 && StateNo >= 0);
	Total = 0;
	FOR(i,StateNo)	{
		for(j=(i * m_iDataChar);j<(i+1)*m_iDataChar;j++)	{
			for(k=j+1;k<(i+1)*m_iDataChar;k++)	{
//				cout << "\n\tAdding ("<<j<<","<<k<<"): " << *Q(j,k) <<  " * " << m_ardEqm[j]<< " == " << *Q(j,k) * m_ardEqm[k];
//				cout << "\n\tAdding ("<<k<<","<<j<<"): " << *Q(k,j) <<  " * " << m_ardEqm[k]<< " == " << *Q(k,j) * m_ardEqm[j];
				Total += *Q(j,k) * m_ardEqm[j];
				Total += *Q(k,j) * m_ardEqm[k];
	}	}	}
//	cout << "\n\nTotal rate: "<< Total << "\n\n";
	return Total;
}
// Calculate the rate of transitions
double CTHMMQMat::OverallTransRate()	{
	int i,j;
	double Total = 0.0;
	FOR(i,Char())	{
		FOR(j,Char())	{
			if(i / m_iDataChar == j / m_iDataChar) { continue; }
			Total += *Q(i,j) * m_ardEqm[i];
	}	}
	return Total;
}
// Calculate the transition rate from a process
double CTHMMQMat::TransRateFrom(int From)	{
	int i,j;
	double Total = 0.0;
	// Check entry conditions
//	cout << "\nCalculating rate from HiddenState " << From;
	assert(InRange(From,0,Char() / m_iDataChar) && InRange(From,0,Char() / m_iDataChar));
	for(i=(From * m_iDataChar);i<(From+1) * m_iDataChar;i++)	{
		FOR(j,Char()) {
			if(From == j / m_iDataChar) { continue; }
//			cout << "\n\tAdding rate from Q(" << i << "," << j << "): " << *Q(i,j) << " * " <<  m_ardEqm[i] << " == " << *Q(i,j) * m_ardEqm[i];
			Total += *Q(i,j) * m_ardEqm[i];
	}	}
//	cout << "\nReturning rate: " << Total;
	return Total;
}
// Calculate the transition rate to a process
double CTHMMQMat::TransRateTo(int To)		{
	int i,j;
	double Total = 0.0;
	// Check entry conditions
//	cout << "\nCalculating rate from HiddenState " << To;
	assert(InRange(To,0,Char() / m_iDataChar) && InRange(To,0,Char() / m_iDataChar));
	FOR(i,Char())		{
		for(j=(To * m_iDataChar);j<(To + 1) * m_iDataChar;j++)	{
			if(To == i / m_iDataChar) { continue; }
//			cout << "\n\tAdding rate from Q(" << i << "," << j << "): " << *Q(i,j) << " * " <<  m_ardEqm[i] << " == " << *Q(i,j) * m_ardEqm[i];
			Total += *Q(i,j) * m_ardEqm[i];
	}	}
//	cout << "\nReturning rate " << Total;
	return Total;
}
// Scale the Q matrix and the eigen roots to give the process the correct rate.
void CTHMMQMat::ScaleQ(double Rate)	{
	m_bAlwaysI = false;
	// Deal with zero rates
	assert(Rate >= 0.0);
//	cout << "\nTHMMQMat::ScaleQ(" << Rate << ")" << flush;
//	cout << "\nQout: "; CTHMMQMat::OutQ
//	cout << "\n\n---\nScaling Q (Rate=" << OverallSubRate(false) << ") -> " << Rate<<": ";
//	cout << "\nm_dScale = " << m_dScale << "\nEqm:\t" << Eqm();
	if(Rate < 1.0E-6 || m_dScale < DBL_EPSILON) {
		if(OverallSubRate(false) > DBL_EPSILON && m_dScale < DBL_EPSILON) { Error("\nUnexpected values in CQMat::ScaleQ.\n"); }
		m_bAlwaysI = true;  // This is true even for covarion processes -- cannot have zero substitutions with changes of state!
		return;
	}
	assert(m_bScaleReady == true);
	int i;

	double NewRate = Rate / m_dScale;
	FOR(i,m_iChar)	{ m_ardRoot[i] *= NewRate; }
	FOR(i,m_iChar2)	{ m_ardQMat[i] *= NewRate; }
//	cout << "\nRescaled Q: "; OutQ();
//	cout << "\nNew rate: " << OverallSubRate(false) << " cf. " << OverallSubRate(true) << "\n~~~";
//	exit(-1);
	m_dScale = Rate;
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////// CCovEqm class //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Simple CCovEqm contructor
CCovEqm::CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <double> eqm, vector <CPar *> *StateProbs)	: CBaseEqm(Char,ProcPar) {
	m_iNoCovState = (int)StateProbs->size();
	ConstructSimpleCovEqm(eqm,StateProbs);
}
CCovEqm::CCovEqm(int Char, vector <CQPar *> *ProcPar, CData *Data, vector <CPar *> *StateProbs)	: CBaseEqm(Char,ProcPar)			{
	m_iNoCovState = (int)StateProbs->size();
	assert(Char == Data->m_iChar * StateProbs->size() && Data->m_vFreq.size() * StateProbs->size() == Char);
	ConstructSimpleCovEqm(Data->m_vFreq,StateProbs);
}
// ConstructSimpleCovEqm
void CCovEqm::ConstructSimpleCovEqm(vector <double> eqm,vector <CPar *> *StateProbs)	{
	int i;
	vector <int> iV;
	EDataType Type;
	CQPar *Par;
	string Name;
	DoBasicNonReversibleCovarion = false;
	// Check entry conditions
	assert(eqm.size() * m_iNoCovState == m_iChar);
	assert(StateProbs->size() == m_iNoCovState);
	assert(fabs(Sum(&eqm) - 1.0) < FLT_EPSILON);
	switch(eqm.size())	{
	case 4:		Type = DNA; m_iDataChar = 4; break;
	case 20:	Type = AA;	m_iDataChar = 20; break;
	case 64:	Type = COD;	m_iDataChar = 64; break;
	default: Error("Unknown type of eqm data");
	};
	// Set dumb things for complex eqm distributions correctly
	m_iNoCovObsEqm = 1;
	FOR(i,m_iNoCovState) { iV.push_back(i); } m_vviEqm2State.push_back(iV);
	iV.clear();
	assert(m_vviEqm2State.size() == 1);
	// Add the eqm parameters
	FOR(i,(int)eqm.size())	{
		Name = "Freq(" + State(Type,i) + ")";
		Par = new CQPar(Name,m_iChar,eqm[i],true,MIN_PROB,1.0,MULTIPLY);
		m_vpEqmPar.push_back(Par);
		m_pProcPar->push_back(Par);
		Par = NULL;
	}
	ProbabilityScale(&m_vpEqmPar,true,true,true);
	// Add the state probabilities
	FOR(i,m_iNoCovState) { m_vpStateProbs.push_back(StateProbs->at(i)); }
	FOR(i,(int)m_vpStateProbs.size()) { m_vpEqmPar.push_back(m_vpStateProbs[i]); }
	ProbabilityScale(&m_vpStateProbs,true,true,true);
}
//////////////////////////////////////////////////////////////////////
// Constructor function for complex covarion equilibrium
CCovEqm::CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <vector <double> > Eqm, vector <vector <int> > EqmMap, vector <CPar *> *StateProbs, bool AllowOpt) : CBaseEqm(Char,ProcPar) {
	int EqmNo, i ,j;
	double Total;
	EDataType Type;
	CQPar *Par;
	string Name;
	vector <bool> ProcCheck;
	DoBasicNonReversibleCovarion = false;
	vector <CPar *> EqmPar;
	// Check entry conditions
	assert(EqmMap.size() == Eqm.size());
	switch(Eqm[0].size())	{
	case 4:		Type = DNA; m_iDataChar = 4; break;
	case 20:	Type = AA;	m_iDataChar = 20; break;
	case 64:	Type = COD;	m_iDataChar = 64; break;
	default: Error("Unknown type of eqm data");
	};
	assert(Char % m_iDataChar == 0);

	// Set up the hidden process probabilities
	m_iNoCovState = Char / m_iDataChar;
	m_iNoCovObsEqm = (int)Eqm.size();
	assert(m_iNoCovState == StateProbs->size());
	FOR(i,m_iNoCovState) { m_vpStateProbs.push_back(StateProbs->at(i)); }
	ProbabilityScale(&m_vpStateProbs,true,true,true);
	ProcCheck.assign((int)m_iNoCovState,(bool)false);
	FOR(i,(int)EqmMap.size()) { FOR(j,(int)EqmMap[i].size())	{ assert(ProcCheck[EqmMap[i][j]] == false); ProcCheck[EqmMap[i][j]] = true; }	}
	FOR(i,(int)m_iNoCovState) { assert(ProcCheck[i] == true); }
	// Set up the complex eqm distribution
	FOR(EqmNo,(int)Eqm.size())	{
		assert(m_iDataChar == Eqm[EqmNo].size());
		// Check eqm sums to 1.0
		Total = 0.0; FOR(i,m_iDataChar) { Total += Eqm[EqmNo][i]; }
		if(diff(Total,1.0)) { cout << "\nCCovEqm::CCovEqm(...): Sum(Eqm["<<EqmNo<<"]) != 1.0\n\t" << Eqm[EqmNo] << " == " << Total; exit(-1); }
		// Build The Eqm parameters
		assert(EqmPar.empty());
		FOR(i,m_iDataChar) {
			// Do name of type "Freq(<Char>[<process_numbers>])"
			Name = "Freq(" + State(Type,i) + "[" + int_to_string(EqmMap[EqmNo][0]);
			FOR(j,(int)EqmMap[EqmNo].size()-1) { Name += "," + int_to_string(EqmMap[EqmNo][j+1]); } Name += "])";
//			cout << "\nName is: " << Name << " == " << Eqm[EqmNo][i];
			Par = new CQPar(Name,m_iChar,Eqm[EqmNo][i],AllowOpt,MIN_PROB,1.0,MULTIPLY);
			EqmPar.push_back(Par); m_vpEqmPar.push_back(Par); m_pProcPar->push_back(Par);
			Par = NULL;
		}
		// Create the probabilities
		ProbabilityScale(&EqmPar,true,true,true);
		// Clean up
		FOR(i,(int)EqmPar.size()) { EqmPar[i] = NULL; } EqmPar.clear();
	}
	// Add the state probabilitys to m_vpEqmPar
	FOR(i,(int)m_vpStateProbs.size()) { m_vpEqmPar.push_back(m_vpStateProbs[i]); }
	// Transfer the mapping function
	m_vviEqm2State = EqmMap;
//	cout << "\nFinished trying to build CCovEqm::CCovEqm for complex covarion style process";
}

// Destructor function
CCovEqm::~CCovEqm()	{
	int i;
	FOR(i,(int)m_vpStateProbs.size()) { m_vpStateProbs[i] = NULL; } m_vpStateProbs.~vector();
}

void CCovEqm::AdapterFunc()	{
	int StatesDone = 0, i, j, k,l,DoState;
	// Check some entry conditions
//	cout << "\nm_iChar: " << m_iChar << ", m_iNoCovState: " << m_iNoCovState << ", m_iDataChar: " << m_iDataChar;
//	cout << "\nm_vviEqm2State.size(): " << m_vviEqm2State.size();
	assert(m_iChar == m_iNoCovState * m_iDataChar);
	assert(m_iNoCovState == m_vpStateProbs.size());
	assert(m_vviEqm2State.size() == m_iNoCovObsEqm);
	assert(m_iNoCovObsEqm * m_iDataChar == m_vpEqmPar.size() - m_vpStateProbs.size());
	// Loop through the equilibrium distributions
	// Do the normal style of temporal HMM model
		FOR(i,m_iNoCovObsEqm)	{
			FOR(j,(int)m_vviEqm2State[i].size()) {
				StatesDone++;			// Move counter on
				DoState = m_vviEqm2State[i][j];
//				cout << "\n\tDoing state " << j << "; m_vviEqm2State["<<i<<"]["<<j<<"]: " << m_vviEqm2State[i][j];
				FOR(k,m_iChar) {
					FOR(l,m_iDataChar)	{
						// m_ardQModifier applies to all states
						m_ardQModifier[(k*m_iChar) + (DoState * m_iDataChar) + l] = m_vpEqmPar[(i*m_iDataChar)+l]->Val();
	}	}	}	}
	if(DoBasicNonReversibleCovarion) {
		// Do the covarion style of HMM model which becomes non-reversible if states have different eqm. distributions.
		// The current implementation is inefficient
		FOR(i,m_iChar)	{
			for(j=i+1;j<m_iChar;j++)	{
				// Skip when the hidden state remains the same
				if(i / m_iDataChar == j / m_iDataChar) { continue; }
				m_ardQModifier[(i*m_iChar)+j] = m_vpStateProbs[j / m_iDataChar]->Val();
				m_ardQModifier[(j*m_iChar)+i] = m_vpStateProbs[i / m_iDataChar]->Val();
	}	}	} else {
		FOR(i,m_iChar)	{
			for(j=i+1;j<m_iChar;j++)	{
				// Skip when the hidden state remains the same
				if(i / m_iDataChar == j / m_iDataChar) { continue; }
				m_ardQModifier[(i*m_iChar)+j] *= m_vpStateProbs[j / m_iDataChar]->Val();
				m_ardQModifier[(j*m_iChar)+i] *= m_vpStateProbs[i / m_iDataChar]->Val();
	}	}	}
#if DEBUG_HMP_MODEL == 1
	vector <double> NewSubEqm;
	NewSubEqm = SubEqm(i);
	FOR(i,m_iNoCovState)	{
		if(tdiff(1.0,Sum(&NewSubEqm),MIN_PROB)) {
			cout << "\nCCovEqm::AdapterFunc -- probabilities for hidden state " << i << " don't sum to 1.0\nVector: " << NewSubEqm << " == " << Sum(&NewSubEqm);
		}

	}
#endif
//	cout << "\nStateProbs: "; FOR(i,(int)m_vpStateProbs.size()) { cout << *m_vpStateProbs[i] << " "; }
//	cout << "\nFreqs:"; FOR(i,m_iNoCovState) { cout << "\nState["<<i<<"]: " << CCovEqm::SubEqm(i); }
//	cout << "\nm_ardQModifier: "; MatOut(m_iChar,m_ardQModifier);  exit(-1);
	assert(StatesDone == m_iNoCovState);
}

///////////////////////////////////////////////////////////////////////////////////
// Functions for returning various equilibrium distribution functions
vector <double> CCovEqm::Eqm() {
	vector <double> eq;
	int i,j,l,DoState;
	// Loop through the equilibrium distributions
	eq.assign((int)m_iNoCovState *m_iDataChar,-1.0);
//	cout << "\nSet size: " << m_iNoCovState * m_iDataChar;
	FOR(i,m_iNoCovObsEqm)	{
		FOR(j,(int)m_vviEqm2State[i].size()) {
			DoState = m_vviEqm2State[i][j];
			FOR(l,m_iDataChar)	{
				eq[(DoState * m_iDataChar) + l] = m_vpEqmPar[(i*m_iDataChar)+l]->Val() * m_vpStateProbs[DoState]->Val();
	}	}	}
	return eq;
}
// Returns the equilibrium distribution for a specific hidden state (size == m_iDataChar);
vector <double> CCovEqm::SubEqm(int HiddenState)	{
	int i,State;
	double Total = 0.0;
	vector <double> eq;
	// Check entry conditions
	assert(InRange((int)HiddenState,(int) 0,(int)m_iNoCovState));
	// Do what needs to be done
	FOR(i,(int)m_vviEqm2State.size())	{
		if(IsIn(HiddenState,m_vviEqm2State[i])) { State = i; break; } }
	for(i=State*m_iDataChar;i<(State+1)*m_iDataChar;i++) { Total += m_vpEqmPar[i]->Val(); eq.push_back(m_vpEqmPar[i]->Val()); }
	// Check exit conditions
	assert(i != (int) m_vviEqm2State.size() && diff(Total, (double) 1.0) == 0);
	return eq;
}
// Returns the equilibrium distribution for the C matrix (describes transitions between hidden states; size == m_iNoCovState)
vector <double> CCovEqm::TransEqm()					{
	int i;
	double Total = 0.0;
	vector <double> eq;
	assert(m_iNoCovState == (int) m_vpStateProbs.size());
	FOR(i,m_iNoCovState)	{
		Total += m_vpStateProbs[i]->Val();
		eq.push_back(m_vpStateProbs[i]->Val());
	}
	assert((int)eq.size() == m_iNoCovState && diff(Total,1.0) == 0 );
	return eq;
}

// Shuffles the order of states that the eqm process applies to
void CCovEqm::Shuffle()	{
	int i,j,k;
	vector <int> Vals;
	vector <double> Freqs;
	if((int) m_vviEqm2State.size() == 1) { return; }
	// Get original order of frequencies
	FOR(i,(int) m_vviEqm2State.size())	{ FOR(j,(int) m_vviEqm2State[i].size())	{ Vals.push_back(m_vviEqm2State[i][j]); } }
	random_shuffle(Vals.begin(),Vals.end());
	// Get the new order of frequency values
	FOR(i,(int) Vals.size())	{
		FOR(j,m_iDataChar) { Freqs.push_back(m_vpEqmPar[(Vals[i] * m_iDataChar) + j]->Val()); }
	}
	// Transfer them to the frequencies
	k = 0; FOR(i,(int) Vals.size())	{
		FOR(j,m_iDataChar-1) { m_vpEqmPar[(i* m_iDataChar) + j]->SetVal(Freqs[k++],false,true,false); }
		m_vpEqmPar[(i* m_iDataChar) + (m_iDataChar-1)]->SetVal(Freqs[k++],true,true);
	}
	AdapterFunc();
}

void CCovEqm::ResetEqm(vector <double> Vec, bool RandomFactor)	{
	int i,j;
	vector <double> New, Temp;
	if((int)Vec.size() != m_iDataChar) { Error("Error: CCovEqm::ResetEqm -- vector <double> New has wrong number of states...\n"); }
	if(((int)m_vpEqmPar.size() - m_iNoCovState) % m_iDataChar != 0) { Error("Error: Weird number of m_vpEqmPar in CCovEqm::ResetEqm...\n"); }
	// Create a vector of things that look like the original frequencies, but just a bit differentcd
	FOR(i,(int) (m_vpEqmPar.size() - m_iNoCovState) / m_iDataChar)	{
		if(RandomFactor)	{ FOR(j,m_iDataChar) { Temp.push_back(Vec[j] + (0.1 * Random())); } }
		else				{ FOR(j,m_iDataChar) { Temp.push_back(Vec[j]); } }
		Temp = NormaliseVector(Temp);
		FOR(j,m_iDataChar) { New.push_back(Temp[j]); }
		Temp.clear();
	}
	// Do the assignment of parameters
	FOR(i,(int)m_vpEqmPar.size() - m_iNoCovState)	{
		if(m_vpEqmPar[i]->Special()) { continue; }
		else { m_vpEqmPar[i]->SetVal(New[i],true,true,false); }
	}
	FOR(i,(int)m_vpEqmPar.size() - m_iNoCovState)	{
		if(!m_vpEqmPar[i]->Special()) { continue; }
		else { m_vpEqmPar[i]->SetVal(New[i],true,true,true); }
	}
}

/* *******************************************************************************************
	*
	*		General THMM process object
	*
   ******************************************************************************************* */

// Generalised constructor
CTHMMProcess::CTHMMProcess(CData *D, CTree *T,string Name,int NoCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE eqm_type, bool OptFreq, EHIDDEN_TYPE DoHidden, ERateTypes RateTypes) : CBaseProcess(D,T) {
	int i,j;
	vector <CPar *> EqmPar;
	vector <double> Probs;
	CTHMMSubProc *SubProc = NULL;
	CQPar *Par = NULL;
	m_iHiddenChar = NoCovStates;
	m_bUseGammaRates = false; m_Alfa = NULL;

	// Debug code to make the THMM the same as a normal model
#if HMP_MODEL_2_SIMPLE == 1
	eqm_type = obs; DoHidden = H_none; SeperateRates = false;
#endif

	// Check input functions
	assert(D != NULL && T != NULL);
	// Get the space for the process
	MakeBasicSpace(D->m_iChar * NoCovStates);
	m_iDataChar = D->m_iChar;
	Add_QMat(Name,m_iChar);
	if(OriProbs.empty()) { FOR(i,m_iHiddenChar) { OriProbs.push_back(1.0/(double)m_iHiddenChar); } }
	if(OriProbs.size() != NoCovStates) { Error("\nInternal error in CTHMMProcess::CTHMMProcess OriProbs.size != NoCovStates"); }
	// Make the alphabet
	m_sABET = "\0"; FOR(i,NoCovStates) { m_sABET += DataStates(D->m_DataType); } m_iABET_length = LenStates(D->m_DataType);

	// Create the probability of specific covariate statse
	CreateCovProbs(OriProbs);
	// Add simple Eqm
	switch(eqm_type)	{
	case equ:		AddSimpleCovEqm(&m_vpCovProbs,true); break;	// Add equiprobable equilibrium distribution
	case obs:		AddSimpleCovEqm(&m_vpCovProbs,false,OptFreq); break;	// Add observed equilibrium distribution
	case complex:	AddComplexCovEqm(m_iChar,&m_vpCovProbs,OptFreq); break; // Add multiple equilibrium distributions based on sequences
	default: Error("Unknown type of THMM_EQU_TYPE in CFullDNATHMMProcess::CFullDNATHMMProcess\n\n");
	};
	// Create the zeros
	m_vpPar.push_back(CreateZeros());
	// Create the parameter describing changes between states
	if(DoHidden == H_diff)	{	// Different probabilities between each state
		FOR(i,m_iHiddenChar)	{
			for(j=i+1;j<m_iHiddenChar;j++)	{
				m_vpPar.push_back(CreateChangePar(0.1,i,j));
	}	}	}
	else if(DoHidden == H_same) { m_vpPar.push_back(CreateChangePar(0.1)); }	// Single set of probabilities between all states
	else { m_vpPar.push_back(CreateChangePar(0.0,-1,-1,false)); }
	// Finish by adding all the substitution processes to the model
	FOR(i,m_iHiddenChar)	{
		if(RateTypes == varyall)	{ m_vpPar.push_back(AddSubProc(i,1.0,true)); }
		else						{ AddSubProc(i,1.0,false); }
	}
	if(RateTypes == gammarates) {
		m_bUseGammaRates = true; m_iNoGamCat = m_iHiddenChar;
		m_Alfa = new CQPar("Alpha",m_iChar,0.5,true,1.0E-4); m_vpPar.push_back(m_Alfa);
		NoOptimiseHiddenProbs();
	}
}

// Generalised destructor
CTHMMProcess::~CTHMMProcess() {
	int i;
	// Clean vectors holding pointers to parameters
	FOR(i,(int)m_vpCovProbs.size())	{ m_vpCovProbs[i] = NULL; } m_vpCovProbs.~vector();
	FOR(i,(int)m_vpRPars.size())		{ m_vpRPars[i] = NULL; } m_vpRPars.~vector();
	FOR(i,(int)m_vpSPars.size())		{ m_vpSPars[i] = NULL; } m_vpSPars.~vector();
}

// Controls which parameters describing P(hidden_state) will be optimised
void CTHMMProcess::SetHiddenProbsOpt(bool Optimise,int ProcNum)	{
	int i;
	assert(m_iHiddenChar == (int) m_vpCovProbs.size());
	if(ProcNum == -1) {
		FOR(i,m_iHiddenChar) {
			if(m_vpCovProbs[i]->Special()) { continue; }
			m_vpCovProbs[i]->SetOptimise(Optimise);
		}
	} else {
		assert(InRange(ProcNum,0,m_iHiddenChar));
		if(!m_vpCovProbs[ProcNum]->Special()) { m_vpCovProbs[ProcNum]->SetOptimise(Optimise); }
	}
}

// Function that sets the process rates as gamma distributed rates
bool CTHMMProcess::PrepareQMats(vector <int> Qs2do, bool DoScale)	{
	if(m_bUseGammaRates) {
		int i;
		double freqK[100], rK[100];
		// Prepare the rates
		DiscreteGamma(freqK,rK,m_Alfa->Val(),m_Alfa->Val(),m_iNoGamCat,0);

/*		cout << "\nThere are " << m_iHiddenChar;
		FOR(i,m_iHiddenChar) {
			cout << "\nRate["<<i<<"] = " << rK[i]<< ", prob: " << freqK[i];
		}
*/		FOR(i,m_iHiddenChar) { SetSubProcRate(i, rK[i]); }
		/*
		j = -1; FOR(i,m_iHiddenChar)	{
			if(m_vpCovProbs[i]->Special()) { j = i; continue; }
			else						{ m_vpCovProbs[i]->SetVal(freqK[i],true,true,false); }
		}
		assert (j != -1);
		m_vpCovProbs[j]->SetVal(freqK[j],true,true,true);
*/
//		cout << "\nDoing gamma rates";
//		exit(-1);
	}
	return CBaseProcess::PrepareQMats(Qs2do,DoScale);
}

///////////////////////////////////////////////////////////////////////////////////
// Function to create the parameter for setting the correct off diagonal
//  elements to zero
CQPar *CTHMMProcess::CreateZeros()	{
	int i,j;
	CQPar *Par = NULL;
	Par = new CQPar("Covarion_zeros",m_iChar,0,false,0,0,REPLACE);
	FOR(i,m_iChar) {
		for(j=i+1;j<m_iChar;j++)	{
			if(i / m_iDataChar == j / m_iDataChar) { continue; }	// Ignore the on-diagonal elements of the matrix
			if(m_sABET[i] == m_sABET[j]) { continue; }
			Par->AddQij(i,j,true);
	}	}
	assert(Par != NULL);
	return Par;
}

///////////////////////////////////////////////////////////////////////////////////
// Function to create a parameters for describing all changes between states
// From and To describe the changes between specific states
// Default is From=-1 && To==-1 which means it does all states
CQPar *CTHMMProcess::CreateChangePar(double Value, int From, int To,bool DoOpt)	{
	CCSCovPar *Par = NULL;
	string Name = "Hidden_change(" + int_to_string(From) + "<->" + int_to_string(To) + ")";
	if(To == -1 && From == -1) { Name = "Hidden_change(all)"; }
	Par = new CCSCovPar(Name,m_iChar,m_iHiddenChar,Value,DoOpt); assert(Par != NULL);
	Par->AddCSstates(From,To,m_sABET,m_iABET_length);
	m_vpSPars.push_back(Par);
	return Par;
}

/////////////////////////////////////////////////////////////////////////////////
// Function to create the probability parameters describing the relative
//  occurrence of each covarion state
// NOTE: These are applied through the eqm distribution function
void CTHMMProcess::CreateCovProbs(vector <double> Probs)	{
	int state;
	string Name;
	CQPar *Par = NULL;
	// Check entry conditions
	assert(Probs.size() == m_iHiddenChar);
	// Make the probabilities
	Probs = NormaliseVector(Probs);
	FOR(state,m_iHiddenChar)	{
		// Create and name the probability
		Name = "P(Hidden[" + int_to_string(state) + "])";
		Par = new CQPar(Name,m_iChar,Probs[state],true,MIN_PROB,1.0,MULTIPLY);
		m_vpPar.push_back(Par);
		m_vpCovProbs.push_back(Par);
		// Clean up
		Par = NULL;
}	}

/////////////////////////////////////////////////////////////////////////
// Functions for adding the substitution processes to the model
CQPar * CTHMMProcess::AddSubProc(int ProcNum, double InitRate, bool Optimise)	{
	string Name = "OverallRate[SubstitutionProcess=" + int_to_string(ProcNum) + "]";
	CQPar *RatPar = NULL;
	CTHMMSubProc *Proc;
	RatPar = new CQPar(Name,m_iChar,InitRate,Optimise,0.0,1000,REPLACE);
	if(ProcNum == 0) { RatPar->SetVal(1.0); RatPar->SetSpecial(true); RatPar->SetOptimise(false); }
//	if(ProcNum == m_iHiddenChar - 1) { RatPar->SetVal(1.0); RatPar->SetSpecial(true); RatPar->SetOptimise(false); }
	Proc = new CTHMMSubProc(ProcNum,m_iDataChar,m_vpQMat[0]->Q(),m_iChar,RatPar,m_vpEqm[0]);
	m_vpSubProcs.push_back(Proc);
	Proc = NULL;
	return RatPar;
}
// Add extra zero values to a 2-state process so that second state is fixed
void CTHMMProcess::EnforceFixedState(int State2Do)	{
	int i,j,k;
	assert(m_iHiddenChar == 2);
	FOR(i,(int)m_vpPar.size())	{
		if(m_vpPar[i]->Name() == "Covarion_zeros") {
			for(j=m_iDataChar * State2Do;j<m_iDataChar * (State2Do+1);j++)	{
				for(k=j+1;k<m_iDataChar * (State2Do+1);k++)	{
					m_vpPar[i]->AddQij(j,k);
				}	}
			break;
	}	}
//	PrepareLikelihood(true);
//	cout << "\n\nQ matrix: "; m_vpQMat[0]->OutQ();
	if(i == m_vpPar.size()) { Error("Trying to add to covarion zeros when they don't exist!\n\n"); }
}
/////////////////////////////////////////////////////////////////////////////////
// Some interaction functions for the constituent bits of covarion style models
double CTHMMProcess::SetSubProcRate(int Num,double Rate)	{
	assert(InRange(Num,0,m_iHiddenChar) && Rate >= 0);
	return m_vpSubProcs[Num]->Rate(Rate);
}
void CTHMMProcess::SetOptRates(int ProcNum, bool Val)		{
	int i;
	assert(m_iHiddenChar == (int) m_vpSubProcs.size());
	if(ProcNum == -1) {
		FOR(i,m_iHiddenChar) {
			if(m_vpSubProcs[i]->SpecialRate()) { continue; }
			m_vpSubProcs[i]->SetOptRate(Val);
		}
	} else if(ProcNum != 0) {
		assert(InRange(ProcNum,0,m_iHiddenChar));
		m_vpSubProcs[ProcNum]->SetOptRate(Val);
	}
}
double CTHMMProcess::OverallSubRate(bool ForceWholeProcess)	{
	if((int) m_vpQMat.size() > 1) { Error("Cannot allow temporal hidden Markov processes to have more than one state\n"); }
	return m_vpQMat[0]->OverallSubRate(ForceWholeProcess);
}
double CTHMMProcess::OverallTransRate() {
	if((int) m_vpQMat.size() > 1) { Error("Cannot allow temporal hidden Markov processes to have more than one state\n"); }
	return m_vpQMat[0]->OverallTransRate();
}
/////////////////////////////////////////////////////////////////////////
// Functions for adding Q matrices to the model
CQMat *CTHMMProcess::Add_QMat(string Name,EDataType Type)	{ Error("Standard constructor is unavailable to covarion style processes\n"); return NULL; }
CQMat *CTHMMProcess::Add_QMat(string Name, int Char)		{ assert(m_iDataChar >= 0); CTHMMQMat *Mat; Mat = new CTHMMQMat(m_iDataChar,Char, OTHER, Name); Mat->InitQ(m_dBaseVal); m_vpQMat.push_back(Mat); return Mat; }

/////////////////////////////////////////////////////////////////////////////////
// The rate scaling for the covarion model
void CTHMMProcess::DoPostHocQMatUpdate()	{
	int i,j;
	double Total = 0.0;
	if(DoRates())	{
//		cout << "\nOriginal rates: "; FOR(i,m_iHiddenChar) { cout << SubProc(i)->Rate() << " "; }
		// Always check if rates are equal -- if so then change them just a little
		FOR(i,m_iHiddenChar)	{
			for(j=i+1;j<m_iHiddenChar;j++)	{
				// Do the change
				if(fabs(SubProc(i)->Rate() - SubProc(j)->Rate()) < FLT_EPSILON) {
//					cout << "\nChanging rate ("<<j<<"): " << SubProc(j)->Rate();
					if(!SubProc(j)->SpecialRate())	{
						SubProc(j)->Rate(SubProc(j)->Rate() + RandDouble(-1.0E-5,+1.0E-5));
//						cout << " -> (j="<<j<<")" << SubProc(j)->Rate();
					} else {
						SubProc(i)->Rate(SubProc(i)->Rate() + RandDouble(-1.0E-5,+1.0E-5));
//						cout << " -> (i="<<i<<") " << SubProc(i)->Rate();
				}	}

		}	}
		// Correct scaling parameters
		FOR(i,m_iHiddenChar)	{
			if(SubProc(i)->Rate() > DBL_EPSILON) { Total += SubProc(i)->Rate() * m_vpCovProbs[i]->Val(); }
		}
		// Safely catch zero rates
		if(Total < DBL_EPSILON)	{
			FOR(i,m_iHiddenChar) { SubProc(i)->Rate(0.0); SubProc(i)->Scale(); }
		} else {
			FOR(i,m_iHiddenChar)	{
				SubProc(i)->Scale(SubProc(i)->Rate() / Total);
	}	}	}
}
/* *******************************************************************************************
	*					Basic DNA covarion process
   ******************************************************************************************* */
CDNATHMMProcess::CDNATHMMProcess(CData *D, CTree *T,string Name,int NumCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE eqm_type, bool OptFreq, EHIDDEN_TYPE DoHidden, ERateTypes SeperateRates) : CTHMMProcess(D,T,Name,NumCovStates,OriProbs,eqm_type,OptFreq,DoHidden,SeperateRates) {
	// Check input functions
	assert(D != NULL && T != NULL);
	if(D->m_DataType != DNA) { Error("\nTrying to initialise DNA covarion process with data that doesn't look like nucleotides...\n\n"); }
}
CDNATHMMProcess::~CDNATHMMProcess() {
	cout << "\nStill to fix CDNATHMMProcess()";
}
void CDNATHMMProcess::CreateHKYModel()	{
	int i,j,k,count = 0;
	double Max = -1;
	string Name;
	char seq_i, seq_j;
	CQPar *Par;
	Par = new CQPar("Kappa", m_iChar,INITIAL_KAPPA,true); Par->SetOptimise(true);
	// Correct S_ij to ensure that none are over the maximum value allowed by a parameter
	FOR(i,m_iDataChar) {
		seq_i = GetPos(m_sABET,i,DNA)[0];
		FOR(j,i)	{
			seq_j = GetPos(m_sABET,j,DNA)[0];
			// Only do transitions
			if( (seq_i == 'A' && seq_j == 'G') ||
				(seq_i == 'G' && seq_j == 'A') ||
				(seq_i == 'C' && seq_j == 'T') ||
				(seq_i == 'T' && seq_j == 'C'))	{
				FOR(k,m_iHiddenChar) { Par->AddQij((k*4)+i,(k*4)+j); }
	}	}	}
	m_vpPar.push_back(Par); Par = NULL;
}

// Note: This is the correct function description to use. The other one removes amino acid frequencies...
CDNAWARSProcess::CDNAWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv) : CDNATHMMProcess(D,T,Name,NoGamCat+InRange(pInv,DBL_EPSILON,1.0),vector <double>(),obs,false,H_diff,varyall) {
	m_iNoGamCat = NoGamCat;
	int i;
	bool DoInv = true, DoAlfa = true;
	vector <double> OriProbs;
	string SigName = "SigAlpha";
	// Do some basic initialisation
	m_sName = Name;
	m_Alfa = NULL; m_pInv = NULL; m_SigAlfa = NULL; m_SigInv = NULL;
	if(!InRange(pInv,DBL_EPSILON,1.0)) { DoInv = false; pInv = 0; }
	m_bDoRateScaling = true;

	if(Alfa <= DBL_EPSILON || NoGamCat == 1)	{ DoAlfa = false; if(NoGamCat > 1) { Error("\nError in CDNAWARSProcess::CDNAWARSProcess. Setting gamma alpha to 0, but want rate categories?\n");  } }
	// Set all standard parameters to fixed. Only the new parameters will be optimised
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->SetOptimise(false); }
	// Set up the parameters
	if(DoInv)			{ m_pInv = new CQPar("P(Inv)",m_iChar,pInv,true,MIN_PROB,1.0 - MIN_PROB); m_vpPar.push_back(m_pInv);	}
	if(SigAlfa >= 0)	{
		if(SigInv < DBL_EPSILON) { SigName = "Sigma"; }	// Simple name correction
		m_SigAlfa = new CQPar(SigName,m_iChar,SigAlfa); m_vpPar.push_back(m_SigAlfa);
	}
	if(DoInv)			{ if(SigInv >= DBL_EPSILON)	{ m_SigInv = new CQPar("SigInv",m_iChar,SigInv); m_vpPar.push_back(m_SigInv); }
						  else				{ assert(m_SigAlfa != NULL); m_SigInv = m_SigAlfa; }
	} else if(SigInv >= DBL_EPSILON) { Error("\nCannot have SigmaInv when there are no invariant sites...\n\n"); }
	else { assert(m_SigAlfa != NULL); m_SigInv = m_SigAlfa; }

	if(DoAlfa)			{ m_Alfa = new CQPar("Alpha",m_iChar,Alfa,true,1.0E-4); m_vpPar.push_back(m_Alfa); }
	// This needs to be uncommented for final version. Currently is omitted for simplicity
	CreateHKYModel();
}
// Destructor
CDNAWARSProcess::~CDNAWARSProcess() {
	m_Alfa = NULL; m_SigAlfa = NULL; m_SigInv = NULL; m_pInv = NULL;
}
// Do the parameterisations
void CDNAWARSProcess::DoParameterisation()	{
	int i,j;

	double freqK[100], rK[100], pInv = 0,value;
	bool DoInv = false;
	vector <double> Probs(m_iHiddenChar,0), Rates(m_iHiddenChar,0);
	vector <int> Num;
	if(m_pInv != NULL) { pInv = m_pInv->Val(); }
	// Prepare the rates
	if(m_iNoGamCat > 1) { DiscreteGamma(freqK,rK,m_Alfa->Val(),m_Alfa->Val(),m_iNoGamCat,0); }
	else { rK[0] = 1; }
	FOR(i,m_iNoGamCat) {
		if(rK[0] > 1.0E-3) { Rates[i] = rK[i] / rK[0]; }
		else { rK[0] = m_iNoGamCat; for(j=1;j<m_iNoGamCat;j++) { rK[j] = 0.0; } }
		Probs[i] = (1.0 - pInv) / m_iNoGamCat;
	}
	if(Rates[m_iNoGamCat-1] > 500) { value = 500 / Rates[m_iNoGamCat-1]; FOR(i,m_iNoGamCat) { Rates[i] = Rates[i] * value; } }
	if(m_iHiddenChar == m_iNoGamCat + 1) { DoInv = true; Rates[m_iNoGamCat] = 0; Probs[m_iNoGamCat]= pInv; }
	// Apply the rates and probs
	FOR(i,m_iHiddenChar) { if(Probs[i] < 1.0E-3) { Probs[i] = 1.0E-3; } }
	Probs = NormaliseVector(Probs);
	FOR(i,m_iHiddenChar) { SetSubProcRate(i, Rates[i]); }
	j = -1; FOR(i,m_iHiddenChar)	{
		if(m_vpCovProbs[i]->Special()) { j = i; continue; }
		else						{ m_vpCovProbs[i]->SetVal(Probs[i],true,true,false); }
	}
	assert (j != -1);
	m_vpCovProbs[j]->SetVal(Probs[j],true,true,true);
	// Sort out the switching rates
	FOR(i,(int) m_vpSPars.size()) {
		assert(m_vpSPars[i]->Name().find("Hidden_change") != string::npos);
		if(DoInv) {
			if(IsIn(m_iHiddenChar-1,m_vpSPars[i]->To()) || IsIn(m_iHiddenChar-1,m_vpSPars[i]->From())) {
				if(m_SigInv != NULL) { m_vpSPars[i]->SetVal(m_SigInv->Val()); } else { m_vpSPars[i]->SetVal(0.0); }
				continue;
		}	}
		if(m_SigAlfa != NULL) { m_vpSPars[i]->SetVal(m_SigAlfa->Val()); } else { m_vpSPars[i]->SetVal(0.0); }
	}
}
bool CDNAWARSProcess::PrepareQMats(vector <int> Qs2do, bool DoScale) {
	DoParameterisation();
	return CBaseProcess::PrepareQMats(Qs2do,DoScale);
}

/* *******************************************************************************************
		The branchwise covarion process
   ******************************************************************************************* */

CDNABranchWARSProcess::CDNABranchWARSProcess(CData *D, CTree *T, string ModelName, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv, bool VarySig, bool VarypInv) : CDNAWARSProcess(D,T,ModelName,Alfa, NoGamCat,SigAlfa, pInv, SigInv) {
	string SigName = "SigAlpha";
	CQPar *Par = NULL;
	m_RootInv = NULL;
	int i;
	string Name;
	bool DoSigInv = true;
	if(fabs(SigInv) < DBL_EPSILON) { DoSigInv = false; }

	// Set flags appropriately
//	m_bDoStandardDecompose = false;
	// Input checking
	if(!VarySig && !VarypInv) { Error("\nCalling CDNABranchWARSProcess::CDNABranchWARSProcess(...) for nothing varying across branches...\n"); }
	// Do the sigma varying per branch
	if(!DoSigInv) { SigName = "Sigma"; }
	if(T == NULL) { Error("Trying to initialise CDNABranchWARSProcess class without a starting tree...\n"); }
	if(m_SigAlfa != NULL && VarySig) {
		m_SigAlfa->Name(SigName + "[0]"); m_BraSigAlfa.push_back(m_SigAlfa);
		FOR(i,m_pTree->NoBra() - 1) {
			Name = SigName + "[" + int_to_string(i+1) + "]";
			 Par = new CQPar(Name,m_iChar,m_SigAlfa->Val() + (Random() * 0.05));
			//Par = new CQPar(Name,m_iChar,m_SigAlfa->Val()); // HAXORED
			m_BraSigAlfa.push_back(Par);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}
	if(DoSigInv && VarySig) {
		assert(m_SigInv != NULL);
		m_SigInv->Name("SigInv[0]"); m_BraSigInv.push_back(m_SigInv);
		FOR(i,m_pTree->NoBra() - 1) {
			Name = "SigInv[" + int_to_string(i+1) + "]";
			Par = new CQPar(Name,m_iChar,m_SigInv->Val() + (Random() * 0.05));
			m_BraSigInv.push_back(Par);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}

	// Do the invariant sites varying across branches
	if(m_pInv != NULL) {
		if(VarypInv) {
			m_RootInv = new CQPar("Root pInv", m_iChar, 0.2,true,1.0E-4,1.0 - 1.0E-4);
			m_vpPar.push_back(m_RootInv);
			m_pInv->Name("pInv[0]"); m_BrapInv.push_back(m_pInv);
			FOR(i,m_pTree->NoBra() - 1) {
				Name = "pInv[" + int_to_string(i+1) + "]";
				Par = new CQPar(Name,m_iChar,m_pInv->Val() + (Random() * 0.05),true,1.0E-4,1.0 - 1.0E-4);
				m_BrapInv.push_back(Par);
				m_vpPar.push_back(Par);
				Par = NULL;
	}	}	} else if(VarypInv) {
		Error("\nTrying to vary invariant sites across branches when no invariant sites allowed...\n");
	}
/*
	cout << "\nThe parameters are: ";
	FOR(i,(int)m_vpPar.size()) {
		cout << "\nPar["<<i<<"]: " << *m_vpPar[i];
	}
*/
}
// Destructor function
CDNABranchWARSProcess::~CDNABranchWARSProcess() {
	int i;
	FOR(i,(int)m_BraSigAlfa.size())	{ m_BraSigAlfa[i] = NULL; }
	FOR(i,(int)m_BraSigInv.size())	{ m_BraSigInv[i] = NULL; }
	FOR(i,(int)m_BrapInv.size())	{ m_BrapInv[i] = NULL; }

}
// Make_PT function that does the branch-by-branch matrix details
bool CDNABranchWARSProcess::Make_PT(int B, bool RedoRate)	{
	int i,UseB = B;
//	cout << "\nMake_PT("<<B<<","<<RedoRate<<")";
	if(m_pTree->BranchLabels()) { // Do grouped branches
		UseB = m_pTree->Labels()[B];
		// Do some checking
		assert(InRange(UseB,0,m_pTree->NoLabels()));
	} else { // Do one parameter per branch
		// Do some checking
		assert(InRange(UseB,0,m_pTree->NoBra()));
	}
	// Do the m_BraSigAlfa
	if(!m_BraSigAlfa.empty()) {
		assert(InRange(UseB,0,(int)m_BraSigAlfa.size()));
		m_SigAlfa = m_BraSigAlfa[UseB];
		if(m_BraSigInv.empty() && m_SigInv != NULL) {
			m_SigInv = m_SigAlfa;
		}
//		cout << "\n\tUsing m_SigAlfa " << *m_SigAlfa;
	}
	// Do the m_BraSigInv
	if(!m_BraSigInv.empty()) {
		assert(InRange(UseB,0,(int)m_BraSigInv.size()));
		m_SigInv = m_BraSigInv[UseB];
//		cout << "\n\tUsing m_SigInv " << *m_SigInv;
	}
	// Do the m_BrapInv
	if(!m_BrapInv.empty())	{
		assert(InRange(UseB,0,(int)m_BrapInv.size()));
		m_pInv = m_BrapInv[UseB];
//		cout << "\n\tUsing m_pInv " << *m_pInv;
	}
	// Recreate the Q matrix
	if(!PrepareQMats()) { Error("CBranchWARSProcess::Make_PT(...) PrepareQMats failed...."); }
	FOR(i,(int)m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); }
//	cout << "\nAnd using QMat: " << GetQMatSums();
//	m_vpQMat[0]->OutQ();

	bool tempo = CDNAWARSProcess::Make_PT(B,RedoRate);
//	cout << "\n --- The rates: OverallSubRate: " << OverallSubRate() << "; OverallTransRate: " << OverallTransRate() << " --- ";
	return tempo;
}
// Prepare Q matrices for Branch Derivatives
void CDNABranchWARSProcess::PrepareBraDer() {
	int i;
	if(!m_pTree->BranchLabels()) { Error("\nTrying to CBranchWARSProcess::PrepareBraDer when the tree hasn't been labelled...\n\n "); }
	// Prepare partial likelihoods
	Likelihood(true);
	FOR(i,m_pTree->NoBra()) {
		// Do the PT matrix
		Make_PT(i);
		// Calculate the QP Matrix
		MulMat(m_vpQMat[QMat4Bra(i)]->Q(),PT(i),QP(i),m_iChar,m_iChar,m_iChar);
	}
	m_bBraDerReady = true;
}

// The output function for the CDNABranchWARSProcess
ostream &CDNABranchWARSProcess::Output(ostream &os)	{
	int i;
	double t;
	m_pTree->OutLabel();
	CDNATHMMProcess::Output(os);
	CTree TestTree; TestTree = *m_pTree;
	cout << "\nDoing output function for CDNABranchWARSProcess";
	os << "\nRates of switching between classes by branch [subrate; switchrate; branchlength; switchlength]";
	FOR(i,m_pTree->NoBra())	{
		Make_PT(i);
		t = OverallTransRate();
		os << "\n\tBranch["<<i<<"]: " << OverallSubRate() << "\t" << t << "\t" << m_pTree->B(i) << "\t" << m_pTree->B(i) * t;
		TestTree.SetB(i,OverallTransRate() * m_pTree->B(i));
	}
	TestTree.OutBra();
	TestTree.OutName();
	os << "\nSwitchTree:\t" << TestTree;
#if ALLOW_BDNA_PEN == 1
	os << "\nPenalty:\t" << Penalty();
#endif
	return os;
}

#define DO_EQM_CHECK 1		// A debugging function
// The class specific function for getting the root eqm
vector <double> CDNABranchWARSProcess::RootEqm()	{
	int i,j;
	vector <double> E, DataFreq;
	double Mul;
	// Checkers
	assert(m_pData != NULL);
	// Make the equilibrium
	DataFreq = m_pData->m_vFreq;
	assert(!diff(Sum(&DataFreq),1.0));
	assert((int) DataFreq.size() == m_iDataChar);
	// Do the Gamma stuff
	if(m_RootInv != NULL) {
		assert(m_iHiddenChar == m_iNoGamCat);
		Mul = (1.0 - m_RootInv->Val()) / m_iNoGamCat;
	} else { Mul = 1.0 / (double) m_iNoGamCat; }
	FOR(i,m_iNoGamCat) { FOR(j,m_iDataChar) { E.push_back(DataFreq[j] * Mul); } }
	// Do the Inv stuff
	if(m_RootInv != NULL) {
		assert(m_iHiddenChar == m_iNoGamCat + 1);
		Mul = m_RootInv->Val();
		FOR(j,m_iDataChar) { E.push_back(DataFreq[j] * Mul); }
	}
	// Exit checker
	if(abs(Sum(&E) - 1.0) < 1.0E-4) { E = NormaliseVector(E); }
	else { Error("\nEquilibrium vector doesn't sum to 1.0 in vector <double> CDNABranchWARSProcess::RootEqm()"); }
//	cout << "\nReturning root eqm: " << E;
	return E;
}

// ----------- Functions for doing penalised likelihood -------------
// Penalty of the form: Penalty = (1 - exp(-abs(a-b)) for each comparison a b.
// Penalty is calculated using a traversal on a tree, with the penalty on each branch
// calculated by the difference between its value and the 4 attached branches. This can
// be efficiently calculated by examining the three branches for each internal node and taking the
// three possible pairwise comparisons.
// Therefore, for an unrooted tree there are n-2 internal branches => 3n -6 pairwise comparisons
// Penalty is calculated for each value in the vectors m_BraSigAlfa, m_BraSigInv, and m_BrapInv
// The penalty is scaled for *EACH* vector independently to MaxPen (defined in front of model).
// This means the total penalty is 3 * MaxPen if all parameter are present
double CDNABranchWARSProcess::Penalty() {
	double Pen = 0.0;
	Pen = PenaltyTraverse(m_pTree->StartCalc(),-1) * CBranchWARSProcess_PenaltyIntensity;
	if(!m_BraSigAlfa.empty()) { Pen /= (3*m_pTree->NoSeq()) - 6; }
	if(!m_BraSigInv.empty())  { Pen /= (3*m_pTree->NoSeq()) - 6; }
	if(!m_BrapInv.empty())    { Pen /= (3*m_pTree->NoSeq()) - 6; }
//	cout << "\nReturning penalty: "<< Pen; exit(-1);
	return Pen;
}
double CDNABranchWARSProcess::PenaltyTraverse(int iNoTo, int iNoFr) {
	int i,j;
	double Penalty = 0.0;
//	cout << "\nVisiting branchnode["<<iNoTo<<"]";
	// Do post order tree traversal
	FOR(i,m_pTree->NoLinks(iNoTo)) {
		if(m_pTree->NodeLink(iNoTo,i) == iNoFr || m_pTree->NodeLink(iNoTo,i) == -1) { continue; }		// Skip the node it came from
		Penalty += PenaltyTraverse(m_pTree->NodeLink(iNoTo,i),iNoTo);	// Traverse the tree
	}
	if(m_pTree->NodeType(iNoTo) == branch) {
//		cout << "\nCalculating branchnode["<<iNoTo<<"]";
		FOR(i,m_pTree->NoLinks(iNoTo)) {
			for(j=i+1;j<m_pTree->NoLinks(iNoTo);j++) {
//				cout << "\n\t["<<i<<"] " << m_pTree->NodeBra(iNoTo,i) << " : ["<<j<<"] " << m_pTree->NodeBra(iNoTo,j) << flush;
//				cout << " == " <<1 - exp(fabs(m_BraSigAlfa[m_pTree->NodeBra(iNoTo,i)]->Val() - m_BraSigAlfa[m_pTree->NodeBra(iNoTo,j)]->Val())) << flush;
				if(!m_BraSigAlfa.empty()) { Penalty += 1 - exp(fabs(m_BraSigAlfa[m_pTree->NodeBra(iNoTo,i)]->Val() - m_BraSigAlfa[m_pTree->NodeBra(iNoTo,j)]->Val())); }
				if(!m_BraSigInv.empty())  { Penalty += 1 - exp(fabs(m_BraSigInv[m_pTree->NodeBra(iNoTo,i)]->Val() - m_BraSigInv[m_pTree->NodeBra(iNoTo,j)]->Val())); }
				if(!m_BrapInv.empty())    { Penalty += 1 - exp(fabs(m_BrapInv[m_pTree->NodeBra(iNoTo,i)]->Val() - m_BrapInv[m_pTree->NodeBra(iNoTo,j)]->Val())); }
//				cout << " *" << flush;
		}	}
	}
	return Penalty;
}


/* *******************************************************************************************
	*					Basic Amino Acid covarion process
   ******************************************************************************************* */
CAATHMMProcess::CAATHMMProcess(CData *D, CTree *T,string Name,int NumCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE eqm_type, bool OptFreq, EHIDDEN_TYPE DoHidden, ERateTypes SeperateRates) : CTHMMProcess(D,T,Name,NumCovStates,OriProbs,eqm_type,OptFreq,DoHidden,SeperateRates) {
	// Check input functions
	assert(D != NULL && T != NULL);
	if(D->m_DataType != AA) { Error("\nTrying to initialise AA covarion process with data that doesn't look like amino acids...\n\n"); }
}
CAATHMMProcess::~CAATHMMProcess() {
	cout << "\nStill to fix CAATHMMProcess()";
}
void CAATHMMProcess::CreateEMPModel(vector <double> S_ij,vector <double> Freq,bool AddF)	{
	int i,j,k,count = 0;
	double Max = -1;
	string Name;
	vector <double> F;
	CQPar *Par;
	assert((int)S_ij.size() == 190 );

	// Correct S_ij to ensure that none are over the maximum value allowed by a parameter
	FOR(i,m_iDataChar) {
		FOR(j,i)	{
			Name = GetPos(m_sABET,i,AA) + " <-> " + GetPos(m_sABET,j,AA);
			Par = new CQPar(Name,m_iChar,S_ij[count++],false,0.0,BIG_NUMBER);
			// Loop through any hidden states to fully set up the matrix
			FOR(k,m_iHiddenChar) { Par->AddQij((k*20)+i,(k*20)+j); }
			Par->SetOptimise(false);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}
}

///////////////////////////////// The WARS Model ///////////////////////////////////////////////

// DEBUG FUNCTION HEADER WITH EQUIPROBABLE AA FREQS
//CWARSProcess::CWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv) : CAATHMMProcess(D,T,Name,NoGamCat+InRange(pInv,DBL_EPSILON,1.0),vector <double>(),equ,false,H_diff,true) {

// Note: This is the correct function description to use. The other one removes amino acid frequencies...
 CWARSProcess::CWARSProcess(CData *D, CTree *T,string Name, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv) : CAATHMMProcess(D,T,Name,NoGamCat+InRange(pInv,DBL_EPSILON,1.0),vector <double>(),obs,false,H_diff,varyall) {
	m_iNoGamCat = NoGamCat;
	int i;
	bool DoInv = true, DoAlfa = true;
	vector <double> OriProbs;
	string SigName = "SigAlpha";
	// Do some basic initialisation
	m_sName = Name;
	m_Alfa = NULL; m_pInv = NULL; m_SigAlfa = NULL; m_SigInv = NULL;
	if(!InRange(pInv,DBL_EPSILON,1.0)) { DoInv = false; pInv = 0; }
	m_bDoRateScaling = true;

	if(Alfa <= DBL_EPSILON || NoGamCat == 1)	{ DoAlfa = false; if(NoGamCat > 1) { Error("\nError in CWARSProcess::CWARSProcess. Setting gamma alpha to 0, but want rate categories?\n");  } }
	// Set all standard parameters to fixed. Only the new parameters will be optimised
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->SetOptimise(false); }
	// Set up the parameters
	if(DoInv)			{ m_pInv = new CQPar("P(Inv)",m_iChar,pInv,true,MIN_PROB,1.0 - MIN_PROB); m_vpPar.push_back(m_pInv);	}
	if(SigAlfa >= 0)	{
		if(SigInv < DBL_EPSILON) { SigName = "Sigma"; }	// Simple name correction
		m_SigAlfa = new CQPar(SigName,m_iChar,SigAlfa); m_vpPar.push_back(m_SigAlfa);
	}
	if(DoInv)			{ if(SigInv >= DBL_EPSILON)	{ m_SigInv = new CQPar("SigInv",m_iChar,SigInv); m_vpPar.push_back(m_SigInv); }
						  else				{ assert(m_SigAlfa != NULL); m_SigInv = m_SigAlfa; }
	} else if(SigInv >= DBL_EPSILON) { Error("\nCannot have SigmaInv when there are no invariant sites...\n\n"); }
	else { assert(m_SigAlfa != NULL); m_SigInv = m_SigAlfa; }

	if(DoAlfa)			{ m_Alfa = new CQPar("Alpha",m_iChar,Alfa,true,1.0E-4); m_vpPar.push_back(m_Alfa); }
	// This needs to be uncommented for final version. Currently is omitted for simplicity
	CreateEMPModel(vWAGVal(),vWAGFreq(),false);
//	CreateEMPModel(vcpREVVal(),vcpREVFreq(),false);
}
// Destructor
CWARSProcess::~CWARSProcess() {
	cout << "\nStill to fix CWARSProcess()";
	m_Alfa = NULL; m_SigAlfa = NULL; m_SigInv = NULL; m_pInv = NULL;
}
// Do the parameterisations
void CWARSProcess::DoParameterisation()	{
	int i,j;
	double freqK[100], rK[100], pInv = 0;
	bool DoInv = false;
	vector <double> Probs(m_iHiddenChar,0), Rates(m_iHiddenChar,0);
	vector <int> Num;
	if(m_pInv != NULL) { pInv = m_pInv->Val(); }
	// Prepare the rates
	if(m_iNoGamCat > 1) { DiscreteGamma(freqK,rK,m_Alfa->Val(),m_Alfa->Val(),m_iNoGamCat,0); }
	else { rK[0] = 1; }
	FOR(i,m_iNoGamCat) {
		if(rK[0] > 1.0E-3) { Rates[i] = rK[i] / rK[0]; }
		else { rK[0] = m_iNoGamCat; for(j=1;j<m_iNoGamCat;j++) { rK[j] = 0.0; } }
		Probs[i] = (1.0 - pInv) / m_iNoGamCat;
	}
	if(m_iHiddenChar == m_iNoGamCat + 1) { DoInv = true; Rates[m_iNoGamCat] = 0; Probs[m_iNoGamCat]= pInv; }
	// Apply the rates and probs
	FOR(i,m_iHiddenChar) { if(Probs[i] < 1.0E-3) { Probs[i] = 1.0E-3; } }
	Probs = NormaliseVector(Probs);
	FOR(i,m_iHiddenChar) { SetSubProcRate(i, Rates[i]); }
	j = -1; FOR(i,m_iHiddenChar)	{
		if(m_vpCovProbs[i]->Special()) { j = i; continue; }
		else						{ m_vpCovProbs[i]->SetVal(Probs[i],true,true,false); }
	}
	assert (j != -1);
	m_vpCovProbs[j]->SetVal(Probs[j],true,true,true);
	// Sort out the switching rates
	FOR(i,(int) m_vpSPars.size()) {
		assert(m_vpSPars[i]->Name().find("Hidden_change") != string::npos);
		if(DoInv) {
			if(IsIn(m_iHiddenChar-1,m_vpSPars[i]->To()) || IsIn(m_iHiddenChar-1,m_vpSPars[i]->From())) {
				if(m_SigInv != NULL) { m_vpSPars[i]->SetVal(m_SigInv->Val()); } else { m_vpSPars[i]->SetVal(0.0); }
				continue;
		}	}
		if(m_SigAlfa != NULL) { m_vpSPars[i]->SetVal(m_SigAlfa->Val()); } else { m_vpSPars[i]->SetVal(0.0); }
	}
}
bool CWARSProcess::PrepareQMats(vector <int> Qs2do, bool DoScale) {
	DoParameterisation();
	return CBaseProcess::PrepareQMats(Qs2do,DoScale);
}

ostream &CWARSProcess::Output(ostream &os)	{
	int i;
	double t;
	CAATHMMProcess::Output(os);
	CTree TestTree; TestTree = *m_pTree;
	cout << "\nDoing output function for CBranchWARSProcess";
	os << "\nRates of switching between classes by branch [subrate; switchrate; branchlength; switchlength]";
	FOR(i,m_pTree->NoBra())	{
		Make_PT(i);
		t = OverallTransRate();
		os << "\n\tBranch["<<i<<"]: " << OverallSubRate() << "\t" << t << "\t" << m_pTree->B(i) << "\t" << m_pTree->B(i) * t;
		TestTree.SetB(i,OverallTransRate() * m_pTree->B(i));
	}
	TestTree.OutBra();
	TestTree.OutName();
	os << "\nSwitchTree:\t" << TestTree;
#if ALLOW_BAA_PEN ==1
	os << "\nPenalty:\t" << Penalty();
#endif
	return os;
}

/* *******************************************************************************************
		The branchwise covarion process
   ******************************************************************************************* */

CBranchWARSProcess::CBranchWARSProcess(CData *D, CTree *T, string ModelName, double Alfa, int NoGamCat, double SigAlfa, double pInv, double SigInv, bool VarySig, bool VarypInv) : CWARSProcess(D,T,ModelName,Alfa, NoGamCat,SigAlfa, pInv, SigInv) {
	string SigName = "SigAlpha";
	CQPar *Par = NULL;
	m_RootInv = NULL;
	int i;
	string Name;
	bool DoSigInv = true;
	if(fabs(SigInv) < DBL_EPSILON) { DoSigInv = false; }
	vector <int> Labels;
	// Set flags appropriately
//	m_bDoStandardDecompose = false;
	// Input checking
	if(!VarySig && !VarypInv) { Error("\nCalling CBranchWARSProcess::CBranchWARSProcess(...) for nothing varying across branches...\n"); }
	// Do the sigma varying per branch
	if(!DoSigInv) { SigName = "Sigma"; }
	if(T == NULL) { Error("Trying to initialise CBranchWARSProcess class without a starting tree...\n"); }
	if(m_SigAlfa != NULL && VarySig) {
		if(!m_pTree->BranchLabels()) { m_pTree->CreateBranchLabels(true); }
		m_SigAlfa->Name(SigName + "[0]"); m_BraSigAlfa.push_back(m_SigAlfa);
		FOR(i,m_pTree->NoLabels() - 1) {
			Name = SigName + "[" + int_to_string(i+1) + "]";
			Par = new CQPar(Name,m_iChar,m_SigAlfa->Val() + (Random() * 0.05));
			//Par = new CQPar(Name,m_iChar,m_SigAlfa->Val()); // HAXORED
			m_BraSigAlfa.push_back(Par);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}
	if(DoSigInv && VarySig) {
		if(m_pTree->BranchLabels()) {
			cout << "\nNeed to do something clever!!!";
			exit(-1);
		} else {
			assert(m_SigInv != NULL);
			m_SigInv->Name("SigInv[0]"); m_BraSigInv.push_back(m_SigInv);
			FOR(i,m_pTree->NoBra() - 1) {
				Name = "SigInv[" + int_to_string(i+1) + "]";
				Par = new CQPar(Name,m_iChar,m_SigInv->Val() + (Random()
					* 0.05));
				m_BraSigInv.push_back(Par);
				m_vpPar.push_back(Par);
				Par = NULL;
			}
	}	}
	// Do the invariant sites varying across branches
	if(m_pInv != NULL) {
		if(VarypInv) {
			m_RootInv = new CQPar("Root pInv", m_iChar, 0.2,true,MIN_PROB,1.0 - MIN_PROB);
			m_vpPar.push_back(m_RootInv);
			m_pInv->Name("pInv[0]"); m_BrapInv.push_back(m_pInv);
			FOR(i,m_pTree->NoBra() - 1) {
				Name = "pInv[" + int_to_string(i+1) + "]";
				Par = new CQPar(Name,m_iChar,m_pInv->Val() + (Random() * 0.05),true,MIN_PROB,1.0 - MIN_PROB);
				m_BrapInv.push_back(Par);
				m_vpPar.push_back(Par);
				Par = NULL;
			}
		} else {
			m_RootInv = m_pInv;
	}	} else if(VarypInv) {
		Error("\nTrying to vary invariant sites across branches when no invariant sites allowed...\n");
	}
/*
	cout << "\nThe parameters are: ";
	FOR(i,(int)m_vpPar.size()) {
		cout << "\nPar["<<i<<"]: " << *m_vpPar[i];
	}
*/

}
// Destructor function
CBranchWARSProcess::~CBranchWARSProcess() {
	int i;
	m_RootInv = NULL;
	FOR(i,(int)m_BraSigAlfa.size())	{ m_BraSigAlfa[i] = NULL; }
	FOR(i,(int)m_BraSigInv.size())	{ m_BraSigInv[i] = NULL; }
	FOR(i,(int)m_BrapInv.size())	{ m_BrapInv[i] = NULL; }
}
// Make_PT function that does the branch-by-branch matrix details
bool CBranchWARSProcess::Make_PT(int B, bool RedoRate)	{
	int i,UseB = B;
//	cout << "\nCBranchWARSProcess::Make_PT("<<B<<","<<RedoRate<<"): Branch = " << m_pTree->B(B);
	if(m_pTree->BranchLabels()) { // Do grouped branches
		UseB = m_pTree->Labels()[B];
		// Do some checking
		assert(InRange(UseB,0,m_pTree->NoLabels()));
	} else { // Do one parameter per branch
		Error("\nCalculations for CBranchWARSProcess::Make_PT(...) requires the tree have labelled branches...\n");
	}
	// Do the m_BraSigAlfa
	if(!m_BraSigAlfa.empty()) {
		assert(InRange(UseB,0,(int)m_BraSigAlfa.size()));
		m_SigAlfa = m_BraSigAlfa[UseB];
		if(m_BraSigInv.empty() && m_SigInv != NULL) {
			m_SigInv = m_SigAlfa;
		}
//		cout << "\n\tUsing m_SigAlfa " << *m_SigAlfa;
	}
	// Do the m_BraSigInv
	if(!m_BraSigInv.empty()) {
		assert(InRange(UseB,0,(int)m_BraSigInv.size()));
		m_SigInv = m_BraSigInv[UseB];
//		cout << "\n\tUsing m_SigInv " << *m_SigInv;
	}
	// Do the m_BrapInv
	if(!m_BrapInv.empty())	{
		assert(InRange(UseB,0,(int)m_BrapInv.size()));
		m_pInv = m_BrapInv[UseB];
//		cout << "\n\tUsing m_pInv " << *m_pInv;
	}
	// Recreate the Q matrix
	if(!PrepareQMats()) { Error("CBranchWARSProcess::Make_PT(...) PrepareQMats failed...."); }
	FOR(i,(int)m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); }
//	cout << "\nAnd using QMat: " << GetQMatSums();
//	if(B == 2) { cout << " ... Sigma: " << *m_SigAlfa; m_vpQMat[0]->OutQ(); }
	bool tempo = CBaseProcess::Make_PT(B,RedoRate);
//	cout << "\n --- The rates: OverallSubRate: " << OverallSubRate() << "; OverallTransRate: " << OverallTransRate() << " --- ";
	return tempo;
}
void CBranchWARSProcess::PrepareBraDer() {
	int i;
	if(!m_pTree->BranchLabels()) { Error("\nTrying to CBranchWARSProcess::PrepareBraDer when the tree hasn't been labelled...\n\n "); }
	// Prepare partial likelihoods
	Likelihood(true);
	FOR(i,m_pTree->NoBra()) {
		// Do the PT matrix
		Make_PT(i);
		// Calculate the QP Matrix
		MulMat(m_vpQMat[QMat4Bra(i)]->Q(),PT(i),QP(i),m_iChar,m_iChar,m_iChar);
	}
	m_bBraDerReady = true;
}
// The output function for the CBranchWARSProcess
ostream &CBranchWARSProcess::Output(ostream &os)	{
	int i;
	double t;
	CAATHMMProcess::Output(os);
	CTree TestTree; TestTree = *m_pTree;
	cout << "\nDoing output function for CBranchWARSProcess";
	os << "\nRates of switching between classes by branch [subrate; switchrate; branchlength; switchlength]";
	FOR(i,m_pTree->NoBra())	{
		Make_PT(i);
		t = OverallTransRate();
		os << "\n\tBranch["<<i<<"]: " << OverallSubRate() << "\t" << t << "\t" << m_pTree->B(i) << "\t" << m_pTree->B(i) * t;
		TestTree.SetB(i,OverallTransRate() * m_pTree->B(i));
	}
	TestTree.OutBra();
	TestTree.OutName();
	os << "\nSwitchTree:\t" << TestTree;
#if ALLOW_BAA_PEN ==1
	os << "\nPenalty:\t" << Penalty();
#endif
	return os;
}

#define DO_EQM_CHECK 1		// A debugging function
// The class specific function for getting the root eqm
vector <double> CBranchWARSProcess::RootEqm()	{
	int i,j;
	vector <double> E, DataFreq;
	double Mul;
	// Checkers
	assert(m_pData != NULL); \
	// Make the equilibrium
	DataFreq = m_pData->m_vFreq;
	assert(tdiff(Sum(&DataFreq),1.0,1.0E-5) == 0);
	assert((int) DataFreq.size() == m_iDataChar);
	// Do the Gamma stuff
	if(m_RootInv != NULL) {
		assert(m_iHiddenChar == m_iNoGamCat + 1);
		Mul = (1.0 - m_RootInv->Val()) / m_iNoGamCat;
	} else { Mul = 1.0 / (double) m_iNoGamCat; }
	FOR(i,m_iNoGamCat) { FOR(j,m_iDataChar) { E.push_back(DataFreq[j] * Mul); } }
	// Do the Inv stuff
	if(m_RootInv != NULL) {
		assert(m_iHiddenChar == m_iNoGamCat + 1);
		Mul = m_RootInv->Val();
		FOR(j,m_iDataChar) { E.push_back(DataFreq[j] * Mul); }
	}
	// Exit checker
	if(abs(Sum(&E) - 1.0) < 1.0E-4) { E = NormaliseVector(E); }
	else { Error("\nEquilibrium vector doesn't sum to 1.0 in vector <double> CBranchWARSProcess::RootEqm()"); }
//	cout << "\nReturning root eqm["<<E.size()<<"]: \n" << E;
	return E;
}

// A function that shuffles the parameters around in the sub processes in the THMM
void CBranchWARSProcess::ShuffleParameters()	{
	int i,j,k;
	vector <double> Par;
	vector <double> New;
	double Val;
	// Do the branch specific values of sigma
	FOR(i,m_BraSigAlfa.size()) { m_BraSigAlfa[i]->SetVal(RandDouble(0.5,1.5),true,true); }
	// Do the other parameters
	FOR(i,m_vpPar.size()) {
		if(m_vpPar[i]->Name().find("P(Inv)") != string::npos) { m_vpPar[i]->SetVal(RandDouble(0.25,0.75),true,true,true); }
		if(m_vpPar[i]->Name().find("Alpha") != string::npos) { m_vpPar[i]->SetVal(RandDouble(0.2,1.5),true,true); }
	}
}

// ----------- Functions for doing penalised likelihood -------------
// Penalty of the form: Penalty = (1 - exp(-abs(a-b)) for each comparison a b.
// Penalty is calculated using a traversal on a tree, with the penalty on each branch
// calculated by the difference between its value and the 4 attached branches. This can
// be efficiently calculated by examining the three branches for each internal node and taking the
// three possible pairwise comparisons.
// Therefore, for an unrooted tree there are n-2 internal branches => 3n -6 pairwise comparisons
// Penalty is calculated for each value in the vectors m_BraSigAlfa, m_BraSigInv, and m_BrapInv
// The penalty is scaled for *EACH* vector independently to MaxPen (defined in front of model).
// This means the total penalty is 3 * MaxPen if all parameter are present
double CBranchWARSProcess::Penalty() {
	double Pen = 0.0;
	Pen = PenaltyTraverse(m_pTree->StartCalc(),-1) * CBranchWARSProcess_PenaltyIntensity;
	if(!m_BraSigAlfa.empty()) { Pen /= (3*m_pTree->NoSeq()) - 6; }
	if(!m_BraSigInv.empty())  { Pen /= (3*m_pTree->NoSeq()) - 6; }
	if(!m_BrapInv.empty())    { Pen /= (3*m_pTree->NoSeq()) - 6; }
//	cout << "\nReturning penalty: "<< Pen; exit(-1);
	return Pen;
}
double CBranchWARSProcess::PenaltyTraverse(int iNoTo, int iNoFr) {
	int i,j,Bra_i,Bra_j;
	double Penalty = 0.0;
//	cout << "\nVisiting branchnode["<<iNoTo<<"]";
	// Do post order tree traversal
	FOR(i,m_pTree->NoLinks(iNoTo)) {
		if(m_pTree->NodeLink(iNoTo,i) == iNoFr || m_pTree->NodeLink(iNoTo,i) == -1) { continue; }		// Skip the node it came from
		Penalty += PenaltyTraverse(m_pTree->NodeLink(iNoTo,i),iNoTo);	// Traverse the tree
	}
	if(m_pTree->NodeType(iNoTo) == branch) {
//		cout << "\nCalculating branchnode["<<iNoTo<<"]";
		FOR(i,m_pTree->NoLinks(iNoTo)) {
			for(j=i+1;j<m_pTree->NoLinks(iNoTo);j++) {
				if(m_pTree->BranchLabels()) { Bra_i = m_pTree->Labels()[m_pTree->NodeBra(iNoTo,i)]; Bra_j = m_pTree->Labels()[m_pTree->NodeBra(iNoTo,j)]; }
				else { Bra_i = m_pTree->NodeBra(iNoTo,i); Bra_j = m_pTree->NodeBra(iNoTo,j); }
//				cout << "\n\t["<<i<<"] " << m_pTree->NodeBra(iNoTo,i) << " : ["<<j<<"] " << m_pTree->NodeBra(iNoTo,j) << flush;
//				cout << " == " <<1 - exp(fabs(m_BraSigAlfa[m_pTree->NodeBra(iNoTo,i)]->Val() - m_BraSigAlfa[m_pTree->NodeBra(iNoTo,j)]->Val())) << flush;
				if(!m_BraSigAlfa.empty()) { Penalty += 1 - exp(fabs(m_BraSigAlfa[Bra_i]->Val() - m_BraSigAlfa[Bra_j]->Val())); }
				if(!m_BraSigInv.empty())  { Penalty += 1 - exp(fabs(m_BraSigInv[Bra_i]->Val() - m_BraSigInv[Bra_j]->Val())); }
				if(!m_BrapInv.empty())    { Penalty += 1 - exp(fabs(m_BrapInv[Bra_i]->Val() - m_BrapInv[Bra_j]->Val())); }
//				cout << " *" << flush;
		}	}
	}
	return Penalty;
}
/* *************************** Basic DNA covarion process ************************************ */
// This is a basic constructor function that all Markov modulated Markov models use
CFullDNATHMMProcess::CFullDNATHMMProcess(CData *D, CTree *T,string Name,int NumCovStates, vector <double> OriProbs, ETHMM_EQM_TYPE eqm_type, bool OptFreq, EHIDDEN_TYPE DoHidden, ERateTypes SeperateRates) : CTHMMProcess(D,T,Name,NumCovStates,OriProbs,eqm_type,OptFreq,DoHidden,SeperateRates) {
	// Check input functions
	assert(D != NULL && T != NULL);
	if(D->m_DataType != DNA) { Error("\nTrying to initialise DNA covarion process with data that doesn't look like nucleotides...\n\n"); }
}
///////////////////////////////////////////////////////////////////////////////////
// Destructor function
CFullDNATHMMProcess::~CFullDNATHMMProcess()	{
}

/////////////////////////////////////////////////////////////////////////////
// Functions for adding simple substitution parameters to the model
// Make the substitution model a normal REV
void CFullDNATHMMProcess::AddREV(int HiddenProc)	{
	int i,j;
	FOR(i,m_iDataChar)	{
		for(j=i+1;j<m_iDataChar;j++)	{
			AddDNASubPar(i,j,HiddenProc);
}	}	}
// Add a normal substitution parameter
CQPar * CFullDNATHMMProcess::AddDNASubPar(int S1, int S2, int HiddenProc)	{
	int i;
	double Val = 1.0;
	CQPar * Par = NULL;
	string Name, temp;
	// Ensure that the order is correct
	if(S1 > S2) { i = S1; S1 = S2; S2 = i; }
	assert(InRange(S1,0,m_iDataChar) && InRange(S2,0,m_iDataChar) && InRange(HiddenProc,-1,m_iHiddenChar));
	// Add the parameter
	Name = "Substitution_parameter(" + m_sABET.substr((S1 * m_iABET_length),m_iABET_length) + "<->" + m_sABET.substr((S2 * m_iABET_length),m_iABET_length) + ")";
	if(HiddenProc != -1) { Name += "(HiddenState=" + int_to_string(HiddenProc) + ")"; } else { Name += "(all)"; }
	if(m_pData->m_DataType == DNA && m_iDataChar == 4)	{
		if(	(m_sABET.substr((S1 * m_iABET_length),m_iABET_length) == "A" && m_sABET.substr((S2 * m_iABET_length),m_iABET_length) == "G")||
			(m_sABET.substr((S1 * m_iABET_length),m_iABET_length) == "G" && m_sABET.substr((S2 * m_iABET_length),m_iABET_length) == "A") ||
			(m_sABET.substr((S1 * m_iABET_length),m_iABET_length) == "C" && m_sABET.substr((S2 * m_iABET_length),m_iABET_length) == "T") ||
			(m_sABET.substr((S1 * m_iABET_length),m_iABET_length) == "T" && m_sABET.substr((S2 * m_iABET_length),m_iABET_length) == "C"))	{
				if(ALLOW_KAPPA_GUESS) { Val = m_pData->GuessKappa(); } else { Val = INITIAL_KAPPA; }
	}	}
	Par = new CQPar(Name,m_iChar,Val);
	if(HiddenProc == -1)	{
		FOR(i,m_iHiddenChar) {
			Par->AddQij((i*m_iDataChar) + S1,(i*m_iDataChar) + S2);
	}	} else {
		Par->AddQij((HiddenProc * m_iDataChar) + S1,(HiddenProc * m_iDataChar) + S2);
	}
	m_vpPar.push_back(Par);
	return Par;
}
// Add a transition/transversion parameter
CQPar *CFullDNATHMMProcess::AddKappa(int HiddenProc)						{
	int i,j,k;
	double Val = 1.0;
	CQPar * Par = NULL;
	string Name, temp;
	assert(InRange(HiddenProc,-1,m_iHiddenChar));
	// Add the parameter
	Name = "Kappa({AG}<->{CT})";
	if(HiddenProc != -1) { Name += "(HiddenState=" + int_to_string(HiddenProc) + ")"; } else { Name += "(all)"; }
	if(ALLOW_KAPPA_GUESS) { Val = m_pData->GuessKappa(); } else { Val = INITIAL_KAPPA; }
	Par = new CQPar(Name,m_iChar,Val);
	// Add the Kappa
	if(HiddenProc == -1)	{
		FOR(k,m_iHiddenChar)	{
			for(i=(k*m_iDataChar);i<(k+1)*m_iDataChar;i++)	{
				for(j=i+1 ;j<((k+1)*m_iDataChar);j++)	{
					if(	(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "A" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "G")||
						(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "G" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "A")||
						(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "C" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "T")||
						(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "T" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "C"))	{
						Par->AddQij(i,j);
	}	}	}	}	}else {
		for(i=(HiddenProc*m_iDataChar);i<(HiddenProc+1)*m_iDataChar;i++)	{
			for(j=i+1 ;j<((HiddenProc+1)*m_iDataChar);j++)	{
				if(	(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "A" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "G")||
					(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "G" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "A")||
					(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "C" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "T")||
					(m_sABET.substr((i * m_iABET_length),m_iABET_length) == "T" && m_sABET.substr((j * m_iABET_length),m_iABET_length) == "C"))	{
					Par->AddQij(i,j);
	}	}	}	}
	m_vpPar.push_back(Par);
	return Par;
}
// Add a transition/transversion ratio parameter for every process
void CFullDNATHMMProcess::AddSeperateKappas()	{
	int i;
#if HMP_MODEL_2_SIMPLE == 1
	AddKappa(); return;
#endif
	FOR(i,m_iHiddenChar)	{ AddKappa(i); }
}
/////////////////////////////////////////////////////////////////////////////////
// Output functions
ostream &CFullDNATHMMProcess::Output(ostream &os)	{
	int i;
	CBaseProcess::Output(os);
	os << "\n\n\t\t === Temporal Hidden Markov Model Consists of " << m_iHiddenChar << " hidden states ===";
	os << "\n\t\tGeneral statistics:\n\t\t\tOverall substitution rate: " << OverallSubRate() << "\n\t\t\tOverall transition rate:   " << OverallTransRate();
	os << "\n\t\t\tRatio trans/subs:          " << OverallTransRate()/OverallSubRate();

	assert((int) m_vpSubProcs.size() == m_iHiddenChar);
	FOR(i,m_iHiddenChar)	{
		os << " \n\t\t---\n\t\t\tHidden state " << i + 1 << ": ";
		os << "\n\t\t\t\tProbability of process:      " << m_vpCovProbs[i]->Val();
		os << "\n\t\t\t\tOverall substitution rate:   " << m_vpSubProcs[i]->Rate() << " cf. " << m_vpSubProcs[i]->CalcRate();
		os << "\n\t\t\t\tSubprocess equilibrium:      " << m_vpSubProcs[i]->Eqm();
		os << "\n\t\t\t\tRates between hidden states: in= " << m_vpQMat[0]->TransRateTo(i) << "; out= " << m_vpQMat[0]->TransRateFrom(i);
		m_vpSubProcs[i]->Output(os,'\t',"\t\t\t\t");
	}
	return os;
}

// A function that shuffles the parameters around in the sub processes in the THMM
void CFullDNATHMMProcess::ShuffleParameters()	{
	int i,j,k;
	vector <double> Par;
	vector <double> New;
	double Val;
	// Do rates
	FOR(i,m_iHiddenChar)	{ if(SubProc(i)->SpecialRate()) { continue; } Par.push_back(SubProc(i)->Rate()); }
//	cout << "\nOriginal rates: " << Par;
	random_shuffle(Par.begin(),Par.end());
//	cout << "\nNew rates:      " << Par;
	j= 0; FOR(i,m_iHiddenChar)	{ if(SubProc(i)->SpecialRate()) { continue; } SubProc(i)->Rate(Par[j++]); }
	Par.clear();

	// Do state probabilities
	FOR(i,m_iHiddenChar)	{ Par.push_back(m_vpCovProbs[i]->Val()); }
	Val = 0.0; FOR(i,m_iHiddenChar) { Val += fabs(Val - 0.25); }
	if(Val > 0.5) { FOR(i,m_iHiddenChar) { Par[i] = 0.5 + (0.25 * Random()); } }
//	cout << "\nOriginal probs: " << Par << " == " << Sum(Par);
	random_shuffle(Par.begin(),Par.end());
//	cout << "\nNew probs:      " << Par << " == " << Sum(Par);
	j = -1; FOR(i,m_iHiddenChar)	{
		if(m_vpCovProbs[i]->Special()) { j = i; continue; }
		else						{ m_vpCovProbs[i]->SetVal(Par[i],true,true,false); }
	}
	assert (j != -1);
	m_vpCovProbs[j]->SetVal(Par[j],true,true,true);
	Par.clear();
	// Do the eqm distributions
	FOR(i,(int)m_vpEqm.size()) {
		Par = m_vpEqm[i]->Eqm(); New = m_pData->GetFreq(-1);
		Val = 0; FOR(j,m_iChar) { if(j % m_iDataChar == 0) { k = 0; } Val += fabs(Par[j] - New[k++]); }
		if(Val < 0.5 * m_iHiddenChar) { m_vpEqm[i]->Shuffle(); }	// If values aren't that different from original, then shuffle
		else {	// Otherwise reset to the data values;
			m_vpEqm[i]->ResetEqm(New,true);
		}
	}
	// Reset all the hidden changes to 0.1
	FOR(i,(int)m_vpPar.size())	{
		if(m_vpPar[i]->Opt() && m_vpPar[i]->Name().find("Hidden_change") != string::npos) { m_vpPar[i]->SetVal(0.1,true,true); }
	}

}

/* ***************************************************************************
		Class implementation for substitution processes in THMMs
   ***************************************************************************/
CTHMMSubProc::CTHMMSubProc(int ProcNum,int DataChar, double *QMat, int QChar, CPar *Rate, CBaseEqm *Eqm)	{
	int i,j;
	// Set some basic values
	assert(QMat != NULL && Rate != NULL && Eqm != NULL);
	m_iHiddenState = ProcNum; m_pRatePar = Rate; m_pEqm = Eqm;
	m_iChar = DataChar;
	// Get the rate matrix
	for(i=ProcNum * m_iChar;i<(ProcNum+1) * m_iChar;i++)	{
		for(j=ProcNum * m_iChar;j<(ProcNum+1) * m_iChar;j++)	{ m_vpQ.push_back(&QMat[(i*QChar)+j]); }
	}
	assert((int)m_vpQ.size() == m_iChar * m_iChar);
}
CTHMMSubProc::~CTHMMSubProc()	{
	int i;
	m_iHiddenState = -1;
	m_pRatePar = NULL;
	m_pEqm = NULL;
	FOR(i,(int)m_vpQ.size()) { m_vpQ[i] = NULL; } m_vpQ.~vector();
}
// Calculate the overall rate of the substitution process
double CTHMMSubProc::CalcRate()	{
	int i,j;
	double Total = 0;
	vector <double> eqm = m_pEqm->SubEqm(m_iHiddenState);
	// Check entry conditions
	assert((int) m_vpQ.size() == m_iChar * m_iChar && m_pRatePar != NULL && m_pEqm != NULL);
	// Calculate the rate
	FOR(i,m_iChar) {
		FOR(j,m_iChar)	{
			if(i==j) { continue; }
			Total += eqm[i] * *m_vpQ[(i*m_iChar)+j];
	}	}
	return Total;
}
// Scale the process so its mean rate is defined by m_pRatePar
void CTHMMSubProc::Scale(double SetRate)	{
	int i;
	double RateAdjuster, CurRate;
	// Check entry conditions
	assert((int) m_vpQ.size() == m_iChar * m_iChar && m_pRatePar != NULL && m_pEqm != NULL);
	// Pick up default entry
	if(SetRate < 0) { SetRate = m_pRatePar->Val(); }

//	cout << "\nScaling process: " << SetRate << " cf. " << m_pRatePar->Val();

	// Initialise rate if required
	if(SetRate > FLT_EPSILON)	{
		// Do the scaling
		CurRate = CalcRate();
		if(CurRate > FLT_EPSILON)	{
			RateAdjuster = SetRate/CurRate; FOR(i,m_iChar*m_iChar) { *m_vpQ[i] *= RateAdjuster; }
		} else if(DEBUG_HMP_MODEL && diff(CurRate,0) == 0)	{
			cout << "\nCurrent process:";
			CTHMMSubProc::Output();
			cout << "\nWarning: trying to rescale a CTHMMSubProc to rate " << m_pRatePar->Val() << " to 0\n"; exit(-1);
	}	} else {
		FOR(i,m_iChar*m_iChar) { *m_vpQ[i] = 0.0; }
	}
}
// Output the process's Q matrix
void CTHMMSubProc::Output(ostream &os, char delim,string lead)	{
	int i,j;
	double Total;
	os << endl << lead << "Substitution process["<< m_iHiddenState << "]:";
	FOR(i,m_iChar)	{
		os << endl << lead << "\t";
		// Get the diagonal
		Total = 0.0;
		FOR(j,m_iChar) { if(i==j) { continue; } Total -= *m_vpQ[(i*m_iChar)+j]; }
		FOR(j,m_iChar) { if(i==j) { os << Total << delim; continue; } os << *m_vpQ[(i*m_iChar)+j] << delim; }
	}
}

/* ***************************************************************************
	Classes describing specific covarion process
   *************************************************************************** */

C2StateCovProc::C2StateCovProc(CData *D, CTree *T, string Name, vector <double> Probs, ETHMM_EQM_TYPE  eqm_type) : CFullDNATHMMProcess(D,T,Name,2,Probs,eqm_type,true,H_same,same)	{
	EnforceFixedState();
	SetSubProcRate(0,1.0); SetSubProcRate(1,0.0);
	NoOptimiseRates();
	OptimiseHiddenProbs();
}

/* ****************************************************************************
 * 		Coevolution models
 * **************************************************************************** */

CPairwiseCoevoProcess::CPairwiseCoevoProcess(CData *D, vector <double> Left, vector <double> Right, CTree *T, vector <double> *R, double init_psi) : CBaseProcess(D,T) {
	int i;
	double Total;
	cout << "\nCPairwiseCoevoProcess::CPairwiseCoevoProcess(...)" << flush;
	// Do some initialisation and checking
	MakeBasicSpace(D->m_DataType); AddSimpleEqm();
	m_pPsi = NULL; OneCharEqmLeft = NULL; OneCharEqmRight = NULL;
	// Create appropriate values
	m_iChar = D->m_iChar;
	switch(D->m_DataType) {
	case RY2: m_iOneChar = 2; break;
	case DNA2: m_iOneChar = 4; break;
	case AA2: m_iOneChar = 20; break;
	default:
		Error("\nPassing an incorrect data format to CPairwiseCoevoProcess constructor...\n");
	}
	m_vRMat.assign(m_iOneChar*m_iOneChar,0);
	if(R!=NULL) {
		assert(R->size() == m_vRMat.size());
		m_vRMat = *R;
	}
	m_pPsi = new CQPar("Psi",m_iChar, init_psi, true,-(MAX_PAR_VALUE),MAX_PAR_VALUE);
	m_vpPar.push_back(m_pPsi);
	// Set up equilbrium distributions
	OneCharEqmLeft = new CSimpleEqm(m_iOneChar,&m_vpPar,Left);
	OneCharEqmRight = new CSimpleEqm(m_iOneChar,&m_vpPar,Right);
}

CPairwiseCoevoProcess::~CPairwiseCoevoProcess() {
	if(m_pPsi != NULL) { delete m_pPsi; }
	delete OneCharEqmLeft; delete OneCharEqmRight;
}

// Create the new equilibrium distribution from the propensity matrix, the psi parameter, and the Left and Right One Character Eqms.
// The formula for creating the the joint probability is pi_{ik} = pi_{i} * pi_{k} * e^(psi * r_{ik}) * c
vector <double> CPairwiseCoevoProcess::MakeCoevoEqm(double Psi) {
	int i,j;
	vector <double> CoevoEqm(m_iChar,0), Left, Right;
	double Total = 0.0;
	// Check the requisite parameters and values are in the model
	assert(m_vRMat.size() == m_iChar && m_pPsi != NULL);
	assert(OneCharEqmLeft != NULL && OneCharEqmRight != NULL);
	// Create the new eqm
	Left = OneCharEqmLeft->Eqm(); Right = OneCharEqmRight->Eqm();
//	cout << "\nSums\nLeft:  " << Left << "\nRight: " << Right;
	if(fabs(Psi + BIG_NUMBER) < 0.01) { Psi = m_pPsi->Val(); }
	FOR(i,m_iOneChar) {
		FOR(j,m_iOneChar) {
			Total += CoevoEqm[(i*m_iOneChar)+j] = Left[i] * Right[j] * exp(Psi * m_vRMat[(i*m_iOneChar)+j]);
	}	}
	// Scale
//	cout << "\nEqm:   " << CoevoEqm;
//	cout << "\nSum = " << Total << " and exp(Psi * R[0][1]) == exp(" << Psi << " * " << m_vRMat[1] << "): " << exp(Psi * m_vRMat[(i*m_iOneChar)+j]);
	FOR(i,m_iChar) { CoevoEqm[i] /= Total; }
//	cout << "\nNeed to sort rate in CPairwiseCoevoProcess::MakeCoevoEqm(...)";
	// Set the rate; Note that it's twice the scaling value because it covers two columns in the alignment
//	Rate(Total * 2);
//	cout << "\nRate set to: " << Rate();
	return CoevoEqm;

}

// Function for creating the stationary distribution and the rate of the process
vector <double> CPairwiseCoevoProcess::RootEqm()	{
	return MakeCoevoEqm(m_pPsi->Val());
}

// Functions controlling the Q matrices
bool CPairwiseCoevoProcess::PrepareQMats(vector <int> Qs2do, bool DoScale) {
		// Make the eqm distribution
	// Get the eqm and scaling factor
//	cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Resetting m_vpEqm[0]";
//	cout << "\nBefore -> m_vpEqm[0]->Eqm: " << m_vpEqm[0]->Eqm();
//	cout << "\nRootEqm() returns this:    " << RootEqm();
	m_vpEqm[0]->ResetEqm(RootEqm(),false);
//	cout << "\nAfter -> m_vpEqm[0]->Eqm:  " << m_vpEqm[0]->Eqm();
	// No build the matrix
	return CBaseProcess::PrepareQMats(Qs2do,DoScale);
}

/* Basic amino acid Coevolution model
 * ---
 */

CAACoevoProcess::CAACoevoProcess(CData *D, vector <double> Left, vector <double> Right, CTree *T, vector<double> *R, double init_psi,bool DoF) : CPairwiseCoevoProcess(D,Left,Right,T,R,init_psi) {
	int i, j, k, l, Counter = 0;
	string sOneABET = AA_ABET, ABET = D->m_sABET, Name;
	CQPar *ZeroPar = NULL, *TempPar = NULL;
	vector <CQPar *> MatPar;
	double WAGVal[] = {0.551571,0.509848,0.635346,0.738998,0.147304,5.429420,1.027040,0.528191,0.265256,0.0302949,0.908598,3.035500,1.543640,0.616783,0.0988179,1.582850,0.439157,0.947198,6.174160,0.021352,5.469470,1.416720,0.584665,1.125560,0.865584,0.306674,0.330052,0.567717,0.316954,2.137150,3.956290,0.930676,0.248972,4.294110,0.570025,0.249410,0.193335,0.186979,0.554236,0.039437,0.170135,0.113917,0.127395,0.0304501,0.138190,0.397915,0.497671,0.131528,0.0848047,0.384287,0.869489,0.154263,0.0613037,0.499462,3.170970,0.906265,5.351420,3.012010,0.479855,0.0740339,3.894900,2.584430,0.373558,0.890432,0.323832,0.257555,0.893496,0.683162,0.198221,0.103754,0.390482,1.545260,0.315124,0.174100,0.404141,4.257460,4.854020,0.934276,0.210494,0.102711,0.0961621,0.0467304,0.398020,0.0999208,0.0811339,0.049931,0.679371,1.059470,2.115170,0.088836,1.190630,1.438550,0.679489,0.195081,0.423984,0.109404,0.933372,0.682355,0.243570,0.696198,0.0999288,0.415844,0.556896,0.171329,0.161444,3.370790,1.224190,3.974230,1.071760,1.407660,1.028870,0.704939,1.341820,0.740169,0.319440,0.344739,0.967130,0.493905,0.545931,1.613280,2.121110,0.554413,2.030060,0.374866,0.512984,0.857928,0.822765,0.225833,0.473307,1.458160,0.326622,1.386980,1.516120,0.171903,0.795384,4.378020,0.113133,1.163920,0.0719167,0.129767,0.717070,0.215737,0.156557,0.336983,0.262569,0.212483,0.665309,0.137505,0.515706,1.529640,0.139405,0.523742,0.110864,0.240735,0.381533,1.086000,0.325711,0.543833,0.227710,0.196303,0.103604,3.873440,0.420170,0.398618,0.133264,0.428437,6.454280,0.216046,0.786993,0.291148,2.485390,2.006010,0.251849,0.196246,0.152335,1.002140,0.301281,0.588731,0.187247,0.118358,7.821300,1.800340,0.305434,2.058450,0.649892,0.314887,0.232739,1.388230,0.365369,0.314730};


	cout << "\nCAACoevoProcess constructor...";
	cout << "\nm_iChar: " << m_iChar;
	cout << "\nsOneABET: " << sOneABET;
	cout << "\nABET:     " << ABET; int count = 0;
	cout << "\nNeed to change so it uses wag...";
	// Create Q Matrix
	Add_QMat("Coevo",AA2);
	// Create WAG-based exchangeability parameters (including zero)
	ZeroPar = new CQPar("DoubleZeros",m_iChar,0,false,0,0,REPLACE);
	FOR(i,m_iOneChar) {
		for(j=i+1;j<m_iOneChar;j++) {
			// Initialise
			Name = ""; Name += sOneABET[i]; Name += "<->"; Name += sOneABET[j];
//			TempPar = new CQPar(Name,m_iChar,WAGVal[Counter++],false);
			TempPar = new CQPar(Name,m_iChar,1.0,false);	// Bit that stops WAG
//			cout << "\nPreparing parameter " << Name; count = 0;
			// Get the appropriate values
			FOR(k,m_iChar) {
				for(l=k+1;l<m_iChar;l++) {
	//				string GetPos(string Seq, int Pos, int PosLength);
						// Add the zero where both characters differ
						if(ABET[k*2] != ABET[l*2] && ABET[(k*2)+1] != ABET[(l*2)+1]) {
							// cout << "\n[" << k<< "," << l << "]" << GetPos(ABET,k,2) << "," << GetPos(ABET,l,2); cout << "\n\tZero : [" << k << "," << l <<"] = [" << ABET[(k*2)] << ABET[(k*2)+1] << "," << ABET[(l*2)] << ABET[(l*2)+1] << "]";
							ZeroPar->AddQij(k,l,true);
							continue;
						}
						// Check the first character
						if(ABET[(k*2)] == sOneABET[i] && ABET[(l*2)] == sOneABET[j]) {
							// cout << "\n[" << k<< "," << l << "]" << GetPos(ABET,k,2) << "," << GetPos(ABET,l,2); cout << "\n\tMatch first " << count++ <<": [" << k << "," << l <<"] = [" << ABET[(k*2)] << ABET[(k*2)+1] << "," << ABET[(l*2)] << ABET[(l*2)+1] << "]"; cout << " --> " << ABET[(k*2)] << " -> " << ABET[(l*2)];
							TempPar->AddQij(k,l,true);
							continue;
						}
						// Check the second character
						if(ABET[(k*2)+1] == sOneABET[i] && ABET[(l*2)+1] == sOneABET[j]) {
							// cout << "\n[" << k<< "," << l << "]" << GetPos(ABET,k,2) << "," << GetPos(ABET,l,2); cout << "\n\tMatch second " << count++ <<": [" << k << "," << l <<"] = [" << ABET[(k*2)] << ABET[(k*2)+1] << "," << ABET[(l*2)] << ABET[(l*2)+1] << "]"; cout << " --> " << ABET[(k*2)+1] << " -> " << ABET[(l*2)+1];
							TempPar->AddQij(k,l,true);
							continue;
			}	}	}
			// Store the parameter
//			cout << "\nValue: " << *TempPar;
			m_vpPar.push_back(TempPar); TempPar = NULL;
	}	}
	m_vpPar.push_back(ZeroPar); ZeroPar = NULL;
	// Copy the propensity matrix to m_vRMat
	m_vRMat = *R;
	// Set up the OneChar equilibrium distributions
	// Exit
	cout << "\nDone with CAACoevoProcess\n\n";
//	exit(-1);

}

CAACoevoProcess::~CAACoevoProcess() {

}


