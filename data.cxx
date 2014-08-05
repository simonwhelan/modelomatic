/////////////////////////////////////////////////
// Data.cxx
////////////////////////////////////////////////

#include "data.h"  // header for class
#include "tree.h"
// Debugging info
#define DEBUG_DATA 0

#if DO_MEMORY_CHECK
extern CMemChecker memory_check;
#endif


//////////////////////////////////////////////
// Constructors for data object
/////////////////////////////////////////////////////////////

CData::CData(int NoSeq, int Size, EDataType Type, vector <string> *Names)	{
#if DO_MEMORY_CHECK
	memory_check.CountCData++;
#endif
	int i;vector <int> viTemp;
	// Some original preperation
    m_bValid = false; m_iGenCode = -1;
	// Do some basic initialisation
	m_iNoSeq = NoSeq; m_iSize = m_iTrueSize = Size; m_DataType = Type;
	m_iChar = NumStates(Type);
	m_sABET = DataStates(Type);
	if(Names == NULL)	{ FOR(i,NoSeq) { m_vsName.push_back("Sequence["+int_to_string(i+1)+"]");  } }
	else				{
		if((int)Names->size() != m_iNoSeq) { Error("\nTrying to create data CData::CData(...) with non-matching tree and names...\n\n"); }
		FOR(i,NoSeq) { m_vsName.push_back(Names->at(i)); }
	}
	FOR(i,m_iChar) { m_vFreq.push_back((double) ( 1.0 / (double) m_iChar)); }
}

CData::CData(string file, EDataType SpecType, bool AllowFail, streampos FilePos) {
#if DO_MEMORY_CHECK
	memory_check.CountCData++;
#endif
    int i,j, CumSize,LastSize;
	vector <string> vName, vData,Toks;
	vector <int> ariPos;
	vector <int> SiteLabels;
	string store;
	bool First = true, DoSiteLabels = false;
	EDataType  Type,Old;

	// Some original preperation
    m_bValid = false; m_iGenCode = -1;

	// Check whether a NEXUS file first
    ifstream input(file.c_str());
   	getline(input,store);
   	input.close();
   	if(input.eof()) { cout << "\nUnexpected end of file at beginning of file......"; if(!AllowFail) { Error(); } }
   	Toks = Tokenise(store);
	if(Toks[0].find("#NEXUS") != string::npos) {
		InputNexus(file);
		return;
	}
    ReadData(file,vName,vData);
	/////////////////////////////////////////////////////////////////////////
	// Do a bit more checking and initialise the data object
    ////////////////////////////////////////////////
	// Guess the type of data
    // New code for guessing data type
    int DNA_count = 0, AA_count = 0, Error_count = 0;
    if(SpecType == NONE) {
    	FOR(i,(int)vData.size()) 	{
    		Old = GuessDataType(vData[i]);
    		if(Old == DNA) { DNA_count++; } else if(Old == AA) { AA_count++; } else { Error_count++; }
    	}
    	if(Error_count > (int)vData.size() *0.75) {
    		cout << "\nError when guessing data type: " << Error_count << "/" << vData.size() << " sequences could not be identified.\nPlease inspect your data for excesses of weird or gap characters.\nThis message will be also triggerred when there are a lot of gaps (>" << double_to_string(100*(1.0-MIN_DATA_PERCENT)) << "%) or frequent non-standard characters.\n\n"; Error(); }
    	if(DNA_count < (int)vData.size() /2 && AA_count < (int)vData.size()) {
    		cout << "\nError when guessing data type. Of " << (int)vData.size() << " sequences: " << DNA_count << " appear as DNA and " << AA_count << " appear as AA. \nPlease inspect your data for excesses of weird or gap characters.\nThis message will be also triggerred when there are a lot of gaps (>" << double_to_string(100*(1.0-MIN_DATA_PERCENT)) << "%) or frequent non-standard characters.\n\n"; Error();}
    	if(DNA_count > AA_count && AA_count > 3) { cout << "\nWARNING: Identified DNA sequences, but " << AA_count << "/" << vData.size() << " look like amino acid sequences. Proceeding with analysis, but be cautious...\n"; }
    	if(AA_count > DNA_count && DNA_count > 3) { cout << "\nWARNING: Identified DNA sequences, but " << DNA_count << "/" << vData.size() << " look like DNA sequences. Proceeding with analysis, but be cautious...\n"; }
    	if(Error_count > (int)vData.size() *0.1) { cout << "\nWARNING: " << Error_count << "/" << vData.size() << " sequences could not be reliably identified.\nInspect your data for excesses of weird or gap characters."; }

    	if(DNA_count > AA_count) { Type = DNA; } else { Type = AA; }
    }

    /* OLD CODE FOR GUESSING DATA TYPE
    int ErrorCount = 0;
	if(SpecType == NONE)	{
		// Find a stable data type
		FOR(i,(int)vData.size())	{
			Old = GuessDataType(vData[i]);
			if(Old != NONE) { break; }
		}
		if(i == (int)vData.size())  { Error("Could not guess data type... Seems to be very little data here... All sequences have >"+double_to_string(100*(1.0-MIN_DATA_PERCENT))+"% gaps...\n\n"); }
		FOR(i,(int)vData.size())	{
			Type = GuessDataType(vData[i]);
			if(Type == NONE) { ErrorCount++; Type = Old; }
			if(Type != Old) { Type = NONE; break; }
	}	} else { Type = SpecType; }
	// Ensure the data is of a real type
	if(ErrorCount > (int)vData.size() / 3) { cout << "\nError: There are " << ErrorCount << " sequences with little or no data..."; exit(-1); }
	*/
	if(Type == NONE) { Error("\nCouldn't guess data type. Please inspect your data for excesses of weird or gap characters.\nThis message will be also triggerred when there are a lot of gaps (>" + double_to_string(100*(1.0-MIN_DATA_PERCENT)) + "%) or frequent non-standard characters.\n\n" ); }
    InputData(Type,vData,vName,SiteLabels,AllowFail);		// Put the sequence data into the object
}

// Inputs data from a set of arrays (for pairwise distances)
CData::CData(int NoSeq,int Size, vector <string> InSeq, vector <string> InName,EDataType SpecType) {
#if DO_MEMORY_CHECK
	memory_check.CountCData++;
#endif
	int i;
	bool flag = false;
	vector <int> SiteLabels;
	EDataType Type, Old;
	m_iNoSeq = NoSeq;
	m_iSize = m_iTrueSize = Size;
	// Some original preperation
    m_bValid = false; m_iGenCode = -1;
	if((int)InSeq.size() != m_iNoSeq) { cout << "\nCData() constructor only passed " << InSeq.size() << " sequences when expecting " << m_iNoSeq; exit(-1); }
	FOR(i,m_iNoSeq) {
		if((int)InSeq[i].size() != Size) { cout << "\nCData() constructor: seq["<<i<<"] of length " << InSeq[i].size() << " when expecting " << Size; }
	}
	// Guess the type of data
	if(SpecType == NONE)	{
		Old = GuessDataType(InSeq[0]);
		FOR(i,m_iNoSeq)	{
			Type = GuessDataType(InSeq[i]);
			if(Type == NONE || Type != Old) { Type = NONE; break; }
	}	} else { Type = SpecType; }
    InputData(Type,InSeq,InName,SiteLabels);		// Put the sequence data into the object
}

// Copy constructor
CData::CData(const CData & CopD)	{
#if DO_MEMORY_CHECK
	memory_check.CountCData++;
#endif
	m_iNoSeq = CopD.m_iNoSeq;
	m_iSize = CopD.m_iSize;
	m_DataType = CopD.m_DataType;
	m_sABET = CopD.m_sABET;
	m_iChar = CopD.m_iChar;
	m_ariSeq = CopD.m_ariSeq;
	m_viNoChange = CopD.m_viNoChange;
	m_ariPatOcc = CopD.m_ariPatOcc;
	m_iName = CopD.m_iName;
	m_vsName = CopD.m_vsName;
	m_vFreq = CopD.m_vFreq;
	m_viSiteLabels = CopD.m_viSiteLabels;
	m_vviCoevoMapping = CopD.m_vviCoevoMapping;
	m_viRealPatOcc = CopD.m_viRealPatOcc;
	m_ariPatMap = CopD.m_ariPatMap;
	m_iTrueSize= CopD.m_iTrueSize;
	m_vsTrueSeq = CopD.m_vsTrueSeq;
	m_bValid = CopD.m_bValid;
	m_iGenCode = CopD.m_iGenCode;
}

CData::~CData() {
#if DO_MEMORY_CHECK
	memory_check.CountCData--;
#endif
	Clean();
}

void CData::Clean() {
	int i;
	CleanBootstrap();
	if(!m_ariSeq.empty()) { FOR(i,m_iNoSeq) { m_ariSeq[i].clear(); } m_ariSeq.clear(); }
	m_ariPatOcc.clear(); m_ariPatMap.clear();

	m_vFreq.clear();
	m_iChar = m_iName = m_iNoSeq = m_iSize = m_iTrueSize= - 1;
	m_sABET = "";
	m_viNoChange.clear();
	m_vsName.clear();
	m_vFreq.clear();
	m_viSiteLabels.clear();
	m_vsTrueSeq.clear();
	m_viRealPatOcc.clear();
	m_vviCoevoMapping.clear();
	m_bValid = false;
}

////////////////////////////////////////////////////////////////////////////////////////
// Function to read data from file
bool ReadData(string File, vector <string > &Names, vector <string > &Seqs, bool AllowFail)	{
    int i,j, CumSize,LastSize;
	vector <string > Toks;
	string store;
	int NoSeq = -1, Size = -1;
	vector <int> ariPos;
	ifstream input(File.c_str());

	// Input checking
	if (!Names.empty()) { Error("\nTrying to fill already assigned Names vector\n\n"); }
	if (!Seqs.empty()) { Error("\nTrying to fill already assigned Seqs vector\n\n"); }
    // All file types start with #sequences #length

	getline(input,store);
	if(input.eof()) { cout << "\nUnexpected end of file at beginning of file..."; if(AllowFail) { return false; } else { Error(); } }
	Toks = Tokenise(store);
	if(Toks.size() < 2) { cout << "\nFirst line of datafile should contain: #seq	#length"; if(AllowFail) { return false; } Error(); }
	NoSeq = atoi(Toks[0].c_str()); Size = atoi(Toks[1].c_str());
	if(NoSeq <= 0 || Size <=0) { return false; }
    // Now initialise memory
	Names.assign(NoSeq,"");
	Seqs.assign(NoSeq,"");
	ariPos.assign(NoSeq,0);

	// If the data is partitioned then input the partitioning information
	if(Toks.size() == 3) {
		cout << "\nData partitioned: entering partition information, this is currently not available...\n\n"; exit(-1);
/*		int NoPart = atoi(Toks[2].c_str());
		Toks = Tokenise(GetDataLine(&input));
		if((int)Toks.size() != Size) { cout << "\nNumber of sites with defined data partitions = " << Toks.size() << "; expecting " << Size; Error("\nFailed\n\n"); }
		FOR(i,(int)Toks.size()) {
			if(atoi(Toks[i].c_str()) >= NoPart) { cout << "\nMaximum number of partitions defined as " << NoPart << "; just found one defined as " << Toks[i]; Error("\nFailed\n\n"); }
			m_viSiteLabels.push_back(atoi(Toks[i].c_str()));
		} */
	}
	// Get the data
	/////////////////////////////////////////////////////////////
	// Currently I can get the following types of data
	// 1.) Sequential	-	PAML or FASTA
	// 2.) Interleaved	-	Phylip

    // Move to the first line of sequence info and use it to guess the data type
	store = GetDataLine(&input);
	Toks = Tokenise(store);
	// PAML and FASTA data have no sequences on the first line
	if(Toks.size() == 1) {
		FOR(i,NoSeq)	{
			// Get the name
			Names[i] = GetName(Toks[0]);
			// Get the sequence
			CumSize = 0; while(CumSize < Size)	{
				store = GetDataLine(&input); store = EatWhiteSpace(store);
				if(input.eof()) { cout << "\nUnexpected end of file... Try adding a return at the end of the data file?\n"; if(AllowFail) { return false; } else { Error(); } }
				Seqs[i] += store; CumSize += (int)store.size();
			}
			if((int) Seqs[i].size() != Size) { cout << "\nSequence["<<i<<"] " << Names[i] << " of length " << Seqs[i].size() << ", expected " << Size << "\n" << Seqs[i] << "\n\n"; }
			if(i < NoSeq - 1) {
				store = GetDataLine(&input); Toks = Tokenise(store);
				if(Toks.size() > 1)  { cout << "Name line should not contain data in sequential or fasta format"; if(AllowFail) { return false; } else { Error(); } }
		}	}
	} else {
	// Assume that its interleaved -- Note that when there is extra white-space then the Tokens will need compressing
		CumSize = 0;
		while(CumSize < Size) {
			FOR(i,NoSeq) {
				if(Names[i].empty()) { Names[i] = Toks[0]; } // If required get the name
				// Condense the remaining Tokens
				if(Toks[0] == Names[i]) { j = 1; } else { j = 0; }
				store = ""; for(j;j<(int)Toks.size();j++) { store += Toks[j]; }
				// Get the data and check
				Seqs[i] += store;
				if(i == 0)	{ CumSize += LastSize = (int)store.size(); }
				else		{ if((int)store.size() != LastSize) { cout << "\nData for species["<<i<<"] " << Names[i] << " of unexpected length... [Obs: " << store.size() << " LastSize: " << LastSize << " cf. Exp: " << Size << ">\n" << store << "\nThis can be caused by strange line end characters between file header and sequences\n"; if(AllowFail) { return false; } Error(); } }
				if(CumSize < Size || i < NoSeq - 1) { Toks = Tokenise(GetDataLine(&input)); }
				if(input.eof()) { cout << "\nUnexpected end of file... Try adding a return at the end of the data file?\n"; if(AllowFail) { return false; } else { Error(); } }
			}
		}
	}
	input.close();
	return true;
}

////////////////////////////////////////////////////////////
// Function to count number of analysable characters
int CData::CountMSAChars()	{
	int i,j, count = 0;
	FOR(i,m_iNoSeq) {
		FOR(j,m_iSize) { if(m_ariSeq[i][j] != m_iChar) { count += m_ariPatOcc[j]; } }
	}
	return count;
}

// Get Name function
////////////////////////////////////////////////////////////
// Arguments are a char string
// returns an integer that signals where the name ends in
// that char string
////////////////////////////////////////////////////////////

char* GetName(char *string, char *name)	{
    int i,j=0;
	if(name != NULL) { DEL_MEM(name); }
	GET_MEM(name,char,(int)strlen(string) + 2);
    // Skip '>' on fasta data
    if(string[0] == '>') { j = 1; }
    char top[] = "cosn";
    for(i=0,j;j<(int)strlen(string)+1;j++)	{ if(!isspace(string[j])) { name[i++] = string[j]; } }
	return name;
}

string GetName(string name)	{
	if(name.empty()) { Error("Empty name in data"); }
	if(name[0] == '>') { name.erase(0,1); }
	return EatWhiteSpace(name);
}



///////////////////////////////////////////////////////////////////
// Bootstrapping routine
void CData::Bootstrap()	{
	int i,count=0;
	vector <double> Probs;
	CleanBootstrap();
	// Store the data
	FOR(i,m_iSize) {
		m_viRealPatOcc.push_back(m_ariPatOcc[i]);
		Probs.push_back((double) m_ariPatOcc[i]);
		count += m_ariPatOcc[i];
		m_ariPatOcc[i] = 0;
	}
	Probs = NormaliseVector(Probs,1.0);
	assert(count == m_iTrueSize);
	FOR(i,m_iTrueSize) { m_ariPatOcc[DiscreteRand(&Probs)]++; }
}

void CData::CleanBootstrap()	{
	int i;
	if(m_viRealPatOcc.empty()) { return; }
	assert((int)m_viRealPatOcc.size() == m_iSize);
	FOR(i,m_iSize) { m_ariPatOcc[i] = m_viRealPatOcc[i]; }
	m_viRealPatOcc.clear();
}
////////////////////////////////////////////////////////////
// Function to get the multinomial lnL of the data
double CData::MultilnL()	{
	int i;
	double ret = 0;
	if(!Valid()) { Error("\nCannot get multinomial lnL on non CData::Valid() data\n\n"); }
	FOR(i,m_iSize)	{ ret += ( (double)m_ariPatOcc[i] * log((double)m_ariPatOcc[i]));	}
	ret -= (m_iTrueSize * log((double) m_iTrueSize));
	return ret;
}

//////////////////////////////////////////////////////////////////////////
// Function that outputs the distribution of character values
//	Pattersns are dealt with as:
//		Value of each pattern = = 1 / 4^{Number of gaps}
//		Absolute pattern (no missing data) == 1
//
void CData::OutDataDist(string File, EDataType Type)	{
	int i,j,k,gaps,Count;
	unsigned int NoPat;
	vector <vector <int> > Pats;
	vector <int> iters;
	vector <double> Vals;
	string ABET = DataStates(Type);
	string PatCheck, temp;
	// Do some initialisation and error checking
	if(Type == AA && m_DataType != AA) { Error("\nTrying to output Amino acids from non-amino acid data\n\n"); }
	if(m_iNoSeq == 0) { Error("\nCannot output patterns without data... m_iNoSeq == 0 in CData::OutDataDist(...)\n"); }
	NoPat = (int) pow((double) NumStates(Type), (double) m_iNoSeq);
	cout << "\nCreating frequency distribution of data points\n\tm_iNoSeq: " << m_iNoSeq << "\tNoPats: " << NoPat;
	if(NoPat == 0 || NoPat > 1.5E+6) { Error("\nTrying to output too many datapatterns..."); }
	Vals.assign(NoPat,0.0);
	iters.assign(m_iNoSeq,0);
	// Create the data patterns
	cout << "\n\nGenerating patterns";
	CreatePatterns(0,&iters,Type,&Pats);
	cout << "\n\nScanning patterns from data" << flush; Count = 0;
	// Grab patterns from data
	FOR(i,m_iSize)	{
		if(Count ++ % 60 == 0) { cout << "\n"; } cout << "." << flush;

//		cout << "\nSite["<<i<<"]:";
		// Get the gaps
		gaps = 0;
		FOR(j,m_iNoSeq) { if(m_ariSeq[j][i] == m_iChar) { gaps++; } }
		// Do the pattern comparison
		FOR(j,(int) Pats.size()) {
			FOR(k,m_iNoSeq) {
				if(m_ariSeq[k][i] == m_iChar) { continue; }
				if(m_ariSeq[k][i] != Pats[j][k]) { break; }
			}
			if(k == m_iNoSeq) {
//				cout << "\n\tComparing: " << Pats[j] << " : "; FOR(k,m_iNoSeq) { cout << m_ariSeq[k][i]; }
				Vals[j] += m_ariPatOcc[i] / pow((double)m_iChar,(double)gaps);
//				cout << " == " << m_ariPatOcc[i] / pow((double)m_iChar,(double)gaps);
				if(gaps == 0) { break; }	// If no gaps then its a unique match
			}
		}
	}
	double Sums = 0;
	// Output the patterns
	ofstream out(File.c_str());
	FOR(i,(int)Pats.size()) { FOR(j,(int) Pats[i].size()) { out << ABET[Pats[i][j]]; } out << "\t"; } out << "Total\n";
	Sums = 0.0; FOR(i,(int)Vals.size()) { out << Vals[i] << "\t"; Sums += (int) Vals[i]; } out << Sums << "\n";
	out.close();
}

void CData::CreatePatterns(int Pos, vector <int> *iters,EDataType Type, vector <vector <int> > *Pats )	{
	int i;
	vector <int> temp;
	static int Count = 0;

	// Do counter
	if(Count % 60000 == 0) { cout << "\n"; }
	if(Count++ % 1000 == 0) { cout << "." << flush; }

	// Input the pattern
	if(Pos == (int)iters->size()) {
		temp.clear();
		FOR(i,(int)iters->size()) { temp.push_back(iters->at(i)); }
		Pats->push_back(temp);
		return;
	}
	// Do the recursion
	FOR(iters->at(Pos), NumStates(Type))	{ CreatePatterns(Pos+1,iters,Type,Pats); }
}

////////////////////////////////////////////////////////////
// Function to get gap mask; note: this uses the m_vsTrueSeq
void CData::GetGapMask(vector <vector <bool> > *Mask)	{
	int i,j,k;
	if(Mask == NULL) { Error("\nNeed to pass vector <vector <bool> > to CData::GetGapMask...\n"); }
	if(m_iNoSeq < 1) { Error("\nCan only create gap mask when CData::GetGapMask has sequences\n"); }
	assert(m_iTrueSize == (int)m_vsTrueSeq[0].size());
	Mask->clear();
	Mask->assign(m_iNoSeq,vector <bool>());
	// Get the mask
	FOR(i,m_iNoSeq) {
		FOR(j,m_iTrueSize)	{
			FOR(k,m_iChar) { if(m_vsTrueSeq[i][j] == m_sABET[k]) { break; } }
			if(k == m_iChar) { Mask->at(i).push_back(true); } else { Mask->at(i).push_back(false); }
	}	}
}
///////////////////////////////////////////////////////////////////
// Input data function
// This is the meaty function that actually puts the
// data into the CData object
// 	TODO:
//		Update to STL
// 		Update to remove all gap columns
//////////////////////////////////////////////////////////////////

void CData::InputData(EDataType Type, vector <string> cInputSeq, vector <string> cInputName, vector <int> SiteLabels, bool AllowFail)
{
    int i,j,k,l;                // Counters
    // Stuff for processing into patterns
	int **ariPat = NULL;				// Preliminary store patterns[m_iNoSeq][m_iSize]
	int *iTempPat = NULL;				// Temporary store for patterns
    int *ariPatOcc = NULL;             // Store for how often patterns have occurred
    double *ardFreqCount = NULL;		// Store for frequency counts for m_ariFreq
    int total=0;              // Counter
    int var_site=0;           // Number of variable sites
    int okay;                 // Flag to check what's going on
    vector <int> vtemp;
	// Set data size and number of sequences
	m_iNoSeq = (int) cInputSeq.size(); m_iTrueSize = m_iSize = (int) cInputSeq[0].size() / LenStates(Type);
	FOR(i,m_iNoSeq){ transform(cInputSeq[i].begin(),cInputSeq[i].end(),cInputSeq[i].begin(),(int(*)(int)) toupper); }
	// Do some checks regarding the sequence data
	assert(Type != NONE);
	assert(m_vsTrueSeq.empty());
	m_iChar = NumStates(Type);
	m_DataType = Type;
	FOR(i,m_iNoSeq)	{
		if((int)cInputSeq[i].size() != m_iSize * LenStates(Type))	{
			cout << "\nSequence " << i << " length of alignment passed to InputData seems incorrect: ";
			cout << "obs = " << cInputSeq[i].size() << " cf. exp=" << m_iSize * LenStates(Type) << "\n\n";
			if(cInputSeq[i][m_iSize] != '\0') { cout << "\nThis may be because there is no termination in the sequence added\n\n"; }
			FOR(j,(int)cInputSeq[i].size())	{ if(j%60 == 0) { cout << endl; } cout << cInputSeq[i][j]; }
			if(AllowFail) { return; } else { Error(); }
	}	}
	if(!m_ariPatMap.empty() || !m_ariPatOcc.empty() || !m_ariSeq.empty()) { Error("\nCData::InputData in a dirty object...\n\n"); }
	// Specify the alphabet
	m_sABET = DataStates(Type);
	// Initial memory allocation and initialise some variables
	m_vFreq.assign(m_iChar,0);
	GET_MEM(ardFreqCount,double,m_iChar+1); FOR(i,m_iChar+1){ ardFreqCount[i] = 0.0; }
	FOR(i,m_iNoSeq)	{ m_vsName.push_back(cInputName[i]); }
	GET_MEM(ariPat,int*,m_iNoSeq*3); FOR(i,m_iNoSeq) { GET_MEM(ariPat[i],int,m_iSize); }
	GET_MEM(iTempPat,int,m_iNoSeq*3);
	GET_MEM(ariPatOcc,int,m_iSize*3); FOR(i,m_iSize) { ariPatOcc[i] = 0; }
	m_ariPatMap.assign(m_iSize,0);

	/////////////////////////////////////////////////////////////////////////
    // Condense the data to its sufficient statistics
	/////////////////////////////////////////////////////////////////////////

    FOR(i,m_iSize)	{	// Go through sequences storing patterns and calculating frequencies
		okay=0;
        FOR(j,m_iNoSeq)	{
			// Identify which character (if unrecognised it is assumed to be a gap)
			k = FindState(Type,GetPos(cInputSeq[j],i, Type));
			ardFreqCount[k]+= 1.0; iTempPat[j] = k;	 // Pattern used is defined here
			if(iTempPat[j] != m_iChar) { okay ++; }
        }
//		FOR(j,m_iNoSeq) { cout << iTempPat[j] << " "<< flush; }
        // Remove pointless sites
#ifndef MATCH_PAML
		if(okay < 2) {  // Skips site if only one or zero sequences have a charactor at it
#else
		if(okay == 0)	{	// Skips site if there are no sequences with a character in it
#endif
//			cout << "\nRemoving site[" << i << "]: "; FOR(j,m_iNoSeq) { cout << cInputSeq[j][i]; }
			FOR(j,m_iNoSeq) { cInputSeq[j].erase(i*LenStates(Type),LenStates(Type)); }
			i--;
			m_iSize --;
			continue;
		}
		FOR(j,m_iNoSeq) { if(iTempPat[j] < m_iChar) { break; } } // Don't store sites which are all gaps...
		if(j == m_iNoSeq) { continue; }
        // See if this pattern already occurs in previous patterns
        okay=0;
		FOR(j,var_site)	{
            // Case when there is a match
            l=0; // Flag whether already seen or not
			FOR(k,m_iNoSeq)	{ if(ariPat[k][j] == iTempPat[k]) { l++; } }
            if(l==m_iNoSeq)  { okay=1;break; }
        }
		// Store the pattern mapping information
		m_ariPatMap[i] = j;
        // Add an occurence of this pattern
        if(okay==1)     { ariPatOcc[j]++; continue; }
		else	{	// If its not been seen before add it to the list
            ariPatOcc[j]=1;
			FOR(j,m_iNoSeq) { ariPat[j][var_site] = iTempPat[j]; }
            var_site++;
	}	}

	// Now the data are processed, store the true sequences
	m_iTrueSize = (int) cInputSeq[0].size() / LenStates(Type);
	FOR(i,m_iNoSeq) { m_vsTrueSeq.push_back(cInputSeq[i]); }

	/////////////////////////////////////////////////////////////////////////////
    // Dataset now condensed and ready to be transferred to the Object.
	/////////////////////////////////////////////////////////////////////////////
    m_iSize=var_site;
    vtemp.assign(m_iSize,0);
    m_ariSeq.assign(m_iNoSeq,vtemp);
	m_viNoChange.assign(m_iSize,-1);
	// Get the sequences
	FOR(i,m_iNoSeq) { FOR(j,m_iSize) { m_ariSeq[i][j] = ariPat[i][j]; } }
	// Get the places where no changes occur
	FOR(i,m_iSize)	{
		for(j=1;j<m_iNoSeq;j++) { if(m_ariSeq[j][i] != m_ariSeq[0][i]) { break; } }
		if(j == m_iNoSeq) {
			if(InRange(m_ariSeq[0][i],0,m_iChar)) { m_viNoChange[i] = m_ariSeq[0][i]; }
		}
	}
	m_ariPatOcc.assign(m_iSize,0);
	FOR(i,m_iSize) { m_ariPatOcc[i] = ariPatOcc[i]; }
    // The last thing we need to do is the base frequencies

	FOR(i,m_iChar)	{ total += (int) ardFreqCount[i]; }
	FOR(i,m_iChar)	{
		if(total == 0) { m_vFreq[i] = 0.0; continue; } // For blank data
		m_vFreq[i]=(ardFreqCount[i]/(double) total);
	}
	// Now to finish up we'll clear up the memory
	DEL_MEM(iTempPat);
	DEL_MEM(ardFreqCount);
	DEL_MEM(ariPatOcc);
	FOR(i,m_iNoSeq) { DEL_MEM(ariPat[i]); } DEL_MEM(ariPat);
	// Tell the data it is all okay
	m_bValid = true;
}

///////////////////////////////////////////////////////
// ofstream overloading for output

ostream& operator<<(ostream& os, const CData& DATA)	{
#if DEBUG_DATA
    os << "\nInternal members are: " << "\nm_iNoSeq: " << DATA.m_iNoSeq;
    os << "\nm_iSize: " << DATA.m_iSize;
#endif
    int i,j,k;
    os << "\n\nPatterned sequences are:\n";
    for(i=0;i<DATA.m_iSize;i++)	{
        os << "\nPattern " << i << " x " << DATA.m_ariPatOcc[i] << " \t  ";
        for(j=0;j<DATA.m_iNoSeq;j++)	{
            if(DATA.m_ariSeq[j][i]<DATA.m_iChar) { // If a member of the standard alphabet
				os << GetPos(DATA.m_sABET,DATA.m_ariSeq[j][i],DATA.m_DataType);;
			}
			else { FOR(k,LenStates(DATA.m_DataType)) { os << "-"; } }
			os << " ";
        }
//		FOR(j,DATA.m_iNoSeq) { os << " " << DATA.m_ariSeq[j][i]; }
    }
    return os;
}

void CData::OutRealData(ostream &os)	{
	int i;
	os << m_iNoSeq << "  " << m_iTrueSize << "\n";
	FOR(i,m_iNoSeq) {
		os << "\n" << m_vsName[i] << "  \t" << m_vsTrueSeq[i];
	}
}

////////////////////////////////////////////////////////////
// Function for getting character frequencies from data
vector <double> CData::GetFreq(int Seq)	{
	int i,j;
	vector <double> Freq(m_iChar,0.0);
	double Total = 0.0;
	if(Seq == -1)	{
		FOR(i,m_iNoSeq)	{
			FOR(j,m_iSize)	{
				if(m_ariSeq[i][j] == m_iChar) { continue; }
				Freq[m_ariSeq[i][j]] += (double) m_ariPatOcc[j]; Total += (double) m_ariPatOcc[j];
		}	}
	} else {
		FOR(j,m_iSize)	{
			if(m_ariSeq[Seq][j] == m_iChar) { continue; }
			Freq[m_ariSeq[Seq][j]] += (double) m_ariPatOcc[j]; Total += (double) m_ariPatOcc[j];
	}	}
	FOR(i,m_iChar) { Freq[i] /= Total; }
	return Freq;
}

//////////////////////////////////////////////////////////////
// Function for permanently removing invariant sites from an alignment
// This is written very inelegantly and inefficiently... Not fully error checked either
void CData::RemoveInvariantSites()	{
	int i,j,c1,NoInv,Pos, Char;
	int** Seqs = NULL;
	char cChar = '\0';
	vector <string> vsSeqs(10,"");
	vector <int> vtemp;
	///////////////// Do the invariant sites for the numerical values
	// Get memory
	GET_MEM(Seqs,int*,m_iNoSeq); FOR(i,m_iNoSeq) { GET_MEM(Seqs[i],int,m_iSize); }
	Pos = 0; NoInv = 0;
	FOR(i,m_iSize) {
		FOR(j,m_iNoSeq)	{ if(m_ariSeq[j][i] != m_iChar) { Char = m_ariSeq[j][i]; } }
		c1 = 0;
		FOR(j,m_iNoSeq)	{
			if(m_ariSeq[j][i] == m_iChar) { continue; }
			if(m_ariSeq[j][i] != Char) { break; }
			c1++;
		}
		if(c1 > 0 && j == m_iNoSeq) { NoInv++; continue; }
		// If not invariant, copy the site across
		FOR(j,m_iNoSeq) { Seqs[j][Pos] = m_ariSeq[j][i]; } Pos++;
	}
	if(Pos + NoInv != m_iSize) { Error("\nFailed CData::RemoveInvariantSites... Mismatch between number of sites...\n"); }
	m_iSize = Pos;
	// Transfer the sequence data across
	FOR(i,m_iNoSeq) { m_ariSeq[i].clear(); } m_ariSeq.clear();
	vtemp.assign(m_iSize, 0); m_ariSeq.assign(m_iNoSeq,vtemp);
	FOR(i,m_iNoSeq) { FOR(j,m_iSize) { m_ariSeq[i][j] = Seqs[i][j]; } }
	cout << "\nDoing to here!" << flush;
	////////////// Do the invariant sites for the strings
	Pos = 0; NoInv = 0;
	FOR(i,m_iTrueSize) {
		cout << "\nDoing i: "<< i << flush;
		FOR(j,m_iNoSeq)	{ if(IsDataType(m_DataType,m_vsTrueSeq[j][i])) { cChar = m_vsTrueSeq[j][i]; } }
		c1 = 0;
		FOR(j,m_iNoSeq)	{
			if(!IsDataType(m_DataType,m_vsTrueSeq[j][i])) { continue; }
			if(m_vsTrueSeq[j][i] != cChar) { break; }
			c1++;
		}
		if(c1 > 0 && j == m_iNoSeq) { NoInv++; continue; }
		// If not invariant, copy the site across
		FOR(j,m_iNoSeq) { vsSeqs[j] = vsSeqs[j] + m_vsTrueSeq[j][i]; } Pos++;
	}
	cout << "\nPos: " << Pos << " + Inv: " << NoInv << " == " << m_iTrueSize;
	if(Pos + NoInv != m_iTrueSize) { Error("\nFailed CData::RemoveInvariantSites... Mismatch between number of sites for real data...\n"); }
	m_iTrueSize = Pos;
	exit(-1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Function for removing sequences with little or no data in them (i.e. all or mostly gaps)
// Can also edit a tree if needed
void CData::RemoveSparseSeqs(bool Sparse,CTree *Tree, bool out) {
	int i,j,count, OriNoSeq = m_iNoSeq;
	double perc;
	bool Flag;
	string name;

/*	cout << "\n\n" << m_iNoSeq << "  " << m_iSize;
	FOR(i,m_iNoSeq) {
		cout << "\n" << m_vsName[i] << "  " << m_vsTrueSeq[i];
	}*/

	if(out) { cout << "\nChecking for "; if(Sparse) { cout << "sparse (>" << (1-MIN_DATA_PERCENT)*100 <<"% gaps)"; } else { cout << "all gap"; } cout << " sequences\n\tStarting with " << m_iNoSeq; }
	FOR(i,m_iNoSeq) {
//		cout << "\nChecking[" << i << "]" << flush;
		count = 0;
		FOR(j,(int)m_vsTrueSeq[i].size()) {
			if(IsGap(m_vsTrueSeq[i][j]) || (m_DataType == DNA && m_vsTrueSeq[i][j] == 'N')) { count++;}
		}
		perc = ((double) count / (double) m_vsTrueSeq[i].size());
//		cout << " ... percent gap: "<< perc;
		if((!Sparse && 1.0 - perc < FLT_EPSILON) || (Sparse && (1.0 - MIN_DATA_PERCENT) - perc  < FLT_EPSILON) ) {
//			cout << "\n\tRemoving["<<i<<"]: " << m_vsName[i] << " has " << perc *100 << "% gaps" << flush;
			name = m_vsName[i];
			if(RemoveSeq(i,Tree)) { Error("\nFailed to remove sequence " + name + " in RemoveSparseSeqs(...)\n"); } i--;
		}
	}
	if(out)	{ cout << " ... after removal there are " << m_iNoSeq << " sequences" << flush; }
	// Removing sequences can leave sites with just gaps
	if(m_iNoSeq != OriNoSeq) {
		Flag = false;
		FOR(j,m_iSize)	{
			FOR(i,m_iNoSeq)	{
				if(!IsGap(m_vsTrueSeq[i][j])) { break; }
			}
			if(i == m_iNoSeq) { Flag = true; break; }
		}
		if(Flag)	{
			vector <string> Seqs, Names;
			vector <int> Blank;
			EDataType Type = m_DataType;
			Seqs = m_vsTrueSeq; Names = m_vsName;
			Clean();
			InputData(Type, Seqs, Names,Blank);
		}
	}

}

// Condense gaps function to remove all sites where there are only gaps
void CData::CondenseGaps()		{
	int i,j,k;
	// Do the patterned bit
	FOR(j,m_iSize)	{		// Loop through the length of the sequence
		FOR(i,m_iNoSeq)	{
			if(m_ariSeq[i][j] == m_iChar) {
				m_iSize--;
				for(i=j;i<m_iSize;i++)	{ FOR(k,m_iNoSeq) { m_ariSeq[k][i] = m_ariSeq[k][i+1]; } m_ariPatOcc[i] = m_ariPatOcc[i+1]; }
				j--;
				break;
	}	}	}
	j=0; FOR(i,m_iSize) { j+= m_ariPatOcc[i];}
	// Do the TrueSequences
	int count = 0;
	FOR(j,(int)m_vsTrueSeq[0].size())	{
		FOR(i,m_iNoSeq)	{
			// Identify which character (if unrecognised it is assumed to be a gap)
			FOR(k,m_iChar)	{ if(toupper(m_vsTrueSeq[i][j])==toupper(m_sABET[k])) { break; } }
			if(k == m_iChar) {
				FOR(k,m_iNoSeq) { m_vsTrueSeq[k].replace(j,1,""); }
				count ++;
				j--;
				break;
			}
		}
	}
	m_iTrueSize = (int)m_vsTrueSeq[0].size();
}

////////////////////////////////////////////////////////////////
// Changes a DNA sequence into ÄDNA form
void CData::DNA2RY()	{
	int Seq,Site;
	if(m_DataType == RY) { return; }
	if(m_DataType != DNA) { Error("\nTrying to translate DNA to RY for wrong EDataType\n\n"); }
	assert(m_iTrueSize == (int) m_vsTrueSeq[0].size());
	vector <string> NewSeq(m_iNoSeq,"");
	vector <string> Names; Names = m_vsName;
	vector <int> Labels; Labels = m_viSiteLabels;
	// Create new sequences
	FOR(Site,m_iTrueSize) {
		FOR(Seq,m_iNoSeq) {
			if(m_vsTrueSeq[Seq][Site] == 'A' || m_vsTrueSeq[Seq][Site] == 'G')		{ NewSeq[Seq]+= 'R'; }
			else if(m_vsTrueSeq[Seq][Site] == 'C' || m_vsTrueSeq[Seq][Site] == 'T')	{ NewSeq[Seq]+= 'Y'; }
			else { NewSeq[Seq] += '-'; }

	}	}

	int Count1 = 0, Count2 =0;
	FOR(Site,m_iTrueSize) {
		FOR(Seq,m_iNoSeq) {
			if(m_vsTrueSeq[Seq][Site] == '-' && NewSeq[Seq][Site] != '-') { cout << "\n\tSeq["<<Seq << "][" << Site << "]"; }
			if(m_vsTrueSeq[Seq][Site] == '-') { Count1++; } if(NewSeq[Seq][Site] == '-') { Count2++; }
	}	}

	// Clean the memory
	Clean();
	// Create the new sequences
	InputData(RY,NewSeq,Names,Labels,false);
}

////////////////////////////////////////////////////////////////
// Cleans DNA data so it produces something matching in codons
// Ensures blocks of 3 alignable ACGTs
void CData::CleanToDNACodon()	{
	int codon,site,seq;
	bool Okay;
	if(m_DataType != DNA)  { Error("\nCannot CleanToDNACodon non DNA data...\n"); }
	if(m_iTrueSize%3 != 0) { Error("Can only build codon sequences from data divisible by 3...\n"); }
	assert(m_iTrueSize == (int) m_vsTrueSeq[0].size());
	vector <string> NewSeq(m_iNoSeq,"");
	vector <string> Names; Names = m_vsName;
	vector <int> Labels; Labels = m_viSiteLabels;
	// Create new sequences
	FOR(codon,m_iTrueSize/3) {
		FOR(seq,m_iNoSeq)	{
			Okay = true;
			// Check if a codon is okay
			FOR(site,3) { if(!(m_vsTrueSeq[seq][(codon*3)+site] == 'A' ||	m_vsTrueSeq[seq][(codon*3)+site] == 'C' || m_vsTrueSeq[seq][(codon*3)+site] == 'G' || m_vsTrueSeq[seq][(codon*3)+site] == 'T')) { Okay = false; break; } }
			if(Okay) 	{ FOR(site,3) { NewSeq[seq] += m_vsTrueSeq[seq][(codon*3)+site]; } }
			else 		{ NewSeq[seq] += "---"; }
		}
	}
	// Clean the memory
	Clean();
	// Create the new sequences
	InputData(DNA,NewSeq,Names,Labels,false);
}


// Remove sequence functions
int CData::RemoveSeq(int RemSeq,CTree *TREE)	{
	int i =0;
	string NewTree;
//	cout <<"\n\nRemoveSeq(" << RemSeq << ")";
	if(TREE!=NULL) {
		// Hash job of trimming the sequence from the tree. It's inefficient, but at least works
		TREE->OutName();
		TREE->SetNames(m_vsName,true);
		i+=TREE->RemoveLeafNode(RemSeq);
		if(RemSeq == 0) { TREE->SetStartCalc(1); }	// Get a new traversal position if the first sequence is removed.
		i+=RemoveSeq(RemSeq);
		ostringstream os; os << *TREE;
		TREE->CleanTree();
		TREE->CreateTree(os.str(),m_iNoSeq,true,false,false,this);
//		cout << "\n\tNew tree: " << *TREE;
	} else {
	// Remove the sequence from the data structure
		i+=RemoveSeq(RemSeq);
	}
	return i;
}

int CData::RemoveSeq(int RemSeq)	{
/* Removes a sequence from an alignments and adjust all required info
	returns 0 if successful or 1 otherwise	*/
	int i,j;
	double *ardFreqCount, total;
	GET_MEM(ardFreqCount,double,m_iChar);
	// Section to remove a sequence
	m_iNoSeq--;
	m_vsTrueSeq.erase(m_vsTrueSeq.begin() + RemSeq);
	m_vsName.erase(m_vsName.begin() + RemSeq);
	m_ariSeq.erase(m_ariSeq.begin() + RemSeq);

	// Now finished adjust the frequency counts
	for(i=0;i<m_iChar;i++) { ardFreqCount[i] = 0.0; } total = 0.0;
	for(i=0;i<m_iNoSeq;i++) {
		for(j=0;j<m_iSize;j++) {
			if(m_ariSeq[i][j] != m_iChar)	{ ardFreqCount[m_ariSeq[i][j]] += 1.0; total = total + 1.0; }
	}	}
	for(i=0;i<m_iChar;i++) { m_vFreq[i] = ardFreqCount[i] / total; }
	DEL_MEM(ardFreqCount);
	return 0;

}

// -------------------------- Function for adding a column to a data set --------------------------------
void CData::AddColumn(vector <int> Pattern, int Occur) {
	int i,j;
	vector <string> RealSeq(m_iNoSeq,"");
	vector <int> viTemp;
	// Data checking
	assert(m_DataType != NONE); assert(Pattern.size() == m_iNoSeq);
	if(m_vsTrueSeq.empty()) { m_vsTrueSeq.assign(m_iNoSeq,""); }
	if(m_ariSeq.empty())	{ m_ariSeq.assign(m_iNoSeq,viTemp); }
	FOR(i,Pattern.size()) {
		// Sort LenStates
		if(Pattern[i] == m_iChar) 	{ FOR(j,LenStates(m_DataType)) { RealSeq[i].push_back('-'); } continue; }
		if(!InRange(Pattern[i],0,m_iChar)) { Error("\nTrying to add pattern to data [CData::AddColumn(...)] that has characters out of range\n\n"); }
		else						{ FOR(j,LenStates(m_DataType)) { RealSeq[i].push_back(m_sABET[(Pattern[i]*LenStates(m_DataType))+j]); } }
	}
	// Add the sequence to the TrueSeq information
//	cout << "\nAdding Pattern: " << RealSeq << flush;
//	cout << "\nOriginal sequences: "; FOR(i,m_iNoSeq) { cout << "\nSeq[" << i<<"]: " << m_vsTrueSeq[i] << flush; }
	m_iTrueSize++;
	FOR(j,m_iNoSeq) { FOR(i,Occur) { m_vsTrueSeq[j] += RealSeq[j]; } }
//	cout << "\n---\n----\nUpdated sequences:  "; FOR(i,m_iNoSeq) { cout << "\nSeq[" << i<<"]: " << m_vsTrueSeq[i]; }
//	cout << "\nPattern: "; FOR(i,m_iNoSeq) { cout << Pattern[i] << " " << flush; }
	// Find data pattern in the current data set: m_ariSeq[Seq][Site]
	FOR(j,m_iSize) {
//		cout << "\n[" << j << "]: "; FOR(i,m_iNoSeq) { cout << m_ariSeq[i][j] << " "; }
		FOR(i,m_iNoSeq) { if(m_ariSeq[i][j] != Pattern[i]) { break; } }
		if(i != m_iNoSeq) { continue; }
		else {
//			cout << "\nFound pattern match: ";
//			cout << "\nOriginal data["<<j<<"]: "; FOR(i,m_iNoSeq) { cout << m_ariSeq[i][j] << " "; }
//			cout << "\nPattern data:           "; FOR(i,m_iNoSeq) { cout << Pattern[i] << " "; }
			m_ariPatOcc[j] += Occur;
			return;
		}
	}
	assert(j==m_iSize);
//	cout << "\nNew pattern...";
	FOR(i,m_iNoSeq) { m_ariSeq[i].push_back(Pattern[i]); }
	m_ariPatOcc.push_back(Occur);
	m_iSize++;

}

// Count the differences between two sequences;
double CData::PropDiff(int S1, int S2, bool IgnoreGaps)	{
	int i,Length = 0,Diff = 0;
	FOR(i,m_iSize)	{
		if(IgnoreGaps && (m_ariSeq[S1][i] == m_iChar || m_ariSeq[S2][i] == m_iChar) ) { continue; }
		if(m_ariSeq[S1][i] != m_ariSeq[S2][i]) { Diff += m_ariPatOcc[i]; }
		Length += m_ariPatOcc[i];
	}
	return (double) Diff / (double) Length;
}
// Some very basic distance calculations
double CData::PoissonDist(int S1, int S2)	{
	double P = PropDiff(S1,S2,true), n1 = ((double) m_iChar - 1.0) / (double) m_iChar;
	if(n1 == 0) { Error("Trying Poisson dist for 1/0 characters...\n"); }
	if(P >= n1) { return MAX_BRANCH; }
	if(1 - (P/n1) >= 1.0 - DBL_EPSILON) { return MAX_BRANCH; }
	return -n1 * log(1 - (P/n1));
}
double CData::PoissonVar(int S1, int S2)	{
	double p = PropDiff(S1,S2,true);
	return (p * (1-p)) / (pow((1 - ((m_iChar / (m_iChar - 1)) * p)),2) * m_iTrueSize);
}

double CData::GetPar(vector <string> s1, vector <string> s2)	{
	double Val = 1.0;
	return Val;
}

double CData::GuessKappa()	{
	// Check entry conditions
	if(!Valid()) { return INITIAL_KAPPA; }
	if(m_DataType != DNA) { return INITIAL_KAPPA; }
	assert(m_iChar == 4);
	// Initialise variables
	int No_Comp = ( m_iNoSeq* (m_iNoSeq - 1) ) / 2;
	int i,j,k,count; double Ts,Tv,UsedSize,c1,c2,c5,c6;
	double *Kappa =NULL,*VarKappa=NULL, Estimate;
	string Alfa = "ACGT";
	double Freq[4];
	// Get memory
	GET_MEM(Kappa,double,No_Comp); GET_MEM(VarKappa,double,No_Comp);
	// Transfer the frequencies into the right order
	FOR(i,4) {
		FOR(j,4) { if(m_sABET[j] == Alfa[i]) { break; } }
		assert(j!=4);
		Freq[i] = m_vFreq[j];
	}
	// Estimate kappas
	count = 0; FOR(i,m_iNoSeq)	{
		for(j=i+1;j<m_iNoSeq;j++)	{
			for(k=0,Ts = Tv = UsedSize = 0.0;k<m_iSize;k++)	{	// Get Value for each family
				if(m_ariSeq[i][k] == m_ariSeq[j][k]) { UsedSize+=m_ariPatOcc[k]; continue; }	// Skip when there is no change
				if(m_ariSeq[i][k] == m_iChar) { continue; }					// Skip gaps
				if(	m_sABET[m_ariSeq[i][k]] == 'A' && m_sABET[m_ariSeq[j][k]] == 'G' ||
					m_sABET[m_ariSeq[i][k]] == 'G' && m_sABET[m_ariSeq[j][k]] == 'A' ||
					m_sABET[m_ariSeq[i][k]] == 'C' && m_sABET[m_ariSeq[j][k]] == 'T' ||
					m_sABET[m_ariSeq[i][k]] == 'T' && m_sABET[m_ariSeq[j][k]] == 'C')	{
					Ts += m_ariPatOcc[k]; UsedSize += m_ariPatOcc[k]; continue;
				}
				Tv+=m_ariPatOcc[k]; UsedSize+=m_ariPatOcc[k];
			}
			Ts /= UsedSize; Tv /= UsedSize;

//			cout << "\nTs: "<< Ts << ", Tv: "<< Tv << " UsedSize: " << UsedSize;
//			cout << "\nTaking log of i): " << 1-(2*Ts)-Tv << " = " <<log(1-(2*Ts)-Tv)<< ", and ii): " <<1-(2*Tv) << " = " << log(1-(2*Tv));

			if(Tv >= 0.5 || Ts < DBL_EPSILON || Tv < DBL_EPSILON) { Kappa[count] = 0; VarKappa[count++] = 0.0; continue; }
			if((1-(2*Tv)) < DBL_EPSILON || (1-(2*Ts)-Tv) < DBL_EPSILON) { Kappa[count] = 0; VarKappa[count++] = 0.0; continue;}
			Kappa[count] = log(1-(2*Ts)-Tv)/log(1-(2*Tv)) - 1/2;
			c1 = 1/(1-(2*Ts)-Tv);
			c2 = 1/(1 - (2*Tv));
			c5 = (2*c1)/log(1 - (2*Tv));
			c6 = (c5 + ((4*c2*log(1 - (2*Ts) - Tv)) / pow(log(1 - (2 * Tv)),2)))/2;
			VarKappa[count] = UsedSize / ( (c5*c5*Ts) + (c6*c6*Tv) - pow((c5*Ts)+(c6*Tv),2) );

//			cout << "\nEst[" << count << "]:\tKappa: " << Kappa[count] << "\tVar: " << 1/VarKappa[count];

			count ++;
	}	}

	// Get estimate of kappa
	for(Ts = 0.0,i=0;i<count;i++)	{ Ts += VarKappa[i]; }
	if(Ts>DBL_EPSILON)	{
		for(Estimate = 0.0,i=0;i<count;i++) { Estimate += ( (VarKappa[i]/Ts) * Kappa[i] ); }
	}

//	cout << "\nScaled Kappa: "<< 2*Estimate;

	// Adjust for eqm frequencies
	Estimate = (Estimate * (Freq[0] + Freq[2]) * (Freq[1] + Freq[3]))/ ((Freq[0]*Freq[2]) + (Freq[1]*Freq[3]));

	// Check in bounds
	if(Estimate > 100) { Estimate = 100; }
	if(Estimate < 1) { Estimate = 1; }

//	cout << "\nAdjusted Kappa: "<< Estimate;

	// Empty memory
	DEL_MEM(Kappa); DEL_MEM(VarKappa);
	return Estimate;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Translate function
// ------------------
// If data is amino acid it returns immediately with warning
// If data is DNA it will do the translation providing it is divisible by 3 -- will throw error if problem occurs
// Genetic codes are specified in data.h

void CData::Translate(int GenCode)	{
	int site,seq,i;
	vector <string> vName;
	vector <string> vSeq;
	vector <int> SiteLabels;
	string Codon;

//	cout << "\nm_iTrueSize = " << m_iTrueSize << " : " << m_iTrueSize % 3;

	// Do some preprocessing
	if(m_DataType == AA) { cout << "\nWarning in CData::Translate(...): already amino acid sequence data\n"; return; }		// Already amino acid
	// Do some error checking
	if(!InRange(GenCode,0,11)) { Error("\nTrying to do DNA CData::Translate() with out of range GenCode\n"); }
	if(m_iSize < 3 || m_iTrueSize % 3 != 0) { Error("Trying CData::Translate(...) when DNA data not divisible by 3\n"); }
	assert(m_DataType == DNA); 	// Only works with DNA

	// Do some assignment
	vName = m_vsName;
	vSeq.assign(m_iNoSeq,string());
	FOR(seq,m_iNoSeq)	{
//		cout << "\nSeq["<<seq<<"]: " << m_vsTrueSeq[seq];
		vName.push_back(m_vsName[seq]);
		FOR(site,m_iTrueSize/3)	{
			Codon = "";
			FOR(i,3) { Codon += m_vsTrueSeq[seq][(site*3)+i]; }
			if(FindState(COD,Codon) == 64) {
				vSeq[seq] += '-';
			} else {
				if(!InRange(GenCodes[GenCode][FindState(COD,Codon)],0,64)) { Error("\nUnrecognised codon: " + Codon + " for data codon position " + int_to_string(site) + " in sequence " + int_to_string(seq) + "..."); }
				vSeq[seq] += State(AA,GenCodes[GenCode][FindState(COD,Codon)]);
			}
	}	}
	// Clean the data
	Clean();
	// Reassign everything
	InputData(AA,vSeq,vName,SiteLabels,false);
}

////////////////////////////////////////////////////////////////
// Function that changes DNA data into Codon data

void CData::MakeCodonData()	{
	vector <string> vName = m_vsName;
	vector <string> vSeq = m_vsTrueSeq;
	vector <int> SiteLabels = CData::m_viSiteLabels;
	// Do some initialisation and checking
	if(m_DataType == COD || m_DataType == COD_RED) { return; }
	if(m_DataType == AA) { Error("\nCannot make Codon data out of AA data in CData::MakeCodonData(...)\n"); }
	if(m_iTrueSize % 3 != 0) { Error("\nCannot make data into codons because DNA sequence length "+int_to_string(m_iTrueSize)+" is not divisible by 3"); }
	// Clean the data
	Clean();
	// Reassign everything
	InputData(COD,vSeq,vName,SiteLabels,false);
}

////////////////////////////////////////////////////////////////
// Function that takes standard codon data and removes all reference to stop codons
// For example changes a 64 state model to a 61 state for the universal code
void CData::ReduceCodonData(int GenCode) {
	int i,j,num; m_iChar = 64;
	string::iterator Current = m_sABET.end() - 3;
	vector<double>::iterator fCurrent = m_vFreq.end();
	if(ALLOW_CODON_REDUCTION == 0) { cout << "\nWarning: ALLOW_CODON_REDUCTION == 0"; }
	if(m_DataType == COD_RED) { if(GenCode != m_iGenCode) { Error("\nChanging genetic code in CData::ReduceCodonData()\n\n"); } return; }
	if(m_DataType != COD) { Error("\nCan not CData::ReduceCodonData() for non-codon data...\n"); }
	rFOR(num,64)	{
//		cout << "\nFor: " << *(Current)  << *(Current+1) << *(Current +2);
		fCurrent--;
		if(GenCodes[GenCode][num] == -1) { // Do the reduction
//			cout << "\nBeing removed...";
			FOR(i,m_iNoSeq) { FOR(j,m_iSize) { if(m_ariSeq[i][j] >= num) { m_ariSeq[i][j]--; } } }
			m_sABET.erase(Current,Current+3);
			m_vFreq.erase(fCurrent);
			m_iChar--;
		}
		Current-=3;
	}
	m_iGenCode = GenCode;
	m_DataType = COD_RED;
	m_vFreq = NormaliseVector(m_vFreq);
	// Also rebuild the m_viNoChange vector so it has the right dimensionality
	FOR(i,m_iSize) {
			if(m_viNoChange[i] != -1) {
			m_viNoChange[i] = m_ariSeq[0][i];
			if(!InRange(m_viNoChange[i],0,m_iChar)) {
				cout << "\nHaving a problem with m_viNoChange["<<i<<"]: set to: " << m_viNoChange[i];
				cout << "\n"; FOR(j,m_iNoSeq) { cout << "[" << j <<"]: " << m_ariSeq[j][i] << " "; }
				cout << "\nHave set it to: " << m_ariSeq[0][i] << " == " << m_viNoChange[i];
				cout << "\nAnd the sequences..."; FOR(j,m_iNoSeq) { cout << "\n["<<j<<"]: " << m_vsTrueSeq[j]; }
			}
			assert(InRange(m_viNoChange[i],0,m_iChar));
		}
	}


//	cout << "\nNew frequency vector:"; FOR(i,m_iChar) { cout << "\n" << m_sABET.substr(i*3,3) << " == " << m_vFreq[i]; }
//	cout << "\nState-space:"; FOR(i,m_iChar) { cout << "\n["<<i<<"]: " << m_sABET.substr(i*3,3); }
//	cout << "\nStates = " << m_iChar;
}

void CData::ExpandCodonData(int GenCode) {
	int i,j,num;
	if(m_DataType == COD_RED && GenCode != m_iGenCode) { Error("\nChanging genetic code in CData::ReduceCodonData()\n\n"); }
	if(m_DataType != COD_RED) { return; }
	if(ALLOW_CODON_REDUCTION == 0) { cout << "\nWarning: ALLOW_CODON_REDUCTION == 0"; }
	FOR(num,64)	{
		if(GenCodes[GenCode][num] == -1) { // Do the reduction
			FOR(i,m_iNoSeq) { FOR(j,m_iSize) { if(m_ariSeq[i][j] >= num) { m_ariSeq[i][j]++; } } }
	}	}

}

//////////////////////////////////////////////////////////////////
// Function that can be used to extract a subset of codon positions
// ---
// Extracts all the bools that are true
bool CData::GetCodonPositions(bool First, bool Second, bool Third)      {
        int i,j,pos;
        vector <string> Names, Sequences;
        vector <int> SiteLabels;
        assert(!m_vsTrueSeq.empty());
        assert(m_vsTrueSeq[0].size() == m_iTrueSize);
        if(m_iTrueSize %3 != 0) { cout << "\nTrying to extract codon positions in GetCodonPositons(...) with data not divisible by 3. " << m_iTrueSize << "% 3 = " << m_iTrueSize % 3 << " ...\n"; Error("Boom"); }
        // Get Names
        Names = m_vsName;
        // DEBUG STUFF IF NEEDED
//        cout << "\nInto CData::GetCodonPositions : Data returned will be of form (" << First << ","<< Second << "," << Third << ");";
        // Get sequences
        Sequences.assign(m_iNoSeq,"");
        for(i=0;i<m_iTrueSize;i+=3) {
                // Do first position
                if(First) {
                        pos = 0;
                        FOR(j,m_iNoSeq) { Sequences[j] += m_vsTrueSeq[j][i+pos]; }
                }
                // Do second position
                if(Second) {
                        pos = 1;
                        FOR(j,m_iNoSeq) { Sequences[j] += m_vsTrueSeq[j][i+pos]; }
                }
                // Do third position
                if(Third) {
                        pos = 2;
                        FOR(j,m_iNoSeq) { Sequences[j] += m_vsTrueSeq[j][i+pos]; }
                }
        }

        Clean();
        InputData(DNA,Sequences,Names,SiteLabels,false);
        return true;
}


//////////////////////////////////////////////////////////////////
// Functions that return data for Tree HMMs

CData *CData::MakeMatchData()	{
	int i,j;
	CData *MatchData = NULL;
	assert(!m_vsTrueSeq.empty());
//	cout << "\n--- Creating Match states ---";
	string NewSeq;
	vector <string > Toks;
	vector <string > NewData;
	FOR(i,m_iNoSeq) {
		NewSeq.clear();
		Toks = Tokenise(GetProfilePath(m_vsTrueSeq[i],m_sABET));
		FOR(j,(int)Toks.size()) {
			assert((int)Toks[j].size() == 2);
			if(Toks[j][0] == 'M') {
				if(Toks[j][1] == 'M') { NewSeq = NewSeq + "R"; }
				if(Toks[j][1] == 'D') { NewSeq = NewSeq + "Y"; }
			} else { NewSeq = NewSeq + "-"; }
		}
		NewData.push_back(NewSeq);
//		cout << "\n\nOriSeq: " << m_vsTrueSeq[i] << "\nTokSeq: ";
//		FOR(j,(int)Toks.size()) { cout << Toks[j] << " "; }
//		cout << "\nMatSeq: " << NewSeq;
	}
	MatchData = new CData(m_iNoSeq,m_iTrueSize-1,NewData,m_vsName,RY);
//	cout << "\n" << m_iNoSeq << " " << m_iTrueSize-1;
//	FOR(i,m_iNoSeq) { cout << "\n\n" << m_vsName[i] << "\t"; FOR(j,m_iTrueSize-1) { if(NewData[i][j]=='R') { cout << "A"; } else if(NewData[i][j]=='Y') { cout << "T"; } else { cout << "-"; } } }
//	cout << "\n\n";
	return MatchData;
}

CData *CData::MakeDeleteData()	{
	int i,j;
	CData *DeleteData = NULL;
	assert(!m_vsTrueSeq.empty());
//	cout << "\n--- Creating Delete states ---";
	string NewSeq;
	vector <string > Toks;
	vector <string > NewData;
	FOR(i,m_iNoSeq) {
		NewSeq.clear();
		Toks = Tokenise(GetProfilePath(m_vsTrueSeq[i],m_sABET));
		FOR(j,(int)Toks.size()) {
			assert((int)Toks[j].size() == 2);
			if(Toks[j][0] == 'D') {
				if(Toks[j][1] == 'D') { NewSeq = NewSeq + "R"; }
				if(Toks[j][1] == 'M') { NewSeq = NewSeq + "Y"; }
			} else { NewSeq = NewSeq + "-"; }
		}
		NewData.push_back(NewSeq);
//		cout << "\n\nOriSeq: " << m_vsTrueSeq[i] << "\nTokSeq: ";
//		FOR(j,(int)Toks.size()) { cout << Toks[j] << " "; }
//		cout << "\nMatSeq: " << NewSeq;
	}
	DeleteData = new CData(m_iNoSeq,m_iTrueSize-1,NewData,m_vsName,RY);
//	cout << "\n" << m_iNoSeq << " " << m_iTrueSize-1;
//	FOR(i,m_iNoSeq) { cout << "\n\n" << m_vsName[i] << "\t"; FOR(j,m_iTrueSize-1) { if(NewData[i][j]=='R') { cout << "A"; } else if(NewData[i][j]=='Y') { cout << "T"; } else { cout << "-"; } } }
//	cout << "\n\n";

	return DeleteData;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions relating to creating pairwise data
// ---
// Takes one (where SecondData == NULL) or two data sets and return a new object with all pairs between those data sets

// Blank data array
//	CData(int NoSeq, int Size, EDataType Type, vector <string> *Names = NULL);

CData * CData::CreateAllPairs(CData *SecondData, bool output) 	{
	int i,j,sp;
	CData *NewData;
	string NewABET;
	vector <string> NewSeq(m_iNoSeq,"");
	vector <int> Pattern(m_iNoSeq,-1) , Mapping(2,-1);
	if(SecondData == NULL) { SecondData = this; }
	// Some basic error checking
	assert(m_iNoSeq == SecondData->m_iNoSeq); assert(m_DataType == SecondData->m_DataType);
	// Create new data structure
	switch(m_DataType) {
	case RY: 	NewData = new CData(m_iNoSeq,0,RY2,&m_vsName);  break;
	case DNA:	NewData = new CData(m_iNoSeq,0,DNA2,&m_vsName); break;
	case AA:	NewData = new CData(m_iNoSeq,0,AA2,&m_vsName);  break;
	default:	Error("\nUnexpected datatype in CData::CreateAllPairs...\n"); exit(-1);
	};
	// Create the new sequences
	// Two options depending on whether one data set or two
	// 1. One data set
	//		- each pairwise comparison should only be made once (assumes no direction to interaction i.e A->B == B->A)
	//		- self comparisons are meaningless

	cout << "\nCreating all pairs... " << flush;

	if(output) { cout << "\nCreating paired sequences"; }
	if(SecondData == this) {
		FOR(i,m_iSize) {
			for(j=i+1;j<m_iSize;j++) {
				FOR(sp,m_iNoSeq) {
//					cout << "\n["<<sp<<"]["<<i<<","<<j<<"]: " << State(AA,m_ariSeq[sp][i]) << ":" << State(AA,SecondData->m_ariSeq[sp][j]) << flush;
					if(m_ariSeq[sp][i] == m_iChar || SecondData->m_ariSeq[sp][j] == m_iChar) { Pattern[sp] = m_iChar * m_iChar; }
					else { Pattern[sp] = (m_iChar*m_ariSeq[sp][i]) + SecondData->m_ariSeq[sp][j]; }
//					cout << " cf. [" << Pattern[sp] << "]: " << State(AA2,Pattern[sp]) << flush;
//					NewSeq[sp].push_back(m_vsTrueSeq[sp][i]); NewSeq[sp].push_back(SecondData->m_vsTrueSeq[sp][j]);
				}
//				cout << "\nAdding Column" << flush;
				NewData->AddColumn(Pattern,m_ariPatOcc[i]*SecondData->m_ariPatOcc[j]);
				Mapping[0] = i; Mapping[1] = j; NewData->m_vviCoevoMapping.push_back(Mapping);
//				cout << " ... done" << flush;
		} 	}
	} else {
		// 2. Two data sets
		cout << "\nWarning two data sets is not checked yet... Current implementation is clearly wrong!!!" << flush; exit(-1);
		FOR(i,m_iTrueSize) {
			FOR(j,SecondData->m_iTrueSize) {
				FOR(sp,m_iNoSeq) {
//					cout << "\n["<<sp<<"]["<<i<<","<<j<<"]: " << m_vsTrueSeq[sp][i] << ":" << SecondData->m_vsTrueSeq[sp][j];
					NewSeq[sp].push_back(m_vsTrueSeq[sp][i]); NewSeq[sp].push_back(SecondData->m_vsTrueSeq[sp][j]);
	}	}	}	}
//	FOR(i,m_iNoSeq) { cout << "\nNewSeq[" << i << "]: " << NewSeq[i]; }
	// Create the new data object
	if(output) { cout << " ... placing them into a new data object"; }
	// Clean up
	SecondData = NULL;
	NewData->m_vFreq = NewData->GetFreq();
	NewData->m_bValid = true;
	return NewData;





/*	OLD VERSION
 * 	---
 *
	int i,j,sp;
	CData *NewData;
	string NewABET;
	vector <string> NewSeq(m_iNoSeq,"");
	if(SecondData == NULL) { SecondData = this; }
	// Some basic error checking
	assert(m_iNoSeq == SecondData->m_iNoSeq); assert(m_DataType == SecondData->m_DataType);
	// Create the new sequences
	// Two options depending on whether one data set or two
	// 1. One data set
	//		- each pairwise comparison should only be made once (assumes no direction to interaction i.e A->B == B->A)
	//		- self comparisons are meaningless
	if(output) { cout << "\nCreating paired sequences"; }
	if(SecondData == this) {
		FOR(i,m_iTrueSize) {
			for(j=i+1;j<m_iTrueSize;j++) {
				FOR(sp,m_iNoSeq) {
//					cout << "\n["<<sp<<"]["<<i<<","<<j<<"]: " << m_vsTrueSeq[sp][i] << ":" << SecondData->m_vsTrueSeq[sp][j];
					NewSeq[sp].push_back(m_vsTrueSeq[sp][i]); NewSeq[sp].push_back(SecondData->m_vsTrueSeq[sp][j]);
		} } }
	} else {
		// 2. Two data sets
		cout << "\nWarning two data sets is not checked yet..." << flush; exit(-1);
		FOR(i,m_iTrueSize) {
			FOR(j,SecondData->m_iTrueSize) {
				FOR(sp,m_iNoSeq) {
//					cout << "\n["<<sp<<"]["<<i<<","<<j<<"]: " << m_vsTrueSeq[sp][i] << ":" << SecondData->m_vsTrueSeq[sp][j];
					NewSeq[sp].push_back(m_vsTrueSeq[sp][i]); NewSeq[sp].push_back(SecondData->m_vsTrueSeq[sp][j]);
	}	}	}	}
//	FOR(i,m_iNoSeq) { cout << "\nNewSeq[" << i << "]: " << NewSeq[i]; }
	// Create the new data object
	if(output) { cout << " ... placing them into a new data object"; }
	switch(m_DataType) {
	case RY: 	NewData = new CData(m_iNoSeq,NewSeq[0].size(),NewSeq,m_vsName,RY2); break;
	case DNA:	NewData = new CData(m_iNoSeq,NewSeq[0].size(),NewSeq,m_vsName,DNA2); break;
	case AA:	NewData = new CData(m_iNoSeq,NewSeq[0].size(),NewSeq,m_vsName,AA2); break;
	default:	Error("\nUnexpected datatype in CData::CreateAllPairs...\n"); exit(-1);
	};
	// Clean up
	SecondData = NULL;
	return NewData;
*/
}

// Function that finds the mapping between a pair of data sets and the currently stored (Joint) data set
bool CData::CreateAllPairMapping(CData *D1, CData *D2) {
	int i,j,site;

	switch(m_DataType) {
	case RY2: assert(D1->m_DataType == RY && D2->m_DataType == RY);
	case AA2: assert(D1->m_DataType == AA && D2->m_DataType == AA);
	case DNA2: assert(D1->m_DataType == DNA && D2->m_DataType == DNA);
	};
	FOR(site,m_iSize) {

	}
}

//////////////////////////////////////////////////////////////////////////////////////
// Get the scaling factor between amino acid models and codon models
// ---
// Idea and much of the implementation is based on Seo and Kishino
// The formula is: Adj = \Sum_{i} \Sum_{j} log(CodonFreq[i][j] / AminoAcidFreq[i][j])
// The novelty w.r.t. Seo and Kishino comes from the calculation of the frequencies and how they're treated as parameters in the model
// There are several possible parameterisations based on how those frequencies are produced
// 1. Amino acid frequencies
//		-- Those fed in from the model (e.g. database frequencies for LG without the +F option); 0 d.f.
//		-- The empirical frequencies estimated from the data; additional 19 d.f.
// 2. Codon frequencies
//		-- No information about codon frequencies. The frequency of a codon CodonFreq[i][j] is AminoAcidFreq[i][j] / #degenerate_codons
//		-- Information about coding frequencies as parameters in the model. There are a very wide range of possible parameterisations here,
//			but we will only consider the case of empirically estimated sense codon frequencies, EmpCodon[i][j], taken direct from the data.
//			Let K = {k_i,k_j} be the set of codons associated with AminoAcidFreq[i][j] so
//				CodonFreq[i][j] = AminoAcidFreq[i][j] * ( EmpCodon[i][j] / \Sum_{k \in K} EmpCodon[k_i][k_j] );
//			Note when the AminoAcidFreq are taken empirically from the data, then this value is equivalent of:
//				CodonFreq[i][j] = EmpFreq[i][j] because \Sum_{k \in K} EmpCodon[k_i][k_j] = AminoAcidFreq[i][j]
//			This option adds (60 df from codon model - 19 df. from amino acid model =) 41 degrees of freedom to the model
// ---
// This approach means that for amino acid models there are 4 possible scaling factors.
// In this function the best scaling factor is chosen according to AIC, and the return value has the AIC penalty associated with it.
// *****************************************************************************
// Update: 31 Oct 2012
// ---
// I have realised that there are actually only 2 possible scaling factors regardless of which model is examined
// Consider equation (6) from Seo and Kishino:
//		Adj = \Sum_{i} \Sum_{j} log(CodonFreq[i][j] / AminoAcidFreq[i][j])
// The factors CodonFreq[i][j] / AminoAcidFreq[i][j] are effectively the proportions by which a codon codes for an amino acid
// So there are two natural ways to treat these quantities when the codon frequencies are not part of the model
// i. Treat them as equiprobable, which means the adjustment factor is:
//		1 / Number_of_codons[i][j] and this approach adds 0 degrees of freedom to the model
// ii. Take the ML estimates of the codon frequencies, which minimises the adjustment factor.
//		The way to maximise the adjustment factor is to find the set of P[k] frequencies that maximise the function
//			adj = \Sum_{i} \Sum{j} log ( P[{site[i][j]}] );
//		This is the equivalent of obtaining the MLEs for frequencies for a multinomial distribution, which is simply taking the empirical frequencies
// 		This approach is the equivalent of using the empirical codon frequencies (EmpFreq) and empirical amino acid frequencies (EmpAA):
//		EmpFrq[i][j] / EmpAA[i][j]
// 		And adds (#CodingCodons - 20)  degrees of freedom to the model
// Note that this approach then makes the correction independent of amino acid model frequencies, allowing the correction to be applied for all models
double CData::GetAminoToCodonlnLScale(int GeneticCode, int *df)	{
	bool AAFrDiff = false;
	int i, count;
	double EqFrq, EmFrq, RetValue = 0.0;
	vector <double> DataAAFreq(20,0), EquCod, EmpCod;
	if(m_DataType != DNA) { Error("\nError: CData::GetAminoToCodonlnLScale() only works for NUCLEOTIDE data\n"); }
	if(!InRange(GeneticCode,0,11)) { Error("\nTrying to do DNA CData::GetAminoToCodonlnLScale() with out of range GenCode\n"); }
	// Initialise objects associated with the function
	CData AA_Data = *this; AA_Data.Translate(GeneticCode); DataAAFreq = AA_Data.m_vFreq;
	CData COD_Data = *this; COD_Data.MakeCodonData();	EmpCod = COD_Data.m_vFreq; // Note the genetic code here is in *64* state-space
	EquCod.assign(COD_Data.m_iChar,0);
	// Do some basic error checking
	assert(AA_Data.m_iTrueSize == COD_Data.m_iTrueSize && AA_Data.m_iNoSeq == COD_Data.m_iNoSeq);
	assert((int)DataAAFreq.size() == 20);
	assert(df != NULL); *df = 0;
	// Do the calculations
	// 1. Set up the frequencies
	count = 0; FOR(i,64) { if(GenCodes[GeneticCode][i] != -1) { count++; } }
	EquCod.assign(64, (double) 1 / (double) count );
	FOR(i,64) { if(GenCodes[GeneticCode][i] == -1) { EquCod[i] = 0.0; } }
	EquCod = EnforceAAFreqOnCodon(EquCod,DataAAFreq,GeneticCode);
	EmpCod = EnforceAAFreqOnCodon(EmpCod,DataAAFreq,GeneticCode);
	// 2. Get the two possible adjustment factors
	RetValue = EqFrq = GetAdjustmentScore(&AA_Data,&COD_Data,DataAAFreq,EquCod,GeneticCode); *df = 0;
	EmFrq = GetAdjustmentScore(&AA_Data,&COD_Data,DataAAFreq,EmpCod,GeneticCode);
	// Organise the return value if AIC(EmFrq) > AIC(EqFrq)
	if(EmFrq - count + 20 > EqFrq) {
		RetValue = EmFrq; *df = count - 20;
	}
//	cout << "\nEmFrq: " << EmFrq  << " (" << EmFrq - count + 20 << ") and EqFrq: " << EqFrq << " df = " << *df << " (0:" << count - 20 << ")";
//	if(EmFrq - count + 20 > EqFrq) { cout << " *EMP*"; } else { cout << " *EQU*"; } cout << " -- return: " << RetValue << endl;
	// Return
	return RetValue;
}


// Function returns the likelihood and df the degrees of freedom for the correction
/*  OLD VERSION THAT PRODUCED THE CORRECTION IN LOTS OF DIFFERENT WAYS
double CData::GetAminoToCodonlnLScale(int GeneticCode, vector <double> ModelAAFreq, int *df)	{
	bool AAFrDiff = false;
	int i, count;
	double EqFrq, EmFrq, RetValue = 0.0;
	double Adj_ModelF_EQUCod, Adj_ModelF_EmpCod, Adj_EmpF_EQUCod, Adj_EmpF_EmpCod;			 // The store for the different adjustments
	vector <double> DataAAFreq(20,0), DataEquCod, DataEmpCod, ModelEquCod, ModelEmpCod;

	cout << "\n----------------- CData::GetAminotoCodonlnLScale(...) -----------------\n";

	if(m_DataType != DNA) { Error("\nError: CData::GetAminoToCodonlnLScale() only works for NUCLEOTIDE data\n"); }
	if(!InRange(GeneticCode,0,11)) { Error("\nTrying to do DNA CData::GetAminoToCodonlnLScale() with out of range GenCode\n"); }
	// Initialise objects associated with the function
	Adj_ModelF_EQUCod = Adj_ModelF_EmpCod = Adj_EmpF_EQUCod = Adj_EmpF_EmpCod = -BIG_NUMBER;
	CData AA_Data = *this; AA_Data.Translate(GeneticCode);
	CData COD_Data = *this; COD_Data.MakeCodonData();	// Note the genetic code here is in *64* state-space
	DataEquCod.assign(COD_Data.m_iChar,0); DataEmpCod = ModelEquCod = ModelEmpCod = DataEquCod;
	// Do some basic error checking
	assert(AA_Data.m_iTrueSize == COD_Data.m_iTrueSize && AA_Data.m_iNoSeq == COD_Data.m_iNoSeq);
	assert((int)DataAAFreq.size() == 20);
	assert((int)ModelAAFreq.size() == 20);
	assert(df != NULL); *df = 0;
	// Make the various sets of frequencies
	// a. Get data amino acid frequencies
	DataAAFreq = AA_Data.m_vFreq;
	FOR(i,20) { cout << "\ntest: [" << DataAAFreq[i] << " - " << ModelAAFreq[i] << "] " << fabs(DataAAFreq[i] - ModelAAFreq[i]); if(fabs(DataAAFreq[i] - ModelAAFreq[i]) > 1.0E-3) { AAFrDiff = true; break; } }
	// b. Get the empirical codon frequencies via function EnforceAAFreqOnCodon(...)
	count = 0; FOR(i,64) { if(GenCodes[GeneticCode][i] != -1) { count++; } }
	DataEquCod.assign(64, (double) 1 / (double) count );
	FOR(i,64) { if(GenCodes[GeneticCode][i] == -1) { DataEquCod[i] = 0; } }
	ModelEquCod = DataEquCod;
	DataEquCod = EnforceAAFreqOnCodon(DataEquCod,DataAAFreq,GeneticCode);
	ModelEquCod = EnforceAAFreqOnCodon(ModelEquCod,ModelAAFreq,GeneticCode);
	// c. Get the other codon frequencies
	DataEmpCod = ModelEmpCod = COD_Data.m_vFreq;
	ModelEmpCod = EnforceAAFreqOnCodon(ModelEmpCod,ModelAAFreq,GeneticCode);

	cout << "\nAA\nData \t" << DataAAFreq << "\nModel\t" << ModelAAFreq;
	cout << "\n>>>\nCodon\nModelEqu\t" << ModelEquCod << "\nModelEmp\t" << ModelEmpCod;
	cout << "\nDataEqu \t" << DataEquCod << "\nDataEmp \t" << DataEmpCod << flush;

	// d. Count the number of non-stop codons
	count = 0; FOR(i,64) { if(GenCodes[GeneticCode][i]!=-1) { count++; } }
	// Compute the different adjustment factors Adj_ModelF_EQUCod, Adj_ModelF_EmpCod, Adj_EmpF_EQUCod, Adj_EmpF_EmpCod;
	// 1. Equiprobable codons, with model's amino acid frequencies. 0 extra degrees of freedom
	RetValue = Adj_ModelF_EQUCod = GetAdjustmentScore(&AA_Data,&COD_Data,ModelAAFreq,ModelEquCod,GeneticCode);
	*df = 0;
	// 2. Empirical codons, with model's amino acid frequencies. (count - 1) - 19 extra degrees of freedom
	Adj_ModelF_EmpCod =  GetAdjustmentScore(&AA_Data,&COD_Data,ModelAAFreq,ModelEmpCod,GeneticCode);
	if(Adj_ModelF_EmpCod < RetValue) { RetValue = Adj_ModelF_EmpCod; *df = (count - 20); }
	// Do the remaining adjustments if needed
	if(AAFrDiff)	{
		// 3. Equiprobobable codons, with empirical amino acid frequencies. 0 extra degrees of freedom
		Adj_EmpF_EQUCod = GetAdjustmentScore(&AA_Data,&COD_Data,DataAAFreq,DataEquCod,GeneticCode);
		if(Adj_EmpF_EQUCod < RetValue) { RetValue = Adj_EmpF_EQUCod; *df = 0;}
		// 4. Empirical codons, with empirical amino acid frequencies. (count - 1) - 19 extra degrees of freedom
		Adj_EmpF_EmpCod = GetAdjustmentScore(&AA_Data,&COD_Data,DataAAFreq,DataEmpCod,GeneticCode);
		if(Adj_EmpF_EmpCod  < RetValue) { RetValue = Adj_EmpF_EmpCod; *df = (count - 20); }
	}
	cout << "\nThe adjustment factors:\n\tCod_equ + Mod_aa = " << Adj_ModelF_EQUCod;
	cout << "\n\tCod_emp + Mod_aa = " << Adj_ModelF_EmpCod << "\n\tCod_equ + Emp_aa = " << Adj_EmpF_EQUCod << "\n\tCod_emp + Emp_aa = " << Adj_EmpF_EmpCod;

	// Return
	return RetValue;
}
*/
#define GetAdjustmentScore_DEBUG 1
double CData::GetAdjustmentScore(CData *AA_Data, CData *COD_Data, vector <double> AAFreq, vector <double> CodFreq, int GenCode) {
	int i, j;
	double RetValue = 0.0;
	// Some basic error checking
	assert(AA_Data->m_DataType == AA); assert(COD_Data->m_DataType == COD);
	assert(AAFreq.size() == 20); assert(CodFreq.size() == 64);
	assert(AA_Data->m_iTrueSize == COD_Data->m_iTrueSize && AA_Data->m_iNoSeq == COD_Data->m_iNoSeq);
	// If required do the hard debug
#if GetAdjustmentScore_DEBUG == 1
	int count;
	double Checker = 0.0;
	FOR(i,20) {	// Loop through the amino acids
		Checker = 0.0;
		FOR(j,64)	{ // Loop through the codons
			if(GenCodes[GenCode][j] == i) { assert(count < CodFreq.size()); Checker += CodFreq[j]; }
		}
		if(fabs(Checker - AAFreq[i]) > 1.0E-5) { Error("Frequencies in GetAdjustmentScore(...) do not match...\n\n"); }
	}
#endif
/*	cout << "\nDEBUG CHECKING";
	cout << "\nCODS:"; FOR(i,64) { cout << "\t" << State(COD,i); }
	cout << "\nCodFrqs1:\t" << CodFreq;
	cout << "\nCodFrqs2:\t" << COD_Data->m_vFreq;
	cout << "\nCodons:"; FOR(i,64) { cout << "\t" << GenCodes[GenCode][i]; }
	cout << "\nAAFrq:\t" << AAFreq;
*/	// Do the computations. Arbitrarily indexing on the amino acid sequence.
	// Note I have to work in the true sequence space because of the redundancy of the amino acid sequence
	RetValue = 0.0;
	FOR(i,AA_Data->m_iNoSeq)	{
		FOR(j,AA_Data->m_iTrueSize)	{
//			cout << "\n[" << i<< "][" << j << "]: ";
			if(AA_Data->m_ariSeq[i][AA_Data->m_ariPatMap[j]] == AA_Data->m_iChar) { continue; } // Skip gaps
//			cout << State(COD,COD_Data->m_ariSeq[i][COD_Data->m_ariPatMap[j]]) << " = " << CodFreq[COD_Data->m_ariSeq[i][COD_Data->m_ariPatMap[j]]] << " : " << AAFreq[AA_Data->m_ariSeq[i][ AA_Data->m_ariPatMap[j]]];
			RetValue += log(CodFreq[COD_Data->m_ariSeq[i][COD_Data->m_ariPatMap[j]]] / AAFreq[AA_Data->m_ariSeq[i][ AA_Data->m_ariPatMap[j]]]);
		}
	}
	return RetValue;
}


double CData::OldGetAminoToCodonlnLScale(int GeneticCode)	{
	int i, j;
	double RetValue = 0.0;
	if(m_DataType != DNA) { Error("\nError: CData::GetAminoToCodonlnLScale() only works for NUCLEOTIDE data\n"); }
	if(!InRange(GeneticCode,0,11)) { Error("\nTrying to do DNA CData::GetAminoToCodonlnLScale() with out of range GenCode\n"); }
	CData AA_Data = *this; AA_Data.Translate(GeneticCode);
	CData COD_Data = *this; COD_Data.MakeCodonData();
	assert(AA_Data.m_iTrueSize == COD_Data.m_iTrueSize && AA_Data.m_iNoSeq == COD_Data.m_iNoSeq);
	assert((int)COD_Data.m_vFreq.size() == 64 && (int)AA_Data.m_vFreq.size() == 20);
	FOR(i,AA_Data.m_iNoSeq) {
		FOR(j,AA_Data.m_iTrueSize) {
			if(AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]] == AA_Data.m_iChar) { continue; } // Skip gaps
			 cout << "\n["<<i<<"]["<<j<<"]: " <<  COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] << " / " << AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]];
			 cout << " = " << COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]] << " / " << AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]];
			 cout << " = " << COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]];
			 cout << " --> " << log(COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]]);
			RetValue += log(COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]]);
		}
	}
	return RetValue;
}

// Simple function for nucleotide sequences
double CData::GetNT2MissinglnLScale()	{
	int site,sp;
	double Adj = 0.0;

	vector <int> FreqCount(4,0);

	assert(m_DataType == DNA);
	assert(fabs(Sum(&m_vFreq) -1.0) < 1.0E-5);
	FOR(site,m_iSize)	{
		FOR(sp,m_iNoSeq) {
			if(m_ariSeq[sp][site] == m_iChar) { continue; } // skip gaps
			FreqCount[m_ariSeq[sp][site]] += m_ariPatOcc[site];
			Adj += log(m_vFreq[ m_ariSeq[sp][site] ]) * m_ariPatOcc[site];
		}
	}

	cout << "\nFreqs  " << m_vFreq;
	cout << "\nCounts " << FreqCount;



	return Adj;
}

////////////////////////////////////////////////////////////////////////////////////
// Function that produces corrects a set of codon frequencies so that they match a set of amino acid frequencies
// For a codon c_{i,j} and the corresponding amino acid a_c_{i,j} we know
//   a_c_{i,j} = \sum_{c_{i,j} \in a_c_{i,j}} c_{i,j}
// For a given set of c_{i,j} they are normalised to produce the correct a_c_{i,j}
// Note: Assumes codon frequencies are passed in 64-state space
vector <double> EnforceAAFreqOnCodon(vector<double> CodFreq, vector <double> AAFreq, int GCode)	{
	int i, j;
	int count, aa;
	double Total;
	vector <double> temp;
	// Check entry conditions
	assert(CodFreq.size() == 64); assert(AAFreq.size() == 20);
	assert(fabs(Sum(&CodFreq) - 1.0) < 1.0E-5); assert(fabs(Sum(&AAFreq) - 1.0) < 1.0E-5);
	// Do the codon frequencies on a per amino acid basis
	FOR(aa,20) { // Loop through the amino acids
		// Count the number of codons for a particular amino acid
		Total = 0; count = 0; temp.clear();
		FOR(i,64) { if(GenCodes[GCode][i] == aa) { count++; Total+=CodFreq[i]; temp.push_back(CodFreq[i]); } }
		Total /= AAFreq[aa];
		temp.clear();
		// Non-zero it's simple
		if(Total > 1.0E-5) { FOR(i,64) { if(GenCodes[GCode][i] == aa) { CodFreq[i] /= Total; temp.push_back(CodFreq[i]); } }assert(fabs(Sum(&temp) - AAFreq[aa]) < 1.0E-5); }
		// If it's zero arbitrarily set to even
		else { FOR(i,64) { if(GenCodes[GCode][i] == aa) { CodFreq[i] = AAFreq[aa]/ (double) count; } } }
	}
	// Error checking on exit
	assert(fabs(Sum(&CodFreq) - 1.0) < 1.0E-5);
	return CodFreq;
}

/////////////////////////////////////////////////////////////////////////////////////
// Get the scaling factor between RY recoding and codon models
// ** This function is based on the premise that the RY to nucleotide projection is adequate **
// Similar to the AA -> Codon projection there are two valid ways of producing an adjustment factor
// 1. Take the empirical frequencies, which is the equivalent of taking the MLEs from a multinomial. Adds 2 df
// 2. Equiprobable frequencies, which are always 1/2. Adds 0 df.
double CData::GetRYToCodonlnLScale(int GeneticCode, int *df)	{
	int i, j, count =0;
	double RetValue = 0.0, EmpVal = 0.0, EquVal = 0.0; *df = 0;
	if(m_DataType != DNA) { Error("\nError: CData::GetRYToCodonlnLScale() only works for NUCLEOTIDE data\n"); }
	if(!InRange(GeneticCode,0,11)) { Error("\nTrying to do DNA CData::GetAminoToCodonlnLScale() with out of range GenCode\n"); }
	CData RY_Data = *this; RY_Data.DNA2RY();
	assert(RY_Data.m_iTrueSize == m_iTrueSize && RY_Data.m_iNoSeq == m_iNoSeq);
	assert((int)RY_Data.m_vFreq.size() == 2 && (int)m_vFreq.size() == 4);
	FOR(i,m_iNoSeq) {
		FOR(j,m_iTrueSize) {
			if(m_ariSeq[i][m_ariPatMap[j]] == m_iChar) { continue; } // Skip gaps
/*			cout << "\n["<<i<<"]["<<j<<"]: " <<  "[" <<m_ariSeq[i][m_ariPatMap[j]] << "]=" << m_sABET[m_ariSeq[i][m_ariPatMap[j]]] << " / ";
			cout << "[" << RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]] << "]=" << RY_Data.m_sABET[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]] << " : ";
			cout <<  m_vFreq[m_ariSeq[i][m_ariPatMap[j]]] << " / " << RY_Data.m_vFreq[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]] << " = " << m_vFreq[m_ariSeq[i][m_ariPatMap[j]]] / RY_Data.m_vFreq[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]];
			cout << " --> " << log(m_vFreq[m_ariSeq[i][m_ariPatMap[j]]] / RY_Data.m_vFreq[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]]);
*/			EmpVal += log(m_vFreq[m_ariSeq[i][m_ariPatMap[j]]] / RY_Data.m_vFreq[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]]);
			count++;
		}
	}
	EquVal = count * log(0.5);
	if(EmpVal -2 > EquVal) { RetValue = EmpVal; *df = 2; } else { RetValue = EquVal; }
//	cout << "\n--- NT---\nEmp: " << EmpVal << " (" << EmpVal - 2 << ") and Equ: " << EquVal; if(EmpVal - 2 > EquVal) { cout << "  *EMP*"; } else { cout << " *EQU*"; } cout << endl;
	return RetValue;
}


/* ******************************* NEXUS FILE SUPPORT ********************************* */
void CData::InputNexus(string File) {
    int i,j;
	vector <string> vName, vData,Toks;
	string store;
	EDataType  Type,Old;

	// Check whether a NEXUS file first
    ifstream input(File.c_str());
   	Toks = Tokenise(store);
	if(Toks[0].find("#NEXUS") == string::npos) { cout << "\nIn CData::InputNexus(string file): expecting nexus file in <"<<File<<">\n\n"; Error(); }

	while(!input.eof()) {
		// Some basic initialisation
	   	getline(input,store);
	   	if(input.eof()) { cout << "\nUnexpected end of file..."; Error(); }
	   	Toks = Tokenise(store); if(Toks[0][0] == '[') { continue; }
	   	// Get blocks
	   	if(Toks[0].find("begin") != string::npos) {

	   	}
	}


	input.close();

	cout << "\nSuccessful read...";
	exit(-1);
}
