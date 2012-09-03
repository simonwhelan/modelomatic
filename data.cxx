/////////////////////////////////////////////////
// Data.cxx
////////////////////////////////////////////////


#include "data.h"  // header for class
#include "tree.h"
// Debugging info
#define DEBUG_DATA 0

//////////////////////////////////////////////
// Constructors for data object
/////////////////////////////////////////////////////////////

CData::CData(int NoSeq, int Size, EDataType Type, vector <string> *Names)	{
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
    int i,j, CumSize,LastSize;
	vector <string> vName, vData,Toks;
	vector <int> ariPos;
	vector <int> SiteLabels;
	string store;
	bool First = true, DoSiteLabels = false;
	EDataType  Type,Old;

	// Some original preperation
    m_bValid = false; m_iGenCode = -1;

    ReadData(file,vName,vData);
	/////////////////////////////////////////////////////////////////////////
	// Do a bit more checking and initialise the data object
	// Guess the type of data
	if(SpecType == NONE)	{
		Old = GuessDataType(vData[0]);
		FOR(i,(int)vData.size())	{
			Type = GuessDataType(vData[i]);
			if(Type == NONE || Type != Old) { Type = NONE; break; }
	}	} else { Type = SpecType; }
	// Ensure the data is of a real type
	if(Type == NONE) { Error("\nCouldn't guess data type. Please inspect your data for excesses of weird or gap characters.\nThis message will be triggerred when there are a lot of gaps (>" + double_to_string(100*(1.0-MIN_DATA_PERCENT)) + "%) in a single sequence.\n\n" ); }
    InputData(Type,vData,vName,SiteLabels,AllowFail);		// Put the sequence data into the object
}

// Inputs data from a set of arrays (for pairwise distances)
CData::CData(int NoSeq,int Size, vector <string> InSeq, vector <string> InName,EDataType SpecType) {
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

CData::~CData() { Clean(); }

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
	if(input.eof()) { cout << "\nUnexpected end of file..."; if(AllowFail) { return false; } else { Error(); } }
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
				if(input.eof()) { cout << "\nUnexpected end of file..."; if(AllowFail) { return false; } else { Error(); } }
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
				else		{ if((int)store.size() != LastSize) { cout << "\nData for species["<<i<<"] " << Names[i] << " of unexpected length...\n" << store << "\n\n"; if(AllowFail) { return false; } Error(); } }
				if(CumSize < Size || i < NoSeq - 1) { Toks = Tokenise(GetDataLine(&input)); }
				if(input.eof()) { cout << "\nUnexpected end of file..."; if(AllowFail) { return false; } else { Error(); } }
			}
		}
	}
	input.close();
	return true;
}

// Get Name function
////////////////////////////////////////////////////////////
// Arguements are a char string
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
// TODO: It was written before I used STL -- Needs updating
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
	// Allocate memory and initialise some variables
	m_vFreq.assign(m_iChar,0);
	GET_MEM(ardFreqCount,double,m_iChar+1); FOR(i,m_iChar+1){ ardFreqCount[i] = 0.0; }
	FOR(i,m_iNoSeq)	{ m_vsName.push_back(cInputName[i]); }
	m_ariPatMap.assign(m_iSize,0);
	GET_MEM(ariPat,int*,m_iNoSeq*3); FOR(i,m_iNoSeq) { GET_MEM(ariPat[i],int,m_iSize); }
	GET_MEM(iTempPat,int,m_iNoSeq*3);
	GET_MEM(ariPatOcc,int,m_iSize*3); FOR(i,m_iSize) { ariPatOcc[i] = 0; }
	FOR(i,m_iNoSeq) { m_vsTrueSeq.push_back(cInputSeq[i]); }

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
#ifndef MATCH_PAML
		if(okay < 2) { continue; } // Skips site if only one or zero sequences have a charactor at it
#endif
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
		if(j == m_iNoSeq) { m_viNoChange[i] = m_ariSeq[0][i]; }
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
// Changes a DNA sequence into RY form
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
	// Clean the memory
	Clean();
	// Create the new sequences
	InputData(RY,NewSeq,Names,Labels,false);
}


// Remove sequence functions
int CData::RemoveSeq(int RemSeq,CTree *TREE)	{
	int i =0;
	i+=RemoveSeq(RemSeq);
	if(TREE!=NULL) { i+=TREE->RemoveLeafNode(RemSeq); }
	return i;
}

int CData::RemoveSeq(int RemSeq)	{
/* Removes a sequence from an alignments and adjust all required info
	returns 0 if successful or 1 otherwise	*/
	int i,j;
	double *ardFreqCount, total;
	GET_MEM(ardFreqCount,double,m_iChar);
	// Section to remove a sequence
	for(i=RemSeq;i<m_iNoSeq-1;i++)		{
		m_vsName[i] = m_vsName[i+1];
		for(j=0;j<m_iSize;j++) { m_ariSeq[i][j] = m_ariSeq[i+1][j]; }
		m_vsTrueSeq[i] = m_vsTrueSeq[i+1];
	}
	m_iNoSeq--;
	m_ariSeq.erase(m_ariSeq.end());
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
				if(!InRange(GenCodes[GenCode][FindState(COD,Codon)],0,64)) { Error("\nUnrecognised codon: " + Codon + " for data..."); }
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

double CData::GetAminoToCodonlnLScale(int GeneticCode)	{
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
			 // cout << "\n["<<i<<"]["<<j<<"]: " <<  COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] << " / " << AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]];
			 // cout << " = " << COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]] << " / " << AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]];
			 // cout << " = " << COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]];
			 // cout << " --> " << log(COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]]);
			RetValue += log(COD_Data.m_vFreq[COD_Data.m_ariSeq[i][COD_Data.m_ariPatMap[j]]] / AA_Data.m_vFreq[AA_Data.m_ariSeq[i][AA_Data.m_ariPatMap[j]]]);
		}
	}
	return RetValue;
}

/////////////////////////////////////////////////////////////////////////////////////
// Get the scaling factor between RY recoding and codon models
// ** This function is based on the premise that the RY to nucleotide projection is adequate **
double CData::GetRYToCodonlnLScale(int GeneticCode)	{
	int i, j;
	double RetValue = 0.0;
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
*/			RetValue += log(m_vFreq[m_ariSeq[i][m_ariPatMap[j]]] / RY_Data.m_vFreq[RY_Data.m_ariSeq[i][RY_Data.m_ariPatMap[j]]]);
		}
	}
	return RetValue;
}


