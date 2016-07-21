/*
 * PfamModel.cxx
 *
 *  Created on: May 12, 2016
 *      Author: simon
 */

#include "ding.h"

/*
 *
 *
class CPfamModel : public CBaseModel	{
public:
	CPfamModel(CData *Data, CTree *Tree, string ProbFileName);		// Data/Tree/File name of frequencies
	~CPfamModel();
	// Variables
	vector <vector <double> > m_vvdSiteFreq;						// The storage of the sitewise probability vectors from the ProbFileName (derived from pfam)

	// Implementation
	double lnL(bool ForceReal = false);	// perform a likelihood calculation (if over-ridden, be careful other likelihood functions are too e.g. DoBralnL(...) )
};

class CPfamProcess : public CBaseProcess {
public:
	CPfamProcess(CData *Data, CTree *Tree, vector <vector <double> > *Freqs);		// Data/Tree/Pointer to a permanent store of the frequencies (to be held in the model above)
	~CPfamProcess();
	// Variables
	vector <CQMat *> m_vpSiteQ;	// Array of Q matrices for the process (better because then it doesn't involve rebuilding the matrix each time)

	// Implementation
	bool Likelihood(bool ForceReal = false);		// Main likelihood function
};
 */

#define ALLOW_BRANCH_OPTIMISE 1 // Specifies whether fast branch optimisation is allowed... (Should always be 1!)
#define OUTPUT_BRANCH_OPTIMISE 0

string MainData = "MitoData.phy";
string SmallName = "mito55_small.phy";
string SmallTree = "mito55_small.tre";

enum RunType { LG,PFAM,HET };

bool SAR11Run = true;		// Whether to run the models an SAR11. Will not accept HET model runs
bool DoShort = false;
RunType DoPfamModel = PFAM;
bool DoGamma = true;
bool DoOpt = true;			// Whether to do the optimisation. Used for debug
bool DoH1 = true;
bool DoH2 = true;
bool DoH3 = true;
bool DoOptHere = false;		// Whether to do optimisation in the H1/H2/H3 and general tree functions

double GammaLG[3] = { 0.636298 , 0.636299 , 0.636508 };
double GammaPfam[3] = { 0.605827 , 0.606925 , 0.606983 };
double GammaHet[3] = { 0.633508 , 0.633970 , 0.633941 };

/////////////////////////////////////// Pfam model main code ////////////////////////////////////////////////////////////////
void PfamModelAnalysis()	{
	int i,j;
	double curlnL;
	CBaseModel *Model= NULL;
	string freq_file = "PfamFreqFile.txt", Name1,Name2, Extra;
	vector <int> EukNames;
	vector <int> RikNames;
	vector <int> OthNames;
	vector <int> ProNames;
	vector <int> RootNames;
	vector <vector <int> > SubNames;
	CPfamModel *PfamLG_H1, *PfamLG_H2,*PfamLG_H3;
	CEMP *LG_H1, *LG_H2, *LG_H3;
	CHeteroEMP *LGHet_H1, *LGHet_H2, *LGHet_H3;
	CPfamModel *PfamM = NULL;
	CEMP *LGM = NULL;
	CTree *T = NULL;
	CData *SarDat;


	// Set model header
	switch(DoPfamModel) {
	case LG: Extra = ".LG"; break;
	case PFAM: Extra = ".pfam"; break;
	case HET: Extra = ".het"; break;
	}
	if(DoGamma) { Extra = Extra + ".dG"; }
	cout << "\nDoing " << Extra << " analysis\n------------------------------";
	cout << "\nInitialising data and hypotheses" << flush; cout.precision(10);

	//////////////////////////////// Create smaller data sets ////////////////////////
	if(false)	{
		CData NewDat("mito55_exRog.phy",AA, false,0,true,false);
//		CTree NewTree("mito55_exRog_H1.tre", true, &DingDat);
		int NewNoSeq = 8;		// Number of new sequences
		int NewSize = 1000;	// Length of new sequences
		string NewName = "mito55_small.phy";

		ofstream newout(SmallName.c_str());
		assert(NewNoSeq <= NewDat.m_iNoSeq && NewSize <= NewDat.m_iSize);
		newout << NewNoSeq << "  " << NewSize;
		FOR(i,NewNoSeq) {
			newout << "\n" << NewDat.m_vsName[i] << "\t";
			FOR(j,NewSize) { newout << NewDat.m_vsTrueSeq[i][j]; }
		}
		newout << endl << flush;
		newout.close();
		exit(-1);
	}

	//////////////////////////////// Test development code ///////////////////////////
	if(DoShort) {
		// Data
		cout << "\nReading data" << flush;
		CData TDat(SmallName.c_str(),AA,false,0,true,false);
		cout << "\nData read successfully ("<<TDat.m_iNoSeq << " x " << TDat.m_iSize << ")" << flush;

//		cout << "\nPattern freqs\n---"; FOR(i,TDat.m_iSize) { cout << "\n["<<i<<"] x " << TDat.m_ariPatOcc[i] << "\t"; FOR(j,TDat.m_iNoSeq) { cout << TDat.m_vsTrueSeq[j][i]; } } exit(-1);

		// Tree
		cout << "\nMaking tree" << flush;
//		CTree TT("((1:0.1,2:0.1):0.1,3:0.1,(4:0.1,5:0.1):0.1);",5);
		CTree TT("(((1:0.1,2:0.1):0.1,(3:0.1,4:0.1):0.1):0.1,(5:0.1,6:0.1):0.1,(7:0.1,8:0.1):0.1);",8);
//		CTree TT("(1:10,2:1.273048672,((3:1.448444937,4:10):0,((5:10,6:10):10,(7:10,8:10):10):10):0);",8);
		cout << " ... done" << flush;
		// Model
		cout << "\nMaking model " << flush;
		CBaseModel *Model;
		CPfamModel TPfam(&TDat,&TT,freq_file);
		CEMP LG_Model(&TDat,&TT,"LG",true,(double*)dLGVal,(double*)dLGFreq);		// Set for +F
		switch(DoPfamModel) {
		case LG: Model = &LG_Model; break;
		case PFAM: Model = &TPfam; break;
		case HET: cout << "\nShort version not appropriate with heterogeneous model"; exit(-1);
		}
		cout << "\nAdding Gamma" << flush;
		Model->MakeGammaModel(0,4,0.75);
		cout << " ... done" << flush;

//		Model->SetFastBranchOpt(false);

//		cout << "\nTree: " << TT;

		// Likelihoods
		cout << "\nLikelihood: " << flush;
		curlnL = Model->lnL(true);
		cout << " = " << curlnL << flush;
//		cout << "\nModel: " << *Model; exit(-1);
/*
		cout << "\nSpace check: ";
		FOR(i,3) {
			cout << "\nSite["<<i<<"]: ";
			double *p_test = Model->m_vpProc[0]->ForceRealFd(2,i);
			FOR(j,20) {
				cout << "\t" << *p_test++;
			}
		}
		exit(-1);
*/
/*
		cout << "\nSitewise likelihoods: ";
		Model->m_vpProc[0]->PrepareBraDer();
		FOR(i,100) {
			cout << "\nSite["<<i<<"]: lnL: " << Model->m_arL[i].LogP() << flush;
//			cout << Model->m_vpProc[0]->ModelL(i).Prob() << flush;
		}
		int NFtest = 5;
		FOR(NFtest,8) {
			cout << "\nScales from " << NFtest << " <<<<<<<<<<<<<<<<< REQUIRES UNPROTECTION OF CBASEPROCESS >>>>>>>>>>>>>>>>>>>>>>>>>>";
			FOR(i,100) { if(i%5 == 0)  { cout << endl; } cout << Model->m_vpProc[0]->ForceRealFdSc(NFtest,i) << " = " << *Model->m_vpProc[0]->ForceRealFdSc(NFtest,i) << "  "; }
		}
		exit(-1);
*/
		// Optimisation
		cout << "\nOptimisation of " << Model->Name() << " ... " << flush;
		FullOpt(Model,true, true, false, -BIG_NUMBER,true,DEFAULT_OPTNUM,curlnL,FULL_LIK_ACC,true,true);
		cout << "\nFinished opt... Likelihood: " << flush;
		curlnL = Model->lnL(true);
		cout << " = " << curlnL << flush;
		cout << "\nFinal tree: " << *Model->Tree();
		// Finish
		cout << "\ndone...\n";
		cout << "\nModel: \n" << *Model;
		ofstream outter("ModelDetails.test.txt");
		outter << *Model;
		outter.close();
		Model->OutputSitelnL("SiteLikelihoods.test.txt");

		Model = NULL;
		exit(-1);
	}
	//////////////////////////////// Do SAR 11 run ///////////////////////////////////
	if(SAR11Run)	{
/*
		cout << "\nReached right place";
		double *P,*R;
		int NoC = 10; double Alpha = 0.6;
		GET_MEM(P,double,NoC); GET_MEM(R,double,NoC);
		// Check some entry conditions
		// Get the gamma distributed rates
		DiscreteGamma(P,R,Alpha,Alpha,NoC,0);

		cout << "\nGetting gamma -- alpha ["<<NoC<<"] = " << Alpha << " == { "; FOR(i,NoC) { cout << "\nCat["<<i<<"]: Prob= " << P[i] << ", Rate= " << R[i] << " "; } cout << "}";
		exit(-1);
*/

		int run;
		const int NumTrees = 4;
		char *Trees[] = {"mito55_H1Sr1.tre","mito55_H1Sr2.tre","mito55_H3Sr1.tre","mito55_H3Sr2.tre" };
		string DataFile = "SarData.phy";
		string OutputName;

		CData SDat(DataFile,AA, false,0,true);
		CTree TT(Trees[0],true,&SDat);
		CEMP TM(&SDat,&TT,"LG",true,(double*)dLGVal,(double*)dLGFreq);

		switch(DoPfamModel) {
		case LG:
			SarDat = new CData(DataFile,AA, false,0,true);
			break;
		case HET:
		case PFAM:
			SarDat = new CData(DataFile,AA, false,0,true,false); break;
		}
		cout << "\nData read successfully ("<<SarDat->m_iNoSeq << " x " << SarDat->m_iSize << ")" << flush;
		switch(DoPfamModel) {
		case LG: cout << " ... LG model"; break;
		case HET:
		case PFAM: cout << " ... HET/PFAM model"; break;
		};

		FOR(run,NumTrees) {
			// Set stuff up
			switch(DoPfamModel) {
			case LG:
				T = new CTree(Trees[run],true,SarDat);
				LGM = new CEMP(SarDat,T,"LG",true,(double*)dLGVal,(double*)dLGFreq); Model = LGM;
				break;
			case PFAM:
				T = new CTree(Trees[run],true,SarDat);
				PfamM = new CPfamModel(SarDat,T,freq_file); Model = PfamM;
				break;
			case HET: cout << "\nError: General trees not accepted by heterogeneous model"; exit(-1);
			};
			if(DoGamma) { Model->MakeGammaModel(0,4,0.630); }
			// Do the analysis
			cout << "\n--------------------------------------- Doing analysis: " << Trees[run] << " ---------------------------------------";
			OutputName = Trees[run];
			if(DoOptHere)	{
				cout << "\nStarting optimisation" << flush;
				if(DoOpt) { curlnL = FullOpt(Model,true, true, false, -BIG_NUMBER,true,DEFAULT_OPTNUM,curlnL,FULL_LIK_ACC,true,true); }
		//		cout << "\nModel:\n" << *Model;
			} else {
				curlnL = Model->lnL(true);
			}
			// Output stuff
			cout << "\n >>>>>>>>>>>>>>>>>>>>>>>>> Done. Final likelihood: " << curlnL << " <<<<<<<<<<<<<<<<<<<<<<<<<<" << flush;
			Name1 = "ModelDetails." + OutputName + Extra + ".txt";
			Name2 = "SiteLikelihoods." + OutputName + Extra + ".txt";
			ofstream poutter(Name1.c_str());
			poutter << *Model;
			poutter.close();
			Model->OutputSitelnL(Name2.c_str());
			// Clean up
			Model = NULL;
			if(LGM != NULL) { delete LGM; LGM = NULL; }
			if(PfamM != NULL) { delete PfamM; PfamM = NULL; }
			delete T; T = NULL;
		}
		cout << "\nFINISHED...\n";
		exit(-1);
	}

	// Input
	//   Data
	cout << "\nReading <"<<MainData<<">" << flush;
	CData *DingDat;
//	cout << "\nTest" << flush;
//	CData *DoDah; DingDat = new CData(MainData,AA,false,0,true,false);
//	CData DoDah(MainData,AA,false,0,true,false);
//	cout << "... okay" << flush;
//	exit(-1);

	CTree *T_H1, *T_H2, *T_H3;
	// Trees
	// 	 Tree H1
	switch(DoPfamModel) {
	case LG:
		DingDat = new CData(MainData,AA, false,0,true);
		T_H1 = new CTree("MitoTree.LG.H1.tre",true,DingDat);
		T_H2 = new CTree("MitoTree.LG.H2.tre",true,DingDat);
		T_H3 = new CTree("MitoTree.LG.H3.tre",true,DingDat);
		break;
	case PFAM:
		DingDat = new CData(MainData,AA, false,0,true,false);
		T_H1 = new CTree("MitoTree.pfam.H1.tre",true,DingDat);
		T_H2 = new CTree("MitoTree.pfam.H2.tre",true,DingDat);
		T_H3 = new CTree("MitoTree.pfam.H3.tre",true,DingDat);
		break;
	case HET:
		DingDat = new CData(MainData,AA,false,0,true,false);
		T_H1 = new CTree("MitoTree.het.H1.tre",true,DingDat);
		T_H2 = new CTree("MitoTree.het.H2.tre",true,DingDat);
		T_H3 = new CTree("MitoTree.het.H3.tre",true,DingDat);
		break;
	};
	cout << "\nData read successfully ("<<DingDat->m_iNoSeq << " x " << DingDat->m_iSize << ")" << flush;
/*	cout << "\nPattern freqs\n---";
	FOR(i,DingDat.m_iSize) { cout << "\n["<<i<<"] x " << DingDat.m_ariPatOcc[i] << "\t"; FOR(j,DingDat.m_iNoSeq) { cout << DingDat.m_vsTrueSeq[j][i]; } } exit(-1);
*/

/*	cout << "\nReading trees: H1" << flush;
	CTree OriT_H1("mito55_exRog_H1.tre", true, &DingDat);
	OriT_H1.OutBra(); OriT_H1.OutName();
	CTree T_H1 = OriT_H1;
	//   Tree H2
	cout << " ... H2" << flush;
	CTree OriT_H2("mito55_exRog_H2.tre", true, &DingDat);
	OriT_H2.OutBra(); OriT_H2.OutName();
	CTree T_H2 = OriT_H2;
	//   Tree H3
	cout << " ... H3" << flush;
	CTree OriT_H3("mito55_exRog_H3.tre", true, &DingDat);
	OriT_H3.OutBra(); OriT_H3.OutName();
	CTree T_H3 = OriT_H3;*/
	// Initialise models

	/*
	 * 	// Name files
	 *
	 */

	//////////////////////////// Standard run ////////////////////////////////////////
	switch(DoPfamModel) {
	case LG:
		LG_H1 = new CEMP(DingDat,T_H1,"LG",true,(double*)dLGVal,(double*)dLGFreq);
		LG_H2 = new CEMP(DingDat,T_H2,"LG",true,(double*)dLGVal,(double*)dLGFreq);
		LG_H3 = new CEMP(DingDat,T_H3,"LG",true,(double*)dLGVal,(double*)dLGFreq);
		break;
	case PFAM:
		PfamLG_H1 = new CPfamModel(DingDat,T_H1,freq_file);
		PfamLG_H2 = new CPfamModel(DingDat,T_H2,freq_file);
		PfamLG_H3 = new CPfamModel(DingDat,T_H3,freq_file);
		break;
	case HET:
		// Organise name files
		EukNames = ReadNameFile("mito55_Euk.names",DingDat);
		RikNames = ReadNameFile("mito55_Rik.names",DingDat);
		OthNames = ReadNameFile("mito55_Oth.names",DingDat);
		ProNames = ReadNameFile("mito55_Pro.names",DingDat);
		RootNames = ReadNameFile("mito55_root.names",DingDat);
		assert(EukNames.size() + RikNames.size() + OthNames.size() + ProNames.size() == DingDat->m_vsName.size());
		SubNames.push_back(EukNames); SubNames.push_back(RikNames); SubNames.push_back(OthNames); SubNames.push_back(ProNames);
		vector <int> Temp;
		// Build pairwise sets of these four => 4C2 = 6
		Temp = VecCon(EukNames,RikNames); sort(Temp.begin(),Temp.end());
		vector <int> EukRik = Temp;
		Temp = VecCon(EukNames,OthNames); sort(Temp.begin(),Temp.end());
		vector <int> EukOth = Temp;
		Temp = VecCon(EukNames,ProNames); sort(Temp.begin(),Temp.end());
		vector <int> EukPro = Temp;
		Temp = VecCon(RikNames,OthNames); sort(Temp.begin(),Temp.end());
		vector <int> RikOth = Temp;
		Temp = VecCon(RikNames,ProNames); sort(Temp.begin(),Temp.end());
		vector <int> RikPro = Temp;
		Temp = VecCon(OthNames,ProNames); sort(Temp.begin(),Temp.end());
		vector <int> OthPro = Temp;
		LGHet_H1 = new CHeteroEMP(DingDat,T_H1,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);
		LGHet_H2 = new CHeteroEMP(DingDat,T_H2,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);
		LGHet_H3 = new CHeteroEMP(DingDat,T_H3,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);
		break;

	}

/*
	CHeteroEMP LGHet_H1(&DingDat,&T_H1,"LG_Hetero",(double*)dLGVal,SubNames,RootNames);
	cout << "\nBuilding model from file: <" << freq_file << ">";
	CPfamModel PfamLG_H1(&DingDat,&T_H1,freq_file);
	CEMP LG_H1(&DingDat,&T_H1,"LG",true,(double*)dLGVal,(double*)dLGFreq);
	CPfamModel PfamLG_H2(&DingDat,&T_H2,freq_file);
	CEMP LG_H2(&DingDat,&T_H2,"LG",true,(double*)dLGVal,(double*)dLGFreq);
	CPfamModel PfamLG_H3(&DingDat,&T_H3,freq_file);
	CEMP LG_H3(&DingDat,&T_H3,"LG",true,(double*)dLGVal,(double*)dLGFreq);
	cout << " ... success" << flush;
*/
//	CEMP LG_H1(&DingDat,&T_H1,"LG",false,(double*)dLGVal,(double*)dLGFreq);
//	cout << "\nLG model" << LG_H1.lnL(true); exit(-1);
	if(DoH1) {
		cout << "\n---------------------------------------- Likelihood H1: " << flush;
		switch(DoPfamModel) {
		case LG:
			Model = LG_H1;
			Model->MakeGammaModel(0,4,GammaLG[0]);
		break;
		case PFAM:
			Model = PfamLG_H1;
			Model->MakeGammaModel(0,4,GammaPfam[0]);
			break;
		case HET:
			Model = LGHet_H1;
			Model->MakeGammaModel(0,4,GammaHet[0]);
		break;
		}
		curlnL = Model->lnL(true);
		cout << curlnL << " ------------------------------------------\n";
		if(DoOptHere)	{
			cout << "\nStarting optimisation" << flush;
			if(DoOpt) { curlnL = FullOpt(Model,true, true, false, -BIG_NUMBER,true,DEFAULT_OPTNUM,curlnL,FULL_LIK_ACC,true,true); }
			cout << "\n >>>>>>>>>>>>>>>>>>>>>>>>> Done. Final likelihood: " << curlnL << " <<<<<<<<<<<<<<<<<<<<<<<<<<" << flush;
	//		cout << "\nModel:\n" << *Model;
			Name1 = "ModelDetails.H1" + Extra + ".txt";
			Name2 = "SiteLikelihoods.H1" + Extra + ".txt";
			ofstream outter(Name1.c_str());
			outter << *Model;
			outter.close();
			Model->OutputSitelnL(Name2.c_str());
		}
		delete Model;	// Clears the space and memory
	//	if(!DoOpt) { exit(-1); }
	}
	if(DoH2) {
		cout << "\n---------------------------------------- Likelihood H2: " << flush;
		switch(DoPfamModel) {
		case LG:
			Model = LG_H2;
			Model->MakeGammaModel(0,4,GammaLG[1]);
		break;
		case PFAM:
			Model = PfamLG_H2;
			Model->MakeGammaModel(0,4,GammaPfam[1]);
			break;
		case HET:
			Model = LGHet_H2;
			Model->MakeGammaModel(0,4,GammaHet[1]);
		break;
		}
		curlnL = Model->lnL(true);
		cout << curlnL << " ------------------------------------------\n";
		if(DoOptHere)	{
			if(DoOpt) { curlnL = FullOpt(Model,true, true, false, -BIG_NUMBER,true,DEFAULT_OPTNUM,curlnL,FULL_LIK_ACC,true,true); }
			cout << "\nDone. Final likelihood: " << curlnL;
			cout << "\n >>>>>>>>>>>>>>>>>>>>>>>>> Done. Final likelihood: " << curlnL << " <<<<<<<<<<<<<<<<<<<<<<<<<<" << flush;
	//		cout << "\nModel:\n" << *Model;
			Name1 = "ModelDetails.H2" + Extra + ".txt";
			Name2 = "SiteLikelihoods.H2" + Extra + ".txt";
			ofstream outter(Name1.c_str());
			outter << *Model;
			outter.close();
			Model->OutputSitelnL(Name2.c_str());
		}
		delete Model;
	}
	if(DoH3) {
		cout << "\n---------------------------------------- Likelihood H3: " << flush;
		switch(DoPfamModel) {
		case LG:
			Model = LG_H3;
			Model->MakeGammaModel(0,4,GammaLG[2]);
		break;
		case PFAM:
			Model = PfamLG_H3;
			Model->MakeGammaModel(0,4,GammaPfam[2]);
			break;
		case HET:
			Model = LGHet_H3;
			Model->MakeGammaModel(0,4,GammaHet[2]);
		break;
		}
		curlnL = Model->lnL(true);
		cout << curlnL << " ------------------------------------------\n";
		if(DoOptHere)	{
			if(DoOpt) { curlnL = FullOpt(Model,true, true, false, -BIG_NUMBER,true,DEFAULT_OPTNUM,curlnL,FULL_LIK_ACC,true,true); }
			cout << "\nDone. Final likelihood: " << curlnL;
			cout << "\n >>>>>>>>>>>>>>>>>>>>>>>>> Done. Final likelihood: " << curlnL << " <<<<<<<<<<<<<<<<<<<<<<<<<<" << flush;
	//		cout << "\nModel:\n" << *Model;
			Name1 = "ModelDetails.H3" + Extra + ".txt";
			Name2 = "SiteLikelihoods.H3" + Extra + ".txt";
			ofstream outter(Name1.c_str());
			outter << *Model;
			outter.close();
			Model->OutputSitelnL(Name2.c_str());
		}
		delete Model;
	}
	delete T_H1, T_H2, T_H3, DingDat;
	cout << "\nDone...\n" << flush;

}

//////////////////////////////////////////// Definition of the Pfam model class ////////////////////////////////////////////////////////////////////////////////////////

CPfamModel::CPfamModel(CData *D, CTree *T, string File) : CBaseModel(D,T)	{
	int i,j;
	string store;
	vector <string> Toks;
	vector <double> OutFreq;
	CPfamProcess *TempProc;
	m_sName = "PfamModel";
	// Check tree
	if(T->IsRooted()) { cout << "\nTree cannot be rooted for Pfam Model" << flush; exit(-1); }
	// Some temporary code to output some frequencies //
	if(false) 	{
		ofstream out(File.c_str());
		if(false) { OutFreq = D->m_vFreq; FOR(i,D->m_iSize) { out << OutFreq << "\n"; } } // Use same freqs per site (identical to LG+dG
		else { // Build per site freqs
			FOR(j,D->m_iSize) {
				double c = 0;
				OutFreq.clear(); FOR(i,20) { OutFreq.push_back(c * D->m_vFreq[i]); }		// La Place counts of c * D->m_vFreq;
				FOR(i,D->m_iNoSeq) {
					if(D->m_ariSeq[i][j] == 20) { continue; }
					OutFreq[D->m_ariSeq[i][j]]+= 1.0;
				}
				OutFreq = NormaliseVector(OutFreq);

//				OutFreq.clear(); FOR(i,20) { OutFreq.push_back(dLGFreq[i]); }	// Outputs the LG freq
				out << OutFreq << "\n";
			}
		}
		out.close();
		cout << "\nOutput some freqs to <"<<File << ">";
		exit(-1);
	}
	// Read the file with the frequencies
	ifstream in(File.c_str()); assert(!in.eof());
	FOR(i,D->m_iSize)	{
		getline(in,store);	// Can add delim here
		assert(!in.eof());
		Toks = Tokenise(store); assert(Toks.size()==20);
		OutFreq.clear(); FOR(j,20) { OutFreq.push_back(atof(Toks[j].c_str())); }
		if(fabs(Sum(&OutFreq) - 1) > 1.0E-4)  { cout << "\nError: Frequency line ["<<i<<"] sums to " << Sum(&OutFreq) << "\n" << Toks; }
		m_vvdSiteFreq.push_back(OutFreq);
	}
	// Initialise the model processes
	TempProc = new CPfamProcess(D,T,&m_vvdSiteFreq);
	m_vpProc.push_back(TempProc);
	TempProc = NULL;
	// Final initialisation
//	SetFastBranchOpt(false);		// Initially won't use fast branch calculations
	FinalInitialisation();
}

CPfamModel::~CPfamModel()	{

}

bool BoolAllowBranchAnalytics = true;

/////////////////////////////////////////////////////////////////////////
// Derivative calculations
vector <double> CPfamModel::GetDerivatives(double CurlnL, bool *pOK)	{
	int i;
	bool OK = true, ForceNumBra = false;
	vector <double> Grads;
	vector <double> temp;
	double temp_lnL;
	if(m_vpAllOptPar.empty())  { return Grads; }
	////////////////////////// Do likelihood derivatives /////////////////
	if(ALLOW_ANALYTICAL_BRANCH_DERIVATIVES && BoolAllowBranchAnalytics)	{
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
			// PreparelnL(); // No PreparelnL() in this function since it's done sitewise in the individual processes
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
	}	} else {
		FOR(i,(int)m_vpAllOptPar.size())	{ m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL)); }
	}
	// Store the gradients
	FOR(i,(int)m_vpAllOptPar.size()) { Grads.push_back(m_vpAllOptPar[i]->grad()); }
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

double CPfamModel::lnL(bool ForceReal){
	int i;
	double logL = 0.0, TotalProb = 0.0;
	static int lnLCount = 0;
//	cout << "\n\n========== Starting new likelihood calc "  << lnLCount << " ================ "<< flush;
	lnLCount ++;
	// New stuff dealing with submodels (Left in, although I'm not sure how useful it is)
	if(!m_vpAssociatedModels.empty()) {
		FOR(i,(int) m_vpAssociatedModels.size()) {
			// Fix cases where I can assign a tree ok
			if(m_vpAssociatedModels[i] == NULL) { Error("Associated model NULL...\nProbably a problem with TreeHMMs...\n"); }
			logL += m_vpAssociatedModels[i]->lnL(ForceReal);
	}	}
	if(!m_pData->Valid()) { Error("\nTrying to do likelihood computation with invalid data\n\n"); }
	assert(m_pTree != NULL && m_pData != NULL);
	if(m_pData->m_iNoSeq != m_pTree->NoSeq()) { Error("\nMismatch between number of sequences in data and tree\n\n"); }
	// Prepare the Q mats
	// PreparelnL(); // No PreparelnL() in this function since it's done sitewise in the individual processes
	// Get all of the processes ready individually and calculate likelihoods
	FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->Likelihood(ForceReal);
	}
	// Get probabilities of the processes; TODO: This is currently really basic
	if(!FormMixtureSitewiseL()) { cout << "\nReturning error..."; return -BIG_NUMBER; }
	// Get the log likelihood and return it
	FOR(i,m_pData->m_iSize) {
			logL += m_pData->m_ariPatOcc[i] * m_arL[i].LogP();
	}
	// Apply extra stuff to the likelihood function if needed
	if(pLikelihood != NULL) {
		logL -= pLikelihood(NULL); // Called as blank. Other arguments intended to allow functionality
	}
//	cout << "\n\tReturning likelihood: " << logL << endl;  // exit(-1);
	return logL;
}

//////////////////////////////////////////// Definition of the Pfam process class ////////////////////////////////////////////////////////////////////////////////////////


CPfamProcess::CPfamProcess(CData *D, CTree *T, vector <vector <double> > *Freqs) : CBaseProcess(D,T) {
	int site,i,j, count;
	string Name;
	CQPar *Par = NULL;
	CQMat *tQMat = NULL;
	CBaseEqm *tEqm = NULL;
	double *S_ij = (double*) dLGVal,*pdQMat;
	// Some basic error checking on input
	assert(D->m_iSize == D->m_iTrueSize && D->m_iSize == (int) Freqs->size());
	m_sName = "PfamProcess";
	// Make the space
	MakeBasicSpace(AA);
	// Make the S_ij (assuming LG)
	// Correct S_ij to ensure that none are over the maximum value allowed by a parameter
	count = 0;
	FOR(i,m_iChar) {
		FOR(j,i)	{
			Name = GetPos(m_sABET,i,m_iABET_length/m_iChar) + " <-> " + GetPos(m_sABET,j,m_iABET_length/m_iChar);
			Par = new CQPar(Name,m_iChar,S_ij[count++],false,0.0,BIG_NUMBER);
			Par->AddQij(i,j);
			Par->SetOptimise(false);
			m_vpPar.push_back(Par);
			Par = NULL;
	}	}
//	cout << "\nBuilt S_ij and now working on QMats. The parameters are:"; FOR(i,m_vpPar.size()) { cout << "\nPar["<<i<<"]" << *m_vpPar[i]; }
	// Make QMats - one for each site
	FOR(site,D->m_iSize) {
		Name = "Qmat[" + int_to_string(site) + "]";
		tQMat = Add_QMat(Name,AA); // Grab a pointer for building the matrix
		tQMat->InitQ(m_dBaseVal);
		FOR(i,m_vpPar.size()) { m_vpPar[i]->UpdatePar(); tQMat->ApplyPar2Q(m_vpPar[i]); }
		// Apply the eqm here
		pdQMat = tQMat->Q();
		FOR(i,m_iChar) { // By *row*
			FOR(j,m_iChar) { // By column
				if(i==j) { continue; }
				pdQMat[(i*m_iChar)+j] *= Freqs->at(site)[j];
			}
		}
		pdQMat = NULL;
		tQMat->DoQDiag();
		tQMat->Decompose(Freqs->at(site),true,true,m_pRate->Val());
		tQMat->Lock();
		m_vpSiteQ.push_back(tQMat);
		//		cout << "\nThe QMat is: " << *tQMat;
		//		exit(-1);
		tQMat = NULL;	// Clear the pointer
	}
	// Clean up
	S_ij = NULL;
	assert(m_vpQMat.size() == m_vpSiteQ.size()); FOR(i,m_vpQMat.size()) { m_vpQMat[i] = NULL; } m_vpQMat.clear(); // Remove the m_vpQMat stuff since this is now stored in m_vpSiteQ.
}

CPfamProcess::~CPfamProcess() {
	int i;
	if(!m_vpSiteQ.empty()) { FOR(i,m_vpSiteQ.size()) { delete m_vpSiteQ[i]; m_vpSiteQ[i] = NULL; } m_vpSiteQ.clear(); }
}

// Functions dealing with rate process copies
CBaseProcess * CPfamProcess::RateProcessCopy()	{
	CPfamProcess *NewProc;
	static int CopyNum = 0;
	string Name;
	int i;
	vector <vector <double> > Freqs;
	FOR(i,m_vpSiteQ.size()) { Freqs.push_back(m_vpSiteQ[i]->Eqm()); }
	// Do naming and initialise
	Name = "Pseudo" + m_sName + "(CopyID=" + int_to_string(CopyNum++) + ")";
	NewProc = new CPfamProcess(m_pData, m_pTree,&Freqs);
	NewProc->m_sName = Name;
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

bool CPfamProcess::Likelihood(bool ForceReal) {
	bool RetVal = true;
	int i,j, site;
	static int Counter =0;
//	cout << "\n----------------------------------------------------------------------------------\nPfam process likelihood comp: " << Counter++;
//	cout << "\n--------------------------------------------------\nLikelihood\nTree: " << *m_pTree << flush;

	// Some basic entry functions
	if(m_vSpace.empty()) { MakeCalcSpace(false); }
	// Otherwise do the usual computations
	bool OldComp = m_bCompressedSpace;
	m_bBraDerReady = false;
	m_bFailedL = false;
	assert(m_bDoingPartial == false);
	// Apply the parameters (necessary for the gamma distribution)
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->GlobalApply(); }
	// If the probability of the process is small then stop here
	if(Prob() < SMALL_PROB) { return true; }
	// Preprare Q matrices for all the sites
	FOR(i,(int)m_vpSiteQ.size()) {	m_vpSiteQ[i]->ScaleQ(m_pRate->Val()); }
	// For each site replace m_vpQMat (standard likelihood stuff) with an element from m_vpSiteQ. The calculations are then pretty simple
	FOR(site,m_pData->m_iSize)	{
//		cout << "\nDoing likelihood for site["<<site<<"]" << flush;
		// Identify to partial likelihood functions what sites being calculated
		m_iCurSite = site;
		// Put the current sites matrix into the appropriate pointer
		assert(m_vpQMat.empty()); m_vpQMat.push_back(m_vpSiteQ[site]);
		// Prepare the P(t) matrices for each site
		FOR(i,Tree()->NoBra()) { Tree()->UpdateB(i); Make_PT(i); }
//		cout << "\nP(t)[0]"; OutPT(cout,0); cout << " done" << flush;
		// Adjust m_bCompressedSpace if required
		if(ForceReal == true) { m_bCompressedSpace = false; }
		// Call the local PartialL functions to do the calculation for a single site
		CleanScale(-1,ForceReal);	// Clean the scaling array in the storage space
		PartialL(Tree(),Tree()->StartCalc(),-1,-1,true);
		// Store the resultant likelihood
//		cout << "\nReady to assign" << flush;
		m_ardL[site].Assign(Lsum(site));
//		cout << "\nAssigned likelihood: " << m_ardL[site] << " : " << m_ardL[site].LogP() << flush;
		// Clean up
		m_vpQMat[0] = NULL; m_vpQMat.clear(); m_iCurSite = -1;

//		if(site > 10) { exit(-1); }
	}
//	cout << "\nHave likelihood";
	// Get the final likelihoods for the process
//	FOR(i,10) { cout << "\nSite["<< i<< "]: " << m_ardL[i].LogP(); }
	// Return if okay
	m_bCompressedSpace = OldComp;
	return RetVal;
}

///////////////////////////////// Full likelihood functions /////////////////////////////////

void CPfamProcess::MakeZeroRateLikelihood()	{
	int i;
	vector <double> eqm;
	FOR(i,m_pData->m_iSize)	{
		eqm = m_vpSiteQ[i]->Eqm();
		if(m_pData->m_viNoChange[i] == -1) { m_ardL[i].Zero(); continue; }
		if(!InRange(m_pData->m_viNoChange[i],0,(int)eqm.size())) {
			cout << "\nGoing to throw assert error... Site["<<i<<"] which is: ";
			int j; FOR(j,m_pData->m_iNoSeq) { cout << m_pData->m_ariSeq[j][i] << " "; }
			cout << "\nAnd m_pData->m_viNoChange[" << i<<"]: " << m_pData->m_viNoChange[i];
			cout << "\nEqm[size=" << eqm.size() << "]: " << m_vdEqm;
			cout << "\nAnd the vector of m_viNoChange:\n" << m_pData->m_viNoChange;
		}
		assert( InRange(m_pData->m_viNoChange[i],0,(int)eqm.size()) );
		m_ardL[i].Assign(eqm[m_pData->m_viNoChange[i]]);
	}
}

void CPfamProcess:: MakeMaxRateLikelihood() {
	int i,j;
	vector <double> eqm;
	FOR(i,m_pData->m_iSize) {
		m_ardL[i].Assign(1.0);
		eqm = m_vpSiteQ[i]->Eqm();
		FOR(j,NoSeq()) {
			if(m_pData->m_ariSeq[j][i] == m_iChar) { continue; }
			m_ardL[i].Multiply(eqm[m_pData->m_ariSeq[j][i]],true);
}	}	}


// Calculates partial likelihoods recursively
// --
// Note BlockWriteback is used for PreparePartialL where you don't what to write to iNoFr for the starting node only
void CPfamProcess::PartialL(CTree *pTree, int iNoTo, int iNoFr, int Branch, bool NodeFirst, bool BlockWriteback)	{
    int i,j,NodePos1, NodePos2;
	bool First = true;
	double *I = NULL, *p_a = NULL, *p_b = NULL;
//	cout << "\nPartialL Site: " << m_iCurSite <<" - NT: " << iNoTo << ", iNoFr" << iNoFr << ", Branch: " << Branch << ", NodeFirst: " << NodeFirst << flush;
	// Entry
	assert(InRange(m_iCurSite,0,m_pData->m_iSize));
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
			NodePos1 = InitNodePos(iNoTo) + m_iCurSite;								// Space from
			NodePos2 = InitNodePos(PartLNode()) + m_iCurSite;						// Space to (the storage node)
			if(m_bCompressedSpace) {
				i = m_iCurSite;
				if(QkFdReal(NodePos2))	{
					p_a = QkFd(NodePos2); p_b = QkFd(NodePos1);
					FOR(j,m_iChar) { *(p_a++) *= *(p_b++); }
					// do the scale
					*LScale(i) += *QkFdSc(NodePos1);
					NodePos1; NodePos2;
			}	} else { // Uncompressed space
				i = m_iCurSite;
				p_a = QkForceRealFd(NodePos2); p_b = QkForceRealFd(NodePos1);
				FOR(j,m_iChar) { *(p_a++) *= *(p_b++); }
				// do the scale
				*LScale(i) += *QkForceRealFdSc(NodePos1);
				NodePos1++; NodePos2++;
			}
//			cout << "\nAnswer:   "; FOR(i,20) { cout << " " << ForceRealFd(PartLNode(),0)[i]; }
			p_a = NULL; p_b = NULL;
	}	} else if(Tree()->NoSeq() == 2) {	// Do the final calculation for 2 species trees
		if(!m_viLeafMap.empty())	{ assert((int)m_viLeafMap.size() >= iNoTo); iNoTo = m_viLeafMap[iNoTo]; }
		NodePos1 = InitNodePos(PartLNode()) + m_iCurSite;
		i = m_iCurSite;
		Vec_by_Data(m_pData->m_ariSeq[iNoTo][i],QkFd(NodePos1++));
	}
}
// Simple calculations for BranchNodes
void CPfamProcess::BranchNodePartialL(int SpFrom, int SpTo, double *PT, bool First)	{
	int i,k,NodePos1,NodePos2;
    double *p_a,*p_b, *p_c;
	// Entry
	assert(InRange(m_iCurSite,0,m_pData->m_iSize));
//	cout << "\n\tBranchNode : " << m_iCurSite << " - " << SpFrom << " ; " << SpTo << flush;
	// If the node comes from last link in StartCalc() then copy information to PartLNode()
	if(SpTo == Tree()->StartCalc() && !m_bDoingPartial) {
		FOR(k,Tree()->NoLinks(SpTo)) { if(Tree()->NodeLink(SpTo,k) == SpFrom) { break; } }
		assert(k != Tree()->NoLinks(SpTo));
		switch(k)	{	// Copy stuff to the internal node if required
		case 1:
			NodePos1 = InitNodePos(SpTo) + m_iCurSite; NodePos2 = InitNodePos(PartLNode()) + m_iCurSite;
			i = m_iCurSite;
			m_vBackSp[NodePos1++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));
			break;
		case 2:
			NodePos1 = InitNodePos(SpTo) + m_iCurSite; NodePos2 = InitNodePos(PartLNode()) + m_iCurSite;
			i = m_iCurSite;
			m_vSpace[NodePos1++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));
			break;
		};
		SpTo = PartLNode();
	}
	// Initialise quick space access
	NodePos1 = InitNodePos(SpFrom) + m_iCurSite;
	NodePos2 = InitNodePos(SpTo) + m_iCurSite;
	if(First == true)	{ // Copies the probabilities
		// Do first round of calculations for the node
		if(m_bCompressedSpace)	{
			k = m_iCurSite;
			if(!QkFdReal(NodePos2)) { NodePos1++; NodePos2++; }
			else { VMat(QkFd(NodePos1++),PT,QkFd(NodePos2++),m_iChar); }
		} else {	// Incompressed space
			k = m_iCurSite;
			VMat(QkForceRealFd(NodePos1++),PT,QkForceRealFd(NodePos2++),m_iChar);
	}	} else { // Replaces the probabilities (This has been hard checked)
		if(m_bCompressedSpace)	{
			// Do the final round of calculations, including the backspace storage
			k = m_iCurSite;
				// Deal with when the node isn't real
				if(!QkFdReal(NodePos2)) { NodePos1++; NodePos2++; }
				else{
					p_a = QkFd(NodePos1++); p_b = QkBk(NodePos2);
					VMat(p_a,PT,p_b,m_iChar);
					p_c = QkFd(NodePos2++);
					FOR(i,m_iChar) { *(p_c++) *= *(p_b++); }
				}
		} else {	// Uncompressed space
			k = m_iCurSite;
			p_a = QkForceRealFd(NodePos1++); p_b = QkForceRealBk(NodePos2);
			VMat(p_a,PT,p_b,m_iChar);
			p_c = QkForceRealFd(NodePos2++);
			FOR(i,m_iChar) { *(p_c++) *= *(p_b++); }
	}	}
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

void CPfamProcess::LeafNodePartialL(CTree *pTree, int LeafNode, int Branch, int Sp, double *PT, bool First)	{
	int site,i,k,NodePos, NodePos2;
	static double TempSp[MAX_SPACE],*p_a, *p_b, *p_c;
	vector <double> teqm;
//	cout << "\n\tLeafNode : " << m_iCurSite << " - " << LeafNode<< " ; " << Sp << " & PT: " << flush; FOR(i,10) { cout << "  " << PT[i]; }
	// If comes from StartCalc(), then transfer info to final node
	// If the node comes from last link in StartCalc() then copy information to PartLNode()
	if(Sp == Tree()->StartCalc() && !m_bDoingPartial) {
		if(Tree()->NoSeq() > 2) {
			FOR(k,Tree()->NoLinks(Sp)) { if(Tree()->NodeLink(Sp,k) == LeafNode) { break; } }
			assert(k != Tree()->NoLinks(Sp));
			switch(k)	{	// Copy stuff to the internal node if required
			case 1:
				NodePos = InitNodePos(Sp) + m_iCurSite; NodePos2 = InitNodePos(PartLNode()) + m_iCurSite;
				i = m_iCurSite;
				m_vBackSp[NodePos++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));
				break;
			case 2:
				NodePos = InitNodePos(Sp) + m_iCurSite; NodePos2 = InitNodePos(PartLNode()) + m_iCurSite;
				i = m_iCurSite;
				m_vSpace[NodePos++].CopyVals(&m_vSpace[NodePos2++],FlipBool(m_bCompressedSpace));
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
		NodePos = InitNodePos(Sp) + m_iCurSite;
	} else		{		// For final likelihood space
		NodePos = InitNodePos(PartLNode()) + m_iCurSite;
	}
	// Loop through sites
	// If TempSp copies to the space
	if(First == true)	{
		if(m_bCompressedSpace)	{
			site = m_iCurSite;
			if(!QkFdReal(NodePos)) { NodePos++; }							// Skip if required
			else {
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);			// Get data vector
				p_a = QkFd(NodePos++); p_b = TempSp;									// Set Pointers
				// *COPY* TempSp to the calculation space; BkSp not needed
				FOR(i,m_iChar) { *(p_a++) = *(p_b++); }
			}
		} else {
			site = m_iCurSite;
			Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);			// Get data vector
			p_a = QkForceRealFd(NodePos++); p_b = TempSp;									// Set Pointers
			// *COPY* TempSp to the calculation space; BkSp not needed
			FOR(i,m_iChar) { *(p_a++) = *(p_b++); }
	}	}	else {
	// Check entry conditions
	// If TempSp multiplies with the space
		if(m_bCompressedSpace)	{
			site = m_iCurSite;
			if(!QkFdReal(NodePos))	{ NodePos++; }					// Skip if required
			else {
				Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);		// Get data vector
				p_a = QkFd(NodePos); p_b = TempSp; p_c = QkBk(NodePos++);			// Set Pointers
				// *MULTIPLY* TempSp and calculation space; BkSp has a copy
				FOR(i,m_iChar) {
					*(p_c++) = *(p_b);			// Copy to BkSp
					*(p_a++) *= *(p_b++);		// Multiply by SpaceTo
		}	}	} else {
			site = m_iCurSite;
			Data2PartL(m_pData->m_ariSeq[LeafNode][site],PT,TempSp,&teqm);		// Get data vector
			p_a = QkForceRealFd(NodePos); p_b = TempSp; p_c = QkForceRealBk(NodePos++);			// Set Pointers
			// *MULTIPLY* TempSp and calculation space; BkSp has a copy
			FOR(i,m_iChar) {
				*(p_c++) = *(p_b);			// Copy to BkSp
				*(p_a++) *= *(p_b++);		// Multiply by SpaceTo
	}	}	}
	p_a = NULL; p_b = NULL; p_c = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
// Space access functions
void CPfamProcess::TransScale(int NT, int NF,bool First, bool Partial,bool ForceReal)	{
#if ALLOW_SCALE == 1
	int i = m_iCurSite, NPosT, NPosF;
	assert(InRange(NT,-1,m_pTree->NoNode()+1) && InRange(NF,0,m_pTree->NoNode()));
	// if NT >= 0 then its a normal tranfer; otherwise transfer to PartLNode
	if(NT >= 0)		{ NPosT = InitNodePos(NT) + m_iCurSite; } else { NPosT = InitNodePos(PartLNode()) + m_iCurSite; assert(First == true); }
	// Set the other pointer
	NPosF = InitNodePos(NF) + m_iCurSite;
	// Do the normal scale
	if(ForceReal == true || m_bCompressedSpace == false)	{
		if(Partial == true)	{
			if(First == true)	{
				*QkForceRealFdSc(NPosT++) = *QkForceRealFdSc(NPosF++);
			}	else {
				*QkForceRealBkSc(NPosT) = *QkForceRealFdSc(NPosF);
				*QkForceRealFdSc(NPosT) += *QkForceRealBkSc(NPosT);
				NPosF++; NPosT++;
		}	} else {
			if(First == true)	{
				*QkForceRealFdSc(NPosT++) = *QkForceRealFdSc(NPosF++);
			} else {
				*QkForceRealFdSc(NPosT++) += *QkForceRealFdSc(NPosF++);
	}	}	} else {
		// Do the column sorted scale
		if(Partial == true)	{
			if(First == true)	{
				if(!QkFdReal(NPosT)) { NPosT++; NPosF++; }
				else { *QkFdSc(NPosT++) = *QkFdSc(NPosF++); }
		} else {
				if(!QkFdReal(NPosT)) { NPosT++; NPosF++; }
				else {
					*QkBkSc(NPosT) = *QkFdSc(NPosF);
					*QkFdSc(NPosT) += *QkBkSc(NPosT);
					NPosF++; NPosT++;
		}	}	} else {
			if(First == true)	{
				if(!QkFdReal(NPosT)) { NPosT++; NPosF++; }
				else { *QkFdSc(NPosT++) = *QkFdSc(NPosF++); }
			} else {
					if(!QkFdReal(NPosT)) { NPosT++; NPosF++; }
					else { *QkFdSc(NPosT++) += *QkFdSc(NPosF++); }
	}	}	}
#endif
}
void CPfamProcess::DoScale(int Node, bool ForceRealScale)	{
#if ALLOW_SCALE == 1
//	cout << "\nTrying Scale: node= "<< Node << " pos= "<< m_iCurSite;
	int j, NodePos = InitNodePos(Node) + m_iCurSite;
	double BkMax,FdMax, *p_f = NULL, *p_b = NULL;
	if(Node == m_pTree->NoNode()) { return; }	// Don't scale if its from the last node
	assert(Node >= 0 && Node < m_pTree->NoNode());
	if(!QkFdReal(NodePos) && ForceRealScale == false && m_bCompressedSpace == true) { return; }
	// Find Max at site
	p_f = QkForceRealFd(NodePos); p_b = QkForceRealBk(NodePos);
	BkMax = FdMax = 0.0; FOR(j,m_iChar)	{
		if(*(p_f) > FdMax) { FdMax = *(p_f); }
		if(*(p_b) > BkMax) { BkMax = *(p_b); }
		p_f++; p_b++;
	}
	// If max are small enough then scale
	if(BkMax < P_SCALE_VAL)	{
//		cout << " ... BkMax scale -> Address: " << QkForceRealBkSc(NodePos) << " = " << QkForceRealBkSc(NodePos)[0];
		// Allow true zero
		if(Double_Zero(BkMax) || *QkForceRealBkSc(NodePos) > 1000000000) { p_b = QkForceRealBk(NodePos); FOR(j,m_iChar) { *(p_b++) = 0.0; *QkForceRealBkSc(NodePos) = 0; } }
		else {
			// Do normal scale
			while(BkMax < 1.0)	{ BkMax *= 10; p_b = QkForceRealBk(NodePos); FOR(j,m_iChar)	{ *(p_b++) *= 10; } QkForceRealBkSc(NodePos)[0]++; }
			while(BkMax > 10.0)	{ BkMax /= 10; p_b = QkForceRealBk(NodePos); FOR(j,m_iChar)	{ *(p_b++) /= 10; } QkForceRealBkSc(NodePos)[0]--; }
		}
//		cout << " to " << QkForceRealBkSc(NodePos)[0];
	}
	if(FdMax < P_SCALE_VAL)	{
//		cout << " ... FdMax scale -> Address: " << QkForceRealFdSc(NodePos) << " = " << QkForceRealFdSc(NodePos)[0];
		// Allow true zero
		if(Double_Zero(FdMax)  || *QkForceRealFdSc(NodePos) > 1000000000) { FOR(j,m_iChar) { QkForceRealFd(NodePos)[j] = 0.0; } *QkForceRealFdSc(NodePos) = 0;  }
		else {
			// Do normal scale
			while(FdMax < 1.0)	{ FdMax *= 10; p_f = QkForceRealFd(NodePos); FOR(j,m_iChar)	{ *(p_f++) *= 10; } QkForceRealFdSc(NodePos)[0]++; }
			while(FdMax > 10.0)	{ FdMax /= 10; p_f = QkForceRealFd(NodePos); FOR(j,m_iChar)	{ *(p_f++) /= 10; } QkForceRealFdSc(NodePos)[0]--; }
		}
//		cout << " to " << QkForceRealFdSc(NodePos)[0];
	}
	p_f = NULL; p_b = NULL;
#endif
}

void CPfamProcess::CleanScale(int NodeNum, bool ForceAll)	{
	int i, NodePos;
	assert(m_bCheckSpace == true);
	m_bBraDerReady = false;
	assert(InRange(NodeNum,-1,m_pTree->NoNode()) || (NodeNum == 2 && m_pData->m_iNoSeq == 2));
	if(NodeNum == -1) { // Clean the base scale
		*LScale(m_iCurSite) = 0;
	} else {			// Clean the node scale
		NodePos = InitNodePos(NodeNum) + m_iCurSite;
		if(ForceAll == false)	{ *QkFdSc(NodePos) = 0; *QkBkSc(NodePos++) = 0; }
		else					{ *QkForceRealFdSc(NodePos) = 0; *QkForceRealBkSc(NodePos++) = 0; }
}	}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Analytic branch derivative functions
///////////////////////////////////////////////////////////////////////////////////////////////////////

// Public function used to prepare calculations
// CPfam -- Changed so QP(t) matrices are not prepared before hand.
void CPfamProcess::PrepareBraDer()	{
	int i;
	// Get the Q.P(t) matrices required for analytic derivatives
	if(m_pRate->Val() > DX)	{
		Likelihood(true);	// this likelihood function is critical and I assume that all the preparation for the P(t) matrices is done correctly here
//		FOR(i,Tree()->NoBra()) {
//			if(Tree()->GoodBra(i) == true) { MulMat(m_vpQMat[QMat4Bra(i)]->Q(),PT(i),QP(i),m_iChar,m_iChar,m_iChar); }
//		}
}
	// NB: For small rates there is no point getting QP because it is zero
	m_bBraDerReady = true;
}

// Public function to get branch derivatives
// CPfam -- Currently same as CBaseProcess
bool CPfamProcess::GetBraDer(CProb *ModelL)     {
	int i, BrError = 1;
//		cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Entered GetBraDer, m_vpQMat: " << m_vpQMat.size() << "\n" << flush;
	// Initialise
	m_arModelL = ModelL;

//	cout << "\nSitewise likelihoods: "; FOR(i,m_pData->m_iSize) { cout << "\nSite["<<i<<"]: lnL: " << m_arModelL[i].LogP() << flush; } exit(-1);

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
//	cout << "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Done GetBraDer\n" << flush;
//	exit(-1);
//	cout << "\nBranch derivatives: "; FOR(i,Tree()->NoBra()) { cout << Tree()->pBra(i)->grad() << "\t"; } exit(-1);
	return FlipBool(m_bFailedL);
}


// In-order tree traversal to calculate the gradients
// CPfam -- Currently same as CBaseProcess
void CPfamProcess::Branch_dT(int NTo, int NFr, int Branch, CTree *pTree, int First, int *BrError)	{
	int i;
//	cout << "\nBranch_dT(NTo= " << NTo << ", NFr= " << NFr << ", Br= " << Branch << "): m_vpQMat = " << m_vpQMat.size();
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

// For the following functions CPfamProcess calculates QP(t) on the fly (otherwise there's just too much memory taken up)

void CPfamProcess::LeafNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError)	{
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n-------------- LeafNode_dT: NTo: " << NTo << "; NFr: " << NFr << "; Br: " << Br << " -------------";
#endif
	bool CharCheck,LeafSeq = true;
	int i,site,Seq, SiteScale = 0;	// Counters
	double Total, Value, *p_a = NULL, *p_b = NULL;
	static double Vec[MAX_CHAR];
	CProb NewL;

//	cout << "\nScales:";FOR(site,m_iSize) { if(site %20 == 0)  { cout << endl; } cout << *ForceRealFdSc(NFr,site) << "  "; }

	// Variables refreshed at every site
	vector <double> eqm;
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
		// Put the QMatrix in the m_vpQMat spot
		m_iCurSite = site;
		assert(m_vpQMat.empty());
		m_vpQMat.push_back(m_vpSiteQ[site]);

		// Make the P(t) matrix for the branch
		m_vpSiteQ[site]->MakePT(Tree()->B(Br),PT(Br));
		// Make the QP(t) matrix
		MulMat(m_vpSiteQ[site]->Q(),PT(Br),QP(Br),m_iChar,m_iChar,m_iChar);
		// Build the QP(t) matrix for the site
		eqm = m_vpQMat[QMat4Bra(Br)]->Eqm();

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
//			if(site < BIG_NUMBER) { cout << "\nSite["<<site<<"]: Value: " << Value << " += " << Total << " ModelL: "<< m_arModelL[site].Prob() << " -> " << m_arModelL[site].LogP() << " scale: " << SiteScale << " - " << ModelL(site).Scale() << " = " << SiteScale - ModelL(site).Scale() << " == " << PartialGrad(site,Total,SiteScale - ModelL(site).Scale()); }
			Value += PartialGrad(site,Total,SiteScale - ModelL(site).Scale());
#if ANALYTIC_DERIVATIVE_DEBUG == 1
			if(site < ADD_SITE_MAX) { cout << "\nTotal = " << Value; }
#endif
		}
		// update the space for this site only
		if(First != 1) { LeafNode_Update(NTo,NFr,Br,pTree,First); }
		// Clean up
		m_vpQMat[0] = NULL; m_vpQMat.clear();
	}
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n\tBranch[" << Br <<"] = " << Tree()->B(Br) << ": Value=" << Value << " Prob()=" << Prob() << " BScale=" << Tree()->BScale(Br) << " -- ret: " << -(Value * Prob() * Tree()->BScale(Br));
#endif
	Tree()->pBra(Br)->grad(-(Value * Prob() * Tree()->BScale(Br)));
	p_a = NULL; p_b = NULL;

}
void CPfamProcess::BranNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError)	{
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\nBranNode_dT: NTo: " << NTo << "; NFr: " << NFr << "; Br: " << Br;
//	cout << "\nP(t): "; OutPT(cout,Br);
#endif
	int i,site,SiteScale;	// Counters
	double Value,Total, *p_a = NULL, *p_c = NULL, *p_d = NULL;
	static double Vec[MAX_CHAR];
	vector <double> eqm;
	assert(InRange(NFr,pTree->NoSeq(),pTree->NoNode()) || First == -1);
	assert(pTree->NodeType(NTo) == branch); assert(m_arModelL != NULL);
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	Value = 0.0;
	FOR(site,m_iSize)	{
		// Put the QMatrix in the m_vpQMat spot
		m_iCurSite = site;
		assert(m_vpQMat.empty());
		m_vpQMat.push_back(m_vpSiteQ[site]);

		// Make the P(t) matrix for the branch
		m_vpSiteQ[site]->MakePT(Tree()->B(Br),PT(Br));
		// Make the QP(t) matrix
		MulMat(m_vpSiteQ[site]->Q(),PT(Br),QP(Br),m_iChar,m_iChar,m_iChar);
		// Build the QP(t) matrix for the site
		eqm = m_vpQMat[QMat4Bra(Br)]->Eqm();

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
//			if(site < 500000) { cout << "\nValue: " << Value << " += " << Total << ", " << SiteScale - ModelL(site).Scale() << " == " << PartialGrad(site,Total,SiteScale - ModelL(site).Scale());;; }
			Value += PartialGrad(site,Total,SiteScale - ModelL(site).Scale());
			// Update the space
			BranNode_Update(NTo, NFr, Br, pTree, First);

		}
		// Clean up
		m_vpQMat[0] = NULL; m_vpQMat.clear();
	}
#if ANALYTIC_DERIVATIVE_DEBUG == 1
	cout << "\n\tBranch[" << Br <<"] = " << Tree()->B(Br) << ": Value=" << Value << " Prob()=" << Prob() << " BScale=" << Tree()->BScale(Br) << " -- ret: " << -(Value * Prob() * Tree()->BScale(Br)) << " = " << PartialGrad(site,Total,SiteScale - ModelL(site).Scale());;
#endif
	Tree()->pBra(Br)->grad(-(Value * Prob() * Tree()->BScale(Br)));
	// Tidy up
	p_a = NULL; p_c = NULL; p_d = NULL;
}

double CPfamProcess::PartialGrad(int site,double Total,int SiteScale)	{

	if(m_pData->IsBootstrap()) { if(m_pData->m_ariPatOcc[site] == 0) { return 0.0; } }
//#if ANALYTIC_DERIVATIVE_DEBUG == 1
//	if(site < ADD_SITE_MAX) {
//		cout << "\n\t\tDoing site: " << site << " [occ=" << m_pData->m_ariPatOcc[site]<<"]: Total= " << Total << "*10^" << SiteScale << "; Partial= " << ModelL(site) << "; return value: " << (Total / ModelL(site).ScalVal() * m_pData->m_ariPatOcc[site] * pow((double)10,-SiteScale) );
//	}
//#endif
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

// Update routines that allow the backward calculations
// Pfam -- This function assumes that there is a single suitable m_vpQmat available for analyses. This may need to be initialised prior to entry into the function
void CPfamProcess::LeafNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoCompleteUpdate)	{
	bool CharCheck,LeafSeq = true;
	int i,site,Seq, SiteScale = 0;	// Counters
	double Vec[MAX_CHAR], Value, *p_a = NULL, *p_b = NULL;
	CProb NewL;
	assert(m_vpQMat.size() == 1 && QMat4Bra(Br) == 0);
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

	// Adjust First so more updating is done if DoCompleteUpdate == true (Only done on entry so m_iCurSite == 0)
//	if(DoCompleteUpdate == true && m_iCurSite == 0 )	{ if(First == 0) { First = -1; } else if(First == 1) { First = 0; } }
	if(DoCompleteUpdate == true)	{ if(First == 0) { First = -1; } else if(First == 1) { First = 0; } }
	// Note: You can't speed up these computations the same way as you can speed up forwards
	site = m_iCurSite;

	if(site==0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << "\nDoing LeafNode_Update(NTo: "<< NTo << ", NFr: " << NFr << ", Br: "<< Br << ")"; }

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
			FOR(i,m_iChar)	{
				*(p_a++) = *(p_b) * Vec[i];
				*(p_b++) = Vec[i];
			}
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

// Pfam -- This function assumes that there is a single suitable m_vpQmat available for analyses. This may need to be initialised prior to entry into the function
void CPfamProcess::BranNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoNTo, bool DoNFr, bool DoCompleteUpdate)	{
	int i,site,SiteScale,SiteScale2;	// Counters
	double *p_a = NULL, *p_b = NULL, *p_c = NULL, *p_d = NULL;
	double Vec[MAX_CHAR], Vec2[MAX_CHAR];
	assert(InRange(NFr,pTree->NoSeq(),pTree->NoNode()) || First == -1);
	assert(pTree->NodeType(NTo) == branch);
	assert(m_vpQMat.size() == 1 && QMat4Bra(Br) == 0);
	site = m_iCurSite;

	if(site==0 && OUTPUT_BRANCH_OPTIMISE == 1) {
		cout << "\nDoing BranNode_Update(NTo: "<< NTo << ", NFr: " << NFr << ", Br: "<< Br << "): First = " << First << " DCU: " << DoCompleteUpdate;
		cout << "\n\tT: " << Tree()->B(Br) << " : P(t): "; FOR(i,3) { cout << "  " << PT(Br)[i]; }
	}
	// Direction makes no difference to calculation so whether first or later nodes doesn't matter
	// Do derivative calculation and updating procedure
	// For Update: if(first == true) {		then Fd(NFr) = vP(t) & Bk(NFr) *= vP(t)
	//										else Fd(NFr) *= vP(t)
	/////////////////////////////////////////////////////////////////////////
	if(DoNTo && DoNFr)	{
		if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " Type 1: "; }
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
		DoScale(NFr,true); DoScale(NTo,true);
	} else if(DoNTo)	{
		if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " Type 2: "; }
		///////////////////////////////////////////////////////
		// Only update NTo
		site = m_iCurSite;
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
		DoScale(NTo,true);
	} else {
		if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " Type 3: "; }
		// Only done on entry
		if(DoCompleteUpdate == true) {  if(First == 0) { First = -1; } else if(First == 1) { First = 0; } }
//		if(m_iCurSite < 10) { cout << "  First = " << First; }
//		if(site < 10) { cout << "\n\tData: "; FOR(i,5) { cout << "  " << ForceRealFd(NTo,site)[i]; } }
		///////////////////////////////////////////////////////
		// Only update NFr
		site = m_iCurSite;
		// Do the updating procedure
		/////////////////////////////////////////////////////////////////
		// i) Update appropriate memory according to First
		switch(First)	{
		case -1:	// Do the update for the StartCalc() node -- Also requires update of BackSp in StartCalc()
			if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " case -1"; }
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
			if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " case 0"; }
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
			if(site == 0 && OUTPUT_BRANCH_OPTIMISE == 1) { cout << " case 1"; }
			break;
		default:
			Error("Unknown first...");
		};
		DoScale(NFr,true);
	}
}

////////////////////////////////////////////////////////////////
// Fast optimisation function for branches
////////////////////////////////////////////////////////////////
// This function will alway require at least 2 passes of the tree
// The first one is to get reasonable sets of branches at a lower tolerance
double CPfamModel::FastBranchOpt(double CurlnL, double tol, bool *Conv, int NoIter, bool CheckPars)	{
	int i,j, Branches = 0;
	double newlnL, BestlnL = 0.0, working_tol;
	assert(CurlnL < 0);
	if(!m_bAllowFastBranchOpt) { return CurlnL; }
//#if FASTBRANCHOPT_DEBUG == 1
	cout << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW ROUND OF FAST BRANCH OPT >>>>>>>>>>>>>>>>>>>>" << flush;
//#endif
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

	// Assumed not neccessary for now

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
	BestlnL = lnL(true);							// 1. Do the first calculation and initialise things
	FOR(i,NoIter)	{
		newlnL = BestlnL;	// new_lnL hold current optimal likelihood
//#if FASTBRANCHOPT_DEBUG == 1
//#if DEVELOPER_BUILD == 1
//		cout << "\n\n--- Round " << i<< ". " << newlnL << " (tol: "<< working_tol << ") ---";;
//		cout << "\nOriginal branches:  "; int j; FOR(j,Tree()->NoBra()) { cout << Tree()->B(j) << " "; }
		Tree()->OutBra(); cout << "\n["<<i<<"]" << *Tree() << " == " << newlnL << flush;
//#endif
		BranchOpt(-1,Tree()->StartCalc(),-1, &BestlnL,working_tol);	// 2. Run the fast optimisation routine
//		double plop = BestlnL; cout << "\nDone branch: BestlnL: " << BestlnL;
		BestlnL = lnL(true);			// Calculation must be done here using true space, otherwise the next round of the optimisation fails.
//		cout << " cf. new: " << BestlnL << " diff = " << fabs(plop - BestlnL);
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
	if(OUTPUT_BRANCH_OPTIMISE == 1) {
		cout << "\nReturning: " << BestlnL << " cf. " << lnL() << " fabs: " << fabs(BestlnL - lnL()); exit(-1);
	}
	assert(BestlnL < 0);
	return BestlnL;
}

// Optimise the set of branches
// PfamModel is more complicated since each site needs its own update
void CPfamModel::BranchOpt(int First,int NTo, int NFr, double *BestlnL,double tol)	{
	int i,j,OriFirst = First;

	if(NFr == -1)	{
		rFOR(i,Tree()->NoLinks(NTo)) { if(Tree()->NodeLink(NTo,i) == -1) { continue; } BranchOpt(First,Tree()->NodeLink(NTo,i),NTo,BestlnL,tol);
	}	} else {
		// Always perform the calculations in the first place
		if(Tree()->NodeType(NFr) == leaf)	{ // Do the leaf calculations for first calc
#if ALLOW_BRANCH_OPTIMISE == 1
			DoBraOpt(First,NTo,NFr,Tree()->NodeBra(NFr,0),true,BestlnL,tol);
#else
			DoBralnL(Tree()->NodeBra(NFr,0),NFr,NTo);
			FOR(i,(int)m_vpAssociatedModels.size()) { m_vpAssociatedModels[i]->Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true); }
			Leaf_update(NTo,NFr,Tree()->NodeBra(NFr,0),Tree(),First,true);

#endif
		} else if(Tree()->NodeType(NTo) == leaf)	{ // Do the leaf calculations for other calcs
#if ALLOW_BRANCH_OPTIMISE == 1
			DoBraOpt(First,NFr,NTo,Tree()->NodeBra(NTo,0),true,BestlnL,tol);
#else
			DoBralnL(Tree()->NodeBra(NTo,0),NTo,NFr);
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
		Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
		FOR(j,(int) m_vpAssociatedModels.size()) {
			m_vpAssociatedModels[j]->Bran_update(NTo,NFr,Tree()->FindBra(NTo,NFr),Tree(),OriFirst,false,true,true);
}	}	}	}










// Function that calculates Node[NTo] * PT(Br), then multiplies elementwise by Node[NFr]
// --
// Note the current implementation works at the double level, whereas normal likelihood function works with CProbs.
// This difference is a potential source of numerical instability.

void CPfamProcess::GetBranchPartL(CProb **arpP, int NT, int NF, int B)	{
//	cout << "\n--> Entered CPfamProcess::GetBranchPartL(CProb **arpP, int NT, int NF, int B)" << flush;
	bool NTreal = true, NFreal = true;
	int i,SiteScale = 0,site;
	double *p_a = NULL, Total = 0.0;
	static double V[MAX_SPACE];
	CProb Pr;
	vector <double> eqm;
//	cout << "\nIn CBaseProcess::GetBranchPartL: eqm: " << eqm;
	// Do Garbage Collector rate
	//////////////////////////////////////////////
	if(MaxRate()) {
		FOR(site,m_pData->m_iSize) {
			eqm = m_vpSiteQ[site]->Eqm();
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
		m_pTree->UpdateB(B);

		// Organise whether NT and NF are real sequences
		if(IsSubTree()) { // Deal with subtrees
			if(NT < Tree()->NoSeq()) { if(m_viLeafMap[NT] == -1) { NTreal = false; } else { NT = m_viLeafMap[NT]; } } else { NTreal = false; }
			if(NF < Tree()->NoSeq()) { if(m_viLeafMap[NF] == -1) { NFreal = false; } else { NF = m_viLeafMap[NF]; } } else { NFreal = false; }
		} else { if(NT >= Tree()->NoSeq()) { NTreal = false; } if(NF >= Tree()->NoSeq()) { NFreal = false; } }
		//////////////////////////////////////////////////////
		// Do the calculations
//		int siteCheck = 3;
		FOR(site,m_pData->m_iSize)	{
//			if(site < siteCheck) { cout << "\nProcSite["<<site<<"] NT=" << NT <<":"; }
			// Put the QMatrix in the m_vpQMat spot
			m_iCurSite = site;
			assert(m_vpQMat.empty());
			m_vpQMat.push_back(m_vpSiteQ[site]);
			// Make the P(t) matrix for the branch
			m_vpSiteQ[site]->MakePT(Tree()->B(B),PT(B));
			// Get the eqm
			eqm = m_vpQMat[QMat4Bra(B)]->Eqm();
//			if(site < siteCheck) { cout << "\n\teqm: " << eqm; }
			SiteScale = 0; Total = 0.0;
			// Get first vector of calc
			if(NTreal)	{
				Data2PartL(m_pData->m_ariSeq[NT][site],PT(B),V,&eqm);
//				FOR(i,m_iChar) { if(i == m_pData->m_ariSeq[NT][site] || m_pData->m_ariSeq[NT][site] == m_iChar) { cout << "\t1"; } else { cout << "\t0"; } }
			} else {
//				if(site < siteCheck) { FOR(i,m_iChar)	{ cout << "\t" << ForceRealFd(NT,site)[i]; } }
				VMat(ForceRealFd(NT,site),PT(B),V,m_iChar); SiteScale += *ForceRealFdSc(NT,site);
			}
//			int j; 	cout << "\n\tRight:\t"; if(NFreal) { FOR(i,m_iChar) { if(i == m_pData->m_ariSeq[NF][site] || m_pData->m_ariSeq[NF][site] == m_iChar) { cout << "\t1"; } else { cout << "\t0"; } } } else { FOR(i,m_iChar)	{ cout << "\t" << ForceRealFd(NF,site)[i]; } } cout << "\n\tLeft * P(t):"; FOR(j,m_iChar) { cout << "\t" << V[j]; }

			// Get calculation of total = sum(Vec[i] = Vec[i] * BranchNode[i] * Eqm[i]);
//			if(site < siteCheck) { cout << "\n"; }
			if(NFreal)	{
//				if(site < siteCheck) { cout << "Sum_Vec"; }
				Total = Sum_Vec(m_pData->m_ariSeq[NF][site],V,eqm);
			} else {	// Do partial likelihoods
				p_a = ForceRealFd(NF,site);
				double *p_test = ForceRealFd(NF,site);
//				if(site < siteCheck) { FOR(i,m_iChar) { cout << "("<<V[i] << "*" << *p_test++ << ")\t"; } }
				FOR(i,m_iChar)	{ Total += V[i] * *(p_a++) * eqm[i]; }
				SiteScale += *ForceRealFdSc(NF,site);
			}
			Total *= Prob();

			// Assign the likelihood
			Pr.Assign(Total,SiteScale);
			arpP[site]->Add(Pr,true);

//			if(site < siteCheck) { cout << "\n\tlogL: " << Pr.LogP(); }

			// Clean up
			m_vpQMat[0] = NULL; m_vpQMat.clear();
		}
	} else {
	/////////////////////////////////////////////////
	// Do zero rate processes
		FOR(site,m_pData->m_iSize)	{
			if(m_pData->m_viNoChange[site] != -1) {
				eqm = m_vpSiteQ[site]->Eqm();
				Pr.Assign(eqm[m_pData->m_viNoChange[site]] * Prob());
				arpP[site]->Add(Pr,true);
	}	}	}
#if DEVELOPER_BUILD == 1
		cout << "\nFinished with Eqm: " << eqm;
#endif

}

///////////////////////////////////////////////////////////////////////
// Space updating functions
void CPfamModel::Leaf_update(int NTo, int NFr, int Br, CTree *T, int First,bool DoFullUpdate)	{
	int i; FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->FullLeafNode_Update(NTo,NFr,Br,T,First,DoFullUpdate);
}	}

void CPfamModel::Bran_update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo, bool DoNFr, bool DoFullUpdate)	{
	int i; FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->FullBranNode_Update(NTo,NFr,Br,T,First,DoNTo,DoNFr,DoFullUpdate);
}	}

void CPfamProcess::FullLeafNode_Update(int NTo, int NFr, int Br, CTree *pTree, int First, bool DoCompleteUpdate) {
	int site;
	m_pTree->UpdateB(Br);
	FOR(site,m_pData->m_iSize) {
		// Put the QMatrix in the m_vpQMat spot
		m_iCurSite = site;
		assert(m_vpQMat.empty());
		m_vpQMat.push_back(m_vpSiteQ[site]);
		// Make the P(t) matrix for the branch
		m_vpSiteQ[site]->MakePT(Tree()->B(Br),PT(Br));
		// Do the update for this site
		LeafNode_Update(NTo,NFr,Br,pTree,First,DoCompleteUpdate);
		// Clean up
		m_vpQMat[0] = NULL; m_vpQMat.clear();
	}
}

void CPfamProcess::FullBranNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo, bool DoNFr, bool DoCompleteUpdate) {
	int site;
	m_pTree->UpdateB(Br);
	FOR(site,m_pData->m_iSize) {
		// Put the QMatrix in the m_vpQMat spot
		m_iCurSite = site;
		assert(m_vpQMat.empty());
		m_vpQMat.push_back(m_vpSiteQ[site]);
		// Make the P(t) matrix for the branch
		m_vpSiteQ[site]->MakePT(Tree()->B(Br),PT(Br));
		// Do the update for this site
		BranNode_Update(NTo,NFr,Br,T,First,DoNTo,DoNFr,DoCompleteUpdate);
		// Clean up
		m_vpQMat[0] = NULL; m_vpQMat.clear();
	}
}





