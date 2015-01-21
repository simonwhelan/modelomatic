/*
 * codon_model.cxx
 * ---
 * File created to house the new set of per site codon models
 * Will be merged with model.cxx in due course...
 *
 *  Created on: May 22, 2014
 *      Author: simon
 */

#include "model.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Site specific codon models
// ---
// Nucleotide models distributed across the 3 codon positions. Can have single or joined models across the codon positions
// The way it's set up BranchPar and ModelPar specify the sets of model and branch parameters and how they're shared between codon positions
// T.ex: {0,0,1} means positions 1 and 2 share parameters and position 3 has its own parameters
//
// The BranchPar and ModelPar can specify -1 at a position and this codon position will not be optimised or included in the calculation,
//		but the likelihood adjustment factor will be calculated and computed in all comparisons.
//
// The model will be set up inefficiently to start with. The main CBaseModel object will have nothing to optimise in it
CSiteCodon::CSiteCodon(CData *D, CTree *T, vector <int> ModelPar, vector <int> BranchPar, EModel CoreModel, bool Gamma) : CBaseModel(D,T)	{
	int i,j,k;
	CData *TempData;
	vector <CTree *> vTempTrees;
	vector <double> FreqCount;
	string t_string;
	// Some error checking on entry
	assert(ModelPar.size() == 3); assert(BranchPar.size() == 3);
	FOR(i,3) {
		assert(ModelPar[i] <= i); assert(InRange(ModelPar[i],-1,3)); assert(BranchPar[i] <= i); assert(InRange(BranchPar[i],-1,3));
		if(ModelPar[i] == -1) { assert(BranchPar[i] == -1); }
	}
	assert(InRange(CoreModel,JC,RY_model)); // Check it's a DNA model
	m_dlnLAdjustment = -BIG_NUMBER;			// Initialise the adjustment factor
//	cout << "\nWorking with core model: " << CoreModel;
//	cout << "\nInitialising okay..." << endl << ModelPar << endl << BranchPar;

	// Create data for first, second, and third codon position
	FOR(i,3) {
		// Create data
		TempData = new CData(0,0,DNA); *TempData = *D;
		vector <bool> WhichSite(3,false);
		WhichSite[i] = true;
		TempData->GetCodonPositions(WhichSite[0],WhichSite[1],WhichSite[2]);
		m_vpDataSites.push_back(TempData);
		TempData = NULL;
	}

	// Fix the frequencies in the data for those data that share a model
	FOR(i,3) { // Do it for each model
		if(ModelPar[i] == -1) { continue; } // Skip sites not included
		FreqCount.clear(); FreqCount.assign(4,0.0);
		FOR(j,3) { if(ModelPar[j] == i) { FOR(k,4) { FreqCount[k] += m_vpDataSites[j]->m_vFreq[k]; } } }
		if(Sum(&FreqCount) < FLT_EPSILON) { continue; }
		FreqCount = NormaliseVector(FreqCount);	// Safe to just average the frequencies because the counts will be the same
		FOR(j,3) { if(ModelPar[j] == i) { FOR(k,4) { m_vpDataSites[j]->m_vFreq[k] = FreqCount[k];  } } }
	}

	// Create the tree objects
	vTempTrees.assign(3,NULL); FOR(i,3) {
		vTempTrees[i] = new CTree(); *vTempTrees[i] = *T;
		if(ModelPar[i] == -1) { continue; }
		t_string = ""; FOR(j,3) { if(BranchPar[j] == i) { t_string = t_string + int_to_string(j+1); } }
		FOR(j,vTempTrees[i]->NoBra()) {
			vTempTrees[i]->pBra(j)->Name("site[" + t_string + "]::" + vTempTrees[i]->pBra(j)->Name());	// Name the branches
			vTempTrees[i]->SetB( j , max( 0.001 , vTempTrees[i]->B(j) + RandDouble(0.0,0.03) - 0.015) ,true,true );								// Slightly randomise the start branch to ensure better optimisation
		}
	}
	m_vpTreeSites.assign(3,NULL);
	FOR(i,3) { if(BranchPar[i]!=-1) { m_vpTreeSites[i] = vTempTrees[BranchPar[i]]; } else { m_vpTreeSites[i] = vTempTrees[2]; } }	// Could potentially cause a memory leak, but it's small fry... (I hope!)
	FOR(i,3) { vTempTrees[i] = NULL; }
	// Create the model objects
	m_viModelMap = ModelPar;
	m_viTreeMap = BranchPar;
	FOR(i,3) {
		m_vpAssociatedModels.push_back(GetMyModel(CoreModel,m_vpDataSites[i],m_vpTreeSites[i]));
		if(m_viModelMap[i] == -1) { continue; }
		t_string = ""; FOR(j,3) { if(m_viModelMap[j] == m_viModelMap[i]) { t_string = t_string + int_to_string(j+1); } }
		FOR(j,(int)m_vpAssociatedModels[i]->m_vpPar.size()) { m_vpAssociatedModels[i]->m_vpPar[j]->Name("site[" + t_string + "]::" + m_vpAssociatedModels[i]->m_vpPar[j]->Name()); }

	}
	NormaliseParameters();
	// Finish by calculating the standard adjustment factor. This could be added to other routines, but kept here for clarity.
	m_dlnLAdjustment = 0;
	FOR(i,3)	{
		if(m_viModelMap[i] == -1) {
			assert(m_viTreeMap[i] == -1);
			m_dlnLAdjustment += m_vpDataSites[i]->GetNT2MissinglnLScale();
	}	}

	cout << "\nAnd the adjustment factor is: " << m_dlnLAdjustment << flush;
	/*
	cout << "\nGrabbing parameters";
	vector <double*> ParVals = GetOptPar(true,true,true,false);
	cout << "\nParameters are: ";
	FOR(i,ParVals.size()) {
		cout << "\nPar["<<i<<"] == " << *ParVals[i] << " : " << *m_vpAllOptPar[i];
	}
*/

//	CBaseModel *Tester = GetMyModel(CoreModel,D,T);

//	cout << "\nComplete model" << Tester->lnL(true);
//	cout << "\nModel: " << *Tester;
//	double Summer = 0;
//	cout << "\nChecking individual submodels...";
//	FOR(i,3) {
//		Summer += m_vpAssociatedModels[i]->lnL(true);
//		cout << "\nModel["<<i<<"]: " << m_vpAssociatedModels[i]->lnL(true);
//		cout << "\n" << *m_vpAssociatedModels[i];
//	}
//	cout << "\nTotal lnL: " << Summer;

	/*
	 * Some issues that need to be addressed when establishing models
	 * ---
	 * 1.The equilibrium distributions from the processes is currently wacky because each instancing of the model gets its own eqm. Need to edit the data objects so the reflect the model assumptions. Can take the sum across codon positions to be consistent
	 * 		General thought: Want to set it up so eqm is never touched so fix here and have an assert during likelihood computation -- Should be done and set in NormaliseParameters(...)
	 * 2. The likelihood function needs at least a wrapper function that ensures consistency of parameters between models that are meant to share parameters -- This function needs to call NormaliseParameters(...)
	 * 3. I need to control which parameters are passed to the optimiser. The ordering of BranchPar and ModelPar mean I have if(i=BranchPar[i]) then I know it's the first instance. Should I be storing this?
	 *
	 */



}

////////////////////////////////////////
// Destructor
CSiteCodon::~CSiteCodon()	{
	int i,j;
	// Clear m_vpDataSites structure
	if(!m_vpDataSites.empty()) {
		assert(m_vpDataSites.size() == 3);
		FOR(i,3) {
			assert(m_vpDataSites[i] != NULL);
			delete m_vpDataSites[i];
		}
		m_vpDataSites.clear();
	}
	// Clear m_vpTreeSites structure
	if(!m_vpTreeSites.empty()) {
		assert(m_vpTreeSites.size() == 3);
		FOR(i,3) {
			// Clean up the hierarchy of pointers
			if(m_vpTreeSites[i] == NULL) { continue; }
			for(j=i+1;j<3;j++) {
				if(m_vpTreeSites[j] == m_vpTreeSites[i]) { m_vpTreeSites[j] = NULL; }
			}
			delete m_vpTreeSites[i]; m_vpTreeSites[i] = NULL;
		}
		m_vpTreeSites.clear();
	}
	// Clear pointer information
	m_viModelMap.clear(); m_viTreeMap.clear();

}

////////////////////////////////////////////////
// Function that normalises model parameters
// (branches should be done automatically)
// Assumes models are all the same
bool CSiteCodon::NormaliseParameters()	{
	int site,i,j;
//	cout << "\nInto NormaliseParameters(...)" << flush;
	vector <int> Checks(3,-1);	// Captures first usage of model and takes a link to it
	// Check input conditions
//	cout << "\nModelMap: " << m_viModelMap;
	assert(m_viModelMap.size() == 3); FOR(i,3) { assert(m_viModelMap[i] <= i); }

//	cout << "\nGoing through sites..." << flush;
	// Loop through the sites and do the model parameters
	FOR(site,3) {
//		cout << "\n --- Site " << site << " --- ";
		// Skip sites where the model is not applied
		if(m_viModelMap[site] == -1) { assert(m_viTreeMap[site] == -1); continue; }
		// Case when the model's not been seen before
		if(Checks[m_viModelMap[site]] == -1) {
//			cout << "\nFirst instance of Model[" << m_viModelMap[site] << "] at site " << site;
			Checks[m_viModelMap[site]] = site; continue;
		} else {
//			cout << "\n\tModel at site[" << site << "] seen before at " << m_viModelMap[site];
			assert(m_vpAssociatedModels[site]->m_vpProc.size() == m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc.size());
			FOR(i,(int)m_vpAssociatedModels[site]->m_vpProc.size()) {
//				cout << "\nProcess["<<i<<"]: ";
				if(m_vpAssociatedModels[site]->m_vpProc[i]->PseudoProcess()) { continue; }
				assert(m_vpAssociatedModels[site]->m_vpProc[i]->NoPar() == m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->NoPar());
				FOR(j,(int)m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->NoPar()) {
					assert(m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j)->Name() == m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->pPar(j)->Name());
					if(m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j)->Name().find("Freq") != string::npos) { assert(fabs(m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j)->Val() - m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->pPar(j)->Val()) < 0.0001); continue; }
					m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j)->SetVal(m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->pPar(j)->Val(),true);
					m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j)->SetOptimise(false);	// Does every time as its just as quick to assign a bool than to test it
//					cout << "\nTaking Parameter " << *m_vpAssociatedModels[m_viModelMap[site]]->m_vpProc[i]->pPar(j) << " --> " << *m_vpAssociatedModels[site]->m_vpProc[i]->pPar(j);
				}
			}
		}
	}
	return true;
}

/////////////////////////////////////////////
// Override of the lnL function
double CSiteCodon::lnL(bool ForceReal) {
	int i;
	double logL = 0.0;
	bool AddAdj = false;
	NormaliseParameters();
	assert(m_vpAssociatedModels.size() == 3);
	if(m_vbUseInBraCalc.empty()) {
		FOR(i,(int) m_vpAssociatedModels.size()) {
			if(m_viModelMap[i] == -1) { AddAdj = true; continue; }
			logL += m_vpAssociatedModels[i]->lnL(ForceReal);
		}
	} else {	// Should only happen for FastBranch calculations. Here for debugging
		assert(m_vbUseInBraCalc.size()==3);
		FOR(i,(int) m_vpAssociatedModels.size()) {

			if(m_vbUseInBraCalc[i]) {
				if(m_viModelMap[i] == -1) { AddAdj = true; continue; }
				logL += m_vpAssociatedModels[i]->lnL(ForceReal);
			}
	}	}
	if(AddAdj) logL += m_dlnLAdjustment;
	return logL;
}

//////////////////////////////////////////////////////////////////////////////////////////
// --- Branch optimisation routines ---
//////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////
// DoFastBraOpt -- Main driver function
// ---
// Goes through the available trees and optimises each of their branches
double CSiteCodon::FastBranchOpt(double CurlnL, double tol, bool *Conv, int NoIter, bool CheckPars) 	{
	int site,i,j,k;
	vector <bool> TreeDone(3,false);
	double working_tol, BestlnL, newlnL, RetVal;
cout << "\nCSiteCodon::FastBranchOpt(...) with " << flush; cout << lnL(true); FOR(i,3) { if(m_viTreeMap[i] != -1) {  cout << " ["<<i<<"]: " << m_vpAssociatedModels[i]->lnL(true); } } cout << " ... done" << flush;

	// Checking everything okay going in
	assert(m_vbUseInBraCalc.empty()); assert(m_vpAssociatedModels.size() == 3);
	assert(m_viModelMap.size() ==3 && m_viTreeMap.size() == 3); FOR(i,3) { assert(InRange(m_viTreeMap[i],-1,3)); }
	// Loop through the sites and see what branches need to be optimised
	NormaliseParameters();
//	cout << "\nTreeMap: " << m_viTreeMap;
	RetVal = 0.0;
	FOR(site,3) {
		if(m_viTreeMap[site] == -1) { continue; }
//		cout << "\n\tDoing branch set["<<site<<"] == " << m_viTreeMap[site];
		// Check whether the trees been done yet
		if(TreeDone[m_viTreeMap[site]]) { continue; } TreeDone[m_viTreeMap[site]] = true;
		// If not, find the other processes associated with it. Also initialises the model ready for computation
		m_vbUseInBraCalc.assign(3,false); CurlnL = 0;
		FOR(j,3) {
			if(m_viTreeMap[site] == -1) { m_vbUseInBraCalc[j] = false; continue; }
			if(m_viTreeMap[site] == m_viTreeMap[j]) {
				m_vbUseInBraCalc[j] = true; CurlnL += m_vpAssociatedModels[j]->lnL(true);
				FOR(k,m_vpAssociatedModels[site]->m_vpProc.size()) { m_vpAssociatedModels[site]->m_vpProc[k]->PrepareBraDer(); }
		}	}
		if(tol > 1.0E-3) { working_tol = tol; } else { working_tol = 1.0E-3; }

		///////////////////////////////////////////////////////////////////
		// Only do cyclical optimisation with multiple branches
		if(Tree()->NoBra() == 1) {
			DoBraOpt(true,0,1,0,true,&CurlnL,tol,false);
			continue;
		}
		FOR(i,NoIter)	{
			BestlnL = newlnL = CurlnL;							// 1. Do the first calculation
	#if FASTBRANCHOPT_DEBUG == 1
	//#if DEVELOPER_BUILD == 1
			cout << "\n\n--- Round " << i<< ". " << newlnL << " (tol: "<< working_tol << ") ---";;
			cout << "\nOriginal branches:  "; int j; FOR(j,m_vpAssociatedModels[site]->Tree()->NoBra()) { cout << m_vpAssociatedModels[site]->Tree()->B(j) << " "; }
			cout << flush;
	#endif
			BranchOpt(-1,m_vpAssociatedModels[site]->Tree()->StartCalc(),-1, &BestlnL,working_tol);	// 2. Run the fast optimisation routine
			BestlnL = lnL(true);	// This is also important as it correctly cleans up the space and reassigns things
			if(working_tol > tol) { working_tol = max(tol,working_tol/10); }
	//#if DEVELOPER_BUILD == 1
	#if FASTBRANCHOPT_DEBUG == 1
			cout << "; " << BestlnL << " == " << lnL() << "; diff = " << BestlnL - lnL();
			cout << "\nTree: " << *Tree();
			if(fabs(BestlnL - lnL()) > tol) { cout << "\nBig Error..."; exit(-1); }
	#endif
//			cout << "\nIter["<<i<< "]: best: " << BestlnL << " cf. " << newlnL << " and diff: " << fabs(BestlnL - newlnL) << " cf. tol=" << tol;
			if(fabs(BestlnL - newlnL) < tol) { break; }				// 3. Control exit
			if(newlnL - BestlnL > tol * 100) { cout << "\nDetected an increase in likelihood for CSiteCodon::FastBranchOpt(...) -- new = " << newlnL << " cf. best = " << BestlnL << " diff = " << newlnL - BestlnL; exit(-1); }
			CurlnL = BestlnL;
		}
		RetVal += BestlnL;
		if(Conv != NULL) { if(i==NoIter) { *Conv = false; } else { *Conv = true; } }
	//	cout << "\nReturning: " << BestlnL << " cf. " << lnL() << " fabs: " << fabs(BestlnL - lnL()); // exit(-1);
		assert(BestlnL < 0);
		// Finish up for this subset of data
		m_vbUseInBraCalc.clear();
	}
	CurlnL = 0; FOR(i,3) { if(m_viTreeMap[i] == -1) { continue; } CurlnL += m_vpAssociatedModels[i]->lnL(true); } 	// There's a better way to get this number without recalculating everything
//	cout << "\n\t ... Done and returning: " << CurlnL << " cf. " << RetVal << " == " << CurlnL - RetVal; FOR(i,3) { cout << " [" << i<< "]: " << m_vpAssociatedModels[i]->lnL(true); }
	return CurlnL + m_dlnLAdjustment;
}

/////////////////////////////////////////////
// Override the DoBraOpt function
// Optimise a single branch providing the partial likelihoods are correctly assigned
// This probably didn't need over-riding so brutally, but I'm doing it to ensure everything works cleanly
#define RETURN_DOBRAOPT Par = NULL; FOR(i,3) { if(m_vbUseInBraCalc[i]) { m_vpAssociatedModels[i]->PreparePT(Br); } } \
	if(AllowUpdate) { if(IsExtBra) { FOR(i,3) { if(m_vbUseInBraCalc[i]) { m_vpAssociatedModels[i]->Leaf_update(NTo,NFr,Br,Tree(),First,true); } } } \
	else { FOR(i,3) { if(m_vbUseInBraCalc[i]) { m_vpAssociatedModels[i]->Bran_update(NTo,NFr,Br,Tree(),First,true,false); } } } }  return;
#define GS_DELTA 5
#define DO_BRENT 0		// Whether to do Brent's algorithm in fast search

const double phi = (1 + sqrt(5)) /2;
const double resphi = 2 - phi;


void CSiteCodon::DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol,bool AllowUpdate) {

	cout << "CSITECODON::DOBRAOPT NEEDS FIXING..."; exit(-1);

	int i;
	double *p_x,x1,x2,x3,x1_lnL = 1.0,x2_lnL = -fabs(*BestlnL),x3_lnL = 1.0, xi,temp, ori_x, ori_lnL;
	tol = max(tol,FULL_LIK_ACC);
	double dx = DX;
	CPar *Par;
	// Some stuff to check at the beginning
	NormaliseParameters();
	assert(m_vpAssociatedModels.size() == 3 && m_viTreeMap.size() == 3 && m_vbUseInBraCalc.size() == 3);
	FOR(i,3) { if(m_vbUseInBraCalc[i]) { p_x = m_vpAssociatedModels[i]->Tree()->OptimiserB(Br); Par = m_vpAssociatedModels[i]->Tree()->pBra(Br); break; } }
	assert(i!=3);
	ori_x = x2 = *p_x;
//	cout << "\n---\nBranch["<<Br<<"] has DoBralnL: "<< DoBralnL(Br,NFr,NTo) << " cf. " << DoBralnL(Br,NFr,NTo) << " cf. best: " << *BestlnL; cout << " Using parameter: " << *Par << " == " << *p_x << " == " << ori_x;
//	if(fabs(DoBralnL(Br,NFr,NTo) - *BestlnL) > 0.001) { cout.precision(12); cout << "\nError: DoBralnL("<<Br<<","<<NFr<< "," << NTo << "): " << DoBralnL(Br,NFr,NTo) << " cf. " << *BestlnL; cout << " ... and full lnL: " << lnL(true); exit(-1); }
	ori_lnL = *BestlnL; // = DoBralnL(Br,NFr,NTo);
	// ------------------------------------- Catch for when at lower bounds branch length and whether to proceed ------------------------------------
	if(x2 < DX) {	// There's some annoying behaviour with low bounds
//		cout << "\n--- Caught small branch length ---";
		x1 = *p_x = ori_x; x1_lnL = ori_lnL; // = DoBralnL(Br,NFr,NTo);	// This extra call should be unnecessary, but for some reason it's needed...
		*p_x +=DX; x2 = *p_x; x2_lnL = DoBralnL(Br,NFr,NTo);
		if(x1_lnL > x2_lnL)	{	// The return point if it's really at the lower bound
			*p_x = ori_x; *BestlnL = ori_lnL;
			RETURN_DOBRAOPT;
		}
		*p_x = x1 = ori_x; x1_lnL = ori_lnL; // DoBralnL(Br,NFr,NTo);
		assert(x2_lnL > x1_lnL);	// Final sanity check
	}

	// ----------------------------------- Bracketing routine -------------------------------------
	// Get left bracketing
	while(x2_lnL < x1_lnL)	{
		dx = DX;
		*p_x = x1 = max(0,x2 - (dx * 100));
		x1_lnL = DoBralnL(Br,NFr,NTo); m_iFastBralnL_Bracket++;
//		cout << "L[" << *p_x << "]" << flush;
		if(x1_lnL > x2_lnL)	{
			x3 = x2; x3_lnL = x2_lnL;
			x2 = x1; x2_lnL = x1_lnL;
			x1_lnL = 1.0;
		}
		dx *= GS_DELTA;
		if(!Par->CheckLowBound()) { x1 = *p_x; break; }
	}
	// Get right bracketing
	if(x3_lnL > 0.0)	{
		while(x2_lnL < x3_lnL)	{
			*p_x = x3 = max(0,max(x2 + (dx* 100),0.01));
			x3_lnL = DoBralnL(Br,NFr,NTo); m_iFastBralnL_Bracket++;
			x3 = *p_x;
//			cout << "R[" << *p_x << "]" << flush;
//			cout << "\nRight: " << x3 << " == " << x3_lnL;
			if(x3_lnL > x2_lnL)	{
				x1 = x2; x1_lnL = x2_lnL;
				x2 = x3; x2_lnL = x3_lnL;
				x3_lnL = 1.0;
			}
			dx *= GS_DELTA;
			if(!Par->CheckUpBound()) { break; }
	}	}
	if(x1 < 0 || x2 < 0 || x3 < 0 || *p_x < 0) { cout << "\nHave detected an error for DoBraOpt: x1: " << x1 << " x2: " <<  x2 << " x3: " <<  x3 <<" x: " <<  *p_x << "!!!"; exit(-1); }
//	cout << "\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
//	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")"; *p_x = x2; cout << " checking " << DoBralnL(Br,NFr,NTo);
//	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
	// Do various checks to make sure it looks like a hill
	if(!Par->CheckBound()) { *p_x = x2; *BestlnL = x2_lnL; RETURN_DOBRAOPT; }
	if(x1_lnL > max(x2_lnL,x3_lnL)) {
		// Deal with the U shape case
		if(x3_lnL > x2_lnL) {
			if(x1_lnL > x3_lnL) { *p_x = x1; DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol); }
			else				{ *p_x = x3; DoBraOpt(First,NTo,NFr,Br,IsExtBra,BestlnL,tol); }
		} else { *p_x = x1; *BestlnL = x1_lnL; RETURN_DOBRAOPT; }
	}
	else if (x3_lnL > x2_lnL)		{ *p_x = x3; *BestlnL = x3_lnL; RETURN_DOBRAOPT; }
	// New Goldensection
	assert(x2_lnL + FLT_EPSILON > x1_lnL && x2_lnL + FLT_EPSILON> x3_lnL);
//	cout << "\nInto Goldensection (" << x1 << "," << x2 << "," << x3 << ")";
	FOR(i,20) {
		// New value
		if(x3 - x2 > x2 - x1) 	{ *p_x = xi = x2 + resphi * (x3-x2); }
		else					{ *p_x = xi = x2 - resphi * (x2-x1); }
		// break condition
		if(fabs(x3_lnL - x1_lnL) < tol) { /* cout << "\nBreaking at tol=" << tol << " fabs(" << x3_lnL << " - " << x1_lnL << ")";  */ *p_x = x2 = (x1+x3)/2; break; }
		// Search
		temp = DoBralnL(Br,NFr,NTo);
//		cout << "\n\t[i="<<i<<"] xi:" << xi << ": " << temp;
		if(temp > x2_lnL) {
			if(x3 - x2 > x2 - x1) 	{ x1 = x2; x1_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
			else					{ x3 = x2; x3_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
		} else {
			if(x3 - x2 > x2 - x1)	{ x3 = xi; x3_lnL = temp; }
			else					{ x1 = xi; x1_lnL = temp; }
		}
//		cout << " [x1:" << x1 << ", x2:" << x2 << ", x3:" << x3 << ","<< max(x1_lnL,max(x2_lnL,x3_lnL)) << "]" << flush;
	}

/*
	cout << "\nFinished search: x: " << *p_x << " = " << temp << " == " << DoBralnL(Br,NFr,NTo);
	cout << ": tol= " << max(x2_lnL - x1_lnL,x2_lnL - x3_lnL);
	cout << "\n---\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
*/
	m_iFastBralnL_Calls++;

	if(x2_lnL > ori_lnL) { *p_x = x2; *BestlnL = x2_lnL; } else { if(fabs(x2_lnL - ori_lnL) > tol) { cout << "\nWeird... optimiser (tol=" << tol << ") made worse likelihood in CSiteCodon::DoBraOpt(...)\n"; cout << " should have: " << ori_lnL << " and have " << DoBralnL(Br,NFr,NTo) << " diff = " << DoBralnL(Br,NFr,NTo) - ori_lnL; } *p_x = ori_x; *BestlnL = ori_lnL;  }
//	*BestlnL = DoBralnL(Br,NFr,NTo); // Can probably avoid this...
//	cout << "\nRETURNING p_x = " << *p_x << ": " << *BestlnL << " cf. " << DoBralnL(Br,NFr,NTo);
//	exit(-1);
	RETURN_DOBRAOPT;

}

/////////////////////////////////////////////
// Override the DoBraOpt function
// Optimise the set of branches
void CSiteCodon::BranchOpt(int First,int NTo, int NFr, double *BestlnL,double tol)	{
	int i,j,OriFirst = First;
	CTree *ProcTree = NULL;
	FOR(i,3) { if(m_vbUseInBraCalc[i]) { ProcTree = m_vpAssociatedModels[i]->Tree(); } }
	if(NFr == -1)	{
		rFOR(i,ProcTree->NoLinks(NTo)) { if(ProcTree->NodeLink(NTo,i) == -1) { continue; } BranchOpt(First,ProcTree->NodeLink(NTo,i),NTo,BestlnL,tol);
	}	} else {
		// Always perform the calculations in the first place
		if(ProcTree->NodeType(NFr) == leaf)	{ // Do the leaf calculations for first calc
			DoBraOpt(First,NTo,NFr,ProcTree->NodeBra(NFr,0),true,BestlnL,tol);
		} else if(ProcTree->NodeType(NTo) == leaf)	{ // Do the leaf calculations for other calcs
			DoBraOpt(First,NFr,NTo,ProcTree->NodeBra(NTo,0),true,BestlnL,tol);
		} else { // Do the internal calculations
			FOR(i,ProcTree->NoLinks(NTo))	{ if(ProcTree->NodeLink(NTo,i) == NFr || ProcTree->NodeLink(NTo,i) == -1) { break; } }
			assert(i != ProcTree->NoLinks(NTo));
			// If the node from isn't a leaf node do the internal calculation (i.e. avoids first node)
			if(ProcTree->NodeType(ProcTree->NodeLink(NTo,i)) != leaf)	{
				DoBraOpt(First,NTo,NFr, ProcTree->NodeBra(NTo,i),false,BestlnL,tol);
		}	}
		// Do the looping
		First = 0;
		FOR(i,ProcTree->NoLinks(NTo))	{
			if(ProcTree->NodeLink(NTo,i) == NFr || ProcTree->NodeLink(NTo,i) == -1) { continue; }
			BranchOpt(First,ProcTree->NodeLink(NTo,i),NTo,BestlnL,tol);
			First = 1;
	}	}
	if(NFr != -1) { if(ProcTree->NodeType(NTo) != leaf && ProcTree->NodeType(NFr) != leaf) {
		// Do the updates needed for the process
		FOR(i,3) {
			if(!m_vbUseInBraCalc[i]) { continue; }
			m_vpAssociatedModels[i]->PreparePT(ProcTree->FindBra(NTo,NFr));
			m_vpAssociatedModels[i]->Bran_update(NTo,NFr,ProcTree->FindBra(NTo,NFr),ProcTree,OriFirst,false,true,true);
	}	}	}
	ProcTree = NULL;
}


/////////////////////////////////////////////
// The DoBralnL routine for CSiteCodon models
// ---
// Excessive overload, but I don't want DoBralnL to become public
double CSiteCodon::DoBralnL(int B, int NF,int NT)	{

	int site,i;
	double RetlogL = 0.0, logL = 0.0; // Having position specific and return log separate for debugging
	CProb **P = NULL;
	int CProbSize = -1;
//	cout << "\nInto DoBralnL(B=" << B << "; NF: " << NF << "; NT: " << NT << ") :: m_vbUseInBraCalc " << m_vbUseInBraCalc;
	// Check going in
	assert(m_vpAssociatedModels.size() == 3 && m_viTreeMap.size() == 3 && m_vbUseInBraCalc.size() == 3);
//	cout << "\nUsing sites ";
	FOR(site,3) {
		if(!m_vbUseInBraCalc[site]) { continue; }	// Skip processes not needed for this likelihood
		// Get the memory
		logL = 0.0;
		GET_MEM(P,CProb *,m_vpAssociatedModels[site]->m_pData->m_iSize);
		FOR(i,m_vpAssociatedModels[site]->m_pData->m_iSize) { P[i] = new CProb(0.0,0); }
		CProbSize = m_vpAssociatedModels[site]->m_pData->m_iSize;
		// Get the partial likelihoods
		FOR(i,(int)m_vpAssociatedModels[site]->m_vpProc.size())	{
			if(m_vpAssociatedModels[site]->m_vpProc[i]->Prob() < MIN_PROB) { continue; }
			m_vpAssociatedModels[site]->m_vpProc[i]->GetBranchPartL(P,NT,NF,B);
		}
		FOR(i,m_vpAssociatedModels[site]->m_pData->m_iSize) { logL += P[i]->LogP() * m_vpAssociatedModels[site]->m_pData->m_ariPatOcc[i]; }
//		cout << "\n\tLikelihood for site["<<site<<"] -- lnL: " << logL; // << " cf. " << lnL(true);
		// Apply extra stuff to the likelihood function if needed
		if(m_vpAssociatedModels[site]->pLikelihood != NULL) {
			logL -= m_vpAssociatedModels[site]->pLikelihood(NULL); // Called as blank. Other arguments intended to allow functionality
		}
	//	cout << "\n>>>>>>>>>>>>>>>>>>>> DoBralnL(Branch="<<B<<") returning: "<< logL;
		m_iFastBralnL++;
//		cout << " ["<< site<<"]+= " << logL;
		RetlogL += logL;
		// Clear up after the process
		FOR(i,CProbSize) { delete P[i]; } delete [] P;
	}
//	cout << " == " << RetlogL;
	return RetlogL;
}


/////////////////////////////////////////////
// Override the DoBraOpt function
vector <double *> CSiteCodon::GetOptPar(bool ExtBranch, bool IntBranch, bool Parameters, bool Eqm)	{
	int site,i,j, grad_pointer = 0;
	bool ProcEBra, ProcIBra, ProcPar, ProcEqm;
	vector <double *> OptVal, PlaceHold;
	vector <bool> TreeDone(3,false), ParDone(3,false);
	// Error checking going in
	assert(m_vpAssociatedModels.size() == 3);

//	cout << "\nInto CSiteCodon::GetOptPar(...)" << flush;

	// Clean the parameter space
	FOR(i,(int)m_vpAllOptPar.size()) { m_vpAllOptPar[i] = NULL; } m_vpAllOptPar.clear();
	m_vbDoBranchDer.clear();
	// Loop through the individual site models collecting the necessary parameters as we go
	FOR(site,3) {
		if(m_viTreeMap[site] == -1) { continue; }
		// Work out what needs to be collected
		ProcEBra = ExtBranch; ProcIBra = IntBranch; ProcPar = Parameters; ProcEqm = Eqm;
		if(TreeDone[m_viTreeMap[site]] == false) {	TreeDone[m_viTreeMap[site]] = true; } else { ProcEBra = false; ProcIBra = false; }	// Only collect first instance of branches/parameters
		if(ParDone[m_viModelMap[site]] == false) {	ParDone[m_viTreeMap[site]] = true; } else { ProcPar = false; ProcEqm = false; }
		// Get the values
		PlaceHold = m_vpAssociatedModels[site]->GetOptPar(ProcEBra,ProcIBra,ProcPar,ProcEqm);
		assert(PlaceHold.size() == m_vpAssociatedModels[site]->m_vpAllOptPar.size());
		FOR(i,(int)PlaceHold.size()) { OptVal.push_back(PlaceHold[i]); PlaceHold[i] = NULL; m_vpAllOptPar.push_back(m_vpAssociatedModels[site]->m_vpAllOptPar[i]); } PlaceHold.clear();
//		cout << "\nProcess["<<site<<"]: totalpar = " << OptVal.size();


	} /*
		cout << "\nGetting Opt Par ["<< m_vpAllOptPar.size();
		cout << ":" << OptVal.size() <<"]";
		FOR(i,m_vpAllOptPar.size()) {
			cout << "\n\tPar["<<i<<"] " << m_vpAllOptPar[i]->Name() << " = " << m_vpAllOptPar[i]->Val() << " == " << *OptVal[i];
		}
		cout << "\n---------" << flush;
*/
	return OptVal;
}

/////////////////////////////////////////////////////////////////////////
// Derivative calculations
vector <double> CSiteCodon::GetDerivatives(double CurlnL, bool *pOK)	{
	int i;
	bool OK = true, ForceNumBra = false;
	vector <double> Grads;
	vector <double> temp;
	double temp_lnL;
	if(m_vpAllOptPar.empty())  { return Grads; }
	assert(m_vpAssociatedModels.size() == 3);
//	cout << "\nInto GetDer";
	if(IsRMSDCalc())	{	////////////////////////// Do RMSD derivatives ///////////////////////
		FOR(i,(int)m_vpAllOptPar.size()) {
			Grads.push_back(m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL)));
			if(Grads[i] > RMSD_GRAD_LIM)		{ OK = false; Grads[i] = RMSD_GRAD_LIM; }
			else if(Grads[i] < -RMSD_GRAD_LIM)	{ OK = false; Grads[i] = -RMSD_GRAD_LIM; }
	}	} else {			////////////////////////// Do likelihood derivatives /////////////////
		FOR(i,(int)m_vpAllOptPar.size())	{
//			cout << "\n\tOptimising: " << m_vpAllOptPar[i]->Name();
			m_vpAllOptPar[i]->grad(GetNumDerivative(m_vpAllOptPar[i]->OptimiserValue(),CurlnL)); }
		// Store the gradients
		FOR(i,(int)m_vpAllOptPar.size()) { Grads.push_back(m_vpAllOptPar[i]->grad()); }
	}
	if(pOK != NULL) { *pOK = OK; }
	return Grads;
}

////////////////////////////////////////////////////////////////////////////
// Empirical codon models
CEMPCodonREST::CEMPCodonREST(CData *D, CTree *Tree, bool PlusFreq, int GenCode) : CBaseModel(D,Tree) {
	int i;
	string m_sName = sModelNames[(int)CodonEMPRest];
	assert(GenCode == 0); 						// Currently only available for universal genetic code
	// Do genetic code
	D->MakeCodonData();
	// Reduce the model and the data to the correct genetic code
	m_pData->ReduceCodonData(GenCode);

	// Sort out the eqm distribution (should probably be dealt with else where, but this works
	if(PlusFreq) {		// The '+F' option comparable to amino acid empirical models. Note only works with F64
		assert(D->m_vFreq.size() == 61);
		FOR(i,61) { D->m_vFreq[i] = dECMrestFreq[i]; }
	}
	// Add the process
	m_vpProc.push_back(AddCodonProcess(D,Tree,pCodonEMPRest,F64,GenCode));
	FinalInitialisation();
}

CEMPCodonUNREST::CEMPCodonUNREST(CData *D, CTree *Tree, bool PlusFreq, int GenCode) : CBaseModel(D,Tree) {
	int i;
	string m_sName = sModelNames[(int)CodonEMPRest];
	assert(GenCode == 0); 						// Currently only available for universal genetic code
	// Do genetic code
	D->MakeCodonData();
	// Reduce the model and the data to the correct genetic code
	m_pData->ReduceCodonData(GenCode);

	// Sort out the eqm distribution (should probably be dealt with else where, but this works
	if(PlusFreq) {		// The '+F' option comparable to amino acid empirical models. Note only works with F64
		assert(D->m_vFreq.size() == 61);
		FOR(i,61) { D->m_vFreq[i] = dECMunrestFreq[i]; }
	}
	// Add the process
	m_vpProc.push_back(AddCodonProcess(D,Tree,pCodonEMPUnrest,F64,GenCode));
	FinalInitialisation();
}

///////////////////////////////////////////////////////////////////////////////////////////
// A codon model built from an empirical amino acid model
CAAEMPCodon::CAAEMPCodon(CData *D, CTree *Tree, ECodonEqm CE,  int GenCode) : CBaseModel(D,Tree) {
	// Maths in development...
	cout << "\nAA empirical model is not ready yet... "; exit(0);
	/* Code below comes from M0 and needs modifying
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
	*/
}

///////////////////////////////////////////////////////////////////////////////////////////
// A codon model for conservative and radical amino acid changes, defined by matrix Radical.mat
CCodonDrDc::CCodonDrDc(CData *Data, CTree *Tree, ECodonEqm CE, string RadicalFile, int GenCode) : CBaseModel(Data,Tree) {
	int i, j;
	m_sName = sModelNames[(int)CodonM0_DrDc];
	// Do genetic code
	Data->MakeCodonData();
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
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,RadicalFile));

	// Store the RadicalFile matrix
	m_sRadicalFile = RadicalFile;
	m_viRadMat.assign(20*20,-1);
	FINOPEN(Radin, RadicalFile.c_str());
		FOR(i,20)	{
			FOR(j,i)	{
				Radin >> m_viRadMat[(i*20)+j];
				if(!InRange(m_viRadMat[(i*20)+j],0,2)) { cout << "\nError reading Radical Matrix from: " << RadicalFile << " at ["<<i << "," << j << "] = " << m_viRadMat[(i*20)+j] << "\nMatrix so far: " << MatOut(20,m_viRadMat); exit(-1); }
				m_viRadMat[(j*20)+i] = m_viRadMat[(i*20)+j];
			}
		}
	Radin.close();

	FinalInitialisation();

}

///////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate how often we would expect to observe particular changes
// Should provide numbers comparable to Kr/Kc
// 1. Gets the expected number of substitutions that actually occur per unit time using the rate matrix
double GetAminoAcidCountFromCodonQ(CQMat *Mat, int GenCode, vector <int> RadMat, int ChangeType) {
	int i,j,i_count, CurChar_i, CurChar_j;
	bool RedData = false;	// Whether data is reduced from 64 codons
	double ExpObs = 0.0;
	vector <double> eqm;

//	cout << "\n------------------------ Calculating overall obs for ChangeType: " << ChangeType << " ------------------- ";

	// Some input checking
	assert(InRange(ChangeType,0,2));
	assert(Mat != NULL); assert(InRange(GenCode,0,NumGenCode));
	i_count = 0.0; FOR(i,64) { if(GenCodes[GenCode][i] != -1) { i_count++; } }
	assert(Mat->Char() == i_count || Mat->Char() == 64);
	if(Mat->Char() != 64) { RedData = true; }
	// Do the counting
	eqm = Mat->Eqm();
	assert(eqm.size() == Mat->Char());
	CurChar_i = 0;
	FOR(i,64)	{
		CurChar_j = CurChar_i+1;
		if(GenCodes[GenCode][i] == -1) { if(!RedData) { CurChar_i++; } continue; }
		for(j=i+1;j<64;j++) {
			if(GenCodes[GenCode][j] == -1) { if(!RedData) { CurChar_j++; } continue; }
			if(GenCodes[GenCode][i] != GenCodes[GenCode][j] && RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] == ChangeType) {
//				cout << "\nGetting " << State(COD,i) << "[" << GenCodes[GenCode][i] << "] -> " << State(COD,j) << "[" << GenCodes[GenCode][j] << "] == RadMat: " <<  RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] << " == " << ChangeType;
				ExpObs += eqm[CurChar_i] * *Mat->Q(CurChar_i,CurChar_j);			// i -> j
				ExpObs += eqm[CurChar_j] * *Mat->Q(CurChar_j,CurChar_i);			// j -> i
			}
			CurChar_j++;
		}
		CurChar_i++;
	}
//	cout << "\nOverall obs for ChangeType[" << ChangeType << "]: " << ExpObs << "\n//";
	return ExpObs;
}
// 2. Gets the expected number of observable substitutions after time t, using a given P(t) matrix. Needs QMat for eqm
double GetAminoAcidCountFromCodonPt(CQMat *QMat, double time, int GenCode, vector <int> RadMat, int ChangeType) {
	int i,j,i_count, CurChar_i, CurChar_j;
	bool RedData = false;	// Whether data is reduced from 64 codons
	double ExpObs = 0.0;
	vector <double> WorkMat;
	vector <double> eqm;
//	cout << "\n------------------------ Calculating overall obs for ChangeType: " << ChangeType << " ------------------- ";

	// Some input checking
	assert(InRange(ChangeType,0,2));
	 assert(InRange(GenCode,0,NumGenCode));
	i_count = 0.0; FOR(i,64) { if(GenCodes[GenCode][i] != -1) { i_count++; } }
	assert(QMat->Char() == i_count || QMat->Char() == 64);
	if(QMat->Char() != 64) { RedData = true; }

	// If time is 0 (or lower for simplicity of -1 pass) pass the Q matrix to WorkMat
	if(time <= FLT_EPSILON) { FOR(i,QMat->Char() ) { FOR(j,QMat->Char() ) { WorkMat.push_back(*QMat->Q(i,j)); } } }
	// Otherwise make a P(t) matrix and transfer that to WorkMat
	else {
		double PT[64*64];
		QMat->MakePT(time,PT);
		FOR(i,QMat->Char() * QMat->Char()) { WorkMat.push_back(PT[i]); }
	}

	// Do the counting
	eqm = QMat->Eqm();
	assert(eqm.size() == QMat->Char());
	CurChar_i = 0;
	FOR(i,64)	{
		CurChar_j = CurChar_i+1;
		if(GenCodes[GenCode][i] == -1) { if(!RedData) { CurChar_i++; } continue; }
		for(j=i+1;j<64;j++) {
			if(GenCodes[GenCode][j] == -1) { if(!RedData) { CurChar_j++; } continue; }
			if(GenCodes[GenCode][i] != GenCodes[GenCode][j] && RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] == ChangeType) {
//				cout << "\nGetting " << State(COD,i) << "[" << GenCodes[GenCode][i] << "] -> " << State(COD,j) << "[" << GenCodes[GenCode][j] << "] == RadMat: " <<  RadMat[(GenCodes[GenCode][i]*20)+GenCodes[GenCode][j]] << " == " << ChangeType;
				ExpObs += eqm[CurChar_i] * WorkMat[(CurChar_i*QMat->Char()) + CurChar_j];	// i -> j
//				ExpObs += eqm[CurChar_i] * *QMat->Q(CurChar_i,CurChar_j);			// i -> j
				ExpObs += eqm[CurChar_j] * WorkMat[(CurChar_j*QMat->Char()) + CurChar_i];	// j -> i
//				ExpObs += eqm[CurChar_j] * *QMat->Q(CurChar_j,CurChar_i);			// j -> i
			}
			CurChar_j++;
		}
		CurChar_i++;
	}
//	cout << "\nOverall obs for ChangeType[" << ChangeType << "]: " << ExpObs << "\n//";
	return ExpObs;

}

///////////////////////////////////////////////////////////////////////////////////////////
// General mixture model class for DrDc models
CCodonDrDcMixer::CCodonDrDcMixer(CData *Data, CTree *Tree, ECodonEqm CE, string File, int GenCode) : CCodonDrDc(Data,Tree,CE,File,GenCode)	{




}

CCodonDrDcMixer::~CCodonDrDcMixer()	{

}

// Function that prepares the codon based processes to be scaled correctly so the synonymous rate of substitution is the same throughout
void CCodonDrDcMixer::PreparelnL(bool ForceRemake)	{
//	cout << "\nNew CCodonDrDcMixer::PreparelnL";
	int i, j;
	double TotalRate = 0.0;
	vector <double> Probs, Rates;
	// Go through and apply all global parameters
	FOR(i,(int)m_vpPar.size()) { m_vpPar[i]->GlobalApply(); }
	// Prepare for likelihood computations
	FOR(i,(int)m_vpProc.size()) { m_vpProc[i]->PrepareLikelihood(true,ForceRemake,false); }		// Prepare the process Q matrices, but don't scale them yet
	// Get process probabilities
	FOR(i,(int) m_vpProc.size()) {
//		cout << "\nMatrix["<<i<<"]\n"; m_vpProc[i]->OutQ(0); cout << "\nEqm: " << m_vpProc[i]->Eqm(0);
		if(m_vpProc[i]->MaxRate()) { Probs.push_back(0.0); } else { Probs.push_back(m_vpProc[i]->Prob()); } }
	Probs = NormaliseVector(Probs);
	// Get the current rate of the individual processes (they currently scaled correctly relative to one another, just not the right mean rate)
	FOR(i,(int)m_vpProc.size()) { if(m_vpProc[i]->MaxRate()) { continue; }  Rates.push_back(m_vpProc[i]->CalcRate()); TotalRate += Rates[i] * Probs[i]; }
	// Do the rate assignment
//	cout << "\nModel: " << Name();
	FOR(i,(int)m_vpProc.size()) {
		m_vpProc[i]->Rate(Rates[i]/TotalRate);
//		cout << "\nProcess["<<i<<"] has rate: " << m_vpProc[i]->Rate() << " == " << Rates[i] / TotalRate;
		m_vpProc[i]->ScaleQ();
//		cout << "\nProcess["<<i<<"] syn rate: " << m_vpProc[i]->SynRate(true) << " and non-syn rate: " << m_vpProc[i]->NonsynRate(true);
	}

	// Output for checking
/*
	FOR(i,(int)m_vpProc.size()) {
		cout << "\n-------- Process " << i << " with prob " << Probs[i] << " = " << m_vpProc[i]->Prob() << "[AND] rate = " << Rates[i] << " == " << m_vpProc[i]->CalcRate() << " -------";
		FOR(j,m_vpProc[i]->NoPar()) { cout << "\nPar[" << j << "] " << *m_vpProc[i]->pPar(j); }
		cout << "\nQMat:";
		m_vpProc[i]->OutQ();
		cout << "\nEqm\n" << m_vpProc[i]->RootEqm();
	}
*/
//	cout << "\nAnd CCodonDrDcMixer::PreparelnL is done..."; // exit(-1);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Specific mixture models for DrDc models
CCodon2Dr1Dc::CCodon2Dr1Dc(CData *Data, CTree *Tree, ECodonEqm CE, string File, int GenCode) : CCodonDrDcMixer(Data,Tree,CE,File,GenCode) {
	int i;
	// Initialisation stuff
	m_sName = m_sName + ".Mix(2Dr1Dc)";

	// Processes
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,m_sRadicalFile));

	// Create matching between processes (conditioned on same data so eqm matches already)
	// Assumes following naming convention
	//  Kappa -> "Kappa"
	//  OmegaDr -> "Omega_Radical"
	//  OmegaDc -> "Omega_Conservative"
	// Remove parameters
	m_vpProc[1]->RemovePar("Omega_Conservative");
	m_vpProc[1]->RemovePar("Kappa");
	// Add them back
	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Kappa"));
	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Omega_Conservative"));


//	m_vpProc[0]->GetPar("Omega_Conservative")->SetVal(0.1);
//	m_vpProc[0]->GetPar("Omega_Radical")->SetVal(0.2);
//	m_vpProc[1]->GetPar("Omega_Radical")->SetVal(0.5);

	PrepareProcessProbs(true);		// Prepare the process probabilities and true == optimise them
	FinalInitialisation();

}

CCodon2Dr1Dc::~CCodon2Dr1Dc() { /* INTENTIONALLY BLANK */ }

CCodon1Dr2Dc::CCodon1Dr2Dc(CData *Data, CTree *Tree, ECodonEqm CE, string File, int GenCode) : CCodonDrDcMixer(Data,Tree,CE,File,GenCode) {
	int i;
	// Initialisation stuff
	m_sName = m_sName + ".Mix(1Dr2Dc)";

	// Processes
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,m_sRadicalFile));

	// Create matching between processes (conditioned on same data so eqm matches already)
	// Assumes following naming convention
	//  Kappa -> "Kappa"
	//  OmegaDr -> "Omega_Radical"
	//  OmegaDc -> "Omega_Conservative"
	// Remove parameters
	m_vpProc[1]->RemovePar("Omega_Radical");
	m_vpProc[1]->RemovePar("Kappa");
	// Add them back
	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Kappa"));
	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Omega_Radical"));


	// Make sure the parameters are seperate so optimisation is easy
/*	m_vpProc[0]->ProbPar()->SetVal(1.0/3.0); m_vpProc[1]->ProbPar()->SetVal(2.0/3.0);
	m_vpProc[0]->GetPar("Omega_Radical")->SetVal(0.2);
	m_vpProc[0]->GetPar("Omega_Conservative")->SetVal(0.4);
	m_vpProc[1]->GetPar("Omega_Conservative")->SetVal(0.8);
	m_vpProc[0]->GetPar("Kappa")->SetVal(2.5);
*/
	PrepareProcessProbs(true);		// Prepare the process probabilities and true == optimise them
	FinalInitialisation();

}

CCodon1Dr2Dc::~CCodon1Dr2Dc() { /* INTENTIONALLY BLANK */ }

CCodon2Dr2Dc::CCodon2Dr2Dc(CData *Data, CTree *Tree, ECodonEqm CE, string File, int GenCode) : CCodonDrDcMixer(Data,Tree,CE,File,GenCode) {
	int i;
	// Initialisation stuff
	m_sName = m_sName + ".Mix(2Dr2Dc)";

	// Processes
	// 0: Dr[0], Dc[0]
	// 1: Dr[0], Dc[1]
	// 2: Dr[1], Dc[0]
	// 3: Dr[1], Dc[1]
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,m_sRadicalFile));
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,m_sRadicalFile));
	m_vpProc.push_back(AddCodonProcess(Data,Tree,pM0DrDc,CE,GenCode,m_sRadicalFile));

	// Create matching between processes (conditioned on same data so eqm matches already)
	// Assumes following naming convention
	//  Kappa -> "Kappa"
	//  OmegaDr -> "Omega_Radical"
	//  OmegaDc -> "Omega_Conservative"
	// Remove parameters
	m_vpProc[1]->RemovePar("Omega_Radical");
	m_vpProc[3]->RemovePar("Omega_Radical");

	m_vpProc[2]->RemovePar("Omega_Conservative");
	m_vpProc[3]->RemovePar("Omega_Conservative");

	m_vpProc[1]->RemovePar("Kappa");
	m_vpProc[2]->RemovePar("Kappa");
	m_vpProc[3]->RemovePar("Kappa");
	// Add them back
	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Omega_Radical"));
	m_vpProc[3]->AddQPar(m_vpProc[2]->GetPar("Omega_Radical"));

	m_vpProc[2]->AddQPar(m_vpProc[0]->GetPar("Omega_Conservative"));
	m_vpProc[3]->AddQPar(m_vpProc[1]->GetPar("Omega_Conservative"));

	m_vpProc[1]->AddQPar(m_vpProc[0]->GetPar("Kappa"));
	m_vpProc[2]->AddQPar(m_vpProc[0]->GetPar("Kappa"));
	m_vpProc[3]->AddQPar(m_vpProc[0]->GetPar("Kappa"));


	PrepareProcessProbs(true);		// Prepare the process probabilities and true == optimise them
	FinalInitialisation();

}

CCodon2Dr2Dc::~CCodon2Dr2Dc() { /* INTENTIONALLY BLANK */ }






