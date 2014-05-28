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

// The model will be set up inefficiently to start with. The main CBaseModel object will have nothing to optimise in it
CSiteCodon::CSiteCodon(CData *D, CTree *T, vector <int> ModelPar, vector <int> BranchPar, EModel CoreModel, bool Gamma) : CBaseModel(D,T)	{
	int i,j,k;
	CData *TempData;
	vector <CTree *> vTempTrees;
	vector <double> FreqCount;
	string t_string;
	// Some error checking on entry
	assert(ModelPar.size() == 3); assert(BranchPar.size() == 3);
	FOR(i,3) { assert(ModelPar[i] <= i); assert(InRange(ModelPar[i],0,3)); assert(BranchPar[i] <= i); assert(InRange(BranchPar[i],0,3)); }
	assert(InRange(CoreModel,JC,RY_model)); // Check it's a DNA model

	cout << "\nWorking with core model: " << CoreModel;
	cout << "\nInitialising okay..." << endl << ModelPar << endl << BranchPar;

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
		FreqCount.clear(); FreqCount.assign(4,0.0);
		FOR(j,3) { if(ModelPar[j] == i) { FOR(k,4) { FreqCount[k] += m_vpDataSites[j]->m_vFreq[k]; } } }
		if(Sum(&FreqCount) < FLT_EPSILON) { continue; }
		FreqCount = NormaliseVector(FreqCount);	// Safe to just average the frequencies because the counts will be the same
		FOR(j,3) { if(ModelPar[j] == i) { FOR(k,4) { m_vpDataSites[j]->m_vFreq[k] = FreqCount[k];  } } }
	}
	// Create the tree objects
	vTempTrees.assign(3,NULL); FOR(i,3) {
		vTempTrees[i] = new CTree(); *vTempTrees[i] = *T;
		t_string = ""; FOR(j,3) { if(BranchPar[j] == i) { t_string = t_string + int_to_string(j+1); } }
		FOR(j,vTempTrees[i]->NoBra()) {  vTempTrees[i]->pBra(j)->Name("site[" + t_string + "]::" + vTempTrees[i]->pBra(j)->Name()); }
	}
	m_vpTreeSites.assign(3,NULL);
	FOR(i,3) { m_vpTreeSites[i] = vTempTrees[BranchPar[i]]; }	// Could potentially cause a memory leak, but it's small fry... (I hope!)
	FOR(i,3) { vTempTrees[i] = NULL; }
	// Create the model objects
	m_viModelMap = ModelPar;
	m_viTreeMap = BranchPar;

	FOR(i,3) {
		m_vpAssociatedModels.push_back(GetMyModel(CoreModel,m_vpDataSites[i],m_vpTreeSites[i]));
		t_string = ""; FOR(j,3) { if(m_viModelMap[j] == m_viModelMap[i]) { t_string = t_string + int_to_string(j+1); } }
		FOR(j,(int)m_vpAssociatedModels[i]->m_vpPar.size()) { m_vpAssociatedModels[i]->m_vpPar[j]->Name("site[" + t_string + "]::" + m_vpAssociatedModels[i]->m_vpPar[j]->Name()); }

	}
	NormaliseParameters();

	cout << "\nGrabbing parameters";
	vector <double*> ParVals = GetOptPar(true,true,true,false);
	cout << "\nParameters are: ";
	FOR(i,ParVals.size()) {
		cout << "\nPar["<<i<<"] == " << *ParVals[i] << " : " << *m_vpAllOptPar[i];
	}


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
	NormaliseParameters();
	assert(m_vpAssociatedModels.size() == 3);
	if(m_vbUseInBraCalc.empty()) {
		FOR(i,(int) m_vpAssociatedModels.size()) { logL += m_vpAssociatedModels[i]->lnL(ForceReal); }
	} else {	// Should only happen for FastBranch calculations. Here for debugging
		assert(m_vbUseInBraCalc.size()==3);
		FOR(i,(int) m_vpAssociatedModels.size()) { if(m_vbUseInBraCalc[i]) { logL += m_vpAssociatedModels[i]->lnL(ForceReal); } }
	}
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
	double working_tol, BestlnL, newlnL;
cout << "\nCSiteCodon::FastBranchOpt(...)";

	// Checking everything okay going in
	assert(m_vbUseInBraCalc.empty()); assert(m_vpAssociatedModels.size() == 3);
	assert(m_viModelMap.size() ==3 && m_viTreeMap.size() == 3); FOR(i,3) { assert(InRange(m_viTreeMap[i],0,3)); }
	// Loop through the sites and see what branches need to be optimised
	NormaliseParameters();
	cout << "\nTreeMap: " << m_viTreeMap;
	FOR(site,3) {
		cout << "\n\tDoing branch set["<<site<<"] == " << m_viTreeMap[site];
		// Check whether the trees been done yet
		if(TreeDone[m_viTreeMap[site]]) { cout << " ... skipping"; continue; } TreeDone[m_viTreeMap[site]] = true;
		// If not, find the other processes associated with it. Also initialises the model ready for computation
		m_vbUseInBraCalc.assign(3,false); CurlnL = 0;
		FOR(j,3) {
			if(m_viTreeMap[site] == m_viTreeMap[j]) {
				m_vbUseInBraCalc[j] = true; CurlnL += m_vpAssociatedModels[j]->lnL(true);
				FOR(k,m_vpAssociatedModels[site]->m_vpProc.size()) { m_vpAssociatedModels[site]->m_vpProc[k]->PrepareBraDer(); }
		}	}
		if(tol > 1.0E-3) { working_tol = tol; } else { working_tol = 1.0E-3; }

		cout << "\n\t" << m_vbUseInBraCalc << " == lnL: " << CurlnL;

		///////////////////////////////////////////////////////////////////
		// Only do cyclical optimisation with multiple branches
		if(Tree()->NoBra() == 1) {
			DoBraOpt(true,0,1,0,true,&CurlnL,tol,false);
			continue;
		}
		FOR(i,NoIter)	{
			BestlnL = newlnL = CurlnL;							// 1. Do the first calculation
//	#if FASTBRANCHOPT_DEBUG == 1
	//#if DEVELOPER_BUILD == 1
			cout << "\n\n--- Round " << i<< ". " << newlnL << " (tol: "<< working_tol << ") ---";;
			cout << "\nOriginal branches:  "; int j; FOR(j,Tree()->NoBra()) { cout << Tree()->B(j) << " "; }
			cout << flush;
//	#endif
			BranchOpt(-1,m_vpAssociatedModels[site]->Tree()->StartCalc(),-1, &BestlnL,working_tol);	// 2. Run the fast optimisation routine
			BestlnL = lnL();
			if(working_tol > tol) { working_tol = max(tol,working_tol/10); }
	//#if DEVELOPER_BUILD == 1
//	#if FASTBRANCHOPT_DEBUG == 1
			cout << "; " << BestlnL << " == " << lnL() << "; diff = " << BestlnL - lnL();
			cout << "\nTree: " << *Tree();
			if(fabs(BestlnL - lnL()) > tol) { cout << "\nBig Error..."; exit(-1); }
//	#endif
			if(fabs(BestlnL - newlnL) < tol) { break; }				// 3. Control exit
		}
		if(Conv != NULL) { if(i==NoIter) { *Conv = false; } else { *Conv = true; } }
	//	cout << "\nReturning: " << BestlnL << " cf. " << lnL() << " fabs: " << fabs(BestlnL - lnL()); // exit(-1);
		assert(BestlnL < 0);
		// Finish up for this subset of data
		m_vbUseInBraCalc.clear();
	}
	CurlnL = 0; FOR(i,3) { CurlnL += m_vpAssociatedModels[i]->lnL(true); } 	// There's a better way to get this number without recalculating everything
	cout << "\nDone with CSiteCodon::FastBranchOpt(...) -- Returning: " << CurlnL;
	return CurlnL;
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
	NormaliseParameters();
	int i;
	double *p_x,x1,x2,x3,x1_lnL = 1.0,x2_lnL = -fabs(*BestlnL),x3_lnL = 1.0, xi,temp;
	tol = max(tol,FULL_LIK_ACC);
	double dx = DX;
	CPar *Par;
	// Some stuff to check at the beginning
	assert(m_vpAssociatedModels.size() == 3 && m_viTreeMap.size() == 3 && m_vbUseInBraCalc.size() == 3);
	FOR(i,3) { if(m_vbUseInBraCalc[i]) { p_x = m_vpAssociatedModels[i]->Tree()->OptimiserB(Br); Par = m_vpAssociatedModels[i]->Tree()->pBra(Br); break; } }
	assert(i!=3);
	cout << "\nBranch["<<Br<<"] has DoBralnL: "<< DoBralnL(Br,NFr,NTo) << " cf. " << DoBralnL(Br,NFr,NTo); cout << " Using parameter: " << *Par << " == " << *p_x;
//	cout << "\nReturning from CBaseModel::DoBraOpt (including branch updates)";
//	RETURN_DOBRAOPT;
//	cout << "\nDoing branch["<<Br<<"]: sent bestlnL: " << *BestlnL << " cf. DoBralnL(Br,NFr,NTo): " << DoBralnL(Br,NFr,NTo); //  << " cf. real " << lnL(); exit(-1);
	*BestlnL = DoBralnL(Br,NFr,NTo);
	x2 = *p_x;
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
		if(!Par->CheckLowBound()) { break; }
	}
	// ------------------------------------- Catch for when lower bounds branch length ------------------------------------
	if(fabs(x1 - x2) < DX) {
		*p_x = x2 = DX; x2_lnL = DoBralnL(Br,NFr,NTo); m_iFastBralnL_Bracket++;
		if(x1_lnL > x2_lnL) { *p_x = x1; 	*BestlnL = x1_lnL; RETURN_DOBRAOPT; }
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
	cout << "\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
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
	cout << "\nInto Goldensection";
	FOR(i,20) {
		// New value
		if(x3 - x2 > x2 - x1) 	{ *p_x = xi = x2 + resphi * (x3-x2); }
		else					{ *p_x = xi = x2 - resphi * (x2-x1); }
		// break condition
		if(fabs(x3_lnL - x1_lnL) < tol) { /* cout << "\nBreaking at tol=" << tol << " fabs(" << x3_lnL << " - " << x1_lnL << ")";  */ *p_x = x2 = (x1+x3)/2; break; }
		// Search
		temp = DoBralnL(Br,NFr,NTo);
		cout << "\n\t[i="<<i<<"] xi:" << xi << ": " << temp;
		if(temp > x2_lnL) {
			if(x3 - x2 > x2 - x1) 	{ x1 = x2; x1_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
			else					{ x3 = x2; x3_lnL = x2_lnL;  x2 = xi; x2_lnL = temp; }
		} else {
			if(x3 - x2 > x2 - x1)	{ x3 = xi; x3_lnL = temp; }
			else					{ x1 = xi; x1_lnL = temp; }
		}
		cout << " [x1:" << x1 << ", x2:" << x2 << ", x3:" << x3 << ","<< max(x1_lnL,max(x2_lnL,x3_lnL)) << "]" << flush;
	}

/*
	cout << "\nFinished search: x: " << *p_x << " = " << temp << " == " << DoBralnL(Br,NFr,NTo);
	cout << ": tol= " << max(x2_lnL - x1_lnL,x2_lnL - x3_lnL);
	cout << "\n---\nx1: " << x1 << " == " << x1_lnL << " (diff="<<x2_lnL - x1_lnL << ")";
	cout << "\nx2: " << x2 << " == " << x2_lnL << " (diff="<<x2_lnL - x2_lnL << ")";
	cout << "\nx3: " << x3 << " == " << x3_lnL << " (diff="<<x2_lnL - x3_lnL << ")";
*/	// Finish by doing the calculation again to correctly update the partial likelihoods
	m_iFastBralnL_Calls++;
	*BestlnL = DoBralnL(Br,NFr,NTo);
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
		RetlogL += logL;
		// Clear up after the process
		FOR(i,CProbSize) { delete P[i]; } delete [] P;
	}
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
	// Clean the parameter space
	FOR(i,(int)m_vpAllOptPar.size()) { m_vpAllOptPar[i] = NULL; } m_vpAllOptPar.clear();
	m_vbDoBranchDer.clear();
	// Loop through the individual site models collecting the necessary parameters as we go
	FOR(site,3) {
		// Work out what needs to be collected
		ProcEBra = ExtBranch; ProcIBra = IntBranch; ProcPar = Parameters; ProcEqm = Eqm;
		if(TreeDone[m_viTreeMap[site]] == false) {	TreeDone[m_viTreeMap[site]] = true; } else { ProcEBra = false; ProcIBra = false; }	// Only collect first instance of branches/parameters
		if(ParDone[m_viModelMap[site]] == false) {	ParDone[m_viTreeMap[site]] = true; } else { ProcPar = false; ProcEqm = false; }
		// Get the values
		PlaceHold = m_vpAssociatedModels[site]->GetOptPar(ProcEBra,ProcIBra,ProcPar,ProcEqm);
		assert(PlaceHold.size() == m_vpAssociatedModels[site]->m_vpAllOptPar.size());
		FOR(i,(int)PlaceHold.size()) { OptVal.push_back(PlaceHold[i]); PlaceHold[i] = NULL; m_vpAllOptPar.push_back(m_vpAssociatedModels[site]->m_vpAllOptPar[i]); } PlaceHold.clear();
		cout << "\nProcess["<<site<<"]: totalpar = " << OptVal.size();


	}
	/*	cout << "\nGetting Opt Par ["<< m_vpAllOptPar.size();
		cout << ":" << OptVal.size() <<"]";
		FOR(i,m_vpAllOptPar.size()) {
			cout << "\n\tPar["<<i<<"] " << m_vpAllOptPar[i]->Name() << " = " << m_vpAllOptPar[i]->Val() << " == " << *OptVal[i];
		}
		cout << "\n---------";
	*/
	return OptVal;
}


