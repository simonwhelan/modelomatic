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
	vTempTrees.assign(3,NULL); FOR(i,3) { vTempTrees[i] = new CTree(); *vTempTrees[i] = *T; }
	m_vpTreeSites.assign(3,NULL);
	FOR(i,3) { m_vpTreeSites[i] = vTempTrees[BranchPar[i]]; }
	FOR(i,3) { vTempTrees[i] = NULL; }
	// Create the model objects
	m_viModelMap = ModelPar;
	FOR(i,3) { m_vpAssociatedModels.push_back(GetMyModel(CoreModel,m_vpDataSites[i],m_vpTreeSites[i])); }
	NormaliseParameters();


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
	FOR(i,(int) m_vpAssociatedModels.size()) {
		logL += m_vpAssociatedModels[i]->lnL(ForceReal);
	}
	return logL;
}

/////////////////////////////////////////////
// Override the DoBraOpt function
void CSiteCodon::DoBraOpt(int First, int NTo, int NFr, int Br, bool IsExtBra,double *BestlnL,double tol,bool AllowUpdate) {
	NormaliseParameters();
	cout << "\nHaven't sorted DoBraOpt just yet... I think I have to be smarter than I am at the moment..."; exit(-1);
}



