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
CSiteCodon::CSiteCodon(CData *D, CTree *T, vector <int> ModelPar, vector <int> BranchPar, EModel CoreModel) : CBaseModel(D,T)	{
	int i;
	CData *TempData;
	vector <CTree *> vTempTrees;
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
	// Create the tree objects
	vTempTrees.assign(3,NULL); FOR(i,3) { vTempTrees[i] = new CTree(); *vTempTrees[i] = *T; }
	m_vpTreeSites.assign(3,NULL);
	FOR(i,3) { m_vpTreeSites[i] = vTempTrees[BranchPar[i]]; }
	FOR(i,3) { vTempTrees[i] = NULL; }
	// Create the model objects
	FOR(i,3) { m_vpAssociatedModels.push_back(GetMyModel(CoreModel,m_vpDataSites[i],m_vpTreeSites[i])); }

	CBaseModel *Tester = GetMyModel(CoreModel,D,T);

	cout << "\nComplete model" << Tester->lnL(true);
	double Summer = 0;
	cout << "\nChecking individual submodels...";
	FOR(i,3) {
		Summer += m_vpAssociatedModels[i]->lnL(true);
		cout << "\nModel["<<i<<"]: " << m_vpAssociatedModels[i]->lnL(true);
	}
	cout << "\nTotal lnL: " << Summer;

	/*
	 * Some issues that need to be addressed when establishing models
	 * ---
	 * 1.The equilibrium distributions from the processes is currently wacky because each instancing of the model gets its own eqm. Need to edit the data objects so the reflect the model assumptions. Can take the sum across codon positions to be consistent
	 * 		General thought: Want to set it up so eqm is never touched so fix here and have an assert during likelihood computation
	 * 2. The likelihood function needs at least a wrapper function that ensures consistency of parameters between models that are meant to share parameters
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

