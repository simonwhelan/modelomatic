/*	///////////////////////////////////////////////////////////////////////
			Definitions of the interface for Leaphy
	/////////////////////////////////////////////////////////////////////// */

#ifndef LEAPHY_INTERFACE
#define LEAPHY_INTERFACE

#include "interface.h"


void ModelMenu(CPhyloDat *PhyDat)	{
	bool DoMenu = true;
	int i,j,k,count;
	// Default Model settings
	EModel Model = PhyDat->Model();
	bool DoGamma = PhyDat->DoGam();
	int GamCat = PhyDat->NoGamCat();
	bool DoInv = PhyDat->DoInv();
	bool DoBioNJ = PhyDat->bionj();
	bool DoSA = PhyDat->bionj();

	// Organise the data before anything else
	PhyDat->GetData();
	PhyDat->GetOutputFile();

	if(Model == UNKNOWN) { if(PhyDat->DataType() == DNA) { Model = DEFAULT_DNA_MODEL; } else if (PhyDat->DataType() == AA) { Model = DEFAULT_AA_MODEL; } else { Error("Unknown data type...\n"); } }

	while(DoMenu)	{
		if(system("clear") == 1) { if(system("cls") == 1) { cout << "\nApologies for messy output -- unknown prompt..."; }; } ;
		cout << "\n\t==================================================================";
		cout << "\n\t                       " << PROG_NAME << " (version " << VERSION_NUM << ")";
		cout << "\n\t==================================================================";
		cout << "\n\n\tWorking with " << PhyDat->DataType() << " data from file <"<<PhyDat->In()<<">:\n\t\t" << PhyDat->NoSeq() << " sequences of length " << PhyDat->TrueSize() << " [unique_columns=" << PhyDat->Size() << "]";
		cout << "\n\n\tOutput file: " << PhyDat->Out();
		cout << "\n\n\tMenu:";
		if(Model == THMM_FULLDNA)	{ cout << "\n\t\t1. Model:                 temporal hidden Markov Model (THMM" << PhyDat->THMM_Details() << ")"; } 
		else { cout << "\n\t\t1. Model:                 " << Model; if(PhyDat->DoF()) { cout << "+F"; } }
		cout << "\n\t\t2. Gamma rate variation:  ";
		if(DoGamma)	{ cout << "Yes (" << GamCat << " categories)"; } else { cout << "No"; }
		cout << "\n\t\t3. Invariant sites:       ";
		if(DoInv) { cout << "Yes"; } else { cout << "No"; }
		cout << "\n\n\t\t4. Starting tree:         ";
		if(PhyDat->IsTree()) { 
			if(PhyDat->TreeFile().empty()) { Error("Problem with treefile input..."); }
			cout << "from file <" << PhyDat->TreeFile() << ">";
		} else {
			if(DoBioNJ && DoSA) { cout << "Best of bionj & stepwise addition"; }
			else if(DoBioNJ) { cout << "bionj"; } else if(DoSA) { cout << "stepwise addition"; }
		}
		cout << "\n\n\t\t0. Perform analysis\n\n";
		

		switch(GetOption("\n\tPlease choose an option",0,4))	{
		case 0: DoMenu = false; break;
		case 1:
			cout << "\n\n\tModels available: ";
			if(PhyDat->DataType() == DNA)			{ i = 0; j = NoDNAModels; }
			else if(PhyDat->DataType() == AA)		{ i = NoDNAModels; j = NoDNAModels + NoAAModels; }
			count = 1; for(k=i;k<j;k++)	{ cout << "\n\t\t" << count++ << ": " << (EModel) k; }
			Model = (EModel) (GetOption("\nPlease choose a model",1,j-i) + i - 1); 
			if(PhyDat->DataType() == AA) { PhyDat->SetDoF(GetYesNoOption("\nUse observed equilibrium distribution ")); }
			if(Model == THMM_FULLDNA)	{
				cout << "\nTemporal HMMs require more information:"; 
				PhyDat->GetTHMMStdIn();
			}
			break;
		case 2:
			DoGamma = GetYesNoOption("Would you like the model to have gamma distributed among site rate variation");
			if(DoGamma) { GamCat = GetOption("How many rate categories?",1,10); }
			break;
		case 3:
			DoInv = GetYesNoOption("Would you like the model to have some invariant sites");
			break; 
		case 4: 
			if(GetYesNoOption("Get tree from file"))	{
				if(PhyDat->IsTree()) { PhyDat->CleanTree(); }
				PhyDat->SetTreeFile(GetInFileName());
				PhyDat->GetTree();
			} else {
				DoSA = DoBioNJ = false;
				while(!DoSA && !DoBioNJ)	{
					DoSA = GetYesNoOption("Do stepwise addition starting tree");
					DoBioNJ = GetYesNoOption("Do bionj starting tree");
					if(!DoSA && !DoBioNJ) { cout << "\nNeed to get a starting tree from somewhere...\n"; }
				}
			}
			break;
		default:	Error("Unknown option in menu");
		}
	}

	PhyDat->SetBioNJ(DoBioNJ); PhyDat->SetSA(DoSA);
	PhyDat->SetModel(Model); PhyDat->SetInv(DoInv);
	PhyDat->SetGamma(DoGamma); PhyDat->SetGamCat(GamCat);
}


#endif
