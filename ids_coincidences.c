#include "ids_coincidences.h"

// Script for producing coincidence matrices, written during
// IS665 run (2021) by J.G. Cubiss.
// Use config file "config_coincidences.cfg" to set prompt
// windows and to choose which histos are created
// [1=create, 0=do not create]

//void ids_coincidences(const char* coinc)
void ids_coincidences(int optionAlphaGamma, int optionGammaGamma, int optionBetaGammaGamma)
{
	// Get conditions for different histos depending on what
	// coincidence options were selected
	if(optionAlphaGamma==1)
	{
		cout << "  Alpha-gamma prompt window?\n";
		//cin  >> promptWindowAlphaGamma;
		promptWindowAlphaGamma=55;
	}


	if(optionGammaGamma==1||optionBetaGammaGamma==1) 
	{
		cout << "  Gamma-gamma prompt window?\n";
		//cin  >> promptWindowGammaGamma;
		promptWindowGammaGamma=55;

		if(optionBetaGammaGamma==1)
		{
			cout << "  Lower limit beta energy cut?\n";
			//cin  >> betaEnergyMin;
			betaEnergyMin=20;
			cout << "  Upper limit beta energy cut?\n";
			//cin  >> betaEnergyMax;
			betaEnergyMax=2000;
		}
	}

	// Output to terminal what will be drawn
	outputConfiguration(optionAlphaGamma,optionGammaGamma,optionBetaGammaGamma);

	// Load tree
	TTree *ids = (TTree*)gROOT->GetFile()->Get("ids");
	int nEntries = ids->GetEntries();

	// Create histogram objects
	if(optionAlphaGamma==1)      ag_mat = new TH2D("ag_mat","ag_mat",5e3,0,10e3,8e3,0,8e3);
	if(optionGammaGamma==1)	     gg_mat = new TH2D("gg_mat","gg_mat",5e3,0,5e3,5e3,0,5e3);
	if(optionBetaGammaGamma==1) bgg_mat = new TH2D("bgg_mat","bgg_mat",5e3,0,5e3,5e3,0,5e3);

	// Set branch addresses of input tree to variables
	ids->SetBranchAddress("Multiplicity",&MULT);
	ids->SetBranchAddress("Time_vs_ref",&TIME_REF);
	ids->SetBranchAddress("Timestamp",&TIMESTAMP);
	ids->SetBranchAddress("Energy_Clov",&E_Clov);
	ids->SetBranchAddress("Time_Clov",&T_Clov);
	ids->SetBranchAddress("Mult_Clov",&M_Clov);
	ids->SetBranchAddress("Energy_Beta",&E_Beta);
	ids->SetBranchAddress("Time_Beta",&T_Beta);
	ids->SetBranchAddress("Mult_Beta",&M_Beta);
	ids->SetBranchAddress("Energy_Si",&E_Si);
	ids->SetBranchAddress("Time_Si",&T_Si);
	ids->SetBranchAddress("Mult_Si",&M_Si);
	ids->SetBranchAddress("Energy_Proton",&E_Proton_T1);
	ids->SetBranchAddress("Time_Proton",&T_Proton_T1);
	ids->SetBranchAddress("Mult_Proton",&M_Proton_T1);

	// Loop over all entries
	for(int i=0; i<nEntries; i++)
	{
		ids->GetEntry(i);
		
		// Alpha-gamma matrix
		if(optionAlphaGamma==1)
		{
			// Check if a-g event
			if(M_Clov>0 && M_Si>0)
			{
				// Loop over all deetector and clover crystal combinations
				// to check if time condition meets prompt condition, if so,
				// add event to matrix
				for(int j=0; j<detectors; j++)
					for(int k=0; k<crystals; k++)
					{
						if(E_Si[j]>10 && E_Clov[k]>10 && abs(T_Si[j]-T_Clov[k])<promptWindowAlphaGamma )
						{
							ag_mat->Fill(E_Si[j],E_Clov[k]);
						}
					}
			}
		}

		// Same as alpha-gamma but for gamma-gamma events, or beta-gamma-gamma
		if(optionGammaGamma==1 || optionBetaGammaGamma==1)
		{
			if(M_Clov>1)
			{
				for(int j=0; j<crystals; j++)
					for(int k=0; k<crystals; k++)
					{
						if(j!=k && E_Clov[j]>10 && E_Clov[k]>10 && abs(T_Clov[j]-T_Clov[k])<promptWindowGammaGamma )
						{
							if(optionGammaGamma==1) gg_mat->Fill(E_Clov[j],E_Clov[k]);
							for(int l=0; l<detectors; l++) 
								if(optionBetaGammaGamma==1 && (M_Beta>0 || (M_Si>0 && E_Si[l]>betaEnergyMin && E_Si[l]<betaEnergyMax )) ) bgg_mat->Fill(E_Clov[j],E_Clov[k]);
						}
					}
			}
		}
	}

	// Draw histos and reset options
	if(optionAlphaGamma==1)
	{
		TCanvas *ag = new TCanvas();
		ag_mat->Draw("colz");
	}
	if(optionGammaGamma==1)
	{
		TCanvas *gg = new TCanvas();
		gg_mat->Draw("colz");
	}
	if(optionBetaGammaGamma==1)
	{
		TCanvas *bgg = new TCanvas();
		bgg_mat->Draw("colz");
	}
}

// Output info on what is being made to terminal
void outputConfiguration(int optionAlphaGamma, int optionGammaGamma, int optionBetaGammaGamma)
{
	string fileName;
	fileName = gROOT->GetFile()->GetName();
	cout << "\n\n  *****************************************************\n";
	cout << "  Building coincidences for " << fileName  << ":\n";
	if(optionAlphaGamma==1) cout << "  Alpha-Gamma coincidences, Prompt window: " << promptWindowAlphaGamma << "\n";
	if(optionGammaGamma==1||optionBetaGammaGamma) cout << "  Gamma-Gamma coincidences, Prompt window: " << promptWindowGammaGamma << "\n";
	if(optionBetaGammaGamma==1)
	{
		cout << "  Beta cut min: " << betaEnergyMin << "\n";
		cout << "  Beta cut max: " << betaEnergyMax << "\n";
	}
	cout << "  *****************************************************\n\n";
}
