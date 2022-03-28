#include <iostream>
#include <fstream>
#include <sstream>
//#include <sys/stat.h>

#include "analyse_hfs.h"

void analyse_hfs(const char *spinState, double eGateMin, double eGateMax)//(int Mass, int nRun1, int nRun2)
{
	/***************************************************
	*   Open files and read in data from laser file    *
	***************************************************/
	// Check if file is open and load trees
	if(!loadTrees()) return;
	loadTrees();

	// Get file namebase
	TFile *f = (TFile*)gROOT->GetFile();
	string nameBase = f->GetName();
	nameBase.erase(nameBase.find_last_of("."), string::npos);

	// Open laser file from RILIS and output file for results
	ostringstream laserFile, datFile;
	laserFile << "../bin/"  << nameBase << ".txt";
	datFile   << "hfs_dat/" << nameBase << "_" << spinState << ".dat";
	std::string::size_type sz;
	int Mass = stoi(nameBase.substr(0,3),&sz);
	int nRun = stoi(nameBase.substr(nameBase.length() - 4),&sz);

	// Check if directory for output file exists, if not, make it
	struct stat buffer;
	if(stat("hfs_dat",&buffer)!=0) system("mkdir hfs_dat");

	// Check if RILIS laser file exists
	ifstream lFile;
	lFile.open(laserFile.str());
	if(!lFile.is_open())
	{
		cout << "Could not open laser file: " << laserFile.str() << "\nExiting program\n";
		return;
	}

	// Check the output file has been made 
	ofstream dFile;
	dFile.open(datFile.str());
	if(!dFile.is_open())
	{
			cout << "Could not open .dat file: " << datFile.str() << "\nExiting program\n";
			return;
	}

	// Output the laser file anme to read in, and the name of the output data file.
	cout << "\n  ***********************************************************************\n";
	cout <<   "  * Reading laser file:  " << laserFile.str() << endl;
	cout <<   "  * Output results to:   " << datFile.str() << endl;
	cout <<   "  * Alpha energy gates:  " << eGateMin << " to " << eGateMax << endl;
	cout <<   "  ***********************************************************************\n\n";

	/***************************************************
	*                   Extract HFS                    *
	***************************************************/
	string headerLine;
	getline(lFile,headerLine);

	// Read in data from laser file, check the data format in file matches that which is trying to be read in
	while(lFile >> dS >> dS >> dI >> dD >> _freq >> dI >> _protons >> _power >> dS >> dS)
	{
		   freq.push_back(_freq);
		protons.push_back(_protons);
		  power.push_back(_power);		
	}

	// Column headers for output file
	dFile << "Wavenumber\tAlpha counts\tLaser power\tProton current" << endl;

	int steps = freq.size();
	double progress = 0.0;
	double pStep = 100000/steps;

	cout << "  Extracting HFS:\n";

	// Go through root file, gate on each laser step with alpha energy gate,
	// output how many alphas at that laser step, the wavenumber setting and the integrated
	// proton current and laser power at that step to .dat file
	for (int i=0; i<(int)freq.size(); i++) 
	{
		ostringstream gateCondition;
		gateCondition << "((SEng[1]>" << eGateMin << "&&SEng[1]<" << eGateMax << ")||(SEng[2]>" << eGateMin << "&&SEng[2]<" << eGateMax << "))&&(SsCycle==" << i+1 << ")";
		counts.push_back(hf->GetEntries(gateCondition.str().c_str()));

		if (freq[i] != 0)
		{
			dFile << fixed << setprecision(6) << freq[i] << "\t" << fixed << setprecision(0) << counts[i] << "\t\t" << fixed << setprecision(6) << power[i] << "\t" << protons[i] << endl;
		}


		//Progress bar!!!
		int barWidth = 50;
		if(i % 5 == 0)
		{
			cout << "   [";
			Float_t pos = barWidth*progress/100000;
			for(int j=0; j<barWidth; j++)
			{
				if(j<pos) {cout << "=";}
				else if(j==pos) {cout << ">";}
				else {cout << " ";}
			}
			cout << left << "] " << fixed << setprecision(0) << right << 100*progress/100000 << " %     \r ";
			cout.flush();
		}
		progress += pStep;

	}
	dFile.close();

	cout << "\n  Finished extracting HFS\n";

	/******************************************************
	*                    Make plots                       *
	******************************************************/
	TCanvas *aHFS = new TCanvas("aHFS","aHFS",825,915);
	aHFS->Divide(1,2);
	aHFS->cd(1);
	TH1D *se = new TH1D("se","se",1.2e3,6.3e3,7.5e3);
	hf->Draw("SEng>>se","SEng[1]>5e3||SEng[2]>5e3");
	TLine *xl = new TLine();
	xl->SetLineColor(2);
	xl->SetLineWidth(2);
	xl->DrawLine(eGateMin,0,eGateMin,se->GetMaximum());
	xl->DrawLine(eGateMax,0,eGateMax,se->GetMaximum());

	aHFS->cd(2);
	vector<double> eCounts;
	for(int i=0; i<(int)freq.size(); i++) eCounts.push_back(sqrt(counts[i]));
	TGraphErrors *hfs = new TGraphErrors(freq.size(),&(freq[0]),&(counts[0]),0,&(eCounts[0]));
	hfs->SetMarkerStyle(7);
	hfs->Draw("AP");

	TCanvas *lp = new TCanvas("lp","lp",825,915);
	TPad *sp[5];
	sp[0] = new TPad("sp0","sp0",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp1","sp1",0.0,2./3,1.0,1.0);
	sp[2] = new TPad("sp2","sp2",0.0,1./3,1.0,2./3);
	sp[3] = new TPad("sp3","sp3",0.0,0.0,1.0,1./3);
	sp[1]->SetBottomMargin(0);
	sp[2]->SetBottomMargin(0);
	sp[2]->SetTopMargin(0);
	sp[3]->SetTopMargin(0);

	lp->cd();
	sp[0]->Draw();
	for(int i=0; i<3; i++)
	{
		sp[0]->cd();
		sp[i+1]->Draw();
	}

	sp[1]->cd();
	hfs->Draw("AP");
	sp[2]->cd();
	TGraph *lPower = new TGraph(freq.size(),&(freq[0]),&(power[0]));
	lPower->SetMarkerStyle(7);
	lPower->GetYaxis()->CenterTitle();
	lPower->GetYaxis()->SetTitle("Laser power [Watts]");
	lPower->GetYaxis()->SetTitleSize(0.06);
	lPower->GetYaxis()->SetTitleOffset(0.8);
	lPower->Draw("AP");
	lPower->SetMinimum(0);
	sp[3]->cd();
	TGraph *nProtons = new TGraph(freq.size(),&(freq[0]),&(protons[0]));
	nProtons->SetMarkerStyle(7);
	nProtons->GetYaxis()->CenterTitle();
	nProtons->GetYaxis()->SetTitle("No. of protons");
	nProtons->GetXaxis()->CenterTitle();
	nProtons->GetXaxis()->SetTitle("Wavenumber [cm^{-1}]");
	nProtons->GetYaxis()->SetTitleSize(0.06);
	nProtons->GetYaxis()->SetTitleOffset(0.8);
	nProtons->Draw("AP");
	nProtons->SetMinimum(0);

	ostringstream ac, pc;
	ac << "hfs_dat/alpha_hfs_" << nRun << "_" << spinState << ".pdf";
	pc << "hfs_dat/power_hfs_" << nRun << "_" << spinState << ".pdf";
	aHFS->SaveAs(ac.str().c_str());
	lp->SaveAs(pc.str().c_str());

	freq.clear();
	counts.clear();
	power.clear();
	protons.clear();
}
