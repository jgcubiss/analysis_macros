#include "TVirtualFitter.h"

// You should be using at least version
// 6.04 of root for this.
//
// You need a config file with two columns,
// 1st column for the energies of your peak,
// the second for a rough estimate of their
// heights. The code will read in the values
// in the config file and fit however many
// peaks there are values for and output the
// individual peak integrals to terminal.

// Declare functions
double CBeff(double *x, double *par);

// Global variable for number of peaks
void effLege(double eGamma)
//void CB_LEGe()
{
//	Double_t eGamma=150;
//	TCanvas *c = new TCanvas("c","c",750,550);

	TGraphErrors *effGraph = new TGraphErrors("/home/jgc/macros/LEGe_eff.dat","%lg %lg %lg");
	effGraph->SetMarkerStyle(20);
//	effGraph->Draw("AP");
	effGraph->GetYaxis()->SetRangeUser(0.05,9.);
	effGraph->GetXaxis()->SetRangeUser(0,3000);
	effGraph->SetTitle(0);

	effGraph->GetXaxis()->SetTitle("E_{#gamma} [keV]");
	effGraph->GetXaxis()->SetLabelSize(0.03);
	effGraph->GetXaxis()->CenterTitle();

	effGraph->GetYaxis()->SetTitle("Absolute Efficiency [%]");
	effGraph->GetYaxis()->CenterTitle();
	effGraph->GetYaxis()->SetRangeUser(0,10);

	TF1 *eff_fit = new TF1("eff_fit",CBeff,0,1500,5);
	eff_fit->SetNpx(10000);

	// Parameter limits for le fit
	double lim_Height[3] = {9, 8, 10};		// height
	double E_tol[3]      = {80,70,90};			// +/- E est
	double lim_n[3]      = {14, 1, 100};	// n
	double lim_Alpha[3]  = {-1.65, -12, 0};	// alpha
	double lim_Sigma[3]  = {108, 20, 150};		// sigma

	eff_fit->SetParName(0,"alpha");
	eff_fit->SetParameter(0,lim_Alpha[0]);
	eff_fit->SetParLimits(0,lim_Alpha[1],lim_Alpha[2]);

	eff_fit->SetParName(1,"n");
	eff_fit->SetParameter(1,lim_n[0]);
	eff_fit->SetParLimits(1,lim_n[1],lim_n[2]);

	eff_fit->SetParName(2,"sig");
	eff_fit->SetParameter(2,lim_Sigma[0]);
	eff_fit->SetParLimits(2,lim_Sigma[1],lim_Sigma[2]);

	eff_fit->SetParName(3,"mu");
	eff_fit->SetParameter(3,E_tol[0]);
	eff_fit->SetParLimits(3,E_tol[1],E_tol[2]);

	eff_fit->SetParName(4,"H");
	eff_fit->SetParameter(4,lim_Height[0]);
	eff_fit->SetParLimits(4,lim_Height[1],lim_Height[2]);

	TFitResultPtr r = effGraph->Fit(eff_fit,"S Q");
//	effGraph->Fit(eff_fit);

	double eff = eff_fit->Eval(eGamma);

	double x[1] = {eGamma};
	double err[1];
	r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false);
/*
	TLine *x1 = new TLine(eGamma,-0.01,eGamma,eff);
	x1->SetLineColor(4);
	x1->SetLineWidth(2);
	x1->Draw();
*/
	cout << "Efficiency @ " << eGamma << " keV = " << eff << " +/- " << err[0] << endl;

	int ngr = 24;
	TH1D *ci = new TH1D("ci","ci",1500,0,1500);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(ci,0.683);
	ci->SetLineColor(-1);
	ci->SetFillColorAlpha(2,0.25);
	ci->SetFillStyle(3002);
	ci->Draw("E3");
	ci->GetXaxis()->CenterTitle();
	ci->GetYaxis()->CenterTitle();
	ci->GetXaxis()->SetTitle("Efficiency [%]");
	ci->GetYaxis()->SetTitle("Energy [keV]");
	ci->GetXaxis()->SetTitleSize(0.05);
	ci->GetYaxis()->SetTitleSize(0.05);
	ci->GetXaxis()->SetTitleOffset(0.95);
	ci->GetYaxis()->SetTitleOffset(0.8);
	effGraph->Draw("p");
//	x1->Draw();
	eff_fit->SetLineColor(4);
	eff_fit->Draw("same");
	effGraph->Draw("psame");
}

/********************************************************
*			Fit functions			*
********************************************************/
// Fit function for multiple peaks
double CBeff(double *x, double *par)
{
	double xcur  = x[0];
	double alpha = par[0]; 
	double n     = par[1]; 
	double sigma = par[2]; 
	double mu    = par[3];
	double H     = par[4];

	return H*ROOT::Math::crystalball_function(xcur, alpha, n, sigma, mu);
}
