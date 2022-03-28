#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "Math/PdfFuncMathCore.h"
#include "TVirtualFitter.h"

#include "fits.h"

void gamFit(const char* histo, double rMin, double rMax, int nOfPeaks, ...)
{
	TH1D *h = (TH1D*)df->Get(histo);
	double maxArea = h->GetEntries();

	double lim_Area[2]  = {10, maxArea};		// Area
	double E_tol        = 1.5;			// +/- E est
	double lim_sigma[3] = {0.5, 0.15, 1.5};		// sigma
//	double lim_sigma[3] = {0.5, 0.15, 100.5};		// sigma

	rangeMin = rMin;
	rangeMax = rMax;

	int x1 = h->GetBin(rangeMin);
	int x2 = h->GetBin(rangeMax);
	double y1 = h->GetBinContent(x1);
	double y2 = h->GetBinContent(x2);

	binFactor = h->GetBinWidth(0);

	double grad = (y2-y1)/(x2-x1);

	TF1 *f1 = new TF1("f1",fGam,rangeMin,rangeMax,3+2*nOfPeaks);
	TF1 *f2 = new TF1("f2",bg,rangeMin,rangeMax,2);
	f1->SetNpx(1000);

	f1->SetParameter(1,grad);
	f1->SetParameter(2,lim_sigma[0]);
	f1->SetParLimits(2,lim_sigma[1],lim_sigma[2]);

	f1->SetParName(0,"BG offset");
	f1->SetParName(1,"BG grad");
	f1->SetParName(2,"sigma");

	// Vectors for storing centroid & height estimates
	vector<double> ePeak, aPeak;

	// Variable list for reading in your N variables
	va_list vl;
	va_start(vl,nOfPeaks);

	// Start reading in your variables, find the min and max energy value
	// and set your histo+fit limits to min-100,max+100
	double val;
	cout << "\n";

	for(int i=0; i<nOfPeaks; i++)
	{
		val = va_arg(vl, double);
		ePeak.push_back(val);
	}

	nPeaks = ePeak.size();

	for(int i=0; i<nOfPeaks;i++)
	{
		ostringstream eName,aName;
		eName << i+1 << ". Energy";
		aName << i+1 << ". Area";

		f1->SetParName(2*i+3,eName.str().c_str());
		f1->SetParameter(2*i+3,ePeak[i]);
		f1->SetParLimits(2*i+3,ePeak[i]-E_tol,ePeak[i]+E_tol);

		f1->SetParName(2*i+4,aName.str().c_str());
		f1->SetParameter(2*i+4,100);
		f1->SetParLimits(2*i+4,lim_Area[0],lim_Area[1]);
	}

	h->Fit(f1,"R");
	f2->FixParameter(0,f1->GetParameter(0));
	f2->FixParameter(1,f1->GetParameter(1));
//	f2->SetLineColor(kAzure+2);
	f2->SetLineColor(4);
	
	TF1 *ip[nOfPeaks];
	for(int i=0; i<nOfPeaks; i++)
	{
		ostringstream ipT;
		ipT << "ip[" << i << "]";
		ip[i] = new TF1(ipT.str().c_str(),fGam,rangeMin,rangeMax,5);
		ip[i]->FixParameter(0,f1->GetParameter(0));
		ip[i]->FixParameter(1,f1->GetParameter(1));
		ip[i]->FixParameter(2,f1->GetParameter(2));
		ip[i]->FixParameter(3,f1->GetParameter(2*i+3));
		ip[i]->FixParameter(4,f1->GetParameter(2*i+4));
		ip[i]->SetLineColor(4);
		ip[i]->SetLineStyle(2);
		ip[i]->SetLineWidth(2);
		ip[i]->SetNpx(1000);
		h->Fit(ip[i],"Q R+");
	}

	h->Fit(f2,"Q R+");
}

void fitPeaks(const char* histo, double rMin, double rMax, int nOfPeaks, ...)
{
	double lim_Area[2]  = {0.1, 180000};		// height
	double E_tol        = 20.;			// +/- E est
	double lim_n[3]     = {0.01, 0.01, 25.2};	// n
	double lim_Alpha[3] = {0.01, 0.01, 1.8};	// alpha
//	double lim_Sigma[3] = {12, 2.5, 1500};		// sigma
	double lim_Sigma[3] = {12, 7, 15};		// sigma

	eMin = rMin;
	eMax = rMax;

	// Vectors for storing centroid & height estimates
	vector<double> ePeak, hPeak;

	// Variable list for reading in your N variables
	va_list vl;
	va_start(vl,nOfPeaks);

	// Start reading in your variables, find the min and max energy value
	// and set your histo+fit limits to min-100,max+100
	double val;
	cout << "\n";
	for(Int_t i=0; i<(2*nOfPeaks); i++)
	{
		val = (double)va_arg(vl,int);
		if(i%2==0)
		{
			ePeak.push_back(val);
		}
		else
		{
			hPeak.push_back(val);
		}
	}

	// Work out how many peaks and how many parameters there are to fit
	nPeaks = ePeak.size();
	int nParams = 3+(2*nPeaks);

	// Open file
	TFile *file1;
//	checkFile();// TFile *file1 = (TFile*)gROOT->GetFile();
	if(checkFile()) file1 = (TFile*)gROOT->GetFile();
	else
	{
		errorFile();
		return;
	}

	// Draw that bad boy
	TCanvas *cbFit = new TCanvas("cbFit","cbFit");
	TPad *sp[4];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[2]",0.0,0.4,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[1]",0.0,0.0,1.0,0.4);
	cbFit->cd();
	sp[0]->Draw();
	sp[0]->cd();
	sp[1]->Draw();
	sp[1]->cd();

	TH1F *alpha;
	if( checkHisto(df,histo) ) alpha = (TH1F*)df->Get(histo);
	else	
	{
		errorHisto(histo);
		return;
	}

	//TH1F* alpha = (TH1F*)df->Get(histo);
	alpha->Draw();
	logy();
	alpha->SetTitle("Si1+Si2 alpha spectrum");
	alpha->GetXaxis()->SetRangeUser(rMin,rMax);

	// Get binFactor
	binFactor = ((alpha->GetXaxis()->GetXmax())-(alpha->GetXaxis()->GetXmin()))/(alpha->GetNbinsX());
	cout << "\n\nTHE BIN FACTOR IS: " << binFactor << "\n";

	// Declare the fit
	TF1 *f1 = new TF1("f1",CB,eMin,eMax,nParams);

	// Peak alpha value (where tail joins gaussian)
	f1->SetParName(0,"alpha");
	f1->SetParameter(0,lim_Alpha[0]);
	f1->SetParLimits(0,lim_Alpha[1],lim_Alpha[2]);

	// Peak n value (how quickly the tail drops off)
	f1->SetParName(1,"n");
	f1->SetParameter(1,lim_n[0]);
	f1->SetParLimits(1,lim_n[1],lim_n[2]);

	// Peak sigma value (width of Gaussian component)
	f1->SetParName(2,"sigma");
	f1->SetParameter(2,lim_Sigma[0]);
	f1->SetParLimits(2,lim_Sigma[1],lim_Sigma[2]);

	// Set parameters and limits for fit
	for(int i=0; i<nPeaks; i++)
	{
		ostringstream ossH;
		ossH << i+1 << ". Height";
		string sArea = ossH.str();
		f1->SetParName(3+2*i,sArea.c_str());	
		f1->SetParameter(3+2*i,hPeak[i]);
		f1->SetParLimits(3+2*i,lim_Area[0],lim_Area[1]);

		ostringstream ossE;
		ossE << i+1 << ". Energy";
		string sEnergy = ossE.str();
		// Peak energy (centroid of peak)
		f1->SetParName(4+2*i,sEnergy.c_str());
		f1->SetParameter(4+2*i,ePeak[i]);
		f1->SetParLimits(4+2*i,ePeak[i]-E_tol,ePeak[i]+E_tol);
	}
	
	// Fit it with binned-likelihood method
	f1->SetNpx(500);
//	TFitResultPtr r = alpha->Fit(f1,"LL S R");
	TFitResultPtr r = alpha->Fit(f1,"W S R");

	// Store fit results as histogram with confidence error intervals in order to calculate
	// residual between fit and histo later on
	TH1D *res = new TH1D("res","res",alpha->GetNbinsX(),alpha->GetXaxis()->GetXmin(),alpha->GetXaxis()->GetXmax());
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(res,conf_int);
	TMatrixDSym covTot = r->GetCovarianceMatrix();

	// Declare array of functions for individual peaks and vectors for their parameters
	TF1 *fP[nPeaks];
	vector<double> pH, pE, pIntegral, eIntegral;

	// Read in alpha, n and sigma values of peaks	
	double pAlpha = f1->GetParameter(0);
	double pn     = f1->GetParameter(1);
	double pSig   = f1->GetParameter(2);

	// Draw individual peaks and integrate are underneath them for intensities
	for(int i=0; i<nPeaks; i++)
	{
		// Declare fits
		ostringstream title;
		title << "fP[" << i << "]";
		fP[i] = new TF1(title.str().c_str(),CB1,eMin,eMax,5);
		fP[i]->SetNpx(500);

		// Store parameters of individual peaks from mother fit
		pH.push_back(f1->GetParameter(3+2*i));
		pE.push_back(f1->GetParameter(4+2*i));

		// Get covariance matric for indiividual fit components
		// based on cov. matrix of mother fit
		TMatrixD covMatrix(5,5);
		for(int j=0; j<5; j++)
		{
			for(Int_t k=0; k<5; k++)
			{
				int s = j;
				int t = k;
				if(j>2) s = j+(2*i);
				if(k>2) t = k+(2*i);
				
				covMatrix[j][k] = covTot[s][t];
			}
		}

		// Set parameters of individual peaks
/*		fP[i]->FixParameter(0,pAlpha);
		fP[i]->FixParameter(1,pn);
		fP[i]->FixParameter(2,pSig);
		fP[i]->FixParameter(3,pH[i]);
		fP[i]->FixParameter(4,pE[i]);
*/
		fP[i]->SetParameter(0,pAlpha);
		fP[i]->SetParameter(1,pn);
		fP[i]->SetParameter(2,pSig);
		fP[i]->SetParameter(3,pH[i]);
		fP[i]->SetParameter(4,pE[i]);

		// Draw and integrate peaks
		fP[i]->SetLineColor(4);
		//alpha->Fit(fP[i],"Q R+");
		pIntegral.push_back((1./binFactor)*fP[i]->Integral(eMin,eMax));
		eIntegral.push_back((1./binFactor)*fP[i]->IntegralError(eMin,eMax,fP[i]->GetParameters(),covMatrix.GetMatrixArray()));
//		alpha->Fit(fP[i],"Q R+");
		fP[i]->Draw("same");
	}

	// Draw mother fit on top of individual peaks	
	f1->Draw("same");

	// Calculate and output reduced chi^2
	double redCHI2 = (f1->GetChisquare())/(f1->GetNDF());
	cout << "\n\nChi2/ndf = " << redCHI2 << endl;

	// Output integral+error of individual peaks
	for(int i=0; i<nPeaks; i++)
	{
		cout << "Integral of peak @ " << pE[i] << " keV = " << pIntegral[i] << " +/- " << eIntegral[i] << endl;
	}

	// Draw residual
	sp[0]->cd();
	sp[2]->Draw();
	sp[2]->cd();
	res->Add(alpha,-1);
	res->Scale(-1);
	res->Draw();
	res->GetXaxis()->SetRangeUser(rMin,rMax);

	sp[1]->cd();

	// clean up
	hPeak.clear();
	ePeak.clear();
	pH.clear();
	pE.clear();
	eIntegral.clear();
	eIntegral.clear();
}

void alphasBG(const char* histo, double rMin, double rMax, int nOfPeaks, ...)
{
	double lim_Area[2]  = {1, 180000};		// height
//	double E_tol        = 10.;			// +/- E est
	double E_tol        = 10.;			// +/- E est
	double lim_n[3]     = {1.4, 0.7, 2.83};	// n
//	double lim_n[3]     = {2.83, 0.9, 3.0};	// n
	double lim_Alpha[3] = {1.65, 1.02, 1.8};	// alpha
//	double lim_Alpha[3] = {0.52, 0.42, 1.8};	// alpha
	double lim_Sigma[3] = {20, 2, 30};		// sigma

	eMin = rMin;
	eMax = rMax;

	// Vectors for storing centroid & height estimates
	vector<double> ePeak, hPeak;

	// Variable list for reading in your N variables
	va_list vl;
	va_start(vl,nOfPeaks);

	// Start reading in your variables, find the min and max energy value
	// and set your histo+fit limits to min-100,max+100
	double val;
	cout << "\n";
	for(int i=0; i<(2*nOfPeaks); i++)
	{
		val = (double)va_arg(vl,int);
		if(i%2==0)
		{
			ePeak.push_back(val);
		}
		else
		{
			hPeak.push_back(val);
		}
	}

	// Work out how many peaks and how many parameters there are to fit
	nPeaks = ePeak.size();
	int nParams = 3+(2*nPeaks);

	// Open file
	TFile *file1 = (TFile*)gROOT->GetFile();

	// Draw that bad boy
	TCanvas *cbFit = new TCanvas("cbFit","cbFit");
	TPad *sp[4];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[2]",0.0,0.4,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[1]",0.0,0.0,1.0,0.4);
	cbFit->cd();
	sp[0]->Draw();
	sp[0]->cd();
	sp[1]->Draw();
	sp[1]->cd();
	TH1F* alpha = (TH1F*)df->Get(histo);
	alpha->Draw();
	logy();
	alpha->SetTitle("Si1+Si2 alpha spectrum");

	alpha->GetXaxis()->SetRangeUser(rMin,rMax);

	// Get binFactor
	binFactor = ((alpha->GetXaxis()->GetXmax())-(alpha->GetXaxis()->GetXmin()))/(alpha->GetNbinsX());
	cout << "\n\nTHE BIN FACTOR IS: " << binFactor << "\n";

	// Declare the fit
	TF1 *f1 = new TF1("f1",CBBG,eMin,eMax,nParams+2);

	// Peak alpha value (where tail joins gaussian)
	f1->SetParName(0,"alpha");
	f1->SetParameter(0,lim_Alpha[0]);
	f1->SetParLimits(0,lim_Alpha[1],lim_Alpha[2]);

	// Peak n value (how quickly the tail drops off)
	f1->SetParName(1,"n");
	f1->SetParameter(1,lim_n[0]);
	f1->SetParLimits(1,lim_n[1],lim_n[2]);

	// Peak sigma value (width of Gaussian component)
	f1->SetParName(2,"sigma");
	f1->SetParameter(2,lim_Sigma[0]);
	f1->SetParLimits(2,lim_Sigma[1],lim_Sigma[2]);

	// Set parameters and limits for fit
	for(int i=0; i<nPeaks; i++)
	{
		ostringstream ossH;
		ossH << i+1 << ". Height";
		string sArea = ossH.str();
		f1->SetParName(3+2*i,sArea.c_str());	
		f1->SetParameter(3+2*i,hPeak[i]);
		f1->SetParLimits(3+2*i,lim_Area[0],lim_Area[1]);

		ostringstream ossE;
		ossE << i+1 << ". Energy";
		string sEnergy = ossE.str();
		// Peak energy (centroid of peak)
		f1->SetParName(4+2*i,sEnergy.c_str());
		f1->SetParameter(4+2*i,ePeak[i]);
		f1->SetParLimits(4+2*i,ePeak[i]-E_tol,ePeak[i]+E_tol);
	}
	
	// Fit it with binned-likelihood method
	f1->SetNpx(500);
	TFitResultPtr r = alpha->Fit(f1,"LL S R");
//	TFitResultPtr r = alpha->Fit(f1,"S R");
//	TFitResultPtr r = alpha->Fit(f1,"WW S R");

	// Store fit results as histogram with confidence error intervals in order to calculate
	// residual between fit and histo later on
	TH1D *res = new TH1D("res","res",alpha->GetNbinsX(),alpha->GetXaxis()->GetXmin(),alpha->GetXaxis()->GetXmax());
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(res,conf_int);
	TMatrixDSym covTot = r->GetCovarianceMatrix();

	// Declare array of functions for individual peaks and vectors for their parameters
	TF1 *fP[nPeaks];
	vector<double> pH, pE, pIntegral, eIntegral;

	// Read in alpha, n and sigma values of peaks	
	double pAlpha = f1->GetParameter(0);
	double pn     = f1->GetParameter(1);
	double pSig   = f1->GetParameter(2);

	// Draw individual peaks and integrate are underneath them for intensities
	for(int i=0; i<nPeaks; i++)
	{
		// Declare fits
		ostringstream title;
		title << "fP[" << i << "]";
		fP[i] = new TF1(title.str().c_str(),CB1,eMin,eMax,5);
		fP[i]->SetNpx(500);

		// Store parameters of individual peaks from mother fit
		pH.push_back(f1->GetParameter(3+2*i));
		pE.push_back(f1->GetParameter(4+2*i));

		// Get covariance matric for indiividual fit components
		// based on cov. matrix of mother fit
		TMatrixD covMatrix(5,5);
		for(int j=0; j<5; j++)
		{
			for(int k=0; k<5; k++)
			{
				int s = j;
				int t = k;
				if(j>2) s = j+(2*i);
				if(k>2) t = k+(2*i);
				
				covMatrix[j][k] = covTot[s][t];
			}
		}

		// Set parameters of individual peaks
		fP[i]->FixParameter(0,pAlpha);
		fP[i]->FixParameter(1,pn);
		fP[i]->FixParameter(2,pSig);
		fP[i]->FixParameter(3,pH[i]);
		fP[i]->FixParameter(4,pE[i]);

		// Draw and integrate peaks
		fP[i]->SetLineColor(4);
		alpha->Fit(fP[i],"Q R+");
		pIntegral.push_back((1./binFactor)*fP[i]->Integral(eMin,eMax));
		eIntegral.push_back((1./binFactor)*fP[i]->IntegralError(eMin,eMax,fP[i]->GetParameters(),covMatrix.GetMatrixArray()));
	}

	// Draw mother fit on top of individual peaks	
	f1->Draw("same");

	// Calculate and output reduced chi^2
	double redCHI2 = (f1->GetChisquare())/(f1->GetNDF());
	cout << "\n\nChi2/ndf = " << redCHI2 << endl;

	// Output integral+error of individual peaks
	for(int i=0; i<nPeaks; i++)
	{
		cout << "Integral of peak @ " << pE[i] << " keV = " << pIntegral[i] << " +/- " << eIntegral[i] << endl;
	}

	// Draw residual
	sp[0]->cd();
	sp[2]->Draw();
	sp[2]->cd();
	res->Add(alpha,-1);
	res->Scale(-1);
	res->Draw();
	res->GetXaxis()->SetRangeUser(rMin,rMax);

	sp[1]->cd();

	// clean up
	hPeak.clear();
	ePeak.clear();
	pH.clear();
	pE.clear();
	eIntegral.clear();
	eIntegral.clear();
}

void graphCB(const char* graph, double rMin, double rMax, int nOfPeaks, ...)
{
	double lim_Area[2]  = {10, 180000};		// height
	double E_tol        = 10.;			// +/- E est
	double lim_n[3]     = {1.4, 0.9, 2.8};	// n
	double lim_Alpha[3] = {1.65, 0.5, 2.8};	// alpha
	double lim_Sigma[3] = {12, 3.9, 1500};		// sigma

	eMin = rMin;
	eMax = rMax;

	// Vectors for storing centroid & height estimates
	vector<double> ePeak, hPeak;

	// Variable list for reading in your N variables
	va_list vl;
	va_start(vl,nOfPeaks);

	// Start reading in your variables, find the min and max energy value
	// and set your histo+fit limits to min-100,max+100
	double val;
	cout << "\n";
	for(Int_t i=0; i<(2*nOfPeaks); i++)
	{
		val = (double)va_arg(vl,int);
		if(i%2==0)
		{
			ePeak.push_back(val);
		}
		else
		{
			hPeak.push_back(val);
		}
	}

	// Work out how many peaks and how many parameters there are to fit
	nPeaks = ePeak.size();
	int nParams = 3+(2*nPeaks);

	// Open file
	TFile *file1 = (TFile*)gROOT->GetFile();

	// Draw that bad boy
	TCanvas *cbFit = new TCanvas("cbFit","cbFit");
	TPad *sp[4];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[2]",0.0,0.4,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[1]",0.0,0.0,1.0,0.4);
	cbFit->cd();
	sp[0]->Draw();
	sp[0]->cd();
	sp[1]->Draw();
	sp[1]->cd();
	TGraph *alpha = (TGraph*)file1->Get(graph);
	alpha->Draw("A*");
	alpha->SetMarkerStyle(20);
	alpha->SetMarkerSize(0.5);
	logy();
	alpha->SetTitle("Si1+Si2 alpha spectrum");
	alpha->GetXaxis()->SetRangeUser(rMin,rMax);

	// Get binFactor
//	binFactor = ((alpha->GetXaxis()->GetXmax())-(alpha->GetXaxis()->GetXmin()))/(alpha->GetNbinsX());
//	cout << "\n\nTHE BIN FACTOR IS: " << binFactor << "\n";

	// Declare the fit
	TF1 *f1 = new TF1("f1",CB,eMin,eMax,nParams);

	// Peak alpha value (where tail joins gaussian)
	f1->SetParName(0,"alpha");
	f1->SetParameter(0,lim_Alpha[0]);
	f1->SetParLimits(0,lim_Alpha[1],lim_Alpha[2]);

	// Peak n value (how quickly the tail drops off)
	f1->SetParName(1,"n");
	f1->SetParameter(1,lim_n[0]);
	f1->SetParLimits(1,lim_n[1],lim_n[2]);

	// Peak sigma value (width of Gaussian component)
	f1->SetParName(2,"sigma");
	f1->SetParameter(2,lim_Sigma[0]);
	f1->SetParLimits(2,lim_Sigma[1],lim_Sigma[2]);

	// Set parameters and limits for fit
	for(int i=0; i<nPeaks; i++)
	{
		ostringstream ossH;
		ossH << i+1 << ". Height";
		string sArea = ossH.str();
		f1->SetParName(3+2*i,sArea.c_str());	
		f1->SetParameter(3+2*i,hPeak[i]);
		f1->SetParLimits(3+2*i,lim_Area[0],lim_Area[1]);

		ostringstream ossE;
		ossE << i+1 << ". Energy";
		string sEnergy = ossE.str();
		// Peak energy (centroid of peak)
		f1->SetParName(4+2*i,sEnergy.c_str());
		f1->SetParameter(4+2*i,ePeak[i]);
		f1->SetParLimits(4+2*i,ePeak[i]-E_tol,ePeak[i]+E_tol);
	}
	
	// Fit it with binned-likelihood method
	f1->SetNpx(500);
	TFitResultPtr r = alpha->Fit(f1,"S R");

	// Store fit results as histogram with confidence error intervals in order to calculate
	// residual between fit and histo later on
/*	TH1D *res = new TH1D("res","res",alpha->GetNbinsX(),alpha->GetXaxis()->GetXmin(),alpha->GetXaxis()->GetXmax());
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(res,conf_int);
	TMatrixDSym covTot = r->GetCovarianceMatrix();
*/
	// Declare array of functions for individual peaks and vectors for their parameters
	TF1 *fP[nPeaks];
	vector<double> pH, pE, pIntegral, eIntegral;

	// Read in alpha, n and sigma values of peaks	
	double pAlpha = f1->GetParameter(0);
	double pn     = f1->GetParameter(1);
	double pSig   = f1->GetParameter(2);

	// Draw individual peaks and integrate are underneath them for intensities
	for(int i=0; i<nPeaks; i++)
	{
		// Declare fits
		ostringstream title;
		title << "fP[" << i << "]";
		fP[i] = new TF1(title.str().c_str(),CB1,eMin,eMax,5);
		fP[i]->SetNpx(500);

		// Store parameters of individual peaks from mother fit
		pH.push_back(f1->GetParameter(3+2*i));
		pE.push_back(f1->GetParameter(4+2*i));

		// Get covariance matric for indiividual fit components
		// based on cov. matrix of mother fit
/*		TMatrixD covMatrix(5,5);
		for(int j=0; j<5; j++)
		{
			for(Int_t k=0; k<5; k++)
			{
				int s = j;
				int t = k;
				if(j>2) s = j+(2*i);
				if(k>2) t = k+(2*i);
				
				covMatrix[j][k] = covTot[s][t];
			}
		}
*/
		// Set parameters of individual peaks
		fP[i]->FixParameter(0,pAlpha);
		fP[i]->FixParameter(1,pn);
		fP[i]->FixParameter(2,pSig);
		fP[i]->FixParameter(3,pH[i]);
		fP[i]->FixParameter(4,pE[i]);

		fP[i]->SetParameter(0,pAlpha);
		fP[i]->SetParameter(1,pn);
		fP[i]->SetParameter(2,pSig);
		fP[i]->SetParameter(3,pH[i]);
		fP[i]->SetParameter(4,pE[i]);

		// Draw and integrate peaks
		fP[i]->SetLineColor(4);
		//alpha->Fit(fP[i],"Q R+");
//		pIntegral.push_back((1./binFactor)*fP[i]->Integral(eMin,eMax));
//		eIntegral.push_back((1./binFactor)*fP[i]->IntegralError(eMin,eMax,fP[i]->GetParameters(),covMatrix.GetMatrixArray()));
//		alpha->Fit(fP[i],"Q R+");
		fP[i]->Draw("same");
	}

	// Draw mother fit on top of individual peaks	
	f1->Draw("same");

	// Calculate and output reduced chi^2
	double redCHI2 = (f1->GetChisquare())/(f1->GetNDF());
	cout << "\n\nChi2/ndf = " << redCHI2 << endl;

/*	// Output integral+error of individual peaks
	for(int i=0; i<nPeaks; i++)
	{
		cout << "Integral of peak @ " << pE[i] << " keV = " << pIntegral[i] << " +/- " << eIntegral[i] << endl;
	}
*/
	// Draw residual
	sp[0]->cd();
	sp[2]->Draw();
	sp[2]->cd();
/*	res->Add(alpha,-1);
	res->Scale(-1);
	res->Draw();
	res->GetXaxis()->SetRangeUser(rMin,rMax);
*/
	sp[1]->cd();

	// clean up
	hPeak.clear();
	ePeak.clear();
	pH.clear();
	pE.clear();
	eIntegral.clear();
	eIntegral.clear();

}

void graphAlphasBG(const char* graph, double rMin, double rMax, int nOfPeaks, ...)
{
	double lim_Area[2]  = {1, 180000};		// height
	double E_tol        = 10.;			// +/- E est
	double lim_n[3]     = {1.4, 1.0, 2.83};	// n
	double lim_Alpha[3] = {1.65, 0.12, 1.8};	// alpha
	double lim_Sigma[3] = {12, 1.3, 15};		// sigma

	eMin = rMin;
	eMax = rMax;

	// Vectors for storing centroid & height estimates
	vector<double> ePeak, hPeak;

	// Variable list for reading in your N variables
	va_list vl;
	va_start(vl,nOfPeaks);

	// Start reading in your variables, find the min and max energy value
	// and set your histo+fit limits to min-100,max+100
	double val;
	cout << "\n";
	for(int i=0; i<(2*nOfPeaks); i++)
	{
		val = (double)va_arg(vl,int);
		if(i%2==0)
		{
			ePeak.push_back(val);
		}
		else
		{
			hPeak.push_back(val);
		}
	}

	// Work out how many peaks and how many parameters there are to fit
	nPeaks = ePeak.size();
	int nParams = 3+(2*nPeaks);

	// Open file
	TFile *file1 = (TFile*)gROOT->GetFile();

	// Draw that bad boy
	TCanvas *cbFit = new TCanvas("cbFit","cbFit");
	TPad *sp[4];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[2]",0.0,0.4,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[1]",0.0,0.0,1.0,0.4);
	cbFit->cd();
	sp[0]->Draw();
	sp[0]->cd();
	sp[1]->Draw();
	sp[1]->cd();
	TGraph* alpha = (TGraph*)df->Get(graph);
	alpha->Draw("A*");
	alpha->SetMarkerStyle(20);
	alpha->SetMarkerSize(0.5);
	logy();
	alpha->SetTitle("Si1+Si2 alpha spectrum");

	alpha->GetXaxis()->SetRangeUser(rMin,rMax);

	// Get binFactor
//	binFactor = ((alpha->GetXaxis()->GetXmax())-(alpha->GetXaxis()->GetXmin()))/(alpha->GetNbinsX());
//	cout << "\n\nTHE BIN FACTOR IS: " << binFactor << "\n";

	// Declare the fit
	TF1 *f1 = new TF1("f1",CBBG,eMin,eMax,nParams+2);

	// Peak alpha value (where tail joins gaussian)
	f1->SetParName(0,"alpha");
	f1->SetParameter(0,lim_Alpha[0]);
	f1->SetParLimits(0,lim_Alpha[1],lim_Alpha[2]);

	// Peak n value (how quickly the tail drops off)
	f1->SetParName(1,"n");
	f1->SetParameter(1,lim_n[0]);
	f1->SetParLimits(1,lim_n[1],lim_n[2]);

	// Peak sigma value (width of Gaussian component)
	f1->SetParName(2,"sigma");
	f1->SetParameter(2,lim_Sigma[0]);
	f1->SetParLimits(2,lim_Sigma[1],lim_Sigma[2]);

	// Set parameters and limits for fit
	for(int i=0; i<nPeaks; i++)
	{
		ostringstream ossH;
		ossH << i+1 << ". Height";
		string sArea = ossH.str();
		f1->SetParName(3+2*i,sArea.c_str());	
		f1->SetParameter(3+2*i,hPeak[i]);
		f1->SetParLimits(3+2*i,lim_Area[0],lim_Area[1]);

		ostringstream ossE;
		ossE << i+1 << ". Energy";
		string sEnergy = ossE.str();
		// Peak energy (centroid of peak)
		f1->SetParName(4+2*i,sEnergy.c_str());
		f1->SetParameter(4+2*i,ePeak[i]);
		f1->SetParLimits(4+2*i,ePeak[i]-E_tol,ePeak[i]+E_tol);
	}
	
	// Fit it with binned-likelihood method
	f1->SetNpx(500);
	TFitResultPtr r = alpha->Fit(f1,"S R");

	// Store fit results as histogram with confidence error intervals in order to calculate
	// residual between fit and histo later on
//	TH1D *res = new TH1D("res","res",alpha->GetNbinsX(),alpha->GetXaxis()->GetXmin(),alpha->GetXaxis()->GetXmax());
//	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(res,conf_int);
	TGraphErrors *gr = new TGraphErrors(alpha->GetN());
	for(int i=0; i<alpha->GetN(); i++)
	{
		double res = alpha->GetX()[i] - f1->Eval(alpha->GetX()[i]);
		gr->SetPoint(i,alpha->GetX()[i],0);
	}

	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr,conf_int);

	TGraphErrors *grint = new TGraphErrors(alpha->GetN());
	for(int i=0; i<alpha->GetN(); i++)
	{
		grint->SetPoint(i,alpha->GetX()[i],alpha->GetY()[i] - gr->GetY()[i]);
		grint->SetPointError(i,0,gr->GetErrorY(i));
	}

	// Get covariance matrix
	TMatrixDSym covTot = r->GetCovarianceMatrix();

	// Declare array of functions for individual peaks and vectors for their parameters
	TF1 *fP[nPeaks];
	vector<double> pH, pE, pIntegral, eIntegral;

	// Read in alpha, n and sigma values of peaks	
	double pAlpha = f1->GetParameter(0);
	double pn     = f1->GetParameter(1);
	double pSig   = f1->GetParameter(2);

	// Draw individual peaks and integrate are underneath them for intensities
	for(int i=0; i<nPeaks; i++)
	{
		// Declare fits
		ostringstream title;
		title << "fP[" << i << "]";
		fP[i] = new TF1(title.str().c_str(),CB1,eMin,eMax,5);
		fP[i]->SetNpx(500);

		// Store parameters of individual peaks from mother fit
		pH.push_back(f1->GetParameter(3+2*i));
		pE.push_back(f1->GetParameter(4+2*i));

		// Get covariance matric for indiividual fit components
		// based on cov. matrix of mother fit
		TMatrixD covMatrix(5,5);
		for(int j=0; j<5; j++)
		{
			for(int k=0; k<5; k++)
			{
				int s = j;
				int t = k;
				if(j>2) s = j+(2*i);
				if(k>2) t = k+(2*i);
				
				covMatrix[j][k] = covTot[s][t];
			}
		}

		// Set parameters of individual peaks
		fP[i]->FixParameter(0,pAlpha);
		fP[i]->FixParameter(1,pn);
		fP[i]->FixParameter(2,pSig);
		fP[i]->FixParameter(3,pH[i]);
		fP[i]->FixParameter(4,pE[i]);

		// Draw and integrate peaks
		fP[i]->SetLineColor(4);
		alpha->Fit(fP[i],"Q R+");
		pIntegral.push_back(fP[i]->Integral(eMin,eMax));
		eIntegral.push_back(fP[i]->IntegralError(eMin,eMax,fP[i]->GetParameters(),covMatrix.GetMatrixArray()));
	}

	// Draw mother fit on top of individual peaks	
	f1->Draw("same");

	// Calculate and output reduced chi^2
	double redCHI2 = (f1->GetChisquare())/(f1->GetNDF());
	cout << "\n\nChi2/ndf = " << redCHI2 << endl;

	// Output integral+error of individual peaks
	for(int i=0; i<nPeaks; i++)
	{
		cout << "Integral of peak @ " << pE[i] << " keV = " << pIntegral[i] << " +/- " << eIntegral[i] << endl;
	}

	// Draw residual
	sp[0]->cd();
	sp[2]->Draw();
	sp[2]->cd();
	grint->Draw("AP");
	grint->SetMarkerStyle(20);
	grint->SetMarkerSize(0.5);
	grint->GetXaxis()->SetRangeUser(rMin,rMax);
	sp[1]->cd();

	// clean up
	hPeak.clear();
	ePeak.clear();
	pH.clear();
	pE.clear();
	eIntegral.clear();
	eIntegral.clear();
}

void HLFit(const char *histo, double rMin, double rMax, double HL)
{
	TH1D *h = (TH1D*)df->Get(histo);
	rangeMin = rMin;
	rangeMax = rMax;

	TF1 *f1 = new TF1("f1","expo",rangeMin,rangeMax);
	f1->SetNpx(1000);
	// 
	f1->SetParameter(1,TMath::Log(2)/HL);

	h->Fit(f1,"R");

	cout << "\n\n  Halflife is: " << -1.*TMath::Log(2)/f1->GetParameter(1) << " +/- " << -1.*TMath::Log(2)*f1->GetParError(1)/f1->GetParameter(1) << " s\n\n";

//	h->Fit(f2,"Q R+");
}

/********************************************************
*                 BG + gammas function                  *
********************************************************/
// Fit for multiple peaks
double fGam(double *x, double *par)
{
	double xcur  = x[0];
 
	double y0  = par[0]; 
	double grd = par[1]; 
	double sig = par[2];
	
	double sum = y0 + grd*xcur;

	for (int i=0; i<nPeaks; i++)
	{
		double C = par[3+2*i];
		double A = binFactor*par[4+2*i];
		sum += A*TMath::Gaus(xcur, C, sig, kTRUE);
	}
	return sum;
}

double bg(double *x, double *par)
{
	return par[0] + par[1]*x[0];
}

/********************************************************
*                   CB functions                        *
********************************************************/
// Fit function for multiple peaks
double CB(double *x, double *par)
{
	double xcur  = x[0];
 
	double alpha = par[0]; 
	double n     = par[1]; 
	double sigma = par[2];
	
	double sum = 0;
	for (int i=0; i<nPeaks; i++)
	{
		double H     = par[3+2*i];
		double mu    = par[4+2*i];

		sum += H*ROOT::Math::crystalball_function(xcur, alpha, n, sigma, mu);
	}
	return sum;
}

// Fit function for multiple peaks + linear background
double CBBG(double *x, double *par)
{
	double xcur  = x[0];
 
	double alpha = par[0]; 
	double n     = par[1]; 
	double sigma = par[2];
	
	double sum = 0;
	for (int i=0; i<nPeaks; i++)
	{
		double H     = par[3+2*i];
		double mu    = par[4+2*i];
		sum += H*ROOT::Math::crystalball_function(xcur, alpha, n, sigma, mu);
	}
	sum+= (par[5+2*(nPeaks-1)]*xcur) + par[6+2*(nPeaks-1)];
	return sum;
}

// Fit function for individual peaks
double CB1(double *x, double *par)
{
	double xcur  = x[0];
 
	double alpha = par[0]; 
	double n     = par[1]; 
	double sigma = par[2]; 
	double H     = par[3];
	double mu    = par[4];

	return H*ROOT::Math::crystalball_function(xcur, alpha, n, sigma, mu);
}

