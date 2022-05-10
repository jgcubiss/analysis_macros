#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TF1.h"

#include "checks.c"

/*****************************************************
*              list available spectra                *
*****************************************************/
void spectra()
{
	cout << "\n\n"
	     << " *************************************************************************\n" 
	     << " *   Commands for drawing spectra:                                       *\n"
	     << " *                                                                       *\n"
	     << " *   Singles:                                                            *\n"
	     << " *   single alpha    |       a([Si no.],[bin width])                     *\n"
	     << " *   single gamma    |       g([Ge no.],[bin width])                     *\n"
	     << " *   all Si spectra  |   fullA([Si no.],[bin width])                     *\n"
	     << " *   Si1+Si2         |    impA([bin width])                              *\n"
	     << " *   Si3+Si4         |    decA([bin width])                              *\n"
	     << " *   compare Si1+Si2 | compS12([bin width])                              *\n"
	     << " *   compare Si3+Si4 | compS34([bin width])                              *\n"
	     << " *   compare Si2+Si3 | compS23([bin width])                              *\n"
	     << " *                                                                       *\n"
	     << " *   Coincidences:                                                       *\n"
	     << " *   Si1:Si2         |    SS12([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *   Si3:Si4         |    SS34([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *   Ge1:Si1+2       |      GS([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *   comp. G1:S1/S2  |     GSC([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *   Ge2:Si1+2       |     G2S([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *   Ge1:Ge2         |      GG([timeWindowLow],[timeWindowHigh])         *\n"
	     << " *                                                                       *\n"
	     << " *   Projection and zooms:                                               *\n"
	     << " *   Y-Projection    |      yP([histoName],[eMin],[eMax])                *\n"
	     << " *   X-Projection    |      xP([histoName],[eMin],[eMax])                *\n"
	     << " *   Y-axis zoom     |      zY([histoName],[eMin],[eMax])                *\n"
	     << " *   Z-axis zoom     |      zX([histoName],[eMin],[eMax])                *\n"
	     << " *   XY-axis zoom    |     zXY([histoName],[xMin],[xMax],[yMin],[yMax])  *\n"
	     << " *************************************************************************\n\n";
}

/*****************************************************
*                 Singles spectra                    *
*****************************************************/
// Single Si spectrum
void a(int i, double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}
	loadTrees();
	if(gROOT->FindObject("a")) delete gROOT->FindObject("a");
	ostringstream drw, gate;
	drw << "sEng[" << i << "]>>a";
	gate << "sEng[" << i << "]>0";
	TH1D *a = new TH1D("a","a",(1/binWidth)*10000,0,10000);
	sd->Draw(drw.str().c_str(),gate.str().c_str());
	a->GetXaxis()->CenterTitle();
	a->GetXaxis()->SetTitle("E_{#alpha} [keV]");
	a->GetYaxis()->CenterTitle();
	ostringstream yTit;
	yTit << "Counts / " << binWidth << " keV";
	a->GetYaxis()->SetTitle(yTit.str().c_str());
}

// Single Si spectrum
void circ(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}
	loadTrees();
	if(gROOT->FindObject("circ")) delete gROOT->FindObject("circ");
	ostringstream drw, gate;
	drw << "E_circ>>circ";
	gate << "E_circ>0";
	TH1D *circ = new TH1D("circ","circ",(1/binWidth)*10000,0,10000);
	g4->Draw(drw.str().c_str(),gate.str().c_str());
	circ->GetXaxis()->CenterTitle();
	circ->GetXaxis()->SetTitle("E_{#alpha} [keV]");
	circ->GetYaxis()->CenterTitle();
	ostringstream yTit;
	yTit << "Counts / " << binWidth << " keV";
	circ->GetYaxis()->SetTitle(yTit.str().c_str());
}

// Single Ge spectrum
void g(int i, double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	if(gROOT->FindObject("g")) delete gROOT->FindObject("g");
	ostringstream drw, gate;
	drw << "gEng[" << i << "]>>g";
	gate << "gEng[" << i << "]>0";
	TH1D *g = new TH1D("g","g",(1/binWidth)*5000,0,5000);
	sd->Draw(drw.str().c_str(),gate.str().c_str());
}

// All 4 Si spectra
void fullA()
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH1D *aH[5];
	TCanvas *c = new TCanvas("c","c",700,900);
	c->Divide(1,4);

	for(int i=0; i<4; i++)
	{
		ostringstream hT, drw, gate;
		hT << "aH[" << i+1 << "]";
		drw << "sEng[" << i+1 << "]>>aH[" << i+1 << "]";
		gate << "sEng[" << i+1 << "]>0";

		aH[i] = new TH1D(hT.str().c_str(),hT.str().c_str(),10000,0,10000);

		c->cd(i+1);
		sd->Draw(drw.str().c_str(),gate.str().c_str());
	}
}

// Sum of Si1+2
void impA(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH1D *a1 = new TH1D("a1","a1",(1/binWidth)*10000,0,10000);
	TH1D *a2 = new TH1D("a2","a2",(1/binWidth)*10000,0,10000);
	TH1D *impA = new TH1D("impA","impA",(1/binWidth)*10000,0,10000);

	sd->Project("a1","sEng[1]","sEng[1]>0");
	sd->Project("a2","sEng[2]","sEng[2]>0");
	impA->Add(a1,1);
	impA->Add(a2,1);
	impA->Draw();
	impA->GetXaxis()->CenterTitle();
	impA->GetXaxis()->SetTitle("E_{#alpha} [keV]");
	impA->GetYaxis()->CenterTitle();
	ostringstream yTit;
	yTit << "Counts / " << binWidth << " keV";
	impA->GetYaxis()->SetTitle(yTit.str().c_str());

}

void impShade(double binWidth, double aMin, double aMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH1D *a1 = new TH1D("a1","a1",(1/binWidth)*10000,0,10000);
	TH1D *s1 = new TH1D("s1","s1",(1/binWidth)*10000,0,10000);
	TH1D *a2 = new TH1D("a2","a2",(1/binWidth)*10000,0,10000);
	TH1D *s2 = new TH1D("s2","s2",(1/binWidth)*10000,0,10000);
	TH1D *impA = new TH1D("impA","impA",(1/binWidth)*10000,0,10000);
	TH1D *impS = new TH1D("impS","impS",(1/binWidth)*10000,0,10000);

	ostringstream g1, g2;
	g1 << "sEng[1]>" << aMin << "&&sEng[1]<" << aMax;
	g2 << "sEng[2]>" << aMin << "&&sEng[2]<" << aMax;

	sd->Project("a1","sEng[1]","sEng[1]>0");
	sd->Project("s1","sEng[1]",g1.str().c_str());
	sd->Project("a2","sEng[2]","sEng[2]>0");
	sd->Project("s2","sEng[2]",g2.str().c_str());
	impA->Add(a1,1);
	impA->Add(a2,1);
	impS->Add(s1,1);
	impS->Add(s2,1);
	impS->SetFillColorAlpha(2,0.35);
	impA->Draw();
	impS->Draw("hist same");
	impA->GetXaxis()->CenterTitle();
	impA->GetXaxis()->SetTitle("E_{#alpha} [keV]");
	impA->GetYaxis()->CenterTitle();
	ostringstream yTit;
	yTit << "Counts / " << binWidth << " keV";
	impA->GetYaxis()->SetTitle(yTit.str().c_str());

}

// Sum of Si3+4
void decA(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH1D *a3 = new TH1D("a3","a3",(1/binWidth)*10000,0,10000);
	TH1D *a4 = new TH1D("a4","a4",(1/binWidth)*10000,0,10000);
	TH1D *decA = new TH1D("decA","decA",(1/binWidth)*10000,0,10000);

	sd->Project("a3","sEng[3]","sEng[3]>0");
	sd->Project("a4","sEng[4]","sEng[4]>0");
	decA->Add(a3,1);
	decA->Add(a4,1);
	decA->Draw();
}

// Compare silicons
void compS12(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	int nBins = (10000./binWidth);
	TH1D *s1 = new TH1D("s1","s1",nBins,0,10000);
	TH1D *s2 = new TH1D("s2","s2",nBins,0,10000);

	sd->Draw("sEng[2]>>s2","sEng[2]>0","hist");
	sd->Draw("sEng[1]>>s1","sEng[1]>0","hist same");
	s1->SetLineColor(2);
}

void ratioS12(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	int nBins = (10000./binWidth);
	TH1D *s1 = new TH1D("s1","s1",nBins,0,10000);
	TH1D *s2 = new TH1D("s2","s2",nBins,0,10000);
	TH1D *s3 = new TH1D("s3","s3",nBins,0,10000);

	sd->Draw("sEng[2]>>s2","sEng[2]>0","hist");
	sd->Draw("sEng[1]>>s1","sEng[1]>0","hist same");
	sd->Draw("sEng[1]>>s3","sEng[1]>0","hist same");
	s3->Divide(s2);
	s3->Draw("same");
	s2->SetLineColor(2);
	s3->SetLineColor(4);
	s3->SetMarkerColor(4);
}

void compS34(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	int nBins = (10000./binWidth);
	TH1D *s3 = new TH1D("s3","s3",nBins,0,10000);
	TH1D *s4 = new TH1D("s4","s4",nBins,0,10000);

	sd->Draw("sEng[3]>>s3","sEng[3]>0");
	sd->Draw("sEng[4]>>s4","sEng[4]>0","same");
	s4->SetLineColor(2);
}
void compS23(double binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	int nBins = (10000./binWidth);
	TH1D *s3 = new TH1D("s3","s3",nBins,0,10000);
	TH1D *s2 = new TH1D("s2","s2",nBins,0,10000);

	sd->Draw("sEng[2]>>s2","sEng[2]>0");
	sd->Draw("sEng[3]>>s3","sEng[3]>0","same");
	s3->SetLineColor(2);
}

void hfs(double aMin, double aMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TCanvas *HFS = new TCanvas();
//	TH1D *lf = new TH1D("lf","lf");
	ostringstream gate;
	gate << "((SEng[1]>" << aMin << "&&SEng[1]<" << aMax << ")||(SEng[2]>" << aMin << "&&SEng[2]<" << aMax << "))&&(Slf>10000)";
	hf->Draw("Slf>>lf",gate.str().c_str());
}


void hfs2D(double aBin, double aMin, double aMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TCanvas *HFS = new TCanvas();
	ostringstream gate;
	gate << "(SEng>" << aMin << ")&&((SEng[1]>" << aMin << "&&SEng[1]<" << aMax << ")||(SEng[2]>" << aMin << "&&SEng[2]<" << aMax << "))&&(Slf>10000)";
//	hf->Draw("SEng:Slf>>lftmp",gate.str().c_str());
//	TH1D *lftmp = (TH1D*)gROOT->GetFile()->Get("lftmp");

	double lfMin   = hf->GetMinimum("Slf");
	double lfMax   = hf->GetMaximum("Slf");
	double lfRange = lfMax - lfMin;
	int  nFreqBins = (int)50*(lfRange+0.2);

//	TH2D *lf = new TH2D("lf","lf",60,lfMin-0.1*lfRange,lfMax+0.1*lfRange,(aMax-aMin+200)/aBin,aMin,aMax);
	TH2D *lf = new TH2D("lf","lf",(aMax-aMin+200)/aBin,aMin,aMax,nFreqBins,lfMin-0.1*lfRange,lfMax+0.1*lfRange);
//	hf->Draw("SEng:Slf>>lf",gate.str().c_str(),"surf1");
	hf->Draw("Slf:SEng>>lf",gate.str().c_str(),"colz");
	lf->GetXaxis()->SetTitle("Wavenumber [cm^{-1}]");
	lf->GetXaxis()->SetTitleOffset(1.34);
	lf->GetXaxis()->CenterTitle();
	lf->GetYaxis()->SetTitle("E_{#alpha} [keV]");
	lf->GetYaxis()->SetTitleOffset(1.54);
	lf->GetYaxis()->CenterTitle();
}



void simTot(int binWidth)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TCanvas *simRes = new TCanvas();
	ostringstream gate;
	gate << "E_ann>1 || E_circ>1";
	TH1D *sim1   = new TH1D("sim1","sim1",(1./binWidth)*10000,0,10000);
	TH1D *sim2   = new TH1D("sim2","sim2",(1./binWidth)*10000,0,10000);
	TH1D *simTot = new TH1D("simTot","simTot",(1./binWidth)*10000,0,10000);

	g4->Project("sim1","E_ann",gate.str().c_str());
	g4->Project("sim2","E_circ",gate.str().c_str());

	simTot->Add(sim1,1);
	simTot->Add(sim2,1);
	simTot->Draw();
}

void deadTime()
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	double timeMax = 1.2*(sd->GetMaximum("pulserTime") / 1e6);

	TH1D *pTime = new TH1D("pTime","pTime",timeMax,0,timeMax);
//	sd->Draw("pulserTime/1e6>>pTime","pulserEvent>0");
	sd->Project("pTime","pulserTime/1e6","pulserEvent>0");

	TF1 *aveLT = new TF1("aveLT","pol0",1,timeMax-1);
	TFitResultPtr r = pTime->Fit(aveLT,"S Q N");
/*
	TCanvas *dt = new TCanvas();
	TH1D *resLT = new TH1D("resLT","resLT",timeMax,0,timeMax);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(resLT,0.683);
	resLT->SetLineColor(-1);
	resLT->SetFillColorAlpha(2,0.25);
	resLT->Draw("E3 same");
	aveLT->Draw("same");
*/
	cout << "\n\n   Average deadtime: " << 100-aveLT->GetParameter(0) << "%\n\n";

}

/*****************************************************
*               Coincidence spectra                  *
*****************************************************/

// Si1-Si2 coincs
void SSI(double tMin, double tMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *SSI = new TH2D("SSI","SSI",2000,0,10000,2000,0,10000);
	sg->Draw("senA[1]:senA[2]>>SSI","","colz");
	SSI->SetMinimum();
}

// Si3-Si4 coincs
void SSD(double tMin, double tMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *SSD = new TH2D("SSD","SSD",2000,0,10000,2000,0,10000);
	sg->Draw("senA[3]:senA[4]>>SSD","","colz");
	SSD->SetMinimum();
}

// Ge1:Si1+2 coincs
void GS(double tMin, double tMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	ostringstream tGate;
	tGate << "(gtd>" << tMin << "&&gtd<" << tMax << ")";

	double gamBinWidth=0.25, alfBinWidth=1.0;
	double gMin=0, gMax=2500;
	double aMin=0, aMax=10000;

	TH2D *gs1 = new TH2D("gs1","gs1",(1/alfBinWidth)*(aMax-aMin),aMin,aMax,(1/gamBinWidth)*(gMax-gMin),gMin,gMax);
	TH2D *gs2 = new TH2D("gs2","gs2",(1/alfBinWidth)*(aMax-aMin),aMin,aMax,(1/gamBinWidth)*(gMax-gMin),gMin,gMax);
	TH2D *GS = new TH2D("GS","GS",   (1/alfBinWidth)*(aMax-aMin),aMin,aMax,(1/gamBinWidth)*(gMax-gMin),gMin,gMax);

	sg->Project("gs1","gen[1]:senA[1]",tGate.str().c_str());
	sg->Project("gs2","gen[1]:senA[2]",tGate.str().c_str());

	GS->Add(gs1,1);
	GS->Add(gs2,1);
	//GS->Draw();
	GS->Draw("colz");
	GS->SetMinimum(1);
//	gPad->SetLogz();
}
// Ge2:Si1+2 coincs
void G2S(double tMin, double tMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	ostringstream tGate;
	tGate << "(gtd>" << tMin << "&&gtd<" << tMax << ")";

	TH2D *gs1 = new TH2D("gs1","gs1",5000,0,10000,10000,0,2500);
	TH2D *gs2 = new TH2D("gs2","gs2",5000,0,10000,10000,0,2500);
	TH2D *G2S = new TH2D("G2S","G2S",5000,0,10000,10000,0,2500);

	sg->Project("gs1","gen[2]:senA[1]",tGate.str().c_str());
	sg->Project("gs2","gen[2]:senA[2]",tGate.str().c_str());

	G2S->Add(gs1,1);
	G2S->Add(gs2,1);
	G2S->Draw("colz");
	G2S->SetMinimum(1);
	gPad->SetLogz();
}

// Ge1:Ge2 coincs
void GG(double wBin, double tMin, double tMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	ostringstream tGate;
	tGate << "(gtd>" << tMin << "&&gtd<" << tMax << ")";

	double eMin = 0;
	double eMax = 2500;
	double nBin = (eMax-eMin)/wBin;

	TCanvas *ggCanv = new TCanvas();
//	TH2D *GG = new TH2D("GG","GG",10000,0,2500,10000,0,2500);
//	TH2D *GG = new TH2D("GG","GG",2500,0,2500,2500,0,2500);
	TH2D *GG = new TH2D("GG","GG",nBin,eMin,eMax,nBin,eMin,eMax);
	gg->Draw("gen[1]:gen[2]>>GG",tGate.str().c_str(),"colz");
	GG->SetMinimum(1);
	ostringstream title;
	title << "#gamma-#gamma for " << tMin << "<#Deltat<" << tMax << " ns";
	GG->SetTitle(title.str().c_str());
	gPad->SetLogz();
}

/*****************************************************
*                Coinc time distrs                   *
*****************************************************/
void agtd(double aMin, double aMax, double gMin, double gMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	ostringstream gate, title;
	gate << "((senA[1]> " << aMin << "&&senA[1]<" << aMax << ")||(senA[2]> " << aMin << "&&senA[2]<" << aMax << ")) && (gen[1]>" << gMin << "&&gen[1]<" << gMax << ")";

	title << aMin << "#leqE_{#alpha}#leq" << aMax << "      " << gMin << "#leqE_{#gamma}#leq" << gMax;

	TH1D *agtd = new TH1D("agtd",title.str().c_str(),20000,-20000,20000);

	TCanvas *agTD = new TCanvas();
	sg->Draw("gtd>>agtd",gate.str().c_str(),"hist");
}

void ggtd(double gMin1, double gMax1, double gMin2, double gMax2)
{
	loadTrees();
	ostringstream gate;
	gate << "(gen[1]> " << gMin1 << "&&gen[1]<" << gMax1 << ") && (gen[2]>" << gMin2<< "&&gen[2]<" << gMax2 << ")";

	TH1D *ggtd = new TH1D("ggtd","ggtd",20000,-20000,20000);

	TCanvas *ggTD = new TCanvas();
	gg->Draw("gtd>>ggtd",gate.str().c_str(),"hist");
}
/*****************************************************
*                Matrix projections                  *
*****************************************************/
void mSlice(const char* histo)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *xS = h->ProjectionX("xS");
	TH1D *yS = h->ProjectionY("yS");

	TCanvas *sl = new TCanvas();
	TPad *sp[3];
	sp[0] = new TPad("sp[0]","sp[0]",0,0,1,1);
	sp[1] = new TPad("sp[1]","sp[1]",0,0.5,1,1);
	sp[2] = new TPad("sp[2]","sp[2]",0,0,1,0.5);
	sp[1]->SetBottomMargin(0);
	sp[2]->SetTopMargin(0);
	sl->cd();
	sp[0]->Draw();
	sp[0]->cd();
	sp[1]->Draw();
	sp[1]->cd();
	xS->Draw("hist");
	sp[0]->cd();
	sp[2]->Draw();
	sp[2]->cd();
	yS->Draw("hist");
}
// y-axis projection
void yP(const char* histo, double eMin, double eMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *yP = h->ProjectionY("yP",h->GetXaxis()->FindBin(eMin)+1,h->GetXaxis()->FindBin(eMax));
	TCanvas *yProjCanv = new TCanvas("yProjCanv","yProjCanv",1030,610);
	yP->Draw("H E0");
}

void compYP(const char* histo, double eMin1, double eMax1, double eMin2, double eMax2)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *h = (TH2D*)df->Get(histo);
	ostringstream name1, name2;
	name1 << "yP_" << eMin1 << "_" << eMax1;
	name2 << "yP_" << eMin2 << "_" << eMax2;

	TH1D *yP1 = h->ProjectionY(name1.str().c_str(),h->GetXaxis()->FindBin(eMin1)+1,h->GetXaxis()->FindBin(eMax1));
	TH1D *yP2 = h->ProjectionY(name2.str().c_str(),h->GetXaxis()->FindBin(eMin2)+1,h->GetXaxis()->FindBin(eMax2));

	TCanvas *yProjCanv = new TCanvas("yProjCanv","yProjCanv",1030,610);
	yP1->Draw("H E0");
	yP2->Draw("H E0 same");
	yP2->SetLineColor(2);

	if(yP2->GetMaximum() > yP1->GetMaximum()) yP1->SetMaximum(1.2*yP2->GetMaximum());
}

// x-axis projection
void xP(const char* histo, double eMin, double eMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *xP = h->ProjectionX("xP",h->GetYaxis()->FindBin(eMin)+1,h->GetYaxis()->FindBin(eMax));
	TCanvas *xProjCanv = new TCanvas("xProjCanv","xProjCanv",1030,610);
	xP->Draw("H E0");
}

// gamma-gamma projections
void ggP(const char* histo, double eMin, double eMax, double bgMin, double bgMax)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	loadTrees();
	TCanvas *ggP = new TCanvas("ggP","ggP",1300,600);
	TPad *sp[5];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[1]",0.0,0.5,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[2]",0.0,0.0,1.0,0.5);

	sp[1]->SetBottomMargin(0);
	sp[2]->SetTopMargin(0);

	ggP->cd();
	sp[0]->Draw();
	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *cp[3], *bg[3];

	for(int i=0; i<2; i++)
	{
		sp[0]->cd();
		sp[i+1]->Draw();
		sp[i+1]->cd();

		ostringstream pT, bgT;
		pT  << "xP_" << i+1;
		bgT << "bg_" << i+1;

		if(i==0)
		{
			cp[i] = h->ProjectionX(pT.str().c_str(), h->GetYaxis()->FindBin(eMin)+1, h->GetYaxis()->FindBin(eMax));
			bg[i] = h->ProjectionX(bgT.str().c_str(),h->GetYaxis()->FindBin(bgMin)+1,h->GetYaxis()->FindBin(bgMax));
		}
		else
		{
			cp[i] = h->ProjectionY(pT.str().c_str(), h->GetYaxis()->FindBin(eMin)+1, h->GetYaxis()->FindBin(eMax));

			bg[i] = h->ProjectionY(bgT.str().c_str(),h->GetYaxis()->FindBin(bgMin)+1,h->GetYaxis()->FindBin(bgMax));
		}


		cp[i]->Add(bg[i],-1*(eMax-eMin)/(bgMax-bgMin));
		cp[i]->Draw("hist");
	}
}



// gamma-gamma or
void ggPor(const char* histo, int nGates, ...)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	vector<double> gateVal;
	double val;

	va_list vl;
	va_start(vl,nGates);

	// Read in values
	for(int i=0; i<4*nGates; i++)
	{
		val = va_arg(vl, double);
		gateVal.push_back(val);
	}

	loadTrees();
	TCanvas *ggP = new TCanvas("ggP","ggP",1300,600);
	TPad *sp[5];
	sp[0] = new TPad("sp[0]","sp[0]",0.0,0.0,1.0,1.0);
	sp[1] = new TPad("sp[1]","sp[1]",0.0,0.5,1.0,1.0);
	sp[2] = new TPad("sp[2]","sp[2]",0.0,0.0,1.0,0.5);

	sp[1]->SetBottomMargin(0);
	sp[2]->SetTopMargin(0);

	ggP->cd();
	sp[0]->Draw();
	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *cp[3], *bg[3];
	
	TH1D *totHist[2];
	totHist[0] = new TH1D("tH0","tH0",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
	totHist[1] = new TH1D("tH1","tH1",h->GetNbinsY(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
	cout << "\n\n";
	for(int i=0; i<nGates; i++)
	{
		for(int j=0; j<2; j++)
		{
			ostringstream pT, bgT;
			pT  << "xP_" << j+1;
			bgT << "bg_" << j+1;

			double eMin  = gateVal[4*i+0];
			double eMax  = gateVal[4*i+1];
			double bgMin = gateVal[4*i+2];
			double bgMax = gateVal[4*i+3];

			if(j==0)
			{
				cout << "   Gate on: " << eMin << " to " << eMax << "\t BG on: " << bgMin << " to " << bgMax << "\n";
				cp[j] = h->ProjectionX(pT.str().c_str(), h->GetYaxis()->FindBin(eMin)+1, h->GetYaxis()->FindBin(eMax));
				bg[j] = h->ProjectionX(bgT.str().c_str(),h->GetYaxis()->FindBin(bgMin)+1,h->GetYaxis()->FindBin(bgMax));
			}
			else
			{
				cp[j] = h->ProjectionY(pT.str().c_str(), h->GetYaxis()->FindBin(eMin)+1, h->GetYaxis()->FindBin(eMax));

				bg[j] = h->ProjectionY(bgT.str().c_str(),h->GetYaxis()->FindBin(bgMin)+1,h->GetYaxis()->FindBin(bgMax));
			}


			cp[j]->Add(bg[j],-1*(eMax-eMin)/(bgMax-bgMin));
			totHist[j]->Add(cp[j],1);
		}
	}
	cout << "\n\n";
	for(int i=0; i<2; i++)
	{
		sp[0]->cd();
		sp[i+1]->Draw();
		sp[i+1]->cd();
		totHist[i]->Draw("hist");
	}
	gateVal.clear();
}




// gamma-gamma or
void pgg(const char* histo, int nGates, ...)
{
	if( !checkFile() )
	{
		errorFile();
		return;
	}

	vector<double> gateVal;
	double val;

	va_list vl;
	va_start(vl,nGates);

	// Read in values
	for(int i=0; i<4*nGates; i++)
	{
		val = va_arg(vl, double);
		gateVal.push_back(val);
	}

	loadTrees();
	TCanvas *pgg = new TCanvas("pgg","pgg",800,600);
	pgg->cd();

	TH2D *h = (TH2D*)df->Get(histo);
	TH1D *cp, *bg;
	
	TH1D *totHist;
	totHist = new TH1D("tH0","tH0",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
	cout << "\n\n";

	for(int i=0; i<nGates; i++)
	{
		ostringstream pT, bgT;
		pT  << "xP";
		bgT << "bg";

		double eMin  = gateVal[4*i+0];
		double eMax  = gateVal[4*i+1];
		double bgMin = gateVal[4*i+2];
		double bgMax = gateVal[4*i+3];

		cout << "   Gate on: " << eMin << " to " << eMax << "\t BG on: " << bgMin << " to " << bgMax << "\n";
		cp = h->ProjectionX(pT.str().c_str(), h->GetYaxis()->FindBin(eMin)+1, h->GetYaxis()->FindBin(eMax));
		bg = h->ProjectionX(bgT.str().c_str(),h->GetYaxis()->FindBin(bgMin)+1,h->GetYaxis()->FindBin(bgMax));

		cp->Add(bg,-1*(eMax-eMin)/(bgMax-bgMin));
		totHist->Add(cp,1);
	}

	totHist->Draw("hist");
	gateVal.clear();
}






/*****************************************************
*                  Zoom functions                    *
*****************************************************/

void zX(const char* histo, double eMin, double eMax)
{
	if(df->Get(histo)->InheritsFrom(TH1D::Class()))
	{
		TH1D *h = (TH1D*)df->Get(histo);
		h->GetXaxis()->SetRangeUser(eMin,eMax);
		h->Draw();
	}
		
	if(df->Get(histo)->InheritsFrom(TH2D::Class()))
	{
		TH2D *h = (TH2D*)df->Get(histo);
		h->GetXaxis()->SetRangeUser(eMin,eMax);
		h->Draw("colz");
	}
}

void zY(const char* histo, double eMin, double eMax)
{
	if(df->Get(histo)->InheritsFrom(TH1D::Class()))
	{
		TH1D *h = (TH1D*)df->Get(histo);
		h->GetYaxis()->SetRangeUser(eMin,eMax);
		h->Draw();
	}
		
	if(df->Get(histo)->InheritsFrom(TH2D::Class()))
	{
		TH2D *h = (TH2D*)df->Get(histo);
		h->GetYaxis()->SetRangeUser(eMin,eMax);
		h->Draw("colz");
	}
}

void zZ(const char* histo, double eMin, double eMax)
{
	if(df->Get(histo)->InheritsFrom(TH1D::Class()))
	{
		TH1D *h = (TH1D*)df->Get(histo);
		h->GetZaxis()->SetRangeUser(eMin,eMax);
		h->Draw();
	}
		
	if(df->Get(histo)->InheritsFrom(TH2D::Class()))
	{
		TH2D *h = (TH2D*)df->Get(histo);
		h->GetZaxis()->SetRangeUser(eMin,eMax);
		h->Draw("colz");
	}
}


void zXY(const char* histo, double xMin, double xMax, double yMin, double yMax)
{
	TFile *df = (TFile*)gROOT->GetFile();
	if(gROOT->GetFile()->Get(histo)->InheritsFrom(TH1D::Class()))
	{
		TH1D *h = (TH1D*)df->Get(histo);
		h->GetXaxis()->SetRangeUser(xMin,xMax);
		h->GetYaxis()->SetRangeUser(yMin,yMax);
		h->Draw();
	}
		
	if(gROOT->GetFile()->Get(histo)->InheritsFrom(TH2D::Class()))
	{
		TH2D *h = (TH2D*)df->Get(histo);
		h->GetXaxis()->SetRangeUser(xMin,xMax);
		h->GetYaxis()->SetRangeUser(yMin,yMax);
		h->Draw("colz");
	}
}
void logx()
{
	if(!gPad->GetLogx()) gPad->SetLogx();
	else gPad->SetLogx(0);
}
void logy()
{
	if(!gPad->GetLogy()) gPad->SetLogy();
	else gPad->SetLogy(0);
}
void logz()
{
	if(!gPad->GetLogz()) gPad->SetLogz();
	else gPad->SetLogz(0);
}
