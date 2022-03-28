// Based on https://doi.org/10.1007/s100500070129
#include "TRandom.h"

// Simulate sum of nOfDecs decay curves. For example, two
// curves with T1/2 10 and 2000 units, and N= 300 and 400
// decays, use command:
// KHS_sim(2,10,300,2000,400)
void KHS_sim(int nOfDecs, ...)
{
	TRandom *r = new TRandom();

	TH1D *dt = new TH1D("dt","dt",4010,-100,40000);
	TH1D *kt = new TH1D("kt","kt",2000,-100,100);

	TH1D *h[nOfDecs+1];
	TH1D *k[nOfDecs+1];

	va_list vl;
	va_start(vl,nOfDecs);

	double val;
	vector<double> par;
	vector<double> col = {2, 4, 6, 5};

	for(int i=0; i<(2*nOfDecs); i++)
	{
		val = (double)va_arg(vl,int);
		par.push_back(val);
	}

	for(int i=0; i<nOfDecs; i++)
	{
		cout << "HL: " << par[i] << "\n";
		cout << "NC: " << par[i+1] << "\n";

		ostringstream title1, title2;
		title1 << "h[" << i+1 << "]";
		title2 << "k[" << i+1 << "]";
		h[i+1] = new TH1D(title1.str().c_str(),title1.str().c_str(),40100,-100,40000);
		k[i+1] = new TH1D(title2.str().c_str(),title2.str().c_str(),2000,-100,100);


		for(int j=0; j<par[2*i+1]; j++)
		{
			double a = r->Exp(par[2*i]/(TMath::Log(2)));
			h[i+1]->Fill(a);
			dt->Fill(a);
			k[i+1]->Fill(TMath::Log10(a));
			kt->Fill(TMath::Log10(a));
		}
	}

	TCanvas *c = new TCanvas();
	c->Divide(1,2);
	c->cd(1);
	dt->SetLineColor(1);
	dt->Draw("hist E1");

	c->cd(2);
	kt->SetLineColor(1);
	kt->Draw("hist E1");

	for(int i=0; i<nOfDecs; i++)
	{
		if(i>3) col.push_back(i);
		c->cd(1);
		h[i+1]->SetLineColor(col[i]);
		h[i+1]->Draw("hist E1 same");

		c->cd(2);
		k[i+1]->SetLineColor(col[i]);
		k[i+1]->Draw("hist E1 same");
	}
}
