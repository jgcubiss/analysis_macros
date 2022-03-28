void old_sim()
{
	TRandom *r = new TRandom();

	TH1D *h1 = new TH1D("h1","h1",401,-100,40000);
	TH1D *h2 = new TH1D("h1","h1",401,-100,40000);
	TH1D *h3 = new TH1D("h1","h1",401,-100,40100);
	TH1D *h4 = new TH1D("h2","h2",2000,-100,100);
	TH1D *h5 = new TH1D("h2","h2",2000,-100,100);
	TH1D *h6 = new TH1D("h2","h2",2000,-100,100);

	for(int i=0; i<3300; i++)
	{
		double a = r->Exp(30000.0/(TMath::Log(2)));
		double b = r->Exp(1700.0/(TMath::Log(2)));

		h1->Fill(a);
		h2->Fill(b);
		h3->Fill(a);
		h3->Fill(b);
		h4->Fill(TMath::Log(a));
		h5->Fill(TMath::Log(b));
		h6->Fill(TMath::Log(a));
		h6->Fill(TMath::Log(b));
	}

	TCanvas *c = new TCanvas();
	c->Divide(1,2);
	c->cd(1);
	h1->SetLineColor(2);
	h2->SetLineColor(4);
	h3->SetLineColor(1);
	h3->Draw();
	h1->Draw("same");
	h2->Draw("same");

	c->cd(2);
	h4->SetLineColor(2);
	h5->SetLineColor(4);
	h6->SetLineColor(1);
	h6->Draw();
	h4->Draw("same");
	h5->Draw("same");


}
