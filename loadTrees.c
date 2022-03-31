#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "loadTrees.h"

#include "spectra.c"
#include "fits.c"
#include "hfsSim.c"
#include "KHS_sim.c"
#include "effLege.c"
//#include "analyse_hfs.c"

bool loadTrees()
{
	if(!gROOT->GetFile())
	{
		cout << "\n\t No root file open\n\n";
		return false;
	}

	df = gROOT->GetFile();

	if( gROOT->GetFile()->Get("single_data") != NULL)  sd = (TTree*)gROOT->GetFile()->Get("single_data");
	if( gROOT->GetFile()->Get("sg")          != NULL)  sg = (TTree*)gROOT->GetFile()->Get("sg");
	if( gROOT->GetFile()->Get("gg")          != NULL)  gg = (TTree*)gROOT->GetFile()->Get("gg");
	if( gROOT->GetFile()->Get("hfs")         != NULL)  hf = (TTree*)gROOT->GetFile()->Get("hfs");
	if( gROOT->GetFile()->Get("Broad")       != NULL)  g4 = (TTree*)gROOT->GetFile()->Get("Broad");
	if( gROOT->GetFile()->Get("ids")         != NULL) ids = (TTree*)gROOT->GetFile()->Get("ids");
	return true;
}

void unloadTrees()
{
	gROOT->GetFile()->cd();
	sd->Delete();
	sg->Delete();
	gg->Delete();
	hf->Delete();
	g4->Delete();
}
