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

#include "checks.h"

bool checkFile()
{
	if(gROOT->GetFile()) return true;
	else
	{
		errorFile();
		return false;
	}
}

void errorFile()
{
	cout << "\n\n\tERROR:\n";
	cout << "\tNo file open\n\n\n";
	return;
}

bool checkHisto(TFile *tmpFile, const char *nameHisto)
{
	if(tmpFile->Get(nameHisto)) return true;
	else 
	{
		errorHisto(nameHisto);
		return false;
	}
}

void errorHisto(const char *nameHisto)
{
	cout << "\n\n\tERROR:\n";
	cout << "\tHisto: \"" << nameHisto << "\" does not exist\n\n\n";
	return;
}

