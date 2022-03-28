#include "hfsSim.h"

#include "Math/SpecFuncMathMore.h"
#include "Math/SpecFunc.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>

//float nSpin  = 1.5;		// Spin of nucleus
//float eSpinL = 0.5;		// Spin lower e- level
//float eSpinU = 0.5;		// Spin upper e- level

float nSpin;
float eSpinL;		// Spin lower e- level
float eSpinU;		// Spin upper e- level

double CoG = 0;	// Estimate centre of gravity

//double A_U = 700;    // Estimate A factor
//double A_L = 300;    // Estimate A factor

double A_U = -52000;    // Estimate A factor
double A_L = -10000;    // Estimate A factor

double B_U = 694;   // Estimate B factor
double B_L = -226;   // Estimate B factor

double wGaus = 800;
double wLore = 20;

/************************************************************************
*			Main part of code				*
************************************************************************/
void hfsSim(float spinEL, float spinEU, float spinN)
{
	nSpin  = spinN;
	eSpinL = spinEL;
	eSpinU = spinEU;
	
	/**********************************************************************
	*                Calculate positions & intensities                    *
	**********************************************************************/
	// How many transitions (peaks) are there?
	nPeaks = peaks(nSpin,eSpinL,eSpinU);

	vector<string> trans;
	
	// Calculate relative intensities of peaks and their positions relative to the estimated CoG
	for (int i=0; i<F_states(eSpinL); i++)
	{
		FL.push_back(fabs(nSpin - eSpinL)+i);
		for (int j=0; j<F_states(eSpinU); j++)
		{
			FU.push_back(fabs(nSpin - eSpinU)+j);
			// Check transition is valid
			if( ( FL.at(i)==0 && FU.at(j)==0 ) || ( abs(FL.at(i)-FU.at(j))>1 )  )
			{}
			else
			{
				double guess = freq_shift(FL.at(i),FU.at(j),A_L,B_L,A_U,B_U)+CoG;

				ostringstream oss;
				string temp;
				double intRel = intensity(FL.at(i),FU.at(j));

				rI.push_back(intRel);
				f_l.push_back(FL.at(i));
				f_u.push_back(FU.at(j));
					
				if(FL.at(i)==FU.at(j)||(FL.at(i)==FU.at(j)+1)||(FL.at(i)==FU.at(j)-1))
				{
					gC.push_back(guess);
					oss << FL.at(i) << "->" << FU.at(j) << " line";
					temp = oss.str();
					trans.push_back(temp);
				}

			} // Close if statement on whether transition is valid
		} // Close for loop over possible F_upper states
	} // Close for loop over possible  F_lower states

	// Find the most intense transition and the maximum point of the HFS
	double maxI = *max_element(rI.begin(),rI.end());
	
	cout << "\n\nNumber of peaks: " << nPeaks << "\n\n";
	
	cout << "Transition\tIntensity\tPosition" << endl;
	for(int i=0; i<(int)rI.size(); i++)
			/*if(valTrans.at(i) == 1)*/cout << setprecision(0) << "   " << f_l.at(i) << "->" << f_u.at(i) << setprecision(2) << fixed << "\t\t   " << rI.at(i)/maxI << "\t\t " << gC.at(i) << endl;
	cout << "\n\n";

	// Calculate width of structure
	double minR = *min_element(gC.begin(),gC.end());
	double maxR = *max_element(gC.begin(),gC.end());
	double mrgR = 0.35*(maxR-minR);
	
	/**********************************************************************
	*                       Fit in frequency space                        *
	**********************************************************************/

	// Calculate the number of peaks and declare a voigt function with nPeaks Voigts
	TF1 *sf = new TF1("",voigt,minR-mrgR,maxR+mrgR,7+nPeaks);

	// Set parameters for function
	sf->SetParameter(0,CoG);
	sf->SetParameter(1,A_U);
	sf->SetParameter(2,A_L);
	sf->SetParameter(3,B_U);
	sf->SetParameter(4,B_L);
	sf->SetParameter(5,wGaus);
	sf->SetParameter(6,wLore);

	// Set parameter names for output
	sf->SetParName(0,"CoG");
	sf->SetParName(1,"A_U");
	sf->SetParName(2,"A_L");
	sf->SetParName(3,"B_U");
	sf->SetParName(4,"B_L");
	sf->SetParName(5,"wGaus");
	sf->SetParName(6,"wLore");

	// Set heigth+limit for each peak, as well as output which F_i -> F_f transitions this is
	vector<double> pAreaFreq, pHeight;
	for(int i=0; i<nPeaks; i++)
	{
		// Set name of parameter to "F_i->F_f" of transitions
		string str = trans.at(i);
		char *name = new char [str.length()+1];
		strcpy(name, str.c_str());
		sf->SetParName(7+i,name);
		sf->SetParameter(7+i,rI.at(i));
	}

	// Draw simulated hfs
	sf->SetNpx(1000);
	sf->Draw();

	/**********************************************************************
	*                             Tidy up                                 *
	**********************************************************************/
	// Clear vectors
	wNum.clear();
	rFreq.clear();
	nAlpha.clear();
	lPower.clear();
	pCurrent.clear();
	eNAlpha.clear();
	FL.clear();
	FU.clear();
	gC.clear();
	rI.clear();
}

/************************************************************************
*			Voigt func. for n peaks				*
************************************************************************/
double voigt(double *x, double *par)
{
	Float_t  xx  = x[0];    // Curent x value
	double CoG = par[0];  // Centre-of-Gravity
	double AU  = par[1];  // A-factor upper
	double AL  = par[2];  // A-factor lower
	double BU  = par[3];  // B-factor upper
	double BL  = par[4];  // B-factor lower
	double sig = par[5];	// Gaussian width
	double gam = par[6];	// Lorentz width

	// Declare function
	double f = 0.;
	
	// Declare variables that will be used for fitting Voigt function
	int _FL, _FU;
	int peak = 0;
	double _rI, mean;

	// Calculate where peaks should appear and their intensity
	// whilst minimising on CoG, A- and B-factor values
	for (int j=0; j<F_states(eSpinL); j++)
	{
		_FL = fabs(nSpin - eSpinL)+j;

		for (int k=0; k<F_states(eSpinU); k++)
		{
			_FU = fabs(nSpin - eSpinU)+k;

			if( ( _FL==0 && _FU==0 ) || ( abs(_FU-_FL)>1 ) )
			{
				f+=0;
			}
			else
			{
				mean = freq_shift(_FL,_FU,AL,BL,AU,BU)+CoG;
				_rI  = par[7+peak];
				f += _rI*TMath::Voigt(xx-mean,sig,gam,4);
				peak++;
			}
		}
	}

	// Return function
	return f;
}

/************************************************************************
*			How many transitions				*
************************************************************************/
int peaks(float nS, float eSL, float eSU)
{
	// Function for calculating how many transitions are possible
	// between the upper and lower states
	float F_lo_max = eSL+nS;
	float F_hi_max = eSU+nS;
	
	vector<float> F_lo, F_hi;
	int q = 0;
	while (fabs(nS-eSL)+q<=F_lo_max)
		{
			F_lo.push_back(fabs(nS-eSL)+q);
			q++;
		}
	q=0;
	while (fabs(nS-eSU)+q<=F_hi_max)
		{
			F_hi.push_back(fabs(nS-eSU)+q);
			q++;
		}
	
	int n = 0;
	for(int i=0; i<(int)F_lo.size(); i++)
		{
			if(F_lo.at(i)!=0)if(std::find(F_hi.begin(), F_hi.end(),F_lo.at(i)) != F_hi.end()) n++;
			if(std::find(F_hi.begin(), F_hi.end(),F_lo.at(i)-1) != F_hi.end()) n++;
			if(std::find(F_hi.begin(), F_hi.end(),F_lo.at(i)+1) != F_hi.end()) n++;
		}
	return n;
}
/************************************************************************
*			Calc. no. of F states				*
************************************************************************/
int F_states(float eS)
{
	// Function to calculate how many F-states there are
	// after hyperfine splitting of a level
	float F_hi = nSpin + eS;
	float F_lo = fabs(nSpin - eS);
	float nF = F_hi-F_lo+1;
	return (int)nF;
}

/************************************************************************
*			Calc. F state E level				*
************************************************************************/
double W_F(float F, float eS, float A, float B)
{
	// Funciton for calculating the shift in energy og an F state,
	// relative to the inital level
	float  I = nSpin;
	float  J = eS;
	double K = F*(F+1) - I*(I+1) - J*(J+1);

	double W_F;
	if(J==0)                    W_F = 0;
	else if((I==0.5)||(J==0.5)) W_F = 0.5*A*K;
	else                        W_F = (0.5*A*K) + (0.75*B*K*(K+1)
					  - I*(I+1)*J*(J+1))/(2*(2*I-1)*(2*J-1)*I*J);
	return W_F;
}

/************************************************************************
*			Calc. shift in centroid				*
************************************************************************/
double freq_shift(float FL, float FU, float AL, float BL, float AU, float BU)
{
	// Function to calculate the frequency shift of an F state, based on the
	// energy shift
	double   E_d = W_F(FL,eSpinL,AL,BL);
	double   E_u = W_F(FU,eSpinU,AU,BU);
	double delta = E_u - E_d;
	return delta;
}

/************************************************************************
*			Calc. relative peak h				*
************************************************************************/
double intensity(float FL, float FU)
{
	// Function for calculating the intensity of a transition between
	// initial and final F states
	double sixJ  = ROOT::Math::wigner_6j(2*FL, 2*FU, 2*1, 2*eSpinU, 2*eSpinL, 2*nSpin);
	double I_rel = ((2*FL)+1) * ((2*FU)+1) * sixJ * sixJ;

	return I_rel;
}
