// Some variables
double C = 299792458;		// Speed of light m/s

// Variables for reading in from file and plotting
double wN, nA, lP, pC;
vector<double> wNum, rFreq, nAlpha, lPower, pCurrent, eNAlpha;

// Vectors for the F values of the lower and upper electron states
vector<double> FL, FU, f_l, f_u;

// Vectors for estimated peak centroids and relative intensities
vector<double> gC, rI;

vector<int> valTrans;

//Double_t voigt(Double_t *x, Double_t *par);
double voigt(double *x, double *par);
double wnVoigt(double *x, double *par);
int peaks(float nS, float eSL, float eSU);
int F_states(float eS);
double W_F(float F, float eS, float A, float B);
double freq_shift(float FL, float FU, float AL, float BL, float AU, float BU);
double intensity(float FL, float FU);
